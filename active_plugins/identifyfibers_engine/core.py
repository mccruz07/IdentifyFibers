"""Algorithmic core for two-dimensional fiber extraction.

The implementation is independent of CellProfiler. It accepts a grayscale
image and, optionally, a label image defining parent objects. Every non-zero
label is analyzed independently, so fibers are clipped at object boundaries.

This is an original Python implementation inspired by the scientific problem
addressed by CT-FIRE. It is not a port and is not expected to reproduce
CT-FIRE fiber identities numerically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from .geometry import partition_regions, coverage_diagnostic, centerline_direction_tensor
from .path_identity import resolve_shared_tails, resolve_shared_segments
from skimage import exposure, filters, morphology, segmentation, measure


@dataclass(frozen=True)
class FiberSettings:
    """Settings controlling enhancement, segmentation and path filtering."""

    tracing_method: str = "skeleton"
    min_width: float = 1.0
    max_width: float = 8.0
    num_scales: int = 4
    threshold_method: str = "otsu"  # otsu, percentile, manual
    threshold_value: float = 0.15
    threshold_percentile: float = 75.0
    threshold_correction: float = 1.0
    threshold_lower: float = 0.0
    threshold_upper: float = 1.0
    minimum_length: float = 15.0
    minimum_area: int = 10
    closing_radius: int = 1
    reconstruct_crossings: bool = True
    crossing_search_radius: int = 4
    crossing_max_angle: float = 35.0
    tangent_points: int = 5
    pixel_size: float = 1.0
    polarity: str = "bright"  # bright or dark
    boundary_direction_weight: float = 1.0

    def validate(self) -> None:
        if not np.isfinite(self.boundary_direction_weight) or self.boundary_direction_weight < 0:
            raise ValueError("Object boundary direction guidance must be finite and nonnegative")
        if self.tracing_method not in {"skeleton", "distance_guided", "fire2d"}:
            raise ValueError("Unknown tracing method")
        if self.min_width <= 0 or self.max_width < self.min_width:
            raise ValueError("Fiber widths must satisfy 0 < min_width <= max_width")
        if self.num_scales < 1:
            raise ValueError("num_scales must be at least 1")
        if self.threshold_method not in {"otsu", "li", "percentile", "manual"}:
            raise ValueError("Unknown threshold method")
        if self.threshold_correction <= 0 or not 0 <= self.threshold_lower <= self.threshold_upper <= 1:
            raise ValueError("Threshold correction must be positive and bounds must satisfy 0 <= lower <= upper <= 1")
        if not 0 <= self.threshold_value <= 1:
            raise ValueError("threshold_value must be in [0, 1]")
        if not 0 <= self.threshold_percentile <= 100:
            raise ValueError("threshold_percentile must be in [0, 100]")
        if self.minimum_length < 0 or self.pixel_size <= 0:
            raise ValueError("minimum_length must be nonnegative and pixel_size positive")
        if self.polarity not in {"bright", "dark"}:
            raise ValueError("polarity must be bright or dark")
        if self.crossing_search_radius < 1 or self.tangent_points < 2:
            raise ValueError("crossing_search_radius must be >= 1 and tangent_points >= 2")
        if not 0 <= self.crossing_max_angle <= 90:
            raise ValueError("crossing_max_angle must be in [0, 90]")


@dataclass
class Fiber:
    """One extracted centerline and its measurements."""

    number: int
    parent_object: int
    coordinates: np.ndarray
    length: float
    centerline_length: float
    end_to_end_length: float
    straightness: float
    tortuosity: float
    orientation: float
    width_mean: float
    width_median: float
    width_maximum: float
    curvature_mean: float
    curvature_maximum: float
    intensity_mean: float
    intensity_maximum: float
    intensity_integrated: float
    center_y: float
    center_x: float
    start_y: float
    start_x: float
    end_y: float
    end_x: float
    nearest_distance: float = float("nan")
    neighbor_alignment: float = float("nan")
    isolated: int = 0


def measure_neighbors(fibers, radius):
    """Centroid distance and axial alignment within radius, per parent region."""
    if radius <= 0:
        raise ValueError("Neighborhood radius must be positive")
    for parent in {fiber.parent_object for fiber in fibers}:
        group = [fiber for fiber in fibers if fiber.parent_object == parent]
        if not group:
            continue
        xy = np.array([(f.center_y, f.center_x) for f in group], dtype=float)
        tree = cKDTree(xy)
        distances, _ = tree.query(xy, k=2) if len(group) > 1 else (None, None)
        for index, fiber in enumerate(group):
            fiber.nearest_distance = float(distances[index, 1]) if distances is not None else float("nan")
            neighbors = [j for j in tree.query_ball_point(xy[index], radius) if j != index]
            fiber.isolated = int(not neighbors)
            if neighbors and np.isfinite(fiber.orientation):
                angles = np.array([group[j].orientation for j in neighbors], dtype=float)
                angles = angles[np.isfinite(angles)]
                if angles.size:
                    fiber.neighbor_alignment = float(np.mean((1 + np.cos(np.deg2rad(2 * (angles - fiber.orientation)))) / 2))


@dataclass
class FiberAnalysisResult:
    fibers: list[Fiber] = field(default_factory=list)
    enhanced: np.ndarray | None = None
    binary_mask: np.ndarray | None = None
    skeleton: np.ndarray | None = None
    label_image: np.ndarray | None = None
    rejected_short: np.ndarray | None = None
    uncovered_skeleton: np.ndarray | None = None
    coverage_percent: float = float("nan")
    reference_skeleton_pixels: int = 0
    resolved_shared_tails: int = 0
    resolved_shared_segments: int = 0


_NEIGHBORS = tuple(
    (dy, dx)
    for dy in (-1, 0, 1)
    for dx in (-1, 0, 1)
    if not (dy == 0 and dx == 0)
)


def centerline_labels(fibers, shape):
    """Label accepted centerline pixels for visual quality control.

    At shared pixels the last fiber wins; such crossings are therefore not
    represented as distinct overlapping labels. This never changes objects.
    """
    labels = np.zeros(shape, dtype=np.int32)
    for fiber in fibers:
        coords = np.asarray(fiber.coordinates, dtype=int)
        if coords.size:
            labels[coords[:, 0], coords[:, 1]] = fiber.number
    return labels


def _normalize(image: np.ndarray, mask: np.ndarray) -> np.ndarray:
    values = image[mask & np.isfinite(image)]
    if values.size == 0:
        return np.zeros(image.shape, dtype=np.float32)
    low, high = np.percentile(values, (1.0, 99.5))
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return np.zeros(image.shape, dtype=np.float32)
    result = exposure.rescale_intensity(image, in_range=(low, high), out_range=(0.0, 1.0))
    result = np.asarray(result, dtype=np.float32)
    result[~mask] = 0
    result[~np.isfinite(result)] = 0
    return result


def _enhance(image: np.ndarray, mask: np.ndarray, settings: FiberSettings) -> np.ndarray:
    normalized = _normalize(image, mask)
    if settings.polarity == "dark":
        normalized = np.where(mask, 1.0 - normalized, 0.0)
    sigmas = np.linspace(
        max(0.5, settings.min_width / 2.0),
        max(0.5, settings.max_width / 2.0),
        settings.num_scales,
    )
    response = filters.sato(
        normalized,
        sigmas=sigmas,
        black_ridges=False,
        mode="constant",
    )
    response = _normalize(np.asarray(response, dtype=np.float32), mask)
    response[~mask] = 0
    return response


def _segment(enhanced: np.ndarray, region: np.ndarray, settings: FiberSettings) -> np.ndarray:
    values = enhanced[region]
    values = values[np.isfinite(values)]
    if values.size == 0 or np.max(values) <= 0:
        return np.zeros(region.shape, dtype=bool)
    if settings.threshold_method == "otsu":
        positive = values[values > 0]
        threshold = filters.threshold_otsu(positive) if positive.size else 1.0
    elif settings.threshold_method == "li":
        positive = values[values > 0]
        threshold = filters.threshold_li(positive) if positive.size and np.ptp(positive) > 0 else (float(positive[0]) if positive.size else 1.0)
    elif settings.threshold_method == "percentile":
        threshold = np.percentile(values, settings.threshold_percentile)
    else:
        threshold = settings.threshold_value
    threshold = float(np.clip(threshold * settings.threshold_correction, settings.threshold_lower, settings.threshold_upper))
    binary = region & (enhanced >= threshold) & (enhanced > 0)
    if settings.closing_radius > 0:
        binary = morphology.closing(binary, morphology.disk(settings.closing_radius))
    # Remove small components without relying on a scikit-image argument whose
    # name and boundary semantics changed after the version in CellProfiler.
    binary &= region
    component_labels, _ = ndi.label(binary)
    component_sizes = np.bincount(component_labels.ravel())
    remove = component_sizes < settings.minimum_area
    remove[0] = False
    binary = binary & ~remove[component_labels]
    return np.asarray(binary, dtype=bool)


def _neighbors(point: tuple[int, int], pixels: set[tuple[int, int]]) -> list[tuple[int, int]]:
    y, x = point
    # Do not add a diagonal shortcut when an orthogonal connection exists.
    # Such triangles turn ordinary pixel-scale bends into false junctions.
    return sorted((y + dy, x + dx) for dy, dx in _NEIGHBORS
                  if (y + dy, x + dx) in pixels
                  and not (dy and dx and ((y + dy, x) in pixels or (y, x + dx) in pixels)))


def _edge_key(a: tuple[int, int], b: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return (a, b) if a <= b else (b, a)


def _trace_paths(skeleton: np.ndarray) -> list[np.ndarray]:
    """Split a skeleton into ordered paths between endpoints and junctions."""

    pixels = {tuple(p) for p in np.argwhere(skeleton)}
    if not pixels:
        return []
    adjacency = {p: _neighbors(p, pixels) for p in pixels}
    nodes = {p for p, ns in adjacency.items() if len(ns) != 2}
    visited: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    paths: list[np.ndarray] = []

    for start in sorted(nodes):
        for nxt in adjacency[start]:
            edge = _edge_key(start, nxt)
            if edge in visited:
                continue
            visited.add(edge)
            path = [start, nxt]
            previous, current = start, nxt
            while current not in nodes:
                candidates = [p for p in adjacency[current] if p != previous]
                if not candidates:
                    break
                following = candidates[0]
                visited.add(_edge_key(current, following))
                path.append(following)
                previous, current = current, following
            paths.append(np.asarray(path, dtype=np.int32))

    # Components that contain no endpoints or junctions are closed loops.
    for start in sorted(pixels):
        remaining = [n for n in adjacency[start] if _edge_key(start, n) not in visited]
        if not remaining:
            continue
        path = [start]
        previous = None
        current = start
        while True:
            candidates = [
                n for n in adjacency[current]
                if n != previous and _edge_key(current, n) not in visited
            ]
            if not candidates:
                break
            following = candidates[0]
            visited.add(_edge_key(current, following))
            path.append(following)
            previous, current = current, following
            if current == start:
                break
        if len(path) > 1:
            paths.append(np.asarray(path, dtype=np.int32))
    return paths


def _endpoint_tangent(path: np.ndarray, endpoint: int, points: int) -> np.ndarray | None:
    """Return the unit vector pointing from a junction into a branch."""

    if len(path) < 2:
        return None
    step = min(points - 1, len(path) - 1)
    if endpoint == 0:
        vector = path[step].astype(float) - path[0].astype(float)
    else:
        vector = path[-step - 1].astype(float) - path[-1].astype(float)
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 0 else None


def _junction_regions(skeleton: np.ndarray, radius: int) -> np.ndarray:
    pixels = {tuple(p) for p in np.argwhere(skeleton)}
    junctions = np.zeros_like(skeleton, dtype=bool)
    for point in pixels:
        if len(_neighbors(point, pixels)) > 2:
            junctions[point] = True
    if not np.any(junctions):
        return np.zeros(skeleton.shape, dtype=np.int32)
    expanded = morphology.dilation(junctions, morphology.disk(radius))
    return ndi.label(expanded)[0].astype(np.int32)


def _skeleton_connector(a, b, skeleton, radius):
    """Shortest local route in skeleton steps, with no interpolation across background."""
    from collections import deque
    start, target = tuple(map(int, a)), tuple(map(int, b))
    if start == target:
        return np.asarray([start], dtype=np.int32)
    margin = radius + 1
    lo = np.maximum(np.minimum(a, b) - margin, 0).astype(int)
    hi = np.minimum(np.maximum(a, b) + margin + 1, skeleton.shape).astype(int)
    pixels = {tuple(p + lo) for p in np.argwhere(skeleton[lo[0]:hi[0], lo[1]:hi[1]])}
    queue, previous = deque([start]), {start: None}
    while queue:
        point = queue.popleft()
        if point == target:
            route = []
            while point is not None:
                route.append(point)
                point = previous[point]
            route = np.asarray(route[::-1], dtype=np.int32)
            return route if _path_length(route, 1) <= 2 * radius + 2 else None
        for neighbor in _neighbors(point, pixels):
            if neighbor not in previous:
                previous[neighbor] = point
                queue.append(neighbor)
    return None


def _merge_paths_at_crossings(
    paths: list[np.ndarray], skeleton: np.ndarray, settings: FiberSettings
) -> list[np.ndarray]:
    """Reconnect branches through junctions using local directional continuity.

    Each branch endpoint can be paired only once. Candidate pairs must terminate
    in the same dilated junction region and their outward tangent vectors must
    form a nearly straight continuation. Connectors must follow a short local
    skeleton route without doubling back through branch interiors.
    Greedy selection is deterministic and
    begins with the smallest angular deflection.
    """

    if not settings.reconstruct_crossings or len(paths) < 2:
        return paths
    regions = _junction_regions(skeleton, settings.crossing_search_radius)
    if regions.max() == 0:
        return paths

    endpoint_region: dict[tuple[int, int], int] = {}
    tangents: dict[tuple[int, int], np.ndarray] = {}
    for index, path in enumerate(paths):
        if len(path) < 2:
            continue
        for endpoint, coordinate in ((0, path[0]), (1, path[-1])):
            y, x = map(int, coordinate)
            region = int(regions[y, x])
            tangent = _endpoint_tangent(path, endpoint, settings.tangent_points)
            if region > 0 and tangent is not None:
                token = (index, endpoint)
                endpoint_region[token] = region
                tangents[token] = tangent

    candidates: list[tuple[float, float, tuple[int, int], tuple[int, int]]] = []
    tokens = sorted(endpoint_region)
    for position, first in enumerate(tokens):
        for second in tokens[position + 1 :]:
            if first[0] == second[0] or endpoint_region[first] != endpoint_region[second]:
                continue
            p1 = paths[first[0]][0 if first[1] == 0 else -1]
            p2 = paths[second[0]][0 if second[1] == 0 else -1]
            distance = float(np.linalg.norm(p1.astype(float) - p2.astype(float)))
            if distance > 2 * settings.crossing_search_radius + 1:
                continue
            dot = float(np.clip(np.dot(tangents[first], tangents[second]), -1.0, 1.0))
            deflection = float(np.degrees(np.arccos(-dot)))
            if deflection <= settings.crossing_max_angle:
                candidates.append((deflection, distance, first, second))

    links: dict[tuple[int, int], tuple[int, int]] = {}
    connectors = {}
    for _, _, first, second in sorted(candidates):
        if first not in links and second not in links:
            p1 = paths[first[0]][0 if first[1] == 0 else -1]
            p2 = paths[second[0]][0 if second[1] == 0 else -1]
            connector = _skeleton_connector(p1, p2, skeleton, settings.crossing_search_radius)
            if connector is None:
                continue
            # A connection must not double back through either branch interior.
            interior = set(map(tuple, connector[1:-1]))
            branch_interior = set(map(tuple, paths[first[0]][1:-1])) | set(map(tuple, paths[second[0]][1:-1]))
            if interior & branch_interior:
                continue
            links[first] = second
            links[second] = first
            connectors[first] = connector
            connectors[second] = connector[::-1]

    merged: list[np.ndarray] = []
    visited: set[int] = set()
    # Start with chain ends; otherwise starting in the middle consumes only
    # one direction and incorrectly fragments a continuous reconstructed fiber.
    order = sorted(range(len(paths)), key=lambda i: (
        0 if any((i, end) not in links for end in (0, 1)) else 1, i
    ))
    for initial in order:
        if initial in visited:
            continue
        unlinked = [end for end in (0, 1) if (initial, end) not in links]
        entry = unlinked[0] if unlinked else 0
        current, current_entry = initial, entry
        assembled: np.ndarray | None = None
        while current not in visited:
            visited.add(current)
            oriented = paths[current] if current_entry == 0 else paths[current][::-1]
            if assembled is None:
                assembled = oriented.copy()
            else:
                connector = connectors[(current, current_entry)][::-1]
                assembled = np.vstack((assembled, connector[1:], oriented[1:]))
            exit_token = (current, 1 - current_entry)
            next_token = links.get(exit_token)
            if next_token is None or next_token[0] in visited:
                break
            current, current_entry = next_token
        if assembled is not None:
            keep = np.ones(len(assembled), dtype=bool)
            keep[1:] = np.any(assembled[1:] != assembled[:-1], axis=1)
            merged.append(assembled[keep])
    return merged


def _trace_distance_guided(skeleton, distance, settings):
    """Experimental 2D distance-guided graph tracker, inspired by FIRE.

    Traverses skeleton edges from endpoints then remaining junctions; direction
    continuity and distance-to-background guide local branch selection. This is
    not a translation of FIRE and has no curvelet preprocessing.
    """
    pixels = {tuple(p) for p in np.argwhere(skeleton)}
    adjacency = {p: _neighbors(p, pixels) for p in pixels}
    edges = {_edge_key(p, q) for p in pixels for q in adjacency[p]}
    remaining = set(edges)
    paths = []
    # Endpoints make useful starting positions; then consume unused cycles or branches.
    order = sorted(pixels, key=lambda p: (len(adjacency[p]) != 1, -distance[p], p))
    for start in order:
        while any(_edge_key(start, n) in remaining for n in adjacency[start]):
            path = [start]
            while True:
                current = path[-1]
                candidates = [n for n in adjacency[current] if _edge_key(current, n) in remaining]
                if not candidates:
                    break
                if len(path) > 1:
                    origin = np.asarray(path[max(0, len(path)-settings.tangent_points)], float)
                    forward = np.asarray(current, float) - origin
                    norm = np.linalg.norm(forward)
                    scored = []
                    for n in candidates:
                        direction = np.asarray(n, float) - current
                        cosine = float(np.dot(forward, direction) / (norm * np.linalg.norm(direction))) if norm else 1.0
                        angle = float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
                        if len(adjacency[current]) <= 2 or angle <= settings.crossing_max_angle:
                            medial = distance[n] / (distance[current] + distance[n] + 1e-12)
                            scored.append((-(cosine + 0.15 * medial), n))
                    if not scored:
                        break
                    nxt = min(scored)[1]
                else:
                    nxt = min(candidates, key=lambda n: (-distance[n], n))
                remaining.remove(_edge_key(current, nxt))
                path.append(nxt)
                if nxt == start:
                    break
            if len(path) >= 2:
                paths.append(np.asarray(path, dtype=np.int32))
    return paths


def _path_length(coords: np.ndarray, pixel_size: float) -> float:
    if len(coords) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(coords.astype(float), axis=0), axis=1).sum() * pixel_size)


def _curvature(coords: np.ndarray, pixel_size: float) -> tuple[float, float]:
    if len(coords) < 3:
        return 0.0, 0.0
    points = coords.astype(float) * pixel_size
    v1 = points[1:-1] - points[:-2]
    v2 = points[2:] - points[1:-1]
    n1 = np.linalg.norm(v1, axis=1)
    n2 = np.linalg.norm(v2, axis=1)
    valid = (n1 > 0) & (n2 > 0)
    if not np.any(valid):
        return 0.0, 0.0
    cosine = np.sum(v1[valid] * v2[valid], axis=1) / (n1[valid] * n2[valid])
    angle = np.arccos(np.clip(cosine, -1.0, 1.0))
    arc = (n1[valid] + n2[valid]) / 2.0
    values = angle / np.maximum(arc, np.finfo(float).eps)
    return float(np.mean(values)), float(np.max(values))


def _orientation(coords: np.ndarray) -> float:
    if len(coords) < 2:
        return float("nan")
    xy = coords[:, ::-1].astype(float)
    centered = xy - xy.mean(axis=0)
    covariance = centered.T @ centered
    values, vectors = np.linalg.eigh(covariance)
    vx, vy = vectors[:, np.argmax(values)]
    return float(np.degrees(np.arctan2(vy, vx)) % 180.0)


def _measure_path(
    number: int,
    parent: int,
    coords: np.ndarray,
    image: np.ndarray,
    distance: np.ndarray,
    settings: FiberSettings,
) -> Fiber:
    y, x = coords[:, 0], coords[:, 1]
    length = _path_length(coords, settings.pixel_size)
    end_to_end = float(np.linalg.norm((coords[-1] - coords[0]).astype(float)) * settings.pixel_size)
    straightness = end_to_end / length if length > 0 else float("nan")
    widths = 2.0 * distance[y, x] * settings.pixel_size
    intensities = np.asarray(image[y, x], dtype=float)
    curvature_mean, curvature_maximum = _curvature(coords, settings.pixel_size)
    center = coords.astype(float).mean(axis=0)
    return Fiber(
        number=number,
        parent_object=parent,
        coordinates=coords,
        length=length,
        centerline_length=length,
        end_to_end_length=end_to_end,
        straightness=straightness,
        tortuosity=(length / end_to_end if end_to_end > 0 else float("nan")),
        orientation=_orientation(coords),
        width_mean=float(np.mean(widths)),
        width_median=float(np.median(widths)),
        width_maximum=float(np.max(widths)),
        curvature_mean=curvature_mean,
        curvature_maximum=curvature_maximum,
        intensity_mean=float(np.mean(intensities)),
        intensity_maximum=float(np.max(intensities)),
        intensity_integrated=float(np.sum(intensities)),
        center_y=float(center[0]), center_x=float(center[1]),
        start_y=float(coords[0, 0]), start_x=float(coords[0, 1]),
        end_y=float(coords[-1, 0]), end_x=float(coords[-1, 1]),
    )


def analyze_fibers(
    image: np.ndarray,
    labels: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    settings: Optional[FiberSettings] = None,
    segmentation_mask: Optional[np.ndarray] = None,
) -> FiberAnalysisResult:
    """Identify and measure fibers in a whole image or within labeled objects.

    When ``labels`` is supplied, every positive label is processed separately.
    Fibers crossing a parent boundary are clipped at that boundary. If labels
    is omitted, the valid mask is treated as one parent region with identifier 0.
    """

    settings = settings or FiberSettings()
    settings.validate()
    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError("IdentifyFibers currently supports two-dimensional grayscale images")
    valid = np.isfinite(image) if mask is None else (np.asarray(mask, dtype=bool) & np.isfinite(image))
    if valid.shape != image.shape:
        raise ValueError("mask and image must have the same shape")
    if segmentation_mask is not None:
        segmentation_mask = np.asarray(segmentation_mask, dtype=bool)
        if segmentation_mask.shape != image.shape:
            raise ValueError("segmentation mask and image must have the same shape")
    if labels is not None:
        labels = np.asarray(labels)
        if labels.shape != image.shape or labels.ndim != 2:
            raise ValueError("labels and image must be two-dimensional and have the same shape")
        parent_ids: Iterable[int] = [int(x) for x in np.unique(labels) if x > 0]
    else:
        parent_ids = [0]

    combined_enhanced = np.zeros(image.shape, dtype=np.float32)
    combined_binary = np.zeros(image.shape, dtype=bool)
    combined_skeleton = np.zeros(image.shape, dtype=bool)
    output_labels = np.zeros(image.shape, dtype=np.int32)
    rejected_short = np.zeros(image.shape, dtype=bool)
    uncovered_skeleton = np.zeros(image.shape, dtype=bool)
    reference_count = 0
    resolved_shared_tails = 0
    resolved_shared_segments = 0
    fibers: list[Fiber] = []

    for parent in parent_ids:
        region = valid if labels is None else (valid & (labels == parent))
        if not np.any(region):
            continue
        enhanced = _enhance(image, region, settings)
        if segmentation_mask is None:
            binary = _segment(enhanced, region, settings)
        else:
            binary = segmentation_mask & region
            if settings.closing_radius:
                binary = morphology.closing(binary, morphology.disk(settings.closing_radius)) & region
            components, _ = ndi.label(binary)
            sizes = np.bincount(components.ravel())
            binary &= (components > 0) & (sizes[components] >= settings.minimum_area)
        skeleton = (np.zeros_like(binary) if settings.tracing_method == "fire2d"
                    else morphology.skeletonize(binary) & region)
        distance = ndi.distance_transform_edt(binary)
        combined_enhanced = np.maximum(combined_enhanced, enhanced)
        combined_binary |= binary
        if settings.tracing_method != "fire2d":
            combined_skeleton |= skeleton
        markers = np.zeros(image.shape, dtype=np.int32)
        first_fiber_index = len(fibers)
        if settings.tracing_method == "fire2d":
            from .fire2d import trace_network
            paths, _, seeds = trace_network(binary, settings.crossing_max_angle, settings.tangent_points)
        elif settings.tracing_method == "distance_guided":
            paths = _trace_distance_guided(skeleton, distance, settings)
        else:
            paths = _merge_paths_at_crossings(_trace_paths(skeleton), skeleton, settings)
        paths, resolved = resolve_shared_tails(paths, distance, binary, settings.minimum_length/settings.pixel_size)
        resolved_shared_tails += resolved
        paths, resolved = resolve_shared_segments(paths, distance, binary, settings.minimum_length/settings.pixel_size, image=image)
        resolved_shared_segments += resolved
        if settings.tracing_method == "fire2d":
            for path in paths:
                combined_skeleton[path[:, 0], path[:, 1]] = True
        for coords in paths:
            if _path_length(coords, settings.pixel_size) < settings.minimum_length:
                rejected_short[coords[:, 0], coords[:, 1]] = True
                continue
            # Marker pixels must remain inside the segmented parent region.
            inside = binary[coords[:, 0], coords[:, 1]] & region[coords[:, 0], coords[:, 1]]
            points = coords[inside]
            available = markers[points[:, 0], points[:, 1]] == 0
            if not np.any(available):
                continue
            number = len(fibers) + 1
            fiber = _measure_path(number, parent, coords, image, distance, settings)
            fibers.append(fiber)
            markers[points[available, 0], points[available, 1]] = number

        reference = morphology.skeletonize(binary) & region
        _, uncovered = coverage_diagnostic(reference, markers > 0)
        uncovered_skeleton |= uncovered
        reference_count += int(reference.sum())

        # Restrict propagation to components containing an accepted fiber. This
        # prevents discarded isolated fragments from becoming fiber objects.
        if len(fibers) > first_fiber_index:
            components, _ = ndi.label(binary, structure=np.ones((3, 3)))
            seeded = np.unique(components[markers > 0])
            accepted_mask = binary & np.isin(components, seeded[seeded > 0]) & region
            directions = (centerline_direction_tensor(fibers[first_fiber_index:], image.shape)
                          if settings.boundary_direction_weight else None)
            grown = partition_regions(accepted_mask, markers, image, directions, settings.boundary_direction_weight)
            output_labels[region] = grown[region]
            sums = np.bincount(grown.ravel(), weights=np.where(region, image, 0).ravel())
            for prop in measure.regionprops(grown, intensity_image=image):
                fiber = fibers[prop.label - 1]
                fiber.length = float(prop.major_axis_length * settings.pixel_size)
                fiber.intensity_mean = float(prop.mean_intensity)
                fiber.intensity_maximum = float(prop.max_intensity)
                fiber.intensity_integrated = float(sums[prop.label])
                fiber.center_y, fiber.center_x = map(float, prop.centroid)

    return FiberAnalysisResult(
        fibers=fibers,
        resolved_shared_tails=resolved_shared_tails,
        resolved_shared_segments=resolved_shared_segments,
        enhanced=combined_enhanced,
        binary_mask=combined_binary,
        skeleton=combined_skeleton,
        label_image=output_labels,
        rejected_short=rejected_short & (output_labels == 0),
        uncovered_skeleton=uncovered_skeleton,
        reference_skeleton_pixels=reference_count,
        coverage_percent=(100.*(reference_count-int(uncovered_skeleton.sum()))/reference_count
                          if reference_count else float("nan")),
    )

