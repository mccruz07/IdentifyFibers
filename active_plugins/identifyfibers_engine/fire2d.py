"""Independent Sato/FIRE-inspired 2D network extraction.

Nucleation and extension principles follow Stein's FIRE (see LICENSE_FIRE.txt).
This is an independent implementation, not a numerically equivalent FIRE port.
Only NumPy, SciPy and scikit-image are required. No skeleton is used to seed or
extend the network. Widths are measured later on the segmentation mask.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.draw import line
from scipy.signal import find_peaks
from functools import lru_cache


def _segment(a, b, mask):
    rr, cc = line(*map(int, a), *map(int, b))
    if np.all(mask[rr, cc]):
        return np.column_stack((rr, cc))
    return None


@lru_cache(maxsize=64)
def _ring_offsets(radius):
    r = int(radius)
    return np.array([(-r, dx) for dx in range(-r, r)] +
                    [(dy, r) for dy in range(-r, r)] +
                    [(r, dx) for dx in range(r, -r, -1)] +
                    [(dy, -r) for dy in range(r, -r, -1)], dtype=int)


def _ring_candidates(point, radius, distance, mask, all_foreground=False):
    """True local maxima on a closed square perimeter, including flat peak tops."""
    coords = _ring_offsets(radius) + point
    valid = ((coords >= 0) & (coords < mask.shape)).all(axis=1)
    scores = np.zeros(len(coords))
    indices = np.flatnonzero(valid)
    scores[indices] = distance[coords[indices, 0], coords[indices, 1]]
    if all_foreground:
        indices = np.flatnonzero(valid & (scores > .1))
    else:
        # Tiling makes plateau handling circular. Compiled peak detection also
        # rejects flat minima, unlike >= comparisons of neighboring pixels.
        n = len(scores)
        _, props = find_peaks(np.tile(np.round(scores, 9), 3), plateau_size=(1, None), prominence=max(.02, .1 * scores.max()))
        middle = (props['left_edges'] + props['right_edges'] + 1) // 2
        indices = middle[(middle >= n) & (middle < 2*n)] - n
        indices = indices[scores[indices] >= (.5 * scores.max() if scores.max() > 4 else .1)]
        indices = sorted(indices, key=lambda i: (-scores[i], i))
    return [coords[i] for i in indices
            if valid[i] and mask[tuple(coords[i])] and _segment(point, coords[i], mask) is not None]


def nucleate(mask, distance, spacing=4):
    """Distance maxima with deterministic suppression and per-component fallback."""
    if not mask.any():
        return np.empty((0,2),int)
    # Coordinate-only perturbation breaks flat ridges reproducibly, no global RNG.
    y,x=np.indices(mask.shape)
    perturbed=distance + 1e-7*((y*73856093 ^ x*19349663) % 997)/997
    maxima=mask & (perturbed == ndi.maximum_filter(perturbed, size=2*spacing+1, mode='constant'))
    seeds=list(map(tuple,np.argwhere(maxima)))
    labels,count=ndi.label(mask, structure=np.ones((3,3)))
    present=set(labels[maxima])
    for index in range(1,count+1):
        if index not in present:
            region=np.argwhere(labels==index)
            p=region[np.argmax(distance[region[:,0],region[:,1]])]
            seeds.append(tuple(p))
    return np.array(sorted(seeds,key=lambda p:(-distance[p],p)),int)


def _tangents(path, window=5):
    indices = np.arange(len(path))
    delta = path[np.minimum(indices + window, len(path) - 1)] - path[np.maximum(indices - window, 0)]
    norm = np.linalg.norm(delta, axis=1)
    return delta / np.maximum(norm[:, None], 1e-12)


def _deduplicate(paths, distance_map=None, mask=None):
    """Suppress re-traces within the local fiber width, checking orientation.

    A fixed pixel tolerance missed opposite sides of the same broad fiber.
    Width-based matching still requires a foreground-supported connector, so
    separated parallel fibers are not collapsed across a background gap.
    """
    retained = []
    trees = []
    directions = []
    boxes = []
    margin = 1.5 if distance_map is None else 2 * float(distance_map.max())
    for path in sorted(paths, key=lambda p: -len(p)):
        tangent = _tangents(path)
        low, high = path.min(axis=0), path.max(axis=0)
        duplicate = False
        for tree, other_direction, (a, b) in zip(trees, directions, boxes):
            if np.any(high + margin < a) or np.any(b + margin < low):
                continue
            distance, index = tree.query(path)
            same_direction = np.abs(np.sum(tangent * other_direction[index], axis=1)) >= np.cos(np.deg2rad(35))
            tolerance = 1.5
            if distance_map is not None:
                nearest = tree.data[index].astype(int)
                tolerance = np.maximum(1.5, .8 * (distance_map[path[:, 0], path[:, 1]] +
                                                   distance_map[nearest[:, 0], nearest[:, 1]]))
            matching = (distance <= tolerance) & same_direction
            if mask is not None and matching.any():
                nearest = tree.data[index].astype(int)
                for j in np.flatnonzero(matching & (distance > 1)):
                    if _segment(path[j], nearest[j], mask) is None:
                        matching[j] = False
            if np.mean(matching) >= .8:
                duplicate = True
                break
        if duplicate:
            continue
        retained.append(path)
        trees.append(cKDTree(path))
        directions.append(tangent)
        boxes.append((low, high))
    return retained


def _join(paths, mask, max_angle, tangent_points, distance_map=None):
    """Join compatible terminal portions, trimming their overlapping tails.

    Nucleation from opposite sides can overshoot the same crossing. Comparing
    only the final pixels left these trajectories split even when they already
    overlapped. Each end may be paired once; branch interiors stay untouched.
    """
    if len(paths) < 2:
        return paths
    points, tokens = [], []
    for i, path in enumerate(paths):
        arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
        for end in (0, 1):
            radius = 2.0 if distance_map is None else float(distance_map[tuple(path[0 if end == 0 else -1])])
            reach = max(3., min(12., 2 * radius + tangent_points))
            terminal = arc if end == 0 else arc[-1] - arc
            for index in np.flatnonzero((terminal <= reach) & (terminal < .4 * arc[-1])):
                span = tangent_points
                if distance_map is not None:
                    span = max(span, int(np.ceil(3 * distance_map[tuple(path[index])])))
                other = min(len(path)-1, index + span) if end == 0 else max(0, index-span)
                delta = (path[other] - path[index]).astype(float)
                norm = np.linalg.norm(delta)
                if norm == 0:
                    continue
                points.append(path[index])
                tokens.append((i, end, int(index), float(terminal[index]), delta/norm))
    if len(points) < 2:
        return paths
    options = []
    for a, b in sorted(cKDTree(points).query_pairs(1.5)):
        i, e, ia, ta, u = tokens[a]
        j, f, ib, tb, v = tokens[b]
        if i == j:
            continue
        angle = float(np.degrees(np.arccos(np.clip(-np.dot(u, v), -1, 1))))
        if angle > max_angle:
            continue
        bridge = _segment(points[a], points[b], mask)
        if bridge is None:
            continue
        options.append((angle + ta + tb, angle, ta + tb, i, e, j, f, ia, ib, bridge))
    links, cuts, bridges = {}, {}, {}
    for _, _, _, i, e, j, f, ia, ib, bridge in sorted(options, key=lambda q: q[:-1]):
        a, b = (i, e), (j, f)
        if a in links or b in links:
            continue
        # Two joins must leave a nonempty ordered body between the cuts.
        if (e == 0 and ia >= cuts.get((i, 1), len(paths[i])-1) or
            e == 1 and ia <= cuts.get((i, 0), 0) or
            f == 0 and ib >= cuts.get((j, 1), len(paths[j])-1) or
            f == 1 and ib <= cuts.get((j, 0), 0)):
            continue
        links[a], links[b] = b, a
        cuts[a], cuts[b] = ia, ib
        bridges[a], bridges[b] = bridge, bridge[::-1]
    visited, result = set(), []
    order = sorted(range(len(paths)), key=lambda i: ((i, 0) in links and (i, 1) in links, i))
    for first in order:
        if first in visited:
            continue
        entry = 0 if (first, 0) not in links else 1
        current, chunks = first, []
        while current not in visited:
            visited.add(current)
            path = paths[current][cuts.get((current, 0), 0):cuts.get((current, 1), len(paths[current])-1)+1]
            chunks.append(path if entry == 0 else path[::-1])
            out = (current, 1-entry)
            if out not in links:
                break
            dest = links[out]
            if dest[0] in visited:
                break
            chunks.append(bridges[out])
            current, entry = dest
        merged = np.vstack(chunks)
        keep = np.r_[True, np.any(np.diff(merged, axis=0) != 0, axis=1)]
        result.append(merged[keep])
    return result


def _extend(seed, first, distance, mask, max_angle):
    """Continue between distance ridges on a scale wider than the local fiber.

    The perimeter must reach beyond the fiber radius; searching inside a wide
    fiber left the old tracer free to follow an arbitrary parallel offset.
    Smaller radii are retried at curves and tips, never across background.
    """
    direction = (first - seed).astype(float)
    direction /= np.linalg.norm(direction)
    initial = _segment(seed, first, mask)
    chunks = [initial]
    seen = set(map(tuple, initial))
    current = first
    for _ in range(mask.size):
        local_radius = float(distance[tuple(current)])
        if local_radius >= 2:
            reach = int(np.ceil(1.6 * local_radius)) + 1
            end = np.rint(current + reach * direction).astype(int)
            rr, cc = line(*current, *end)
            samples = np.column_stack((rr, cc))
            valid = ((samples >= 0) & (samples < mask.shape)).all(axis=1)
            supported = np.zeros(len(samples), bool)
            supported[valid] = mask[rr[valid], cc[valid]]
            outside = np.flatnonzero(~supported)
            if len(outside) and outside[0] > 1:
                ray = samples[:outside[0]]
                widths = distance[ray[:, 0], ray[:, 1]]
                if (np.all(np.diff(widths) <= .15) and widths[-1] < .65 * local_radius
                        and not any(tuple(q) in seen for q in ray[1:])):
                    chunks.append(ray[1:])
                    break
        target = int(np.clip(np.ceil(2 * distance[tuple(current)]), 2, 16))
        candidates = []
        for radius in range(target, 0, -1):
            peaks = _ring_candidates(current, radius, distance, mask,
                                     radius == 1 and distance[tuple(current)] < 1.5)
            for nxt in peaks:
                unit = (nxt - current).astype(float)
                unit /= np.linalg.norm(unit)
                cosine = float(np.dot(direction, unit))
                limit = max(max_angle, 70 if radius == 1 else np.degrees(np.arctan2(1.0, radius)))
                if cosine < np.cos(np.deg2rad(min(limit, 85))):
                    continue
                segment = _segment(current, nxt, mask)
                if any(tuple(q) in seen for q in segment[1:]):
                    continue
                candidates.append((-cosine, -distance[tuple(nxt)], tuple(nxt), unit, segment))
            if candidates:
                break
        if not candidates:
            break
        _, _, coord, unit, segment = min(candidates, key=lambda q: q[:3])
        chunks.append(segment[1:])
        seen.update(map(tuple, segment[1:]))
        current = np.array(coord)
        direction = direction + unit
        direction /= np.linalg.norm(direction)
    return np.vstack(chunks)


def trace_network(mask, max_angle=35, tangent_points=5, seed_spacing=4):
    """Euclidean distance seeds -> adaptive ridge tracing -> terminal stitching.

    Width-aware duplicate removal precedes and follows terminal stitching.
    A trace can pass through a junction;
    centerline pixels at true crossings may belong to more than one fiber.
    All raster connectors are required to remain inside the foreground mask.
    """
    mask = np.asarray(mask, bool)
    if mask.ndim != 2:
        raise ValueError('Expected a 2D mask')
    d = ndi.distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1]
    distance = ndi.gaussian_filter(d, .5, mode='constant')
    seeds = nucleate(mask, distance, seed_spacing)
    paths = []
    for seed in seeds:
        radius = int(np.clip(np.ceil(2 * distance[tuple(seed)]), 2, 16))
        first_points = _ring_candidates(seed, radius, distance, mask)
        if distance[tuple(seed)] < 1.5:
            # A one-pixel digital staircase has distance maxima orthogonal to
            # its actual axis. Estimate the initial axis from nearby foreground.
            lo = np.maximum(seed - 4, 0); hi = np.minimum(seed + 5, mask.shape)
            local = np.argwhere(mask[lo[0]:hi[0], lo[1]:hi[1]]) + lo - seed
            _, vectors = np.linalg.eigh(local.T @ local)
            axis = vectors[:, -1]
            candidates = _ring_candidates(seed, 1, distance, mask, True)
            first_points = []
            for sign in (-1, 1):
                forward = [q for q in candidates if np.dot(q-seed, sign*axis) > .2]
                if forward:
                    first_points.append(max(forward, key=lambda q: np.dot(q-seed, sign*axis)/np.linalg.norm(q-seed)))
        arms = [_extend(seed, first, distance, mask, max_angle) for first in first_points]
        # Opposite arms from the same seed form a single full fiber. Joining
        # only after global deduplication used to discard half of a trace.
        paths.extend(_join(_deduplicate(arms, distance, mask), mask, max_angle, tangent_points, distance))
    paths = _deduplicate(paths, distance, mask)
    paths = _join(paths, mask, max_angle, tangent_points, distance)
    return _deduplicate(paths, distance, mask), distance, seeds
