"""Calibrated Sato enhancement, 3D centerlines and exclusive voxel ownership."""

from dataclasses import dataclass
from heapq import heapify, heappop, heappush

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from .volume_trace import OFFSETS, traces


@dataclass(frozen=True)
class Settings3D:
    width_min: float = 0.5
    width_max: float = 2.0
    min_length: float = 3.0
    correction: float = 1.0
    max_angle: float = 35.0


@dataclass
class VolumeResult:
    labels: np.ndarray
    centers: np.ndarray
    enhanced: np.ndarray
    fibers: list


def normalize(image, valid):
    out = np.zeros(image.shape, np.float32)
    values = image[valid]
    if values.size:
        lo, hi = np.percentile(values, [1, 99.5])
        if hi <= lo:
            lo, hi = values.min(), values.max()
        if hi > lo:
            out[valid] = np.clip((values.astype(float) - lo) / (hi - lo), 0, 1)
    return out


def _derivative(a, axis, kernel):
    return ndi.correlate1d(a, kernel, axis, mode="constant")


def tubeness(image, valid, spacing, settings):
    """Multiscale Sato tubular response from a physical-coordinate Hessian."""
    original = normalize(image, valid)
    out = np.zeros_like(original)
    if not original.any():
        return out
    for sigma in np.unique(
        np.geomspace(settings.width_min / 2, settings.width_max / 2, 3)
    ):
        scales = sigma / spacing
        # Limit temporary Hessian/eigenvalue arrays by processing slabs with halos.
        depth = max(1, 180000 // (image.shape[1] * image.shape[2]))
        halo = int(4 * scales[0] + 0.5) + 2
        for z in range(0, len(image), depth):
            stop = min(len(image), z + depth)
            lo, hi = max(0, z - halo), min(len(image), stop + halo)
            crop = slice(z - lo, stop - lo)
            smooth = ndi.gaussian_filter(original[lo:hi], scales, mode="constant")
            h = np.empty((stop - z,) + image.shape[1:] + (3, 3), np.float32)
            for i in range(3):
                h[..., i, i] = (
                    _derivative(smooth, i, [1, -2, 1])[crop] * (sigma / spacing[i]) ** 2
                )
                for j in range(i):
                    v = _derivative(
                        _derivative(smooth, i, [-0.5, 0, 0.5]), j, [-0.5, 0, 0.5]
                    )[crop]
                    h[..., i, j] = h[..., j, i] = (
                        v * sigma**2 / (spacing[i] * spacing[j])
                    )
            eigen = np.linalg.eigvalsh(h)
            value = np.sqrt(
                np.maximum(-eigen[..., 0], 0) * np.maximum(-eigen[..., 1], 0)
            )
            np.maximum(out[z:stop], value, out=out[z:stop])
    return normalize(out, valid)


def _contrast_evidence(raw, valid, spacing, cfg):
    """Estimate local bright-ridge contrast and the noise of its small-scale average."""
    signal = normalize(raw, valid)
    estimates = []
    for axis in range(3):
        a = np.moveaxis(signal, axis, 0)
        domain = np.moveaxis(valid, axis, 0)
        differences = (a[1:] - a[:-1])[domain[1:] & domain[:-1]]
        if differences.size >= 16:
            median = np.median(differences)
            estimates.append(np.median(np.abs(differences - median)) / 0.95387255)
    if not estimates or min(estimates) <= 0:
        return None
    small = (cfg.width_min / 2) / spacing
    large = (2 * cfg.width_max) / spacing

    def mean(sigma):
        weight = ndi.gaussian_filter(valid.astype(np.float32), sigma, mode="constant")
        return ndi.gaussian_filter(signal, sigma, mode="constant") / np.maximum(
            weight, 1e-7
        )

    contrast = mean(small) - mean(large)
    gain = 1.0
    for sigma in small:
        radius = int(4 * sigma + 0.5)
        x = np.arange(-radius, radius + 1)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        gain *= np.linalg.norm(kernel) / kernel.sum()
    return contrast, 3 * min(estimates) * gain


def _assign(mask, paths, image, spacing):
    labels = np.zeros(mask.shape, np.int32)
    for i, path in enumerate(paths, 1):
        points = tuple(path.T)
        labels[points] = np.where(labels[points] == 0, i, labels[points])
    if len(paths) == 1:
        labels[mask] = 1
        return labels
    fixed = labels > 0
    signal = normalize(image, mask)
    costs = np.full(mask.shape, np.inf)
    costs[fixed] = 0
    queue = [(0.0, int(labels[tuple(p)]), *map(int, p)) for p in np.argwhere(fixed)]
    heapify(queue)
    steps = [(d, np.linalg.norm(spacing * d)) for d in OFFSETS]
    nz, ny, nx = mask.shape
    while queue:
        cost, owner, z, y, x = heappop(queue)
        if cost != costs[z, y, x] or labels[z, y, x] != owner:
            continue
        for (dz, dy, dx), length in steps:
            zz, yy, xx = z + dz, y + dy, x + dx
            if not (0 <= zz < nz and 0 <= yy < ny and 0 <= xx < nx):
                continue
            q = (zz, yy, xx)
            if not mask[q] or fixed[q]:
                continue
            value = cost + length * (1 + 2 * abs(float(signal[q]) - signal[z, y, x]))
            if value < costs[q] or (value == costs[q] and owner < labels[q]):
                labels[q], costs[q] = owner, value
                heappush(queue, (value, owner, *q))
    return labels


def _curve_metrics(path, spacing):
    """Suppress voxel stair steps on a physical arc-length grid; retain both endpoints."""
    xyz = path * spacing
    steps = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    arc = np.r_[0.0, np.cumsum(steps)]
    if len(path) < 3 or arc[-1] == 0:
        return xyz, arc[-1], arc[-1]
    grid = np.linspace(0, arc[-1], len(path))
    sampled = np.column_stack([np.interp(grid, arc, v) for v in xyz.T])
    # One traversed voxel step in physical units sets the smoothing scale.
    sigma = min(steps.max() / grid[1], (len(path) - 1) / 4)
    if np.array_equal(path[0], path[-1]):
        smoothed = ndi.gaussian_filter1d(sampled[:-1], sigma, axis=0, mode="wrap")
        smoothed = np.vstack((smoothed, smoothed[0]))
    else:
        pad = int(np.ceil(4 * sigma)) + 1
        extended = np.pad(
            sampled, ((pad, pad), (0, 0)), mode="reflect", reflect_type="odd"
        )
        smoothed = ndi.gaussian_filter1d(extended, sigma, axis=0)[pad:-pad]
        smoothed[[0, -1]] = xyz[[0, -1]]
    return smoothed, np.linalg.norm(np.diff(smoothed, axis=0), axis=1).sum(), arc[-1]


def _orientation(curve):
    xyz = curve[:, ::-1].copy()
    xyz -= xyz.mean(axis=0)
    values, axes = np.linalg.eigh(xyz.T @ xyz)
    if values[-1] <= 0 or values[-1] - values[-2] <= 1e-6 * values[-1]:
        return np.nan, np.nan
    x, y, z = axes[:, -1]
    if np.hypot(x, y) < 1e-8:
        return np.nan, 90.0
    if y < -1e-12 or (abs(y) <= 1e-12 and x < 0):
        x, y, z = -x, -y, -z
    return np.degrees(np.arctan2(y, x)) % 180, np.degrees(np.arctan2(z, np.hypot(x, y)))


FEATURES = (
    "PathLength_um",
    "VoxelPathLength_um",
    "EndToEnd_um",
    "Straightness",
    "Volume_um3",
    "ApparentWidth_um",
    "Azimuth_deg",
    "Inclination_deg",
    "IntensityMean",
    "IntensityIntegrated",
    "TouchesImageEdge",
    "FallbackTrace",
)


def analyze_volume(image, spacing, settings=None, mask=None, parents=None, binary=None):
    image = np.asarray(image)
    settings = settings or Settings3D()
    spacing = np.asarray(spacing, float)
    if image.ndim != 3 or min(image.shape) < 1:
        raise ValueError("Expected a grayscale volume with axes Z,Y,X.")
    if image.shape[0] < 2:
        raise ValueError(
            "Volume mode requires at least two Z planes. Use 2D image for one plane."
        )
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or min(spacing) <= 0:
        raise ValueError("Specify three positive finite voxel sizes in Z,Y,X order.")
    if not np.isfinite(list(vars(settings).values())).all() or not (
        0 < settings.width_min <= settings.width_max
        and settings.min_length >= 0
        and settings.correction > 0
        and 0 <= settings.max_angle <= 90
    ):
        raise ValueError("Check fiber widths, length, threshold factor and angle.")
    for item in (mask, parents, binary):
        if item is not None and np.shape(item) != image.shape:
            raise ValueError(
                "All masks and parent objects must match the input Z,Y,X dimensions."
            )
    if parents is not None:
        parents = np.asarray(parents)
        if not np.issubdtype(parents.dtype, np.integer) or np.any(parents < 0):
            raise ValueError("Parent labels must be nonnegative integers.")
    valid = np.isfinite(image)
    if mask is not None:
        valid &= np.asarray(mask, bool)
    if binary is not None:
        binary = np.asarray(binary)
        if binary.dtype != bool:
            raise ValueError("External segmentation must be binary.")
    result = VolumeResult(
        np.zeros(image.shape, np.int32),
        np.zeros(image.shape, bool),
        np.zeros(image.shape, np.float32),
        [],
    )
    if parents is None:
        regions = [(0, tuple(slice(0, n) for n in image.shape))]
    else:
        ids, inverse = np.unique(parents[valid & (parents > 0)], return_inverse=True)
        dense = np.zeros(image.shape, np.int32)
        dense[valid & (parents > 0)] = inverse + 1
        regions = zip(ids, ndi.find_objects(dense))
    for parent, box in regions:
        roi = valid[box].copy()
        if parents is not None:
            roi &= parents[box] == parent
        raw = image[box]
        evidence = None
        if binary is None:
            response = tubeness(raw, roi, spacing, settings)
            result.enhanced[box][roi] = response[roi]
            positives = response[roi & (response > 0)]
            if not positives.size:
                continue
            threshold = threshold_otsu(positives) * settings.correction
            foreground = roi & (response > 0) & (response >= threshold)
            if foreground.any():
                evidence = _contrast_evidence(raw, roi, spacing, settings)
        else:
            foreground = roi & binary[box]
        components, n = ndi.label(foreground, structure=np.ones((3, 3, 3)))
        for number, sub in enumerate(ndi.find_objects(components), 1):
            comp = components[sub] == number
            if comp.sum() < 2:
                continue
            offset = np.array([b.start + s.start for b, s in zip(box, sub)])
            distance = ndi.distance_transform_edt(np.pad(comp, 1), sampling=spacing)[
                1:-1, 1:-1, 1:-1
            ]
            candidates, fallback = traces(
                comp,
                distance,
                spacing,
                settings.max_angle,
                max(settings.width_max, 2 * min(spacing)),
            )
            paths = []
            lengths = []
            for path in candidates:
                curve, length, voxel_length = _curve_metrics(path, spacing)
                if evidence is not None:
                    contrast, floor = evidence
                    if np.median(contrast[sub][tuple(path.T)]) < floor:
                        continue
                if length >= settings.min_length and length > 0:
                    paths.append(path)
                    lengths.append((curve, length, voxel_length))
            if not paths:
                continue
            grown = _assign(comp, paths, raw[sub], spacing)
            bounds = ndi.find_objects(grown, max_label=len(paths))
            for i, (path, (curve, length, voxel_length)) in enumerate(
                zip(paths, lengths), 1
            ):
                box_i = bounds[i - 1]
                if box_i is None:
                    continue
                vox = np.argwhere(grown[box_i] == i) + [s.start for s in box_i]
                if not len(vox):
                    continue
                original = raw[sub][tuple(vox.T)]
                end = np.linalg.norm((path[-1] - path[0]) * spacing)
                az, inc = _orientation(curve)
                widths = 2 * distance[tuple(path.T)]
                points = path + offset
                vox += offset
                label = len(result.fibers) + 1
                result.labels[tuple(vox.T)] = label
                result.centers[tuple(points.T)] = True
                border = np.any((vox == 0) | (vox == np.asarray(image.shape) - 1))
                values = (
                    length,
                    voxel_length,
                    end,
                    min(1.0, end / length),
                    len(vox) * spacing.prod(),
                    np.mean(widths),
                    az,
                    inc,
                    np.mean(original, dtype=float),
                    np.sum(original, dtype=float),
                    int(border),
                    int(fallback),
                )
                result.fibers.append(
                    dict(
                        parent=int(parent),
                        path=points,
                        centroid=vox.mean(axis=0),
                        **dict(zip(FEATURES, values)),
                    )
                )
    return result
