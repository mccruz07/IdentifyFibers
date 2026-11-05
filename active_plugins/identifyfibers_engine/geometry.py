"""Exclusive, centerline-seeded regions and geometric quality diagnostics.

Independent region partition. Separate path-identity checks use intensity upstream.
"""
import heapq
import colorsys
import numpy as np
from scipy import ndimage as ndi


def centerline_direction_tensor(fibers, shape):
    """Continuous local tangent information, averaged as axial tensors.

    The tensor is built from entire accepted centerlines, never their centroids.
    Equal competing orthogonal directions cancel the anisotropic term, so an
    unresolved crossing does not arbitrarily inherit the last fiber's direction.
    Nearest-line extension is only a direction prior; propagation stays masked.
    """
    yy, yx, xx = [np.zeros(shape, float) for _ in range(3)]
    count = np.zeros(shape, float)
    for fiber in fibers:
        p = np.asarray(fiber.coordinates, dtype=int)
        if len(p) < 2:
            continue
        index = np.arange(len(p))
        tangent = p[np.minimum(index+3, len(p)-1)] - p[np.maximum(index-3, 0)]
        length = np.linalg.norm(tangent, axis=1)
        good = length > 0
        p = p[good]; tangent = tangent[good]/length[good, None]
        indices = (p[:, 0], p[:, 1])
        np.add.at(yy, indices, tangent[:, 0]**2)
        np.add.at(yx, indices, tangent[:, 0]*tangent[:, 1])
        np.add.at(xx, indices, tangent[:, 1]**2)
        np.add.at(count, indices, 1)
    support = count > 0
    if not support.any():
        return yy, yx, xx
    # Smooth only locally around the paths, then extend their tangent estimates.
    weight = ndi.gaussian_filter(count, 1., mode="constant")
    tensors = [ndi.gaussian_filter(a, 1., mode="constant")/np.maximum(weight, 1e-12)
               for a in (yy, yx, xx)]
    nearest = ndi.distance_transform_edt(~support, return_distances=False, return_indices=True)
    a, b, c = [t[tuple(nearest)] for t in tensors]
    trace = a+c
    contrast = np.sqrt((a-c)**2+4*b*b)
    minor = (trace-contrast)/2
    # PSD tensor with eigenvalues [coherence, 0]. Straight axes are strong;
    # isotropic crossings have negligible directional preference.
    return ((a-minor)/np.maximum(trace, 1e-12), b/np.maximum(trace, 1e-12),
            (c-minor)/np.maximum(trace, 1e-12))


def partition_regions(mask, markers, intensity, direction_tensor=None, direction_weight=1.):
    """Multi-source geodesic Voronoi partition of a foreground mask.

    Every centerline marker is fixed. Eight-connected path cost is Euclidean
    step length multiplied by (1 + 2*absolute normalized intensity change +
    direction_weight*axial tangent penalty). The shared direction field makes
    this a label-independent geodesic metric, so ordinary multi-source Dijkstra
    is valid. Growth along tangents is penalized relative to transverse growth;
    full centerlines already seed the longitudinal dimension of each fiber.
    Ties use object ID deterministically. Unseeded components remain zero.
    No propagation into background; diagonal foreground neighbors are connected.
    """
    if not np.isfinite(direction_weight) or direction_weight < 0:
        raise ValueError("Direction guidance must be finite and nonnegative")
    mask = np.asarray(mask, bool)
    labels = np.where(mask, markers, 0).astype(np.int32)
    values = np.asarray(intensity, float)
    finite = mask & np.isfinite(values)
    normalized = np.zeros(mask.shape, float)
    if finite.any():
        lo, hi = np.percentile(values[finite], (1, 99))
        if hi > lo:
            normalized[finite] = np.clip((values[finite]-lo)/(hi-lo), 0, 1)
    costs = np.full(mask.shape, np.inf)
    seeds = np.argwhere(labels > 0)
    queue = [(0., int(labels[y, x]), int(y), int(x)) for y, x in seeds]
    costs[labels > 0] = 0
    heapq.heapify(queue)
    height, width = mask.shape
    steps = [(dy, dx, np.hypot(dy, dx)) for dy in (-1, 0, 1)
             for dx in (-1, 0, 1) if dy or dx]
    while queue:
        cost, owner, y, x = heapq.heappop(queue)
        if cost != costs[y, x] or owner != labels[y, x]:
            continue
        for dy, dx, length in steps:
            yy, xx = y + dy, x + dx
            if not (0 <= yy < height and 0 <= xx < width) or not mask[yy, xx]:
                continue
            if markers[yy, xx] > 0:
                continue
            penalty = 0.
            if direction_tensor is not None and direction_weight:
                a, b, c = direction_tensor
                # Average the tensor at both ends to keep edge costs symmetric.
                penalty = ((a[y,x]+a[yy,xx])*dy*dy +
                           2*(b[y,x]+b[yy,xx])*dy*dx +
                           (c[y,x]+c[yy,xx])*dx*dx)/(2*length*length)
            candidate = cost + length * (1 + 2*abs(normalized[y, x]-normalized[yy, xx]) +
                                          direction_weight*max(0., penalty))
            previous = costs[yy, xx]
            if candidate < previous or (candidate == previous and owner < labels[yy, xx]):
                costs[yy, xx] = candidate
                labels[yy, xx] = owner
                heapq.heappush(queue, (candidate, owner, yy, xx))
    return labels


def coverage_diagnostic(reference_skeleton, accepted_paths, tolerance=2.):
    """Percent of reference skeleton pixels <= tolerance pixels from a path.

    Empty reference is undefined (NaN), not perfect coverage. The reference is
    a skeleton of the segmentation, not biological ground truth.
    """
    reference = np.asarray(reference_skeleton, bool)
    paths = np.asarray(accepted_paths, bool)
    near = ndi.distance_transform_edt(~paths) <= tolerance if paths.any() else np.zeros_like(paths)
    uncovered = reference & ~near
    total = int(reference.sum())
    percentage = 100.*(total-int(uncovered.sum()))/total if total else float('nan')
    return percentage, uncovered


def label_overlay(labels, image=None, alpha=1.):
    """Stable object-ID palette, even when different subsets are displayed."""
    labels = np.asarray(labels, dtype=int)
    palette = np.zeros((int(labels.max(initial=0))+1, 3))
    for number in range(1, len(palette)):
        palette[number] = colorsys.hsv_to_rgb((number*.61803398875) % 1, .85, 1.)
    base = np.zeros(labels.shape) if image is None else np.clip(np.asarray(image, float), 0, 1)
    rgb = np.repeat(base[..., None], 3, axis=2)
    foreground = labels > 0
    rgb[foreground] = (1-alpha)*rgb[foreground] + alpha*palette[labels[foreground]]
    return rgb
