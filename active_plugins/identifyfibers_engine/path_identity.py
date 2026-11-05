"""Resolve long shared terminal and internal traces before object assignment.

This post-tracing stage changes centerline ownership, not the foreground mask.
Short crossing overlaps are retained. Internal identity ties use an explicit
geometric convention rather than allowing duplicate longitudinal ownership.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
from skimage.draw import line


def _arc(p):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]


def _tangent(p, span=5):
    k = np.arange(len(p))
    v = p[np.minimum(k+span, len(p)-1)] - p[np.maximum(k-span, 0)]
    return v/np.maximum(np.linalg.norm(v, axis=1)[:, None], 1e-12)


def _turn(p, index, span):
    if index < 3 or len(p)-1-index < 3:
        return None
    incoming = p[index]-p[max(0, index-span)]
    outgoing = p[min(len(p)-1, index+span)]-p[index]
    norm = np.linalg.norm(incoming)*np.linalg.norm(outgoing)
    if norm == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(incoming, outgoing)/norm, -1, 1))))


def _tail(path, other, distance, mask):
    """Candidate shared tail of path (caller also checks reversed path)."""
    separation, indices = cKDTree(other).query(path)
    near = separation <= 1.5
    # A second neighboring fiber across a background gap is not a duplicate.
    for k in np.flatnonzero(near & (separation > 0)):
        rr, cc = line(*path[k], *other[indices[k]])
        if not np.all(mask[rr, cc]):
            near[k] = False
    parallel = np.abs(np.sum(_tangent(path)*_tangent(other)[indices], axis=1)) >= np.cos(np.deg2rad(25))
    matched = near & parallel
    end = len(path)-1
    while end >= 0 and not matched[end] and len(path)-1-end <= 2:
        end -= 1
    if end < 0 or len(path)-1-end > 2 or not matched[end]:
        return None
    start = end
    while start > 0 and matched[start-1]:
        start -= 1
    arc = _arc(path)
    radius = float(np.median(distance[path[start:end+1, 0], path[start:end+1, 1]]))
    if arc[end]-arc[start] < max(8., 4*radius):
        return None
    # Tangent windows straddle the junction; recover the earlier meeting point
    # before assessing the incoming continuation, rather than cutting downstream.
    cut = start
    while cut > 0 and near[cut-1]:
        cut -= 1
    if cut < 3:
        return None  # whole-path duplicates belong to the existing deduplicator
    other_index = int(indices[cut])
    if indices[end] < other_index:
        other = other[::-1]
        other_index = len(other)-1-other_index
    span = max(8, int(np.ceil(4*radius)))
    turn = _turn(path, cut, span)
    other_turn = _turn(other, other_index, span)
    if turn is None or other_turn is None or turn-other_turn < 5.:
        return None
    # Preserve at least two points and never remove the entire incoming branch.
    return cut, turn-other_turn, arc[-1]-arc[cut]


def resolve_shared_tails(paths, distance, mask, minimum_length=0.):
    """Keep the smoother continuation when two traces reuse a long terminal arc.

    Matching requires <=1.5 px separation, axial agreement <=25 degrees,
    foreground support and shared length >=max(8 px, twice local full width).
    A >=5 degree continuation advantage is required; ties are not forced.
    Each accepted edit strictly shortens a tail, so iteration terminates.
    The winning path must satisfy minimum_length (in pixel arclength units).
    Returns paths and the number of resolved tail conflicts.
    """
    paths = [np.asarray(p, dtype=int).copy() for p in paths]
    resolved = 0
    while True:
        choices = []
        boxes = [(p.min(axis=0), p.max(axis=0)) for p in paths]
        lengths = [_arc(p)[-1] for p in paths]
        for i, p in enumerate(paths):
            if len(p) < 12:
                continue
            for j, q in enumerate(paths):
                if i == j or len(q) < 12 or lengths[j] < minimum_length:
                    continue
                if np.any(boxes[i][1]+2 < boxes[j][0]) or np.any(boxes[j][1]+2 < boxes[i][0]):
                    continue
                for reverse in (False, True):
                    oriented = p[::-1] if reverse else p
                    candidate = _tail(oriented, q, distance, mask)
                    if candidate is not None:
                        cut, advantage, removed_length = candidate
                        # Coordinate-based tie-breaking avoids object-ID priority.
                        key = tuple(oriented[0])+tuple(oriented[-1])
                        choices.append((-advantage, -removed_length, key, i, reverse, cut))
        if not choices:
            break
        _, _, _, i, reverse, cut = min(choices)
        p = paths[i][::-1] if reverse else paths[i]
        shortened = p[:cut+1]
        paths[i] = shortened[::-1] if reverse else shortened
        resolved += 1
    return paths, resolved


def _canonical(path):
    """Make geometric tie breaks independent of tracing order and direction."""
    p = np.asarray(path, dtype=int)
    return p if tuple(p.ravel()) <= tuple(p[::-1].ravel()) else p[::-1].copy()


def _prepare_signal(image, mask):
    """Light denoising for transverse valley checks, never for measurements."""
    if image is None:
        return None
    image = np.asarray(image, float)
    if image.shape != mask.shape:
        raise ValueError("Identity image and mask must have the same shape")
    valid = np.isfinite(image)
    if not np.any(valid & mask):
        return None
    low, high = np.percentile(image[valid & mask], (1, 99.5))
    if high <= low:
        return None  # no photometric evidence for a wider matching tolerance
    # Normalize without upper clipping: a bright peak must not erase a valley.
    signal = np.maximum(np.where(valid, image, low)-low, 0)/(high-low)
    return ndi.gaussian_filter(signal, .6, mode="nearest")


def _corridor_matches(p, q, distance, mask, signal=None, expand=True):
    """Return nearby supported samples and their nearest points on the other path.

    Wide matching is enabled only with a usable original-intensity image.
    Cross-profile valleys veto matching even if thresholding/closing connected
    two bright fibers. This is an image-resolution convention, not a proof of
    biological identity.
    """
    separation, indices = cKDTree(q).query(p)
    radius_p = distance[p[:, 0], p[:, 1]]
    radius_q = distance[q[indices, 0], q[indices, 1]]
    tolerance = (np.maximum(1.5, np.minimum(4., .75*np.minimum(radius_p, radius_q)))
                 if signal is not None and expand else np.full(len(p), 1.5))
    near = separation <= tolerance
    for k in np.flatnonzero(near & (separation > 0)):
        rr, cc = line(*p[k], *q[indices[k]])
        if not np.all(mask[rr, cc]):
            near[k] = False
            continue
        if signal is not None and len(rr) > 2:
            profile = signal[rr, cc]
            peak = min(profile[0], profile[-1])
            valley = float(np.min(profile[1:-1]))
            # Use relative contrast so dim but resolved fibers are protected
            # even when a separate bright bundle sets the image dynamic range.
            # Light smoothing and repeated longitudinal support reject speckle.
            if peak > 1e-6 and peak-valley > .20*peak:
                near[k] = False
    return near, indices, separation


def _candidate_runs(p, q, distance, mask, signal=None, expand=False):
    """Long, supported, parallel corridors, including internal overlaps.

    This is a centerline identity test, not an object-distance or dilation rule.
    Nearby fibers separated by background are never treated as one corridor.
    """
    near, indices, separation = _corridor_matches(p, q, distance, mask, signal, expand)
    parallel = np.abs(np.sum(_tangent(p)*_tangent(q)[indices], axis=1)) >= np.cos(np.deg2rad(25))
    edges = np.diff(np.r_[False, near, False].astype(int))
    ap, aq = _arc(p), _arc(q)
    for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        end = stop-1
        if end-start < 4 or np.mean(parallel[start:stop]) < .8:
            continue
        j0, j1 = sorted((int(indices[start]), int(indices[end])))
        # Require reciprocal support across the entire other interval: nearest
        # indices can jump over a real loop/detour and must not erase that arm.
        reverse_near, _, _ = _corridor_matches(q[j0:j1+1], p[start:stop], distance, mask, signal, expand)
        if not np.all(reverse_near):
            continue
        radius = float(np.median(distance[p[start:stop, 0], p[start:stop, 1]]))
        shared = min(ap[end]-ap[start], aq[j1]-aq[j0])
        # Between nearby junctions a genuine duplicate can span only one
        # local width. Requiring two widths left these short trunks with two
        # owners. Keep the absolute 8 px floor and all direction/support tests
        # so an ordinary point crossing is not classified as a shared trunk.
        if shared < max(8., 2*radius):
            continue
        wide = np.any(separation[start:stop] > 1.5)
        if wide:
            # Longer evidence is required when the paths are not coincident:
            # at least three local full widths, and 12 pixels of common arc.
            if shared < max(12., 6*radius):
                continue
        if signal is not None:
            tangent = _tangent(p)[start:stop]
            delta = q[indices[start:stop]]-p[start:stop]
            signed = tangent[:, 0]*delta[:, 1]-tangent[:, 1]*delta[:, 0]
            # Passing from one side to the other is a crossing, not a retrace.
            signed = signed[parallel[start:stop]]
            edge = max(3, len(signed)//4)
            first, last = np.median(signed[:edge]), np.median(signed[-edge:])
            # A sustained side change distinguishes a shallow crossing from
            # one-pixel raster jitter along an otherwise coincident run.
            if first*last < 0 and min(abs(first), abs(last)) > .75:
                continue
        # A folded path visiting the other fiber twice is not a shared arc.
        traveled = np.sum(np.abs(np.diff(aq[indices[start:stop]])))
        if traveled > 1.25*max(aq[j1]-aq[j0], 1.):
            continue
        yield start, end, j0, j1, shared, radius


def _shared_runs(p, q, distance, mask, signal=None):
    # Preserve established coincident-run repairs. Expanding the search alone
    # can otherwise join a valid narrow run to a nearby diverging arm, causing
    # the stricter wide-run tests to reject the entire original candidate.
    seen = set()
    for expand in (False, True) if signal is not None else (False,):
        for run in _candidate_runs(p, q, distance, mask, signal, expand):
            key = run[:4]
            if key not in seen:
                seen.add(key)
                yield run


def resolve_shared_segments(paths, distance, mask, minimum_length=0., image=None):
    """Give each long shared corridor one centerline owner before propagation.

    Retain the path with the smoother continuation through the corridor. If
    turns differ by less than 5 degrees, retain the longer path, then use a
    geometric tie break. This tie break is computational, not biological proof.
    Split the competing path at the shared interval; keep its exclusive arms
    as separate candidate paths, subject to the normal minimum-length filter.
    Junction endpoints may coincide; accepted shared runs get one owner.
    With an original image, nearby noncoincident paths can also match within
    0.75 times the smaller local radius (capped at 4 px). Such runs require
    >=max(12 px, three local widths), no resolved intensity valley and no
    side-switch crossing. Without image evidence, the original 1.5 px limit
    remains. Neither the mask nor the measured input intensity is modified.
    """
    signal = _prepare_signal(image, mask)
    margin = 4. if signal is not None else 2.
    paths = sorted((_canonical(p).copy() for p in paths if len(p) >= 2),
                   key=lambda p: tuple(p.ravel()))
    resolved = 0
    while True:
        choices = []
        boxes = [(p.min(axis=0), p.max(axis=0)) for p in paths]
        lengths = [_arc(p)[-1] for p in paths]
        for i, p in enumerate(paths):
            if len(p) < 6:
                continue
            for j in range(i+1, len(paths)):
                q = paths[j]
                if len(q) < 6 or max(lengths[i], lengths[j]) < minimum_length:
                    continue
                if np.any(boxes[i][1]+margin < boxes[j][0]) or np.any(boxes[j][1]+margin < boxes[i][0]):
                    continue
                for a, b, c, d, shared, radius in _shared_runs(p, q, distance, mask, signal):
                    span = max(8, int(np.ceil(4*radius)))
                    turns = []
                    for path, lo, hi in ((p, a, b), (q, c, d)):
                        values = [_turn(path, k, span) for k in (lo, hi)]
                        values = [v for v in values if v is not None]
                        turns.append(float(np.mean(values)) if values else 0.)
                    if lengths[i] < minimum_length:
                        winner = j
                    elif lengths[j] < minimum_length:
                        winner = i
                    elif abs(turns[0]-turns[1]) >= 5.:
                        winner = i if turns[0] < turns[1] else j
                    else:
                        winner = min((i, j), key=lambda k: (-lengths[k], tuple(paths[k].ravel())))
                    loser, lo, hi = (j, c, d) if winner == i else (i, a, b)
                    key = tuple(paths[winner].ravel()), tuple(paths[loser].ravel())
                    choices.append((-shared, key, loser, lo, hi))
        if not choices:
            break
        _, _, loser, lo, hi = min(choices)
        p = paths.pop(loser)
        # Keep junction points as endpoints, never reconnect the two exclusive
        # arms through the just-removed common arc.
        paths.extend(_canonical(arm).copy() for arm in (p[:lo+1], p[hi:]) if len(arm) >= 2)
        resolved += 1
    return paths, resolved
