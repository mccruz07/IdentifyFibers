"""Distance-guided 3D centerlines. Each non-junction branch is consumed once."""

from itertools import product
from heapq import heappop, heappush
import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import skeletonize_3d

OFFSETS = tuple(p for p in product((-1, 0, 1), repeat=3) if any(p))


def _graph(points):
    points = set(map(tuple, points))
    graph = {p: set() for p in points}
    for p in points:
        for d in OFFSETS:
            q = tuple(a + b for a, b in zip(p, d))
            if q <= p or q not in points:
                continue
            # Suppress corner shortcuts when a route through a face neighbor exists.
            box = product(*[(a, b) if a != b else (a,) for a, b in zip(p, q)])
            if any(r != p and r != q and r in points for r in box):
                continue
            graph[p].add(q)
            graph[q].add(p)
    return graph


def _components(nodes, graph):
    todo = set(nodes)
    groups = []
    while todo:
        start = min(todo)
        todo.remove(start)
        group, queue = {start}, [start]
        while queue:
            for q in graph[queue.pop()]:
                if q in todo:
                    todo.remove(q)
                    queue.append(q)
                    group.add(q)
        groups.append(group)
    return groups


def _connector(a, b, zone, graph, spacing):
    if a == b:
        return [a]
    queue, costs, previous = [(0.0, a)], {a: 0.0}, {}
    while queue:
        cost, p = heappop(queue)
        if p == b:
            route = [b]
            while route[-1] != a:
                route.append(previous[route[-1]])
            return route[::-1]
        if cost != costs[p]:
            continue
        for q in graph[p] & zone:
            value = cost + np.linalg.norm((np.subtract(q, p)) * spacing)
            if value < costs.get(q, np.inf):
                costs[q] = value
                previous[q] = p
                heappush(queue, (value, q))
    raise RuntimeError("Disconnected junction")


def _ordered_paths(graph, spacing, angle_limit, tangent_range):
    # Collapse a cluster of adjacent junction voxels to one decision.
    zones = _components((p for p, v in graph.items() if len(v) > 2), graph)
    zone_of = {p: i for i, z in enumerate(zones) for p in z}
    branches, arms = [], [[] for z in zones]
    junction_voxels = set(zone_of)
    for chain in _components(set(graph) - junction_voxels, graph):
        tips = sorted(p for p in chain if len(graph[p] & chain) < 2)
        start = tips[0] if tips else min(chain)
        path, used = [start], {start}
        while True:
            options = (graph[path[-1]] & chain) - used
            if not options:
                break
            q = min(options)
            path.append(q)
            used.add(q)
        if not tips and len(path) > 2:
            path.append(start)
        contacts = [sorted(graph[p] & junction_voxels) for p in (path[0], path[-1])]
        if len(path) == 1:
            # A singleton may connect two decision zones.
            contacts = [contacts[0][:1], contacts[0][1:2]]
        for side, choices in enumerate(contacts):
            if choices:
                anchor = choices[0]
                arms[zone_of[anchor]].append((len(branches), side))
                if side:
                    path.append(anchor)
                else:
                    path.insert(0, anchor)
        branches.append(path)
    links = {}
    for zone, ends in zip(zones, arms):
        vectors = {}
        for edge, side in ends:
            path = branches[edge][::-1] if side else branches[edge]
            xyz = np.asarray(path) * spacing
            arc = np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))
            j = min(np.searchsorted(arc, tangent_range) + 1, len(path) - 1)
            # The first voxel after a broad junction can be a diagonal shortcut.
            # Estimate the continuation along the arm beyond that voxel.
            base = 1 if j >= 3 else 0
            v = xyz[j] - xyz[base]
            if np.linalg.norm(v) > 0:
                vectors[edge, side] = v / np.linalg.norm(v)
        options = []
        for a in sorted(vectors):
            for b in sorted(vectors):
                if a >= b or a[0] == b[0]:
                    continue
                angle = np.degrees(np.arccos(np.clip(-vectors[a] @ vectors[b], -1, 1)))
                if angle <= angle_limit + 1e-8:
                    options.append((angle, a, b))
        for _, a, b in sorted(options):
            if a in links or b in links:
                continue
            links[a], links[b] = (b, zone), (a, zone)
    visited = set()
    starts = [(i, s) for i in range(len(branches)) for s in (0, 1)]
    for start in sorted(starts, key=lambda k: (k in links, k)):
        edge, side = start
        if edge in visited:
            continue
        result = []
        while edge not in visited:
            visited.add(edge)
            segment = branches[edge][::-1] if side else branches[edge]
            result.extend(segment if not result else segment[1:])
            link = links.get((edge, 1 - side))
            if link is None:
                break
            (edge, side), zone = link
            next_anchor = branches[edge][-1 if side else 0]
            result.extend(_connector(result[-1], next_anchor, zone, graph, spacing)[1:])
        if len(result) > 1:
            yield np.asarray(result, dtype=int)


def _fallback(mask, distance, spacing):
    """A center-favoring geodesic diameter when thinning erases an entire tube."""
    from skimage.graph import MCP_Geometric

    cost = np.full(mask.shape, np.inf)
    cost[mask] = 1 / np.maximum(distance[mask], 0.5 * min(spacing))
    mcp = MCP_Geometric(cost, sampling=tuple(spacing))
    seed = np.unravel_index(np.argmax(distance), mask.shape)
    ridge = mask & (distance == ndi.maximum_filter(distance, size=3))
    targets = ridge if ridge.sum() > 1 else mask
    for _ in range(2):
        costs, _ = mcp.find_costs([seed])
        seed = np.unravel_index(
            np.argmax(np.where(targets & np.isfinite(costs), costs, -1)), mask.shape
        )
    return np.asarray(mcp.traceback(seed), dtype=int)


def traces(mask, distance, spacing, angle_limit, tangent_range):
    # Preserve digital lines that already have no filled voxel faces.
    # In particular, Lee thinning can displace the meeting point of a thin T.
    has_face = False
    for a, b in ((0, 1), (0, 2), (1, 2)):
        grid = np.moveaxis(mask, (a, b), (0, 1))
        if np.any(grid[:-1, :-1] & grid[1:, :-1] & grid[:-1, 1:] & grid[1:, 1:]):
            has_face = True
            break
    sk = skeletonize_3d(np.ascontiguousarray(mask)).astype(bool) if has_face else mask
    graph = _graph(np.argwhere(sk))
    paths = list(_ordered_paths(graph, spacing, angle_limit, tangent_range))
    if not paths and mask.sum() > 1:
        # Thinning can erase an even-width tube. A sheet or compact blob
        # also has a geodesic diameter, but it gives no unique longitudinal axis.
        xyz = np.argwhere(mask) * spacing
        xyz -= xyz.mean(axis=0)
        eigen = np.linalg.eigvalsh(xyz.T @ xyz)
        if eigen[-1] >= 4 * eigen[-2] and eigen[-1] > 0:
            path = _fallback(mask, distance, spacing)
            if len(path) > 1:
                return [path], True
    return paths, False
