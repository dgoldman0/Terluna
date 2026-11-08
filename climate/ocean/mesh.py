"""Triangle meshes of the seas for SCHISM, from the terrain's pixels.

The mesh is a quadtree on the 16 px/deg terrain: each base cell of 8 x 8
pixels (half a degree) splits into leaves of 8, 4, 2 or 1 pixels as a level
map asks, with neighbouring base cells at most one level apart. A leaf is
water where at least half its pixels belong to the chosen sea, with their
mean depth. A leaf whose edges meet no finer neighbour splits along a
diagonal into two triangles; one with finer neighbours fans from a node at
its centre through the neighbours' corner nodes, so the mesh has no hanging
nodes. The coastline follows the leaves' edges, as the MITgcm grid's does,
so the two models share the geometry at equal resolution and differ in their
numerics.

Where leaves touch only at a corner, each side gets its own node there, as
water cannot pass through a corner in either model; SCHISM then sees each
boundary node on one boundary. Node depths are the mean depth of the leaves
that share the node.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BASE = 8          # pixels in a base cell's side (half a degree at 16 px/deg)


@dataclass
class Mesh:
    lon: np.ndarray          # node longitudes (deg east)
    lat: np.ndarray          # node latitudes
    depth: np.ndarray        # node depths (m, positive)
    triangles: np.ndarray    # (ne, 3) node indices, anticlockwise in lon/lat

    @property
    def nodes(self):
        return len(self.lon)

    @property
    def elements(self):
        return len(self.triangles)

    def element_areas(self, radius):
        """Element areas (m^2) by the local metric at each element's mean latitude."""
        la = np.radians(self.lat[self.triangles])
        lo = np.radians(self.lon[self.triangles])
        c = np.cos(la.mean(axis=1))
        x = radius * lo * c[:, None]
        y = radius * la
        return 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))

    def boundary_loops(self):
        """Closed loops of boundary nodes, each in order; the first is the longest (the exterior)."""
        edges = np.sort(np.concatenate([self.triangles[:, [0, 1]], self.triangles[:, [1, 2]],
                                        self.triangles[:, [2, 0]]]), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        boundary = unique[counts == 1]
        neighbours = {}
        for a, b in boundary:
            neighbours.setdefault(int(a), []).append(int(b))
            neighbours.setdefault(int(b), []).append(int(a))
        if any(len(v) != 2 for v in neighbours.values()):
            raise ValueError("mesh boundary is not a set of simple loops")
        loops, seen = [], set()
        for start in neighbours:
            if start in seen:
                continue
            loop, prev, node = [start], start, neighbours[start][0]
            seen.add(start)
            while node != start:
                loop.append(node)
                seen.add(node)
                a, b = neighbours[node]
                prev, node = node, (b if a == prev else a)
            loops.append(loop)
        loops.sort(key=len, reverse=True)
        return loops


def balance(level):
    """Raise levels until face neighbours differ by at most one."""
    level = level.copy()
    while True:
        padded = np.pad(level, 1, mode="edge")
        neighbour = np.maximum.reduce([padded[:-2, 1:-1], padded[2:, 1:-1], padded[1:-1, :-2], padded[1:-1, 2:]])
        raised = np.maximum(level, neighbour - 1)
        if np.array_equal(raised, level):
            return level
        level = raised


def build(depth_px, sea_px, level, lat_px, lon_px, min_depth=5.0):
    """A mesh from pixel depths (m) and the sea's pixels, with a level (0-3) per base cell.

    depth_px, sea_px: (ny, nx) on the terrain's pixels (latitudes ascending); level: (ny // 8, nx // 8).
    """
    level = balance(np.asarray(level, int))
    water = sea_px & (depth_px > 0)
    cw = np.pad(np.cumsum(np.cumsum(water, 0), 1), ((1, 0), (1, 0)))
    cd = np.pad(np.cumsum(np.cumsum(np.where(water, depth_px, 0.0), 0), 1), ((1, 0), (1, 0)))
    def box_sum(c, r0, r1, c0, c1):
        return c[r1, c1] - c[r0, c1] - c[r1, c0] + c[r0, c0]
    leaves = []                          # (row0, col0, size, depth) in pixels
    for bj in range(level.shape[0]):
        for bi in range(level.shape[1]):
            size = BASE >> int(level[bj, bi])
            for r0 in range(bj * BASE, (bj + 1) * BASE, size):
                for c0 in range(bi * BASE, (bi + 1) * BASE, size):
                    n = box_sum(cw, r0, r0 + size, c0, c0 + size)
                    if n * 2 < size * size:
                        continue
                    d = box_sum(cd, r0, r0 + size, c0, c0 + size) / n
                    if d >= min_depth:
                        leaves.append((r0, c0, size, float(d)))
    leaves = largest_component(leaves)
    step = abs(lat_px[1] - lat_px[0])
    return triangulate(leaves, lat_px[0] - step / 2, lon_px[0] - step / 2, step)


def from_cells(depth, lat_c, lon_c):
    """A mesh of a regular grid's wet cells, two triangles each (the checks on the atlas's cells)."""
    leaves = [(int(j), int(i), 1, float(depth[j, i])) for j, i in zip(*np.nonzero(depth > 0))]
    step = lat_c[1] - lat_c[0]
    return triangulate(leaves, lat_c[0] - step / 2, lon_c[0] - (lon_c[1] - lon_c[0]) / 2, step,
                       lon_step=lon_c[1] - lon_c[0])


def edge_pairs(leaves):
    """Pairs of leaves sharing (part of) an edge, from their unit edges."""
    edges = {}
    for k, (r0, c0, s, _) in enumerate(leaves):
        for r in range(r0, r0 + s):
            edges.setdefault((0, r, c0), []).append(k)
            edges.setdefault((0, r, c0 + s), []).append(k)
        for c in range(c0, c0 + s):
            edges.setdefault((1, r0, c), []).append(k)
            edges.setdefault((1, r0 + s, c), []).append(k)
    return {(min(ks), max(ks)) for ks in edges.values() if len(ks) == 2}


def largest_component(leaves):
    parent = list(range(len(leaves)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edge_pairs(leaves):
        parent[find(a)] = find(b)
    roots = [find(k) for k in range(len(leaves))]
    weight = {}
    for k, root in enumerate(roots):
        weight[root] = weight.get(root, 0) + leaves[k][2] ** 2
    keep = max(weight, key=weight.get)
    return [leaf for k, leaf in enumerate(leaves) if roots[k] == keep]


def triangulate(leaves, lat0, lon0, step, lon_step=None):
    """Triangles over square leaves on a lattice of `step` degrees from the corner (lat0, lon0)."""
    lon_step = step if lon_step is None else lon_step
    corners = set()
    for r0, c0, s, _ in leaves:
        corners.update(((r0, c0), (r0, c0 + s), (r0 + s, c0), (r0 + s, c0 + s)))
    rings = []
    touching = {}                        # lattice point -> leaves whose ring holds it
    for k, (r0, c0, s, _) in enumerate(leaves):
        ring = []
        for (ra, ca), (rb, cb) in (((r0, c0), (r0, c0 + s)), ((r0, c0 + s), (r0 + s, c0 + s)),
                                   ((r0 + s, c0 + s), (r0 + s, c0)), ((r0 + s, c0), (r0, c0))):
            ring.append((ra, ca))
            if s > 1 and s % 2 == 0:
                mid = ((ra + rb) // 2, (ca + cb) // 2)
                if mid in corners:
                    ring.append(mid)
        rings.append(ring)
        for key in ring:
            touching.setdefault(key, []).append(k)
    adjacency = {}
    for a, b in edge_pairs(leaves):
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    node_of = {}                         # (point, leaf) -> node index
    points, depths = [], []
    for key, ks in touching.items():
        members = set(ks)
        while members:
            seed = members.pop()
            group, todo = [seed], [seed]
            while todo:
                k = todo.pop()
                for m in adjacency.get(k, ()):
                    if m in members:
                        members.remove(m)
                        group.append(m)
                        todo.append(m)
            index = len(points)
            points.append(key)
            depths.append(float(np.mean([leaves[k][3] for k in group])))
            for k in group:
                node_of[(key, k)] = index
    triangles = []
    for k, ((r0, c0, s, d), ring) in enumerate(zip(leaves, rings)):
        ids = [node_of[(key, k)] for key in ring]
        if len(ids) == 4:
            triangles += [(ids[0], ids[1], ids[2]), (ids[0], ids[2], ids[3])]
        else:
            centre = len(points)
            points.append((r0 + s / 2.0, c0 + s / 2.0))
            depths.append(d)
            triangles += [(centre, ids[j], ids[(j + 1) % len(ids)]) for j in range(len(ids))]
    rc = np.array(points, dtype=float)
    return Mesh(lon=lon0 + rc[:, 1] * lon_step, lat=lat0 + rc[:, 0] * step, depth=np.array(depths),
                triangles=np.array(triangles, dtype=np.int64))
