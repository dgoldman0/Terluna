"""Flat pattern pieces as triangle meshes: every triangle keeps its own corner coordinates in the cut fabric.

A pattern piece is cut flat, so its rest state is the plane. Each triangle stores the pattern (u, v) coordinates of its
three corners in metres, which lets a piece close on itself (a sleeve sewn into a tube) without a seam in the rest
geometry. The fabric's grain is a direction in the pattern plane: the warp, with the weft at right angles to it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Piece:
    """A pattern piece: V vertices, T triangles, corner coordinates (T, 3, 2) in the flat fabric, and a grain."""
    F: np.ndarray                    # (T, 3) int64 vertex indices, counter-clockwise in the pattern
    UV: np.ndarray                   # (T, 3, 2) pattern coordinates of each corner, metres
    n_vertices: int
    warp: tuple = (1.0, 0.0)         # warp direction in the pattern plane
    rings: dict = field(default_factory=dict)   # named vertex index lists (waist, hem, cuff, end, ...)
    vertex_uv: np.ndarray | None = None          # (V, 2) one pattern position per vertex (seam side arbitrary)

    def edge_length(self) -> float:
        e = np.concatenate([self.UV[:, 1] - self.UV[:, 0], self.UV[:, 2] - self.UV[:, 1], self.UV[:, 0] - self.UV[:, 2]])
        return float(np.linalg.norm(e, axis=1).mean())

    def area(self) -> float:
        e1 = self.UV[:, 1] - self.UV[:, 0]
        e2 = self.UV[:, 2] - self.UV[:, 0]
        return float(0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]).sum())


def _orient(F, UV):
    e1 = UV[:, 1] - UV[:, 0]
    e2 = UV[:, 2] - UV[:, 0]
    flip = (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) < 0
    F[flip] = F[flip][:, [0, 2, 1]]
    UV[flip] = UV[flip][:, [0, 2, 1]]
    return F, UV


def rectangle(width: float, length: float, h: float, wrap: bool = False, warp=(0.0, 1.0)) -> Piece:
    """A width x length rectangle (u across, v along) on a square grid of spacing about h, each cell split along
    alternating diagonals. With wrap=True the u edges are sewn together into a tube of circumference `width`.
    Rings: 'start' (v = 0) and 'end' (v = length); for an open strip also 'left' (u = 0) and 'right'."""
    nu = max(2, int(round(width / h)))
    if wrap and nu % 2:
        nu += 1
    nv = max(1, int(round(length / h)))
    du, dv = width / nu, length / nv
    cols = nu if wrap else nu + 1
    vid = lambda i, j: (i % nu if wrap else i) + cols * j
    F, UV = [], []
    for j in range(nv):
        for i in range(nu):
            c00, c10, c01, c11 = (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)
            tris = [(c00, c10, c11), (c00, c11, c01)] if (i + j) % 2 == 0 else [(c00, c10, c01), (c10, c11, c01)]
            for t in tris:
                F.append([vid(*c) for c in t])
                UV.append([(c[0] * du, c[1] * dv) for c in t])
    F = np.array(F, np.int64)
    UV = np.array(UV, float)
    F, UV = _orient(F, UV)
    V = cols * (nv + 1)
    vuv = np.array([((k % cols) * du, (k // cols) * dv) for k in range(V)])
    rings = dict(start=[vid(i, 0) for i in range(cols)], end=[vid(i, nv) for i in range(cols)])
    if not wrap:
        rings['left'] = [vid(0, j) for j in range(nv + 1)]
        rings['right'] = [vid(nu, j) for j in range(nv + 1)]
    return Piece(F, UV, V, warp=tuple(warp), rings=rings, vertex_uv=vuv)


def annulus(r_in: float, r_out: float, h: float, seed: int = 0, warp=(1.0, 0.0)) -> Piece:
    """A flat ring (a full circle skirt's pattern): concentric rows of points about h apart, each row offset by half a
    step, joined by a Delaunay triangulation; the hole is removed. Rings: 'inner' (the waist), 'outer' (the hem)."""
    from scipy.spatial import Delaunay
    nr = max(1, int(round((r_out - r_in) / (h * np.sqrt(3) / 2))))
    radii = np.linspace(r_in, r_out, nr + 1)
    rng = np.random.default_rng(seed)
    pts, ring_of = [], []
    for k, r in enumerate(radii):
        n = max(6, int(round(2 * np.pi * r / h)))
        phase = (0.5 * (k % 2) + 0.1 * rng.uniform(-1, 1)) * 2 * np.pi / n
        a = phase + 2 * np.pi * np.arange(n) / n
        pts.append(np.c_[r * np.cos(a), r * np.sin(a)])
        ring_of += [k] * n
    P = np.concatenate(pts)
    ring_of = np.array(ring_of)
    tri = Delaunay(P).simplices.astype(np.int64)
    cen = P[tri].mean(axis=1)
    keep = np.hypot(cen[:, 0], cen[:, 1]) > r_in
    # also drop slivers spanning the hole between two inner-ring points
    tri = tri[keep]
    UV = P[tri]
    tri, UV = _orient(tri, UV)
    # drop degenerate triangles (area far below the mean)
    e1 = UV[:, 1] - UV[:, 0]
    e2 = UV[:, 2] - UV[:, 0]
    ar = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    good = ar > 1e-3 * ar.mean()
    tri, UV = tri[good], UV[good]
    rings = dict(inner=list(np.where(ring_of == 0)[0]), outer=list(np.where(ring_of == nr)[0]))
    return Piece(tri, UV, len(P), warp=tuple(warp), rings=rings, vertex_uv=P)


def mesh_edges(F: np.ndarray) -> np.ndarray:
    """Unique undirected edges (E, 2), i < j."""
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.sort(e, axis=1)
    return np.unique(e, axis=0)


def boundary_vertices(F: np.ndarray) -> np.ndarray:
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.sort(e, axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return np.unique(u[c == 1].ravel())
