"""Flat pattern pieces as triangle meshes: every triangle keeps its own corner coordinates in the cut fabric.

A pattern piece is cut flat, so its rest state is the plane. Each triangle stores the pattern (u, v) coordinates of its
three corners in metres, which lets a piece close on itself (a sleeve sewn into a tube) without a seam in the rest
geometry. The fabric's grain is a direction in the pattern plane: the warp, with the weft at right angles to it.
"""
from __future__ import annotations

import math
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
    corner_rot: np.ndarray | None = None         # (T, 3) turn of each corner's pattern position from its vertex's
                                                 # vertex_uv, for seams that turn the pattern (None: no such seam)

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


def gore(r_in: float, r_out: float, angle: float, h: float, warp=None) -> Piece:
    """An open ring sector of `angle` radians: a flared panel (a gore, or a petal left unsewn) whose straight edges
    stay free. Rows of points about h apart from r_in to r_out run from one edge to the other, neighbouring rows
    joined by strips. The warp defaults to the panel's centre line; a bias cut turns it by 45 degrees. Rings: 'inner'
    (the top edge) and 'outer' (the hem), each in order of angle."""
    if warp is None:
        warp = (math.cos(angle / 2), math.sin(angle / 2))
    nr = max(1, int(round((r_out - r_in) / (h * math.sqrt(3) / 2))))
    radii = np.linspace(r_in, r_out, nr + 1)
    rows, off = [], 0
    for r in radii:
        n = max(2, int(round(angle * r / h)) + 1)
        rows.append((off, n, np.linspace(0.0, angle, n), r))
        off += n
    vuv = np.concatenate([np.c_[r * np.cos(psi), r * np.sin(psi)] for _, _, psi, r in rows])
    F = []
    for k in range(nr):
        (oa, na, pa, _), (ob, nb, pb, _) = rows[k], rows[k + 1]
        i = j = 0
        while i < na - 1 or j < nb - 1:
            if j >= nb - 1 or (i < na - 1 and pa[i + 1] <= pb[j + 1]):
                F.append([oa + i, oa + i + 1, ob + j])
                i += 1
            else:
                F.append([oa + i, ob + j + 1, ob + j])
                j += 1
    F = np.array(F, np.int64)
    UV = vuv[F]
    e1 = UV[:, 1] - UV[:, 0]
    e2 = UV[:, 2] - UV[:, 0]
    flip = (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) < 0
    F[flip] = F[flip][:, [0, 2, 1]]
    UV = vuv[F]
    rings = dict(inner=list(range(rows[0][0], rows[0][0] + rows[0][1])),
                 outer=list(range(rows[-1][0], rows[-1][0] + rows[-1][1])))
    return Piece(F, UV, off, warp=tuple(warp), rings=rings, vertex_uv=vuv)


def sector(r_in: float, r_out: float, angle: float, h: float, warp=None, hem=None) -> Piece:
    """A ring sector of `angle` radians whose two straight edges are sewn together: a skirt cut from part of a
    circle (angle pi for a half circle, 2 pi for a full one) or from more than one (4 pi for a double circle, whose
    pattern covers the plane twice). Rows of points about h apart from r_in to r_out are spaced evenly in angle, each
    row's points offset by half a step from the row before, and neighbouring rows are joined by a strip that wraps
    across the seam. A corner reached across the seam lies at its vertex's pattern position turned by `angle`
    (corner_rot), so every triangle keeps its own flat shape and the sewn piece rests as the cone it makes. The warp
    defaults to the straight grain at the middle of the arc, which is the centre front when the seam is at the
    centre back. Rings: 'inner' (the waist) and 'outer' (the hem), each in order of angle.

    hem: an optional function of the fraction u (0 to 1) along the arc that gives the hem's distance from r_in there,
    for a hem that is not an arc (a hem cut on a slant). Each row then lies at a fixed fraction of its radial line, and
    r_out (the longest line) sets the number of rows."""
    if warp is None:
        warp = (math.cos(angle / 2), math.sin(angle / 2))
    nr = max(1, int(round((r_out - r_in) / (h * math.sqrt(3) / 2))))
    radii = np.linspace(r_in, r_out, nr + 1)
    rows, off = [], 0
    for k, r in enumerate(radii):
        n = max(6, int(round(angle * r / h)))
        psi = (np.arange(n) + 0.5 * (k % 2)) * angle / n
        rr = np.full(n, r) if hem is None else r_in + (k / nr) * np.asarray(hem(psi / angle), float)
        rows.append((off, n, psi, rr))
        off += n
    V = off
    vuv = np.concatenate([np.c_[rr * np.cos(psi), rr * np.sin(psi)] for _, _, psi, rr in rows])
    F, UV, ROT = [], [], []

    def corner(row, l):
        o, n, psi, rr = row
        lap = l // n                                   # times this corner has been carried across the seam
        a = psi[l % n] + angle * lap
        r = rr[l % n]
        return o + l % n, (r * math.cos(a), r * math.sin(a)), angle * lap

    for k in range(nr):
        ra, rb = rows[k], rows[k + 1]
        na, nb = ra[1], rb[1]
        i = j = 0
        while i < na or j < nb:
            next_a = ra[2][(i + 1) % na] + angle * ((i + 1) // na)
            next_b = rb[2][(j + 1) % nb] + angle * ((j + 1) // nb)
            if j >= nb or (i < na and next_a <= next_b):
                cs = [corner(ra, i), corner(ra, i + 1), corner(rb, j)]
                i += 1
            else:
                cs = [corner(ra, i), corner(rb, j + 1), corner(rb, j)]
                j += 1
            F.append([c[0] for c in cs])
            UV.append([c[1] for c in cs])
            ROT.append([c[2] for c in cs])
    F = np.array(F, np.int64)
    UV = np.array(UV, float)
    ROT = np.array(ROT, float)
    e1 = UV[:, 1] - UV[:, 0]
    e2 = UV[:, 2] - UV[:, 0]
    flip = (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) < 0
    F[flip] = F[flip][:, [0, 2, 1]]
    UV[flip] = UV[flip][:, [0, 2, 1]]
    ROT[flip] = ROT[flip][:, [0, 2, 1]]
    rings = dict(inner=list(range(rows[0][0], rows[0][0] + rows[0][1])),
                 outer=list(range(rows[-1][0], rows[-1][0] + rows[-1][1])))
    return Piece(F, UV, V, warp=tuple(warp), rings=rings, vertex_uv=vuv, corner_rot=ROT)


def grid(P, warp=(1.0, 0.0), wrap: bool = False) -> Piece:
    """A pattern piece on a structured grid of flat positions P (ni, nj, 2), i along u and j along v, each cell split
    along alternating diagonals as in rectangle(). With wrap=True the last column is the first one again across a
    seam (its flat positions those of column 0 carried round the tube), so the piece is a tube of ni - 1 columns.
    Rings: 'start' (j = 0), 'end' (j = nj - 1), and for an open piece 'left' (i = 0) and 'right' (i = ni - 1)."""
    P = np.asarray(P, float)
    ni, nj = P.shape[:2]
    cols = ni - 1 if wrap else ni
    vid = lambda i, j: (i % cols) + cols * j
    F, UV = [], []
    for j in range(nj - 1):
        for i in range(ni - 1):
            c00, c10, c01, c11 = (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)
            tris = [(c00, c10, c11), (c00, c11, c01)] if (i + j) % 2 == 0 else [(c00, c10, c01), (c10, c11, c01)]
            for t in tris:
                F.append([vid(*c) for c in t])
                UV.append([P[c] for c in t])
    F, UV = _orient(np.array(F, np.int64), np.array(UV, float))
    V = cols * nj
    vuv = np.array([P[k % cols, k // cols] for k in range(V)])
    rings = dict(start=[vid(i, 0) for i in range(cols)], end=[vid(i, nj - 1) for i in range(cols)])
    if not wrap:
        rings['left'] = [vid(0, j) for j in range(nj)]
        rings['right'] = [vid(ni - 1, j) for j in range(nj)]
    return Piece(F, UV, V, warp=tuple(warp), rings=rings, vertex_uv=vuv)


def develop(X, common_rows: bool = False) -> np.ndarray:
    """Flat positions (ni, nj, 2) for a structured grid of points X (ni, nj, 3) that keep the lengths of its columns,
    each column i laid along v by its own arc length, and of its rows along u: each row j by its own arc length,
    centred, or with common_rows every row by the longest, step by step. A piece cut so starts on the grid without
    strain along its columns. Where rows of one surface differ in length (a band gathered into folds round a body,
    whose crests run further round than its hollows) rows of their own lengths slant the columns between them, so such
    a piece is cut to its longest rows and its shorter rows start compressed, which cloth takes by buckling."""
    X = np.asarray(X, float)
    du = np.linalg.norm(np.diff(X, axis=0), axis=2)
    if common_rows:
        du = np.broadcast_to(du.max(1, keepdims=True), du.shape)
    u = np.concatenate([np.zeros((1, X.shape[1])), np.cumsum(du, 0)])
    v = np.concatenate([np.zeros((X.shape[0], 1)), np.cumsum(np.linalg.norm(np.diff(X, axis=1), axis=2), 1)], 1)
    return np.stack([u - u[-1:] / 2, v], -1)


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
