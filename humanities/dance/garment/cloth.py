"""XPBD cloth with quasi-steady aerodynamics, body capsules, the floor and self-contact, on the GPU in float32.

Each substep is captured once as a CUDA graph and replayed; the self-contact pair search runs outside the graph and
writes into fixed-capacity buffers.

Method. Extended position-based dynamics (Macklin, Müller & Chentanez 2016) in its small-steps form (Macklin et al.
2019): many substeps per frame, one constraint pass each, Jacobi-parallel over the mesh with mass splitting (each
vertex's inverse mass is multiplied by the number of constraints that touch it in the denominators only), so a stiff
constraint converges like an average and a soft one acts with its true stiffness. XPBD was chosen over projective
dynamics because its per-substep work is a handful of gathers and scatters with no global solve, which suits a
GPU and float32, and because substepping with one iteration converges to the implicit solution with little
artificial damping.

Fabric. Woven cloth is close to inextensible along warp and weft and shears easily. Each triangle carries three
strain constraints in the grain directions of the flat pattern (strain-based dynamics, Müller et al. 2014): warp
stretch |F d_w| - 1, weft stretch |F d_f| - 1 and the shear cos(angle between the deformed threads), with compliances
1 / (k A). Bending is the isometric quadratic energy (Bergou et al. 2006; Wardetzky et al. 2007), written per vertex
as the cotangent Laplacian of the rest pattern: E = sum_i B / (2 A_i) |(M x)_i|^2, which equals B/2 * integral of
(k1 + k2)^2 dA for a bent sheet. At boundary vertices the rest pattern's own in-plane Laplacian is carried along by
the local deformation gradient, so M x vanishes for every isometric placement of the flat piece, edges included. M is
constant, so the bending constraint is linear with a constant gradient. Fabric damping (bending and shear
hysteresis) is not modelled; the optional Rayleigh damping of bending is off by default, because with one pass per
substep its (1 + gamma) denominator softens the static stiffness (by 30-45% in the fine cantilever tests).

Air. Each triangle feels pressure along its normal, -1/2 rho C_N |v_n| v_n A n, and skin friction along its plane,
-1/2 rho C_T |v_t| v_t A, relative to still air, applied semi-implicitly (the velocity factor k dt / (1 + k dt)) so
that very light fabric stays stable. Stacked layers in contact (a folded skirt, a flattened sleeve) move together and
present one face to the air: a vertex with another layer within 1.5 contact distances and a parallel or opposed normal
(|cos| > 0.7) takes half the force.

Contact. Vertices are kept outside the body's capsules (moved with the motion each substep) by the contact offset,
with Coulomb friction relative to the body's surface; above the floor; and apart from one another by the
self-contact distance, found with a sorted spatial hash rebuilt every few substeps. Contacts are inelastic: the
approach made within a substep is undone at once and enters the velocity; older penetration (cloth caught where two
capsules overlap) is undone at the depenetration speed in position only. Pinned vertices follow their attachment
targets.

Precision. Positions and velocities are float32. The velocity is updated from the projections' increments,
v + D / dt, never from position differences, whose float32 rounding (ulp / dt, 0.3-2 mm/s at 1-10 m and
dt = 1/2400 s) exceeds lunar gravity's gain per substep (0.7 mm/s) and froze a falling lunar sheet in testing.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, asdict, field

import numpy as np
import torch

from .patterns import Piece, mesh_edges


@dataclass(frozen=True)
class World:
    name: str
    gravity: float          # m/s^2, downward
    air_density: float      # kg/m^3


MOON = World('moon', 1.6242, 1.40)
EARTH = World('earth', 9.80665, 1.225)


@dataclass(frozen=True)
class Material:
    name: str
    areal_density: float        # kg/m^2
    bending_rigidity: float     # N m: bending moment per unit width per unit curvature (KES B, Peirce's G)
    shear_stiffness: float      # N/m per radian of in-plane shear
    tensile_stiffness: float    # N/m along warp and weft
    friction: float = 0.3       # Coulomb coefficient against the body and the floor
    normal_factor: float = 1.0  # share of a closed sheet's normal pressure an open weave or net takes (1: closed)


@dataclass
class Settings:
    dt: float = 1.0 / 2400.0
    c_normal: float = 1.17          # pressure coefficient for flow normal to a plate
    c_tangent: float = 0.01         # skin friction, both faces together, on the one-sided area
    shielding: bool = True
    bend_damping_time: float = 0.0     # s, Rayleigh damping of bending (0: air is the only damping)
    stretch_iterations: int = 2
    body_offset: float = 0.004      # m, vertices kept this far outside the capsules
    floor_offset: float = 0.002
    self_collision: bool = True
    self_distance: float = 0.0      # m; 0 means 0.75 x the mean edge length
    self_rebuild_every: int = 4
    self_max_per_cell: int = 12
    depenetration_speed: float = 0.5   # m/s: rate at which older penetration is undone (positions only)

    def to_dict(self):
        return asdict(self)


# ---------------------------------------------------------------------------------------------------------- operators

def triangle_rest(piece: Piece):
    """Rest areas, strain coefficient rows for warp and weft (T, 3) each, and Dm^-1."""
    UV = piece.UV
    Dm = np.stack([UV[:, 1] - UV[:, 0], UV[:, 2] - UV[:, 0]], axis=2)       # (T, 2, 2), columns are edges
    det = Dm[:, 0, 0] * Dm[:, 1, 1] - Dm[:, 0, 1] * Dm[:, 1, 0]
    area = 0.5 * np.abs(det)
    Dinv = np.linalg.inv(Dm)
    w = np.asarray(piece.warp, float)
    w = w / np.linalg.norm(w)
    f = np.array([-w[1], w[0]])
    a = Dinv @ w
    b = Dinv @ f
    c1 = np.stack([-a[:, 0] - a[:, 1], a[:, 0], a[:, 1]], 1)
    c2 = np.stack([-b[:, 0] - b[:, 1], b[:, 0], b[:, 1]], 1)
    return area, c1, c2, Dinv


def bending_operator(piece: Piece):
    """M = L - P (scipy CSR, V x V) and the barycentric vertex areas. L is the cotangent Laplacian of the flat
    pattern; P carries each vertex's in-plane rest Laplacian through the angle-weighted deformation gradient of its
    triangles, so M x = 0 for every isometric placement of the flat piece. A piece sewn across a seam that turns the
    pattern (a skirt cut from part of a circle) gives each corner the turn from its vertex's own pattern frame to its
    triangle's (piece.corner_rot), so each vertex's rest Laplacian is summed in one frame and read in each triangle's."""
    import scipy.sparse as sp
    F, UV, V = piece.F, piece.UV, piece.n_vertices
    area, _, _, Dinv = triangle_rest(piece)
    rot = getattr(piece, 'corner_rot', None)
    turn = (lambda vec, a: vec) if rot is None else _turn
    rows, cols, vals = [], [], []
    r = np.zeros((V, 2))
    ang = np.zeros((len(F), 3))
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        u = UV[:, j] - UV[:, i]
        v = UV[:, k] - UV[:, i]
        cr = np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
        dot = (u * v).sum(1)
        ang[:, i] = np.arctan2(cr, dot)
        wgt = 0.5 * dot / cr                                   # half the cotangent at corner i, opposite edge jk
        J, K = F[:, j], F[:, k]
        rows += [J, K, J, K]
        cols += [K, J, J, K]
        vals += [wgt, wgt, -wgt, -wgt]
        np.add.at(r, J, turn(wgt[:, None] * (UV[:, k] - UV[:, j]), None if rot is None else -rot[:, j]))
        np.add.at(r, K, turn(wgt[:, None] * (UV[:, j] - UV[:, k]), None if rot is None else -rot[:, k]))
    L = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(V, V)).tocsr()
    angsum = np.zeros(V)
    np.add.at(angsum, F, ang)
    prow, pcol, pval = [], [], []
    for q in range(3):
        vq = F[:, q]
        rq = turn(r[vq], None if rot is None else rot[:, q])   # the vertex's rest Laplacian in this triangle's frame
        d = np.einsum('tij,tj->ti', Dinv, rq)                  # Dm^-1 r for the corner's vertex
        coef = np.stack([-d[:, 0] - d[:, 1], d[:, 0], d[:, 1]], 1) * (ang[:, q] / angsum[vq])[:, None]
        for c in range(3):
            prow.append(vq)
            pcol.append(F[:, c])
            pval.append(coef[:, c])
    P = sp.coo_matrix((np.concatenate(pval), (np.concatenate(prow), np.concatenate(pcol))), shape=(V, V)).tocsr()
    M = (L - P).tocsr()
    M.eliminate_zeros()
    return M, mixed_areas(piece)


def _turn(vec, a):
    """Rotate 2-vectors (n, 2) by angles a (n,) in the pattern plane."""
    c, s = np.cos(a), np.sin(a)
    return np.stack([c * vec[:, 0] - s * vec[:, 1], s * vec[:, 0] + c * vec[:, 1]], 1)


def mixed_areas(piece: Piece) -> np.ndarray:
    """Vertex areas matched to the cotangent Laplacian: Voronoi areas for non-obtuse triangles, and half (at the obtuse
    corner) or a quarter of the triangle otherwise (Meyer, Desbrun, Schröder & Barr 2003)."""
    F, UV, V = piece.F, piece.UV, piece.n_vertices
    Av = np.zeros(V)
    e = [UV[:, 1] - UV[:, 0], UV[:, 2] - UV[:, 1], UV[:, 0] - UV[:, 2]]       # edge opposite corner 2, 0, 1
    area = 0.5 * np.abs(e[0][:, 0] * (-e[2][:, 1]) - e[0][:, 1] * (-e[2][:, 0]))
    cot = np.zeros((len(F), 3))
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        u = UV[:, j] - UV[:, i]
        v = UV[:, k] - UV[:, i]
        cot[:, i] = (u * v).sum(1) / np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
    obtuse = cot < 0
    any_obt = obtuse.any(1)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        lij = ((UV[:, j] - UV[:, i]) ** 2).sum(1)
        lik = ((UV[:, k] - UV[:, i]) ** 2).sum(1)
        vor = (lij * cot[:, k] + lik * cot[:, j]) / 8.0
        a = np.where(~any_obt, vor, np.where(obtuse[:, i], area / 2, area / 4))
        np.add.at(Av, F[:, i], a)
    return Av


def bending_energy(piece: Piece, X: np.ndarray, B: float) -> float:
    M, Av = bending_operator(piece)
    C = M @ X
    return float((B / (2 * Av) * (C * C).sum(1)).sum())


def _csr(m, dev):
    m = m.tocsr()
    return torch.sparse_csr_tensor(torch.as_tensor(m.indptr, dtype=torch.int64), torch.as_tensor(m.indices, dtype=torch.int64),
                                   torch.as_tensor(m.data, dtype=torch.float32), size=m.shape).to(dev)


def quat_rotate(q, v):
    """Rotate vectors v (..., 3) by unit quaternions q (..., 4), (w, x, y, z)."""
    w = q[..., :1]
    u = q[..., 1:]
    t = 2.0 * torch.cross(u, v, dim=-1)
    return v + w * t + torch.cross(u, t, dim=-1)


# ------------------------------------------------------------------------------------------------------------ solver

class Cloth:
    """One garment piece. `pinned` vertices follow targets given by the caller each substep."""

    def __init__(self, piece: Piece, X0: np.ndarray, material: Material, world: World, settings: Settings | None = None,
                 pinned=None, device=None, extra_mass=None, tethers=None):
        s = settings or Settings()
        self.s, self.material, self.world, self.piece = s, material, world, piece
        dev = torch.device(device or ('cuda' if torch.cuda.is_available() else 'cpu'))
        self.dev = dev
        V = piece.n_vertices
        F = np.asarray(piece.F, np.int64)
        area, c1, c2, _ = triangle_rest(piece)
        m = np.zeros(V)
        np.add.at(m, F, (material.areal_density * area / 3)[:, None])
        if extra_mass is not None:                     # mass sewn to the cloth (a weighted hem) that the air does not see
            m = m + np.asarray(extra_mass, float)
        pin = np.zeros(V, bool)
        if pinned is not None:
            pin[np.asarray(pinned, np.int64)] = True
        w = np.where(pin, 0.0, 1.0 / m)
        tri_n = np.bincount(F.ravel(), minlength=V).astype(float)
        nw = w * tri_n
        W11 = (nw[F] * c1 * c1).sum(1)
        W22 = (nw[F] * c2 * c2).sum(1)
        W12 = (nw[F] * c1 * c2).sum(1)
        M, Av = bending_operator(piece)
        Mabs2 = M.copy()
        Mabs2.data = Mabs2.data ** 2
        bend_n = np.asarray((M != 0).sum(axis=0)).ravel().astype(float)       # constraints touching each vertex
        S_bend = Mabs2 @ (w * bend_n)
        f32 = dict(device=dev, dtype=torch.float32)
        self.V, self.T = V, len(F)
        self.F = torch.as_tensor(F, device=dev)
        self.c1 = torch.as_tensor(c1, **f32)
        self.c2 = torch.as_tensor(c2, **f32)
        self.area0 = torch.as_tensor(area, **f32)
        self.mtri = torch.as_tensor(material.areal_density * area, **f32)
        self.m = torch.as_tensor(m, **f32)
        self.w = torch.as_tensor(w, **f32)
        self.pin = torch.as_tensor(pin, device=dev)
        self.free = ~self.pin
        self.W11 = torch.as_tensor(W11, **f32)
        self.W22 = torch.as_tensor(W22, **f32)
        self.W12 = torch.as_tensor(W12, **f32)
        self.M = _csr(M, dev)
        self.MT = _csr(M.T.tocsr(), dev)
        self.S_bend = torch.as_tensor(S_bend, **f32)
        self.Av = torch.as_tensor(Av, **f32)
        self.bend_active = torch.as_tensor(S_bend > 0, device=dev)
        self.x = torch.as_tensor(X0, **f32).clone()
        self.v = torch.zeros_like(self.x)
        self.edge_len = piece.edge_length()
        self.d_self = s.self_distance or 0.75 * self.edge_len
        edges = mesh_edges(F)
        self.excl = torch.as_tensor(np.sort(edges[:, 0] * V + edges[:, 1]), device=dev)
        self.pairs = None
        self.stacked = torch.zeros(V, dtype=torch.bool, device=dev)
        self.n_sub = 0
        self.gvec = torch.tensor([0.0, 0.0, -world.gravity], **f32)
        self.pin_idx = torch.nonzero(self.pin).ravel()
        self.has_pins = bool(pin.any())
        # long-range attachments (Kim, Chentanez & Mueller 2012): each tethered vertex stays within a set distance of
        # a pinned anchor, so a heavy hem's pull cannot stretch a light fabric that one stretch pass cannot stiffen
        self.has_tethers = tethers is not None
        if self.has_tethers:
            anchor, length = (np.asarray(a) for a in tethers)
            keep = (anchor >= 0) & ~pin
            self.t_idx = torch.as_tensor(np.nonzero(keep)[0], device=dev)
            self.t_anchor = torch.as_tensor(anchor[keep], dtype=torch.int64, device=dev)
            self.t_len = torch.as_tensor(length[keep], **f32)
        self.has_body = False
        # static inputs for the substep (filled by set_frame / set_fraction), so the substep can be a CUDA graph
        self.capA0 = self.capdA = self.capB0 = self.capdB = self.cap_r = None
        self.tg0 = torch.zeros((len(self.pin_idx), 3), **f32)
        self.tgd = torch.zeros_like(self.tg0)
        self.frac0 = torch.zeros((), **f32)
        self.frac1 = torch.zeros((), **f32)
        self.damp = torch.zeros((), **f32)
        cap = max(1024, 8 * V)
        self.pair_cap = cap
        self.pi = torch.zeros(cap, dtype=torch.int64, device=dev)
        self.pj = torch.zeros(cap, dtype=torch.int64, device=dev)
        self.pvalid = torch.zeros(cap, **f32)
        self.n_pairs = 0
        self.graph = None
        self.graph_dt = None
        self.use_graph = dev.type == 'cuda'
        self.stats = dict(pairs_max=0, pairs_truncated=0)

    def set_body(self, cap_r, cap_rb=None):
        """Capsule radii at their a ends, and at their b ends for tapered pieces (default: the same)."""
        n = len(cap_r)
        f32 = dict(device=self.dev, dtype=torch.float32)
        self.cap_r = torch.as_tensor(cap_r, **f32)
        self.cap_dr = torch.as_tensor(cap_rb, **f32) - self.cap_r if cap_rb is not None else torch.zeros_like(self.cap_r)
        self.capA0 = torch.zeros((n, 3), **f32)
        self.capdA = torch.zeros((n, 3), **f32)
        self.capB0 = torch.zeros((n, 3), **f32)
        self.capdB = torch.zeros((n, 3), **f32)
        self.has_body = True
        self.graph = None

    def set_frame(self, A0=None, A1=None, B0=None, B1=None, T0=None, T1=None, damping=None):
        """Pose for one frame: capsule end points and pinned targets at its start and end (device tensors)."""
        if A0 is not None:
            self.capA0.copy_(A0)
            self.capdA.copy_(A1 - A0)
            self.capB0.copy_(B0)
            self.capdB.copy_(B1 - B0)
        if T0 is not None:
            self.tg0.copy_(T0)
            self.tgd.copy_(T1 - T0)
        if damping is not None:
            self.damp.fill_(float(damping))

    def set_targets(self, T):
        """Hold pinned vertices at T (no motion within the frame)."""
        self.tg0.copy_(T)
        self.tgd.zero_()

    # ------------------------------------------------------------------------------------------------ pieces of a step
    def tri_geometry(self, x):
        F = self.F
        x0, x1, x2 = x[F[:, 0]], x[F[:, 1]], x[F[:, 2]]
        n = torch.cross(x1 - x0, x2 - x0, dim=1)
        A2 = n.norm(dim=1).clamp_min(1e-12)
        return x0, x1, x2, n / A2[:, None], 0.5 * A2

    def air(self, dt):
        rho = self.world.air_density
        if rho <= 0:
            return
        s = self.s
        F = self.F
        _, _, _, n, A = self.tri_geometry(self.x)
        vt = (self.v[F[:, 0]] + self.v[F[:, 1]] + self.v[F[:, 2]]) / 3.0
        vn = (vt * n).sum(1)
        vtan = vt - vn[:, None] * n
        vtn = vtan.norm(dim=1)
        kn = 0.5 * rho * s.c_normal * self.material.normal_factor * vn.abs() * A / self.mtri
        kt = 0.5 * rho * s.c_tangent * vtn * A / self.mtri
        dv = -(kn * dt / (1 + kn * dt) * vn)[:, None] * n - (kt * dt / (1 + kt * dt))[:, None] * vtan
        imp = (self.mtri / 3.0)[:, None] * dv
        acc = torch.zeros_like(self.v)
        for c in range(3):
            acc.index_add_(0, F[:, c], imp)
        dvv = acc / self.m[:, None]
        if s.shielding:
            dvv = torch.where(self.stacked[:, None], 0.5 * dvv, dvv)
        self.v += dvv * self.free[:, None]

    def stretch(self, x, dt):
        F = self.F
        x0, x1, x2 = x[F[:, 0]], x[F[:, 1]], x[F[:, 2]]
        c1, c2 = self.c1, self.c2
        f1 = c1[:, 0:1] * x0 + c1[:, 1:2] * x1 + c1[:, 2:3] * x2
        f2 = c2[:, 0:1] * x0 + c2[:, 1:2] * x1 + c2[:, 2:3] * x2
        l1 = f1.norm(dim=1).clamp_min(1e-9)
        l2 = f2.norm(dim=1).clamp_min(1e-9)
        g1 = f1 / l1[:, None]
        g2 = f2 / l2[:, None]
        C1 = l1 - 1.0
        C2 = l2 - 1.0
        C3 = (g1 * g2).sum(1)
        p = (g2 - C3[:, None] * g1) / l1[:, None]
        q = (g1 - C3[:, None] * g2) / l2[:, None]
        mat = self.material
        at = 1.0 / (mat.tensile_stiffness * self.area0 * dt * dt)
        ash = 1.0 / (mat.shear_stiffness * self.area0 * dt * dt)
        d1 = -C1 / (self.W11 + at)
        d2 = -C2 / (self.W22 + at)
        den3 = (p * p).sum(1) * self.W11 + 2 * (p * q).sum(1) * self.W12 + (q * q).sum(1) * self.W22
        d3 = -C3 / (den3 + ash)
        dx = torch.zeros_like(x)
        for k in range(3):
            corr = (c1[:, k:k + 1] * (g1 * d1[:, None] + p * d3[:, None]) + c2[:, k:k + 1] * (g2 * d2[:, None] + q * d3[:, None]))
            dx.index_add_(0, F[:, k], corr)
        return dx * self.w[:, None]

    def tether(self, x):
        """Pull each tethered vertex back to its tether's length from its anchor, where it has gone beyond."""
        v = x[self.t_idx] - x[self.t_anchor]
        n = v.norm(dim=1).clamp_min(1e-12)
        ex = (n - self.t_len).clamp_min(0.0)
        out = torch.zeros_like(x)
        out.index_copy_(0, self.t_idx, -(v / n[:, None]) * ex[:, None])
        return out

    def bend(self, x, x_start, dt):
        B = self.material.bending_rigidity
        if B <= 0:
            return torch.zeros_like(x)
        gam = self.s.bend_damping_time / dt
        C = self.M @ ((1 + gam) * x - gam * x_start)
        alpha = self.Av / (B * dt * dt)
        dl = -C / ((1 + gam) * self.S_bend + alpha)[:, None]
        dl = dl * self.bend_active[:, None]
        return self.w[:, None] * (self.MT @ dl)

    def build_pairs(self):
        x = self.x
        V = self.V
        rs = self.d_self + 0.5 * self.edge_len
        c = torch.floor(x / rs).to(torch.int64) + 1024
        key = (c[:, 0] * 2048 + c[:, 1]) * 2048 + c[:, 2]
        skey, order = torch.sort(key)
        offs = torch.tensor([(dx * 2048 + dy) * 2048 + dz for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)],
                            device=self.dev)
        nk = key[:, None] + offs[None, :]
        start = torch.searchsorted(skey, nk)
        end = torch.searchsorted(skey, nk, right=True)
        K = self.s.self_max_per_cell
        cnt = (end - start).clamp(max=K)
        kk = torch.arange(K, device=self.dev)
        idx = (start[:, :, None] + kk).clamp(max=V - 1)
        valid = kk[None, None, :] < cnt[:, :, None]
        j = order[idx]
        i = torch.arange(V, device=self.dev)[:, None, None].expand_as(j)
        msk = valid & (j > i)
        i, j = i[msk], j[msk]
        d = (x[i] - x[j]).norm(dim=1)
        near = d < rs
        i, j = i[near], j[near]
        pk = i * V + j
        pos = torch.searchsorted(self.excl, pk).clamp(max=len(self.excl) - 1)
        keep = self.excl[pos] != pk
        i, j = i[keep], j[keep]
        n = len(i)
        self.stats['pairs_max'] = max(self.stats['pairs_max'], n)
        if n > self.pair_cap:
            self.stats['pairs_truncated'] += 1
            i, j, n = i[:self.pair_cap], j[:self.pair_cap], self.pair_cap
        self.pi[:n] = i
        self.pj[:n] = j
        self.pi[n:] = 0
        self.pj[n:] = 0
        self.pvalid[:n] = 1.0
        self.pvalid[n:] = 0.0
        self.n_pairs = n
        self.pairs = (i, j)
        if self.s.shielding:
            _, _, _, n, A = self.tri_geometry(x)
            vn = torch.zeros_like(x)
            for k in range(3):
                vn.index_add_(0, self.F[:, k], n * A[:, None])
            vn = vn / vn.norm(dim=1, keepdim=True).clamp_min(1e-12)
            pi, pj = self.pairs
            close = ((x[pi] - x[pj]).norm(dim=1) < 1.5 * self.d_self) & ((vn[pi] * vn[pj]).sum(1).abs() > 0.7)
            st = torch.zeros(V, dtype=torch.float32, device=self.dev)
            st.index_add_(0, pi[close], torch.ones_like(pi[close], dtype=torch.float32))
            st.index_add_(0, pj[close], torch.ones_like(pj[close], dtype=torch.float32))
            self.stacked.copy_(st > 0)

    def self_contact(self, x, disp):
        """Pairs closer than the contact distance are pushed apart: the approach made in this substep is undone at
        once (an inelastic contact), older overlap at the depenetration speed in position only."""
        pi, pj = self.pi, self.pj
        d = x[pi] - x[pj]
        dist = d.norm(dim=1).clamp_min(1e-9)
        n = d / dist[:, None]
        pen = (self.d_self - dist).clamp_min(0.0) * self.pvalid
        approach = (-((disp[pi] - disp[pj]) * n).sum(1)).clamp_min(0.0)
        p_pos = torch.minimum(pen, approach + self.s.depenetration_speed * self.dt_now)
        p_vel = torch.minimum(pen, approach)
        wi, wj = self.w[pi], self.w[pj]
        ws = (wi + wj).clamp_min(1e-12)
        cnt = torch.zeros(self.V, device=self.dev)
        one = (pen > 0).to(x.dtype)
        cnt.index_add_(0, pi, one)
        cnt.index_add_(0, pj, one)
        cnt = cnt.clamp_min(1.0)[:, None]
        out = []
        for p in (p_pos, p_vel):
            corr = (p / ws)[:, None] * n
            acc = torch.zeros_like(x)
            acc.index_add_(0, pi, corr * wi[:, None])
            acc.index_add_(0, pj, -corr * wj[:, None])
            out.append(acc / cnt)
        return out[0], out[1]

    def body_contact(self, x, disp):
        """x: working positions; disp: their displacement over this substep so far (dt v + corrections)."""
        if not self.has_body:
            z = torch.zeros_like(x)
            return z, z
        a = self.capA0 + self.frac1 * self.capdA
        b = self.capB0 + self.frac1 * self.capdB
        a0 = self.capA0 + self.frac0 * self.capdA
        b0 = self.capB0 + self.frac0 * self.capdB
        r = self.cap_r
        mu = self.material.friction
        ab = b - a
        L2 = (ab * ab).sum(1).clamp_min(1e-12)
        xa = x[:, None, :] - a[None]
        t = ((xa * ab[None]).sum(2) / L2[None]).clamp(0.0, 1.0)
        cpt = a[None] + t[..., None] * ab[None]
        d = x[:, None, :] - cpt
        dist = d.norm(dim=2)
        pen = r[None] + t * self.cap_dr[None] + self.s.body_offset - dist
        pmax, best = pen.max(dim=1)
        act = (pmax > 0) & self.free
        ar = torch.arange(self.V, device=self.dev)
        db = d[ar, best]
        nb = db / dist[ar, best].clamp_min(1e-9)[:, None]
        tb = t[ar, best]
        c_now = cpt[ar, best]
        c_was = a0[best] + tb[:, None] * (b0[best] - a0[best])
        # resolve the penetration made in this substep at once, and older penetration (a vertex caught where two
        # capsules overlap, as an arm resting on the hip does) at the depenetration speed only
        approach = (((c_now - c_was) * nb).sum(1) - (disp * nb).sum(1)).clamp_min(0)
        depth = torch.minimum(pmax.clamp_min(0), approach + self.s.depenetration_speed * self.dt_now) * act
        depth_v = torch.minimum(pmax.clamp_min(0), approach) * act
        push = nb * depth[:, None]
        # Coulomb friction relative to the moving surface (the axis point's motion over the substep)
        rel = (disp + push) - (c_now - c_was)
        rel_t = rel - (rel * nb).sum(1, keepdim=True) * nb
        rt = rel_t.norm(dim=1).clamp_min(1e-12)
        k = torch.clamp(mu * depth / rt, max=1.0) * act
        fr = -rel_t * k[:, None]
        return push + fr, nb * depth_v[:, None] + fr

    def floor_contact(self, x, disp):
        z0 = self.s.floor_offset
        pen0 = (z0 - x[:, 2]).clamp_min(0) * self.free
        approach = (-disp[:, 2]).clamp_min(0)
        pen = torch.minimum(pen0, approach + self.s.depenetration_speed * self.dt_now)
        rel = disp[:, :2]
        rt = rel.norm(dim=1).clamp_min(1e-12)
        k = torch.clamp(self.material.friction * pen / rt, max=1.0)
        d = torch.cat([-rel * k[:, None], pen[:, None]], dim=1)
        dv = torch.cat([-rel * k[:, None], torch.minimum(pen0, approach)[:, None]], dim=1)
        return d, dv

    # -------------------------------------------------------------------------------------------------------- the step
    dt_now = 1.0 / 2400.0

    def _step(self, dt):
        """One substep on static tensors (capturable as a CUDA graph). Constraint projections return increments D;
        the velocity is updated as v + D / dt, so that it never inherits the rounding of float32 positions (a position
        difference at 1-10 m carries an error of ulp / dt, 0.3-2 mm/s at dt = 1/2400 s, larger than lunar gravity's
        increment g dt = 0.7 mm/s)."""
        s = self.s
        self.dt_now = dt
        self.air(dt)
        v = (self.v + dt * self.gvec * self.free[:, None]) / (1.0 + self.damp * dt)
        x0 = self.x
        step = dt * v
        D = torch.zeros_like(x0)
        y = x0 + step
        tg = self.tg0 + self.frac1 * self.tgd
        if self.has_pins:
            y[self.pin_idx] = tg
        for _ in range(s.stretch_iterations):
            d = self.stretch(y, dt)
            D = D + d
            y = y + d
        if self.has_tethers:
            d = self.tether(y)
            D = D + d
            y = y + d
        d = self.bend(y, x0, dt)
        D = D + d
        y = y + d
        Dv = D.clone()
        if s.self_collision:
            d, dv = self.self_contact(y, step + D)
            D = D + d
            Dv = Dv + dv
            y = y + d
        d, dv = self.body_contact(y, step + D)
        D = D + d
        Dv = Dv + dv
        y = y + d
        d, dv = self.floor_contact(y, step + D)
        D = D + d
        Dv = Dv + dv
        y = y + d
        v = v + Dv / dt
        x = y
        if self.has_pins:
            v[self.pin_idx] = (tg - x0[self.pin_idx]) / dt
            x[self.pin_idx] = tg
        self.x.copy_(x)
        self.v.copy_(v)

    def substep(self, dt, s0=0.0, s1=1.0):
        """Advance one substep of length dt; s0, s1 locate it within the current frame (for body and pins)."""
        s = self.s
        if s.self_collision and self.n_sub % s.self_rebuild_every == 0:
            self.build_pairs()
        self.frac0.fill_(s0)
        self.frac1.fill_(s1)
        if self.use_graph:
            if self.graph is None or self.graph_dt != dt:
                self._capture(dt)
            self.graph.replay()
        else:
            self._step(dt)
        self.n_sub += 1

    def _capture(self, dt):
        x, v = self.x.clone(), self.v.clone()
        st = torch.cuda.Stream()
        st.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(st):
            for _ in range(2):
                self._step(dt)
        torch.cuda.current_stream().wait_stream(st)
        g = torch.cuda.CUDAGraph()
        with torch.cuda.graph(g):
            self._step(dt)
        self.x.copy_(x)
        self.v.copy_(v)
        self.graph, self.graph_dt = g, dt

    # ---------------------------------------------------------------------------------------------------- diagnostics
    def strain_stats(self):
        F = self.F
        x = self.x
        x0, x1, x2 = x[F[:, 0]], x[F[:, 1]], x[F[:, 2]]
        f1 = self.c1[:, 0:1] * x0 + self.c1[:, 1:2] * x1 + self.c1[:, 2:3] * x2
        f2 = self.c2[:, 0:1] * x0 + self.c2[:, 1:2] * x1 + self.c2[:, 2:3] * x2
        e = torch.cat([f1.norm(dim=1) - 1, f2.norm(dim=1) - 1])
        return dict(max_stretch=float(e.max()), mean_abs_strain=float(e.abs().mean()),
                    p99_stretch=float(torch.quantile(e, 0.99)))

    def energy(self):
        """Kinetic, gravitational and bending energy (J)."""
        ke = 0.5 * float((self.m * (self.v * self.v).sum(1) * self.free).sum())
        pe = float((self.m * self.x[:, 2]).sum()) * self.world.gravity
        C = self.M @ self.x
        be = float((self.material.bending_rigidity / (2 * self.Av) * (C * C).sum(1)).sum())
        return dict(kinetic=ke, potential=pe, bending=be)
