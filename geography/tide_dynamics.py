"""How closely each sea follows its equilibrium tide, and its seiche periods.

The equilibrium tide assumes the water settles at once. A sea settles in the
time a long wave takes to cross it, sqrt(g h) for depth h, which at lunar
gravity is hours to days; the tide changes over weeks. This module solves the
linear shallow-water equations on each sea's 1-degree cells for a single
frequency of forcing, with the Moon's slow rotation (f = 2 Omega sin lat) and
a weak linear bed friction, and compares the result with the equilibrium
response. It also finds each sea's slowest free oscillations (seiches).

Discretisation: an Arakawa C grid on the atlas's latitude-longitude cells.
Levels sit at cell centres and flows at the faces between neighbouring water
cells; faces against land carry no flow. Face depth is the shallower of the
two cells, at least 1 m. The Coriolis term averages the four nearest
perpendicular flows. The forced problem is solved directly in complex form;
the seiches come from the generalised eigenproblem of the same operator
without rotation or friction.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, SIDEREAL_MONTH_DAYS

FRICTION_M_S = 2.5e-5      # linear bed friction r = Cd u, with Cd 0.0025 and tidal flow of 1 cm/s
MIN_DEPTH_M = 1.0
OMEGA = 2 * math.pi / (SIDEREAL_MONTH_DAYS * 86400)


def largest_component(mask):
    """The largest region of `mask` joined through shared cell faces, wrapping in longitude."""
    from scipy.sparse.csgraph import connected_components
    rows, cols = np.nonzero(mask)
    index = -np.ones(mask.shape, dtype=int)
    index[rows, cols] = np.arange(len(rows))
    east = index[rows, (cols + 1) % mask.shape[1]]
    north = np.where(rows > 0, index[np.maximum(rows - 1, 0), cols], -1)
    a = np.r_[np.arange(len(rows))[east >= 0], np.arange(len(rows))[north >= 0]]
    b = np.r_[east[east >= 0], north[north >= 0]]
    graph = sparse.csr_matrix((np.ones(len(a)), (a, b)), shape=(len(rows), len(rows)))
    count, label = connected_components(graph, directed=False)
    keep = np.bincount(label).argmax()
    out = np.zeros_like(mask)
    out[rows[label == keep], cols[label == keep]] = True
    return out


class SeaGrid:
    """One sea's cells and faces on a global latitude-longitude grid (latitudes descending)."""

    def __init__(self, labels, depth, lat_deg, lon_deg, body, *, main_body=True):
        mask = labels == body
        self.body_cells = int(mask.sum())
        if main_body:
            mask = largest_component(mask)
        self.mask = mask
        ny, nx = mask.shape
        self.rows, self.cols = np.nonzero(mask)
        index = -np.ones(mask.shape, dtype=int)
        index[self.rows, self.cols] = np.arange(len(self.rows))
        self.index = index
        dlat = math.radians(abs(lat_deg[1] - lat_deg[0]))
        dlon = math.radians(abs(lon_deg[1] - lon_deg[0]))
        lat = np.radians(lat_deg)
        r = MOON_RADIUS
        self.lat = lat[self.rows]
        self.lon = np.radians(lon_deg)[self.cols]
        self.depth = np.maximum(depth[self.rows, self.cols], MIN_DEPTH_M)
        self.area = r * r * dlon * (np.sin(self.lat + dlat / 2) - np.sin(self.lat - dlat / 2))
        # East faces: (j, i) to (j, i+1), wrapping in longitude.
        east = index[self.rows, (self.cols + 1) % nx]
        has = east >= 0
        self.u_a, self.u_b = np.arange(len(self.rows))[has], east[has]
        self.u_len = np.full(has.sum(), r * dlat)
        self.u_dist = r * np.cos(self.lat[has]) * dlon
        self.u_lat = self.lat[has]
        # North faces: (j, i) to (j-1, i), the row above in a descending-latitude grid.
        north = np.where(self.rows > 0, index[np.maximum(self.rows - 1, 0), self.cols], -1)
        has = north >= 0
        self.v_a, self.v_b = np.arange(len(self.rows))[has], north[has]
        self.v_lat = self.lat[has] + dlat / 2
        self.v_len = r * np.cos(self.v_lat) * dlon
        self.v_dist = np.full(has.sum(), r * dlat)
        self.u_depth = np.maximum(np.minimum(self.depth[self.u_a], self.depth[self.u_b]), MIN_DEPTH_M)
        self.v_depth = np.maximum(np.minimum(self.depth[self.v_a], self.depth[self.v_b]), MIN_DEPTH_M)

    @property
    def cells(self):
        return len(self.rows)

    def gradient(self):
        """Face differences (b - a) of a cell field, over the face distance: (faces, cells)."""
        n = self.cells
        def block(a, b, dist):
            k = len(a)
            data = np.r_[1 / dist, -1 / dist]
            return sparse.csr_matrix((data, (np.r_[np.arange(k), np.arange(k)], np.r_[b, a])), shape=(k, n))
        return block(self.u_a, self.u_b, self.u_dist), block(self.v_a, self.v_b, self.v_dist)

    def divergence(self):
        """Net outflow of face transports (H u times face length), per cell: (cells, faces)."""
        n = self.cells
        def block(a, b, depth, length):
            k = len(a)
            data = np.r_[depth * length, -depth * length]
            return sparse.csr_matrix((data, (np.r_[a, b], np.r_[np.arange(k), np.arange(k)])), shape=(n, k))
        return block(self.u_a, self.u_b, self.u_depth, self.u_len), block(self.v_a, self.v_b, self.v_depth, self.v_len)

    def coriolis_averages(self):
        """Four-point averages of v at u faces and of u at v faces."""
        n = self.cells
        def faces_of(cells_a, cells_b, count):
            # Each cell's list of adjacent faces of the other family, as a (cells, faces) incidence.
            rows = np.r_[cells_a, cells_b]
            cols = np.r_[np.arange(count), np.arange(count)]
            return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, count))
        v_of_cell = faces_of(self.v_a, self.v_b, len(self.v_a))
        u_of_cell = faces_of(self.u_a, self.u_b, len(self.u_a))
        def pick(a, b, k):
            return sparse.csr_matrix((np.ones(2 * k), (np.r_[np.arange(k), np.arange(k)], np.r_[a, b])), shape=(k, n))
        v_at_u = 0.25 * (pick(self.u_a, self.u_b, len(self.u_a)) @ v_of_cell)
        u_at_v = 0.25 * (pick(self.v_a, self.v_b, len(self.v_a)) @ u_of_cell)
        # A face shared by both cells is counted once per cell; keep the classic four points.
        v_at_u.data = np.minimum(v_at_u.data, 0.25)
        u_at_v.data = np.minimum(u_at_v.data, 0.25)
        return v_at_u, u_at_v

    def equilibrium(self, forcing):
        return forcing - np.sum(forcing * self.area) / np.sum(self.area)


def forced_response(sea, forcing, frequency_cpd, *, gravity=MOON_SURFACE_GRAVITY, friction=FRICTION_M_S, rotation=True):
    """Complex level at each cell for forcing height `forcing` (complex, per cell) at one frequency."""
    omega = 2 * math.pi * frequency_cpd / 86400
    gu, gv = sea.gradient()
    du, dv = sea.divergence()
    nu, nv, n = len(sea.u_a), len(sea.v_a), sea.cells
    fu = 2 * OMEGA * np.sin(sea.u_lat) if rotation else np.zeros(nu)
    fv = 2 * OMEGA * np.sin(sea.v_lat) if rotation else np.zeros(nv)
    v_at_u, u_at_v = sea.coriolis_averages()
    iw = 1j * omega
    blocks = [[sparse.diags(iw * sea.area), du, dv],
              [gravity * gu, sparse.diags(iw + friction / sea.u_depth), -sparse.diags(fu) @ v_at_u],
              [gravity * gv, sparse.diags(fv) @ u_at_v, sparse.diags(iw + friction / sea.v_depth)]]
    matrix = sparse.bmat(blocks, format="csc").astype(complex)
    rhs = np.r_[np.zeros(n), gravity * (gu @ forcing), gravity * (gv @ forcing)].astype(complex)
    solution = linalg.spsolve(matrix, rhs)
    return solution[:n]


def seiche_periods(sea, count=5, *, gravity=MOON_SURFACE_GRAVITY):
    """The slowest free oscillations of the sea, without rotation or friction: periods in hours."""
    gu, gv = sea.gradient()
    # Energy-consistent stiffness: K = G^T diag(g H L d) G, with mass A (cell areas).
    wu = gravity * sea.u_depth * sea.u_len * sea.u_dist
    wv = gravity * sea.v_depth * sea.v_len * sea.v_dist
    stiffness = (gu.T @ sparse.diags(wu) @ gu + gv.T @ sparse.diags(wv) @ gv).tocsc()
    mass = sparse.diags(sea.area).tocsc()
    shift = -1e-14 * float(stiffness.diagonal().max() / sea.area.min())
    values, vectors = linalg.eigsh(stiffness, k=count + 1, M=mass, sigma=shift, which="LM")
    order = np.argsort(values)
    values = values[order][1:]            # drop the uniform level, a zero mode
    return 2 * math.pi / np.sqrt(values) / 3600.0, vectors[:, order[1:]]


def compare(sea, forcing, frequency_cpd, **kwargs):
    """Dynamic against equilibrium response for one constituent's forcing pattern."""
    dynamic = forced_response(sea, forcing, frequency_cpd, **kwargs)
    equilibrium = sea.equilibrium(forcing)
    scale = np.max(np.abs(equilibrium))
    weights = sea.area / sea.area.sum()
    peak = int(np.argmax(np.abs(equilibrium)))
    return dict(frequency_cpd=frequency_cpd, equilibrium_max_m=float(scale),
                departure_max_relative=float(np.max(np.abs(dynamic - equilibrium)) / scale),
                rms_amplitude_ratio=float(np.sqrt(np.sum(weights * np.abs(dynamic)**2) / np.sum(weights * np.abs(equilibrium)**2))),
                peak_amplitude_ratio=float(np.abs(dynamic[peak]) / np.abs(equilibrium[peak])),
                peak_phase_lag_deg=float(np.degrees(np.angle(dynamic[peak] / equilibrium[peak]))))
