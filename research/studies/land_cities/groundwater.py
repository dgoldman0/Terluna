"""The water table in two dimensions around a lava-tube candidate: a steady Dupuit mound on the atlas grid, capped by
the ground, and how it rises once the rain begins.

Unconfined flow over a flat base Zb with uniform conductivity K obeys div(K H grad H) + R = 0, H = Z - Zb; with
Phi = H^2 / 2 that is K lap(Phi) + R = 0, a linear Poisson equation (the Dupuit mound of Woessner & Poeter 2020, Box 6,
in two dimensions). The atlas's seas and lakes hold the water table at their own levels; where the mound would rise
above the ground, the ground holds it (a seepage face: springs and wet ground), found by pinning those cells and
releasing any that would draw water, until neither changes. The rise in time is the Boussinesq equation
S dZ/dt = div(K H grad Z) + R, stepped implicitly from the water table at the sea's level under the land (the basalt
below it saturated), with the ground holding the water table the same way within every step.

Finite volumes on the atlas's 4 px/deg cells over a window about the candidate, with no flow across the window's
edge: north and south faces carry cos(latitude) of their length, east and west faces 1/cos(latitude).
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS

from research.studies.land_cities import grid

SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY
WINDOW_KM = 400.0               # half-width of the window about a candidate
MAX_PINNING_ROUNDS = 200
# The rise's implicit time steps: the first step, each step's growth over the last, and the run's length (years).
FIRST_STEP_YEARS = 0.5
STEP_GROWTH = 1.1
HORIZON_YEARS = 20000.0
MAX_SET_ROUNDS = 50             # passes per time step to settle which cells stand at the ground


class Window:
    """Cells about a candidate, with the finite-volume geometry of the Laplacian on them."""

    def __init__(self, wet, level, ground, area_m2, lat_deg, centre, ppd):
        self.wet, self.level, self.ground, self.area = wet, level, ground, area_m2
        self.centre = centre
        self.shape = wet.shape
        lat_c = np.radians(np.asarray(lat_deg, float))
        half = math.radians(0.5 / ppd)
        ones = np.ones((1, self.shape[1]))
        self.g_north = np.cos(lat_c + half)[:, None] * ones
        self.g_south = np.cos(lat_c - half)[:, None] * ones
        self.g_east = (1.0 / np.cos(lat_c))[:, None] * ones

    @classmethod
    def about(cls, at: grid.Atlas, lat, lon, half_km=WINDOW_KM):
        """The atlas cells within about half_km of a point each way (longitude periodic)."""
        half_lat = half_km / (math.radians(1.0) * grid.RADIUS_KM)
        half_lon = min(half_lat / max(math.cos(math.radians(lat)), 0.2), 90.0)
        i0, j0 = at.cell(lat, lon)
        n_lat, n_lon = int(round(half_lat * at.ppd)), int(round(half_lon * at.ppd))
        rows = np.arange(max(i0 - n_lat, 0), min(i0 + n_lat, at.lat.size - 1) + 1)
        cols = np.arange(j0 - n_lon, j0 + n_lon + 1) % at.lon.size
        ix = np.ix_(rows, cols)
        win = cls(at.wet()[ix], at.water_level_m()[ix], at.height_m[ix], at.area[ix] * 1e6, at.lat[rows],
                  (i0 - rows[0], n_lon), at.ppd)
        win.ix = ix                 # to take other atlas fields over the same cells
        return win

    def neighbours(self, cells):
        """For each face direction: the cells' rows and columns, their neighbours' and the face geometry."""
        ny, nx = self.shape
        for di, dj, g in ((-1, 0, self.g_north), (1, 0, self.g_south), (0, -1, self.g_east), (0, 1, self.g_east)):
            ys, xs = np.nonzero(cells)
            yy, xx = ys + di, xs + dj
            ok = (yy >= 0) & (yy < ny) & (xx >= 0) & (xx < nx)
            yield ys[ok], xs[ok], yy[ok], xx[ok], g[ys[ok], xs[ok]]

    def operator(self, fixed, weight=None):
        """-div(weight grad .) on the free cells as a sparse matrix, its index map, and each free cell's coupling to
        its fixed neighbours (free index, fixed flat index, coefficient). Faces take the mean weight of their cells."""
        free = ~fixed
        n = int(free.sum())
        idx = -np.ones(self.shape, int)
        idx[free] = np.arange(n)
        w = np.ones(self.shape) if weight is None else weight
        diag = np.zeros(self.shape)
        rows, cols, vals, couple = [], [], [], []
        for ys, xs, yy, xx, g in self.neighbours(free):
            c = g * 0.5 * (w[ys, xs] + w[yy, xx])
            np.add.at(diag, (ys, xs), c)
            inner = free[yy, xx]
            rows.append(idx[ys[inner], xs[inner]])
            cols.append(idx[yy[inner], xx[inner]])
            vals.append(-c[inner])
            couple.append((idx[ys[~inner], xs[~inner]], yy[~inner] * self.shape[1] + xx[~inner], c[~inner]))
        rows.append(np.arange(n))
        cols.append(np.arange(n))
        vals.append(diag[free])
        A = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
        return A, idx, couple


def coupled(rhs, couple, values):
    out = rhs.copy()
    for fi, fj, c in couple:
        np.add.at(out, fi, c * values[fj])
    return out


def net_inflow(win: Window, field, cells, weight=None):
    """Water each of the given cells takes from its neighbours (per unit conductivity): the sum over its faces of the
    face geometry times the face's mean weight times the difference of the field across it."""
    w = np.ones(win.shape) if weight is None else weight
    net = np.zeros(win.shape)
    for ys, xs, yy, xx, g in win.neighbours(cells):
        np.add.at(net, (ys, xs), g * 0.5 * (w[ys, xs] + w[yy, xx]) * (field[yy, xx] - field[ys, xs]))
    return net


def seepage(win: Window, phi, source, cells):
    """Net water leaving each of the given cells by seepage at steady state: its recharge and the inflow from its
    neighbours."""
    return np.where(cells, source, 0.0) + net_inflow(win, phi, cells)


def steady(win: Window, recharge, k, base, cap=True):
    """The steady water table (m above the geoid) on the window, with the seas and lakes at their levels and a flat
    base `base` below the window's lowest water; the ground caps the mound where cap is set. Returns the field, the
    pinned cells, the base and the pinning rounds."""
    lowest = float(win.level[win.wet].min()) if win.wet.any() else float(win.ground.min())
    zb = lowest - base
    phi_fixed = np.where(win.wet, 0.5 * np.maximum(np.nan_to_num(win.level, nan=lowest) - zb, 0.0) ** 2, 0.0)
    phi_ground = 0.5 * np.maximum(win.ground - zb, 0.0) ** 2
    source = recharge * win.area / k
    scale = max(float(np.abs(source).max()), 1e-30)
    pinned = np.zeros(win.shape, bool)
    rounds = 0
    while True:
        fixed = win.wet | pinned
        A, _, couple = win.operator(fixed)
        values = np.where(pinned, phi_ground, phi_fixed)
        phi = values.copy()
        phi[~fixed] = spl.spsolve(A, coupled(source[~fixed].astype(float), couple, values.ravel()))
        z = zb + np.sqrt(2.0 * np.maximum(phi, 0.0))
        rounds += 1
        if not cap:
            break
        new = (~fixed) & (z > win.ground + 1e-6)
        release = pinned & (seepage(win, phi, source, pinned) < -1e-9 * scale)
        if not new.any() and not release.any():
            break
        if rounds >= MAX_PINNING_ROUNDS:
            raise RuntimeError('the wet ground did not settle in the steady solve')
        pinned = (pinned | new) & ~release
    return dict(z=z, pinned=pinned, base=zb, rounds=rounds)


def rise(win: Window, recharge, k, specific_yield, start_level, end, first_step=None, growth=None, horizon=None):
    """The water table at the window's centre through time, from start_level under the land (the ground where that
    lies lower; the seas and lakes at their own levels) towards the end (steady) field, by the Boussinesq equation
    S dZ/dt = div(K H grad Z) + R stepped implicitly on growing steps, with K H taken from the step's start. Cells
    whose water would stand above the ground hold it at the ground and shed the excess (springs and wet ground),
    and let go once they would draw water from the surface, settled within each step as in the steady solve.
    Returns the times (years) and the centre's level after each step, and the factorisations taken."""
    ci, cj = win.centre
    base = end['base']
    z = np.where(win.wet, win.level, np.minimum(start_level, win.ground))
    pinned = ~win.wet & (z >= win.ground)
    times, levels = [0.0], [float(z[ci, cj])]
    growth = STEP_GROWTH if growth is None else growth
    horizon = HORIZON_YEARS if horizon is None else horizon
    t, dt_years, solves = 0.0, FIRST_STEP_YEARS if first_step is None else first_step, 0
    while t < horizon:
        dt = dt_years * SECONDS_PER_YEAR
        weight = k * np.maximum(z - base, 1.0)
        storage = specific_yield * win.area / dt
        for attempt in range(MAX_SET_ROUNDS + 1):
            if attempt == MAX_SET_ROUNDS:
                raise RuntimeError('the wet ground did not settle within a step of the rise')
            fixed = win.wet | pinned
            free = ~fixed
            L, _, couple = win.operator(fixed, weight=weight)
            values = np.where(pinned, win.ground, z)
            new = values.copy()
            rhs = coupled(storage[free] * z[free] + recharge[free] * win.area[free], couple, values.ravel())
            new[free] = spl.spsolve((sp.diags(storage[free]) + L).tocsc(), rhs)
            solves += 1
            above = free & ~win.wet & (new > win.ground + 1e-6)
            # A held cell sheds what reaches it beyond what refills its own pores; below zero it would draw water.
            shed = (recharge * win.area + net_inflow(win, new, pinned, weight)
                    - storage * (win.ground - z))
            release = pinned & (shed < 0.0)
            if not above.any() and not release.any():
                break
            pinned = (pinned | above) & ~release
        z = new
        t += dt_years
        times.append(t)
        levels.append(float(z[ci, cj]))
        dt_years *= growth
    return dict(times=np.array(times), levels=np.array(levels), solves=solves)


def first_reached(series, level):
    """The first time (years) the centre reaches a level, interpolated between steps; 0 if there from the start and
    None if never within the run."""
    t, z = series['times'], series['levels']
    if z[0] >= level:
        return 0.0
    hit = np.flatnonzero(z >= level)
    if hit.size == 0:
        return None
    i = int(hit[0])
    return float(t[i - 1] + (level - z[i - 1]) / (z[i] - z[i - 1]) * (t[i] - t[i - 1]))
