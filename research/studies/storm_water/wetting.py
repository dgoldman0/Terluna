"""The wetting between storms (storm-water study, phase 2a): how the water a storm leaves in the topsoil spreads
downward and dries before the next storm, at any gravity, and how the next storm infiltrates into what is left.

The wetted layer is Smith, Corradini and Melone's (1993) rectangular profile, as Ogden and Saghafian (1997) carry it
between Green and Ampt storms: I mm of water above the soil's initial content, at a uniform effective saturation se0
down to a sharp front at depth Z = I / (theta_e (se0 - se_i)). Mass is I; the front advances under gravity and under
the capillary pull of the dry soil below it,

    theta_e Z dse0/dt = -E - K(se0) - p Ks G(se_i, se0) / Z,        dI/dt = -E,

with Brooks and Corey's curves, K(se) = Ks se^(3 + 2/lam) and the capillary drive
G(a, b) = integral of K dh / Ks = hb (b^m - a^m) / (1 + 3 lam), m = 3 + 1/lam, between the initial and the layer's
contents, and Smith et al.'s profile factor p. The surface dries at the smaller of the demand and the desorptivity's
rate S_d / (2 sqrt(t)), with S_d^2 = 2 theta_e Ks hb se0^(4 + 1/lam) / (1 + 4 lam), the sharp-front analogue of the
sorptivity for drying to the residual content (an estimate). Gravity enters only through Ks (times g) and hb (over g):
the capillary products Ks hb, the diffusivity K dh/dtheta, the drive and the desorptivity are the same at any gravity,
and only the gravity drainage K(se0) falls to a sixth on the Moon.

The next storm starts a new Green and Ampt front into the layer, with the deficit theta_e (1 - se0) and the front
suction phase 1 uses, scaled by the Brooks and Corey drive from se0 to saturation; when that front reaches the old
front's depth the two merge into one front carrying all the water into the initial soil. After the storm the profile
collapses back to one layer at the deeper of the two fronts, keeping its water. The infiltration itself is phase 1's
(infiltration.py), so with no water left from earlier storms each storm meets the mean's soil exactly as before.
"""
from __future__ import annotations

import numpy as np

from research.studies.storm_water import infiltration as IN

# Brooks and Corey's curves by texture class at Earth's gravity: Rawls, Brakensiek and Saxton (1982) class means as the
# GSSHA and SWMM5 tables reproduce them (residual content, bubbling pressure, lam the arithmetic mean), read at that
# secondary. Ks is twice the Green and Ampt conductivity phase 1 uses (Rawls et al. 1983, eq. 6, K = Ks / 2, read at the
# source), so both models describe one soil; the same table lists 1.32 cm/h for a loam (Rawls et al. 1982), twice this.
# The bubbling pressure's spread is carried as the class suction's (phase2a.SUCTION_SD); lam's spread was not found.
CURVES = {
    'sandy_loam': dict(porosity=0.453, residual=0.041, bubbling_mm=146.6, lam=0.378, ks_mm_h=21.8),
    'loam': dict(porosity=0.463, residual=0.027, bubbling_mm=111.5, lam=0.252, ks_mm_h=6.8),
    'sand': dict(porosity=0.437, residual=0.020, bubbling_mm=72.6, lam=0.694, ks_mm_h=235.6),
}
P_SHAPE = 1.7                    # Smith et al.'s profile factor (recalled); Richards' equation agrees (test_phase2a.py)
SUBSTEPS = 48                    # log-spaced substeps per gap (convergence tested against four times as many)
FIRST_SUBSTEP_H = 0.01


def curve(soil, gravity=IN.MOON):
    """A soil's Brooks and Corey curve at a gravity (in Earth's): Ks (mm/h) times g, hb (mm) over g."""
    c = CURVES[soil]
    lam, hb = c['lam'], c['bubbling_mm']
    # phase 1's effective porosity where the soil is one of its surfaces', so both models share one deficit
    theta_e = IN.SOILS[soil]['effective_porosity'] if soil in IN.SOILS else c['porosity'] - c['residual']
    return dict(soil=soil, lam=lam, n=3.0 + 2.0 / lam, m=3.0 + 1.0 / lam, hb=hb / gravity, ks=c['ks_mm_h'] * gravity,
                theta_e=theta_e, residual=c['residual'], porosity=c['porosity'], gravity=gravity)


def conductivity(c, se):
    return c['ks'] * np.power(np.clip(se, 0.0, 1.0), c['n'])


def drive(c, se_a, se_b):
    """Ks G(se_a, se_b) in mm2/h: the unsaturated capillary drive between two contents (the same at any gravity)."""
    return c['ks'] * c['hb'] / (1.0 + 3.0 * c['lam']) * (np.power(se_b, c['m']) - np.power(se_a, c['m']))


def desorptivity_sq(c, se0):
    """S_d^2 (mm2/h) for drying a layer at se0 to the residual content at the surface (sharp-front estimate)."""
    return 2.0 * c['theta_e'] * c['ks'] * c['hb'] * np.power(se0, 4.0 + 1.0 / c['lam']) / (1.0 + 4.0 * c['lam'])


def front_suction_factor(c, se0):
    """The share of the front suction (phase 1's, at a dry soil) left for a saturated front entering soil at se0:
    the Brooks and Corey drive from se0 to saturation, the saturated segment hb included, over that from dry soil."""
    return (2.0 + 3.0 * c['lam'] - np.power(se0, c['m'])) / (2.0 + 3.0 * c['lam'])


class Profile:
    """The soil water of a set of columns: background saturation, the layer's water and saturation, the current storm's
    new front, the surface store, and running totals (mm)."""

    def __init__(self, se_i):
        self.se_i = np.array(se_i, float)
        z = np.zeros_like(self.se_i)
        self.layer = z.copy()            # I, mm above the background
        self.se0 = self.se_i.copy()
        self.f2 = z.copy()               # the current storm's front, mm
        self.store = z.copy()
        self.runoff = z.copy()
        self.energy = z.copy()
        self.evaporated = z.copy()
        self.dry_h = z.copy()            # hours since the surface began to dry

    def depth(self, c):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(self.layer > 0.0, self.layer / (c['theta_e'] * np.maximum(self.se0 - self.se_i, 1e-12)), 0.0)


def _capacity(p, c, pr):
    """Infiltration capacity (mm/h) of the current storm's front on the profile, before this step's intake."""
    th = c['theta_e']
    d_i = IN.deficit(p, pr.se_i)
    z1 = pr.depth(c)
    d2 = np.maximum(th * (1.0 - pr.se0), 1e-9)
    old = pr.layer > 0.0
    merged = ~old | (pr.f2 / d2 >= z1)
    q = dict(p)
    q['psi'] = p['psi'] * np.where(old, front_suction_factor(c, pr.se0), 1.0)
    cap_new = IN.capacity(q, pr.f2, d2, pr.energy) if p['forming'] else IN.capacity(q, pr.f2, d2)
    cap_all = IN.capacity(p, pr.layer + pr.f2, d_i, pr.energy) if p['forming'] else IN.capacity(p, pr.layer + pr.f2, d_i)
    return np.where(merged, cap_all, cap_new), merged


def storm_step(p, c, pr, rate, dt_h):
    """One step of rain (mm/h): phase 1's midpoint rule on the profile's capacity, the store and the runoff."""
    supply = rate * dt_h + pr.store
    de = p['energy_factor'] * IN.earth_rain_energy(rate) * rate * dt_h if p['forming'] else 0.0
    cap, _ = _capacity(p, c, pr)
    first = np.minimum(supply, cap * dt_h)
    f2 = pr.f2
    pr.f2 = f2 + 0.5 * first
    e0 = pr.energy
    pr.energy = e0 + 0.5 * de
    cap, _ = _capacity(p, c, pr)
    pr.f2, pr.energy = f2, e0
    take = np.minimum(supply, cap * dt_h)
    pr.f2 = pr.f2 + take
    held = pr.store + rate * dt_h - take
    off = np.maximum(held - p['retention'], 0.0)
    pr.runoff = pr.runoff + off
    pr.store = held - off
    pr.energy = pr.energy + de


def end_storm(p, c, pr):
    """The store soaks in and the two fronts collapse to one layer at the deeper front, keeping their water."""
    pr.f2 = pr.f2 + pr.store
    pr.store = np.zeros_like(pr.store)
    th = c['theta_e']
    _, merged = _capacity(p, c, pr)
    total = pr.layer + pr.f2
    z_old = pr.depth(c)
    d2 = np.maximum(th * (1.0 - pr.se0), 1e-9)
    z_new = np.where(merged, total / np.maximum(th * (1.0 - pr.se_i), 1e-9), np.maximum(z_old, pr.f2 / d2))
    with np.errstate(divide='ignore', invalid='ignore'):
        se0 = np.where(total > 0.0, pr.se_i + total / (th * np.maximum(z_new, 1e-12)), pr.se_i)
    pr.se0 = np.minimum(se0, 1.0)
    pr.layer = total
    pr.f2 = np.zeros_like(pr.f2)
    pr.dry_h = np.zeros_like(pr.dry_h)
    if p['forming']:
        pr.energy = np.zeros_like(pr.energy)


def redistribute(c, pr, hours, demand_mm_h, p_shape=P_SHAPE, substeps=None):
    """Carry the layer through a gap of `hours` (per column) with drying at the demand (mm/h) where the soil allows."""
    substeps = SUBSTEPS if substeps is None else substeps
    hours = np.broadcast_to(np.asarray(hours, float), pr.se_i.shape)
    edges = np.expm1(np.linspace(0.0, 1.0, substeps + 1)[:, None] * np.log1p(hours / FIRST_SUBSTEP_H)[None, :])
    edges = edges * FIRST_SUBSTEP_H
    th = c['theta_e']
    for k in range(substeps):
        dt = edges[k + 1] - edges[k]
        live = (pr.layer > 1e-9) & (dt > 0.0)
        if not live.any():
            break
        excess = np.maximum(pr.se0 - pr.se_i, 1e-12)
        z = pr.layer / (th * excess)
        tmid = pr.dry_h + 0.5 * dt
        z = np.where(live, z, 1.0)
        e_soil = np.sqrt(desorptivity_sq(c, pr.se0)) / (2.0 * np.sqrt(np.maximum(tmid, 1e-6)))
        e = np.minimum(np.minimum(demand_mm_h, e_soil), pr.layer / np.maximum(dt, 1e-12))
        loss = e + conductivity(c, pr.se0) - conductivity(c, pr.se_i) + p_shape * drive(c, pr.se_i, pr.se0) / z
        dse = -loss * dt / (th * z)
        se_new = np.maximum(pr.se0 + dse, pr.se_i + 0.02 * excess)          # an Euler step never crosses the background
        layer_new = pr.layer - e * dt
        pr.evaporated = np.where(live, pr.evaporated + e * dt, pr.evaporated)
        pr.se0 = np.where(live, se_new, pr.se0)
        pr.layer = np.where(live, layer_new, pr.layer)
        pr.dry_h = np.where(live, pr.dry_h + dt, pr.dry_h)
    gone = pr.layer <= 1e-9
    pr.layer = np.where(gone, 0.0, pr.layer)
    pr.se0 = np.where(gone, pr.se_i, pr.se0)


def mean_storms(p, c, lib, saturations, n_max, mean_hours, demand_mm_h, copies=2, seed=20261101, p_shape=P_SHAPE,
                substeps=None):
    """Expected runoff, water left in the layer and evaporation (mm) of a mean holding exactly N = 0..n_max storms at
    uniformly random times within mean_hours, the layer redistributing and drying between them: (n_max + 1, sat)."""
    rng = np.random.default_rng(seed)
    k = len(lib['storms'])
    w = np.asarray(lib['weights'], float)
    if not np.allclose(w, w[0]):
        raise ValueError('mean_storms expects a library of equally weighted storms')
    dt = lib['step_min'] / 60.0
    reals = copies * k
    sat = np.asarray(saturations, float)
    ns = np.arange(1, n_max + 1)
    # Columns: every N, every realization, every saturation.
    n_col = np.repeat(ns, reals * sat.size)
    s_col = np.tile(sat, ns.size * reals)
    pr = Profile(s_col)
    picks = np.stack([np.concatenate([rng.permutation(k) for _ in range(copies)]) for _ in range(n_max)])  # (pos, real)
    pick_col = np.tile(np.repeat(np.arange(reals), sat.size), ns.size)
    lengths = np.array([len(s) * dt for s in lib['storms']])
    starts = rng.uniform(0.0, mean_hours, size=(ns.size, reals, n_max))
    starts = np.sort(np.where(np.arange(n_max)[None, None, :] < ns[:, None, None], starts, np.inf), axis=2)
    start_col = np.repeat(starts.reshape(ns.size * reals, n_max), sat.size, axis=0)        # (col, pos)
    for pos in range(n_max):
        active = n_col > pos
        chosen = picks[pos][pick_col]
        if pos > 0:
            prev_end = start_col[:, pos - 1] + lengths[picks[pos - 1][pick_col]]
            gap = np.where(active, np.maximum(start_col[:, pos] - prev_end, 0.0), 0.0)
            redistribute(c, pr, gap, demand_mm_h, p_shape, substeps)
        n_steps = int(max(len(lib['storms'][i]) for i in np.unique(chosen[active])))
        rates = np.zeros((n_steps, pr.se_i.size))
        for i in np.unique(chosen[active]):
            s = np.asarray(lib['storms'][i], float)
            m = active & (chosen == i)
            rates[: len(s), m] = s[:, None]
        for t in range(n_steps):
            storm_step(p, c, pr, rates[t], dt)
        end_storm(p, c, pr)
    out = []
    for arr in (pr.runoff, pr.layer, pr.evaporated):
        a = arr.reshape(ns.size, reals, sat.size).mean(axis=1)
        out.append(np.vstack([np.zeros((1, sat.size)), a]))
    return tuple(out)


def below_root_zone(c, water_mm, se_i, hours, demand_mm_h, root_zone_mm, p_shape=P_SHAPE):
    """Water (mm) of an infiltrated depth that lies below the root zone after a period of redistribution and drying,
    the profile starting saturated at the surface: the rectangular layer's share beyond root_zone_mm."""
    water = np.asarray(water_mm, float)
    pr = Profile(np.broadcast_to(se_i, water.shape))
    pr.layer = water.copy()
    pr.se0 = np.ones_like(water)
    redistribute(c, pr, hours, demand_mm_h, p_shape)
    z = pr.depth(c)
    with np.errstate(divide='ignore', invalid='ignore'):
        share = np.where(z > root_zone_mm, 1.0 - root_zone_mm / np.maximum(z, 1e-12), 0.0)
    return pr.layer * share, pr.layer, z


def richards_redistribution(c, water_mm, se_i, hours, dz_mm=2.0, depth_mm=1200.0):
    """A numerical control: Richards' equation (Celia et al.'s mixed form, implicit, Picard iteration) for a column of
    the same Brooks and Corey soil at the same gravity, a saturated layer holding water_mm above se_i at the start,
    no flux through the surface and free drainage at the base. Returns the surface layer's effective saturation and the
    depth (mm) above which nine tenths of the added water lies, at each of `hours`."""
    from scipy.linalg import solve_banded
    n = int(depth_mm / dz_mm)
    th_r, th_e, hb, lam, ks, nn = c['residual'], c['theta_e'], c['hb'], c['lam'], c['ks'], c['n']

    def theta(psi):
        se = np.where(psi < -hb, (hb / np.maximum(-psi, hb)) ** lam, 1.0)
        return th_r + th_e * se, se

    def capacity_(psi):
        return np.where(psi < -hb, th_e * lam / np.maximum(-psi, hb) * (hb / np.maximum(-psi, hb)) ** lam, 1e-7)

    psi = np.full(n, -hb * se_i ** (-1.0 / lam))
    z0 = water_mm / (th_e * (1.0 - se_i))
    wet = (np.arange(n) + 0.5) * dz_mm < z0
    psi[wet] = -hb
    th_bg = theta(np.array([psi[-1]]))[0][0]
    # Scale the wet cells so the column holds exactly the added water (the last wet cell partly filled).
    th0, _ = theta(psi)
    excess = (th0 - th_bg).sum() * dz_mm
    if excess > water_mm and wet.any():
        k = np.flatnonzero(wet)[-1]
        frac = 1.0 - (excess - water_mm) / ((th0[k] - th_bg) * dz_mm)
        se_k = se_i + frac * (1.0 - se_i)
        psi[k] = -hb * max(se_k, 1e-9) ** (-1.0 / lam) if se_k < 1.0 else -hb
    t, dt, out = 0.0, 1e-5, []
    for t_out in np.atleast_1d(hours):
        while t < t_out - 1e-12:
            step = min(dt, t_out - t)
            th_old, _ = theta(psi)
            new = psi.copy()
            ok = False
            for it in range(40):
                th_m, se_m = theta(new)
                k_m = ks * se_m ** nn
                kf = 0.5 * (k_m[1:] + k_m[:-1])                          # conductivity between cells
                cm = capacity_(new)
                # q_down(i+1/2) = -kf ((psi_{i+1} - psi_i)/dz - 1); bottom free drainage q = K_n.
                a = np.zeros(n); b = cm / step; cc = np.zeros(n)
                b[:-1] += kf / dz_mm ** 2; cc[:-1] = -kf / dz_mm ** 2
                b[1:] += kf / dz_mm ** 2; a[1:] = -kf / dz_mm ** 2
                rhs = cm * new / step - (th_m - th_old) / step
                rhs[:-1] -= kf / dz_mm
                rhs[1:] += kf / dz_mm
                rhs[-1] -= k_m[-1] / dz_mm
                band = np.vstack([np.r_[0.0, cc[:-1]], b, np.r_[a[1:], 0.0]])
                sol = solve_banded((1, 1), band, rhs)
                change = np.abs(sol - new).max()
                new = sol
                if change < 1e-3 * hb:
                    ok = True
                    break
            if not ok:
                dt *= 0.3
                continue
            psi, t = new, t + step
            dt = min(dt * (1.3 if it < 6 else 0.9), 0.5)
        th, se = theta(psi)
        added = np.cumsum((th - th_bg) * dz_mm)
        z90 = float((np.searchsorted(added, 0.9 * added[-1]) + 1) * dz_mm)
        out.append((float(se[0]), z90))
    return out
