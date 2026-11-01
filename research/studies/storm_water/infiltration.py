"""Green and Ampt infiltration on the storm-water study's four designed surfaces, at any gravity.

Each surface is a soil with an optional crust on top. The capacity at wetted depth L = F / dtheta (F the water taken in
since the wetting began, dtheta the soil's deficit) is K (psi + L) / L on a bare soil and
(psi + L) / (R + (L - Zc) / K) under a crust of thickness Zc and resistance R = Zc / Kc, with K (psi + L) / L inside the
crust's own depth using its conductivity Kc. This is Green and Ampt's two-layer form (Hillel and Gardner's crust-topped
profile), and with Rawls et al.'s (1990) seal resistance it reproduces their crust factor exactly (test_run.py).

Gravity enters through two terms and nothing else. A conductivity is k rho g / mu, so every K and every Kc scales with
g, and a resistance with 1/g; a suction head is a capillary pressure over rho g, so psi scales with 1/g. The products
K psi, and with them the sorptivity and the capillary part of the capacity, are the same at any gravity; the gravity
part, K, falls to a sixth on the Moon. A seal forming under a storm's drops follows WEPP's exponential approach to its
final conductance in the cumulative kinetic energy of the rain, and a raindrop's kinetic energy per millimetre of
rain scales with v^2, so with g / rho_air (the plan's scaling, fall speed ~ sqrt(g / rho_air) at a fixed drop size).

Time steps are one minute; within a step the capacity is taken at the step's midpoint in F. Rain fills a surface
store of retention_mm first; what the store cannot hold runs off, and what it holds infiltrates after the rain.
"""
from __future__ import annotations

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_ATMOSPHERE, STANDARD_GRAVITY

MOON = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY            # the Moon's gravity in Earth's
WSMAX_M = 0.5            # PlaSim's bucket capacity, WSMAX_EARTH (landmod.f90 line 14), unchanged by the design run
DRHSFULL = 0.4           # PlaSim's wetness factor reaches 1 at this share of the bucket (landmod.f90, drhsfull)

# Rawls, Brakensiek and Miller (1983) texture-class means (Chow, Maidment and Mays 1988, Table 4.3.1; the conductivities
# and suctions as HEC-HMS reproduces them); field capacity (33 kPa) and wilting point (1500 kPa) from Rawls,
# Brakensiek and Saxton (1982), recalled, not re-read.
SOILS = {
    'sandy_loam': dict(ks_mm_h=10.9, suction_mm=110.1, porosity=0.453, effective_porosity=0.412,
                       field_capacity=0.207, wilting_point=0.095),
    'loam': dict(ks_mm_h=3.4, suction_mm=88.9, porosity=0.463, effective_porosity=0.434,
                 field_capacity=0.270, wilting_point=0.117),
}
# The surfaces are designed soils (design choices); the crusts' Earth parameters carry their sources.
# Rawls et al. (1990) Table 1 for a loam: the matric drop across the seal and the subcrust factor; WEPP's 5 mm crust;
# WEPP's stability regression (eq. 7.9.8) for a loam of 40% sand, 20% clay and 15 meq/100 g gives 0.0034 m2/J.
SEAL = dict(kind='seal', thickness_mm=5.0, matric_drop_mm=70.0, subcrust_factor=0.82)
SURFACES = {
    'sandy_loam': dict(soil='sandy_loam', crust=None, retention_mm=1.0),
    'loam': dict(soil='loam', crust=None, retention_mm=1.0),
    'sealed_loam': dict(soil='loam', crust=dict(SEAL, forming=False), retention_mm=1.0),
    'sealed_loam_forming': dict(soil='loam', crust=dict(SEAL, forming=True, stability_m2_j=0.0034), retention_mm=1.0),
    'biocrust': dict(soil='sandy_loam', crust=dict(kind='biocrust', thickness_mm=5.0, conductivity_ratio=0.5),
                     retention_mm=1.0),
}
MAIN = ('sandy_loam', 'loam', 'sealed_loam', 'biocrust')


def drop_energy_factor(pressure_pa, gravity=MOON):
    """Kinetic energy per mm of rain against Earth's at the same drop size and air temperature: v^2 ~ g / rho_air."""
    return gravity * STANDARD_ATMOSPHERE / np.asarray(pressure_pa, float)


def earth_rain_energy(rate_mm_h):
    """Kinetic energy of Earth's natural rain (J/m2 per mm) at an intensity: Brown and Foster (1987), recalled."""
    return 29.0 * (1.0 - 0.72 * np.exp(-0.05 * np.asarray(rate_mm_h, float)))


def parameters(name, gravity=MOON, ks_factor=1.0, suction_factor=1.0, retention_mm=None, ratio=None,
               energy_factor=None, pressure_pa=None):
    """A surface's Green-Ampt parameters at a gravity (in Earth's): mm, mm/h, h."""
    s = SURFACES[name]
    soil = SOILS[s['soil']]
    ks = soil['ks_mm_h'] * ks_factor * gravity
    p = dict(name=name, gravity=gravity, psi=soil['suction_mm'] * suction_factor / gravity, ks=ks,
             k_sub=ks, zc=0.0, r=0.0, kc=ks, forming=False,
             retention=s['retention_mm'] if retention_mm is None else retention_mm, soil=soil)
    c = s['crust']
    if c is not None:
        zc = c['thickness_mm']
        if c['kind'] == 'seal':
            sc = c['subcrust_factor']
            # Rawls et al.'s Earth resistance (their eq. 4, from the matric drop measured on Earth), then 1/g
            r_final = (c['matric_drop_mm'] + zc) / (sc * soil['ks_mm_h'] * ks_factor) / gravity
            p.update(zc=zc, r=r_final, kc=zc / r_final, k_sub=sc * ks, subcrust_factor=sc)
            if c['forming']:
                if energy_factor is None:
                    energy_factor = 1.0 if gravity == 1.0 else float(drop_energy_factor(pressure_pa, gravity))
                p.update(forming=True, stability=c['stability_m2_j'], energy_factor=energy_factor, r_final=r_final)
        else:
            kc = (c['conductivity_ratio'] if ratio is None else ratio) * ks
            p.update(zc=zc, kc=kc, r=zc / kc)
    return p


def deficit(p, relative_saturation):
    """The soil's moisture deficit dtheta for an initial effective saturation (0 at residual, 1 at saturation)."""
    return (1.0 - np.asarray(relative_saturation, float)) * p['soil']['effective_porosity']


def initial_saturation(soil, mrso_m, mapping='plant_available'):
    """The topsoil's effective saturation at storm time from PlaSim's bucket water (m). plant_available reads the
    bucket as the water between the wilting point (empty) and field capacity (full); wetness_factor fills that range
    at PlaSim's wetness factor min(1, W / (0.4 wsmax)); air_dry leaves the topsoil at its residual water."""
    sb = np.clip(np.asarray(mrso_m, float) / WSMAX_M, 0.0, 1.0)
    if mapping == 'wetness_factor':
        sb = np.clip(np.asarray(mrso_m, float) / (DRHSFULL * WSMAX_M), 0.0, 1.0)
    elif mapping == 'air_dry':
        return np.zeros_like(sb)
    theta = soil['wilting_point'] + sb * (soil['field_capacity'] - soil['wilting_point'])
    residual = soil['porosity'] - soil['effective_porosity']
    return (theta - residual) / soil['effective_porosity']


def capacity(p, f_mm, dtheta, energy=None):
    """Infiltration capacity (mm/h) after f_mm of infiltration at deficit dtheta."""
    depth = np.asarray(f_mm, float) / dtheta
    psi, zc = p['psi'], p['zc']
    kc, r, k_sub = p['kc'], p['r'], p['k_sub']
    if p['forming']:
        w = np.exp(-p['stability'] * energy)
        g = 1.0 / p['r_final'] + (p['ks'] / zc - 1.0 / p['r_final']) * w      # the seal's conductance, 1/h
        r, kc = 1.0 / g, zc * g
        k_sub = p['ks'] * (p['subcrust_factor'] + (1.0 - p['subcrust_factor']) * w)
    with np.errstate(divide='ignore', invalid='ignore'):
        inside = np.where(depth > 0.0, kc * (psi + depth) / depth, np.inf)
        if zc == 0.0:
            return inside
        below = (psi + depth) / (r + (depth - zc) / k_sub)
    return np.where(depth <= zc, inside, below)


class State:
    """Water taken in (mm), water held at the surface (mm), runoff (mm) and the storm's rain energy (J/m2)."""

    def __init__(self, shape, f0=0.0):
        self.f = np.zeros(shape) + f0
        self.store = np.zeros(shape)
        self.runoff = np.zeros(shape)
        self.energy = np.zeros(shape)


def step(p, st, dtheta, rate, dt_h):
    """One step of rain (mm/h) over dt_h hours: the midpoint capacity takes what it can of the rain and the store."""
    supply = rate * dt_h + st.store
    de = 0.0
    if p['forming']:
        de = p['energy_factor'] * earth_rain_energy(rate) * rate * dt_h
    first = np.minimum(supply, capacity(p, st.f, dtheta, st.energy) * dt_h)
    take = np.minimum(supply, capacity(p, st.f + 0.5 * first, dtheta, st.energy + 0.5 * de) * dt_h)
    st.f = st.f + take
    held = st.store + rate * dt_h - take
    off = np.maximum(held - p['retention'], 0.0)
    st.runoff = st.runoff + off
    st.store = held - off
    st.energy = st.energy + de


def end_storm(p, st):
    """Between storms the held water soaks in, and a forming seal starts the next storm fresh."""
    st.f = st.f + st.store
    st.store = np.zeros_like(st.store)
    if p['forming']:
        st.energy = np.zeros_like(st.energy)


def single_storms(p, lib, saturations):
    """Runoff (mm) of each library storm on fresh soil at each initial saturation: (storm, saturation)."""
    k, dt = len(lib['storms']), lib['step_min'] / 60.0
    n = max(len(s) for s in lib['storms'])
    rates = np.zeros((n, 1, k))
    for j, s in enumerate(lib['storms']):
        rates[: len(s), 0, j] = s
    dth = deficit(p, saturations)[:, None] * np.ones((1, k))
    st = State(dth.shape)
    for t in range(n):
        step(p, st, dth, rates[t], dt)
    end_storm(p, st)
    return st.runoff.T, st.f.T


def sequences(p, lib, saturations, n_max, copies=4, seed=20261101):
    """Expected runoff and water taken in (mm) after n = 0..n_max storms drawn from the library, the wetting front
    carried from storm to storm: arrays (n_max + 1, saturation). The libraries weight their storms equally, so each
    position draws every storm the same number of times (copies shuffled permutations), which makes the first storm's
    expectation exact and keeps every position's draw balanced."""
    rng = np.random.default_rng(seed)
    w = np.asarray(lib['weights'], float)
    k = len(lib['storms'])
    if not np.allclose(w, w[0]):
        raise ValueError('sequences expects a library of equally weighted storms')
    realizations = copies * k
    picks = np.stack([np.concatenate([rng.permutation(k) for _ in range(copies)]) for _ in range(n_max)])
    dt = lib['step_min'] / 60.0
    dth = deficit(p, saturations)[:, None] * np.ones((1, realizations))
    st = State(dth.shape)
    runoff, taken = [np.zeros(dth.shape[0])], [np.zeros(dth.shape[0])]
    for pos in range(n_max):
        chosen = [lib['storms'][i] for i in picks[pos]]
        n = max(len(s) for s in chosen)
        rates = np.zeros((n, realizations))
        for m, s in enumerate(chosen):
            rates[: len(s), m] = s
        for t in range(n):
            step(p, st, dth, rates[t][None, :], dt)
        end_storm(p, st)
        runoff.append(st.runoff.mean(axis=1))
        taken.append(st.f.mean(axis=1))
    return np.array(runoff), np.array(taken)


def steady_rain(p, depth_mm, rate_mm_h, dtheta, f0, dt_h=0.25):
    """Runoff (mm) of rain at a constant rate (mm/h) to a given depth on a front holding f0 mm: arrays alike."""
    depth = np.asarray(depth_mm, float)
    st = State(depth.shape, f0)
    left = depth.copy()
    while np.any(left > 1e-9):
        r = np.minimum(left, rate_mm_h * dt_h) / dt_h
        step(p, st, dtheta, r, dt_h)
        left = left - r * dt_h
    end_storm(p, st)
    return st.runoff
