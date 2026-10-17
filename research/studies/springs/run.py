"""The springs' load: how long groundwater stays underground, what it dissolves on the way, how hard the springs are,
where they build travertine, and what the groundwater carries of the weathering and carbon the carbon cycle counts.

    python -m research.studies.springs.run          # -> results/springs.json

Groundwater's third stage, on the author's word of 14 October 2026. Stage 1 (geography/groundwater_uptake.py) follows
how fast the crust takes the water down; stage 2 (geography/groundwater_table.py) solves the water table on land with
its springs, seeps and the water it trades with the seas and lakes. This study re-solves stage 2's states and follows
the water through them:

- transit times: each cell's spread of times from rain to outflow, through its soil, the unsaturated rock above the
  water table and the saturated crust below it, where the permeability falls with depth so that a little of the water
  goes deep and stays long; the spread at each spring, lake and sea is the mixture of the recharge that reaches it along
  stage 2's flow;
- what the water dissolves: the soil water carries the carbon cycle's catchment load at the soils' age; above a deep
  water table the water keeps reacting while open to the CO2 in the rock's pores; below it the water consumes the CO2
  it carries, at a rate that slows as the wetted rock ages, and lays calcite in the rock where it passes saturation;
- the springs' hardness and travertine: the water at each spring, after it meets the air, against the saturations at
  which Earth's streams and lake barriers build tufa;
- the carbon cycle: the weathering the groundwater carries against the catchment law's total, and the carbon the filling
  crust takes down with the water it draws from the land, the lakes and the seas;
- warm springs: how much of the springs' water circulates deep enough to be warmed.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from geography import groundwater_table as gt
from geography import groundwater_uptake as gu
from research.studies.carbon_cycle import run as cc
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'springs.json'
GRIDS = ROOT / 'research' / 'runs' / 'resources' / 'groundwater' / 'springs_2ppd.npz'
SCHEMA = 'terluna.research.springs/1'
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
ROCK = 'research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json'
AIR = 'research/studies/atmospheric_co2/results/atmospheric_co2.json'
INPUTS = {'geography/results/atlas.json': 'terluna.geography.atlas/1',
          'geography/results/groundwater_uptake.json': 'terluna.geography.groundwater-uptake/1',
          'geography/results/groundwater_table.json': 'terluna.geography.groundwater-table/1',
          'geography/results/drainage.json': 'terluna.geography.drainage/1',
          ROCK: 'terluna.research.lunar-cycle-ecology/1', AIR: 'terluna.research.atmospheric-co2/1'}
IGNORED = gt.IGNORED
# The carbon cycle's weathering law, soils' factor and constants are imported from its code and recorded in the
# setting, where the tests hold them to the carbon cycle's own; its code is left out of the binding so that the carbon
# cycle, which reads this product, can change without making it stale.
CODE = ('research/studies/springs/run.py', 'geography/groundwater_table.py', 'geography/groundwater_uptake.py')
ATM_PA = 101325.0
CACO3_G_MOL = 100.087
# Travertine thresholds, as the terrain study uses them: Earth's streams build tufa conspicuously above a mean saturation
# index of 0.75 (Merz-Preiss & Riding 1999); Plitvice's barriers grow above three times saturation with pH above 8;
# engineered builders are taken to build from twice saturation (design guess).
THRESHOLDS = dict(tufa=0.75, plitvice=math.log10(3.0), designed=math.log10(2.0))

# --- Transit -------------------------------------------------------------------------------------------------------
# The unsaturated rock drains at unit gradient, so its water content is the one at which its conductivity equals the
# recharge: theta = phi (R/K)^(1/(2b + 3)) (Campbell 1974), with b for loamy sand (Clapp & Hornberger 1978) in the
# shattered rock. The developed soil above it: a root zone of 1.5 m (design guess) of loam, porosity 0.45, saturated
# conductivity 7e-6 m/s and b 5.4 (Clapp & Hornberger 1978).
B_ROCK = 4.38
SOIL_DEPTH_M = 1.5
SOIL_POROSITY = 0.45
SOIL_K_M_S = 7.0e-6
B_SOIL = 5.39
# Exceedance levels for the transit-time spread: P(T > t) runs from 1 to 0; the nodes sit at the middles of these
# bins, finer in the long tail.
P_EDGES = np.concatenate([np.linspace(1.0, 0.1, 19), np.geomspace(0.1, 1e-4, 7)[1:], [0.0]])
AGE_THRESHOLDS_YR = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6)
CDF_YEARS = 10.0 ** np.arange(-1.0, 7.01, 0.25)   # the spread carried to the springs, every quarter decade
STATE_AGE_YR = dict(young_1000_yr=1000.0, young_10000_yr=10000.0, mature=1.0e6)
DEPTH_STEP_M = 2.0


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def log(msg, t0=[time.time()]):
    print(f'{time.time() - t0[0]:7.1f} s  {msg}', file=sys.stderr, flush=True)


def profile(k_1km, phi0, phi_e, bottom_m=gu.PORE_BOTTOM_M):
    """The crust's tables on a depth grid: conductivity K (m/s) at stage 2's water temperature, porosity, KI(z), the
    conductivity integrated from z to the bottom (m2/s), and PHI(z), the porosity over KI integrated from the top."""
    z = np.arange(0.0, bottom_m + DEPTH_STEP_M, DEPTH_STEP_M)
    z[-1] = min(z[-1], bottom_m)
    k = gu.permeability(z, k_1km) * gt.RHO_WATER * MOON_SURFACE_GRAVITY / gt.viscosity(gt.WATER_K)
    phi = phi0 * np.exp(-z / phi_e)
    seg = 0.5 * (k[1:] + k[:-1]) * np.diff(z)
    ki = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    # PHI diverges where KI reaches zero at the bottom; integrate in midpoints and leave the last node at the limit.
    mid = 0.5 * (ki[1:] + ki[:-1])
    f = 0.5 * (phi[1:] + phi[:-1]) / np.maximum(mid, 1e-300)
    big = np.concatenate([[0.0], np.cumsum(f * np.diff(z))])
    return dict(z=z, k=k, phi=phi, ki=ki, big=big, pores=np.concatenate([[0.0], np.cumsum(
        0.5 * (phi[1:] + phi[:-1]) * np.diff(z))]))


def at_depth(tab, key, depth):
    return np.interp(depth, tab['z'], tab[key])


def depth_where_ki(tab, ki):
    """The depth at which KI falls to the given values (KI decreases with depth)."""
    return np.interp(-ki, -tab['ki'], tab['z'])


def unsaturated_years(tab, depth, rech_m_yr):
    """Years for recharge to cross the soil, and then the rock between the soil and the water table, at unit gradient.
    The rock's water content is phi (R/K)^(1/n) where the recharge R is below the conductivity K and phi where it is
    not, and K falls with depth, so the crossing time is R^(1/n - 1) J(z') + (P(d) - P(z'))/R with z' the shallower of
    the water table d and the depth where K falls to R: J and P are the profile's running integrals of phi K^(-1/n) and
    phi below the soil, one-dimensional tables shared by every cell."""
    r = np.maximum(rech_m_yr, 1e-9)
    n = 2.0 * B_ROCK + 3.0
    soil_k = SOIL_K_M_S * YEAR_S
    theta_soil = SOIL_POROSITY * np.minimum(1.0, (r / soil_k) ** (1.0 / (2.0 * B_SOIL + 3.0)))
    soil = theta_soil * SOIL_DEPTH_M / r
    z, k, phi = tab['z'], tab['k'] * YEAR_S, tab['phi']
    run = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(z))])
    j_tab, p_tab = run(phi * k ** (-1.0 / n)), run(phi)
    top = SOIL_DEPTH_M
    j0, p0 = np.interp(top, z, j_tab), np.interp(top, z, p_tab)
    d = np.maximum(depth, top)
    # The depth where K falls to R: K holds its cap near the surface, then falls; search only the falling part.
    falling = k < k[0]
    if falling.any():
        z_star = np.where(r >= k[0], top, np.interp(-r, -k[falling], z[falling], left=z[falling][0], right=z[-1]))
    else:
        z_star = np.where(r >= k[0], top, z[-1])
    z_lim = np.clip(np.minimum(d, z_star), top, None)
    wet = (r ** (1.0 / n - 1.0)) * (np.interp(z_lim, z, j_tab) - j0)
    full = (np.interp(d, z, p_tab) - np.interp(z_lim, z, p_tab)) / r
    return soil, np.where(depth > top, wet + full, 0.0)


def saturated_spread(tab, depth, rech_m_yr):
    """The saturated crust's transit times, a generalized Vogel model: in Dupuit flow with uniform recharge R over a
    strip, the share of the outflow older than t is G(z), the share of the flow passing below the depth z that the
    streamline of age t reaches, with t(z) = (1/R) * integral from the water table to z of phi/G, whatever the strip's
    length (Vogel 1967 for constant properties). Returns the node depths and ages for the exceedance bins, and the mean
    age, the saturated pores over the recharge."""
    r = np.maximum(rech_m_yr, 1e-9)
    ki_d = at_depth(tab, 'ki', depth)
    big_d = at_depth(tab, 'big', depth)
    nodes = 0.5 * (P_EDGES[:-1] + P_EDGES[1:])
    zp = depth_where_ki(tab, nodes[None, :] * ki_d[:, None])
    age = (ki_d[:, None] / (r[:, None] / YEAR_S)) * (at_depth(tab, 'big', zp) - big_d[:, None]) / YEAR_S
    mean = (tab['pores'][-1] - at_depth(tab, 'pores', depth)) / r
    return zp, np.maximum(age, 0.0), mean


def older_than(tab, depth, rech_m_yr, delay, years):
    """P(T > t) at each of the years, for transit through the delay and then the saturated spread."""
    r = np.maximum(rech_m_yr, 1e-9) / YEAR_S
    ki_d = at_depth(tab, 'ki', depth)
    big_d = at_depth(tab, 'big', depth)
    t = np.maximum(np.asarray(years)[None, :] - delay[:, None], 0.0) * YEAR_S
    target = big_d[:, None] + t * r[:, None] / np.maximum(ki_d[:, None], 1e-300)
    z = np.interp(target, tab['big'], tab['z'])
    p = at_depth(tab, 'ki', z) / np.maximum(ki_d[:, None], 1e-300)
    return np.where(np.asarray(years)[None, :] <= delay[:, None], 1.0, np.clip(p, 0.0, 1.0))


# --- Flow and mixing -----------------------------------------------------------------------------------------------

def faces(s, h, k_1km):
    """Stage 2's flows between neighbouring cells (m3/s, from a to b), from its heads."""
    rows, cols = s['rows'], s['cols']
    dx, dy, dx_face, area = gt.geometry(s)
    idx = np.arange(rows * cols).reshape(rows, cols)
    a = np.concatenate([idx.ravel(), idx[:-1].ravel()])
    b = np.concatenate([np.roll(idx, -1, axis=1).ravel(), idx[1:].ravel()])
    g = np.concatenate([np.repeat(dy / dx, cols), np.repeat(dx_face / dy, cols)])
    hh = h.ravel()
    t = gt.transmissivity(s['ground'].ravel() - hh, k_1km, gu.PERM_CAP_M2, gu.PERM_EXPONENT)
    c = 2.0 * t[a] * t[b] / (t[a] + t[b]) * g
    return a, b, c * (hh[a] - hh[b])


def mixer(s, r, k_1km):
    """The mixing of recharge along stage 2's steady flow: each land cell's water is its recharge and the flows into it,
    well mixed; it leaves as seepage and as flow to its neighbours, the seas and the lakes. Returns a solver that turns
    each cell's recharge property into the property of each cell's outflow, with the flows themselves."""
    n = s['rows'] * s['cols']
    fixed = (s['sea'] | s['lake']).ravel()
    a, b, f = faces(s, r['h'], k_1km)
    src = np.where(f > 0, a, b)
    dst = np.where(f > 0, b, a)
    q = np.abs(f)
    rech = r['recharge'].ravel()
    # Water the seas or lakes give the land enters it as new water, mixed in as recharge of the receiving cell's own.
    gift = np.bincount(dst[fixed[src] & ~fixed[dst]], weights=q[fixed[src] & ~fixed[dst]], minlength=n)
    inflow = np.bincount(dst[~fixed[src] & ~fixed[dst]], weights=q[~fixed[src] & ~fixed[dst]], minlength=n)
    through = rech + gift + inflow
    free = np.flatnonzero(~fixed)
    pos = np.full(n, -1, np.int64)
    pos[free] = np.arange(free.size)
    keep = ~fixed[src] & ~fixed[dst]
    w = q[keep] / np.maximum(through[dst[keep]], 1e-300)
    m = free.size
    A = sp.identity(m, format='csc') - sp.csc_matrix((w, (pos[dst[keep]], pos[src[keep]])), shape=(m, m))
    lu = splu(A.tocsc(), permc_spec='COLAMD')
    new = (rech + gift)[free] / np.maximum(through[free], 1e-300)

    def solve(prop):
        """prop: (cells, k) recharge property on all cells; returns the (cells, k) outflow property of each cell."""
        rhs = prop[free] * new[:, None]
        out = np.zeros((n, prop.shape[1]))
        out[free] = lu.solve(rhs)
        return out
    to_fixed = fixed[dst] & ~fixed[src]
    return solve, dict(src=src[to_fixed], dst=dst[to_fixed], q=q[to_fixed], through=through, rech=rech, gift=gift)


# --- Chemistry -----------------------------------------------------------------------------------------------------
# The carbonate system in fresh water, on arrays: the carbon cycle's constants (Plummer & Busenberg 1982 for CO2's
# solubility, carbonic acid and calcite; Harned & Owen's water product) and its Davies activity coefficients.
U_LO, U_HI = math.log(1e-13), math.log(1e-2)     # [H+] brackets, pH 2 to 13
BISECTIONS = 46
SPECIES = ('ta', 'ca', 'mg', 'mono', 'dic')      # eq/kg, then mol/kg (mono in eq), per kg of water


def constants(t_k):
    t = np.asarray(t_k, dtype=float)
    lt = np.log10(t)
    return dict(
        kh=10 ** (108.3865 + 0.01985076 * t - 6919.53 / t - 40.45154 * lt + 669365.0 / t ** 2),
        k1=10 ** (-356.3094 - 0.06091964 * t + 21834.37 / t + 126.8339 * lt - 1684915.0 / t ** 2),
        k2=10 ** (-107.8871 - 0.03252849 * t + 5151.79 / t + 38.92561 * lt - 563713.9 / t ** 2),
        ksp=10 ** (-171.9065 - 0.077993 * t + 2839.319 / t + 71.595 * lt),
        kw=10 ** -(4470.99 / t - 6.0875 + 0.01706 * t),
        davies=0.4883 + 8.074e-4 * (t - 273.15))


def effective(c, w):
    """Concentration constants at the water's ionic strength: K1/g1^2, K2/g2, Kw/g1^2 and Ksp/g2^2, with g1."""
    i = 0.5 * (4.0 * w['ca'] + 4.0 * w['mg'] + w['mono'] + np.maximum(w['ta'], 0.0))
    r = np.sqrt(i)
    f = r / (1.0 + r) - 0.3 * i
    g1, g2 = 10 ** (-c['davies'] * f), 10 ** (-4.0 * c['davies'] * f)
    return dict(k1=c['k1'] / g1 ** 2, k2=c['k2'] / g2, kw=c['kw'] / g1 ** 2, ksp=c['ksp'] / g2 ** 2, kh=c['kh'], g1=g1)


def bisect(g, lo, hi):
    """The root of g, decreasing in u, between lo and hi, elementwise."""
    for _ in range(BISECTIONS):
        mid = 0.5 * (lo + hi)
        up = g(mid) > 0.0
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return 0.5 * (lo + hi)


def full(x, like):
    return np.broadcast_to(np.asarray(x, dtype=float), like.shape).copy()


def open_h(ta, co2, k):
    """[H+] of water at alkalinity ta in equilibrium with dissolved CO2 co2."""
    a, b = k['k1'] * co2 + k['kw'], 2.0 * k['k1'] * k['k2'] * co2
    g = lambda u: a * np.exp(-u) + b * np.exp(-2.0 * u) - np.exp(u) - ta
    return np.exp(bisect(g, full(U_LO, ta), full(U_HI, ta)))


def fractions(h, k):
    d = h * h + k['k1'] * h + k['k1'] * k['k2']
    return k['k1'] * h / d, k['k1'] * k['k2'] / d


def closed_h(ta, dic, k):
    """[H+] of water at alkalinity ta and dissolved inorganic carbon dic, sealed from any gas."""
    def g(u):
        h = np.exp(u)
        a1, a2 = fractions(h, k)
        return dic * (a1 + 2.0 * a2) + k['kw'] / h - h - ta
    return np.exp(bisect(g, full(U_LO, ta), full(U_HI, ta)))


def add(w, xi, shares):
    """Weathering: xi eq/kg of the rock's cations, in its shares by equivalent, each balanced by alkalinity."""
    return dict(ta=w['ta'] + xi, ca=w['ca'] + shares['ca'] * xi / 2.0, mg=w['mg'] + shares['mg'] * xi / 2.0,
                mono=w['mono'] + shares['mono'] * xi, dic=w['dic'] + 0.0 * xi)


def open_water(w, p_atm, c, omega_c, passes=3):
    """Water held open to a gas at p_atm of CO2, with calcite coming out above omega_c. Returns the water, the calcite
    laid (mol/kg) and its pH. The calcite is solved again at the activities of the water it leaves."""
    co2 = c['kh'] * p_atm
    x = np.zeros_like(w['ta'])
    k = effective(c, w)
    h0 = open_h(w['ta'], co2, k)
    over = w['ca'] * k['k1'] * k['k2'] * co2 / h0 ** 2 / k['ksp'] > omega_c
    for _ in range(passes if over.any() else 0):
        k = effective(c, dict(w, ta=w['ta'] - 2.0 * x, ca=w['ca'] - x))
        sel = lambda a: a[over] if np.ndim(a) else a
        kk = {key: sel(v) for key, v in k.items()}
        ta, ca, co2s = w['ta'][over], w['ca'][over], sel(full(co2, w['ta']))
        a, b = kk['k1'] * co2s, kk['k1'] * kk['k2'] * co2s

        def g(u):
            h = np.exp(u)
            xx = 0.5 * (ta - (a / h + 2.0 * b / h ** 2 + kk['kw'] / h - h))
            om = (ca - xx) * b / h ** 2 / kk['ksp']
            return np.where(xx < 0.0, 1.0, np.where(xx > np.minimum(ca, 0.5 * ta), -1.0, om - sel(full(omega_c,
                                                                                                        w['ta']))))
        u = bisect(g, full(U_LO, ta), full(U_HI, ta))
        h = np.exp(u)
        x[over] = np.clip(0.5 * (ta - (a / h + 2.0 * b / h ** 2 + kk['kw'] / h - h)), 0.0, np.minimum(ca, 0.5 * ta))
    out = dict(ta=w['ta'] - 2.0 * x, ca=w['ca'] - x, mg=w['mg'], mono=w['mono'])
    k = effective(c, out)
    h = open_h(out['ta'], co2, k)
    a1, a2 = fractions(h, k)
    out['dic'] = co2 + co2 * k['k1'] / h * (1.0 + k['k2'] / h)
    return out, x, -np.log10(h * k['g1'])


def closed_water(w, c, omega_c):
    """Sealed water with calcite coming out above omega_c: the water, the calcite laid (mol/kg), its pH, CO2 (atm) and
    saturation."""
    k = effective(c, w)
    h0 = closed_h(w['ta'], w['dic'], k)
    a1, a2 = fractions(h0, k)
    over = w['ca'] * w['dic'] * a2 / k['ksp'] > omega_c
    x = np.zeros_like(w['ta'])
    for _ in range(3 if over.any() else 0):
        k = effective(c, dict(w, ta=w['ta'] - 2.0 * x, ca=w['ca'] - x))
        kk = {key: (v[over] if np.ndim(v) else v) for key, v in k.items()}
        ta, ca, dic = w['ta'][over], w['ca'][over], w['dic'][over]
        oc = full(omega_c, w['ta'])[over]

        def g(u):
            h = np.exp(u)
            f1, f2 = fractions(h, kk)
            aa = f1 + 2.0 * f2
            xx = (ta - dic * aa - kk['kw'] / h + h) / (2.0 - aa)
            om = (ca - xx) * (dic - xx) * f2 / kk['ksp']
            return np.where(xx < 0.0, 1.0, np.where(xx > np.minimum(ca, dic), -1.0, om - oc))
        u = bisect(g, full(U_LO, ta), full(U_HI, ta))
        h = np.exp(u)
        f1, f2 = fractions(h, kk)
        aa = f1 + 2.0 * f2
        x[over] = np.clip((ta - dic * aa - kk['kw'] / h + h) / (2.0 - aa), 0.0, np.minimum(ca, dic))
    out = dict(ta=w['ta'] - 2.0 * x, ca=w['ca'] - x, mg=w['mg'], mono=w['mono'], dic=w['dic'] - x)
    k = effective(c, out)
    h = closed_h(out['ta'], out['dic'], k)
    a1, a2 = fractions(h, k)
    co2 = out['dic'] * (1.0 - a1 - a2)
    return out, x, -np.log10(h * k['g1']), co2 / c['kh'], out['ca'] * out['dic'] * a2 / k['ksp']


def reaction_room(w, c, shares, ph_max, omega_c):
    """The weathering (eq/kg) that takes sealed water to ph_max, calcite coming out above omega_c on the way: beyond it
    the secondary minerals hold the pH and the rock stops dissolving. Exact at the starting water's activities."""
    k = effective(c, w)
    h = 10.0 ** -ph_max / k['g1']
    a1, a2 = fractions(h, k)
    am, wm = a1 + 2.0 * a2, k['kw'] / h - h
    xi_np = w['dic'] * am + wm - w['ta']
    omega_np = (w['ca'] + shares['ca'] * np.maximum(xi_np, 0.0) / 2.0) * w['dic'] * a2 / k['ksp']
    c1 = w['ca'] + shares['ca'] / 2.0 * xi_np
    c2 = shares['ca'] / 2.0 * (2.0 - am) - 1.0
    q = omega_c * k['ksp'] / a2
    lo, hi = np.zeros_like(xi_np), w['dic'].copy()
    for _ in range(BISECTIONS):
        mid = 0.5 * (lo + hi)
        up = (c1 + c2 * mid) * (w['dic'] - mid) - q > 0.0
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    x = 0.5 * (lo + hi)
    xi_p = (w['dic'] - x) * am + wm - w['ta'] + 2.0 * x
    return np.maximum(np.where(omega_np <= omega_c, xi_np, xi_p), 0.0)


# --- The water's path ----------------------------------------------------------------------------------------------
# Soil CO2 against actual evapotranspiration, log10 pCO2 (atm) = a + b (1 - exp(-c AET)) with AET in mm a year (Brook,
# Folkoff & Box 1983).
BROOK = (-3.47, 2.09, 0.00172)
OMEGA_ROCK = 1.0                # calcite comes out of soil water and groundwater at saturation
PH_MAX = 9.5                    # design guess: basaltic groundwaters stay below about pH 9.5 (secondary minerals hold it)
TA_CEILING = 406.0 / 61.017e3   # eq/kg: the Columbia Plateau's highest bicarbonate, 406 mg/L (Steinkampf & Hearn 1996)
# Weathering in the wetted rock slows as White & Brantley's (2003) field rates do, r = r0 (t / t0)^-0.61, from an anchor
# rate r0 at a rock age t0 and temperature, carried in temperature by Dessert et al.'s (2003) exponent. The anchors are
# set in KINETICS below.
AGEING = 0.61


def soil_co2(aet_mm_yr, p_air, shift=0.0):
    a, b, c = BROOK
    return np.maximum(10.0 ** (a + shift + b * (1.0 - np.exp(-c * np.maximum(aet_mm_yr, 0.0)))), p_air)


def rate(t_wet_yr, t_k, anchor):
    """eq/kg/yr of weathering in rock wetted t_wet_yr ago, at t_k."""
    return (anchor['rate_eq_kg_yr'] * (np.maximum(t_wet_yr, 1.0) / anchor['rock_age_yr']) ** -AGEING
            * np.exp(cc.DESSERT_B * (t_k - 273.15 - anchor['temperature_c'])))


def zero(n):
    return {key: np.zeros(n) for key in SPECIES}


def path(cells, chem, state_age, ttd):
    """Each land cell's water from rain to outflow. The soil water carries the carbon cycle's catchment load at the
    soils' age, spread by temperature, open to the soil's CO2; calcite comes out above saturation (soil carbonate).
    Below the soil the water is sealed: it keeps dissolving rock at the wetted rock's rate through the unsaturated rock
    and the saturated crust, consuming the CO2 it carries, laying calcite in the rock above saturation, until pH_MAX.
    In the open-rock variant the unsaturated rock is open to pore CO2 between the soil's and the air's. Returns the soil
    water, the water at the water table, the mean outflow and the calcite laid per kg of recharge."""
    n = cells['t_k'].size
    c = constants(cells['t_k'])
    shares = chem['shares']
    p_s = soil_co2(cells['aet_mm_yr'], chem['p_air'], chem['brook_shift'])
    load = chem['soil_load_eq_kg'] * cells['soil_weight']
    soil, x_soil, _ = open_water(add(zero(n), load, shares), p_s, c, OMEGA_ROCK)
    r = rate(state_age, cells['t_k'], chem['kinetics'])
    entry, x_vadose = soil, np.zeros(n)
    sealed_years = ttd['vadose']
    if chem['open_rock_share'] > 0.0:
        p_v = chem['p_air'] + chem['open_rock_share'] * (p_s - chem['p_air'])
        deep = ttd['vadose'] > 0.0
        # Past calcite saturation each equivalent of weathering keeps at least its magnesium and alkali share of
        # alkalinity, so the ceiling is reached within TA_CEILING / (1 - calcium share) of weathering.
        xi = np.minimum(r * ttd['vadose'], 1.01 * TA_CEILING / (1.0 - shares['ca']))
        trial, x_vadose, _ = open_water(add(soil, xi, shares), p_v, c, OMEGA_ROCK)
        # Hold the water at the ceiling: find the weathering that brings it there where it would pass it.
        over = deep & (trial['ta'] > TA_CEILING)
        if over.any():
            lo, hi = np.zeros(over.sum()), xi[over]
            sub = {key: v[over] for key, v in soil.items()}
            cs = {key: (v[over] if np.ndim(v) else v) for key, v in c.items()}
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                t_ = open_water(add(sub, mid, shares), p_v[over], cs, OMEGA_ROCK)[0]['ta']
                up = t_ < TA_CEILING
                lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
            xi[over] = lo
            trial, x_vadose, _ = open_water(add(soil, xi, shares), p_v, c, OMEGA_ROCK)
        entry = {key: np.where(deep, trial[key], soil[key]) for key in SPECIES}
        x_vadose = np.where(deep, x_vadose, 0.0)
        sealed_years = np.zeros(n)
    room = reaction_room(entry, c, shares, PH_MAX, OMEGA_ROCK)
    ages = np.minimum(ttd['ages'], np.maximum(state_age - ttd['delay'], 0.0)[:, None])
    reacting = np.minimum(sealed_years[:, None] + ages, state_age)
    xi = np.minimum(r[:, None] * reacting, room[:, None])
    wide = {key: np.repeat(v[:, None], ages.shape[1], axis=1) for key, v in entry.items()}
    cw = {key: (np.repeat(v[:, None], ages.shape[1], axis=1) if np.ndim(v) else v) for key, v in c.items()}
    out, x_rock, _, _, _ = closed_water(add(wide, xi, shares), cw, OMEGA_ROCK)
    weights = -np.diff(P_EDGES)
    mean = {key: (out[key] * weights[None, :]).sum(1) for key in SPECIES}
    return dict(soil=soil, entry=entry, outflow=mean, p_soil=p_s, rate=r,
                calcite_mol_kg=dict(soil=x_soil, unsaturated_rock=x_vadose,
                                    saturated_rock=(x_rock * weights[None, :]).sum(1)))


def surface(w, t_k, p_air):
    """Spring water as it comes out (sealed) and after it meets the air (open, before any calcite)."""
    c = constants(t_k)
    _, _, ph, pco2, omega = closed_water(dict(w), c, 1e9)
    air, _, ph_air = open_water(dict(w), p_air, c, 1e9)
    k = effective(c, air)
    h = 10.0 ** -ph_air / k['g1']
    omega_air = air['ca'] * k['k1'] * k['k2'] * (c['kh'] * p_air) / h ** 2 / k['ksp']
    return dict(ph=ph, pco2_atm=pco2, omega=omega, ph_air=ph_air, omega_air=omega_air, c=c)


# --- One state -----------------------------------------------------------------------------------------------------
DEEP_WARMING_K = (10.0, 30.0)   # spring water that circulated where the mature crust is this much above the ground


def state_recharge(s, stage1, case, state, bands):
    """Stage 2's recharge for the state: the full recharge in the mature Moon; in the young Moon less stage 1's deep
    uptake under each band's land over the span ending at the state's year."""
    p = stage1['cases'][case]['params']
    full_r = gt.recharge(s, p['recharge_share'], gu.ARID_RECHARGE_SHARE)
    if state == 'mature':
        return full_r
    names = [f'{lo:g}-{hi:g}' for lo, hi in gu.BANDS_DEG]
    frac = np.zeros(len(names))
    year = gt.SPANS[state].split('-')[1]
    for i, name in enumerate(names):
        col = stage1['cases'][case]['bands'][name]['columns']['land']
        frac[i] = min(col['by_year'][year]['uptake_m_yr_since_last'] / max(col['recharge_mm_yr'] / 1e3, 1e-12), 1.0)
    return full_r * (1.0 - frac[bands])


def hydrology(s, base, stage1, case, state, h0=None):
    """Stage 2's state re-solved, with each land cell's transit spread and the mixing along its flow."""
    p = stage1['cases'][case]['params']
    k_1km = p['perm_1km']
    rech = state_recharge(s, stage1, case, state, base['bands'])
    r = gt.solve(s, rech, k_1km, gu.PERM_CAP_M2, gu.PERM_EXPONENT, h0)
    solve_mix, fx = mixer(s, r, k_1km)
    area = base['area'].ravel()
    li = np.flatnonzero(base['land'].ravel())
    depth = np.where(r['seeping'], 0.0, np.maximum(s['ground'] - r['h'], 0.0)).ravel()[li]
    r_eff = ((fx['rech'] + fx['gift']) / area * YEAR_S)[li]                 # m/yr into each land cell
    tab = profile(k_1km, base['phi0'], base['phi_e'])
    soil_t, rock_t = unsaturated_years(tab, depth, r_eff)
    zp, ages, sat_mean = saturated_spread(tab, depth, r_eff)
    delay = soil_t + rock_t
    age_yr = STATE_AGE_YR[state]
    cdf = 1.0 - older_than(tab, depth, r_eff, delay, CDF_YEARS)
    cdf[:, CDF_YEARS >= age_yr] = 1.0
    ki_d = at_depth(tab, 'ki', depth)
    gradient = p['heat_flow'] / p['k_wet']
    deep = [np.where(depth >= dk / gradient, 1.0, at_depth(tab, 'ki', np.full_like(depth, dk / gradient))
                     / np.maximum(ki_d, 1e-300)) for dk in DEEP_WARMING_K]
    mean_age = np.minimum(delay + sat_mean, age_yr)
    cols = [np.ones_like(depth)] + [cdf[:, j] for j in range(cdf.shape[1])] + deep + [mean_age]
    prop = np.zeros((area.size, len(cols)))
    prop[li] = np.stack(cols, axis=1)
    rech_full = gt.recharge(s, p['recharge_share'], gu.ARID_RECHARGE_SHARE)
    return dict(r=r, fx=fx, solve_mix=solve_mix, li=li, depth=depth, r_eff=r_eff, soil_t=soil_t, rock_t=rock_t,
                rech_full=rech_full,
                ttd=dict(vadose=rock_t, delay=delay, ages=ages), mixed=solve_mix(prop), ncdf=cdf.shape[1],
                rech=rech, age_yr=age_yr, gradient_k_m=gradient, k_1km=k_1km)


def chemistry(hy, base, chem):
    """The water's path on a solved state, mixed along its flow to each spring, lake and sea."""
    li = hy['li']
    cells = dict(t_k=base['t_k'].ravel()[li], aet_mm_yr=base['aet_mm_yr'].ravel()[li],
                 soil_weight=base['soil_weight'].ravel()[li])
    w = path(cells, dict(chem, soil_load_eq_kg=chem['soil_load_eq_kg_at'](hy['age_yr'])), hy['age_yr'], hy['ttd'])
    prop = np.zeros((base['area'].size, len(SPECIES)))
    prop[li] = np.stack([w['outflow'][key] for key in SPECIES], axis=1)
    return dict(water=w, mixed=hy['solve_mix'](prop))


# --- The setting ---------------------------------------------------------------------------------------------------

def base_setting():
    """The 2 pixel/degree grid of stage 2 with the ground's temperature and evaporation, the soil load's spread by
    temperature, and stage 1's porosity."""
    s4 = gt.setting()
    s = gt.coarsen(s4)
    clim = dict(np.load(ROOT / gt.IGNORED[0]))
    clim.pop('metadata')
    ts4 = gt.bilinear(clim['ts_k'], clim['lat'], clim['lon'], s4['lat'], s4['lon'])
    rows, cols = s['rows'], s['cols']
    t_k = ts4.reshape(rows, 2, cols, 2).mean((1, 3))
    area = gt.geometry(s)[3]
    land = ~(s['sea'] | s['lake'])
    lat2 = s['lat'][:, None] * np.ones((1, cols))
    runoff = np.maximum(s['rain'] - s['evap'], 0.0) * area * land
    b = cc.DESSERT_B
    t_ref = math.log(float((runoff * np.exp(b * (t_k - 300.0))).sum() / runoff.sum())) / b + 300.0
    stage1 = json.loads((ROOT / gt.INPUTS[1]).read_text())
    drainage = json.loads((ROOT / 'geography/results/drainage.json').read_text())
    return s, dict(discharge_m3_yr=drainage['rivers']['discharge_to_sea_m3s'] * YEAR_S, area=area, land=land, bands=gt.band_of(lat2), t_k=t_k, aet_mm_yr=s['evap'] * 1e3,
                   soil_weight=np.exp(b * (t_k - t_ref)), runoff_m3_yr=runoff, t_ref_k=t_ref,
                   phi0=stage1['setting']['porosity_phi0'], phi_e=stage1['setting']['porosity_e_folding_m']), stage1


# --- Cases ---------------------------------------------------------------------------------------------------------
# The wetted rock's weathering rate. Central: White & Brantley's (2003) field rate, R = 3.1e-13 t^-0.61 mol/m2/s with t
# in years, on Maher & Chamberlain's (2014) specific surface of 100 m2/kg, over the rock a kilogram of pore water wets at
# stage 1's surface porosity and 2,800 kg/m3, at two equivalents a mole as anorthite gives, taken at 15 C (design
# guess for the field rates' temperature). Fast: ten times the surface, about that of fine lunar soil. Slow: the
# Columbia Plateau's basalt, where groundwater gains about 2 meq/L over its mean corrected radiocarbon age of 3,830
# years in 6-17 million-year-old flows near 16 C (Steinkampf & Hearn 1996; the gain assumes recharge near 1 meq/L).
SPECIFIC_SURFACE_M2_KG = 100.0
ROCK_DENSITY = 2800.0
WHITE_BRANTLEY = 3.1e-13
COLUMBIA = dict(gain_eq_kg=2.0e-3, age_yr=3830.0, rock_age_yr=1.5e7, temperature_c=16.0)
VARIANTS = dict(central=dict(kinetics='central', open_rock=0.0), slow_rock=dict(kinetics='slow', open_rock=0.0),
                fast_rock=dict(kinetics='fast', open_rock=0.0), open_rock=dict(kinetics='central', open_rock=1.0))
STATES = gt.STATES


def kinetics(phi0):
    per_water = (1.0 - phi0) * ROCK_DENSITY / (phi0 * gt.RHO_WATER)
    field = WHITE_BRANTLEY * SPECIFIC_SURFACE_M2_KG * per_water * 2.0 * YEAR_S
    return dict(central=dict(rate_eq_kg_yr=field, rock_age_yr=1.0, temperature_c=15.0),
                fast=dict(rate_eq_kg_yr=10.0 * field, rock_age_yr=1.0, temperature_c=15.0),
                slow=dict(rate_eq_kg_yr=COLUMBIA['gain_eq_kg'] / COLUMBIA['age_yr'],
                          rock_age_yr=COLUMBIA['rock_age_yr'], temperature_c=COLUMBIA['temperature_c']),
                rock_per_kg_water=per_water)


def rock_shares(data):
    """The highland soil's weatherable cations by equivalent, as the carbon cycle takes them (Apollo 16 soil)."""
    rock = data[ROCK]['return_limb']['rock_stock_top_metre']['highland_soil']
    mol = dict(ca=rock['Ca_kg_m2'] / cc.M_CA, mg=rock['Mg_kg_m2'] / cc.M_MG,
               mono=rock['Na_kg_m2'] / cc.M_NA + rock['K_kg_m2'] / cc.M_K)
    eq = 2.0 * mol['ca'] + 2.0 * mol['mg'] + mol['mono']
    return dict(ca=2.0 * mol['ca'] / eq, mg=2.0 * mol['mg'] / eq, mono=mol['mono'] / eq)


# --- Summaries -----------------------------------------------------------------------------------------------------

def weighted(x, f):
    return float((x * f).sum() / f.sum()) if f.sum() > 0 else 0.0


def median_age(cdf_rows, f):
    cdf = (cdf_rows * f[:, None]).sum(0) / max(f.sum(), 1e-300)
    return float(10 ** np.interp(0.5, cdf, np.log10(CDF_YEARS))), cdf


def mixture(w, f):
    return {key: np.array([weighted(w[key], f)]) for key in SPECIES}


def summarise(s, base, hy, ch, chem):
    """Springs, lakes, rivers and calcite by band and Moon-wide, for one state and chemistry."""
    km3 = YEAR_S / 1e9
    kg_yr = gt.RHO_WATER * YEAR_S                       # kg a year per m3/s
    area = base['area'].ravel()
    seep = hy['r']['seep'].ravel()
    spring = hy['r']['seeping'].ravel() & (seep > 0)
    mixed, x = hy['mixed'], ch['mixed']
    ncdf = hy['ncdf']
    w_all = {key: x[:, i] for i, key in enumerate(SPECIES)}
    tk = base['t_k'].ravel()
    bands = base['bands'].ravel()
    names = [f'{lo:g}-{hi:g}' for lo, hi in gu.BANDS_DEG]
    sp = np.flatnonzero(spring)
    ws = {key: v[sp] for key, v in w_all.items()}
    surf = surface(ws, tk[sp], chem['p_air'])
    si = np.log10(np.maximum(surf['omega_air'], 1e-30))
    build = dict(tufa=si >= THRESHOLDS['tufa'], plitvice=(si >= THRESHOLDS['plitvice']) & (surf['ph_air'] > 8.0),
                 designed=si >= THRESHOLDS['designed'])
    c_sp = surf['c']
    dep_tufa = open_water(dict(ws), chem['p_air'], c_sp, 10.0 ** THRESHOLDS['tufa'])[1]
    dep_sat = open_water(dict(ws), chem['p_air'], c_sp, 1.0)[1]
    # Rivers: quick flow off the soil at the soil water's load, and the groundwater that comes back up. What soaks in
    # leaves the quick flow whether it reaches the water table or, in the young Moon, the filling crust below.
    li = hy['li']
    water = ch['water']
    quick = np.zeros(area.size)
    quick[li] = np.maximum(base['runoff_m3_yr'].ravel()[li] / YEAR_S - hy['rech_full'].ravel()[li] * area[li], 0.0)
    soil_w = {key: np.zeros(area.size) for key in SPECIES}
    for key in SPECIES:
        soil_w[key][li] = water['soil'][key]
    soil_ta = soil_w['ta']
    soil_x = np.zeros(area.size)
    soil_x[li] = water['calcite_mol_kg']['soil']
    rech_m3s = hy['fx']['rech']
    calc = {key: np.zeros(area.size) for key in water['calcite_mol_kg']}
    for key, v in water['calcite_mol_kg'].items():
        calc[key][li] = v
    src, q_fix, dst = hy['fx']['src'], hy['fx']['q'], hy['fx']['dst']
    to_sea = s['sea'].ravel()[dst]

    def block(mask_cells, mask_springs):
        f = seep[sp][mask_springs]
        if f.sum() <= 0:
            return dict(seepage_km3_yr=0.0)
        med, cdf = median_age(mixed[sp][mask_springs, 1:1 + ncdf], f)
        mix = mixture({key: v[mask_springs] for key, v in ws.items()}, f)
        ms = surface(mix, np.array([weighted(tk[sp][mask_springs], f)]), chem['p_air'])
        land_a = area[mask_cells].sum()
        out = dict(
            seepage_km3_yr=float(f.sum() * km3),
            spring_land_share=float(area[sp][mask_springs].sum() / max(land_a, 1.0)),
            median_age_yr=med,
            younger_than={f'{y:g}_yr': float(np.interp(math.log10(y), np.log10(CDF_YEARS), cdf))
                          for y in AGE_THRESHOLDS_YR},
            deeper_than_warming={f'{dk:g}_k': weighted(mixed[sp][mask_springs, 1 + ncdf + j], f)
                                 for j, dk in enumerate(DEEP_WARMING_K)},
            alkalinity_meq_kg=float(mix['ta'][0] * 1e3), calcium_mmol_kg=float(mix['ca'][0] * 1e3),
            magnesium_mmol_kg=float(mix['mg'][0] * 1e3),
            hardness_mg_l_caco3=float((mix['ca'][0] + mix['mg'][0]) * CACO3_G_MOL * 1e3),
            ph=float(ms['ph'][0]), pco2_uatm=float(ms['pco2_atm'][0] * 1e6),
            saturation_index_in_air=float(math.log10(max(ms['omega_air'][0], 1e-30))),
            building_share_of_flow={k: float(f[v[mask_springs]].sum() / f.sum()) for k, v in build.items()},
            building_share_of_spring_land={k: float(area[sp][mask_springs][v[mask_springs]].sum()
                                                    / max(area[sp][mask_springs].sum(), 1.0))
                                           for k, v in build.items()},
            deposit_kt_caco3_yr=dict(
                to_tufa_threshold=float((dep_tufa[mask_springs] * f).sum() * kg_yr * CACO3_G_MOL / 1e9),
                to_saturation=float((dep_sat[mask_springs] * f).sum() * kg_yr * CACO3_G_MOL / 1e9)))
        return out

    land = base['land'].ravel()
    by_band = {}
    for i, name in enumerate(names):
        cells = land & (bands == i)
        if not cells.any():
            continue
        row = block(cells, bands[sp] == i)
        lb = bands[li] == i
        fl = hy['r_eff'][lb] * area[li][lb]
        row.update(
            recharge_km3_yr=float((rech_m3s[li][lb]).sum() * km3),
            soil_co2_percent=float(np.median(water['p_soil'][lb]) * 100.0),
            soil_water_alkalinity_meq_kg=weighted(water['soil']['ta'][lb], fl) * 1e3,
            soil_years=weighted(hy['soil_t'][lb], fl), unsaturated_rock_years=weighted(hy['rock_t'][lb], fl),
            rock_rate_meq_kg_yr=float(np.median(water['rate'][lb]) * 1e3),
            water_table_entry=dict(alkalinity_meq_kg=weighted(water['entry']['ta'][lb], fl) * 1e3,
                                   dic_mmol_kg=weighted(water['entry']['dic'][lb], fl) * 1e3))
        by_band[name] = row
    moon = block(land, np.ones(sp.size, bool))
    # Weathering carried off the land: quick flow at the soil water's load, groundwater at its own.
    eq_quick = float((quick * soil_ta).sum() * kg_yr)
    eq_seep = float((seep[sp] * ws['ta']).sum() * kg_yr)
    eq_fixed = float((q_fix * w_all['ta'][src]).sum() * kg_yr)
    river_ta = (eq_quick + eq_seep) / max((quick.sum() + seep[sp].sum()) * kg_yr, 1.0)
    soil_mean = float((quick * soil_ta).sum() / max(quick.sum(), 1e-300))
    laid = {key: float(((rech_m3s + (quick if key == 'soil' else 0.0)) * v).sum() * kg_yr) for key, v in calc.items()}
    cations = {key: float(((quick * soil_w[key]).sum() + (seep[sp] * ws[key]).sum()
                           + (q_fix * w_all[key][src]).sum()) * kg_yr) for key in ('ca', 'mg', 'mono')}
    discharge = base['discharge_m3_yr']
    law_total = chem['soil_load_eq_kg_at'](hy['age_yr']) * discharge * gt.RHO_WATER
    moon.update(
        drainage_discharge_km3_yr=discharge / 1e9, land_runoff_km3_yr=float(base['runoff_m3_yr'].sum() / 1e9),
        recharge_km3_yr=float(rech_m3s[land].sum() * km3), quick_flow_km3_yr=float(quick.sum() * km3),
        into_filling_crust_km3_yr=float(((hy['rech_full'] - hy['rech']).ravel() * area)[land].sum() * km3),
        into_lakes_km3_yr=float(q_fix[~to_sea].sum() * km3), into_seas_km3_yr=float(q_fix[to_sea].sum() * km3),
        weathering=dict(
            quick_flow_mol_yr=eq_quick, springs_mol_yr=eq_seep, into_lakes_and_seas_mol_yr=eq_fixed,
            total_gt_c_yr=(eq_quick + eq_seep + eq_fixed) * cc.M_C / 1e15,
            groundwater_share=(eq_seep + eq_fixed) / max(eq_quick + eq_seep + eq_fixed, 1e-300),
            carbon_cycle_catchment_mol_yr=law_total,
            over_carbon_cycle_catchment=(eq_quick + eq_seep + eq_fixed) / max(law_total, 1e-300),
            rivers_to_sea_alkalinity_meq_kg=(eq_quick + eq_seep + eq_fixed) / (discharge * gt.RHO_WATER) * 1e3,
            land_water_alkalinity_meq_kg=river_ta * 1e3, soil_water_alkalinity_meq_kg=soil_mean * 1e3,
            springs_over_quick_flow=float(moon['alkalinity_meq_kg'] / max(soil_mean * 1e3, 1e-12)),
            cations_mol_yr=cations),
        calcite_laid_gt_c_yr={key: v * cc.M_C / 1e15 for key, v in laid.items()},
        calcite_laid_mt_caco3_yr={key: v * CACO3_G_MOL / 1e12 for key, v in laid.items()})
    return dict(by_band=by_band, moon=moon)


# --- The product ---------------------------------------------------------------------------------------------------
EVIDENCE = ('A calculation on committed products: stage 2\'s water table re-solved on the atlas at 15 km, a generalized '
            'Vogel model of each cell\'s transit times through stage 1\'s permeability and porosity, unit-gradient '
            'drainage through the soil and the rock above the water table, and the recharge mixed along the steady '
            'flow; the carbonate system at each cell\'s ground temperature with the carbon cycle\'s constants, its '
            'catchment load in the soil water, Brook et al.\'s soil CO2, field weathering rates that age as White & '
            'Brantley\'s and calcite at saturation. A runnable calculation checked against analytic solutions and the '
            'carbon cycle\'s solver. Its reaction rates, the rock\'s reactive surface and the pH at which the rock stops '
            'dissolving are Earth measurements and design guesses; no lunar groundwater has been measured.')
READING_RULE = ('Years from the moment the planned climate and the water are in place; the mature Moon is taken at a '
                'million years, when the deep crust has filled. Spring values are means over the springs\' flow, the '
                'water mixed before it meets the air; saturation indices are of calcite once the water has come to the '
                'air\'s CO2, before any calcite comes out. Alkalinity in meq/kg, calcium and magnesium in mmol/kg, '
                'hardness as mg/L of CaCO3, CO2 in microatmospheres. Ages are the spread of times from rain to outflow '
                'at the springs: the median, and the shares younger than each age. Building shares are of the springs\' '
                'flow and of the land where groundwater comes out. Deposits are the calcite the springs\' water would lay '
                'after it meets the air, down to the tufa threshold or to saturation. Calcite laid is what comes out '
                'in the soil, the rock above the water table and the saturated rock a year, as carbon and as CaCO3.')


def cell_median_age(rows):
    """Each cell's median age from its spread at the quarter decades."""
    k = np.clip((rows >= 0.5).argmax(1), 1, rows.shape[1] - 1)
    i = np.arange(rows.shape[0])
    lo, hi = rows[i, k - 1], rows[i, k]
    t = np.clip((0.5 - lo) / np.maximum(hi - lo, 1e-12), 0.0, 1.0)
    return 10 ** (np.log10(CDF_YEARS)[k - 1] + t * 0.25)


def results():
    data = {}
    for path, schema in INPUTS.items():
        data[path] = json.loads((ROOT / path).read_text())
        if schema and data[path].get('schema') != schema:
            raise ValueError(f'{path}: unexpected schema {data[path].get("schema")}')
    s, base, stage1 = base_setting()
    log('setting')
    kin = kinetics(base['phi0'])
    shares = rock_shares(data)
    p_air = data[AIR]['design_air']['co2_pa'] / ATM_PA
    law = cc.DESSERT_A * math.exp(cc.DESSERT_B * (cc.LAND_AIR_K - 273.15)) / 1e6
    corner = dict(fresh_factor=cc.FRESH_FACTOR[1], soil_age_yr=cc.SOIL_AGE_AT_START_YR[1])

    def chem_for(variant):
        v = VARIANTS[variant]
        return dict(shares=shares, p_air=p_air, brook_shift=0.0,
                    soil_load_eq_kg_at=lambda t: law * cc.weathering_factor(t, corner),
                    kinetics=kin[v['kinetics']], open_rock_share=v['open_rock'])
    cases, grids = {}, {}
    for case in gu.CASES:
        h0, rows = None, {}
        for state in STATES:
            hy = hydrology(s, base, stage1, case, state, h0)
            if state == 'mature':
                h0 = hy['r']['h']
            rows[state] = dict(soil_load_factor=cc.weathering_factor(STATE_AGE_YR[state], corner),
                               rock_age_yr=STATE_AGE_YR[state], variants={})
            for variant in (VARIANTS if case == 'central' else ('central',)):
                chem = chem_for(variant)
                ch = chemistry(hy, base, chem)
                rows[state]['variants'][variant] = summarise(s, base, hy, ch, chem)
                log(f'{case} {state} {variant}')
                if case == 'central' and variant == 'central':
                    seep = hy['r']['seep'].ravel()
                    spring = hy['r']['seeping'].ravel() & (seep > 0)
                    w = {key: ch['mixed'][:, i] for i, key in enumerate(SPECIES)}
                    sp = np.flatnonzero(spring)
                    surf = surface({key: v[sp] for key, v in w.items()}, base['t_k'].ravel()[sp], p_air)
                    shape = (s['rows'], s['cols'])
                    age = np.full(seep.size, np.nan)
                    age[sp] = cell_median_age(hy['mixed'][sp, 1:1 + hy['ncdf']])
                    hard = np.full(seep.size, np.nan)
                    hard[sp] = (w['ca'][sp] + w['mg'][sp]) * CACO3_G_MOL * 1e3
                    si = np.full(seep.size, np.nan)
                    si[sp] = np.log10(np.maximum(surf['omega_air'], 1e-30))
                    ph = np.full(seep.size, np.nan)
                    ph[sp] = surf['ph']
                    grids.update({f'spring_median_age_yr_{state}': age.reshape(shape).astype(np.float32),
                                  f'spring_hardness_mg_l_{state}': hard.reshape(shape).astype(np.float32),
                                  f'spring_saturation_index_in_air_{state}': si.reshape(shape).astype(np.float32),
                                  f'spring_ph_{state}': ph.reshape(shape).astype(np.float32)})
        cases[case] = dict(perm_1km_m2=stage1['cases'][case]['params']['perm_1km'],
                           recharge_share=stage1['cases'][case]['params']['recharge_share'], states=rows)
    grids.update(lat_deg=s['lat'].astype(np.float32), lon_deg=s['lon'].astype(np.float32))
    setting = dict(
        grid_pixels_per_degree=gt.PPD, air_co2_atm=p_air, ground_reference_k=base['t_ref_k'], cation_shares=shares,
        catchment_load_meq_kg=law * 1e3, soil_co2_brook=list(BROOK), calcite_saturation_in_rock=OMEGA_ROCK,
        carbon_cycle=dict(dessert_a=cc.DESSERT_A, dessert_b=cc.DESSERT_B, land_air_k=cc.LAND_AIR_K,
                          fresh_factor=corner['fresh_factor'], soil_age_yr=corner['soil_age_yr'],
                          age_exponent=cc.AGE_EXPONENT),
        ph_where_rock_stops=PH_MAX, alkalinity_ceiling_meq_kg=TA_CEILING * 1e3,
        rock_per_kg_of_water_kg=kin['rock_per_kg_water'],
        rates_eq_kg_yr={k: v for k, v in kin.items() if k != 'rock_per_kg_water'}, ageing_exponent=AGEING,
        soil=dict(depth_m=SOIL_DEPTH_M, porosity=SOIL_POROSITY, k_m_s=SOIL_K_M_S, b=B_SOIL), rock_b=B_ROCK,
        thresholds=THRESHOLDS, deep_warming_k=list(DEEP_WARMING_K), state_ages_yr=STATE_AGE_YR,
        variants=VARIANTS)
    return dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, setting=setting, cases=cases), grids


def main(argv=None) -> int:
    from geography.groundwater_uptake import rnd
    product, grids = results()
    product = rnd(product)
    product['producer'] = dict(
        study='research/studies/springs', files={p: sha256(ROOT / p) for p in CODE},
        inputs={p: sha256(ROOT / p) for p in INPUTS}, ignored_inputs={p: sha256(ROOT / p) for p in IGNORED})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    GRIDS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GRIDS, metadata=json.dumps(dict(schema=SCHEMA + '-grids', reading_rule=READING_RULE)), **grids)
    for state, row in product['cases']['central']['states'].items():
        m = row['variants']['central']['moon']
        print(state, 'median age', m['median_age_yr'], 'hardness', m['hardness_mg_l_caco3'], 'SI air',
              m['saturation_index_in_air'], 'groundwater share', m['weathering']['groundwater_share'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
