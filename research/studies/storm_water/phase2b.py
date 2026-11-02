"""Phase 2b of the storm-water study: the pools through the lunar night, the leverage on PlaSim's land scheme, the
biome maps' wet-patch layer, and what sustains or wears each surface.

The hollows of phase 2a hold the routed water: geography's closed depressions on the 1.9 km grid, each with an
area-volume curve from the depression tree, and the crater bowls below the grid as paraboloids, Pike's fresh form at
half its depth. Each pool runs through the ten design years in half-mean steps with its routed inflow, the rain on its
own surface, its evaporation, its seepage and its spill. The evaporation comes from an energy balance of the pool and
its floor under the design years' radiation, air, humidity and wind by tenth: by day it is Penman's evaporation at
water's albedo, at night the pool's stored heat and the night air's sensible heat set the limit. Seepage passes a clay
seal over the subsoil at lunar gravity, either everywhere (clay-floored) or below a stage of 1 m with sandy banks
above it (sandy-banked, the HAPEX-Sahel pools' form). A depression that fills spills into the one below it in the
tree. The slopes' soil loss follows WEPP's interrill and rill detachment with the lunar drops' energy and the flow's
shear at lunar gravity, against the project's soil formation; the salts follow each pool's solutes to its floor.
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2a as P2
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R1
from research.studies.storm_water import wetting as W
from shared.constants import (JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, STEFAN_BOLTZMANN,
                              SYNODIC_MONTH_DAYS)

STATES = P2.STATES
SEEPAGE = ('clay_floored', 'sandy_banked')
MEAN_H = R1.MEAN_DAYS * 24.0
LATENT, CP = EX.LATENT, EX.CP
R_DRY, R_VAP = 287.05, 461.5
RHO_W, C_W = 1000.0, 4186.0
KARMAN = 0.4
# The pool and its air (each value a design choice, an Earth value or a physical property, as named).
WATER_ALBEDO = 0.07                 # open water, the plan's reading
WATER_EMISSIVITY = 0.97             # a property of water
FLOOR_HEAT = dict(depth_m=0.3, rho_c_j_m3_k=2.0e6)   # design: the wet floor's soil joins the pool's heat store
CHARNOCK = 0.011                    # Earth value; the roughness z0 = a u*^2 / g + 0.11 nu / u* carries the gravity
NU_AIR = 1.4e-5                     # m2/s, air near 300 K at the drylands' 1.16 atm
Z0_SCALAR_M = 1.0e-4                # Earth value: heat and vapour roughness over water, set by the viscous sublayer
Z_WIND_M, Z_AIR_M = 10.0, 2.0       # PlaSim's 10 m wind and 2 m air; its lowest level's humidity (about 309 m) is
                                    # taken as the surface layer's (a limit)
D_REF = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])       # mean depths (m) at which the energy balance runs
STEPS_PER_MEAN = 2
WET_M, DEEP_M = 0.01, 0.10          # a pool holds water at 1 cm at its deepest point; the plan's measure is 10 cm
DAY_FRACTION = np.array([0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0])    # share of each tenth in daylight
CHECKPOINTS = dict(sunset=7, sunrise=1, before_storms=3)      # the tenth at whose end each check falls
SPINUP_MEANS = 120                  # the first design year fills the pools; statistics use years 2-10
# Seepage (design choices at stated Earth values; conductivities times g, suction heads over g).
SEAL = dict(k_earth_mm_h=0.01, thickness_m=0.05, subsoil='sandy_loam')
# A 5 cm clay seal at 2.8e-9 m/s over a sandy loam passes about 5.5 mm a day under 0.5 m of water on Earth, the
# order of the HAPEX-Sahel plateau pools, which lose under half their water to seepage against a few mm a day of
# evaporation (Desconnets et al. 1997, abstract, via the plan); a design choice, not a measured lunar floor.
BANK = dict(k_earth_mm_h=RT.BED['k_eff_in_h'] * RT.INCH_MM, clogged_stage_m=1.0)
# Sandy banks: the channels' sand (Lane's class 2, 3 in/h) above a clogged floor 1 m deep, the lower end of the
# HAPEX-Sahel pools' 1-2 m threshold above which most stored water infiltrates (Desconnets et al. 1997).
# Solutes (stated design inputs).
SALTS = dict(rain_mg_l=2.0, dust_soluble_g_m2_yr=0.5, soda_share=0.02, saline_kg_m2=0.31)
# Rain from fresh seas carries no sea salt: 2 mg/L of continental solutes (design); the dust's soluble part, 0.5 g/m2
# a year, dissolved in the year's rain (design); sodium and potassium 2% of the solutes, the weathering share of the
# sea chemistry study, kept as soda once calcite, magnesium clays and opal leave; 0.31 kg/m2 of soluble salt in the
# top 0.3 m (4 dS/m in a saturation extract at 0.4 water content, 0.64 g/L per dS/m), the threshold of a saline soil.
# The slopes' wear: WEPP's baseline erodibilities for soils with 30% sand or more, Ki = 2,728,000 + 192,100 VFS
# (kg s/m4), Kr = 0.00197 + 0.00030 VFS + 0.03863 exp(-1.84 OM) (s/m), tau_c = 2.67 + 0.065 clay - 0.058 VFS (Pa),
# percentages (Alberts et al. 1995, read in a search summary of WEPP's documentation).
TEXTURES = dict(sandy_loam=dict(sand=65.0, clay=10.0, vfs=15.0, om=1.0), loam=dict(sand=40.0, clay=20.0, vfs=10.0, om=1.0))
BIOCRUST = dict(erodibility_factor=0.1, critical_shear_factor=3.0)    # design: crusts bind the surface (recalled)
RILL = dict(spacing_m=1.0, width_a=1.13, width_b=0.303)               # WEPP's rill spacing and width (recalled)
BULK_DENSITY = 1400.0
SOIL_FORMATION_MM_YR = (0.017, 0.036)  # natural soil formation (Montgomery 2007), as biosphere/ecology/people.py has it
# The Earth checks' rules, set before each was run.
MEAD_RULE = dict(tolerance=0.25, note='the bulk transfer at Earth\'s gravity with the measured water-surface temperature, '
                 'air vapour pressure and the period\'s mean wind gives the period\'s latent heat within 25% of the '
                 'Bowen-ratio energy budget\'s 164 W/m2 (Moreo and Swancar 2013, Table 4)')
NIGHT_RULE = dict(day_over_pt=(0.8, 2.0), night_mm=(2.5, 30.0), depth_m=0.3,
                  note='over the semi-arid cells, the day\'s evaporation of a 0.3 m pool lies within 0.8-2.0 times '
                       'Priestley and Taylor\'s demand at water\'s albedo (the plan: small pools in dry air up to about '
                       'twice it), and the night\'s within 2.5-30 mm, a factor of two around the plan\'s 5-15 mm')
# Lake Mead, March 2010-February 2012 (Moreo and Swancar 2013, USGS SIR 2013-5229, Table 4, read at the source):
# monthly air and water-surface temperature (C), vapour pressure at air temperature (kPa) and the Bowen-ratio energy
# budget's latent heat (W/m2); the daily mean wind over the record, 3.74 m/s, from sensors 3.18-9.66 m above water.
MEAD = dict(ta=[15.6, 18.8, 22.8, 31.0, 35.3, 33.0, 29.6, 23.3, 16.3, 12.9, 11.1, 12.0,
                16.4, 19.4, 22.9, 29.6, 32.7, 34.5, 29.6, 22.5, 15.7, 11.1, 12.3, 13.6],
            to=[13.6, 15.4, 18.5, 23.1, 27.6, 27.5, 26.6, 23.8, 19.2, 15.4, 13.4, 12.6,
                13.8, 16.0, 18.8, 22.8, 26.9, 28.6, 28.0, 23.3, 18.7, 15.6, 13.9, 13.3],
            ea=[0.79, 0.80, 0.84, 1.09, 1.55, 1.55, 1.32, 1.28, 0.79, 0.92, 0.72, 0.64,
                0.82, 0.82, 0.89, 1.02, 1.65, 1.56, 1.60, 1.14, 0.81, 0.64, 0.64, 0.68],
            le=[150, 131, 220, 189, 291, 235, 165, 189, 193, 105, 77, 101,
                91, 166, 210, 219, 161, 230, 156, 167, 148, 134, 102, 109],
            le_period=164.0, wind_m_s=3.74, height_m=6.0, pressure_pa=97000.0)


def esat_pa(t_k):
    t = t_k - 273.15
    return 610.8 * np.exp(17.27 * t / (t + 237.3))


def qsat(t_k, p_pa):
    e = esat_pa(t_k)
    return 0.622 * e / (p_pa - 0.378 * e)


def transfer(wind, gravity, z_wind=Z_WIND_M, z_air=Z_AIR_M):
    """Neutral bulk transfer coefficient for heat and vapour over water at a gravity (in Earth's): the log law with
    Charnock's roughness at that gravity plus the smooth-flow term, and a viscous scalar roughness."""
    u = np.maximum(np.asarray(wind, float), 0.3)
    g = gravity * STANDARD_GRAVITY
    z0 = np.full_like(u, 2e-4)
    for _ in range(6):
        ustar = KARMAN * u / np.log(z_wind / z0)
        z0 = CHARNOCK * ustar ** 2 / g + 0.11 * NU_AIR / ustar
    return KARMAN ** 2 / (np.log(z_wind / z0) * np.log(z_air / Z0_SCALAR_M)), z0


def mean_by_tenth(a, tenth):
    out = np.zeros((10,) + a.shape[1:])
    for k in range(10):
        m = tenth == k
        out[k] = (a * m).sum(axis=0) / np.maximum(m.sum(axis=0), 1)
    return out


def forcing(x):
    """Each dryland cell's air and radiation by tenth (tenth, cell): sunlight down from the net and the albedo,
    longwave down from the net and the ground's emission, 2 m air, the lowest level's humidity, the pressure and the
    10 m wind speed of each 3-day mean."""
    t = x['tenth']
    alb = np.clip(x['alb'].astype(float), 0.0, 0.95)
    sw = x['rss'].astype(float) / (1.0 - alb)
    lw = x['rls'].astype(float) + STEFAN_BOLTZMANN * x['ts'].astype(float) ** 4
    wind = np.hypot(x['uas'].astype(float), x['vas'].astype(float))
    f = dict(sw_in=sw, lw_in=lw, ta=x['tas'].astype(float), qa=x['hus_low'].astype(float),
             p=x['ps'].astype(float) * 100.0, wind=wind, rls=x['rls'].astype(float))
    return {k: mean_by_tenth(v, t) for k, v in f.items()}


def fluxes(tw, f, ce, k=None):
    """Net radiation, latent and sensible heat (W/m2, upward positive for the turbulent fluxes) of water at tw."""
    sl = slice(None) if k is None else k
    p, ta, qa, u = f['p'][sl], f['ta'][sl], f['qa'][sl], f['wind'][sl]
    rho = p / (R_DRY * ta * (1.0 + 0.61 * qa))
    rn = (1.0 - WATER_ALBEDO) * f['sw_in'][sl] + WATER_EMISSIVITY * (f['lw_in'][sl] - STEFAN_BOLTZMANN * tw ** 4)
    le = rho * LATENT * ce * u * (qsat(tw, p) - qa)
    h = rho * CP * ce * u * (tw - ta)
    return rn, le, h


def penman(f, gravity, factor=1.0):
    """Penman's evaporation in its exact form: the water temperature at which net radiation balances the latent and
    sensible heat with no storage, by tenth and cell (mm/h), with that temperature."""
    ce, _ = transfer(f['wind'], gravity)
    ce = ce * factor
    tw = f['ta'].copy()
    for _ in range(40):
        rn, le, h = fluxes(tw, f, ce)
        r = rn - le - h
        d = 0.05
        rn2, le2, h2 = fluxes(tw + d, f, ce)
        tw = tw + r / np.maximum((r - (rn2 - le2 - h2)) / d, 1.0)
    rn, le, h = fluxes(tw, f, ce)
    return le / LATENT * 3600.0, tw


def energy_balance(f, gravity, depths=D_REF, cycles=3, dt_h=1.0, factor=1.0):
    """The pool's temperature through the lunar cycle at each mean depth, the forcing constant within each tenth,
    with the floor's soil in the heat store: mean evaporation (mm/h) and water temperature (K) by depth, tenth and
    cell from the last cycle."""
    ce, _ = transfer(f['wind'], gravity)
    ce = ce * factor
    nd, nc = depths.size, f['ta'].shape[1]
    heat = (RHO_W * C_W * depths + FLOOR_HEAT['depth_m'] * FLOOR_HEAT['rho_c_j_m3_k'])[:, None]
    steps = int(round(MEAN_H / dt_h))
    tw = np.repeat(f['ta'][4][None, :], nd, axis=0)
    e_out, t_out = np.zeros((nd, 10, nc)), np.zeros((nd, 10, nc))
    dt = dt_h * 3600.0
    for cyc in range(cycles):
        for k in range(10):
            fk = {q: v[k][None, :] for q, v in f.items()}
            cek = ce[k][None, :]
            esum, tsum = np.zeros((nd, nc)), np.zeros((nd, nc))
            for _ in range(steps):
                rn, le, h = fluxes(tw, fk, cek)
                r = rn - le - h
                rn2, le2, h2 = fluxes(tw + 0.05, fk, cek)
                slope = np.maximum((r - (rn2 - le2 - h2)) / 0.05, 0.0)
                tw = tw + dt * r / (heat + dt * slope)
                esum += le
                tsum += tw
            if cyc == cycles - 1:
                e_out[:, k] = esum / steps / LATENT * 3600.0
                t_out[:, k] = tsum / steps
    return e_out.transpose(2, 0, 1), t_out.transpose(2, 0, 1)                 # (cell, depth, tenth)


def pt_water(f):
    """Priestley and Taylor's demand at water's albedo by tenth and cell (mm/h): the sunlight absorbed at water's
    albedo with the land's net longwave, as the plan reads it, in the extract's formula."""
    net = (1.0 - WATER_ALBEDO) * f['sw_in'] + f['rls']
    return np.maximum(EX.demand(f['ta'], f['p'] / 100.0, net), 0.0) * 1000.0 * 3600.0


def seal_seepage(depth_m, gravity, seal=SEAL):
    """Seepage (m/s) through a clay seal under depth_m of water over unsaturated subsoil draining under gravity:
    K_seal (d + L + psi) / L = K_soil(psi), Brooks and Corey's K(psi) = Ks (hb / psi)^(2 + 3 lam), solved for the
    suction head psi below the seal."""
    c = W.curve(seal['subsoil'], gravity)
    kc = seal['k_earth_mm_h'] * gravity / 1000.0 / 3600.0
    ks, hb, eta = c['ks'] / 1000.0 / 3600.0, c['hb'] / 1000.0, 2.0 + 3.0 * c['lam']
    d = np.asarray(depth_m, float)
    lo, hi = np.full_like(d, hb), np.full_like(d, 1e4)
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        f = kc * (d + seal['thickness_m'] + mid) / seal['thickness_m'] - ks * (hb / mid) ** eta
        lo, hi = np.where(f < 0, mid, lo), np.where(f < 0, hi, mid)
    psi = np.sqrt(lo * hi)
    return kc * (d + seal['thickness_m'] + psi) / seal['thickness_m']


def bank_seepage(gravity):
    return BANK['k_earth_mm_h'] * gravity / 1000.0 / 3600.0


# ---------------------------------------------------------------------------------------------------------------------
# The pools' geometry: every pool as a piecewise curve of area against stage, linear between breakpoints.

def depression_curves(z, a, n_below, pass_m, vmax):
    """Breakpoints (stage, area, volume) of each depression from its lowest cells: the area rises from nothing at the
    floor cell's height to each cell's full area at the next cell's height, and to all of them at the pass."""
    n, k = z.shape
    s = np.full((n, k + 1), np.nan)
    A = np.zeros((n, k + 1))
    z0 = z[:, :1]
    s[:, :k] = z - z0
    cum = np.cumsum(np.nan_to_num(a), axis=1)
    A[:, 1:k] = cum[:, :k - 1]
    nb = np.minimum(n_below, k)
    rows = np.arange(n)
    full = n_below <= k
    # the pass closes the curve where every cell below it is kept; beyond k cells the walls stand vertical
    s[rows, nb] = np.where(full, pass_m - z0[:, 0], s[rows, np.maximum(nb - 1, 0)] + 1e-3)
    A[rows, nb] = cum[rows, np.maximum(nb - 1, 0)]
    for j in range(k + 1):
        pad = j > nb
        s[pad, j] = s[pad, nb[pad]] if pad.any() else s[pad, j]
        A[pad, j] = A[pad, nb[pad]] if pad.any() else A[pad, j]
    s = np.maximum.accumulate(np.nan_to_num(s, nan=0.0), axis=1)
    ds = np.diff(s, axis=1)
    V = np.zeros_like(A)
    V[:, 1:] = np.cumsum(0.5 * (A[:, 1:] + A[:, :-1]) * ds, axis=1)
    top_v = np.where(full, np.maximum(vmax, V[rows, nb]), V[rows, nb] + A[rows, nb] * 1e3)
    return s, A, V, top_v


def bowl_curves(diameter_m):
    """A paraboloid at half Pike's fresh depth: one segment from nothing at the floor to the rim's area at the rim."""
    d = np.asarray(diameter_m, float)
    depth = RT.BOWL_DEPTH_RATIO * d
    area = math.pi * d ** 2 / 4.0
    s = np.stack([np.zeros_like(d), depth], axis=1)
    A = np.stack([np.zeros_like(d), area], axis=1)
    V = np.stack([np.zeros_like(d), 0.5 * area * depth], axis=1)
    return s, A, V, V[:, 1].copy()


def pad_curves(curves, width):
    out = []
    for s, A, V, top in curves:
        k = s.shape[1]
        if k < width:
            s = np.concatenate([s, np.repeat(s[:, -1:], width - k, axis=1)], axis=1)
            A = np.concatenate([A, np.repeat(A[:, -1:], width - k, axis=1)], axis=1)
            V = np.concatenate([V, np.repeat(V[:, -1:], width - k, axis=1)], axis=1)
        out.append((s, A, V, top))
    return [np.concatenate(v) for v in zip(*out)]


def on_curve(vol, s, A, V):
    """Stage and area of each pool at its volume, from its breakpoints (area linear in stage within a segment)."""
    j = np.clip((vol[:, None] > V).sum(axis=1) - 1, 0, V.shape[1] - 2)
    r = np.arange(vol.size)
    s0, a0, v0 = s[r, j], A[r, j], V[r, j]
    ds = np.maximum(s[r, j + 1] - s0, 1e-9)
    g = (A[r, j + 1] - a0) / ds
    dv = np.maximum(vol - v0, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        u = np.where(g > 1e-12, (np.sqrt(a0 ** 2 + 2.0 * g * dv) - a0) / np.maximum(g, 1e-12), dv / np.maximum(a0, 1e-9))
    u = np.minimum(u, ds)
    return s0 + u, a0 + g * u


def area_at_stage(stage, s, A):
    j = np.clip((stage[:, None] >= s).sum(axis=1) - 1, 0, s.shape[1] - 2)
    r = np.arange(stage.size)
    ds = np.maximum(s[r, j + 1] - s[r, j], 1e-9)
    w = np.clip((stage - s[r, j]) / ds, 0.0, 1.0)
    return A[r, j] + w * (A[r, j + 1] - A[r, j])


# ---------------------------------------------------------------------------------------------------------------------
# The pools through the design years.

def simulate(pools, x, evap_tab, seep_case, gravity, off, solute_mg_l):
    """Every pool through the means in STEPS_PER_MEAN steps each: inflow spread through each mean, the rain on the
    pool's surface, evaporation at the pool's mean depth in its cell's tenth, seepage, and spill to the depression it
    drains into. The volume advances implicitly in sqrt(V), exact for a paraboloid drying under constant losses."""
    n = pools['cell'].size
    cell, s, A, V, top = pools['cell'], pools['s'], pools['A'], pools['V'], pools['top']
    c1 = np.sqrt(2.0 * np.maximum(A[:, 1], 1e-9) / np.maximum(s[:, 1], 1e-9))      # A = c1 sqrt(V) in the first segment
    clog = area_at_stage(np.full(n, BANK['clogged_stage_m']), s, A) if seep_case == 'sandy_banked' else None
    q_bank = bank_seepage(gravity)
    depth_tab = np.linspace(0.0, 30.0, 601)
    seal_tab = seal_seepage(depth_tab, gravity)
    tenth = x['tenth'].astype(np.int64)
    nm = tenth.shape[0]
    years = nm * R1.MEAN_DAYS / JULIAN_YEAR_DAYS                                   # Julian years in the record
    tot = np.stack([(off * (tenth == k)).sum(axis=0) for k in range(10)])
    share = off / np.maximum(np.take_along_axis(tot, tenth, axis=0), 1e-30)        # each mean's part of its tenth
    rain = x['pr'].astype(float) * R1.MM_PER_MEAN / 1000.0                         # m per mean
    with np.errstate(divide='ignore', invalid='ignore'):
        direct = np.clip(1.0 - np.where(rain > 0, off / 1000.0 / rain, 0.0), 0.0, 1.0)   # rain on a pool not yet counted
    parent = pools['parent_pos']
    dt = MEAN_H / STEPS_PER_MEAN * 3600.0
    ln_d = np.log(D_REF)
    idx = np.arange(n)
    vol, stage, area, spill_in, salt = (np.zeros(n) for _ in range(5))
    nd = D_REF.size
    ev = np.ascontiguousarray(evap_tab.transpose(0, 2, 1)).reshape(-1) / 1000.0 / 3600.0    # (cell, tenth, depth) m/s
    g1 = A[:, 1] / np.maximum(s[:, 1], 1e-9)                                       # area per metre of stage, first segment
    v1 = V[:, 1]

    def geometry(v):
        st_ = np.sqrt(2.0 * v / np.maximum(g1, 1e-12))
        ar = g1 * st_
        big = v > v1
        if big.any():
            st_[big], ar[big] = on_curve(v[big], s[big], A[big], V[big])
        return st_, ar
    z = lambda: np.zeros(n)
    st = dict(evap=np.zeros((n, 10)), seep=z(), spill_out=z(), inflow=z(), rain=z(), wet=z(), deep=z(), stage_sum=z(),
              stage_max=z(), area=np.zeros((n, 10)), steps_by_tenth=np.zeros(10), dry_ever=np.zeros(n, bool),
              deposit=z(), salt_in=z(), salt_out=z(), checks={c: z() for c in CHECKPOINTS},
              n_at={c: z() for c in CHECKPOINTS})
    for m in range(nm):
        kk = tenth[m, cell]
        q_step = pools['inflow'][idx, kk] * years * share[m, cell] / STEPS_PER_MEAN   # m3 per step
        r_step = rain[m, cell] * direct[m, cell] / STEPS_PER_MEAN                      # m per step
        rec = m >= SPINUP_MEANS
        if m == SPINUP_MEANS:
            st['storage_start'] = vol.copy()
        base = (cell * 10 + kk) * nd
        for _ in range(STEPS_PER_MEAN):
            dmean = np.where(area > 0, vol / np.maximum(area, 1e-9), 0.0)
            ld = np.clip(np.log(np.maximum(dmean, D_REF[0])), ln_d[0], ln_d[-1])
            j = np.clip(np.searchsorted(ln_d, ld) - 1, 0, D_REF.size - 2)
            w = (ld - ln_d[j]) / (ln_d[j + 1] - ln_d[j])
            e = (1.0 - w) * ev[base + j] + w * ev[base + j + 1]
            q_seal = np.interp(dmean, depth_tab, seal_tab)
            if clog is None:
                seep = q_seal
            else:
                seep = np.where(area > clog, (q_seal * clog + q_bank * (area - clog)) / np.maximum(area, 1e-9), q_seal)
            loss = e + seep
            water_in = q_step + r_step * area
            gain = water_in + spill_in
            spill_in = np.zeros(n)
            c = np.where(vol > 0, area / np.sqrt(np.maximum(vol, 1e-12)), c1)
            b = dt * loss * c - 2.0 * np.sqrt(vol)
            new = ((-b + np.sqrt(b * b + 8.0 * gain)) / 4.0) ** 2
            lost = vol + gain - new
            with np.errstate(divide='ignore', invalid='ignore'):
                fe = np.where(np.abs(loss) > 0, e / loss, 0.0)
            over = np.maximum(new - top, 0.0)
            new -= over
            has = over > 0
            if has.any():
                tgt, amt = parent[has], over[has]
                np.add.at(spill_in, tgt[tgt >= 0], amt[tgt >= 0])
            # solutes come with the water, leave with the seepage and the spill at the pool's concentration, and
            # settle on the floor when the pool dries
            sin = water_in * solute_mg_l[cell]                                         # g
            salt += sin
            conc = salt / np.maximum(vol + gain, 1e-12)
            sout = np.minimum(conc * (np.maximum(lost * (1.0 - fe), 0.0) + over), salt)
            salt -= sout
            vol = new
            stage, area = geometry(vol)
            wet = stage >= WET_M
            dry = vol <= 1e-9
            if rec:
                st['evap'][idx, kk] += lost * fe
                st['seep'] += lost * (1.0 - fe)
                st['spill_out'] += np.where(parent >= 0, 0.0, over)
                st['inflow'] += q_step
                st['rain'] += r_step * area
                st['salt_in'] += sin
                st['salt_out'] += sout
                st['deposit'] += np.where(dry, salt, 0.0)
                st['wet'] += wet
                st['deep'] += stage >= DEEP_M
                st['stage_sum'] += np.where(wet, stage, 0.0)
                np.maximum(st['stage_max'], stage, out=st['stage_max'])
                st['area'][idx, kk] += np.where(wet, area, 0.0)
                st['dry_ever'] |= ~wet
            salt = np.where(dry, 0.0, salt)
        if rec:
            for cname, ck in CHECKPOINTS.items():
                at = kk == ck
                st['n_at'][cname] += at
                st['checks'][cname] += at & wet
    st['salt_left'] = salt
    st['storage_end'] = vol.copy()
    st['steps'] = (nm - SPINUP_MEANS) * STEPS_PER_MEAN
    st['record_years'] = (nm - SPINUP_MEANS) * R1.MEAN_DAYS / JULIAN_YEAR_DAYS
    # steps per tenth for each pool's cell (the same for every pool of a cell): the means after spin-up in each tenth
    st['steps_in_tenth'] = np.stack([(tenth[SPINUP_MEANS:, cell] == k).sum(axis=0) for k in range(10)], axis=1) * STEPS_PER_MEAN
    return st


def persistence(st):
    """Each pool's persistence class from its checks over the record: never wet; dry within the storm period; lasting
    the night (wet at sunrise in at least half the lunar days); lasting from one storm period to the next (wet at the
    end of the morning, before the storms, in at least half); standing (wet at every step)."""
    share = {c: st['checks'][c] / np.maximum(st['n_at'][c], 1) for c in CHECKPOINTS}
    cls = np.where(st['wet'] > 0, 1, 0)
    cls = np.where(share['sunrise'] >= 0.5, 2, cls)
    cls = np.where(share['before_storms'] >= 0.5, 3, cls)
    cls = np.where((st['wet'] > 0) & ~st['dry_ever'], 4, cls)
    return cls.astype(np.int8), share


CLASS_NAMES = ('never_wet', 'dry_within_the_storm_period', 'lasting_the_night', 'storm_period_to_the_next', 'standing')


# ---------------------------------------------------------------------------------------------------------------------
# The pools of one surface under one library: the grid's depressions and the bowls below the grid, both hollow ages.

def build_pools(inp, b, name):
    """The pools as one set: per hollow age, the dryland depressions (fed ones and those they spill into) and, per
    T21 cell and diameter bin, a wall-fed bowl (its own walls' runoff and what its rim catches from the slope) and a
    channel-fed bowl (with a share of the channels' capture), each standing for its number of bowls."""
    keep = inp['floor_t21'] >= 0
    nodes = inp['node'][keep]
    curves_d = depression_curves(inp['hyps_z'][keep].astype(float), inp['hyps_a'][keep].astype(float),
                                 inp['n_below'][keep], inp['pass_m'][keep].astype(float), inp['vmax_m3'][keep])
    par = inp['parent'][keep]
    pos = np.searchsorted(nodes, par)
    ok = (pos < nodes.size) & (nodes[np.minimum(pos, nodes.size - 1)] == par)
    parent_local = np.where(ok, pos, -1)
    parts, curves, offset = [], [], 0
    dry, reach = inp['t21_dry_km2'], inp['t21_reach_km']
    ncell = dry.size
    for si, s in enumerate(STATES):
        inflow = inp[f'inflow__{b}__{name}__{s}'].astype(float)
        parts.append(dict(kind=np.zeros(nodes.size, np.int8), state=np.full(nodes.size, si, np.int8),
                          cell=inp['floor_t21'][keep].astype(np.int64), inflow=inflow[keep],
                          parent_pos=np.where(parent_local >= 0, parent_local + offset, -1),
                          weight=np.ones(nodes.size), bin=np.full(nodes.size, -1, np.int16), node=nodes,
                          diameter_m=np.zeros(nodes.size), lost_outside=np.full(nodes.size, 0.0)))
        parts[-1]['lost_outside'][0] = float(inflow[~keep].sum())
        curves.append(curves_d)
        offset += nodes.size
    for si, s in enumerate(STATES):
        bc = RT.bowl_classes(s)
        nb = bc['number_km2'].size
        a_i = bc['number_km2'] * math.pi * bc['diameter_km'] ** 2 / 4.0
        r_i = bc['number_km2'] * bc['diameter_km']
        number = bc['number_km2'][None, :] * dry[:, None]                       # bowls per cell and bin
        crossed = np.minimum(number, r_i[None, :] * reach[:, None])
        inside = inp[f'bowls_inside__{b}__{name}__{s}'].astype(float)
        slope = inp[f'bowls_slope__{b}__{name}__{s}'].astype(float)
        chan = inp[f'bins__{b}__{name}__{s}'].astype(float)                    # (tenth, cell, bin)
        vol = dry[None, :, None] * 1000.0                                       # mm a year over the cell to m3 a year
        with np.errstate(divide='ignore', invalid='ignore'):
            own = (inside[:, :, None] * (a_i / a_i.sum())[None, None, :] + slope[:, :, None] * (r_i / r_i.sum())[None, None, :]) \
                * vol / np.maximum(number[None], 1e-30)
            ch = chan * vol / np.maximum(crossed[None], 1e-30)
        cells, bins = np.meshgrid(np.arange(ncell), np.arange(nb), indexing='ij')
        for kind, wgt, q in ((1, number - crossed, own), (2, crossed, own + ch)):
            use = (wgt > 0) & (q.sum(axis=0) > 0)
            c, k = cells[use], bins[use]
            parts.append(dict(kind=np.full(c.size, kind, np.int8), state=np.full(c.size, si, np.int8), cell=c.astype(np.int64),
                              inflow=q[:, use].T, parent_pos=np.full(c.size, -1), weight=wgt[use], bin=k.astype(np.int16),
                              node=np.zeros(c.size, np.int64), diameter_m=bc['diameter_km'][k] * 1000.0,
                              lost_outside=np.zeros(c.size)))
            curves.append(bowl_curves(bc['diameter_km'][k] * 1000.0))
            offset += c.size
    pools = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    width = max(c[0].shape[1] for c in curves)
    pools['s'], pools['A'], pools['V'], pools['top'] = pad_curves(curves, width)
    return pools


def case_summary(pools, st, inp, x):
    """Per hollow age: the depressions' and bowls' persistence, depths and fates, and per T21 cell the pools'
    evaporation, seepage and wet area by tenth and by persistence class."""
    cls, share = persistence(st)
    w = pools['weight']
    yrs = st['record_years']
    ncell = x['cls'].size
    dry_m2 = inp['t21_dry_km2'] * 1e6
    steps_k = np.maximum(st['steps_in_tenth'], 1)
    mean_area = st['area'] / steps_k                                              # (pool, tenth) m2 while wet
    out, cells = {}, {}
    for si, s in enumerate(STATES):
        rows = {}
        for kind, label in ((0, 'depressions'), (1, 'bowls_wall_fed'), (2, 'bowls_channel_fed')):
            m = (pools['state'] == si) & (pools['kind'] == kind)
            fed = m & ((pools['inflow'].sum(axis=1) > 0) | (st['wet'] > 0))
            ww = w[fed]
            wet_share = st['wet'][fed] / max(st['steps'], 1)
            with np.errstate(invalid='ignore', divide='ignore'):
                mean_depth = np.where(st['wet'][fed] > 0, st['stage_sum'][fed] / np.maximum(st['wet'][fed], 1), 0.0)
            counts = {CLASS_NAMES[c]: float(ww[cls[fed] == c].sum()) for c in range(len(CLASS_NAMES))}
            infl = float((ww * st['inflow'][fed]).sum() / yrs)
            row = dict(number=float(ww.sum()), classes=counts,
                       share_in_class={k: v / max(ww.sum(), 1e-30) for k, v in counts.items()},
                       inflow_km3_yr=infl * 1e-9, rain_on_pools_km3_yr=float((ww * st['rain'][fed]).sum() / yrs) * 1e-9,
                       evaporated_km3_yr=float((ww[:, None] * st['evap'][fed]).sum() / yrs) * 1e-9,
                       seeped_km3_yr=float((ww * st['seep'][fed]).sum() / yrs) * 1e-9,
                       spilled_out_km3_yr=float((ww * st['spill_out'][fed]).sum() / yrs) * 1e-9,
                       time_wet_share_mean=float((ww * wet_share).sum() / max(ww.sum(), 1e-30)),
                       time_deep_share_mean=float((ww * st['deep'][fed] / max(st['steps'], 1)).sum() / max(ww.sum(), 1e-30)),
                       carried_to_next_storms_share=float((ww * share['before_storms'][fed]).sum() / max(ww.sum(), 1e-30)),
                       mean_depth_m_median=float(np.median(mean_depth[mean_depth > 0])) if (mean_depth > 0).any() else 0.0,
                       largest_depth_m_quantiles=[float(v) for v in np.quantile(st['stage_max'][fed], [0.1, 0.5, 0.9, 0.99])]
                       if fed.any() else [0.0] * 4,
                       wet_area_km2_by_tenth=[float(v) for v in (ww[:, None] * mean_area[fed]).sum(axis=0) * 1e-6])
            if kind == 0:
                row['inflow_to_floors_outside_drylands_km3_yr'] = float(pools['lost_outside'][m].sum()) * 1e-9
            rows[label] = row
        out[s] = rows
        m = pools['state'] == si
        c = pools['cell'][m]
        ev = np.stack([np.bincount(c, weights=w[m] * st['evap'][m, k], minlength=ncell) for k in range(10)]) / yrs
        wet = np.stack([np.bincount(c, weights=w[m] * mean_area[m, k], minlength=ncell) for k in range(10)])
        by_cls = np.stack([np.bincount(c, weights=w[m] * mean_area[m].mean(axis=1) * (cls[m] == q), minlength=ncell)
                           for q in range(1, len(CLASS_NAMES))])
        cells[s] = dict(evap_mm_yr=ev / dry_m2 * 1000.0, seep_mm_yr=np.bincount(c, weights=w[m] * st['seep'][m], minlength=ncell)
                        / yrs / dry_m2 * 1000.0, rain_on_pools_mm_yr=np.bincount(c, weights=w[m] * st['rain'][m], minlength=ncell)
                        / yrs / dry_m2 * 1000.0, spill_out_mm_yr=np.bincount(c, weights=w[m] * st['spill_out'][m], minlength=ncell)
                        / yrs / dry_m2 * 1000.0, wet_fraction=wet / dry_m2, wet_fraction_by_class=by_cls / dry_m2,
                        salt_retained_g_yr=np.bincount(c, weights=w[m] * (st['salt_in'][m] - st['salt_out'][m]), minlength=ncell) / yrs)
    return out, cells, cls, share, mean_area


def salts(pools, st, cls, mean_area):
    """The solutes each pool keeps (left on its floor or held in its water), per unit of its mean wet area, and the
    years until its sodium and potassium carbonate make the top 0.3 m of its floor a saline soil."""
    yrs = st['record_years']
    kept = (st['salt_in'] - st['salt_out']) / yrs                                   # g a year
    area = np.maximum(mean_area.mean(axis=1), 1e-9)
    soda = SALTS['soda_share'] * kept / area                                        # g/m2 a year
    years = np.where(soda > 0, SALTS['saline_kg_m2'] * 1000.0 / np.maximum(soda, 1e-30), np.inf)
    hardpan = (1.0 - SALTS['soda_share']) * kept / area / (2700.0 * 0.7) * 1e6     # mm a millennium of 30%-porous cement
    return years, hardpan


def salt_summary(pools, st, cls, mean_area):
    years, hardpan = salts(pools, st, cls, mean_area)
    w = pools['weight']
    out = {}
    for si, s in enumerate(STATES):
        rows = {}
        for kind, label in ((0, 'depressions'), (1, 'bowls_wall_fed'), (2, 'bowls_channel_fed')):
            m = (pools['state'] == si) & (pools['kind'] == kind) & (st['wet'] > 0)
            ww = w[m]
            tot = max(ww.sum(), 1e-30)
            rows[label] = dict(saline_within_years={str(t): float(ww[years[m] <= t].sum() / tot) for t in (100, 300, 1000, 10000)},
                               saline_within_years_number={str(t): float(ww[years[m] <= t].sum()) for t in (100, 300, 1000)},
                               hardpan_mm_per_millennium_median=float(np.median(hardpan[m])) if m.any() else 0.0,
                               years_to_saline_median=(float(np.median(years[m])) if m.any() and np.isfinite(np.median(years[m]))
                                                       else None))
        out[s] = rows
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The leverage on PlaSim's land scheme.

def leverage_2b(x, inp, b, name, cells_case, dry_km2):
    """The drylands' evaporation with the pools against the bucket's, by year and by tenth, per class and hollow age:
    the water leaving the slopes' soil (phase 2a's fates, with the rain falling on pools) comes off the bucket's
    evaporation, either in the tenth it runs off (an upper bound of a tenth's change) or spread as the bucket's
    evaporation is (its timing), and the pools give back what they evaporate in the tenth they do."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    evap = R1.by_tenth(-x['evap'].astype(float) * R1.MM_PER_MEAN, x['tenth']) / years * yr     # (tenth, cell) mm a year
    rain = R1.by_tenth(x['pr'].astype(float) * R1.MM_PER_MEAN, x['tenth']) / years * yr
    rows = {}
    for s in STATES:
        f = {q: inp[f'{q}__{b}__{name}__{s}'].astype(float) for q in ('bowls_inside', 'bowls_slope', 'bowls_channel',
             'depressions', 'slope_recharge', 'bed_recharge', 'sea', 'humid_land')}
        hollows = f['bowls_inside'] + f['bowls_slope'] + f['bowls_channel'] + f['depressions']
        cc = cells_case[s]
        scale = dry_km2 / x['area_km2']                                           # pools' mm over the grid's dryland area
        pool_rain = cc['rain_on_pools_mm_yr'] * scale
        with np.errstate(invalid='ignore', divide='ignore'):
            rain_k = hollows / np.maximum(hollows.sum(axis=0), 1e-30) * pool_rain[None, :]
        leaves = hollows + f['slope_recharge'] + f['bed_recharge'] + f['sea'] + f['humid_land'] + rain_k
        pools = cc['evap_mm_yr'] * scale[None, :]
        out = {}
        for cname in R1.CLASSES + ('drylands',):
            m = np.ones(x['cls'].size, bool) if cname == 'drylands' else (x['cls'] == cname)
            a = x['area_km2'][m]
            e_k, l_k, p_k, r_k = ((v[:, m] * a).sum(axis=1) for v in (evap, leaves, pools, rain))
            e_y = e_k.sum()
            storm = p_k - l_k
            timed = p_k - l_k.sum() * e_k / e_y
            out[cname] = dict(change_share_of_evaporation_yearly=float((p_k.sum() - l_k.sum()) / e_y),
                              pools_evaporation_share_of_evaporation=float(p_k.sum() / e_y),
                              leaving_soil_share_of_evaporation=float(l_k.sum() / e_y),
                              largest_tenth_change_storm_timing=float(np.max(np.abs(storm) / np.maximum(e_k, 1e-30))),
                              largest_tenth_change_bucket_timing=float(np.max(np.abs(timed) / np.maximum(e_k, 1e-30))),
                              change_by_tenth_storm_timing=[float(v) for v in storm / np.maximum(e_k, 1e-30)],
                              change_by_tenth_bucket_timing=[float(v) for v in timed / np.maximum(e_k, 1e-30)],
                              pools_evaporation_by_tenth_share=[float(v) for v in p_k / np.maximum(e_k, 1e-30)],
                              recharge_share_of_rain=float(((f['slope_recharge'] + f['bed_recharge'])[:, m] * a).sum() / r_k.sum()
                                                           + (cc['seep_mm_yr'][m] * scale[m] * a).sum() / r_k.sum()),
                              passes_yearly=bool(abs(p_k.sum() - l_k.sum()) / e_y > R1.THRESHOLDS['yearly']),
                              passes_tenth_bucket_timing=bool(np.max(np.abs(timed) / np.maximum(e_k, 1e-30)) > R1.THRESHOLDS['tenth']),
                              passes_tenth_storm_timing=bool(np.max(np.abs(storm) / np.maximum(e_k, 1e-30)) > R1.THRESHOLDS['tenth']))
        rows[s] = out
    return rows


# ---------------------------------------------------------------------------------------------------------------------
# What sustains or wears each surface.

def wepp(texture, crust=False):
    vfs, clay, om = texture['vfs'], texture['clay'], texture['om']
    ki = 2.728e6 + 1.921e5 * vfs
    kr = 0.00197 + 0.00030 * vfs + 0.03863 * math.exp(-1.84 * om)
    tc = 2.67 + 0.065 * clay - 0.058 * vfs
    if crust:
        ki, kr, tc = ki * BIOCRUST['erodibility_factor'], kr * BIOCRUST['erodibility_factor'], tc * BIOCRUST['critical_shear_factor']
    return ki, kr, tc


def event_wear(c, texture, crust, gravity, kappa):
    """Soil detached by one runoff event on the design hillslope (kg/m2): interrill detachment Ki I sigma with the
    drops' energy factor kappa, and rill detachment Kr (tau - tau_c) in rills 1 m apart whose flow is Manning's with
    gravity explicit (deeper and slower on the Moon), tau = rho g h S."""
    ki, kr, tc = wepp(texture, crust)
    t = max(c['duration_h'], 1e-3) * 3600.0
    sigma = c['runoff_mm'] / 1000.0 / t
    rate = sigma + c['capacity_mm_h'] / 1000.0 / 3600.0
    length, slope = RT.HILLSLOPE['length_m'], RT.HILLSLOPE['slope']
    xx = (np.arange(200) + 0.5) / 200.0 * length
    q = sigma * xx * RILL['spacing_m']
    width = RILL['width_a'] * q ** RILL['width_b']
    n_eff = RT.HILLSLOPE['manning_n_earth'] / math.sqrt(gravity)
    h = (q * n_eff / (width * math.sqrt(slope))) ** 0.6
    tau = RHO_W * gravity * STANDARD_GRAVITY * h * slope
    rill_share = np.minimum(width / RILL['spacing_m'], 1.0)
    rill = float((kr * np.maximum(tau - tc, 0.0) * rill_share).mean()) * t
    inter = ki * kappa * rate * sigma * t * float((1.0 - rill_share).mean())
    return inter, rill


def slope_wear(x, off, name, lib, pressure):
    """Per dryland cell, the yearly soil loss (mm) of one surface under one library, with the same events at Earth's
    gravity and drops for comparison."""
    p = IN.parameters(name, pressure_pa=pressure)
    soil = IN.SURFACES[name]['soil']
    se = P2.ref_saturation(x, soil)
    classes = RT.event_classes(p, lib, se)
    kappa = float(IN.drop_energy_factor(pressure))
    texture = TEXTURES['loam' if soil == 'loam' else 'sandy_loam']
    crust = name == 'biocrust'
    years = off.shape[0] * R1.MEAN_DAYS / JULIAN_YEAR_DAYS
    runoff = off.sum(axis=0) / years                                                # mm a year
    tot = sum(c['probability'] * c['runoff_mm'] for c in classes)
    moon, earth, parts = np.zeros_like(runoff), np.zeros_like(runoff), []
    for c in classes:
        events = runoff * (c['probability'] * c['runoff_mm'] / tot) / c['runoff_mm']
        im, rm = event_wear(c, texture, crust, IN.MOON, kappa)
        ie, re = event_wear(c, texture, crust, 1.0, 1.0)
        moon += events * (im + rm)
        earth += events * (ie + re)
        parts.append(dict(runoff_mm=c['runoff_mm'], duration_h=c['duration_h'], interrill_kg_m2=im, rill_kg_m2=rm,
                          earth_interrill_kg_m2=ie, earth_rill_kg_m2=re))
    to_mm = 1000.0 / BULK_DENSITY
    return moon * to_mm, earth * to_mm, parts


def wear_summary(x, loss_mm):
    out = {}
    for cname in R1.CLASSES + ('drylands',):
        m = np.ones(x['cls'].size, bool) if cname == 'drylands' else (x['cls'] == cname)
        a = x['area_km2'][m]
        v = loss_mm[m]
        out[cname] = dict(loss_mm_yr=float((v * a).sum() / a.sum()), loss_mm_yr_max=float(v.max()),
                          share_of_area_above_formation={str(f): float(a[v > f].sum() / a.sum()) for f in SOIL_FORMATION_MM_YR})
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Checks.

def mead_check():
    """The bulk transfer at Earth's gravity against Lake Mead's measured latent heat (rule set before running)."""
    ce, z0 = transfer(np.array([MEAD['wind_m_s']]), 1.0, z_wind=MEAD['height_m'], z_air=MEAD['height_m'])
    p = MEAD['pressure_pa']
    ta = np.array(MEAD['ta']) + 273.15
    tw = np.array(MEAD['to']) + 273.15
    ea = np.array(MEAD['ea']) * 1000.0
    qa = 0.622 * ea / (p - 0.378 * ea)
    rho = p / (R_DRY * ta * (1.0 + 0.61 * qa))
    lam = 2.501e6 - 2370.0 * (tw - 273.15)
    le = rho * lam * ce[0] * MEAD['wind_m_s'] * (qsat(tw, p) - qa)
    ratio = float(le.mean() / MEAD['le_period'])
    monthly = le / np.array(MEAD['le'], float)
    return dict(rule=MEAD_RULE, transfer_coefficient=float(ce[0]), roughness_m=float(z0[0]), modelled_w_m2=float(le.mean()),
                measured_w_m2=MEAD['le_period'], ratio=ratio, monthly_ratio_range=[float(monthly.min()), float(monthly.max())],
                monthly_correlation=float(np.corrcoef(le, MEAD['le'])[0, 1]),
                modelled_mm_yr=float(le.mean() / 2.45e6 * 86400 * JULIAN_YEAR_DAYS), measured_mm_yr=[1881.0, 2074.0],
                passes=bool(abs(ratio - 1.0) <= MEAD_RULE['tolerance']))


def night_check(x, f, e_tab, pen):
    """The pools' day and night against the plan's hand estimates (rules set before running)."""
    semi = x['cls'] == 'semi-arid'
    a = x['area_km2'][semi]
    j = int(np.argmin(np.abs(D_REF - NIGHT_RULE['depth_m'])))
    e = e_tab[semi, j, :] * MEAN_H                                                   # mm in each tenth, a lunar day
    pt = pt_water(f)[:, semi].T * MEAN_H
    pn = pen[:, semi].T * MEAN_H
    wm = lambda v: float((v * a).sum() / a.sum())
    day, night = wm((e * DAY_FRACTION).sum(axis=1)), wm((e * (1.0 - DAY_FRACTION)).sum(axis=1))
    pt_day, pen_day = wm((pt * DAY_FRACTION).sum(axis=1)), wm((pn * DAY_FRACTION).sum(axis=1))
    shallow = wm((e_tab[semi, 0, :] * MEAN_H * DAY_FRACTION).sum(axis=1))
    by_depth = {f'{d:g}': dict(day_mm=wm((e_tab[semi, i, :] * MEAN_H * DAY_FRACTION).sum(axis=1)),
                               night_mm=wm((e_tab[semi, i, :] * MEAN_H * (1.0 - DAY_FRACTION)).sum(axis=1)))
                for i, d in enumerate(D_REF)}
    return dict(rule=NIGHT_RULE, day_mm=day, night_mm=night, lunar_day_mm=day + night, pt_water_day_mm=pt_day,
                penman_steady_day_mm=pen_day, shallowest_day_mm=shallow, day_over_pt=day / pt_day, by_depth=by_depth,
                hand_estimates=dict(register_lunar_day_mm=60.0, plan_day_mm=75.0, plan_night_mm=[5.0, 15.0],
                                    plan_lunar_day_mm=[80.0, 165.0]),
                passes=bool(NIGHT_RULE['day_over_pt'][0] <= day / pt_day <= NIGHT_RULE['day_over_pt'][1]
                            and NIGHT_RULE['night_mm'][0] <= night <= NIGHT_RULE['night_mm'][1]))


# ---------------------------------------------------------------------------------------------------------------------
# The biome maps' wet-patch layer.

WET_PATCHES = R1.HERE / 'results' / 'wet_patches.json'
WET_CELLS = P2.RUNS / 'wet_patches_2b.npz'
WET_SCHEMA = 'terluna.research.storm-water-wet-patches/1'
SENSITIVITY_TRANSFER = 1.45          # the factor that would close Lake Mead's latent heat: a calibration, kept out of
                                     # the check and used only as a sensitivity case


def rel(path):
    try:
        return str(path.relative_to(R1.ROOT))
    except ValueError:
        return str(path)


def sig(v, n=3):
    v = np.asarray(v, float)
    out = np.zeros_like(v)
    nz = v != 0
    mag = np.floor(np.log10(np.abs(v[nz])))
    out[nz] = np.round(v[nz] / 10 ** (mag - n + 1)) * 10 ** (mag - n + 1)
    return [float(a) for a in out.ravel()]


def wet_cells(pools, mean_area, cls, inp, si):
    """The grid depressions' pools on the 1.9 km cells: each pool's mean wet area by tenth spread over its lowest
    cells in height order (the floor cell first), with the pool's persistence class."""
    keep = inp['floor_t21'] >= 0
    m = np.flatnonzero((pools['kind'] == 0) & (pools['state'] == si))
    cell_ids = inp['hyps_cell'][keep]
    a = inp['hyps_a'][keep].astype(np.float32)
    before = np.cumsum(a, axis=1) - a
    area = mean_area[m].astype(np.float32)                                         # (dep, tenth)
    cover = np.clip(area[:, None, :] - before[:, :, None], 0.0, a[:, :, None])
    with np.errstate(invalid='ignore', divide='ignore'):
        frac = np.where(a[:, :, None] > 0, cover / np.maximum(a[:, :, None], 1e-9), 0.0)
    on = (frac.max(axis=2) > 0) & (cell_ids >= 0)
    rows = np.repeat(np.arange(m.size)[:, None], a.shape[1], axis=1)[on]
    return (cell_ids[on].astype(np.int64), np.round(frac[on] * 250.0).astype(np.uint8), cls[m][rows].astype(np.int8))


# ---------------------------------------------------------------------------------------------------------------------
# The whole phase.

def run_phase2b(x, libs, pressure):
    """Everything phase 2b adds: the product's section, the per-cell additions and the wet-patch layer."""
    with np.load(P2.PHASE2B_INPUTS) as z:
        meta_in = json.loads(str(z['metadata']))
        inp = {k: z[k] for k in z.files if k != 'metadata'}
    f = forcing(x)
    e_tab, t_tab = energy_balance(f, IN.MOON)
    pen, pen_t = penman(f, IN.MOON)
    e_hi, _ = energy_balance(f, IN.MOON, factor=SENSITIVITY_TRANSFER)
    pen_hi, _ = penman(f, IN.MOON, factor=SENSITIVITY_TRANSFER)
    years = x['pr'].shape[0] * R1.MEAN_DAYS / JULIAN_YEAR_DAYS
    rain_m_yr = x['pr'].astype(float).sum(axis=0) * R1.MM_PER_MEAN / 1000.0 / years
    solute = SALTS['rain_mg_l'] + SALTS['dust_soluble_g_m2_yr'] / np.maximum(rain_m_yr, 1e-3)    # g/m3 in the runoff
    ncell = x['cls'].size
    dry_km2 = inp['t21_dry_km2']
    pools_out, lev_out, salt_out, wear_out, layer, layer_cells, cells_rows = {}, {}, {}, {}, {}, {}, {}
    for b in R1.BRACKETS:
        for name in IN.MAIN:
            pools = build_pools(inp, b, name)
            off = inp[f'off__{b}__{name}'].astype(float)
            for seep in SEEPAGE:
                st = simulate(pools, x, e_tab, seep, IN.MOON, off, solute)
                summ, cc, cls, share, mean_area = case_summary(pools, st, inp, x)
                w = pools['weight']
                for si, s in enumerate(STATES):
                    m = pools['state'] == si
                    summ[s]['storage_change_km3_yr'] = float((w[m] * (st['storage_end'][m] - st['storage_start'][m])).sum()
                                                             / st['record_years']) * 1e-9
                pools_out.setdefault(b, {}).setdefault(name, {})[seep] = summ
                lev_out.setdefault(b, {}).setdefault(name, {})[seep] = leverage_2b(x, inp, b, name, cc, dry_km2)
                salt_out.setdefault(b, {}).setdefault(name, {})[seep] = salt_summary(pools, st, cls, mean_area)
                for si, s in enumerate(STATES):
                    key = f'{b}__{name}__{s}__{seep}'
                    layer[key] = dict(wet_fraction_by_tenth=cc[s]['wet_fraction'], wet_fraction_by_class=cc[s]['wet_fraction_by_class'])
                    layer_cells[key] = wet_cells(pools, mean_area, cls, inp, si)
                    if b == 'walnut_gulch' and seep == 'clay_floored':
                        cells_rows.setdefault(name, {})[s] = cc[s]
                d = lev_out[b][name][seep]['young']['drylands']
                print(f"pools {b} {name} {seep}: change {d['change_share_of_evaporation_yearly']:+.3f} a year, "
                      f"largest tenth {d['largest_tenth_change_bucket_timing']:.3f} (young)", file=sys.stderr)
            loss, earth, parts = slope_wear(x, off, name, libs[b], pressure)
            wear_out.setdefault(b, {})[name] = dict(classes=wear_summary(x, loss), earth_gravity_and_drops=wear_summary(x, earth),
                                                    events=parts, cell_mm_yr=loss)
    # Sensitivity: the transfer raised by Lake Mead's factor for the loam under Walnut Gulch's storms.
    pools = build_pools(inp, 'walnut_gulch', 'loam')
    st = simulate(pools, x, e_hi, 'clay_floored', IN.MOON, inp['off__walnut_gulch__loam'].astype(float), solute)
    summ, cc, cls, share, mean_area = case_summary(pools, st, inp, x)
    sens = dict(transfer_factor=SENSITIVITY_TRANSFER, case='walnut_gulch, loam, clay_floored',
                pools=summ, leverage=leverage_2b(x, inp, 'walnut_gulch', 'loam', cc, dry_km2),
                night_check=night_check(x, f, e_hi, pen_hi))
    tight = dict(SEAL, k_earth_mm_h=SEAL['k_earth_mm_h'] / 10.0)
    keep_seal = dict(SEAL)
    SEAL.update(tight)
    try:
        st = simulate(pools, x, e_tab, 'clay_floored', IN.MOON, inp['off__walnut_gulch__loam'].astype(float), solute)
    finally:
        SEAL.update(keep_seal)
    summ, cc, cls, share, mean_area = case_summary(pools, st, inp, x)
    sens_seal = dict(seal=tight, case='walnut_gulch, loam, clay_floored', pools=summ,
                     leverage=leverage_2b(x, inp, 'walnut_gulch', 'loam', cc, dry_km2))
    semi = x['cls'] == 'semi-arid'
    a = x['area_km2'][semi]
    energy = dict(depths_m=[float(v) for v in D_REF],
                  semi_arid_evaporation_mm_h_by_tenth={f'{d:g}': sig((e_tab[semi, i, :] * a[:, None]).sum(axis=0) / a.sum())
                                                       for i, d in enumerate(D_REF)},
                  semi_arid_water_minus_air_k_by_tenth={f'{d:g}': sig(((t_tab[semi, i, :] - f['ta'][:, semi].T) * a[:, None]).sum(axis=0) / a.sum())
                                                        for i, d in enumerate(D_REF)},
                  semi_arid_penman_mm_h_by_tenth=sig((pen[:, semi] * a).sum(axis=1) / a.sum()),
                  semi_arid_pt_water_mm_h_by_tenth=sig((pt_water(f)[:, semi] * a).sum(axis=1) / a.sum()),
                  semi_arid_wind_m_s_by_tenth=sig((f['wind'][:, semi] * a).sum(axis=1) / a.sum()),
                  transfer_coefficient_mean=float(transfer(f['wind'], IN.MOON)[0].mean()),
                  seal_seepage_mm_day=dict(depth_m=[0.1, 0.5, 1.0, 3.0],
                                           moon=sig(seal_seepage(np.array([0.1, 0.5, 1.0, 3.0]), IN.MOON) * 8.64e7),
                                           earth=sig(seal_seepage(np.array([0.1, 0.5, 1.0, 3.0]), 1.0) * 8.64e7)),
                  bank_seepage_mm_day=dict(moon=bank_seepage(IN.MOON) * 8.64e7, earth=bank_seepage(1.0) * 8.64e7))
    checks = dict(lake_mead=mead_check(), night=night_check(x, f, e_tab, pen))
    # The wet-patch layer: the per-T21 product and the 1.9 km companion.
    keys = sorted(layer)
    arrays = {}
    for k in keys:
        c, frac, cl = layer_cells[k]
        arrays[f'{k}__cell'], arrays[f'{k}__wet_250ths'], arrays[f'{k}__class'] = c, frac, cl
    lmeta = dict(schema=WET_SCHEMA.replace('/1', '-cells/1'), variants=keys, classes=list(CLASS_NAMES),
                 rule='per variant (library__surface__age__seepage): the 1.9 km cells (row * 5760 + column) a grid '
                      'depression\'s pool covers, the covered share by tenth in 250ths, and the pool\'s persistence '
                      'class; the bowls below the grid are in the per-T21 product and cover every 1.9 km cell of their T21 '
                      'cell alike', drainage_sources=meta_in['drainage_sources'])
    np.savez_compressed(WET_CELLS, metadata=json.dumps(lmeta), **arrays)
    product_layer = dict(
        schema=WET_SCHEMA,
        evidence=('A model result of the offline storm-water chain on the design years (A28R_301k_design_l15_snow_s934, '
                  'years 10-19), with stand-in storms (Walnut Gulch\'s measured family and CM1\'s 6 km window over wet '
                  'ground), designed surfaces and designed pool floors; no climate feedback. Its numbers become final '
                  'when the GCM and CM1 agree and the drylands\' own storms replace the stand-ins.'),
        reading_rule=('wet_ppm_by_tenth (millionths, per cell in the order of cells): the mean share of each T21 cell\'s dryland area under standing water (1 cm '
                      'or more at a pool\'s deepest point) in each tenth of the lunar day, from the tenth centred 162 '
                      'degrees before noon to the one 162 after; wet_ppm_by_class: the day-mean wet share in pools of '
                      'each persistence class (' + ', '.join(CLASS_NAMES[1:]) + '); grid depressions and the crater bowls '
                      'below the grid together. Variants are library__surface__age__seepage.'),
        variants=keys, classes=list(CLASS_NAMES[1:]),
        cells=[dict(lat=R1.r3(x['lat'][c]), lon_e=R1.r3(x['lon'][c]), unep=str(x['cls'][c]), dryland_km2=R1.r3(dry_km2[c]))
               for c in range(ncell)],
        layers={k: dict(wet_ppm_by_tenth=np.rint(layer[k]['wet_fraction_by_tenth'].T * 1e6).astype(int).tolist(),
                        wet_ppm_by_class=np.rint(layer[k]['wet_fraction_by_class'].T * 1e6).astype(int).tolist()) for k in keys},
        cells_1_9km=dict(path=rel(WET_CELLS), schema=lmeta['schema'], sha256=R1.sha256(WET_CELLS)))
    section = dict(
        design=dict(water_albedo=WATER_ALBEDO, water_emissivity=WATER_EMISSIVITY, floor_heat=FLOOR_HEAT, charnock=CHARNOCK,
                    scalar_roughness_m=Z0_SCALAR_M, depths_m=[float(v) for v in D_REF], steps_per_mean=STEPS_PER_MEAN,
                    wet_m=WET_M, deep_m=DEEP_M, checkpoints=CHECKPOINTS, spinup_means=SPINUP_MEANS, seal=SEAL, bank=BANK,
                    salts=SALTS, textures=TEXTURES, biocrust=BIOCRUST, rill=RILL, bulk_density=BULK_DENSITY,
                    soil_formation_mm_yr=SOIL_FORMATION_MM_YR, bowl_shape='paraboloid, Pike\'s fresh form at half its depth',
                    classes=list(CLASS_NAMES)),
        energy=energy, checks=checks, pools=pools_out, leverage=lev_out, salts=salt_out,
        wear={b: {n: {k: v for k, v in d.items() if k != 'cell_mm_yr'} for n, d in v_.items()} for b, v_ in wear_out.items()},
        sensitivity_transfer=sens,
        sensitivity_tight_seal=sens_seal,
        grid_dryland_over_t21_area=float(dry_km2.sum() / x['area_km2'].sum()),
        wet_patches=dict(path=rel(WET_PATCHES), schema=WET_SCHEMA),
        inputs=dict(path=rel(P2.PHASE2B_INPUTS), schema=meta_in['schema']))
    cells = {}
    for c in range(ncell):
        row = {}
        for n in IN.MAIN:
            row[n] = {s: dict(wet_fraction=R1.r3(cells_rows[n][s]['wet_fraction'][:, c].mean()),
                              wet_fraction_by_tenth=R1.r3(cells_rows[n][s]['wet_fraction'][:, c]),
                              pool_evaporation_mm_yr=R1.r3(cells_rows[n][s]['evap_mm_yr'][:, c].sum()),
                              pool_seepage_mm_yr=R1.r3(cells_rows[n][s]['seep_mm_yr'][c])) for s in STATES}
            row[n]['slope_loss_mm_yr'] = R1.r3(wear_out['walnut_gulch'][n]['cell_mm_yr'][c])
        cells[c] = dict(pools_2b_walnut_gulch_clay_floored=row)
    return section, cells, product_layer
