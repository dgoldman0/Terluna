"""How far below the freezing level lunar snowflakes fall before they melt, and where their water and heat go.

    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting earth       # the Earth check, seconds
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting moon        # the lunar cases, minutes
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting assemble    # results/snow_melting.json

A flake leaves the freezing level with the mass, span and fall speed McSnow gives it (Mitchell's (1996) unrimed
aggregates, m = 0.0028 D^2.1 and A = 0.2285 D^1.88 in cgs units, falling by Boehm's law; boehm.py) and falls at its
terminal speed through a measured column of air.  It passes through three stages.

* Dry ice.  Its surface sits at the ice-bulb temperature, where the heat the air conducts to it balances the latent
  heat of the vapour it exchanges with the air (Pruppacher and Klett 2010, ch. 13); in air subsaturated over ice it
  sublimates, cooling itself and the air.  It cannot melt until that temperature reaches 0 C, which is the ice-bulb
  criterion of Heymsfield and co-authors (2021) and lies a little below the level where the wet-bulb temperature
  passes 0 C, by the ratio of the latent heats of sublimation and evaporation.
* Melting.  Its surface stays at 0 C.  The heat the air conducts to it, less the latent heat of the vapour it
  exchanges, melts ice (Mitra and co-authors 1990, as Frick, Seifert and Wernli 2013 give the model, their eqs. 1 and
  7-17): the flake is an oblate spheroid whose axis ratio rises from 0.3 to 1 with its liquid fraction l, its density
  is interpolated from the dry flake's to water's, its capacitance carries Mitra's factor 0.8 + 0.2 l, its fall speed
  moves from the dry flake's to the drop's by Mitra's measured function 0.246 l + 0.754 l^7, and its ventilation is
  Hall and Pruppacher's (1976) at a Reynolds number interpolated between the dry flake and the drop.  The vapour it
  exchanges comes from exposed dry ice, as sublimation, in proportion to the ice's share of its mass, and from its
  meltwater, as evaporation, in the rest; that split changes neither its masses nor the heat it takes from the air.
* Drop.  Once its ice has gone it is a raindrop, followed by the rain solver's single-drop physics
  (atmosphere/rain/drops.py: drag, shape, ventilated evaporation and the drop's own temperature) to the ground, to
  complete evaporation, or past the stability screen, where it breaks.

The air below the freezing level comes from the design-climate CM1 boxes (melting_environment.py): every column that
carries snow through its freezing level, weighted by that snow, each column's temperature, pressure and vapour taken
as they were saved.  Still air: the boxes' updraughts and downdraughts are left out.  The Earth check runs the same
code with Earth's gravity and air through a straight winter sounding.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from shared.constants import GAS_CONSTANT, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.radiative_convective.thermodynamics import DryAir, MOLAR_MASS, earthlike_air
from atmosphere.rain import drops
from atmosphere.snow import boehm, flakes

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = HERE / 'results/snow_melting.json'
RAW = ROOT / 'research/runs/snow_agent_3/melting'
SCHEMA = 'terluna.atmosphere.snow-melting/1'

T0 = 273.15
LF = flakes.LATENT_FUSION                     # J/kg at 0 C
RV = GAS_CONSTANT / MOLAR_MASS['H2O']
RHO_W = drops.WATER_DENSITY
LV0 = drops.latent_heat(T0)
LS0 = LV0 + LF
# McSnow's unrimed aggregates (mo_mass2diam.f90, after Mitchell 1996), SI
AGG_ALPHA, AGG_BETA = 0.0028 * 10 ** (2 * 2.1 - 3), 2.1
AGG_SIGMA, AGG_GAMMA = 0.2285 * 10 ** (2 * 1.88 - 4), 1.88
# Mitra and co-authors (1990) as Frick and co-authors (2013) give them
AXIS_DRY = 0.3
CAPACITANCE_DRY = 0.8
PSI_LINEAR = 0.246
SNOW_DENSITY_CAP = 500.
ICE_DONE = 1e-4                                # the flake is a drop once its ice is below this share of its mass
DEPTH_STEP = 50.                               # m: the grid the budgets are recorded on


# --- the air -------------------------------------------------------------------------------------------------

class Air:
    """Air in still columns given by temperature, pressure and vapour pressure on heights, tabulated with the rain
    solver's transport properties every 20 m and read linearly.  state() serves both the snow and the rain naming."""

    def __init__(self, z, T, p, e, *, gravity, air, updraft=0.):
        order = np.argsort(z)
        z, T, p, e = (np.asarray(a, float)[order] for a in (z, T, p, e))
        self.gravity_m_s2, self.updraft_m_s, self.air = float(gravity), float(updraft), air
        self.top = float(z[-1])
        self.grid = np.arange(max(z[0], 0.), z[-1] + 1e-6, 20.)
        Tg = np.interp(self.grid, z, T)
        pg = np.exp(np.interp(self.grid, z, np.log(p)))
        eg = np.interp(self.grid, z, e)
        props = np.array([drops.transport(a, b, c / b, air.fractions) for a, b, c in zip(Tg, pg, eg)])
        rho = (pg - eg) * air.molar_mass / (GAS_CONSTANT * Tg) + eg / (RV * Tg)
        self.table = dict(T=Tg, p=pg, e=eg, eta=props[:, 0], dv=props[:, 1], ka=props[:, 2], cp=props[:, 3], rho=rho)

    def state(self, z):
        z = float(min(max(z, self.grid[0]), self.grid[-1]))
        i = min(int((z - self.grid[0]) / 20.), len(self.grid) - 2)
        f = (z - self.grid[i]) / 20.
        s = {k: v[i] + (v[i + 1] - v[i]) * f for k, v in self.table.items()}
        return dict(t=s['T'], p=s['p'], e=s['e'], rho=s['rho'], eta=s['eta'], dv=s['dv'], ka=s['ka'],
                    altitude_m=z, temperature_k=s['T'], pressure_pa=s['p'], vapour_pressure_pa=s['e'],
                    relative_humidity=s['e'] / flakes.liquid_saturation_pressure(s['T']), density_kg_m3=s['rho'],
                    viscosity_pa_s=s['eta'], vapour_diffusivity_m2_s=s['dv'], thermal_conductivity_w_m_k=s['ka'],
                    heat_capacity_j_kg_k=s['cp'])

    def freezing_level(self):
        T = self.table['T']
        cross = np.flatnonzero((T[:-1] >= T0) & (T[1:] < T0))
        if T[0] < T0 or not cross.size:
            return None
        k = cross[0]
        return float(self.grid[k] + (T[k] - T0) / (T[k] - T[k + 1]) * 20.)


def straight_sounding(gravity, air, z_freeze, p_freeze, lapse, rh, top_above=500.):
    """A straight temperature profile through the freezing level, hydrostatic, at one relative humidity over water
    below it and over ice above it: the Earth check's winter sounding."""
    z = np.arange(0., z_freeze + top_above + 1., 25.)
    T = T0 - lapse * (z - z_freeze)
    expo = gravity * air.molar_mass / (GAS_CONSTANT * lapse)
    p = p_freeze * (T / T0) ** expo
    sat = np.where(T >= T0, [flakes.liquid_saturation_pressure(x) for x in T], [flakes.ice_saturation_pressure(x) for x in T])
    return Air(z, T, p, rh * sat, gravity=gravity, air=air)


MOON_AIR = earthlike_air(1.2 * 101325.0, 400.0)
EARTH_AIR = DryAir(101325.0)


# --- the flake -------------------------------------------------------------------------------------------------

def aggregate_span(m):
    return (m / AGG_ALPHA) ** (1 / AGG_BETA)


def aggregate_mass(d):
    return AGG_ALPHA * d ** AGG_BETA


def dry_speed(m, s, gravity):
    """Boehm's terminal speed of the McSnow aggregate of mass m (boehm.py), m/s, and its span."""
    d = aggregate_span(m)
    q = min(1., AGG_SIGMA * d ** AGG_GAMMA / (math.pi / 4 * d * d))
    v = boehm.particle_properties(d, m, gravity, dict(rho=s['rho'], eta=s['eta']), q)['v']
    return float(v), d


class DropSpeeds:
    """The rain solver's terminal speeds in one column, tabulated over diameter and height and read bilinearly in
    log diameter."""

    def __init__(self, column, heights, diameters=np.geomspace(1e-4, 0.025, 40)):
        self.lnd, self.h = np.log(diameters), np.asarray(heights, float)
        self.v = np.array([[drops.terminal_speed(d, column, z) for d in diameters] for z in self.h])

    def __call__(self, d, z):
        x = float(np.clip(math.log(d), self.lnd[0], self.lnd[-1]))
        zz = float(np.clip(z, self.h[0], self.h[-1]))
        j = min(int(np.searchsorted(self.lnd, x)), len(self.lnd) - 1)
        j = max(j, 1)
        fx = (x - self.lnd[j - 1]) / (self.lnd[j] - self.lnd[j - 1])
        i = max(min(int(np.searchsorted(self.h, zz)), len(self.h) - 1), 1)
        fz = (zz - self.h[i - 1]) / (self.h[i] - self.h[i - 1])
        a = self.v[i - 1, j - 1] * (1 - fx) + self.v[i - 1, j] * fx
        b = self.v[i, j - 1] * (1 - fx) + self.v[i, j] * fx
        return a * (1 - fz) + b * fz


def ventilation(x):
    """Hall and Pruppacher (1976) for ice, as Mitra's model and Frick and co-authors take it."""
    return 1. + 0.14 * x * x if x < 1. else 0.86 + 0.28 * x


def spheroid_factor(a):
    if a >= 0.9999:
        return 1.
    e = math.sqrt(1. - a * a)
    return e / math.asin(e)


def flake_state(m_i, m_w, z, column, speeds, band):
    """Liquid fraction, span, capacitance, fall speed and the ventilation factors (vapour, heat) of a flake."""
    s = column.state(z)
    m = m_i + m_w
    ell = min(max(m_w / m, 0.), 1.)
    v_s, d_s = dry_speed(m, s, column.gravity_m_s2)
    v_s *= band['speed']
    d_eq = (6 * m / (math.pi * RHO_W)) ** (1 / 3)
    v_r = speeds(d_eq, z)
    a = AXIS_DRY + (1 - AXIS_DRY) * ell
    rho_s = min(6 * m / (math.pi * d_s ** 3 * a), SNOW_DENSITY_CAP)
    rho_m = rho_s + (RHO_W - rho_s) * ell
    d_m = (6 * m / (math.pi * a * rho_m)) ** (1 / 3)
    c = band['capacitance'] * (CAPACITANCE_DRY + (1 - CAPACITANCE_DRY) * ell) * d_m / 2 * spheroid_factor(a)
    psi = PSI_LINEAR * ell + (1 - PSI_LINEAR) * ell ** 7
    v = v_s + (v_r - v_s) * psi
    ratio = drops.axis_ratio(d_eq, column.gravity_m_s2, s['rho'], T0) if ell > 0 else 1.
    re = (1 - ell) * s['rho'] * v_s * d_s / s['eta'] + ell * s['rho'] * v_r * d_eq * ratio ** (-1 / 3) / s['eta']
    sc = s['eta'] / (s['rho'] * s['dv'])
    pr = s['eta'] * s['heat_capacity_j_kg_k'] / s['ka']
    return dict(s=s, ell=ell, d_s=d_s, d_m=d_m, d_eq=d_eq, c=c, v=v, v_s=v_s, v_r=v_r,
                fv=ventilation(sc ** (1 / 3) * math.sqrt(re)), fh=ventilation(pr ** (1 / 3) * math.sqrt(re)))


def ice_bulb(s, fv, fh):
    """Surface temperature of dry ice in steady exchange with the air: conduction balances the latent heat of
    sublimation (or deposition)."""
    rho_v = s['e'] / (RV * s['t'])

    def balance(ts):
        return s['ka'] * fh * (s['t'] - ts) - LS0 * s['dv'] * fv * (flakes.ice_saturation_pressure(ts) / (RV * ts) - rho_v)

    if balance(T0) >= 0:
        return T0
    return brentq(balance, s['t'] - 60., T0, xtol=1e-6)


def rates(m_i, m_w, z, column, speeds, band):
    """Rates (kg/s) of sublimation from ice, melting and evaporation from water, positive as losses to vapour and
    to water, and the fall speed; heat is what the air conducts to the flake (W) and vapour its whole exchange (kg/s).
    Negative sublimation or evaporation is deposition or condensation."""
    st = flake_state(m_i, m_w, z, column, speeds, band)
    s, c = st['s'], st['c']
    rho_v = s['e'] / (RV * s['t'])
    ts = ice_bulb(s, st['fv'], st['fh']) if m_w <= 0 else T0
    if ts < T0:                                                       # dry ice below 0 C at its surface
        sub = 4 * math.pi * c * st['fv'] * s['dv'] * (flakes.ice_saturation_pressure(ts) / (RV * ts) - rho_v)
        q = 4 * math.pi * c * st['fh'] * s['ka'] * (s['t'] - ts)
        return dict(sub=sub, melt=0., evap=0., v=st['v'], ts=ts, st=st, heat=q, vapour=sub)
    q = 4 * math.pi * c * st['fh'] * s['ka'] * (s['t'] - T0)
    e = 4 * math.pi * c * st['fv'] * s['dv'] * (flakes.liquid_saturation_pressure(T0) / (RV * T0) - rho_v)
    ice_loss = (q - LV0 * e) / LF
    if ice_loss < 0 and m_w <= 0:                                     # no heat to melt: stays dry at 0 C
        ice_loss = 0.
    share = m_i / (m_i + m_w)                                         # exposed ice takes this share of the vapour
    sub = share * e
    return dict(sub=sub, melt=ice_loss - sub, evap=(1 - share) * e, v=st['v'], ts=T0, st=st, heat=q, vapour=e)


def follow(d0, column, z_start, speeds, *, band=None, m0=None, z_ground=0., record_step=DEPTH_STEP):
    """One flake from z_start: budgets by depth below z_start every record_step, the depths of the stages and the
    drop's fate.  Masses are fractions of the flake's starting mass."""
    band = band or dict(capacitance=1., speed=1.)
    m0 = aggregate_mass(d0) if m0 is None else m0
    depth_end = z_start - z_ground
    grid = np.arange(0., depth_end + 1e-6, record_step)
    out = dict(span_m=float(aggregate_span(m0)), mass_kg=float(m0),
               melted_diameter_m=float((6 * m0 / (math.pi * RHO_W)) ** (1 / 3)))

    # flake stages, integrated in depth: y = (ice, water, sublimated, melted, evaporated, time) as fractions / s
    def f(depth, y):
        mi, mw = max(y[0], 1e-12) * m0, max(y[1], 0.) * m0
        r = rates(mi, mw, z_start - depth, column, speeds, band)
        v = max(r['v'], 1e-3)
        return [-(r['sub'] + r['melt']) / m0 / v, (r['melt'] - r['evap']) / m0 / v,
                r['sub'] / m0 / v, r['melt'] / m0 / v, r['evap'] / m0 / v, 1. / v]

    def melted(depth, y):
        return y[0] - ICE_DONE * (y[0] + y[1])

    def gone(depth, y):
        return y[0] + y[1] - 1e-6

    def onset(depth, y):
        return y[1] - 1e-6

    def half(depth, y):
        return y[1] - y[0]

    melted.terminal, melted.direction = True, -1
    gone.terminal, gone.direction = True, -1
    onset.direction, half.direction = 1, 1
    sol = solve_ivp(f, (0., depth_end), [1., 0., 0., 0., 0., 0.], method='LSODA', rtol=1e-5,
                    atol=[1e-8, 1e-8, 1e-9, 1e-9, 1e-9, 1e-2], max_step=50., dense_output=True,
                    events=(melted, gone, onset, half))
    if not sol.success:
        raise RuntimeError(sol.message)
    end = float(sol.t[-1])
    first = lambda k: float(sol.t_events[k][0]) if len(sol.t_events[k]) else None
    out.update(onset_depth_m=first(2), half_melted_depth_m=first(3), melted_depth_m=first(0), flake_gone_depth_m=first(1))
    s_on = out['onset_depth_m']
    if s_on is not None:
        st = column.state(z_start - s_on)
        out['onset_temperature_c'] = float(st['t'] - T0)
    if out['melted_depth_m'] is not None:
        out['melted_temperature_c'] = float(column.state(z_start - end)['t'] - T0)
    flake_grid = grid[grid <= end]
    y = sol.sol(flake_grid) if len(flake_grid) else np.zeros((6, 0))
    rec = dict(ice=list(np.clip(y[0], 0, None)), water=list(np.clip(y[1], 0, None)), sub=list(y[2]), melt=list(y[3]),
               evap=list(y[4]), time=list(y[5]))
    y_end = sol.y[:, -1]
    fate = 'ground' if out['melted_depth_m'] is None and out['flake_gone_depth_m'] is None else None
    if out['flake_gone_depth_m'] is not None:
        fate = 'sublimated_or_evaporated_as_flake'
    out['ice_at_end'] = float(max(y_end[0], 0.))
    # the drop
    drop = None
    if out['melted_depth_m'] is not None and end < depth_end - 1.:
        drop = fall_drop(max(y_end[1], 0.) * m0, column, z_start - end, z_ground)
        fate = drop['status']
        rest = grid[grid > end]
        mass = np.interp(rest, drop['depth_below_melt_m'] + end, drop['mass_kg'], right=np.nan) / m0
        mass = np.where(np.isnan(mass), 0., mass)
        lost = max(y_end[1], 0.) - mass
        rec['ice'] += [0.] * len(rest)
        rec['water'] += list(mass)
        rec['sub'] += [float(y_end[2])] * len(rest)
        rec['melt'] += [float(y_end[3])] * len(rest)
        rec['evap'] += list(float(y_end[4]) + lost)
        rec['time'] += list(float(y_end[5]) + np.interp(rest, drop['depth_below_melt_m'] + end, drop['time_s'],
                                                         right=drop['time_s'][-1]))
        out['drop'] = {k: v for k, v in drop.items() if k not in ('depth_below_melt_m', 'mass_kg', 'time_s')}
    elif end < depth_end - 1.:                                       # the flake itself vanished aloft
        rest = grid[grid > end]
        for k in ('ice', 'water'):
            rec[k] += [0.] * len(rest)
        for k in ('sub', 'melt', 'evap', 'time'):
            rec[k] += [rec[k][-1] if rec[k] else 0.] * len(rest)
    out['fate'] = fate
    out['record'] = {k: np.asarray(v, float) for k, v in rec.items()}
    out['depth_grid_m'] = grid
    return out


def fall_drop(m, column, z_start, z_ground):
    """The rain solver's drop from z_start at its terminal speed and 0 C, to the ground or to evaporation; mass and
    time by depth below z_start.  A drop larger than its stability screen breaks: it is halved in diameter into eight
    equal drops, again until they hold together, and one of them is followed for all."""
    d = (6 * m / (math.pi * RHO_W)) ** (1 / 3)
    screen = drops.stable_diameter(column, z_start)
    pieces = 1
    while (6 * m / pieces / (math.pi * RHO_W)) ** (1 / 3) > screen:   # a drop above its screen breaks: halve it
        pieces *= 8
    whole, m = m, m / pieces
    d_piece = (6 * m / (math.pi * RHO_W)) ** (1 / 3)
    t_start = drops.TRIPLE_POINT_K + 1e-3                            # the solver's liquid starts at the triple point
    u0 = drops.terminal_speed(d_piece, column, z_start, t_start)

    def rhs(t, y):
        return drops.drop_rates(y, column)

    def ground(t, y):
        return y[0] - z_ground

    def evaporated(t, y):
        return y[2] - drops.CUTOFF_MASS

    def unstable(t, y):
        return drops.stable_diameter(column, max(y[0], 0.)) - (max(y[2], 0.) * 6 / (RHO_W * math.pi)) ** (1 / 3)

    for ev in (ground, evaporated, unstable):
        ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0., 2e5), [z_start, u0, m, t_start], method='Radau', rtol=1e-5,
                    atol=[1e-3, 1e-6, max(m * 1e-9, 1e-20), 1e-5], max_step=60., events=(ground, evaporated, unstable))
    if not sol.success:
        raise RuntimeError(sol.message)
    names = ('ground', 'evaporated', 'unstable')
    status = next((n for n, hits in zip(names, sol.t_events) if len(hits)), 'time_limit')
    depth = z_start - sol.y[0]
    keep = np.concatenate([[True], np.diff(depth) > 0])
    final = max(sol.y[2, -1], 0.) if status == 'ground' else 0.
    return dict(status=status, initial_diameter_m=float(d), screen_m=float(screen), pieces=pieces,
                depth_below_melt_m=depth[keep], mass_kg=np.clip(sol.y[2][keep], 0., None) * pieces, time_s=sol.t[keep],
                survived_fraction=float(final / m), ground_diameter_m=float((6 * final / (math.pi * RHO_W)) ** (1 / 3))
                if final > 0 else 0., ground_speed_m_s=float(sol.y[1, -1]) if status == 'ground' else None)


# --- the lunar columns -------------------------------------------------------------------------------------------

def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    c = np.cumsum(weights[order])
    return np.interp(np.asarray(q) * c[-1], c, values[order])


def column_air(data, i):
    """One saved CM1 column, ground to about 32 km."""
    return Air(data['zh_m'], data['T'][i], data['p'][i], data['e'][i], gravity=MOON_SURFACE_GRAVITY, air=MOON_AIR)


def composite_air(data, select=None, q=.5, step=100.):
    """The snow-weighted quantile q of temperature and of vapour pressure at each depth below the freezing level over
    the selected columns, hung from their weighted median freezing level; pressure takes the median.  Depths reach the
    deepest ground among the selected columns' freezing levels."""
    sel = np.ones(len(data['z0_m']), bool) if select is None else select
    z0, w = data['z0_m'][sel], data['snow_kg_m3'][sel]
    zh = data['zh_m']
    zf = float(weighted_quantile(z0, w, .5))
    depths = np.arange(0., zf + 1e-6, step)
    out = {k: [] for k in ('T', 'p', 'e')}
    for d in depths:
        z = z0 - d
        ok = z >= zh[0]
        if ok.sum() < 20:
            z = np.maximum(z, zh[0])
            ok = np.ones_like(ok)
        idx = np.clip(np.searchsorted(zh, z), 1, len(zh) - 1)
        f = (z - zh[idx - 1]) / (zh[idx] - zh[idx - 1])
        rows = np.flatnonzero(sel)
        for k in out:
            a = data[k][rows, idx - 1] * (1 - f) + data[k][rows, idx] * f
            out[k].append(float(weighted_quantile(a[ok], w[ok], q if k == 'e' else (.5 if k == 'p' else q))))
    z = zf - depths
    return Air(z, np.array(out['T']), np.array(out['p']), np.array(out['e']), gravity=MOON_SURFACE_GRAVITY,
               air=MOON_AIR), zf


TABLE_SPANS = np.geomspace(1e-3, 0.4, 28)       # m: the flake spans each column's response table holds
BANDS = {'central': dict(capacitance=1., speed=1.),
         # the crystal model's Earth error band (earth_error_band.json) at the edges that make melting deepest and
         # shallowest: less capacitance with more speed, and more capacitance with less speed
         'deep_edge': dict(capacitance=0.8357441864385768, speed=1.1590500041389329),
         'shallow_edge': dict(capacitance=1.1089363415628504, speed=0.5552538482374103)}


def response_table(column, z_start, spans=TABLE_SPANS, band=None, z_ground=0.):
    """follow() for each span, with the budgets on one depth grid."""
    speeds = DropSpeeds(column, np.linspace(z_ground, z_start, 14))
    rows = [follow(d, column, z_start, speeds, band=band, z_ground=z_ground) for d in spans]
    grid = rows[0]['depth_grid_m']
    stack = {k: np.array([r['record'][k] for r in rows]) for k in ('ice', 'water', 'sub', 'melt', 'evap', 'time')}
    return dict(spans=np.asarray(spans), masses=np.array([r['mass_kg'] for r in rows]), grid=grid, stack=stack,
                rows=[{k: v for k, v in r.items() if k not in ('record', 'depth_grid_m')} for r in rows],
                z_start=z_start, z_ground=z_ground)


def population(table, masses, weights):
    """Budgets of a population of flakes (masses at the freezing level, weights their mass flux), as fractions of the
    mass flux arriving at the freezing level, by depth: each flake read from the table in log mass."""
    lm = np.log(table['masses'])
    x = np.clip(np.log(masses), lm[0], lm[-1])
    j = np.clip(np.searchsorted(lm, x), 1, len(lm) - 1)
    f = ((x - lm[j - 1]) / (lm[j] - lm[j - 1]))[:, None]
    w = (weights / weights.sum())[:, None]
    out = {k: (w * (a[j - 1] * (1 - f) + a[j] * f)).sum(axis=0) for k, a in table['stack'].items() if k != 'time'}
    out['grid'] = table['grid']
    out['beyond_table'] = float(weights[(np.log(masses) > lm[-1]) | (np.log(masses) < lm[0])].sum() / weights.sum())
    return out


def depth_where(grid, values, level):
    """First depth at which a falling budget passes below level, linear between grid points; None if never."""
    below = np.flatnonzero(values <= level)
    if not below.size:
        return None
    k = below[0]
    if k == 0:
        return 0.
    return float(grid[k - 1] + (values[k - 1] - level) / (values[k - 1] - values[k]) * (grid[k] - grid[k - 1]))


def layers(pop, column, z_start):
    """Where the snow's ice goes and its latent heat, per mm/h of snowfall at the freezing level: the dry layer down to
    the first melting, the melting layer down to where 99% of the ice is gone, and the rain below; with the depths at
    which 10, 50, 90 and 99% of the ice has gone and the water reaching the ground at each depth."""
    g, ice, water = pop['grid'], pop['ice'], pop['water']
    gone = 1 - ice
    marks = {f'{int(q * 100)}': depth_where(g, ice, 1 - q) for q in (.1, .5, .9, .99)}
    onset = depth_where(g, -pop['melt'], -1e-4)
    bottom = marks['99'] if marks['99'] is not None else float(g[-1])
    edges = dict(dry=(0., onset if onset is not None else bottom), melting=(onset if onset is not None else bottom, bottom),
                 rain=(bottom, float(g[-1])))
    t_air = np.array([column.state(z_start - d)['t'] for d in g])
    lv = np.array([drops.latent_heat(max(t, drops.TRIPLE_POINT_K)) for t in t_air])
    per_mm_h = 1 / 3600.
    out = dict(depths_ice_gone_m=marks, onset_depth_m=onset)
    for name, (a, b) in edges.items():
        ia, ib = np.searchsorted(g, a), np.searchsorted(g, b)
        ib = max(ib, ia)
        def change(k):
            return float(np.interp(b, g, pop[k]) - np.interp(a, g, pop[k]))
        sub, melt, evap = change('sub'), change('melt'), change('evap')
        lv_mean = float(np.mean(lv[ia:ib + 1])) if ib > ia else float(lv[min(ia, len(lv) - 1)])
        out[name] = dict(top_depth_m=a, bottom_depth_m=b, top_height_m=z_start - a, bottom_height_m=z_start - b,
                         sublimated=sub, melted=melt, evaporated=evap,
                         heat_w_m2_per_mm_h=dict(sublimation=LS0 * sub * per_mm_h, melting=LF * melt * per_mm_h,
                                                 evaporation=lv_mean * evap * per_mm_h))
    out['water_at_depth'] = dict(grid_m=g, water=water, ice=ice)
    out['totals'] = dict(sublimated=float(pop['sub'][-1]), melted=float(pop['melt'][-1]), evaporated=float(pop['evap'][-1]),
                         ground_water=float(water[-1]), ground_ice=float(ice[-1]))
    return out


# --- the 20 October column --------------------------------------------------------------------------------------

LADDER = dict(height_m=29245.11370855616, pressure_pa=65925.26459686145, lapse_k_m=0.921346439230548e-3)  # results/snow.json
LADDER_SPANS = (0.002, 0.005, 0.01, 0.02, 0.045, 0.08, 0.11, 0.15, 0.19)


def ladder():
    """The flakes of 20 October and the large flakes in the same straight column the flake model of that day used
    (flakes.py: the resort's sky ladder about its freezing level, at 90% and 70% humidity below it)."""
    rows = []
    for rh in (0.9, 0.7):
        col = straight_sounding(MOON_SURFACE_GRAVITY, MOON_AIR, LADDER['height_m'], LADDER['pressure_pa'],
                                LADDER['lapse_k_m'], rh)
        speeds = DropSpeeds(col, np.linspace(0., LADDER['height_m'], 12))
        for d in LADDER_SPANS:
            r = follow(d, col, LADDER['height_m'], speeds)
            rows.append(dict(relative_humidity=rh, span_m=d, melted_diameter_m=r['melted_diameter_m'],
                             onset_depth_m=r['onset_depth_m'], half_melted_depth_m=r['half_melted_depth_m'],
                             melted_depth_m=r['melted_depth_m'], melted_temperature_c=r.get('melted_temperature_c'),
                             fate=r['fate'], sublimated=float(r['record']['sub'][-1]),
                             reaching_ground=(r.get('drop') or {}).get('survived_fraction')))
    return rows


# --- the lunar cases -------------------------------------------------------------------------------------------

# McSnow columns of growth and aggregation (mcsnow_growth_aggregation.py) whose arriving snow the records hold as
# converged or as the stated single runs; their flakes melt here.
POPULATIONS = {
    'light 0.3 mm/h': ('col_q0.3_f1_sp4', 'col_q0.3_f1_base'),
    'light 0.3 mm/h, a tenth of the crystals': ('col_q0.3_f0.1_base',),
    'light 0.3 mm/h from -30 C': ('col_deep_q0.3_f1_half', 'col_deep_q0.3_f1_base'),
    'moderate 1 mm/h': ('col_q1_f1_base',),
    'heavy 3 mm/h, a hundredth of the crystals': ('col_q3_f0.01_sp4', 'col_q3_f0.01_base'),
    'saturated, a hundredth of the crystals': ('col_sat_f0.01_base',),
}
POOL_H, ARRIVING_M = 12., 2000.
QUARTERS = {'sunrise to noon': (-90., 0.), 'noon to sunset': (0., 90.), 'sunset to midnight': (90., 180.),
            'midnight to sunrise': (-180., -90.)}
SAMPLE_SPANS = (0.005, 0.02, 0.045, 0.08, 0.11, 0.15, 0.19)
SAMPLE_COLUMNS = 120
GROUND_HEIGHTS_M = (0., 1000., 2000., 3000., 5000., 8000., 11600.)


def mcsnow_population(name):
    """Masses and mass-flux weights of the flakes in a McSnow column's lowest 2 km over its last 12 hours."""
    from atmosphere.snow import mcsnow_growth_aggregation as growagg
    snaps = growagg.read_snapshots(name)
    t_end = snaps[-1][0]
    chosen = [a for s, a in snaps if t_end - POOL_H * 3600. - 1. <= s <= t_end]
    a = np.concatenate(chosen)
    z, m, d, mono, xi, v = a.T
    ok = (z < ARRIVING_M) & (m > 0) & np.isfinite(v) & (v > 1e-6) & (v < 50.)
    return m[ok], xi[ok] * m[ok] * v[ok], d[ok]


def quarter_of(hour):
    for name, (a, b) in QUARTERS.items():
        if a <= hour < b:
            return name
    return 'midnight to sunrise'


def summarise_rows(table):
    keep = ('span_m', 'mass_kg', 'melted_diameter_m', 'onset_depth_m', 'onset_temperature_c', 'half_melted_depth_m',
            'melted_depth_m', 'melted_temperature_c', 'fate')
    out = []
    for r in table['rows']:
        row = {k: r.get(k) for k in keep}
        if 'drop' in r:
            row.update(drop_status=r['drop']['status'], drop_survived=r['drop']['survived_fraction'],
                       drop_pieces=r['drop'].get('pieces'),
                       drop_screen_m=r['drop']['screen_m'], drop_ground_diameter_m=r['drop'].get('ground_diameter_m'))
        out.append(row)
    return out


def water_by_ground(pop, z_start, heights=GROUND_HEIGHTS_M):
    """Water and ice reaching a ground at each height (m above the box's surface), as fractions of the snowfall."""
    g = pop['grid']
    out = {}
    for h in heights:
        d = z_start - h
        if d < 0:
            continue
        out[f'{h / 1000:g}_km'] = dict(water=float(np.interp(d, g, pop['water'])), ice=float(np.interp(d, g, pop['ice'])))
    return out


def vapour_depths(pop):
    """Depths by which a quarter, half, three-quarters and nine-tenths of the vapour the snow gives the air (its
    sublimation and evaporation together) has been released, and the same for the drops' evaporation alone."""
    out = {}
    for label, series in (('all_vapour', pop['sub'] + pop['evap']), ('evaporation', pop['evap'])):
        final = series[-1]
        out[label] = ({f'{int(q * 100)}': depth_where(pop['grid'], -series / final, -q) for q in (.25, .5, .75, .9)}
                      if final > 1e-6 else None)
    return out


def moon(samples=SAMPLE_COLUMNS, seed=20261025):
    from atmosphere.snow import melting_environment as env
    RAW.mkdir(parents=True, exist_ok=True)
    from atmosphere.snow import mcsnow_growth_aggregation as growagg
    pops = {}
    for label, names in POPULATIONS.items():
        for name in names:
            if not (growagg.CACHE / 'runs' / name / 'run.json').exists():
                raise FileNotFoundError(f'McSnow run {name} has not finished')
            m, w, d = mcsnow_population(name)
            pops[(label, name)] = (m, w, d)
    out = dict(boxes={}, populations={})
    for (label, name), (m, w, d) in pops.items():
        d_eq = (6 * m / (math.pi * RHO_W)) ** (1 / 3)
        out['populations'][name] = dict(label=label, superparticles=int(len(m)),
                                        span_mm={f'p{q}': float(weighted_quantile(d, w, q / 100) * 1e3) for q in (50, 90, 99)},
                                        largest_span_mm=float(d.max() * 1e3),
                                        melted_diameter_mm={f'p{q}': float(weighted_quantile(d_eq, w, q / 100) * 1e3) for q in (50, 90, 99)})
    rng = np.random.default_rng(seed)
    for box in ('land', 'sea'):
        data, record = env.load(box)
        quarters = np.array([quarter_of(h) for h in data['hour_deg']])
        composites = {'median': composite_air(data)}
        composites['moist'] = composite_air(data, q=.9)
        composites['dry'] = composite_air(data, q=.1)
        for qname in QUARTERS:
            sel = quarters == qname
            if data['snow_kg_m3'][sel].sum() > 0 and sel.sum() > 200:
                composites[qname] = composite_air(data, select=sel)
        bx = dict(environment=dict(record=record, columns=int(len(data['z0_m'])),
                                   snow_share_by_quarter={q: float(data['snow_kg_m3'][quarters == q].sum() / data['snow_kg_m3'].sum())
                                                          for q in QUARTERS},
                                   freezing_level_m={f'p{q}': float(weighted_quantile(data['z0_m'], data['snow_kg_m3'], q / 100))
                                                     for q in (10, 50, 90)}),
                  composites={}, sampled=None)
        for cname, (air, zf) in composites.items():
            bands = BANDS if cname == 'median' else {'central': BANDS['central']}
            entry = dict(freezing_level_m=zf, profile=[dict(depth_m=dd, temperature_c=float(air.state(zf - dd)['t'] - T0),
                                                            relative_humidity=float(air.state(zf - dd)['relative_humidity']),
                                                            pressure_hpa=float(air.state(zf - dd)['p'] / 100))
                                                       for dd in (0., 500., 1000., 2000., 3000., 4000., 6000., 8000., 12000., 16000., 20000.)
                                                       if dd <= zf],
                         bands={})
            for bname, band in bands.items():
                table = response_table(air, zf, band=band)
                np.savez(RAW / f'table_{box}_{cname.replace(" ", "_")}_{bname}.npz', spans=table['spans'], grid=table['grid'],
                         **{k: v for k, v in table['stack'].items()})
                b = dict(flakes=summarise_rows(table), populations={})
                for (label, name), (m, w, d) in pops.items():
                    pop = population(table, m, w)
                    lay = layers(pop, air, zf)
                    b['populations'][name] = dict(label=label, beyond_table=pop['beyond_table'],
                                                  layers={k: lay[k] for k in ('dry', 'melting', 'rain')},
                                                  depths_ice_gone_m=lay['depths_ice_gone_m'], onset_depth_m=lay['onset_depth_m'],
                                                  totals=lay['totals'], ground=water_by_ground(pop, zf),
                                                  vapour_release_depths_m=vapour_depths(pop),
                                                  profile=dict(depth_m=pop['grid'][::10].tolist(), ice=pop['ice'][::10].round(5).tolist(),
                                                               water=pop['water'][::10].round(5).tolist(),
                                                               sublimated=pop['sub'][::10].round(5).tolist(),
                                                               melted=pop['melt'][::10].round(5).tolist(),
                                                               evaporated=pop['evap'][::10].round(5).tolist()))
                entry['bands'][bname] = b
                print(box, cname, bname, 'done', flush=True)
            bx['composites'][cname] = entry
        # columns drawn by their snow, each melting a few spans
        p = data['snow_kg_m3'] / data['snow_kg_m3'].sum()
        picks = rng.choice(len(p), size=samples, replace=True, p=p)
        rows = []
        for i in picks:
            air = column_air(data, i)
            zf = float(data['z0_m'][i])
            speeds = DropSpeeds(air, np.linspace(0., zf, 12))
            for dspan in SAMPLE_SPANS:
                r = follow(dspan, air, zf, speeds)
                rows.append(dict(column=int(i), span_m=dspan, hour_deg=float(data['hour_deg'][i]), freezing_level_m=zf,
                                 onset_depth_m=r['onset_depth_m'], melted_depth_m=r['melted_depth_m'],
                                 melted_temperature_c=r.get('melted_temperature_c'), fate=r['fate'],
                                 survived=(r.get('drop') or {}).get('survived_fraction'),
                                 evaporated_by_m=_evaporated_at(r)))
        (RAW / f'sampled_{box}.json').write_text(json.dumps(rows) + '\n')
        bx['sampled'] = summarise_sample(rows)
        out['boxes'][box] = bx
        print(box, 'sampled columns done', flush=True)
    (RAW / 'moon.json').write_text(json.dumps(out, default=_plain) + '\n')
    return out


def _evaporated_at(r):
    """Depth below the freezing level where the drop vanished, if it evaporated aloft."""
    if r['fate'] != 'evaporated':
        return None
    water = r['record']['water']
    k = np.flatnonzero((water <= 0) & (np.arange(len(water)) > 0) & (r['record']['ice'] <= 0))
    return float(r['depth_grid_m'][k[0]]) if k.size else None


def summarise_sample(rows):
    out = {}
    for dspan in SAMPLE_SPANS:
        sel = [r for r in rows if r['span_m'] == dspan]
        def q(key):
            v = np.array([r[key] for r in sel if r[key] is not None], float)
            return dict(p10=float(np.percentile(v, 10)), p50=float(np.percentile(v, 50)), p90=float(np.percentile(v, 90)),
                        n=int(v.size)) if v.size else None
        fates = {f: sum(r['fate'] == f for r in sel) for f in {r['fate'] for r in sel}}
        surv = np.array([r['survived'] if r['survived'] is not None else 0. for r in sel])
        out[f'{dspan * 1000:g}_mm'] = dict(onset_depth_m=q('onset_depth_m'), melted_depth_m=q('melted_depth_m'),
                                            melted_temperature_c=q('melted_temperature_c'), fates=fates,
                                            water_reaching_ground=dict(mean=float(surv.mean()), p10=float(np.percentile(surv, 10)),
                                                                       p50=float(np.percentile(surv, 50)), p90=float(np.percentile(surv, 90))),
                                            evaporated_depth_m=q('evaporated_by_m'))
    return out


def _plain(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    raise TypeError(type(x))


# --- above the freezing level: the storm paths ----------------------------------------------------------------

PATH_RUNS = ('path_archived_c0.01_base', 'path_archived_c0.01_sp4', 'path_archived_c0.1_base', 'path_archived_c0.1_sp4',
             'path_archived_c1_base', 'path_saturated_c0.01_base', 'path_saturated_c0.01_sp4',
             'path_saturated_c0.1_base', 'path_saturated_c1_base')


def path_budget(name):
    """Ice the storm-path columns of mcsnow_growth_aggregation.py grew, delivered to 0 C, held aloft at the end and
    lost to sublimation before the freezing level.  Each path column's ice aloft is read every 30 minutes; its net
    change plus the snow it delivered in between is its net growth, and the sums of the positive and negative parts
    bound the ice grown and the ice sublimated from below (growth and sublimation in one column within one interval
    cancel)."""
    from atmosphere.snow import mcsnow_growth_aggregation as growagg
    record = json.loads((growagg.CACHE / 'runs' / name / 'run.json').read_text())
    grid = record['configuration']['grid_nml']
    n = grid['nz'] // growagg.N_SUB
    top = grid['dom_top']
    snaps = growagg.read_snapshots(name)
    arr = growagg.read_arrivals(name)
    t_arr, k_arr, _, m_arr, _, _, xi_arr, reason = arr.T
    warm = reason == 1
    times = np.array([s for s, _ in snaps])
    aloft = np.zeros((len(snaps), n))
    for i, (_, a) in enumerate(snaps):
        z, m, d, mono, xi, v = a.T
        c = np.clip(((top - z) // growagg.COLUMN_DZ).astype(int) // growagg.N_SUB, 0, n - 1)
        aloft[i] = np.bincount(c, weights=xi * m, minlength=n)
    left = np.zeros_like(aloft)                       # every superparticle that left a column, by any reason
    for i, t in enumerate(times):
        s = t_arr <= t + 1e-6
        left[i] = np.bincount((k_arr[s] - 1).astype(int), weights=xi_arr[s] * m_arr[s], minlength=n)[:n]
    net = aloft + left
    step = np.diff(np.vstack([np.zeros(n), net]), axis=0)
    grown, sublimated = np.clip(step, 0, None).sum(), np.clip(-step, 0, None).sum()
    sinks = {label: float(xi_arr[reason == r] @ m_arr[reason == r]) for r, label in
             ((1, 'delivered_kg'), (2, 'below_column_kg'), (3, 'censored_kg'))}
    columns_delivering = int(sum((warm & (k_arr == c + 1)).sum() >= 3 for c in range(n)))
    share = lambda x: float(x / grown) if grown > 0 else None
    return dict(paths=n, paths_delivering=columns_delivering, grown_kg=float(grown), sublimated_kg=float(sublimated),
                aloft_at_end_kg=float(aloft[-1].sum()), **sinks,
                sublimated_share_of_grown=share(sublimated), delivered_share_of_grown=share(sinks['delivered_kg']),
                censored_share_of_grown=share(sinks['censored_kg']), aloft_share_of_grown=share(aloft[-1].sum()),
                closure_kg=float(grown - sublimated - sum(sinks.values()) - aloft[-1].sum()),
                snapshot_interval_s=float(np.median(np.diff(times))), duration_h=float(times[-1] / 3600.))


# --- the Earth check ---------------------------------------------------------------------------------------------

EARTH_FREEZING = dict(height_m=1500., pressure_pa=85000.)        # the winter column of flakes.py
EARTH_SPANS = (0.002, 0.005, 0.01, 0.02)
DENSE_FACTOR = 0.0061e-3 * 100 ** 2.1 / AGG_ALPHA                # Heymsfield and co-authors' (2021) m-D, cgs 0.0061 D^2.1
# What the check compares with (melting_sources.json): Matsuo and Sasyo's (1981) computed depths below the freezing
# level where melting begins at 6 K/km; Heymsfield and co-authors' (2021) ice-bulb fits; the observed temperatures
# by which melting completes.
MATSUO_SASYO_ONSET_M = {0.9: 120., 0.8: 250., 0.7: 390., 0.5: 700.}


def heymsfield_onset_c(rh, pressure_hpa):
    """Air temperature where the ice-bulb temperature is 0 C, from Heymsfield and co-authors' (2021) linear fits at
    1000 and 500 hPa (their eqs. 2a, 2b), linear in pressure between them."""
    slope = 9.3 + (12.67 - 9.3) * (pressure_hpa - 500.) / 500.
    return (100. * (1 - rh)) / slope


def earth():
    rows = []
    for lapse in (6.0e-3, 6.5e-3):
        for rh in (1.0, 0.95, 0.9, 0.8, 0.7, 0.5):
            col = straight_sounding(STANDARD_GRAVITY, EARTH_AIR, EARTH_FREEZING['height_m'], EARTH_FREEZING['pressure_pa'],
                                    lapse, rh)
            speeds = DropSpeeds(col, np.linspace(0., EARTH_FREEZING['height_m'], 6))
            for dense in (1., DENSE_FACTOR):
                for d in EARTH_SPANS:
                    m0 = aggregate_mass(d) * dense
                    r = follow(d, col, EARTH_FREEZING['height_m'], speeds, m0=m0)
                    rows.append(dict(lapse_k_km=lapse * 1e3, relative_humidity=rh, mass_factor=dense, span_m=d,
                                     melted_diameter_m=r['melted_diameter_m'], onset_depth_m=r['onset_depth_m'],
                                     onset_temperature_c=r.get('onset_temperature_c'), half_melted_depth_m=r['half_melted_depth_m'],
                                     melted_depth_m=r['melted_depth_m'], melted_temperature_c=r.get('melted_temperature_c'),
                                     fate=r['fate'], fall_speed_m_s=dry_speed(m0, col.state(EARTH_FREEZING['height_m']),
                                                                               STANDARD_GRAVITY)[0]))
    p_hpa = EARTH_FREEZING['pressure_pa'] / 100.
    onset = []
    for rh, ms in MATSUO_SASYO_ONSET_M.items():
        r = next(x for x in rows if x['lapse_k_km'] == 6.0 and x['relative_humidity'] == rh and x['mass_factor'] == 1.
                 and x['span_m'] == 0.01)
        onset.append(dict(relative_humidity=rh, model_onset_depth_m=r['onset_depth_m'], matsuo_sasyo_m=ms,
                          model_onset_c=r['onset_temperature_c'], heymsfield_fit_c=heymsfield_onset_c(rh, p_hpa)))
    return dict(freezing_level=EARTH_FREEZING, rows=rows, onset_comparison=onset)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('earth', 'ladder', 'moon', 'paths', 'assemble'))
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    if args.action == 'earth':
        out = earth()
        (RAW / 'earth.json').write_text(json.dumps(out, default=_plain) + '\n')
        for r in out['onset_comparison']:
            print(r)
    elif args.action == 'ladder':
        rows = ladder()
        (RAW / 'ladder.json').write_text(json.dumps(rows, default=_plain) + '\n')
        for r in rows:
            print(r)
    elif args.action == 'moon':
        moon()
    elif args.action == 'paths':
        out = {name: path_budget(name) for name in PATH_RUNS}
        (RAW / 'paths.json').write_text(json.dumps(out) + '\n')
        for k, v in out.items():
            print(k, v['paths_delivering'], round(v['sublimated_share_of_grown'], 3))
    elif args.action == 'assemble':
        product = assemble()
        print('wrote', PRODUCT.relative_to(ROOT), len(json.dumps(product)), 'bytes')



# --- the product -----------------------------------------------------------------------------------------------

def cm1_processes(box):
    """CM1's own melting and rain evaporation in the snowy columns by depth below each column's freezing level: the
    snow-weighted mean cooling rate (K/day of potential temperature), and the depths above which half and nine-tenths
    of the snow-weighted, density-weighted cooling lies (Morrison's snow and graupel melting, its rain's evaporation)."""
    from atmosphere.snow import melting_environment as env
    from shared.constants import CM1_DRY_AIR_GAS_CONSTANT
    data, _ = env.load(box)
    zh, z0, w = data['zh_m'], data['z0_m'], data['snow_kg_m3']
    depths = np.arange(0., 22000.1, 250.)
    rows = np.arange(len(z0))
    rho = data['p'].astype(float) / (CM1_DRY_AIR_GAS_CONSTANT * data['T'].astype(float))
    prof = {}
    for key, label in (('melt', 'melting'), ('evar', 'rain_evaporation')):
        a = data[key].astype(float)
        rate, heat = [], []
        for d in depths:
            z = z0 - d
            ok = z >= zh[0]
            idx = np.clip(np.searchsorted(zh, z), 1, len(zh) - 1)
            f = (z - zh[idx - 1]) / (zh[idx] - zh[idx - 1])
            v = np.where(ok, a[rows, idx - 1] * (1 - f) + a[rows, idx] * f, 0.)
            r = np.where(ok, rho[rows, idx - 1] * (1 - f) + rho[rows, idx] * f, 0.)
            rate.append(float((v * w).sum() / w.sum()))
            heat.append(float((v * r * w).sum() / w.sum()))
        cum = np.cumsum(np.clip(-np.array(heat), 0., None))          # cooling, made positive
        prof[label] = dict(depth_m=depths[::4].tolist(), mean_k_day=(np.array(rate[::4]) * 86400.).round(4).tolist(),
                           half_depth_m=float(np.interp(.5 * cum[-1], cum, depths)) if cum[-1] != 0 else None,
                           nine_tenths_depth_m=float(np.interp(.9 * cum[-1], cum, depths)) if cum[-1] != 0 else None)
    return prof


def assemble():
    from atmosphere.snow import melting_environment as env
    moon_raw = json.loads((RAW / 'moon.json').read_text())
    keep = ('label', 'totals', 'depths_ice_gone_m', 'onset_depth_m', 'ground')
    for box in moon_raw['boxes'].values():
        # full detail for the median column; its band edges without depth profiles; the humidity and time-of-day
        # composites with each population's totals, depths and ground shares, and every fourth flake span
        for cname, comp in box['composites'].items():
            for bname, band in comp['bands'].items():
                for name, pop in band['populations'].items():
                    if cname != 'median':
                        band['populations'][name] = {k: pop[k] for k in keep}
                    elif bname != 'central':
                        pop.pop('profile', None)
                if cname != 'median':
                    band['flakes'] = band['flakes'][1::4]
                elif bname != 'central':
                    band['flakes'] = band['flakes'][1::2]
    ladder_raw = json.loads((RAW / 'ladder.json').read_text())
    earth_raw = json.loads((RAW / 'earth.json').read_text())
    paths_raw = json.loads((RAW / 'paths.json').read_text())
    files = ('atmosphere/snow/melting.py', 'atmosphere/snow/melting_environment.py', 'atmosphere/snow/boehm.py',
             'atmosphere/rain/drops.py')
    producer = dict(domain='atmosphere', files={f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest()[:16] for f in files})
    from atmosphere.snow import mcsnow_growth_aggregation as growagg
    runs = {}
    for names in POPULATIONS.values():
        for n in names:
            runs[n] = hashlib.sha256((growagg.CACHE / 'runs' / n / 'run.json').read_bytes()).hexdigest()[:16]
    for n in PATH_RUNS:
        runs[n] = hashlib.sha256((growagg.CACHE / 'runs' / n / 'run.json').read_bytes()).hexdigest()[:16]
    product = dict(
        schema=SCHEMA, producer=producer,
        evidence=('Runnable single-flake calculations: Mitra and co-authors\' (1990) wind-tunnel-based melting model as '
                  'Frick and co-authors (2013) give it, McSnow\'s aggregates and fall speeds, the rain solver\'s drops, '
                  'through the archived air of the design-climate CM1 boxes; checked on Earth against observed and '
                  'computed melting depths; no lunar validation.'),
        reading_rule=('Depths are metres below each column\'s freezing level; heights are metres above the box\'s '
                      'surface (sea level for the sea box, the box\'s own ground for the land box). Budgets are fractions '
                      'of the snow mass flux arriving at the freezing level; heat is W/m2 for each mm/h of that snowfall '
                      '(water equivalent). Composites take the snow-weighted median of temperature and vapour pressure '
                      'at each depth below the freezing level; moist and dry take the 90th and 10th percentiles of '
                      'vapour pressure. Sampled columns are drawn by their snow at the freezing level.'),
        assumptions=['Still air: the boxes\' vertical motions are left out.',
                     'Flakes fall alone: no aggregation, collision or breakup of flakes below the freezing level.',
                     'A flake melts to one drop and drops fall alone; a drop above its stability screen breaks into '
                     'eight drops of half its diameter, again until they hold together.',
                     'The melting flake\'s vapour comes from exposed ice in proportion to the ice\'s share of its mass; '
                     'that split moves the record between sublimation, melting and evaporation and leaves the masses '
                     'and the heat taken from the air unchanged.',
                     'Vapour diffusivity is the rain solver\'s, Massman\'s (1998) measured fit since 25 October.',
                     'Gravity holds its surface value through the column, as in CM1.'],
        constants=dict(latent_heat_fusion_j_kg=LF, latent_heat_sublimation_0c_j_kg=LS0, latent_heat_evaporation_0c_j_kg=LV0,
                       moon_gravity_m_s2=MOON_SURFACE_GRAVITY, earth_gravity_m_s2=STANDARD_GRAVITY,
                       aggregate_mass_law_si=[AGG_ALPHA, AGG_BETA], aggregate_area_law_si=[AGG_SIGMA, AGG_GAMMA],
                       mitra=dict(axis_ratio_dry=AXIS_DRY, capacitance_factor_dry=CAPACITANCE_DRY, psi_linear=PSI_LINEAR,
                                  snow_density_cap=SNOW_DENSITY_CAP),
                       vapour_diffusivity_massman=[drops.MASSMAN_DIFFUSIVITY, drops.MASSMAN_EXPONENT], bands=BANDS),
        inputs=dict(environment={b: json.loads((env.ENV_DIR / f'{b}.json').read_text()) for b in ('land', 'sea')},
                    mcsnow_runs=runs),
        earth_check=earth_raw,
        ladder_column=ladder_raw,
        above_freezing_level=paths_raw,
        moon=moon_raw,
        cm1_processes={b: cm1_processes(b) for b in ('land', 'sea')},
        sources='atmosphere/snow/melting_sources.json')
    product = _rounded(json.loads(json.dumps(product, default=_plain)))
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def _rounded(x, digits=4):
    """Floats to four significant figures, the product's precision; strings, integers and hashes as they are."""
    if isinstance(x, dict):
        return {k: _rounded(v, digits) for k, v in x.items()}
    if isinstance(x, list):
        return [_rounded(v, digits) for v in x]
    if isinstance(x, float) and math.isfinite(x) and x != 0.:
        return float(f'{x:.{digits}g}')
    return x


if __name__ == '__main__':
    main()
