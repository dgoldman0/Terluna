"""Atomic oxygen and hydrogen made in the Open Moon's upper air, and what of them escapes.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.oxygen    # writes results/oxygen_escape.json

R2 lists atomic oxygen, photochemistry below the exobase and hydrogen from water among the channels the loss response
leaves out. Behind the titania stack nothing splits O2 below the 0.3 Pa base, but above it the sky's Lyman-alpha glow
and the light that passes gaps do. The glow alone breaks enough O2 to make several kilograms a second of oxygen atoms,
and an oxygen atom, half the mass of the molecules, escapes far more readily, while no magnetosphere holds a neutral.
This step follows the atoms, and the hydrogen that water and H2 give up, through the whole column at a state of the
loss response:

- the column runs from the ground through the middle atmosphere's levels (its stored profile and composition) to the
  0.3 Pa base, and on through the thermal column's solution at the state to the exobase;
- its chemistry is the middle atmosphere's (chemistry.py: fourteen O, H and N species with the JPL rates), with the
  photolysis of sunlight above 202 nm through the window (photolysis.py) and, added here, O2 and water split by the
  glow (isotropic Lyman-alpha absorbed by O2) and by the light the limb tables trace above the base (the far
  ultraviolet and Lyman-alpha one photon to one O2, the extreme ultraviolet and X-rays one ion pair per 35 eV, each
  ending as two oxygen atoms through O2+ recombining);
- species move by eddy diffusion (the middle atmosphere's mixing scaled for the Moon, or above the base that of
  breaking gravity waves, gravity_waves.py, among the setups bounding it) and by molecular diffusion,
  which carries each species toward its own scale height above the homopause (binary diffusion in N2: Banks and
  Kockarts 1973 for O, H and H2, a hard-sphere scaling for the rest), in spherical shells, with the middle
  atmosphere's conditions at the ground;
- O escapes at the exobase by Jeans effusion raised by Earth's tide (tides.py), H and H2 by Jeans effusion; the
  oxygen atoms in the sunlit exosphere outside the shadow are photoionized and charge-exchanged by the solar wind as in
  the exosphere step, and with the September magnets their ions are kept in the share the step finds for the
  molecules;
- atoms made within a collision length or so of the exobase leave hot. Photolysis gives each atom up to 2.5 eV and O2+
  recombining up to 3.5 eV, mostly above the 0.47 eV x (lunar radius / r) an oxygen atom needs to escape, so they
  leave on straight paths. In the unattenuated top a light breaks a molecule at its cross-section times the actinic flux,
  twice the flux onto the column (isotropic glow, or sunlight over the sphere), and an atom made at collision depth s
  escapes with chance E2(s) / 2, which sums to 1 / (4 sigma_c) of collision column; so of the atoms a light makes in a
  column it crosses whole, a share sigma x / (2 sigma_c) leaves hot, with sigma the largest cross-section in the
  light's group of bands, x its absorber's share of the gas at the exobase and sigma_c the thermal column's collision
  cross-section. The atoms all counted above the escape energy and the largest cross-sections make it an upper
  estimate; letting an atom keep its escape energy through two collisions triples it (the upper value).

The atoms do not change the column they move in: the thermal column stays molecular and radiates with the base's
atomic oxygen only, so this step also reports the oxygen atoms' share at the exobase, which shows where that holds.
Left out: O(1D) from the glow's photolysis (all of it is counted as ground-state atoms, which survive longer), the
traced light below the base, ions in the chemistry, the infrared cooling by the atoms made above the base (their
63-um line and their collisions with CO2) and the day-night cycle (rates are diurnal means).
"""
from __future__ import annotations
import argparse
import csv
import dataclasses
import hashlib
import json
import math
from multiprocessing import Pool
import os
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.special import expn

from shared.constants import AVOGADRO, BOLTZMANN, MOON_GM, MOON_RADIUS
from shared.provenance import constants_used
from atmosphere.thermal_column import solve_column
from atmosphere.radiative_convective import thermodynamics as th, climate as cl
from atmosphere.middle_atmosphere import chemistry as ch, photolysis as ph, escape, run as middle_run
from atmosphere.loss_response import model as lr, traced, tides, absorption as ab, exosphere as ex, cycle
from atmosphere.loss_response import gravity_waves as gw

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'oxygen_escape.json'
SCHEMA = 'terluna.atmosphere.oxygen-escape/1'
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
AMU_G = 1.66054e-24
MASS_AMU = dict(O=15.999, O1D=15.999, O3=47.998, H=1.008, OH=17.007, HO2=33.006, H2O2=34.014, H2=2.016,
                NO=30.006, NO2=46.005, NO3=62.004, N2O5=108.01, HNO3=63.012, N2O=44.013)
# Binary diffusion in N2, b = n D = A T^s (cm^-1 s^-1): Banks and Kockarts (1973) for O, H and H2; the other species
# scale from O as hard spheres, b ~ (reduced mass)^-1/2 / (collision diameter)^2, with the diameters below (angstrom).
DIFFUSION = dict(O=(9.69e16, 0.774), H=(4.87e17, 0.698), H2=(2.73e17, 0.741))
DIAMETER_A = dict(O=3.0, O1D=3.0, O3=4.0, OH=3.1, HO2=3.6, H2O2=4.0, NO=3.5, NO2=4.0, NO3=4.4, N2O5=5.0, HNO3=4.5,
                  N2O=3.9)
N2_DIAMETER_A, N2_AMU = 3.75, 28.013
ESCAPING = ('O', 'H', 'H2')
# The largest cross-section (cm^2) in each of the limb tables' groups of bands (0-10, 10-102.7, 102.7-121, 121-122 and
# 122-175 nm) and what absorbs it: the air's photoabsorption at 10 nm and its extreme-ultraviolet peak, O2's bands
# below Lyman-alpha, Lyman-alpha, and the Schumann-Runge continuum's peak at 142 nm.
HOT_CROSS_SECTION_CM2 = (3e-18, 3e-17, 3e-17, escape.O2_LYMAN_ALPHA_CM2, 1.5e-17)
HOT_ABSORBER = ('air', 'air', 'O2', 'O2', 'O2')
HOT_COLLISIONS = 3                 # the upper value: an atom keeps its escape energy through two collisions
LYMAN_ALPHA_H2O_CM2 = 1.59e-17     # water at 121.6 nm (Kley 1984, J. Atmos. Chem. 2, 203), taken to split into H + OH
EV_PER_ION_PAIR = 35.0             # the extreme ultraviolet and X-rays: one ion pair per 35 eV absorbed
PHOTON_GROUPS = (2, 3, 4)          # the limb tables' groups split O2 one photon at a time (102.7-121 nm, Lyman-alpha, FUV)
LAYERS_PER_EFOLD = 4               # levels per e-fold of pressure above the base
MAX_RADIUS_R = 4.0                 # the ring fleet's protected radius
LEVEL = 'standard'
J = dict(O2='O2 -> O + O')


@dataclasses.dataclass(frozen=True)
class Setup:
    """A run's choices: the eddy mixing's scale over the middle atmosphere's (1 keeps it); above the base, the
    middle atmosphere's mixing ('middle') or the breaking gravity waves' (gravity_waves.py, 'gravity_waves', with
    the waves' horizontal wavelengths stretched by wave_stretch and the calibration times wave_factor_scale); the
    water above the base ('middle': the middle atmosphere's top value; 'dry': none) and the levels per e-fold above
    the base."""
    mixing_scale: float = 1.0
    mixing_above_base: str = 'middle'
    wave_stretch: float = gw.STRETCH
    wave_factor_scale: float = 1.0
    water_above_base: str = 'middle'
    layers_per_efold: int = LAYERS_PER_EFOLD


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def middle_case(shield, treatment):
    name = lr.SHIELDS[shield] + lr.TREATMENTS[treatment]
    return middle_run.cases()[name]


def middle_layers(case):
    """The middle atmosphere's stored layers (surface first): pressure, temperature, water and each species' mixing
    ratio; and its summary."""
    rows = list(csv.DictReader(open(MIDDLE / 'profiles' / f'{case.name}.csv', newline='')))
    get = lambda k: np.array([float(r[k]) for r in rows])
    summary = json.loads((MIDDLE / 'cases' / f'{case.name}.json').read_text())
    return dict(p=get('p_pa'), t=get('t_k'), h2o=get('h2o'), mixing={s: get(s) for s in ch.SPECIES}), summary


def levels(case, layers, profile, setup):
    """Level pressures, temperatures and water (surface first): the middle atmosphere's levels below the base, then
    the thermal column's solution from the base to the exobase at LAYERS_PER_EFOLD levels per e-fold."""
    base = lr.BASE_PA
    p_mid = middle_run.eq.Profile(case).p
    keep = p_mid > base * 1.0001
    lp = np.log(layers['p'][::-1])
    t_mid = np.interp(np.log(p_mid[keep]), lp, layers['t'][::-1])
    x_mid = np.interp(np.log(p_mid[keep]), lp, layers['h2o'][::-1])
    p_col, t_col = profile['pressure_Pa'], profile['temperature_K']
    depth = math.log(base / p_col[-1])
    steps = np.r_[np.arange(0.0, depth, 1.0 / setup.layers_per_efold), depth]
    p_up = base * np.exp(-steps)
    t_up = np.interp(np.log(p_up)[::-1], np.log(p_col)[::-1], t_col[::-1])[::-1]
    top_water = float(layers['h2o'][-1]) if setup.water_above_base == 'middle' else 0.0
    x_up = np.full(p_up.size, top_water)
    return np.r_[p_mid[keep], p_up], np.r_[t_mid, t_up], np.r_[x_mid, x_up]


def column(case, layers, summary, profile, setup):
    p, t, x = levels(case, layers, profile, setup)
    col = th.column_from_levels(case.planet, case.air(), p, t, x)
    col['tropopause_pa'] = float(summary['tropopause_pa'])
    col['above_base'] = col['layer_p_pa'] < lr.BASE_PA
    return col


def initial_density(col, layers):
    """The middle atmosphere's mixing ratios at each layer (its top values above it), as densities (cm^-3)."""
    m = th.layer_number_density(col)
    lp = np.log(layers['p'][::-1])
    n = np.zeros((len(ch.SPECIES), m.size))
    for i, s in enumerate(ch.SPECIES):
        f = np.interp(np.log(col['layer_p_pa']), lp, layers['mixing'][s][::-1])
        n[i] = np.maximum(f, 1e-30) * m
    return n


def glow_rates(col, activity):
    """O2 and H2O photolysis (s^-1 per layer) by the sky's Lyman-alpha glow: isotropic light from above absorbed by O2,
    the flux at vertical optical depth tau 2 pi I E3(tau) and the actinic flux 2 pi I E2(tau)."""
    intensity = lr.GLOW_RAYLEIGH * 1e10 / (4 * math.pi) * 1e-4 * lr.activity_of(activity)['glow']   # cm^-2 s^-1 sr^-1
    n_tot = th.layer_number_density(col)
    o2 = col['dry_fractions']['O2'] * (1 - col['layer_x_h2o']) * n_tot
    tau_layer = escape.O2_LYMAN_ALPHA_CM2 * o2 * col['layer_dz_cm']
    above = np.r_[np.cumsum(tau_layer[::-1])[::-1], 0.0]          # at levels, surface first
    flux = 2 * math.pi * intensity * expn(3, above)
    absorbed = flux[1:] - flux[:-1]                                # photons cm^-2 s^-1 absorbed in each layer
    j_o2 = absorbed / np.maximum(o2 * col['layer_dz_cm'], 1e-300)
    middle = 0.5 * (above[1:] + above[:-1])
    j_h2o = LYMAN_ALPHA_H2O_CM2 * 2 * math.pi * intensity * expn(2, np.maximum(middle, 1e-12))
    return j_o2, j_h2o


def group_photon_energy(design):
    """Mean photon energy (J) of the quiet Sun's light in each group of bands the limb tables follow in height."""
    w, f = escape.whi_quiet_sun()
    edges = design['shape_groups_nm']
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (w >= lo) & (w < hi)
        photons = (f[sel] * w[sel] * 1e-9 / (lr.PLANCK * lr.LIGHT)).sum()
        out.append(float(f[sel].sum() / photons) if photons > 0 else float('nan'))
    return out


def traced_production(shield, treatment, activity, gaps, heat, protected_R):
    """Oxygen atoms the traced light makes above the base (per s per m^2 of lunar surface) against log pressure above
    the base (traced.SHAPE_X), the share of them from Lyman-alpha, and all of them by group of bands: each part of the
    light at the state's heat, between the tabulated heats, by group of bands."""
    design = traced.limb_product()['design']
    rows = [e for e in traced.entries(f'{shield}_{treatment}') if 'radii' in e and not e['outflow_limit']]
    q = np.array([e['heat_W_m2'] for e in rows])
    bands, groups = design['bands_nm'], design['shape_groups_nm']
    levels_q = np.asarray(design['shape_quantiles'], float)
    s = traced.scales(activity, bands)
    mids = .5 * (np.asarray(bands[:-1]) + np.asarray(bands[1:]))
    group = np.searchsorted(groups, mids, side='right') - 1
    energy = group_photon_energy(design)
    key = f'{float(protected_R):g}'
    j = int(np.clip(np.searchsorted(q, heat), 1, len(q) - 1))
    w = float(np.clip((heat - q[j - 1]) / (q[j] - q[j - 1]), 0.0, 1.0))
    pieces = (('window', 'window', 'window', 1.0), (traced.FILM, 'annulus', f'annulus_{key}', 1.0),
              ('open', 'outside', f'outside_{key}', 1.0), ('open', 'aperture', f'aperture_{key}', gaps))
    below = np.zeros_like(traced.SHAPE_X)
    lyman = np.zeros_like(traced.SHAPE_X)
    by_group = np.zeros(len(groups) - 1)
    for i, weight in ((j - 1, 1.0 - w), (j, w)):
        for source, zone, place, factor in pieces:
            heats = np.bincount(group, weights=np.asarray(rows[i]['radii'][key][source][zone]) * s,
                                minlength=len(groups) - 1) * factor
            for g, quantiles in enumerate(rows[i]['shapes'][place][source]):
                if heats[g] <= 0 or quantiles is None:
                    continue
                absorbed = heats[g] / escape.HEATING_EFFICIENCY                     # W per m^2 of lunar surface
                per_joule = (2.0 / energy[g]) if g in PHOTON_GROUPS else 2.0 / (EV_PER_ION_PAIR * 1.602177e-19)
                share = traced._share_below(quantiles, levels_q, traced.SHAPE_X)
                made = weight * absorbed * per_joule
                below += made * share
                by_group[g] += made
                if groups[g] == 121.0:
                    lyman += made * share
    return below, lyman, by_group


def traced_rates(col, shield, treatment, activity, gaps, heat, protected_R):
    """O2 and H2O photolysis (s^-1 per layer) by the traced light above the base: the atoms each layer's span of log
    pressure holds, over its volume, the whole of a layer's O2 split; water by the same Lyman-alpha. Also the atoms
    the light makes by group of bands (per s per m^2 of lunar surface)."""
    below, lyman, by_group = traced_production(shield, treatment, activity, gaps, heat, protected_R)
    x_levels = np.log(lr.BASE_PA / col['p_pa'])
    cum = np.interp(np.maximum(x_levels, 0.0), traced.SHAPE_X, below)
    cum_l = np.interp(np.maximum(x_levels, 0.0), traced.SHAPE_X, lyman)
    r_layer = col['planet'].radius_m + col['layer_z_m']
    per_volume = (cum[1:] - cum[:-1]) * (MOON_RADIUS / r_layer) ** 2 / col['layer_dz_cm'] * 1e-4   # cm^-3 s^-1
    per_volume_l = (cum_l[1:] - cum_l[:-1]) * (MOON_RADIUS / r_layer) ** 2 / col['layer_dz_cm'] * 1e-4
    n_tot = th.layer_number_density(col)
    o2 = col['dry_fractions']['O2'] * (1 - col['layer_x_h2o']) * n_tot
    j_o2 = np.where(col['above_base'], per_volume / 2.0 / np.maximum(o2, 1e-300), 0.0)
    j_lyman = np.where(col['above_base'], per_volume_l / 2.0 / np.maximum(o2, 1e-300), 0.0)
    return j_o2, j_lyman * LYMAN_ALPHA_H2O_CM2 / escape.O2_LYMAN_ALPHA_CM2, by_group


def hot_escape(made, sigma_cm2, share, collision_cm2):
    """Atoms per s that leave hot, of `made` per s that a light makes in a column it crosses whole: in the unattenuated
    top a molecule of the absorber (`share` of the gas) breaks at sigma times twice the flux onto the column, and an
    atom made at collision depth s escapes on a straight path with chance E2(s) / 2."""
    return made * sigma_cm2 * share / (2.0 * collision_cm2)


def diffusion_parameter(species, t):
    """n D (cm^-1 s^-1) of a species in N2."""
    if species in DIFFUSION:
        a, s = DIFFUSION[species]
        return a * t ** s
    a, s = DIFFUSION['O']
    mu = lambda m: m * N2_AMU / (m + N2_AMU)
    scale = math.sqrt(mu(MASS_AMU['O']) / mu(MASS_AMU[species]))
    scale *= ((DIAMETER_A['O'] + N2_DIAMETER_A) / (DIAMETER_A[species] + N2_DIAMETER_A)) ** 2
    return a * t ** s * scale


def bernoulli(x):
    """x / (e^x - 1): the weight of exponential fitting (Scharfetter and Gummel) for a drift across a layer."""
    x = np.asarray(x, float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 1.0 - x[small] / 2.0 + x[small] ** 2 / 12.0     # its series, so B(x) / B(-x) stays exp(-x)
    with np.errstate(over='ignore'):
        out[~small] = x[~small] / np.expm1(x[~small])
    return out


class Transport:
    """Eddy and molecular diffusion of each species through spherical shells, with escape at the top."""

    def __init__(self, col, case, setup, escape_velocity):
        self.nl = col['layer_p_pa'].size
        r = (col['planet'].radius_m + col['z_m']) * 100.0                 # cm, levels
        self.r = r
        self.area = r ** 2
        self.volume = (r[1:] ** 3 - r[:-1] ** 3) / 3.0
        zc = col['layer_z_m'] * 100.0
        dz = zc[1:] - zc[:-1]
        p_lev, t_lev = col['p_pa'], col['t_k']
        n_lev = p_lev / (BOLTZMANN * t_lev) * 1e-6
        k_lev = case.mixing().profile(p_lev, col['tropopause_pa']) * setup.mixing_scale
        if setup.mixing_above_base == 'gravity_waves':
            waves = gw.eddy_diffusion(gw.lunar_column(col, MOON_GM, col['planet'].radius_m, case.air().molar_mass),
                                      col['tropopause_pa'], setup.wave_stretch,
                                      gw.calibration()[0] * setup.wave_factor_scale)
            k_lev = np.where(p_lev < lr.BASE_PA, waves, k_lev)
        molar = case.air().molar_mass * (1 - col['x_h2o']) + 0.018015 * col['x_h2o']        # kg/mol at levels
        g = MOON_GM * 1e6 / r ** 2                                                           # cm s^-2
        h_air = BOLTZMANN * 1e7 * t_lev / (molar * 1e3 / AVOGADRO * g)                      # cm
        self.m_layer = th.layer_number_density(col)
        inner = slice(1, -1)
        self.coefficients = {}
        self.diffusion = {}
        self.weights = {}
        for s in ch.SPECIES:
            d = diffusion_parameter(s, t_lev) / n_lev
            h_s = BOLTZMANN * 1e7 * t_lev / (MASS_AMU[s] * AMU_G * g)
            a = ((k_lev + d) * n_lev)[inner]
            c = (d * n_lev * (1.0 / h_s - 1.0 / h_air))[inner]
            pe = c * dz / a
            up, down = a / dz * bernoulli(pe), a / dz * bernoulli(-pe)
            self.weights[s] = (up, down)
            rows, cols, vals = [], [], []
            for k in range(self.nl - 1):
                area = self.area[k + 1]
                # Flux upward through level k+1: up f_k - down f_{k+1}, f = n / M.
                for row, sign in ((k, -1.0), (k + 1, 1.0)):
                    rows += [row, row]
                    cols += [k, k + 1]
                    vals += [sign * area * up[k] / (self.m_layer[k] * self.volume[row]),
                             -sign * area * down[k] / (self.m_layer[k + 1] * self.volume[row])]
            if s in escape_velocity:
                rows.append(self.nl - 1)
                cols.append(self.nl - 1)
                vals.append(-self.area[-1] * escape_velocity[s] / self.volume[-1])
            self.coefficients[s] = (np.array(rows), np.array(cols), np.array(vals))
            self.diffusion[s] = d
        self.eddy = k_lev

    def flux_through(self, n, species, level):
        """Upward flux (cm^-2 s^-1) of a species through a level (1 .. nl-1) between two layers."""
        up, down = self.weights[species]
        k = level - 1
        return float(up[k] * n[k] / self.m_layer[k] - down[k] * n[k + 1] / self.m_layer[k + 1])


def escape_velocities(col):
    """Jeans effusion velocity (cm/s) at the top level for the escaping species, O's raised by Earth's tide."""
    r_m = col['planet'].radius_m + col['z_m'][-1]
    t = float(col['t_k'][-1])
    out, lam = {}, {}
    for s in ESCAPING:
        m = MASS_AMU[s] * AMU_G * 1e-3
        lam[s] = MOON_GM * m / (BOLTZMANN * t * r_m)
        w = math.sqrt(BOLTZMANN * t / (2 * math.pi * m)) * (1 + lam[s]) * math.exp(-lam[s])
        if s == 'O':
            w *= tides.multiplier(lr.tidal_table(), r_m / MOON_RADIUS, lam[s])
        out[s] = w * 100.0
    return out, lam


def solve(col, case, setup, extra_o2, j_h2o, n0, max_steps=4000, tol=1e-7, verbose=False):
    """Steady state of the column's composition: backward-Euler steps of growing length with Newton iterations, as in
    chemistry.solve, with each species' own transport and escape, the photolysis above 202 nm refreshed from the
    composition, and the O2 and water photolysis added above it."""
    fixed = ch.fixed_densities(col)
    t = col['layer_t_k']
    nl, ns = col['layer_p_pa'].size, len(ch.SPECIES)
    velocity, lam = escape_velocities(col)
    geo = Transport(col, case, setup, velocity)
    boundary = case.boundary()
    fixed_surface = {ch.INDEX[s]: f for s, f in boundary.fixed}
    deposition = {ch.INDEX[s]: v for s, v in boundary.deposition}
    soluble = [ch.INDEX[s] for s in boundary.soluble]
    tropo = col['layer_p_pa'] >= col['tropopause_pa']
    water = np.sum(col['layer_mass_kg_m2'] * col['layer_x_h2o'] * 0.018015
                   / (0.018015 * col['layer_x_h2o'] + 0.02897 * (1 - col['layer_x_h2o'])))
    washout = (boundary.precipitation_m_per_year * 1000.0 / 3.156e7) / max(water, 1e-9)
    sun = ph.Sunlight(cl.load_shield(case.shield))
    h2o_source = j_h2o * fixed['H2O']
    n = np.maximum(n0.copy(), 1e-30 * fixed['M'])

    def rates(nn):
        dens = {s: nn[ch.INDEX[s]] for s in ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')}
        j = dict(sun.evaluate(col, dens, case.surface_albedo)['rates'])
        j[J['O2']] = np.asarray(j.get(J['O2'], np.zeros(nl)), float) + extra_o2
        return j

    j_rates = rates(n)

    def residual_and_matrix(nn, dt, old):
        dndt, jac = ch._rates_and_jacobian(nn, fixed, t, j_rates)
        for i, s in enumerate(ch.SPECIES):
            rows, cols, vals = geo.coefficients[s]
            dndt[i] += np.bincount(rows, weights=vals * nn[i][cols], minlength=nl)
            if i in deposition:
                dndt[i, 0] -= deposition[i] * nn[i, 0] * geo.area[0] / geo.volume[0]
            if i in soluble:
                dndt[i] -= np.where(tropo, washout, 0.0) * nn[i]
        dndt[ch.INDEX['H']] += h2o_source
        dndt[ch.INDEX['OH']] += h2o_source
        res = (nn - old) / dt - dndt
        rows_m, cols_m, vals_m = [], [], []
        base = np.arange(nl) * ns
        for i in range(ns):
            for k2 in range(ns):
                v = -jac[:, i, k2]
                if i == k2:
                    v = v + 1.0 / dt
                    if i in deposition:
                        v = v.copy()
                        v[0] += deposition[i] * geo.area[0] / geo.volume[0]
                    if i in soluble:
                        v = v + np.where(tropo, washout, 0.0)
                nz = v != 0
                rows_m.append(base[nz] + i)
                cols_m.append(base[nz] + k2)
                vals_m.append(v[nz])
            tr, tc, tv = geo.coefficients[ch.SPECIES[i]]
            rows_m.append(tr * ns + i)
            cols_m.append(tc * ns + i)
            vals_m.append(-tv)
        for i, f in fixed_surface.items():
            res[i, 0] = nn[i, 0] - f * fixed['M'][0]
        rows_m, cols_m, vals_m = np.concatenate(rows_m), np.concatenate(cols_m), np.concatenate(vals_m)
        keep = np.ones(rows_m.size, dtype=bool)
        for i in fixed_surface:
            keep &= rows_m != i
        rows_m = np.concatenate([rows_m[keep], list(fixed_surface)])
        cols_m = np.concatenate([cols_m[keep], list(fixed_surface)])
        vals_m = np.concatenate([vals_m[keep], np.ones(len(fixed_surface))])
        mat = coo_matrix((vals_m, (rows_m, cols_m)), shape=(nl * ns, nl * ns)).tocsc()
        return res.T.ravel(), mat

    dt, steps, floor = 1e-4, 0, 1e-30 * fixed['M']
    while steps < max_steps:
        steps += 1
        old, trial, ok = n.copy(), n.copy(), False
        for _ in range(8):
            res, mat = residual_and_matrix(trial, dt, old)
            try:
                delta = splu(mat).solve(-res).reshape(nl, ns).T
            except RuntimeError:
                break
            trial = trial + delta
            if not np.all(np.isfinite(trial)):
                break
            trial = np.maximum(trial, floor)
            change = np.max(np.abs(delta) / np.maximum(np.abs(trial), 1e-18 * fixed['M'] + 1.0))
            if change < 1e-3:
                ok = True
                break
        if not ok:
            dt /= 4.0
            if dt < 1e-10:
                raise RuntimeError('Chemistry step size collapsed')
            continue
        rel = np.max(np.abs(trial - old) / np.maximum(trial, 1e-10 * fixed['M']))
        n = trial
        if dt < 1e9 or steps % 5 == 0:
            j_rates = rates(n)
        if verbose and steps % 20 == 0:
            print(f'  step {steps} dt {dt:.2e} rel {rel:.2e}', flush=True)
        if dt > 1e15 and rel < tol:
            j_rates = rates(n)
            break
        dt = min(dt * (2.0 if rel < 0.1 else 1.2), 1e17)
    return n, dict(geo=geo, velocity=velocity, lam=lam, j_rates=j_rates, steps=steps, converged=steps < max_steps,
                   fixed=fixed)


def exosphere_oxygen(profile, n_o_exo, activity, protected_R, outer_R=ex.OUTER_R, points=6000):
    """The sunlit exosphere's oxygen atoms outside the shadow (kg/s): photoions made from the ballistic atoms, and the
    ions the solar wind makes from them by charge exchange on straight paths at full exposure."""
    _, t_c, r_c = ab.exobase_state(profile)
    n_c = n_o_exo * 1e6
    kg = MASS_AMU['O'] * AMU_G * 1e-3
    r = np.geomspace(r_c, outer_R, points)
    bal, esc = ab.exosphere_density(r, n_c, t_c, r_c, mass_kg=kg)
    lit = np.sqrt(np.clip(1.0 - (protected_R / r) ** 2, 0.0, None))
    ions = ab.ionization_rate(activity) * np.trapezoid(lit * kg * bal * 4 * math.pi * r ** 2 * MOON_RADIUS ** 3, r)
    b = np.unique(np.r_[np.linspace(0.0, r_c, 250, endpoint=False), r_c + (outer_R - r_c) * np.linspace(0.0, 1.0, 250) ** 2])
    z0 = np.sqrt(np.clip(r_c ** 2 - b ** 2, 0.0, None))
    z1 = np.sqrt(np.clip(outer_R ** 2 - b ** 2, 0.0, None))
    u = np.linspace(0.0, 1.0, 600) ** 2
    z = z0[:, None] + (z1 - z0)[:, None] * u[None, :]
    rr = np.clip(np.sqrt(b[:, None] ** 2 + z ** 2), r_c, outer_R)
    sides = np.where(b < r_c, 1.0, 2.0)
    ln_bal = np.log(np.maximum(bal, 1e-300))
    ln_all = np.log(np.maximum(bal + esc, 1e-300))
    col_bal = sides * np.trapezoid(np.exp(np.interp(rr, r, ln_bal)), z, axis=1) * MOON_RADIUS
    col_all = sides * np.trapezoid(np.exp(np.interp(rr, r, ln_all)), z, axis=1) * MOON_RADIUS
    flux = lr.SOLAR_WIND['density_m3'] * lr.SOLAR_WIND['speed_m_s']
    exchanged = -np.expm1(-ex.CHARGE_EXCHANGE_M2 * col_all)
    per_path = np.where(col_all > 0, exchanged * col_bal * kg / np.maximum(col_all, 1e-300), 0.0)
    exchange = float(flux * np.trapezoid(per_path * 2 * math.pi * b, b) * MOON_RADIUS ** 2)
    return float(ions), exchange


def state(shield, treatment, activity, transmission, protected_R=MAX_RADIUS_R):
    """The loss response's traced state at a transmission: its heat, the thermal column's summary and profile there,
    and the column's configuration."""
    grid, parts = traced.parts(traced.entries(f'{shield}_{treatment}'), traced.FILM, activity,
                               lr.glow_heat(shield, treatment, activity), protected_R)
    q = traced.first_state(grid, parts, transmission)
    if q is None:
        return None
    cfg = lr.state_config(shield, treatment, activity, protected_R, transmission, q)
    summary, profile, _ = solve_column(q, cfg)
    return dict(heat_w_m2=q, summary=summary, profile=profile, config=cfg)


def run(shield, treatment, activity, transmission, setup=Setup(), protected_R=MAX_RADIUS_R, verbose=False):
    """Oxygen atoms and hydrogen through the column at one state: their making, their profiles, and what escapes."""
    st = state(shield, treatment, activity, transmission, protected_R)
    if st is None:
        return dict(status='runaway')
    case = middle_case(shield, treatment)
    layers, summary = middle_layers(case)
    col = column(case, layers, summary, st['profile'], setup)
    glow_o2, glow_h2o = glow_rates(col, activity)
    trace_o2, trace_h2o, by_group = traced_rates(col, shield, treatment, activity, transmission, st['heat_w_m2'],
                                                 protected_R)
    n, info = solve(col, case, setup, glow_o2 + trace_o2, glow_h2o + trace_h2o, initial_density(col, layers),
                    verbose=verbose)
    geo, fixed = info['geo'], info['fixed']
    i_o, i_o3, i_h, i_h2 = (ch.INDEX[s] for s in ('O', 'O3', 'H', 'H2'))
    shell = 4 * math.pi * geo.r ** 2                                    # cm^2 at levels
    volume = 4 * math.pi * (geo.r[1:] ** 3 - geo.r[:-1] ** 3) / 3.0     # cm^3 per layer
    above = col['above_base']
    kg = lambda species: MASS_AMU[species] * AMU_G * 1e-3
    o2 = fixed['O2']
    made = {name: float(np.sum(2 * j * o2 * volume * above)) for name, j in (('glow', glow_o2), ('traced', trace_o2))}
    made['glow_below_base'] = float(np.sum(2 * glow_o2 * o2 * volume * ~above))
    top = col['layer_p_pa'].size - 1
    escape_rate = {s: float(info['velocity'][s] * n[ch.INDEX[s], top] * shell[-1]) for s in ESCAPING}
    base_level = int(np.argmax(col['p_pa'] < lr.BASE_PA * 1.0001))
    down = -(geo.flux_through(n[i_o], 'O', base_level) + geo.flux_through(n[i_o3], 'O3', base_level)) * shell[base_level]
    homopause = np.nonzero(geo.diffusion['O'] > geo.eddy)[0]
    m_tot = fixed['M']
    x_o = n[i_o] / m_tot
    ions, exchange = exosphere_oxygen(st['profile'], float(n[i_o, top]), activity, protected_R)
    exo_cfg = ex.Exosphere(st['profile'], activity)
    fates = exo_cfg.fates(protected_R, ex.SEPTEMBER_MOMENT, 'central')
    kept = (fates['drained'] + fates['recombined_staying']) / max(sum(fates[k] for k in (
        'outside', 'open', 'convected', 'recombined_escaping', 'recombined_staying', 'drained')), 1e-300)
    wind = ex.no_magnetosphere(ions, exchange, 'central')
    jeans = escape_rate['O'] * kg('O')
    sigma_c = st['config'].collision_cross_section_m2 * 1e4
    o2_share = float(o2[top] / m_tot[top])
    shares = [1.0 if a == 'air' else o2_share for a in HOT_ABSORBER]
    hot_traced = sum(hot_escape(m * 4 * math.pi * MOON_RADIUS ** 2, sigma, x, sigma_c)
                     for m, sigma, x in zip(by_group, HOT_CROSS_SECTION_CM2, shares)) * kg('O')
    hot_glow = hot_escape(made['glow'] + made['glow_below_base'], escape.O2_LYMAN_ALPHA_CM2, o2_share, sigma_c) * kg('O')
    hot = hot_glow + hot_traced
    out = dict(
        status='solved' if info['converged'] else 'not_converged', steps=info['steps'],
        heat_w_m2=st['heat_w_m2'], exobase_radius_R=st['summary']['exobase_radius_R'],
        exobase_temperature_k=st['summary']['exobase_temperature_k'], column_top_radius_R=float(geo.r[-1] / 1e2 / MOON_RADIUS),
        jeans_lambda=info['lam'], escape_velocity_cm_s=info['velocity'],
        oxygen_made_kg_s={k: v * kg('O') for k, v in made.items()},
        oxygen_carried_down_kg_s=float(down * kg('O')),
        oxygen_jeans_escape_kg_s=jeans,
        oxygen_exosphere_kg_s=dict(ions_made=ions, charge_exchange_upper_bound=exchange,
                                   no_magnetosphere=wind['plasma_kg_s'], september_magnetosphere=ions * (1 - kept)),
        oxygen_hot_escape_kg_s=dict(glow=hot_glow, traced=hot_traced, upper=HOT_COLLISIONS * hot),
        oxygen_loss_kg_s=dict(no_magnetosphere=jeans + wind['plasma_kg_s'] + hot,
                              september_magnetosphere=jeans + ions * (1 - kept) + hot),
        oxygen_escaping_share_of_made=(jeans + wind['plasma_kg_s'] + hot) / max(made['glow'] + made['traced'], 1e-300) / kg('O'),
        hydrogen_escape_kg_s=float((escape_rate['H'] + 2 * escape_rate['H2']) * kg('H')),
        oxygen_mixing_ratio=dict(base=float(np.interp(math.log(lr.BASE_PA), np.log(col['layer_p_pa'][::-1]), x_o[::-1])),
                                 homopause=float(x_o[homopause[0]]) if homopause.size else None,
                                 exobase=float(x_o[top])),
        oxygen_share_at_exobase=float(n[i_o, top] / (m_tot[top] + n[i_o, top])),
        homopause_pa=float(col['p_pa'][homopause[0]]) if homopause.size else None,
        ozone_column_du=float(np.sum(n[i_o3] * col['layer_dz_cm']) / ch.DU),
        profile=dict(p_levels_pa=col['p_pa'].tolist(), eddy_cm2_s=geo.eddy.tolist(),
                     p_pa=col['layer_p_pa'].tolist(), o=x_o.tolist(), o3=(n[i_o3] / m_tot).tolist(),
                     h=(n[i_h] / m_tot).tolist(), h2=(n[i_h2] / m_tot).tolist(),
                     oh=(n[ch.INDEX['OH']] / m_tot).tolist(), ho2=(n[ch.INDEX['HO2']] / m_tot).tolist()))
    return out


ACTIVITIES = ('quiet', 'solar_maximum', 'solar_maximum_cycle_19', 'solar_maximum_stress', 'cycles_23_24_mean')
# Titan's measured homopause eddy diffusion, 2e7-1e8 cm^2/s since Cassini (argon: Yelle et al. 2008, JGR 113,
# E10003; Bell et al. 2014, doi:10.1002/2014JA019781), carried to the Moon as Lindzen's mixing scales with the same waves,
# T^(1/2) g^-2: about 0.8 for the Moon's warmer, slightly stronger-pulled thermosphere, so 1.5e7-8e7 against the
# middle atmosphere's Moon-scaled 3.6e8. The analogue takes the middle of that range, a tenth of the Moon-scaled
# mixing; Titan's troposphere is driven by about a hundredth of the Moon's sunlight, so its waves may be the weaker.
TITAN_ANALOGUE = 0.1
SETUPS = {'moon_mixing': Setup(), 'earth_mixing': Setup(mixing_scale=1.0 / middle_run.MOON_KZZ),
          'gravity_waves': Setup(mixing_above_base='gravity_waves'),
          'gravity_waves_weak': Setup(mixing_above_base='gravity_waves', wave_stretch=1.0, wave_factor_scale=1.0 / 3.0),
          'titan_analogue': Setup(mixing_scale=TITAN_ANALOGUE)}
YEARLY = (('titania_stack', 'all_heats'),)          # the case whose oxygen matters, followed year by year
WORKERS = 3                                          # single-threaded processes; the machine is shared
FILES = ('atmosphere/loss_response/oxygen.py', 'atmosphere/loss_response/gravity_waves.py',
         'atmosphere/loss_response/model.py', 'atmosphere/loss_response/traced.py',
         'atmosphere/loss_response/absorption.py', 'atmosphere/loss_response/exosphere.py',
         'atmosphere/loss_response/tides.py', 'atmosphere/thermal_column.py', 'atmosphere/middle_atmosphere/chemistry.py',
         'atmosphere/middle_atmosphere/photolysis.py', 'atmosphere/middle_atmosphere/escape.py',
         'atmosphere/radiative_convective/thermodynamics.py')
PRODUCTS = ('atmosphere/middle_atmosphere/results/limb_heat.json', 'atmosphere/loss_response/results/tidal_escape.json',
            'atmosphere/loss_response/results/solar_cycle.json', 'protection/transmission/results/uv_transmission.json')
# The photolysis reads the middle atmosphere's and the radiative model's external inputs, pinned by their manifests.
MANIFESTS = ('atmosphere/middle_atmosphere/inputs.json', 'atmosphere/radiative_convective/inputs.json')


def activity_named(name):
    if name == 'cycles_23_24_mean':
        return cycle.product()['spans']['cycles_23_24']['mean_spectrum']
    years = {a['label']: a for a in cycle.product()['spans']['cycles_23_24']['years']}
    return years.get(name, name)


def task(job):
    """One run for the pool: (setup name, shield, treatment, activity name, transmission, keep the profile)."""
    key, shield, treatment, name, transmission, keep = job
    result = run(shield, treatment, activity_named(name), transmission, SETUPS[key])
    if not keep:
        result.pop('profile', None)
    return dict(setup=key, shield=shield, treatment=treatment, activity=name, transmission=transmission,
                protected_radius_R=MAX_RADIUS_R, result=result)


def standard_transmission():
    product = json.loads((ROOT / PRODUCTS[3]).read_text())
    return product['levels'][LEVEL]['transmission_by_hole_factor']['1']


def summarise(runs):
    held = [r for r in runs if r['result'].get('status') == 'solved']
    span = lambda get: [min(map(get, held)), max(map(get, held))] if held else None
    out = {}
    for key in SETUPS:
        for shield in lr.SHIELDS:
            rows = [r['result'] for r in held if r['setup'] == key and r['shield'] == shield]
            if not rows:
                continue
            get = lambda f: [min(map(f, rows)), max(map(f, rows))]
            out.setdefault(key, {})[shield] = dict(
                runs=len(rows),
                oxygen_made_above_base_kg_s=get(lambda r: r['oxygen_made_kg_s']['glow'] + r['oxygen_made_kg_s']['traced']),
                oxygen_loss_no_magnetosphere_kg_s=get(lambda r: r['oxygen_loss_kg_s']['no_magnetosphere']),
                oxygen_loss_september_kg_s=get(lambda r: r['oxygen_loss_kg_s']['september_magnetosphere']),
                oxygen_jeans_escape_kg_s=get(lambda r: r['oxygen_jeans_escape_kg_s']),
                oxygen_hot_escape_kg_s=get(lambda r: r['oxygen_hot_escape_kg_s']['glow'] + r['oxygen_hot_escape_kg_s']['traced']),
                oxygen_hot_escape_upper_kg_s=get(lambda r: r['oxygen_hot_escape_kg_s']['upper']),
                oxygen_share_at_exobase=get(lambda r: r['oxygen_share_at_exobase']),
                homopause_pa=get(lambda r: r['homopause_pa']),
                hydrogen_escape_kg_s=get(lambda r: r['hydrogen_escape_kg_s']),
                ozone_column_du=get(lambda r: r['ozone_column_du']))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(argv)
    transmission = standard_transmission()
    jobs = [(key, shield, treatment, name, transmission, key == 'moon_mixing')
            for key in SETUPS for shield in lr.SHIELDS for treatment in lr.TREATMENTS for name in ACTIVITIES
            if key == 'moon_mixing' or shield == 'titania_stack']
    years = [a['label'] for a in cycle.product()['spans']['cycles_23_24']['years']]
    yearly_jobs = [(key, shield, treatment, y, transmission, False)
                   for key in SETUPS for shield, treatment in YEARLY for y in years]
    with Pool(WORKERS, initializer=os.nice, initargs=(10,)) as pool:
        runs = pool.map(task, jobs, chunksize=1)
        yearly = pool.map(task, yearly_jobs, chunksize=1)
    for r in runs:
        loss = r['result'].get('oxygen_loss_kg_s', {})
        print(f"{r['setup']:12s} {r['shield']:13s} {r['treatment']:11s} {r['activity']:18s} {r['result']['status']:8s} "
              f"O loss {loss.get('no_magnetosphere', float('nan')):.2e} / {loss.get('september_magnetosphere', float('nan')):.2e} kg/s, "
              f"H {r['result'].get('hydrogen_escape_kg_s', float('nan')):.2e}", flush=True)
    cycle_mean = {}
    for key in SETUPS:
        for shield, treatment in YEARLY:
            rows = [r['result'] for r in yearly if r['setup'] == key and r['shield'] == shield
                    and r['treatment'] == treatment]
            solved = [r for r in rows if r['status'] == 'solved']
            mean = lambda f: float(np.mean([f(r) for r in solved])) if len(solved) == len(rows) else None
            cycle_mean.setdefault(key, {})[f'{shield}_{treatment}'] = dict(
                years=len(rows), years_without_a_state=len(rows) - len(solved),
                oxygen_loss_no_magnetosphere_kg_s=mean(lambda r: r['oxygen_loss_kg_s']['no_magnetosphere']),
                oxygen_loss_september_kg_s=mean(lambda r: r['oxygen_loss_kg_s']['september_magnetosphere']),
                oxygen_hot_escape_kg_s=mean(lambda r: r['oxygen_hot_escape_kg_s']['glow'] + r['oxygen_hot_escape_kg_s']['traced']),
                largest_yearly_loss_no_magnetosphere_kg_s=max((r['oxygen_loss_kg_s']['no_magnetosphere'] for r in solved),
                                                              default=None),
                oxygen_share_at_exobase=[min(r['oxygen_share_at_exobase'] for r in solved),
                                         max(r['oxygen_share_at_exobase'] for r in solved)] if solved else None)
            print(key, shield, treatment, cycle_mean[key][f'{shield}_{treatment}'], flush=True)
    files = [ROOT / f for f in FILES]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={f: digest(ROOT / f) for f in FILES}, constants=constants_used(files),
                      products={f: digest(ROOT / f) for f in PRODUCTS},
                      inputs=dict({f'atmosphere/middle_atmosphere/results/profiles/{c}.csv': digest(MIDDLE / 'profiles' / f'{c}.csv')
                                   for c in sorted({middle_case(r['shield'], r['treatment']).name for r in runs})},
                                  **{m: digest(ROOT / m) for m in MANIFESTS})),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('Loss rates are kg/s of the gas named. oxygen_loss_kg_s is the oxygen atoms\' Jeans escape with '
                      'Earth\'s tide and the hot atoms made near the exobase (oxygen_hot_escape_kg_s, glow and traced '
                      'light; its upper value lets them keep their escape energy through two collisions) plus, with no '
                      'magnetosphere, the central plasma loss of their photoions outside the '
                      'protected radius and of the ions the solar wind makes from them by charge exchange (the exosphere '
                      'step\'s no_magnetosphere), or with the September magnets the photoions not kept, in the share the '
                      'exosphere step finds for the molecules. hydrogen_escape_kg_s is H and H2 escaping, counted as '
                      'hydrogen. oxygen_made_kg_s counts the atoms the glow and the traced light make above the 0.3 Pa base '
                      '(and the glow below it); oxygen_carried_down_kg_s the odd oxygen going down through the base. Runs '
                      'are at the swarm\'s standard level and the ring fleet\'s 4 lunar radii, the column\'s heat placed '
                      'where the tracing puts it; gravity_waves takes the breaking waves\' mixing above the base '
                      '(gravity_waves.py; gravity_waves_weak with Earth\'s wavelengths and a third of the calibration), '
                      'titan_analogue a tenth of the Moon-scaled mixing, as Titan\'s measured mixing would give, and '
                      'earth_mixing divides the middle atmosphere\'s eddy mixing by its 36-fold '
                      'scaling for the Moon, a bound on weak mixing. cycles_23_24_mean is the state of the two cycles\' '
                      'mean spectrum; cycle_mean averages the yearly states of the warmest case over 1997-2019, which '
                      'leans high.'),
        assumptions=dict(lyman_alpha_h2o_cm2=LYMAN_ALPHA_H2O_CM2, ev_per_ion_pair=EV_PER_ION_PAIR,
                         hot_cross_section_cm2=HOT_CROSS_SECTION_CM2, hot_absorber=HOT_ABSORBER,
                         hot_collisions=HOT_COLLISIONS,
                         layers_per_efold=LAYERS_PER_EFOLD, diffusion_parameters=DIFFUSION, diameters_angstrom=DIAMETER_A,
                         setups={k: dataclasses.asdict(v) for k, v in SETUPS.items()}),
        summary=summarise(runs), cycle_mean=cycle_mean, runs=runs, yearly=yearly)
    args.out.write_text(json.dumps(product, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
