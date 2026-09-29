"""Sunlight at the ground by wavelength under a clear sky: the light plants, eyes and skin receive.

    python -m illumination.surface_light.model                      # writes results/surface_light.json
    python -m illumination.surface_light.model --atlas-dir DIR      # also compares with the sky atlases in DIR

Direct and diffuse downward irradiance at the surface from 202 to 1000 nm in 1-nm bins, at
Sun heights from 1 to 90 degrees, for three columns:

- the design Moon: 1.2 atm (121,590 Pa of dry air) and a 294.9 K surface, the chosen climate's
  global mean (climate/gcm run A28_dim5_moon, years 20-29), behind the titania film passing 5% less sunlight
  at every wavelength (research/decisions.md);
- the same air at 1.0 atm;
- the middle atmosphere's Earth control (1 atm, 288 K, its solved ozone) in unfiltered sunlight.

The design Moon is also given in unfiltered sunlight, which separates what the air does from
what the shield does and matches the sky atlas's configuration.

The radiation is the atmosphere domain's own. Below 500 nm it is the middle atmosphere's
ultraviolet-visible calculation (photolysis.py: 1-nm bins; Rayleigh scattering; O2, O3, NO2
and the other photolysed absorbers; O2 collision-induced bands). From 500 nm it is that
model's near-infrared step: the radiative-convective line-by-line absorption (H2O and O2
lines, the MT_CKD water continuum, collision-induced bands) with ozone's Chappuis band, on
every second level. Both use the delta-two-stream solver with a pseudo-spherical direct beam
(radiative_convective/shortwave.py). Each column is the moist adiabat from its surface
temperature with Manabe-Wetherald humidity up to the tropopause, and above it the middle
atmosphere's solved temperatures; trace gases are its solved mixing ratios. For the Moon at
294.9 K both are interpolated between its stored cases at 288 and 298 K.

The ground is Lambertian with albedo A. The downward light at the ground is the light over a
black ground divided by (1 - A R), where R is the atmosphere's reflectance for diffuse light
arriving from below; the direct beam does not depend on A, and the light returned to space
gains A F T, with F that downward light and T the atmosphere's transmittance for light from
below. This holds exactly in the two-stream adding method (the tests compare it with direct
solutions), so any ground, including one whose albedo varies with wavelength, follows from
the black-ground solution and R.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

from shared.constants import AVOGADRO, PLANCK, SPEED_OF_LIGHT, SYNODIC_MONTH_DAYS
from atmosphere.radiative_convective import climate as cl, optics, shortwave as sw, thermodynamics as th
from atmosphere.radiative_convective import fetch_inputs as rc_inputs
from atmosphere.middle_atmosphere import equilibrium as eq, photolysis as ph
from illumination.sky.colour_matching import cmf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.illumination.surface-light/1'
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
MIDDLE_SCHEMA = 'terluna.atmosphere.middle-atmosphere/1'

JOIN_NM, TOP_NM = 500.0, 1000.0          # the middle atmosphere's split between 1-nm and line-by-line sunlight
EDGES_NM = np.arange(ph.EDGES_NM[0], TOP_NM + 0.5, 1.0)
CENTRES_NM = 0.5 * (EDGES_NM[1:] + EDGES_NM[:-1])
UV_BINS = ph.CENTRES_NM < JOIN_NM
N_UV = int(UV_BINS.sum())
N_LBL = int(TOP_NM - JOIN_NM)
# The middle atmosphere's line-by-line settings (optics.SHORTWAVE) over 500-1000 nm. CO2 has no lines in
# this range, and O2 is computed line by line in each layer instead of from the cached node table.
LBL = optics.OpticsConfig(nu_min=1e7 / TOP_NM, nu_max=1e7 / JOIN_NM, step=optics.SHORTWAVE.step,
                          molecules=('H2O', 'O2'), coarse_factor=optics.SHORTWAVE.coarse_factor, tabulated=())
STRIDE = eq.Case.sw_stride
CHUNK = 40000
MICROMOL = 1e6 / AVOGADRO

SUN_DEG = (90.0, 75.0, 60.0, 45.0, 30.0, 20.0, 15.0, 10.0, 7.0, 5.0, 3.0, 2.0, 1.0)
ALBEDOS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
LATITUDES = (0.0, 30.0, 60.0)
# Photon bands (nm). red_660 and far_red_730 are the 10-nm windows of the red:far-red ratio (Smith 1982).
BANDS_NM = dict(uv_b=(280, 315), uv_a=(315, 400), blue=(400, 500), green=(500, 600), red=(600, 700),
                par=(400, 700), far_red=(700, 750), near_infrared=(750, 1000),
                red_660=(655, 665), far_red_730=(725, 735))
TRACE = ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')
DOBSON_CM2 = 2.6867e16

DESIGN_TS = 294.9                 # the chosen climate's global mean (A28_dim5's 294 K before 2026-09-29)
W298 = (DESIGN_TS - 288.0) / 10.0


@dataclass(frozen=True)
class Column:
    name: str
    planet: th.Planet
    dry_pressure_pa: float
    surface_temperature_k: float
    stored: tuple                        # (middle-atmosphere case, weight): its solved upper air and trace gases
    stratospheric_h2o: float | None = None

    def case(self):
        return eq.Case(self.name, self.planet, self.dry_pressure_pa, surface_temperature_k=self.surface_temperature_k,
                       stratospheric_h2o=self.stratospheric_h2o)


COLUMNS = {c.name: c for c in (
    Column('moon_1.2atm', th.MOON, 121590.0, DESIGN_TS,
           (('moon_1.2atm_titania_stack', 1 - W298), ('moon_1.2atm_titania_stack_ts298', W298))),
    Column('moon_1.0atm', th.MOON, 101325.0, DESIGN_TS,
           (('moon_1.0atm_titania_stack', 1 - W298), ('moon_1.0atm_titania_stack_ts298', W298))),
    Column('earth_control', th.EARTH, 101325.0, 288.0, (('earth_control', 1.0),), stratospheric_h2o=5e-6),
)}
# Sunlight at the top of the air: a shield of the protection product and the share of its light passed.
SUNLIGHTS = dict(design=('titania_stack', 0.95), unfiltered=('none', 1.0))
CASES = (('moon_1.2atm', 'design'), ('moon_1.2atm', 'unfiltered'), ('moon_1.0atm', 'design'),
         ('earth_control', 'unfiltered'))


def stored_case(name):
    """A middle-atmosphere case: its profile table, its summary and the table's hash."""
    raw = (MIDDLE / 'profiles' / f'{name}.csv').read_bytes()
    rows = list(csv.DictReader(raw.decode().splitlines()))
    table = {k: np.array([float(r[k]) if r[k] else np.nan for r in rows]) for k in rows[0]}
    summary = json.loads((MIDDLE / 'cases' / f'{name}.json').read_text())
    return table, summary, hashlib.sha256(raw).hexdigest()[:16]


def build(column: Column):
    """Levels (surface first), dry air, trace-gas densities (cm^-3) and the stored cases used, by hash."""
    head = json.loads((MIDDLE / 'middle_atmosphere.json').read_text())
    if head['schema'] != MIDDLE_SCHEMA:
        raise ValueError('Unexpected middle-atmosphere product schema')
    case = column.case()
    prof = eq.Profile(case)
    loaded = [(name, w, *stored_case(name)) for name, w in column.stored]

    def mix(key, log_p):
        return sum(w * np.interp(-log_p, -np.log(table['p_pa']), table[key]) for _, w, table, _, _ in loaded)
    lp = np.log(prof.p)
    trop_pa = math.exp(sum(w * math.log(summary['tropopause_pa']) for _, w, _, summary, _ in loaded))
    trop = int(np.argmin(np.abs(lp - math.log(trop_pa))))
    t = np.where(np.arange(lp.size) <= trop, prof.t_ad, mix('t_k', lp))
    col = prof.column(t, trop)
    n_tot = th.layer_number_density(col)
    densities = {s: np.maximum(mix(s, np.log(col['layer_p_pa'])), 0.0) * n_tot for s in TRACE}
    return col, prof.air, densities, {name: sha for name, _, _, _, sha in loaded}


def column_facts(col, air, densities):
    x = col['layer_x_h2o']
    m_h2o = th.MOLAR_MASS['H2O']
    water = float(np.sum(col['layer_mass_kg_m2'] * x * m_h2o / (m_h2o * x + air.molar_mass * (1 - x))))
    rayleigh = optics.rayleigh_optical_depth(col, np.array([1e7 / 550.0]))
    return dict(surface_pressure_pa=float(col['p_pa'][0]), surface_temperature_k=float(col['t_k'][0]),
                air_column_kg_m2=float(col['layer_mass_kg_m2'].sum()), water_column_kg_m2=water,
                ozone_du=float(np.sum(densities['O3'] * col['layer_dz_cm']) / DOBSON_CM2),
                rayleigh_optical_depth_550nm=float(rayleigh.sum()),
                tropopause_pa=float(col['tropopause_pa']))


def from_below(tau_abs, tau_sca):
    """Reflectance and transmittance of the whole atmosphere for diffuse light arriving from below.

    Layers top-down; each layer is added beneath the stack above it with the same delta-scaled
    two-stream layer solutions that shortwave.column_fluxes uses (their diffuse terms do not
    depend on the Sun's direction).
    """
    tau = tau_abs + tau_sca
    omega = np.where(tau > 0, tau_sca / np.maximum(tau, 1e-300), 0.0)
    r, t, *_ = sw.layer_properties(tau, omega, np.zeros_like(tau), np.ones_like(tau))
    refl, trans = np.zeros(tau.shape[1]), np.ones(tau.shape[1])
    for k in range(tau.shape[0]):
        d = 1 - r[k] * refl
        trans = t[k] * trans / d
        refl = r[k] + t[k] ** 2 * refl / d
    return refl, trans


def with_ground(direct, diffuse, up, refl, trans, albedo):
    """Surface direct and diffuse and the light returned to space over a Lambertian ground, from black-ground fluxes."""
    total = (direct + diffuse) / (1 - albedo * refl)
    return direct, total - direct, up + albedo * total * trans


class Response:
    """A column's black-ground fluxes per unit sunlight normal to the beam, on its two spectral grids."""

    def __init__(self, col, air, densities, processes=2):
        tau_abs, tau_sca, _ = ph.Sunlight().optics(col, densities)
        self.uv = (tau_abs[::-1][:, UV_BINS], tau_sca[::-1][:, UV_BINS])
        self.uv_radius = (col['planet'].radius_m + col['z_m'])[::-1]
        idx = np.arange(0, col['p_pa'].size, STRIDE)
        if idx[-1] != col['p_pa'].size - 1:
            idx = np.append(idx, col['p_pa'].size - 1)
        coarse = th.column_from_levels(col['planet'], air, col['p_pa'][idx], col['t_k'][idx], col['x_h2o'][idx])
        nu, tau = optics.optical_depth(coarse, replace(LBL, processes=processes))
        # Ozone's Chappuis and Wulf bands, as in the middle atmosphere's near-infrared step.
        o3 = np.add.reduceat(densities['O3'] * col['layer_dz_cm'], idx[:-1])
        sigma = np.interp(1e7 / nu, ph.CENTRES_NM, ph.CrossSections().o3_room, left=0.0, right=0.0)
        tau = tau + o3[:, None] * sigma[None, :]
        self.nu = nu
        self.lbl = (tau[::-1], optics.rayleigh_optical_depth(coarse, nu)[::-1])
        self.lbl_radius = (coarse['planet'].radius_m + coarse['z_m'])[::-1]
        self.lbl_bin = np.clip(np.floor(1e7 / nu).astype(int) - int(JOIN_NM), 0, N_LBL - 1)
        self.below = dict(uv=from_below(*self.uv), lbl=self._chunked(lambda sl: from_below(*(a[:, sl] for a in self.lbl))))

    def _chunked(self, func):
        parts = [func(slice(a, a + CHUNK)) for a in range(0, self.nu.size, CHUNK)]
        return tuple(np.concatenate(p) for p in zip(*parts))

    def fluxes(self, mu0):
        """Black-ground surface direct, surface diffuse and top upward fluxes (fractions of the normal-incidence flux)."""
        def solve(tau_abs, tau_sca, radius):
            direct, diffuse, up = sw.column_fluxes(tau_abs, tau_sca, np.zeros_like(tau_sca), 0.0, mu0, radius)
            return direct[-1], diffuse[-1], up[0]
        uv = solve(*self.uv, self.uv_radius)
        lbl = self._chunked(lambda sl: solve(self.lbl[0][:, sl], self.lbl[1][:, sl], self.lbl_radius))
        return dict(uv=uv, lbl=lbl)


def sunlight(name):
    """Shield transmission function for a named sunlight, and the protection product's hash."""
    shield, share = SUNLIGHTS[name]
    base = cl.load_shield(shield)
    if base is None and share == 1.0:
        return None, ''

    def transmission(lam_nm):
        return share * (base(lam_nm) if base is not None else np.ones_like(np.asarray(lam_nm, dtype=float)))
    return transmission, getattr(base, 'product_sha256', '')


def top_of_air(nu, shield):
    """Sunlight normal to the beam: energy (W/m^2) and photons (umol m^-2 s^-1) per 1-nm bin below 500 nm and
    per line-by-line grid point from 500 nm."""
    energy, count = ph.solar_bins(shield)
    uv = (energy[UV_BINS], count[UV_BINS] * 1e4 * MICROMOL)
    spectrum, _ = sw.solar_spectrum(nu, shield=shield)
    weights = np.full(nu.size, nu[1] - nu[0])
    weights[0] *= 0.5
    weights[-1] *= 0.5
    e = spectrum * weights
    return dict(uv=uv, lbl=(e, e / (PLANCK * SPEED_OF_LIGHT * nu * 100.0) * MICROMOL))


def spectra(resp: Response, light, mu0, albedos=ALBEDOS):
    """1-nm spectra at the ground and at the top of the air for each albedo.

    Returns arrays (albedos, bins): direct and diffuse energy (W/m^2 per nm) and photons
    (umol m^-2 s^-1 per nm) on a horizontal surface, photons returned to space and incident.
    """
    out = {k: np.zeros((len(albedos), CENTRES_NM.size)) for k in
           ('direct_w', 'diffuse_w', 'direct_umol', 'diffuse_umol', 'up_umol', 'incident_umol')}
    fl = resp.fluxes(mu0)
    for grid, (lo, hi) in (('uv', (0, N_UV)), ('lbl', (N_UV, N_UV + N_LBL))):
        e, p = light[grid]
        refl, trans = resp.below[grid]
        for i, a in enumerate(albedos):
            d, f, u = with_ground(*fl[grid], refl, trans, a)
            for key, values in (('direct_w', d * e), ('diffuse_w', f * e), ('direct_umol', d * p),
                                ('diffuse_umol', f * p), ('up_umol', u * p), ('incident_umol', mu0 * p)):
                if grid == 'uv':
                    out[key][i, lo:hi] = values
                else:
                    out[key][i, lo:hi] = np.bincount(resp.lbl_bin, weights=values, minlength=N_LBL)
    return out


def black_ground_reflectance(resp: Response, light):
    """R in 1-nm bins, weighted within each bin by the black-ground light at the ground with the Sun overhead.

    Exact for the downward light to first order in A R; the tests bound the error at A = 0.8.
    """
    fl = resp.fluxes(1.0)
    uv = resp.below['uv'][0]
    d, f, _ = fl['lbl']
    w = (d + f) * light['lbl'][0]
    lbl = (np.bincount(resp.lbl_bin, weights=w * resp.below['lbl'][0], minlength=N_LBL)
           / np.maximum(np.bincount(resp.lbl_bin, weights=w, minlength=N_LBL), 1e-300))
    return np.concatenate([uv, lbl])


def _band(lo, hi):
    return (CENTRES_NM > lo) & (CENTRES_NM < hi)


VBAR = cmf(CENTRES_NM)[:, 1]
ERYTHEMAL = ph.erythemal_weight(CENTRES_NM)


def summarise(s, i):
    """Band photon fluxes (umol m^-2 s^-1), diffuse shares, ratios, illuminance and UV index for albedo row i."""
    pd, pf = s['direct_umol'][i], s['diffuse_umol'][i]
    out = dict(photons_umol_m2_s={}, diffuse_share={}, top_of_air_umol_m2_s={}, returned_to_space={})
    for name, (lo, hi) in BANDS_NM.items():
        m = _band(lo, hi)
        total = float(pd[m].sum() + pf[m].sum())
        out['photons_umol_m2_s'][name] = total
        out['diffuse_share'][name] = float(pf[m].sum()) / total if total > 0 else None
        inc = float(s['incident_umol'][i][m].sum())
        out['top_of_air_umol_m2_s'][name] = inc
        out['returned_to_space'][name] = float(s['up_umol'][i][m].sum()) / inc if inc > 0 else None
    ph_ = out['photons_umol_m2_s']
    energy = s['direct_w'][i] + s['diffuse_w'][i]
    out.update(par_w_m2=float(energy[_band(*BANDS_NM['par'])].sum()),
               illuminance_lux=float(683.0 * energy @ VBAR),
               direct_illuminance_lux=float(683.0 * s['direct_w'][i] @ VBAR),
               uv_index=float(40.0 * energy @ ERYTHEMAL),
               red_far_red=ph_['red_660'] / ph_['far_red_730'] if ph_['far_red_730'] > 0 else None,
               blue_red=ph_['blue'] / ph_['red'] if ph_['red'] > 0 else None)
    return out


def day_integral(heights_deg, flux, latitude_deg, day_s):
    """Integral of a flux over one solar day (flux units x s) with the Sun on the celestial equator.

    The flux is interpolated linearly in the sine of the Sun's height and taken to fall linearly
    to zero between the lowest height given and the horizon; twilight is not counted.
    """
    x = np.sin(np.radians(np.asarray(heights_deg, dtype=float)))
    order = np.argsort(x)
    xp = np.concatenate([[0.0], x[order]])
    fp = np.concatenate([[0.0], np.asarray(flux, dtype=float)[order]])
    hour = np.linspace(-math.pi / 2, math.pi / 2, 20001)
    sin_h = math.cos(math.radians(latitude_deg)) * np.cos(hour)
    return day_s / (2 * math.pi) * float(np.trapezoid(np.interp(sin_h, xp, fp), hour))


def atlas_comparison(results, atlas_dir):
    """Illuminance against the sky solver's spherical multiple-scattering atlases (unfiltered sunlight, albedo 0.1).

    The atlases use exponential optical profiles without water vapour, and the Moon's without
    ozone; the direct beams agree within a few per cent, so the difference in diffuse light at
    low Sun is the two-stream's plane-parallel diffuse field. daylight_spherical_scaling
    integrates the photosynthetic photons again with each Sun height's light scaled by the
    atlas's total illuminance over this model's: an estimate, taking illuminance for photons.
    """
    atlas_dir = Path(atlas_dir)
    y = np.array([0.2126, 0.7152, 0.0722])
    out = {}
    for key, atlas, scaled in (('moon_1.2atm/unfiltered', 'moon_no_ozone', ('moon_1.2atm/design', 'moon_1.2atm/unfiltered')),
                               ('earth_control/unfiltered', 'earth', ('earth_control/unfiltered',))):
        path = atlas_dir / f'{atlas}_atlas.npz'
        if not path.is_file():
            continue
        d = np.load(path)
        rows = []
        for h in SUN_DEG:
            j = int(np.argmin(np.abs(d['suns'] - h)))
            mine = results[key]['summaries'][str(h)][str(0.1)]
            rows.append(dict(sun_deg=h, atlas_sun_deg=float(d['suns'][j]),
                             direct_lux=mine['direct_illuminance_lux'],
                             atlas_direct_lux=float(d['direct_horizontal'][j] @ y),
                             diffuse_lux=mine['illuminance_lux'] - mine['direct_illuminance_lux'],
                             atlas_diffuse_lux=float(d['diffuse'][j] @ y),
                             total_ratio=float((d['direct_horizontal'][j] + d['diffuse'][j]) @ y) / mine['illuminance_lux']))
        daylight = {}
        for case in scaled:
            res = results[case]
            par = [res['summaries'][str(h)][str(0.1)]['photons_umol_m2_s']['par'] * r['total_ratio']
                   for h, r in zip(SUN_DEG, rows)]
            daylight[case] = {str(lat): dict(par_mol_m2_per_day=day_integral(SUN_DEG, par, lat, res['solar_day_s']) * 1e-6,
                                             par_mol_m2_per_24h=day_integral(SUN_DEG, par, lat, res['solar_day_s'])
                                             * 1e-6 * 86400.0 / res['solar_day_s'])
                              for lat in LATITUDES}
        out[key] = dict(atlas=atlas, atlas_sha256=hashlib.sha256(path.read_bytes()).hexdigest()[:16], rows=rows,
                        daylight_spherical_scaling=daylight)
    return out


EVIDENCE = ('Clear-sky calculation for one-dimensional columns: no clouds, aerosols, terrain or vegetation. The '
            'columns are global means (moist adiabat with Manabe-Wetherald humidity, solved upper air and trace '
            'gases from the middle atmosphere); real humidity varies with place and time. Delta-two-stream '
            'radiative transfer with a pseudo-spherical direct beam; the diffuse light is plane-parallel. Against '
            'the spherical sky solver the Moon\'s total light agrees within 5% down to 30 degrees of Sun height '
            'and falls short below it (by about a quarter at 10 degrees and half at 5 degrees; checks.'
            'atlas_comparison), so values at low Sun understate the light, most of all its diffuse part. Ground '
            'albedos are spectrally neutral and Lambertian. The shield is the titania film at normal incidence '
            'passing 5% less at every wavelength; light through gaps in the swarm (a few 1e-4 of sunlight) is left '
            'out. The unfiltered sunlight on the Moon column is for comparison: that air has no ozone, which an '
            'unshielded Moon would grow, so its ultraviolet is not an unshielded Moon\'s. Physical light only: '
            'what plants do with it (action spectra, saturation, photoinhibition, canopies, the month-long day) '
            'is not modelled here.')
READING_RULE = ('centres_nm are 1-nm bins from 202 to 1000 nm. cases[case]["spectra"][sun_deg] holds direct and '
                'black-ground diffuse irradiance on a horizontal surface (W m^-2 nm^-1). For a Lambertian ground of '
                'albedo A (which may vary with wavelength) the downward light is (direct + diffuse_black) / '
                '(1 - A R), R = columns[column]["reflectance_from_below"], and its diffuse part is that minus '
                'direct; photons per bin are energy x wavelength / (h c). cases[case]["summaries"][sun_deg][albedo] '
                'hold band photon fluxes (umol m^-2 s^-1) and shares computed on the line-by-line grid, exact for '
                'those albedos; returned_to_space is the share of the band\'s sunlight at the top of the air that '
                'the column and ground send back to space. daylight integrates the photosynthetic photons over one '
                'solar day with the Sun on the celestial equator: the Moon\'s day is the synodic month, Earth\'s '
                'is an equinox day.')


def run(processes=2, atlas_dir=None, verbose=True):
    say = print if verbose else (lambda *a, **k: None)
    tsis = {f['name']: f['sha256'][:16] for f in rc_inputs.manifest()['files']}['TSIS1_HSRS_stride100.csv']
    columns, results, shields = {}, {}, {}
    responses = {}
    for name, column in COLUMNS.items():
        t0 = time.time()
        col, air, densities, used = build(column)
        resp = Response(col, air, densities, processes=processes)
        responses[name] = resp
        columns[name] = dict(column_facts(col, air, densities), middle_atmosphere_cases=used)
        say(f'{name}: column and optics in {time.time() - t0:.0f} s', flush=True)
    for name, light_name in CASES:
        t0 = time.time()
        resp = responses[name]
        shield, sha = sunlight(light_name)
        shields[light_name] = dict(shield=SUNLIGHTS[light_name][0], share=SUNLIGHTS[light_name][1],
                                   protection_product_sha256=sha)
        light = top_of_air(resp.nu, shield)
        if 'reflectance_from_below' not in columns[name]:
            columns[name]['reflectance_from_below'] = black_ground_reflectance(resp, light).tolist()
        spec, summ = {}, {}
        for h in SUN_DEG:
            s = spectra(resp, light, math.sin(math.radians(h)))
            spec[str(h)] = dict(direct_w_m2_nm=s['direct_w'][0].tolist(), diffuse_black_w_m2_nm=s['diffuse_w'][0].tolist())
            summ[str(h)] = {str(a): summarise(s, i) for i, a in enumerate(ALBEDOS)}
        day_s = (SYNODIC_MONTH_DAYS if COLUMNS[name].planet is th.MOON else 1.0) * 86400.0
        day = {}
        for a in ALBEDOS:
            par = [summ[str(h)][str(a)]['photons_umol_m2_s']['par'] for h in SUN_DEG]
            day[str(a)] = {str(lat): dict(
                par_mol_m2_per_day=day_integral(SUN_DEG, par, lat, day_s) * 1e-6,
                par_mol_m2_per_24h=day_integral(SUN_DEG, par, lat, day_s) * 1e-6 * 86400.0 / day_s,
                noon_par_umol_m2_s=float(np.interp(math.cos(math.radians(lat)),
                                                   np.sin(np.radians(SUN_DEG[::-1])), par[::-1])))
                for lat in LATITUDES}
        results[f'{name}/{light_name}'] = dict(column=name, sunlight=light_name, solar_day_s=day_s,
                                               spectra=spec, summaries=summ, daylight=day)
        say(f'{name}/{light_name}: {len(SUN_DEG)} Sun heights in {time.time() - t0:.0f} s', flush=True)
    earth = results['earth_control/unfiltered']
    check_case = stored_case('earth_control')[1]
    resp = responses['earth_control']
    s = spectra(resp, top_of_air(resp.nu, None), 1.0, albedos=(0.13,))
    checks = dict(earth_subsolar_uv_index=dict(
        this_model=summarise(s, 0)['uv_index'], middle_atmosphere=check_case['uv_index_subsolar'],
        note='Earth control, Sun overhead, albedo 0.13, as the middle atmosphere reports it'))
    if atlas_dir is not None:
        checks['atlas_comparison'] = atlas_comparison(results, atlas_dir)
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='illumination', files={p: digest(p) for p in (
            'illumination/surface_light/model.py', 'illumination/sky/colour_matching.py',
            'atmosphere/middle_atmosphere/photolysis.py', 'atmosphere/middle_atmosphere/equilibrium.py',
            'atmosphere/radiative_convective/shortwave.py', 'atmosphere/radiative_convective/optics.py',
            'atmosphere/radiative_convective/spectroscopy.py', 'atmosphere/radiative_convective/thermodynamics.py',
            'atmosphere/radiative_convective/climate.py')},
            inputs=dict(tsis1_hsrs_sha256=tsis, middle_atmosphere_schema=MIDDLE_SCHEMA)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        units=dict(spectra='W m^-2 nm^-1 on a horizontal surface', photons='umol m^-2 s^-1',
                   daylight='mol m^-2', illuminance='lux (CIE photopic, 683 lm/W)', uv_index='WHO/CIE'),
        settings=dict(sun_deg=SUN_DEG, albedos=ALBEDOS, latitudes_deg=LATITUDES, bands_nm=BANDS_NM,
                      join_nm=JOIN_NM, line_by_line=dict(nu_min=LBL.nu_min, nu_max=LBL.nu_max, step=LBL.step,
                                                         molecules=LBL.molecules, stride=STRIDE),
                      sunlights=shields),
        centres_nm=CENTRES_NM.tolist(), columns=columns, cases=results, checks=checks)
    return product


def _round(obj, digits=6):
    if isinstance(obj, float):
        return float(f'{obj:.{digits}g}')
    if isinstance(obj, dict):
        return {k: _round(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_round(v, digits) for v in obj]
    return obj


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--processes', type=int, default=2, help='worker processes for the line-by-line optics')
    parser.add_argument('--atlas-dir', type=Path, help='folder holding the sky atlases (illumination/sky/data)')
    args = parser.parse_args(argv)
    product = _round(run(args.processes, args.atlas_dir))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'surface_light.json').write_text(json.dumps(product, separators=(',', ':')) + '\n')
    print('wrote', RESULTS / 'surface_light.json')


if __name__ == '__main__':
    main()
