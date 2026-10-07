"""Heat that short-wave sunlight passing the solar shield deposits in the Open Moon's air along slant paths.

    python -m atmosphere.middle_atmosphere.limb_heat      # writes results/limb_heat.json

The escape model once counted the heat of sunlight below 175 nm as a global mean: a quarter of what reaches the lunar
disk, all of it absorbed above the base (escape.film_heat and escape.leakage_heat). The shield's aperture is wider
than the disk. Its climate window covers the disk; its annulus, out to the protected radius, covers rays that pass the
limb at tangent altitudes of thousands of kilometres, and the Open Moon's air is tall, its exobase two to four lunar
radii from the centre. Those rays cross long slant paths through the upper air and give up their energy where the
column along them reaches the inverse of the absorption cross section. This module traces parallel sunlight through
a spherically symmetric atmosphere and tallies the energy each ray leaves below the base (0.3 Pa), in the thermosphere
between the base and the exobase and in the exosphere above it, for the light that passes each zone of the aperture:

- the climate window, 1 um of titania on 10 um of silica with its coating layers, over the disk and its margin;
- the annulus films, 0.1 um of titania on 2 or 4 um of silica with the same coating layers, or the window stack;
- light that passes gaps, a grey share of the band over the whole protected region;
- unfiltered sunlight beyond the aperture.

A ray's zone follows from where it crosses the shield, the ring fleet's screen near 20,000 km: the solar disk's
angular radius smears that crossing by the shield's distance times that radius, the margin the window is sized for.
Each ray's deposits are kept, so the zones of any protected radius can be summed; the product gives them for protected
radii of 2 to 10 lunar radii, which the loss response reads (atmosphere/loss_response/traced.py), and works through
the ring fleet's 4 lunar radii itself. Rays are packed across the smeared edges of 3-6 lunar radii; nearer the limb the
logarithmic spacing is already finer than the smear.

The atmosphere and its heat follow the loss response (atmosphere/loss_response/model.py): its six cases (the
titania-stack and 200-nm-edge middle atmospheres at 1.2 atm under three treatments of the upper air), its column
configuration, the sky's Lyman-alpha glow and Earth's tide on the molecular loss. Solar maximum is 2.5 times the quiet
Sun's ultraviolet and, below 10 nm, FISM2's measured rise of the X-rays by band (escape.xray_cycle), the mean over the
year around each of the last three solar maxima against the WHI 2008 quiet week, as the author chose on 2026-10-07;
each of the three maxima is also traced. The films are the design's in all six cases; the 200-nm-edge atmospheres
stand for a warmer middle atmosphere. The middle atmosphere's profile (temperature, pressure and atomic oxygen, with
the dry air of radiative_convective.thermodynamics) runs from the surface to the base, the thermal column's solution
(molecular N2 and O2 at their fixed ratio) from there to the exobase, shifted to start at the profile's base height,
and the loss response's exosphere (Chamberlain's ballistic and escaping populations, absorption.exosphere_density)
beyond it, out past the aperture. Absorption cross sections: CXRO atomic scattering factors below 25 nm (independent
atoms), and the Leiden database's photoabsorption continua from 25 nm (Heays, Bosman and van Dishoeck 2017); argon,
CO2 and the discrete N2 bands are left out. The spectrum is the WHI 2008 quiet Sun.

The thermal column takes the heat deposited between the base and the exobase, at the escape model's heating
efficiency, as a global mean and with its own heating shape. Light absorbed in the exosphere makes ions and
photoelectrons in collisionless gas, which the exosphere step counts, and does not heat the column. The heat and the
atmosphere it swells are solved together: a scenario's state is the first heat, counting up from zero, at which the
heat the swollen atmosphere takes equals the heat that swells it. The traced light lands between the column's 'low'
and 'middle' heating shapes in mean log-pressure depth, at about the depth of the disk's own light, so each state's
loss is also given with the 'low' shape as a lower bound. The limb's light heats a ring of upper air over the
terminator, which the one-dimensional column spreads over the globe; the circulation that would carry it is not
modelled.
"""
from __future__ import annotations
import csv
import dataclasses
import hashlib
import inspect
import json
import math
import sys
from pathlib import Path
import numpy as np

from shared.constants import AU, AVOGADRO, BOLTZMANN, CLASSICAL_ELECTRON_RADIUS, MOON_RADIUS, SUN_RADIUS
from shared.provenance import constants_used
from atmosphere.thermal_column import ColumnConfig, MOLAR, solve_column
from atmosphere.radiative_convective import thermodynamics as th
from atmosphere.middle_atmosphere import escape, fetch_inputs, fetch_limb_inputs
from atmosphere.loss_response import absorption, traced, model as loss
from protection.spectra import annulus_film

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
OUT = RESULTS / 'limb_heat.json'
CACHE = HERE / 'cache' / 'limb_heat'
SCHEMA = traced.LIMB_SCHEMA
CASES = {f'{shield}_{treatment}': (shield, treatment) for shield in loss.SHIELDS for treatment in loss.TREATMENTS}
DESIGN = dict(shield_distance_m=20.0e6, protected_radii=4, radii_R=[3, 4, 5, 6, 8, 10],
              summary_radii_R=[2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 10], formation_margin_m=50.0,
              base_pa=loss.BASE_PA, band_nm=175.0, crossover_nm=25.0, gap_transmissions=[3e-5, 2e-4],
              films=dict(window=[1.0, 10.0], annulus_2_um=[0.1, 2.0], annulus_4_um=[0.1, 4.0]),
              disk_rays=64, limb_rays=360, edge_rays=32, limb_first_km=0.5, exosphere_points=240, outer_radii=7.0,
              profile_heats_W_m2=[0., 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6.5e-6, 8e-6, 1e-5, 1.25e-5, 1.5e-5, 2e-5,
                                  2.5e-5],
              largest_exobase_radii=6.0, co2_ppm=400.0, budgets_kg_s=list(loss.BUDGETS_KG_S))
# The heat tables depend on these settings only; the rest steer the analysis of the tables.
TABLE_SETTINGS = ('shield_distance_m', 'protected_radii', 'radii_R', 'formation_margin_m', 'base_pa', 'band_nm',
                  'crossover_nm', 'films', 'disk_rays', 'limb_rays', 'edge_rays', 'limb_first_km', 'exosphere_points',
                  'outer_radii', 'profile_heats_W_m2', 'largest_exobase_radii', 'co2_ppm')
# The loss response's quiet Sun and solar maximum, and each of the three solar maxima FISM2's mean is taken over.
ACTIVITIES = dict(loss.ACTIVITY, **{period: dict(loss.ACTIVITY['solar_maximum'], xray=period)
                                    for period in escape.FISM2_MAXIMA})
MAXIMA = tuple(escape.FISM2_MAXIMA)
MOLECULE = dict(N2=(2, 'n'), O2=(2, 'o'), O=(1, 'o'))
MASS_KG = dict(N2=MOLAR[0] / AVOGADRO, O2=MOLAR[1] / AVOGADRO, O=0.0159994 / AVOGADRO)
LIMB_INPUTS = ('n.nff', 'o.nff', 'N2_0.1nm.txt', 'O2_0.1nm.txt', 'O_0.1nm.txt')
CYCLE_INPUTS = (escape.FISM2_QUIET, *escape.FISM2_MAXIMA.values())
WHI = 'whi2008_ref_solar_irradiance_ver2.dat'
SOURCES = ('window', 'annulus_2_um', 'annulus_4_um', 'open')
# The films an annulus could carry: the two light films and the climate window's own stack.
FILMS = {'annulus_2_um': 'annulus_2_um', 'annulus_4_um': 'annulus_4_um', 'window_stack': 'window'}
ZONES = ('window', 'annulus', 'outside', 'aperture', 'disk')
FILES = {name: ROOT / name for name in (
    'atmosphere/middle_atmosphere/limb_heat.py', 'atmosphere/middle_atmosphere/escape.py',
    'atmosphere/thermal_column.py', 'atmosphere/loss_response/model.py', 'atmosphere/loss_response/absorption.py',
    'atmosphere/loss_response/tides.py', 'atmosphere/loss_response/traced.py', 'protection/spectra/annulus_film.py',
    'protection/spectra/short_wave.py', 'protection/model.py')}
# Code outside this module that the heat tables run through, kept in their key by its source.
EXTERNAL_CODE = (loss.base_conditions, escape._profile, escape.base_temperature, escape.whi_quiet_sun,
                 escape.xray_cycle, escape.xray_scale, absorption.exosphere_density, absorption.partition,
                 annulus_film.indices, annulus_film.film)
TABLE_FILES = ('atmosphere/thermal_column.py', 'atmosphere/radiative_convective/thermodynamics.py',
               'protection/spectra/short_wave.py', 'protection/model.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def cross_sections(wavelength_nm):
    """Photoabsorption cross section per molecule (m^2) of N2, O2 and O: CXRO atoms below the crossover, the Leiden
    continua from it (zero beyond each table's range)."""
    w = np.asarray(wavelength_nm, float)
    out = {}
    for species, (count, element) in MOLECULE.items():
        e, _, f2 = np.loadtxt(fetch_limb_inputs.path(f'{element}.nff'), skiprows=1).T
        energy = 1239.84198 / w
        atom = 2 * CLASSICAL_ELECTRON_RADIUS * w * 1e-9 * np.exp(np.interp(np.log(energy), np.log(e), np.log(f2)))
        table = np.loadtxt(fetch_limb_inputs.path(f'{species}_0.1nm.txt'), comments='#')
        leiden = np.interp(w, table[:, 0], table[:, 1], left=0., right=0.) * 1e-4
        out[species] = np.where(w < DESIGN['crossover_nm'], count * atom, leiden)
    return out


def transmissions(wavelength_nm):
    """Transmission of the window stack and the annulus films (protection domain), normal incidence."""
    silica, titania = annulus_film.indices(wavelength_nm)
    return {name: annulus_film.film(wavelength_nm, silica, titania, ti, si, True)[1]
            for name, (ti, si) in DESIGN['films'].items()}


def case_row(case):
    """The middle-atmosphere row behind a case (shield and treatment)."""
    shield, treatment = CASES[case]
    name = loss.SHIELDS[shield] + loss.TREATMENTS[treatment]
    rows = {r['case']: r for r in csv.DictReader(open(RESULTS / 'middle_atmosphere.csv', newline=''))}
    return rows[name]


def profile_path(case):
    return RESULTS / 'profiles' / f"{case_row(case)['case']}.csv"


def middle_profile(case):
    """The middle atmosphere's layers: radius (m), pressure (Pa), temperature (K), O mixing ratio."""
    rows = list(csv.DictReader(open(profile_path(case), newline='')))
    p = np.array([float(r['p_pa']) for r in rows])
    z = np.array([float(r['z_km']) for r in rows]) * 1e3
    t = np.array([float(r['t_k']) for r in rows])
    o = np.array([float(r['O']) for r in rows])
    return MOON_RADIUS + z, p, t, o


def case_config(case):
    """The thermal column's configuration for a case, as the loss response builds it."""
    base = loss.base_conditions(*CASES[case])
    return ColumnConfig(surface_pressure_pa=base['surface_pressure_pa'],
                        surface_temperature_k=base['surface_temperature_k'],
                        lower_temperature_k=base['base_temperature_k'], lower_pressure_pa=loss.BASE_PA,
                        oxygen_mole_fraction=base['oxygen_mole_fraction'])


def glow_heat(case):
    """The loss response's heat from the sky's Lyman-alpha glow at quiet Sun (W/m^2), at the column's base radius."""
    return loss.lyman_glow_heat(solve_column(0.0, case_config(case))[0]['lower_radius_R'])


def atmosphere(case, heat, previous=None):
    """Number densities (m^-3) of N2, O2 and O on a radial grid from the surface to beyond the aperture, with the
    thermal column's solution above the base for the given heat; the column's summary and warm start."""
    cfg = case_config(case)
    row = case_row(case)
    air = th.earthlike_air(float(row['dry_pressure_pa']), DESIGN['co2_ppm']).fractions
    r_mid, p_mid, t_mid, o_mid = middle_profile(case)
    base = DESIGN['base_pa']
    lp = np.log(p_mid)
    r_base = float(np.interp(math.log(base), lp[::-1], r_mid[::-1]))
    o_base = float(np.exp(np.interp(math.log(base), lp[::-1], np.log(np.maximum(o_mid, 1e-30))[::-1])))
    below = p_mid > base
    summary, column, previous = solve_column(heat, cfg, previous)
    column_base = summary['lower_radius_R'] * MOON_RADIUS
    r_col = column['radius_R'] * MOON_RADIUS - column_base + r_base
    p_col, t_col = column['pressure_Pa'], column['temperature_K']
    # Below the base the profile's dry air (N2 and O2 at their shares of the total) and its atomic oxygen; above it
    # the column's molecular mixture with the base's atomic oxygen carried well mixed.
    o2 = cfg.oxygen_mole_fraction
    r = np.r_[r_mid[below], r_col]
    n = np.r_[p_mid[below] / (BOLTZMANN * t_mid[below]), p_col / (BOLTZMANN * t_col)]
    x_o = np.r_[o_mid[below], np.full(len(r_col), o_base)]
    x_n2 = np.r_[np.full(below.sum(), air['N2']), np.full(len(r_col), 1 - o2)]
    x_o2 = np.r_[np.full(below.sum(), air['O2']), np.full(len(r_col), o2)]
    dens = dict(N2=n * x_n2, O2=n * x_o2, O=n * x_o)
    # The loss response's exosphere beyond the column's top: Chamberlain's ballistic and escaping populations of each
    # species from its exobase density, no satellite particles.
    r_exo, t_exo = r_col[-1], t_col[-1]
    outer = max(DESIGN['outer_radii'] * MOON_RADIUS, 1.5 * r_exo)
    r_out = np.geomspace(r_exo, outer, DESIGN['exosphere_points'])[1:]
    for species, values in dens.items():
        ballistic, escaping = absorption.exosphere_density(r_out / MOON_RADIUS, values[-1], t_exo,
                                                           r_exo / MOON_RADIUS, MASS_KG[species])
        dens[species] = np.r_[values, ballistic + escaping]
    r = np.r_[r, r_out]
    order = np.argsort(r)
    keep = np.r_[True, np.diff(r[order]) > 1.]
    r = r[order][keep]
    dens = {k: np.maximum(v[order][keep], 1e-300) for k, v in dens.items()}
    pressure = np.r_[p_mid[below], p_col, np.zeros(len(r_out))][order][keep]
    return dict(radius_m=r, density=dens, pressure_pa=pressure, base_radius_m=r_base,
                base_temperature_k=cfg.lower_temperature_k, column_base_radius_m=column_base, exobase_radius_m=r_exo,
                summary=summary), previous


def escape_covering(distance, protected):
    """The covering radius of protection/dynamics/optical: the protected radius narrowed by the shield's distance and
    widened by the Sun's disk."""
    return protected * (1 - distance / AU) + SUN_RADIUS * distance / AU


def aperture_edge(protected_R=None):
    """The aperture's edge (m) for a protected radius in lunar radii (default the ring fleet's)."""
    protected_R = DESIGN['protected_radii'] if protected_R is None else protected_R
    return escape_covering(DESIGN['shield_distance_m'], protected_R * MOON_RADIUS) + DESIGN['formation_margin_m']


def impact_parameters():
    """Ray impact parameters (m): equal-area rings over the disk, then rays above the limb out past the aperture,
    with extra rays across the smeared edge of each protected radius the rays reach."""
    k = DESIGN['disk_rays']
    disk = MOON_RADIUS * np.sqrt((np.arange(k) + .5) / k)
    height = np.geomspace(DESIGN['limb_first_km'] * 1e3, (DESIGN['outer_radii'] - 1) * MOON_RADIUS * .995,
                          DESIGN['limb_rays'])
    smear = DESIGN['shield_distance_m'] * SUN_RADIUS / AU
    edges = [aperture_edge(x) + np.linspace(-1.5 * smear, 1.5 * smear, DESIGN['edge_rays'])
             for x in DESIGN['radii_R'] if aperture_edge(x) + 1.5 * smear < MOON_RADIUS * DESIGN['outer_radii'] * .995]
    return np.r_[disk, np.unique(np.r_[MOON_RADIUS + height, np.concatenate(edges)])]


def ring_areas(b):
    """Aperture-plane area (m^2) each ray stands for: the equal-area disk rings, then midpoints above the limb."""
    k = DESIGN['disk_rays']
    disk = np.full(k, math.pi * MOON_RADIUS**2 / k)
    limb = b[k:]
    edges = np.r_[MOON_RADIUS, (limb[1:] + limb[:-1]) / 2, limb[-1] + (limb[-1] - limb[-2]) / 2]
    return np.r_[disk, math.pi * (edges[1:]**2 - edges[:-1]**2)]


def zone_weights(b, protected_R=None):
    """Share of each ray's sunlight that crosses the window, the annulus and the open sky beyond the aperture of a
    protected radius: the solar disk smears each crossing radially by up to the shield's distance times the Sun's
    angular radius."""
    d = DESIGN['shield_distance_m']
    a = d * SUN_RADIUS / AU

    def below(edge):
        c = np.clip(edge - b, -a, a)
        return .5 + (c * np.sqrt(a * a - c * c) + a * a * np.arcsin(c / a)) / (math.pi * a * a)
    w = below(escape_covering(d, MOON_RADIUS))
    ap = below(aperture_edge(protected_R))
    return dict(window=w, annulus=ap - w, outside=1 - ap, aperture=ap, disk=(b < MOON_RADIUS).astype(float))


def log_mean(a, b):
    """Integral weight of an exponential segment between densities a and b, per unit length."""
    ratio = a / b
    close = np.abs(ratio - 1) < 1e-6
    return np.where(close, (a + b) / 2, (a - b) / np.log(np.where(close, 2., ratio)))


def deposit(b, profile, sigma, weights):
    """Energy a ray of impact parameter b leaves in each shell (between consecutive profile radii), on its way in and,
    above the limb, on its way out; the rest reaches the ground (b below the radius) or leaves.

    weights: (sources, wavelengths) spectral power per wavelength step (W/m^2 at normal incidence). Returns
    (sources, shells) deposited power per unit aperture area and (sources,) power reaching the ground."""
    r = profile['radius_m']
    start = max(b, MOON_RADIUS)
    inside = r > start
    rr = np.r_[start, r[inside]]
    s = np.sqrt(np.maximum(rr**2 - b**2, 0.))
    dens = {k: np.r_[np.exp(np.interp(start, r, np.log(v))), v[inside]] for k, v in profile['density'].items()}
    # Column (m^-2) of each species from infinity in to each node, by exponential segments.
    tau = np.zeros((sigma['N2'].size, len(rr)))
    for species, n in dens.items():
        segment = log_mean(n[:-1], n[1:]) * np.diff(s)
        column = np.r_[np.cumsum(segment[::-1])[::-1], 0.]
        tau += np.outer(sigma[species], column)
    inward = np.exp(-tau)                      # intensity on the way in at each node
    gained = inward[:, 1:] - inward[:, :-1]    # absorbed between consecutive nodes going in
    if b < MOON_RADIUS:
        shells = gained
        ground = inward[:, 0]
    else:
        half = tau[:, :1]
        outward = np.exp(-(2 * half - tau))    # intensity on the way out at each node
        shells = gained + (outward[:, :-1] - outward[:, 1:])
        ground = np.zeros(sigma['N2'].size)
    # Each segment of the ray goes to the profile shell that holds its middle; a ray that starts below the profile's
    # first layer gives that stretch to the first shell.
    power = weights @ shells
    shell = np.clip(np.searchsorted(r, .5 * (rr[:-1] + rr[1:]), side='right') - 1, 0, len(r) - 2)
    full = np.zeros((len(r) - 1, weights.shape[0]))
    np.add.at(full, shell, power.T)
    return full.T, weights @ ground


def spectra(w, f, films):
    """Spectral power per wavelength step (W/m^2) of each source at each activity: quiet Sun, the solar maximum the
    escape model takes, and each of the three maxima behind it."""
    step = np.gradient(w)
    transmission = dict(films, open=np.ones_like(w))
    names, rows = [], []
    for name in SOURCES:
        for activity, act in ACTIVITIES.items():
            names.append((name, activity))
            rows.append(f * step * transmission[name] *
                        np.where(w < 10.0, escape.xray_scale(w, act['xray']), act['uv']))
    return names, np.array(rows)


def ray_table(profile, b, weights, sigma, shell_weights=None):
    """Each ray's deposit per unit aperture area (W/m^2 at normal incidence): in the thermosphere between the base and
    the exobase, in the exosphere, below the base and on the ground. With shell_weights (zones, rays: each ray's area
    in each zone) also the power each zone leaves in every shell (W), for the pressure profile."""
    r = profile['radius_m'][:-1]
    thermo = (r >= profile['base_radius_m']) & (r < profile['exobase_radius_m'])
    exo = r >= profile['exobase_radius_m']
    out = {k: np.zeros((weights.shape[0], len(b))) for k in ('thermosphere', 'exosphere', 'below_base', 'ground')}
    shells = None if shell_weights is None else np.zeros((weights.shape[0], len(shell_weights), len(r)))
    for i, bi in enumerate(b):
        power, floor = deposit(bi, profile, sigma, weights)
        out['thermosphere'][:, i] = power[:, thermo].sum(axis=1)
        out['exosphere'][:, i] = power[:, exo].sum(axis=1)
        out['below_base'][:, i] = power[:, ~(thermo | exo)].sum(axis=1)
        out['ground'][:, i] = floor
        if shells is not None:
            shells += power[:, None, :] * shell_weights[None, :, i, None]
    return out, shells


def profile_by_pressure(profile, shells, names, zone_names, picks):
    """Heat between the base and the exobase per decade of pressure (W), where each part of the light heats."""
    r = profile['radius_m']
    p = profile['pressure_pa']
    mid = np.sqrt(np.maximum(p[:-1], 1e-30) * np.maximum(p[1:], 1e-30))
    edges = 10.0**np.arange(-16, 0)
    rows = {}
    for source, activity, zone in picks:
        s = names.index((source, activity))
        z = zone_names.index(zone)
        hist = [float(shells[s, z][(mid >= lo) & (mid < hi) & (r[:-1] >= profile['base_radius_m'])].sum())
                for lo, hi in zip(edges[:-1], edges[1:])]
        rows[f'{source}_{zone}_{activity}'] = hist
    return dict(pressure_decades_pa=[[float(lo), float(hi)] for lo, hi in zip(edges[:-1], edges[1:])], watts=rows)


def tables(case, b, sigma, names, weights, key):
    """The heat tables at the profile heats, from the quiescent air up until the exobase passes the largest radius
    or the column leaves its domain: each heat's state and each ray's deposits, kept in the cache under the key of the
    code, settings, inputs and products they come from."""
    meta_path, ray_path = CACHE / f'{case}.json', CACHE / f'{case}.npz'
    if meta_path.is_file() and ray_path.is_file():
        saved = json.loads(meta_path.read_text())
        if saved.get('key') == key and saved.get('complete'):
            with np.load(ray_path) as data:
                return saved['entries'], {k: data[k] for k in data.files}
    areas = ring_areas(b)
    design = zone_weights(b)
    shell_weights = np.array([design[z] * areas for z in ZONES])
    entries, previous = [], None
    arrays = {k: [] for k in ('thermosphere', 'exosphere', 'below_base', 'ground')}
    for q in DESIGN['profile_heats_W_m2']:
        try:
            profile, previous = atmosphere(case, q, previous)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            entries.append(dict(heat_W_m2=q, failed=str(exc)))
            break
        summary = profile['summary']
        rays, shells = ray_table(profile, b, weights, sigma, shell_weights if q == 0. else None)
        for k, v in rays.items():
            arrays[k].append(v)
        entry = dict(heat_W_m2=q, exobase_radius_R=profile['exobase_radius_m'] / MOON_RADIUS,
                     column_exobase_radius_R=summary['exobase_radius_R'],
                     exobase_temperature_k=summary['exobase_temperature_k'],
                     base_radius_km=profile['base_radius_m'] / 1e3,
                     column_base_shift_km=(profile['base_radius_m'] - profile['column_base_radius_m']) / 1e3,
                     outflow_limit='HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED' in summary['domain_flags'])
        if q == 0.:
            entry['by_pressure'] = profile_by_pressure(
                profile, shells, names, list(ZONES),
                [('annulus_4_um', 'solar_maximum', 'annulus'), ('annulus_2_um', 'solar_maximum', 'annulus'),
                 ('window', 'solar_maximum', 'window'), ('open', 'quiet', 'aperture'), ('open', 'quiet', 'disk')])
        entries.append(entry)
        print(f"  {case} heat {q:.2e}: exobase {entry['exobase_radius_R']:.2f} R", flush=True)
        if entry['outflow_limit'] or entry['exobase_radius_R'] > DESIGN['largest_exobase_radii']:
            break
    CACHE.mkdir(parents=True, exist_ok=True)
    stacked = {k: np.array(v) for k, v in arrays.items()}
    np.savez_compressed(ray_path, **stacked)
    meta_path.write_text(json.dumps(dict(key=key, complete=True, entries=entries)))
    return entries, stacked


def zone_summary(rays, areas, zones, incident, names, detailed):
    """Heat in the thermosphere for each source and zone, as a global mean at the escape model's heating efficiency
    (W/m^2 of lunar surface). In detail, for quiet Sun and solar maximum, also the heat above the base, the light
    absorbed in the exosphere and the fate of the light arriving."""
    area = 4 * math.pi * MOON_RADIUS**2
    eff = escape.HEATING_EFFICIENCY
    out = {}
    for zone in ZONES:
        weight = areas * zones[zone]
        parts = {k: v @ weight for k, v in rays.items()}
        arriving = incident * weight.sum()
        for s, (source, activity) in enumerate(names):
            thermo = float(eff * parts['thermosphere'][s] / area)
            if not detailed:
                out.setdefault(source, {}).setdefault(zone, {})[activity] = thermo
                continue
            if activity not in loss.ACTIVITY:
                continue
            a = arriving[s]
            share = (lambda x: float(x / a)) if a > 0 else (lambda x: 0.)
            deposited = parts['thermosphere'][s] + parts['exosphere'][s] + parts['below_base'][s]
            out.setdefault(source, {}).setdefault(zone, {})[activity] = dict(
                thermosphere_W_m2=thermo,
                above_base_W_m2=float(eff * (parts['thermosphere'][s] + parts['exosphere'][s]) / area),
                exosphere_absorbed_W_m2=float(parts['exosphere'][s] / area), arriving_W_m2=float(a / area),
                shares=dict(thermosphere=share(parts['thermosphere'][s]), exosphere=share(parts['exosphere'][s]),
                            below_base=share(parts['below_base'][s]), ground=share(parts['ground'][s]),
                            passing=share(a - deposited - parts['ground'][s])))
    return out


def build_entries(meta, arrays, b, incident, names):
    """Each heat's state with its zone summaries: compact for every protected radius, in detail for the ring
    fleet's."""
    areas = ring_areas(b)
    zones = {x: zone_weights(b, x) for x in DESIGN['summary_radii_R']}
    entries, i = [], 0
    for entry in meta:
        entry = dict(entry)
        if 'failed' not in entry:
            rays = {k: v[i] for k, v in arrays.items()}
            entry['radii'] = {f'{x:g}': zone_summary(rays, areas, zones[x], incident, names, False)
                              for x in DESIGN['summary_radii_R']}
            entry['table'] = zone_summary(rays, areas, zones[DESIGN['protected_radii']], incident, names, True)
            i += 1
        entries.append(entry)
    return entries


def escape_count(w, f, films):
    """The escape model's count before tracing: a quarter of the band's light through each film, all of it above the
    base (escape.film_heat's rule); escape.film_heat itself for the window; and escape.leakage_heat for the whole band
    over the disk."""
    step = np.gradient(w)
    eff = escape.HEATING_EFFICIENCY
    out = {name: {activity: float(.25 * eff * np.sum(f * t * step *
                                                       np.where(w < 10.0, escape.xray_scale(w, act['xray']),
                                                                act['uv'])))
                  for activity, act in loss.ACTIVITY.items()} for name, t in films.items()}
    out['window_film_heat'] = {activity: escape.film_heat(act['uv'], act['xray'])
                               for activity, act in loss.ACTIVITY.items()}
    out['band_over_disk'] = {activity: escape.leakage_heat(1.0, activity=act['uv'], xray_activity=act['xray'])
                             for activity, act in loss.ACTIVITY.items()}
    return out


def cycle_factors():
    """FISM2's rise of the X-rays below 10 nm by band (escape.xray_cycle), and weighted by the light each film passes
    and by the whole band below 10 nm, for the mean the escape model takes and for each maximum."""
    w, f = escape.whi_quiet_sun()
    sel = w < 10.0
    w, f = w[sel], f[sel]
    weights = dict(transmissions(w), open=np.ones_like(w))
    sources = {}
    for period in ('fism2', *MAXIMA):
        scale = escape.xray_scale(w, period)
        sources[period] = {k: float((t * f * scale).sum() / (t * f).sum()) for k, t in weights.items()}
    return dict(bands=escape.xray_cycle(), sources=sources)


def molecular_loss(cfg, heat):
    """The thermal column at a heat: exobase and molecular loss with Earth's tide, or why it has none."""
    try:
        summary, _, _ = solve_column(heat, cfg)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        return dict(heat_W_m2=heat, status='failed', note=str(exc))
    if 'HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED' in summary['domain_flags']:
        return dict(heat_W_m2=heat, status='outside_thermal_column_domain')
    tidal, multiplier = loss.with_tides(summary)
    return dict(heat_W_m2=heat, status='thermal_column', exobase_temperature_k=summary['exobase_temperature_k'],
                exobase_radius_R=summary['exobase_radius_R'], molecular_loss_kg_s=tidal,
                molecular_loss_two_body_kg_s=summary['molecular_loss_kg_s'], tidal_multiplier_N2=multiplier)


def heat_for_loss(cfg, budget, ceiling):
    """The heat (W/m^2) at which the molecular loss with Earth's tide reaches the budget; None above the ceiling."""
    def lost(q):
        row = molecular_loss(cfg, q)
        return row.get('molecular_loss_kg_s', math.inf)
    if lost(ceiling) < budget:
        return None
    lo, hi = 0., ceiling
    for _ in range(60):
        mid = .5 * (lo + hi)
        lo, hi = (mid, hi) if lost(mid) < budget else (lo, mid)
    return .5 * (lo + hi)


def disk_count(case, activity, counts, glow):
    """The fixed heat of the escape model's count before tracing: the window film behind the titania stack (nothing
    behind the 200-nm edge) and the glow."""
    shield, _ = CASES[case]
    film = counts['window_film_heat'][activity] if shield == 'titania_stack' else 0.
    return film + glow * loss.ACTIVITY[activity]['glow']


def state_row(cfg, low, q, p, gaps):
    """A scenario's traced state, its parts and its loss, with the loss if the heat lay as low as the column's 'low'
    shape puts it."""
    state = traced.first_state(q, p, gaps)
    if state is None:
        return dict(status='runaway_beyond_profile_heats', beyond_W_m2=float(q[-1]))
    bound = molecular_loss(low, state)
    kept = ('status', 'exobase_temperature_k', 'molecular_loss_kg_s')
    return dict(molecular_loss(cfg, state), parts_W_m2=traced.shares(q, p, gaps, state),
                low_shape={k: bound[k] for k in kept if k in bound})


def analyse(case, entries, counts):
    """Scenario states, their losses, and the gap transmission each budget allows at the ring fleet's protected
    radius, traced and by the escape model's count before tracing; for each of the three maxima as well."""
    cfg = case_config(case)
    low = dataclasses.replace(cfg, heating_shape='low')
    glow = glow_heat(case)
    unit = counts['band_over_disk']
    ceiling = max(e['heat_W_m2'] for e in entries if 'radii' in e and not e['outflow_limit'])
    heats = {budget: heat_for_loss(cfg, budget, ceiling) for budget in DESIGN['budgets_kg_s']}
    scenarios, budgets = {}, []
    for film, source in FILMS.items():
        for activity, act in ACTIVITIES.items():
            q, p = traced.parts(entries, source, activity, glow * act['glow'], DESIGN['protected_radii'])
            counted = activity in loss.ACTIVITY
            for gaps in [0.] + DESIGN['gap_transmissions']:
                row = dict(film=film, activity=activity, gap_transmission=gaps,
                           traced=state_row(cfg, low, q, p, gaps))
                if counted:
                    row['disk_count'] = molecular_loss(cfg, disk_count(case, activity, counts, glow)
                                                       + gaps * unit[activity])
                scenarios[f'{film}_{activity}_gaps_{gaps:g}'] = row
            for budget, heat in heats.items():
                row = dict(film=film, activity=activity, budget_kg_s=budget, heat_W_m2=heat)
                if heat is not None:
                    row['traced'] = traced.allowed(q, p, heat)
                    if counted:
                        f = (heat - disk_count(case, activity, counts, glow)) / unit[activity]
                        row['disk_count'] = f if f > 0 else None
                        if row['traced'] and row['disk_count']:
                            row['traced_over_disk_count'] = row['traced'] / row['disk_count']
                        row['gap_heat_over_disk_count'] = float(np.interp(heat, q, p['gap_unit']) / unit[activity])
                budgets.append(row)
    return dict(glow_heat_quiet_W_m2=glow, scenarios=scenarios, budgets=budgets)


TABLE_CODE = (cross_sections, transmissions, case_row, profile_path, middle_profile, case_config, atmosphere,
              escape_covering, aperture_edge, impact_parameters, ring_areas, zone_weights, log_mean, deposit, spectra,
              ray_table, profile_by_pressure, tables)


def table_key(inputs, products, constants):
    """Digest of everything the heat tables depend on: the functions that make them, here and outside, their
    settings and values, the other code files they run, the shared constants, the inputs and the products they
    read."""
    source = [inspect.getsource(f) for f in TABLE_CODE + EXTERNAL_CODE]
    settings = {k: DESIGN[k] for k in TABLE_SETTINGS}
    values = dict(shields=loss.SHIELDS, treatments=loss.TREATMENTS, base_pa=loss.BASE_PA, activities=ACTIVITIES,
                  solar_maximum=escape.SOLAR_MAXIMUM, xray_bands=escape.XRAY_BANDS_NM,
                  fism2=[escape.FISM2_QUIET, escape.FISM2_MAXIMA], film_split_nm=annulus_film.SPLIT_NM)
    files = {name: digest(ROOT / name) for name in TABLE_FILES}
    blob = json.dumps([source, settings, values, CASES, MOLECULE, MASS_KG, SOURCES, ZONES, files, constants, inputs,
                       products], sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def provenance():
    """Digests of the code, inputs and products, the shared constants read, and the heat tables' key."""
    code = {name: digest(path) for name, path in FILES.items()}
    inputs = {n: digest(fetch_limb_inputs.path(n)) for n in LIMB_INPUTS + CYCLE_INPUTS}
    inputs[WHI] = digest(fetch_inputs.path(WHI))
    profiles = {'atmosphere/middle_atmosphere/results/middle_atmosphere.csv': digest(RESULTS / 'middle_atmosphere.csv')}
    for case in CASES:
        profiles[str(profile_path(case).relative_to(ROOT))] = digest(profile_path(case))
    products = dict(profiles, **{'atmosphere/loss_response/results/tidal_escape.json': digest(loss.TIDES),
                                 'protection/spectra/stack_short_wave.json': digest(escape.FILM_PRODUCT)})
    constants = constants_used(list(FILES.values()))
    return code, inputs, products, constants, table_key(inputs, profiles, constants)


def main(argv=None) -> int:
    whi_w, whi_f = escape.whi_quiet_sun()
    sel = (whi_w > 0.1) & (whi_w < DESIGN['band_nm'])
    w, f = whi_w[sel], whi_f[sel]
    films = transmissions(w)
    sigma = cross_sections(w)
    names, weights = spectra(w, f, films)
    b = impact_parameters()
    incident = weights.sum(axis=1)
    counts = escape_count(w, f, films)
    code, inputs, products, constants, key = provenance()
    cases = {}
    for case in CASES:
        meta, arrays = tables(case, b, sigma, names, weights, key)
        entries = build_entries(meta, arrays, b, incident, names)
        result = analyse(case, entries, counts)
        cases[case] = dict(by_profile_heat=entries, **result)
        print(case, flush=True)
        for row in result['budgets']:
            if row['film'] == 'annulus_4_um':
                print('   ', {k: (f'{v:.3g}' if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    product = dict(
        schema=SCHEMA, producer=dict(domain='atmosphere', files=code, constants=constants, inputs=inputs,
                                     products=products),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=(
            'Heats are W per m^2 of lunar surface, global means at the escape model\'s heating efficiency, as the '
            'thermal column takes them. In each case\'s by_profile_heat, "radii"[R][source][zone][activity] is the '
            'heat in the thermosphere, between the 0.3 Pa base and the exobase, at a protected radius of R lunar '
            'radii; "table" gives it in detail at the ring fleet\'s 4 lunar radii, with the light absorbed in the '
            'exosphere ("exosphere_absorbed", which the column does not take) and the fate of the light arriving. '
            'Activities: quiet, solar_maximum (2.5 on the ultraviolet and FISM2\'s mean rise below 10 nm) and each '
            'maximum behind that mean (cycle_23, cycle_24, cycle_25). Transmissions through gaps are grey shares of '
            'the band below 175 nm. A traced state puts the window stack on the window, a film on the annulus, the '
            'gaps\' light over the whole aperture, unfiltered sunlight beyond it and the sky\'s glow on the air, and '
            'is the first heat from zero that the air it swells returns. "disk_count" is the escape model\'s count '
            'before tracing (a quarter of the light that reaches the disk, all of it above the base). "low_shape" '
            'is the loss at the same heat with the column\'s low heating shape. "window_stack" puts the climate '
            'window\'s stack on the annulus too. A budget row\'s transmissions are the largest grey transmissions '
            'through gaps whose molecular loss with Earth\'s tide stays within the budget at the ring fleet\'s 4 '
            'lunar radii, the solar wind and the exosphere step\'s losses not deducted.'),
        design=DESIGN, activities=ACTIVITIES, escape_count=counts, xray_cycle=cycle_factors(), cases=cases)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
