"""Heat that short-wave sunlight passing the solar shield deposits in the Open Moon's air along slant paths.

    python -m atmosphere.middle_atmosphere.limb_heat      # writes results/limb_heat.json

The escape model counts the heat of sunlight below 175 nm as a global mean: a quarter of what reaches the lunar
disk, all of it absorbed above the base (escape.film_heat and escape.leakage_heat), and the loss response derives
O1's allowed transmissions with that count. The shield's aperture is wider than the disk. Its climate window covers
the disk; its annulus, out to four lunar radii, covers rays that pass the limb at tangent altitudes up to about
5,200 km, and the Open Moon's air is tall, its exobase two to four lunar radii from the centre. Those rays cross long
slant paths through the upper air and give up their energy where the column along them reaches the inverse of the
absorption cross section. This module traces parallel sunlight through a spherically symmetric atmosphere and
tallies the energy each ray leaves in every shell, below the base (0.3 Pa), in the thermosphere between the base and
the exobase and in the exosphere above it, for the light that passes each zone of the aperture:

- the climate window, 1 um of titania on 10 um of silica with its coating layers, over the disk and its margin;
- the annulus films, 0.1 um of titania on 2 or 4 um of silica with the same coating layers, or the window stack;
- light that bypasses the films through gaps, a grey share of the band over the whole protected region;
- unfiltered sunlight beyond the aperture.

A ray's zone follows from where it crosses the shield (a ring screen near 20,000 km): the solar disk's angular
radius smears that crossing by the shield's distance times that radius, the margin the window is sized for.

The atmosphere and its heat follow the loss response (atmosphere/loss_response/model.py): its six cases (the
titania-stack and 200-nm-edge middle atmospheres at 1.2 atm under three treatments of the upper air), its column
configuration, the sky's Lyman-alpha glow, solar maximum as 2.5 times the quiet Sun's ultraviolet and 100 times its
X-rays below 10 nm, and Earth's tide on the molecular loss. The films are the design's in all six; the 200-nm-edge
atmospheres stand for a warmer middle atmosphere. The films pass only X-rays near 1 nm, so their heat at solar
maximum follows the X-rays' solar-cycle factor. FISM2's daily spectra (Chamberlin et al. 2020) give that factor for
each film's band and for the whole band below 10 nm, as the year around each of the last three solar maxima over the
WHI 2008 quiet week, and the states are also given with those factors. The middle atmosphere's profile (temperature,
pressure and atomic oxygen, with the dry air of radiative_convective.thermodynamics) runs from the surface to the
base, the thermal column's solution (molecular N2 and O2 at their fixed ratio) from there to the exobase, shifted to
start at the profile's base height, and the loss response's exosphere (Chamberlain's ballistic and escaping
populations, absorption.exosphere_density) beyond it, out past the aperture. Absorption cross sections: CXRO atomic
scattering factors below 25 nm (independent atoms), and the Leiden database's photoabsorption continua from 25 nm
(Heays, Bosman and van Dishoeck 2017); argon, CO2 and the discrete N2 bands are left out. The spectrum is the WHI
2008 quiet Sun.

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
from atmosphere.loss_response import absorption, model as loss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
OUT = RESULTS / 'limb_heat.json'
CACHE = HERE / 'cache' / 'limb_heat'
SCHEMA = 'terluna.atmosphere.limb-heat/1'
CASES = {f'{shield}_{treatment}': (shield, treatment) for shield in loss.SHIELDS for treatment in loss.TREATMENTS}
DESIGN = dict(shield_distance_m=20.0e6, protected_radii=4, formation_margin_m=50.0, base_pa=loss.BASE_PA,
              band_nm=175.0, crossover_nm=25.0, gap_transmissions=[3e-5, 2e-4],
              films=dict(window=[1.0, 10.0], annulus_2_um=[0.1, 2.0], annulus_4_um=[0.1, 4.0]),
              disk_rays=64, limb_rays=360, edge_rays=48, limb_first_km=0.5, exosphere_points=240, outer_radii=7.0,
              profile_heats_W_m2=[0., 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6.5e-6, 8e-6, 1e-5, 1.25e-5, 1.5e-5, 2e-5,
                                  2.5e-5],
              largest_exobase_radii=6.0, co2_ppm=400.0, budgets_kg_s=list(loss.BUDGETS_KG_S))
# The heat tables depend on these settings only; the rest steer the analysis of the tables.
TABLE_SETTINGS = ('shield_distance_m', 'protected_radii', 'formation_margin_m', 'base_pa', 'band_nm', 'crossover_nm',
                  'films', 'disk_rays', 'limb_rays', 'edge_rays', 'limb_first_km', 'exosphere_points', 'outer_radii',
                  'profile_heats_W_m2', 'largest_exobase_radii', 'co2_ppm')
MOLECULE = dict(N2=(2, 'n'), O2=(2, 'o'), O=(1, 'o'))
MASS_KG = dict(N2=MOLAR[0] / AVOGADRO, O2=MOLAR[1] / AVOGADRO, O=0.0159994 / AVOGADRO)
LIMB_INPUTS = ('n.nff', 'o.nff', 'N2_0.1nm.txt', 'O2_0.1nm.txt', 'O_0.1nm.txt')
# FISM2 daily X-ray spectra: the WHI 2008 quiet-Sun week, and the year around each of the last three solar maxima.
CYCLE_QUIET = 'fism2_xray_2008-04-10_2008-04-16.csv'
CYCLE_MAXIMA = {'cycle_23': 'fism2_xray_2001-05_2002-04.csv', 'cycle_24': 'fism2_xray_2013-10_2014-09.csv',
                'cycle_25': 'fism2_xray_2024-04_2025-03.csv'}
WHI = 'whi2008_ref_solar_irradiance_ver2.dat'
# 'open_flat' is the open band with its X-rays at the ultraviolet's solar-maximum factor, as escape.leakage_heat
# scales the band a filter passes; 'open' gives them their own factor, as escape.film_heat does.
SOURCES = ('window', 'annulus_2_um', 'annulus_4_um', 'open', 'open_flat')
# The films an annulus could carry: the two light films and the climate window's own stack.
FILMS = {'annulus_2_um': 'annulus_2_um', 'annulus_4_um': 'annulus_4_um', 'window_stack': 'window'}
FILES = {name: ROOT / name for name in (
    'atmosphere/middle_atmosphere/limb_heat.py', 'atmosphere/middle_atmosphere/escape.py',
    'atmosphere/thermal_column.py', 'atmosphere/loss_response/model.py', 'atmosphere/loss_response/absorption.py',
    'atmosphere/loss_response/tides.py', 'protection/spectra/annulus_film.py', 'protection/spectra/short_wave.py',
    'protection/model.py')}


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
    from protection.spectra import annulus_film as af
    silica, titania = af.indices(wavelength_nm)
    return {name: af.film(wavelength_nm, silica, titania, ti, si, True)[1]
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


def aperture_edge():
    return (escape_covering(DESIGN['shield_distance_m'], DESIGN['protected_radii'] * MOON_RADIUS)
            + DESIGN['formation_margin_m'])


def impact_parameters():
    """Ray impact parameters (m): equal-area rings over the disk, then rays above the limb out past the aperture,
    with extra rays across the aperture's smeared edge."""
    k = DESIGN['disk_rays']
    disk = MOON_RADIUS * np.sqrt((np.arange(k) + .5) / k)
    height = np.geomspace(DESIGN['limb_first_km'] * 1e3, (DESIGN['outer_radii'] - 1) * MOON_RADIUS * .995,
                          DESIGN['limb_rays'])
    smear = DESIGN['shield_distance_m'] * SUN_RADIUS / AU
    edge = aperture_edge() + np.linspace(-1.5 * smear, 1.5 * smear, DESIGN['edge_rays'])
    return np.r_[disk, np.unique(np.r_[MOON_RADIUS + height, edge])]


def ring_areas(b):
    """Aperture-plane area (m^2) each ray stands for: the equal-area disk rings, then midpoints above the limb."""
    k = DESIGN['disk_rays']
    disk = np.full(k, math.pi * MOON_RADIUS**2 / k)
    limb = b[k:]
    edges = np.r_[MOON_RADIUS, (limb[1:] + limb[:-1]) / 2, limb[-1] + (limb[-1] - limb[-2]) / 2]
    return np.r_[disk, math.pi * (edges[1:]**2 - edges[:-1]**2)]


def zone_weights(b):
    """Share of each ray's sunlight that crosses the window, the annulus and the open sky beyond the aperture: the
    solar disk smears each crossing radially by up to the shield's distance times the Sun's angular radius."""
    d = DESIGN['shield_distance_m']
    a = d * SUN_RADIUS / AU

    def below(edge):
        c = np.clip(edge - b, -a, a)
        return .5 + (c * np.sqrt(a * a - c * c) + a * a * np.arcsin(c / a)) / (math.pi * a * a)
    w = below(escape_covering(d, MOON_RADIUS))
    ap = below(aperture_edge())
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
    """Spectral power per wavelength step (W/m^2) of each source at quiet Sun and solar maximum."""
    step = np.gradient(w)
    transmission = dict(films, open=np.ones_like(w), open_flat=np.ones_like(w))
    names, rows = [], []
    for name in SOURCES:
        for activity, act in loss.ACTIVITY.items():
            xray = act['uv'] if name == 'open_flat' else act['xray']
            names.append((name, activity))
            rows.append(f * step * transmission[name] * np.where(w < 10.0, xray, act['uv']))
    return names, np.array(rows)


def heat_table(profile, rays, weights, sigma):
    """Deposited power (W) in each shell and on the ground, and the power arriving, from every source over every zone
    of the aperture."""
    b, areas, zones = rays
    zone_names = list(zones)
    shells = np.zeros((weights.shape[0], len(zone_names), len(profile['radius_m']) - 1))
    ground = np.zeros((weights.shape[0], len(zone_names)))
    incident = np.zeros_like(ground)
    supply = weights.sum(axis=1)
    for i, bi in enumerate(b):
        power, floor = deposit(bi, profile, sigma, weights)
        for z, name in enumerate(zone_names):
            share = zones[name][i] * areas[i]
            if share > 0:
                shells[:, z] += power * share
                ground[:, z] += floor * share
                incident[:, z] += supply * share
    return shells, ground, incident, zone_names


def summarise(profile, shells, ground, incident, names, zone_names):
    """Heat in the thermosphere and above the base for each source and zone, as a global mean at the escape model's
    heating efficiency; light absorbed in the exosphere as a global mean; and the fate of the light arriving."""
    r = profile['radius_m'][:-1]
    thermo = (r >= profile['base_radius_m']) & (r < profile['exobase_radius_m'])
    exo = r >= profile['exobase_radius_m']
    area = 4 * math.pi * MOON_RADIUS**2
    eff = escape.HEATING_EFFICIENCY
    out = {}
    for s, (source, activity) in enumerate(names):
        for z, zone in enumerate(zone_names):
            row, arriving = shells[s, z], incident[s, z]
            share = (lambda x: float(x / arriving)) if arriving > 0 else (lambda x: 0.)
            out.setdefault(source, {}).setdefault(zone, {})[activity] = dict(
                thermosphere_W_m2=float(eff * row[thermo].sum() / area),
                above_base_W_m2=float(eff * row[thermo | exo].sum() / area),
                exosphere_absorbed_W_m2=float(row[exo].sum() / area),
                arriving_W_m2=float(arriving / area),
                shares=dict(thermosphere=share(row[thermo].sum()), exosphere=share(row[exo].sum()),
                            below_base=share(row[~(thermo | exo)].sum()), ground=share(ground[s, z]),
                            passing=share(arriving - row.sum() - ground[s, z])))
    return out


def escape_count(w, f, films):
    """The escape model's count: a quarter of the band's light through each film, all of it above the base
    (escape.film_heat's rule); escape.film_heat itself for the window; and escape.leakage_heat for the whole band
    over the disk."""
    step = np.gradient(w)
    eff = escape.HEATING_EFFICIENCY
    out = {name: {activity: float(.25 * eff * np.sum(f * t * step * np.where(w < 10.0, act['xray'], act['uv'])))
                  for activity, act in loss.ACTIVITY.items()} for name, t in films.items()}
    out['window_film_heat'] = {activity: escape.film_heat(act['uv'], act['xray'])
                               for activity, act in loss.ACTIVITY.items()}
    out['band_over_disk'] = {activity: escape.leakage_heat(1.0, activity=act['uv'])
                             for activity, act in loss.ACTIVITY.items()}
    return out


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


def tables(case, rays, sigma, names, weights, key):
    """The heat tables at the profile heats, from the quiescent air up until the exobase passes the largest radius
    or the column leaves its domain; kept in the cache under the code's and inputs' digest."""
    path = CACHE / f'{case}.json'
    if path.is_file():
        saved = json.loads(path.read_text())
        if saved.get('key') == key and saved.get('complete'):
            return saved['entries']
    entries, previous = [], None
    for q in DESIGN['profile_heats_W_m2']:
        try:
            profile, previous = atmosphere(case, q, previous)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            entries.append(dict(heat_W_m2=q, failed=str(exc)))
            break
        summary = profile['summary']
        shells, ground, incident, zone_names = heat_table(profile, rays, weights, sigma)
        entry = dict(heat_W_m2=q, exobase_radius_R=profile['exobase_radius_m'] / MOON_RADIUS,
                     column_exobase_radius_R=summary['exobase_radius_R'],
                     exobase_temperature_k=summary['exobase_temperature_k'],
                     base_radius_km=profile['base_radius_m'] / 1e3,
                     column_base_shift_km=(profile['base_radius_m'] - profile['column_base_radius_m']) / 1e3,
                     outflow_limit='HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED' in summary['domain_flags'],
                     table=summarise(profile, shells, ground, incident, names, zone_names))
        if q == 0.:
            entry['by_pressure'] = profile_by_pressure(
                profile, shells, names, zone_names,
                [('annulus_4_um', 'solar_maximum', 'annulus'), ('annulus_2_um', 'solar_maximum', 'annulus'),
                 ('window', 'solar_maximum', 'window'), ('open', 'quiet', 'aperture'), ('open', 'quiet', 'disk')])
        entries.append(entry)
        print(f"  {case} heat {q:.2e}: exobase {entry['exobase_radius_R']:.2f} R", flush=True)
        if entry['outflow_limit'] or entry['exobase_radius_R'] > DESIGN['largest_exobase_radii']:
            break
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(key=key, complete=True, entries=entries)))
    return entries


def parts(entries, film, activity, glow, gap_source='open', xray=None):
    """Heat (W/m^2) each part of the light leaves in the thermosphere against the heat that swells the air: the
    window stack, the annulus film, sunlight beyond the aperture, the sky's glow, and the whole band through gaps
    over the aperture, per unit of transmission.

    xray: at solar maximum, a factor on the X-rays below 10 nm in place of the convention's, or one for each source
    (window, the films, open). Each table is linear in the spectrum, so the X-rays' part of a solar-maximum value is
    (maximum - uv * quiet) / (convention - uv)."""
    rows = [e for e in entries if 'table' in e and not e['outflow_limit']]
    q = np.array([e['heat_W_m2'] for e in rows])
    act = loss.ACTIVITY[activity]

    def get(source, zone):
        value = np.array([e['table'][source][zone][activity]['thermosphere_W_m2'] for e in rows])
        if xray is None:
            return value
        quiet = np.array([e['table'][source][zone]['quiet']['thermosphere_W_m2'] for e in rows])
        xrays = (value - act['uv'] * quiet) / (act['xray'] - act['uv'])
        factor = xray[source] if isinstance(xray, dict) else xray
        return act['uv'] * quiet + (factor - act['uv']) * xrays
    return q, dict(window=get('window', 'window'), annulus=get(FILMS[film], 'annulus'),
                   beyond_aperture=get('open', 'outside'), glow=np.full(len(q), glow * act['glow']),
                   gap_unit=get(gap_source, 'aperture'))


def cycle_factors():
    """Solar-cycle factors of the X-rays below 10 nm from FISM2: each maximum year's mean over the WHI quiet week's,
    for the light each film passes (weighted by its transmission) and for the whole band below 10 nm (open)."""
    def mean_spectrum(name):
        rows = np.genfromtxt(fetch_limb_inputs.path(name), delimiter=',', skip_header=1)
        w = np.unique(rows[:, 1])
        return w, rows[:, 2].reshape(-1, len(w)).mean(axis=0)
    w, quiet = mean_spectrum(CYCLE_QUIET)
    weights = dict(transmissions(w), open=np.ones_like(w))
    bands = dict(hard_0p1_2p5_nm=(0.1, 2.5), soft_2p5_10_nm=(2.5, 10.0))
    out = {}
    for period, name in CYCLE_MAXIMA.items():
        w2, peak = mean_spectrum(name)
        assert np.allclose(w2, w)
        out[period] = dict(sources={k: float((peak * t).sum() / (quiet * t).sum()) for k, t in weights.items()},
                           bands={k: float(peak[(w >= lo) & (w < hi)].sum() / quiet[(w >= lo) & (w < hi)].sum())
                                  for k, (lo, hi) in bands.items()})
    return out


def first_state(q, p, gaps):
    """The first heat, counting up from zero, at which the heat the swollen air takes equals the heat that swells it;
    None if the heat runs past the largest profile heat."""
    fine = np.linspace(0., q[-1], 4001)
    taken = sum(np.interp(fine, q, p[k]) for k in ('window', 'annulus', 'beyond_aperture', 'glow'))
    gap = gaps * np.interp(fine, q, p['gap_unit'])
    excess = taken + gap - fine
    settled = np.nonzero(excess <= 0)[0]
    if not len(settled):
        return None
    i = settled[0]
    return float(fine[i - 1] + excess[i - 1] * (fine[i] - fine[i - 1]) / (excess[i - 1] - excess[i]))


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


def allowed_traced(q, p, heat):
    """Largest transmission through gaps whose state stays at or below the heat; None if the films, beyond-aperture
    light and glow alone
    pass it."""
    state = first_state(q, p, 0.)
    if state is None or state > heat:
        return None
    lo, hi = math.log(1e-9), math.log(.3)
    if (first_state(q, p, .3) or math.inf) <= heat:
        return .3
    for _ in range(60):
        mid = .5 * (lo + hi)
        s = first_state(q, p, math.exp(mid))
        lo, hi = (mid, hi) if s is not None and s <= heat else (lo, mid)
    return math.exp(.5 * (lo + hi))


def counted_heats(case, film, activity, counts, glow):
    """The fixed heats as the escape model counts them: the loss response's own (the window film behind the titania
    stack, nothing behind the 200-nm edge, and the glow), and its rule applied to this module's films."""
    shield, _ = CASES[case]
    g = glow * loss.ACTIVITY[activity]['glow']
    o1 = (counts['window_film_heat'][activity] if shield == 'titania_stack' else 0.) + g
    rule = counts['window_film_heat'][activity] + counts[FILMS[film]][activity] + g
    return dict(o1_count=o1, escape_rule=rule)


def state_row(cfg, low, q, p, gaps):
    """A scenario's traced state, its parts and its loss, with the loss if the heat lay as low as the column's 'low'
    shape puts it."""
    state = first_state(q, p, gaps)
    if state is None:
        return dict(status='runaway_beyond_profile_heats', beyond_W_m2=float(q[-1]))
    shares = {k: float(np.interp(state, q, v)) for k, v in p.items() if k != 'gap_unit'}
    shares['gaps'] = float(gaps * np.interp(state, q, p['gap_unit']))
    bound = molecular_loss(low, state)
    kept = ('status', 'exobase_temperature_k', 'molecular_loss_kg_s')
    return dict(molecular_loss(cfg, state), parts_W_m2=shares, low_shape={k: bound[k] for k in kept if k in bound})


def analyse(case, entries, counts, cycle):
    """Scenario states, their losses, and the gap transmission each budget allows, traced and counted; at solar
    maximum also with
    each source's X-rays at FISM2's factor for each of the last three maxima."""
    cfg = case_config(case)
    low = dataclasses.replace(cfg, heating_shape='low')
    glow = glow_heat(case)
    unit = counts['band_over_disk']
    scenarios, budgets = {}, []
    for film in FILMS:
        for activity in loss.ACTIVITY:
            q, p = parts(entries, film, activity, glow)
            fixed = counted_heats(case, film, activity, counts, glow)
            periods = list(cycle) if activity == 'solar_maximum' else []
            scaled = {k: parts(entries, film, activity, glow, xray=cycle[k]['sources'])[1] for k in periods}
            for gaps in [0.] + DESIGN['gap_transmissions']:
                row = dict(film=film, activity=activity, gap_transmission=gaps, traced=state_row(cfg, low, q, p, gaps))
                for count, heat in fixed.items():
                    row[count] = molecular_loss(cfg, heat + gaps * unit[activity])
                scenarios[f'{film}_{activity}_gaps_{gaps:g}'] = row
                for period, scaled_parts in scaled.items():
                    scenarios[f'{film}_{activity}_fism2_{period}_gaps_{gaps:g}'] = dict(
                        film=film, activity=activity, gap_transmission=gaps, xrays='fism2_' + period,
                        traced=state_row(cfg, low, q, scaled_parts, gaps))
            flat = parts(entries, film, activity, glow, gap_source='open_flat')[1]
            for budget in DESIGN['budgets_kg_s']:
                heat = heat_for_loss(cfg, budget, float(q[-1]))
                row = dict(film=film, activity=activity, budget_kg_s=budget, heat_W_m2=heat)
                if heat is not None:
                    for count, fixed_heat in fixed.items():
                        f = (heat - fixed_heat) / unit[activity]
                        row[count] = f if f > 0 else None
                    row['traced'] = allowed_traced(q, p, heat)
                    if row['traced'] and row['o1_count']:
                        row['traced_over_o1_count'] = row['traced'] / row['o1_count']
                    row['traced_fism2'] = {k: allowed_traced(q, scaled[k], heat) for k in periods}
                    row['gap_heat_over_disk_count'] = float(np.interp(heat, q, p['gap_unit']) / unit[activity])
                    row['gap_heat_over_disk_count_xrays_at_uv_factor'] = float(np.interp(heat, q, flat['gap_unit'])
                                                                               / unit[activity])
                budgets.append(row)
    return dict(glow_heat_quiet_W_m2=glow, scenarios=scenarios, budgets=budgets)


TABLE_CODE = (cross_sections, transmissions, case_row, profile_path, middle_profile, case_config, atmosphere,
              escape_covering, aperture_edge, impact_parameters, ring_areas, zone_weights, log_mean, deposit, spectra,
              heat_table, summarise, profile_by_pressure, tables)


def table_key(code, inputs, products, constants):
    """Digest of everything the heat tables depend on: the functions that make them, their settings, the other
    code they call, the shared constants, the inputs and the products they read."""
    source = [inspect.getsource(f) for f in TABLE_CODE]
    settings = {k: DESIGN[k] for k in TABLE_SETTINGS}
    other = {k: v for k, v in code.items() if not k.endswith('limb_heat.py')}
    blob = json.dumps([source, settings, CASES, MOLECULE, MASS_KG, SOURCES, other, constants, inputs, products],
                      sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def provenance():
    """Digests of the code, inputs and products, the shared constants read, and the heat tables' key."""
    code = {name: digest(path) for name, path in FILES.items()}
    inputs = {n: digest(fetch_limb_inputs.path(n)) for n in LIMB_INPUTS}
    inputs[WHI] = digest(fetch_inputs.path(WHI))
    cycle_inputs = {n: digest(fetch_limb_inputs.path(n)) for n in [CYCLE_QUIET, *CYCLE_MAXIMA.values()]}
    products = {'atmosphere/middle_atmosphere/results/middle_atmosphere.csv': digest(RESULTS / 'middle_atmosphere.csv'),
                'atmosphere/loss_response/results/tidal_escape.json': digest(loss.TIDES),
                'protection/spectra/stack_short_wave.json': digest(escape.FILM_PRODUCT)}
    for case in CASES:
        products[str(profile_path(case).relative_to(ROOT))] = digest(profile_path(case))
    constants = constants_used(list(FILES.values()))
    return code, inputs, cycle_inputs, products, constants, table_key(code, inputs, products, constants)


def main(argv=None) -> int:
    whi_w, whi_f = escape.whi_quiet_sun()
    sel = (whi_w > 0.1) & (whi_w < DESIGN['band_nm'])
    w, f = whi_w[sel], whi_f[sel]
    films = transmissions(w)
    sigma = cross_sections(w)
    names, weights = spectra(w, f, films)
    b = impact_parameters()
    rays = (b, ring_areas(b), zone_weights(b))
    counts = escape_count(w, f, films)
    code, inputs, cycle_inputs, products, constants, key = provenance()
    cycle = cycle_factors()
    cases = {}
    for case in CASES:
        entries = tables(case, rays, sigma, names, weights, key)
        result = analyse(case, entries, counts, cycle)
        cases[case] = dict(by_profile_heat=entries, **result)
        print(case, flush=True)
        for row in result['budgets']:
            print('   ', {k: (f'{v:.3g}' if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    product = dict(
        schema=SCHEMA, producer=dict(domain='atmosphere', files=code, constants=constants,
                                     inputs=dict(inputs, **cycle_inputs), products=products),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=(
            'Heats are W per m^2 of lunar surface, global means at the escape model\'s heating efficiency, as the '
            'thermal column takes them; "thermosphere" is the heat between the 0.3 Pa base and the exobase, which the '
            'column takes, and "exosphere_absorbed" the light absorbed above the exobase, which it does not. '
            'Transmissions through gaps are grey shares of the band below 175 nm. A traced state puts the window '
            'stack on the window, a film on the annulus, the gaps\' light over the whole aperture, unfiltered sunlight '
            'beyond it and the sky\'s glow on '
            'the air, and is the first heat from zero that the air it swells returns. "o1_count" is the loss '
            'response\'s count, from which O1\'s allowed transmissions were derived; "escape_rule" applies that '
            'count to this module\'s films; "gap_heat_over_disk_count_xrays_at_uv_factor" scales the gaps\' X-rays '
            'as the loss response does. "fism2" rows and "traced_fism2" give each source\'s X-rays below 10 nm '
            'the factor FISM2 finds between the WHI quiet week and the year around each of the last three solar '
            'maxima ("xray_cycle"), in place of the convention\'s 100. "low_shape" is the loss at the same heat with '
            'the column\'s low heating shape. "window_stack" puts the climate window\'s stack on the annulus too. A '
            'budget row\'s transmissions are the largest grey transmissions through gaps whose molecular loss with '
            'Earth\'s tide stays within '
            'the budget, the solar wind and the exosphere step\'s losses not deducted.'),
        design=DESIGN, escape_count=counts, xray_cycle=cycle, cases=cases)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
