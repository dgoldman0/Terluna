"""The rain product: fall speed and survival by drop size in the Moon's air at 2 m in the three target climates,
beside Earth's air, and the check of CM1's rain fall speeds against the single-drop solver.

    python -m atmosphere.rain.product        # about three and a half minutes on one thread
    # -> atmosphere/rain/results/rain_drops.json

The air at 2 m is each target's median over land from the climate targets' product
(research/studies/climate_targets): its 2 m air temperature and 2 m wet-bulb temperature, from which the specific
humidity follows by inverting that study's own psychrometric relation. The Moon's air is the repository's design air,
1.2 atm with Earth's oxygen partial pressure and 400 ppm of CO2 (atmosphere/radiative_convective/thermodynamics.py,
earthlike_air), under lunar gravity. Earth's air is Earth's dry air at 1 atm under standard gravity, at the same
temperature and relative humidity. A third air, the Moon's under Earth's gravity, isolates gravity for the CM1 check.

Survival starts each drop at the parcel cloud base at its terminal speed, in still air, and follows it to the ground;
the rain work's bounded cloud-base spectra weigh the classes into the share of the water that arrives. The CM1 check
sets the solver's slowing by gravity alone beside the one factor CM1's Morrison scheme applies at every size,
(g / 9.81)^((B + 1) / 3) with B = 0.8, and sets CM1's whole rain law in the Moon's air beside the solver.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from shared.constants import (CM1_DRY_AIR_GAS_CONSTANT, CM1_GRAVITY, CM1_RAIN_DENSITY_EXPONENT,
                              CM1_RAIN_DENSITY_REFERENCE_PRESSURE, CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE,
                              CM1_RAIN_FALL_CAP, CM1_RAIN_FALL_COEFFICIENT, CM1_RAIN_FALL_EXPONENT,
                              MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)
from shared.provenance import constants_used
from atmosphere.radiative_convective.thermodynamics import DryAir, earthlike_air
from atmosphere.rain import distributions, drops

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.rain-drops/1'
OUT = HERE / 'results' / 'rain_drops.json'
TARGETS_PRODUCT = ROOT / 'research' / 'studies' / 'climate_targets' / 'results' / 'climate_targets.json'
TARGETS_SCHEMA = 'terluna.research.climate-targets/2'
CODE = ('atmosphere/rain/drops.py', 'atmosphere/rain/distributions.py', 'atmosphere/rain/product.py',
        'atmosphere/radiative_convective/thermodynamics.py')
# Each target as this product names it, and its name in the climate targets' product.
TARGETS = (('hot_world', 'no further dimming'), ('resort', '300 K'), ('lock', '295 K'))
ATMOSPHERE_PA = 101325.0
MOON_PRESSURE_PA = 1.2 * ATMOSPHERE_PA       # the design pressure at the water level (research/decisions.md)
LAND_MEAN_PRESSURE_PA = 1.16 * ATMOSPHERE_PA  # the land's mean surface pressure in the GCM (climate/gcm/README.md)
CO2_PPM = 400.0
SPEED_MM = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)
CLOUD_BASE_MM = (0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0)
# The bounded cloud-base spectra (the rain work's scenario set of 11 October): gamma concentrations with mass-weighted
# mean diameters of 0.75, 1.5 and 2.5 mm and shapes 0 and 3, on 64 classes from 0.1 to 6 mm, entering at 5 mm/h. With
# independent drops the share reaching the ground is the same at any rain rate, and the encounters grow with it.
SPECTRUM_GRID_MM = (0.1, 6.0, 64)
SPECTRUM_LEVELS = 33
SPECTRA = ((0.75, 0.0), (0.75, 3.0), (1.5, 0.0), (1.5, 3.0), (2.5, 0.0), (2.5, 3.0))
CENTRAL_SPECTRUM = (1.5, 3.0)
RAIN_MM_H = 5.0
# The climate targets' wet-bulb temperature (research/studies/climate_targets/extract.py, wet_bulb): Bolton's (1980)
# saturation over water and the psychrometric balance with the GCM's heat capacity of dry air.
GCM_HEAT_CAPACITY = 1004.64

EVIDENCE = ('A runnable calculation of independent single drops: the solver reproduces 34 Gunn-Kinzer fall speeds '
            'on Earth to 0.20 m/s (root mean square) and 7.4% at most, and its lunar values are that physics carried '
            'to lunar gravity and air. The air at 2 m pairs each target\'s median land temperature with its median '
            'land wet-bulb temperature. Survival holds the parcel column below the cloud fixed and the drops '
            'independent, so it filters sizes under that air; the spectra are bounded scenario boundary '
            'conditions. The CM1 check compares fall-speed laws; CM1 was not rerun.')
READING_RULE = ('Diameters are volume-equivalent, in mm. fall_speed gives the still-air terminal speed at 2 m in m/s, '
                'null above the stability screen. survival starts each drop at the cloud base and gives what '
                'reaches the ground; its rain section gives the share of the water leaving the cloud base that '
                'reaches the ground for each spectrum, the same at any rain rate, and the encounters '
                'along the central spectrum\'s paths at 5 mm/h, which grow in proportion to it. cm1_check gives, by '
                'diameter, the solver\'s ratio of Earth-gravity to lunar-gravity speed in the Moon\'s air beside '
                'CM1\'s single factor, CM1\'s whole lunar law in that air over the solver\'s speed, and the water '
                'lost per metre of fall at CM1\'s speed over that at the solver\'s, with the same ventilation. '
                'Interpolate in log diameter.')
ASSUMPTIONS = [
    'Each target\'s air at 2 m is its median land temperature at 2 m with its median land wet-bulb temperature at 2 m, '
    'paired; the specific humidity follows from the targets\' psychrometric relation at 1.2 atm, where the GCM '
    'computed it at each cell\'s own surface pressure, which moves the vapour pressure by 0.6%.',
    'The Moon\'s air is at 1.2 atm, the design pressure at the water level; the land\'s mean surface pressure in the '
    'GCM is 1.16 atm, and pressure_sensitivity gives its effect on the fall speeds.',
    'The column below the cloud cools at the dry adiabat g/cp and holds the 2 m vapour mole fraction, so its cloud '
    'base is the parcel saturation height; the rain neither cools nor moistens it, and the air is still.',
    'Earth\'s air carries the same temperature and relative humidity at 2 m as the Moon\'s, at 1 atm.',
    'Drops are pure liquid water in their equilibrium shape, each falling independently of the others.',
]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _bolton(t):
    return 611.2 * math.exp(17.67 * (t - 273.15) / (t - 29.65))


def humidity_from_wet_bulb(air_k, wet_bulb_k, pressure_pa):
    """Specific humidity (kg/kg) and vapour pressure (Pa) whose psychrometric wet-bulb temperature, as the climate
    targets compute it, is wet_bulb_k: the inverse of qsat(Tw) - q = (cp + 1850 q) (T - Tw) / L(Tw)."""
    e_w = _bolton(wet_bulb_k)
    q_sat = 0.622 * e_w / (pressure_pa - 0.378 * e_w)
    latent = 2.501e6 - 2370.0 * (wet_bulb_k - 273.15)
    depression = air_k - wet_bulb_k
    q = (q_sat - GCM_HEAT_CAPACITY * depression / latent) / (1.0 + 1850.0 * depression / latent)
    return q, q * pressure_pa / (0.622 + 0.378 * q)


def wet_bulb(air_k, q, pressure_pa):
    """The climate targets' wet-bulb temperature (K) by bisection, for checking the inverse."""
    lo, hi = 200.0, air_k
    for _ in range(60):
        tw = 0.5 * (lo + hi)
        e = _bolton(tw)
        f = 0.622 * e / (pressure_pa - 0.378 * e) - q - (GCM_HEAT_CAPACITY + 1850.0 * q) / (
            2.501e6 - 2370.0 * (tw - 273.15)) * (air_k - tw)
        lo, hi = (tw, hi) if f < 0 else (lo, tw)
    return 0.5 * (lo + hi)


def target_air(targets):
    """The 2 m air of each target over land: temperature, wet bulb, humidity at 1.2 atm."""
    if targets.get('schema') != TARGETS_SCHEMA:
        raise ValueError(f'Expected {TARGETS_SCHEMA}, found {targets.get("schema")}')
    by_name = {t['name']: t for t in targets['targets']}
    out = {}
    for key, name in TARGETS:
        climate = by_name[name]['climate']
        t = climate['land_air_2m_median_c'] + 273.15
        tw = climate['wet_bulb_2m_median_c'] + 273.15
        q, e = humidity_from_wet_bulb(t, tw, MOON_PRESSURE_PA)
        out[key] = dict(name_in_targets=name, ground_k=by_name[name]['ground_k'], air_2m_k=t, wet_bulb_2m_k=tw,
                        specific_humidity=q, vapour_pressure_pa=e,
                        relative_humidity=e / drops.liquid_saturation_pressure(t))
    return out


def columns(air):
    """The three airs at a target's 2 m state: the Moon's, Earth's, and the Moon's under Earth's gravity."""
    t, rh = air['air_2m_k'], air['relative_humidity']
    moon_air, earth_air = earthlike_air(MOON_PRESSURE_PA, CO2_PPM), DryAir(ATMOSPHERE_PA)
    moon = drops.Column(gravity_m_s2=MOON_SURFACE_GRAVITY, surface_pressure_pa=MOON_PRESSURE_PA,
                        surface_temperature_k=t, surface_relative_humidity=rh, air=moon_air)
    earth = drops.Column(gravity_m_s2=STANDARD_GRAVITY, surface_pressure_pa=ATMOSPHERE_PA,
                         surface_temperature_k=t, surface_relative_humidity=rh, air=earth_air)
    moon_earth_g = replace(moon, gravity_m_s2=STANDARD_GRAVITY, lapse_k_m=None, cloud_base_m=None)
    return dict(moon=moon, earth=earth, moon_air_earth_gravity=moon_earth_g)


def speeds(column):
    """Terminal speed and axis ratio at 2 m by diameter, null above the stability screen."""
    state = column.state(0.0)
    screen = drops.stable_diameter(column)
    v, ratio = [], []
    for mm in SPEED_MM:
        d = mm * 1e-3
        stable = d <= screen
        v.append(drops.terminal_speed(d, column) if stable else None)
        ratio.append(drops.axis_ratio(d, column.gravity_m_s2, state['density_kg_m3'], state['temperature_k'])
                     if stable else None)
    return dict(terminal_speed_m_s=v, axis_ratio=ratio, stable_diameter_mm=screen * 1e3)


def _survives(mm, column):
    return drops.fall(mm * 1e-3, column, sample_heights_m=[column.cloud_base_m, 0.0])['survived']


def smallest_surviving(column, precision=1.005):
    """The smallest drop leaving the cloud base that reaches the ground (mm), by bisection in log diameter."""
    lo, hi = 0.02, 6.0
    if _survives(lo, column) or not _survives(hi, column):
        raise ValueError('The survival threshold lies outside 0.02-6 mm')
    while hi / lo > precision:
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if _survives(mid, column) else (mid, hi)
    return math.sqrt(lo * hi)


def spectra(column, bins=SPECTRUM_GRID_MM[2]):
    """The share of the water leaving the cloud base that reaches the ground, for the bounded gamma spectra, and the
    collision warning along the paths of the central one."""
    d, width = distributions.diameter_grid(SPECTRUM_GRID_MM[0], SPECTRUM_GRID_MM[1], bins)
    heights = np.linspace(column.cloud_base_m, 0.0, SPECTRUM_LEVELS)
    runs = [drops.fall(x, column, sample_heights_m=heights) for x in d]
    entry = np.array([r['initial_speed_m_s'] for r in runs])
    alive = np.array([r['survived'] for r in runs])
    surface_d = np.array([r['surface_diameter_m'] if r['survived'] else np.nan for r in runs])
    surface_v = np.array([r['surface_speed_m_s'] if r['survived'] else np.nan for r in runs])
    path = {k: np.array([[np.nan if x is None else x for x in r['trajectory'][k]] for r in runs]).T
            for k in ('diameter_m', 'speed_m_s', 'terminal_speed_m_s')}
    rows, central = [], None
    for dm, mu in SPECTRA:
        source = distributions.gamma_flux(d, width, entry, dm_mm=dm, shape_mu=mu, cloud_base_rain_mm_h=RAIN_MM_H)
        ground = distributions.surface_summary(d, source['flux_m2_s'], surface_d, surface_v, alive)
        rows.append(dict(dm_mm=dm, shape_mu=mu, realized_dm_mm=source['realized_dm_mm'],
                         water_survival_fraction=ground['water_survival_fraction'],
                         number_survival_fraction=ground['number_survival_fraction'],
                         arrival_diameter_mm_by_number=(ground['arrival_diameter_mm_by_number'] or {}).get('p50'),
                         arrival_speed_m_s_by_number=(ground['arrival_speed_m_s_by_number'] or {}).get('p50')))
        if (dm, mu) == CENTRAL_SPECTRUM:
            meet = distributions.encounters(heights, path['diameter_m'], path['speed_m_s'],
                                            path['terminal_speed_m_s'], source['flux_m2_s'])
            central = dict(dm_mm=dm, shape_mu=mu, rain_mm_h=RAIN_MM_H,
                           **{k: meet[k] for k in ('tau_by_surface_water', 'volume_exposure_by_surface_water',
                                                   'surface_water_share_tau_above_1')})
    return dict(grid_mm=[SPECTRUM_GRID_MM[0], SPECTRUM_GRID_MM[1], bins], levels=SPECTRUM_LEVELS, spectra=rows,
                central_encounters=central)


def grid_check(column):
    """The spectra's shares on twice the classes, which resolve the survival threshold twice as finely."""
    fine = spectra(column, bins=2 * SPECTRUM_GRID_MM[2])
    return dict(classes=2 * SPECTRUM_GRID_MM[2],
                spectra=[dict(dm_mm=r['dm_mm'], shape_mu=r['shape_mu'],
                              water_survival_fraction=r['water_survival_fraction']) for r in fine['spectra']])


def survival(column):
    """What reaches the ground from each size leaving the parcel cloud base."""
    rows = []
    for mm in CLOUD_BASE_MM:
        if mm * 1e-3 > drops.stable_diameter(column, column.cloud_base_m):
            rows.append(dict(initial_diameter_mm=mm, status='above the stability screen'))
            continue
        r = drops.fall(mm * 1e-3, column, sample_heights_m=[column.cloud_base_m, 0.0])
        rows.append(dict(initial_diameter_mm=mm, status=r['status'],
                         surviving_water_fraction=1.0 - r['evaporated_mass_fraction'] if r['survived'] else 0.0,
                         surface_diameter_mm=r['surface_diameter_m'] * 1e3 if r['survived'] else None,
                         impact_speed_m_s=r['surface_speed_m_s'] if r['survived'] else None,
                         fall_time_s=r['fall_time_s']))
    base = column.state(column.cloud_base_m)
    return dict(cloud_base_m=column.cloud_base_m, cloud_base_k=base['temperature_k'],
                cloud_base_pa=base['pressure_pa'], lapse_k_km=column.lapse_k_m * 1e3,
                smallest_surviving_mm=smallest_surviving(column), drops=rows, rain=spectra(column))


def cm1_density_factor(density):
    """Morrison's (RHOSU / rho)^0.54, with RHOSU the air at 850 hPa and 0 C."""
    reference = (CM1_RAIN_DENSITY_REFERENCE_PRESSURE
                 / (CM1_DRY_AIR_GAS_CONSTANT * CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE))
    return (reference / density)**CM1_RAIN_DENSITY_EXPONENT


def cm1_gravity_factor(gravity):
    """The Terluna builds' scaling of Morrison's rain law and cap, (g / 9.81)^((B + 1) / 3)."""
    return (gravity / CM1_GRAVITY)**((CM1_RAIN_FALL_EXPONENT + 1) / 3)


def cm1_rain_speed(diameter_m, density, gravity):
    """Morrison's raindrop fall speed (m/s) at a diameter, as CM1's build for that gravity sets it."""
    return (CM1_RAIN_FALL_COEFFICIENT * cm1_gravity_factor(gravity) * cm1_density_factor(density)
            * diameter_m**CM1_RAIN_FALL_EXPONENT)


def per_metre_loss_ratio(diameter_m, speed_a, speed_b, state):
    """Water lost per metre of fall at speed_a over that at speed_b: f(Re_a) / v_a over f(Re_b) / v_b, with the same
    ventilation law and air."""
    schmidt = state['viscosity_pa_s'] / (state['vapour_diffusivity_m2_s'] * state['density_kg_m3'])

    def per_metre(v):
        reynolds = state['density_kg_m3'] * v * diameter_m / state['viscosity_pa_s']
        return drops.ventilation(reynolds, schmidt) / v

    return per_metre(speed_a) / per_metre(speed_b)


def cm1_check(cols):
    """The solver's slowing by gravity beside CM1's factor, and CM1's whole lunar law in the Moon's air at 2 m."""
    moon, earth_g, earth = cols['moon'], cols['moon_air_earth_gravity'], cols['earth']
    state, earth_state = moon.state(0.0), earth.state(0.0)
    factor = 1.0 / cm1_gravity_factor(MOON_SURFACE_GRAVITY)

    def slowing(mm):
        d = mm * 1e-3
        return drops.terminal_speed(d, earth_g) / drops.terminal_speed(d, moon)

    rows = dict(diameter_mm=[], slowing_by_gravity=[], rigid_sphere_slowing_by_gravity=[], solver_speed_m_s=[],
                cm1_speed_m_s=[], cm1_over_solver=[], earth_law_over_solver_on_earth=[],
                water_lost_per_metre_cm1_over_solver=[])
    for mm in SPEED_MM:
        d = mm * 1e-3
        if d > drops.stable_diameter(earth_g):
            continue
        v_moon = drops.terminal_speed(d, moon)
        v_cm1 = cm1_rain_speed(d, state['density_kg_m3'], MOON_SURFACE_GRAVITY)
        rows['diameter_mm'].append(mm)
        rows['slowing_by_gravity'].append(slowing(mm))
        rows['rigid_sphere_slowing_by_gravity'].append(drops.terminal_speed(d, earth_g, deform=False)
                                                       / drops.terminal_speed(d, moon, deform=False))
        rows['solver_speed_m_s'].append(v_moon)
        rows['cm1_speed_m_s'].append(v_cm1)
        rows['cm1_over_solver'].append(v_cm1 / v_moon)
        rows['earth_law_over_solver_on_earth'].append(
            cm1_rain_speed(d, earth_state['density_kg_m3'], CM1_GRAVITY) / drops.terminal_speed(d, earth))
        rows['water_lost_per_metre_cm1_over_solver'].append(per_metre_loss_ratio(d, v_cm1, v_moon, state))
    crossing = brentq(lambda mm: slowing(mm) - factor, 0.5, 3.0, xtol=1e-4)
    equal_speed = brentq(lambda mm: cm1_rain_speed(mm * 1e-3, state['density_kg_m3'], MOON_SURFACE_GRAVITY)
                         - drops.terminal_speed(mm * 1e-3, moon), 0.2, 6.0, xtol=1e-4)
    cap = CM1_RAIN_FALL_CAP * cm1_gravity_factor(MOON_SURFACE_GRAVITY) * cm1_density_factor(state['density_kg_m3'])
    at_cap = brentq(lambda mm: drops.terminal_speed(mm * 1e-3, moon) - cap, 1.0, 8.0, xtol=1e-4)
    return dict(cm1_slowing_at_every_size=factor, crossing_diameter_mm=crossing, equal_speed_diameter_mm=equal_speed,
                cm1_mean_fall_speed_cap_m_s=cap, solver_diameter_at_the_cap_mm=at_cap, **rows)


def pressure_sensitivity(air):
    """Fall speeds at the land's mean surface pressure over those at 1.2 atm, at a target's 2 m air (the same
    composition, which is uniform with height)."""
    composition = earthlike_air(MOON_PRESSURE_PA, CO2_PPM)

    def at(pressure):
        return drops.Column(gravity_m_s2=MOON_SURFACE_GRAVITY, surface_pressure_pa=pressure,
                            surface_temperature_k=air['air_2m_k'], surface_relative_humidity=air['relative_humidity'],
                            air=composition)

    high, low = at(MOON_PRESSURE_PA), at(LAND_MEAN_PRESSURE_PA)
    return dict(pressure_pa=LAND_MEAN_PRESSURE_PA, diameter_mm=list(SPEED_MM[:12]),
                speed_ratio=[drops.terminal_speed(mm * 1e-3, low) / drops.terminal_speed(mm * 1e-3, high)
                             for mm in SPEED_MM[:12]])


def describe(column):
    state = column.state(0.0)
    return dict(gravity_m_s2=column.gravity_m_s2, pressure_pa=column.surface_pressure_pa,
                temperature_k=column.surface_temperature_k, relative_humidity=column.surface_relative_humidity,
                dry_fractions=dict(column.air.fractions), dry_molar_mass_kg_mol=column.air.molar_mass,
                density_kg_m3=state['density_kg_m3'], viscosity_pa_s=state['viscosity_pa_s'],
                vapour_diffusivity_m2_s=state['vapour_diffusivity_m2_s'])


def build():
    targets = json.loads(TARGETS_PRODUCT.read_text())
    air = target_air(targets)
    product = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS,
                   targets=air, diameters_mm=list(SPEED_MM), cloud_base_diameters_mm=list(CLOUD_BASE_MM),
                   airs={}, fall_speed={}, survival={}, cm1_check={})
    for key, _ in TARGETS:
        cols = columns(air[key])
        product['airs'][key] = {world: describe(c) for world, c in cols.items()}
        product['fall_speed'][key] = {world: speeds(c) for world, c in cols.items()}
        product['survival'][key] = {world: survival(cols[world]) for world in ('moon', 'earth')}
        product['cm1_check'][key] = cm1_check(cols)
    product['pressure_sensitivity'] = dict(resort=pressure_sensitivity(air['resort']))
    resort = columns(air['resort'])
    product['grid_check'] = dict(resort={world: grid_check(resort[world]) for world in ('moon', 'earth')})
    product['producer'] = dict(domain='atmosphere', files={p: sha256(ROOT / p) for p in CODE},
                               inputs={str(TARGETS_PRODUCT.relative_to(ROOT)): sha256(TARGETS_PRODUCT)},
                               constants=constants_used([ROOT / p for p in CODE]))
    return product


def main(argv=None) -> int:
    product = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    for key, _ in TARGETS:
        moon, check, air = product['survival'][key]['moon'], product['cm1_check'][key], product['targets'][key]
        print(f"{key:9s} 2 m {air['air_2m_k'] - 273.15:5.2f} C RH {air['relative_humidity']:.3f}; "
              f"cloud base {moon['cloud_base_m'] / 1e3:5.2f} km; "
              f"smallest surviving {moon['smallest_surviving_mm']:.2f} mm; "
              f"CM1 {check['cm1_slowing_at_every_size']:.2f} crosses the solver at "
              f"{check['crossing_diameter_mm']:.2f} mm")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
