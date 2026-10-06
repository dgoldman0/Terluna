"""Integrated ledger of the shield candidates (step 2 of the integrated comparison).

The author's criteria of 6 October 2026 (research/decisions.md): the shield system meets the shield
requirements, gives reasonable net-positive electricity for habitat, balances overall resource requirements,
and is stable and safe, with every mass flow judged by direction as well as size. This runner gathers what
the products establish for two candidates, adds the arithmetic that follows from stated assumptions, and
marks what no model yet supplies. It runs no new dynamics.

- held_zoned: the held screen with a zoned aperture, 26 g/m2 in the window and 5 g/m2 in the annulus, on the
  selected moving trajectory (zoned_aperture.json, exhaust_isolation.json).
- ring_fleet: a global screen of nested ring bundles in lunar orbit, kept with photon forces
  (ring_bundle.json, ring_keeping.json, ring_screen.json). Its inventory follows from the bundle's geometry at
  a ring radius of 15,000 or 19,000 km; no full-fleet run exists.

    python -m research.studies.solar_shield_array.integrated_ledger
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.optical import covering_radius
from .run import CORE_DIVERTED, SCENARIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/integrated_ledger.json'
FILES = ['research/studies/solar_shield_array/integrated_ledger.py', 'research/studies/solar_shield_array/run.py',
         'protection/dynamics/optical.py']
PRODUCTS = {name: HERE/f'results/{name}.json' for name in
            ('zoned_aperture', 'exhaust_isolation', 'ring_bundle', 'ring_keeping', 'ring_screen')}
# Values fixed by the requirements (research/studies/protection_architecture/requirements.md: R1, S6, O8, E2)
# and the tiles of every study (10 km squares with 9.89 km clear sides).
REQUIREMENTS = dict(loss_budgets_kg_s=[1., 10., 100.], unshielded_loss_kg_s=[300., 80000.],
                    replacement_years=[10., 20.], climate_top_of_air_W_m2=1173.25, dimmer_share=.05,
                    tile_side_m=10e3, clear_side_m=9.89e3, protected_radii=4)
# Stated assumptions: the design brief's illustrative conversion range, demand cases around the 100 TW reference,
# and the two ring radii the plan names.
ASSUMPTIONS = dict(conversion=[.2, .3], demand_TW=[100., 300., 1000.], ring_radii_km=[15000., 19000., 20000.])
STATES = ('met', 'conditional', 'open', 'not met')
YEAR = K.JULIAN_YEAR_DAYS*K.JULIAN_DAY


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def held(zoned):
    prop = SCENARIO['propulsion']
    cant = np.cos(np.radians(prop['cant_deg']))
    radius, window = zoned['aperture_radius_km']*1e3, zoned['core_radius_km']*1e3
    area = np.pi*radius**2
    flow = zoned['propellant_kg_s']
    return dict(aperture_radius_km=radius/1e3, window_radius_km=window/1e3, area_km2=area/1e6,
                panel_equivalents=dict(physical=area/REQUIREMENTS['tile_side_m']**2,
                                       clear=area/REQUIREMENTS['clear_side_m']**2),
                optical_mass_Gt=zoned['base_optical_mass_kg']/1e12, held_mass_Gt=zoned['total_mass_kg']/1e12,
                power_hardware_Gt=zoned['power_hardware_mass_kg']/1e12,
                holding_power_mean_TW=zoned['mean_power_TW'], installed_TW=zoned['installed_design_peak_power_TW'],
                jet_power_TW=prop['efficiency']*zoned['mean_power_TW'],
                propellant_kg_s=flow, propellant_Gt_per_year=flow*YEAR/1e12,
                propellant_per_held_kg_per_year=flow*YEAR/zoned['total_mass_kg'],
                propellant_billion_years_over_moon_mass=flow*1e9*YEAR/(K.MOON_GM/K.GRAVITATIONAL_CONSTANT),
                mean_propulsive_acceleration_m_s2=flow*prop['exhaust_velocity_m_s']*cant/zoned['total_mass_kg'],
                film_replacement_Gt_per_year=[zoned['base_optical_mass_kg']/y/1e12
                                              for y in REQUIREMENTS['replacement_years']],
                unpowered_drift_first_day_km=(.5*flow*prop['exhaust_velocity_m_s']*cant/zoned['total_mass_kg']
                                              *K.JULIAN_DAY**2/1e3),
                s6_propellant_over_unshielded_loss=[flow/x for x in REQUIREMENTS['unshielded_loss_kg_s']])


def ring_fleet(bundle, keeping, screen, zoned):
    """Inventory, mass and control of a full ring screen, from the bundle's geometry."""
    pitch, height = bundle['design']['pitch_m'], bundle['built']['height_pitch_m']
    core, annulus = zoned['core_sigma_kg_m2'], zoned['annulus_sigma_kg_m2']
    held_panels = np.pi*(zoned['aperture_radius_km']*1e3)**2/REQUIREMENTS['tile_side_m']**2
    peak = keeping['keeping']['peak_acceleration_m_s2']
    per_orbit = keeping['keeping']['rings_mean_m_s_per_orbit']
    ve = SCENARIO['propulsion']['exhaust_velocity_m_s']
    rows = []
    for r_km in ASSUMPTIONS['ring_radii_km']:
        r = r_km*1e3
        aperture = float(covering_radius(r, K.AU, REQUIREMENTS['protected_radii']*K.MOON_RADIUS)
                         + SCENARIO['formation_margin_m'])
        window = float(covering_radius(r, K.AU, K.MOON_RADIUS))
        rings = int(np.ceil(2*aperture/height))
        per_ring = int(np.ceil(2*np.pi*r/pitch))
        tiles = rings*per_ring
        # Every tile of a ring that crosses the window carries the window's areal mass.
        share = window/aperture
        sigma = share*core+(1-share)*annulus
        area = tiles*REQUIREMENTS['tile_side_m']**2
        mass = area*sigma
        cases = [c for c in screen['cases'] if np.isclose(c['radius_km'], r_km)]
        steering = [c['steering_m_s_per_day'] for c in cases]
        period = 2*np.pi*np.sqrt(r**3/K.MOON_GM)
        orbits_per_year = YEAR/period
        rows.append(dict(
            radius_km=r_km, aperture_radius_km=aperture/1e3, window_radius_km=window/1e3, rings=rings,
            tiles_per_ring=per_ring, tiles=tiles, tiles_over_held_panels=tiles/held_panels,
            area_over_aperture=area/(np.pi*aperture**2), window_ring_share=share, mean_sigma_g_m2=sigma*1e3,
            optical_mass_Gt=mass/1e12,
            film_replacement_Gt_per_year=[mass/y/1e12 for y in REQUIREMENTS['replacement_years']],
            steering_m_s_per_day=[min(steering), max(steering)] if steering else None,
            crossing_height_error_max_km=([min(c['height_error_max_km'] for c in cases),
                                           max(c['height_error_max_km'] for c in cases)] if cases else None),
            # For comparison only: the same control by electric thrust would release this exhaust in lunar orbit.
            electric_alternative_kg_s=dict(
                keeping=mass*per_orbit*orbits_per_year/ve/YEAR,
                steering=[mass*s*K.JULIAN_YEAR_DAYS/ve/YEAR for s in (min(steering), max(steering))]
                if steering else None)))
    sail = {f'{s*1e3:g}_g_m2': CORE_DIVERTED*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*s) for s in (annulus, core)}
    return dict(by_radius=rows,
                electric_alternative_note=('Exhaust if the kept run\'s mean keeping (mostly the fixed keeping frame '
                                           'working against the common regression) and the year-long screen\'s '
                                           'steering of the steady drift were done by 30 km/s electric thrust.'),
                photon_keeping=dict(
        sail_acceleration_m_s2=sail, keeping_peak_m_s2=peak,
        keeping_peak_over_sail={k: peak/v for k, v in sail.items()},
        steering_m_s2={f"{row['radius_km']:g}_km": [s/K.JULIAN_DAY for s in row['steering_m_s_per_day']]
                       for row in rows if row['steering_m_s_per_day']},
        note=('The redirected band gives the sail acceleration shown; tilting the tile turns part of it along-track '
              'or normal. Whether its direction can follow the demand around each orbit is the open check.')))


def electricity(held_row, zoned):
    s = K.SOLAR_CONSTANT
    prop = SCENARIO['propulsion']
    thrust_speed = prop['exhaust_velocity_m_s']*np.cos(np.radians(prop['cant_deg']))
    dimmer = (REQUIREMENTS['dimmer_share']*REQUIREMENTS['climate_top_of_air_W_m2']/(1-REQUIREMENTS['dimmer_share'])
              *np.pi*K.MOON_RADIUS**2)
    annulus = np.pi*((held_row['aperture_radius_km']*1e3)**2-(held_row['window_radius_km']*1e3)**2)
    cases = []
    for demand in ASSUMPTIONS['demand_TW']:
        for conversion in ASSUMPTIONS['conversion']:
            held_area = (demand+held_row['holding_power_mean_TW'])*1e12/(conversion*s)
            fleet_area = demand*1e12/(conversion*s)
            cases.append(dict(
                demand_TW=demand, conversion=conversion,
                held=dict(generation_TW=demand+held_row['holding_power_mean_TW'], collector_km2=held_area/1e6,
                          share_of_annulus=held_area/annulus,
                          # Held collectors burn propellant for their mass and for the sunlight they absorb.
                          propellant_kg_s_per_g_m2=held_area*1e-3*held_row['mean_propulsive_acceleration_m_s2']
                          /thrust_speed,
                          propellant_kg_s_absorbed_light=held_area*s/K.SPEED_OF_LIGHT/thrust_speed),
                ring_fleet=dict(generation_TW=None, collector_km2_for_demand=fleet_area/1e6,
                                note='own load open (attitude); photon keeping uses no thrust power')))
    return dict(dimmer_optical_TW=dimmer/1e12,
                dimmer_electric_TW=[dimmer*c/1e12 for c in ASSUMPTIONS['conversion']],
                annulus_km2=annulus/1e6,
                annulus_one_percent_electric_TW=[.01*annulus*s*c/1e12 for c in ASSUMPTIONS['conversion']],
                cases=cases, open=('Device yield, transfer losses, storage through eclipses and reserves; delivered '
                                   'power to habitat and industry follows once they are modelled.'))


def gates(h, fleet, exhaust, keeping):
    allow = {row['budget_kg_s']: row['share_of_jet_power'] for row in exhaust['allowance']['per_budget']}
    clear = {row['cant_deg']: row for row in exhaust['direct_path']['clear']}
    first = fleet['by_radius'][0]
    slow = [mw for species in exhaust['slow_gas_implied']['by_species'].values() for mw in species['30_days']['deposited_MW']]
    budget_mw = {row['budget_kg_s']: row['deposited_power_MW'] for row in exhaust['allowance']['per_budget']}
    breach = keeping['checks']['clearance']['first_below_150_m_h']/24.
    interior = keeping['checks']['interior_clearance']['min_m']
    days = keeping['design']['orbits']*keeping['period_h']/24.
    return [
        dict(gate='UV transmission (O1, O8)',
             held_zoned=dict(state='open', evidence='The 5 g/m2 annulus film is unproven: the 24.8-120 nm opacity '
                             'proof rests on 22 g/m2 of silica. Continuous area; no panel seams or packing.'),
             ring_fleet=dict(state='open', evidence='Projected overlaps stay positive in the kept bundle; receiver-ray '
                             'coverage, handovers and the annulus film are unproven.')),
        dict(gate='Protected radius (O2)',
             held_zoned=dict(state='met', evidence=f"The {h['aperture_radius_km']:,.0f} km aperture covers four lunar "
                             'radii with the finite Sun on the selected trajectory.'),
             ring_fleet=dict(state='open', evidence='The strip pattern\'s crossing height drifts by '
                             f"{first['crossing_height_error_max_km'][0]:,.0f}-{first['crossing_height_error_max_km'][1]:,.0f}"
                             ' km over a year at 15,000 km; the frozen common orbit is not yet found.')),
        dict(gate='Window spectrum (E2, O4)',
             held_zoned=dict(state='conditional', evidence='The titania stack in the window at 26 g/m2; the dimmer\'s '
                             'form is open, and a band-selective dimmer needs the climate model.'),
             ring_fleet=dict(state='conditional', evidence='As for the held screen, on the rings that cross the window.')),
        dict(gate='Earth\'s shadow and rejected light (O5, O7)',
             held_zoned=dict(state='open', evidence='The optical bound lets redirected light leave in any direction, '
                             'including toward the Moon; the shadow on Earth at eclipse-season new moons is unmapped.'),
             ring_fleet=dict(state='open', evidence='Night-side tiles reflect back past the Moon with light the screen '
                             'already filtered (geometric argument, to check); the shadow on Earth is unmapped.')),
        dict(gate='Stability',
             held_zoned=dict(state='conditional', evidence=f"Held by continuous thrust, {h['installed_TW']:.1f} TW "
                             'installed; a tile that loses thrust drifts about '
                             f"{h['unpowered_drift_first_day_km']:,.0f} km in its first day."),
             ring_fleet=dict(state='open', evidence=f"The kept bundle's edges first come within 150 m at day {breach:.1f}; "
                             f"interior rings stay at least {interior:.0f} m apart for {days:.0f} days in the "
                             'one-tile-per-ring model; independently propagated tiles and the frozen common orbit are '
                             'pending.')),
        dict(gate='Outflows (S7, proposed)',
             held_zoned=dict(state='not met', evidence='At the 45 degree cant the measured plumes send '
                             f"{clear[45.]['gridded_ion_NEXT_share']:.1e} (gridded ion) to "
                             f"{clear[45.]['hall_BPT4000_share']:.1e} (Hall) of their fast ions straight into the "
                             f"protected sphere, against an allowance of {allow[10.]:.1e} of the jet power for "
                             f"10 kg/s and {allow[1.]:.1e} for 1 kg/s. The slow unionized gas alone deposits "
                             f"{min(slow):.0f}-{max(slow):.0f} MW against {budget_mw[1.]:.1f} MW for 1 kg/s and "
                             f"{budget_mw[10.]:.0f} MW for 10 kg/s; the solar wind's pickup of the charged exhaust is "
                             'open and adds to both.'),
             ring_fleet=dict(state='open', evidence='Photon keeping releases nothing; its feasibility check is pending. '
                             'The same control by electric thrust would release '
                             f"{first['electric_alternative_kg_s']['keeping']:,.0f} kg/s at 15,000 km, at or inside "
                             'the planned magnetosphere.')),
        dict(gate='Holding costs less than the atmosphere it saves (S6)',
             held_zoned=dict(state='not met', evidence=f"{h['propellant_kg_s']:,.0f} kg/s of propellant against "
                             f"{REQUIREMENTS['unshielded_loss_kg_s'][0]:,.0f}-{REQUIREMENTS['unshielded_loss_kg_s'][1]:,.0f}"
                             ' kg/s of loss without a shield; met only near the top of that range.'),
             ring_fleet=dict(state='conditional', evidence='No holding propellant if photon keeping suffices; film '
                             'replacement counts as recycled throughput.')),
        dict(gate='Net-positive electricity for habitat',
             held_zoned=dict(state='open', evidence=f"Its own load is {h['holding_power_mean_TW']:.1f} TW mean; "
                             'collectors that ride the screen add propellant; device and transfer are open.'),
             ring_fleet=dict(state='open', evidence='Its own load is the attitude energy, still open; collectors on '
                             'tiles or free-flying carry no holding cost.')),
    ]


def main():
    data = {name: json.loads(path.read_text()) for name, path in PRODUCTS.items()}
    exhaust = data['exhaust_isolation']
    zoned = {c['label']: c for c in data['zoned_aperture']['cases']}[exhaust['design']['case']]
    h = held(zoned)
    fleet = ring_fleet(data['ring_bundle'], data['ring_keeping'], data['ring_screen'], zoned)
    out = dict(schema='terluna.research.integrated-ledger/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={f'results/{name}.json': digest(path) for name, path in PRODUCTS.items()}),
               evidence=('Arithmetic on the named products and on stated requirements and assumptions; no new dynamics. '
                         'The ring fleet\'s inventory is an estimate from the bundle\'s height step, pitch and the '
                         'covering radius at the ring radius, before any oversizing for the strip pattern\'s annual '
                         'motion.'),
               reading_rule=('A candidate must meet every gate; candidates that do are weighed on resources and '
                             'electricity. Gate states are assessments of the cited evidence for the author: met, '
                             'conditional (met if a stated condition holds), open (no model yet decides it) or not met.'),
               requirements=REQUIREMENTS, assumptions=ASSUMPTIONS, states=list(STATES),
               held_zoned=h, ring_fleet=fleet, electricity=electricity(h, zoned),
               outflows=dict(held_zoned=dict(allowance=exhaust['allowance'], direct_path=exhaust['direct_path']['clear'],
                                             slow_gas={k: v for k, v in exhaust['slow_gas'].items() if isinstance(v, dict)
                                                       and 'delivered_by_lifetime' in v},
                                             slow_gas_implied=exhaust['slow_gas_implied'], open=exhaust['open']),
                             ring_fleet=dict(photon_keeping='no exhaust',
                                             electric_alternative_kg_s=[r['electric_alternative_kg_s']
                                                                        for r in fleet['by_radius']])),
               gates=gates(h, fleet, exhaust, data['ring_keeping']))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(held=h, fleet=fleet, electricity={k: v for k, v in out['electricity'].items() if k != 'cases'},
                          gates=[(g['gate'], g['held_zoned']['state'], g['ring_fleet']['state']) for g in out['gates']]),
                     indent=1))


if __name__ == '__main__':
    main()
