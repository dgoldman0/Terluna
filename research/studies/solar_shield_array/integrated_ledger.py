"""Integrated ledger of the shield candidates (step 2 of the integrated comparison).

The author's criteria of 6 October 2026 (research/decisions.md): the shield system meets the shield
requirements, gives reasonable net-positive electricity for habitat, balances overall resource requirements,
and is stable and safe, with every mass flow judged by direction as well as size. This runner gathers what
the products establish for two candidates, adds the arithmetic that follows from stated assumptions, and
marks what no model yet supplies. It runs no new dynamics.

- held_zoned: the held screen with a zoned aperture, 26 g/m2 in the window and 5 g/m2 in the annulus, on the
  selected moving trajectory (zoned_aperture.json, exhaust_isolation.json).
- ring_fleet: a global screen of nested ring bundles in lunar orbit, kept with photon forces
  (ring_bundle.json, ring_keeping.json, ring_screen.json, photon_control.json, frozen_rings.json,
  attitude_schemes.json, plane_motion.json, bundle_validation.json). Its inventory follows from the bundle's
  geometry at ring radii of 15,000-20,000 km; no full-fleet run exists.

Both use the annulus films of the protection domain (protection/spectra/annulus_film.json).

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
            ('zoned_aperture', 'exhaust_isolation', 'ring_bundle', 'ring_keeping', 'ring_screen', 'photon_control',
             'frozen_rings', 'attitude_schemes', 'plane_motion', 'bundle_validation')}
ANNULUS = ROOT/'protection/spectra/annulus_film.json'
# The annulus films the comparison reads: 0.1 um of titania on 2 or 4 um of silica with the stored coating layers.
FILMS = dict(silica_2_um=(0.1, 2.0), silica_4_um=(0.1, 4.0))
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


def annulus_films(annulus, zoned):
    """The chosen annulus films, and the areal mass at which an annulus tile carries the window tiles' sail loading.

    The dynamics models carry the window's redirected band as reflected (pressure coefficient twice the band's share);
    an annulus tile with the film's own coefficient matches that loading at the window's areal mass times the ratio."""
    rows = {(f['titania_um'], f['silica_um'], f['coated']): f for f in annulus['films']}
    window_pressure = 2*CORE_DIVERTED
    core = zoned['core_sigma_kg_m2']
    films = {}
    for name, (ti, si) in FILMS.items():
        f = rows[(ti, si, True)]
        films[name] = dict(titania_um=ti, silica_um=si, areal_mass_g_m2=f['areal_mass_g_m2'],
                           uv_transmission=f['uv_transmission'], upper_air_heat_W_m2=f['upper_air_heat_W_m2'],
                           visible_transmission=f['band_transmission']['visible_400_700'],
                           pressure_coefficient=f['pressure_coefficient'],
                           matched_tile_g_m2=core*1e3*f['pressure_coefficient']/window_pressure)
    w = annulus['window_stack']
    return dict(films=films, window_pressure_coefficient=window_pressure,
                window_stack=dict(uv_transmission=w['uv_transmission'], upper_air_heat_W_m2=w['upper_air_heat_W_m2']),
                gap_heat_W_m2=annulus['gap_heat_W_m2'],
                lightest_euv_fuv=min((f for f in annulus['films'] if f['coated'] and
                                      f['uv_transmission']['max_10_to_175'] < 1e-6),
                                     key=lambda f: f['areal_mass_g_m2'])['areal_mass_g_m2'])


def ring_fleet(bundle, keeping, screen, zoned, planes=None, films=None):
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
        extra = {}
        cases_pm = [c for c in (planes or {}).get('cases', []) if np.isclose(c['radius_km'], r_km)]
        if cases_pm and films:
            # Oversizing for the pattern's common shift (centred on its mean) and a narrower height step that keeps
            # the overlap through the year's opening, taking the better reference plane; without steering, the edge
            # strips' recession over the year (the largest residual) adds to the oversizing.
            best = min(cases_pm, key=lambda c: c['common']['centred_oversizing_share_of_rings'])
            opening = best['differential']['overlap_change_per_step_m']['opening_max']
            step = bundle['built']['height_pitch_m']
            widened = step/min(step, REQUIREMENTS['clear_side_m']-opening)
            common = best['common']['centred_oversizing_share_of_rings']
            edges = best['differential']['residual_km_max']/best['aperture_km']
            factor = (1+common)*widened
            extra = dict(plane_motion=dict(reference_plane=best['plane'], centred_oversizing_share=common,
                                           edge_recession_share=edges, overlap_opening_m=opening,
                                           narrower_step_factor=widened, ring_factor=factor,
                                           ring_factor_without_steering=(1+common+edges)*widened,
                                           steering_needed_over_available_max=max(
                                               best['photon_steering']['needed_over_available_by_ring'])),
                         matched_annulus={name: dict(annulus_tile_g_m2=f['matched_tile_g_m2'],
                                                     mean_sigma_g_m2=share*core*1e3+(1-share)*f['matched_tile_g_m2'],
                                                     optical_mass_Gt=area*(share*core+(1-share)*f['matched_tile_g_m2']*1e-3)
                                                     /1e12,
                                                     oversized_optical_mass_Gt=factor*area*(share*core+(1-share)*
                                                                                            f['matched_tile_g_m2']*1e-3)/1e12,
                                                     oversized_without_steering_Gt=(1+common+edges)*widened*area*(
                                                         share*core+(1-share)*f['matched_tile_g_m2']*1e-3)/1e12)
                                          for name, f in films['films'].items()})
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
                if steering else None), **extra))
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


def frozen_summary(frozen):
    """Ranges from the last pass of the frozen-orbit iteration."""
    rows = frozen['passes'][-1]

    def span(values):
        return [float(min(values)), float(max(values))]
    by_mass = {}
    for sigma in sorted({r['sigma_kg_m2'] for r in rows}):
        group = [r for r in rows if r['sigma_kg_m2'] == sigma]
        by_mass[f'{sigma*1e3:g}_g_m2'] = dict(
            eccentricity=span([np.hypot(*r['centre']) for r in group]),
            crossing_radius_spread_km=span([r['crossing_radius_km'][1]-r['crossing_radius_km'][0] for r in group]),
            periapsis_min_km=float(min(r['periapsis_min_km'] for r in group)))
    settled = [r for r in rows if r['radius_km'] >= 19000. and r['sigma_kg_m2'] >= .026]
    lagging = [r for r in rows if r['radius_km'] == 15000. and r['tilt_deg'] > 0]
    return dict(by_areal_mass=by_mass,
                height_error_19000_20000_km=span([r['crossing_height_error_max_km'] for r in settled]),
                height_drift_15000_tilted_km_per_day=span([abs(r['height_drift_km_per_day']) for r in lagging]))


def attitude_summary(att):
    """What attitude_schemes establishes for the fleet's attitude, push and plane steering."""
    schemes = att['schemes']
    pick = lambda name: dict(radiation_pitch_orbit_mean_N_m=schemes[name]['radiation_orbit_mean_N_m']['pitch'],
                             store_trim_25_N_m_s=schemes[name]['store_N_m_s']['trim_0.25']['median'],
                             store_by_areal_mass_trim_25_N_m_s={k: v['store_N_m_s']['trim_0.25']
                                                                for k, v in schemes[name]['by_areal_mass'].items()},
                             eccentricity_drive_over_as_built=schemes[name]['eccentricity_drive']['over_as_built'])
    spans = [s for row in att['clearance']['envelope'] for s in row['pitch_allowed_deg']]
    return dict(pitch_allowed_deg=[min(s[0] for s in spans), max(s[1] for s in spans)],
                edge_on_pitch=dict(clearance=att['clearance']['laws']['edge_on_pitch'],
                                   eccentricity_drive_over_as_built=att['edge_on_pitch']['eccentricity_drive_over_as_built']),
                as_built=pick('as_built'),
                centred_filter={name: pick(name) for name in ('balanced_240_m', 'balanced_400_m', 'balanced_560_m')},
                rolled={name: dict(pick(name), attitude_torque_p95_N_m=schemes[name]['attitude_torque_N_m']['p95'])
                        for name in ('rolled_30_deg', 'rolled_60_deg', 'rolled_80_deg')},
                moving_mass={k: dict(ballast_over_tile_mass=v['ballast_over_tile_mass'],
                                     store_trim_25_N_m_s=v['store_N_m_s']['trim_0.25'])
                             for k, v in att['moving_mass']['62.7_g_m2']['by_reach'].items()},
                steering_crossing_height_km_per_day=schemes['steering']['crossing_height_rate_km_per_day']['median'],
                steering_strip_turn_deg_per_day=schemes['steering']['strip_turn_deg_per_day']['median'],
                tile_areal_mass_kg_m2=att['tile_areal_mass_kg_m2'], forced_eccentricity=att['forced_eccentricity_as_built'])


def gates(h, fleet, exhaust, keeping, photon, frozen, films, planes):
    allow = {row['budget_kg_s']: row['share_of_jet_power'] for row in exhaust['allowance']['per_budget']}
    held_window = frozen['by_areal_mass']['26_g_m2']
    light = frozen['by_areal_mass']['5_g_m2']
    move = photon['translation']
    pitch = photon['attitude']['by_axis']['pitch_about_orbit_normal']['radiation_orbit_mean_N_m_median']
    centre = photon['attitude']['pressure_centre_offset_m']['median']/1e3
    clear = {row['cant_deg']: row for row in exhaust['direct_path']['clear']}
    first = fleet['by_radius'][0]
    slow = [mw for species in exhaust['slow_gas_implied']['by_species'].values() for mw in species['30_days']['deposited_MW']]
    budget_mw = {row['budget_kg_s']: row['deposited_power_MW'] for row in exhaust['allowance']['per_budget']}
    magnets = {(r['budget_kg_s'], r['thruster'], r['cant_deg']): r['slow_gas_capture_needed']
               for r in exhaust['magnetosphere']['with_magnets']}
    capture = magnets[(1., 'gridded_ion_NEXT', 60.)]
    hall_fails = all(v is None for (b, t, _), v in magnets.items() if b == 1. and t == 'hall_BPT4000')
    breach = keeping['checks']['clearance']['first_below_150_m_h']/24.
    interior = keeping['checks']['interior_clearance']['min_m']
    days = keeping['design']['orbits']*keeping['period_h']/24.
    thin = films['films']['silica_2_um']
    thick = films['films']['silica_4_um']
    window_uv = films['window_stack']['uv_transmission']['below_175']
    gap_max = films['gap_heat_W_m2']['0.0002']['maximum']
    att = fleet['attitude']
    centred = att['centred_filter']['balanced_400_m']
    val = fleet['validation']
    # The plane motion at each radius, from the better reference plane.
    by_radius = {}
    for c in planes['cases']:
        best = by_radius.get(c['radius_km'])
        if best is None or c['common']['centred_oversizing_share_of_rings'] < best['common']['centred_oversizing_share_of_rings']:
            by_radius[c['radius_km']] = c
    near, mid = by_radius[15000.], by_radius[20000.]
    steer = lambda c: max(c['photon_steering']['needed_over_available_by_ring'])
    return [
        dict(gate='UV transmission (O1, O8)',
             held_zoned=dict(state='open', evidence='A tenth of a micrometre of titania with the stack\'s coating '
                             f"layers stops every wavelength from 10 to 175 nm below 1e-6 from {films['lightest_euv_fuv']:.1f} "
                             'g/m2, but light films pass solar X-rays: on 2 um of silica '
                             f"({thin['areal_mass_g_m2']:.1f} g/m2, the 5 g/m2 annulus's class) the film passes "
                             f"{thin['uv_transmission']['below_175']:.1e} of the sunlight below 175 nm, on 4 um "
                             f"{thick['uv_transmission']['below_175']:.1e}, against {window_uv:.1e} for the window "
                             'stack and 3e-5 for O1\'s tight swarm level. Counting every transmitted X-ray as heating '
                             f"the upper air, the thinner film adds {thin['upper_air_heat_W_m2']['maximum']:.1e} W/m2 at "
                             f"solar maximum, against {gap_max:.1e} W/m2 for gaps at the standard level; slant paths "
                             'above the limb and X-ray depths are the atmosphere domain\'s to settle. Continuous area; no '
                             'panel seams or packing.'),
             ring_fleet=dict(state='open', evidence=f"With every tile of a {val['tiles']}-tile patch propagated and kept "
                             f"for {val['orbits']} orbits, none of {val['rays']['tested']:,} receiver rays through its "
                             'interior toward the finite Sun went uncovered (least margin '
                             f"{val['rays']['margin_min_m']:.0f} m). The annulus film as for the held screen; "
                             'coverage of a full fleet, its handovers and its edges are unchecked.')),
        dict(gate='Protected radius (O2)',
             held_zoned=dict(state='met', evidence=f"The {h['aperture_radius_km']:,.0f} km aperture covers four lunar "
                             'radii with the finite Sun on the selected trajectory.'),
             ring_fleet=dict(state='open', evidence='On a Sun-tracking eccentricity (apolune toward the Sun, e '
                             f"{held_window['eccentricity'][0]:.2f}-{held_window['eccentricity'][1]:.2f} for 26 g/m2 tiles) "
                             'the sunward crossing radius stays within '
                             f"{held_window['crossing_radius_spread_km'][0]:,.0f}-{held_window['crossing_radius_spread_km'][1]:,.0f}"
                             ' km over a year. The tide holds the ring planes near the Moon\'s orbit about Earth while the '
                             'Sun moves 5.1 degrees above and below it, so the whole strip pattern shifts each way about '
                             f"its mean by {mid['common']['centred_oversizing_km_per_edge']:,.0f} km at 20,000 km, which "
                             f"oversizing covers with {mid['common']['centred_oversizing_share_of_rings']:.0%} more rings. "
                             'The planes also differ by tilt: at 15,000 km the strips cross within weeks, beyond photon '
                             f"steering by up to {steer(near):.0f} times; at 20,000 km, with rings set from the Moon\'s "
                             'orbit plane, the interior overlaps hold (the year\'s largest opening is '
                             f"{mid['differential']['overlap_change_per_step_m']['opening_max']:.0f} m per height step "
                             'against a 604 m overlap), while the edge strips recede by up to '
                             f"{mid['differential']['residual_km_max']:,.0f} km over the year, beyond photon steering by "
                             f"{steer(mid):.1f} times for the outermost rings. 5 g/m2 tiles are forced to e "
                             f"{light['eccentricity'][0]:.2f}-{light['eccentricity'][1]:.2f}, with perilune as low as "
                             f"{light['periapsis_min_km']:,.0f} km.")),
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
             ring_fleet=dict(state='conditional', evidence=f"The kept bundle's edges first come within 150 m at day "
                             f"{breach:.1f}; interior rings stay at least {interior:.0f} m apart for {days:.0f} days in "
                             f"the one-tile-per-ring model, and with every tile of a {val['tiles']}-tile patch propagated "
                             f"the interior keeps {val['interior_clearance_min_m']:.0f} m for {val['orbits']} orbits. "
                             f"In a bundle of sliding rings a tile may pitch only {att['pitch_allowed_deg'][0]:.0f} to "
                             f"{att['pitch_allowed_deg'][1]:+.0f} degrees from radial-facing, so edge-on turns are ruled "
                             'out. A filter centred on each tile (the pitch plus 400 m) cuts the steady pitch torque from '
                             f"{att['as_built']['radiation_pitch_orbit_mean_N_m']:,.0f} to "
                             f"{centred['radiation_pitch_orbit_mean_N_m']:,.0f} N m, and 25% reflectivity trim then holds "
                             f"it with a momentum store of {centred['store_trim_25_N_m_s']:.1e} N m s per tile: met if the "
                             'centred filter, the trim and the store are built and the along-track keeping holds each '
                             'tile within its overlap.')),
        dict(gate='Outflows (S7, proposed)',
             held_zoned=dict(state='not met', evidence='At the 45 degree cant the measured plumes send '
                             f"{clear[45.]['gridded_ion_NEXT_share']:.1e} (gridded ion) to "
                             f"{clear[45.]['hall_BPT4000_share']:.1e} (Hall) of their fast ions straight into the "
                             f"protected sphere, against an allowance of {allow[10.]:.1e} of the jet power for "
                             f"10 kg/s and {allow[1.]:.1e} for 1 kg/s. The slow unionized gas alone deposits "
                             f"{min(slow):.0f}-{max(slow):.0f} MW against {budget_mw[1.]:.1f} MW for 1 kg/s and "
                             f"{budget_mw[10.]:.0f} MW for 10 kg/s. Without the magnets the solar wind's pickup of "
                             'the charged exhaust adds to both. The September magnets hold off the exhaust ions but '
                             'pass the charge-exchanged fast atoms and the slow gas: for 1 kg/s, gridded ion thrusters '
                             f"canted 60 degrees leave room only if {capture[0]:.0%}-{capture[1]:.0%} of the unionized "
                             'gas is captured' + (', and Hall thrusters exceed the allowance at every cant.'
                                                  if hall_fails else '.')),
             ring_fleet=dict(state='conditional', evidence='Photon keeping releases nothing. Along the kept run, tilts '
                             'up to 2 degrees with 25% reflectivity trim reach '
                             f"{move['2_deg_0.25']['interior']['reach_over_demand']:.1f} times the keeping demand over "
                             f"each orbit, {move['2_deg_0.25']['interior']['instant_share']:.0%} of it at the instant "
                             'asked, and the edge rings need 5 degrees; met if a phase-scheduled law delivers it. The '
                             'attitude is held by trim and a momentum store, which release nothing. The same control by '
                             'electric thrust would '
                             f"release {first['electric_alternative_kg_s']['keeping']:,.0f} kg/s at 15,000 km, at or "
                             'inside the planned magnetosphere.')),
        dict(gate='Holding costs less than the atmosphere it saves (S6)',
             held_zoned=dict(state='not met', evidence=f"{h['propellant_kg_s']:,.0f} kg/s of propellant against "
                             f"{REQUIREMENTS['unshielded_loss_kg_s'][0]:,.0f}-{REQUIREMENTS['unshielded_loss_kg_s'][1]:,.0f}"
                             ' kg/s of loss without a shield; met only near the top of that range.'),
             ring_fleet=dict(state='conditional', evidence='No holding propellant if photon keeping suffices; film '
                             'replacement counts as recycled throughput.')),
        dict(gate='Net-positive electricity for habitat',
             held_zoned=dict(state='open', evidence=f"Its own load is {h['holding_power_mean_TW']:.1f} TW mean; "
                             'collectors that ride the screen add propellant; device and transfer are open.'),
             ring_fleet=dict(state='open', evidence='Its own loads are reflectivity trim and a momentum store near '
                             f"{centred['store_trim_25_N_m_s']:.0e} N m s per tile, neither sized for power yet; "
                             'collectors on tiles or free-flying carry no holding cost.')),
    ]


def main():
    data = {name: json.loads(path.read_text()) for name, path in PRODUCTS.items()}
    exhaust = data['exhaust_isolation']
    zoned = {c['label']: c for c in data['zoned_aperture']['cases']}[exhaust['design']['case']]
    h = held(zoned)
    annulus = json.loads(ANNULUS.read_text())
    films = annulus_films(annulus, zoned)
    fleet = ring_fleet(data['ring_bundle'], data['ring_keeping'], data['ring_screen'], zoned, data['plane_motion'], films)
    frozen = frozen_summary(data['frozen_rings'])
    fleet['frozen_orbit'] = frozen
    photon = data['photon_control']
    fleet['photon_control'] = dict(translation=photon['translation'],
                                   pitch_torque_orbit_mean_N_m=photon['attitude']['by_axis']['pitch_about_orbit_normal']
                                   ['radiation_orbit_mean_N_m_median'],
                                   pressure_centre_offset_m=photon['attitude']['pressure_centre_offset_m'],
                                   trim_needed_lit_median=photon['attitude']['trim_needed_lit_median'])
    fleet['attitude'] = attitude_summary(data['attitude_schemes'])
    val = data['bundle_validation']
    fleet['validation'] = dict(orbits=val['design']['orbits'], tiles=589, interior_clearance_min_m=val['interior_clearance']['min_m'],
                               patch_clearance_min_m=val['clearance']['min_m'],
                               patch_first_below_150_m_h=val['clearance']['first_below_150_m_h'],
                               service_overlaps_min_m=dict(along=val['service']['along_ring_overlap_min_m'],
                                                           across=val['service']['adjacent_ring_overlap_min_m']),
                               rays=dict(tested=val['rays']['tested'], uncovered=val['rays']['uncovered'],
                                         margin_min_m=val['rays']['margin_min_m']),
                               keeping_m_s2=dict(interior_p95=val['keeping']['interior_tiles_m_s2']['p95'],
                                                 ends_and_edges_p95=val['keeping']['patch_ends_and_edges_m_s2']['p95'],
                                                 sail=val['keeping']['sail_acceleration_normal_incidence_m_s2']))
    out = dict(schema='terluna.research.integrated-ledger/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={**{f'results/{name}.json': digest(path) for name, path in PRODUCTS.items()},
                                     'protection/spectra/annulus_film.json': digest(ANNULUS)}),
               evidence=('Arithmetic on the named products and on stated requirements and assumptions; no new dynamics. '
                         'The ring fleet\'s inventory is an estimate from the bundle\'s height step, pitch and the '
                         'covering radius at the ring radius, before any oversizing for the strip pattern\'s annual '
                         'motion.'),
               reading_rule=('A candidate must meet every gate; candidates that do are weighed on resources and '
                             'electricity. Gate states are assessments of the cited evidence for the author: met, '
                             'conditional (met if a stated condition holds), open (no model yet decides it) or not met.'),
               requirements=REQUIREMENTS, assumptions=ASSUMPTIONS, states=list(STATES),
               held_zoned=h, ring_fleet=fleet, annulus=films, electricity=electricity(h, zoned),
               outflows=dict(held_zoned=dict(allowance=exhaust['allowance'], direct_path=exhaust['direct_path']['clear'],
                                             slow_gas={k: v for k, v in exhaust['slow_gas'].items() if isinstance(v, dict)
                                                       and 'delivered_by_lifetime' in v},
                                             slow_gas_implied=exhaust['slow_gas_implied'],
                                             magnetosphere=exhaust['magnetosphere'], open=exhaust['open']),
                             ring_fleet=dict(photon_keeping='no exhaust',
                                             electric_alternative_kg_s=[r['electric_alternative_kg_s']
                                                                        for r in fleet['by_radius']])),
               gates=gates(h, fleet, exhaust, data['ring_keeping'], photon, frozen, films, data['plane_motion']))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(held=h, fleet=fleet, electricity={k: v for k, v in out['electricity'].items() if k != 'cases'},
                          gates=[(g['gate'], g['held_zoned']['state'], g['ring_fleet']['state']) for g in out['gates']]),
                     indent=1))


if __name__ == '__main__':
    main()
