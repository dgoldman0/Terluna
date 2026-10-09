"""The array's place in industry and the material flows of the near horizons.

    python -m resources.industry
    # -> resources/results/industry.json

Seven results, read from committed products with the literature values named in industry_sources.json:
- **Places.** The Moon's surface, platforms near the ring fleet, the Sun-Earth L1/L2 hubs, Earth and three source
  populations: the well to climb, sunlight, where released heat goes, light time, and the transfer to 1 AU.
- **Flying with the fleet.** The bundle's single sail loading, the annulus tile's mass beyond its film, and what
  hardware aboard the tiles costs against platforms of their own.
- **The plant's traffic.** Tile exchanges, the velocity change between the plant and the rings, the time on the
  tiles' own sails and the exhaust tugs would release.
- **Computing.** Hardware mass, renewal and heat for the integrated ledger's 100 and 1,000 TW cases.
- **Flows after the build.** The air's losses and the sinks that draw on it, the water, the nitrogen and the film,
  each with its cycle time, over years 500-1,000, years 1,000-10,000 and the 10^9-year operation.
- **People's stocks.** Population scenarios for the Moon and the array: stocks per person, the habitats' shielding,
  turnover, copper and zinc, energy and its heat.
- **The operation.** What the Solar-System-wide operation supplies during and after the build, its packets and its
  braking energy.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from shared.constants import (AU, EARTH_GM, EARTH_MOON_DISTANCE, EARTH_RADIUS, JULIAN_DAY, JULIAN_YEAR_DAYS,
                              MOON_GM, MOON_RADIUS, SOLAR_CONSTANT, SPEED_OF_LIGHT, SUN_GM)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.resources.industry/1'
SCREEN = HERE / 'results' / 'first_screen.json'
HEAT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'array_heat.json'
LEDGER = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'integrated_ledger.json'
LAYOUT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'ring_layout.json'
JOINT = ROOT / 'research' / 'studies' / 'joint_synthesis' / 'results' / 'ledgers.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
SEPTEMBER = ROOT / 'engineering' / 'reference' / 'industrial_architecture' / 'results' / 'derived_industry_calculations.json'
CO2 = ROOT / 'research' / 'studies' / 'atmospheric_co2' / 'results' / 'atmospheric_co2.json'
SOURCES = HERE / 'industry_sources.json'
OUT = HERE / 'results' / 'industry.json'
INPUTS = (SCREEN, HEAT, LEDGER, LAYOUT, JOINT, ATLAS, DRAINAGE, SEPTEMBER, CO2, SOURCES)

YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
MOON_AREA_M2 = 4 * np.pi * MOON_RADIUS ** 2
HUB_M = AU * ((EARTH_GM + MOON_GM) / (3 * SUN_GM)) ** (1 / 3)   # the Sun-Earth L1/L2 points' distance from Earth
BUILD_YEARS = 500.0                      # decisions register: construction within 500 years of 2026
WINDOWS_YR = {'years_500_1000': 500.0, 'years_1000_10000': 9000.0, 'operation_1e9': 1e9}
RING_RADIUS_M = 20000e3                  # the ring fleet's middle radius (decisions register, ring_layout design)
TILE_AREA_M2 = 1e8                       # 10 km tiles (integrated ledger requirements, tile_side_m)

# --- Places -----------------------------------------------------------------------------------------------------
LEO_ALTITUDE_M = 400e3                   # a reference low Earth orbit, the space station's height
SOURCE_AU = {'main_belt': 2.7, 'jupiter_trojans': 5.2, 'kuiper_belt': 40.0}   # representative distances (assumed)
PHASING_FROM_HALO_M_S = (14.0, 45.0)     # Sun-Earth halo orbit to a lunar encounter (Chen, Kawakatsu & Hanada 2016)
HALO_KEEPING_M_S_YR = 1.0                # station keeping of ACE and WIND at Sun-Earth L1 (Roberts 2011)

# --- Flying with the fleet --------------------------------------------------------------------------------------
KEEPING_REACH = 0.25                     # photon keeping reaches about a quarter of the sail force (integrated comparison)
GLOW_PER_PERCENT_W_M2 = (1.66, 2.91)     # glow on the Moon per percent of sunlight absorbed fleet-wide (provisioning
                                         # branch 52dcdb6, README finding 1)
MOON_SHARE_OF_TILE_HEAT = (1 / 303, 1 / 221)   # the Moon receives 0.33-0.45% of the heat the films absorb (same)
TW_PER_KELVIN = (38.0, 61.0)             # TW used on the Moon that warm it 1 K (provisioning branch, finding 2)

# --- The plant's traffic ----------------------------------------------------------------------------------------
PLANTS = (1, 8)                          # one plant in the stack's middle plane, or eight spread over its tilts (assumed)
TUG_EXHAUST_M_S = 30e3                   # the integrated ledger's electric alternative, 30 km/s
HELD_SCREEN_EXHAUST_KG_S = 51002.0       # the zoned held screen's propellant (integrated ledger, held_zoned)

# --- Computing --------------------------------------------------------------------------------------------------
COMPUTE_TW = (100.0, 1000.0)             # the integrated ledger's reference and high electricity cases
SERVER_KG_PER_W = 130.45 / 10.2e3        # NVIDIA DGX H100 datasheet: 130.45 kg, 10.2 kW
SERVER_LIFE_YR = (4.0, 6.0)              # Microsoft's server life, raised from four to six years (The Register 2022)

# --- Flows after the build --------------------------------------------------------------------------------------
LOSS_BUDGETS_KG_S = (1.0, 10.0, 100.0)   # decisions register: designed across 1-100 kg/s
WATER_OVER_H2 = 18.015 / 2.016           # mass of water per mass of hydrogen it holds
CRUST_FILL_YEARS = (1000.0, 10000.0)     # if the crust's pores fill after the seas (assumed span)
LAKE_BURIAL_G_C_M2_YR = (19.0, 48.0)     # lakes and reservoirs (Mendonca et al. 2017), applied to the seas (finding 6)
CN_ATOMIC = (4.0, 20.0)                  # algae 4-10, vascular land plants 20 and more (Meyers 1994)
N_OVER_C = 14.007 / 12.011
CO2_OVER_C = 44.009 / 12.011
CO2_CASE = dict(runoff_mm_yr=349, temperature_c=26.0)    # the CO2 study's design-like weathering case, its fresh-glass
                                                         # factors 1 and 10 bracketing the net sink
RECOVERY = ('0.99', '0.999')             # the film plant's recovery cases in first_screen.json
STOCKPILE_YEARS = 1000.0                 # a pre-air reserve of fresh feed for the first thousand years (assumed)
LIFT_YEARS = 100.0                       # the airless window used for lifting (protection report: shield years 150-350)

# --- People's stocks --------------------------------------------------------------------------------------------
SCENARIOS = {                            # (Moon, array) people: assumptions of this screen within the author's frame
    'first_millennium': (1e9, 2e9),
    'frame_low': (3e9, 7e9),
    'frame_high': (6e9, 14e9),
}
EARTH_PEAK = 10.3e9                      # UN World Population Prospects 2024, the mid-2080s peak
STOCK_T = (32.0 + 11.0, 335.0)           # t/person: decent living, direct + indirect (Velez-Henao & Pauliuk 2023);
                                         # industrial countries in 2010 (Krausmann et al. 2017)
LIFETIME_YR = (65.0, 81.0)               # mean lifetime of demolished buildings, European and US cities
                                         # (Berglund-Brown et al. 2025)
RECYCLED_SHARE_OF_INFLOW = 0.12          # recycling's share of inflows to stocks in 2010 (Krausmann et al. 2017)
SOIL_DENSITY = 1500.0                    # kg/m3, a developed soil (resources/screen.py)
SP413 = dict(people=1e4, shield_t=9.9e6, structure_t=1.5e5, projected_m2=6.8e5, shield_t_m2=4.5)
                                         # NASA SP-413 Table 4-1(a), single torus at 1 rpm; 4.5 t/m2 for 0.5 rem/yr
GROWTH_YEARS = (1000.0, 5000.0)          # years over which the array reaches its scenario population (assumed)
ENERGY_KWH = (33183.0, 75051.0)          # primary energy per person in 2025, Germany and the US (Energy Institute
                                         # Statistical Review 2026, via Our World in Data)
COPPER_KG = (140.0, 300.0)               # in-use copper per person, more-developed countries (UNEP 2010, Table 1)
ZINC_KG = (80.0, 200.0)                  # in-use zinc per person, more-developed countries (UNEP 2010, Table 1)
COPPER_LIFE_YR = (25.0, 40.0)            # copper's lifetime in a building (UNEP 2010)
EOL_RECYCLING = 0.5                      # copper and zinc are among 18 metals above 50% end-of-life (UNEP 2011)
CU_PPM = dict(highland_soil=4.6, ci_chondrite=127.0)    # Apollo 16 soil 67461 (Meyer 2010); CI mean (Lodders 2003)
ZN_PPM = dict(highland_soil=10.4, ci_chondrite=310.0)   # same sources
COPPER_EARTH_T = dict(identified_incl_past=2.1e9, reserves=9.8e8, mine_2024=2.3e7)   # USGS 2025

# --- The operation ----------------------------------------------------------------------------------------------
PACKET_KG = (1e6, 1e8)                   # traffic-safety rule: mature dense traffic in packets of 1e6-1e8 kg

EVIDENCE = ('A screen on committed products (the first resources screen, the array heat screen, the integrated ledger '
            'and ring layout, the joint ledgers, the geography atlas and drainage, the 19 September calculations) with '
            'the literature values named in industry_sources.json and three values from the provisioning branch '
            '(52dcdb6). Closed-form arithmetic: coplanar Hohmann transfers, heat released evenly in every direction, '
            'population, growth and lifetime scenarios. No orbital, industrial, economic or population model is run.')
READING_RULE = ('Masses in kg unless named (t, Mt, Gt = 1e12 kg); flows in kg/s or per year as named; energies in '
                'MJ/kg; powers in W, TW or PW as named; heat in W/m2 averaged over the Moon; times in years or days '
                'as named. Ranges are low-high and pair the low ends and the high ends of their inputs. cycle_years is '
                'the inventory over the flow that drains it. Scenarios (population, growth period, crust filling, '
                'plant count, stockpile) are assumptions of this screen.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def sig(x, n=3):
    """Round to n significant figures, element-wise for lists."""
    if isinstance(x, (list, tuple)):
        return [sig(v, n) for v in x]
    return float(f'{x:.{n}g}')


def span(values):
    return [min(values), max(values)]


def heat_per_tw(distance_m):
    """W/m2 averaged over the Moon from a terawatt radiated evenly in every direction at this distance."""
    return 1e12 * MOON_RADIUS ** 2 / (4 * distance_m ** 2) / MOON_AREA_M2


def hohmann(r_au):
    """A coplanar Hohmann transfer between circular orbits at r_au and 1 AU."""
    r1, r2 = AU, r_au * AU
    a = (r1 + r2) / 2
    v_p, v_a = np.sqrt(SUN_GM * (2 / r1 - 1 / a)), np.sqrt(SUN_GM * (2 / r2 - 1 / a))
    v_inf = abs(v_p - np.sqrt(SUN_GM / r1))
    return dict(departure_dv_km_s=sig(abs(np.sqrt(SUN_GM / r2) - v_a) / 1e3), arrival_v_inf_km_s=sig(v_inf / 1e3),
                braking_mj_kg=sig(0.5 * v_inf ** 2 / 1e6), transit_years=sig(np.pi * np.sqrt(a ** 3 / SUN_GM) / YEAR_S))


def places(screen, heat, joint):
    hub_moon_m = [HUB_M - EARTH_MOON_DISTANCE, HUB_M + EARTH_MOON_DISTANCE]
    leo = EARTH_RADIUS + LEO_ALTITUDE_M
    tile_share = [s * 1e12 / MOON_AREA_M2 for s in MOON_SHARE_OF_TILE_HEAT]
    sources = {}
    for name, r in SOURCE_AU.items():
        sources[name] = dict(distance_au=r, sunlight_w_m2=sig(SOLAR_CONSTANT / r ** 2),
                             light_time_hours=sig([(r - 1) * AU / SPEED_OF_LIGHT / 3600, (r + 1) * AU / SPEED_OF_LIGHT / 3600]),
                             **hohmann(r))
    return dict(
        moon_surface=dict(escape_mj_kg=sig(MOON_GM / MOON_RADIUS / 1e6), to_ring_orbit_mj_kg=screen['film']['lift_to_ring_mj_kg'],
                          sunlight_w_m2=heat['placement']['design_sunlight_global_mean_W_m2'],
                          heat_w_m2_per_tw_used=sig(heat['placement']['moon_W_m2_per_TW_used_on_moon']),
                          note='Counts energy brought from outside the Moon\'s own sunlight, wind and rivers (fusion, '
                               'fission, beams); those move heat the sunlight already brings.'),
        ring_platforms=dict(radius_km=RING_RADIUS_M / 1e3, to_lunar_escape_mj_kg=sig(MOON_GM / (2 * RING_RADIUS_M) / 1e6),
                            sunlight_w_m2=SOLAR_CONSTANT, light_time_s=sig(RING_RADIUS_M / SPEED_OF_LIGHT),
                            heat_w_m2_per_tw_even=sig(heat_per_tw(RING_RADIUS_M)),
                            heat_w_m2_per_tw_from_tile_faces=sig(tile_share),
                            night_half_intercepts_pw=sig(joint['energy']['fleet_night_half_intercepts_w'] / 1e15)),
        sun_earth_l1_l2=dict(distance_from_earth_km=sig(HUB_M / 1e3), distance_from_moon_km=sig([d / 1e3 for d in hub_moon_m]),
                             light_time_to_moon_s=sig([d / SPEED_OF_LIGHT for d in hub_moon_m]),
                             heat_w_m2_per_tw_even=sig([heat_per_tw(d) for d in hub_moon_m[::-1]]),
                             phasing_to_lunar_encounter_m_s=list(PHASING_FROM_HALO_M_S),
                             station_keeping_m_s_per_year=HALO_KEEPING_M_S_YR, sunlight_w_m2=SOLAR_CONSTANT),
        earth=dict(escape_mj_kg=sig(EARTH_GM / EARTH_RADIUS / 1e6),
                   to_low_orbit_mj_kg=sig((EARTH_GM / (2 * leo) + EARTH_GM * (1 / EARTH_RADIUS - 1 / leo)) / 1e6),
                   light_time_to_moon_s=sig(EARTH_MOON_DISTANCE / SPEED_OF_LIGHT)),
        sources=sources,
        ratios=dict(moon_over_ring=sig(heat['placement']['moon_W_m2_per_TW_used_on_moon'] / heat_per_tw(RING_RADIUS_M)),
                    ring_over_hub=sig([heat_per_tw(RING_RADIUS_M) / heat_per_tw(d) for d in hub_moon_m])))


def fleet_carriage(ledger, heat):
    films = ledger['annulus']['films']['silica_4_um']
    loading = films['matched_tile_g_m2'] / films['pressure_coefficient']       # g/m2 per unit pressure coefficient
    window_c = ledger['annulus']['window_pressure_coefficient']
    allowance = films['matched_tile_g_m2'] - films['areal_mass_g_m2']          # g/m2 beyond the film
    fleet = heat['fleet']
    annulus_m2 = [a * 1e6 * (1 - fleet['window_ring_share']) for a in fleet['area_km2']]
    allowance_gt = [allowance * 1e-3 * a / 1e12 for a in annulus_m2]
    rotors = span([g for r in heat['store']['rotors'] for g in r['fleet_rotor_Gt']])
    holding = heat['trim']['eccentricity_holding_share']
    percent_gt = [0.01 * m for m in fleet['optical_mass_Gt']]
    renewal_years = [m / u for m, u in zip(fleet['optical_mass_Gt'], fleet['film_upkeep_Gt_per_year'])]
    used_on_moon = heat['placement']['moon_W_m2_per_TW_used_on_moon']
    return dict(
        loading_g_m2_per_unit_pressure=sig(loading, 4), window_pressure_coefficient=sig(window_c, 4),
        window_tile_g_m2=sig(loading * window_c, 4), annulus_pressure_coefficient=sig(films['pressure_coefficient'], 4),
        annulus_tile_g_m2=sig(films['matched_tile_g_m2'], 4), annulus_film_g_m2=sig(films['areal_mass_g_m2'], 4),
        allowance_g_m2=sig(allowance, 3), allowance_share_of_tile=sig(allowance / films['matched_tile_g_m2'], 2),
        allowance_per_annulus_tile_t=sig(allowance * TILE_AREA_M2 / 1e6), allowance_fleet_gt=sig(allowance_gt),
        store_rotors_gt=sig(rotors, 2), store_share_of_allowance=sig([rotors[0] / allowance_gt[1], rotors[1] / allowance_gt[0]], 2),
        mean_tile_t=sig(fleet['tile_mass_kg'] / 1e3, 4),
        cargo_at_full_keeping_t=sig(KEEPING_REACH * fleet['tile_mass_kg'] / 1e3),
        eccentricity_holding_share_outer_rings=[sig(holding['600']['max'], 2), sig(holding['1000']['max'], 2)],
        uniform_percent=dict(fleet_gt=sig(percent_gt), plant_gt_per_year=sig([percent_gt[0] / renewal_years[0], percent_gt[1] / renewal_years[1]]),
                             glow_w_m2_per_percent_absorbed=list(GLOW_PER_PERCENT_W_M2),
                             equivalent_tw_used_on_moon=sig([g / used_on_moon for g in GLOW_PER_PERCENT_W_M2]),
                             share_of_design_sunlight=sig([g / heat['placement']['design_sunlight_global_mean_W_m2'] for g in GLOW_PER_PERCENT_W_M2], 2)),
        device_heat_w_m2_per_tw=dict(on_tiles=sig([s * 1e12 / MOON_AREA_M2 for s in MOON_SHARE_OF_TILE_HEAT]),
                                     platform_even=sig(heat_per_tw(RING_RADIUS_M))))


def plant_traffic(heat, ledger, layout):
    fleet = heat['fleet']
    tile_kg = fleet['tile_mass_kg']
    upkeep_kg_s = [g * 1e12 / YEAR_S for g in fleet['film_upkeep_Gt_per_year']]
    per_day = [g * 1e12 / tile_kg / JULIAN_YEAR_DAYS for g in fleet['film_upkeep_Gt_per_year']]
    stack = layout['stack']
    rows = stack['ideal_radius_km']['by_height']
    h_km, r_km = np.array([b['height_km'] for b in rows]), np.array([b['radius_km'] for b in rows])
    order = np.argsort(h_km)
    heights = np.linspace(h_km.min(), h_km.max(), stack['rings'])            # rings evenly spread in height
    radii = np.interp(heights, h_km[order], r_km[order]) * 1e3
    tilts = np.degrees(np.arcsin(heights * 1e3 / radii))
    speed = np.sqrt(MOON_GM / radii)
    sail = ledger['ring_fleet']['photon_keeping']['sail_acceleration_m_s2']['26_g_m2'] * KEEPING_REACH
    plants = {}
    for n in PLANTS:
        bands = np.array_split(np.argsort(tilts), n)
        dv = np.empty_like(tilts)
        for band in bands:
            i_plant = band[len(band) // 2]                                       # the plant sits at its band's median
            plane = 2 * speed[band] * np.sin(np.radians(np.abs(tilts[band] - tilts[i_plant])) / 2)
            radius = np.abs(speed[band] - speed[i_plant])                        # small radius change, Hohmann total
            dv[band] = plane + radius
        median, top = float(np.median(dv)), float(dv.max())
        days = median / sail / JULIAN_DAY
        transit = [2 * p * days for p in per_day]
        exhaust = [2 * u * median / TUG_EXHAUST_M_S for u in upkeep_kg_s]
        plants[str(n)] = dict(dv_median_m_s=sig(median), dv_max_m_s=sig(top), sail_days_per_leg=sig(days),
                              tiles_in_transit=sig(transit), share_of_fleet_in_transit=sig([t / f for t, f in zip(transit, fleet['tiles'])], 2),
                              tug_exhaust_kg_s=sig(exhaust), tug_exhaust_share_of_held_screen=sig([e / HELD_SCREEN_EXHAUST_KG_S for e in exhaust], 2),
                              tug_jet_power_tw=sig([0.5 * e * TUG_EXHAUST_M_S ** 2 / 1e12 for e in exhaust]))
    return dict(tile_t=sig(tile_kg / 1e3, 4), upkeep_gt_yr=sig(fleet['film_upkeep_Gt_per_year'], 4), exchanges_per_day=sig(per_day),
                ring_radius_km=sig(span(list(radii / 1e3))), ring_tilt_deg=sig([float(tilts.min()), float(tilts.max())]),
                sail_acceleration_m_s2=sig(sail), plants=plants)


def computing(heat):
    c = heat['computing']
    collectors = heat['collectors']
    heat_per_electric = span([k['heat_per_electric'] for k in collectors])
    server_mt = SERVER_KG_PER_W * 1e12 / 1e9
    radiator_mt = c['radiator_Mt_per_TW']
    cases = {}
    for tw in COMPUTE_TW:
        released = [tw * (1 + h) for h in heat_per_electric]
        cases[f'{tw:g}_tw'] = dict(
            servers_gt=sig(tw * server_mt / 1e3), radiators_gt=sig([tw * m / 1e3 for m in radiator_mt]),
            hardware_gt=sig([tw * (server_mt + m) / 1e3 for m in radiator_mt]),
            server_renewal_gt_yr=sig([tw * server_mt / 1e3 / y for y in SERVER_LIFE_YR[::-1]]),
            collectors_km2=sig([tw * a for a in c['collector_km2_per_TW']]),
            heat_released_tw=sig(released),
            moon_w_m2_near_fleet=sig([r * heat_per_tw(RING_RADIUS_M) for r in released]),
            moon_w_m2_at_hubs=sig([r * heat_per_tw(HUB_M - EARTH_MOON_DISTANCE) for r in released]))
    return dict(server_mt_per_tw=sig(server_mt), radiator_mt_per_tw=sig(radiator_mt), server_life_years=list(SERVER_LIFE_YR),
                heat_per_electric=sig(heat_per_electric), cases=cases)


def windows(kg_s, inventory_kg, label='inventories', operation_kg_s=None):
    """Mass a steady flow moves in each window, and how many times it takes the inventory it drains."""
    out = {}
    for name, years in WINDOWS_YR.items():
        rate = operation_kg_s if (operation_kg_s is not None and name == 'operation_1e9') else kg_s
        mass = [k * years * YEAR_S for k in rate]
        out[name] = {'kg': sig(mass), label: sig([m / inventory_kg for m in mass])}
    return out


def carbon(screen, co2, carbon_kg_s, used_on_moon):
    """The air's CO2 against weathering and burial: the near horizons' one fast drain on the air."""
    cases = [c for c in co2['weathering']['cases']
             if c['runoff_mm_yr'] == CO2_CASE['runoff_mm_yr'] and c['temperature_c'] == CO2_CASE['temperature_c']]
    net = sorted(c['net_gt_yr'] for c in cases)                                   # Gt CO2 a year, the study's runoff
    scale = screen['oxygen_sink']['runoff_rescale']                               # to the corrected runoff (first screen)
    weathering = [n * scale * 1e12 / YEAR_S for n in net]                         # kg CO2/s
    burial = [c * CO2_OVER_C for c in carbon_kg_s]
    total = [w + b for w, b in zip(weathering, burial)]
    inventory = screen['inventories']['co2_kg']
    kilns = [w * YEAR_S / 1e12 * co2['return_step']['heat_gw_per_gt_co2_per_year'] for w in weathering]   # GW
    return dict(inventory_kg=inventory, runoff_rescale=scale, weathering_gt_co2_yr=sig([w * YEAR_S / 1e12 for w in weathering]),
                burial_gt_co2_yr=sig([b * YEAR_S / 1e12 for b in burial]),
                cycle_years=dict(weathering=sig([inventory / (w * YEAR_S) for w in weathering[::-1]]),
                                 weathering_and_burial=sig([inventory / (t * YEAR_S) for t in total[::-1]])),
                kilns_gw=sig(kilns), kilns_heat_w_m2=sig([k / 1e3 * used_on_moon for k in kilns]),
                import_as_carbon_kg_s=sig([t / CO2_OVER_C for t in total]),
                burial_import_as_carbon_kg_s=sig([b / CO2_OVER_C for b in burial]),
                burial_o2_release_gt_yr=sig([b / CO2_OVER_C * 31.998 / 12.011 * YEAR_S / 1e12 for b in burial]),
                note='Kilns return the weathered carbon on the Moon; buried organic carbon needs its own return path or an import.')


def after_build(screen, joint, atlas, heat, co2):
    inv = screen['inventories']
    air, n2, o2 = inv['air_kg'], inv['n2_kg'], inv['o2_kg']
    rates = dict(warm_case_no_magnets=joint['air']['losses']['warm_case_no_magnets_kg_s']['kg_s'],
                 warm_case_september_magnets=joint['air']['losses']['warm_case_september_magnets_kg_s']['kg_s'],
                 cooler_cases_low=joint['air']['losses']['cooler_cases_kg_s']['kg_s'][0],
                 cooler_cases_high=joint['air']['losses']['cooler_cases_kg_s']['kg_s'][1])
    rates.update({f'budget_{b:g}_kg_s': b for b in LOSS_BUDGETS_KG_S})
    losses = {k: dict(kg_s=v, cycle_years=sig(air / (v * YEAR_S)), n2_kg_s=sig(v * n2 / air)) for k, v in rates.items()}
    computed = span([rates[k] for k in ('warm_case_no_magnets', 'warm_case_september_magnets', 'cooler_cases_low', 'cooler_cases_high')])
    o2_kg_s = [m * 1e9 / YEAR_S for m in screen['oxygen_sink']['o2_sink_mt_yr']]
    o2_erosion_kg_s = [m * 1e9 / YEAR_S for m in screen['oxygen_sink']['erosion_limited_o2_mt_yr']]   # over 1e9 years
    water = inv['water']['seas_kg'] + inv['water']['lakes_kg']
    water_loss = [h * WATER_OVER_H2 for h in joint['air']['hydrogen_escape_kg_s']]
    crust = inv['water']['crust_kg']
    crust_kg_s = [crust[0] / (CRUST_FILL_YEARS[1] * YEAR_S), crust[1] / (CRUST_FILL_YEARS[0] * YEAR_S)]
    build_water_kg_s = water / (BUILD_YEARS * YEAR_S)
    sea_m2 = atlas['water_share'] * MOON_AREA_M2
    burial_c = [g * sea_m2 / 1e3 / YEAR_S for g in LAKE_BURIAL_G_C_M2_YR]                  # kg C/s
    nitrogen = [burial_c[0] * N_OVER_C / CN_ATOMIC[1], burial_c[1] * N_OVER_C / CN_ATOMIC[0]]
    film = screen['film']
    fleet_kg = [m * 1e12 for m in film['optical_mass_gt']]
    fresh = {r: film['fresh_feed_kg_s'][r] for r in RECOVERY}
    stock_kg = [STOCKPILE_YEARS * YEAR_S * f for f in fresh['0.99']]
    lift = film['lift_to_ring_mj_kg'] * 1e6
    return dict(
        air=dict(inventory_kg=sig(air, 4), losses=losses, computed_kg_s=computed, windows=windows(computed, air, 'airs')),
        oxygen_to_rock=dict(kg_s=sig(o2_kg_s), cycle_years=sig([o2 / (k * YEAR_S) for k in o2_kg_s[::-1]]),
                            times_warm_case_escape=sig([k / computed[1] for k in o2_kg_s]),
                            erosion_limited_kg_s=sig(o2_erosion_kg_s),
                            windows=windows(o2_kg_s, o2, 'air_oxygens', o2_erosion_kg_s),
                            note='The operation window takes the erosion-limited sink, as first_screen.json does.'),
        water=dict(inventory_kg=sig(water, 4), hydrogen_escape_as_water_kg_s=sig(water_loss),
                   cycle_years=sig([water / (w * YEAR_S) for w in water_loss[::-1]]),
                   crust_kg=crust, crust_fill_years=list(CRUST_FILL_YEARS), crust_kg_s_if_after_the_seas=sig(crust_kg_s),
                   build_water_kg_s=sig(build_water_kg_s), crust_over_build_stream=sig([c / build_water_kg_s for c in crust_kg_s], 2)),
        nitrogen=dict(escape_kg_s=sig([computed[0] * n2 / air, computed[1] * n2 / air]),
                      burial_at_lake_rates_kg_s=sig(nitrogen), sea_carbon_burial_gt_yr=sig([c * YEAR_S / 1e12 for c in burial_c]),
                      burial_cycle_years=sig([n2 / (k * YEAR_S) for k in nitrogen[::-1]]),
                      windows=windows(nitrogen, n2, 'air_nitrogens')),
        carbon=carbon(screen, co2, burial_c, heat['placement']['moon_W_m2_per_TW_used_on_moon']),
        film=dict(upkeep_gt_yr=sig(film['upkeep_gt_yr'], 4), fresh_feed_kg_s=fresh,
                  fleet_replaced_by_fresh_feed_years={r: sig([fleet_kg[0] / (fresh[r][1] * YEAR_S), fleet_kg[1] / (fresh[r][0] * YEAR_S)])
                                                      for r in RECOVERY},
                  windows={r: windows(fresh[r], fleet_kg[0], 'fleets') for r in RECOVERY}),
        pre_air_lift=dict(lift_mj_kg=film['lift_to_ring_mj_kg'], fleet_gt=sig(film['optical_mass_gt']),
                          reserve_years=STOCKPILE_YEARS, reserve_gt_at_99=sig([s / 1e12 for s in stock_kg]),
                          power_gw_over_100_years=dict(fleet=sig([m * lift / (LIFT_YEARS * YEAR_S) / 1e9 for m in fleet_kg]),
                                                       reserve=sig([s * lift / (LIFT_YEARS * YEAR_S) / 1e9 for s in stock_kg]))))


def people(screen, atlas, drainage, heat):
    land_m2 = (atlas['main_land_share'] + atlas['island_share'] - drainage['lakes']['area_share']) * MOON_AREA_M2
    used_on_moon = heat['placement']['moon_W_m2_per_TW_used_on_moon']
    heat_per_electric = span([k['heat_per_electric'] for k in heat['collectors']])
    collector_km2 = heat['computing']['collector_km2_per_TW']
    fleet_kg = [m * 1e12 for m in screen['film']['optical_mass_gt']]
    water = screen['inventories']['water']['total_kg']
    shield_t = SP413['shield_t'] / SP413['people']
    structure_t = SP413['structure_t'] / SP413['people']
    energy_w = [e * 3.6e6 / YEAR_S for e in ENERGY_KWH]                     # W per person
    out = {}
    for name, (moon, array) in SCENARIOS.items():
        stocks = [moon * s * 1e3 for s in STOCK_T]
        turnover = [stocks[0] / (LIFETIME_YR[1] * YEAR_S), stocks[1] / (LIFETIME_YR[0] * YEAR_S)]   # kg/s
        moon_tw = [moon * e / 1e12 for e in energy_w]
        copper = [moon * c for c in COPPER_KG]
        zinc = [moon * z for z in ZINC_KG]
        cu_top = [copper[0] / COPPER_LIFE_YR[1] * EOL_RECYCLING, copper[1] / COPPER_LIFE_YR[0] * EOL_RECYCLING]   # kg/yr
        habitat = [array * (shield_t + structure_t) * 1e3] * 2
        goods = [array * s * 1e3 for s in STOCK_T]
        array_total = [habitat[0] + goods[0], habitat[1] + goods[1]]
        array_tw = [array * e / 1e12 for e in energy_w]
        released = [array_tw[0] * (1 + heat_per_electric[0]), array_tw[1] * (1 + heat_per_electric[1])]
        out[name] = dict(
            moon_people=moon, array_people=array, total_with_earth_peak=sig(moon + array + EARTH_PEAK),
            moon=dict(stocks_gt=sig([s / 1e12 for s in stocks]), soil_depth_over_dry_land_m=sig([s / (land_m2 * SOIL_DENSITY) for s in stocks]),
                      turnover_gt_yr=sig([t * YEAR_S / 1e12 for t in turnover]),
                      fresh_gt_yr_at_earths_recycling=sig([t * YEAR_S / 1e12 * (1 - RECYCLED_SHARE_OF_INFLOW) for t in turnover]),
                      energy_tw=sig(moon_tw), heat_w_m2_if_brought_in=sig([t * used_on_moon for t in moon_tw]),
                      warming_k_if_brought_in=sig([moon_tw[0] / TW_PER_KELVIN[1], moon_tw[1] / TW_PER_KELVIN[0]], 2),
                      copper=dict(stock_gt=sig([c / 1e12 for c in copper]),
                                  over_earths_identified=sig([c / 1e3 / COPPER_EARTH_T['identified_incl_past'] for c in copper], 2),
                                  ore_gt={k: sig([c / (p * 1e-6) / 1e12 for c in copper]) for k, p in CU_PPM.items()},
                                  dry_land_soil_m=sig([c / (CU_PPM['highland_soil'] * 1e-6) / (land_m2 * SOIL_DENSITY) for c in copper]),
                                  top_up_mt_yr=sig([t / 1e9 for t in cu_top]), top_up_kg_s=sig([t / YEAR_S for t in cu_top]),
                                  top_up_over_earths_mine_output=sig([t / 1e3 / COPPER_EARTH_T['mine_2024'] for t in cu_top], 2)),
                      zinc=dict(stock_gt=sig([z / 1e12 for z in zinc]), ore_gt={k: sig([z / (p * 1e-6) / 1e12 for z in zinc]) for k, p in ZN_PPM.items()})),
            array=dict(shield_t_per_person=sig(shield_t), habitat_gt=sig([h / 1e12 for h in habitat]),
                       goods_gt=sig([g / 1e12 for g in goods]), total_gt=sig([t / 1e12 for t in array_total]),
                       fleets=sig([array_total[0] / fleet_kg[1], array_total[1] / fleet_kg[0]]),
                       share_of_build_water=sig([habitat[0] / water[1], habitat[1] / water[0]], 2),
                       habitat_kg_s_over_growth={f'{y:g}_years': sig([h / (y * YEAR_S) for h in habitat]) for y in GROWTH_YEARS},
                       total_kg_s_over_growth={f'{y:g}_years': sig([t / (y * YEAR_S) for t in array_total]) for y in GROWTH_YEARS},
                       energy_tw=sig(array_tw), collectors_km2=sig([array_tw[0] * collector_km2[0], array_tw[1] * collector_km2[1]]),
                       heat_released_tw=sig(released), moon_w_m2=sig([r * heat_per_tw(RING_RADIUS_M) for r in released]),
                       pre_air_lift_tw_over_100_years=sig([h * screen['film']['lift_to_ring_mj_kg'] * 1e6 / (LIFT_YEARS * YEAR_S) / 1e12 for h in habitat])))
    return dict(dry_land_km2=sig(land_m2 / 1e6), stock_t_per_person=list(STOCK_T), lifetime_years=list(LIFETIME_YR),
                energy_gj_per_person=sig([e * 3.6e-3 for e in ENERGY_KWH]), sp413=SP413,
                shield_t_per_m2_projected=sig(SP413['shield_t'] / SP413['projected_m2']), scenarios=out)


def operation(screen, after, stocks, september, places_):
    build = screen['delivery']['flow_kg_s']
    scen = stocks['scenarios']
    array_kg_s = span([v for s in scen.values() for k in s['array']['habitat_kg_s_over_growth'].values() for v in k])
    copper = span([v for s in scen.values() for v in s['moon']['copper']['top_up_kg_s']])
    streams = dict(air_escape=after['air']['computed_kg_s'], film_fresh_feed=[after['film']['fresh_feed_kg_s']['0.999'][0], after['film']['fresh_feed_kg_s']['0.99'][1]],
                   array_habitats=array_kg_s, moon_copper_top_up=copper)
    total = [sum(v[0] for v in streams.values()), sum(v[1] for v in streams.values())]
    braking = {name: src['braking_mj_kg'] for name, src in places_['sources'].items()}
    lo, hi = min(braking.values()), max(braking.values())
    per_day = lambda kg_s: sig([kg_s[0] * JULIAN_DAY / PACKET_KG[1], kg_s[1] * JULIAN_DAY / PACKET_KG[0]])
    transit = {name: sig([b * src['transit_years'] * YEAR_S for b in build]) for name, src in places_['sources'].items()}
    k117 = september['summary']['civilization_power_W']
    return dict(
        build_kg_s=build, after_build_streams_kg_s=streams, after_build_total_kg_s=sig(total),
        build_over_after=sig([build[0] / total[1], build[1] / total[0]]),
        conditional_streams_kg_s=dict(crust_water_if_after_the_seas=after['water']['crust_kg_s_if_after_the_seas'],
                                      carbon_if_imported_for_weathering_and_burial=after['carbon']['import_as_carbon_kg_s']),
        packets_per_day=dict(build=per_day(build), after_build=per_day(total)),
        braking_mj_kg=braking,
        braking_power=dict(build_pw=sig([build[0] * lo * 1e6 / 1e15, build[1] * hi * 1e6 / 1e15]),
                           after_build_tw=sig([total[0] * lo * 1e6 / 1e12, total[1] * hi * 1e6 / 1e12]),
                           build_share_of_k1_17=sig([build[0] * lo * 1e6 / k117, build[1] * hi * 1e6 / k117], 2)),
        build_in_transit_kg=transit, k1_17_w=sig(k117))


def main(argv=None) -> int:
    screen, heat, ledger, layout, joint, atlas, drainage, september, co2 = (
        json.loads(p.read_text()) for p in (SCREEN, HEAT, LEDGER, LAYOUT, JOINT, ATLAS, DRAINAGE, SEPTEMBER, CO2))
    files = ['resources/industry.py']
    places_ = places(screen, heat, joint)
    after = after_build(screen, joint, atlas, heat, co2)
    stocks = people(screen, atlas, drainage, heat)
    comp = computing(heat)
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='resources', files={f: digest(ROOT / f) for f in files},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS}, constants=constants_used(files)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        places=places_, fleet_carriage=fleet_carriage(ledger, heat), plant_traffic=plant_traffic(heat, ledger, layout),
        computing=comp, after_build=after, people=stocks,
        operation=operation(screen, after, stocks, september, places_))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: v for k, v in product.items() if k not in ('producer', 'evidence', 'reading_rule')}, indent=1)[:12000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
