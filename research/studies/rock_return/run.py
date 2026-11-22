"""The living return of rock: can animals carry back what erosion takes, and does the loop hold itself?

    python -m research.studies.rock_return.run        # seconds; reads results/rock.json (rock.py, about four minutes)

The terrain study (13 October 2026) found that keeping the Open Moon's lakes and relief is one loop, the sediment
dredged and carried back uphill, and priced it for machines. Two habits that are each common on Earth, carrying earth
and stones for nests and migrating far as part of a life, would do the same work if one animal joined them. The
register's row of 22 November 2026 sets the aim: the hybrid and used land may have designed and engineered upkeep, and
the wild should if possible keep a dynamically stable equilibrium on its own. This study sizes that living return on
the design surface. It asks five things:

- the job: how much rock leaves whose land, where it first settles, and how far and how high the way back is (rock.py);
- the reach: how much of it lies where a flyer, a wader or a diver can get at it, from the river deltas screen's shares
  for the delta plain and front and the floors' depths;
- the carriers: what flying it back costs at lunar gravity against Earth's, how many flyers of the sky ecology study's
  albatross scale it takes, and what they eat;
- staying put: how thin the returned rock lies, how close to its source it must be set down for soil creep to pass it
  on without oversteepening, and what cement the rock that leaves in solution and the carriers' own nitrogen can give;
- the loop: whether carriers that multiply where sediment gathers return exactly what arrives, and how fast they
  recover against how slowly the waters fill.

Writes results/rock_return.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from research.studies.rock_return import carriers as CA
from research.studies.sky_ecology import food_web as FW
from shared.constants import JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'rock_return.json'
SCHEMA = 'terluna.research.rock-return/1'
ROCK = HERE / 'results' / 'rock.json'
SKY = ROOT / 'research' / 'studies' / 'sky_ecology' / 'results' / 'sky_ecology.json'
DELTAS = ROOT / 'research' / 'studies' / 'river_deltas' / 'results' / 'river_deltas.json'
TERRAIN = ROOT / 'research' / 'studies' / 'terrain_through_time' / 'results' / 'terrain_through_time.json'
BIOMASS = ROOT / 'research' / 'studies' / 'light_limitation' / 'results' / 'biomass_comparison.json'
PHOSPHORUS = ROOT / 'resources' / 'results' / 'sea_phosphorus.json'
INPUTS = (ROCK, SKY, DELTAS, TERRAIN, BIOMASS, PHOSPHORUS)
YEAR_S = 86400.0 * JULIAN_YEAR_DAYS

FLYER = 'wandering albatross (Diomedea exulans)'     # the sky ecology study's flyer six times the albatross in every length
PAYLOADS = (0.25, 0.6, 1.0)          # load over body mass (design guesses: a quarter is a heavy load for a bird on Earth, where the same wing bears a sixth of the mass)
TRIPS = (300, 1000, 3000)            # loaded trips a year (design guesses)
SIZES_KG = dict(hornero=0.06, crow=1.0, goose=10.0)   # smaller carriers, by Earth birds of those weights (design guesses); upkeep by the sky ecology study's allometry
EARTH_PAYLOADS = (0.05, 0.25)        # loads Earth's birds carry: a hornero's pellet of mud, a raptor's prey (recalled)
CRUISE_M_S = 15.0                    # design guess
HANDLING_S = 1200.0                  # design guess: loading and setting down
TEXTURES = ('clay', 'clay_loam', 'loam', 'sandy_loam', 'sand')      # the river deltas screen's mature soils
DIVE_M = 200.0                       # a diving bird's reach by breath; Earth's deepest dive to about 200 m (recalled)
EARTH_DIVES_M = dict(diving_bird=200.0, seal=1500.0, whale=2000.0)  # recalled orders of Earth's deepest divers
DELTA_ORGANIC_C = (0.007, 0.015)     # organic carbon of deltaic sediment (Hedges & Keil 1995, as the carbon cycle study takes it)
KJ_PER_G_C = 40.0                    # energy of organic matter per gram of its carbon (18 kJ/g dry at 45% carbon)
DIGESTIBLE = (0.05, 0.3)             # share of a sediment's organic matter an animal uses (design guess on deposit feeders)
FOOD_KJ_PER_G = 15.0                 # usable energy of dry food, fruit and prey (design guess)
FOOD_N = (0.01, 0.03)                # nitrogen in dry food (design guess: fruit to prey)
CARBON_IN_DRY_BIOMASS = 0.45         # the biomass comparison's
LOOSE_SOIL_KG_M3 = 1300.0            # returned rock as loose soil (design guess)
REPOSE_SLOPE = math.tan(math.radians(30.0))   # about where loose soil slides (design guess)
LOOPS = dict(slow=dict(growth=0.05, death=0.02, half_years=100.0), central=dict(growth=0.10, death=0.03, half_years=30.0),
             fast=dict(growth=0.30, death=0.10, half_years=10.0))   # design guesses: a great flyer's yearly increase and loss, and the years of supply standing where carriers find half their fill
SPANS_YR = (1e4, 1e5, 1e6, 1e7, 1e8)
UNFILLED_TO_THE_WILD = 0.5           # the register's row of 22 November 2026: half of the open land people do not fill goes to the wild, half to the hybrid land; its rock is taken to divide as its area does


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def flyer_of(sky):
    row = next(f for f in sky['food_web']['flyers'] if f['flyer'] == FLYER and f['edge'] == 'low')
    return dict(mass_kg=row['mass_kg'], span_m=row['span_m'], upkeep_j_day=row['fmr_j_day'], upkeep_range_j_day=row['fmr_range_j_day'])


def delta_shares(deltas):
    """The share of a river's load the delta plain keeps, and the plain and front together, for each mature soil."""
    out = {}
    for t in TEXTURES:
        a = deltas['states'][f'mature/{t}']['central']['all']
        out[t] = dict(plain=a['keep_share_by_load'], plain_and_front=a['keep_with_underflow_by_load'])
    return out


def the_job(rock, phosphorus):
    central = rock['rock']['kg_yr_central']
    slow, fast = rock['setting']['slow_and_fast_against_central']
    land_m2 = rock['rock']['land_km2'] * 1e6
    lowering = rock['rock']['mean_lowering_m_per_myr']
    grade = phosphorus['land_loss']['land_soil_grade']
    by_kind = {}
    for key, kk in rock['kinds'].items():
        by_kind[key] = dict(what=kk['what'], **{name: dict(land_km2=kk[name]['land_km2'], gt_yr=kk[name]['rock_share'] * central / 1e12,
                                                           share=kk[name]['rock_share']) for name in ('wild', 'hybrid', 'peoples', 'open_beyond_need')})
    return dict(gt_yr=dict(slow=slow * central / 1e12, central=central / 1e12, fast=fast * central / 1e12),
                kg_m2_yr_central=central / land_m2, mean_lowering_m_per_myr=lowering,
                returned_where_it_left_mm_per_century=dict(as_rock=lowering * 0.1, as_loose_soil=lowering * 0.1 * rock['setting']['rock_density_kg_m3'] / LOOSE_SOIL_KG_M3),
                phosphorus_with_it_mt_yr=[central * g / 1e9 for g in grade], by_kind_of_land=by_kind,
                against_the_terrain_study=dict(rock_gt_yr=rock['setting']['terrain_study']['rock_gt_yr_central'],
                                               mean_lift_m=rock['setting']['terrain_study']['mean_lift_m'],
                                               mean_distance_km=rock['setting']['terrain_study']['mean_distance_km'],
                                               here_mean_lift_m=rock['the_way_back']['all']['lift_m']['mean'],
                                               here_mean_distance_km=rock['the_way_back']['all']['haul_m']['mean'] / 1e3))


def the_reach(rock, shares):
    """Who can reach what: of all the rock and of the wild land's, by the deltas' texture from clay to sand."""
    floors = rock['floors']
    within = dict(lakes=floors['lakes_that_spill'][f'under {DIVE_M:g} m'], closed=floors['lakes_that_pass_nothing_on'][f'under {DIVE_M:g} m'],
                  seas=floors['seas'][f'under {DIVE_M:g} m'])

    def one(by_place):
        total = sum(by_place.values())
        out = {}
        for t in ('clay', 'sand'):
            s = shares[t]
            parts = [CA.reach(0.0, by_place['lake_that_spills'] / total, s['plain'], s['plain_and_front'], within['lakes']),
                     CA.reach(0.0, by_place['lake_that_passes_nothing_on'] / total, s['plain'], s['plain_and_front'], within['closed']),
                     CA.reach(0.0, by_place['sea'] / total, s['plain'], s['plain_and_front'], within['seas'])]
            row = {k: by_place['dry_hollow'] / total * (0.0 if k == 'left_on_the_deeper_floors' else 1.0) + sum(p[k] for p in parts) for k in parts[0]}
            out[f'deltas of {t}'] = row
        return dict(first_settles={k: v / total for k, v in by_place.items()}, reach=out)

    everything = {k: rock['the_way_back'][k]['rock_share'] for k in ('sea', 'lake_that_spills', 'lake_that_passes_nothing_on', 'dry_hollow')}
    out = dict(all_the_land=one(everything), floor_within_a_dive=within, dive_m=DIVE_M,
               pressure_at_the_dive_atm=CA.pressure_atm(DIVE_M),
               moons_depth_at_the_pressure_of_earths_deepest_dives_m={k: CA.depth_at_earths_pressure(v) for k, v in EARTH_DIVES_M.items()},
               raising_a_kg_of_grains_through_water_j={f'{d:g} m': [CA.lift_in_water_j_per_kg(d, e) for e in (0.25, 0.2)] for d in (200.0, 1000.0, 3000.0)})
    for key, kk in rock['kinds'].items():
        out[f'the_wild_land_{key}'] = one(kk['wild']['by_place'])
    return out


def the_carriers(rock, sky, job, biomass):
    glide, muscle = sky['food_web']['ranges']['glide_ratio'], sky['food_web']['ranges']['muscle_efficiency']
    flyer = flyer_of(sky)
    cheap, dear = (glide[1], muscle[1]), (glide[0], muscle[0])

    def costs(way, gravity, payloads):
        haul, lift = way['haul_m']['mean'], way['lift_m']['mean']
        home = way['glides_home_share']
        rows = {}
        for p in payloads:
            lo = CA.delivered_cost_j_per_kg(haul, lift, *cheap, p, home[f'glide {cheap[0]:g}'], gravity)
            hi = CA.delivered_cost_j_per_kg(haul, lift, *dear, p, home[f'glide {dear[0]:g}'], gravity)
            rows[f'load {p:g} of body'] = [lo / 1e3, hi / 1e3]
        return dict(the_load_alone=[CA.payload_cost_j_per_kg(haul, lift, *cheap, gravity) / 1e3, CA.payload_cost_j_per_kg(haul, lift, *dear, gravity) / 1e3],
                    with_the_carrier=rows)

    way = rock['the_way_back']
    cost = dict(moon={k: costs(way[k], MOON_SURFACE_GRAVITY, PAYLOADS) for k in ('all', 'dry_hollow', 'lake_that_spills')},
                earth=costs(way['all'], STANDARD_GRAVITY, EARTH_PAYLOADS))
    usable = [1e3 * c * KJ_PER_G_C * d for c, d in zip(DELTA_ORGANIC_C, DIGESTIBLE)]
    cost['usable_energy_in_a_kg_of_sediment_kj'] = usable
    cost['sediments_energy_over_the_flight'] = dict(moon=[usable[0] / cost['moon']['all']['with_the_carrier']['load 0.25 of body'][1],
                                                          usable[1] / cost['moon']['all']['with_the_carrier']['load 1 of body'][0]],
                                                    earth=[usable[0] / cost['earth']['with_the_carrier']['load 0.05 of body'][1],
                                                           usable[1] / cost['earth']['with_the_carrier']['load 0.25 of body'][0]])

    # how many flyers of the albatross scale, and what the hauling adds to their upkeep
    central = job['gt_yr']['central'] * 1e12
    land_m2 = rock['rock']['land_km2'] * 1e6
    haul = way['all']['haul_m']['mean']
    numbers = {}
    for trips in TRIPS:
        for p in PAYLOADS:
            n = CA.carriers_for(central, flyer['mass_kg'], p, trips)
            per_kg = [1e3 * x for x in cost['moon']['all']['with_the_carrier'][f'load {p:g} of body']]
            hauling_day = [trips * p * flyer['mass_kg'] * x / JULIAN_YEAR_DAYS for x in per_kg]
            whole = [n * (flyer['upkeep_j_day'] + h) * JULIAN_YEAR_DAYS for h in hauling_day]           # J a year
            numbers[f'{trips} trips, load {p:g} of body'] = dict(
                rock_t_a_flyer_a_year=trips * p * flyer['mass_kg'] / 1e3, body_weights_a_year=trips * p,
                flyers_for_all_the_land_central=n, live_mass_mt=n * flyer['mass_kg'] / 1e9,
                hours_a_year=CA.hours_a_year(trips, haul, CRUISE_M_S, HANDLING_S), share_of_the_year=CA.hours_a_year(trips, haul, CRUISE_M_S, HANDLING_S) / (24.0 * JULIAN_YEAR_DAYS),
                hauling_over_upkeep=[h / flyer['upkeep_j_day'] for h in hauling_day],
                metabolism_gw=[w / YEAR_S / 1e9 for w in whole], dry_food_mt_yr=[w / (FOOD_KJ_PER_G * 1e6) / 1e9 for w in whole],
                dry_food_g_per_m2_of_land_yr=[w / (FOOD_KJ_PER_G * 1e3) / land_m2 for w in whole])
    plant = {k: v['land_g_c_m2_year'] / CARBON_IN_DRY_BIOMASS for k, v in biomass['per_area']['cases'].items()}
    middle = numbers['1000 trips, load 0.6 of body']
    by_kind = {key: {name: middle['flyers_for_all_the_land_central'] * kk[name]['share'] for name in ('wild', 'hybrid', 'peoples', 'open_beyond_need')}
               for key, kk in job['by_kind_of_land'].items()}
    # the same job by smaller carriers: more of them, closer together, dearer to keep
    sizes = {}
    for name, mass in dict(SIZES_KG, albatross_scale=flyer['mass_kg']).items():
        upkeep = FW.bird_fmr_j_day(mass)
        n = CA.carriers_for(central, mass, 0.6, 1000)
        hauling_day = [1000 * 0.6 * mass * 1e3 * x / JULIAN_YEAR_DAYS for x in cost['moon']['all']['with_the_carrier']['load 0.6 of body']]
        whole = [n * (upkeep + h) * JULIAN_YEAR_DAYS for h in hauling_day]
        food = [w / (FOOD_KJ_PER_G * 1e3) / land_m2 for w in whole]
        sizes[name] = dict(mass_kg=mass, upkeep_j_day=upkeep, carriers=n, live_mass_mt=n * mass / 1e9, a_km2_of_land=n / (land_m2 / 1e6),
                           metres_between_carriers=math.sqrt(land_m2 / n), hauling_over_upkeep=[h / upkeep for h in hauling_day],
                           dry_food_mt_yr=[w / (FOOD_KJ_PER_G * 1e6) / 1e9 for w in whole], dry_food_g_per_m2_of_land_yr=food,
                           food_over_plant_growth=[food[0] / max(plant.values()), food[1] / min(plant.values())])
    fed = [r['as_grazers'] for r in sky['food_web']['flyers_fed'] if r['flyer'] == FLYER and r['edge'] == 'low' and abs(r['cover'] - 0.01) < 1e-9][0]
    return dict(glide_ratio=glide, muscle_efficiency=muscle, flyer=flyer, cost_kj_per_kg_delivered=cost, numbers=numbers,
                the_lands_plant_growth_g_dry_m2_yr=plant,
                food_over_plant_growth=[middle['dry_food_g_per_m2_of_land_yr'][0] / max(plant.values()), middle['dry_food_g_per_m2_of_land_yr'][1] / min(plant.values())],
                flyers_by_kind_of_land_at_1000_trips_and_0p6=by_kind, by_size_at_1000_trips_and_0p6=sizes,
                slow_and_fast_erosion_scale_the_numbers_by=rock['setting']['slow_and_fast_against_central'],
                sky_ecology_for_comparison=dict(flyers_the_sky_feeds_at_1_percent_cover=fed,
                                                flyers_to_return_the_lands_dissolved_phosphorus=next(t['flyers'] for t in sky['phosphorus']['sea_to_land']['targets']
                                                                                                    if t['target'] == 'dissolved_loss' and t['flyer'] == FLYER and t['mass_kg'] == flyer['mass_kg'])))


def staying_put(rock, terrain, job, carriers):
    lowering = job['mean_lowering_m_per_myr'] * 1e-6
    slow, fast = rock['setting']['slow_and_fast_against_central']
    kappa = terrain['small_craters']['kappa_m2_yr']
    reach = {f'creep {k:g} m2/yr': {name: CA.reach_for_slope(REPOSE_SLOPE, lowering * f, k) for name, f in (('slow', slow), ('central', 1.0), ('fast', fast))}
             for k in kappa}
    stone = terrain['upkeep']['built_land']['seas_stone_mt_yr']
    rock_mt = [job['gt_yr'][k] * 1e3 for k in ('fast', 'central', 'slow')]
    lime = stone['calcite'] + stone['magnesite']
    middle = carriers['numbers']['1000 trips, load 0.6 of body']
    nitrogen_mt = [middle['dry_food_mt_yr'][0] * FOOD_N[0], middle['dry_food_mt_yr'][1] * FOOD_N[1]]
    calcite_mt = [CA.calcite_by_urea_kg(n) for n in nitrogen_mt]
    return dict(repose_slope=REPOSE_SLOPE,
                how_far_each_side_a_line_of_returned_rock_can_feed_m=reach,
                slope_needed_if_only_main_crests_are_fed={f'{half:g} m to each side': {f'creep {k:g} m2/yr': CA.crest_slope(lowering, half, k) for k in kappa}
                                                         for half in (100.0, 1000.0, 5000.0)},
                the_rock_that_leaves_in_solution_mt_yr=stone,
                lime_over_the_particles_by_mass=dict(zip(('fast', 'central', 'slow'), (lime / r for r in rock_mt))),
                opal_over_the_particles_by_mass=dict(zip(('fast', 'central', 'slow'), (stone['opal'] / r for r in rock_mt))),
                the_carriers_nitrogen_mt_yr=nitrogen_mt, calcite_it_could_lay_by_the_urea_route_mt_yr=calcite_mt,
                that_calcite_over_the_rock_returned=[c / rock_mt[1] for c in calcite_mt])


def the_loop(rock, carriers):
    """The resting state and the recovery of carriers that multiply where sediment gathers, with the supply the central
    rock and the duty a thousand trips at 0.6 of body weight, against how slowly the lakes fill."""
    supply = rock['rock']['kg_yr_central']
    duty = 1000 * 0.6
    fill = rock['lakes']['years_to_fill_by_share_of_area']
    out = {}
    for name, s in LOOPS.items():
        half = s['half_years'] * supply
        deposit, flyers = CA.loop_equilibrium(supply, duty, s['growth'], s['death'], half)
        eig = CA.loop_eigenvalues(supply, duty, s['growth'], s['death'], half)
        years = 40.0 / max(-eig.real.max(), 1e-9)
        t, dep, car, ret = CA.loop_run(supply, duty, s['growth'], s['death'], half, deposit, 0.5 * flyers, min(years, 20000.0), dt=0.1)
        t2, dep2, car2, ret2 = CA.loop_run(2.0 * supply, duty, s['growth'], s['death'], half, deposit, flyers, min(years, 20000.0), dt=0.1)
        out[name] = dict(setting=s, standing_deposit_years_of_supply=deposit / supply, carriers_live_mass_mt=flyers / 1e9,
                         carriers_work_at_share_of_full_duty=s['death'] / s['growth'],
                         eigenvalues_per_yr=[[float(e.real), float(e.imag)] for e in eig], settles_in_years=1.0 / max(-eig.real.max(), 1e-9),
                         after_half_the_carriers_are_lost=dict(years_until_the_return_is_within_5_percent=CA.years_until_within(t, ret, supply),
                                                               most_sediment_left_standing_years_of_supply=float(dep.max() / supply)),
                         after_the_erosion_doubles=dict(years_until_the_return_is_within_5_percent=CA.years_until_within(t2, ret2, 2.0 * supply),
                                                        carriers_end_at_times_the_start=float(car2[-1] / flyers),
                                                        standing_deposit_ends_at_years_of_the_old_supply=float(dep2[-1] / supply)))
    out['the_lakes_for_scale'] = dict(years_for_a_hundredth_of_their_area_to_fill=fill['0.01'], years_for_a_tenth=fill['0.1'], years_for_half=fill['0.5'])
    return out


def the_waters(rock, reach, terrain):
    fill = rock['lakes']['years_to_fill_by_share_of_area']
    lakes = {'nothing returned': {f'{t:g}': CA.open_after(t, fill, 0.0) for t in SPANS_YR}}
    for texture, row in reach['all_the_land']['reach'].items():
        # of what the lakes keep, the share each carrier reaches: the dry hollows' part is not the lakes'
        s = row
        for who in ('flyers_and_waders', 'with_divers_at_the_delta_fronts', 'with_divers_on_the_floors_within_a_dive'):
            hollow = reach['all_the_land']['first_settles']['dry_hollow']
            share = (s[who] - hollow) / (1.0 - hollow)
            lakes[f'{texture}, {who}'] = dict(share_of_the_lakes_rock_returned=share, **{f'{t:g}': CA.open_after(t, fill, share) for t in SPANS_YR})
    stone_km3 = terrain['upkeep']['built_land']['seas_stone_total_km3_yr']
    volume = sum(s['volume_km3'] for s in rock['seas'].values())
    rock_km3 = sum(s['rock_share_central'] for s in rock['seas'].values()) * rock['rock']['kg_yr_central'] / 1000.0 / 1e9
    seas = dict(volume_km3=volume, the_seas_own_stone_km3_yr=stone_km3, years_for_their_own_stone_to_fill_them=volume / stone_km3,
                the_lands_rock_km3_yr_at_1000_kg_m3=rock_km3, years_on_rock_and_stone_together=volume / (stone_km3 + rock_km3),
                each={name: dict(years_on_the_lands_rock=s['years_to_fill_on_the_lands_rock']['1000 kg/m3'], floor_within_a_dive=s['floor'][f'under {DIVE_M:g} m'])
                      for name, s in rock['seas'].items()})
    return dict(lake_area_still_open=lakes, seas=seas)


def wild_and_managed(job, terrain):
    """What the designed upkeep of the hybrid and used land would carry if machines did it at the terrain study's
    costs, beside the wild land's rock that the carriers take."""
    r = terrain['upkeep']['rates']['50 m/Myr']
    per_gt = [(r['lift_gw_in_place'] + r['haulage_gw_at_rail_energy'] + f) / r['rock_gt_yr'] for f in r['lake_dredging_fuel_gw']]
    out = {}
    for key, kk in job['by_kind_of_land'].items():
        managed = kk['hybrid']['gt_yr'] + kk['peoples']['gt_yr']
        # the register's row of 22 November: half of the open land people do not fill goes to the hybrid land, half to the wild
        spare, spare_km2 = kk['open_beyond_need']['gt_yr'], kk['open_beyond_need']['land_km2']
        land = sum(kk[n]['land_km2'] for n in ('wild', 'hybrid', 'peoples', 'open_beyond_need'))
        shared = dict(wild_gt_yr=kk['wild']['gt_yr'] + UNFILLED_TO_THE_WILD * spare,
                      hybrid_gt_yr=kk['hybrid']['gt_yr'] + (1.0 - UNFILLED_TO_THE_WILD) * spare, peoples_gt_yr=kk['peoples']['gt_yr'],
                      share_of_land=dict(wild=(kk['wild']['land_km2'] + UNFILLED_TO_THE_WILD * spare_km2) / land,
                                         hybrid=(kk['hybrid']['land_km2'] + (1.0 - UNFILLED_TO_THE_WILD) * spare_km2) / land,
                                         peoples=kk['peoples']['land_km2'] / land))
        shared['machines_for_the_hybrid_and_peoples_land_gw'] = [(shared['hybrid_gt_yr'] + shared['peoples_gt_yr']) * x for x in per_gt]
        out[key] = dict(what=kk['what'], wild_gt_yr=kk['wild']['gt_yr'], open_beyond_need_gt_yr=spare,
                        hybrid_and_peoples_gt_yr=managed, machines_for_the_hybrid_and_peoples_land_gw=[managed * x for x in per_gt],
                        with_the_unfilled_land_shared=shared)
    out['machines_gw_per_gt'] = per_gt
    return out


def build():
    rock, sky, deltas, terrain, biomass, phosphorus = (load(p) for p in INPUTS)
    shares = delta_shares(deltas)
    job = the_job(rock, phosphorus)
    reach = the_reach(rock, shares)
    carriers = the_carriers(rock, sky, job, biomass)
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.run',
                      code_sha256={p.name: sha256(p) for p in (Path(__file__), HERE / 'carriers.py')}),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in INPUTS},
        evidence=('A sizing on committed products: the rock and its way back on the design surface (rock.py), the sky '
                  'ecology study\'s flyer and flight rule, the river deltas screen\'s shares, the terrain study\'s creep '
                  'and stone and the biomass comparison\'s plant growth. No animal carries sediment kilometres uphill; '
                  'loads, trips, the reach of a dive, what an animal digests of a sediment and the loop\'s rates are '
                  'design guesses, named in the setting. Flight with a load at lunar gravity is not measured. Not validated.'),
        reading_rule=('job: the rock to return and whose land it leaves. reach: the shares of the rock a flyer or '
                      'wader, a diver at the delta fronts and a diver on the floors within a dive can get at, with the '
                      'deltas\' soil from clay to sand; the dry hollows\' part is all within a flyer\'s reach. carriers: '
                      'kJ of metabolism for each kg delivered, pairs from the best glide and muscle to the worst, on '
                      'the Moon and for the same way on Earth; numbers of flyers of the albatross scale for all the '
                      'land\'s central rock by trips and load, with what the hauling adds to their upkeep and what '
                      'they eat. staying_put: how far each side a line of returned rock can feed ground losing rock at '
                      'each rate before creep needs more than the repose slope; lime and opal against the particles by '
                      'mass; the calcite the carriers\' own nitrogen could lay. loop: the resting deposit and carriers, '
                      'the rates at which a disturbance dies away, and recoveries. waters: the share of the lakes\' '
                      'area still open after each span for each reach.'),
        setting=dict(flyer=FLYER, smaller_carriers_kg=SIZES_KG, payloads_over_body=list(PAYLOADS), trips_a_year=list(TRIPS), earth_payloads=list(EARTH_PAYLOADS),
                     cruise_m_s=CRUISE_M_S, handling_s=HANDLING_S, dive_m=DIVE_M, earth_dives_m=EARTH_DIVES_M,
                     delta_organic_carbon=list(DELTA_ORGANIC_C), kj_per_g_c=KJ_PER_G_C, digestible=list(DIGESTIBLE),
                     food_kj_per_g=FOOD_KJ_PER_G, food_nitrogen=list(FOOD_N), loose_soil_kg_m3=LOOSE_SOIL_KG_M3,
                     repose_slope=REPOSE_SLOPE, loops=LOOPS, delta_shares=shares),
        job=job, reach=reach, carriers=carriers, staying_put=staying_put(rock, terrain, job, carriers),
        loop=the_loop(rock, carriers), waters=the_waters(rock, reach, terrain), wild_and_managed=wild_and_managed(job, terrain))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    j, c = out['job'], out['carriers']
    print(f'{j["gt_yr"]["central"]:.2f} Gt of rock a year; {c["numbers"]["1000 trips, load 0.6 of body"]["flyers_for_all_the_land_central"] / 1e6:.2f} '
          f'million flyers at a thousand trips and 0.6 of body; wrote {OUT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
