"""Where a living return of rock could live and how far it could carry: a first look by the land's dryness.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.rock_return.dryness   # about three minutes on one core, 6 GB; reads the research drive

The runs of this study set the returned rock down within a range of the place it settled, 20 km or 50 km, and take
carriers to live wherever rock is shed. Neither is derived from anything. A carrier's range has two senses, how far
a load is carried and where the carriers can live, and this module asks what would limit each before any ecology of
carriers and their partners is stated:

1. By dryness, in the sand study's classes of rain over demand on the 1.9 km grid: the land, the rock it sheds and
   where that rock settles, the plants' growth there to feed carriers, the lasting sand that gathers at the ends of
   drainage, and the wind's own sand. Plant growth by dryness is no project product: the land's mean, the biomass
   comparison's, is spread here by two shapes as a bracket, the globe's plant cover curve (a display rule) and the
   rain branch of Earth's Miami rule (for its shape only).
2. What a longer haul costs a carrier in energy, flying hours and food, by the study's flight rule, for one typical
   trip at the mean lift.
3. How dry the landscape model's regions are, so that a run's "dry" hollows are not read as an arid climate.

Writes results/dryness.json.
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import sand as S
from research.studies.atlas_round import sand_seas as SS
from research.studies.rock_return import carriers as CA
from research.studies.rock_return import evolve as EV
from research.studies.rock_return import rock as RK
from research.studies.rock_return import run as RU
from shared.constants import JULIAN_YEAR_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'dryness.json'
SCHEMA = 'terluna.research.rock-return-dryness/1'
SAND_SEAS = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'sand_seas.json'
HAULS_M = (9400.0, 12500.0, 31000.0, 100000.0)   # the straight way back; the mean way of a bell with 86% inside 20 km and inside 50 km (1.25 of its width); a long commute
LIFT_M = 1000.0                                  # the rock's mean climb back (rock.py: 1.0 km)
LOAD = 0.6                                       # of the body, the sizing's middle load
TRIPS = 1000.0                                   # loaded trips a year, the sizing's middle
BUILD_YEARS = '5000'                             # the made seas' span of building in the sand seas study's plan
MIAMI_G_M2 = 3000.0                              # Lieth's Miami rule, rain branch: g of dry matter per m2 a year on Earth, 3000 (1 - exp(-0.000664 P)), P in mm (recalled)
MIAMI_PER_MM = 0.000664
SEDIMENT_KINDS = ('lake deltas', 'closed hollows', 'sea shores')


def miami(rain_mm):
    """Plant growth by rain on Earth (g of dry matter per m2 a year); used here for its shape alone."""
    return MIAMI_G_M2 * (1.0 - np.exp(-MIAMI_PER_MM * np.asarray(rain_mm, float)))


def over_the_mean(shape, area):
    """A field scaled so that its mean over the area is one: the land's mean plant growth spread by a shape."""
    shape, area = np.asarray(shape, float), np.asarray(area, float)
    return shape * area.sum() / (shape * area).sum()


def food_multiple(rock_per_m2, mean_rock_per_m2, growth_over_mean):
    """The return's food as a share of the plants' growth on some ground, over the same share on the land's mean:
    the rock to return there against the mean, over the growth there against the mean."""
    return rock_per_m2 / mean_rock_per_m2 / growth_over_mean


def haul_rows(total_kg, land_m2, plant_g, glide, muscle, sizes, upkeep_j_day, hauls_m=HAULS_M, lift_m=LIFT_M, load=LOAD, trips=TRIPS,
              speed=RU.CRUISE_M_S, handling=RU.HANDLING_S, food_kj_per_g=RU.FOOD_KJ_PER_G):
    """For each mean haul: the energy for each kilogram delivered (cheap and dear ends of the flight rule, the way
    home glided as far as the height gained allows), the flying hours for the year's trips, and for each size of
    carrier its number, live mass and food, both at the same trips and at the flying hours of the first haul, where
    longer trips mean fewer of them and more carriers."""
    first = CA.hours_a_year(trips, hauls_m[0], speed, handling)
    out = {}
    for haul in hauls_m:
        cost = [CA.delivered_cost_j_per_kg(haul, lift_m, g, e, load, min(1.0, lift_m * g / haul)) for g, e in ((glide[1], muscle[1]), (glide[0], muscle[0]))]
        hours = CA.hours_a_year(trips, haul, speed, handling)
        fewer = trips * first / hours
        row = dict(kj_per_kg_delivered=[x / 1e3 for x in cost], hours_a_year=hours, share_of_the_year=hours / (24.0 * JULIAN_YEAR_DAYS),
                   trips_in_the_first_hauls_hours=fewer, by_size={})
        for name, mass in sizes.items():
            sized = {}
            for label, n_trips in (('at_the_same_trips', trips), ('in_the_same_hours', fewer)):
                count = CA.carriers_for(total_kg, mass, load, n_trips)
                whole = [count * (upkeep_j_day[name] * JULIAN_YEAR_DAYS + n_trips * load * mass * x) for x in cost]            # J a year
                food = [w / (food_kj_per_g * 1e6) / 1e9 for w in whole]                                                       # Mt of dry food a year
                sized[label] = dict(carriers=count, live_mass_mt=count * mass / 1e9, dry_food_mt_yr=food,
                                    of_the_lands_plant_growth=[food[0] * 1e12 / (land_m2 * plant_g[1]), food[1] * 1e12 / (land_m2 * plant_g[0])])
            row['by_size'][name] = sized
        out[f'{haul / 1e3:g} km'] = row
    return out


def regions_dryness(aridity, rain, cover, water):
    """The dry ground of each region placement.py cut, by the sand study's dryness: its classes and its middle values."""
    with np.load(EV.REGIONS, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
    n = meta['cells']
    out = {}
    for name, r in meta['regions'].items():
        rr = np.arange(r['row0'], r['row0'] + n)
        cc = np.arange(r['col0'], r['col0'] + n) % aridity.shape[1]
        take = lambda a: a[np.ix_(rr, cc)]                                               # noqa: E731
        land = ~take(water)
        a = take(aridity)[land]
        k = S.class_of(a)
        out[name] = dict(dry_ground_share=float(land.mean()), aridity=[float(np.percentile(a, q)) for q in (10, 50, 90)],
                         rain_mm_yr=[float(np.percentile(take(rain)[land], q)) for q in (10, 50, 90)],
                         plant_cover=[float(np.percentile(take(cover)[land], q)) for q in (10, 50, 90)],
                         classes={S.CLASS_NAMES[q]: float((k == q).mean()) for q in range(6)})
    return out


def build(log=print):
    started = time.time()
    seas = json.loads(SAND_SEAS.read_text())
    placement = json.loads(EV.PLACEMENT.read_text())
    F = S.fields(log)
    src = SS.sources(F)
    shape = F.aridity.shape
    aridity = F.aridity.astype(np.float32).copy()
    cover = F.cover.astype(np.float32).copy()
    sand_water = F.water.copy()
    typical = np.interp(np.abs(F.lat_f), S.EROSION_CENTRES_DEG, F.band_rain['mean_mm_yr'])
    rain_rule = F.rain_rule.astype(np.float64)
    rain = (typical[:, None] * rain_rule ** (1.0 / S.EROSION_RAIN_POWER)).astype(np.float32)     # the rain behind the sand study's rule, mm a year
    del F
    gc.collect()
    regions = regions_dryness(aridity, rain, cover, sand_water)
    aridity, cover, rain, rain_rule, sand_water = (a.ravel() for a in (aridity, cover, rain, rain_rule, sand_water))
    cls = S.class_of(aridity)
    names = S.CLASS_NAMES
    log(f'the sand study\'s fields and the ends of drainage, {time.time() - started:.0f} s')

    c = RK.LO.prepare(log)
    f = RK.fields(c, log, started)
    rock, eroding, first, area, n = f['rock'], f['eroding'], f['first'], f['area'], f['n']
    donor, body = c['donor'], c['body']
    spills_at = np.zeros(n, bool)
    spills_at[body >= 0] = c['spills'][body[body >= 0]]
    place = np.where(~donor[first], 0, np.where(body[first] >= 0, np.where(spills_at[first], 1, 2), 3)).astype(np.int8)
    agree = float((np.isfinite(f['water']) == sand_water).mean())
    total, land = f['total'], float(area[eroding].sum())
    cells = np.flatnonzero(eroding)
    a, m = area[cells], rock[cells]
    k_shed, k_settles, where = cls[cells], cls[first[cells]], place[cells]
    with_rain = m * rain_rule[cells]
    with_rain *= total / with_rain.sum()                          # the same total, spread also by the sand study's rule of rain
    mean_rock = total / land
    growth = dict(by_cover=over_the_mean(cover[cells], a), by_rain=over_the_mean(miami(rain[cells]), a))
    plant = RU.load(RU.BIOMASS)['per_area']['cases']
    plant_g = [plant[key]['land_g_c_m2_year'] / RU.CARBON_IN_DRY_BIOMASS for key in ('low', 'high')]
    flows = np.zeros((6, 6))
    np.add.at(flows, (k_shed, k_settles), m)

    by = {}
    for k, name in enumerate(names):
        s, e = k_shed == k, k_settles == k
        ground = float(a[s].sum())
        if ground == 0:
            continue
        over = {key: float((g[s] * a[s]).sum() / ground) for key, g in growth.items()}
        shed, shed_rain = float(m[s].sum()), float(with_rain[s].sum())
        ends = cls[src['cell']] == k
        by[name] = dict(
            land_share=ground / land, aridity=float((aridity[cells][s] * a[s]).sum() / ground), rain_mm_yr=float((rain[cells][s] * a[s]).sum() / ground),
            plant_cover=float((cover[cells][s] * a[s]).sum() / ground),
            rock_shed_gt_yr=shed / 1e12, rock_shed_share=shed / total, rock_shed_g_m2_yr=shed / ground * 1e3,
            rock_shed_with_the_rule_of_rain_g_m2_yr=shed_rain / ground * 1e3,
            rock_settling_here_gt_yr={RK.PLACES[p]: float(m[e & (where == p)].sum() / 1e12) for p in range(4)},
            of_the_rock_settling_here_shed_in_this_class=float(flows[k, k] / max(flows[:, k].sum(), 1e-30)),
            plant_growth_over_the_lands_mean=over,
            plant_growth_g_dry_m2_yr={key: [p * v for p in plant_g] for key, v in over.items()},
            food_over_the_lands_mean_share=dict(
                {key: food_multiple(shed / ground, mean_rock, v) for key, v in over.items()},
                **{f'{key}_with_the_rule_of_rain': food_multiple(shed_rain / ground, mean_rock, v) for key, v in over.items()}),
            lasting_sand_at_ends_km3_yr={SEDIMENT_KINDS[q]: float(src['m3_yr'][ends & (src['kind'] == q)].sum() / 1e9) for q in range(3)},
            winds_sand_km2=seas['map']['open']['by_class_km2'][name])
        log(f'{name}: {100 * by[name]["land_share"]:.1f}% of the land, {by[name]["rock_shed_g_m2_yr"]:.0f} g of rock a m2, growth {over}')

    # the made seas: the lasting sand at the ends of drainage within each sea's radius of building, by the dryness of the end
    made = {}
    for name, sea in seas['seas'].items():
        plan = sea['haul']['within_years'][BUILD_YEARS]
        near = (SS.distance_km(src['lat'], src['lon'], *sea['centre_today']) <= plan['radius_km']) & (src['kind'] < 2)
        made[name] = dict(sand_km3=sea['sand_km3'], radius_km=plan['radius_km'], need_km3_yr=sea['sand_km3'] / plan['years'],
                          at_ends_within_the_radius_km3_yr={names[k]: float(src['m3_yr'][near & (cls[src['cell']] == k)].sum() / 1e9) for k in range(6)})
    sand_all = sum(sea['sand_km3'] for sea in seas['seas'].values())
    wind = seas['supply']['sand_study_capped_rule']
    sand = dict(made_seas=made,
                made_seas_rock_over_the_return=sand_all * 1e9 * S.BULK_DENSITY / float(BUILD_YEARS) / total,
                winds_take_km3_yr=[wind['wind_offer_central']['own_ground_km3'] / seas['years'], wind['wind_offer_opensand']['own_ground_km3'] / seas['years']],
                winds_sand_km2=seas['map']['open']['sand_km2'])

    sky = RU.load(RU.SKY)
    flyer = RU.flyer_of(sky)
    sizes = dict(RU.SIZES_KG, albatross_scale=flyer['mass_kg'])
    upkeep = {name: RU.FW.bird_fmr_j_day(mass) for name, mass in sizes.items()}
    hauls = haul_rows(total, land, plant_g, sky['food_web']['ranges']['glide_ratio'], sky['food_web']['ranges']['muscle_efficiency'], sizes, upkeep)

    code = ('dryness.py', 'rock.py', 'carriers.py', 'run.py')
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.dryness', code_sha256={p: EV.sha256(HERE / p) for p in code},
                      sand_study=seas['producer']['files']),
        inputs={str(p.relative_to(ROOT)): EV.sha256(p) for p in (SAND_SEAS, RU.BIOMASS, RU.SKY, EV.PLACEMENT, ROOT / RK.LO.ROUTING_RECORD)},
        ignored_inputs=dict(seas['ignored_inputs'], **{str(EV.REGIONS): placement['product']['sha256']}),
        evidence=('A first look, nothing in it observed. The dryness, the rain and the plant cover are the sand study\'s '
                  'fields on the 1.9 km grid, from one climate run; the rock is this study\'s, spread by relief and '
                  'not by rain; the lasting sand at the ends of drainage is the sand seas study\'s. Plant growth by '
                  'dryness is no project product: the biomass comparison\'s mean for the land is spread by the '
                  'globe\'s plant cover curve, a display rule, and by the rain branch of Earth\'s Miami rule, as a '
                  'bracket. A haul is one typical trip at the rock\'s mean climb; loads, trips and food are design '
                  'guesses. No carrier, plant or cover is designed here, so this says what energy and food allow and '
                  'nothing of what an ecology of them would do.'),
        reading_rule=('by_dryness: for each class of rain over demand, its share of the dry land, the rock it sheds '
                      '(by this study\'s rule, and with the sand study\'s rule of rain beside it), the rock that '
                      'settles in it by kind of place, the plants\' growth there over the land\'s mean by each shape, '
                      'food_over_the_lands_mean_share (the return\'s food as a share of the plants\' growth there, '
                      'over that share on the land\'s mean; multiply by the sizing\'s shares, 0.05-1.2% by the '
                      'carrier\'s size), the lasting sand its ends of drainage gather and the wind\'s own sand in the '
                      'sand study\'s open case. sand: each made sea\'s need a year against the lasting sand at the '
                      'ends within its radius of building, by the dryness of the end; the made seas\' rock over the '
                      'whole return while they are built; what the wind takes a year. hauls: for each mean haul, the '
                      'energy for a kilogram delivered, the flying hours for a thousand trips, and by size the '
                      'carriers and their food at the same trips and in the same flying hours. regions: the classes '
                      'of each region\'s dry ground and its middle dryness, rain and plant cover (10th, 50th, 90th '
                      'percentile).'),
        setting=dict(classes=list(names), class_edges=list(S.CLASS_EDGES), hauls_m=list(HAULS_M), lift_m=LIFT_M, load_over_body=LOAD, trips_a_year=TRIPS,
                     cruise_m_s=RU.CRUISE_M_S, handling_s=RU.HANDLING_S, food_kj_per_g=RU.FOOD_KJ_PER_G, sizes_kg=sizes,
                     miami=dict(g_m2=MIAMI_G_M2, per_mm=MIAMI_PER_MM), build_years=float(BUILD_YEARS)),
        land_m2=land, rock_kg_yr=total, mean_rock_g_m2_yr=mean_rock * 1e3, plant_growth_land_mean_g_dry_m2_yr=plant_g,
        the_two_studies_water_agrees=agree, by_dryness=by, rock_from_class_to_class_gt_yr=(flows / 1e12).tolist(), sand=sand, hauls=hauls,
        regions=regions, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build(log=lambda *a: print(*a, flush=True))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
