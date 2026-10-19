"""The resort at 301 K across the domains: where the branches' work meets, set side by side and worked out together.

    python -m research.studies.resort_301k.run        # under a second

Arithmetic on committed products; no model runs. It binds the lakes' and seas' methane at 301 K, the lakes' anoxia
(the phosphorus lines), the domain synthesis's gas ledger, the sea phosphorus screen, the food systems, the
heat-levers screen, the target climates, the coast and island water study and the power screen, and writes
results/resort_301k.json. Four questions: what each gas pair asks of the waters' methane; what phosphorus the deep
lakes can take from people and still stay nutrient-poor, and what the seas get; the resort's dimming with the levers
worth pursuing; and what watering the deserts would take, as a scale for the author's choice.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'resort_301k.json'
SCHEMA = 'terluna.research.resort-301k/1'
CODE = ('research/studies/resort_301k/run.py',)
METHANE = 'research/studies/lake_methane/results/lake_methane.json'
ANOXIA = 'research/studies/lake_anoxia/results/lake_anoxia.json'
LEDGERS = 'research/studies/domain_synthesis/results/ledgers.json'
PHOSPHORUS = 'resources/results/sea_phosphorus.json'
FOOD = 'research/studies/food_systems/results/food_systems.json'
LEVERS = 'research/studies/heat_levers/results/heat_levers.json'
TARGETS = 'research/studies/climate_targets/results/climate_targets.json'
COAST = 'research/studies/coast_water/results/coast_water.json'
POWER = 'provisioning/results/power.json'
SEA_OXYGEN = 'research/studies/sea_oxygen/results/sea_oxygen.json'
INPUTS = {METHANE: 'terluna.research.lake-methane/1', ANOXIA: 'terluna.research.lake-anoxia/2',
          LEDGERS: 'terluna.research.domain-synthesis-ledgers/1', PHOSPHORUS: 'terluna.resources.sea-phosphorus/1',
          FOOD: None, LEVERS: None, TARGETS: 'terluna.research.climate-targets/5',
          COAST: 'terluna.research.coast-water/1', POWER: 'terluna.provisioning.power/1',
          SEA_OXYGEN: 'terluna.research.sea-oxygen/1'}

TARGET_K = 301.0
PAIRS_PPM = dict(held=1.87, low=0.3, grid_lowest=0.1)      # methane: the held pair, the levers' low pair, the grid's
SETTLING = ('earth', 'aggregates', 'stokes')                # lake_anoxia's settling: Earth's, sqrt(g/6), g/6
MOON_AREA_M2 = 4.0 * math.pi * 1737.4e3 ** 2
POINT_W_M2 = 2.254                                          # a point of dimming: 1% of the hot world's 225.4 W/m2
WATERED_SHARES = (0.01, 0.05, 0.10, 0.25)                   # shares of the desert land watered, for the scale
# design guesses for the watering scale, each named in the README
GREENING_ALBEDO_DROP = (0.10, 0.20)      # a watered garden or park against bare desert ground
LIFT_EFFICIENCY = 0.7                    # pumps
CONVEYANCE_KWH_M3 = (0.1, 0.3)           # pipes and channels to the gardens, beyond the lift
GRAVITY = 1.625                          # m/s2 at the lunar surface

EVIDENCE = ('Arithmetic on committed products, no model run: the waters\' methane traced at 301 K and the domain '
            'synthesis\'s gas ledger; the lakes\' anoxia study\'s phosphorus lines and the sea phosphorus screen\'s '
            'loads; the targets\' resort band and the heat levers\' ranking; the coast water study\'s deserts, '
            'demand and water table with the power screen\'s heat and the levers\' albedo kernels. Its own design '
            'guesses are named in design_guesses.')
READING_RULE = ('Pairs are [low, high]. methane gives, for each trophic state and methanotroph community, the waters\' '
                'Tg a year, the sources with the ledger\'s others cut 90%, the air the soils alone leave, and at '
                'each pair the removal and the share of the air its contactors pass a year. phosphorus gives the '
                'oligotrophic line by settling and the room under it for people\'s dissolved phosphorus, in kg a year '
                'and g a person a year at 4 and 8 billion, with the share of the design loss that may arrive dissolved '
                '[least, most]. '
                'dimming gives the resort\'s at 301 K with the levers\' first four taken off. watering gives, for each '
                'share of the deserts watered and each interpolation form, the land, the water a year, the lift, the '
                'power, its heat in points and the points the greener ground costs.')
DESIGN_GUESSES = dict(
    greening_albedo_drop=dict(range=list(GREENING_ALBEDO_DROP), basis='a watered garden or park against bare desert '
                                                                       'ground'),
    lift_efficiency=dict(value=LIFT_EFFICIENCY, basis='pumps'),
    conveyance_kwh_m3=dict(range=list(CONVEYANCE_KWH_M3), basis='pipes and channels to the gardens, beyond the lift'))


def digest(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load():
    data = {p: json.loads((ROOT / p).read_text()) for p in INPUTS}
    for p, schema in INPUTS.items():
        if schema and data[p].get('schema') != schema:
            raise SystemExit(f'{p} is {data[p].get("schema")}, expected {schema}')
    return data


def methane(data):
    """Each pair's methane removal with the waters traced at 301 K: the ledger's other sources cut by biology's 90%,
    the waters at their traced flux, the soils' sink and the air's methane per ppm, as lake_methane reckons the held
    pair."""
    air = data[METHANE]['air']['lakes_and_seas']
    gases = data[LEDGERS]['gases']['methane']
    tg_per_ppm = gases['tg_per_ppm']
    sink = air['soil_sink_tg_yr_per_ppm']
    others = [0.1 * (c['sources_tg_yr']['total'] - air['ledger_lakes_tg_yr']) for c in gases['central']]
    out = dict(soil_sink_tg_yr_per_ppm=sink, tg_per_ppm=tg_per_ppm, other_sources_after_cut_tg_yr=others,
               pairs_ppm=PAIRS_PPM, waters={})
    for trophic, by_community in data[METHANE]['air']['lakes_and_seas_tg_yr'].items():
        rows = {}
        for community, waters in by_community.items():
            sources = [o + waters for o in others]
            row = dict(waters_tg_yr=waters, sources_tg_yr=sources, steady_ppm_by_soils=[s / sink for s in sources])
            for name, ppm in PAIRS_PPM.items():
                removal = [max(s - sink * ppm, 0.0) for s in sources]
                row[name] = dict(removal_tg_yr=removal,
                                 air_passed_per_year_share=[r / (ppm * tg_per_ppm) for r in removal])
            rows[community] = row
        out['waters'][trophic] = rows
    return out


def phosphorus(data):
    """The deep lakes' room for people's dissolved phosphorus under the oligotrophic line, by settling, against the
    plan's losses; the seas' phosphorus as the lakes leave it."""
    sp, la, food = data[PHOSPHORUS], data[ANOXIA], data[FOOD]
    area = sp['partition']['lake_area_m2']
    land_kg = sp['routes']['bioavailable_input_kg_yr']['lakes']['base']
    land_g = [x * 1e3 / area for x in land_kg]
    people = sp['people']['working_range_people']
    lost = sp['people']['ecology_design_target']['lost_kg_per_person']
    perm = la['phosphorus']['permissible_load_g_m2_yr']
    glass = food['systems']['fog_desert_glass_vented']['phosphorus_mt_year']['discharged']   # Mt a year a billion
    room = {}
    for s in SETTLING:
        line = perm[s]['oligotrophic']
        r = [(line - g) * area / 1e3 for g in land_g]                      # kg a year of bioavailable P, by land loss
        room[s] = dict(oligotrophic_line_g_m2_yr=line, room_kg_yr=r,
                       per_person_g_yr={f'{int(n / 1e9)}_billion': [x / n * 1e3 for x in r] for n in people},
                       dissolved_share_of_the_loss_allowed={
                           f'{int(n / 1e9)}_billion': [r[1] / (n * lost[1]), r[0] / (n * lost[0])] for n in people})
    seas = sp['seas']['at_1000_years_mg_m3']
    return dict(lake_area_m2=area, land_bioavailable_g_m2_yr=land_g, people=people, lost_kg_per_person=lost,
                glass_discharge_kg_per_person=[x for x in glass], room=room,
                earth_practice_lost_kg_per_person=sp['people']['earth_practice_lost_kg_per_person'],
                seas_at_1000_years_mg_m3=seas, earth_ocean_mean_mg_m3=sp['seas']['earth_ocean']['mean_mg_m3'])


def dimming(data):
    """The resort's dimming at 301 K with the levers that ease the shield, and what the low pair's methane half asks
    of the contactors once the waters are traced."""
    row = next(r for r in data[TARGETS]['resort_band']['rows'] if r['ground_k'] == TARGET_K)
    lv = data[LEVERS]
    first = lv['ranking']['first_four_points']
    total = row['dimming']['with_the_heat_beyond_sunlight_percent']
    own = row['dimming']['own_percent']
    held = lv['gases']['costs']['held']['methane cut 0.9, N2O cut 0.95']
    tw_per_share = [f / s for f, s in zip(held['methane_fan_tw'], held['methane_air_share_per_year'])]
    low = lv['gases']['options']['low pair, 0.3 and 0.03 ppm']
    return dict(resort_own_percent=own, resort_in_all_percent=total, levers_first_four_points=first,
                levers=[dict(rank=x['rank'], lever=x['lever'], points=x['points']) for x in lv['ranking']['levers']],
                resort_in_all_with_the_first_four_percent=[total[0] - first[1], total[1] - first[0]],
                low_pair_points=low['points'], methane_fan_tw_per_share_of_air_a_year=tw_per_share)


def watering(data):
    """What watering the deserts at 301 K would take, per share of the desert land: water, the lift and its power, the
    heat it leaves and the sunlight a greener ground keeps."""
    coast, lv, power = data[COAST], data[LEVERS], data[POWER]
    demand = coast['demand']['desert_uptake_mm_day']
    deserts = {}
    for form, mask in (('additive', 'rain_mask'), ('log_linear', 'log_rain_mask')):
        cell = coast['water_table'][mask]['every_cell']
        deserts[form] = dict(land_km2=cell['land_km2'], median_depth_m=cell['median_depth_m'],
                             bands={b: dict(land_km2=v['land_km2'], median_depth_m=v['median_depth_m'])
                                    for b, v in cell['bands'].items() if v['land_km2']})
    heat_per_w = power['absorption']['link']['ground_heat_per_w']
    # the levers' desert albedo: W/m2 for each unit of albedo on each unit of the Moon's sunlit desert share
    per_unit = []
    for name in ('hot world', 'lock'):
        d = lv['albedo']['deserts'][name]
        share = d['sunlit_desert_share']['every_cell']
        per_unit += [d['w_m2'][0] / (0.10 * share), d['w_m2'][1] / (0.30 * share)]
    per_unit = [min(per_unit), max(per_unit)]
    sunlit_over_area = []
    for name in ('hot world', 'lock'):
        d = lv['albedo']['deserts'][name]
        sunlit_over_area.append(d['sunlit_desert_share']['every_cell'] / d['desert_share_of_moon']['every_cell'])
    rows = []
    for share in WATERED_SHARES:
        row = dict(share_of_deserts=share)
        for form, d in deserts.items():
            km2 = d['land_km2'] * share
            water_km3 = km2 * 1e6 * demand * 1e-3 * 365.25 / 1e9
            lifts = [d['bands'].get('30-45', {}).get('median_depth_m') or d['median_depth_m'], d['median_depth_m']]
            kwh_m3 = [1000.0 * GRAVITY * h / LIFT_EFFICIENCY / 3.6e6 + c for h, c in zip(lifts, CONVEYANCE_KWH_M3)]
            gw = [water_km3 * 1e9 * k * 3.6e6 / (365.25 * 86400) / 1e9 for k in kwh_m3]  # m3/yr * kWh/m3 -> GW
            heat_points = [g * 1e9 * heat_per_w / MOON_AREA_M2 / POINT_W_M2 for g in gw]
            sunlit = [km2 * 1e6 / MOON_AREA_M2 * s for s in sorted(sunlit_over_area)]
            albedo_w = [sunlit[0] * GREENING_ALBEDO_DROP[0] * per_unit[0],
                        sunlit[1] * GREENING_ALBEDO_DROP[1] * per_unit[1]]
            row[form] = dict(land_km2=km2, water_km3_yr=water_km3, lift_m=lifts, kwh_per_m3=kwh_m3, power_gw=gw,
                             heat_points=heat_points, albedo_w_m2=albedo_w,
                             albedo_points=[w / POINT_W_M2 for w in albedo_w])
        rows.append(row)
    return dict(demand_mm_day=demand, deserts=deserts, heat_per_w=heat_per_w,
                albedo_w_m2_per_unit_albedo_per_unit_sunlit_share=per_unit,
                sunlit_over_area=sorted(sunlit_over_area),
                rivers_km3_yr=data[PHOSPHORUS]['seas']['discharge_m3_yr'] / 1e9,
                rows=rows)


def build():
    data = load()
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, design_guesses=DESIGN_GUESSES,
               target_k=TARGET_K, methane=methane(data), phosphorus=phosphorus(data), dimming=dimming(data),
               watering=watering(data))
    out['producer'] = dict(files={p: digest(p) for p in CODE if (ROOT / p).exists()},
                           inputs={p: digest(p) for p in INPUTS})
    return out


def main() -> int:
    res = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1) + '\n')
    poor = res['methane']['waters']['oligotrophic']['earth']
    ppm = poor['steady_ppm_by_soils']
    resort = res['dimming']['resort_in_all_with_the_first_four_percent']
    print(f"wrote {OUT.relative_to(ROOT)}: nutrient-poor waters {poor['waters_tg_yr']:.2f} Tg a year, air "
          f"{ppm[0]:.2f}-{ppm[1]:.2f} ppm; resort in all with the levers {resort[0]:.1f}-{resort[1]:.1f}%")
    return 0


if __name__ == '__main__':
    sys.exit(main())
