"""The resort at 301 K across the domains: where the branches' work meets, set side by side and worked out together.

    python -m research.studies.resort_301k.run        # under a second

Arithmetic on committed products; no model runs. It binds the lakes' and seas' methane at 301 K, the lakes by kind
and their anoxia, the seas' oxygen, the domain synthesis's gas ledger, the sea phosphorus screen, the food systems,
the heat-levers screen, the target climates, the coast and island water study and the power screen, and writes
results/resort_301k.json. Four questions: what each gas pair asks of the waters' methane; what phosphorus the lakes can
take from people and still stay nutrient-poor, and what the seas get; the resort's dimming with the levers worth
pursuing; and what watering the deserts would take, as a scale for the author's choice.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'resort_301k.json'
SCHEMA = 'terluna.research.resort-301k/1'
CODE = ('research/studies/resort_301k/run.py',)
METHANE = 'research/studies/lake_methane/results/lake_methane.json'
ANOXIA = 'research/studies/lake_anoxia/results/lake_anoxia.json'
KINDS = 'research/studies/lake_kinds/results/lake_kinds.json'
SEA_OXYGEN = 'research/studies/sea_oxygen/results/sea_oxygen.json'
LEDGERS = 'research/studies/domain_synthesis/results/ledgers.json'
PHOSPHORUS = 'resources/results/sea_phosphorus.json'
FOOD = 'research/studies/food_systems/results/food_systems.json'
LEVERS = 'research/studies/heat_levers/results/heat_levers.json'
TARGETS = 'research/studies/climate_targets/results/climate_targets.json'
COAST = 'research/studies/coast_water/results/coast_water.json'
POWER = 'provisioning/results/power.json'
INPUTS = {METHANE: 'terluna.research.lake-methane/1', ANOXIA: 'terluna.research.lake-anoxia/2',
          KINDS: 'terluna.research.lake-kinds/1', SEA_OXYGEN: 'terluna.research.sea-oxygen/1',
          LEDGERS: 'terluna.research.domain-synthesis-ledgers/1', PHOSPHORUS: 'terluna.resources.sea-phosphorus/1',
          FOOD: 'terluna.research.food-systems/1', LEVERS: 'terluna.research.heat-levers/2',
          TARGETS: 'terluna.research.climate-targets/5', COAST: 'terluna.research.coast-water/1',
          POWER: 'terluna.provisioning.power/1'}

TARGET_K = 301.0
PAIRS_PPM = dict(held=1.87, low=0.3, grid_lowest=0.1)      # methane: the held pair, the levers' low pair, the grid's
SETTLING = ('earth', 'aggregates', 'stokes')                # lake_anoxia's settling: Earth's, sqrt(g/6), g/6
KINDS_BY_AREA = ('self_ventilated', 'deep', 'giant')        # the lake kinds' census, which together cover the lakes
MOON_AREA_M2 = 4.0 * math.pi * MOON_RADIUS ** 2
WATERED_SHARES = (0.01, 0.05, 0.10, 0.25)                   # shares of the desert land watered, for the scale
GARDEN_KM2 = 1000.0                                         # the per-band rows' unit of watered land
# design guesses for the watering scale, each named in the README
GREENING_ALBEDO_DROP = (0.10, 0.20)      # a watered garden or park against bare desert ground
LIFT_EFFICIENCY = 0.7                    # pumps
CONVEYANCE_KWH_M3 = (0.1, 0.3)           # pipes and channels to the gardens, beyond the lift

EVIDENCE = ('Arithmetic on committed products, no model run: the waters\' methane traced at 301 K and the domain '
            'synthesis\'s gas ledger; the lakes\' oligotrophic lines by kind and the sea phosphorus screen\'s loads; '
            'the seas\' oxygen by trophic state; the heat levers\' first four by case; the coast water study\'s '
            'deserts, demand and water table with the power screen\'s heat and the levers\' albedo kernels. Its own '
            'design guesses are named in design_guesses.')
READING_RULE = ('Pairs are [low, high] unless named [larger, smaller]. methane gives, for each trophic state and '
                'methanotroph community, and for the methane share\'s 0.8 and Shih and Bouffard\'s mixing beside the '
                'central case, the waters\' Tg a year, the sources with the ledger\'s others cut 90%, the air the '
                'soils alone leave, and at each pair the removal, the share of the air its contactors pass a year and '
                'the years they take to pass the air\'s methane. phosphorus gives the oligotrophic line by lake kind '
                'and settling, the room under it for people\'s dissolved phosphorus [larger, smaller] at the land\'s '
                'low and high load, in kg a year and g a person a year at 4 and 8 billion, with the giants\' floors '
                'oxic and with them anoxic, and the share of the design loss that may arrive dissolved [least, most]. '
                'dimming gives the resort\'s at 301 K with the levers\' first four taken off, each case paired as the '
                'heat levers pair it. watering gives the scale by share of the deserts and, per 1,000 km2 of garden, '
                'by latitude band: the water a year, the lift, the power, its heat in points and the points the '
                'greener ground costs.')
DESIGN_GUESSES = dict(
    greening_albedo_drop=dict(range=list(GREENING_ALBEDO_DROP),
                              basis='a watered garden or park against bare desert ground'),
    lift_efficiency=dict(value=LIFT_EFFICIENCY, basis='pumps'),
    conveyance_kwh_m3=dict(range=list(CONVEYANCE_KWH_M3), basis='pipes and channels to the gardens, beyond the lift'))


def digest(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load():
    data = {p: json.loads((ROOT / p).read_text()) for p in INPUTS}
    for p, schema in INPUTS.items():
        if data[p].get('schema') != schema:
            raise SystemExit(f'{p} is {data[p].get("schema")}, expected {schema}')
    return data


def point_w_m2(data):
    """A point of dimming: 1% of the hot world's absorbed sunlight, as the targets count it."""
    return data[TARGETS]['runs']['test']['top_of_air_w_m2']['absorbed'] / 100.0


def budget(waters, others, sink, tg_per_ppm):
    """The air's methane with the waters at a flux: the sources, the steady air the soils alone leave, and at each pair
    the removal, the air through the contactors a year and the years they take to pass the air's methane."""
    sources = [o + waters for o in others]
    row = dict(waters_tg_yr=waters, sources_tg_yr=sources, steady_ppm_by_soils=[s / sink for s in sources])
    for name, ppm in PAIRS_PPM.items():
        removal = [max(s - sink * ppm, 0.0) for s in sources]
        held = ppm * tg_per_ppm
        row[name] = dict(removal_tg_yr=removal, air_passed_per_year_share=[r / held for r in removal],
                         contactors_turnover_yr=[held / r if r > 0 else None for r in removal],
                         turnover_with_the_soils_yr=[held / (r + sink * ppm) for r in removal])
    return row


def methane(data):
    """Each pair's methane removal with the waters traced at 301 K: the ledger's other sources cut by biology's 90%,
    the waters at their traced flux, the soils' sink and the air's methane per ppm, as lake_methane reckons the held
    pair; at the methane share's centre and 0.8 under Osborn's rule, and under Shih and Bouffard's mixing."""
    air = data[METHANE]['air']
    gases = data[LEDGERS]['gases']['methane']
    tg_per_ppm = gases['tg_per_ppm']
    sink = air['lakes_and_seas']['soil_sink_tg_yr_per_ppm']
    ledger_lakes = air['lakes_and_seas']['ledger_lakes_tg_yr']
    others = [0.1 * (c['sources_tg_yr']['total'] - ledger_lakes) for c in gases['central']]
    out = dict(soil_sink_tg_yr_per_ppm=sink, tg_per_ppm=tg_per_ppm, other_sources_after_cut_tg_yr=others,
               pairs_ppm=PAIRS_PPM, waters={}, share_0_8={}, shih_bouffard={})
    for trophic, by_community in air['lakes_and_seas_tg_yr'].items():
        out['waters'][trophic] = {c: budget(w, others, sink, tg_per_ppm) for c, w in by_community.items()}
        out['share_0_8'][trophic] = budget(air['lakes_and_seas_share_tg_yr'][trophic]['high'], others, sink,
                                           tg_per_ppm)
        out['shih_bouffard'][trophic] = {c: budget(w, others, sink, tg_per_ppm)
                                         for c, w in air['lakes_and_seas_shih_tg_yr'][trophic].items()}
    costs = data[LEVERS]['gases']['costs']['low_pair']['methane cut 0.9, N2O cut 0.95']
    out['low_pair_n2o_decomposition_tw'] = costs['n2o_decomposition_tw']
    out['seas_mesotrophic_tg_yr'] = air['seas_tg_yr']['mesotrophic']['earth']
    return out


def room(lines, areas, land_g, people):
    """People's dissolved phosphorus the lakes can take under their lines: each kind's room at the land's load, none
    where the land alone reaches its line, in kg a year [larger, smaller] at the land's low and high load."""
    kg = [sum(max(lines[k] - g, 0.0) * areas[k] * 1e6 / 1e3 for k in KINDS_BY_AREA) for g in land_g]
    return dict(room_kg_yr=kg, per_person_g_yr={f'{int(n / 1e9)}_billion': [x / n * 1e3 for x in kg] for n in people})


def phosphorus(data):
    """The lakes' room for people's dissolved phosphorus under the oligotrophic line, kind by kind and by settling,
    with the giants' floors oxic (Osborn's rule) and anoxic (Shih and Bouffard's form), against the plan's losses; and
    the seas' phosphorus as the lakes and the dust leave it."""
    sp, food, kinds = data[PHOSPHORUS], data[FOOD], data[KINDS]['needs']
    area = sp['partition']['lake_area_m2']
    land_g = [x * 1e3 / area for x in sp['routes']['bioavailable_input_kg_yr']['lakes']['base']]
    people = sp['people']['working_range_people']
    lost = sp['people']['ecology_design_target']['lost_kg_per_person']
    areas = {k: kinds[k]['area_km2'] for k in KINDS_BY_AREA}
    glass = food['systems']['fog_desert_glass_vented']['phosphorus_mt_year']['discharged']   # Mt a year a billion
    out = dict(lake_area_m2=area, kinds_area_km2=areas,
               giants_share_of_lake_area=areas['giant'] / sum(areas.values()), land_bioavailable_g_m2_yr=land_g,
               people=people, lost_kg_per_person=lost, glass_discharge_kg_per_person=list(glass),
               earth_practice_lost_kg_per_person=sp['people']['earth_practice_lost_kg_per_person'], by_settling={})
    for s in SETTLING:
        lines = {k: kinds[k]['oligotrophic_load_g_p_m2_yr'][s] for k in KINDS_BY_AREA}
        anoxic = dict(lines, giant=kinds['giant']['oligotrophic_load_g_p_m2_yr']['anoxic_floors'])
        cases = {}
        for name, ls in (('giants_oxic', lines), ('giants_anoxic', anoxic)):
            r = room(ls, areas, land_g, people)
            r['area_weighted_line_g_m2_yr'] = sum(ls[k] * areas[k] for k in KINDS_BY_AREA) / sum(areas.values())
            r['dissolved_share_of_the_loss_allowed'] = {
                f'{int(n / 1e9)}_billion': [r['room_kg_yr'][1] / (n * lost[1]), r['room_kg_yr'][0] / (n * lost[0])]
                for n in people}
            cases[name] = r
        out['by_settling'][s] = dict(lines_g_m2_yr=lines, anoxic_line_g_m2_yr=anoxic['giant'], **cases)
    seas = sp['seas']['at_1000_years_mg_m3']
    out['seas'] = dict(at_1000_years_mg_m3=dict(rivers=seas['dissolved_after_lakes'],
                                                rivers_and_dust=seas['dissolved_after_lakes_and_dust'],
                                                rivers_and_dust_top=seas['dissolved_after_lakes_and_dust_top'],
                                                lakes_pass_dissolved=seas['lakes_pass_dissolved']),
                       earth_ocean_mean_mg_m3=sp['seas']['earth_ocean']['mean_mg_m3'],
                       lakes_keep_share_of_the_land_loss=sp['partition']['kept_by_lakes_share_of_loss'])
    wide = data[SEA_OXYGEN]['seas_wide']
    out['seas']['deepest_tenth_mg_l'] = {
        form: {t: wide[key][t]['mean'] for t in ('oligotrophic', 'mesotrophic', 'eutrophic')}
        for form, key in (('osborn', 'typical'), ('shih_bouffard_at_the_two_great_basins', 'typical_shih_bouffard'))}
    return out


def dimming(data):
    """The resort's dimming at 301 K with the levers' first four taken off, each case paired as the heat levers pair
    it (the glow's answers at the glow's efficacy, food made in orbit at the beams' end)."""
    row = next(r for r in data[TARGETS]['resort_band']['rows'] if r['ground_k'] == TARGET_K)
    rank = data[LEVERS]['ranking']
    held = data[LEVERS]['gases']['costs']['held']['methane cut 0.9, N2O cut 0.95']
    tw_per_share = [f / s for f, s in zip(held['methane_fan_tw'], held['methane_air_share_per_year'])]
    low = data[LEVERS]['gases']['options']['low pair, 0.3 and 0.03 ppm']
    return dict(resort_own_percent=row['dimming']['own_percent'],
                resort_in_all_percent=row['dimming']['with_the_heat_beyond_sunlight_percent'],
                levers_first_four_points=rank['first_four_points'],
                levers=[dict(rank=x['rank'], lever=x['lever'], points=x['points']) for x in rank['levers']],
                resort_in_all_with_the_first_four_percent=rank['first_four_left_percent']['resort, 301 K'],
                first_four_share_of_the_resorts_dimming=rank['first_four_share_of_dimming_in_all']['resort, 301 K'],
                low_pair_points=low['points'], methane_fan_tw_per_share_of_air_a_year=tw_per_share)


def watering(data):
    """What watering the deserts at 301 K would take: by latitude band, per 1,000 km2 of garden, the water, the lift
    from the band's own groundwater, the power, its heat and the sunlight a greener ground keeps at the band's own
    sunlight; and the scale by share of the deserts, the watered land spread over them as their area lies."""
    coast, lv, power = data[COAST], data[LEVERS], data[POWER]
    point = point_w_m2(data)
    demand = coast['demand']['desert_uptake_mm_day']
    per_km2_m3 = demand * 1e-3 * 365.25 * 1e6
    heat_per_w = power['absorption']['link']['ground_heat_per_w']
    # the levers' desert albedo: W/m2 for each unit of albedo on each unit of the Moon's sunlit desert share
    per_unit = []
    for name in ('hot world', 'lock'):
        d = lv['albedo']['deserts'][name]
        share = d['sunlit_desert_share']['every_cell']
        per_unit += [d['w_m2'][0] / (0.10 * share), d['w_m2'][1] / (0.30 * share)]
    per_unit = [min(per_unit), max(per_unit)]

    def cost(km2, lift_m, sun):
        water = km2 * per_km2_m3
        kwh = [1000.0 * MOON_SURFACE_GRAVITY * lift_m / LIFT_EFFICIENCY / 3.6e6 + c for c in CONVEYANCE_KWH_M3]
        gw = [water * k * 3.6e6 / (365.25 * 86400) / 1e9 for k in kwh]
        sunlit = km2 * 1e6 / MOON_AREA_M2 * sun
        albedo = [sunlit * GREENING_ALBEDO_DROP[0] * per_unit[0], sunlit * GREENING_ALBEDO_DROP[1] * per_unit[1]]
        return dict(land_km2=km2, water_km3_yr=water / 1e9, lift_m=lift_m, kwh_per_m3=kwh, power_gw=gw,
                    heat_points=[g * 1e9 * heat_per_w / MOON_AREA_M2 / point for g in gw],
                    albedo_points=[w / point for w in albedo])

    deserts, bands = {}, {}
    for form, mask in (('additive', 'rain_mask'), ('log_linear', 'log_rain_mask')):
        cell = coast['water_table'][mask]['every_cell']
        rows = {}
        for b, v in cell['bands'].items():
            if not v['land_km2']:
                continue
            lo, hi = (float(x) for x in b.split('-'))
            sun = 4.0 / math.pi * math.cos(math.radians(0.5 * (lo + hi)))     # the band's diurnal mean sunlight
            rows[b] = dict(desert_km2=v['land_km2'], median_depth_m=v['median_depth_m'], sunlight_over_mean=sun,
                           per_1000_km2=cost(GARDEN_KM2, v['median_depth_m'], sun))
        total = sum(r['desert_km2'] for r in rows.values())
        deserts[form] = dict(land_km2=cell['land_km2'], median_depth_m=cell['median_depth_m'],
                             area_weighted_sunlight=sum(r['desert_km2'] * r['sunlight_over_mean']
                                                        for r in rows.values()) / total,
                             area_weighted_depth_m=sum(r['desert_km2'] * r['median_depth_m']
                                                       for r in rows.values()) / total,
                             young_median_depth_m=[cell['young_1000_yr']['median_depth_m'],
                                                   cell['young_10000_yr']['median_depth_m']])
        bands[form] = rows
    scale = []
    for share in WATERED_SHARES:
        row = dict(share_of_deserts=share)
        for form, d in deserts.items():
            row[form] = cost(d['land_km2'] * share, d['area_weighted_depth_m'], d['area_weighted_sunlight'])
        scale.append(row)
    return dict(demand_mm_day=demand, water_m3_yr_per_km2=per_km2_m3, deserts=deserts, bands=bands,
                heat_per_w=heat_per_w, point_w_m2=point,
                albedo_w_m2_per_unit_albedo_per_unit_sunlit_share=per_unit,
                rivers_km3_yr=data[PHOSPHORUS]['seas']['discharge_m3_yr'] / 1e9, rows=scale)


def build():
    data = load()
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, design_guesses=DESIGN_GUESSES,
               target_k=TARGET_K, methane=methane(data), phosphorus=phosphorus(data), dimming=dimming(data),
               watering=watering(data))
    out['producer'] = dict(files={p: digest(p) for p in CODE}, inputs={p: digest(p) for p in INPUTS},
                           constants=dict(MOON_RADIUS=MOON_RADIUS, MOON_SURFACE_GRAVITY=MOON_SURFACE_GRAVITY))
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
