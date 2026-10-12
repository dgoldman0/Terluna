"""The author's blend of the three food mixes that fit, and heritage meat from the homesteads.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.food_systems.blend   # -> results/food_blend.json

On 12 October 2026 the author chose a blend of the three mixes that fit under the wild share at 4-8 billion,
power-heavy, power-made and fully lamp-lit, and left its weights to an optimum, with small amounts of heritage meat from
the homesteads and hobby farms as a luxury (decisions register, "Automation, the network and food"). This module builds
each blend as one case from the food study's own mix rules and basis:
- power-heavy's share of the diet is grown by its rule, its fields filling the rain-fed bands around the towns of
  every resident, so the fields of a blend fill the wettest land first;
- power-made's and fully lamp-lit's shares add their power, gases, phosphorus, materials and night energy, and no
  land beyond the towns;
- heat, dimming, N2O and methane follow from the totals by the food study's own functions.
A share of the diet is grown for that share of the residents, so a blend of one mix alone repeats that mix's case,
towns aside for the two without fields.

The optimum takes the weights in steps of 5%, each mix carrying at least a tenth of the diet, so all three stay in use
and can grow when another falters, and none more than half, so that losing any one leaves at least half the food in
place. Among those it takes the blend that takes the least sunlight from wild life at 8 billion under the author's
compromise, at the high-cost end and summed over both ends: the fields' area, plus the area-equivalent of the added
dimming that answers the food's heat under fusion, the planning assumption for local power (the dimming takes the same
share of sunlight from every plant on the Moon). Ties go to the blend with less heat.

The heritage meat is a basket of the author's examples, priced per kilogram with Earth's farms' land use and methane:
nothing here runs a livestock or grassland model.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS, MOON_RADIUS
from shared.provenance import constants_used

from research.studies.food_systems import mixes as mixes_mod
from research.studies.food_systems import run as food
from research.studies.food_systems import systems as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.food-blend/1'
OUT = HERE / 'results' / 'food_blend.json'
CODE = ['blend.py'] + food.CODE
FOOD_PRODUCT = 'research/studies/food_systems/results/food_systems.json'

THREE = ('power_heavy', 'power_made', 'fully_lamp_lit')
STEP = 0.05                      # weight grid
FLOOR, CAP = 0.10, 0.50          # each mix's least and greatest share of the diet
PLANNED_N = 8.0                  # billion lunar residents, the author's C'
PLANNED_RULE = 'wild_30_hybrid_20'
SUPPLY = 'fusion'                # the planning assumption for local power (provisioning README)
MOON_KM2 = 4.0 * math.pi * MOON_RADIUS ** 2 / 1e6

# Heritage meat: the author's examples as a basket, kg of each per kg of heritage meat (design guess).
BASKET = dict(poultry=0.40, pork=0.30, beef_dairy_herd=0.25, beef_beef_herd=0.05)
BASKET_NOTE = ('The author\'s examples: homestead dairy and vaca vieja (beef from dairy herds), some wagyu (beef from '
               'beef herds), heritage pork and chicken. The shares are a design guess.')
# Land a kg of product takes a year, m2, global means of Poore & Nemecek (2018), Science 360, 987, as given by Our
# World in Data's "Land use per kilogram of food product" (Poore, 2010 data).
LAND_M2_PER_KG = dict(poultry=12.22, pork=17.36, beef_dairy_herd=43.24, beef_beef_herd=326.21)
# Cattle's enteric methane: FAO GLEAM's global intensities of 67.8 (beef herds) and 18.4 (dairy herds) kg CO2-eq per kg
# carcass weight, of which enteric fermentation is 43% in beef and 46% in dairy supply chains (Gerber et al. 2013,
# p. 23; Opio et al. 2013 give methane as 44% and 50% of the totals, over 97% and 92% of it enteric), at the GWP of 25
# the assessment used. Per kg of fat- and bone-free meat, Poore & Nemecek's basis, the carcass value is divided by the
# carcass's yield of boneless (0.667) or retail (0.705) beef, USDA ERS (1992), AH-697 Table 7, 1988-1991. Pigs' and
# poultry's enteric methane is small: IPCC (2006) Tier 1 gives swine 1.5 kg a head a year against 27-128 for cattle.
GLEAM_CO2EQ_PER_KG_CW = dict(beef_beef_herd=67.8, beef_dairy_herd=18.4)
ENTERIC_SHARE = dict(beef_beef_herd=0.43, beef_dairy_herd=0.46)
GWP_CH4 = 25.0
CARCASS_YIELD = (0.705, 0.667)   # kg of meat a kg of carcass: retail (low-cost end), boneless (high-cost end)
PROTEIN_FRACTION = (0.20, 0.25)  # protein in lean raw meat, an assumption for the share of the diet's animal part
HERITAGE_KG = (0.5, 1.0, 2.0)    # kg of heritage meat a lunar resident a year


def part(basis, mix, n_part, n_all, end, kept):
    """One mix's share of a blend: its rule for n_part billion residents' diet. Power-heavy's case carries the towns
    of all n_all billion residents, so its fields fill around them; the other two carry no towns."""
    b = copy.deepcopy(basis)
    factor = n_all / n_part if mix == 'power_heavy' else 0.0
    b['diet']['settlement_m2'] = [x * factor for x in b['diet']['settlement_m2']]
    return mixes_mod.mix_case(b, mix, n_part, end, kept)


def towns_km2(basis, n, end):
    return n * (1.0 - mixes_mod.ALOFT_SHARE) * S.pick(basis['diet']['settlement_m2'], end) * S.B / S.KM2


def blend_case(basis, weights, n, end, kept):
    """A blend of the three mixes for n billion residents at one end under one land rule."""
    parts = {m: part(basis, m, w * n, n, end, kept) for m, w in zip(THREE, weights) if w > 0.0}
    towns = towns_km2(basis, n, end)
    ph = parts.get('power_heavy')
    land = dict(ph['land_km2']) if ph else dict(rain_belt=0.0, middle_wet=0.0, fog_desert=0.0, highlands=towns,
                                               polar=0.0, total=towns, glass=0.0)
    power = {}
    for p in parts.values():
        for key, value in p['power_tw'].items():
            power[key] = power.get(key, 0.0) + value
    total_power = power['total']
    heat = S.heat(total_power, basis, end)
    fields = land['total'] - towns
    dimming = heat[SUPPLY]['added_dimming_percent']
    kt_n = sum(p['n2o']['kt_n_year'] for p in parts.values())
    ch4 = sum(p['methane']['mt_year'] for p in parts.values())
    phosphorus = {}
    for p in parts.values():
        for key in ('eaten', 'lost', 'in_microbial_biomass', 'glass_discharge'):
            phosphorus[key] = phosphorus.get(key, 0.0) + p['phosphorus_mt_year'][key]
    structure = {}
    for p in parts.values():
        for key, value in p['structure_gt'].items():
            structure[key] = structure.get(key, 0.0) + value
    night = sum(p['night_energy_twh'] for p in parts.values())
    kept_dry = basis['zones']['dry_land_km2'] * (1.0 - kept)
    return dict(
        weights=dict(zip(THREE, weights)),
        people_billion=n,
        land_km2=land, towns_km2=towns, fields_km2=fields,
        dry_land_wild_share_outside_kept=1.0 - land['total'] / kept_dry,
        power_tw=power,
        heat=heat,
        light_taken_km2=fields + dimming / 100.0 * MOON_KM2,
        light_taken_share_of_moon=(fields + dimming / 100.0 * MOON_KM2) / MOON_KM2,
        water_km3_year=ph['water_km3_year'] if ph else dict(rain_used=0.0, withdrawn_from_seas=0.0, rain_on_glass=0.0),
        phosphorus_mt_year=phosphorus,
        n2o=S.n2o(kt_n, basis),
        methane=S.ch4(ch4, basis, end),
        structure_gt=structure,
        night_energy_twh=night,
        night_energy_over_european_night=night / (n * basis['night']['european_final_kwh']),
        middle_wet_fields_km2=land['middle_wet'],
        feasible=ph['feasible'] if ph else True,
        overflowing=ph['overflowing'] if ph else [])


def grid():
    """Every weighting in steps of STEP with each share between FLOOR and CAP."""
    k = round(1.0 / STEP)
    lo, hi = round(FLOOR / STEP), round(CAP / STEP)
    return [(a * STEP, b * STEP, (k - a - b) * STEP) for a in range(lo, hi + 1) for b in range(lo, hi + 1)
            if lo <= k - a - b <= hi]


def rule_by_key(basis, key):
    return next(r for r in basis['land_rules'] if food.sys_mod.kept_key(r) == key)


def optimum(basis):
    """The least light taken at the planned case's high-cost end, then over both ends, then the least heat."""
    kept = rule_by_key(basis, PLANNED_RULE)
    scored = []
    for w in grid():
        lo = blend_case(basis, w, PLANNED_N, 'low', kept)
        hi = blend_case(basis, w, PLANNED_N, 'high', kept)
        scored.append((round(hi['light_taken_km2'], -3), round(lo['light_taken_km2'] + hi['light_taken_km2'], -3),
                       hi['heat'][SUPPLY]['heat_tw'], w, lo, hi))
    scored.sort(key=lambda s: s[:3])
    least_heat = min(scored, key=lambda s: (s[2], s[0]))
    return scored, least_heat


def heritage_meat(basis, kept):
    """Land and cattle methane for heritage meat at each level and population, against the homestead region: the
    middle-wet zone's land outside the kept share, which the blend's fields leave free."""
    land_per_kg = sum(BASKET[k] * LAND_M2_PER_KG[k] for k in BASKET)
    beef_ch4 = {k: [GLEAM_CO2EQ_PER_KG_CW[k] * ENTERIC_SHARE[k] / GWP_CH4 / y for y in CARCASS_YIELD]
                for k in GLEAM_CO2EQ_PER_KG_CW}
    ch4_per_kg = [sum(BASKET[k] * beef_ch4[k][i] for k in beef_ch4) for i in (0, 1)]
    region = basis['zones']['zones']['middle_wet']['dry_land_km2'] * (1.0 - kept)
    animal_g = basis['diet']['animal_source_protein_g_day']
    out = dict(land_m2_per_kg=land_per_kg, cattle_ch4_kg_per_kg=ch4_per_kg, beef_ch4_kg_per_kg_by_herd=beef_ch4,
               homestead_region_km2=region, animal_part_protein_g_day=animal_g, cases={})
    for n in mixes_mod.POPULATIONS:
        for kg in HERITAGE_KG:
            mt = n * kg                                   # Mt of heritage meat a year (billion people x kg)
            land = mt * 1e9 * land_per_kg / S.KM2
            ch4 = [mt * c for c in ch4_per_kg]
            out['cases'][f'{n:g}_billion_{kg:g}_kg'] = dict(
                meat_mt_year=mt, land_km2=land, share_of_homestead_region=land / region,
                share_of_animal_part=[kg * 1e3 * f / JULIAN_YEAR_DAYS / animal_g for f in PROTEIN_FRACTION],
                cattle_ch4_mt_year=ch4,
                steady_ppm=[S.ch4(ch4[0], basis, 'low')['steady_ppm'], S.ch4(ch4[1], basis, 'high')['steady_ppm']])
    return out


def results():
    data = {p: food.load_json(p, s) for p, s in food.INPUTS.items()}
    eco = food.load_bound('ecology_people')
    screen = food.load_bound('ecology_first_screen')
    sea_p = food.load_bound('resources_sea_phosphorus')
    zones = food.zones_mod.zone_screen(ROOT)
    derived = food.derived_parameters(data, eco)
    basis = food.build_basis(data, eco, screen, zones, derived, sea_p)
    scored, least_heat = optimum(basis)
    best = scored[0]
    w = best[3]
    chosen = {}
    for n in mixes_mod.POPULATIONS:
        chosen[f'{n:g}'] = {food.sys_mod.kept_key(k): S.both(lambda end, n=n, k=k: blend_case(basis, w, n, end, k))
                            for k in basis['land_rules']}
    kept = rule_by_key(basis, PLANNED_RULE)
    reference = {name: S.both(lambda end, ww=ww: blend_case(basis, ww, PLANNED_N, end, kept))
                 for name, ww in (('power_heavy_alone', (1.0, 0.0, 0.0)), ('power_made_alone', (0.0, 1.0, 0.0)),
                                  ('fully_lamp_lit_alone', (0.0, 0.0, 1.0)), ('even', (1 / 3, 1 / 3, 1 / 3)),
                                  ('least_heat', least_heat[3]))}
    per_tenth = {}
    for i, m in enumerate(THREE):
        ww = [0.0, 0.0, 0.0]
        ww[i] = 0.1
        per_tenth[m] = S.both(lambda end, ww=tuple(ww): blend_case(basis, ww, PLANNED_N, end, kept))
    files = [HERE / f for f in CODE]
    product = dict(
        schema=SCHEMA,
        producer=dict(study='food_systems', module='blend',
                      files={str(f.relative_to(ROOT)): food.digest(f) for f in files},
                      constants=constants_used(files),
                      inputs={p: food.digest(ROOT / p) for p in food.INPUTS},
                      food_product=dict(path=FOOD_PRODUCT, sha256=food.sha256(ROOT / FOOD_PRODUCT)),
                      ignored_inputs={p: food.digest(ROOT / p) for p in food.zones_mod.GRID_SCHEMAS},
                      bound_inputs={k: dict(path_at_integration=v['path'], branch=v['branch'], commit=v['commit'],
                                            sha256=v['sha256']) for k, v in food.BOUND.items()}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        criterion=dict(step=STEP, floor=FLOOR, cap=CAP, planned_population_billion=PLANNED_N,
                       planned_land_rule=PLANNED_RULE, supply=SUPPLY, moon_km2=MOON_KM2,
                       rule=('least light_taken_km2 at the high-cost end of the planned case, then summed over both '
                             'ends, then least heat; light_taken_km2 is the fields\' area plus the added dimming '
                             'under the supply as a share of the Moon\'s surface')),
        weights=dict(zip(THREE, w)),
        candidates=[dict(weights=dict(zip(THREE, s[3])), light_taken_km2=[s[4]['light_taken_km2'],
                                                                           s[5]['light_taken_km2']],
                         heat_tw=[s[4]['heat'][SUPPLY]['heat_tw'], s[5]['heat'][SUPPLY]['heat_tw']],
                         fields_km2=[s[4]['fields_km2'], s[5]['fields_km2']]) for s in scored[:8]],
        blend=chosen,
        reference_planned_case=reference,
        per_tenth_planned_case=per_tenth,
        heritage_meat=heritage_meat(basis, kept),
        assumptions=dict(basket=BASKET, basket_note=BASKET_NOTE, land_m2_per_kg=LAND_M2_PER_KG,
                         gleam_co2eq_per_kg_carcass=GLEAM_CO2EQ_PER_KG_CW, enteric_share=ENTERIC_SHARE,
                         gwp_ch4=GWP_CH4, carcass_yield=list(CARCASS_YIELD), protein_fraction=list(PROTEIN_FRACTION),
                         heritage_kg_per_person_year=list(HERITAGE_KG),
                         heritage_in_the_diet=('The heritage meat stands within the diet\'s animal part, which the '
                                               'blend makes as microbial protein; share_of_animal_part gives its '
                                               'protein against that part, so the blend\'s power is left as '
                                               'computed.')),
        sources=SOURCES)
    return food.rd(product)


EVIDENCE = ('A screen on the food study\'s basis and mix rules, recombined: each blend is one case built from the '
            'three mixes\' rules, with heat, dimming, N2O and methane from the totals. The weights come from a stated '
            'rule over a 5% grid. The heritage meat uses Earth\'s global per-kilogram land use and FAO GLEAM\'s '
            'methane intensities for a basket of the author\'s examples. No crop, livestock, grassland, climate or '
            'ecosystem model runs here.')
READING_RULE = ('As the food study: areas km2, power and heat TW, phosphorus Mt a year, N2O kt N and ppb, methane Mt a '
                'year and ppm at steady state, pairs [low-cost end, high-cost end]. weights are shares of the diet. '
                'blend[population][land rule] is the chosen blend; reference_planned_case and per_tenth_planned_case '
                'are at 8 billion under the compromise. light_taken_km2 is the fields\' area plus the area-equivalent '
                'of the added dimming under fusion. heritage_meat cases are named population_level and give Mt of meat '
                'a year, its land, its share of the homestead region and of the diet\'s animal protein [20%, 25% '
                'protein], and the cattle\'s methane per kg of fat- and bone-free meat [retail, boneless carcass '
                'yield], with steady ppm at the soil sink\'s fast and slow ends.')
SOURCES = [
    'Poore, J. & Nemecek, T. 2018. Reducing food\'s environmental impacts through producers and consumers. Science '
    '360, 987-992; per-kg land use as given by Our World in Data (ourworldindata.org/grapher/land-use-per-kg-poore).',
    'Gerber, P.J. et al. 2013. Tackling climate change through livestock. FAO, Rome.',
    'Opio, C. et al. 2013. Greenhouse gas emissions from ruminant supply chains: a global life cycle assessment. FAO, '
    'Rome.',
    'IPCC 2006. Guidelines for National Greenhouse Gas Inventories, volume 4, chapter 10 (enteric fermentation, '
    'Tier 1).',
    'USDA ERS 1992. Weights, Measures, and Conversion Factors for Agricultural Commodities and Their Products. '
    'Agricultural Handbook 697, Table 7.',
]


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(food.dumps(results()) + '\n')
    print(OUT.relative_to(ROOT))


if __name__ == '__main__':
    main()
