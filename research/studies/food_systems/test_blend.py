"""The food blend: its rule, its agreement with the food study's mixes, and the heritage meat's arithmetic."""
import json
import math

import pytest

from research.studies.food_systems import blend as bl
from research.studies.food_systems import run as food
from research.studies.food_systems import zones as zo


@pytest.fixture(scope='module')
def product():
    return json.loads(bl.OUT.read_text())


@pytest.fixture(scope='module')
def food_product():
    return json.loads(food.OUT.read_text())


def missing_inputs():
    missing = [p for p in zo.GRID_SCHEMAS if not (bl.ROOT / p).exists()]
    for key in food.BOUND:
        try:
            food.bound_bytes(key)
        except zo.MissingInput as error:
            missing.append(str(error))
    return missing


@pytest.fixture(scope='module')
def basis():
    missing = missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    data = {p: food.load_json(p, s) for p, s in food.INPUTS.items()}
    eco = food.load_bound('ecology_people')
    zones = food.zones_mod.zone_screen(bl.ROOT)
    return food.build_basis(data, eco, food.load_bound('ecology_first_screen'), zones,
                            food.derived_parameters(data, eco), food.load_bound('resources_sea_phosphorus'))


@pytest.fixture(scope='module')
def fresh(basis):
    return json.loads(food.dumps(bl.results()))


def close(a, b, rel=2e-3):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-9)


def test_weights_obey_the_rule(product):
    w = product['weights']
    assert close(sum(w.values()), 1.0, rel=1e-9)
    assert all(bl.FLOOR - 1e-9 <= v <= bl.CAP + 1e-9 for v in w.values())
    assert product['candidates'][0]['weights'] == w


def test_the_grid_minimum_is_the_chosen_blend(basis, product):
    kept = bl.rule_by_key(basis, bl.PLANNED_RULE)
    grid = bl.grid()
    assert len(grid) == 57
    best = min(grid, key=lambda w: (round(bl.blend_case(basis, w, bl.PLANNED_N, 'high', kept)['light_taken_km2'], -3),
                                    round(sum(bl.blend_case(basis, w, bl.PLANNED_N, e, kept)['light_taken_km2']
                                              for e in ('low', 'high')), -3)))
    assert dict(zip(bl.THREE, best)) == product['weights']


def test_power_heavys_part_carries_every_residents_towns(basis):
    kept = bl.rule_by_key(basis, bl.PLANNED_RULE)
    for end in ('low', 'high'):
        ph = bl.part(basis, 'power_heavy', 0.4 * bl.PLANNED_N, bl.PLANNED_N, end, kept)
        towns = ph['land_km2']['highlands'] - ph['fields']['field_km2_by_zone']['highlands']
        assert close(towns, bl.towns_km2(basis, bl.PLANNED_N, end), rel=1e-9)
        other = bl.part(basis, 'power_made', 0.5 * bl.PLANNED_N, bl.PLANNED_N, end, kept)
        assert other['land_km2']['total'] == 0.0


def test_single_mixes_repeat_the_food_study(product, food_product):
    for name, mix in (('power_heavy_alone', 'power_heavy'), ('power_made_alone', 'power_made'),
                      ('fully_lamp_lit_alone', 'fully_lamp_lit')):
        ours = product['reference_planned_case'][name]
        theirs = food_product['mixes'][mix]['cases']['8'][bl.PLANNED_RULE]
        for end in (0, 1):
            assert close(ours['power_tw']['total'][end], theirs['power_tw']['total'][end])
            assert close(ours['heat']['fusion']['heat_tw'][end], theirs['heat']['fusion']['heat_tw'][end])
            assert close(ours['land_km2']['total'][end], theirs['land_km2']['total'][end])
            assert close(ours['n2o']['kt_n_year'][end], theirs['n2o']['kt_n_year'][end])
            assert close(ours['water_km3_year']['rain_used'][end], theirs['water_km3_year']['rain_used'][end])


def test_heritage_meat_arithmetic(product):
    h = product['heritage_meat']
    per_kg = sum(bl.BASKET[k] * bl.LAND_M2_PER_KG[k] for k in bl.BASKET)
    assert close(h['land_m2_per_kg'], per_kg)
    for i, y in enumerate(bl.CARCASS_YIELD):
        ch4 = sum(bl.BASKET[k] * bl.GLEAM_CO2EQ_PER_KG_CW[k] * bl.ENTERIC_SHARE[k] / bl.GWP_CH4 / y
                  for k in bl.GLEAM_CO2EQ_PER_KG_CW)
        assert close(h['cattle_ch4_kg_per_kg'][i], ch4)
    case = h['cases']['8_billion_1_kg']
    assert close(case['meat_mt_year'], 8.0)
    assert close(case['land_km2'], 8e9 * per_kg / 1e6)
    assert close(case['share_of_homestead_region'], case['land_km2'] / h['homestead_region_km2'])
    assert close(case['cattle_ch4_mt_year'][1], 8.0 * h['cattle_ch4_kg_per_kg'][1])


def test_fresh_run_matches_the_product(product, fresh):
    for key in ('weights', 'candidates', 'blend', 'reference_planned_case', 'per_tenth_planned_case', 'heritage_meat'):
        assert fresh[key] == product[key]
