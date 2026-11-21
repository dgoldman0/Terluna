"""The compact model of the vacuum economy: its technology paths, a resident's needs and costs, and the account."""
import json
import math
from pathlib import Path
import pytest
from research.studies.vacuum_economy_compact import compact as c

HERE = Path(__file__).resolve().parents[1]


def test_a_path_holds_its_stated_values_and_interpolates_between_them():
    t = c.technology()
    years = t['years']
    for name, q in t['quantities'].items():
        if 'central' not in q:
            continue
        for case in ('slow', 'central', 'fast'):
            for year, value in zip(years, q[case]):
                assert c.path(name, case, year) == pytest.approx(value)
            assert c.path(name, case, 2300) == q[case][-1] and c.path(name, case, 2000) == q[case][0]
        assert all(c.path(name, 'static', y) == q['slow'][0] for y in (2041, 2100, 2241))
    # Geometric between stated years: halfway between two stated years the value is their geometric mean.
    q = t['quantities']['array_W_kg']
    a, b = q['slow'][2], q['slow'][3]
    assert c.path('array_W_kg', 'slow', (years[2] + years[3]) / 2) == pytest.approx(math.sqrt(a * b))


def test_the_static_case_returns_the_founding_design():
    d = c.design()
    n = c.needs('static', 2100)
    assert n['food_W'] == pytest.approx(d['food_W'], rel=0.01)
    assert n['mean_W'] == pytest.approx(d['food_W'] + d['basic_W'] + d['gardens_W'] + d['industry_W'] + d['cooling_W'], rel=0.01)
    assert n['local_built_kg'] == pytest.approx(d['local_built_kg']) and n['water_kg_year'] == pytest.approx(d['water_makeup_kg_year'])
    assert n['floor_m2'] == pytest.approx(d['floor_m2'])


@pytest.mark.parametrize('case', ['slow', 'central', 'fast'])
def test_improving_technology_lowers_what_a_resident_needs_and_costs(case):
    years = (2041, 2070, 2100, 2150, 2241)
    food = [c.food_W(case, y) for y in years]
    assert all(b < a for a, b in zip(food, food[1:]))
    build = [sum(c.resident_costs(case, y)['build_USD2041'].values()) for y in years]
    run = [sum(c.resident_costs(case, y)['run_USD2041_year'].values()) for y in years]
    assert all(b < a for a, b in zip(build, build[1:])) and all(b < a for a, b in zip(run, run[1:]))
    static = sum(c.resident_costs('static', 2100)['build_USD2041'].values())
    assert build[2] < static


def test_food_respects_the_limits_the_lookup_found():
    q = c.technology()['quantities']
    assert max(q['led_umol_J']['fast']) < 5.27                       # 10% blue, 90% red: the physical limit
    assert max(q['microbial_efficiency']['fast']) <= 0.27            # the theoretical limit for amino acids


@pytest.mark.parametrize('scenario,case', [('central', 'central'), ('high', 'fast'), ('low', 'static')])
def test_the_account_adds_up_and_visits_keep_their_bound(scenario, case):
    a = c.account('constant_growth', scenario, case)
    rows = a['rows']
    total = 0.0
    for r in rows:
        assert r['net'] == pytest.approx(r['receipts'] - r['line_costs'] - r['build'] - r['running'])
        total += r['net']
        assert r['cumulative'] == pytest.approx(total)
        if r['checks']:
            assert r['checks']['visitors'] <= c.VISITOR_CAP_SHARE[case] * 10.2e9
    assert a['funding_needed_USD2041'] == pytest.approx(-min(0.0, min(r['cumulative'] for r in rows)))


def test_the_saved_product_is_current():
    saved = json.loads((HERE / 'results/compact.json').read_text())
    assert saved['schema'] == c.SCHEMA and saved['producer']['sha256'] == c.digest(HERE / 'compact.py')
    assert saved['dependencies'] == c.dependencies()
    assert len(saved['accounts']) == 2 * 3 * len(c.CASES)
