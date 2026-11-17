"""The vacuum economy period by period: its technology paths, the settlements from their designs, the engineering chosen
by cost, the account with Earth, residents as a result and the saved product."""
import json
import math
from pathlib import Path

import pytest

from research.studies.vacuum_economy_periods import designs, economy, engineering, run, sales, technology as tech

HERE = Path(__file__).resolve().parents[1]


def test_a_path_holds_its_stated_values_interpolates_and_the_static_case_holds_2041():
    t = tech.roadmap()
    years = t['years']
    for name, q in t['quantities'].items():
        for case in ('slow', 'central', 'fast'):
            for year, value in zip(years, q[case]):
                assert tech.value(name, case, year) == pytest.approx(value)
        assert all(tech.value(name, 'static', y) == pytest.approx(q['slow'][0]) for y in (2041, 2100, 2241))
    q = t['quantities']['array_W_kg']
    mid = (years[1] + years[2]) / 2
    assert tech.value('array_W_kg', 'central', mid) == pytest.approx(math.sqrt(q['central'][1] * q['central'][2]))
    for name in t['arrivals']:
        assert tech.arrival(name, 'static') is None


def test_with_2041_technology_the_founding_system_returns_its_design():
    """Four cities of 50,000 and a hundred towns of 500, at the design's own engineering, take what the design's
    production schedule says a resident takes, less the power hardware the model sizes itself."""
    city = designs.build('city', 'static', 2041, True)
    town = designs.build('town', 'static', 2041, True)
    system = json.loads(designs.CITIES.read_text())['cases']['reference']
    d = system['schedule']['production_demand']
    n = system['population']
    mix = {k: 0.8 * city[k] + 0.2 * town[k] for k in designs.STREAMS}
    assert mix['handling_kg'] == pytest.approx(d['handling'] / n, rel=1e-9)
    assert mix['bulk_kg'] == pytest.approx(d['local_bulk'] / n, rel=1e-9)
    assert mix['metal_kg'] == pytest.approx(d['local_metal'] / n, rel=1e-9)
    assert mix['pressure_kg'] == pytest.approx(d['pressure_components'] / n, rel=1e-9)
    assert mix['equipment_kg'] < d['other_imports'] / n
    need = designs.needs('city', 'static', 2041)
    assert need['end_use_W']['food'] == pytest.approx(system['utilities']['mean_load_components_W']['food'] / n, rel=1e-9)


@pytest.mark.parametrize('case', ['slow', 'central', 'fast'])
def test_improving_technology_lowers_what_a_resident_needs(case):
    years = (2041, 2070, 2100, 2150, 2241)
    food = [designs.food_W(case, y) for y in years]
    assert all(b < a for a, b in zip(food, food[1:]))
    built = [designs.build('town', case, y, False)['metal_kg'] for y in years]
    assert all(b <= a for a, b in zip(built, built[1:]))


def test_the_power_supply_is_the_cheapest_option_offered_and_the_polar_grid_waits_for_its_arrival():
    for case, year in (('central', 2060), ('central', 2120), ('fast', 2100), ('static', 2200)):
        p = engineering.Prices(case, year, 155797.0)
        s = p.firm_supply(site_km=1600)
        assert s['chosen'] == min(s['options'], key=lambda k: s['options'][k]['USD_kW_year'])
        assert ('polar_grid' in s['options']) == tech.arrived('polar_grid', case, year)
        assert ('fusion' in s['options']) == tech.arrived('fusion', case, year)


def test_a_chemical_tug_needs_more_water_the_farther_it_flies_and_stops_closing():
    p = engineering.Prices('central', 2049, 155797.0)
    drawn = [engineering.tug('chemical', leg, p)['drawn_kg'] for leg in (1000, 2000, 3000, 4000)]
    assert all(b > a for a, b in zip(drawn, drawn[1:]))
    assert engineering.tug('chemical', 9000, p) is None


@pytest.mark.parametrize('scenario,case,ramp', [('central', 'central', 'constant_growth'), ('low', 'static', 'front_loaded'),
                                                ('high', 'fast', 'constant_growth')])
def test_the_account_adds_up_and_the_bounds_hold(scenario, case, ramp):
    rows = economy.Run(scenario, case, ramp).simulate()
    total = 0.0
    bound = json.loads((HERE / 'period_inputs.json').read_text())['sales']['visitor_cap_share'][case]
    for r in rows:
        assert r['net'] == pytest.approx(r['receipts'] - r['paid_to_Earth'])
        assert r['paid_to_Earth'] == pytest.approx(sum(r['to_Earth'].values()))
        total += r['net']
        assert r['cumulative'] == pytest.approx(total)
        assert r['visitor_share_of_world'] <= bound
        assert 0.0 <= r['lines_served_share'] <= 1.0
    s = economy.summary(rows)
    assert s['funding_needed_USD2041'] == pytest.approx(-min(0.0, min(r['cumulative'] for r in rows)))
    assert s['funding_needed_USD2041'] >= s['opening_funding_USD2041']


def test_better_technology_needs_less_funding_for_the_plan():
    funding = [economy.summary(economy.Run('central', c, 'constant_growth').simulate())['funding_needed_by_year_USD2041']['200']
               for c in ('static', 'slow', 'central', 'fast')]
    assert funding[0] > funding[1] > funding[3]


def test_residents_as_a_result_take_the_highest_rates_the_test_allows_to_the_plan_and_past_it():
    rates, sched = run.solve('central', 'central', 'steady')
    assert rates['sustained'] and rates['plan_year'] is not None and rates['past_plan'] > 0
    assert sched == run.steady_schedule(rates['to_plan'], rates['past_plan'])
    assert run.sustained('central', 'central', sched, ())[0]
    assert not run.sustained('central', 'central', run.steady_schedule(rates['to_plan'] + 0.002, 0.0), ())[0]
    assert not run.sustained('central', 'central', run.steady_schedule(rates['to_plan'], rates['past_plan'] + 0.002), ())[0]
    rows = economy.Run('central', 'central', schedule=sched).simulate()
    s = economy.summary(rows)
    assert s['NPV_USD2041']['0.03'][str(run.HORIZON)] >= 0 and s['deficit_decades_after_repayment'] == 0
    assert rows[run.HORIZON]['residents'] > run.planned()


def test_the_test_binds_year_by_year_where_the_present_value_alone_would_not():
    """Growth carried on at the programme's pace passes a test on the present value at 3%, since discounting shrinks
    its late deficits, and fails the test that once repaid the account pays its way decade by decade."""
    rows = economy.Run('central', 'central', 'constant_growth').simulate()
    s = economy.summary(rows)
    assert s['NPV_USD2041']['0.03'][str(run.HORIZON)] > 0
    assert s['deficit_decades_after_repayment'] > 0
    assert not run.sustained('central', 'central', [rows[y]['residents'] for y in range(run.HORIZON + 1)], ())[0]


def test_a_ramp_carries_on_past_the_plan_at_its_final_rate_until_the_sites_are_full():
    r = economy.Run('central', 'central', 'constant_growth', last_year=run.HORIZON + 20)
    rows = r.simulate()
    n = [x['residents'] for x in rows]
    assert n[200] == pytest.approx(run.planned())
    assert n[201] / n[200] == pytest.approx(n[200] / n[199], rel=1e-9)
    assert max(n) <= r.capacity() * (1 + 1e-12)
    assert economy.summary(rows)['sites_full_year'] is not None
    assert all(x['housing_short'] == 0 for x in rows)
    held = economy.Run('central', 'central', 'constant_growth', past_plan='hold').simulate()
    assert held[run.HORIZON]['residents'] == pytest.approx(run.planned())
    assert [x['net'] for x in held[:200]] == pytest.approx([x['net'] for x in rows[:200]])


def test_newcomers_past_the_initial_quota_live_on_the_eventual_set_of_sites():
    r = economy.Run('central', 'central', 'constant_growth')
    rows = r.simulate()
    inputs = json.loads((HERE / 'period_inputs.json').read_text())['settlements']
    assert r.capacity() == pytest.approx(4 * 50000 + inputs['eventual_crater_sites']['value'] * 750 + 50000)
    assert rows[200]['housed_beyond_initial_quota'] == 0
    assert rows[run.HORIZON]['housed_beyond_initial_quota'] > 0
    assert rows[run.HORIZON]['housed_town'] <= inputs['eventual_crater_sites']['value'] * 750 + 50000 + 1


def test_water_set_aside_lengthens_the_legs_and_brings_exhaustion_sooner():
    assert engineering.usable_water_kg('p50', 0.9) == pytest.approx(0.1 * engineering.usable_water_kg('p50'))
    for drawn in (1e12, 1e13, 1e14):
        assert engineering.leg_for(drawn, 'p50', 0.9) >= engineering.leg_for(drawn, 'p50')
    whole = run.exhaustion('central', 'central', 'ramp', (('reach', 'p50'), ('water_set_aside', 0.0)))
    kept = run.exhaustion('central', 'central', 'ramp', (('reach', 'p50'), ('water_set_aside', 0.9)))
    assert kept['water_exhausted_year'] < whole['water_exhausted_year']
    for e in (whole, kept):
        assert e['last_year'] == e['water_exhausted_year']
        if e['water_out_because'] == 'all drawn':
            assert e['contained_drawn_kg_at_end'] >= e['usable_water_kg']
        else:
            assert e['water_out_because'] == 'the rest beyond every tug' and e['contained_drawn_kg_at_end'] < e['usable_water_kg']
    assert kept['water_out_because'] == 'all drawn'
    assert whole['residents_when_sites_full'] == pytest.approx(whole['capacity_residents'])


def test_the_saved_product_is_current():
    saved = json.loads((HERE / 'results/periods.json').read_text())
    assert saved['schema'] == run.SCHEMA
    assert saved['producer']['sha256'] == run.digest(HERE / 'run.py')
    assert saved['dependencies'] == run.dependencies()
    assert [c['name'] for c in saved['cases']] == [c[0] for c in run.cases()]
