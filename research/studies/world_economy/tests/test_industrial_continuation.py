"""Economic controls independent of a pressure or pumping cost assumption."""
import math

import pytest

from research.studies.world_economy.industrial_continuation import evaluate_continuation
from research.studies.world_economy import industrial_continuation as m


def sources():
    # Invest100 at0; earn30 at1,2; close at2 for40 salvage less5 removal.
    # Or earn30 at3,4 with20 renewal at3, then20 salvage less5 removal at4.
    rate = .1
    common = -100+30/1.1+30/1.1**2
    close = common+35/1.1**2
    long = common+(30-20)/1.1**3+(30+15)/1.1**4
    return dict(years=2, unlevered_project_npv_U=close), dict(years=4, unlevered_project_npv_U=long)


def evaluate(**changes):
    close, long = sources()
    options = dict(real_discount_rate=.1, future_output_kg_year=10)
    options.update(changes)
    return evaluate_continuation(close, long, **options)


def test_increment_includes_foregone_and_later_terminal_credits_once():
    result = evaluate()['results']
    # Decision-year2 value: forgo35 now; earn net10 at3, net45 at4.
    expected = -35+10/1.1+45/1.1**2
    assert result['base_incremental_value_at_decision_U'] == pytest.approx(expected)
    assert result['incremental_npv_vs_close_U'] == pytest.approx(expected/1.1**2)


def test_zero_added_cost_reproduces_original_long_case():
    _, long = sources()
    result = evaluate()['results']
    assert result['continuation_project_npv_after_added_cost_U'] == long['unlevered_project_npv_U']
    assert result['added_cost_pv_at_initial_year_U'] == 0


def test_costs_are_charged_at_decision_and_following_end_years():
    close, long = sources()
    result = evaluate(retrofit_capex_U=2, added_annual_cost_U=1)['results']
    expected = long['unlevered_project_npv_U']-2/1.1**2-1/1.1**3-1/1.1**4
    assert result['continuation_project_npv_after_added_cost_U'] == pytest.approx(expected)
    assert result['incremental_npv_vs_close_U'] == pytest.approx(expected-close['unlevered_project_npv_U'])
    assert result['added_annual_cost_U_per_kg'] == .1


def test_continuation_and_original_investment_recovery_differ():
    result = evaluate()['results']
    assert result['continuation_preferred_at_given_costs']
    assert not result['whole_project_recovers_at_given_costs']
    assert result['continuation_retrofit_limit_U'] > 0
    assert result['whole_project_retrofit_limit_U'] < 0
    result_at_limit = evaluate(retrofit_capex_U=result['continuation_retrofit_limit_U'])['results']
    assert result_at_limit['incremental_npv_vs_close_U'] == pytest.approx(0, abs=1e-12)
    assert result_at_limit['continuation_project_npv_after_added_cost_U'] < 0


def test_annual_cost_limit_closes_incremental_value():
    limit = evaluate(retrofit_capex_U=1)['results']['continuation_annual_cost_limit_U']
    at_limit = evaluate(retrofit_capex_U=1, added_annual_cost_U=limit)['results']
    assert at_limit['incremental_npv_vs_close_U'] == pytest.approx(0, abs=1e-12)


def test_zero_rate_zero_output_and_signed_allowances():
    close = dict(years=1, unlevered_project_npv_U=10)
    long = dict(years=4, unlevered_project_npv_U=5)
    value = evaluate_continuation(close, long, real_discount_rate=0,
        future_output_kg_year=0, retrofit_capex_U=1, added_annual_cost_U=2)
    assert value['settings']['annuity_factor_at_decision'] == 3
    assert value['results']['incremental_npv_vs_close_U'] == -12
    assert value['results']['continuation_annual_cost_limit_U'] == -2
    assert value['results']['continuation_annual_cost_limit_U_per_kg'] is None


@pytest.mark.parametrize('values', [dict(retrofit_capex_U=-1),
    dict(added_annual_cost_U=float('nan')), dict(retrofit_capex_U=float('inf')),
    dict(real_discount_rate=-.01), dict(future_output_kg_year=-1)])
def test_nonfinite_or_invalid_cost_inputs_rejected(values):
    with pytest.raises(ValueError):
        evaluate(**values)


@pytest.mark.parametrize('year', [1, 2, 2.5, True])
def test_invalid_order_or_year_type_rejected(year):
    close, long = sources(); long['years'] = year
    with pytest.raises(ValueError):
        evaluate_continuation(close, long, real_discount_rate=.1, future_output_kg_year=10)


def test_mismatched_initial_plant_rejected():
    close, long = sources(); close['initial_capital_U']=100;long['initial_capital_U']=110
    with pytest.raises(ValueError):
        evaluate_continuation(close, long, real_discount_rate=.1, future_output_kg_year=10)


def test_nonfinite_source_npv_rejected():
    close, long = sources();long['unlevered_project_npv_U']=float('nan')
    with pytest.raises(ValueError):
        evaluate_continuation(close, long, real_discount_rate=.1, future_output_kg_year=10)


def test_product_reproduces():
    import json
    assert m.results() == json.loads((m.HERE/'results/industrial_continuation.json').read_text())


def test_bound_product_change_rejected(tmp_path):
    import json
    s, _ = m.load_inputs()
    s['products'][next(iter(s['products']))]['sha256'] = '0'*64
    p=tmp_path/'changed.json';p.write_text(json.dumps(s))
    with pytest.raises(ValueError,match='Input changed'):
        m.load_inputs(p)


def test_zero_pressure_difference_still_requires_process_gas_pumping():
    r=m.pumped_chamber(.1,.1,150,1,.02)
    assert r['ideal_reversible_work_J'] == r['loaded_pumpdown_s'] == 0
    assert r['maintenance_pumping_speed_threshold_m3_s'] == pytest.approx(.2)


def test_gas_balance_floor_and_transient_solve_ode():
    ambient,target,V,S,Q,C=100,.1,150,1,.02,1e-5
    r=m.pumped_chamber(ambient,target,V,S,Q,C)
    floor=r['chosen_pump_asymptotic_pressure_Pa']
    assert S*floor == pytest.approx(Q+C*(ambient-floor))
    t=r['loaded_pumpdown_s']
    assert floor+(ambient-floor)*math.exp(-(S+C)*t/V) == pytest.approx(target)
    assert t > r['no_gas_load_pumpdown_s']


def test_internal_gas_can_prevent_target_even_without_external_leak():
    r=m.pumped_chamber(1,.1,150,1,.2)
    assert not r['target_reachable_with_chosen_pump']
    assert r['loaded_pumpdown_s'] is None


def test_lower_pressure_reduces_work_much_more_than_ideal_pumpdown_time():
    low=m.pumped_chamber(100,.1,150,1)
    earth=m.pumped_chamber(101325,.1,150,1)
    assert low['ideal_reversible_work_J']/earth['ideal_reversible_work_J'] < .001
    assert .49 < low['no_gas_load_pumpdown_s']/earth['no_gas_load_pumpdown_s'] < .51


@pytest.mark.parametrize('values', [(1,.1,0,1),(1,.1,150,0),(.01,.1,150,1),(1,.1,150,float('nan'))])
def test_invalid_chamber_conditions(values):
    with pytest.raises(ValueError):m.pumped_chamber(*values)


def test_source_histories_share_operations_before_decision_and_constant_later_output():
    from research.studies.world_economy import industrial_markets as old
    s,p=old.load_inputs();bound=old.anchors(p)
    plant=next(x for x in s['plants'] if x['id']==s['reference']['plant'])
    cost=s['cost_cases'][s['reference']['cost_case']]
    short=old.run_case(s,bound,plant,cost,50,keep_annual=True)
    long=old.run_case(s,bound,plant,cost,100,keep_annual=True)
    assert short['annual'] == long['annual'][:50]
    assert all(r['produced_kg']==pytest.approx(63.7e9) for r in long['annual'][50:])
    now,_=m.load_inputs()
    assert now['continuation']['real_discount_rate']==s['capital']['economic_real_rate']


def test_linear_loading_does_not_assign_full_pressure_at_first_air():
    rows=m.results()['loading_example']['samples']
    assert [r['reference_column_pressure_Pa'] for r in rows] == pytest.approx([270.2,2702,13510])
