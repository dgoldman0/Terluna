import copy
import json
import pytest
from research.studies.world_economy import early_exports as model

@pytest.fixture(scope='module')
def product():return model.results()


def test_committed_screen_reproduces(product):
    assert product==json.loads((model.HERE/'results/early_exports.json').read_text())


def test_source_changes_require_rebinding(tmp_path):
    settings,_=model.load_inputs()
    settings['products'][next(iter(settings['products']))]['sha256']='0'*64
    path=tmp_path/'bad.json';path.write_text(json.dumps(settings))
    with pytest.raises(ValueError,match='Input changed'):model.load_inputs(path)


def test_scenarios_respect_capacity_and_customers(product):
    for rows in product['scenarios'].values():
        for r in rows:
            assert r['realized_visitor_nights']<=r['capacity_visitor_nights']+1e-5
            assert r['factory_sold_kg']<=r['factory_saleable_capacity_kg']+1e-5
            assert r['manufacturing']['sales_U']<=r['factory_customer_orders_budget_U']+1
            for o in r['origin'].values():
                assert o['realized_trips']<=o['trips_before_bed_cap']+1e-5
                assert o['destination_receipts_U']+o['access_fares_excluded_U']<=o['budget_U']+1


def test_cash_and_production_bridges_close(product):
    for rows in product['scenarios'].values():
        for r in rows:
            scale=max(1,r['first_year_total_facility_capital_U'])
            assert abs(r['allocation_closure_U'])<scale*1e-12
            for s in ('tourism','manufacturing'):
                f=r[s]
                assert abs(f['cash_closure_residual_U'])<scale*1e-12
                assert abs(f['cash_bridge_residual_U'])<scale*1e-12
                assert f['initial_borrowing_U']+f['initial_sponsor_or_equity_U']==pytest.approx(f['initial_capital_U'])
                assert f['interest_U']+f['principal_U']==pytest.approx(f['facility_debt_service_U'])


def test_negative_cases_are_explicit_and_do_not_finance_development(product):
    for r in product['scenarios']['weak_demand_high_cost']:
        assert r['assigned_development_cash_U']==0
        assert r['unfunded_cash_gap_U']>0
    high=product['scenarios']['visitor_focused_expansion'][-1]
    assert high['actual_occupancy']<.8
    assert high['origin']['earth']['realized_trips']<high['origin']['earth']['requested_trips']


def test_no_customers_leave_fixed_costs_and_debt():
    settings,_=model.load_inputs();stage=copy.deepcopy(settings['stages'][0]);case=settings['cases']['reference']
    stage['annual_origin_participation']={'earth':0,'array':0}
    settings['factory_customer_budget_shares_of_gwp']={'earth':0,'array':0}
    r=model.evaluate(settings,stage,case)
    assert r['realized_trips']==r['factory_sold_kg']==0
    assert r['unfunded_cash_gap_U']>0
    assert r['tourism']['fixed_operations_U']>0
    assert r['total_closing_facility_debt_U']>0


def test_spare_life_support_can_limit_visitors():
    settings,_=model.load_inputs();settings['spare_life_support_guest_places_per_bed']=.01
    r=model.evaluate(settings,settings['stages'][0],settings['cases']['reference'])
    assert r['actual_occupancy']==pytest.approx(.01)


def test_factory_absence_removes_factory_capital_and_sales(product):
    for r in product['scenarios']['no_manufacturing_business']:
        f=r['manufacturing']
        assert f['sales_U']==f['initial_capital_U']==f['closing_facility_debt_U']==0


def test_old_loans_retire_but_renewal_continues(product):
    for r in product['financing_scale']:
        late=r['after_facility_loan_retirement']
        assert late['total_closing_facility_debt_U']==0
        assert late['tourism']['cfc_and_replacement_U']>0
        assert late['tourism']['facility_debt_service_U']==0
        assert r['outstanding_new_facility_debt_after_year1_U']>0


def test_fares_reduce_affordability_without_becoming_destination_receipts():
    settings,_=model.load_inputs();stage=copy.deepcopy(settings['stages'][-1]);case=settings['cases']['visitor_focused_expansion']
    reference=model.evaluate(settings,stage,case)
    stage['access_fare_U']={w:v*100 for w,v in stage['access_fare_U'].items()}
    costly=model.evaluate(settings,stage,case)
    assert costly['tourism']['sales_U']<reference['tourism']['sales_U']
    assert costly['realized_trips']<reference['realized_trips']
