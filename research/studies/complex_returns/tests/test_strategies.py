from copy import deepcopy
import json
import pytest
from research.studies.complex_returns import strategies as m


@pytest.fixture(scope='module')
def model():return m.inputs()


def test_source_identity_and_current_account(model):
    q,d,b,e,t=model
    assert d['producer']['export_script_sha256']==m.run.digest(m.HERE/'provenance/export_profiles.py')
    assert d['profiles']['current_large']['external']==e['high_value_capacity_case']['external']
    for p in d['profiles'].values():
        assert p['engineering']['feasible']
        assert p['joint_replication']['at_limit']['engineering']['feasible']
        assert p['repricing']['delivered_lunar_goods_kg_year']+p['repricing']['delivered_source_goods_kg_year']==pytest.approx(p['engineering']['saleable_export_kg_year'])


def test_individual_maxima_cannot_be_added_and_inventory_is_not_revenue(model):
    _,d,_,_,_=model
    assert not d['rejected_controls']['all_individual_maxima']['engineering']['feasible']
    assert d['rejected_controls']['all_individual_maxima']['engineering']['capacity_use']['source_heat']>1.7
    assert d['profiles']['source_refinery_inventory']['external']['USD2041']['gross_external_sales']==d['baseline']['external']['USD2041']['gross_external_sales']


@pytest.mark.parametrize('name',['specialty','premium_lunar','propellant_logistics','compute_digital','diversified','services'])
def test_cash_and_earned_capital_conservation(model,name):
    q,d,b,e,t=model;x=m.simulate(q,d,b,e,t,name,'breakthrough',100);earned=0
    for r in x['annual']:
        uses=r['industrial_orders_USD2041']+r['service_orders_USD2041']+r['dividend_USD2041']
        earned+=r['pooled_USD2041']
        if uses:assert uses<=max(0,earned)+.1
        earned-=uses
        assert r['earned_USD2041']==pytest.approx(earned,abs=.1)
        assert r['closing_cash_USD2041']==pytest.approx(r['opening_cash_USD2041']+r['pooled_USD2041']-uses,abs=.1)
        assert r['operation']['resource_screen']['feasible']
        assert r['starts_scale']<=r['scale']*r['growth_fraction_ceiling']+1e-9
    assert x['summary']['external_funding_after_first_production_USD2041']==0
    assert -sum(x['capital_provider_flows_USD2041'][:8])==pytest.approx(x['summary']['external_seed_USD2041'])


def test_services_and_premium_share_industrial_construction_capacity(model):
    q,d,b,_,_=model
    a=m.resource_screen(q,d,b,'premium_lunar',1,0)
    z=m.resource_screen(q,d,b,'premium_lunar',1,20)
    assert z['capacity_use']['ground_firm_power']>a['capacity_use']['ground_firm_power']
    assert z['growth_fraction_ceiling']<a['growth_fraction_ceiling']
    assert m.resource_screen(q,d,b,'premium_lunar',.05,20)['feasible'] is False


def test_commercial_replication_is_paid_and_built(model):
    q,d,b,e,t=model;p=m.profile(q,d,'compute_digital');a=m.capital(q,d,b,e,t,'compute_digital')
    g=p['remaining_balanced_growth'];v=g['at_limit']['external']['USD2041'];base=p['external']['USD2041']
    host=(v['expansion_external_purchases']+v['resource_provision']-base['resource_provision'])/g['fraction_of_residents_year']
    assert a['marginal_post_founding_USD2041_per_scale']>host+4e11
    assert p['joint_replication']['fraction_of_residents_year']<g['fraction_of_residents_year']


def test_joint_customer_budgets_and_compatible_missions(model):
    q,d,b,_,_=model
    for name in q['strategies']:
        a=m.operation(q,d,b,name,'global',5,5,20,30)
        r=a['revenue_by_sector_USD2041'];market=q['markets']['global']
        assert r['materials']+r['source_products']<=a['industrial_budget_USD2041']+1
        assert r['sports']+r['entertainment']+r['licenses']<=market['household_budget_USD2041_year']+1
        assert a['mission_research_USD2041']+r['propellant']+r['mission_support']<=market['mission_budget_USD2041_year']+1
        assert r['research']-a['mission_research_USD2041']+r['professional']+r['computing']<=market['institutional_services_budget_USD2041_year']+1
        assert a['delivered_propellant_kg']<=market['funded_compatible_missions_year']*q['logistics']['delivered_fuel_kg_per_mission']+1
        assert a['external_compute_W']<=market['compute_contracted_capacity_W']+1


def test_research_displaces_mission_sales_and_internal_compute_cannot_be_sold(model):
    q,d,b,_,_=model
    a=m.operation(q,d,b,'propellant_logistics','global',1,1,0,30)
    z=m.operation(q,d,b,'propellant_logistics','global',1,1,5,30)
    assert z['revenue_by_sector_USD2041']['propellant']<a['revenue_by_sector_USD2041']['propellant']
    a=m.operation(q,d,b,'compute_digital','breakthrough',1,1,5,30)
    z=m.operation(q,d,b,'compute_digital','breakthrough',1,1,5,30,synergies=False)
    assert a['external_compute_W']+a['internal_compute_W']==pytest.approx(5e9)
    assert a['revenue_by_sector_USD2041']['computing']<z['revenue_by_sector_USD2041']['computing']
    assert a['avoided_external_digital_cost_USD2041']>0


def test_premium_failure_keeps_costs_and_research_slots_are_reserved(model):
    q,d,b,_,_=model;q=deepcopy(q)
    a=m.operation(q,d,b,'premium_lunar','breakthrough',1,1,1,30)
    q['premium']['accepted_fraction']=0
    z=m.operation(q,d,b,'premium_lunar','breakthrough',1,1,1,30)
    assert a['premium_sold_kg']>0 and z['premium_sold_kg']==0
    assert a['premium_cost_USD2041']==z['premium_cost_USD2041']>0
    assert a['internal_research_slots']>0
    assert a['revenue_by_sector_USD2041']['research']<(1000*b['scenarios']['breakthrough']['research_USD2041_contract'])


def test_linked_competition_reduces_launch_purchases_and_fuel_price(model):
    q,d,b,_,_=model;v=deepcopy(q);v['competition']['Earth_launch_cost_factor']=.25
    a=m.operation(q,d,b,'propellant_logistics','breakthrough',.05,.05,1,30)
    z=m.operation(v,d,b,'propellant_logistics','breakthrough',.05,.05,1,30)
    assert z['core_imports_USD2041']<a['core_imports_USD2041']
    assert z['propellant_price_USD2041_kg']<a['propellant_price_USD2041_kg']


def test_shared_services_can_finance_their_host_without_free_construction(model):
    q,d,b,e,t=model;x=m.simulate(q,d,b,e,t,'services','breakthrough',100)
    assert x['summary']['final_scale']>b['seed_scale']
    assert any(r['joint_programme_reserve_USD2041']>0 for r in x['annual'])
    assert x['summary']['total_industrial_orders_USD2041']>0


def test_failure_and_continuing_operating_reserve(model):
    q,d,b,e,t=model;x=m.simulate(q,d,b,e,t,'commodity_refined','limited',100)
    assert x['summary']['failure']['year']==x['annual'][-1]['year']
    assert x['summary']['terminal_liquid_USD2041']==0
    y=m.simulate(q,d,b,e,t,'propellant_logistics','breakthrough',100)
    assert y['summary']['terminal_liquid_USD2041']==pytest.approx(max(0,y['annual'][-1]['closing_cash_USD2041']-y['annual'][-1]['buffer_USD2041']))


def test_reviewed_product_reproduces():
    assert m.results()==json.loads((m.HERE/'results/strategies.json').read_text())


def test_imposed_cap_is_distinct_from_financed_growth():
    r=json.loads((m.HERE/'results/strategies.json').read_text())
    fixed=next(x for x in r['cap_operation_checks'] if x['strategy']=='diversified' and x['market']=='breakthrough')
    financed=next(x for x in r['horizons'] if x['strategy']=='diversified' and x['market']=='breakthrough' and x['horizon_elapsed_years']==800)
    assert fixed['feasible_under_assumed_replication']
    assert fixed['installed_scale']>financed['final_scale']
    assert fixed['available_USD2041_year']<0<financed['final_available_USD2041']
