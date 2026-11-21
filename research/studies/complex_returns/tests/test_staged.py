from copy import deepcopy
import json
import pytest
from research.studies.complex_returns import staged


@pytest.fixture
def model():
    return staged.inputs()


def test_matched_campaign_decomposition(model):
    _,e,t=model;c=staged.capital_account(e,t)
    assert (c['establishment_USD2041']+c['initial_stocks_USD2041']+
        c['original_supplier_operation_and_renewal_USD2041']+
        c['original_construction_makeup_USD2041']) == pytest.approx(c['original_campaign_USD2041'])
    assert c['one_time_bill_USD2041']==pytest.approx(17.566349012275703e12)


def test_core_replaces_overlapping_sales_and_identified_costs(model):
    q,e,_=model;a=e['high_value_capacity_case']['external']
    margin=sum(a['sales_by_activity_USD2041'][k]*(1-v) for k,v in q['replaced_baseline_Earth_service_fractions'].items())
    x=staged.core_account(q,e,1,1)
    assert x['available_USD2041']==pytest.approx(a['USD2041']['after_renewal_and_resource']-margin)
    assert staged.core_account(q,e,2,2)['available_USD2041']==pytest.approx(2*x['available_USD2041'])


def test_bootstrap_funding_has_no_free_seed_or_omitted_carry(model):
    q,e,t=model;x=staged.simulate(q,e,t,'global',30);s=x['summary']
    assert -sum(x['capital_provider_flows_USD2041'][:8])==pytest.approx(s['external_seed_USD2041'])
    assert s['external_seed_USD2041']==pytest.approx(s['seed_hardware_USD2041']+s['seed_liquidity_USD2041']+s['bootstrap_carry_USD2041'])
    assert s['bootstrap_carry_USD2041']>0
    assert s['external_funding_after_first_production_USD2041']==0


@pytest.mark.parametrize('scenario',['limited','global','breakthrough'])
def test_cash_conservation_and_only_earned_money_expands(model,scenario):
    q,e,t=model;x=staged.simulate(q,e,t,scenario,800);earned=0
    for row in x['annual']:
        uses=row['industrial_orders_USD2041']+row['service_orders_USD2041']+row['dividend_USD2041']
        assert row['closing_cash_USD2041']==pytest.approx(row['opening_cash_USD2041']+row['pooled_USD2041']-uses,abs=.01)
        earned+=row['pooled_USD2041']
        if uses: assert uses <= max(0,earned)+.1
        earned-=uses
        assert row['retained_earned_cash_USD2041']==pytest.approx(earned,abs=.1)
        assert row['pooled_USD2041']+row['other_stakeholder_cash_USD2041']==pytest.approx(row['available_USD2041'])


def test_commissioning_ramp_capacity_staff_and_cap(model):
    q,e,t=model;x=staged.simulate(q,e,t,'breakthrough',800)
    cohorts=[(8,q['seed_scale'])];cap=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']/250000
    frontier=e['high_value_capacity_with_growth']['engineering']['growth_people_year']/250000
    for row in x['annual']:
        year=row['year']
        installed=sum(n for due,n in cohorts if due<=year)
        productive=sum(n*min(1,(year-due+1)/5) for due,n in cohorts if due<=year)
        assert row['installed_scale']==pytest.approx(installed)
        assert row['productive_scale']==pytest.approx(productive)
        assert row['starts_scale'] <= installed*frontier+1e-9
        assert row['installed_scale']+row['pipeline_scale']<=cap+1e-9
        assert row['services']['power_W']<=q['reserved_service_power_W_per_founding_scale']*installed+1
        assert row['services']['staff']<=q['maximum_new_service_staff_fraction_of_residents']*250000*installed+1e-6
        if row['starts_scale']:cohorts.append((year+q['commissioning_years'],row['starts_scale']))
    assert x['summary']['cap_year'] is not None
    assert x['annual'][-1]['starts_scale']==0


def test_cultural_reach_independent_of_staffed_b2b_capacity(model):
    q,_,_=model;a=staged.service_account(q,'breakthrough',20,0)
    b=staged.service_account(q,'breakthrough',20,20)
    assert a['audience']==b['audience']==2e9
    assert a['revenue_by_activity_USD2041']['culture']==b['revenue_by_activity_USD2041']['culture']
    assert a['research_contracts']==a['professional_clients']==0
    assert b['research_contracts']==20000


def test_acquisition_charges_initial_cohort_and_replacement(model):
    q,_,_=model;a=staged.service_account(q,'global',20,1)
    b=staged.service_account(q,'global',21,1,a['paid_accounts'],a['licenses'])
    assert a['acquisition_USD2041']==pytest.approx(a['paid_accounts']*30+a['licenses']*8)
    assert b['acquisition_USD2041']==pytest.approx(a['paid_accounts']*.25*30+a['licenses']*.2*8)


def test_insolvency_is_reported_without_refinancing(model):
    q,e,t=model;x=staged.simulate(q,e,t,'limited',800)
    assert x['summary']['failure']['year']==x['annual'][-1]['year']
    assert x['summary']['failure']['shortfall_USD2041']==-x['annual'][-1]['closing_cash_USD2041']
    assert x['summary']['terminal_liquid_USD2041']==0
    assert x['summary']['total_industrial_orders_USD2041']==0


def test_terminal_cash_preserves_continuing_buffer(model):
    q,e,t=model;x=staged.simulate(q,e,t,'breakthrough',100)
    a=x['annual'][-1];s=x['summary']
    assert s['terminal_liquid_USD2041']==pytest.approx(max(0,a['closing_cash_USD2041']-a['liquidity_buffer_USD2041']))
    assert s['terminal_liquid_USD2041']+s['retained_operating_liquidity_USD2041']==pytest.approx(a['closing_cash_USD2041'])


def test_market_gate_prevents_unfunded_overbuilding_despite_large_service_income(model):
    q,e,t=model;q=deepcopy(q);q['material_compute_demand_cap_multiple']=q['seed_scale']
    x=staged.simulate(q,e,t,'breakthrough',100)
    assert x['summary']['final_scale']==q['seed_scale']
    assert x['summary']['total_industrial_orders_USD2041']==0
    assert x['summary']['total_dividends_USD2041']>1e12


def test_doubling_growth_bill_does_not_reprice_original_seed(model):
    q,e,t=model;a=staged.simulate(q,e,t,'breakthrough',100)
    q=deepcopy(q);q['marginal_capital_multiplier']=2
    b=staged.simulate(q,e,t,'breakthrough',100)
    assert a['summary']['external_seed_USD2041']==b['summary']['external_seed_USD2041']
    assert b['summary']['total_industrial_orders_USD2041']>=a['summary']['total_industrial_orders_USD2041']
    assert b['summary']['NPV_USD2041']['0.03']<a['summary']['NPV_USD2041']['0.03']


def test_reviewed_product_reproduces_with_input_and_producer_hashes():
    assert staged.results()==json.loads((staged.HERE/'results/staged.json').read_text())


def test_earnings_accumulate_for_an_indivisible_service_facility(model):
    q,e,t=model;q=deepcopy(q);q['seed_scale']=.1;q['service_kit']['capital_USD2041']=100e9
    rows=staged.simulate(q,e,t,'global',100)['annual']
    saved=[a for a in rows if a['planned_service_saving_USD2041']>0]
    assert len(saved)>2
    assert any(a['service_orders_USD2041']==100e9 for a in rows)
    assert saved[1]['planned_service_saving_USD2041']>saved[0]['planned_service_saving_USD2041']


def test_fixed_commitment_comparison_uses_same_elapsed_horizon():
    rows=staged.results()['fixed_commitment_comparisons']
    assert len(rows)==4
    for row in rows:
        assert row['fixed_supported_residents_after_year33']==250000
        assert row['staged_supported_residents']<250000
        assert row['fixed_campaign_NPV_at_3pct_USD2041']<row['staged_NPV_at_3pct_USD2041']
