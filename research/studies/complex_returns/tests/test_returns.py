from copy import deepcopy
import json
import math
import pytest
from research.studies.complex_returns import run


def test_snapshot_identity_and_saved_reproduction():
    assert run.results()==json.loads((run.HERE/'results/returns.json').read_text())


@pytest.mark.parametrize('r',[0,.01,.03,.07])
def test_constant_cash_matches_independent_annuity(r):
    flow=[-100]+[8]*30
    expected=-100+8*(30 if r==0 else (1-(1+r)**-30)/r)
    assert run.pv(flow,r)==pytest.approx(expected)


def test_foundation_cash_matches_pinned_account_and_excludes_gdp():
    q,e,_=run.inputs()
    a=run.evaluate(q,e,'large',10,.03)
    assert a['operating'][0]['external_cash_USD2041']==pytest.approx(e['high_value_capacity_case']['external']['USD2041']['after_renewal_and_resource'])
    e['high_value_capacity_case']['domestic']['USD2041']['GDP_GVA']*=20
    b=run.evaluate(q,e,'large',10,.03)
    assert a['cash_NPV_USD2041']==b['cash_NPV_USD2041']
    assert a['proxy_economic_NPV_USD2041']==b['proxy_economic_NPV_USD2041']


@pytest.mark.parametrize('r',[0,.005,.01])
def test_growth_pipeline_capacity_and_plateau(r):
    q,e,_=run.inputs();p=run.growth_path(q,e,r,800)
    assert p[0]['population']==p[1]['population']==250000
    assert p[2]['population']==pytest.approx(250000+p[0]['starts'])
    assert all(0<=v['capacity_fraction']<=1+1e-12 for v in p)
    cap=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']
    assert all(x['population']+x['pipeline']<=cap+1e-6 for x in p)
    if r:
        assert p[-1]['population']==pytest.approx(cap)
        assert p[-1]['starts']==0
        assert p[-1]['pipeline']==0
    else:assert p[-1]['population']==250000


def test_growth_reconciles_pinned_marginal_expansion():
    q,e,_=run.inputs()
    a=run.evaluate(q,e,'large',10,0,growth=.01)
    rows=a['operating'];base=e['high_value_capacity_case'];g=e['high_value_capacity_with_growth']
    coef=(g['external']['USD2041']['after_renewal_and_resource']-g['external']['USD2041']['after_expansion_and_resource'])/g['engineering']['growth_people_year']
    assert coef==pytest.approx(22873645.899,rel=1e-6)
    extra=(g['domestic']['USD2041']['GDP_GVA']-base['domestic']['USD2041']['GDP_GVA'])/g['engineering']['growth_people_year']
    assert rows[0]['domestic_GVA_USD2041']==pytest.approx(base['domestic']['USD2041']['GDP_GVA']+extra*rows[0]['starts'])


def test_demand_loss_reduces_trade_and_value_added_without_free_cost_reduction():
    q,e,_=run.inputs()
    a=run.evaluate(q,e,'large',100,.03,growth=.01)
    b=run.evaluate(q,e,'large',100,.03,growth=.01,demand_cap=1)
    x,y=a['operating'][-1],b['operating'][-1]
    lost=x['exports_USD2041']-y['exports_USD2041']
    assert x['imports_USD2041']==y['imports_USD2041']
    assert x['domestic_GVA_USD2041']-y['domestic_GVA_USD2041']==pytest.approx(lost)
    assert x['external_cash_USD2041']-y['external_cash_USD2041']==pytest.approx(lost)


def test_social_valuation_does_not_generate_cash():
    q,e,_=run.inputs()
    a=run.evaluate(q,e,'large',30,.03,service_value=0,export_value=0)
    b=run.evaluate(q,e,'large',30,.03,service_value=1,spillover=1e12)
    assert a['cash_NPV_USD2041']==b['cash_NPV_USD2041']
    assert b['proxy_economic_NPV_USD2041']>a['proxy_economic_NPV_USD2041']


def test_reserve_is_cash_restriction_not_economic_loss():
    q,e,_=run.inputs();a=run.evaluate(q,e,'large',10,0)
    e['high_value_capacity_commissioning']['additional_resource_provision_USD2041']+=1e9
    b=run.evaluate(q,e,'large',10,0)
    assert b['cash_NPV_USD2041']==pytest.approx(a['cash_NPV_USD2041']-1e9)
    assert b['proxy_economic_NPV_USD2041']==pytest.approx(a['proxy_economic_NPV_USD2041'])


def test_financial_and_economic_switching_values_close_their_own_gaps():
    q,e,_=run.inputs();a=run.evaluate(q,e,'large',30,.03)
    factors=run.pv([0]*q['campaign_years']+[1]*30,.03)
    assert a['cash_NPV_USD2041']+factors*a['required_extra_annual_cash_USD2041']==pytest.approx(0,abs=.01)
    assert a['proxy_economic_NPV_USD2041']+factors*a['required_extra_annual_net_benefit_USD2041']==pytest.approx(0,abs=.01)


def test_financing_claims_redistribute_without_creating_value():
    for row in run.results()['illustrative_capital_claims']:
        assert row['investor_NPV_USD2041']+row['other_stakeholder_cash_NPV_USD2041']==pytest.approx(row['combined_NPV_USD2041'])


def test_incremental_upgrade_keeps_whole_campaign_sunk_in_both_options():
    q,e,_=run.inputs();result=run.results()['incremental_commercial_programme'][1]
    annual=e['high_value_capacity_case']['external']['USD2041']['after_renewal_and_resource']-e['baseline_profiles']['local_structures']['external']['USD2041']['after_renewal_and_resource']
    assert result['annual_external_balance_improvement_USD2041']==pytest.approx(annual)
    assert result['incremental_cash_NPV_USD2041']>0
    assert run.evaluate(q,e,'large',30,.03)['cash_NPV_USD2041']<0
