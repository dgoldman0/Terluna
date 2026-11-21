"""Independent controls for revised returns and resource/wealth interfaces."""
from copy import deepcopy
import gzip,json,math
from pathlib import Path
import numpy as np
import pytest
from research.studies.enclosed_returns import finance,revised as r,revised_orders
from research.studies.enclosed_returns.tests.test_finance import path,config


def test_rollout_orderbook_preserves_residential_schedule_and_customer_limits():
    new=revised_orders.product();old=json.loads((r.HERE/'results/orders.json').read_text());lookup={c['name']:c for c in old['cases']}
    assert len(new['cases'])==6
    for case in new['cases']:
        base=lookup[case['experiment']['baseline_name']]
        for a,b in zip(case['annual'],base['annual']):
            assert all(a[k]==b[k] for k in ('year','population','new_residential_people','household_budget_USD2041','institutional_budget_USD2041','mission_budget_USD2041','industrial_budget_USD2041'))
    assert new==json.loads((r.HERE/'results/revised_orders.json').read_text())


def test_sponsor_cash_counts_seed_working_reserve_once():
    p=path(end=40);cfg=config();p.update(seed_direct_cost=p['seed_hardware_and_carry'],seed_capital_cost_allowance=0)
    sim=finance.finance(p,cfg,'equity');cash=r.sponsor_cash(p,sim,40);m=r.measures(p,sim,cfg,40)
    assert sum(cash)==pytest.approx(m['owner_payouts']-m['all_funder_contributions'])
    assert finance.valuation.pv(cash,.03)==pytest.approx(m['discounted']['0.03']['sponsor_distributed_cash_NPV'])
    assert m['initial_full_funding']==p['seed']+sim['initial_reserve']


def test_resource_probability_uses_water_once_and_does_not_reprice_cash():
    physical={'annual':[{'year':y,'physical':{'source_feed_kg_year':10.,'lunar_raw_kg_year':2.,'source_water_kg_year':1.}} for y in range(8,11)],'seed':{'source_feed_kg':100.,'lunar_raw_kg':20.}}
    stock={'body_mass_kg':[1000.,2000.],'water_equivalent_kg':[10.,20.]}
    population={'scenarios':{'central':{'pools':{key:{'total':stock} for key in ('wider','all_no_don_quixote')}}}}
    spec={'horizons':[10],'water_recovery':.5,'availability_shares':[.5,1.]}
    before=deepcopy(physical);result=r.resource_account(physical,None,population,spec)
    assert physical==before
    d=result['checkpoints'][0]['demand'];assert d['feed_kg']==130 and d['water_kg_bounds']==[3.,10.]
    trial=result['checkpoints'][0]['trials'][0]
    assert trial['availability'][1]['joint_stock_sufficiency_frequency_bounds']==[.5,1.]
    assert result['qualified_source_supply'] is None and not result['extraction_authorized']


def test_exact_seed_loop_demand_and_failed_schedule():
    audit={'annual':[{'year':y,'source_feed_kg':10.,'lunar_raw_kg':2.,'source_water_equivalent_kg':1.,'source_nitrogen_kg':.1,'source_virgin_metal_kg':.2} for y in range(11)]}
    d=r.cumulative_demand({},audit,10);assert d['feed_kg']==110 and d['water_kg_bounds']==[11.,11.]
    with pytest.raises(ValueError):r.cumulative_demand({}, {'annual':audit['annual'][1:]},10)
    failed={'convergence':{'converged':True},'first_capacity_shortfall':8}
    x=r.evaluate({'name':'failed'},failed,{},[],'combined')
    assert x['status']=='physical_requirements_unresolved' and 'checkpoints' not in x


def test_flow_export_debt_and_payout_reconciliation():
    p=path(capital={8:20e6},end=40);p.update(seed_direct_cost=p['seed_hardware_and_carry'],seed_capital_cost_allowance=0,seed_cost_scope='synthetic',unpriced_deployment_requirements=[])
    cfg=config();cfg['horizons']=[40];eq=finance.finance(p,cfg,'equity');sim=finance.finance(p,cfg,'U_credit',eq)
    physical={'seed':{'total_resource_provision_USD2041':0}}
    case={'name':'synthetic'};flow=r.export_flow(case,p,sim,cfg,physical)
    assert [x['year'] for x in flow['annual']]==list(range(41))
    for row in flow['annual']:
        assert row['debt']==pytest.approx(row['local_U_debt']+row['supplier_U_debt']+row['external_debt'])
    assert sum(x['outside_equity_contribution'] for x in flow['annual'])==pytest.approx(p['seed']+sim['initial_reserve']+sum(x['contribution'] for x in sim['annual']))
    assert sum(x['owner_payout'] for x in flow['annual'])==sum(x['owner_payout'] for x in sim['annual'])


@pytest.fixture(scope='module')
def saved():return r.read(r.HERE/'results/revised.json.gz')


def test_saved_complete_scope_and_original_controls(saved):
    assert len(saved['primary_cases'])==162
    assert saved['original_baseline_reproduction']=={'cases':81,'modes':3,'horizons':[100,200,220],'exact':True}
    old={(x['name'],x['mode'],x['horizon']):x for x in r.read(r.HERE/'results/returns.json')['checkpoints']}
    for case in saved['primary_cases']:
        if case['status']!='conditional_account':assert 'checkpoints' not in case;continue
        assert len(case['checkpoints'])==9
        if case['variant']=='baseline':
            for row in case['checkpoints']:
                before=old[row['name'],row['mode'],row['horizon']]
                assert row['discounted']==before['discounted']
                assert row['owner_payouts']==before['owner_payouts']
    for path,digest in saved['dependencies'].items():assert r.sha(r.HERE/path)==digest


def test_saved_comparison_has_no_new_revenue_and_no_money_as_social_benefit(saved):
    assert len(saved['matched_differences'])==729
    for d in saved['matched_differences']:assert d['delta']['cumulative_gross_receipts']==0
    for case in saved['primary_cases']:
        if case['status']!='conditional_account':continue
        for h in (100,200,220):
            rows=[v for v in case['checkpoints'] if v['horizon']==h]
            for rate in ('0.0','0.01','0.03','0.07'):
                assert len({v['discounted'][rate]['full_price_socioeconomic_NPV'] for v in rows})==1


def test_saved_stratification_and_wealth_flows_reconcile(saved):
    for case in saved['primary_cases']:
        if 'stratification' not in case:continue
        for h in (100,200,220):
            periods=[p for p in case['stratification']['periods'] if p['last_year']<=h]
            row=next(v for v in case['checkpoints'] if v['horizon']==h and v['mode']=='U_credit')
            assert sum(p['values']['0.03']['net_real_cash'] for p in periods)==pytest.approx(row['discounted']['0.03']['unlevered_cash_NPV'],rel=1e-12,abs=1.)
            assert abs(case['stratification']['financing_bridges'][str(h)]['reconciliation_residual_USD2041'])<.25
    flows=r.read(r.HERE/'results/revised_flows.json.gz');assert len(flows['cases'])==18
    cases={c['name']:c for c in saved['primary_cases'] if c['variant']=='combined'}
    for flow in flows['cases']:
        assert [x['year'] for x in flow['annual']]==list(range(221))
        for h in (100,200,220):
            row=next(x for x in cases[flow['name']]['checkpoints'] if x['mode']=='U_credit' and x['horizon']==h)
            assert sum(x['owner_payout'] for x in flow['annual'] if x['year']<=h)==pytest.approx(row['owner_payouts'])
            assert sum(x['outside_equity_contribution'] for x in flow['annual'] if x['year']<=h)==pytest.approx(row['all_funder_contributions'])


def test_saved_failed_interventions_keep_failure_and_resources_separate(saved):
    for case in saved['physical_sensitivities']:
        assert case['resources']['qualified_source_supply'] is None
        if case['status']!='conditional_account':assert 'checkpoints' not in case
    for condition in saved['resource_conditions']:
        assert condition['qualified_source_supply'] is None and not condition['extraction_authorized']


def test_poor_source_case_uses_its_own_recovery_fraction():
    audit={'annual':[{'year':y,'source_feed_kg':10.,'lunar_raw_kg':2.,'source_water_equivalent_kg':1.,'source_nitrogen_kg':.1,'source_virgin_metal_kg':.2} for y in range(11)],'source_coefficients':{'water_recovery':.5}}
    stock={'body_mass_kg':[1000.,2000.],'water_equivalent_kg':[20.,40.]}
    population={'scenarios':{'central':{'pools':{key:{'total':stock} for key in ('wider','all_no_don_quixote')}}}}
    spec={'horizons':[10],'water_recovery':.7,'availability_shares':[1.]}
    result=r.resource_account({'annual':[{}]},audit,population,spec)
    assert result['water_recovery_used']==.5
    assert result['checkpoints'][0]['trials'][0]['availability'][0]['joint_stock_sufficiency_frequency_bounds']==[.5,.5]


def test_saved_payment_and_financing_stresses_reconcile(saved):
    assert len(saved['financial_risks'])==18
    for case in saved['financial_risks']:
        assert len(case['comparisons'])==8
        assert len(case['payment_episodes']['cases'])==48
        for episode in case['payment_episodes']['cases']:
            result=episode['result']
            assert result['max_accounting_residual']<1.
            if result['status']=='paid_and_repaid':
                assert result['terminal_bridge']<=1.
                assert abs(result['terminal_receivables'])<1.
                assert result['closing_cash']>=-1.
            else:assert result['failure']['reasons']
