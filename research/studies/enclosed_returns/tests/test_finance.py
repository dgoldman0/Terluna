"""Accounting controls independent of the production solver's equations."""
from copy import deepcopy
import json
import math
from research.studies.enclosed_returns import finance as f


def config():
    q=json.loads((f.HERE/'finance_inputs.json').read_text())
    q.update(distribution_policy='retain_to_completion',initial_runway_years=1,liquidity_buffer_years=1,retain_through_year=200)
    return q


def path(capital=None,earnings=None,end=40):
    capital=capital or {8:20e6};earnings=earnings or {}
    rows=[]
    for y in range(8,end+1):
        cost=1e6;available=earnings.get(y,10e6);cap=capital.get(y,0.)
        rows.append(dict(year=y,population=100,service_kits=1,capital=cap,external_capital=cap/2,local_capital=cap/2,gross=cost+available,
            operating_cost=cost,operating_available=available,external_operating=cost,resource_provision=0.,resident_payroll=0.,resident_value=1e9,domestic_transaction_basket=1e9,
            programme_NPV_screen=1e9,revenue_by_sector={'services':cost+available},maximum_capacity_use=1.))
    return dict(name='synthetic',portfolio='services',market='test',growth='static',seed=1000e6,seed_hardware_and_carry=999e6,
        initial_working=1e6,seed_resource_provision=0.,seed_resident_payroll=0.,seed_population=100,annual=rows)


def test_cash_identity_and_credit_is_not_operating_income():
    p=path();q=config();sim=f.finance(p,q,'funded_credit')
    assert sim['annual'][0]['draw']>0
    for r,a in zip(sim['annual'],p['annual']):
        assert abs(r['cash_residual'])<1e-5
        assert r['operating_available']==a['gross']-a['operating_cost']
    debt=0.
    for r in sim['annual']:
        debt+=r['draw']-r['principal'];assert math.isclose(debt,r['debt'],abs_tol=1e-6)
    assert sim['annual'][-1]['debt']==0


def test_retained_earnings_pay_later_orders_without_dividend_topup_churn():
    p=path(capital={10:20e6,15:50e6});q=config();sim=f.finance(p,q,'equity')
    assert all(r['dividend']==0 for r in sim['annual'])
    assert all(r['contribution']==0 for r in sim['annual'])
    assert sim['annual'][-1]['cash']==p['initial_working']+sum(r['operating_available']-r['capital'] for r in p['annual'])


def test_gap_is_a_requirement_and_disqualifies_equity_self_funding():
    p=path(capital={8:60e6});sim=f.finance(p,config(),'equity')
    assert sim['first_unfunded_year']==8
    assert not sim['self_financing_prescribed_schedule']
    assert sim['annual'][0]['contribution']==50e6
    assert math.isclose(sum(sim['annual'][0]['topup_by_purpose'].values()),50e6)


def test_negative_operation_cannot_borrow_for_loss():
    p=path(earnings={8:-10e6});sim=f.finance(p,config(),'U_credit')
    assert sim['annual'][0]['draw']==0
    assert sim['annual'][0]['topup_by_purpose']['operating']==10e6


def test_original_seed_debt_ceiling_survives_large_equity_topups():
    p=path(capital={y:1000e6 for y in range(8,41)},earnings={y:100e6 for y in range(8,41)})
    sim=f.finance(p,config(),'funded_credit')
    assert sum(r['contribution'] for r in sim['annual'])>p['seed']
    assert max(r['debt'] for r in sim['annual'])<=config()['underwriting']['max_debt_fraction_of_seed']*p['seed']+1e-6


def test_refused_local_credit_is_not_replaced_by_external_loan():
    p=path();q=config()
    for r in p['annual']:r.update(resident_value=1.,domestic_transaction_basket=1.,external_capital=0.,local_capital=r['capital'],gross=1e6,operating_cost=1.,operating_available=1e6-1.,external_operating=1.)
    q['monetary']['annual_relative_price_limit']=0.;q['monetary']['cumulative_relative_price_limit']=0.
    control=f.finance(p,q,'equity');sim=f.finance(p,q,'U_credit',control)
    r=sim['annual'][0]
    assert r['external_draw']==0
    assert r['contribution']>0
    assert r['local_U_draw']<.75*p['annual'][0]['capital']


def test_retained_cash_is_terminal_value_not_distributed_return():
    p=path(end=20);q=config();sim=f.finance(p,q,'equity');m=f.metrics(p,sim,q,20)
    assert m['owner_payouts']==0
    assert m['terminal_cash']>0
    v=m['discounted']['0.0']
    assert math.isclose(v['sponsor_cash_and_terminal_NPV']-v['sponsor_distributed_cash_NPV'],v['terminal_financial_claim_PV'])


def test_socioeconomic_credit_is_neither_output_nor_new_resource_cost():
    p=path(end=40);q=config();control=f.finance(p,q,'equity');credit=f.finance(p,q,'U_credit',control)
    a=f.metrics(p,control,q,40);b=f.metrics(p,credit,q,40)
    for rate in a['discounted']:
        assert a['discounted'][rate]['full_price_socioeconomic_NPV']==b['discounted'][rate]['full_price_socioeconomic_NPV']


def test_seed_working_cash_not_counted_as_consumed_resource():
    p=path(end=20);q=config();sim=f.finance(p,q,'equity');a=f.metrics(p,sim,q,20)
    p2=deepcopy(p);p2['initial_working']+=100e6;p2['seed']+=100e6
    sim2=f.finance(p2,q,'equity');b=f.metrics(p2,sim2,q,20)
    assert a['discounted']['0.0']['discounted_resource_cost']==b['discounted']['0.0']['discounted_resource_cost']
    assert math.isclose(a['discounted']['0.0']['sponsor_cash_and_terminal_NPV'],b['discounted']['0.0']['sponsor_cash_and_terminal_NPV'],abs_tol=1e-5)


def test_completion_reserve_allows_early_payout_without_future_equity():
    p=path(capital={15:70e6},end=20);q=config();q['distribution_policy']='completion_funded'
    sim=f.finance(p,q,'equity')
    assert sim['annual'][0]['dividend']>0
    assert all(r['contribution']==0 for r in sim['annual'])
    for i,r in enumerate(sim['annual']):
        if r['dividend']>1:assert not any(x['contribution']>1 for x in sim['annual'][i+1:])


def test_social_benefit_weight_does_not_create_monetary_headroom():
    p=path();q=config();a=f.finance(p,q,'U_credit')
    other=deepcopy(p)
    for r in other['annual']:r['resident_value']*=100
    b=f.finance(other,q,'U_credit')
    assert [r['price'] for r in a['annual']]==[r['price'] for r in b['annual']]


def test_announced_future_shock_is_retained_before_distribution():
    p=path(capital={15:40e6},end=20);q=config();q['distribution_policy']='completion_funded'
    shock=dict(start=17,end=20,factor=0.)
    sim=f.finance(p,q,'equity',receipt_shock=shock)
    for i,r in enumerate(sim['annual']):
        if r['dividend']>1:assert not any(x['contribution']>1 for x in sim['annual'][i+1:])
    m=f.metrics(p,sim,q,20)
    assert sum(m['revenue_by_sector'].values())==m['annual_gross']==0.


def tiny_products():
    request=[];annual=[]
    for y,stock,known,actual_replacement in ((8,100.,20.,6.),(9,200.,40.,1.),(10,200.,40.,1.)):
        r=dict(year=y,population=1.,service_kit_count=1,visitor_beds=0.,service_installed_capital_USD2041=stock,service_new_capital_USD2041=stock if y==8 else (100. if y==9 else 0.),service_replacement_allowance_USD2041=10.,commercial_operating_USD2041=0.,commercial_payroll_USD2041=0.,gross_USD2041=100.,revenue_by_sector_USD2041={'services':100.})
        c=dict(identified_service_stock_value_USD2041=known,identified_service_mean_renewal_allowance_USD2041=2.,capital_external_USD2041=0.,capital_local_USD2041=0.,operating_external_replacement_USD2041=actual_replacement,operating_external_logistics_materials_USD2041=0.,operating_external_services_USD2041=0.,operating_local_payroll_USD2041=0.,operating_total_USD2041=actual_replacement,capital_resource_provision_USD2041=0.,resource_provision_USD2041=0.,resident_payroll_USD2041=0.)
        request.append(r);annual.append(dict(year=y,cash=c,maximum_capacity_use=1.))
    case=dict(name='fixture',portfolio='services',market='test',market_growth_name='static',annual=request)
    physical=dict(convergence={'converged':True},first_capacity_shortfall=None,seed={'cost':{'total_USD2041':100.},'total_resource_provision_USD2041':0.,'resident_payroll_USD2041':0.},annual=annual)
    cfg=config();cfg['capital_cost_allowance_fraction']=0.;cfg['preopening_commercial_fixed_USD2041_year']=0.;cfg['preopening_commercial_payroll_USD2041_year']=0.
    return case,physical,cfg


def test_preopening_facility_order_enters_seed_once():
    case,physical,cfg=tiny_products();p=f.plan(case,physical,cfg)
    # Initial known equipment is already in the physical100. Residual80 for the
    # initial package plus40 of the next package precede commercial operations.
    assert p['seed_hardware_and_carry']==220.
    assert p['annual'][0]['service_residual_capital']==40.
    assert p['annual'][1]['service_residual_capital']==0.


def test_smooth_facility_allowance_does_not_duplicate_lumpy_equipment_renewal():
    case,physical,cfg=tiny_products();p=f.plan(case,physical,cfg)
    assert [r['service_residual_renewal'] for r in p['annual']]==[8.,8.,8.]
    assert sum(r['operating_cost'] for r in p['annual'])==6.+1.+1.+3*8.


def test_unfulfilled_physical_plan_has_no_completed_return_account():
    case,physical,cfg=tiny_products();physical['first_capacity_shortfall']=9
    import pytest
    with pytest.raises(ValueError,match='unresolved requirements'):f.plan(case,physical,cfg)


def test_capital_allowance_changes_budget_without_replicating_operating_assets():
    case,physical,cfg=tiny_products();direct=f.plan(case,physical,cfg)
    cfg['capital_cost_allowance_fraction']=.35;budget=f.plan(case,physical,cfg)
    assert budget['seed_hardware_and_carry']==1.35*direct['seed_hardware_and_carry']
    for a,b in zip(direct['annual'],budget['annual']):
        assert b['operating_cost']==a['operating_cost']
        assert b['capital']==1.35*a['capital']


def test_unanticipated_receipt_shock_cannot_change_earlier_payouts():
    p=path(capital={15:40e6},end=20);q=config();q['distribution_policy']='completion_funded'
    baseline=f.finance(p,q,'equity')
    shock=f.finance(p,q,'equity',receipt_shock=dict(start=17,end=20,factor=0.,known_from=17))
    assert [r['owner_payout'] for r in shock['annual'] if r['year']<17]==[r['owner_payout'] for r in baseline['annual'] if r['year']<17]
    assert sum(r['contribution'] for r in shock['annual'])>sum(r['contribution'] for r in baseline['annual'])
