import json
import math
from copy import deepcopy

import pytest

from research.studies.complex_returns import credit as c


@pytest.fixture(scope='module')
def inputs():
    return (*c.s.inputs(),json.loads((c.HERE/'credit_inputs.json').read_text()))


def test_restricted_proceeds_cannot_support_dividends(inputs):
    q,d,b,e,t,config=inputs
    cap=c.s.capital(q,d,b,e,t,'diversified')
    coeff=c.external_capital_coefficients(b,t,c.s.profile(q,d,'diversified'))
    plan=dict(committed=.05,limit=0,kit_cost=0,budget=0,joint=False,eligible_local=0)
    x=c.allocate(b,cap,coeff,plan,200.,100.,100.,0.,100.,config['monetary'])
    assert x['dividend']==0
    assert x['fx_cash']-x['fx_restricted']==100


def test_principal_interest_and_grace_are_separate():
    loans=[dict(year=10,original=100.,remaining=100.,rate=.03,term=10,grace=2,kind='U')]
    assert c.debt_due(loans,10)['interest']==0
    assert c.debt_due(loans,11)['interest']==3
    assert c.debt_due(loans,12)['principal']==0
    assert c.repay(loans,13)['principal']==10
    assert loans[0]['remaining']==90
    assert c.debt_due(loans,14)['interest']==pytest.approx(2.7)
    for year in range(14,23):c.repay(loans,year)
    assert loans[0]['remaining']==0


def test_concession_values_contract_against_benchmark():
    assert c.concession(100,.03,20,2,.03)==pytest.approx(0,abs=1e-10)
    assert 0<c.concession(100,.015,40,2,.03)<100


def test_external_capital_floor_respects_founding_boundary(inputs):
    q,d,b,e,t,config=inputs
    for name in q['strategies']:
        p=c.s.profile(q,d,name);cap=c.s.capital(q,d,b,e,t,name)
        coeff=c.external_capital_coefficients(b,t,p)
        x=c.capital_split(b,cap,coeff,.9,.2,0.,config['monetary'])
        assert x['identified_external_floor']==pytest.approx(.1*(coeff[0]+coeff[1]))
        assert x['eligible_local']+x['other_externally_funded']==pytest.approx(x['cost'])
        assert x['other_externally_funded']>=x['identified_external_floor']


def test_loan_money_displaces_other_issuance_before_becoming_excess(inputs):
    *_,config=inputs;m=config['monetary']
    x=c.price_screen(100,0,0,100,10,0,0,0,10,m)
    y=c.price_screen(100,0,0,100,0,0,0,0,10,m)
    assert x['price']==y['price']==1
    assert x['other_money_creation']+x['circulating_credit']==pytest.approx(x['total_money_creation'])
    z=c.price_screen(100,0,0,100,1000,0,0,0,10,m)
    assert z['price']>1
    escrow=c.price_screen(100,0,0,100,1000,1000,0,0,10,m)
    assert escrow['price']==1


@pytest.mark.parametrize('market',['global','breakthrough'])
def test_financed_paths_conserve_cash_and_debt_and_respect_joint_build(inputs,market):
    q,d,b,e,t,config=inputs
    for name in ('diversified','propellant_logistics'):
        control=c.simulate(q,d,b,e,t,name,market,100,config,'retained_earnings')
        for mode in config['financing']:
            sim=control if mode=='retained_earnings' else c.simulate(q,d,b,e,t,name,market,100,config,mode,control)
            debt=reserve=custody=0.;U_debt=0.;last_kits=1
            for r in sim['annual']:
                if 'failure' in r:break
                debt+=r['external_draw']+r['U_draw']-r['principal']
                U_debt+=r['U_draw']-r['U_principal']
                reserve+=r['reserve_contribution']
                custody+=config['monetary']['incremental_distribution_spending_fraction']*r['dividend']-r['local_payout']
                assert r['debt']==pytest.approx(debt,abs=.01)
                assert r['U_debt']==pytest.approx(U_debt,abs=.01)
                assert r['reserves']==pytest.approx(sim['initial_reserve']+reserve)
                assert r['custody']==pytest.approx(custody,abs=.001)
                assert abs(r['cash_residual'])<1
                assert r['cash']+1>=r['FX_restricted']
                assert r['starts']<=r['scale']*min(b['desired_growth_rate'],r['physical_growth_ceiling'])+1e-10
                if r['covenant_lock']:assert r['dividend']==0
                if mode=='retained_earnings':assert r['external_draw']==r['U_draw']==0
            assert all(not r['price_limit_breach'] for r in sim['annual'] if 'failure' not in r)


def test_no_local_classification_means_no_U_lending(inputs):
    q,d,b,e,t,config=inputs;cfg=deepcopy(config)
    cfg['monetary']['unclassified_capital_local_fraction']=0
    cfg['monetary']['service_capital_local_fraction']=0
    a=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,cfg,'retained_earnings')
    z=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,cfg,'U_local_credit',a)
    assert sum(r.get('U_draw',0) for r in z['annual'])==0


def test_default_recovery_is_dated_and_shared_consistently(inputs):
    q,d,b,e,t,config=inputs;cfg=deepcopy(config)
    cfg['stress'].update(start_year=20,end_year=100,buyer_budget_factor=0.)
    a=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,cfg,'retained_earnings',stress=True)
    z=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,cfg,'mixed_credit',a,stress=True)
    r=c.summary(z,100,cfg,a)
    assert r['failure'] is not None
    last=z['annual'][-1]
    claims=last['debt']+z['failure']['accrued_interest']
    assert r['default_cash_recovery']+r['lender_loss_under_cash_only_recovery']==pytest.approx(claims)
    assert c.summary(z,800,cfg,a)['external_lender_cash_NPV_at_3_percent']==r['external_lender_cash_NPV_at_3_percent']


def test_supplier_U_claims_match_loan_terms_and_stay_outside_circulation(inputs):
    q,d,b,e,t,config=inputs
    a=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,config,'retained_earnings')
    z=c.simulate(q,d,b,e,t,'diversified','breakthrough',100,config,'international_U_mix',a)
    claims=prepayments=0.
    for r in z['annual']:
        if 'failure' in r:break
        claims+=r['supplier_credit_origination']-r['supplier_principal']
        prepayments+=r['supplier_credit_origination']-r['U_financed_external_capital']
        assert r['outside_U']==pytest.approx(claims,abs=.01)
        assert r['supplier_prepay']==pytest.approx(prepayments,abs=.01)
        assert r['supplier_prepay']<=r['U_escrow']+1
        assert r['outside_U']<=r['U_debt']+1
        assert r['circulating_credit']<=max(0.,r['U_debt']-r['outside_U'])+1
    summary=c.summary(z,100,config,a)
    assert summary['cumulative_supplier_U_finance']>0
    assert summary['terminal_supplier_U_claims']==0
    assert summary['price_limit_breach_years']==[]


def test_saved_comparison_reproduces_and_has_valid_reference_identity():
    saved=json.loads((c.HERE/'results/credit.json').read_text())
    assert c.product()==saved
    assert len(saved['horizons'])==9*3*6*4
    assert len(saved['stress'])==9*3*5
    assert len(saved['severe_stress'])==9*5
    for row in saved['horizons']:
        assert row['max_cash_identity_residual']<1
        assert not row['price_limit_breach_years']
        assert row['max_relative_price_change_percent']<=5+1e-7
        assert row['max_annual_positive_premium_growth_percent']<=.5+1e-7
        assert row['default_cash_recovery']>=0
        assert row['reserve_capital_loss']<=row['terminal_reserves']
        if row['mode']=='retained_earnings':
            assert row['cumulative_external_loan_draws']==0
            assert row['cumulative_U_credit_draws']==0
        if row['mode']!='international_U_mix':assert row['cumulative_supplier_U_finance']==0
