"""Accounting counterexamples and probability controls, not implementation mirrors."""
from copy import deepcopy
import numpy as np
import pytest
from research.studies.production_income import model as m
from research.studies.production_income.extract import HERE, admitted, digest, read

@pytest.fixture(scope='module')
def data():
    return read(HERE/'inputs/year200.json.gz')

@pytest.fixture(scope='module')
def spec():
    return read(HERE/'income_inputs.json')

@pytest.fixture
def example(data,spec):
    rows={r['name']:r for r in data['cases']}
    p={k:np.array([v[1]]) for k,v in spec['triangular_priors'].items()}
    return deepcopy(rows['propellant_logistics:breakthrough:wider']), rows[spec['base_resident_case']], p


def test_imported_capital_has_no_direct_gdp_or_income_gain(example,spec):
    r,b,p=example
    before=m.accounts(r,b,p,spec)
    r['imported_growth_capital']+=1e12
    after=m.accounts(r,b,p,spec)
    np.testing.assert_allclose(after['gross_capital_formation']-before['gross_capital_formation'],1e12,rtol=1e-12,atol=.01)
    np.testing.assert_allclose(after['total_imports']-before['total_imports'],1e12,rtol=1e-12,atol=.01)
    assert after['gdp']==before['gdp']
    assert m.income(after,.75)['gni']==m.income(before,.75)['gni']


def test_in_kind_services_add_compensation_without_new_profit(example,spec):
    r,b,p=example
    a=m.accounts(r,b,p,spec)
    p['cost_valuation_return']=np.array([.04])
    z=m.accounts(r,b,p,spec)
    np.testing.assert_allclose(z['gdp']-a['gdp'],z['compensation']-a['compensation'],atol=.01)
    np.testing.assert_allclose(z['net_profit'],a['net_profit'],atol=.01)
    assert z['cash_compensation']==a['cash_compensation']


def test_gross_income_keeps_depreciation(example,spec):
    r,b,p=example
    a=m.accounts(r,b,p,spec)
    r['depreciation']+=1e12
    z=m.accounts(r,b,p,spec)
    assert z['gdp']==a['gdp']
    np.testing.assert_allclose(m.income(z,0)['gni']-m.income(a,0)['gni'],1e12,atol=.01)
    np.testing.assert_allclose(m.income(z,1)['gni'],m.income(a,1)['gni'],atol=.01)


def test_distributions_and_core_restoration_are_not_earned_income(example,spec):
    r,b,p=example
    a=m.accounts(r,b,p,spec)
    r['annual_owner_payout']+=100e12
    r['retained_core_service_restoration']=100e12
    z=m.accounts(r,b,p,spec)
    for k in a:
        np.testing.assert_array_equal(a[k],z[k])


def test_signed_foreign_losses_are_not_floored():
    a=dict(gdp=np.array([3e12]),net_profit=np.array([-2e12]),external_interest=np.array([.1e12]))
    assert m.income(a,0)['gni']==4.9e12
    assert m.income(a,1)['gni']==2.9e12
    np.testing.assert_array_equal(m.income(a,.5)['gnp'],m.income(a,.5)['gni'])
    with pytest.raises(ValueError,match='Portfolio-equity'):
        m.income(a,.5,foreign_direct_investment_share=.5)


def test_outside_primary_income_is_separate_from_domestic_output(example,spec):
    r,b,p=example
    a=m.accounts(r,b,p,spec)
    z=m.income(a,.5,outside_income=1e12)
    np.testing.assert_allclose(z['gni']-m.income(a,.5)['gni'],1e12,atol=.01)
    np.testing.assert_allclose(m.income(a,1)['gni'],a['gdp']-a['external_interest'],atol=.01)


def test_prior_support_and_weighted_quantile_gap():
    np.testing.assert_array_equal(m.weighted_quantile([1,2,100,101],[.25]*4,[.1,.5,.5001,.9]),[1,2,100,101])
    for mode in (0,.4,1):
        x=m.triangular_ppf(np.linspace(0,1,101),0,mode,1)
        assert x[0]==0 and x[-1]==1 and np.all(np.diff(x)>=0)


def test_probability_inputs_and_source_identity(data,spec):
    m.verify_inputs(spec,data)
    s=deepcopy(spec)
    s['products']['orders']['sha256']='0'*64
    with pytest.raises(ValueError,match='Source hash'):
        admitted(s,'orders')
    bad=deepcopy(data)
    bad['cases'].pop()
    with pytest.raises(ValueError,match='combinations'):
        m.verify_inputs(spec,bad)
    bad=deepcopy(data)
    bad['cases'][0]['revenue_by_sector']['sports']+=1e12
    with pytest.raises(ValueError,match='receipts'):
        m.verify_inputs(spec,bad)


def test_case_coverage_asset_boundary_and_actual_ownership(data):
    assert len(data['cases'])==81
    assert sum(r['asset_method']=='exported' for r in data['cases'])==18
    assert sum(bool(r['local_shares']) for r in data['cases'])==18
    assert all(r['service_capital']==0 for r in data['cases'])
    assert max(abs(e) for e in data['calibration']['leave_one_out_relative_error_range'])<.14
    assert all(0<=s<=1 for r in data['cases'] for s in r['local_shares'].values())
    assert all(r['external_interest']>=0 and r['local_interest']>=0 for r in data['cases'])


def test_all_accounts_close_on_different_draws(data,spec):
    p=m.draws(spec,257,19)
    base=next(r for r in data['cases'] if r['name']==spec['base_resident_case'])
    for r in data['cases']:
        a=m.accounts(r,base,p,spec)
        assert np.max(np.abs(a['gdp']-a['expenditure_gdp']))<.1
        assert np.max(np.abs(a['gdp']-a['income_gdp']))<.1
        for share in (0,.5,.75,1):
            expected=a['compensation']+a['depreciation']+a['local_interest']+share*a['net_profit']
            np.testing.assert_allclose(m.income(a,share)['gni'],expected,rtol=1e-12,atol=.1)
        assert np.all(a['intermediate_imports']>=0)
        assert np.all(a['imported_capital']>=0)
        assert np.all(a['exports']<=r['exports'])


def test_saved_result_and_sampling_properties(data,spec):
    x=read(HERE/'results/year200.json.gz')
    assert x['producer']['sha256']==digest(HERE/'model.py')
    for path,value in x['dependencies'].items():
        assert digest(HERE/path)==value
    assert max(x['controls'].values())<.1
    for prior in x['distributions'].values():
        assert len(prior)==9
        for portfolio in prior.values():
            for case in portfolio.values():
                assert abs(sum(r['probability'] for r in case['cases'])-1)<1e-12
                for d in [case['gdp'],*case['gni'].values()]:
                    curve=d['inverse_cdf']
                    assert np.all(np.diff(curve['value'])>=0)
                    assert curve['probability'][0]==0 and curve['probability'][-1]==1
                    assert all(0<=v<=1.0000000001 for v in d['probability_below'].values())
                assert 0<=case['probability_year200_operating_deficit']<=1.0000000001
    for case in x['actual_ownership_conditional_breakthrough'].values():
        assert len(case['cases'])==2


def test_portfolio_choices_are_never_pooled(data,spec):
    x=read(HERE/'results/year200.json.gz')
    for prior in x['distributions'].values():
        for portfolio,result in prior.items():
            for dist in result.values():
                assert all(r['name'].split(':')[0]==portfolio for r in dist['cases'])
