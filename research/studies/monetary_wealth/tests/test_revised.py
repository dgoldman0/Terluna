"""Revised physical/financial linkage and conserved monetary ownership accounts."""
from copy import deepcopy
import gzip
import json

import pytest

from research.studies.monetary_wealth import model, revised


@pytest.fixture(scope='module')
def settings():
    result = json.loads((revised.HERE/'wealth_inputs.json').read_text())
    result['growth_cases'] = ['reference']
    return result


def fixture_products():
    name = 'synthetic:breakthrough:reference'
    rows = []
    assets = []
    for year in range(221):
        row = {k: 0. for k in revised.STOCK_FIELDS+revised.FLOW_FIELDS}
        row.update(year=year, population=100., cash=20., reserves=10.,
                   gross=100. if year>=8 else 0.,
                   resident_payroll=10. if year>=8 else 0.,
                   owner_payout=5. if year>=8 else 0.,
                   dividend=5. if year>=8 else 0.,
                   outside_equity_contribution=1030. if year==0 else 0.)
        rows.append(row)
        assets.append({'year':year,
                       'totals':{'net_remaining_service_value_USD2041':1000.,
                                 'construction_WIP_value_USD2041':0.},
                       'by_functional_sector':{'industry':{'net_remaining_service_value_USD2041':1000.}}})
    finance_case = {'name':name, 'baseline_case':name, 'variant':'combined',
                    'market':'breakthrough', 'growth':'reference',
                    'portfolio':'synthetic', 'seed_cost_scope':'Synthetic accounting fixture',
                    'source_case_sha256':'c'*64,'source_loop_case_sha256':'d'*64,
                    'annual':rows}
    asset_case = {'name':name+':combined', 'baseline_case':name,
                  'variant':'combined', 'source_case_sha256':'c'*64,
                  'source_loop_case_sha256':'d'*64,'annual':assets}
    shared = 'a'*64
    return ({'schema':'terluna.humanities.wealth-finance-flows/1',
             'source_physical_sha256':shared,'cases':[finance_case]},
            {'schema':'terluna.provisioning.enclosed-asset-stock/1',
             'source_physical_sha256':shared,'cases':[asset_case],
             'unvalued_assets':['Unpriced intellectual property']})


def small_config():
    return {'variant':'combined','expected_financial_cases':1,
            'expected_portfolios':1,'rules':['Synthetic accounting fixture']}


def test_admission_pins_exact_bytes_schema_and_has_no_unadmitted_fallback(tmp_path):
    path=tmp_path/'x.json.gz'
    path.write_bytes(gzip.compress(b'{"schema":"expected","cases":[]}\n',mtime=0))
    spec={'path':path.name,'schema':'expected','sha256':revised.sha(path)}
    assert revised.read_product(spec,tmp_path)['cases']==[]
    with pytest.raises(ValueError,match='Unadmitted'):
        revised.read_product({**spec,'sha256':None},tmp_path)
    with pytest.raises(ValueError,match='schema'):
        revised.read_product({**spec,'schema':'other'},tmp_path)
    path.write_bytes(path.read_bytes()+b' ')
    with pytest.raises(ValueError,match='Changed'):
        revised.read_product(spec,tmp_path)


def test_linkage_rejects_mixed_physical_sources_duplicate_cases_and_missing_dates(settings):
    finance,assets=fixture_products()
    assert len(revised.matched_cases(finance,assets,small_config(),settings))==1
    with pytest.raises(ValueError,match='physical sources'):
        revised.matched_cases(finance,{**assets,'source_physical_sha256':'b'*64},small_config(),settings)
    with pytest.raises(ValueError,match='Duplicate'):
        revised.matched_cases(finance,{**assets,'cases':assets['cases']*2},small_config(),settings)
    broken=deepcopy(assets);broken['cases'][0]['annual'].pop(20)
    with pytest.raises(ValueError,match='asset date'):
        revised.matched_cases(finance,broken,small_config(),settings)
    broken=deepcopy(finance);broken['cases'][0]['variant']='baseline'
    with pytest.raises(ValueError,match='variant'):
        revised.matched_cases(broken,assets,small_config(),settings)


def test_sector_debt_and_custody_stock_flows_are_verified(settings):
    finance,_=fixture_products();case=finance['cases'][0]
    for y in range(20,30):
        case['annual'][y].update(debt=20.,local_U_debt=20.)
    case['annual'][20]['local_U_draw']=20.
    case['annual'][30]['local_U_principal']=20.
    revised.validate_finance(case,settings)
    case['annual'][30]['local_U_principal']=18.
    with pytest.raises(ValueError,match='stock/flow'):
        revised.validate_finance(case,settings)
    case['annual'][30]['local_U_principal']=20.
    case['annual'][40]['owner_custody']=3.
    with pytest.raises(ValueError,match='custody'):
        revised.validate_finance(case,settings)


@pytest.mark.parametrize('identity',('source_case_sha256','source_loop_case_sha256'))
def test_correct_product_with_wrong_case_identity_is_rejected(settings,identity):
    finance,assets=fixture_products()
    finance['cases'][0][identity]='f'*64
    with pytest.raises(ValueError,match='case identities differ'):
        revised.matched_cases(finance,assets,small_config(),settings)
    del finance['cases'][0][identity]
    with pytest.raises(ValueError,match='case identities differ'):
        revised.matched_cases(finance,assets,small_config(),settings)


@pytest.fixture(scope='module')
def built(settings):
    finance,assets=fixture_products();case=finance['cases'][0];asset=assets['cases'][0]
    baseline={'cases':[model.path(case,asset,settings,ownership)
                       for ownership in settings['ownership_scenarios']]}
    summaries=revised.baseline_checkpoints(baseline)
    original=json.dumps((finance,assets,settings),sort_keys=True)
    result=revised.build(small_config(),settings,finance,assets,summaries,
                         {'schema':'baseline','producer':{'sha256':'b'*64},
                          'unvalued_assets':['Unpriced intellectual property']})
    assert json.dumps((finance,assets,settings),sort_keys=True)==original
    return finance,assets,result


def test_same_inputs_have_zero_matched_baseline_changes_and_preserve_all_payouts(built):
    finance,_,result=built
    assert len(result['cases'])==3 and len(result['baseline_comparisons'])==15
    for comparison in result['baseline_comparisons']:
        assert all(v is None or v==0 for v in comparison['change'].values())
    for path in result['cases']:
        for flow,row in zip(finance['cases'][0]['annual'],path['annual']):
            assert row['flows']['local_payouts']+row['flows']['outside_payouts']==pytest.approx(flow['owner_payout'])
            assert row['book_balance']['controls']['passed']
            assert row['programme_cash']==flow['cash']
            assert row['programme_reserves']==flow['reserves']
            assert abs(row['money_stock_flow']['reconciliation_residual'])<1e-6
    # Detailed valuation is compacted only between the retained checkpoints.
    assert 'cash_claim_balance' not in result['cases'][0]['annual'][50]
    assert 'cash_claim_balance' in result['cases'][0]['checkpoints']['100']


def test_terminal_equity_includes_programme_cash_and_reserves_only_once(built):
    _,_,result=built
    path=next(p for p in result['cases'] if p['ownership']=='external_equity')
    row=path['checkpoints']['200'];summary=revised.checkpoint_summary(row)
    assert summary['enterprise_own_funds_equity']==1030.
    assert summary['outside_enterprise_equity_at_own_funds']==1030.
    cash=summary['outside_cash_NPV_at_3percent']
    assert summary['outside_cash_plus_terminal_claim_NPV_at_3percent']==pytest.approx(cash+1030/1.03**200)
    assert not result['controls']['terminal_equity_added_to_distributed_cash_returns']
    assert summary['monetary_base']==summary['bank_settlement_reserves']
    assert summary['consolidated_local_net_wealth']==pytest.approx(summary['household_net_wealth'])


def test_unreleased_owner_custody_is_separate_from_enterprise_equity(settings):
    finance,assets=fixture_products();case=finance['cases'][0]
    case['annual'][200].update(owner_payout=4.96,owner_custody=.04)
    case['annual'][201].update(owner_payout=5.04)
    revised.validate_finance(case,settings)
    path=model.path(case,assets['cases'][0],settings,'external_equity')
    summary=revised.checkpoint_summary(path['checkpoints']['200'])
    assert summary['outside_custody_claim']==pytest.approx(.04)
    assert summary['enterprise_own_funds_equity']==1030.
    terminal_pv=summary['outside_cash_plus_terminal_claim_NPV_at_3percent']-summary['outside_cash_NPV_at_3percent']
    assert terminal_pv==pytest.approx(1030.04/1.03**200,rel=1e-12,abs=1e-10)
    assert path['annual'][201]['flows']['outside_payouts']==pytest.approx(5.04)


def test_ownership_changes_claims_without_inventing_programme_income(built):
    _,_,result=built
    external=next(p for p in result['cases'] if p['ownership']=='external_equity')
    gifted=next(p for p in result['cases'] if p['ownership']=='negotiated_vesting')
    for y in ('32','100','200','220'):
        a=revised.checkpoint_summary(external['checkpoints'][y])
        b=revised.checkpoint_summary(gifted['checkpoints'][y])
        assert a['included_real_assets']==b['included_real_assets']
        assert a['gross_receipts_cumulative']==b['gross_receipts_cumulative']
        assert a['programme_payouts_cumulative']==pytest.approx(b['programme_payouts_cumulative'])
        assert b['outside_payouts_cumulative']+b['local_payouts_cumulative']==pytest.approx(a['programme_payouts_cumulative'])
    assert gifted['checkpoints']['100']['local_economic_ownership_fraction']==pytest.approx(.75)


def test_restored_service_value_changes_equity_not_cash_receipts(settings):
    finance,assets=fixture_products();case=finance['cases'][0];asset=assets['cases'][0]
    original=model.path(case,asset,settings,'negotiated_vesting')
    revised_asset=deepcopy(asset)
    for row in revised_asset['annual'][100:]:
        row['totals']['net_remaining_service_value_USD2041']+=100.
    changed=model.path(case,revised_asset,settings,'negotiated_vesting')
    a=revised.checkpoint_summary(original['checkpoints']['200'])
    b=revised.checkpoint_summary(changed['checkpoints']['200'])
    assert b['included_real_assets']-a['included_real_assets']==100.
    assert b['local_enterprise_equity_at_own_funds']-a['local_enterprise_equity_at_own_funds']==pytest.approx(75.)
    for key in ('gross_receipts_cumulative','programme_payouts_cumulative','outside_cash_NPV_at_3percent'):
        assert b[key]==a[key]


def test_historical_files_are_pinned_and_remain_unchanged():
    cfg=json.loads((revised.HERE/'revised_inputs.json').read_text())
    for key in ('settings_product','baseline_product'):
        spec=cfg[key]
        assert revised.sha(revised.HERE/spec['path'])==spec['sha256']
    for name,expected in cfg['implementation_dependencies'].items():
        assert revised.sha(revised.HERE/name)==expected
