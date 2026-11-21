"""Cross-product controls for the stratified account and bounded interventions."""
from copy import deepcopy
import gzip
import json
import math
import pytest
from research.studies.enclosed_returns import bottlenecks as b, bottleneck_orders as o, markets, finance


@pytest.fixture(scope='module')
def result():
    return json.loads((b.HERE/'results/bottlenecks.json').read_text())


@pytest.fixture(scope='module')
def experiment_orders():
    return json.loads((b.HERE/'results/bottleneck_orders.json').read_text())


def test_experimental_controls_match_every_baseline_request_and_population(experiment_orders):
    reference=json.loads((b.HERE/'results/orders.json').read_text())
    by={c['name']:c for c in reference['cases']}
    for case in experiment_orders['cases']:
        base=by[case['experiment']['baseline_name']]
        if case['experiment']['variant'] in ('baseline','tug_life_30'):
            assert case['annual']==base['annual']
        for row,ref in zip(case['annual'],base['annual']):
            for key in ('year','population','commissioned_people','new_residential_people'):
                assert row[key]==ref[key]
    assert o.product()==experiment_orders


def test_demand_changes_preserve_audience_and_unit_quotes():
    cfg,q,d,staged,_,_=markets.inputs();before=deepcopy((q,staged))
    qq,bb=o.demand_variant(q,staged,'breakthrough',1.25)
    assert (q,staged)==before
    x=bb['scenarios']['breakthrough'];old=staged['scenarios']['breakthrough']
    assert x['audience']==old['audience']==2e9
    assert x['paid_fraction']==.5
    for key in ('paid_USD2041_year','license_USD2041_year','research_USD2041_contract','professional_USD2041_client_year'):
        assert x[key]==old[key]
    assert qq['strategies']==q['strategies']
    with pytest.raises(ValueError,match='participation'):
        o.demand_variant(q,staged,'breakthrough',3.)


def test_all_baseline_financial_accounts_stay_equal_to_committed_product(result):
    reference=json.loads((b.HERE/'results/returns.json').read_text())
    by={(r['name'],r['horizon']):r for r in reference['checkpoints'] if r['mode']=='U_credit'}
    assert len(result['baseline_stratification'])==81
    for case in result['baseline_stratification']:
        for r in case['checkpoints']:
            original={k:v for k,v in by[case['name'],r['horizon']].items() if k!='resources'}
            assert r['checkpoint']==original
            for d in r['discounted'].values():
                assert abs(d['baseline_record_residual'])<.1


def test_every_period_cost_stratum_reconciles_to_cash_before_financing(result):
    for c in result['baseline_stratification']:
        for h in (100,200,220):
            periods=[p for p in c['periods'] if p['last_year']<=h]
            expected=next(r['checkpoint'] for r in c['checkpoints'] if r['horizon']==h)
            for rate in ('0.0','0.03'):
                actual=math.fsum(p['values'][rate]['net_real_cash'] for p in periods)
                assert actual==pytest.approx(expected['discounted'][rate]['unlevered_cash_NPV'],rel=1e-12,abs=.1)
            bridge=c['financing_bridges'][str(h)]
            assert abs(bridge['reconciliation_residual_USD2041'])<.25


def test_capacity_failures_never_receive_completed_return_metrics(result):
    failed=[r for r in result['interventions'] if r['status']=='physical_requirements_unresolved']
    assert len(failed)==7
    for r in failed:
        assert 'checkpoints' not in r
        assert r['workforce_failures'] and r['first_capacity_shortfall']==r['workforce_failures'][0]['year']
        assert all(x['required_jobs']>x['available_jobs'] for x in r['workforce_failures'])
    for r in result['interventions']:
        assert all(v['qualified_source_supply'] is None for v in r['resources'])
    by={r['portfolio']:r for r in result['intervention_differences'] if r['horizon']==200}
    assert by['diversified']['demand_rollout_interaction'] is None
    assert by['services']['demand_rollout_interaction'] is not None


def test_intervention_baselines_preserve_all_financial_numbers(result):
    by={c['name']:c for c in result['baseline_stratification']}
    for c in result['interventions']:
        if c['experiment']['variant']!='baseline':continue
        for row,reference in zip(c['checkpoints'],by[c['experiment']['baseline_name']]['checkpoints']):
            actual=deepcopy(row);actual['checkpoint']['name']=reference['checkpoint']['name']
            assert actual==reference


def test_irr_roundoff_screen_does_not_modify_cash_accounts(result):
    x=next(c for c in result['baseline_stratification'] if c['name']=='prior_manufacturing_mix:breakthrough:reference')
    r=next(r for r in x['checkpoints'] if r['horizon']==200)
    assert r['checkpoint']['pooled_funder_IRR_if_unique'] is None
    assert r['real_IRR_roundoff_screened']==pytest.approx(.05746931,abs=1e-8)
    assert r['checkpoint']['additional_equity_required']<1.


def test_margin_target_has_no_free_output_or_unpriced_seed_change(result):
    bases={c['name']:c for c in result['baseline_stratification']}
    for case in result['commercial_margin_targets']:
        baseline=bases[case['name']]
        assert case['socioeconomic_change'] is None
        assert all(x['cost_reduction']<=x['covered_variable_purchases'] for x in case['annual_cost_reductions'])
        for row,ref in zip(case['checkpoints'],baseline['checkpoints']):
            actual,old=row['checkpoint'],ref['checkpoint']
            assert all('full_price_socioeconomic_NPV' not in d for d in actual['discounted'].values())
            assert row['initial_full_funding']==ref['initial_full_funding']
            assert actual['annual_gross']==old['annual_gross']
            covered=sum(old['revenue_by_sector'][k] for k in ('sports','entertainment','licenses','research','professional'))
            assert actual['annual_operating_available']-old['annual_operating_available']==pytest.approx(.05*covered)
            assert actual['annual_capital']==old['annual_capital']


def test_saved_stratification_reproduces_exactly(result):
    assert b.product()==result
