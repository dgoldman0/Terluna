"""Reproduction and cross-account checks on admitted financial results."""
import json
from research.studies.enclosed_returns import finance as f


def saved():return json.loads((f.HERE/'results/returns.json').read_text())


def test_saved_result_reproduces_from_admitted_products():
    expected=saved();actual=f.product(*f.inputs())
    assert actual==expected


def test_coverage_and_resource_conditions_are_explicit():
    d=saved();assert not d['unresolved'];assert len(d['checkpoints'])==87*3*3
    keys={(r['name'],r['mode'],r['horizon']) for r in d['checkpoints']}
    assert len(keys)==len(d['checkpoints'])
    for r in d['checkpoints']:
        assert r['resources']['qualified_source_supply'] is None
        assert len(r['resources']['conditional_mass_allowances'])==7
        assert r['seed_capital_cost_allowance']>0
        assert r['unpriced_deployment_requirements']


def test_financing_does_not_create_output_or_social_benefit():
    d=saved();by={(r['name'],r['mode'],r['horizon']):r for r in d['checkpoints']}
    for r in d['checkpoints']:
        equity=by[r['name'],'equity',r['horizon']]
        assert r['annual_gross']==equity['annual_gross']
        assert r['annual_operating_available']==equity['annual_operating_available']
        for rate in r['discounted']:
            assert r['discounted'][rate]['full_price_socioeconomic_NPV']==equity['discounted'][rate]['full_price_socioeconomic_NPV']


def test_planned_distributions_never_cause_later_unfunded_commitments():
    for path in saved()['reference_finance_paths']:
        seen_payout=False
        for r in path['annual']:
            if seen_payout:assert r['contribution']<=1.
            if r['dividend']>1.:seen_payout=True
