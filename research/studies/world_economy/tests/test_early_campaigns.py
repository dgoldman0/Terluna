import json
import pytest
from research.studies.world_economy import early_campaigns as m


def test_product_reproduces():
    assert m.results()==json.loads((m.HERE/'results/early_campaigns.json').read_text())


def test_bound_source_change_rejected(tmp_path):
    s,_=m.load_inputs();s['products'][next(iter(s['products']))]['sha256']='0'*64
    p=tmp_path/'changed.json';p.write_text(json.dumps(s))
    with pytest.raises(ValueError,match='Input changed'):m.load_inputs(p)


def test_reference_contract_capital_funding_and_independence():
    r=m.results()['reference']
    assert r['delivered_qualified_mass_kg']==pytest.approx(14.805e12)
    assert r['delivered_qualified_mass_kg']+r['remaining_orders_kg']==pytest.approx(15e12)
    assert r['initial_capital_U']>19e12
    assert r['peak_funding_with_customer_advance_U']<7e12
    assert r['unlevered_project_npv_U']>0
    assert r['sustained_capital_recovery_year']==23
    assert r['illustrative_independent_fraction_of_positive_final_surplus_U']==0


def test_ten_year_and_modest_advantage_losses_are_preserved():
    rows=m.results()['scenarios']
    assert all(r['unlevered_project_npv_U']<0 for r in rows if r['years']==10 or r['cost_case']=='modest_site_advantage')


def test_financing_never_changes_production_or_project_npv_in_grid():
    rows=m.results()['scenarios']
    paired={}
    for r in rows:
        key=tuple(r[k] for k in ('raw_capacity_kg_year','cost_case','years','buyer_budget_U_year'))
        if key in paired:
            a=paired[key]
            for field in ('delivered_qualified_mass_kg','unlevered_project_npv_U','final_cash_after_all_capital_U'):
                assert r[field]==a[field]
        else:paired[key]=r


def test_reference_cash_includes_factory_recovery_and_wind_down_once():
    r=m.results()['reference']
    value=-r['initial_capital_U']+sum(x['operating_cash_before_capital_U']-x['actual_replacement_U'] for x in r['annual'])+r['terminal_recovery_U']-r['decommission_U']
    assert r['final_cash_after_all_capital_U']==pytest.approx(value)
    assert r['customer_finance_cash_closure_U']==pytest.approx(0,abs=.1)
