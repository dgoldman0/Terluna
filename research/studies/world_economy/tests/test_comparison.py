import json
import math
import pytest
from research.studies.world_economy import compare


@pytest.fixture(scope="module")
def product():
    return compare.results()


def test_committed_comparison_reproduces(product):
    assert product == json.loads((compare.HERE/'results/comparison.json').read_text())


def test_each_financing_case_conserves_three_world_claims(product):
    for row in product['reference_financing'].values():
        ledger=row['accounts']
        scale=abs(ledger['combined']['gross_world_product'])
        for residual in ledger['closure_residuals'].values():
            assert abs(residual)<scale*1e-12
        assert abs(row['required_net_exports_sum_U'])<scale*1e-12
        for flow in row['trade_requirements'].values():
            assert abs(flow['resource_pool_closure_residual_U'])<scale*1e-12


def test_debt_changes_income_without_changing_territorial_product(product):
    base=product['reference_financing']['legacy_retired']['accounts']['worlds']
    for scenario in product['reference_financing'].values():
        for name, row in scenario['accounts']['worlds'].items():
            assert row['gross_world_product']==base[name]['gross_world_product']
            assert row['net_world_product']==base[name]['net_world_product']
    stress=product['reference_financing']['deferred_stress']['accounts']['worlds']
    assert stress['moon']['cash_resources_after_principal']<0
    assert stress['array']['cash_resources_after_principal']<0


def test_water_changes_debt_but_not_maintained_endowment_depreciation(product):
    rows=[r for r in product['sensitivities']['cost_and_water_under_patient_finance'] if r['cost_case']=='advanced_industrial']
    assert len(set(r['worlds']['moon']['net_world_product'] for r in rows))==1
    cash=[r['worlds']['moon']['cash_resources_after_principal'] for r in rows]
    assert cash[0]>cash[1]>cash[2]


def test_gas_power_subset_and_film_replacement_are_not_added_to_gva(product):
    for case in product['operating_funding_memo']['cases'].values():
        assert all(a<b for a,b in zip(case['lifting_hydrogen_energy_subset_U_year'],case['delivered_electricity_U_year']))
    sensitivity=product['sensitivities']['collective_operating_allocation']
    for world in ('earth','moon','array'):
        assert len(set(r['worlds'][world]['net_world_product'] for r in sensitivity))==1
        assert len(set(r['worlds'][world]['cash_resources_after_principal'] for r in sensitivity))==1
