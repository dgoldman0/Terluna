import json
import math

import pytest
from research.studies.world_economy import amenities, finance


@pytest.fixture(scope='module')
def product():
    return amenities.results()


def test_committed_extension_reproduces_without_changing_baseline(product):
    assert product == json.loads((amenities.HERE/'results/amenities.json').read_text())
    _, inputs = amenities.load_inputs()
    baseline = inputs['research/studies/world_economy/results/comparison.json']
    assert baseline['production']['close_parity']['worlds']['moon']['gross_product'] == pytest.approx(839.2e12)
    for case in product['tourism'].values():
        assert all(w['baseline_gwp_adjustment_U'] == 0 for w in case['worlds'].values())


def test_changed_product_requires_rebinding(tmp_path):
    settings, _ = amenities.load_inputs()
    first = next(iter(settings['products']))
    settings['products'][first]['sha256'] = '0'*64
    path = tmp_path/'inputs.json'
    path.write_text(json.dumps(settings))
    with pytest.raises(ValueError, match='Input changed'):
        amenities.load_inputs(path)


def test_tourism_trade_closes_and_capacities_bind(product):
    for case in product['tourism'].values():
        assert len(case['routes']) == 6
        assert abs(case['closure']['three_world_combined_net_external_cash_U']) < 1
        for world in case['worlds'].values():
            assert world['origin_budget']['resident_trip_spend_U'] <= world['origin_budget']['outbound_budget_U']+1
            assert world['inbound_capacity']['annual_capacity_utilization'] <= 1
    strong = product['tourism']['lunar_destination']
    assert strong['worlds']['earth']['origin_budget']['unspent_outbound_budget_U'] > 20e12
    assert product['tourism']['balanced']['worlds']['moon']['trade']['tourism_net_external_cash_U'] < 0
    assert product['tourism']['lunar_destination_double_capacity']['worlds']['moon']['financing_comparison']['tourism_fraction_of_positive_requirement'] < 1


def test_development_conserves_cash_claims_and_reproduces_prior_baseline(product):
    _, inputs = amenities.load_inputs()
    baseline = inputs['research/studies/world_economy/results/comparison.json']
    for case in product['development_scenarios']:
        assert abs(case['cross_world_transfer_closure_U']) < 1
        for world in case['worlds'].values():
            scale = world['debt_at_full_population_without_early_cash_U']
            assert abs(world['debt_avoidance_identity_residual_U']) < scale*1e-12
            assert world['annual_cash_closure_max_U'] < scale*1e-12
            assert 0 <= world['debt_reduction_share'] <= 1+1e-12
        # Loan principal/interest move resident resources between worlds only.
        assert sum(case['mature_cash_resources_U'].values()) == pytest.approx(
            sum(w['net_world_product'] for w in baseline['reference_financing']['legacy_retired']['accounts']['worlds'].values()))
    ref = product['development_reference']
    for name in ('moon','array'):
        history = baseline['reference_financing']['patient']['finance_history'][name]
        assert ref['worlds'][name]['debt_at_full_population_without_early_cash_U'] == pytest.approx(history['debt_at_full_population'])


def test_array_transfers_are_new_cash_released_once(product):
    rows = [r for r in product['development_scenarios'] if r['activity_case']=='phased' and r['cash_share']==.05 and r['annual_rate']==.001]
    own = next(r for r in rows if r['array_transfer_share']==0)
    shared = next(r for r in rows if r['array_transfer_share']==1)
    transfer = shared['worlds']['array']['cumulative_paid_transfer_U']
    assert transfer == pytest.approx(own['worlds']['array']['cumulative_released_cash_retained_U'])
    assert shared['worlds']['array']['cumulative_released_cash_retained_U'] == 0
    assert shared['worlds']['moon']['cumulative_received_transfer_U'] == transfer
    assert shared['worlds']['moon']['debt_at_full_population_with_early_cash_U'] < own['worlds']['moon']['debt_at_full_population_with_early_cash_U']
    assert shared['worlds']['array']['debt_at_full_population_with_early_cash_U'] == own['worlds']['array']['debt_at_full_population_with_early_cash_U'] == 0


def test_late_population_has_matching_baseline_age(product):
    ref = product['development_reference']
    late = next(r for r in product['development_scenarios'] if r['activity_case']=='enclosed_then_open' and r['cash_share']==.05 and r['annual_rate']==.001 and r['array_transfer_share']==0)
    assert late['full_population_year'] == 750
    for w in ('moon','array'):
        assert late['worlds'][w]['debt_at_full_population_without_early_cash_U'] == pytest.approx(ref['worlds'][w]['debt_at_full_population_without_early_cash_U']*1.001**250)
    assert late['worlds']['moon']['checkpoints'][-2]['activity_fraction'] == .005


def test_idle_stock_care_has_funding_and_does_not_raise_fixed_sponsor_grants(product):
    reference, stress = product['development_reference'], product['idle_stock_funding_stress']
    for w in ('moon','array'):
        r,s = reference['worlds'][w],stress['worlds'][w]
        assert s['cumulative_sponsor_funding_U'] == r['cumulative_sponsor_funding_U']
        assert s['cumulative_extra_idle_stock_care_U'] > 0
        assert s['debt_at_full_population_without_early_cash_U'] > r['debt_at_full_population_without_early_cash_U']
    assert stress['worlds']['moon']['debt_reduction_share'] < reference['worlds']['moon']['debt_reduction_share']


@pytest.mark.parametrize('knots', [[[0,0],[10,.2],[9,1]], [[1,0],[10,1]], [[0,0],[10,1.1]]])
def test_invalid_activity_profile_rejected(knots):
    with pytest.raises(ValueError):
        amenities.activity_curve(knots,10)


def test_mining_reservations_nutrients_and_cash_are_separate(product):
    assert product['mining']['geological_export_authorization_established'] is False
    for row in product['mining']['rows']:
        x = row['feed_equivalent_account']
        assert x['remaining_inventory_kg'] >= x['remaining_combined_reservation_kg']
        assert row['exported_phosphorus_kg'] <= 1e9*(1+1e-12)
        assert abs(row['feed_mass_closure_kg']) <= max(x['quantity_kg'],1)*1e-12
        assert x['independent_project_cash_U']+x['excluded_programme_funded_cash_U'] == pytest.approx(x['project_cash_before_buyer_filter_U'])
        if row['independent_buyer_fraction'] == 0:
            assert row['standalone_debt_avoided_at_year500_U'] == 0
        assert row['standalone_reduction_share_of_patient_moon_debt'] < .001
