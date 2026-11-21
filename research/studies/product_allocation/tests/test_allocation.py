import json
import math

import pytest

from research.studies.product_allocation import allocation as model


@pytest.fixture(scope='module')
def result():
    return model.product()


def test_result_reproduces_with_all_admitted_identities(result):
    assert result == model.read(model.HERE / 'results/allocation.json')
    assert len(result['market_products']) == 21
    assert len(result['refining_location']) == 54


def test_refining_threshold_reconciles_independent_route_cash_flows():
    for grade in (.001, .1, .8):
        for earth_recovery, space_recovery in ((.95, .9), (.95, .95), (.9, .98)):
            for freight in (10, 150, 1500):
                price, survival, earth_cost = 40000., .99, 25.
                threshold = model.maximum_refining_cost(grade, earth_recovery, space_recovery,
                                                        price, freight, survival, earth_cost)
                on_earth = survival * (grade * earth_recovery * price - earth_cost) - freight
                off_earth = survival * grade * space_recovery * price - freight * grade * space_recovery - threshold
                assert on_earth == pytest.approx(off_earth)


def test_better_recovery_can_outweigh_returned_mass_saving():
    # Losing five percent of the contained valuable material can reverse the route.
    assert model.maximum_refining_cost(.1, .95, .9, 40000, 150, .99, 0) < 0
    assert model.maximum_refining_cost(.1, .95, .95, 40000, 150, .99, 0) > 0


def test_mass_added_during_fabrication_has_correct_opportunity_cost():
    assert model.fabrication_mass_penalty(.005, 150) == pytest.approx(29850)
    assert model.fabrication_mass_penalty(1, 150) == 0
    with pytest.raises(ValueError):
        model.fabrication_mass_penalty(0, 150)


def test_missing_assays_and_oxide_units_do_not_create_invented_sales(result):
    lunar = next(r for r in result['normalized_recovery'] if r['profile'] == 'kreep15386' and r['scenario'] == 'middle')
    ci = next(r for r in result['normalized_recovery'] if r['profile'] == 'ci' and r['scenario'] == 'middle')
    assert lunar['groups_kg']['PGM6'] is None
    assert ci['groups_kg']['PGM6'] == pytest.approx(1715.7)
    assert lunar['groups_kg']['NdPr'] < lunar['groups_kg']['REY']
    assert 'Os' in ci['element_kg']
    for r in result['joint_circuit_budgets']:
        assert 'Os' not in r['individual_contributions_USD2026_per_feed_kg']
        assert not any('oxide' in k for k in r['individual_contributions_USD2026_per_feed_kg'])


def test_commodity_handling_failure_and_lower_freight_sensitivity_remain_distinct(result):
    nickel = next(r for r in result['market_products'] if r['id'] == 'Ni')
    costs = {r['freight_USD2026_kg']: r for r in nickel['transport_screens']}
    assert costs[150]['remaining_cost_budget_USD2026_per_kg_shipped'] < 0
    assert costs[10]['remaining_cost_budget_USD2026_per_kg_shipped'] > 0
    assert costs[150]['at_or_above_inherited_ground_floor']
    assert not costs[10]['at_or_above_inherited_ground_floor']


def test_joint_budget_is_once_per_feed_and_declines_with_freight(result):
    for price_factor in (.5, 1, 2):
        cases = [r for r in result['joint_circuit_budgets'] if r['price_multiplier'] == price_factor]
        budgets = [r['maximum_joint_incremental_cost_USD2026_per_feed_kg'] for r in cases]
        assert budgets == sorted(budgets, reverse=True)
        for row in cases:
            assert row['maximum_joint_incremental_cost_USD2026_per_feed_kg'] == pytest.approx(
                sum(row['individual_contributions_USD2026_per_feed_kg'].values()))
            assert all(v >= 0 for v in row['individual_contributions_USD2026_per_feed_kg'].values())


def test_checkpoint_opportunity_uses_existing_feed_and_books_no_revenue(result):
    cfg, assays, flows, recovery, _, _ = model.load()
    cases = {c['name']: c for c in flows['cases']}
    process = math.prod(recovery['recovery_scenarios']['middle'].values())
    assert len(result['annual_opportunity']) == 8
    for row in result['annual_opportunity']:
        source = next(p for p in cases[row['case']]['checkpoints'] if p['year'] == row['year'])
        expected = source['source_feed_kg_year'] * row['treated_share'] * assays['profiles']['ci']['elements']['Pt']['ppm'] * 1e-6 * process
        assert row['refined_element_kg_year']['Pt'] == pytest.approx(expected)
        assert row['receipts_booked_USD2041'] == 0
        # Return survival belongs only in the comparison with a delivered Earth market.
        assert row['ratio_to_2025_world_primary_output']['Pt'] == pytest.approx(expected * cfg['return_survival'] / 170000)


def test_price_unit_conversion_preserves_the_usgs_product_basis():
    register = model.read(model.HERE / 'market_benchmarks.json')
    p = {p['id']: p for p in register['products']}
    assert p['Ni']['USD2025_per_product_kg'] == 15
    assert p['Pt']['USD2025_per_product_kg'] * register['unit_conversions']['kg_per_troy_ounce'] == pytest.approx(1200)
    assert p['NdPr_oxide']['USD2025_per_product_kg'] == 69
    assert p['Sc_metal']['estimated_world_primary_production_2025_t'] is None
    assert p['Sc_oxide']['production_basis'] == 'scandium oxide'


def test_customer_targets_preserve_the_existing_purchase_account(result):
    previous = model.read(model.HERE / '../enclosed_returns/results/earth_purchases.json.gz')
    cases = {c['name']: c for c in previous['cases']}
    for row in result['existing_customer_requirements']:
        p = next(p for p in cases[row['case']]['checkpoints'] if p['horizon'] == row['year'])
        assert row['Earth_goods_kg_year'] == p['annual_at_checkpoint']['mass_kg']['Earth_purchased_goods']
        assert row['Earth_goods_invoices_USD2041_year'] == sum(p['annual_at_checkpoint']['revenue_by_sector_USD2041'][k] for k in ('materials', 'source_products'))
