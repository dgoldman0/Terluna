"""Physical and cash controls for the inventory comparison."""
import math

import pytest

from research.studies.world_economy import industrial_stocks as stocks


def case(**changes):
    settings = dict(initial_stock=50, annual_demand=1, horizon_years=50,
        acquisition_cost_per_kg=1, future_buy_cost_per_kg=1.5)
    settings.update(changes)
    return stocks.evaluate_stock(**settings)


def test_loss_free_stock_and_cash_closure():
    result = case()
    totals = result['totals']
    assert totals['stock_withdrawals_kg'] == 50
    assert totals['terminal_physical_kg'] == 0
    assert totals['physical_loss_kg'] == 0
    assert totals['pv_stock_plus_shortfall_U'] == 50
    assert totals['pv_all_just_in_time_U'] == 75
    assert totals['break_even_acquisition_cost_U_per_kg'] == 1.5
    assert totals['original_manufacture_gva_recounted_on_withdrawal_U'] == 0


def test_compatibility_loss_preserves_mass_and_requires_new_purchases():
    result = case(initial_stock=10, annual_demand=6, horizon_years=1,
                  physical_attrition=.1, compatibility_decay=.5)
    row = result['annual'][0]
    assert row['physical_loss_kg'] == 1
    assert row['newly_incompatible_kg'] == 4.5
    assert row['withdrawal_kg'] == 4.5
    assert row['just_in_time_shortfall_purchase_kg'] == 1.5
    assert row['closing_physical_kg'] == row['closing_incompatible_kg'] == 4.5
    assert result['totals']['physical_closure_residual_kg'] == 0


@pytest.mark.parametrize('decay', [0, .005, .05])
def test_required_stock_serves_each_order_and_conserves_mass(decay):
    demands = [1, 3, 2, 4, 7] * 10
    required = stocks.required_initial_stock(demands, 50, .0005, decay, .8)
    result = case(initial_stock=required, annual_demand=demands,
                  physical_attrition=.0005, compatibility_decay=decay,
                  initial_compatible_fraction=.8)
    totals = result['totals']
    assert totals['stock_withdrawals_kg'] == pytest.approx(sum(demands))
    assert totals['just_in_time_shortfall_purchases_kg'] == pytest.approx(0, abs=1e-10)
    assert totals['terminal_compatible_kg'] == pytest.approx(0, abs=1e-10)
    assert totals['physical_closure_residual_kg'] == pytest.approx(0, abs=1e-10)
    for row in result['annual']:
        assert row['physical_closure_residual_kg'] == pytest.approx(0, abs=1e-10)
        assert row['demand_closure_residual_kg'] == pytest.approx(0)
        assert row['withdrawal_kg'] <= row['compatible_before_withdrawal_kg']


def test_required_stock_matches_closed_form():
    q = (1-.0005)*(1-.005)
    assert stocks.required_initial_stock(1, 50, .0005, .005) == pytest.approx(sum(q**-year for year in range(1, 51)))
    assert stocks.required_initial_stock(3, 20) == 60


def test_storage_on_opening_stock_and_discounted_future_prices():
    result = case(initial_stock=2, horizon_years=2,
        storage_cost_per_kg_year=.2, real_discount_rate=.1,
        annual_price_decline=.1, preposition_cost_per_kg=.3)
    totals = result['totals']
    expected_storage = .4/1.1 + .2/1.1**2
    expected_jit = 1.5*.9/1.1 + 1.5*.9**2/1.1**2
    assert totals['pv_storage_U'] == pytest.approx(expected_storage)
    assert totals['pv_all_just_in_time_U'] == pytest.approx(expected_jit)
    assert totals['pv_stock_plus_shortfall_U'] == pytest.approx(2+.6+expected_storage)


@pytest.mark.parametrize('storage_settings', [dict(storage_cost_per_kg_year=.02), dict(storage_fraction_initial_acquisition_value=.02)])
def test_break_even_acquisition_reproduces_equal_cost(storage_settings):
    settings = dict(initial_stock=10, horizon_years=12, annual_demand=1,
        physical_attrition=.001, compatibility_decay=.005,
        real_discount_rate=.005, preposition_cost_per_kg=.1,
        annual_price_decline=.003, **storage_settings)
    first = case(**settings)
    threshold = first['totals']['break_even_acquisition_cost_U_per_kg']
    assert threshold > 0
    at_threshold = case(**settings, acquisition_cost_per_kg=threshold)
    assert at_threshold['totals']['pv_saving_vs_just_in_time_U'] == pytest.approx(0, abs=1e-12)


def test_terminal_reference_value_does_not_fund_current_costs():
    result = case(initial_stock=100, annual_demand=0, horizon_years=1,
                  compatibility_decay=.5, storage_cost_per_kg_year=.002)
    totals = result['totals']
    assert totals['terminal_physical_kg'] == 100
    assert totals['terminal_compatible_kg'] == 50
    assert totals['terminal_incompatible_kg'] == 50
    assert totals['terminal_compatible_market_reference_value_U'] == 75
    assert totals['terminal_liquidation_credit_in_pv_U'] == 0
    assert totals['pv_stock_plus_shortfall_U'] == 100.2
    assert totals['break_even_acquisition_cost_U_per_kg'] < 0


def test_equal_service_horizon_comparison_includes_no_reserve():
    compared = stocks.evaluate_reserve_horizons(1, (5, 20, 50), 50,
        acquisition_cost_per_kg=1, future_buy_cost_per_kg=1,
        storage_cost_per_kg_year=.002, real_discount_rate=.005)
    assert compared['cheapest_tested_reserve_years'] == 0
    assert [row['reserve_years'] for row in compared['cases']] == [0, 5, 20, 50]
    for row in compared['cases']:
        totals = row['result']['totals']
        assert totals['stock_withdrawals_kg'] + totals['just_in_time_shortfall_purchases_kg'] == pytest.approx(50)
    assert compared['cases'][0]['result']['totals']['break_even_acquisition_cost_U_per_kg'] is None


@pytest.mark.parametrize('settings', [dict(physical_attrition=1.1), dict(compatibility_decay=-.1),
    dict(initial_stock=float('nan')), dict(real_discount_rate=-.01),
    dict(horizon_years=1.5), dict(annual_price_decline=1.1),
    dict(storage_cost_per_kg_year=.1, storage_fraction_initial_acquisition_value=.01)])
def test_reject_invalid_settings(settings):
    with pytest.raises(ValueError):
        case(**settings)


def test_zero_survival_and_zero_demand_are_distinct():
    with pytest.raises(ValueError):
        stocks.required_initial_stock(1, 50, compatibility_decay=1)
    assert stocks.required_initial_stock(0, 50, compatibility_decay=1) == 0
    result = case(physical_attrition=1)
    assert result['totals']['physical_loss_kg'] == 50
    assert result['totals']['just_in_time_shortfall_purchases_kg'] == 50


def test_nonrepresentable_stock_requirement_rejected():
    with pytest.raises(ValueError):
        stocks.required_initial_stock(1, 1000, compatibility_decay=.99)
