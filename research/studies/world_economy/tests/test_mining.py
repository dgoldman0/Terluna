"""Finite inventories, accounting boundaries and independent financing controls."""

import math

import pytest

from research.studies.world_economy.mining import finite_export_screen


def _case(**changes):
    values = dict(
        recoverable_inventory_kg=1000,
        reserved_future_kg=200,
        protected_kg=100,
        uncertainty_buffer_kg=100,
        buyer_demand_cap_kg=450,
        years=5,
        annual_export_capacity_kg=100,
        price_U_per_kg=20,
        intermediate_cost_U_per_kg=5,
        cfc_U_per_kg=2,
        depletion_U_per_kg=3,
        other_senior_cash_U_per_kg=1,
        retained_growth_capex_U_per_kg=2,
        depletion_reserve_U_per_kg=4,
        independent_buyer_fraction=0.75,
    )
    values.update(changes)
    return finite_export_screen(**values)


def test_export_preserves_disjoint_reservations_and_inventory_identity():
    r = _case(buyer_demand_cap_kg=2000, annual_export_capacity_kg=1000)
    assert r["quantity_kg"] == 600
    assert r["remaining_inventory_kg"] == 400
    assert r["remaining_reserved_future_kg"] == 200
    assert r["remaining_protected_kg"] == 100
    assert r["remaining_uncertainty_buffer_kg"] == 100
    assert r["remaining_unreserved_inventory_kg"] == 0
    assert r["remaining_inventory_kg"] + r["quantity_kg"] == 1000


def test_full_reservation_leaves_zero_sales_and_no_cash():
    r = _case(reserved_future_kg=800)
    assert r["quantity_kg"] == r["sales_U"] == r["independent_project_cash_U"] == 0
    assert r["remaining_inventory_kg"] == r["remaining_combined_reservation_kg"] == 1000
    assert r["annual_uniform_independent_project_cash_U"] == [0] * 5


def test_overlapping_or_excess_reservations_are_rejected_not_silently_clipped():
    with pytest.raises(ValueError, match="combined reservations"):
        _case(reserved_future_kg=801)


def test_demand_and_capacity_can_independently_limit_exports():
    demand_limited = _case()
    assert demand_limited["quantity_kg"] == 450
    assert demand_limited["remaining_unreserved_inventory_kg"] == 150
    capacity_limited = _case(annual_export_capacity_kg=50)
    assert capacity_limited["quantity_kg"] == 250
    assert capacity_limited["annual_uniform_export_kg"] == [50] * 5


def test_nwp_and_cash_use_different_depletion_and_reinvestment_deductions():
    r = _case()
    assert r["sales_U"] == 450 * 20
    assert r["gross_value_added_U"] == 450 * (20 - 5)
    assert r["net_world_product_U"] == 450 * (20 - 5 - 2 - 3)
    assert r["project_cash_before_buyer_filter_U"] == 450 * (20 - 5 - 2 - 1 - 2 - 4)
    assert r["independent_project_cash_U"] == 2025
    assert r["excluded_programme_funded_cash_U"] == 675
    assert sum(r["annual_uniform_independent_project_cash_U"]) == 2025
    assert sum(r["annual_uniform_excluded_programme_funded_cash_U"]) == 675


def test_depletion_changes_nwp_without_subtracting_a_second_cash_reserve():
    base, changed = _case(), _case(depletion_U_per_kg=8)
    assert changed["net_world_product_U"] == base["net_world_product_U"] - 450 * 5
    assert changed["independent_project_cash_U"] == base["independent_project_cash_U"]
    reserve_changed = _case(depletion_reserve_U_per_kg=6)
    assert reserve_changed["net_world_product_U"] == base["net_world_product_U"]
    assert reserve_changed["project_cash_before_buyer_filter_U"] == 450 * 4


def test_programme_funded_purchases_produce_output_but_no_independent_finance():
    base, dependent = _case(), _case(independent_buyer_fraction=0)
    for key in ("quantity_kg", "sales_U", "gross_value_added_U", "net_world_product_U"):
        assert dependent[key] == base[key]
    assert dependent["independent_project_cash_U"] == 0
    assert dependent["excluded_programme_funded_cash_U"] == base["project_cash_before_buyer_filter_U"]


def test_negative_cash_remains_a_funding_requirement_in_both_allocations():
    r = _case(other_senior_cash_U_per_kg=10)
    assert r["net_world_product_U"] == 4500
    assert r["project_cash_before_buyer_filter_U"] == -1350
    assert r["independent_project_cash_U"] == -1012.5
    assert r["excluded_programme_funded_cash_U"] == -337.5
    assert all(c < 0 for c in r["annual_uniform_independent_project_cash_U"])


@pytest.mark.parametrize("changes", [
    {"buyer_demand_cap_kg": 0},
    {"annual_export_capacity_kg": 0},
    {"recoverable_inventory_kg": 0, "reserved_future_kg": 0,
     "protected_kg": 0, "uncertainty_buffer_kg": 0},
])
def test_zero_resource_demand_or_capacity_produces_zero_output(changes):
    r = _case(**changes)
    assert r["quantity_kg"] == r["gross_value_added_U"] == r["net_world_product_U"] == 0
    assert r["project_cash_before_buyer_filter_U"] == 0


@pytest.mark.parametrize("name,value", [
    ("recoverable_inventory_kg", math.inf),
    ("price_U_per_kg", math.nan),
    ("protected_kg", -1),
    ("depletion_reserve_U_per_kg", -1),
    ("independent_buyer_fraction", 1.01),
    ("independent_buyer_fraction", -0.01),
    ("years", 0),
    ("years", 2.5),
    ("years", True),
    ("price_U_per_kg", "20"),
])
def test_invalid_inputs_are_rejected(name, value):
    with pytest.raises(ValueError):
        _case(**{name: value})


def test_finite_inputs_that_overflow_totals_are_rejected():
    with pytest.raises(ValueError, match="finite arithmetic"):
        _case(price_U_per_kg=1e308)


def test_uniform_cash_series_and_stock_reservations_close_with_fractional_values():
    r = _case(years=7, independent_buyer_fraction=1 / 3,
              buyer_demand_cap_kg=413.2)
    assert math.fsum(r["annual_uniform_export_kg"]) == pytest.approx(r["quantity_kg"])
    assert math.fsum(r["annual_uniform_independent_project_cash_U"]) == pytest.approx(r["independent_project_cash_U"])
    assert r["independent_project_cash_U"] + r["excluded_programme_funded_cash_U"] == pytest.approx(r["project_cash_before_buyer_filter_U"])
    assert r["remaining_inventory_kg"] == pytest.approx(
        r["remaining_combined_reservation_kg"] + r["remaining_unreserved_inventory_kg"])
