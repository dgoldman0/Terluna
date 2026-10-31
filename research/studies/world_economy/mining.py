"""Finite mineral-export accounting under explicit inventory and demand assumptions.

This screen does not measure reserves, establish extraction feasibility, infer a
market price, or allocate ownership. Reserved future use, protection and an
uncertainty buffer are disjoint quantities within recoverable inventory.
"""

from math import fsum, isfinite
from numbers import Real


def _nonnegative_finite(name, value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite nonnegative number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite nonnegative number") from exc
    if not isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite nonnegative number")
    return result


def _finite_result(name, value):
    if not isfinite(value):
        raise ValueError(f"{name} exceeds finite arithmetic")
    return value


def finite_export_screen(
    recoverable_inventory_kg,
    reserved_future_kg,
    protected_kg,
    uncertainty_buffer_kg,
    buyer_demand_cap_kg,
    years,
    annual_export_capacity_kg,
    price_U_per_kg,
    intermediate_cost_U_per_kg,
    cfc_U_per_kg,
    depletion_U_per_kg,
    other_senior_cash_U_per_kg,
    retained_growth_capex_U_per_kg,
    depletion_reserve_U_per_kg,
    independent_buyer_fraction,
):
    """Return finite exports, production and cash with uniform annual delivery.

    All masses refer to the same saleable material and inventory boundary; prices
    and costs use the same real unit U. The demand cap is cumulative over the
    positive integer number of years. The three reservations must sum to no more
    than inventory. Delivery is limited by available inventory, buyer demand and
    annual capacity times years.

    GVA is sales less intermediate inputs. NWP additionally deducts CFC and
    depletion. Cash instead deducts intermediate inputs, replacement equal to
    CFC, senior obligations, retained growth capital and the explicit depletion
    reserve. Depletion is not also subtracted from cash. Replacement equal to CFC
    is a simplifying assumption, not a general accounting identity.

    Independent buyers supply the specified fraction of project cash; buyers
    financed by the same construction programme supply the excluded fraction.
    This filter changes neither real production nor its associated NWP. Negative
    cash is retained as a funding requirement, including for excluded sales.
    Uniform delivery and cash timing are assumptions, not a commissioning plan.
    """
    raw = {
        "recoverable_inventory_kg": recoverable_inventory_kg,
        "reserved_future_kg": reserved_future_kg,
        "protected_kg": protected_kg,
        "uncertainty_buffer_kg": uncertainty_buffer_kg,
        "buyer_demand_cap_kg": buyer_demand_cap_kg,
        "annual_export_capacity_kg": annual_export_capacity_kg,
        "price_U_per_kg": price_U_per_kg,
        "intermediate_cost_U_per_kg": intermediate_cost_U_per_kg,
        "cfc_U_per_kg": cfc_U_per_kg,
        "depletion_U_per_kg": depletion_U_per_kg,
        "other_senior_cash_U_per_kg": other_senior_cash_U_per_kg,
        "retained_growth_capex_U_per_kg": retained_growth_capex_U_per_kg,
        "depletion_reserve_U_per_kg": depletion_reserve_U_per_kg,
        "independent_buyer_fraction": independent_buyer_fraction,
    }
    p = {name: _nonnegative_finite(name, value) for name, value in raw.items()}
    if isinstance(years, bool) or not isinstance(years, int) or years <= 0:
        raise ValueError("years must be a positive integer")
    fraction = p["independent_buyer_fraction"]
    if fraction > 1:
        raise ValueError("independent_buyer_fraction must lie between zero and one")

    inventory = p["recoverable_inventory_kg"]
    # Check reservations before summing; this also rejects an overflowing sum.
    reservations = (
        p["reserved_future_kg"], p["protected_kg"], p["uncertainty_buffer_kg"]
    )
    try:
        reserved = fsum(reservations)
    except OverflowError as exc:
        raise ValueError("combined reservations exceed recoverable inventory") from exc
    if reserved > inventory:
        raise ValueError("combined reservations exceed recoverable inventory")
    available = max(inventory - reserved, 0.0)
    try:
        capacity = years * p["annual_export_capacity_kg"]
    except OverflowError as exc:
        raise ValueError("cumulative capacity exceeds finite arithmetic") from exc
    _finite_result("cumulative capacity", capacity)
    quantity = min(available, p["buyer_demand_cap_kg"], capacity)

    per_kg_gva = p["price_U_per_kg"] - p["intermediate_cost_U_per_kg"]
    per_kg_nwp = per_kg_gva - p["cfc_U_per_kg"] - p["depletion_U_per_kg"]
    per_kg_cash = (
        per_kg_gva
        - p["cfc_U_per_kg"]
        - p["other_senior_cash_U_per_kg"]
        - p["retained_growth_capex_U_per_kg"]
        - p["depletion_reserve_U_per_kg"]
    )
    for name, value in (
        ("per-kg GVA", per_kg_gva),
        ("per-kg NWP", per_kg_nwp),
        ("per-kg cash", per_kg_cash),
    ):
        _finite_result(name, value)

    sales = _finite_result("sales", quantity * p["price_U_per_kg"])
    gva = _finite_result("GVA", quantity * per_kg_gva)
    nwp = _finite_result("NWP", quantity * per_kg_nwp)
    cash = _finite_result("project cash", quantity * per_kg_cash)
    independent_cash = cash * fraction
    excluded_cash = cash * (1 - fraction)
    return {
        "schema": "terluna.world-economy.finite-export-screen/1",
        "evidence": (
            "Conditional inventory, reservation, demand, capacity, price and cost "
            "scenario. No mineral-resource measurement, market feasibility, "
            "extraction feasibility or ownership allocation is established."
        ),
        "quantity_kg": quantity,
        "sales_U": sales,
        "gross_value_added_U": gva,
        "net_world_product_U": nwp,
        "project_cash_before_buyer_filter_U": cash,
        "independent_project_cash_U": independent_cash,
        "excluded_programme_funded_cash_U": excluded_cash,
        "remaining_inventory_kg": inventory - quantity,
        "remaining_reserved_future_kg": p["reserved_future_kg"],
        "remaining_protected_kg": p["protected_kg"],
        "remaining_uncertainty_buffer_kg": p["uncertainty_buffer_kg"],
        "remaining_combined_reservation_kg": reserved,
        "remaining_unreserved_inventory_kg": max(available - quantity, 0.0),
        "export_limits_kg": {
            "unreserved_inventory": available,
            "buyer_demand": p["buyer_demand_cap_kg"],
            "cumulative_capacity": capacity,
        },
        "per_kg": {
            "gross_value_added_U": per_kg_gva,
            "net_world_product_U": per_kg_nwp,
            "project_cash_before_buyer_filter_U": per_kg_cash,
        },
        "years": years,
        "annual_uniform_export_kg": [quantity / years] * years,
        "annual_uniform_independent_project_cash_U": [independent_cash / years] * years,
        "annual_uniform_excluded_programme_funded_cash_U": [excluded_cash / years] * years,
    }
