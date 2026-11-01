"""Accounting identity, classification and cross-world transfer controls."""

import copy

import pytest

from research.studies.world_economy import accounts, finance


def _worlds():
    return {
        "moon": {"sector_gva": {"services": 80, "industry": 40},
                 "ordinary_cfc": 15, "engineered_cfc": 5,
                 "natural_depletion": 2, "essential_final_operations": 8},
        "earth": {"sector_gva": {"services": 100, "industry": 50},
                  "ordinary_cfc": 20, "engineered_cfc": 3,
                  "natural_depletion": 4, "essential_final_operations": 6},
    }


def test_real_costs_have_one_deduction_and_operating_allocation_does_not_change_nwp():
    result = accounts.consolidated_accounts(_worlds())
    moon = result["worlds"]["moon"]
    assert moon["gross_world_product"] == 120
    assert moon["net_world_product"] == 98
    assert moon["resources_after_operations_and_principal"] == 90
    inputs = _worlds()
    inputs["moon"]["essential_final_operations"] += 7
    changed = accounts.consolidated_accounts(inputs)["worlds"]["moon"]
    assert changed["net_world_product"] == 98
    assert changed["net_resident_income"] == 98
    assert changed["resources_after_operations_and_principal"] == 83


def test_interest_profits_and_principal_have_distinct_accounts_and_creditor_counterparts():
    transfers = [
        {"payer": "moon", "recipient": "earth", "kind": "interest", "amount": 4},
        {"payer": "moon", "recipient": "earth", "kind": "principal", "amount": 6},
        {"payer": "earth", "recipient": "moon", "kind": "profit", "amount": 2},
    ]
    result = accounts.consolidated_accounts(_worlds(), transfers)
    moon, earth = result["worlds"]["moon"], result["worlds"]["earth"]
    assert (moon["net_world_product"], earth["net_world_product"]) == (98, 123)
    assert (moon["net_resident_income"], earth["net_resident_income"]) == (96, 125)
    assert (moon["cash_resources_after_principal"], earth["cash_resources_after_principal"]) == (90, 131)
    assert result["combined"]["net_world_product"] == 221
    assert result["combined"]["net_resident_income"] == 221
    assert result["combined"]["cash_resources_after_principal"] == 221
    assert set(result["closure_residuals"].values()) == {0}


def test_internal_debt_does_not_reduce_world_product_income_or_resources():
    baseline = accounts.consolidated_accounts(_worlds())
    flows = accounts.debt_transfers("moon", "moon", 7, 11)
    result = accounts.consolidated_accounts(_worlds(), flows)
    assert result["combined"] == baseline["combined"]
    assert result["worlds"]["moon"]["internal_interest_payments"] == 7
    assert result["worlds"]["moon"]["internal_principal_payments"] == 11
    assert result["worlds"]["moon"]["external_interest_payments"] == 0


def test_outside_system_flows_remain_and_close_without_an_invented_internal_creditor():
    flows = accounts.debt_transfers("moon", None, 10, 20, 0.5)
    flows.append({"payer": None, "recipient": "earth", "kind": "profit", "amount": 3})
    result = accounts.consolidated_accounts(_worlds(), flows)
    assert result["combined"]["net_resident_income"] == 219
    assert result["combined"]["cash_resources_after_principal"] == 209
    assert result["net_receipts_from_outside_system"]["interest"] == -5
    assert result["net_receipts_from_outside_system"]["principal"] == -10
    assert set(result["closure_residuals"].values()) == {0}


def test_full_population_age_preserves_remaining_term_and_excludes_repaid_loans():
    aged = accounts.mature_service_transfers("moon", "earth", 100, 0, 100, 75)
    assert aged["debt_at_maturity"] == 25
    assert aged["payment_first_mature_year"] == 1
    retired = accounts.mature_service_transfers("moon", "earth", 100, 0.01, 100, 150)
    assert retired["debt_at_maturity"] == 0
    assert retired["transfers"] == []


def test_staged_service_enters_world_accounts_without_restarting_amortization():
    staged = finance.staged_construction_finance(1000, 500, 50, 0.01, 100)
    flows = accounts.debt_transfers(
        "moon", "earth", staged["interest_first_mature_year"],
        staged["principal_repaid_first_mature_year"], 0.4,
    )
    baseline = accounts.consolidated_accounts(_worlds())
    result = accounts.consolidated_accounts(_worlds(), flows)
    service = staged["payment_first_mature_year"] * 0.4
    assert result["worlds"]["moon"]["cash_resources_after_principal"] == pytest.approx(98 - service)
    assert result["worlds"]["earth"]["cash_resources_after_principal"] == pytest.approx(123 + service)
    assert result["combined"]["cash_resources_after_principal"] == pytest.approx(
        baseline["combined"]["cash_resources_after_principal"]
    )


def test_matched_worlds_have_identical_accounts_and_negative_stress_results_are_preserved():
    identical = {"sector_gva": {"production": 100}, "ordinary_cfc": 15,
                 "engineered_cfc": 3, "natural_depletion": 2,
                 "essential_final_operations": 4}
    result = accounts.consolidated_accounts({"moon": identical, "earth": copy.deepcopy(identical)})
    assert result["worlds"]["moon"] == result["worlds"]["earth"]
    identical["engineered_cfc"] = 120
    result = accounts.consolidated_accounts({"moon": identical})
    assert result["worlds"]["moon"]["net_world_product"] == -37
    assert result["worlds"]["moon"]["resources_after_operations_and_principal"] == -41


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -1])
def test_invalid_costs_fail(invalid):
    inputs = _worlds()
    inputs["moon"]["ordinary_cfc"] = invalid
    with pytest.raises(ValueError):
        accounts.consolidated_accounts(inputs)


def test_generic_operating_bills_and_unnamed_creditors_are_rejected():
    inputs = _worlds()
    inputs["moon"]["producer_operating_costs"] = 10
    with pytest.raises(ValueError, match="Unknown world-account inputs"):
        accounts.consolidated_accounts(inputs)
    with pytest.raises(ValueError, match="endpoint"):
        accounts.consolidated_accounts(_worlds(), [
            {"payer": "moon", "recipient": "unknown", "kind": "interest", "amount": 10}
        ])
