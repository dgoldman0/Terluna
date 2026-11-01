"""Numerical controls for the study's construction and debt cash-flow cases."""

import math

import pytest

from research.studies.world_economy import finance


def test_zero_rate_and_tiny_rate_have_correct_limits():
    assert finance.annuity_payment(500, 0, 1000) == 0.5
    assert finance.annuity_payment(500, 1e-14, 1000) == pytest.approx(0.5, rel=1e-10)
    result = finance.construction_finance(1000, 500, 0, 0.7, 1, 100)
    assert result["debt_at_full_occupancy"] == pytest.approx(700)
    assert result["total_accrued_interest"] == 0


@pytest.mark.parametrize("rate,term", [(0, 100), (0.001, 1000), (0.005, 250), (0.02, 500)])
def test_discounted_installments_close_loan_and_balances_conserve_debt(rate, term):
    principal = 123.45
    payment = finance.annuity_payment(principal, rate, term)
    pv = math.fsum(payment / (1 + rate) ** year for year in range(1, term + 1))
    assert pv == pytest.approx(principal, rel=2e-12)
    for paid in (0, term // 2, term - 1):
        before = finance.loan_balance_after_payments(principal, rate, term, paid)
        after = finance.loan_balance_after_payments(principal, rate, term, paid + 1)
        assert before * (1 + rate) - payment == pytest.approx(after, abs=2e-12)
    assert finance.loan_balance_after_payments(principal, rate, term, term) == 0
    assert finance.loan_balance_after_payments(principal, rate, term, term + 100) == 0


def test_partial_interest_capitalization_conserves_draws_and_interest():
    result = finance.construction_finance(1000, 250, 0.01, 0.6, 0.4, 40)
    assert result["debt_at_full_occupancy"] == pytest.approx(
        result["financed_construction"] + result["total_capitalized_interest"]
    )
    assert result["total_accrued_interest"] == pytest.approx(
        result["total_capitalized_interest"] + result["total_current_interest"]
    )
    assert result["debt_at_full_occupancy"] == pytest.approx(
        result["debt_at_completion"] * (1 + 0.4 * 0.01) ** 40
    )


def test_current_interest_only_financing_has_no_hidden_capitalization():
    result = finance.construction_finance(1000, 500, 0.02, 1, 0, 25)
    assert result["debt_at_full_occupancy"] == pytest.approx(1000)
    assert result["construction_current_interest"] == pytest.approx(0.02 * 1000 * 499 / 2)
    assert result["gap_current_interest"] == pytest.approx(0.02 * 1000 * 25)


def test_five_century_capitalization_matches_geometric_sum():
    result = finance.construction_finance(1, 500, 0.02)
    expected = math.expm1(500 * math.log1p(0.02)) / (500 * 0.02)
    assert expected > 1900
    assert result["debt_at_completion"] == pytest.approx(expected, rel=2e-13)
    service = finance.debt_service_first_year(expected, 0.02, 1000)
    assert service["payment"] >= 0.02 * expected
    assert service["interest"] + service["principal_repaid"] == pytest.approx(service["payment"])


def test_maturity_excludes_retired_tranches_and_preserves_surviving_loans():
    result = finance.rolling_tranches([
        {"principal": 100, "completion_year": 0, "annual_rate": 0, "repayment_years": 100},
        {"principal": 100, "completion_year": 50, "annual_rate": 0, "repayment_years": 100},
        {"principal": 100, "completion_year": 100, "annual_rate": 0, "repayment_years": 100},
    ], 100)
    assert result["debt_at_maturity"] == 150
    assert result["payment_first_mature_year"] == 2
    assert result["cumulative_payment"] == 150
    assert result["tranches"][0]["payment_first_mature_year"] == 0


def test_staging_reduces_final_debt_by_exposing_payments_during_construction():
    result = finance.staged_construction_finance(1000, 500, 100, 0, 100)
    assert result["debt_at_maturity"] == 200
    assert result["cumulative_principal_repaid"] == 800
    assert result["cash_paid_before_maturity"] == 800
    assert result["payment_first_mature_year"] == 2
    assert result["financed_construction"] + result["construction_capitalized_interest"] == pytest.approx(
        result["debt_at_maturity"] + result["cumulative_principal_repaid"]
    )
    assert sum(row["payment_first_mature_year"] > 0 for row in result["tranches"]) == 1


def test_uneven_last_tranche_and_interest_are_accounted_for():
    result = finance.staged_construction_finance(1000, 250, 100, 0.005, 150, 0.8, 0.5, 20)
    assert [row["completion_year"] for row in result["tranches"]] == [100, 200, 250]
    assert result["financed_construction"] + result["construction_capitalized_interest"] == pytest.approx(
        result["debt_at_maturity"] + result["cumulative_principal_repaid"]
    )
    assert result["cumulative_payment"] == pytest.approx(
        result["cumulative_interest_paid"] + result["cumulative_principal_repaid"]
    )


@pytest.mark.parametrize("args", [(-1, 0.01, 100), (1, -0.01, 100), (1, 0.01, 0), (1, 0.01, 1.5)])
def test_invalid_debt_inputs_fail_explicitly(args):
    with pytest.raises(ValueError):
        finance.annuity_payment(*args)
