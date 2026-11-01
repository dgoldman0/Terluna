import math

import pytest

from research.studies.world_economy.development import (compare_with_deferred, development_cashflows,
                         independent_business_cash)


def test_deferred_matches_end_year_draws_with_fixed_grants():
    x = development_cashflows([100, 100, 100], [0, 0, 0], .1, annual_grants=[20, 20, 20])
    assert x["closing_debt"] == pytest.approx(80 * (1.1**2 + 1.1 + 1))
    assert x["totals"]["construction_grant"] == 60


def test_early_profit_avoids_borrowing_and_its_future_interest():
    x = compare_with_deferred([100, 100, 0], [50, 0, 0], .1)
    assert x["closing_debt_reduction"] == pytest.approx(50 * 1.1**2)
    assert x["interest_avoided"] == pytest.approx(10.5)


def test_receipts_funded_by_the_project_are_not_independent_cash():
    x = compare_with_deferred([100, 100], [0, 0], .1, annual_debt_financed_receipts=[100, 100])
    assert x["closing_debt_reduction"] == 0
    assert x["early"]["totals"]["excluded_project_debt_funded_receipts"] == 200


def test_business_keeps_operating_and_replacement_and_growth_costs_first():
    x = independent_business_cash(100, 45, 20, 15,
                                  other_senior_cash_obligations=5,
                                  own_project_debt_funded_receipts=1000)
    assert x["independent_net_cash"] == 15


def test_priority_changes_classification_without_creating_cash():
    a = development_cashflows([100], [50], .1, initial_debt=100, allocation="interest_first")
    b = development_cashflows([100], [50], .1, initial_debt=100, allocation="draws_first")
    assert a["annual"][0]["interest_paid_from_early_cash"] == 10
    assert b["annual"][0]["interest_paid_from_early_cash"] == 0
    assert a["closing_debt"] == b["closing_debt"] == 160


def test_surplus_prepays_debt_and_cannot_make_debt_negative():
    x = development_cashflows([0, 0], [1000, 0], .1, initial_debt=100)
    assert x["annual"][0]["principal_prepayment"] == 100
    assert x["closing_debt"] == 0
    assert x["closing_reserve"] == 890


def test_retained_surplus_can_fund_later_construction():
    x = development_cashflows([0, 100], [150, 0], 0)
    assert x["closing_debt"] == 0
    assert x["closing_reserve"] == 50
    assert x["annual"][1]["construction_paid_from_early_cash"] == 100


def test_no_prepayment_keeps_cash_and_continues_interest():
    x = development_cashflows([0, 0], [200, 0], .1, initial_debt=100, allow_prepayment=False)
    assert x["closing_debt"] == 100
    assert x["closing_reserve"] == 180
    assert x["totals"]["accrued_interest"] == 20


def test_operating_losses_need_explicit_support_outside_construction_loan():
    with pytest.raises(ValueError):
        development_cashflows([100], [-20], .1)


@pytest.mark.parametrize("allocation", ["interest_first", "draws_first"])
def test_multiyear_cash_and_debt_close(allocation):
    costs = [100, 120, 90, 0, 20]
    x = development_cashflows(costs, [10, 80, 300, 100, 0], .03,
                              annual_grants=[20, 10, 0, 0, 5], allocation=allocation)
    for row in x["annual"]:
        assert math.isclose(row["cash_sources"], row["cash_uses"], rel_tol=1e-13, abs_tol=1e-12)
        assert abs(row["debt_closure_residual"]) < 1e-10
    assert x["closing_debt"] + x["totals"]["principal_prepayment"] == pytest.approx(
        x["totals"]["new_construction_borrowing"] + x["totals"]["capitalized_interest"])


def test_no_extra_grant_assumed_for_early_income_case():
    x = compare_with_deferred([100, 100], [90, 90], .1, annual_grants=[20, 20])
    assert x["early"]["totals"]["construction_grant"] == x["deferred"]["totals"]["construction_grant"] == 40


def test_released_array_cash_can_fund_moon_without_remaining_on_array():
    array = development_cashflows([100, 100, 0], [0, 300, 50], .1, release_excess_cash=True)
    released = [r["cash_released"] for r in array["annual"]]
    assert released == pytest.approx([0, 90, 50])
    assert array["closing_reserve"] == 0
    assert array["debt_free_year_ends"] == [2, 3]
    assert array["first_year_remaining_debt_free_through_horizon"] == 2
    moon = compare_with_deferred([100, 100, 0], released, .1)
    assert moon["closing_debt_reduction"] == pytest.approx(90*1.1+50)
    assert array["totals"]["external_net_cash"] == pytest.approx(
        array["totals"]["early_cash_applied_to_project"] + sum(released))


def test_compounded_cash_identity_and_final_annuity():
    x = compare_with_deferred([100, 100, 100, 0], [5, 30, 60, 15], .02,
                             annual_grants=[50, 50, 50, 0], mature_repayment_years=100)
    assert x["debt_avoidance_identity_residual"] == pytest.approx(0, abs=1e-12)
    d = x["early"]["closing_debt"]
    assert x["early"]["final_annual_service"]["payment"] == pytest.approx(d*.02/(1-1.02**-100))
    assert x["early"]["totals"]["external_net_cash"] == pytest.approx(
        x["early"]["totals"]["early_cash_applied_to_project"] + x["early"]["closing_reserve"])


def test_century_scale_and_large_currency_values_close_relatively():
    x = compare_with_deferred([1e15]*500, [0]*100+[2e14]*400, .01, annual_grants=[5e14]*500)
    assert x["early"]["max_relative_cash_closure_residual"] < 1e-14
    assert x["early"]["max_relative_debt_closure_residual"] < 1e-14
    assert abs(x["debt_avoidance_identity_relative_residual"]) < 1e-12


@pytest.mark.parametrize("kwargs", [
    {"annual_grants": [101]}, {"annual_grants": [-1]}, {"allocation": "unknown"},
    {"annual_grants": [0, 0]}, {"annual_rate": -1}, {"initial_debt": float("nan")},
])
def test_invalid_inputs_rejected(kwargs):
    arguments = dict(annual_construction=[100], annual_external_net_cash=[0], annual_rate=.1)
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        development_cashflows(**arguments)
