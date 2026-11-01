"""Transparent annual debt scenarios; all monetary inputs use the same units.

Rates are effective annual real rates. Draws occur at each year end; interest
accrues on the opening balance. A completed tranche's first repayment falls
one year after completion. Debt at full occupancy is measured immediately
after that year's payments and construction draws; reported mature service
is the following year's payment. These are cash-flow calculations, not
deductions from GDP or net domestic product.
"""

from __future__ import annotations

import math
from numbers import Integral


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _years(value: int, name: str, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        raise ValueError(f"{name} is outside its allowed range")
    return int(value)


def _fraction(value: float, name: str) -> float:
    value = _nonnegative(value, name)
    if value > 1:
        raise ValueError(f"{name} must not exceed one")
    return value


def annuity_payment(principal: float, annual_rate: float, years: int) -> float:
    """Equal end-year payment that amortizes a loan over ``years`` payments."""
    principal = _nonnegative(principal, "principal")
    rate = _nonnegative(annual_rate, "annual_rate")
    years = _years(years, "years")
    if rate == 0:
        return principal / years
    return principal * rate / -math.expm1(-years * math.log1p(rate))


def loan_balance_after_payments(
    principal: float, annual_rate: float, repayment_years: int, payments_made: int
) -> float:
    """Balance after completed annuity payments, with fully paid loans at zero.

    The remaining-payment present-value ratio avoids cancellation between a
    compounded principal and a nearly equal accumulated payment balance.
    """
    principal = _nonnegative(principal, "principal")
    rate = _nonnegative(annual_rate, "annual_rate")
    term = _years(repayment_years, "repayment_years")
    paid = _years(payments_made, "payments_made", allow_zero=True)
    remaining = max(0, term - paid)
    if rate == 0:
        return principal * remaining / term
    log_rate = math.log1p(rate)
    return principal * (
        -math.expm1(-remaining * log_rate)
        / -math.expm1(-term * log_rate)
    )


def debt_service_first_year(
    principal: float, annual_rate: float, repayment_years: int
) -> dict[str, float]:
    """Split the first payment into interest and repayment of principal."""
    principal = _nonnegative(principal, "principal")
    rate = _nonnegative(annual_rate, "annual_rate")
    payment = annuity_payment(principal, rate, repayment_years)
    interest = principal * rate
    principal_repaid = max(0.0, payment - interest)
    return {
        "payment": payment,
        "interest": interest,
        "principal_repaid": principal_repaid,
        "closing_debt": loan_balance_after_payments(
            principal, rate, repayment_years, 1
        ),
    }


def construction_finance(
    total_cost: float,
    construction_years: int,
    annual_rate: float,
    financed_fraction: float = 1.0,
    capitalized_interest_fraction: float = 1.0,
    occupancy_gap_years: int = 0,
) -> dict[str, float | int]:
    """Finance equal annual end-year construction draws, then an optional gap.

    The financed fraction applies to construction cost only. The specified
    fraction of every year's accrued interest is added to debt; the remainder
    is paid from contemporaneous resources. No principal is repaid before full
    occupancy in this single-project scenario. Gap years continue that same
    interest policy without new construction draws. A zero capitalization
    fraction therefore leaves principal unchanged during the gap.
    """
    cost = _nonnegative(total_cost, "total_cost")
    years = _years(construction_years, "construction_years")
    gap = _years(occupancy_gap_years, "occupancy_gap_years", allow_zero=True)
    rate = _nonnegative(annual_rate, "annual_rate")
    fraction = _fraction(financed_fraction, "financed_fraction")
    capitalized_fraction = _fraction(
        capitalized_interest_fraction, "capitalized_interest_fraction"
    )
    financed_cost = cost * fraction
    annual_draw = financed_cost / years
    balance = 0.0
    construction_interest = []
    gap_interest = []
    completion_debt = 0.0
    for year in range(1, years + gap + 1):
        interest = balance * rate
        (construction_interest if year <= years else gap_interest).append(interest)
        balance += interest * capitalized_fraction
        if year <= years:
            balance += annual_draw
        if not math.isfinite(balance):
            raise OverflowError("Debt exceeds finite floating-point range")
        if year == years:
            completion_debt = balance
    construction_accrued = math.fsum(construction_interest)
    gap_accrued = math.fsum(gap_interest)
    accrued = construction_accrued + gap_accrued
    return {
        "construction_years": years,
        "occupancy_gap_years": gap,
        "financed_construction": financed_cost,
        "unfinanced_construction": cost - financed_cost,
        "annual_construction_draw": annual_draw,
        "debt_at_completion": completion_debt,
        "debt_at_full_occupancy": balance,
        "construction_accrued_interest": construction_accrued,
        "construction_capitalized_interest": construction_accrued * capitalized_fraction,
        "construction_current_interest": construction_accrued * (1 - capitalized_fraction),
        "gap_accrued_interest": gap_accrued,
        "gap_capitalized_interest": gap_accrued * capitalized_fraction,
        "gap_current_interest": gap_accrued * (1 - capitalized_fraction),
        "total_accrued_interest": accrued,
        "total_capitalized_interest": accrued * capitalized_fraction,
        "total_current_interest": accrued * (1 - capitalized_fraction),
    }


def rolling_tranches(
    tranches: list[dict[str, float | int]], maturity_year: int
) -> dict[str, object]:
    """Combine loans completed at different dates without reviving retired debt.

    Each item specifies ``principal``, ``completion_year``, ``annual_rate`` and
    ``repayment_years``. Principal is the liability at completion, including
    any construction interest already capitalized. Completion year zero is
    allowed; every completion must precede or equal the maturity snapshot.
    The first payment occurs at completion_year + 1. The output includes
    repayments made by the snapshot and debt service in its following year.
    """
    maturity = _years(maturity_year, "maturity_year", allow_zero=True)
    rows = []
    for tranche in tranches:
        principal = _nonnegative(tranche["principal"], "principal")
        completion = _years(tranche["completion_year"], "completion_year", allow_zero=True)
        if completion > maturity:
            raise ValueError("A tranche completes after maturity_year")
        rate = _nonnegative(tranche["annual_rate"], "annual_rate")
        term = _years(tranche["repayment_years"], "repayment_years")
        paid = min(maturity - completion, term)
        balance = loan_balance_after_payments(principal, rate, term, paid)
        original_payment = annuity_payment(principal, rate, term)
        next_payment = original_payment if paid < term else 0.0
        next_interest = balance * rate
        repaid = principal - balance
        cumulative_payment = original_payment * paid
        rows.append({
            "principal_at_completion": principal,
            "completion_year": completion,
            "payments_made": paid,
            "debt_at_maturity": balance,
            "cumulative_payment": cumulative_payment,
            "cumulative_principal_repaid": repaid,
            "cumulative_interest_paid": max(0.0, cumulative_payment - repaid),
            "payment_first_mature_year": next_payment,
            "interest_first_mature_year": next_interest,
            "principal_repaid_first_mature_year": max(0.0, next_payment - next_interest),
        })
    sum_fields = (
        "principal_at_completion", "debt_at_maturity", "cumulative_payment",
        "cumulative_principal_repaid", "cumulative_interest_paid",
        "payment_first_mature_year", "interest_first_mature_year",
        "principal_repaid_first_mature_year",
    )
    result = {field: math.fsum(row[field] for row in rows) for field in sum_fields}
    result.update({"maturity_year": maturity, "tranches": rows})
    return result


def staged_construction_finance(
    total_cost: float,
    construction_years: int,
    tranche_years: int,
    annual_rate: float,
    repayment_years: int,
    financed_fraction: float = 1.0,
    capitalized_interest_fraction: float = 1.0,
    occupancy_gap_years: int = 0,
) -> dict[str, object]:
    """Construct equal-rate tranches and amortize each after its own completion.

    Annual real construction expenditure is constant; a shorter last tranche
    costs proportionally less. These projects must earn revenue or receive
    transfers before final world occupancy to fund their early repayments.
    This scenario never treats those early payments as free resources.
    """
    cost = _nonnegative(total_cost, "total_cost")
    years = _years(construction_years, "construction_years")
    tranche_length = _years(tranche_years, "tranche_years")
    gap = _years(occupancy_gap_years, "occupancy_gap_years", allow_zero=True)
    term = _years(repayment_years, "repayment_years")
    rate = _nonnegative(annual_rate, "annual_rate")
    fraction = _fraction(financed_fraction, "financed_fraction")
    cap_fraction = _fraction(capitalized_interest_fraction, "capitalized_interest_fraction")
    loans = []
    construction_rows = []
    for start in range(0, years, tranche_length):
        duration = min(tranche_length, years - start)
        finance = construction_finance(
            cost * duration / years, duration, rate, fraction, cap_fraction
        )
        construction_rows.append(finance)
        loans.append({
            "principal": finance["debt_at_completion"],
            "completion_year": start + duration,
            "annual_rate": rate,
            "repayment_years": term,
        })
    result = rolling_tranches(loans, years + gap)
    result.update({
        "financed_construction": cost * fraction,
        "unfinanced_construction": cost * (1 - fraction),
        "construction_capitalized_interest": math.fsum(
            row["total_capitalized_interest"] for row in construction_rows
        ),
        "construction_current_interest": math.fsum(
            row["total_current_interest"] for row in construction_rows
        ),
    })
    result["cash_paid_before_maturity"] = (
        result["unfinanced_construction"]
        + result["construction_current_interest"]
        + result["cumulative_payment"]
    )
    return result
