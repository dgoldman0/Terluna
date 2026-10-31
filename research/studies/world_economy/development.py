"""Annual development cash flows with independent early operating cash.

Construction occurs at year end. Interest accrues on opening debt; all unpaid
interest is capitalized. Construction grants are fixed, restricted to that
year's construction, and identical in the early-income and deferred cases.
External net cash is already after operating costs, replacement expenditure,
tax/other senior cash obligations and investment retained by the operating
business. It is independently earned outside the ring-fenced project whose
debt is compared. Receipts funded by that project's own borrowing are a memo
only and never add to available cash. Outside buyers' ordinary financing is
not automatically excluded: the prohibition concerns circular project funding.

Cash can first pay interest or current construction, then prepay opening debt.
Unused cash stays in a zero-return reserve for later years, or is released
explicitly once current obligations and all debt have been paid. Negative
operating cash is rejected: its funding must be resolved before this loan model.
No GDP, CFC, market demand, asset availability or repayment contract is inferred.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from research.studies.world_economy.finance import debt_service_first_year


def _number(value, name, *, signed=False):
    value = float(value)
    if not math.isfinite(value) or (not signed and value < 0):
        raise ValueError(f"{name} must be finite" + ("" if signed else " and nonnegative"))
    return value


def _series(values, name, length=None, *, signed=False):
    result = [_number(v, name, signed=signed) for v in values]
    if not result or (length is not None and len(result) != length):
        raise ValueError(f"{name} must have the same nonzero length as construction")
    return result


def independent_business_cash(
    external_receipts: float,
    external_operating_costs: float,
    replacement_capex: float,
    retained_growth_capex: float,
    *,
    other_senior_cash_obligations: float = 0.0,
    own_project_debt_funded_receipts: float = 0.0,
) -> dict[str, float]:
    """Cash from the independent business, with circular receipts excluded.

    All costs refer to the external-service business represented by receipts.
    The own-project receipt is a separate memo; it is not inside
    ``external_receipts`` and must not appear in that argument as well. Spending
    on project construction and its suppliers is accounted in construction,
    not a second time among this business's external-service operating costs.
    """
    values = {k: _number(v, k) for k, v in {
        "external_receipts": external_receipts,
        "external_operating_costs": external_operating_costs,
        "replacement_capex": replacement_capex,
        "retained_growth_capex": retained_growth_capex,
        "other_senior_cash_obligations": other_senior_cash_obligations,
        "own_project_debt_funded_receipts": own_project_debt_funded_receipts,
    }.items()}
    values["independent_net_cash"] = (
        values["external_receipts"] - values["external_operating_costs"]
        - values["replacement_capex"] - values["retained_growth_capex"]
        - values["other_senior_cash_obligations"]
    )
    return values


def development_cashflows(
    annual_construction: Sequence[float],
    annual_external_net_cash: Sequence[float],
    annual_rate: float,
    *,
    annual_grants: Sequence[float] | None = None,
    initial_debt: float = 0.0,
    allocation: str = "interest_first",
    allow_prepayment: bool = True,
    annual_debt_financed_receipts: Sequence[float] | None = None,
    release_excess_cash: bool = False,
    mature_repayment_years: int | None = None,
) -> dict[str, object]:
    """Run one ring-fenced loan with fixed annual construction and grants.

    ``allocation`` is ``interest_first`` or ``draws_first``. Surplus cash then
    prepays debt if permitted, otherwise remains reserved. Project-debt-funded
    receipts are excluded memos; their corresponding borrowing must already
    be represented by construction or separately reported outside this loan.
    Negative operating cash is rejected. ``release_excess_cash`` pays out
    remaining cash only when closing debt is zero. Its annual ``cash_released``
    series may fund another project's cash flow once, without also leaving that
    cash in this project's reserves. ``excess_cash_retained`` is the positive
    annual reserve addition; ``closing_reserve`` is a stock, not annual income.
    An optional mature term reports next-year annuity service on closing debt.
    """
    construction = _series(annual_construction, "construction")
    n = len(construction)
    earnings = _series(annual_external_net_cash, "external net cash", n)
    grants = _series([0.0] * n if annual_grants is None else annual_grants, "grants", n)
    excluded = _series([0.0] * n if annual_debt_financed_receipts is None else annual_debt_financed_receipts,
                       "project-debt-funded receipts", n)
    rate = _number(annual_rate, "annual_rate")
    debt = _number(initial_debt, "initial_debt")
    if allocation not in {"interest_first", "draws_first"}:
        raise ValueError("allocation must be interest_first or draws_first")
    if not isinstance(allow_prepayment, bool):
        raise ValueError("allow_prepayment must be boolean")
    if not isinstance(release_excess_cash, bool):
        raise ValueError("release_excess_cash must be boolean")
    reserve = 0.0
    rows = []
    for year, (cost, net, grant, circular) in enumerate(zip(construction, earnings, grants, excluded), 1):
        if grant > cost:
            raise ValueError("A construction grant exceeds that year's construction")
        opening_debt, opening_reserve = debt, reserve
        interest = opening_debt * rate
        cash = opening_reserve + net
        required_construction = cost - grant
        if allocation == "interest_first":
            interest_paid = min(cash, interest)
            cash -= interest_paid
            construction_paid = min(cash, required_construction)
            cash -= construction_paid
        else:
            construction_paid = min(cash, required_construction)
            cash -= construction_paid
            interest_paid = min(cash, interest)
            cash -= interest_paid
        new_draw = required_construction - construction_paid
        capitalized_interest = interest - interest_paid
        principal_prepayment = min(cash, opening_debt) if allow_prepayment else 0.0
        cash -= principal_prepayment
        debt = opening_debt + new_draw + capitalized_interest - principal_prepayment
        cash_available_for_release = cash if debt == 0 else 0.0
        cash_released = cash_available_for_release if release_excess_cash else 0.0
        cash -= cash_released
        reserve = cash
        if not all(math.isfinite(v) for v in (debt, reserve, interest)):
            raise OverflowError("Development balance exceeds finite floating-point range")
        cash_sources = math.fsum([opening_reserve, net, grant, new_draw])
        cash_uses = math.fsum([cost, interest_paid, principal_prepayment, reserve, cash_released])
        debt_sources = math.fsum([opening_debt, new_draw, capitalized_interest])
        debt_uses = math.fsum([debt, principal_prepayment])
        rows.append(dict(
            year=year, opening_debt=opening_debt, opening_reserve=opening_reserve,
            construction=cost, construction_grant=grant, external_net_cash=net,
            excluded_project_debt_funded_receipts=circular,
            accrued_interest=interest, interest_paid_from_early_cash=interest_paid,
            construction_paid_from_early_cash=construction_paid,
            new_construction_borrowing=new_draw, capitalized_interest=capitalized_interest,
            principal_prepayment=principal_prepayment, closing_debt=debt,
            closing_reserve=reserve, net_debt=debt-reserve,
            early_cash_applied_to_project=interest_paid+construction_paid+principal_prepayment,
            cash_available_for_release=cash_available_for_release, cash_released=cash_released,
            excess_cash_retained=max(reserve-opening_reserve, 0.0),
            reserve_drawdown=max(opening_reserve-reserve, 0.0),
            cash_sources=cash_sources, cash_uses=cash_uses,
            cash_closure_residual=cash_sources-cash_uses,
            debt_closure_residual=debt_sources-debt_uses,
            cash_closure_relative_residual=(cash_sources-cash_uses)/max(1.0, abs(cash_sources), abs(cash_uses)),
            debt_closure_relative_residual=(debt_sources-debt_uses)/max(1.0, abs(debt_sources), abs(debt_uses)),
        ))
    summed = ("construction", "construction_grant", "external_net_cash",
              "excluded_project_debt_funded_receipts",
              "accrued_interest", "interest_paid_from_early_cash", "construction_paid_from_early_cash",
              "new_construction_borrowing", "capitalized_interest", "principal_prepayment",
              "early_cash_applied_to_project", "cash_released", "excess_cash_retained", "reserve_drawdown")
    totals = {k: math.fsum(row[k] for row in rows) for k in summed}
    debt_free_years = [row["year"] for row in rows if row["closing_debt"] == 0]
    last_debt_year = max((row["year"] for row in rows if row["closing_debt"] > 0), default=0)
    compounded_cash = 0.0
    for row in rows:
        compounded_cash = compounded_cash * (1+rate) + row["early_cash_applied_to_project"]
    final_service = None if mature_repayment_years is None else debt_service_first_year(debt, rate, mature_repayment_years)
    return dict(years=n, annual_rate=rate, initial_debt=initial_debt, allocation=allocation,
                allow_prepayment=allow_prepayment, annual=rows, totals=totals,
                release_excess_cash=release_excess_cash,
                closing_debt=debt, closing_reserve=reserve, net_debt=debt-reserve,
                debt_free_year_ends=debt_free_years,
                first_year_remaining_debt_free_through_horizon=last_debt_year+1 if last_debt_year<n else None,
                compounded_early_cash_applied=compounded_cash,
                final_annual_service=final_service,
                max_abs_cash_closure_residual=max(abs(r["cash_closure_residual"]) for r in rows),
                max_abs_debt_closure_residual=max(abs(r["debt_closure_residual"]) for r in rows),
                max_relative_cash_closure_residual=max(abs(r["cash_closure_relative_residual"]) for r in rows),
                max_relative_debt_closure_residual=max(abs(r["debt_closure_relative_residual"]) for r in rows))


def compare_with_deferred(
    annual_construction: Sequence[float],
    annual_external_net_cash: Sequence[float],
    annual_rate: float,
    **kwargs,
) -> dict[str, object]:
    """Compare identical construction, grants and rate; only early cash changes."""
    costs = list(annual_construction)
    early = development_cashflows(costs, annual_external_net_cash, annual_rate, **kwargs)
    baseline_kwargs = dict(kwargs)
    baseline_kwargs["annual_debt_financed_receipts"] = [0.0] * len(costs)
    deferred = development_cashflows(costs, [0.0] * len(costs), annual_rate, **baseline_kwargs)
    return dict(early=early, deferred=deferred,
                closing_debt_reduction=deferred["closing_debt"]-early["closing_debt"],
                net_debt_reduction=deferred["net_debt"]-early["net_debt"],
                interest_avoided=deferred["totals"]["accrued_interest"]-early["totals"]["accrued_interest"],
                debt_avoidance_identity_residual=deferred["closing_debt"]-early["closing_debt"]-early["compounded_early_cash_applied"],
                debt_avoidance_identity_relative_residual=(deferred["closing_debt"]-early["closing_debt"]-early["compounded_early_cash_applied"])/max(1.0, deferred["closing_debt"], early["compounded_early_cash_applied"]),
                comparison_rule="same construction schedule, construction grants, initial debt, annual rate and cash policy; deferred case has zero independent early operating cash")
