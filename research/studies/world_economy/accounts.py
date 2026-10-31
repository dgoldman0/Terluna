"""Additive mature-world accounts and a separate resource-allocation budget.

Study conventions
-----------------
* Every number is an annual flow in the same common-price unit. Sector GVA
  already excludes intermediate consumption. Net product taxes, if included,
  must be supplied as their own GVA-account adjustment in ``sector_gva``.
* ``ordinary_cfc`` and ``engineered_cfc`` are disjoint economic depreciation
  allowances. The latter contains only additional engineered assets omitted
  from the ordinary allowance. There is no second replacement-cost deduction.
* ``natural_depletion`` is separate from depreciation and from delivered
  material purchases. It represents resource value lost in this territory.
* ``essential_final_operations`` is an explicitly chosen operating component
  of essential collective final consumption. It excludes intermediate inputs
  of market producers, depreciation, replacement/new capital expenditure,
  income transfers, and debt service. It is not the full SNA valuation of all
  collective services. Subtracting it yields a study-defined allocation
  residual, not a second definition of net world product.
* Interest and profits affect resident income, not territorial product.
  Principal changes the cash resource balance, not product or income. This
  cash interpretation assumes replacement expenditure equals the depreciation
  reserve and no unlisted borrowing, asset sales, transfers or additional
  investment. It is not a complete balance-of-payments forecast.
* Flows have explicit payer and recipient worlds. ``None`` means outside the
  modeled system. Same-world payments are recorded but cancel on consolidation;
  this does not cancel the underlying contracts or distributional obligations.
* The same equations and asset boundaries apply to future Earth and the Moon.
  An Earth debt-free case must be a named scenario, not an accounting default
  inferred from its inherited infrastructure. Age and staged debt payments are
  calculated in ``finance`` before their flows enter these accounts.

These conventions prevent mechanical duplicate deductions; the caller must
still classify assets and expenditure consistently. No prices, costs, rates,
asset lives, ownership shares or productivity scenarios are imposed here.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

from . import finance


_COST_FIELDS = (
    "ordinary_cfc", "engineered_cfc", "natural_depletion",
    "essential_final_operations",
)
_FLOW_KINDS = ("interest", "profit", "other_primary_income", "principal")
_PRIMARY_KINDS = ("interest", "profit", "other_primary_income")


def _finite(value: float, name: str, *, nonnegative: bool = False) -> float:
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"{name} must be finite" + (" and nonnegative" if nonnegative else ""))
    return value


def _sum(values: Iterable[float], name: str) -> float:
    return _finite(math.fsum(values), name)


def consolidated_accounts(
    worlds: Mapping[str, Mapping[str, object]],
    transfers: Iterable[Mapping[str, object]] = (),
) -> dict[str, object]:
    """Compute world accounts, their sum, and independent consolidation checks.

    Each world has ``sector_gva`` (a mapping of sector name to value added)
    and optional zero-default cost fields named in the module documentation.
    Sector GVA can be negative, allowing loss-making sectors and net subsidies;
    all cost allocations must be nonnegative. Negative net/residual outcomes
    are preserved as scenario failures rather than clipped to zero.

    A transfer specifies ``payer``, ``recipient``, ``kind`` and ``amount``.
    Kinds are interest, profit, other_primary_income and principal. Positive
    amounts flow from payer to recipient; endpoints must name supplied worlds
    or be None for the unmodeled rest of the system.
    """
    if not worlds:
        raise ValueError("At least one world is required")
    rows = {}
    transfer_terms = {}
    for world, inputs in worlds.items():
        if not isinstance(world, str) or not world:
            raise ValueError("World identifiers must be nonempty strings")
        unknown = set(inputs) - {"sector_gva", *_COST_FIELDS}
        if unknown:
            raise ValueError(f"Unknown world-account inputs: {sorted(unknown)}")
        gva = inputs.get("sector_gva")
        if not isinstance(gva, Mapping):
            raise ValueError("sector_gva must map sectors to value added")
        sectors = {str(k): _finite(v, f"{world}.{k}") for k, v in gva.items()}
        row = {field: _finite(inputs.get(field, 0), field, nonnegative=True) for field in _COST_FIELDS}
        row["sector_gva"] = sectors
        row["gross_world_product"] = _sum(sectors.values(), "gross_world_product")
        row["total_cfc"] = _sum((row["ordinary_cfc"], row["engineered_cfc"]), "total_cfc")
        row["net_world_product"] = _sum(
            (row["gross_world_product"], -row["total_cfc"], -row["natural_depletion"]),
            "net_world_product",
        )
        rows[world] = row
        transfer_terms[world] = {
            f"{scope}_{kind}_{direction}": []
            for scope in ("external", "internal")
            for kind in _FLOW_KINDS
            for direction in ("receipts", "payments")
        }

    flows = []
    outside_net = {kind: [] for kind in _FLOW_KINDS}
    for transfer in transfers:
        if set(transfer) != {"payer", "recipient", "kind", "amount"}:
            raise ValueError("A transfer needs payer, recipient, kind and amount")
        payer, recipient, kind = transfer["payer"], transfer["recipient"], transfer["kind"]
        if payer is None and recipient is None:
            raise ValueError("A transfer must involve at least one modeled world")
        if any(endpoint is not None and endpoint not in rows for endpoint in (payer, recipient)):
            raise ValueError("A transfer endpoint does not name a modeled world")
        if kind not in _FLOW_KINDS:
            raise ValueError(f"Unknown transfer kind: {kind}")
        amount = _finite(transfer["amount"], "transfer amount", nonnegative=True)
        scope = "internal" if payer == recipient else "external"
        if payer is not None:
            transfer_terms[payer][f"{scope}_{kind}_payments"].append(amount)
        if recipient is not None:
            transfer_terms[recipient][f"{scope}_{kind}_receipts"].append(amount)
        if payer is None:
            outside_net[kind].append(amount)
        if recipient is None:
            outside_net[kind].append(-amount)
        flows.append({"payer": payer, "recipient": recipient, "kind": kind, "amount": amount})

    for world, row in rows.items():
        row.update({field: _sum(values, field) for field, values in transfer_terms[world].items()})
        row["net_primary_income_from_other_worlds"] = _sum(
            (
                sign * row[f"external_{kind}_{direction}"]
                for kind in _PRIMARY_KINDS
                for direction, sign in (("receipts", 1), ("payments", -1))
            ),
            "net_primary_income_from_other_worlds",
        )
        row["net_resident_income"] = _sum(
            (row["net_world_product"], row["net_primary_income_from_other_worlds"]),
            "net_resident_income",
        )
        row["net_principal_receipts"] = _sum(
            (row["external_principal_receipts"], -row["external_principal_payments"]),
            "net_principal_receipts",
        )
        row["cash_resources_after_principal"] = _sum(
            (row["net_resident_income"], row["net_principal_receipts"]),
            "cash_resources_after_principal",
        )
        row["resources_after_operations_and_principal"] = _sum(
            (row["cash_resources_after_principal"], -row["essential_final_operations"]),
            "resources_after_operations_and_principal",
        )

    fields = (
        "gross_world_product", "ordinary_cfc", "engineered_cfc", "total_cfc",
        "natural_depletion", "net_world_product", "net_primary_income_from_other_worlds",
        "net_resident_income", "net_principal_receipts", "cash_resources_after_principal",
        "essential_final_operations", "resources_after_operations_and_principal",
    )
    combined = {field: _sum((row[field] for row in rows.values()), field) for field in fields}
    outside = {kind: _sum(values, kind) for kind, values in outside_net.items()}
    outside_primary = _sum((outside[kind] for kind in _PRIMARY_KINDS), "outside primary income")
    checks = {
        "product_identity": _sum((
            combined["gross_world_product"], -combined["total_cfc"],
            -combined["natural_depletion"], -combined["net_world_product"],
        ), "product closure"),
        "income_consolidation": _sum((
            combined["net_resident_income"], -combined["net_world_product"],
            -outside_primary,
        ), "income closure"),
        "principal_consolidation": _sum((
            combined["cash_resources_after_principal"], -combined["net_resident_income"],
            -outside["principal"],
        ), "principal closure"),
        "allocation_identity": _sum((
            combined["resources_after_operations_and_principal"],
            -combined["cash_resources_after_principal"], combined["essential_final_operations"],
        ), "allocation closure"),
    }
    return {
        "worlds": rows, "combined": combined, "transfers": flows,
        "net_receipts_from_outside_system": outside,
        "closure_residuals": checks,
    }


def debt_transfers(
    debtor: str,
    creditor: str | None,
    interest: float,
    principal_repaid: float,
    creditor_share: float = 1.0,
) -> list[dict[str, object]]:
    """Assign annual debt service to one creditor and domestic residual owners.

    The given creditor owns ``creditor_share``; the remaining share belongs to
    residents of the debtor world. Set creditor equal to debtor for wholly
    domestic debt. For staged loans, pass finance's first-mature-year interest
    and principal fields here; do not re-amortize their remaining stock.
    """
    share = _finite(creditor_share, "creditor_share", nonnegative=True)
    if share > 1:
        raise ValueError("creditor_share must not exceed one")
    interest = _finite(interest, "interest", nonnegative=True)
    principal = _finite(principal_repaid, "principal_repaid", nonnegative=True)
    return [
        {"payer": debtor, "recipient": recipient, "kind": kind, "amount": value * fraction}
        for recipient, fraction in ((creditor, share), (debtor, 1 - share))
        for kind, value in (("interest", interest), ("principal", principal))
        if value * fraction != 0
    ]


def mature_service_transfers(
    debtor: str,
    creditor: str | None,
    original_principal: float,
    annual_rate: float,
    repayment_years: int,
    payments_already_made: int = 0,
    creditor_share: float = 1.0,
) -> dict[str, object]:
    """Service an existing loan at full population without resetting its age.

    ``original_principal`` is the debt when repayments began, not its current
    remaining balance. The snapshot follows ``payments_already_made`` annual
    payments. The returned flows occur in the next year. Loans already paid
    off contribute zero mature debt service.
    """
    service = finance.rolling_tranches([{
        "principal": original_principal, "completion_year": 0,
        "annual_rate": annual_rate, "repayment_years": repayment_years,
    }], payments_already_made)
    service["transfers"] = debt_transfers(
        debtor, creditor, service["interest_first_mature_year"],
        service["principal_repaid_first_mature_year"], creditor_share,
    )
    return service
