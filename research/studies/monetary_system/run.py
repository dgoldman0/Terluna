"""Reproduce bounded monetary transactions; no GWP/NWP or FX forecast."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import ledger as m

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_context_identities(source_record: dict, root: Path = ROOT):
    """Preserve exact baseline context; these files do not supply model numbers."""
    for row in source_record.get("local_context", []):
        path = root / row["path"]
        if digest(path) != row["sha256"]:
            raise ValueError(f"Local context changed; review its source identity: {row['path']}")


def load_inputs(path: Path = HERE / "inputs.json") -> dict:
    spec = json.loads(Path(path).read_text())
    if spec["schema"] != "terluna.research.monetary-system-inputs/1":
        raise ValueError("Unexpected monetary input schema")
    return spec


def opening(spec: dict) -> m.Ledger:
    claims = {}
    for issuer, holder, instrument, currency, amount in spec["opening_claims"]:
        key = m.claim(issuer, holder, instrument, currency)
        if key in claims:
            raise ValueError("Duplicate opening claim")
        claims[key] = amount
    return m.Ledger(spec["actors"], spec["prices_U"], claims,
                    spec["opening_real_assets_U"])


def monetary_totals(book: m.Ledger) -> dict:
    """Separate resident deposits, domestic reserves, foreign settlement assets."""
    def total(instrument, currency, *, holders=None):
        return sum(v for k, v in book.claims.items()
                   if k.instrument == instrument and k.currency == currency
                   and (holders is None or k.holder in holders))
    moon = {a for a, place in book.actors.items() if place == "Moon"}
    return {"resident_L_deposits": total("deposit", "L", holders=moon),
            "bank_L_reserves": total("reserve", "L"),
            "moon_total_S_settlement_assets": total("settlement", "S", holders=moon),
            "cb_S_settlement_assets": book.amount(m.claim("settlement_issuer", "moon_cb", "settlement", "S")),
            "bank_a_S_settlement_assets": book.amount(m.claim("settlement_issuer", "moon_bank_a", "settlement", "S")),
            "bank_a_capital_U": book.net_worth()["moon_bank_a"],
            "treasury_bonds_L": total("bond", "L")}


def capture(book: m.Ledger) -> dict:
    return {"monetary_totals": monetary_totals(book), **book.snapshot()}


def domestic_scenario(spec: dict) -> dict:
    p = spec["domestic_liquidity"]
    book = opening(spec)
    stages = {"opening": capture(book)}
    m.bank_loan(book, "moon_bank_a", "moon_firms", p["new_loan_L"])
    stages["loan_created"] = capture(book)
    reserves = book.amount(m.claim("moon_cb", "moon_bank_a", "reserve"))
    gap = max(0, p["interbank_payment_L"] - reserves)
    try:
        m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_households",
                           "moon_bank_b", p["interbank_payment_L"])
    except m.LiquidityError:
        blocked = True
    else:
        blocked = False
    if not blocked:
        raise ValueError("Reference stress must expose an initial settlement shortage")
    m.central_bank_loan(book, "moon_bank_a", gap,
                        p["eligible_collateral_after_haircut_L"])
    m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_households",
                       "moon_bank_b", p["interbank_payment_L"])
    stages["after_liquidity_and_settlement"] = capture(book)
    m.destroy_produced_capital(book, "moon_firms", p["damaged_produced_capital_U"])
    stages["after_physical_capital_damage"] = capture(book)
    m.cancel_loan(book, "moon_bank_a", "moon_firms", p["agreed_debt_cancellation_L"])
    stages["after_agreed_debt_cancellation"] = capture(book)
    try:
        m.central_bank_loan(book, "moon_bank_a", p["counterfactual_additional_liquidity_L"],
                            p["eligible_collateral_after_haircut_L"])
    except ValueError:
        ineligible = True
    else:
        ineligible = False
    if not ineligible:
        raise ValueError("Reference loss must expose insolvency under the chosen access rule")
    m.central_bank_loan(book, "moon_bank_a", p["counterfactual_additional_liquidity_L"],
                        p["eligible_collateral_after_haircut_L"], require_solvent=False)
    stages["counterfactual_liquidity_does_not_restore_capital"] = capture(book)
    return {"settlement_shortfall_L": gap, "initial_payment_blocked": blocked,
            "ordinary_liquidity_request_rejected_after_loss": ineligible,
            "stages": stages, "journal": book.journal,
            "boundary": "Current transfers isolate settlement. No production is inferred. Assumed capital damage reduces real wealth before debt forgiveness. Agreed debt cancellation is a capital transfer that legally releases the borrower; it is not an impairment provision. Last credit is an expressly ineligible accounting counterfactual."}


def fiscal_scenarios(spec: dict) -> dict:
    p = spec["fiscal_transactions"]
    outcomes = {}
    for name, buyer in [("nonbank_buyer", "moon_households"),
                        ("commercial_bank_buyer", "moon_bank_b")]:
        book = opening(spec)
        stages = {"opening": capture(book)}
        m.issue_treasury_bond(book, buyer, "moon_bank_b", p["bond_issue_L"])
        stages["bond_issued"] = capture(book)
        m.treasury_transfer(book, "moon_firms", "moon_bank_a", p["transfer_L"])
        stages["treasury_spent"] = capture(book)
        m.central_bank_bond_purchase(book, buyer, "moon_bank_b", p["secondary_purchase_L"])
        stages["central_bank_secondary_purchase"] = capture(book)
        outcomes[name] = {"stages": stages, "journal": book.journal}
    return {"alternatives": outcomes,
            "boundary": "The unconditional fiscal current transfer is a financial example with no assumed production; it is not an investment grant. A procurement scenario would also require real delivery and cost/value-added accounts. CB secondary purchases transfer the existing bond; they supply no new proceeds to the Treasury. Public consolidation changes the form, not the outstanding amount, of public liabilities to the private sector in the purchase step."}


def external_debt_scenarios(spec: dict) -> dict:
    """Same creditor, different denomination; residence and currency stay distinct."""
    cases = {}
    amount = spec["fiscal_transactions"]["bond_issue_L"]
    for currency in ("L", "S"):
        book = opening(spec)
        before = capture(book)
        m.external_treasury_bond(book, "earth_counterparty", amount, currency)
        settlement = amount * book.prices_U[currency] / book.prices_U["S"]
        m.treasury_sells_foreign_receipt(book, settlement)
        m.treasury_transfer(book, "moon_firms", "moon_bank_a",
                            settlement * book.prices_U["S"] / book.prices_U["L"])
        borrowed = capture(book)
        book.revalue("L", spec["foreign_liquidity"]["depreciated_U_per_L"])
        after = capture(book)
        change = {w: after["net_foreign_assets_U"][w] - borrowed["net_foreign_assets_U"][w]
                  for w in before["net_foreign_assets_U"]}
        cases[currency] = {"stages": {"opening": before, "funded_and_spent": borrowed,
                                      "after_L_depreciation": after},
                           "external_debt_face": amount, "denomination": currency,
                           "debt_value_U_after": amount * book.prices_U[currency],
                           "debt_value_L_after": amount * book.prices_U[currency] / book.prices_U["L"],
                           "current_account_U": {w: 0.0 for w in change},
                           "NFA_valuation_U": change,
                           "journal": book.journal}
    return {"alternatives": cases,
            "boundary": "External purchase is assumed, not predicted. Borrowing creates a financial inflow and external liability; it is not export income. Funds retained in the lunar CB remain an external asset. Local fiscal transfers create no current-account entry. Contract face is unchanged by imposed L depreciation. No foreign creditor loss/repricing response, hedging or inflation response is modeled."}


def foreign_scenario(spec: dict) -> dict:
    p = spec["foreign_liquidity"]
    book = opening(spec)
    stages = {"opening": capture(book)}
    m.bank_loan(book, "moon_bank_a", "moon_firms", p["working_capital_loan_L"])
    m.central_bank_loan(book, "moon_bank_a", p["initial_cb_loan_L"],
                        p["eligible_collateral_after_haircut_L"])
    stages["domestic_liquidity_available"] = capture(book)
    reserve_key = m.claim("settlement_issuer", "moon_cb", "settlement", "S")
    foreign_available = book.amount(reserve_key)
    initial_gap = max(0, p["requested_import_S"] - foreign_available)
    try:
        m.buy_foreign_reserves(book, "moon_bank_a", p["requested_import_S"])
    except m.LiquidityError:
        fx_blocked = True
    else:
        fx_blocked = False
    if not fx_blocked:
        raise ValueError("Reference must expose insufficient external settlement assets")
    first_paid = min(p["requested_import_S"], foreign_available)
    m.buy_foreign_reserves(book, "moon_bank_a", first_paid)
    m.import_asset(book, "moon_firms", "moon_bank_a", "earth_counterparty", first_paid)
    stages["available_S_spent_on_import"] = capture(book)
    m.central_bank_loan(book, "moon_bank_a", p["additional_local_liquidity_L"],
                        p["eligible_collateral_after_haircut_L"])
    stages["additional_L_liquidity"] = capture(book)
    gap_after_local_credit = max(0, initial_gap - book.amount(reserve_key))
    m.export_service_receipt(book, "moon_firms", "moon_bank_a", "array_counterparty",
                             p["array_service_receipt_S"])
    stages["external_receipt_owned_by_bank"] = capture(book)
    m.sell_foreign_reserves(book, "moon_bank_a", p["array_service_receipt_S"])
    stages["cb_purchases_external_receipt"] = capture(book)
    remaining = min(initial_gap, book.amount(reserve_key))
    m.buy_foreign_reserves(book, "moon_bank_a", remaining)
    m.import_asset(book, "moon_firms", "moon_bank_a", "earth_counterparty", remaining)
    stages["remaining_import_settled"] = capture(book)
    current_account = {"Moon": p["array_service_receipt_S"] - first_paid - remaining,
                       "Earth": first_paid + remaining,
                       "Array": -p["array_service_receipt_S"]}
    nfa0 = stages["opening"]["net_foreign_assets_U"]
    nfa1 = stages["remaining_import_settled"]["net_foreign_assets_U"]
    ca_U = {w: x * spec["prices_U"]["S"] for w, x in current_account.items()}
    residual = {w: nfa1[w] - nfa0[w] - ca_U[w] for w in ca_U}
    # Valuation after all real transactions; the U values of real assets are held fixed.
    book.revalue("L", p["depreciated_U_per_L"])
    stages["exogenous_L_depreciation"] = capture(book)
    nfa2 = stages["exogenous_L_depreciation"]["net_foreign_assets_U"]
    valuation = {w: nfa2[w] - nfa1[w] for w in nfa1}
    return {"initial_foreign_funding_gap_S": initial_gap,
            "gap_after_additional_local_credit_S": gap_after_local_credit,
            "import_paid_S": first_paid + remaining,
            "remaining_unpaid_import_S": initial_gap - remaining,
            "stages": stages, "journal": book.journal,
            "current_account_U": ca_U,
            "NFA_transaction_reconciliation_residual_U": residual,
            "NFA_valuation_U": valuation,
            "boundary": "Earth delivers existing movable produced equipment to the Moon and transfers economic ownership; array consumes a lunar service. Costs, value added, freight limits and prices are not modeled. No foreign debt/interest/current transfers/capital transfers occur in this sequence. External receipt first belongs to the bank, not the CB. Depreciation is imposed after settlement, with fixed U real-asset values and no behavioural feedback; its NFA effect is zero here because cross-border claims are all S-denominated."}


def results(input_path: Path = HERE / "inputs.json") -> dict:
    spec = load_inputs(input_path)
    source = HERE / "sources.json"
    if not source.is_file():
        raise FileNotFoundError("The monetary source record must exist before producing results")
    verify_context_identities(json.loads(source.read_text()))
    model_hash = hashlib.sha256(b"".join(p.read_bytes() for p in
                                      (HERE / "ledger.py", HERE / "run.py"))).hexdigest()
    return {"schema": "terluna.research.monetary-system/1",
            "producer": {"module": "research.studies.monetary_system.run",
                         "model_sha256": model_hash,
                         "inputs_sha256": digest(Path(input_path)),
                         "sources_sha256": digest(source)},
            "evidence": spec["evidence"], "reading_rule": spec["units"],
            "assumptions": spec,
            "domestic_liquidity_and_solvency": domestic_scenario(spec),
            "treasury_and_central_bank": fiscal_scenarios(spec),
            "external_debt_denomination": external_debt_scenarios(spec),
            "foreign_settlement": foreign_scenario(spec),
            "boundaries": ["No GWP/NWP baseline or existing economic product is changed.",
                           "Accounting consistency is not evidence of sustainable fiscal or monetary policy.",
                           "No exchange-rate equilibrium, inflation, credit demand, interest-rate response, default probability, regulatory capital rule or liquidity-premium forecast is made.",
                           "Money creation, real resources, liquidity, solvency and financing capacity require separate evidence."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "research/runs/monetary_system/reference.json")
    args = parser.parse_args()
    product = results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
