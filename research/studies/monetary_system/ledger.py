"""A bounded financial transaction ledger, not a macroeconomic simulator.

Each claim names its issuer, holder, instrument and denomination. An edge is
simultaneously the holder's asset and issuer's liability, preventing missing
financial counterparties. Net worth is carried forward independently through
explicit transaction or valuation entries, then reconciled to every balance
sheet after each posting. Bank capital here is assets less debt liabilities;
shares and market capitalisation are outside this small reference.

U only values the balance sheets. L and S are hypothetical contract currencies.
No monetary regime, prices, production response or capital requirement is adopted.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping


class LiquidityError(ValueError):
    """An attempted transaction exceeds a specified liquid claim."""


def finite(value: float, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("A finite numerical quantity is required")
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError("Quantity must be finite and, where required, nonnegative")
    return value


@dataclass(frozen=True, order=True)
class Claim:
    issuer: str
    holder: str
    instrument: str
    currency: str


def claim(issuer: str, holder: str, instrument: str, currency: str = "L") -> Claim:
    return Claim(issuer, holder, instrument, currency)


class Ledger:
    """Explicit claims and real assets with atomic, reconciled postings.

    Real assets use fixed U valuations unless an event explicitly changes them.
    Opening net worth may be negative for a borrower. Later bank solvency policy
    is checked by the relevant operation, not by the accounting identity itself.
    """

    def __init__(self, actors: Mapping[str, str], prices_U: Mapping[str, float],
                 claims: Mapping[Claim, float], real_assets_U: Mapping[str, float]):
        self.actors = dict(actors)
        if not self.actors or any(not name or not place for name, place in self.actors.items()):
            raise ValueError("Named actors and territories required")
        self.prices_U = {k: finite(v, nonnegative=True) for k, v in prices_U.items()}
        if any(v == 0 for v in self.prices_U.values()):
            raise ValueError("Currency valuations must be positive")
        self.claims = {}
        for key, value in claims.items():
            self._validate_claim(key)
            self.claims[key] = finite(value, nonnegative=True)
        if set(real_assets_U) - self.actors.keys():
            raise ValueError("Unknown real-asset owner")
        self.real_assets_U = {a: finite(real_assets_U.get(a, 0), nonnegative=True)
                              for a in self.actors}
        self.carried_net_worth_U = self.net_worth()
        self.journal: list[dict] = []

    def _validate_claim(self, key: Claim):
        if key.issuer not in self.actors or key.holder not in self.actors:
            raise ValueError("Unknown claim counterparty")
        if key.issuer == key.holder or not key.instrument:
            raise ValueError("Claims require different counterparties and an instrument")
        if key.currency not in self.prices_U:
            raise ValueError("Unknown contract currency")

    def amount(self, key: Claim) -> float:
        self._validate_claim(key)
        return self.claims.get(key, 0.0)

    def net_worth(self) -> dict[str, float]:
        values = self.real_assets_U.copy()
        for key, amount in self.claims.items():
            value = amount * self.prices_U[key.currency]
            values[key.holder] += value
            values[key.issuer] -= value
        return values

    def post(self, label: str, changes: Mapping[Claim, float], *,
             real_changes_U: Mapping[str, float] | None = None,
             net_worth_changes_U: Mapping[str, float] | None = None,
             flow_class: str = "financial_transaction"):
        """Post all legs or none; income/loss/real-asset changes are explicit.

        Missing net-worth counterparts therefore fail, rather than becoming an
        unexplained balancing residual. No negative cash/claim balance is hidden.
        """
        real_changes_U = dict(real_changes_U or {})
        net_worth_changes_U = dict(net_worth_changes_U or {})
        if (set(real_changes_U) | set(net_worth_changes_U)) - self.actors.keys():
            raise ValueError("Unknown actor in flow account")
        next_claims = self.claims.copy()
        for key, delta in changes.items():
            self._validate_claim(key)
            updated = next_claims.get(key, 0) + finite(delta)
            if updated < -1e-9:
                raise LiquidityError(f"Insufficient {key.currency} {key.instrument}: {key.issuer}->{key.holder}")
            next_claims[key] = max(0.0, updated)
        next_real = self.real_assets_U.copy()
        for actor, delta in real_changes_U.items():
            next_real[actor] += finite(delta)
            if next_real[actor] < -1e-9:
                raise ValueError("Real assets cannot be overdrawn")
        expected = self.carried_net_worth_U.copy()
        for actor, delta in net_worth_changes_U.items():
            expected[actor] += finite(delta)
        observed = next_real.copy()
        for key, amount in next_claims.items():
            value = amount * self.prices_U[key.currency]
            observed[key.holder] += value
            observed[key.issuer] -= value
        residual = max(abs(observed[a] - expected[a]) for a in self.actors)
        if residual > 1e-8:
            raise ValueError(f"Unexplained change in net worth: {residual:g} U")
        self.claims, self.real_assets_U = next_claims, next_real
        self.carried_net_worth_U = expected
        self.journal.append({"event": label, "flow_class": flow_class, "claim_changes": [
            {**vars(key), "amount": delta} for key, delta in sorted(changes.items())],
            "real_asset_changes_U": real_changes_U,
            "net_worth_changes_U": net_worth_changes_U,
            "max_balance_sheet_residual_U": residual})

    def revalue(self, currency: str, price_U: float):
        """Exogenous valuation only; no assumed price/output/FX-market response."""
        price = finite(price_U, nonnegative=True)
        if currency not in self.prices_U or price <= 0:
            raise ValueError("Known currency and positive valuation required")
        change = price - self.prices_U[currency]
        valuation = {a: 0.0 for a in self.actors}
        for key, amount in self.claims.items():
            if key.currency == currency:
                valuation[key.holder] += amount * change
                valuation[key.issuer] -= amount * change
        self.prices_U[currency] = price
        self.carried_net_worth_U = {a: v + valuation[a]
                                    for a, v in self.carried_net_worth_U.items()}
        residual = max(abs(self.net_worth()[a] - self.carried_net_worth_U[a])
                       for a in self.actors)
        if residual > 1e-8:
            raise ValueError("Currency valuation failed reconciliation")
        self.journal.append({"event": f"Revalue {currency}", "flow_class": "revaluation", "valuation_changes_U": valuation,
                             "max_balance_sheet_residual_U": residual})

    def snapshot(self) -> dict:
        balances = {a: {"financial_assets_U": 0.0, "debt_liabilities_U": 0.0,
                        "real_assets_U": self.real_assets_U[a]} for a in self.actors}
        nfa = {territory: 0.0 for territory in set(self.actors.values())}
        for key, amount in self.claims.items():
            value = amount * self.prices_U[key.currency]
            balances[key.holder]["financial_assets_U"] += value
            balances[key.issuer]["debt_liabilities_U"] += value
            ih, hh = self.actors[key.issuer], self.actors[key.holder]
            if ih != hh:
                nfa[hh] += value
                nfa[ih] -= value
        for actor, row in balances.items():
            row["net_worth_U"] = self.carried_net_worth_U[actor]
        total_nw = sum(self.carried_net_worth_U.values())
        return {"prices_U": self.prices_U.copy(), "balance_sheets": balances,
                "claims": [{**vars(k), "amount": v} for k, v in sorted(self.claims.items()) if v],
                "net_foreign_assets_U": dict(sorted(nfa.items())),
                "aggregate_net_worth_U": total_nw,
                "aggregate_real_assets_U": sum(self.real_assets_U.values()),
                "global_financial_residual_U": sum(nfa.values()),
                "net_worth_minus_real_assets_U": total_nw - sum(self.real_assets_U.values()),
                "posting_count": len(self.journal)}


def bank_loan(book: Ledger, bank: str, borrower: str, amount: float):
    amount = finite(amount, nonnegative=True)
    book.post("Commercial bank loan creates matching deposit", {
        claim(borrower, bank, "loan"): amount,
        claim(bank, borrower, "deposit"): amount})


def loan_repayment(book: Ledger, bank: str, borrower: str, amount: float):
    amount = finite(amount, nonnegative=True)
    book.post("Principal repayment cancels loan and deposit", {
        claim(borrower, bank, "loan"): -amount,
        claim(bank, borrower, "deposit"): -amount})


def domestic_payment(book: Ledger, payer: str, payer_bank: str, payee: str,
                     payee_bank: str, amount: float):
    """Current transfer between actors, or a financial transfer between own accounts.

    The reference uses transfers to isolate monetary mechanics. Output is not
    inferred from the payment. Within-bank payments need no reserve transfer.
    """
    amount = finite(amount, nonnegative=True)
    changes = {claim(payer_bank, payer, "deposit"): -amount,
               claim(payee_bank, payee, "deposit"): amount}
    if payer == payee and payer_bank == payee_bank:
        raise ValueError("Payment needs distinct accounts")
    if payer_bank != payee_bank:
        changes.update({claim("moon_cb", payer_bank, "reserve"): -amount,
                        claim("moon_cb", payee_bank, "reserve"): amount})
    nw = {payer: -amount * book.prices_U["L"], payee: amount * book.prices_U["L"]}
    if payer == payee:
        nw = {}
    book.post("Domestic transfer and reserve settlement", changes,
              net_worth_changes_U=nw,
              flow_class="financial_transaction" if payer == payee else "current_transfer")


def central_bank_loan(book: Ledger, bank: str, amount: float, collateral_limit_L: float,
                      *, require_solvent: bool = True):
    amount = finite(amount, nonnegative=True)
    limit = finite(collateral_limit_L, nonnegative=True)
    if require_solvent and book.net_worth()[bank] < 0:
        raise ValueError("Insolvent bank does not qualify for the reference liquidity facility")
    debt = claim(bank, "moon_cb", "central_bank_loan")
    if book.amount(debt) + amount > limit:
        raise ValueError("Credit exceeds assumed eligible collateral after haircut")
    book.post("Own-currency central-bank lending", {
        debt: amount, claim("moon_cb", bank, "reserve"): amount})


def cancel_loan(book: Ledger, bank: str, borrower: str, amount: float):
    """Agreed debt cancellation: deposits survive; legal debt is released."""
    amount = finite(amount, nonnegative=True)
    value = amount * book.prices_U["L"]
    book.post("Agreed cancellation of borrower debt", {claim(borrower, bank, "loan"): -amount},
              net_worth_changes_U={bank: -value, borrower: value},
              flow_class="agreed_debt_forgiveness_capital_transfer")


def issue_treasury_bond(book: Ledger, buyer: str, bank: str, amount: float):
    amount = finite(amount, nonnegative=True)
    changes = {claim("moon_treasury", buyer, "bond"): amount,
               claim("moon_cb", bank, "reserve"): -amount,
               claim("moon_cb", "moon_treasury", "treasury_account"): amount}
    if buyer != bank:
        changes[claim(bank, buyer, "deposit")] = -amount
    book.post("Treasury bond issue", changes)


def treasury_transfer(book: Ledger, recipient: str, bank: str, amount: float):
    """Unconditional current transfer; no purchased output or investment grant."""
    amount = finite(amount, nonnegative=True)
    value = amount * book.prices_U["L"]
    book.post("Treasury transfer from its funded account", {
        claim("moon_cb", "moon_treasury", "treasury_account"): -amount,
        claim("moon_cb", bank, "reserve"): amount,
        claim(bank, recipient, "deposit"): amount},
        net_worth_changes_U={"moon_treasury": -value, recipient: value},
        flow_class="current_transfer")


def central_bank_bond_purchase(book: Ledger, seller: str, bank: str, amount: float):
    amount = finite(amount, nonnegative=True)
    changes = {claim("moon_treasury", seller, "bond"): -amount,
               claim("moon_treasury", "moon_cb", "bond"): amount,
               claim("moon_cb", bank, "reserve"): amount}
    if seller != bank:
        changes[claim(bank, seller, "deposit")] = amount
    book.post("Central bank secondary-market bond purchase", changes)


def buy_foreign_reserves(book: Ledger, bank: str, amount_S: float):
    """Bank purchases S from the CB against L reserves; no outside S is created."""
    amount = finite(amount_S, nonnegative=True)
    cost_L = amount * book.prices_U["S"] / book.prices_U["L"]
    book.post("Bank buys existing foreign settlement assets from lunar CB", {
        claim("settlement_issuer", "moon_cb", "settlement", "S"): -amount,
        claim("settlement_issuer", bank, "settlement", "S"): amount,
        claim("moon_cb", bank, "reserve"): -cost_L})


def sell_foreign_reserves(book: Ledger, bank: str, amount_S: float):
    amount = finite(amount_S, nonnegative=True)
    proceeds_L = amount * book.prices_U["S"] / book.prices_U["L"]
    book.post("Lunar CB buys foreign settlement assets from bank", {
        claim("settlement_issuer", bank, "settlement", "S"): -amount,
        claim("settlement_issuer", "moon_cb", "settlement", "S"): amount,
        claim("moon_cb", bank, "reserve"): proceeds_L})


def import_asset(book: Ledger, buyer: str, bank: str, seller: str, amount_S: float):
    """Import movable produced equipment, changing location and economic owner.

    Equipment is delivered from the seller's territory to the buyer's territory.
    Existing assets move at the same U valuation; no new production is inferred.
    Land, financial claims and purchases left in the source territory are excluded.
    """
    if book.actors[buyer] == book.actors[seller]:
        raise ValueError("Imports require different buyer and seller territories")
    amount = finite(amount_S, nonnegative=True)
    value = amount * book.prices_U["S"]
    book.post("Import existing asset; settle in S", {
        claim(bank, buyer, "deposit"): -value / book.prices_U["L"],
        claim("settlement_issuer", bank, "settlement", "S"): -amount,
        claim("settlement_issuer", seller, "settlement", "S"): amount},
        real_changes_U={buyer: value, seller: -value},
        flow_class="goods_trade_existing_produced_equipment")


def export_service_receipt(book: Ledger, seller: str, bank: str, buyer: str, amount_S: float):
    """Receipt for a consumed service; no claim about its costs or value added."""
    if book.actors[seller] == book.actors[buyer]:
        raise ValueError("Exports require different buyer and seller territories")
    amount = finite(amount_S, nonnegative=True)
    value = amount * book.prices_U["S"]
    book.post("External customer pays lunar service provider", {
        claim("settlement_issuer", buyer, "settlement", "S"): -amount,
        claim("settlement_issuer", bank, "settlement", "S"): amount,
        claim(bank, seller, "deposit"): value / book.prices_U["L"]},
        net_worth_changes_U={seller: value, buyer: -value},
        flow_class="current_service_trade")


def external_treasury_bond(book: Ledger, buyer: str, amount: float, currency: str):
    """Foreign investor accepts the specified denomination; settles in S."""
    if book.actors[buyer] == book.actors["moon_treasury"]:
        raise ValueError("External debt requires a nonresident creditor")
    amount = finite(amount, nonnegative=True)
    if currency not in book.prices_U:
        raise ValueError("Unknown debt denomination")
    settlement = amount * book.prices_U[currency] / book.prices_U["S"]
    book.post("External investor purchases newly issued Treasury debt", {
        claim("moon_treasury", buyer, "bond", currency): amount,
        claim("settlement_issuer", buyer, "settlement", "S"): -settlement,
        claim("settlement_issuer", "moon_treasury", "settlement", "S"): settlement})


def treasury_sells_foreign_receipt(book: Ledger, amount_S: float):
    amount = finite(amount_S, nonnegative=True)
    proceeds = amount * book.prices_U["S"] / book.prices_U["L"]
    book.post("Treasury exchanges its S proceeds with lunar CB", {
        claim("settlement_issuer", "moon_treasury", "settlement", "S"): -amount,
        claim("settlement_issuer", "moon_cb", "settlement", "S"): amount,
        claim("moon_cb", "moon_treasury", "treasury_account"): proceeds})


def destroy_produced_capital(book: Ledger, owner: str, amount_U: float):
    """Assumed physical damage; an other-volume change with no financial payment."""
    amount = finite(amount_U, nonnegative=True)
    book.post("Assumed destruction of produced capital", {},
              real_changes_U={owner: -amount}, net_worth_changes_U={owner: -amount},
              flow_class="other_volume_change")
