# Money, banking and interworld finance

31 October 2026. **The Open Moon needs a monetary analysis at the scale of an economy: credit creation, settlement, public finance, external payments and the institutions that absorb losses.** This study compares arrangements for early settlement and a mature Moon integrated with Earth, the array and the wider Solar System. The [sports and media account](../world_economy/sports_finance.md) remains a separate preliminary business scenario.

The purpose is to establish which monetary mechanisms could assist development and what each would require. The reader is assumed to know the three-world production comparison, but not its banking arrangements. Institutional descriptions draw on primary sources; lunar arrangements are conditional proposals. The reference ledger tests specified transactions and balance-sheet identities. It neither estimates a future money supply nor selects a currency regime.

## Current findings

**Currency use, access to central-bank facilities and fiscal risk sharing are separate decisions.** An early settlement can use a supplier's currency without becoming entitled to its lender of last resort. A monetary union could supply negotiated membership and support, but must decide representation, supervision, resolution and loss sharing. Three worlds need not imply three currencies or three sovereign governments.

**Monetary autonomy can have a material production benefit when finance is the constraint.** Credit or public spending can help available lunar facilities, inputs and useful work remain active through a demand or payment crisis. Its effect depends on what expenditure buys. A missing external input requires an accepting supplier and delivery, and an obligation in outside currency requires access to that currency. This distinction becomes more consequential as domestic supply broadens over the settlement's development.

**A bank loan, settlement liquidity and solvency capital do different work.** In the reference example a new loan of 20 L creates 20 L of deposits. An interbank payment then exposes a 12 L reserve shortage at Bank A, which eligible central-bank lending can bridge. Bank B initially holds 12 L: private interbank funding could redistribute those reserves, but is excluded from this chosen facility example. After assumed physical damage to the borrower’s assets, an agreed debt cancellation reduces the bank's capital from 10 U to −5 U; another liquidity loan leaves it at −5 U. Recapitalization or resolution would require a separate loss-allocation decision. The illustrative commercial-bank access rule is distinct from the accounting and operational position of a currency-issuing central bank.

**Foreign receipts and reserve backing require ownership accounts.** An export initially pays its supplier; it does not automatically replenish official reserves. A currency board's base-money backing also differs from backing every customer deposit. In the external example another 50 L of central-bank liquidity cannot supply a missing 3 S of external settlement assets. A subsequent export and actual reserve purchase can. This is a settlement constraint in a fixed-price example, not a predicted exchange-rate path.

**Construction finance and monetary policy interact without becoming the same account.** Fiscal spending, bond issue, central-bank purchases, foreign financing and valuation changes have distinct effects. The inherited construction bill and mature product estimates remain unchanged. A prospective benefit from lower financing risk or better stabilization would need a behavioral model before being added to income or deducted from debt.

**No automatic progression toward a floating lunar currency is established.** External-currency use, negotiated union, a board or peg, and an independent float remain alternatives whose merits depend on invoicing, domestic supply, shocks, mobility, fiscal capacity and political arrangements. A mature Moon could still prefer a shared currency. A distinct local unit could also become useful before maturity if its institutions and economic conditions support it.

## Analysis and reproducible account

| Component | What it establishes | Limit |
|---|---|---|
| [Institutions](institutions.md) | Four monetary architectures; banking, fiscal roles, stabilization, resolution, seigniorage and payment continuity | Conditional institutional comparison; no adopted constitution, law or monetary mandate |
| [External accounts and construction](external_accounts.md) | Current and financial accounts, currency mismatch, public purchases, reserve requirements and connection to inherited debt scenarios | No equilibrium trade matrix, inflation forecast or estimate of monetary-policy gains |
| [Reference ledger](ledger.py) and [runner](run.py) | Paired financial claims, explicit transaction effects, interbank settlement, credit loss, fiscal operations and external valuation | Arbitrary normalized balances; constant valuations except the stated currency revaluation |
| [Inputs](inputs.json) and [result](results/reference.json) | Reproducible opening balance sheets, event sequences and residual checks | No mapping of illustrative monetary stocks onto mature world output |
| [Sources](sources.json) and [checks](checks.json) | Targeted reading scopes, local baseline identities, verification and unresolved work | Targeted institutional research, not a completed empirical calibration or full literature review |

The model identifies two lunar banks, a monetary authority, a Treasury, firms, households, Earth and array counterparties, and an outside S issuer. L and S are hypothetical currencies; U only values their claims. Every scenario starts afresh from the same opening balances. The accounts include a matching outside issuer and real assets so that financial claims do not manufacture aggregate net worth. Bank net worth is unallocated capital: shareholder securities and distributions are outside the model.

Four examples expose different mechanisms:

1. A loan creates deposits; payment requires reserves; a credit loss requires separate treatment from liquidity.
2. Treasury issuance and transfers are traced through bank and nonbank bond buyers, followed by central-bank secondary purchases from each.
3. An importer exhausts external settlement assets; an array customer pays for a lunar service; a subsequent official reserve purchase permits further imports. Territorial current-account balances reconcile with external financial positions.
4. The same outside creditor holds lunar- or external-currency debt. A specified depreciation changes either the creditor's valuation or the borrower's local repayment requirement, without inventing an export receipt.

The code validates transaction feasibility before changing balances, records expected net-worth effects and tests failed operations as well as successful settlement. The loss sequence first destroys 15 U of real assets, reducing aggregate wealth from 365 U to 350 U. Agreed cancellation then reduces the borrower's obligation by 15 L and transfers that financial loss to its bank without a second loss of aggregate wealth. It does not attempt bank provisioning or a unilateral accounting write-down while a legal claim remains outstanding. Collateral ceilings, access rules and exchange prices are illustrative assumptions. There is no calibrated bank regulation, interest accrual, endogenous exchange rate, price formation or behavioral credit demand.

Reproduce with:

```sh
python -m research.studies.monetary_system.run
python -m pytest research/studies/monetary_system/tests -q
```

The default runner writes under ignored `research/runs/monetary_system/`. Use `--output research/studies/monetary_system/results/reference.json` to regenerate the committed data product. Its schema, producer hashes, source identity, evidence statement and reading rule travel with the result. It is a standalone reference product; no domain imports a second domain's implementation.

## Evidence needed for a quantitative policy comparison

A monetary regime cannot be ranked solely from population, low gravity or the scale of construction. The next comparison needs staged domestic supply and import accounts; currencies of invoices, deposits and debt; creditor residence and maturities; fiscal revenue and guarantees; and possible joint shocks to the Moon, Earth and array. Those inputs would reveal where monetary flexibility preserves output, where a common currency reduces risk, and how much external liquidity remains necessary.

The political decisions also remain substantive. Banking support requires a bearer of losses; a monetary union requires representation and common rules; an independent issuer needs a mandate and a fiscal relationship. Automation can improve administration while leaving those choices intact. The current result is an institutional and accounting foundation for evaluating them, with the construction, commercial and monetary evidence kept at their respective scales.
