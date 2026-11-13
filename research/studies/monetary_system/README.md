# Money, banking and interworld finance

The current [U trade and currency transition plan](trade_transition.md) concerns the planned off-Earth economic zone before the Open Moon, planetary array and subsequent outer Solar System buildout. It develops a substantial transition from USD toward either a single U currency or a basket of U-related currencies, with separately funded partial reserves. Its transaction, basket and reserve calculations preserve the older institutional reference below. The [nine-portfolio pressure study](pressure.md) now couples actual distributions and reserve contributions to conditional monetary-stock transmission. It distinguishes domestic price levels, outside balance demand, issuance shortfalls/excesses, import costs and withdrawal stress. It does not forecast floating exchange rates or calibrate U inflation.

31 October 2026. **Any monetary arrangement for the Open Moon has to issue the claims that circulate, settle payments, extend credit, bear losses and risk, and take its collective decisions through whoever holds the decision rights.** This study compares one family of arrangements, the currency areas, banks, fiscal authorities and central banks of present economies, for early settlement and for a mature Moon integrated with Earth, the array and the wider Solar System. The [sports and media account](../world_economy/sports_finance.md) remains a separate preliminary business scenario.

It sets out which mechanisms could help development and what each would require, for a reader who knows the three-world production comparison. The descriptions of present institutions draw on primary sources, and the lunar arrangements are proposals for comparison; the choice among them belongs to the communities' collective decisions ([register](../../decisions.md#automation-the-network-and-food-12-october-2026)). The reference ledger tests specified transactions and balance-sheet identities at normalized sizes, leaving the money supply and the currency arrangement to those decisions.

## Current findings

**Currency use, access to central-bank facilities and fiscal risk sharing are separate decisions.** An early settlement can use a supplier's currency without becoming entitled to its lender of last resort. A monetary union could supply negotiated membership and support, but must decide representation, supervision, resolution and loss sharing. Three worlds need not imply three currencies or three separate issuers.

**Monetary autonomy can have a material production benefit when finance is the constraint.** Credit or public spending can help available lunar facilities, inputs and useful work remain active through a demand or payment crisis. Its effect depends on what expenditure buys. A missing external input requires an accepting supplier and delivery, and an obligation in outside currency requires access to that currency. This distinction becomes more consequential as domestic supply broadens over the settlement's development.

**A bank loan, settlement liquidity and solvency capital do different work.** In the reference example a new loan of 20 L creates 20 L of deposits. An interbank payment then exposes a 12 L reserve shortage at Bank A, which eligible central-bank lending can bridge. Bank B initially holds 12 L: private interbank funding could redistribute those reserves, but is excluded from this chosen facility example. After assumed physical damage to the borrower’s assets, an agreed debt cancellation reduces the bank's capital from 10 U to −5 U; another liquidity loan leaves it at −5 U. Recapitalization or resolution would require a separate loss-allocation decision. The illustrative commercial-bank access rule is distinct from the accounting and operational position of a currency-issuing central bank.

**Foreign receipts and reserve backing require ownership accounts.** An export initially pays its supplier, and official reserves rise only when the authority buys the receipts. A currency board's base-money backing also differs from backing every customer deposit. In the external example another 50 L of central-bank liquidity cannot supply a missing 3 S of external settlement assets. A subsequent export and actual reserve purchase can. This is a settlement constraint in a fixed-price example.

**Construction finance and monetary policy interact without becoming the same account.** Fiscal spending, bond issue, central-bank purchases, foreign financing and valuation changes have distinct effects. The inherited construction bill and mature product estimates remain unchanged. A prospective benefit from lower financing risk or better stabilization would need a behavioral model before being added to income or deducted from debt.

**Each currency arrangement can last.** External-currency use, a negotiated union, a board or peg, and an independent float are alternatives whose merits depend on invoicing, domestic supply, shocks, mobility, collective finance and the arrangements the communities choose. A mature Moon could keep a shared currency, and a distinct local unit could serve before maturity where its institutions and economic conditions support it.

## Analysis and reproducible account

| Component | What it establishes | Limit |
|---|---|---|
| [Institutions](institutions.md) | Four currency arrangements; banking, fiscal roles, stabilization, resolution, seigniorage and payment continuity | Conditional institutional comparison; no adopted constitution, law or monetary mandate |
| [External accounts and construction](external_accounts.md) | Current and financial accounts, currency mismatch, public purchases, reserve requirements and connection to inherited debt scenarios | No equilibrium trade matrix, inflation forecast or estimate of monetary-policy gains |
| [Reference ledger](ledger.py) and [runner](run.py) | Paired financial claims, explicit transaction effects, interbank settlement, credit loss, fiscal operations and external valuation | Arbitrary normalized balances; constant valuations except the stated currency revaluation |
| [Inputs](inputs.json) and [result](results/reference.json) | Reproducible opening balance sheets, event sequences and residual checks | No mapping of illustrative monetary stocks onto mature world output |
| [Sources](sources.json) and [checks](checks.json) | Targeted reading scopes, local baseline identities, verification and unresolved work | Targeted institutional research; empirical calibration and a full literature review outstanding |
| [U trade transition](trade_transition.md), [calculation](trade_transition.py) and [results](results/trade_transition.json) | Single U versus constituent-currency basket, phased trade adoption, fixed-stock FX transactions, funded partial reserves and transaction-balance scenarios for nine portfolios | Proposed institutions and normalized demonstrations; no adopted constituents, issuer, peg or monetary forecast |

The ledger's entities are two lunar banks, a monetary authority, a Treasury standing for the collective budget, firms, households, Earth and array counterparties, and an outside S issuer. L and S are hypothetical currencies; U only values their claims. Every scenario starts afresh from the same opening balances. The accounts include a matching outside issuer and real assets, so financial claims add no aggregate net worth of their own. Bank net worth is unallocated capital; shareholder securities and distributions sit outside the model.

Four examples expose different mechanisms:

1. A loan creates deposits; payment requires reserves; a credit loss requires separate treatment from liquidity.
2. Treasury issuance and transfers are traced through bank and nonbank bond buyers, followed by central-bank secondary purchases from each.
3. An importer exhausts external settlement assets; an array customer pays for a lunar service; a subsequent official reserve purchase permits further imports. Territorial current-account balances reconcile with external financial positions.
4. The same outside creditor holds lunar- or external-currency debt. A specified depreciation changes either the creditor's valuation or the borrower's local repayment requirement, without inventing an export receipt.

The code validates transaction feasibility before changing balances, records expected net-worth effects and tests failed operations as well as successful settlement. The loss sequence first destroys 15 U of real assets, reducing aggregate wealth from 365 U to 350 U. Agreed cancellation then reduces the borrower's obligation by 15 L and transfers that financial loss to its bank without a second loss of aggregate wealth. Bank provisioning and unilateral write-downs while a claim stays outstanding are outside the ledger, as are calibrated bank regulation, interest accrual, endogenous exchange rates, price formation and behavioural credit demand; collateral ceilings, access rules and exchange prices are illustrative assumptions.

Reproduce with:

```sh
python -m research.studies.monetary_system.run
python -m pytest research/studies/monetary_system/tests -q
```

The default runner writes under ignored `research/runs/monetary_system/`. Use `--output research/studies/monetary_system/results/reference.json` to regenerate the committed data product. Its schema, producer hashes, source identity, evidence statement and reading rule travel with the result. It is a standalone reference product.

## Evidence needed for a quantitative policy comparison

Ranking monetary arrangements needs more than population, low gravity or the scale of construction. The next comparison needs staged domestic supply and import accounts; the currencies of invoices, deposits and debt; creditor residence and maturities; collective revenue and guarantees; and possible joint shocks to the Moon, Earth and array. Those inputs would show where monetary flexibility preserves output, where a common currency reduces risk, and how much external liquidity remains necessary.

The collective decisions are substantive too, and they belong to the communities' decision-making ([register](../../decisions.md#automation-the-network-and-food-12-october-2026)): banking support needs a bearer of losses, a monetary union representation and common rules, and an independent issuer a mandate and a relationship with the collective budget. Automation can improve the administration of each while those choices stay with the people who make them. The present result is an institutional and accounting foundation for evaluating them, with the construction, commercial and monetary evidence kept at their respective scales.

The [credit-supported financial comparison](../complex_returns/credit.md) applies
this monetary work to six investment-financing policies, including outside
suppliers accepting U term claims. It maintains separate external purchasing
power, loan-created deposits, restricted prepayments, loss capital and funded
payout deferral. Its continuous matched price-premium rule is a policy experiment,
with actual consumer inflation and floating exchange rates still uncalibrated.

The [U monetary-policy comparison](../complex_returns/policy.md)
separates modest local loan-rate relief, wider supplier acceptance and temporary
receivables support. In the strong year-100 propellant/logistics case, a six-month
payment delay requires 144.1 billion USD2041-equivalent of bridge funding. A
conditional facility repays by month eight using foreign reserves and voluntary
term holdings, with only 91 million in additional circulating local balances.
Its benefit depends on collateral, saving contracts and reserve access; neither
low inflation nor a U-related currency basket supplies those automatically.

The [integrated returns synthesis](../complex_returns/consolidated.md)
brings staged capital, worldwide services, U financing and resilience into one
set of dated accounts. Strong-demand sponsor cash NPV is 2.64–4.47 trillion
USD2041-equivalent at 3% over 100 years, including initial monetary-reserve funding;
all portfolios remain negative at 7%. It reports payback, internal return rates,
socioeconomic switching values and payment stress separately, and distinguishes
profitable staged capacity from completing the enclosed architectural build-out.
