# U trade transition verification

12 November 2026. The transition plan is a proposed monetary and commercial development sequence for the currently planned off-Earth economic zone. It compares a single U currency with a basket of U-related currencies and includes separately funded partial reserve backing. Institutions, constituents, monetary mandates, adoption targets and funding shares remain proposals or illustrative assumptions.

## Numerical and structural checks

`python3 /tmp/check-u-transition.py --plan` selected the monetary-system layer and test scope. This wrapper imports the current main `shared/check_changed.py` with its root set to the isolated humanities checkout; it does not alter branch tooling. `python3 /tmp/check-u-transition.py` passed **41 tests**, comprising the 24 existing monetary tests and 17 transition tests. The layer check examined **1,655 files with zero violations**.

The new tests establish paired financial claims and net-worth reconciliation through all trade steps; identical real trade and fixed U stock across holding choices; unchanged total foreign U holdings when a customer draws dealer inventory; common valuation-unit rescaling; the basket's constituent identities and external valuation arithmetic; reserve source/use conservation, withdrawal shortfalls and valuation losses; transaction stocks proportional to holding time; source-hash rejection; exact dated receipt agreement; failed-path exclusion; missing-checkpoint rejection; and reproduction of the saved product.

The consumer reads the pinned `complex_returns/results/strategies.json` product, schema `terluna.research.export-strategy-returns/1`. Its full hash and producer/input identities are retained in the transition inputs and result. The 108 portfolio/date rows comprise 73 funded comparisons and 35 excluded failed paths. Each funded row contains nine settlement-share/holding-period alternatives. These are source-value equivalents, not a USD/U conversion or an addition to receipts. The financial producer, portfolio results and original monetary `results/reference.json` were not changed.

JSON parsing, the new study's local Markdown links and `git diff --check` passed. All source/dependency paths in the new product resolve. The full repository suite was not run because the change is bounded to a new read-only consumer, synthetic transactions and research prose. No new figure or manuscript PDF was produced.

## Source and substantive review

`trade_transition_sources.json` records targeted primary and institutional reading, with the external-basket comparisons kept separate from the U-family basket proposal. The unavailable IMF private-ECU article remains a recorded retrieval failure and supplies no detailed claim.

Two bounded internal reviews examined the monetary architecture, distinction between invoice/settlement/held balances, constituency and basket rules, payment-claim credit risk, reserve ownership, capital allocation and partial coverage. The final numerical reading verified the fixed-stock trade table, +40 external-position change, 130 = 95 + 11 + 24 reserve funding account, withdrawal and haircut shortfalls, the $1 trillion flow example and all nine year-100 transaction-stock comparisons. The completed prose was read for the specified economic boundary, explicit U-related constituents, absence of an implied Earth-currency peg, proposed adoption milestones and the distinction between the existing price screen and a monetary forecast.

The review led to explicit treatment of basket reference value versus the executable value of a particular issuer's claim, constituent entry/exit/default, and contract continuity when basket quantities change. These are internal accounting and editorial checks, not empirical validation, a full literature review or independent peer review.

## Limits

The fixed-stock transaction cases take an illustrative FX quote as given. They do not price dealer inventory, solve a currency market or estimate new issuance. The partial-reserve example is separate from those opening balances and from the existing industrial seed budget. Its allocations show where funds come from and which other uses they displace; they do not establish reserve adequacy or revise financial returns. The basket's synthetic constituents do not designate actual issuers.

An integrated monetary forecast still requires the selected constitution and policy framework, measured contract and payment currencies, reserve ownership and access, credit issuance, customer adoption and balance demand, and endogenous price/production/investment responses. No forecast of U appreciation, depreciation, purchasing power or annual inflation is claimed.
