# Revised enclosed-industry accounts

The [revised producer](revised_industry.py) applies the combined recovery and mechanical-refurbishment scenario to every one of the 81 original commercial orderbooks. It also evaluates six cases with poorer source grades and recovery, and six cases combining reclamation with faster commercial rollout. The [earlier 36-case reclamation study](reclamation.md), its inputs and outputs, the original physical baseline and the first asset export remain unchanged.

This is a physical and cost-accounting comparison for the planned enclosed communities and their supporting industry and transport. It retains the prescribed population path, reaching 6.468 million residents in year 200, and reports operation through year 220. Customer demand, resource access, repair qualification and the chemical freight routes remain scenario conditions. The financial and monetary consumers read these products separately; neither resource uncertainty nor a successful numerical solve is a probability-weighted financial return.

## Comparison scope and source identity

The primary cases preserve all original combinations of nine portfolios, three market strengths and three real market-growth assumptions. Their water, goods, structural recovery and delayed mechanical refurbishment follow the combined scenario already defined in the reclamation study, including its processing equipment, energy, spare parts, consumables and direct labour. The legacy 5 kg/person/day landscape/agricultural water allowance remains in every primary case. Its removal remains a separate endpoint in the earlier study.

The source sensitivity retains the services, propellant/logistics and diversified orderbooks under the strongest market, at 0.5% and 1% real annual market growth. Source water grade and recovery become 5% and 50%; nitrogen grade and recovery become 0.12% and 50%; recoverable-metal grade and recovery become 4% and 50%. These are the existing adverse engineering coefficients, not measurements of the selected source population. The solver recalculates ore handling, processing, power, freight, replacement and staffing instead of multiplying the entire economy's costs by an arbitrary factor.

The faster-rollout comparison uses six separately frozen humanities orderbooks for those same portfolios and growth assumptions. Annual commercial offer growth rises from 4% to 8%, the service-staffing allowance from 10% to 20% of residents, and the service-package order limit from one to four per year. Worldwide buyer budgets, prices, qualification rules, growth assumptions and the residential schedule stay unchanged. These are planning interventions, not demonstrated labour productivity. The physical producer consumes their dated requests without rewriting their commercial assumptions. Each case has an explicit `order_case_name`, canonical `baseline_case` and distinct `variant`: `combined`, `combined_poor_source` or `combined_rollout`. The nested `physical_case` retains the original orderbook metadata. Failed capacity or solver cases retain their requirements or error record and provide no qualified completed production path.

The [input record](revised_industry_inputs.json) pins the original engineering and orders, the prior reclamation producer and record, the asset valuation methods, and the additional rollout orders. The frozen baseline is used directly as the control for all 81 cases. The six combined cases already published in the reclamation study must agree exactly before entering this wider comparison; the original 81 physical cases are not solved a second time merely to repeat the control.

## Assets after refurbishment

The [revised asset product](results/revised_asset_stock.json.gz) covers the 18 primary cases with the strongest market and either reference or wider market growth. Its source is the same physical solve used for each cash and material account. It includes all registered component families, including the added water-treatment, sorting and repair facilities, and preserves functional sectors rather than assigning legal owners.

The gross commissioned stock is surviving component mass at its inherited replacement quote. Net stock uses the existing undiscounted remaining-service proxy. After successful refurbishment, a component receives the original expected service life under the stated engineering hypothesis. The old component retires once and the requalified component enters once. No additional asset is entered for the replacement invoice avoided by repairing it.

A repaired unit nevertheless contains both newly purchased parts and a retained core whose service potential has been restored. The export separates them. Suppose an accepted 100 kg unit requires 20 kg of imported repair parts and retains 80 kg of its original material. At a replacement quote of USD2041 2/kg, its restored complete-unit value is USD2041 200; the new-input quote is USD2041 40, and retained-core service restoration is USD2041 160. The latter is a conditional change in service-potential value, not a receipt, realized profit, or measured resale value. Failure to qualify the repaired unit's full assumed life would invalidate that full restoration.

During the two-year repair period, work in progress includes only the new-unit and repair-part quotes. The retained core receives no separate salvage or resale valuation while awaiting inspection or repair, and its restored service value is recognized only on commissioning. Repair labour, energy and consumables remain in the physical operating account; this export does not capitalize those expenses a second time. Consequently, the remaining-service stock is a valuation proxy rather than an acquisition-cost balance sheet or an estimate of accrued profit.

For each functional sector and year, the exported identities are

\[
G_t=G_{t-1}+C_t-R_t,\qquad
N_t=N_{t-1}+C_t-D_t,
\]

\[
C_t=A_t+S_t,\qquad
W_t=W_{t-1}+F_t-A_t.
\]

Here, \(G\) is gross commissioned replacement value, \(N\) is net remaining-service value, \(C\) is complete commissioned value, \(R\) is retired gross value and \(D\) is the inherited remaining-service consumption proxy. \(A\) is the quote value of new inputs embodied in commissions, \(S\) is retained-core service restoration, \(W\) is new-input work in progress and \(F\) is its dated manufacture accrual. Each quantity uses constant USD2041 replacement-quote equivalents. These identities reconcile stocks without treating internal structure quotes, restoration or avoided expenditure as additional cash flows.

The product retains the ordinary `net_remaining_service_value_USD2041` and `construction_WIP_value_USD2041` fields. Added fields identify new-input commissions, restored retained cores, repaired-unit commissions, imported repair-parts purchases and their masses. Detailed component stocks appear at years 8, 32, 100, 200 and 220; annual totals and functional-sector accounts cover every year from 0 to 220. Unitemized commercial packages, intellectual property, goodwill, mineral rights, human capital and ownership assignments retain their earlier unvalued or consumer-side status.

## Seed resource requirements

Revised cases retain the exact eight-year seed-average products and late-preproduction physical totals supplied by the reclamation producer. For the 81 frozen baseline controls, the original compact product contains total seed feed but lacks the late procurement's water, nitrogen and metal split. The revised producer therefore calculates the exact original seed average and subtracts its feed from the frozen seed total. That residual is the late feed, not an invented product allocation.

For a product with grade \(g\) and recovery \(r\), the late product amount is bounded between zero and \(gr\) times the residual feed. The corresponding complete seed-product interval adds the exact eight-year average product total to both bounds. These are physical identification bounds, not statistical confidence intervals. They preserve seed demand once, avoid substituting the water bound for all products, and make the remaining evidence limitation available to the resource consumer.

## Results and verification

All 81 primary combined cases and all six poorer-source cases converge and pass the dated equipment and resident-workforce screens. Four of the six faster-rollout cases also pass. Both diversified faster-rollout cases first exceed the available workforce in year 24 and remain unresolved completed paths; their requested sales are not admitted as qualified production.

| Comparison | Cases | Pass the complete dated screens | First failure, if any |
|---|---:|---:|---|
| Primary combined recovery and refurbishment | 81 | 81 | None |
| Combined with poorer source grades/recovery | 6 | 6 | None |
| Combined with faster commercial rollout | 6 | 4 | Diversified, both growth paths: year 24 |

Relative to each matching frozen baseline, the combined scenario reduces cumulative source feed through year 200 by 11.05–18.55%. Year-200 physical operating outlays fall by 9.37–26.64%. This operating account includes its represented resource provision and excludes customer receipts, commercial-account expenses, financing and consumer-side capital allowances; selected cases appear below. These results are conditional physical requirements and costs, not a finding that every portfolio is profitable or that its source inventory is accessible.

| Primary case, year 200 | Source feed, trillion kg/year | Physical operating outlay, trillion USD2041/year |
|---|---:|---:|
| services:breakthrough:reference — baseline / combined | 3.008 / 2.335 | 4.128 / 3.029 |
| propellant_logistics:breakthrough:wider — baseline / combined | 3.364 / 2.681 | 4.410 / 3.235 |
| diversified:breakthrough:wider — baseline / combined | 5.974 / 5.245 | 16.812 / 15.132 |

| Poorer-source case, year 200 | Source feed, trillion kg/year | Physical operating outlay, trillion USD2041/year |
|---|---:|---:|
| propellant_logistics:breakthrough:reference | 7.073 | 4.226 |
| services:breakthrough:reference | 6.694 | 4.088 |
| diversified:breakthrough:reference | 10.360 | 9.998 |
| propellant_logistics:breakthrough:wider | 7.692 | 4.458 |
| services:breakthrough:wider | 6.698 | 4.098 |
| diversified:breakthrough:wider | 15.039 | 17.509 |

All 18 revised asset cases reproduce the available component stocks from their corresponding physical solve. In the services/reference case, year-200 net remaining-service value falls from USD2041 32.430 trillion to 27.962 trillion, while work in progress falls from 3.740 trillion to 2.473 trillion. The revised path requires less transport and mining equipment for the same commercial requests, and its repair WIP excludes retained cores. A lower cost-based stock therefore does not establish lower delivered service or welfare. The separate year-200 retained-core service-restoration entry is 0.484 trillion; it is already part of the single requalified unit's value.

The full producer uses one numerical worker and integrity-checked case caches. Each cache requires identical producer bytes, scenario configuration, complete requests and dependency identities. Existing baseline and source products are immutable inputs. The [verification record](revised-industry-checks.md) distinguishes the complete numerical sweep, fresh-case reproduction, stock and material identities, source checks and prose review.

The physical product has schema `terluna.provisioning.revised-industry/1`. The asset product retains `terluna.provisioning.enclosed-asset-stock/1` and adds `source_physical_sha256`, which must match the compressed physical product used by the financial consumer. Case-level identities use canonical sorted compact JSON with a trailing newline. Consumer linkage uses the explicit original case and variant fields rather than inferring them from a name suffix.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.revised_industry
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest provisioning/enclosed_industry/tests/test_revised_industry.py -q
```
