# Commercial and replacement constraints

This diagnostic comparison tests which assumed constraints suppress early receipts or raise later costs in four representative portfolios: services, diversified exports, propellant/logistics, and the prior manufacturing mix. It preserves the accepted enclosed-population path, the physical model and the baseline products. The new [runner](bottlenecks.py) consumes a separately pinned [experimental orderbook](inputs/bottleneck_orders.json); commercial scenarios remain assumptions, and a requested sale receives no completed return estimate if its dated physical requirements fail.

The experiment holds the breakthrough customer scenario and 0.5 percent annual real market growth constant except for each stated intervention. The humanities order producer owns customer budgets, commercial capture, service demand and installation requests. Provisioning reruns the resulting machinery, replacement, freight, source feed and power requirements; it does not count a changed request as delivered output.

| Variant | Changed assumption |
|---|---|
| `baseline` | Original orders and technology. |
| `demand_25` | Buyer budgets and specified customer volumes increase by 25 percent. The audience ceiling remains two billion; paid conversion, licence volume and advertising/rights intensity change as recorded in the orderbook. Most unit prices remain fixed. |
| `capture_8` | Industrial/computing/fuel supplier offers grow at 8 rather than 4 percent per year, with the same initial scale, five-year opening ramp and customer caps. |
| `service_rollout` | Up to four service-package orders per year and a 20 percent service staffing allocation, replacing one order and 10 percent. Package requirements and the two-year lead remain unchanged. |
| `combined_rollout` | Both capture and service-rollout changes. |
| `demand_and_rollout` | Demand change plus both rollout changes. |
| `tug_life_30` | Original orders; source-tug active-component mean life increases from 20 to 30 years. Structural life remains 50 years and purchase price remains USD2026 1,500/kg. |

The first six variants separate demand from delivery timing and permit a demand–rollout interaction comparison. Their effects cannot be added as independent benefits. The last variant is a durability target whose qualification and incremental development costs remain unpriced. It does not change propellant consumption per trip, source grades, structural life or hardware prices. No technical improvement is adopted into the baseline.

The service staffing fraction is a commercial planning allocation; every variant must still pass the unchanged all-sector screen against 40 percent of residents. All added facilities, service power, workers, imported components and commercial transport enter the joint physical calculation. Faster rollout can therefore fail even when customers would pay for the requested output.

The [compressed result](results/bottlenecks.json.gz) preserves exact case names, experiment settings, original and applied tug life, convergence, failures, dated physical/cash rows and initial fuel procurement. Case caches require identical requests, engineering bytes and model/helper hashes, and each cached result has its own integrity hash. The bootstrap fuel addendum is applied once after the case calculations. Baseline source files and their earlier products are unchanged.

From the existing provisioning checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.bottlenecks
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest provisioning/enclosed_industry/tests -q
```

The runner uses one numerical worker. `--case` selects exact case names; `--no-cache` forces physical recomputation. Its uncompressed working result and cache remain under ignored `research/runs/enclosed_industry/bottlenecks/`. The manifest records both compressed and decompressed SHA-256 identities. Source access and mission-by-mission deployment remain unqualified, so a numerical capacity pass retains the baseline evidence boundary.

## Numerical results

All 28 calculations converge. Twenty-one pass the dated equipment and workforce screens; seven exceed the resident-workforce allowance. Faster capture first fails in year 34 for diversified exports and year 26 for the prior manufacturing mix. Combined rollout and demand-plus-rollout first fail in year 26 for diversified exports and year 21 for the prior mix. Service rollout alone first fails in year 41 for the prior mix. Each failure remains in the product with the required workload and zero admitted sales capacity in failing years. The services and propellant/logistics portfolios pass all seven variants.

The four diagnostic baselines reproduce their original seed accounts, every annual row, convergence record and first-shortfall field exactly. The longer-tug-life cases all pass. At year 200, their direct physical operating outlays are:

| Portfolio | Original active life, trillion USD2041/year | Thirty-year active life, trillion USD2041/year |
|---|---:|---:|
| Services | 4.128 | 3.592 |
| Diversified | 9.735 | 8.968 |
| Propellant/logistics | 4.236 | 3.676 |
| Prior manufacturing mix | 6.965 | 6.323 |

These are physical operating accounts before the separate commercial costs, capital allowances and financing. The changed replacement schedule also changes annual construction and seed requirements; a financial comparison must consume the full dated result. Source feed falls by about 1.4 percent at year 200, so the larger financial saving chiefly reflects expensive replacement components rather than a similarly large reduction in extraction. The unchanged lifetime-price assumption is a target sensitivity, with qualification costs still to be estimated.
