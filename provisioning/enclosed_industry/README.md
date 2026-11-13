# Industry sized for enclosed communities

This model sizes equipment from the work requested each year. It separates the settlement's pressure enclosures, ordinary equipment and utilities from the mines, factories, power systems and freight vehicles that provide them. Population determines resident requirements and occupied residential stock; it does not multiply a founding industrial fleet.

The reference retains the chemical freight routes, payloads, source grades, material recovery, process intensities, purchase prices and mean component lives of the pinned initial-industrial-system study at `146a558275f486c9fd9a5cf79908b7c6fdf21888`. It adds no technological learning. Its quantities are conditional engineering requirements, not demonstrated plant performance, qualified mineral reserves or a funded settlement programme. It covers enclosed communities and their actual supporting facilities, before Open Moon or planetary-array construction.

## Inputs and outputs

[The source exporter](../../research/studies/initial_industrial_system/export_enclosed.py) freezes the inherited inputs as [engineering.json](inputs/engineering.json), preserving source file identities. The domain implementation reads this product and imports no research implementation. The separately pinned [orderbook](inputs/orders.json) supplies occupied population, construction orders and actual commercial dispatch. It also supplies service and hotel activity, equipment-package allowances and the market receipts later used by the financial consumer; the physical model does not treat those receipts as productive capacity.

[The runner](run.py) writes a versioned product containing source/model/input hashes, evidence and reading rules. Each case has a separate preproduction seed account and annual physical, procurement and operating accounts. Full component orders and asset inventories can be written to ignored run files. The compact product keeps the important asset stocks, imported-component purchase drivers, workforce and capacity screens. A failed or unconverged case retains its failure; it supplies no qualified sales capacity. Final annual products use deterministic gzip compression with compressed and decompressed hashes in a companion manifest. Compression preserves every annual row.

All cash results use constant USD2041. The inherited USD2026 prices are multiplied once by `1.03**15`. These are fixed purchasing-power accounting units, not a future exchange-rate forecast.

## Stocks, replacement and construction

Each asset has separate active and structural components. Active hardware retains the inherited external-supplier classification; locally fabricated structures use the inherited component partitions. A component enters service in its commissioning year, supplies capacity while it survives, and retires without a sale value. Surviving spare equipment can serve later work. Equipment that is no longer required receives no automatic replacement.

The main calculation uses a Weibull survival function with shape three and the inherited mean life. Its scale is solved so that the sum of annual surviving-stock fractions equals the inherited mean life. This annual calibration avoids the extra half-year of service that beginning-of-year sampling would otherwise add to a continuous distribution. This distribution spreads expected retirements across ages while retaining the reference mean; it is an uncalibrated scheduling assumption. A deterministic-age comparison is available. Expected fractional component/fleet stocks do not establish stochastic reliability.

Orders precede commissioning by two years. Half their fabrication and freight occurs in each construction year. The seed's procurement precedes opening: its component cohorts commission in years zero through seven, and their gross purchased masses account for survival losses to year eight. Additional orders placed in years six and seven are charged once in the seed account. The seed industrial stock is sized against an eight-year construction-throughput screen and opening operations, then retained and reused. The common inherited research/qualification budget is charged once. A separate [bootstrap addendum](bootstrap.py) also prices the inherited starting Earth-route fuel fraction against the actual opening fleet, including Earth purchase and delivery to low Earth orbit. Remaining fuel/storage filling and first-source deployment requirements remain explicitly unpriced; the reported seed total is the costed deployment scope, not a complete executable bootstrap quote. This average bootstrap screen does not establish an executable sequence of mine qualification, staff accommodation, first plant deployment or passenger migration.

The seed and annual accounts price the identified equipment, freight, materials, resource provisions, core payroll and external engineering support directly. They do not apply the former establishment account's 35 percent programme contingency or its separate stock-installation percentage. The once-only research/qualification allowance is USD2026 100 billion, or USD2041 155.797 billion. Installation and commissioning beyond the represented workforce and supplies remain unqualified. The smaller seed account therefore cannot be read as a complete replacement quote for the former full founding campaign.

A conservative fixed-point calculation closes the annual plan. Imported new and replacement hardware determines freight; freight determines propellant, source products and process power; those flows determine the required mines, factories, power, heat rejection and vehicles; their dated component requirements return to the procurement account. Each dated capacity target may increase during solution. It can differ from its neighbour and decline after construction, so this iteration is a feasible capacity envelope rather than an optimizer proving minimum capital cost. The reported end years have forty additional constant-demand continuation years to close construction feedback near the boundary.

Machinery and vehicles retain the inherited 25 percent spare-capacity rule and availability factors. The reported capacity ratio is **required equipment including that margin divided by available equipment**: a ratio of one already includes the margin. Workforce is a separate screen against 40 percent of residents, retaining the prior 60 percent working fraction and 20 percent civic/other reservation.

## Materials, transport and power

Local retired structures enter explicit bulk and metal recycling inventories. Virgin demand is the fabrication requirement that those inventories cannot supply. The source processes the maximum of its water, nitrogen and virgin-metal feed requirements, rather than adding three coproduct feed bills. Separate fields report source cargo, depot-water and source-water components, so a large source fleet can be traced back through its actual payload and propellant demand.

The ferry calculation is algebraically equivalent to the inherited tank-limited chemical-trip solver. It retains payload, return cargo, dry mass, delta-v, specific impulse, propellant margin, transfer losses, cycle time and availability. Reused export capsule buses are launched once per sortie; purchasing a bus does not add another launch charge to that sortie. Passenger transport remains outside this cargo calculation. Cargo and recycled-material arrivals retain the inherited annual-rate route convention; transit inventories and a mission-by-mission commissioning sequence remain to be qualified.

Source ore preparation uses the inherited energy per kilogram of feed, while nitrogen and metal conversion retain their distinct product coefficients. A poorer source grade therefore increases feed processing, power, heat rejection and machinery. Lunar metal processing retains the inherited aggregate 35 kWh/kg-metal scenario. Lower lunar grade increases raw feed and handling; this calculation does not invent a separately calibrated refining-energy law for that grade change.

Service equipment and floor area use inherited research-facility intensities per staff member as an explicit proxy. Hotel stocks, visitor water/air, consumables, power and staff use the inherited per-bed and guest-year coefficients. These provide priced physical demands rather than a complete qualification of sports venues, studios, premium-production processes or hotels. Their aggregate package quotes remain consumer-side allowances. Core technical staffing follows an inherited asset-value intensity, alongside explicit resident food/ordinary-production, factory, computing, service, premium, logistics and hotel jobs; that staffing law remains a planning proxy.

## Consolidated cash accounts

The physical calculation charges imported equipment, freight, consumables, external services and resident core payroll. Its construction-versus-operation allocation evaluates a counterfactual with the same occupied residents, actual exports, replacement work and available recycling but no new-capital manufacture. The incremental freight/material bill is assigned to construction. Capital and operating accounts sum to the directly computed physical cash bill.

Local component values remain capital-account diagnostics. Adding those internal transfer values to the cash spent on the same owned machines, supplies and workers would count domestic production twice. Only actual industrial resident payroll allocated to construction enters `capital_local_USD2041`; it is removed from operating payroll. The remaining core resident payroll is charged once in operation. Service, premium, logistics and hotel jobs enter the labour screen, while their wages remain in the separate commercial account.

`capital_total_USD2041` and `operating_total_USD2041` each include their own resource provision. Their provision subfields permit a financial consumer to retain that funding while a socioeconomic account removes the transfer from consumed-resource costs. The combined `resource_provision_USD2041` is a diagnostic, not another charge. Resource provisions do not buy a qualified source reserve or establish extraction rights.

`identified_service_stock_value_USD2041` and `identified_service_mean_renewal_allowance_USD2041` are noncash credits for the facilities already represented physically. They allow a consumer to define a residual unitemized service-package allowance without subtracting one year's smooth quote from another year's lumpy component replacement. Neither valuation is another physical expenditure. Hospitality is reported separately.

## Published comparison

The [base product](results/industry.json.gz) contains 81 combinations of nine commercial portfolios, three market strengths and three real market-growth assumptions, from year 8 through year 220. All 81 calculations converge and pass their dated equipment and resident-workforce screens. The [poorer lunar grade](results/poor_lunar.json.gz), [poorer source grade and recovery](results/poor_source.json.gz), and [deterministic component-life](results/lifetimes.json.gz) products each contain the services and diversified strong-market reference-growth cases; all six also pass those numerical screens. These results remain conditional on annual-rate transport, source access and the unqualified deployment requirements described above.

The services case illustrates why throughput sizing can still produce a large transport bill. At year 200, its 6.468 million residents and commercial guests require 21.94 billion kg of source water annually. The depots require another 76.13 billion kg, mostly for lunar and Earth ferries, while source-ferry fuel requires 112.50 billion kg. Producing their combined 210.57 billion kg of water at the inherited 7 percent recovered yield requires 3.008 trillion kg of source feed each year. Cargo of 98.93 billion kg per year requires 46,440 source tugs under the inherited three-year cycle, availability and spare-capacity assumptions. The fleet is a consequence of these assumed routes and flows, not a population multiplier.

The same case has USD2041 4.128 trillion of year-200 physical operating outlays: 2.244 trillion for imported replacement equipment, 1.318 trillion for operating freight and materials, 0.054 trillion for external services, 0.046 trillion for resident core payroll, and 0.466 trillion for resource provisions. Source-tug component purchases alone account for 1.354 trillion. The separate new-capital bill is 0.079 trillion that year. These are conditional costs before commercial receipts, package allowances or financing; the financial consumer supplies those accounts.

| Physical assumption, services strong-market reference case | Costed seed, trillion USD2041 | Year-200 operation, trillion USD2041/year | Year-200 capital, trillion USD2041/year | Year-200 source feed, trillion kg/year |
|---|---:|---:|---:|---:|
| Inherited grades, distributed component lives | 2.674 | 4.128 | 0.079 | 3.008 |
| Poorer source grade and recovery | 7.728 | 5.711 | 0.076 | 8.687 |
| Deterministic component lives | 2.618 | 4.008 | 0.711 | 3.058 |

The deterministic-life comparison retains large dated construction and replacement peaks. Its year-200 capital figure is a point on that different schedule, not a general lifetime-cost ranking. Poorer lunar yield is reported separately because it chiefly changes handling and material supply here; the aggregate lunar refining-energy coefficient remains unchanged.

## Reproduction and checks

From the existing provisioning checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m research.studies.initial_industrial_system.export_enclosed
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.run --output research/runs/enclosed_industry/raw.json
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.bootstrap research/runs/enclosed_industry/raw.json research/runs/enclosed_industry/priced.json
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.package research/runs/enclosed_industry/priced.json provisioning/enclosed_industry/results/industry.json.gz
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.run --case services:breakthrough:reference --scenario poor_source --output research/runs/enclosed_industry/poor_source.json
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest provisioning/enclosed_industry/tests -q
```

Updating either frozen input requires an explicit identity update before the runner accepts it. `--details research/runs/enclosed_industry/details` writes the full component ledger; `--case` restricts the requested cases. The other physical sensitivities are `poor_lunar` and `deterministic_lifetimes`. Completed case checkpoints are reused only under identical model, runner, engineering input and complete case identities. The model retains a single numerical worker.

Transport controls reproduce independently evaluated inherited ferry cases. Other focused tests check cohort timing, idle retirement, mean lifetime, seed gross purchases, coproduct feed, feed-bound heat and recycled-inventory conservation. Numerical convergence, provenance checks, accounting identities and engineering qualification remain distinct evidence states.

The separate [commercial and replacement constraint experiment](bottlenecks.md) preserves these baseline products while testing 28 bounded interventions. Its outputs retain failed workforce screens and explicit durability assumptions. [The diagnostic check record](bottlenecks-checks.md) records exact reproduction of the four baseline controls and the completed affected tests.
