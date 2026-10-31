# Construction of the Open Moon and its orbital array

The construction screen makes the cost of obtaining a habitable environment explicit beside the capital that supports ordinary production. Under its middle water inventory, combined lunar and array construction costs **11.2, 96.7 or 2,963.5 years of future Earth's output** across three deliberately different cost cases. Spread over 500 years, those amounts require **2.24%, 19.33% or 592.71%** of Earth's annual output at 9 billion residents and 100,000 common real units of annual gross product per resident. The last case exposes an industrial programme that exceeds that financing scale. The total covers environmental delivery and the initial or end-state capital inventory; renewal and operating expenditure during construction require separate funding. These are conditional scenarios, not forecasts or a newly selected design.

`construction.py` applies the assumptions in `construction_inputs.json`. The headline uses a common zero net-product-tax convention, so gross value added equals GDP. Monetary values use **U**, the same transferable common-price real unit as the world accounts. The physical inputs retain their product hashes in the input record and are checked against the bound inventory. The Moon has 8 billion residents, including 4 billion aloft; the array has 5 billion. The climate target remains 301 K, while revised climate fields and their effects on habitation and food requirements remain unresolved.

Air imports total 3.141×10¹⁸ kg, including 2.529×10¹⁸ kg nitrogen and 6.127×10¹⁷ kg oxygen. Seas and lakes require 1.382×10¹⁹ kg water. The three inventory cases add zero, 1.5×10¹⁹ or 5.93×10¹⁹ kg crust water, giving **1.696, 3.196 or 7.626×10¹⁹ kg delivered**. The middle crust quantity is an illustrative choice within the physical bracket.

| Import assumption | Patient automation | Advanced industrial | Expensive energy, dilute ore |
|---|---:|---:|---:|
| Energy, U/MWh | 0.01 | 0.1 | 1 |
| Recoverable nitrogen in source feed | 3% | 1% | 0.3% |
| Ore processing, MJ/kg feed | 0.2 | 0.7 | 2 |
| Ore equipment/service, U/kg feed | 0.000005 | 0.00002 | 0.0001 |
| Nitrogen chemistry, MJ/kg nitrogen | 5 | 10 | 20 |
| Oxygen chemistry, MJ/kg oxygen | 20 | 22 | 30 |
| Chemical equipment/service, U/kg air | 0.00002 | 0.0001 | 0.0005 |
| Departure propulsion, MJ/kg delivered | 20 | 40 | 120 |
| Departure equipment/service, U/kg delivered | 0.00004 | 0.0002 | 0.001 |
| Capture controls, MJ/kg delivered | 0.5 | 1 | 3 |
| Capture equipment/service, U/kg delivered | 0.00004 | 0.0002 | 0.001 |
| Placement, MJ/kg delivered | 1 | 1 | 3 |
| Placement equipment/service, U/kg delivered | 0.00001 | 0.00005 | 0.0002 |
| Import contingency | 10% | 25% | 50% |

Source material is assumed to contain 50% recoverable water. Joint feed equals the larger of `(water + 1.125 × oxygen)/0.5` and `nitrogen/nitrogen_grade`. Water and nitrogen therefore share extraction costs when recovered from the same material; separate deposits would increase throughput. The middle inventory requires **8.43×10¹⁹–8.43×10²⁰ kg raw feed** across the cost cases. The resulting delivered prices are **0.000185, 0.002465 and 0.082709 U/kg**. All equipment/service prices include the supplier's machinery, renewal, operation and capital recovery. Electricity includes generation capital service; neither the source plants nor their power systems are added again.

The process energies have physical reference points. The project’s departure velocity changes, 3.67–5.64 km/s, correspond to approximately 7–16 MJ/kg payload kinetic energy before inefficiency. Reaction propulsion also accelerates exhaust, increasing energy demand according to its mass ratio; the [NASA rocket-equation derivation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/) supplies that distinction. Reaction enthalpies derived from NIST data are approximately **17.9 MJ/kg oxygen** from water and **3.28 MJ/kg nitrogen** from ammonia. Actual ore separation, salts and organic nitrogen require additional processing, and reaction enthalpy is distinct from electrical work. [NIST water](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=2), [NIST ammonia](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417&Mask=1).

Incoming material carries another **19–70 MJ/kg of braking energy** in the repository's transfer examples. The calculation charges capture hardware and control energy while assuming that hardware rejects or recovers the incoming energy. It charges no purchased electricity for the kinetic energy itself and credits no sale of recovered energy. A catcher network able to perform this operation remains an unresolved engineering requirement.

| Cost case | Crust water | Combined construction, Earth-output years | Annual build / Earth output |
|---|---|---:|---:|
| Patient | None | 8.45 | 1.69% |
| Patient | Middle | 11.19 | 2.24% |
| Patient | Deep | 19.73 | 3.95% |
| Advanced | None | 62.98 | 12.60% |
| Advanced | Middle | 96.66 | 19.33% |
| Advanced | Deep | 196.13 | 39.23% |
| Expensive/dilute | None | 2,033.53 | 406.71% |
| Expensive/dilute | Middle | 2,963.53 | 592.71% |
| Expensive/dilute | Deep | 5,710.13 | 1,142.03% |

The denominator is **900 trillion U/year**. At 50,000 U annual gross product per Earth resident these ratios double; at 250,000 U they become 40% of the values shown. Prices and capital allowances remain fixed while the output benchmark changes. The calculation freezes each case's technology and prices across all 500 construction years; learning and historical price paths require a separate time-dependent model.

Ordinary installed capital is an explicit **300,000/500,000/1,000,000 U per resident** allowance, identical on Earth, Moon and array. It covers normal housing and habitat structure, productive equipment, food systems, ordinary computing, intangibles and R&D. The physical mass-times-price proxies—21,500/150,000/1,000,000 U per resident—are included diagnostics, never additional charges. Their broad correspondence to total capital is a judgment. Equal ordinary capital on the array is optimistic about pressure hulls, rotation and redundancy.

For the middle advanced case, the partition is:

| Capital or development expenditure | Earth | Moon | Array |
|---|---:|---:|---:|
| Ordinary installed capital | 4,500 | 4,000 | 2,500 |
| Imported planetary air and water | — | 78,786.3 | — |
| Aerial carriers/control | — | 420 | — |
| Ground receivers and critical reserves | — | 9.2 | — |
| Environmental establishment allowance | — | 450 | — |
| Optical shield | — | — | 231.5 |
| Orbital generation/transmission | — | — | 10 |
| Magnetic protection allowance | — | — | 90 |
| Passive habitat shielding | — | — | 495 |

Values are **trillion U**. Earth's ordinary capital is a matched reference stock, excluded from the lunar-plus-array construction total. Array shielding is 990 tonnes per resident; its ordinary 15-tonne habitat structure is already within ordinary capital. The aerial carrier allowance is 35 tonnes per aerial resident, excluding district payload already within ordinary capital and reserve lift that is not material.

Additional hardware prices are scenarios: aerial carriers cost **0.5/3/15 U/kg**, optical films **0.5/5/30 U/kg**, passive array shielding **0.01/0.1/1 U/kg**. Orbital power costs **0.01/0.1/1 U per operating watt**, including installed spare generation and transmitter hardware. Ground receivers cost **1/5/20 U/m²** and critical reserves **0.001/0.01/0.1 U/Wh**. Environmental establishment allows **90/450/1,800 trillion U**, including initial aerial hydrogen manufacture/fill, ecosystem seeding and soil engineering. Gas fill has no independent cost estimate; volatile import mass is counted only once. Magnetic protection allows **18/90/450 trillion U**. These two allowances represent unquantified work.

Construction expenditure follows the asset's territory. The totals include development and supplier services; they are not established market asset values or a depreciable stock. The array produces shielding and power services even when lunar residents benefit. Funding shares, ownership and service contracts are separate assumptions; lunar benefit does not transfer the array's depreciation into lunar net product. The annual ordinary-CFC sensitivity also remains independent of this construction screen. At the central capital allowance and 100,000 U annual output per resident, CFC of 15% corresponds arithmetically to a roughly 33-year average replacement period. That comparison does not establish asset lives or a complete depreciation schedule.

**Construction requires a far larger power industry than mature lunar living.** The middle inventory needs **46/99/366 PW** of purchased process power across the three cases. At 100-year equipment life, a 0.5% real capital rate, perfect availability and zero fuel or operating expense, their electricity prices allow generation capital of at most **0.00689/0.0689/0.689 U per continuous watt**. This is a feasibility diagnostic, not another expenditure. Solar flux at 40 AU is 1/1,600 of its 1-AU value; cheap inner-system generation alone cannot supply the assumed cheap source-site energy.

The calculation tests accounting conservation, coproduct feed, inventory and grade effects, territory, capital inclusion and schedule scaling. The thirteen construction tests include source binding and committed-product reproduction. This verifies the arithmetic and classifications. Source-site generation, bulk capture, materials lifetimes and the unquantified establishment work remain conditions of the scenarios.

Reproduce with `python -m research.studies.world_economy.construction`; use `--output research/studies/world_economy/results/construction.json` to regenerate the committed product. Full verification is recorded in [checks.json](checks.json).
