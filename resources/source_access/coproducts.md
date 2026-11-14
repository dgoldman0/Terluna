# Precious metals, rare earths and other extraction coproducts

The earlier source-access calculation quantified water, native iron and major lunar oxides. It established a much broader geological supply than the characterized asteroid sample alone, but left precious metals and rare-earth coproducts largely unquantified. Those elements can materially improve the industrial opportunity: the same feed brought through water extraction or bulk-material processing also contains constituents useful for catalysts, magnets, electronics and specialty alloys. Their contribution needs an element account and a recovery circuit alongside the bulk account.

This assessment adds measured composition analogues, explicit recovery ranges and a comparison with the current industrial feed volumes. At the middle assumptions, one million tonnes of CI-like asteroid feed yields **1.72 tonnes of the six platinum-group metals**, **5,948 tonnes of nickel** and **273 tonnes of cobalt**. One million tonnes of feed resembling the measured lunar KREEP basalt yields **480 tonnes of rare-earth elements plus yttrium**, including **83 tonnes of neodymium and praseodymium**. These are conditional refined-element outputs; the low and high scenarios below retain geological and processing uncertainty. They require neither a claim of proven reserves nor an assumption that all extracted material is already saleable.

The boundary is the planned off-Earth economy served by **U, the project's general economic unit**, before Open Moon, the habitat array and outer-system expansion. U may be a currency or a basket of U-related currencies. The annual comparison uses the existing 81 combined-recovery industrial cases at years 100, 200 and 220. It retains their approximately 6.47-million-person year-200 build and supporting industry. Near-Earth asteroids remain heliocentric sources supplying the system, not material already parked in Earth orbit.

## Composition and the strength of the evidence

The [assay register](coproduct_assays.json) transcribes 40 elements from the [CI-chondrite compilation by Lodders, Bergemann and Palme (2025)](https://ntrs.nasa.gov/citations/20250008646) and 40 from [Neal and Kramer's analysis of Apollo 15 basalt 15386](https://www.lpi.usra.edu/meetings/lpsc2003/pdf/2023.pdf). CI is a meteorite compositional group used here as an analogue for primitive, water-bearing asteroid material. KREEP denotes lunar material enriched in potassium, rare-earth elements and phosphorus. Neither profile describes every asteroid or every lunar excavation site.

Concentrations are elemental mass parts per million: **1 ppm equals 1 gram per tonne, or 1 tonne in a million tonnes of feed**. The six platinum-group metals are platinum, palladium, rhodium, ruthenium, iridium and osmium. The rare-earth total used here, abbreviated REY, includes the fourteen naturally occurring long-lived lanthanides and yttrium; scandium is listed separately. Neodymium and praseodymium are already part of REY, so their subtotal must not be added again.

| Composition analogue | Selected measured or compiled concentration | What this establishes |
|---|---:|---|
| CI chondrites | Six platinum-group metals 3.225 ppm; gold 0.150 ppm; silver 0.206 ppm | Quantified precious-metal content in a bulk primitive-material analogue |
| CI chondrites | Nickel 11,180 ppm; cobalt 514 ppm; copper 133 ppm | Substantial alloy and industrial-element inventories; their chemical forms still require processing |
| CI chondrites | Gallium 9.54 ppm; germanium 33.5 ppm; indium 0.0781 ppm | Smaller specialty-element streams worth retaining and assaying |
| Lunar basalt 15386 | REY 902.45 ppm; neodymium plus praseodymium 155.6 ppm | A directly measured enriched lunar composition, not its spatial extent |
| Lunar basalt 15386 | Scandium 21.7 ppm; zirconium 916.5 ppm; niobium 72.4 ppm | Additional constituents potentially recovered from the same selected rock |
| Returned lunar soils in the USGS-hosted compilation | Lanthanum 5–70 ppm; scandium 8–63 ppm | Strong site and material variation; lanthanum alone does not determine total rare earths |

The soil comparison comes from [Lunar Sourcebook Table 7.15](https://lunar-consortium.astrogeology.usgs.gov/data/chem_min/tab7_15.pdf). Only its selected trace-element rows enter this interpretation. Regional thorium and potassium maps help identify KREEP-related prospecting areas, but do not directly assay each rare earth or establish grade at mining depth. The [NASA sample compendium](https://www.lpi.usra.edu/lunar/samples/atlas/compendium/15386.pdf) places much of basalt 15386's rare-earth inventory in phosphate-bearing material. Its very low siderophile abundances also show why a rare-earth-rich lunar rock cannot simply receive an asteroid precious-metal grade.

Measured meteorites support a wider compositional range than CI alone. The group means compiled by [Day, Brandon and Walker (2016), Table 3](https://msaweb.org/wp-content/uploads/2022/05/REV081C04.pdf) contain about **3.52–4.25 ppm** of the six platinum-group metals across the Rumuruti, carbonaceous, ordinary and enstatite chondrite means. Three measured iron specimens in that table contain **49.9, 88.5 and 136.5 ppm**. Those enriched specimens demonstrate possible compositions; they were not selected as a representative distribution of accessible asteroids.

[Harries (2024)](https://curation.isas.jaxa.jp/symposium/abstract/2024/P-1_Harries.pdf) reports an iron-meteorite platinum mode near **10 ppm**, with concentrations commonly reaching about **40 ppm**, while iridium spans roughly **0.01–60 ppm** at specimen scale. Differentiation therefore makes a single platinum-group ratio unreliable for metal-rich bodies. The model retains a separate platinum sensitivity of **1, 10 and 40 ppm**, giving **0.53, 5.32 and 21.28 tonnes of refined platinum per million tonnes** at the middle recovery. The 1 ppm value is a poor-grade stress; these three points have no assigned probabilities and do not bound all iron meteorites.

### Stronger named evidence, short of a reserve

The stronger evidence combines direct sample chemistry on the Moon and meteorites with observations of particular asteroids. [Sanchez et al. (2021)](https://arxiv.org/pdf/2109.13950) combine radar and spectral evidence for **6178 (1986 DA)** and identify **2016 ED85** as another metal-rich candidate. Their approximately 85% metal mixture is an interpretation of remote observations. It is not an assay showing 85% freely recoverable metal, nor a measurement of precious-metal concentrations.

Both candidates fail the current wider orbital screen's semimajor-axis and eccentricity limits. The retained catalogue gives approximately **2.81 AU and 0.585** for 1986 DA, and **2.98 AU and 0.695** for 2016 ED85, against limits of **1.8 AU and 0.5**. Their evidence belongs in the opportunity and prospecting account; this calculation assigns neither body tonnage to the existing feed schedule. The [candidate extract](inputs/coproduct_candidates.json) preserves catalogue diameter differences from the paper; neither diameter is used to value the body here.

Apollo assays establish composition at sampled locations much more firmly than remote estimates do. They still leave continuity, worked depth, mineral liberation and the protected or scientifically retained fraction open. No economically demonstrated precious-metal or rare-earth reserve is established by the sources examined here. Sample-return targets such as Bennu and Ryugu remain scientific analogues under the existing holds, with no mining credit.

## From contained material to recoverable output

For an element with feed grade *g*, the calculation follows:

**Refined element = original fresh feed × selected treated share × grade × stream retention × separation recovery × finishing recovery.**

Stream retention is the share of that element entering a recoverable concentrate, residue or condensate. Separation then acts on that inventory, and finishing measures elemental mass surviving purification. Finishing recovery is not a product-purity specification. Water loss can enrich a residue's concentration without creating extra precious metal; consequently every yield here uses the original feed denominator once.

| Engineering sensitivity | Retained in recoverable stream | Separated from that stream | Retained through finishing | Combined recovery |
|---|---:|---:|---:|---:|
| Low | 30% | 50% | 80% | 12.0% |
| Middle | 70% | 80% | 95% | 53.2% |
| High | 95% | 95% | 99% | 89.3% |

These are disclosed engineering assumptions for comparing possible circuits. Low-gravity plants have not demonstrated these overall recoveries for the chosen feeds. A plant that is absent, fails to capture the element or cannot separate it produces zero refined output. Using the same sensitivity levels across element families provides a common comparison; it does not establish that one circuit can recover every listed element simultaneously.

The geological sensitivities multiply the CI composition by **0.5, 1 or 2**, and the selected KREEP composition by **0.5, 1 or 1.5**. The former spans more than the differences among the compiled chondritic means; the latter tests dilution or enrichment around a particular rich lunar lithology. They are planning envelopes, not measured 5th–95th percentiles. Published assay deviations remain in the register, separately from transfer to another source and from engineering uncertainty. Whole-profile factors preserve the compositional association instead of assembling an ore from each element's independent maximum.

| Feed and product | Low | Middle | High |
|---|---:|---:|---:|
| CI-like: six platinum-group metals | 0.194 t | 1.716 t | 5.763 t |
| CI-like: gold | 0.009 t | 0.080 t | 0.268 t |
| CI-like: nickel | 671 t | 5,948 t | 19,978 t |
| CI-like: cobalt | 30.8 t | 273 t | 918 t |
| CI-like: gallium | 0.57 t | 5.08 t | 17.05 t |
| CI-like: germanium | 2.01 t | 17.82 t | 59.86 t |
| KREEP-like: REY | 54 t | 480 t | 1,209 t |
| KREEP-like: neodymium and praseodymium, within REY | 9.3 t | 82.8 t | 208.5 t |
| KREEP-like: scandium, additional to REY | 1.3 t | 11.5 t | 29.1 t |

All entries are **refined elemental tonnes per million tonnes of original feed treated**. High-case output can exceed the reference contained quantity because the high geological case has a higher grade. Product specification, plant throughput and access remain additional conditions.

## Scale within the planned industrial system

The reference case is the propellant-and-logistics portfolio with breakthrough demand and wider market growth. Its fresh asteroid feed reaches **2.68 billion tonnes per year at year 200**. Treating 1% of that stream means a **26.8-million-tonne-per-year** coproduct operation, so even this small share requires a substantial plant. At the middle assumptions it produces approximately **46 tonnes of platinum-group metals, 2.1 tonnes of gold, 159,000 tonnes of nickel and 7,300 tonnes of cobalt annually**. These are separate elemental inventories; some nickel or cobalt could already be embedded in the baseline's generic metal products.

| Reference-case checkpoint | Fresh asteroid feed | Platinum-group output from 1% treated | Fresh lunar feed | REY output if 10% is selected KREEP-like feed and treated |
|---|---:|---:|---:|---:|
| Year 100 | 205 million t/year | 0.40–11.8 t/year; middle 3.51 | 6.82 million t/year | 37–825 t/year; middle 327 |
| Year 200 | 2,681 million t/year | 5.19–155 t/year; middle 46.0 | 2.00 million t/year | 10.8–242 t/year; middle 95.9 |
| Year 220 | 2,722 million t/year | 5.27–157 t/year; middle 46.7 | 2.15 million t/year | 11.7–261 t/year; middle 103 |

The brackets combine the low and high geological and processing scenarios. The **1% and 10% shares are selectable treatment scenarios**, not estimates of how much qualifying material has been located. Across all 81 physical cases, middle-recovery platinum-group output at 1% treatment spans **40–99 tonnes per year at year 200**. The result file retains every portfolio and checkpoint, including inherited feasibility flags; none of those flags qualifies the new refining plant.

Increasing the reference asteroid treatment share from 1% to 10% raises the middle estimate from **46 to 460 tonnes per year**. Treating all feed gives about **4,601 tonnes per year** at the same assumptions. These sensitivities show how significant the opportunity could become while also exposing the need to size refining capacity and market demand. They are alternative plant scopes, not outputs to add together. Selecting a different asteroid grade can also alter water yield and the main feed requirement; a coupled industrial rerun must account for that feedback.

For geological scale, consider a hypothetical **100 km² lunar area worked to 5 m**, using the existing density-depth calculation. If it contained the selected KREEP-like material throughout, the combined density and grade scenarios would contain **0.35–1.25 million tonnes of REY**, with a middle value of **0.78 million tonnes**. Applying the processing scenarios gives **0.042–1.11 million tonnes refined**, middle **0.416 million tonnes**, before site restrictions and retention for future use. This is area-times-grade arithmetic for an unlocated deposit, not a mapped inventory or a confidence interval. It demonstrates why small sampled quantities need not imply a small geological resource.

## Cost, product value and the joint operation

Bulk extraction can pay for access, excavation and much of the initial handling. A coproduct account then charges the additional work it causes: selective concentration, mineral breakdown, reagents and their recovery, purification, equipment renewal, electricity, heat rejection, workers, records, residue containment and any freight. A change that diverts an essential alloy constituent or lowers bulk-product yield also carries that opportunity cost. Refining a residue is therefore potentially much cheaper than establishing a separate mine, without being free.

A useful early economic test asks what average realized product value would cover a proposed incremental circuit cost. At the reference grades and middle recovery, the following values would cover the **entire incremental circuit using the named product alone**:

| Product supporting the circuit | At $0.01/kg treated feed | At $0.10/kg treated feed | At $1/kg treated feed |
|---|---:|---:|---:|
| CI-like six-metal platinum group | $5,829/kg product | $58,285/kg | $582,852/kg |
| CI-like nickel | $1.68/kg product | $16.81/kg | $168.13/kg |
| KREEP-like REY mixture | $20.83/kg product | $208.29/kg | $2,082.88/kg |
| KREEP-like neodymium and praseodymium alone | $120.80/kg product | $1,208.03/kg | $12,080.33/kg |

All monetary values are constant **USD2041 purchasing-power equivalents**. The feed-cost levels are sensitivity inputs, not estimated plant costs or U exchange rates. Joint revenues can support the circuit together; the same cost must not be charged in full to every coproduct and then summed. Conversely, mixed output cannot all receive the price of its most valuable element. The platinum group contains six distinct markets, and neodymium plus praseodymium accounts for only about 17% of the selected KREEP REY mass.

The account should allocate refined output among domestic use, strategic inventory and exports. Domestic substitution avoids purchases only where it displaces an actual import or another production cost. Export proceeds require product quality, delivery and customers at the resulting supply volume. Existing generic metal and manufactured exports need their embedded elements reconciled before a new coproduct revenue line is added. **This study establishes additional recovery opportunities and break-even conditions; it does not yet revise financial returns.**

## Preservation, reclamation and process design

The [established protection protocols](protocols.md) apply to these high-value elements as they do to bulk extraction. Survey and representative sampling precede disturbance, with provenance, a native-condition sealed or unopened portion, distributed physical and data archives, and site-specific scientific and heritage review. Geological potential remains visible while access conditions are unresolved. Neither a valuable composition nor a completed accounting calculation releases a protected source.

Scarce-element recovery requires designing the residue streams deliberately. Metal and sulfide concentrates are candidates for precious-metal recovery; lunar phosphate and silicate residues require a different rare-earth route. Volatile constituents need capture before they leave a heated process. Each stream should retain its source identity, assay, mass and location so later processing can find the remaining material. Cleanly segregated, contained residues preserve more future options than indiscriminately mixed spoil or uncontrolled dispersal.

For KREEP processing, thorium and uranium accompany some valuable constituents: the selected sample contains **15.4 and 4 ppm**, respectively. Their separation and containment belong in the process and protection costs; the model assigns them no sales or reactor-fuel credit. Their inventory cannot be ignored while valuing the associated rare earths.

The element balance distinguishes untreated inventory, material outside the selected recovery stream, separation residue, finishing residue and refined output. It conserves the element across those compartments, but **does not assume every residue remains physically recoverable**. Engineering must determine capture, storage and eventual loss. Later recovery from indexed tailings is a second pass through the same inventory, not newly discovered primary supply.

The operating economy also needs return paths for catalysts, magnets, electronics, tools and fabrication scrap. Collection, disassembly, alloy sorting and purification determine how much of each element comes back. The existing bulk-reclamation benefits remain in the inherited feed calculation; no additional perfect-closure or trace-metal recycling credit is added here. A future stock-and-flow account can reduce replacement demand using measured element-specific lifetimes and recovery, while charging collection and treatment once.

## Earth receipts during the initial eight-year build

The accepted slow-growth schedule opens with **12,500 residents at year 8**. Its eight-year founding account processes **298.94 million tonnes of asteroid feed** and **7.703 million tonnes of lunar raw material** in the propellant/logistics reference case. These totals include late founding procurement once and exclude the operating year beginning at time 8. Across the 81 portfolios and market combinations, asteroid feed spans 291.90–307.11 million tonnes. This short-period estimate uses those existing totals; it changes neither the opening population nor the long-run growth schedule.

At the reference CI composition, the asteroid feed contains approximately **964 tonnes of the six platinum-group metals**. If every tonne of feed received the middle recovery process, it could yield **513 tonnes refined**. Earth receipts within eight years depend on when that refining capacity becomes available and when its cargo arrives. The [startup-delivery calculation](startup_deliveries.py) therefore separates production, retention within the off-Earth economy, successful arrival, return loss and cargo still awaiting delivery.

The illustrative programme starts coproduct refining in build year 3. Its source-feed treatment fractions during build years 1–8 are **0%, 0%, 1%, 2%, 4%, 6%, 8% and 10%**. Refinery alternatives use one tenth or ten times those fractions, reaching 1% or 100% of feed by year 8. Each assumes 80% of refined output is allocated to Earth, with 20% retained off Earth. A 1.5-year source-to-Earth delay and the existing capsule model's 1% return loss apply. The delay is a planning proxy alongside the inherited three-year source round trip, not a trajectory calculation.

| Source refinery coverage reached in year 8 | Platinum-group metals received on Earth within eight years, middle assumptions | Geological/recovery scenario range |
|---|---:|---:|
| 1% of asteroid feed | 1.07 tonnes | 0.12–3.60 tonnes |
| 10% of asteroid feed | 10.72 tonnes | 1.21–36.00 tonnes |
| 100% of asteroid feed | 107.18 tonnes | 12.09–360.03 tonnes |

These are capacity alternatives, not probability bands or an estimate of the geological fraction available. With no added coproduct refinery, this programme delivers zero. Across the existing 81 industrial cases, the 10%-coverage middle estimate varies only from **10.47 to 11.01 tonnes**: refinery scale, recovery and arrival timing matter much more than the portfolio differences during startup.

The middle 10%-coverage case treats **16.86 million tonnes of asteroid feed** cumulatively and refines **28.93 tonnes of platinum-group metals**. Approximately **5.79 tonnes stay off Earth**, **12.32 tonnes remain awaiting delivery or in transit at year 8**, **0.11 tonnes are lost on completed returns**, and **10.72 tonnes reach Earth**. Expected annual Earth receipts are zero in build years 1–3, followed by **0.17, 0.62, 1.58, 3.16 and 5.19 tonnes** in years 4–8. These smooth flows describe an assumed schedule; actual dispatch would batch cargo.

| Additional product, same 10%-coverage source case | Cumulative physical receipt on Earth by year 8, middle assumptions |
|---|---:|
| Gold | 0.50 tonnes |
| Silver | 0.68 tonnes |
| Rare earths plus yttrium from selected lunar feed | 148.92 tonnes |
| Rare earths plus yttrium from asteroid feed | 13.79 tonnes |
| Neodymium and praseodymium, already within the two rare-earth totals | 27.60 tonnes |
| Gallium from asteroid feed | 31.71 tonnes |
| Germanium from asteroid feed | 111.34 tonnes |

The lunar comparison separately assumes selected KREEP-like feed can supply the treatment ramp, reaching **10% of total lunar excavation in year 8**. It treats 0.434 million tonnes cumulatively, uses a quarter-year refining-to-Earth handling and delivery lag, and retains the same Earth allocation and return-loss assumptions. Its source selection remains hypothetical. The lunar total does not increase when asteroid refinery coverage changes. The combined rare-earth delivery is about **163 tonnes**, with **18–421 tonnes** across the fixed-capacity geological/recovery scenarios; neodymium and praseodymium must not be added a second time.

These quantities are potential physical receipts, conditional on product quality, funded refining and transport, and buyers. Bulk nickel and cobalt illustrate why recoverability alone is insufficient: the source circuit also refines approximately **100,300 tonnes of nickel** and **4,610 tonnes of cobalt** in this case. Applying the same export allocation and timing mechanically would deliver 37,157 and 1,708 tonnes. Those bulk flows are not included as expected Earth sales. The inherited return model already charges **$150/kg in 2026 purchasing units for Earth ground service alone**, before the other transport and refining costs. For scale, the USGS estimates 2025 nickel at about **$15/kg**; lower return costs or higher-value products would be needed to support that bulk export comparison. This contemporary price is context, not a future quotation. [USGS, Mineral Commodity Summaries 2026, nickel table](https://pubs.usgs.gov/periodicals/mcs2026/mcs2026.pdf).

The original model supplies an eight-year average, without a qualified mining deployment sequence. For this estimate, base feed rises in annual weights 1:2:3:4:5:6:7:8, while the separately recorded late procurement is split between build years 7 and 8. This preserves total extraction, but the last year's base flow is 1.78 times the inherited average. Peak mining capacity, the new refinery's equipment, power and staff, and their extra construction feed have not been solved together. The 10%-coverage case reaches about **6.65 million tonnes/year of source refining** and **0.171 million tonnes/year of selected lunar refining** in its final year. These are required capacities, not existing assets.

Holding composition and recovery at their middle values, alternative feed schedules and source-delivery lags of 0.5–3 years give approximately **3.3–21.5 tonnes of platinum-group metals received within eight years** for the 10%-coverage programme. That sensitivity is separate from the geological/recovery ranges. Cargo arriving after year 8 remains valuable inventory; it is excluded only from this receipt window. The calculation makes no revenue or investment-return credit, and the separate local-power comparison remains outside the inherited feed totals.

## Reproduction and integration

The [calculation](coproducts.py) reads the assay register, disclosed [scenario inputs](coproduct_inputs.json), pinned [industrial flow extract](inputs/coproduct_flows.json), existing lunar density product and protection evaluator. Run `python -m resources.source_access.coproducts` from the resources checkout. The [result](results/coproducts.json.gz) contains per-feed yields, 243 case checkpoints, treatment-share sensitivities, the hypothetical lunar patch and economic thresholds. [Sources and retrieval limits](coproduct_sources.json) identify the exact PDFs, hashes, tables and reading scope; full papers are not redistributed. [Checks](coproduct_checks.md) distinguish numerical accounting from geological, engineering and commercial qualification.

The industrial input is an exact extract of the committed combined-recovery product available at provisioning revision `adeaa8ea`. The extraction script records the original producer, product hash and its earlier source revision separately. The subsequent local-power comparison has not been recoupled into those physical cases. Integrating coproducts into returns will require a selected refining route and capacity, incremental capital and operating costs, an element account for existing exports, and demand-sensitive domestic and external prices. The present result makes that opportunity explicit without treating contained material as booked income.

For the startup extension, run `python -m resources.source_access.startup_deliveries`. Its [pinned founding-feed extract](inputs/startup_feed.json), [assumptions](startup_delivery_inputs.json) and [result](results/startup_deliveries.json.gz) retain source identities and all 81 case comparisons. The feed exporter reads the same committed provisioning product as the long-run coproduct study.
