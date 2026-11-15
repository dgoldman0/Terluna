# Mining requirements for the enclosed build

For current availability estimates, read the [population-potential and reclamation assessment](../source_access/potential.md). It expands beyond the characterized source samples used below and quantifies fresh demand after recovery and refurbishment. The calculations here remain the reproducible physical-demand baseline and finite-sample comparison; their sample shortfalls do not establish a shortage in the wider near-Earth population.

This assessment follows the raw material required by an enclosed lunar settlement and its industrial support. The physical producer determines what must be built, replaced, processed and transported each year; the resources account compares that dated feed with explicitly hypothetical source allowances. Global service sales enter the demand study through their own machinery and operating needs. A larger audience does not automatically require a proportionate increase in mining or a replicated construction fleet.

Earth-side ore extraction is outside this feed comparison.

The saved assessment consumes 81 base cases and six physical sensitivities from provisioning commit `3e04591148223b079ecd255416555b8d4a4d4f7e`. It includes seed extraction in years 0–7 and annual operations from year 8 through year 220, with the prescribed 6.468 million residents reached in year 200. [inputs.json](inputs.json) pins all four physical products and their manifests. Deterministic gzip snapshots preserve the complete source bytes, and both compressed and decompressed hashes are checked. [sources.json](sources.json) distinguishes primary evidence from these conditional calculations; the earlier replicated-industry results supply no physical demand or costs here.

## Dated feed requirements

Every base case exceeds the entire mass of the original 17-body geometry screen before year 100. That finding establishes a shortfall in this particular known-body set even before recovery losses; it does not establish exhaustion of undiscovered or differently screened sources. By year 200, the broader geometry screen also cannot supply any case at the assumed 10% accessible fraction. The physical build therefore requires additional qualified sources or a changed material and transport design, even where the conditional financial account can fund the machinery.

| Elapsed year | Cumulative source feed across 81 base cases, kg | Lunar area excavated, km² | Prospective lunar area at 25% access, km² |
|---|---:|---:|---:|
| 100 | 3.925–7.027×10¹² | 18.8–32.9 | 75.1–131.7 |
| 200 | 1.643–2.915×10¹⁴ | 493.7–1,092.4 | 1,974.6–4,369.6 |
| 220 | 2.249–4.216×10¹⁴ | 505.6–1,411.4 | 2,022.3–5,645.7 |

The areas use the density and depth assumptions below and count virgin lunar feed. The 20-year extension holds the planned population cap while the producer continues replacing equipment and supplying operations. Market and production cases remain scenarios, rather than probabilities or forecasts.

The following entries count cases that fit each hypothetical mass allowance, out of 81. None represents a qualified reserve or delivery schedule.

| Catalogue allowance | Year 100 | Year 200 | Year 220 |
|---|---:|---:|---:|
| Narrow screen, 10% | 0 | 0 | 0 |
| Narrow screen, entire mass ceiling | 0 | 0 | 0 |
| Wider geometry screen, 10% | 81 | 0 | 0 |
| Wider geometry screen, entire mass ceiling | 81 | 73 | 0 |
| Narrow 10% plus broad remainder 1% | 81 | 0 | 0 |
| Broad catalogue, 10% | 81 | 81 | 81 |
| Broad catalogue excluding Don Quixote, 10% | 81 | 81 | 75 |

The last two rows supply a useful upper comparison but no admitted extension of the near-Earth transport architecture. Qualifying orbital windows, composition, extraction and conservation conditions remains necessary even when the row passes. Conversely, a failed mass ceiling identifies a precise limitation of the listed inventory without inventing replacement sources.

The grade sensitivities were run for the services and diversified portfolios under the breakthrough customer scenario and reference market growth. Their year-200 requirements show why lunar yield and source recovery need separate treatment:

| Physical assumption | Services source feed, 10¹⁴ kg | Diversified source feed, 10¹⁴ kg | Services lunar area, km² | Diversified lunar area, km² |
|---|---:|---:|---:|---:|
| Reference grades and survival distribution | 1.646 | 2.310 | 494.8 | 720.3 |
| Lower lunar yields | 1.646 | 2.311 | 818.8 | 1,266.0 |
| Lower source grades and recoveries | 4.795 | 6.724 | 494.8 | 720.3 |
| Deterministic component lifetimes | 1.631 | 2.297 | 492.2 | 718.0 |

Lower source recovery raises cumulative feed by about 2.91 times in these two cases after the physical producer closes equipment and transport feedback. Both then exceed the broad 10% allowance excluding Don Quixote by year 200. Under the broad 10% allowance including Don Quixote, the poor-source diversified case also exhausts the assumed stock before year 220. These comparisons bound particular changes; they do not predict when geology will deteriorate or validate the component lifetime distribution.

Chemical transport explains much of the source requirement. In the reference services case at year 200, the producer requires 210.6 billion kg/year of water-equivalent output: 21.94 billion kg for residents, guests and initial stocks, 76.13 billion kg for depot supply, and 112.50 billion kg for the source ferries themselves. Depot supply includes approximately 51.87 billion kg of lunar-ferry fuel and 24.21 billion kg of Earth-ferry fuel, plus oxygen makeup. At the assumed 7% recovered water yield, that becomes 3.008 trillion kg/year of source feed. The associated 98.93 billion kg/year cargo stream requires 46,440 modeled tugs and 579.5 million kg/year of imported tug components. These are outputs of the inherited chemical-route and replacement assumptions, not a demonstrated transport network. The feed budget and the financial account therefore need the same route qualification; ample aggregate mass cannot make the transport and replacement bill disappear.

## Material distribution and process evidence

The lunar sample analyses support regional sourcing choices. The Apollo 11 mare-soil example is richer in iron and titanium oxides than the unusually aluminous Apollo 16 example, while both contain material relevant to bulk construction. These are sampled compositions, with oxide mass distinct from recoverable, specified metal. The model must separately establish sorting, reduction, alloying and product acceptance. [NASA sample 10084](https://www.lpi.usra.edu/lunar/samples/atlas/compendium/10084.pdf), [NASA sample 67461](https://www.lpi.usra.edu/lunar/samples/atlas/compendium/67461.pdf).

Volatiles require source qualification. Bennu samples provide a measured nitrogen range of 0.23–0.25 wt%; their hydrated chemistry supplies a reason to investigate extraction rather than a universal asteroid grade. Lunar polar ice detections are patchy and concern the optical surface, so they supply no assumed excavatable depth or committed water inventory. Near-Earth asteroids generally remain in heliocentric orbits; departure windows, travel time and velocity requirements determine a source's connection to the Earth–Moon transport system. Material already parked in Earth orbit would require its own capture account. [Glavin et al.](https://ntrs.nasa.gov/api/citations/20250001355/downloads/GlavinAbundantSTI.pdf), [Li et al.](https://www.pnas.org/doi/full/10.1073/pnas.1802345115), [NASA/JPL trajectory documentation](https://ssd-api.jpl.nasa.gov/doc/nhats.html).

The process-energy literature supports calculating mass flow before energy. Leger and colleagues model a particular hydrogen-reduction chain at 24.3 ± 5.8 kWh per kilogram of liquid oxygen for 10 wt% ilmenite feed. Each process receives its actual inlet mass; beneficiation changes how much unreactive material must be heated. This oxygen result does not calibrate metal production or hydrated-asteroid extraction. The useful method is to distinguish feed-dependent heating from product-dependent conversion and finishing, including every recovery loss. [Leger et al.](https://publikationen.bibliothek.kit.edu/1000180152/157902327).

Demonstrations establish narrower capabilities than continuous industrial operation. NASA's May 2025 review reports molten-regolith processing of 25 kg of highland simulant over 36 hours and configuration-dependent water extraction in icy-regolith tests. The optical-mining account reports an 8 kW solar demonstration on asteroid simulant in vacuum. Neither establishes the mature plant's lifetime, maintenance labor, spare-part consumption or production learning rate. Those remain explicit physical scenarios; no automatic learning or repair saving is applied here. [NASA lunar ISRU review](https://ntrs.nasa.gov/api/citations/20250003730/downloads/Progress%20Review%20NASA%20Lunar%20ISRU_Sanders.pdf), [NASA-hosted optical-mining account](https://www.nasa.gov/general/optical-mining-of-asteroids-moons-and-planets-to-enable-sustainable-human-exploration-and-space-industrialization/).

## Grade, recovery and renewal sensitivities

The physical comparisons preserve the reference assumptions while varying the usable yield:

| Material | Reference usable yield | Adverse usable yield |
|---|---:|---:|
| Lunar bulk product | 0.90 kg/kg feed | 0.80 kg/kg feed |
| Lunar metal product | 0.10 kg/kg feed | 0.05 kg/kg feed |
| Asteroid water-equivalent product | 0.10 grade × 0.70 recovery | 0.05 grade × 0.50 recovery |
| Asteroid nitrogen | 0.0024 grade × 0.60 recovery | 0.0012 grade × 0.50 recovery |
| Asteroid metal | 0.08 grade × 0.60 recovery | 0.04 grade × 0.50 recovery |

These adverse values are sensitivities, not a predicted temporal decline. The physical producer recalculates extraction, handling, energy, equipment, replacement and transport for each admitted sensitivity. Co-products share feed: the binding product requirement determines joint extraction, while unused products acquire no automatic sale value. Additional feed heat must not duplicate heating already included in an aggregate refining coefficient. Recycling reduces virgin material requirements only after collection and recovery losses; returning material still passes through inspection, reprocessing and fabrication. The resource consumer charges neither a second recovery factor nor an economy-wide cost multiplier.

The source feed serves several products together. Its required mass follows the largest feed requirement for water, nitrogen or virgin metal at the declared grades and recoveries; adding those three feed requirements would count the same processed material repeatedly. Water includes resident and visitor makeup, depot propellant supply and the source ferry’s own fuel. That transport demand also determines ferry hardware and dated replacement imports. The physical producer closes this feedback before the resources account consumes its feed trajectory.

A broad catalogue mass allowance therefore does not establish the affordability or reliability of chemical transport. Water extraction, propellant preparation, delivery cycles and fleet renewal can remain expensive even when a hypothetical feed allowance is ample. Local lunar water and different freight systems would change those physical requirements, but each needs a separately sized and priced variant. Neither alternative is credited here as an unmodeled saving.

## Finite access allowances and lunar area

The retained catalogue calculation assigns spherical volumes and 1,500 kg/m³ density to reported diameters. Its narrow orbital-geometry screen, after exclusion of named scientific examples, contains 17 bodies with an assumed total mass of approximately \(3.1188\times10^{12}\) kg. The broader dark subset contains 498 bodies and approximately \(8.5112\times10^{15}\) kg. Their overlap is removed before defining a second tier. Missing diameters, unknown composition, conservation restrictions and unresolved routes remain outside this volume calculation.

The broad set is highly concentrated: Don Quixote accounts for 63.3% of its assumed mass. Its semimajor axis is 4.27 AU, eccentricity approximately 0.707 and inclination 31 degrees; its aphelion is approximately 7.29 AU. Those catalogue elements place it outside the narrow screen and supply no basis for using the inherited three-year cargo cycle. A broad-inventory comparison does not add an outer Solar System buildout to the planned near-Earth system.

An additional proposed screen admits dark, diameter-known objects with semimajor axis no greater than 1.8 AU, eccentricity no greater than 0.5 and inclination no greater than 30 degrees, retaining the named-science exclusions. It contains 131 bodies and an assumed mass of 2.1432×10¹⁴ kg. These orbital bounds define a reproducible comparison; transfer dates, velocity changes, navigation, composition and delivery cycles still need qualification. [catalogue.py](catalogue.py) reproduces this screen from the exact archived JPL source, and [the geometry product](inputs/geometry_screen.json) records its identity and selected object IDs.

Seven allowance comparisons use the narrow set alone, a wider geometry extension or the unscreened broad remainder. The narrow and wider screens each receive a 10% allowance and an all-mass ceiling; broad extensions use 1% or 10%, with a separate 10% comparison excluding Don Quixote. Removing that object leaves 3.1242×10¹⁵ kg in the broad catalogue before any extraction fraction. The all-mass ceilings ask whether even complete use of the known screen could cover the request; they are not feasible recovery assumptions. Narrow mass is subtracted from every extension before feed is allocated.

These are hypothetical raw-feed budgets. They assign no body to extraction, and every candidate still requires qualification and permission. Assuming an allowance available from the first modeled year provides a stock comparison; annual prospecting and delivery capability remain unestablished. The allocation function also supports explicit release dates and annual supply limits when qualified inputs become available. Missing bodies and unknown diameters remain unknown, rather than replenishing an exhausted allowance automatically.

The physical producer reports aggregate seed-construction feed separately from operating years. The consumer includes that feed once and allocates it evenly over the declared bootstrap duration immediately before operation. This dating preserves the producer’s average construction screen; it does not establish a detailed deployment schedule. Each row identifies bootstrap averages or operating demand.

Each annual request draws down the tiers in order. When their allowance is exhausted, the output records unallocated feed and its first year; the calculation neither invents another source nor catches up earlier shortages silently. A passing result means only that requested mass fits the chosen hypothetical allowance. `qualified_reserve_sufficiency` remains unknown even when that comparison passes. Each year-100, year-200 and year-220 comparison also carries a resource condition for financial consumers: completing the stated schedule depends on qualified source access and delivery, regardless of its conditional cash return.

For cumulative lunar virgin feed \(F\), the geometric comparison uses assumed bulk density \(\rho=1{,}500\) kg/m³, worked depth \(d=5\) m and accessible fraction \(f=0.25\):

\[
A_{\rm excavated}=F/(\rho d),\qquad
A_{\rm prospective}=F/(\rho d f).
\]

Thus \(10^{12}\) kg of feed requires approximately 133.3 km² actually excavated and 533.3 km² of prospective terrain under these assumptions. Neither number is a surveyed land allocation; processing sites, roads, storage and protection areas require additional accounting.

## Files and checks

[run.py](run.py) implements finite allocation and area accounting. [inputs.json](inputs.json) supplies the geometric and access assumptions; [the catalogue product](inputs/resource_inventory.json) retains its original bytes. [tests/test_run.py](tests/test_run.py) checks conservation, tier overlap, exhausted and delayed supply, seed inclusion and timing, invalid inputs, area conversion and source identity. The saved output is [results/assessment.json.gz](results/assessment.json.gz), schema `terluna.resources.enclosed-mining/1`, with producer and consumed-product hashes.

Run `python -m pytest resources/enclosed_mining/tests` for these numerical controls. `python -m resources.enclosed_mining.run` reproduces the assessment from the pinned products. Numerical reproduction verifies the calculation; empirical mine qualification remains a separate requirement.

The complete resources suite passed 64 tests with one optional, unrelated phosphorus-grid regeneration skipped. The assessment test reproduces all 87 cases, checks source identities, independently sums seed and annual feed, and verifies coproduct and allowance conservation. [checks.json](checks.json) records the commands and scope. Numerical checks establish these accounting and reproduction properties; geological surveys, process qualification and a dated mission schedule remain separate work.

The separate [commercial and replacement intervention assessment](bottlenecks.md) preserves this baseline account while consuming the diagnostic physical comparison. Its [check record](bottlenecks-checks.md) distinguishes exact baseline reproduction, numerical feed allocation and unresolved workforce or source-access conditions.
