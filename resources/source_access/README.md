# Lunar and near-Earth source access

The enclosed build has a plausible lunar bulk-material base, but its present asteroid catalogue is not yet a demonstrated supply base for two centuries of operation. Regional lunar composition is measurable; the depth, recovery and permitted footprint of a mine remain to be established. For asteroids, measured sizes and orbits are much better represented than composition, recurring cargo access or extraction conditions. The useful confidence statement therefore has several parts: measurement uncertainty where a survey supplies it, observed spatial variation, and explicit engineering scenarios where measurements are absent.

This bounded study covers the off-Earth economic zone served by U through the planned 6.47-million-person build at year 200, with checkpoints at years 100 and 220. Near-Earth sources here are asteroids in heliocentric orbits supplying cislunar facilities. Material is not assumed to be waiting in Earth orbit. Open Moon, the array and outer-Solar-System expansion remain outside this boundary. Costs retain constant USD2041 as a common valuation; they do not predict U's purchasing power or exchange rate.

## What the surveys establish

The [Lunar Prospector archive](https://pds-geosciences.wustl.edu/missions/lunarp/grs_elem_abundance.html) supplies 1,790 regional cells with six oxide abundances, potassium, thorium, uranium and a full nine-variable formal covariance matrix. The selected product covers the sphere in approximately five-degree bins. It measures the upper few tens of centimetres, not the five-metre worked layer used in the earlier excavation-area scenario. Its regional estimates retain systematic calibration uncertainty and imperfect oxide closure. [Prettyman et al. (2006)](https://doi.org/10.1029/2005JE002656) describe these limitations.

| Mapped oxide | Area-weighted p05 | Median | p95 | Median formal pixel sigma |
|---|---:|---:|---:|---:|
| FeO | 3.90 wt% | 5.58 wt% | 19.05 wt% | 0.40 percentage points |
| Al₂O₃ | 11.34 wt% | 24.66 wt% | 29.82 wt% | 2.30 percentage points |
| TiO₂ | 0.019 wt% | 0.476 wt% | 4.017 wt% | 0.209 percentage points |

The p05–p95 columns describe spatial variation across the mapped surface. They are neither confidence limits on the Moon's mean nor a distribution of future mine grades. Cell areas, including latitude-boundary intersections, provide the weights. The six oxide estimates sum to 88.69–107.09 wt%; the calculation preserves that nonclosure rather than normalizing it into an apparently complete assay. Eighty cells have zero reported TiO₂ and remain in the analysis. Spatial correlations and common calibration errors prevent averaging formal errors down as though every cell were an independent sample.

These contrasts make prospecting consequential. At an illustrative 70% recovery, the same oxide-equivalent quantity requires 7.50–36.61 kg of feed per kg FeO equivalent, 4.79–12.60 for Al₂O₃ equivalent and 35.56–7,435 for TiO₂ equivalent across the p95-to-p05 grade range. An oxide-equivalent quantity is composition accounting: these numbers establish neither a separable mineral nor a finished oxide or metal product. They show why an undifferentiated lunar “metal yield” cannot describe every industrial product. Titanium-rich sites deserve targeted characterization; the data do not show that scarce elements or nutrients can be substituted with common rock.

## Asteroid inventory and confidence

The retained [JPL Small-Body Database](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html) snapshot contains 43,050 near-Earth asteroids. Only 1,246 have diameters; 498 have both a diameter and reported albedo at or below 0.10. The calculation applies the existing illustrative science exclusions consistently, leaving 496 in the dark catalogue, 17 in the narrow geometry screen and 131 in the wider screen. These selections use reported values without pretending that a hard albedo boundary is physically exact.

| Candidate screen | Objects | Reference whole-body mass | Reported diameter-error coverage | Common-scale/density scenario |
|---|---:|---:|---|---:|
| Narrow geometry | 17 | 3.119 Gt | All 17; propagated sigma 0.055 Gt | 1.065–7.904 Gt |
| Wider geometry | 131 | 214.315 Gt | 127; propagated sigma 9.124 Gt on 180.712 Gt covered mass | 73.153–543.160 Gt |

One gigatonne (Gt) is one billion tonnes, or 10¹² kg. Reference mass assumes spherical volume and density 1,500 kg/m³. Narrow limits are semimajor axis 0.8–1.3 AU, eccentricity at most 0.3 and inclination at most 15°. Wider limits are semimajor axis at most 1.8 AU, eccentricity at most 0.5 and inclination at most 30°. Neither is a trajectory solution.

The formal calculation propagates reported diameter sigma through cubic volume to first order, assuming independent diameter errors and fixed density and membership. Four wider-screen objects, representing 33.603 Gt, lack diameter sigma; their mass is reported separately rather than assigned zero uncertainty. The other column applies a common ±20% diameter scale and densities of 1,000–2,200 kg/m³. It is an engineering stress range, with no assigned probability and no claim to enclose every possible density. A ±30% scale sensitivity widens the wider-screen range to 49.007–690.581 Gt.

This separation matters because catalogue fit errors can be much smaller than method uncertainty. [Nugent et al. (2016)](https://arxiv.org/abs/1606.08923) found approximately 20% diameter and 40% albedo scatter against independent measurements; [Masiero et al. (2021)](https://arxiv.org/abs/2107.07481) found approximately 30% relative diameter uncertainty in another NEOWISE comparison. Those studies motivate stress tests, not a universal statistical interval for this mixed catalogue. Sample-return density measurements near 1,190 kg/m³ for Bennu and Ryugu likewise inform scenarios without measuring the candidate bodies. Both scientific targets remain excluded.

Thirteen detailed JPL records were inspected, preserving physical-parameter references, uncertainties, rotation notes and orbit identifiers. The largest wider-screen object, 256412 (2007 BT2), contributes 32.2% of its assumed mass; the largest five contribute 72.3%. Its 4.444 ± 0.039 km catalogue diameter illustrates why a small formal error cannot settle composition or operational risk. Cacus is dark in this catalogue but has an S-type classification, demonstrating that low albedo alone cannot be treated as hydrated feed. The entire dark catalogue is dominated by Don Quixote, whose orbit reaches well beyond this source architecture; its mass receives no supply credit in the enclosed comparison.

### Independent composition and light-curve catalogues

The independent [MITHNEOS spectral survey, Marsset et al. (2022)](https://doi.org/10.3847/1538-3881/ac532f), adds 491 published spectra. The journal table and its CDS mirror agree on all identifiers and class labels. Number/designation matching identifies 419 distinct table objects, of which 417 match the frozen near-Earth catalogue. The paper reports 420 objects; one object appears under numbered and provisional filename prefixes. This discrepancy remains recorded rather than counting the filenames as separate bodies.

Combining JPL taxonomy with MITHNEOS gives classification records for 691 distinct near-Earth asteroids, compared with 327 in JPL alone. Fifty-three objects occur in both sources and count once. The records are not assumed statistically independent: JPL can cite the same original observations. Low-signal and unusual spectra remain unclassified, conflicting families remain conflicts, and single-colon uncertain classes retain their flags. Near-infrared C/X ambiguity remains a separate group.

| Wider dark screen: integrated spectral group | Classified objects | Share of classified objects | Share of classified reference mass |
|---|---:|---:|---:|
| C/B family | 6 | 35.3% | 74.0% |
| C or X, unresolved between them | 6 | 35.3% | 6.8% |
| S/Q family | 3 | 17.6% | 17.5% |
| X family | 2 | 11.8% | 1.7% |

The 17 objects cover only 13.0% of this 131-object pool and 20.1% of its assumed spherical mass. As a **conditional composition estimate**, the product assigns unclassified objects the observed mixture, once by object count and separately by known-size mass. The two estimates differ because classified bodies differ in size. The mass scenario would assign 74.0% of the pool to C/B; the observations themselves identify only 14.9% of total reference mass as C/B. Allowing the unclassified mass to belong anywhere gives a C/B identification range of 14.9–94.8%. This range assumes the observed labels and reference masses are correct; it does not incorporate diameter, density or taxonomy error and is not a sampling confidence interval. The observed-mixture allocation is useful for scenario planning, but discovery and observing biases prevent treating it as a measured population distribution.

A C/B spectrum is not a water assay, and X is not proof of a metal-rich body. The elemental or mineralogical split within these groups remains unresolved. Their principal immediate value is to select follow-up spectroscopy and sampling; the water-grade cost scenarios below remain independent assumptions.

The dark cut was inherited from water-source prospecting. Removing it while retaining the wider geometry, measured-diameter requirement and scientific holds produces **382 candidates with 667.3 Gt of reference whole-body mass**. Eighty-four have classification records; 75 have usable, non-conflicting groups, covering 59.0% of reference mass. The observed usable count mixture is 61.3% S/Q, 12.0% C/B, 10.7% C-or-X, 8.0% X and 8.0% other groups. This broader comparison matters for rock and metal opportunities, but adds no qualified water inventory. The 131-object screen is not a count of every potential near-Earth resource source.

For an internally consistent estimate across the 382-object pool, the model then separates 131 dark objects, 237 brighter objects and 14 without albedo. It allocates missing classifications within each group using that group's observed mixture. This prevents the overlapping dark and all-albedo pools from assigning incompatible masses when combined. The resulting **conditional distribution** is:

| Spectral family | Estimated share by object count | Estimated share of reference mass |
|---|---:|---:|
| S/Q | 54.8% | 60.4% |
| C/B | 15.3% | 23.9% |
| C or X | 14.4% | 2.2% |
| X | 8.6% | 13.0% |
| Other usable families | 6.9% | 0.6% |

These estimates correct only the choice to keep the albedo groups separate. They still assume that observed and unobserved objects within each group have the same weighted mixture; brightness, size and observing-selection biases remain. No object's missing label is replaced in the evidence record. The mass percentages classify the assumed whole-body inventory by spectral resemblance, not the fraction made of carbon, water or metal. They provide a research-informed planning distribution while preserving the much wider missing-data bounds in the product.

The independent [Asteroid Lightcurve Database](https://minorplanet.info/php/lcdb.php), credited to Warner, Harris and Pravec and using its 1 October 2023 release, contributes observations from the professional–amateur community. It matches 34 of the 131 wider candidates, representing 87.1% of their reference mass; 32 meet its period-quality threshold of 2− or better. Two binary flags, one suspected-binary flag and five tumbling-related notes remain attached for review. Those notes are not five confirmed tumblers, and blank fields do not establish operational simplicity. LCDB's size, albedo and taxonomy fields may contain assumptions and therefore do not replace the independent measurements used here.

All 17 wider candidates with spectral classifications also have light-curve records; only three additionally match the selected NASA human-trajectory query. This overlap identifies a small set for closer operational study, not three qualified mines. [ALCDEF](https://sbn.psi.edu/pds/resource/alcdef.html) also provided actual observations for the largest wider candidate, 2007 BT2: four blocks with 36 photometric measurements. Its catalogued period remains only greater than 20 hours, without a quality rating. No spin solution was fitted. ALCDEF's scientific-use restrictions are recorded, its photometry is not redistributed, and no commercial extraction design is claimed from those data.

### Access evidence

NASA's [Near-Earth Object Human Space Flight Accessible Targets Study](https://ssd-api.jpl.nasa.gov/doc/nhats.html), queried for 2020–2045 launches, round trips of at most 450 days, stays of at least eight days and total delta-v at most 12 km/s, contains 7,208 objects. Exact normalized-name joins find 23 of the wider-screen objects, totaling 8.510 Gt at the same assumed density. Three fall within its 8 km/s limit and one within 6 km/s. These are evidence-coverage comparisons for human Earth round trips. An absent object may fall outside the query conditions; a present object still lacks the project's recurring source-to-depot cargo solution. The inherited 2.5 km/s leg and three-year cargo cycle remain unqualified across the candidate pool.

## Extraction and processing cost ranges

The process screen uses the existing provisioning equipment rates and prices, freezes their identity in [the coefficient product](inputs/cost_coefficients.json), and varies utilization, electricity, equipment valuation and water yield explicitly. It does not import the provisioning model's implementation or replace the financial account. Capital is annualized at 3% real over the relevant equipment life, with the existing 2.5% annual supplier-upkeep provision and 25% capacity margin. This prices a standalone equipment service; a self-fabricating programme must reconcile internal transfers before treating it as external expenditure.

| Conditional cost component, constant USD2041 | Favorable | Reference | Adverse |
|---|---:|---:|---:|
| Lunar excavation, haul, handling energy and restoration provision, per kg feed | $0.111 | $0.200 | $0.559 |
| Lunar bulk material, including inherited bulk plant and process energy, per kg product | $0.324 | $0.750 | $2.378 |
| Asteroid feed extraction/thermal-processing screen and transition provision, per kg feed | $0.331 | $0.639 | $2.138 |
| Water from that feed at the process gate, per kg recovered water | $4.13 | $13.35 | $142.53 |

These are **partial conditional unit costs**, not quotations or confidence intervals for delivered products. They include neither a measured continuous-production lifetime nor the full cost of qualification, native-condition archives, freight, site construction, labour, common infrastructure, water cleanup, source-to-depot delivery or site-specific dust and waste containment. An omitted requirement retains an unknown cost. The 5% preservation capital allowance elsewhere in the project is not a substitute for these duties or an approved protection budget.

The favorable/reference/adverse cases use equipment-value multipliers of 0.5/1/2, delivered-electricity assumptions of $0.05/$0.15/$0.50 per kWh and availability multipliers of 1.2/1/0.6, capped at full availability. Source thermal energy is 0.25/0.50/1.00 kWh per kg feed; thermal-processor capital uses the inherited general chemistry cost per watt as an explicitly provisional proxy. It needs replacement with a tested extraction train. No second energy charge is added to the inherited aggregate metal-processing coefficient.

The water calculation is `C_water = C_feed / (grade × recovery)`, where feed cost sums annualized machinery, process energy and the stated transition provision. Grade and recovery change the tonnes handled per tonne of water; they do not create revenue or remove transport needs.

Recovered water yield is respectively 10 wt% × 80%, 6.84 wt% × 70%, and 3 wt% × 50%. The middle composition is a returned-sample analogue, not the candidate population's mean. [Yokoyama et al.](https://doi.org/10.1126/science.abn7850) measured 6.84 ± 0.34 wt% H₂O equivalent in Ryugu material, predominantly structural hydroxyl; release from those minerals occurred over roughly 300–800°C. This supports separating grade from release conditions. It does not establish 70% collection efficiency or authorize extraction from Ryugu. The three process cases require 12.5, 20.9 and 66.7 kg feed per kg recovered water. Holding process technology fixed instead of varying all assumptions, grade and recovery alone change that feed burden by a factor of 5.33 across the stated corners.

Lunar handling cost also inherits the mapped grade contrast. With reference handling cost and 70% illustrative recovery, p95-to-p05 feed grades imply handling charges of $1.50–7.34 per kg FeO equivalent, $0.96–2.53 per kg Al₂O₃ equivalent and $7.13–1,491 per kg TiO₂ equivalent. These exclude mineral separation, reduction, finishing and yield losses downstream. They must not be presented as iron, aluminium or titanium sale prices. In 39.2% of mapped surface area, TiO₂'s formal mean minus 1.96 sigma is nonpositive; that assay screen provides no finite upper cost bound there. The multiplier is a formal diagnostic, not an empirically established 95% coverage claim.

The inspected [NASA molten-regolith-electrolysis model](https://ntrs.nasa.gov/citations/20240013999) omits downstream metal handling/refinement and beneficiation. The [Leger et al. oxygen-chain calculation](https://doi.org/10.1073/pnas.2306146122) usefully illustrates the sensitivity of feed heating and beneficiation, but its modeled oxygen-energy spread is not a measured metal-production cost distribution. Lunar oxygen remains a by-product of the accepted materials architecture. Finished-metal cost confidence therefore needs a process-specific campaign, rather than a tighter-looking number derived from oxide abundance alone.

## Sustainability, physical preservation and data preservation

[The protocol assessment](protocols.md) maps the actual conservation decisions, proposed principles, source holds and transport requirements to independent evidence conditions. It preserves their different authority states. The study applies proposed requirements conservatively where relevant; it creates no new land allocation, extraction quota or permission.

Before disturbance, a candidate needs scientific and heritage review, representative survey and sampling, protection of affected pristine records, and a plan for preserving scientific functions. Samples require appropriate native conditions, unopened retained material, multiple locations and institutional continuity. Data preservation separately requires recorded provenance, scientific review and open access. A catalogue entry or hash satisfies neither the sample programme nor the archive's long-term stewardship. Future accessible stocks, scarce-element recovery and nutrient retention precede any surplus claim.

The six standing cargo-safety conditions remain distinct: passive failure behavior, Earth-clear catcher geometry, remote braking, packet-energy limits, tracking/interception/corridors and verified capture before a planetary intercept. Dust, exhaust and later outflow destinations require assessment. Open Moon atmospheric numerical budgets are not transferred into this enclosed-phase screen.

The previous broad-minus-narrow stock construction could reintroduce Bennu and Ryugu after excluding them from the narrow pool. The new calculation applies all four retained science holds before every pool operation, removing another 0.654 Gt from that broad comparison. Itokawa and Eros remain held as well. This list is illustrative; the remaining objects have no supplied site-specific protection attestations and remain **unresolved**, not approved. An unknown qualified quantity is not evidence of zero physical resources.

Polar ice receives no commodity-supply credit. The conservation decision calls for scientific excavation, analysis and storage before alteration. [LCROSS](https://doi.org/10.1126/science.1186986) measured an excavated plume at one site; [M³](https://doi.org/10.1073/pnas.1802345115) detected patchy exposed surface ice. Combining their different depths, footprints and definitions into a regional ore grade would manufacture confidence the surveys do not supply.

## Implications for the slow build

| Checkpoint | Cumulative asteroid feed, 81 scenarios | Cumulative lunar feed, 81 scenarios |
|---|---:|---:|
| Year 100 | 3.925–7.027 Gt | 0.141–0.247 Gt |
| Year 200 | 164.298–291.470 Gt | 3.702–8.193 Gt |
| Year 220 | 224.874–421.593 Gt | 3.792–10.586 Gt |

These totals include the one-time seed feed once, followed by the frozen producer's annual demand; the bootstrap is not replicated as population grows. At year 200 the reference wider-screen whole-body mass would have to supply approximately 77–136% of itself as usable feed. Combining the ±20%/density inventory envelope with demand extremes gives approximately 30.2–398.4%. Neither calculation reserves future stocks, establishes composition or removes site exclusions. The 17-object narrow screen remains insufficient even at its entire-mass ceiling before year 100.

The demand comparison retains the previous 10 wt% water and 70% recovery physics. It is not silently recalculated using the new sample analogue. Applying a different water yield to the complete programme requires provisioning to rerun source power, fleet replacement, fuel and recovery feedback together; multiplying aggregate demand by a grade ratio would not do that accounting.

At the previous illustrative density of 1,500 kg/m³ and worked depth of five metres, year-200 lunar feed corresponds to roughly 494–1,092 km² of excavation. This is a geometric translation, not a surveyed footprint or an approved portion of the Moon. The regional maps show opportunities for selecting better feed, but do not establish a grade-depth curve or a specific permitted site. There is therefore no evidence-based calendar at which common lunar bulk material must become more expensive. Costs rise when extraction moves to poorer grades, longer haul distances, harder material or more demanding preservation conditions; site selection and recovery can counter some of those effects.

For asteroid supply the more immediate concern is qualification and recurring access. The known wider inventory is concentrated in a few bodies, and even its reference all-mass total is comparable to year-200 demand. Completing the planned build therefore requires a larger **qualified** source pool, lower source throughput, different cargo economics, or some combination. Additional discoveries are possible; this calculation does not price or count undiscovered supply. Nor does it infer that the entire near-Earth population is insufficient from this deliberately limited screen.

The next evidence programme should prioritize sources whose independent spectra, size constraints and accessible trajectories overlap, obtain representative mineralogical and depth information before disturbing them, and demonstrate extraction and collection at measured energy and maintenance rates. Lunar work should compare contrasting candidate geological units at site scale while resolving their scientific and heritage conditions. Those observations can tighten recovery and cost ranges; enlarging a catalogue without resolving these links cannot by itself turn candidate mass into a sustainable production schedule.

## Reproduction and evidence records

[access_inputs.json](access_inputs.json) pins the original JPL snapshot, the new NHATS query, thirteen detailed JPL records, engineering coefficients and the existing 81-case physical-demand product. [lunar_inputs.json](lunar_inputs.json) pins the PDS table, label and archive documentation. [sources.json](sources.json), [lunar_sources.json](lunar_sources.json), [composition_inputs.json](composition_inputs.json), [lightcurve_inputs.json](lightcurve_inputs.json) and [the protocol mapping](protocols.md) record source roles and reading limits. External archive bytes remain ignored; derived numerical products retain source identity and credit. Full papers are not redistributed.

```sh
python -m resources.source_access.fetch_inputs
python -m resources.source_access.lunar
python -m resources.source_access.composition --download
python -m resources.source_access.lightcurves --fetch
python -m resources.source_access.run
python -m pytest resources/source_access/tests -q
```

The restoration commands refuse changed live API and catalogue bytes. The historical catalogue must be restored from its retained archive if unavailable locally; a fresh query is not a substitute for that snapshot. [results/access.json.gz](results/access.json.gz) carries the numerical assessment, [results/lunar.json.gz](results/lunar.json.gz) the spatial composition and full covariance, and [checks.md](checks.md) the actual verification record. These products support conditional comparison; they establish no mine, released quantity, probabilistic project success rate or revised financial return.
