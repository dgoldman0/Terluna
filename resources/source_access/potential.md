# Material potential and reclamation before the Open Moon

The planned enclosed economy has a much larger potential material base than the previously characterized asteroid sample suggested. Extending the calculation to the full catalogue, inferred missing sizes and an undiscovered population gives **25–73 Gt of contained water equivalent in the wider near-Earth orbital screen**, with a median of **39 Gt**. The broader near-Earth population, excluding Don Quixote, gives **504–2,931 Gt**, with a median of **961 Gt**. These are conditional 5th–95th percentile intervals before extraction recovery or restrictions on access. They support a substantial geological opportunity; the smaller screen alone does not establish the supply limit.

Reclamation changes the demand side as well. Water, household goods and structural materials already receive considerable recovery credit in the physical baseline. Explicitly improving those loops and refurbishing selected active machinery reduces cumulative fresh asteroid feed through year 200 by **12–19%** across six tested demand paths. The associated year-200 physical operating outlay falls **10–27%**, after the experiment charges repair parts, workshops, workers, electricity and their transport. This is a conditional industrial comparison, with the inherited commercial orderbooks and chemical freight architecture retained.

The assessment covers the economic zone served by **U, the project's general economic unit**, through the enclosed build of approximately 6.47 million residents at year 200 and its first 20 years of continued operation. It includes the supporting lunar and near-Earth industrial system; population is a demand input, not the definition of the economic unit or its stakeholders. Open Moon, the habitat array and outer-system buildout are outside this boundary. Near-Earth asteroids are heliocentric bodies; the calculation does not place their material conveniently in Earth orbit. Monetary quantities retain constant USD2041 purchasing-power equivalents without prescribing a U exchange rate.

## From catalogue measurements to population distributions

The [population calculation](population.py) integrates the retained JPL catalogue with MITHNEOS spectroscopy and the published NEOMOD population models. It quantifies 43,039 catalogued objects: 1,242 with measured diameters and 41,797 whose diameters are inferred. Four scientific targets remain excluded, and seven records lacking both absolute magnitude and diameter remain unquantified. Missing characterization therefore widens the inference instead of erasing an object's potential mass.

For an object without a diameter, absolute magnitude constrains size jointly with albedo. The calculation samples albedo using the size-dependent two-component distribution in [NEOMOD3](https://arxiv.org/html/2404.18805v1), weighted by its size-distribution model conditional on magnitude. Reported diameters take precedence. Catalogue albedo, JPL classifications and [MITHNEOS](https://doi.org/10.3847/1538-3881/ac532f) spectra inform material-family assignments. The broad debiased prior of 42% stony, 48% carbonaceous and 10% other/ambiguous objects is an analogue assumption for missing composition. It is not a measured family split in every magnitude or orbital bin. An X classification retains ambiguity and does not make an object a metal deposit.

The undiscovered component uses [NEOMOD2's magnitude-count distribution](https://arxiv.org/html/2312.09406v1) and the published NEOMOD orbital grid. Every known object in the relevant bin, including retained scientific targets, is subtracted before a missing population is inferred. This prevents both double counting and the reappearance of a protected object as an unknown source. Where the newer catalogue exceeds a model bin, the missing component is clipped at zero and the mismatch is recorded. No main-belt stock or future replenishment is added. Counts outside the inferred magnitude range are not extrapolated into an unlimited small-body inventory.

Density and water/native-iron composition use explicitly bounded asteroid and meteorite analogue distributions. Water means hydrogen-bearing material expressed as water equivalent, including bound forms; it does not imply free ice or uniform processing chemistry. All sampled fractions are converted from weight percent to fractions, and contained water plus native iron cannot exceed the body's mass. Native iron is tracked separately from generic structural-metal demand because a mass surplus cannot supply a missing alloy recipe.

| Population pool | Body mass, p05 / median / p95, Gt | Contained water equivalent, p05 / median / p95, Gt |
|---|---:|---:|
| Narrow orbital screen | 20.82 / 37.02 / 118 | 0.34 / 0.82 / 4.31 |
| Wider orbital screen | 1,449 / 1,777 / 2,289 | 24.92 / 39.09 / 72.94 |
| Broader near-Earth population, excluding Don Quixote | 66,403 / 97,126 / 129,311 | 504 / 961 / 2,931 |

One Gt is one billion tonnes, or one trillion kilograms. The narrow geometry is 0.8–1.3 AU semimajor axis, eccentricity at most 0.3 and inclination at most 15°. The wider screen allows semimajor axis up to 1.8 AU, eccentricity up to 0.5 and inclination up to 30°. Neither screen solves recurring cargo trajectories. The whole-population comparison is geological context; using its stock would require a broader, demonstrated access strategy within the intended near-Earth programme. Nested pools overlap and must never be added.

The wider pool's median contained-water contribution is approximately 35.1 Gt from catalogued bodies and 2.61 Gt from inferred undiscovered bodies; the median of their sum is 39.1 Gt. Medians need not add. Most of the increase over the earlier screen comes from using all albedos and inferring missing characterization, rather than from asserting an enormous unseen inventory. Amateur and professional light-curve catalogues continue to constrain rotation and binary follow-up, as documented in the [survey study](README.md#independent-composition-and-light-curve-catalogues); light curves are not converted into water assays.

### What confidence means here

Each central interval summarizes 2,048 coupled Monte Carlo realizations under the stated priors. It includes individual diameter uncertainty, latent albedo/composition, density and grade variability, and the missing-population component. Three population-count branches and common diameter-scale factors of 0.8/1/1.2 or grade factors of 0.5/1/1.5 expose structural uncertainty separately. No probability is assigned to these branches. Their envelopes are sensitivities, not a second statistical confidence interval.

Common calibration errors, regional spectral selection, meteorite representativeness and unmeasured mineralogy remain material uncertainties. The model does not independently fit selection bias in every orbital and size bin. Dense undiscovered populations use a compound-moment approximation; the worst estimated numerical standard error of a single-body water mean is 5.57%, recorded alongside the product. Thus a statement such as “95% of model draws contain enough water” is conditional on this model. It is not a 95% probability of a mine obtaining approval, meeting its delivery schedule or earning a return.

The [source audit](population_sources.json) also examines the historical [USGS asteroid-resource feasibility study](https://doi.org/10.3133/ofr20171041). Its printed program appears to multiply mass by numerical weight-percent grades without dividing by 100, and its carbonaceous-density expression conflicts with the accompanying interval. These are apparent defects in the printed code; the actual execution behind its published table has not been reconstructed. Its large published resource totals are preserved as historical source claims and are not used to calibrate this calculation. The independently implemented model uses explicit kilogram and fraction accounting.

## Fresh demand after recovery

The [provisioning experiment’s pinned export](inputs/reclamation.json.gz) recalculates the physical system for services, propellant/logistics and diversified portfolios, each under strong markets and either 0.5% or 1% annual real market growth. Six variants produce 36 cases. This bounded comparison covers representative orderbooks, rather than claiming a fresh rerun of all nine financial portfolios. The resource consumer also retains the original 81-case demand ranges for context. The producer and its full interpretation are `provisioning/enclosed_industry/reclamation.py` and `reclamation.md` at provisioning commit `210a20bd2cf16ec951f02a01f6ee3758b497df03`.

| Checkpoint | Original 81 cases: fresh source feed, Gt | Six matching controls: fresh source feed, Gt | Combined recovery/refurbishment: fresh source feed, Gt | Combined: source water product, Gt |
|---|---:|---:|---:|---:|
| Year 100 | 3.92–7.03 | 3.96–6.10 | 3.38–5.50 | 0.24–0.38 |
| Year 200 | 164.30–291.47 | 164.59–264.62 | 134.07–232.62 | 9.39–16.28 |
| Year 220 | 224.87–421.59 | 225.27–385.08 | 181.28–338.39 | 12.69–23.69 |

Ranges across orderbooks are scenario ranges, not statistical confidence intervals. Seed production is counted once. The improved-loop case raises domestic/industrial water recovery from 98% to 99.5%, crop recovery from 99.5% to 99.7%, fish-purge recovery from 95% to 98%, structural usable return from 95% to 98.01%, and household-goods return from 90% to 95.06%. These proposed coefficients exceed demonstrated system-wide lunar performance and therefore remain engineering hypotheses.

Refurbishment collects 80% of retiring eligible machinery and accepts 60% of collected units after inspection. Each accepted unit retains 80% of its mass and receives 20% new repair parts; discarded old parts remain recorded. Inspection delays availability by one year, followed by the existing two-year work period. Local workshops at source, depot and surface locations avoid assuming free movement of retired equipment. Successfully requalified units receive the original expected service life, including later refurbishment cycles. That lifetime reset and the purchase price of specialized spares are consequential qualification and cost assumptions.

The services comparison at 0.5% market growth illustrates the mechanisms. Through year 200, fresh source feed falls from 164.59 Gt to 136.48 Gt with improved loops, to 161.73 Gt with refurbishment alone, and to 134.07 Gt with both. In year 200, physical operating outlay falls from 4.128 trillion USD2041 to 3.029 trillion USD2041 in the combined case. Refurbishment saves much more money than its direct feed reduction suggests because it avoids purchasing expensive active equipment. Reduced freight then saves propellant and the additional mines, power systems and vehicles that would have supplied that propellant.

Transport still dominates the water account. In the combined services case at year 200, source-ferry fuel consumes 87.38 billion kg of water equivalent and depot demand consumes 60.27 billion kg, against 15.78 billion kg for residents and visitors. The two transport-related categories account for about 90% of source water output, with a small oxygen provision included in depot demand. Better household closure therefore helps, but the larger long-term opportunity also involves freight mass, equipment renewal and propulsion. Chemical exhaust remains an expended resource in this architecture.

An inherited allowance of 5 kg of water per person-day originally represented landscape and agricultural losses. Later models also carry explicit crop and fish losses, leaving possible overlap. Every primary variant retains that allowance. A separate maximum-credit endpoint removes it while preserving all explicit loop losses; its services feed requirement falls to 81.22 Gt through year 200. This result identifies a valuable accounting and engineering investigation. It does not establish that all of the allowance is redundant or recoverable.

### Which loops are represented

| Material or function | Current treatment | Remaining development |
|---|---|---|
| Domestic, industrial, crop and fish water | Separate processed flow, recovery and fresh loss/makeup; extra treatment hardware and energy in the improved case | Stream-specific contaminant limits, consumables, maintenance and measured achievable recovery |
| Household goods | Existing 90% return, compared with 95.06%; residual carbon feed and local material demand retained | Product lifetimes, polymer degradation and element-specific quality |
| Local bulk and structural metal | Dated retirements, collected material, usable returns, carried inventories and virgin remainder | Alloy separation, contamination, collection routes and measured remanufacturing yields |
| Selected active machinery | Inspection, accepted stock, retained parts, imported repair parts, rejected/discarded mass and delayed return to service | Repeated-life qualification, electronics and other unmodelled assets, spares premiums |
| Ore residues and rejected hardware | No automatic second recovery credit | Retained, assayed and traceable secondary stocks; costed reprocessing before any future credit |
| Food nutrients, carbon and metabolic gases | Existing makeup inputs remain; no complete linked elemental circulation account | Food, waste, biomass, atmosphere and water balances for C, N, P, K and trace nutrients; oxygen/redox accounting |
| Atmosphere leakage and construction stocks | Makeup and initial stocks remain physical requirements | Leakage control and stock recovery during decommissioning |
| Chemical rocket exhaust | Zero collection/recovery credit | Transport architectures that reduce expelled mass; any exhaust recovery would need an explicit capture design |

A loop's fresh demand is its initial stock increase plus irreversible loss, quality rejection and unrecovered withdrawals. Reusing material several times does not erase the first stock requirement. At a stable stock, with no exports of that material and no accumulating waste, a usable return fraction *r* leaves fresh replacement of *(1−r)* times gross turnover. For example, 98% usable return still requires 2% of each turnover, and the material may turn over many times in 200 years. This relationship explains why near-closure is valuable without turning it into perpetual zero-cost supply.

The current [phosphorus-return work](../phosphorus.md) supplies useful conservation principles, but its Open Moon reservoirs cannot be imported as measurements of the enclosed food system. A full enclosed nutrient account needs its own inventories, waste collection, biological uptake, contamination constraints and losses. Similarly, NASA's [oxygen-recovery research](https://www.nasa.gov/spacecraft-oxygen-recovery-scor/) addresses metabolic carbon dioxide, whereas the inherited physical air account represents leakage and new atmospheric stocks. Adding an oxygen-recovery credit to the latter would mix different flows. Complete elemental closure remains an explicit development requirement.

## Does the wider potential cover the dated demand?

The comparison applies recovery once to contained water, then tests 1%, 5%, 10%, 25%, 50% and 100% availability of each population pool. These fractions describe combinations of suitable composition, practical access, preservation, future reserves and material that can actually be allocated. They are planning conditions, not adopted extraction allowances or probabilities of permission. The 100% branch is a theoretical ceiling with no allowance for further exclusions or future-use reserves. Body-mass and contained-water tests are paired within each Monte Carlo draw; raw feed is not multiplied by a second ore-grade factor.

The table reports the fraction of paired draws that cover both source feed and water for the services portfolio at 0.5% market growth. A displayed 100% is a simulation frequency, not certainty; with 2,048 successful draws its 95% Wilson lower bound is approximately 99.8%, even before model uncertainty.

| Central population model; available share | Combined year 100 | Baseline year 200 | Combined year 200 | Combined year 220 |
|---|---:|---:|---:|---:|
| Wider screen; 10% | 100.0% | 0.3% | 0.6% | 0.2% |
| Wider screen; 25% | 100.0% | 7.6% | 17.4% | 5.2% |
| Wider screen; 50% | 100.0% | 72.2% | 92.0% | 60.4% |
| Wider screen; 100% ceiling | 100.0% | 100.0% | 100.0% | 100.0% |
| Broader NEAs without Don Quixote; 5% | 100.0% | 100.0% | 100.0% | 100.0% |
| Broader NEAs without Don Quixote; 10% | 100.0% | 100.0% | 100.0% | 100.0% |

For that combined services case at year 200, the median body-mass requirement is 7.5% of the wider stock; the median water requirement is 34.3% after 70% recovery. The paired required-share distribution reaches its 95th percentile at 53.8%. Across all six combined demand paths, this joint 95th-percentile availability condition ranges from 53.8% to 93.3%. These figures quantify a requirement on useful access, without predicting that the required fraction will be available. A source strategy can improve the fraction of water reached by targeting hydrated bodies rather than processing a representative mixture of all rock.

These are stock tests. All physical variants retain the inherited selective source grade of 10% water and 70% recovery. The central wider-pool mixture has a much lower average grade, so indiscriminate processing would need more feed and energy than the physical fleet is sized to handle. Selective prospecting can concentrate production in better sources, but geological stock sufficiency alone does not demonstrate the required distribution of mineable grades or an executable sequence of mines.

Under the central composition model, the wider pool's water-equivalent mass fraction is 1.41–3.84%, with a median of 2.22%. Processing that representative mixture would require 2.6–7.1 times the feed of a selected 10% source, with a median multiplier of 4.5, at unchanged recovery. The corresponding partial source-process cost is 23.8–64.8 USD2041 per kilogram of recovered water, with a median of 41.2, under the reference equipment/energy case. These intervals concern the pool-wide mixture, not the grade distribution of individual chosen mines.

The processing calculation uses the existing equipment and energy coefficients at the source gate. Freight, labour, shared infrastructure, site construction, native-condition archives and several preservation obligations remain outside those partial costs. Bound water can also require different chemistry. The grade sensitivity therefore identifies a risk to be solved by source selection and processing trials; it is not an integrated delivered-water price or a revised return estimate.

For a longer sustaining comparison, allocating 5% of the broader near-Earth pool without Don Quixote leaves a median of 127 additional stock-equivalent years after year 220 in the combined services case; the conditional 5th–95th percentile interval is 30–544 years. This divides remaining feed and water by their respective mean fresh demand in years 201–220 and retains the smaller span. It holds demand constant, credits no new discoveries or replenishment, and does not establish fleet capacity beyond the published horizon. Different population/grade branches and demand paths remain separate in the product.

## Lunar material distributions and depletion pressure

The [lunar calculation](lunar_potential.py) retains the measured covariance of six Lunar Prospector oxide estimates, samples regional variation by spherical area, and adds explicit density, recovery and process-cost priors. The density profile follows [Hayne et al.](https://doi.org/10.1002/2017JE005387), while its triangular parameter distributions are planning choices around the published thermal and sample constraints. Bulk inventory is calculated at 1, 5 and 10 metres as separate scenarios. Composition-derived inventories use the five-metre case only, with independent grade-depth transfer factors of 0.5, 1 and 1.5; no drilling observation is implied at those depths.

At a prospective footprint of 1,000 km² and a five-metre worked layer, estimated bulk mass is **7.84–9.21 Gt**, with a median of **8.66 Gt**. Across the original 81 cases, cumulative lunar raw demand is **3.70–8.19 Gt by year 200** and **3.79–10.59 Gt by year 220**. This puts bulk demand at the scale of hundreds to roughly a thousand square kilometres under that depth assumption, rather than suggesting Moon-wide bulk depletion. The mapped ±80° band alone has an estimated five-metre inventory of approximately 293,000–344,000 Gt; most of that area would never need to become a mine.

The global ratio cannot establish access to a particular metal. Targeting the mapped iron-rich quarter of the latitude band raises median FeO from 5.70% to 14.03%; targeting the aluminium-rich quarter raises median Al₂O₃ from 24.39% to 28.64%; targeting the titanium-rich tenth raises median TiO₂ from 0.51% to 4.12%. These overlapping prospecting pools are ranking scenarios at coarse regional resolution. Their areas are not independent stocks and must not be added.

With the declared recovery and cost priors, lunar feed handling has a 5th–95th percentile range of **0.16–0.35 USD2041/kg**, and inherited bulk processing gives **0.65–1.42 USD2041/kg of bulk product**. Handling alone per kilogram of titanium-oxide equivalent falls from a median of about 66 to 8.2 USD2041 when sampling the titanium-rich prospecting pool. This shows the value of choosing geology. It does not price beneficiation, oxide separation, reduction or a specified finished metal.

The resulting depletion assessment is asymmetric. The enclosed build is small against lunar bulk potential, while higher-value grades, nutrients and accessible asteroidal water can still become progressively harder to obtain. A defensible marginal-cost path requires a ranked sequence of sites, grade-depth continuity, yields, haul routes and preservation exclusions. Neither a constant cost forever nor an inevitable global shortage follows from the available distribution data. Existing local tailings, retired structures and machinery should enter that sequence as secondary sources before their material is counted as a fresh import.

## Preservation, evidence and use in the economic studies

The existing [physical and data-preservation protocols](protocols.md) apply to each candidate and remain unchanged. Survey must precede irreversible disturbance, including landing, exhaust and dust effects. Representative material must retain native conditions and an unopened fraction, with provenance, distributed archives and long-term access. Physical heritage and continuing natural functions cannot be replaced simply by recording data. Future accessible stocks and scarce-element recovery remain obligations; no universal fraction is approved here.

The four retained asteroid science holds are removed consistently from every relevant pool and from the known-versus-undiscovered calculation. The lunar latitude cut is only a geometric boundary: it does not establish a complete map of heritage, scientific sites or cold-trap records. Cargo still requires the project's remote braking, Earth exclusion, capture and other protection conditions. Unresolved qualification leaves the accessible fraction unknown; it does not turn estimated geological potential into zero.

The principal economic implication is that material abundance, fresh-material efficiency and delivery capability need separate accounts. The wider geological evidence supports a much larger opportunity than the narrow characterized sample. The refurbishment experiment supplies a concrete pathway to lower recurring expenditure, while selective-source grade and chemical freight remain major conditions on that saving. The financial studies can consume the new physical product in a subsequent explicit comparison; the present revision does not silently substitute new costs into their published returns.

## Products and reproduction

[Population inputs](population_inputs.json) and [source records](population_sources.json) retain the assumptions, source identities and reading scopes. [Lunar inputs](lunar_potential_inputs.json) and [source records](lunar_potential_sources.json) distinguish formal measurement covariance from engineering priors. The [comparison inputs](potential_inputs.json) pin these products and the committed reclamation export; [its source record](potential_sources.json) states the cross-domain reading rules. Annual reclamation audits run from years 0–220, with checkpoint comparisons at 100, 200 and 220. The frozen source catalogue and external model grid must be restored by their recorded hashes before rerunning the population producer.

```sh
PYTHONDONTWRITEBYTECODE=1 python -m resources.source_access.population
PYTHONDONTWRITEBYTECODE=1 python -m resources.source_access.lunar_potential
PYTHONDONTWRITEBYTECODE=1 python -m resources.source_access.potential
PYTHONDONTWRITEBYTECODE=1 python visualization/source-access/potential.py
```

The compressed outputs use schemas `terluna.resources.nea-population/1`, `terluna.resources.lunar-potential/1` and `terluna.resources.potential-and-reclamation/1`. They retain source/model hashes, paired draws or summaries, exact demand identities and explicit unknown qualification fields. The [verification record](potential_checks.md) separates accounting tests, source reading, numerical approximations and visual inspection from empirical validation.
