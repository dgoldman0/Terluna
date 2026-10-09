# Returning phosphorus through the aerial biosphere

The author's 9 October discussion requires a functioning phosphorus cycle, explicit aerial life and a common population
comparison. The new [calculation](phosphorus.py) and [saved result](results/phosphorus.json) establish the return rates
that selected productive skies would require. They also follow a finite inventory through six reservoirs. The numerical
cases are requirements and sensitivities; they establish neither an actual closed ecology nor a sustainable population.

## Stock and annual use are different quantities

The common example occupies a projected footprint equal to 1% of the **lunar surface area**, 379,323 km², with assumed
annual net production of 1,000 g carbon per square metre. This denominator is 4πR², not the Moon's solar cross-section.
Carbon is 45% of dry biomass and 5% of annual production is harvested. The example therefore produces 379.3 Mt carbon
or 842.9 Mt dry biomass annually. Its photosynthesis, water supply, lift, light competition and edible fraction remain
biological requirements to establish; this resources calculation does not grant them.

At the canonical molar C:P ratio of 106:1, the **mass** ratio is 41.104:1. That gives:

| Quantity | Elemental phosphorus |
|---|---:|
| P passing through annual new production | 9.228 Mt/year |
| Standing P at 1–10 kg dry biomass per m² of occupied footprint | 4.153–41.527 Mt |
| P carried out in the assumed 5% harvest | 0.4614 Mt/year |
| Stock divided by production throughput | 0.45–4.5 years |

Stock divided by throughput is an accounting turnover time, not organism life expectancy. Internal remobilization can
reuse the same phosphorus repeatedly. Conversely, respiration of carbon does not destroy phosphorus; excretion,
precipitation, dead bodies and harvest redistribute it. The calculation assumes that harvested and unharvested biomass
have the same C:P; P-rich seeds or P-poor storage tissues require their own compositions.

Redfield is a comparator. Martiny et al. (2013) observed marine C:P ratios from 78 to 195 among the regions they studied.
The saved sensitivity covers 78, 106, 137 and 195, plus a hypothetical 500:1 P-poor structural-biomass case. At 195:1 the
annual P throughput falls to 5.016 Mt; at 500:1 it is 1.956 Mt. Long-lived structural tissue, active plankton and food
products should ultimately have separate elemental budgets. Species richness cannot be derived from any of these totals.

## Aerial retention and food return must be counted separately

Let F be annual P throughput, h the harvested share, r the fraction of **unharvested** P recycled within the sky,
q the fraction of harvested P actually delivered back to aerial producers, and U other upward return from the surface.
The sky's annual balance is

\[
\Delta S = U-\{(1-h)(1-r)+h(1-q)\}F.
\]

The food-return fraction includes collection, treatment, recovery, transport and bioavailability. A treatment plant's
recovery efficiency alone is insufficient. The compound q of several imperfect stages is their product. Nutrients
returned to a terrestrial farm remain available to that farm; they have not yet returned to the aerial source biome.

The inherited deposition scenario is 0.0007–0.028 g P/m²/year. Applied across the Moon, with an assumed 10% captured
and made available to the aerial system, it supplies **2,655–106,211 tonnes/year**. These Earth-derived rates are
sensitivity inputs, not a lunar aerosol calculation. Capturing 10% of global flux with 1% projected biological coverage
also requires concentration, extended collectors, motion or repeated interception; that capture has not been demonstrated.

Mahowald et al. (2008) explicitly balance atmospheric sources against deposition. **A captured dust or biological
particle transfers phosphorus already on the Moon. It is not a fresh lunar import.** Their total and soluble P differ,
so availability belongs in the capture assumption. Tipping et al. (2014) also find that larger biological particles
commonly redeposit near their source. Neither study supplies an indefinitely replenished aerial nutrient source.

For the reference biomass, the minimum unharvested-P recycling required is:

| Harvest P delivered back, q | Low upward supply | High upward supply |
|---|---:|---:|
| 0% | Impossible: r would be 105.23% | Impossible: 104.05% |
| 90% | Impossible: 100.50% | 99.315% |
| 99% | Impossible: 100.022% | 98.841% |
| 99.9% | 99.97498% | 98.7938% |

Even perfect natural retention cannot support a 5% harvest without returning at least 99.425% of harvest P at the
low supply, or 76.982% at the high supply. With natural retention of 99.9% and harvest return of 99%, downward export
is 13,381 tonnes/year; the low supply leaves **10,726 tonnes/year of additional upward return** to provide.
The high supply covers that local balance, provided its donor reservoirs continue supplying it.

The full result also tests natural retention of 99%, 99.9% and 99.99% against every q above. The high efficiencies are
requirements, not demonstrated biology or industry. Egle et al.'s (2016) comparative study reports 60–90% wastewater-P
recovery for sludge-ash routes; it does not demonstrate a 99–99.9% collection-to-aerial-delivery chain.

## An explicit return architecture

| Route | Required process and account |
|---|---|
| Inside an aerial organism or community | Resorb nutrients from senescing tissue; retain wastes and dead material in wet microhabitats; decomposers release available phosphate. Count this in r. |
| Sky to people to sky | Collect harvest residues and excreta, recover phosphate, remove contaminants, formulate a biologically available product, and deliver it to the source aerial communities. Count the entire route in q. |
| Sky to land | Rain, settling and dead organisms enter soils; dust, biological dispersal, migrating animals and deliberate collectors can return P. Subtract every return from the terrestrial donor inventory. |
| Land to water | Runoff carries dissolved and particle-bound P into lakes and seas. Retention, river-mouth recovery and catchment return can intercept the flow; they compete with aquatic nutrient needs. |
| Water to sky | Spray and biological transport provide natural paths; concentrated recovery followed by controlled nutrient delivery is the engineered alternative. Fresh lunar seas invalidate unexamined seawater-spray analogues. |
| Water to sediment to circulation | Remineralization returns part of settled P; burial, adsorption and apatite formation retain part. Recovering that stock requires mapped sediment extraction, separation and return, with habitat impacts and contaminant handling. |
| Beyond the recovery system | Trace true escape, unrecoverable dispersal and isolated burial separately. Only these require new reserve extraction or off-world makeup for whole-system closure. |

Ruttenberg and Berner (1993) provide field evidence for sedimentary apatite formation. The Moon's sediment rates,
chemical fractions and extractable inventories remain unknown. A controlled transfer of concentrated nutrient material
to aerial organisms is a candidate for testing. Broadcasting soluble phosphate into the atmosphere would require its
own uptake, aerosol, health and ecological evaluation and is not an adopted delivery design.

## Finite reservoirs change the long-term answer

The six-box model tracks sky, land, sea, sediment, human return and material outside the recovery system. Each transfer
is removed from its donor and added to its recipient. Immediate internal recycling stays inside the sky box. There are
**no imports**. Production throughput scales with remaining sky stock, so a depleted sky cannot keep producing at its
initial rate. This linear model diagnoses nutrient limitation only; it has no photosynthetic ceiling or ecological feedback.

The demonstration starts with ten times the sky's P stock split 72:28 between land and sea; this is an assumed source
inventory, not a lunar geochemical assessment. The initial upward flux uses either end of the deposition scenario and
subsequently declines or rises with donor stock. Land drains to sea on a 100-year timescale and sea to sediment on a
1,000-year timescale. The tested engineered sediment return is 10,000 years. Natural aerial retention is 99.9%, food
return 99%, and human handling time 0.1 year. These deliberately exposed parameters are not fitted lunar rates.

For the 1 kg dry/m² standing-stock example:

| Treatment of sediment | Sky stock after 10,000 years, relative to initial | Long-term outcome in this model |
|---|---:|---|
| No sediment return | 0.00002–0.0698 | Accessible P accumulates in sediment; sky productivity tends to zero |
| Sediment returned, no irreversible loss | 0.3197–1.1523 | Finite equilibrium sky stock; asymptotic range 0.3197–1.1343 |
| Same return, 0.1% lost per processed upward/engineered pass | 0.3191–1.1430 | Loss accumulates; essentially all P reaches the unrecoverable box by 10⁹ years |

The wider initial-stock sensitivity, 1–10 kg dry/m², is in the machine-readable result. High initial upward flux can
increase the sky stock by drawing down land and sea; it does not create nutrients or establish that light can support
the extra stock. All boxes including the unrecoverable box conserve P. Matrix-exponential roundoff remains below
0.00003% of total P over the billion-year runs; splitting shorter steps gives the same trajectories within test tolerance.

A second calculation holds target production fixed and asks what reserve would supply a persistent return deficit.
At 99.9% natural retention, 99% food return and low upward supply, a **hypothetical 1 Gt elemental-P reserve lasts
about 93,200 years**; continuing that deficit for 10⁹ years would require about 10,726 Gt P. This is not a lunar reserve
estimate. It demonstrates why a good annual recycling percentage and a billion-year closed cycle are different tests.
Where the table gives zero deficit, its assumed upward source still needs its own long-term return path.

## Lift the nutrient concentrate, and account for its carrier

Returning the entire reference harvest's 461,415 tonnes P/year to 10 km costs the following **gravity-only** power at
an assumed 70% lifting efficiency:

| Delivered P mass fraction | Material lifted annually | Gravity-only power |
|---|---:|---:|
| 10%, concentrate | 4.61 Mt | 3.37 MW |
| 0.1%, slurry or dilute product | 461 Mt | 337 MW |
| 0.001%, approximately 10 mg P/L water | 46.1 Gt | 33.7 GW |

At 20 km the respective powers are 6.71 MW, 671 MW and 67.1 GW. The calculation uses the lunar gravitational potential
change, rather than Earth g or a constant-g extrapolation. It omits treatment, extraction, concentration, water handling,
pipe friction, drag, vehicles, buoyancy hardware, navigation and storm avoidance, maintenance, or any regenerative
lowering credit. These are transport lower bounds. Passive circulation can move nutrients but still needs a quantified source
of work; its existence does not remove the carrier mass. Concentration and contamination control may dominate a practical
system's costs even when lifting concentrated P is inexpensive.

## Population accounts use the shared cases

New joint work reads [shared/scenarios/population.json](../shared/scenarios/population.json). No allocation is selected.
The independent historical scenarios in [industry.md](industry.md) remain reproducible and are superseded for joint
planning. Applying their generic 43–335 t/person material-stock coefficients gives:

| Shared case | Lunar surface + aerial residents, billions | Generic lunar material proxy, Gt | Orbital residents, billions | Orbital shield + structure comparator, Gt |
|---|---:|---:|---:|---:|
| total20 | 3 + 3 | 258–2,010 | 2 | 1,980 + 30 |
| total25 | 5 + 5 | 430–3,350 | 3 | 2,970 + 45 |
| total30 | 7.5 + 7.5 | 645–5,025 | 3 | 2,970 + 45 |
| total25_earth_heavier | 3 + 3 | 258–2,010 | 3 | 2,970 + 45 |
| total25_moon_heavier | 6 + 6 | 516–4,020 | 3 | 2,970 + 45 |

The lunar quantities are economy-wide Earth-derived material proxies, **not supported aerial payloads**. Aerial district
mass, lifting gas, envelopes and flight structures must be sized separately in [provisioning's inhabited-volume study](https://github.com/dgoldman0/Terluna/blob/domain/provisioning/research/studies/inhabited_volume/README.md);
its district inventories overlap these proxies and must not simply be added to them. The historical stock/lifetime
combination implies generic renewal of 3.19–77.31 Gt/year across these cases, without establishing actual lunar service
lives. Material grades, energy and recovery rates still require individual accounts.

The orbital comparator explicitly retains SP-413's 990 t shield and 15 t structure per person; it is not an optimized
future habitat design and does not apply to lunar atmospheric residents. Earth and elsewhere infrastructure are outside
this screen. P demand for all residents' complete diets remains a provisioning input; the aerial 5% harvest is one
conditional contribution, not permission to multiply the wild ecosystem's yield by population.

## Reproduction and what remains

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m resources.phosphorus`, then
`python -m pytest resources/tests/test_phosphorus.py`. The product binds its producer, shared scenario, historical industry
product and [source register](phosphorus_sources.json), with named shared constants and explicit reading rules.

The next discriminating work is a spatial P-inventory and return model coupled to aerial production: dissolved versus
particle-bound P, organism-specific C:P and tissue turnover, droplet/particle capture and bioavailability, rainfall loss,
food-return contamination and retention, sediment speciation and extraction costs. The source register distinguishes
full-text passages, abstracts and publisher search excerpts. Full source admission, engineering demonstrations and
biological validation remain separate from these conservation tests.

Checks on 9 October: all 24 resource tests passed, including eight phosphorus tests, and an independent review found
no numerical or conservation defect in the saved cases. `make -k check` ran before commit: 1,242 Python tests passed,
111 skipped and 30 failed. Twenty-nine failures require ignored external inputs/products absent from this checkout
(WHI/FISM2 spectra, GCM configuration, Earthshine, water/soil optics or geography arrays); the remaining pre-existing
process-lookup test fails in the sandbox. Layer checks, the provenance/baseline checker and the ensemble validator
passed. Calendar tests passed 13 with one skipped because the spherical optical atlas is absent. The immersion bake
succeeded, but Node 24 rejects the Makefile's directory test argument; explicit test-file diagnostics passed 120 and
failed five in immersion with its `three` dependency absent. These repository-wide issues remain separate from the
new resource results; full `make check` is not green in this checkout.
