# TERLUNA COMPLETE KNOWLEDGE DUMP

This file concatenates the narrative bundle for one-file import. Individual source files remain authoritative for editing.


---

<!-- BEGIN FILE: README.md -->

# Terluna Moon Terraforming Knowledge Bundle

Version 1.0 — 2026-07-12

## Purpose

This ZIP is a portable reconstruction and engineering expansion of the **Terluna** Moon-terraforming concept developed in the Moon World project. It is designed for import into an external AI workspace, wiki, notebook, research repository or knowledge base.

The accessible project context did not include a verbatim export of every earlier message. The bundle therefore distinguishes:

- **Recovered canon** — details clearly preserved from accessible project context.
- **Reconstructed baseline** — quantitative values introduced to make the concept modelable.
- **Engineering inference** — consequences calculated from lunar physics.
- **Speculative worldbuilding** — plausible extensions that require research.

## Fastest import paths

- Import `TERLUNA_COMPLETE_DUMP.md` when the external system accepts one large document.
- Import `13_EXTERNAL_WORKSPACE_CONTEXT.md` for a compact ready-to-paste context.
- Import `data/terluna_context.json` for structured machine context.
- Import `data/qa_handoff.jsonl` for retrieval-style question and answer chunks.
- Import all numbered Markdown files when the workspace supports folders.

## Core concept

Terluna is a deliberately transformed version of Earth’s Moon: a low-gravity world with a globally collisional breathable atmosphere, regional seas, manufactured soils and a managed biosphere adapted to a 29.53-Earth-day solar cycle. The roughly five-hour practical dawn or dusk transition is canonical and lies inside the natural cycle. The concept requires a mature interplanetary civilization able to move and process order-10^18 kg volatiles, sustain enormous power and industry, actively maintain the atmosphere, and resolve major uncertainties in escape, climate, fire, ecology, human partial-gravity health and governance.

## Reading order

1. `00_EXECUTIVE_SUMMARY.md`
2. `TERLUNA_MASTER_SPEC.md`
3. `01_RECOVERED_CANON_AND_CORRECTIONS.md`
4. `02_FEASIBILITY_CASE.md`
5. `03_PHYSICAL_PARAMETERS_AND_CALCULATIONS.md`
6. Domain files `04` through `10`
7. `11_RESEARCH_AND_VALIDATION_PROGRAM.md`
8. `12_QA_HANDOFF.md`
9. `13_EXTERNAL_WORKSPACE_CONTEXT.md`
10. `14_GLOSSARY.md` and `15_SOURCE_NOTES.md`

## Narrative files

- `00_EXECUTIVE_SUMMARY.md` — compact integrated overview.
- `TERLUNA_MASTER_SPEC.md` — principal reference specification.
- `01_RECOVERED_CANON_AND_CORRECTIONS.md` — protected canon, corrected twilight and terminology.
- `02_FEASIBILITY_CASE.md` — why the concept is not automatically forbidden and why it remains extremely difficult.
- `03_PHYSICAL_PARAMETERS_AND_CALCULATIONS.md` — atmosphere, water, energy, flight and artificial-gravity calculations.
- `04_ATMOSPHERE_CLIMATE_AND_WEATHER.md` — composition, circulation, weather, twilight, aviation and escape.
- `05_HYDROSPHERE_GEOSPHERE_AND_SOIL.md` — water, rivers, basins, regolith, soil, dust and nutrients.
- `06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md` — succession, long-cycle biology, aerial ecology and biosafety.
- `07_HUMAN_SETTLEMENT_HEALTH_AND_CULTURE.md` — artificial gravity, architecture, time, food, medicine and culture.
- `08_INFRASTRUCTURE_INDUSTRY_AND_LOGISTICS.md` — volatiles, mining, power, transport and space access.
- `09_IMPLEMENTATION_ROADMAP.md` — staged pressure and ecology program with stop gates.
- `10_RISKS_ETHICS_GOVERNANCE_AND_FAILURES.md` — risk, irreversible loss, institutions and failure responses.
- `11_RESEARCH_AND_VALIDATION_PROGRAM.md` — models, experiments and acceptance metrics.
- `12_QA_HANDOFF.md` — 85 self-contained one- or two-sentence answers.
- `13_EXTERNAL_WORKSPACE_CONTEXT.md` — ready-to-paste context and suggested workspace structure.
- `14_GLOSSARY.md` — terminology registry.
- `15_SOURCE_NOTES.md` — source-by-source interpretation limits.
- `TERLUNA_COMPLETE_DUMP.md` — concatenated single-file version of the narrative corpus.

## Structured data

The `data/` directory contains canon, assumptions, ranges, mass and energy budgets, atmospheric bands, biosphere traits, requirements, risks, dependencies, decisions, work packages, roadmap gates, open questions, sources, glossary data, Q&A JSONL and the central context JSON.

## Figures and models

The `figures/` directory contains atmosphere-mass, atmospheric-band, illumination-cycle and roadmap diagrams. The `models/` directory contains reproducible screening calculations and a limitations note.

## Reference design status

The 80 kPa atmosphere, 24% dry oxygen, 285 K temperature target and 3 × 10^17 kg hydrosphere are **editable reference values**, not recovered canon. The bundle preserves alternative ranges and identifies the evidence required to change them.

## Scientific stance

The bundle does not claim present feasibility. Its narrower conclusion is that a maintained lunar atmosphere is not ruled out by a simple law of physics, while the loss lifetime of a dense modern atmosphere, nitrogen sourcing, long-cycle climate, fire, low-gravity biology and legitimacy remain hard gates.

## Citation convention

Narrative files cite `[S01]` through `[S20]`. Full bibliographic notes and URLs are in `15_SOURCE_NOTES.md`, `data/source_manifest.csv` and `references/references.bib`.

## Integrity

`BUNDLE_MANIFEST.json` lists file sizes and SHA-256 hashes. `SHA256SUMS.txt` provides line-oriented checksums for verification.

<!-- END FILE: README.md -->


---

<!-- BEGIN FILE: 00_EXECUTIVE_SUMMARY.md -->

# Terluna Executive Summary

## Concept

Terluna is a far-future transformation of Earth’s Moon into a maintained low-gravity open world with breathable air, regional seas, manufactured soil, engineered ecosystems, permanent settlements and extensive atmospheric flight. It keeps the Moon’s natural orbit, 0.16 g gravity, synchronous rotation and 29.53-Earth-day solar cycle.

The landscape remains visibly lunar: maria, highlands, craters, scarps and polar basins shape seas, climates and cities. Nearside and farside retain distinct identities. Biology and culture adapt to approximately two weeks of illumination followed by two weeks of darkness.

## Recovered canon

The accessible project context establishes:

- The name **Terluna**.
- A long lunar light-dark cycle as a defining feature.
- Biological resource storage and adaptation across that cycle.
- Extensive aerial or gliding life.
- A canonical practical dawn or dusk transition of about **five hours**, replacing an earlier incorrect 30-hour value.
- A provisional upper-troposphere/lower-stratosphere layer around **35–45 km**.
- Atmospheric gravity waves from terrain and terminator fronts.
- Sparse high-altitude microbes or engineered radiation-hard algal films.
- High-altitude platforms for farside astronomy and nearside Earthlight photometry.
- Descriptive atmospheric layer names rather than L1–L5, which remain reserved for Lagrange points.

## Reference design introduced by this bundle

For quantitative analysis, the bundle uses an editable reference case:

| Parameter | Reference value |
|---|---:|
| Surface pressure | 80 kPa |
| Dry atmosphere | 75% N2, 24% O2, 0.9% Ar, 800 ppm CO2, 200 ppm other trace gases |
| Oxygen partial pressure | 19.2 kPa |
| Mean surface-temperature target | 285 K |
| Hydrosphere | 3 × 10^17 kg |
| Practical twilight | 5 h per dawn or dusk transition |
| Main tropopause region | 35–45 km provisional |

These values are modeling baselines, not recovered canon.

## Quantitative scale

- An 80 kPa atmosphere requires about **1.87 × 10^18 kg** of gas.
- Nitrogen accounts for about **1.35 × 10^18 kg**, the largest material bottleneck.
- Oxygen accounts for about **4.93 × 10^17 kg**.
- Producing that oxygen at 45% regolith oxygen content and 75% recovery requires processing about **1.46 × 10^18 kg of regolith**.
- At an illustrative 20–50 MJ/kg O2, oxygen production alone averages about **313–782 TW over 1,000 years**.
- The atmospheric scale height is about **50 km** at 285 K, producing a deep weather and orbital-drag environment.
- The 80 kPa column mass is about **49,300 kg/m²**, potentially providing strong radiation and small-meteoroid shielding.
- The reference water inventory gives a **7.9 m global equivalent depth** or about **79 m over 10%** of the surface.
- The equatorial sunrise/sunset line moves at about **4.3 m/s**.

## Why the concept might work

No established physical law forbids a maintained collisional atmosphere on the Moon. Ancient volcanic outgassing probably produced a transient lunar atmosphere, and modern escape studies show that lifetime depends strongly on upper-atmosphere temperature and solar forcing [S05][S06].

The Moon receives Earthlike solar flux. A dense atmosphere and substantial regional hydrosphere have heat capacities comparable in scale to the energy absorbed during a half-cycle, offering a plausible route to surviving the long night if circulation and clouds behave favorably. Lunar rock supplies oxygen, metals and construction mass, while low gravity favors aviation and large structures.

## Why the concept might fail

The decisive unknown is dense-atmosphere loss. If pressure requires rapid continuous replacement, the project is untenable regardless of lower-atmosphere habitability.

Independent hard problems remain:

- Sourcing roughly 10^18 kg of nitrogen without damaging Earth or another protected world.
- Importing water, hydrogen, carbon and nutrient volatiles.
- Maintaining habitable regional climates through the long day and night.
- Preventing global fire and oxygen-cycle instability.
- Converting toxic, abrasive regolith into soil.
- Protecting humans across generations at 0.16 g.
- Replacing current low lunar orbit architecture.
- Preserving lunar science, heritage and international legitimacy.

## Climate and weather

Terluna’s slow rotation favors broad overturning circulation, weak familiar Coriolis organization and large weather systems. The moving terminator can generate dawn and dusk fronts, fog, frost, thawing and pressure waves. Regional seas provide thermal storage but can support persistent storms.

Low gravity produces a deep convective atmosphere. Rain and snow microphysics, cloud depth, superrotation and 35–45 km wind shear remain research questions rather than settled canon.

## Biosphere

The biosphere develops by controlled succession:

1. Geochemical conditioning.
2. Microbes and biofilms.
3. Soil crusts, lichens and pioneer plants.
4. Fungi, wetlands, grasslands and shrubs.
5. Forest and aquatic complexity.
6. Pollinators and small animals.
7. Large animals only after long-term stability.

Plants use stored starch, oils, bulbs, rhizomes, dormancy, leaf folding and adjustable pigments. Animals use fat, torpor, caches, burrows, communal roosts and migration. Deep water and under-ice systems act as night refugia. Dense air and low weight create major niches for large soaring, gliding, ballooning and terminator-following organisms.

## Human settlement

Breathable air does not solve lunar gravity. Cities require rotating sleep habitats, clinics, schools, maternity facilities and exercise systems until multigenerational evidence establishes a safe partial-gravity dose.

Buildings must be anchored against winds because they weigh one sixth as much while aerodynamic loads remain substantial. Residents keep approximately 24-hour civil schedules inside a 29.53-day environmental month. Every city retains sealed emergency air, water, power and food reserves.

## Implementation strategy

Terluna advances through stages:

0. Global baseline, archives, protected zones and legitimate authority.
1. Robotic power, mining, oxygen and transport industry.
2. Closed ecological cities and artificial-gravity research.
3. Paraterraformed enclosed basins.
4. Thin global atmosphere for escape and climate experiments.
5. Pressure plateaus at approximately 10, 20, 40 and 60 kPa.
6. Hydrosphere and soil expansion under a mostly inert atmosphere.
7. Controlled oxygenation.
8. Managed open biosphere.
9. Permanent planetary stewardship.

Each stage has a stop gate. Paraterraforming remains a successful fallback and may prove preferable even if full terraforming is technically possible.

## Final verdict

Terluna is scientifically useful and internally coherent as far-future planetary engineering when its resource scale and unknowns remain explicit. It is not a current construction proposal and should never be presented as proven feasible. The bundle’s purpose is to preserve the idea, expose its dependencies, and give an external workspace enough structure to extend it without losing canon or mistaking assumptions for facts.

<!-- END FILE: 00_EXECUTIVE_SUMMARY.md -->


---

<!-- BEGIN FILE: TERLUNA_MASTER_SPEC.md -->

# Terluna Master Specification

## 1. Definition

**Terluna** is a hypothetical transformation of Earth’s Moon into a maintained open world with breathable air, liquid water, biologically active soils, engineered ecosystems, permanent human settlements, and extensive atmospheric flight.

The project keeps the Moon’s existing orbit, mass, gravity, synchronous rotation, and long solar cycle. It does not assume magical gravity alteration, a planetary shell, or a forced 24-hour rotation. Technology is allowed to be extremely advanced, but every major subsystem must close its mass, energy, thermal, ecological, and maintenance budgets.

## 2. Design intent

Terluna should feel recognizably lunar after terraforming. The world retains:

- 0.16 g surface gravity.
- Nearside and farside identity.
- A 29.53-Earth-day sunrise-to-sunrise cycle.
- Roughly 14.765 days from sunrise to sunset and 14.765 days from sunset to sunrise near the equator.
- A bright operational dawn and dusk transition of about five hours inside each corresponding half-cycle, not added to it.
- Craters, maria, highlands, scarps, massifs, and polar basins as dominant geography.
- Earth fixed in the nearside sky apart from libration.
- Weak seasons because the lunar axis is only slightly tilted to the Sun.

The biosphere, economy, and culture adapt to those conditions instead of erasing them.

## 3. Feasibility verdict

The concept is not achievable by current civilization. It is a planetary engineering program for a mature interplanetary economy.

No established physical law makes a maintained atmosphere impossible. Evidence and modeling indicate that the ancient Moon could sustain a transient atmosphere around 1 kPa and that atmospheric loss depends strongly on composition, solar forcing, and upper-atmosphere temperature [S05][S06]. That precedent establishes only that lunar atmospheres can exist. It does not prove that an 80 kPa nitrogen–oxygen atmosphere would remain economical to maintain.

The decisive question is whether the dense atmosphere’s loss lifetime is long enough that replenishment becomes routine rather than civilization-consuming. No global release occurs until a coupled lower-atmosphere, thermosphere, exosphere, photochemistry, solar-wind, and Earth-magnetotail model closes that question.

## 4. Reference environmental envelope

These values are a reconstructed modeling baseline.

| Parameter | Reference | Design range | Status |
|---|---:|---:|---|
| Surface pressure at datum | 80 kPa | 60–100 kPa | Reconstructed baseline |
| Oxygen partial pressure | 19.2 kPa | 18–22 kPa | Engineering target |
| Dry O2 fraction | 24% | 20–30% | Must be fire-tested |
| Dry N2 fraction | 75% | 65–79% | Primary buffer gas |
| Dry CO2 fraction | 800 ppm | 300–3,000 ppm | Actively controlled |
| Mean surface temperature | 285 K | 275–295 K | GCM target, not prediction |
| Reference water inventory | 3 × 10^17 kg | 10^17–10^18 kg | Regional hydrosphere |
| Practical twilight | 5 h | Canonical value | Conversation-derived |
| Main tropopause region | 35–45 km | Model-dependent | Project-derived range |
| Human gravity exposure | 0.16 g plus artificial gravity | Unknown dose | Critical health requirement |

At 80 kPa, the atmosphere has a screening mass of about **1.87 × 10^18 kg**. Its dry reference composition contains approximately:

- 1.35 × 10^18 kg nitrogen.
- 4.93 × 10^17 kg oxygen.
- 2.31 × 10^16 kg argon.
- 2.26 × 10^15 kg carbon dioxide.

This nitrogen requirement is the single largest material problem.

## 5. Why 80 kPa is a useful reference

An 80 kPa atmosphere can supply nearly terrestrial oxygen partial pressure with a 24% oxygen fraction. It reduces total gas mass relative to sea-level Earth while avoiding the extreme oxygen fractions required by very low-pressure open-air designs.

The choice is not settled. A lower pressure saves hundreds of quadrillions of kilograms of gas but increases fire risk, dehydration, pressure sensitivity, weathering differences, and the fraction of oxygen needed. A higher pressure increases resource and escape burdens. The correct value requires integrated fire, biology, climate, and logistics modeling.

## 6. Atmospheric geometry

Lunar gravity gives an Earthlike-composition atmosphere at 285 K a first-order scale height near **50 km**, about six times Earth’s. The dry adiabatic lapse rate is only about **1.62 K per kilometer**. Consequences include:

- A physically deep troposphere.
- Clouds and convection extending much higher than on Earth.
- Significant air pressure at 35–45 km.
- Atmospheric drag extending far above current low lunar orbits.
- A large upper-atmospheric volume exposed to solar ultraviolet radiation and plasma loss.

At 80 kPa, the atmospheric column mass is about 49,300 kg/m², roughly 4.8 times Earth’s sea-level atmospheric mass column. Once established, that column would provide major radiation and meteoroid shielding, although detailed secondary-particle transport must be modeled.

## 7. Climate architecture

The average solar energy available per square meter is similar to Earth’s because the Moon shares Earth’s orbit around the Sun. The difficulty lies in timing and transport: the surface receives sunlight continuously for about two weeks and then none for about two weeks.

The reference climate uses five coupled buffers:

1. **Atmospheric heat transport** — broad day-to-night overturning flow.
2. **Regional seas and deep lakes** — high thermal inertia and night refugia.
3. **Wet soils, aquifers, ice and engineered thermal reservoirs** — distributed heat storage.
4. **Cloud and albedo management** — limit day-side overheating.
5. **Trace-gas control and optional orbital energy support** — prevent night collapse.

Slowly rotating terrestrial climate studies show that broad overturning cells and strong cloud feedbacks can moderate large-scale temperature contrasts, but those studies do not reproduce lunar gravity, topography, or the 29.5-day moving day-night pattern [S11][S12][S13]. A dedicated Terluna GCM is mandatory.

## 8. Weather regime

The world is expected to have weaker Coriolis organization than Earth, broader cells, deep vertical motion, and a moving terminator circulation.

Named atmospheric bands replace ambiguous `L1–L5` labels:

| Band | Approximate altitude | Function and conditions |
|---|---:|---|
| Surface boundary layer | 0–3 km | Terrain winds, dust, fog, convection, settlement weather |
| Lower troposphere | 3–12 km | Most clouds, regional aviation, rainfall systems |
| Middle troposphere | 12–25 km | Long-range flight, broad return flow, layered clouds |
| Upper troposphere | 25–35 km | Deep convection outflow, storm avoidance routes |
| Tropopause/lower stratosphere | 35–45 km | Weak-shear windows, atmospheric gravity waves, research platforms |
| Extended atmosphere | 45 km upward | Stratosphere, mesosphere, thermosphere, escape interface; not an ordinary flight layer |

“Gravity waves” in this context means atmospheric buoyancy waves launched by mountains, convection, and terminator fronts. It does not mean gravitational waves in spacetime.

## 9. Hydrosphere

The reference inventory of 3 × 10^17 kg of water corresponds to a global equivalent depth of about 7.9 m. Concentrated over 10% of the surface, it yields an average depth near 79 m. This supports regional seas and lake districts rather than a deep global ocean.

Water is placed selectively in mare and impact basins after geodesy, heritage protection, and crustal-loading analysis. Deep reservoirs serve as night thermal batteries. Some seas freeze seasonally at their edges while retaining liquid cores beneath ice.

Lunar polar ice and sunlit hydrated material prove that indigenous water exists [S03][S07][S08]. Known lunar water is a strategic seed resource, not a demonstrated source for a global hydrosphere. Most Terluna water must be imported or manufactured from imported hydrogen and lunar oxygen.

## 10. Geosphere and soil

Raw lunar regolith is not soil. It is sharp, abrasive, organic-free, chemically unweathered material [S16]. Plants have germinated in Apollo regolith, but growth was slow and highly stressed [S09].

Terluna soil creation therefore requires:

- Dust capture and particle-size control.
- Washing and leaching.
- Oxidation and controlled weathering.
- Addition of nitrogen, carbon, phosphorus, sulfur and trace nutrients.
- Biochar and organic matter.
- Microbial and fungal inoculation.
- Repeated crop and fallow cycles.
- Containment of toxic metal mobilization.

Large areas remain mineral desert for centuries while soil islands expand outward.

## 11. Biosphere design

The biosphere is assembled as a managed succession:

1. Sterile physical weathering.
2. Microbes and biofilms.
3. Lichens, moss analogues and soil crusts.
4. Fungi and detrital networks.
5. Shrubs, grasses, wetland plants and engineered crops.
6. Pollinators and decomposer invertebrates.
7. Small vertebrates or engineered analogues.
8. Large animals only after centuries of stable nutrient cycles.

Life adapts to two clocks:

- A roughly 24-hour cellular and behavioral rhythm retained for physiology.
- A 708.7-hour environmental rhythm governing growth, dormancy, migration and reproduction.

Plants use storage organs, oils, starches, antifreeze compounds, adjustable pigments, leaf folding, deep roots and dormancy. Animals use fat storage, torpor, burrows, migration, communal heat retention and cached food. Ecosystems shift from high daytime productivity to a long night dominated by respiration, detritivory, fungi, stored energy and aquatic refugia.

## 12. Aerial ecology

Low gravity reduces weight to one sixth while atmospheric density at a given pressure and temperature remains comparable to terrestrial air. This creates a large niche for flight and gliding.

Plausible organisms include:

- Very large soaring animals with low wing loading.
- Cliff- and canopy-launching gliders.
- Aerial plankton and spore clouds in stable layers.
- Migratory fliers tracking the terminator or weather corridors.
- Balloon-assisted organisms using warmed gas or carefully contained light gases.
- Engineered photosynthetic films on research aerostats.

Aerial life must store water, food, nitrogen, minerals and reproductive material. It cannot assume frequent landing during long crossings. Hydrogen or methane lift cells introduce ignition hazards in an oxygenated atmosphere, so such organisms require fire-resistant compartmentalization or alternative buoyancy strategies.

## 13. Human settlement

A breathable atmosphere removes spacesuits from ordinary outdoor life and provides radiation shielding, but it does not solve low gravity.

Long-term health, pregnancy and childhood at 0.16 g remain unknown. NASA identifies altered gravity as a major human-system hazard, and partial-gravity experiments already show cellular changes [S14][S15]. Terluna settlement therefore assumes:

- Rotating residential districts or daily artificial-gravity facilities.
- 1 g medical, pregnancy and childhood environments until evidence supports otherwise.
- Resistance exercise and cardiovascular countermeasures.
- Continuous bone, vision, immune and developmental monitoring.

Cities favor crater walls, lava tubes, basin shores and mountain terraces. Structures can be tall because weight loads are low, but wind, fire, pressure, corrosion and seismic design remain important.

## 14. Industry and resources

### 14.1 Oxygen

Lunar regolith is about 45% oxygen by mass, chemically bound in minerals, and oxygen-extraction processes also yield useful metals [S04]. Producing the reference atmospheric oxygen by regolith processing at 75% recovery would require processing roughly 1.5 × 10^18 kg of regolith, equivalent to about 25 m of regolith averaged over the whole Moon. Concentrated mining reduces the affected area but creates kilometer-scale industrial excavation zones.

### 14.2 Nitrogen

The Moon lacks a demonstrated accessible nitrogen reservoir remotely close to the 10^18 kg requirement. Nitrogen must be imported as N2, ammonia, nitriles, hydrated minerals, or mixed outer-system volatiles. Removing it from Earth would be environmentally unacceptable at this scale. Venus and Titan contain enough nitrogen, but harvesting another world’s atmosphere introduces enormous technical and ethical problems. Small-body resources are the preferred conceptual source, subject to a real inventory.

### 14.3 Water and hydrogen

Water comes from polar ice, imported ice, hydrated asteroids, ammonia cracking, and chemical combination of imported hydrogen with lunar oxygen. Indigenous water is conserved as scientific evidence before large-scale use.

### 14.4 Power

Extracting around 5 × 10^17 kg of oxygen at an illustrative 20–50 MJ/kg requires roughly 10^25 joules. Spread over a millennium, oxygen production alone averages hundreds of terawatts. Transport, mining, nitrogen acquisition, water, construction and losses add substantially.

Terluna needs a mixed system of surface solar, polar solar, fission, fusion if available, orbital solar power, thermal storage, chemical storage and continent-scale grids.

## 15. Space access consequences

A deep lunar atmosphere eliminates the current concept of 50–100 km “low lunar orbit.” First-order pressure remains significant at those altitudes. Existing orbital infrastructure must migrate before atmosphere build-up.

Launch architecture shifts toward:

- High-altitude airports and aerostat ports.
- Mountain mass drivers.
- Electromagnetic launch tracks.
- Tethers or lunar elevators outside the densest atmosphere.
- High circular orbits and Earth–Moon Lagrange logistics.
- Carefully controlled rocket corridors.

The labels L1 through L5 should be used only for actual Lagrange points.

## 16. Implementation sequence

1. **Baseline and consent** — map, archive and legally protect the pre-terraforming Moon.
2. **Robotic industry** — power, mining, oxygen, metals and sealed settlements.
3. **Closed ecological pilots** — run full lunar-cycle habitats for decades.
4. **Thin global atmosphere** — several kilopascals, initially unbreathable, for dust suppression and climate experiments.
5. **Hydrosphere and soil** — fill selected basins and establish nutrient cycles.
6. **Pressure plateau tests** — 10, 20, 40 and 60 kPa with escape verification at each stage.
7. **Oxygenation** — add oxygen only after fire, oxidation and ecology controls mature.
8. **Open biosphere** — staged release with large no-release reserves.
9. **Permanent stewardship** — atmosphere and ecology remain actively managed.

Every stage has a stop condition. The project can remain a network of enclosed biospheres if global atmosphere economics fail.

## 17. Critical no-go gates

The program does not proceed to a breathable global atmosphere unless all of the following are demonstrated:

- Atmospheric replenishment is affordable for at least millennial planning horizons.
- No unacceptable effect on Earth, lunar orbit, or protected external activities is identified.
- Nitrogen supply is secured without destructive appropriation of another inhabited or scientifically protected world.
- A 3-D GCM maintains habitable refugia through repeated lunar cycles.
- Fire behavior is controllable.
- 0.16 g health is mitigated by proven artificial-gravity architecture.
- International consent, liability, preservation and rights systems exist.
- Reversible pilots have operated for multiple decades without hidden ecological collapse.

## 18. Overall conclusion

Terluna works best as a layered idea rather than a single heroic act. The Moon first becomes industrially inhabited, then ecologically enclosed, then climatically modified, and only later potentially opened. Its strongest scientific advantages are Earthlike solar input, abundant mineral oxygen, low-gravity construction and flight, and proximity to Earth. Its strongest objections are volatile mass, atmospheric escape, long-night climate, low-gravity biology, destruction of the pristine lunar environment, and planetary-scale governance.

The concept is credible as disciplined far-future engineering and worldbuilding when those objections remain visible and quantified.

<!-- END FILE: TERLUNA_MASTER_SPEC.md -->


---

<!-- BEGIN FILE: 01_RECOVERED_CANON_AND_CORRECTIONS.md -->

# Recovered Canon and Corrections

## Purpose

This file distinguishes the Terluna details recovered from accessible Moon World project context from the numerical architecture introduced for this bundle.

## Recovered project canon

### Identity

- The terraformed Moon is called **Terluna**.
- The target is Earth’s Moon, not a constructed duplicate or relocated body.
- The project imagines a world with an atmosphere, weather, biology, atmospheric flight and human use.

### Illumination

- Terluna retains the Moon’s long day–night cycle.
- Organisms must adapt to prolonged daylight and prolonged darkness.
- Biological design includes storage of food, water, chemical energy and other resources across the long cycle.
- Practical twilight is about **five hours**.

### Atmosphere and flight

- The atmosphere is vertically structured for practical aviation and aeroship routing.
- The top of the main troposphere or lower-stratosphere transition was placed around **35–45 km** in the project discussion.
- Upper atmospheric flight can encounter weak shear windows and atmospheric gravity waves launched by mountain chains and terminator fronts.
- High-altitude life is sparse, dominated by microbes or engineered films.
- Long-duration high-altitude platforms can support farside astronomy and nearside Earthlight photometry.

### Biology

- Biological notes should cover adaptation to aerial life.
- Biological notes should cover the long day–night cycle.
- Organisms require mechanisms for storing resources over long intervals.
- The world should be treated as a functioning hypothetical biosphere rather than only a collection of settlements.

## Corrections that supersede earlier wording

### 1. Five-hour twilight

Use **about five hours** for the canonical practical twilight transition. Do not repeat a thirty-hour figure.

A scientifically complete model can distinguish operational, civil, nautical and astronomical twilight. Faint skyglow might extend beyond the canonical interval in a deep lunar atmosphere. The project-facing value remains five hours unless deliberately revised.

### 2. Atmospheric layers must not be called L1–L5

The labels `L1` through `L5` are strongly associated with Lagrange points. Reusing them for atmospheric flight layers creates avoidable confusion.

Use named bands:

- Surface boundary layer.
- Lower troposphere.
- Middle troposphere.
- Upper troposphere.
- Tropopause/lower stratosphere.
- Extended atmosphere.

Use L1–L5 only when discussing gravitational equilibrium regions in an orbital system.

### 3. “Gravity waves” are atmospheric waves

In meteorology, a gravity wave is a buoyancy-restored wave in a stably stratified atmosphere. Mountain ranges, convection and moving fronts can launch them.

To prevent confusion with gravitational waves in spacetime, use **atmospheric gravity wave**, **buoyancy wave**, or **mountain wave** on first reference.

## Reconstructed choices introduced in this bundle

The following are not known to have been fixed in the accessible project discussion:

- 80 kPa reference surface pressure.
- 24% dry oxygen mole fraction.
- 75% dry nitrogen mole fraction.
- 800 ppm dry carbon dioxide target.
- 285 K mean surface-temperature target.
- 3 × 10^17 kg reference water inventory.
- Exact altitude boundaries below 35 km.
- Century-by-century implementation schedule.
- Specific volatile sources.
- Requirement for routine artificial-gravity exposure.

These values should be treated as editable assumptions.

## Canon protection rules for an external workspace

1. Do not silently replace recovered canon with a more Earthlike design.
2. Preserve the five-hour practical twilight unless a formal revision is recorded.
3. Preserve the long lunar solar cycle as a central biological and cultural feature.
4. Do not use atmospheric L1–L5 notation.
5. Keep known science, engineering inference and speculative biology visibly separated.
6. Do not claim a settled atmospheric lifetime until a lunar-specific dense-atmosphere model exists.
7. Do not describe Terluna as near-term.
8. Do not erase the ethical cost of destroying the Moon’s present vacuum environment.

## Unknown prior decisions to recover if older notes become available

- Final atmospheric pressure and composition.
- Intended ocean coverage and basin map.
- Whether orbital mirrors are routine, emergency-only, or absent.
- Whether the climate permits widespread night freezing.
- Desired population scale.
- Preferred political structure.
- How much of the Moon remains a protected scientific reserve.
- Whether aerial organisms are primarily evolved, engineered, or both.
- Whether the atmosphere has an artificial magnetic shield.
- Exact names of ecological regions and settlements.

When older material is found, add it to a new canon layer rather than deleting this uncertainty ledger.

<!-- END FILE: 01_RECOVERED_CANON_AND_CORRECTIONS.md -->


---

<!-- BEGIN FILE: 02_FEASIBILITY_CASE.md -->

# Feasibility Case

## 1. Meaning of “can work”

For Terluna, “can work” means that no known conservation law or planetary constraint automatically prohibits a maintained atmosphere, liquid surface water, a biosphere and human settlement on the Moon. It does not mean that the project is economical, politically legitimate, or possible with current technology.

A useful feasibility claim must survive five tests:

1. **Atmospheric binding** — gas remains long enough to be maintained.
2. **Mass closure** — the required gases and water can be sourced.
3. **Energy closure** — extraction, transport and climate operations can be powered.
4. **Biological closure** — ecosystems and humans can function at 0.16 g and on the long solar cycle.
5. **Governance closure** — an irreversible transformation of a shared celestial body can be authorized and managed.

Terluna presently passes only the weakest version of the first statement: the Moon is capable of having a substantial atmosphere for nonzero periods.

## 2. Evidence that lunar atmospheres are possible

The present Moon has only a tenuous exosphere [S01]. That condition results from a low volatile inventory, low gravity, solar ultraviolet radiation, solar-wind interaction and billions of years without sustained outgassing.

Needham and Kring estimated that intense mare volcanism could have produced a lunar atmosphere with surface pressure around 1 kPa and a dissipation time around 70 million years [S05]. Later escape modeling found that the lifetime of a millibar carbon-monoxide atmosphere depends sharply on exobase temperature and solar EUV heating, ranging from hundreds of years in hot energy-limited cases to around a million years in cooler cases [S06].

Those studies differ because atmospheric escape is sensitive to assumptions. That disagreement is useful: it shows why simple statements such as “the Moon cannot hold air” or “heavy gases will last for millions of years” are both inadequate.

A dense N2–O2 atmosphere under the modern Sun could behave very differently from an ancient millibar CO atmosphere. Increased mass can move the exobase, alter radiative cooling, create an ionosphere, and change escape channels. Only a dedicated model can determine whether the reference atmosphere lasts 10^3, 10^5, 10^7, or more years.

## 3. Why low gravity does not imply instant loss

Individual molecules have thermal velocity distributions. Atmospheric loss depends on molecular mass, upper-atmosphere temperature, altitude, chemistry, ionization and plasma interaction, not merely surface gravity.

Nitrogen and oxygen molecules at temperate lower-atmosphere conditions move much slower than the Moon’s 2.38 km/s escape velocity [S02]. The critical region is the hot, rarefied upper atmosphere, where collisions are infrequent and the local escape speed is lower. Hydrogen and helium escape readily; nitrogen and oxygen are harder to lose thermally but can still be removed through photochemical and ion processes.

A thick atmosphere can therefore persist for a meaningful interval while still being a severe maintenance problem.

## 4. Solar input is favorable

The Moon receives essentially the same solar flux as Earth. Its global time-averaged incoming solar energy is about 340 W/m² before albedo, the same first-order value used for Earth.

This means Terluna does not need to move into a habitable orbit. The mean energy budget can support liquid water. The challenge is distribution: two weeks of sunlight followed by two weeks of darkness.

Slow-rotator climate models demonstrate mechanisms that can moderate extremes:

- Global-scale overturning circulation.
- Thick reflective clouds over strongly heated regions.
- Efficient wave adjustment of free-atmosphere temperatures.
- Oceanic heat transport.

These mechanisms support the plausibility of a controlled climate but cannot substitute for a lunar model [S11][S12][S13].

## 5. Low gravity creates genuine advantages

### Atmospheric flight

At the same temperature and pressure, air density is determined by gas composition, not surface gravity. An aircraft or animal weighs only one sixth as much on Terluna while flying in roughly terrestrial-density air. Wing loading and takeoff requirements fall dramatically.

### Construction

Weight-bearing loads are lower. Towers, bridges, canopies and large rotating structures can use less material for gravity loads, although wind and pressure loads remain.

### Launch energy

The vacuum orbital speed of the Moon is only about 1.68 km/s and escape speed about 2.38 km/s. Space access is energetically easier than from Earth, even though a thick atmosphere creates drag and moves useful orbits upward.

### Atmospheric radiation shielding

Surface pressure is weight per unit area. At 80 kPa and 0.16 g, the atmospheric mass above each square meter is about 49 tonnes, almost five times Earth’s sea-level column. That is potentially excellent radiation and meteoroid shielding after the atmosphere is built.

## 6. Local resources help, but do not close the project

### Oxygen

Lunar regolith contains about 45% oxygen by mass, bound in mineral oxides [S04]. Oxygen extraction also produces iron, aluminum, titanium, silicon-bearing material and other industrial feedstocks.

### Metals and ceramics

The Moon provides bulk construction mass, glass, ceramics, shielding, roads, pipes, tanks and electrical materials. These are essential because importing structural mass from Earth would be impossible at Terluna scale.

### Water

Water ice is confirmed in polar permanently shadowed regions, and molecular water or hydroxyl occurs elsewhere [S03][S07][S08]. It can bootstrap settlement and industry.

### Missing volatile inventory

The Moon does not appear to contain enough accessible nitrogen, hydrogen or carbon for an open global biosphere. Local resources reduce imports; they do not eliminate them.

## 7. The material bottleneck

An 80 kPa atmosphere requires 1.87 × 10^18 kg of gas. The reference nitrogen component is about 1.35 × 10^18 kg.

For scale:

- 1.35 × 10^18 kg is the mass of a pure-ice sphere about 136 km in diameter.
- If a source body is only 10% useful nitrogen by mass, the equivalent total processed body is several hundred kilometers across.
- A 3 × 10^17 kg hydrosphere is equivalent to a pure-water sphere about 83 km in diameter.

No single impact can deliver this safely. A direct high-speed impact would release catastrophic energy, eject material, alter the orbit locally and destroy infrastructure. Volatiles must be captured, processed and delivered in controlled packets over centuries.

## 8. The energy bottleneck

The reference atmosphere contains about 4.9 × 10^17 kg of oxygen. At an illustrative 20–50 MJ per kilogram of product oxygen, extraction requires roughly 10^25 joules.

Over 1,000 years, that component alone corresponds to about 300–800 TW of average power. Real total demand is higher because of:

- Mining and crushing.
- Heating molten salts.
- Metal refining.
- Volatile transport and capture.
- Water production.
- Construction.
- Atmosphere compression and distribution.
- Energy storage through the lunar night.
- Loss replacement.

This is an industrial civilization larger than present Earth’s energy system.

## 9. The strongest climate concern

A dense atmosphere can move heat, but the long night still allows land and shallow water to cool for hundreds of hours. If circulation and thermal storage are inadequate, enormous regions freeze. Freeze-thaw can be part of the ecology, but atmosphere-wide water collapse or permanent ocean ice is unacceptable.

The opposite risk is daytime overheating, particularly if water vapor and clouds create positive rather than stabilizing feedback. Both outcomes depend on ocean fraction, cloud microphysics, surface albedo, greenhouse gases and topography.

The concept works only if there are persistent liquid-water refugia and habitable temperature corridors through every cycle.

## 10. The strongest human concern

An atmosphere does not change lunar gravity. There is no evidence that humans can live, reproduce and develop normally for generations at 0.16 g.

Microgravity data show bone, muscle, vision, cardiovascular and immune problems. Partial gravity is probably better than microgravity, but the safe threshold and required exposure dose are unknown [S14][S15].

A viable Terluna therefore includes artificial gravity from the beginning. Rotating homes, medical centers, schools, nurseries, pregnancy habitats and exercise districts are baseline infrastructure rather than optional luxuries.

## 11. The strongest ethical concern

Global terraforming destroys the Moon’s current environmental state. It would:

- Contaminate ancient ice and regolith records.
- Erode or chemically alter geological surfaces.
- End most surface vacuum science.
- Change impact, dust and exosphere processes.
- Threaten heritage sites.
- Affect every existing and future lunar operator.

The Outer Space Treaty requires due regard, consultation over potentially harmful interference, responsibility for national activities and avoidance of harmful contamination [S19]. Present governance is nowhere near sufficient for an irreversible atmosphere project.

## 12. Full terraforming versus alternatives

### Open global Terluna

Benefits:

- Outdoor breathable environment.
- Planet-scale hydrology and ecology.
- Atmospheric transport and weather.
- Natural radiation shielding.

Costs:

- Order-10^18 kg gas imports.
- Global irreversible contamination.
- High atmospheric loss uncertainty.
- Rebuilt space-access architecture.

### Paraterraforming

Large roofed craters, valleys or basin networks use far less gas and preserve most of the lunar vacuum.

Benefits:

- Reversible and incremental.
- Lower volatile mass.
- Easier climate control.
- Protected low-gravity ecology experiments.

Costs:

- Vast structural membranes.
- Local failure modes.
- No global weather system.

### Hybrid Terluna

A thin global atmosphere provides dust suppression, heat transport and meteoroid protection while high-pressure ecological regions remain enclosed.

This architecture may deliver most practical benefits if dense-atmosphere escape or nitrogen logistics fail.

## 13. Required standard of proof

Terluna should be considered feasible only after:

- A validated dense-atmosphere lifetime range is established.
- The annual volatile replacement burden is under a small fraction of industrial throughput.
- Closed habitats survive many full lunar cycles without emergency imports.
- Fire, biology and human development are tested at 0.16 g.
- A nitrogen source and transport system are demonstrated.
- International institutions authorize the irreversible phase.

Until then, Terluna is a coherent research program and worldbuilding system, not a construction proposal.

<!-- END FILE: 02_FEASIBILITY_CASE.md -->


---

<!-- BEGIN FILE: 03_PHYSICAL_PARAMETERS_AND_CALCULATIONS.md -->

# Physical Parameters and Screening Calculations

## 1. Purpose and limits

This file makes the Terluna reference design numerically explicit. The calculations are order-of-magnitude checks intended to expose the controlling variables; they are not substitutes for a general circulation model, kinetic escape model, combustion program, ecological model, or final engineering design.

The reference case uses an **80 kPa** dry surface atmosphere, **285 K** reference temperature, and **3 × 10^17 kg** water inventory. Those values are reconstructed modeling choices rather than recovered canon.

## 2. Fixed lunar quantities

| Quantity | Value | Use |
|---|---:|---|
| Mean radius | 1,737.4 km | Surface area and orbital scales |
| Surface area | 3.793 × 10^13 m² | Atmosphere, water and land inventories |
| Surface gravity | 1.624 m/s² | Pressure-column relation, lapse rate and structures |
| Escape speed | about 2.38 km/s | Atmospheric escape and transport screening |
| Vacuum circular speed near the surface | about 1.68 km/s | Space-access reference before atmosphere |
| Synodic solar cycle | 29.530588 Earth days | Sunrise-to-sunrise environmental cycle |
| Illumination half-cycle | 14.765294 Earth days | Approximate bright and dark halves near the equator |
| Solar angular rate | 0.50795 degrees/hour | Twilight and moving-terminator geometry |
| Equatorial terminator ground speed | about 4.3 m/s | Mobile ecology, weather and operations |

The Moon is synchronously rotating relative to Earth but still rotates once per sidereal month. Atmospheric dynamics therefore feel a rotation rate about one twenty-seventh of Earth’s, while the Sun crosses the sky on the 29.53-day synodic period [S01][S02].

## 3. Atmospheric mass

For a hydrostatic global atmosphere with surface pressure `P`, surface area `A`, and surface gravity `g`, the first-order total mass is:

`M_atm = P A / g`

At 80,000 Pa:

`M_atm ≈ 80,000 × 3.793 × 10^13 / 1.624 ≈ 1.87 × 10^18 kg`

### Pressure sensitivity

| Surface pressure | Atmosphere mass | Main interpretation |
|---:|---:|---|
| 5 kPa | 1.17 × 10^17 kg | Thin global experimental atmosphere |
| 10 kPa | 2.34 × 10^17 kg | Climate and escape plateau; not open-air breathable |
| 20 kPa | 4.67 × 10^17 kg | Very low-pressure habitat design space |
| 40 kPa | 9.34 × 10^17 kg | Intermediate plateau |
| 60 kPa | 1.40 × 10^18 kg | Lower edge of reference open-world range |
| 80 kPa | 1.87 × 10^18 kg | Reference case |
| 101.325 kPa | 2.37 × 10^18 kg | Earth sea-level pressure on the Moon |

Every additional kilopascal requires about **2.34 × 10^16 kg** of gas. Pressure optimization therefore has civilization-scale consequences.

## 4. Reference gas inventory

The reference dry mole fractions are 75% N2, 24% O2, 0.9% Ar, 800 ppm CO2, and 200 ppm other trace gases. The resulting mean molar mass is about 29.09 g/mol.

| Gas | Partial pressure | Approximate mass | Function |
|---|---:|---:|---|
| Nitrogen | 60.0 kPa | 1.35 × 10^18 kg | Inert buffer, nitrogen cycle and pressure reserve |
| Oxygen | 19.2 kPa | 4.93 × 10^17 kg | Respiration, oxidation and ozone precursor |
| Argon | 0.72 kPa | 2.31 × 10^16 kg | Optional inert buffer |
| Carbon dioxide | 0.064 kPa | 2.26 × 10^15 kg | Photosynthesis and greenhouse control |
| Other trace gases | 0.016 kPa | about 3.7 × 10^14 kg | Water vapor, ozone, methane and controlled species |

Mole fraction and mass fraction are different. Oxygen is 24% by molecule count but about 26.4% by mass; nitrogen is 75% by molecule count and about 72.2% by mass.

The reference oxygen partial pressure is close to Earth sea-level oxygen availability, but this does not establish fire safety. Ignition, flame spread and material compatibility depend on oxygen fraction, oxygen partial pressure, total pressure, humidity, flow and low-gravity combustion behavior.

## 5. Column mass and shielding

The column mass above one square meter is:

`m_column = P / g ≈ 49,261 kg/m² at 80 kPa`

This is about 4.8 times the mass column of Earth’s sea-level atmosphere because the same pressure requires more overlying mass in lower gravity. That large column is potentially valuable for cosmic-ray, solar-particle, ultraviolet and meteoroid protection, although shielding quality must be computed with particle-transport models because secondary radiation can matter [S10].

## 6. Scale height and vertical depth

For an isothermal atmosphere:

`H = R T / (M g)`

At 285 K and mean molar mass 29.09 g/mol, `H ≈ 50.2 km`. The reference dry adiabatic lapse rate is approximately:

`Γ_d = g / c_p ≈ 1.62 K/km`

Both values are radically different from Earth’s because lunar gravity is weak. A convective atmosphere can therefore be tens of kilometers deep without exhausting its temperature margin.

### Isothermal pressure screen

| Altitude | Screening pressure from 80 kPa surface |
|---:|---:|
| 10 km | 65.6 kPa |
| 20 km | 53.7 kPa |
| 35 km | 39.8 kPa |
| 40 km | 36.1 kPa |
| 45 km | 32.6 kPa |
| 100 km | 10.9 kPa |
| 200 km | 1.49 kPa |
| 300 km | 0.202 kPa |
| 500 km | 0.00375 kPa, or 3.75 Pa |

These numbers deliberately over-simplify the real atmosphere. Temperature, changing gravity, condensation, photochemistry, molecular diffusion and a hot thermosphere all alter the profile. The calculation is still enough to show that present-day 50–100 km low lunar orbits would lie deep inside a thick atmosphere.

## 7. Twilight geometry

The Sun moves across the lunar sky at about `360° / 708.734 h = 0.50795°/h`. A five-hour practical transition corresponds to about **2.54 degrees** of solar elevation change.

Terluna therefore defines the canonical five hours as the operational brightening or dimming interval around local sunrise or sunset. Very faint astronomical scattering can begin earlier, and mountains can delay or advance local illumination. The five-hour transition is part of the 708.7-hour cycle; it is not extra time added to the bright and dark halves.

A schematic equatorial cycle is:

- Practical dawn: 5 h.
- Full daylight: about 349.37 h.
- Practical dusk: 5 h.
- Full night: about 349.37 h.

## 8. Hydrosphere scales

Global equivalent depth is:

`GED = M_water / (ρ_water A)`

| Water mass | Global equivalent depth | Mean depth over 10% of surface |
|---:|---:|---:|
| 1 × 10^17 kg | 2.64 m | 26.4 m |
| 3 × 10^17 kg | 7.91 m | 79.1 m |
| 1 × 10^18 kg | 26.36 m | 263.6 m |

The 3 × 10^17 kg reference inventory is equivalent to a pure-water sphere about 83 km in diameter. That comparison illustrates the transport problem without implying delivery as one impactor.

If all reference water were manufactured using lunar oxygen, the imported hydrogen alone would be about **3.36 × 10^16 kg**. Importing water or hydrated material may be operationally easier, while indigenous polar ice should first be treated as a limited scientific and strategic resource [S03][S07][S08].

## 9. Oxygen production from regolith

Lunar regolith is roughly 40–45% oxygen by mass bound in oxides, and oxygen-extraction processes have been demonstrated at laboratory scale [S04]. If feedstock contains 45% oxygen and a planetary process recovers 75% of it, producing the reference 4.93 × 10^17 kg atmospheric oxygen requires processing:

`M_regolith = M_O2 / (0.45 × 0.75) ≈ 1.46 × 10^18 kg`

At a bulk density of 1,600 kg/m³, that equals a global average layer near 24 m, although real mining would be concentrated into selected industrial provinces. The operation yields metals or metal alloys as coproducts, so the atmosphere program and construction economy should be designed together.

## 10. Oxygen-production energy

Using an illustrative system energy intensity of 20–50 MJ per kilogram of product oxygen gives:

- Total process energy: about **9.9 × 10^24 to 2.5 × 10^25 J**.
- Average power over 1,000 years: about **313 to 782 TW**.

This excludes much of the mine development, comminution, transport, volatile import, nitrogen processing, water circulation, habitat construction and loss replacement. It is a scale marker, not a cost estimate.

## 11. Planetary thermal-inertia check

At a screening Bond albedo of 0.30, Terluna would absorb roughly **9 petawatts** averaged over the illuminated disk and time. Over one 14.765-day half-cycle, the absorbed energy is of order **1.1 × 10^22 J**.

The heat capacity of 3 × 10^17 kg of liquid water is about 1.25 × 10^21 J/K, so a 10 K participating temperature swing stores about **1.25 × 10^22 J**. The 80 kPa atmosphere has a comparable 10 K sensible-heat scale of about **1.9 × 10^22 J**.

This is one reason the long night is not automatically fatal: the reference atmosphere and hydrosphere contain enough thermal capacity in principle to buffer energy on the same order as a half-cycle’s absorbed sunlight. Distribution, clouds, ice, infrared cooling and heat-transfer rates determine whether that capacity is actually usable.

## 12. Flight scaling

For the same aircraft, atmospheric density, wing area and lift coefficient, stall speed scales with the square root of weight. Reducing gravitational acceleration to 0.16 g gives an idealized stall-speed factor of:

`sqrt(0.16) ≈ 0.40`

A vehicle that stalls at 25 m/s on Earth could, ignoring redesign and atmospheric differences, stall near 10 m/s on Terluna. Aircraft can trade this benefit for larger payload, smaller wings, slower flight or higher operating altitude.

Buoyant lift in kilograms does not gain the same one-sixth factor because both buoyant force and payload weight scale with gravity. Airships still benefit from lower structural loads and a deep atmosphere.

## 13. Artificial-gravity geometry

For a rotating habitat, `a = ω²r`. Representative 1 g rotations are:

| Radius | Rotation rate for 1 g |
|---:|---:|
| 50 m | 4.23 rpm |
| 100 m | 2.99 rpm |
| 250 m | 1.89 rpm |
| 500 m | 1.34 rpm |

Large radii reduce rotation rate and head-to-foot gravity gradients. Terluna settlements should reserve structural corridors for rotating sleep quarters, clinics, maternity facilities and exercise habitats until partial-gravity thresholds are known [S14][S15].

## 14. Imported-material momentum and impact energy

A volatile program must never treat a giant comet impact as an ordinary delivery method. A 10^15 kg packet arriving at 2 km/s carries `2 × 10^21 J`; a 10^18 kg body at the same speed carries `2 × 10^24 J`.

Resources must be subdivided, processed and decelerated through high-orbit depots, electromagnetic systems, tethers, propulsion or controlled aerocapture after a suitable atmosphere exists. Momentum disposal is part of the mass budget.

## 15. Mass fraction of the Moon

The 80 kPa atmosphere is only about **2.5 × 10^-5** of lunar mass. It does not meaningfully alter the Moon’s orbit or surface gravity. The engineering difficulty comes from acquiring, processing and retaining the material, not from its effect on lunar bulk dynamics.

## 16. Reproducibility

The script `models/terluna_screening_calculations.py` reproduces core formulas. The CSV files in `data/` preserve pressure, altitude, water, mass and energy sensitivities. All outputs should be replaced by higher-fidelity models as the concept develops.

<!-- END FILE: 03_PHYSICAL_PARAMETERS_AND_CALCULATIONS.md -->


---

<!-- BEGIN FILE: 04_ATMOSPHERE_CLIMATE_AND_WEATHER.md -->

# Atmosphere, Climate and Weather

## 1. Climate problem statement

Terluna receives nearly the same solar flux as Earth because it shares Earth’s orbit, but it distributes that energy through a 29.53-day sunrise-to-sunrise cycle, low gravity, weak rotation, extreme topography and a physically deep atmosphere. The climate system must keep substantial regions habitable through about 14.8 days of continuous illumination and about 14.8 days without direct sunlight.

The design objective is not uniform weather. It is a bounded climate with survivable refugia, no global atmospheric collapse, no uncontrolled greenhouse excursion, no recurring biosphere-wide freeze, and predictable seasonal and terminator hazards.

## 2. Why a dense atmosphere could help

A dense atmosphere provides four essential functions:

1. **Pressure:** liquid water and unprotected human respiration become possible.
2. **Heat transport:** winds carry day-side energy toward night-side terrain.
3. **Thermal storage:** the atmosphere itself stores heat across the long cycle.
4. **Shielding:** the large atmospheric column attenuates radiation, ultraviolet light and small impactors.

The ancient Moon probably carried a transient volcanogenic atmosphere near the kilopascal scale, showing that lunar atmospheres can exist for geologically nonzero periods [S05]. Later modeling shows that lifetime is highly sensitive to upper-atmosphere temperature and solar forcing [S06]. Neither result settles the lifetime of a modern 60–100 kPa nitrogen–oxygen atmosphere.

## 3. Reference composition and alternatives

### Reference composition

- 75% N2.
- 24% O2.
- 0.9% Ar.
- 800 ppm CO2.
- 200 ppm controlled trace gases before variable water vapor.

### Alternative architectures to model

**Lower-pressure oxygen-rich:** 40–60 kPa with 28–35% O2 reduces total imported gas, but raises combustion and material-compatibility concerns.

**Nitrogen-heavier:** 80–100 kPa with 20–23% O2 provides familiar oxygen fraction and better fire margins, but greatly increases imported nitrogen.

**CO2-assisted early atmosphere:** an initially CO2-rich nonbreathable atmosphere can aid greenhouse warming and industrial carbon storage, then be drawn down before oxygenation. Toxicity, condensation and carbonate reactions constrain this path.

**Argon or other inert supplementation:** argon is chemically useful but resource-constrained. Helium and hydrogen are poor long-term buffer gases because their thermal escape is more severe.

No composition is selected solely by breathing requirements. The optimizer must include escape, greenhouse effect, ozone, fire, plant productivity, nitrogen cycling, gas supply and upper-atmosphere cooling.

## 4. Vertical structure

Low gravity gives an Earthlike gas mixture a scale height around 50 km at 285 K and a dry adiabatic lapse rate near 1.62 K/km. A 35–45 km main tropopause is therefore a reasonable worldbuilding range, but it is not a derived prediction.

The accessible Terluna concept includes a top-of-troposphere/lower-stratosphere operating region near 35–45 km with possible weak-shear windows, atmospheric gravity waves, rare microbial life and engineered radiation-hard algal films. Long-endurance platforms there can support farside astronomy and nearside Earthlight photometry.

### Named operational bands

- **Surface boundary layer, 0–3 km:** settlement weather, dust, fog, local convection and slope flows.
- **Lower troposphere, 3–12 km:** most precipitation, regional aviation and low clouds.
- **Middle troposphere, 12–25 km:** long-range aircraft, layered clouds and broad return flow.
- **Upper troposphere, 25–35 km:** convective outflow, frontal structure and storm avoidance.
- **Tropopause/lower stratosphere, 35–45 km:** stratified research and endurance operations.
- **Extended atmosphere, above 45 km:** atmospheric chemistry, thermosphere, exosphere and escape interface.

The labels L1–L5 are never used for these layers; they refer only to Lagrange points.

## 5. Rotation and circulation

The Moon’s weak rotational rate reduces Coriolis deflection. Expected consequences include:

- Very broad overturning cells rather than many narrow latitude bands.
- Large Rossby deformation lengths and planetary-scale weather structures.
- Weaker familiar jet-stream organization, although strong jets or superrotation can still emerge.
- Deep day-to-night circulation coupled to the moving solar heating pattern.
- Strong topographic control from crater rims, maria, massifs and basin walls.

Slow-rotation terrestrial climate studies show that broad overturning and cloud feedbacks can stabilize climates in some regimes [S11][S12][S13]. Terluna differs in gravity, surface pressure, topography, month-long forcing and hydrosphere geometry, so those studies provide mechanisms rather than answers.

## 6. The moving terminator system

At the equator the sunrise or sunset line advances across the ground near 4.3 m/s, or about 15 km/h. This creates a persistent moving climate zone:

- Dawn fronts warm cold surfaces and release frost, fog and stored gases.
- Dusk fronts shut down convection and can create drainage winds into basins.
- Organisms and vehicles can track the terminator locally, although oceans, mountains and protected areas interrupt routes.
- Industrial operations can schedule heat-intensive work in dawn sectors and maintenance in dusk sectors.

The term “terminator front” means a weather response to moving solar forcing, not a physical wall. Its strength depends on atmosphere and surface thermal inertia.

## 7. Canonical five-hour twilight

The practical dawn and dusk interval is canonically about five hours. Since the Sun moves only 0.508 degrees per hour, five hours spans about 2.54 degrees of solar elevation.

The project should define practical twilight photometrically: for example, the interval between a chosen outdoor-work illuminance threshold and direct sunrise. Faint atmospheric twilight may last longer, and terrain can create locally abrupt or delayed transitions. The five-hour value is an operational canon rather than a claim that all scattered light begins or ends within exactly five hours.

## 8. Day-side climate

Potential day-side conditions include prolonged surface heating, thick convective boundary layers, persistent cloud shields over wet regions, and strong evaporation from seas. A slowly moving subsolar region can favor high, reflective clouds that limit absorbed solar energy, but cloud microphysics at low gravity may behave differently from terrestrial cases.

Critical day-side controls are:

- Cloud albedo and droplet size.
- Soil moisture and evaporative cooling.
- Ocean mixing and heat uptake.
- Vegetation transpiration.
- Dust and aerosol optical depth.
- Greenhouse-gas concentration.
- Topographic ventilation of basins.

Desert surfaces can still become extremely hot even when the global mean is acceptable. Settlement codes must use local extrema, not global averages.

## 9. Night-side climate

During the long night, the system depends on heat stored in water, atmosphere, soil, subsurface reservoirs and engineered infrastructure. Expected features include:

- Strong radiative cooling of exposed highlands and dry deserts.
- Fog, frost and low cloud in basins.
- Katabatic flows descending crater walls.
- Partial ice cover on shallow lakes and sea margins.
- Deep-water and geothermal/industrial refugia remaining liquid.
- Long-lived night storms where warm water feeds cold air.

A thick nitrogen–oxygen atmosphere will not condense at plausible habitable temperatures, but water and CO2 can phase-change locally. Preventing a water-cycle lockup in cold traps is a design problem.

## 10. Hydrological weather

Low gravity changes precipitation in several competing ways:

- Drops fall more slowly for the same size.
- Drops may grow larger before aerodynamic breakup.
- Cloud particles remain suspended longer.
- Convective towers can be physically deeper.
- Snow settles slowly and can drift for long distances.
- Reduced hydrostatic gradients alter cloud pressure-depth relationships.

Rain may be slower, larger-drop and more spatially persistent than Earth rain, but this requires microphysical simulation and parabolic-flight or orbital experiments. Hail, icing and mixed-phase clouds are major aviation risks.

## 11. Winds and storms

Low gravity does not imply weak winds. Atmospheric pressure gradients accelerate air independent of an object’s weight, and a deep atmosphere can organize enormous flows. At the same time, buildings and trees weigh less, so anchoring against wind becomes more important.

Likely hazards include:

- Basin windstorms and crater-rim downslope jets.
- Terminator pressure surges.
- Deep convective systems above warm seas.
- Planetary-scale waves.
- Mountain-launched atmospheric gravity waves.
- Dust storms over unvegetated regolith.
- Long-period oscillations tied to the solar cycle.

“Atmospheric gravity waves” are oscillations restored by buoyancy in stratified air. They can cause turbulence, cloud bands and platform loading; they are unrelated to relativistic gravitational waves.

## 12. Atmospheric superrotation

A slowly rotating atmosphere may develop winds that circulate faster than the solid surface. Superrotation could improve global heat transport but also create persistent high-altitude jets and difficult aviation corridors. The direction and strength depend on wave–mean-flow interactions, topography, atmospheric depth and day-night heating.

Superrotation must be treated as a model outcome, not assumed absent because the Moon rotates slowly.

## 13. Atmospheric chemistry

A mature atmosphere requires controlled cycles for:

- O2 and O3.
- CO2, CO and methane.
- Water vapor and hydrogen escape.
- Nitrogen species, including N2O and NOx.
- Sulfur and halogen compounds.
- Mineral aerosols and sea salt.
- Industrial contaminants.

Ozone can provide ultraviolet shielding after oxygenation, but its vertical distribution in a 50-km-scale-height atmosphere is unknown. Photochemistry must couple to the thermosphere because upper-atmosphere heating strongly affects escape.

## 14. Escape mechanisms

Relevant loss processes include:

- Jeans escape from the high-energy tail of molecular velocities.
- Hydrodynamic or energy-limited outflow under strong heating.
- Photochemical escape.
- Sputtering by solar-wind particles.
- Ion pickup and plasma interaction.
- Impact erosion.
- Preferential loss of hydrogen from water and methane.
- Leakage induced by industrial operations or high-altitude launch systems.

A dense atmosphere can cool its thermosphere through radiatively active species, reducing escape, while ultraviolet absorption can heat it. The outcome is composition- and state-dependent. Terluna’s decisive model must extend from surface weather through the exobase and include solar cycles, flares and passages through Earth’s magnetotail [S06].

## 15. Climate control system

Terluna is designed as an actively managed world. The control system includes:

- Global pressure, composition and isotope monitoring.
- Weather satellites and high-altitude platforms.
- Ocean heat-content and ice sensors.
- Carbon capture and release plants.
- Atmospheric water management.
- Aerosol restrictions and emergency removal capacity.
- Orbital mirrors or shades only as contingency tools, not routine assumptions.
- Protected no-intervention climate reference zones.

Interventions should prefer slow, reversible changes. Large aerosol injections or global albedo manipulation can create cross-border and ecological harms even when technically effective.

## 16. Aviation and aeroship routing

Aircraft exploit low weight, but routes must account for a much deeper weather column. A practical network uses:

- Surface and low-altitude rotorcraft for local transport.
- Fixed-wing aircraft for regional and intercontinental movement.
- Heavy cargo airships where calm corridors exist.
- High-altitude aircraft for fast long-range travel.
- Stratospheric platforms for communications and science.
- Terminator-following routes for energy and thermal advantages.

Aviation weather maps require vertical resolution through at least 45 km. High-altitude “calm” layers are probabilistic windows, not permanent global highways.

## 17. Required climate simulations

The minimum credible model hierarchy is:

1. One-dimensional radiative-convective columns.
2. Two-dimensional terminator circulation and basin cross-sections.
3. Global 3-D GCM with lunar topography and moving solar forcing.
4. Coupled ocean, lake, sea-ice, soil and groundwater models.
5. Interactive clouds, aerosols, vegetation and atmospheric chemistry.
6. Thermosphere–exosphere kinetic escape model.
7. Ensemble runs across pressure, composition, water placement and solar activity.
8. Hardware-in-the-loop control simulations for atmospheric intervention.

The concept passes the climate gate only when multiple independent models and staged real-world pressure experiments converge.

<!-- END FILE: 04_ATMOSPHERE_CLIMATE_AND_WEATHER.md -->


---

<!-- BEGIN FILE: 05_HYDROSPHERE_GEOSPHERE_AND_SOIL.md -->

# Hydrosphere, Geosphere, Soil and Dust

## 1. Hydrosphere design objective

Terluna does not need an Earth-depth global ocean. The reference design uses **3 × 10^17 kg** of water, equivalent to a 7.9 m global layer or an average depth near 79 m if concentrated over 10% of the surface. The intended landscape is a network of regional seas, crater lakes, wetlands, rivers, aquifers, glaciers and controlled polar reservoirs.

Water performs ecological, climatic, industrial and cultural work. Its placement matters as much as its total mass.

## 2. Indigenous water

The Moon contains water ice in permanently shadowed polar regions and molecular water or hydroxyl in other materials [S03][S07][S08]. These discoveries establish local resources but do not establish a global hydrosphere inventory.

Before extraction, the program must:

- Map ice distribution and isotopes at high resolution.
- Preserve representative deposits as records of Solar System history.
- Separate scientific reserves from industrial zones.
- Establish ownership, access and contamination rules.
- Avoid assuming every spectral water signal is economically recoverable.

Indigenous lunar water is best treated as strategic seed stock, emergency reserve and early-industrial feedstock.

## 3. Imported water and hydrogen

Large hydrosphere inventories probably require imported ice, hydrated minerals, ammonia-rich material or hydrogen combined with lunar oxygen. Candidate sources require a real census rather than generic references to “asteroids.”

Source-selection criteria include:

- Total accessible hydrogen, nitrogen, carbon, phosphorus and sulfur.
- Delta-v and transfer time.
- Ability to divide material into safe packets.
- Scientific and ethical protection status.
- Processing contaminants.
- Momentum and heat disposal.
- Long-term replenishment capability.

A single uncontrolled impact is unacceptable. Volatiles are refined or packaged away from the surface, moved through high-orbit depots, and delivered gradually.

## 4. Basin selection

Candidate seas are placed in maria, impact basins and crater chains only after evaluating:

- Elevation and spill points.
- Crustal permeability and faulting.
- Basin volume as water level rises.
- Ground-ice interaction.
- Heritage and science sites.
- Settlement access and wind exposure.
- Seismic loading and mass redistribution.
- Night thermal behavior.
- Evaporation and salinity.

Some iconic basins may remain dry as cultural landscapes or planetary archives.

## 5. Water at low gravity

Hydrostatic pressure increases with depth as `ρgh`, only about one sixth as rapidly as on Earth. A 60 m water column on Terluna produces pressure comparable to roughly 10 m on Earth. Consequences include:

- Organisms can reach greater depths before pressure becomes limiting.
- Gas solubility and bubble behavior differ with depth.
- Submersibles and dams face lower pressure loads for a given depth.
- Waves and sloshing have longer natural periods.
- Deep reservoirs can store large water volumes with modest structural pressure.

Low gravity also changes capillary length, sediment settling, plume buoyancy and shoreline behavior. Every terrestrial hydraulic formula must be revalidated.

## 6. Rivers and erosion

Rivers can flow because slopes and pressure gradients still exist. Expected differences include:

- Slower sediment settling and longer suspended transport.
- Broader or deeper channels for some discharge regimes.
- Lower particle weight but unchanged chemical weathering needs.
- Greater importance of cohesive banks and vegetation roots.
- Long-lived waterfalls and mist plumes.
- Modified meandering, delta formation and flood-wave speed.

Fresh lunar terrain has no mature drainage network. Watersheds are excavated, graded or allowed to evolve over centuries under controlled flooding.

## 7. Tides and long-period oscillations

Earth remains nearly fixed in the nearside sky, so the Earth-raised tidal potential is mostly stationary rather than producing a familiar twice-daily tide. Solar tides, physical libration, orbital eccentricity, basin resonance and atmospheric pressure can still move water on multi-day or monthly periods.

Terluna seas may therefore have long, slow tides and basin seiches rather than Earthlike coastal cycles. Local geometry can amplify even weak forcing.

## 8. Freezing and the long night

The hydrosphere is designed to permit controlled surface freezing without total loss of liquid water. Techniques include:

- Deep basins with large heat capacity.
- Salinity gradients where ecologically acceptable.
- Dark heat-absorbing surfaces and daytime mixing.
- Insulating ice lids.
- Geothermal or industrial waste-heat refugia.
- Pumped circulation between warm and cold layers.
- Covered reservoirs for emergency water supply.

Ice is useful infrastructure: it insulates water, stores seasonal heat, supports transport and creates night habitats. The risk is cumulative migration of water into permanent cold traps.

## 9. Evaporation and atmospheric water

At 60–100 kPa, liquid water is stable across intended surface temperatures. Evaporation during the long day can load the deep atmosphere with water vapor, driving clouds and precipitation far from source basins.

The water cycle must avoid:

- Excess greenhouse warming from vapor.
- Chronic desiccation of continental interiors.
- Ice accumulation in protected or inaccessible basins.
- Salinity concentration.
- Atmospheric escape of photodissociated hydrogen.
- Uncontrolled storm belts downwind of seas.

Desalination, canal networks, groundwater recharge and atmospheric water capture become planetary utilities.

## 10. Raw lunar regolith

Lunar regolith is pulverized rock and impact glass, not soil. It contains sharp, abrasive particles, agglutinates, reactive surfaces and almost no organic carbon or biological structure [S16][S17]. Apollo regolith supported seed germination in a small experiment, but plants grew slowly and showed strong stress responses [S09].

The result demonstrates biological interaction, not agricultural readiness.

## 11. Dust transition after atmosphere

A global atmosphere changes dust behavior:

- Electrostatic lofting and ballistic transport are reduced near the surface.
- Wind erosion and dust storms become possible.
- Fine particles can enter lungs and machinery through ordinary aerosol pathways.
- Rain removes dust but transfers contaminants into soil and water.
- Vegetation and crusts stabilize mature surfaces.

Dust management shifts from vacuum seals and electrostatics toward terrestrial-style filtration, paving, wet suppression, ground cover and land restoration.

## 12. Soil manufacturing pipeline

A scalable soil program includes:

1. **Screening:** remove glass shards and undesirable size fractions.
2. **Washing and leaching:** reduce soluble toxins and process salts.
3. **Mechanical weathering:** crush or round grains where needed.
4. **Chemical conditioning:** set pH, redox state and cation balance.
5. **Nutrient addition:** nitrogen, carbon, phosphorus, sulfur, potassium and trace elements.
6. **Organic structure:** biochar, compost, polymers and root fibers.
7. **Biological inoculation:** bacteria, archaea, fungi, algae and soil fauna.
8. **Cover-crop cycles:** build aggregates and organic matter.
9. **Monitoring:** metals, perchlorates if introduced from imported material, salts and pathogens.
10. **Landscape release:** expand only after runoff and food-chain tests.

Soil is an industrial-biological product measured in centuries and cubic kilometers.

## 13. Nutrient budgets

Oxygen and metals are abundant in rock. The limiting biosphere elements are volatile or concentrated nutrients:

- Nitrogen for proteins and nucleic acids.
- Carbon for biomass and soils.
- Hydrogen for water and organics.
- Phosphorus for energy metabolism and membranes.
- Sulfur for proteins and redox chemistry.
- Biologically available iron, magnesium, calcium and trace metals.

Atmospheric nitrogen is not automatically plant-available. Nitrogen fixation, ammonium production, nitrate cycling and denitrification need engineered reservoirs and controls.

## 14. Carbonate and oxidation sinks

Fresh regolith can consume oxygen and carbon dioxide through oxidation and mineral reactions. Before atmospheric oxygenation, the program must quantify:

- Reduced iron and sulfur sinks.
- Glass-surface reactivity.
- Carbonate formation potential.
- Oxidation heat.
- Release of trapped solar-wind volatiles.
- Changes in mechanical strength during weathering.

An early atmosphere may be chemically drawn down by the ground faster than expected. Passivation tests at large scale are mandatory.

## 15. Groundwater and aquifers

Natural lunar crust is fractured by impacts, but its permeability and sealing behavior under water vary. Managed aquifers can be created by:

- Lining fractured basins.
- Sintering or vitrifying barriers.
- Injecting mineral grouts.
- Using ice barriers in cold zones.
- Excavating sealed reservoirs.

Groundwater is a heat and drought buffer, but uncontrolled leakage can move water into deep cold or vacuum-connected fractures.

## 16. Seismic and geotechnical design

The Moon experiences moonquakes and long seismic ringing. Adding seas, cities and mines changes local loading but remains a tiny fraction of lunar mass. Engineering priorities are:

- Map active faults and deep moonquake zones.
- Avoid reservoir-induced slip in susceptible structures.
- Design flexible dams and pipes.
- Monitor subsidence above mines and lava tubes.
- Anchor low-weight structures against wind and buoyancy.
- Account for weakly compacted regolith under cyclic wetting.

## 17. Protected terrain classes

A mature land-use map should include:

- Pristine vacuum geology reserves preserved before atmosphere.
- Polar ice archives.
- Apollo and other heritage sites.
- Farside radio-quiet zones.
- Dry lunar landscape parks.
- Industrial excavation provinces.
- Hydrological basins.
- Soil-development corridors.
- Settlement and transport zones.
- Biological exclusion areas.

Some reserves may require sealed domes or orbital sample archives once the global atmosphere makes true surface vacuum impossible.

## 18. Hydrosphere success criteria

A candidate design passes only if it demonstrates:

- Liquid refugia throughout repeated nights.
- No progressive transfer of most water to inaccessible cold traps.
- Stable salinity and nutrient balances.
- Manageable storms and floods.
- No unacceptable crustal or heritage impacts.
- Closed hydrogen-loss and replenishment budgets.
- Recoverable failure modes during regional pilots.

<!-- END FILE: 05_HYDROSPHERE_GEOSPHERE_AND_SOIL.md -->


---

<!-- BEGIN FILE: 06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md -->

# Biosphere and Evolutionary Design

## 1. Biosphere premise

Terluna’s biosphere is deliberately assembled. It cannot be created by scattering Earth seeds over wet regolith. Every trophic layer must be matched to the long solar cycle, 0.16 g, deep atmosphere, immature soils, engineered hydrology and active planetary control.

The design goal is a resilient ecology that can survive infrastructure outages and natural variability without becoming uncontrollable.

## 2. Two biological clocks

Organisms must respond to two distinct timing systems:

- **Short internal clocks:** roughly 20–30-hour cellular, sleep, feeding and metabolic rhythms inherited from Earth biology.
- **Long environmental clock:** the 708.7-hour light-dark cycle controlling seasonal-style growth, migration, dormancy and reproduction.

A Terlunan plant might cycle leaf chemistry daily while treating each two-week light period as a growth season and each dark period as winter. Human communities keep 24-hour civil days even as the landscape passes through month-long environmental phases.

## 3. Ecological succession

### Stage 0: sterile geochemical conditioning

Robots weather, wash and passivate regolith. Atmospheric composition is controlled without a free biosphere.

### Stage 1: microbial films

Selected bacteria, archaea, algae and fungi establish sealed mineral-processing communities. Their tasks are rock weathering, nitrogen transformation, carbon fixation, toxin immobilization and organic-matter production.

### Stage 2: crusts and pioneer plants

Lichen analogues, mosses, cyanobacterial crusts and low vascular plants stabilize dust and create soil structure. These systems remain easy to quarantine.

### Stage 3: grassland, shrubland and wetlands

Robust root networks, detrital loops and pollinators expand. Water and nutrient cycling become regional rather than facility-scale.

### Stage 4: forests and complex aquatic systems

Trees, large fungi, fish, amphibious organisms and complex food webs are introduced after repeated long-cycle stability.

### Stage 5: large animals

Large vertebrates or engineered analogues arrive only when food, disease, reproductive and fire models are reliable. Size is constrained by ecology and physiology, not merely by low gravity.

## 4. Plant adaptation to prolonged light

A two-week day risks photodamage, water loss and carbohydrate imbalance. Candidate traits include:

- Movable or folding leaves.
- Adjustable chlorophyll and protective pigments.
- Reflective hairs or waxes.
- Daytime growth pauses despite continuing sunlight.
- Large underground starch, sugar and oil stores.
- Deep or distributed roots.
- High-capacity vascular transport.
- Heat-shock and oxidative-stress systems.
- Symbiosis with fungi and nitrogen fixers.
- Reproduction timed to dawn, dusk or specific weather windows.

Some crops grow continuously under managed shade; others use automated canopies to preserve familiar photoperiods.

## 5. Plant adaptation to prolonged night

Night survival strategies include:

- Metabolic suppression and dormancy.
- Stored carbohydrates and lipids.
- Antifreeze proteins and compatible solutes.
- Deciduous leaf shedding before dusk.
- Insulated stems, bulbs and rhizomes.
- Mycorrhizal nutrient exchange.
- Heat sharing in dense plant mats.
- Low-level artificial light in managed agriculture.

Night ecosystems shift from photosynthesis toward decomposition, fungi, stored-energy consumption and aquatic production under ice or artificial illumination.

## 6. Low-gravity plant form

Reduced weight permits tall stems and broad canopies, but wind loads do not shrink with plant weight. A very tall plant can be easy to hold up in still air yet easy to uproot in a storm.

Useful adaptations include wide anchoring mats, deep cable-like roots, flexible trunks, porous canopies and active posture control. Gravitropism may weaken, so directional light, magnetic cues or engineered developmental pathways help orient growth.

## 7. Primary-production geography

The world does not photosynthesize uniformly. Productivity concentrates in:

- Dawn belts with moderate light and thawing water.
- Wet day-side regions with cloud protection.
- Shallow seas and wetlands.
- Managed agricultural zones.
- Terminator-tracking aerial or mobile systems.
- Night greenhouses powered by stored or nuclear energy.

Dusk initiates harvest and storage, while late night is a period of ecological austerity.

## 8. Microbial governance

Microbes are the foundation and the greatest contamination risk. The program maintains:

- Fully sequenced strain libraries.
- Genetic dependency systems for early releases.
- Competing consortia rather than single monocultures.
- Environmental DNA monitoring.
- Phage and antimicrobial contingency banks.
- Protected wild-type and no-life reference zones.
- Controls against methane, nitrous oxide and toxin overproduction.

No kill switch is assumed perfectly reliable. Physical containment and ecological redundancy remain necessary.

## 9. Fungi and detrital systems

The long night favors fungi and decomposers because stored biomass must be recycled when photosynthesis stops. Fungal networks can move carbon, water and minerals across root communities, stabilize soil and preserve nutrients.

The danger is an oxygen-consuming decomposition pulse late in the night. Atmospheric and soil oxygen sensors must detect regional respiration surges before they become fire or anoxia problems.

## 10. Aquatic biosphere

Seas and deep lakes provide the most reliable night refugia. Candidate systems include:

- Cold-tolerant algae under ice.
- Vertical migration between light and nutrient layers.
- Large dissolved-oxygen reservoirs.
- Artificial mixing and aeration.
- Detritus-based food webs during darkness.
- Thermal refuges near deep water or industrial heat sources.

Low gravity changes swimming, buoyancy, gas exchange and sedimentation. Aquatic organisms require multigenerational testing.

## 11. Aerial ecology

The combination of dense air and one-sixth weight creates extensive flight niches. Possible organisms include:

- Large low-wing-loading soarers.
- Cliff-launching gliders.
- Long-endurance migrators following the terminator.
- Ballooning organisms using heated or light gases.
- Aerial plankton, spores and small filter-feeders.
- Canopy-to-canopy gliders.
- Engineered pollinators serving broad landscapes.

Low gravity does not remove muscle-power, oxygen-delivery or structural constraints. Very large fliers may be possible, but ecological productivity, launch mechanics and cardiovascular design still set limits.

## 12. High-altitude life

Above most conventional ecosystems, the 35–45 km candidate platform layer can host:

- Radiation-tolerant microbes in droplets or particles.
- Engineered algal films attached to research platforms.
- Spore transport corridors.
- Sensors that distinguish natural dispersal from contamination.

The accessible project notes describe life there as sparse. Free-floating complex ecosystems at that altitude are not assumed.

## 13. Animal night strategies

Animals can survive the long night through combinations of:

- Fat and glycogen storage.
- Torpor or hibernation.
- Food caching.
- Communal roosts.
- Burrows and cave habitats.
- Migration toward water, cities or the advancing dawn.
- Diet switching from fresh vegetation to detritus, fungi or stored seeds.

Predator–prey cycles may synchronize with the lunar environmental month. Reproduction could occur around dawn when future food availability is greatest.

## 14. Locomotion and body form

Low gravity enables long jumps, slow falls and lower skeletal loading. Animals may evolve or be designed with:

- Long limbs and elastic tendons.
- Large stabilizing tails or membranes.
- Gripping feet for wind safety.
- Lower bone density balanced against impact and muscle attachment.
- Vestibular systems adapted to slow ballistic motion.
- Strong cardiovascular pumps despite reduced hydrostatic gradients.

Bodies still have inertia. A massive animal moving quickly remains dangerous even if it weighs less.

## 15. Pollination and seed dispersal

Wind can carry pollen and seeds much farther in a deep, low-gravity atmosphere. This helps colonization but threatens protected zones. Strategies include:

- Heavy or sticky pollen near reserves.
- Sterile buffer landscapes.
- Pollinator-specific flowers.
- Seed dormancy tied to soil chemistry.
- Geofenced engineered dependencies.
- Atmospheric eDNA and spore monitoring.

## 16. Fire ecology

Wildfire behavior is a central biosphere gate. Low gravity changes buoyant plume rise, oxygen transport, flame shape, ember lofting and smoke residence. Reduced plant weight may create delicate, wind-thrown fuel structures.

The landscape uses wet firebreaks, mineral corridors, low-flammability plants, compartmented watersheds, atmospheric monitoring and rapid aerial suppression. Oxygenation proceeds only after full-scale 0.16 g combustion research.

## 17. Atmospheric composition and biology

Biology is not allowed to set atmosphere passively. Planetary reservoirs and machines buffer:

- Oxygen against photosynthetic overshoot and fire.
- CO2 against climate swings and plant starvation.
- Fixed nitrogen against eutrophication and denitrification losses.
- Methane and N2O against greenhouse accumulation.
- Aerosols against unintended albedo change.

Terluna is closer to a managed planetary bioreactor than a self-regulating Earth analogue during its first millennia.

## 18. Genetic-engineering philosophy

The program favors conservative, modular changes:

- Modify stress tolerance and timing before inventing new metabolisms.
- Keep early organisms dependent on supplied nutrients or cofactors.
- Maintain unmodified Earth relatives in archives.
- Avoid horizontal-gene-transfer-prone constructs where possible.
- Test ecological interactions across whole lunar cycles.
- Record every released genome and location.

Evolution will eventually escape design intent. Governance must plan for adaptation rather than promise permanent genetic control.

## 19. Disease and invasive species

Dense settlements and engineered ecosystems create novel disease pathways. Risk controls include quarantine between basins, sentinel species, environmental sequencing, vector management, vaccinated animal stocks, and rapid local habitat isolation.

An organism harmless in a sealed habitat may become invasive under low gravity or long photoperiods. Release decisions require landscape-scale trials.

## 20. Biosphere success criteria

The open biosphere passes only after it can:

- Survive multiple full lunar cycles without emergency feeding or heating outside designated managed zones.
- Maintain bounded oxygen, carbon, nitrogen and water cycles.
- Recover from regional fire, freeze, drought and pathogen disturbances.
- Preserve protected sterile reserves.
- Support human food systems without displacing all wild ecology.
- Remain monitorable and governable despite ongoing evolution.

<!-- END FILE: 06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md -->


---

<!-- BEGIN FILE: 07_HUMAN_SETTLEMENT_HEALTH_AND_CULTURE.md -->

# Human Settlement, Health and Culture

## 1. Human-habitation premise

A breathable landscape solves pressure suits and immediate life support; it does not make the Moon physiologically Earthlike. Residents still live in 0.16 g, under a 29.53-day environmental cycle, in an engineered atmosphere whose composition and climate require permanent stewardship.

Settlements therefore combine open-air life with protected buildings, artificial gravity, emergency pressure integrity and closed-loop reserves.

## 2. Breathing environment

The 80 kPa reference case provides 19.2 kPa oxygen partial pressure. Habitability standards must also regulate:

- CO2 accumulation in cities and valleys.
- Carbon monoxide and industrial emissions.
- Ozone at the surface.
- Humidity and dehydration.
- Allergens, spores and dust.
- Radon or other crustal gases where relevant.
- Pressure variation with altitude and weather.

Because the atmosphere is deep, highlands may remain comfortably pressurized compared with equivalent terrestrial altitudes, but local pressure profiles require real models.

## 3. Partial-gravity health

Known spaceflight experience shows that altered gravity affects bone, muscle, cardiovascular function, vision, balance and other systems [S14]. Research also points to cellular and immune effects under lunar and Martian partial gravity, while the lifetime dose-response remains unknown [S15].

Critical unknowns are:

- Childhood skeletal and neurological development.
- Pregnancy, placentation and birth.
- Cardiovascular maturation.
- Immune function and infection.
- Vision and intracranial pressure.
- Fertility and multigenerational epigenetics.
- Minimum effective artificial-gravity dose.

No responsible design assumes 0.16 g is safe for an entire human life.

## 4. Artificial-gravity network

Artificial gravity is a public utility, not a luxury. Cities include:

- Rotating sleep residences.
- 1 g clinics and maternity centers.
- Rotating schools and child-development facilities.
- Exercise centrifuges.
- High-g rehabilitation and emergency treatment.
- Large rotating hotels for visitors adapting between worlds.

A 250 m radius habitat produces 1 g at about 1.9 rpm. Larger rings reduce vestibular effects and head-to-foot gradients. Some structures rotate internally while their outer shell remains fixed.

## 5. Settlement forms

### Crater-rim cities

Rims provide drainage, views, wind exposure and transport access. Buildings need deep anchors because aerodynamic loads can exceed their reduced weight.

### Basin cities

Basins offer water, agriculture and shelter but can trap cold air, fog, pollution and floodwater. Multiple elevated evacuation corridors are mandatory.

### Lava-tube districts

Tubes provide radiation, weather and impact protection, plus stable temperatures. They remain valuable even after atmospheric shielding develops.

### Floating or shoreline cities

Low hydrostatic loads and large seas support floating infrastructure, but long-period waves and ice require flexible moorings.

### High-altitude platform settlements

These remain specialized research and transport nodes. The 35–45 km layer may support long-endurance platforms, but human communities need pressure, radiation and evacuation margins.

## 6. Buildings under low gravity and wind

Structural dead loads are one sixth of Earth values, enabling large spans and tall towers. Wind pressure at a given air density and speed does not scale down with gravity. Consequences include:

- Strong foundation anchoring.
- Lower center-of-pressure designs.
- Permeable facades and wind gaps.
- Active damping for long-period structures.
- Conservative cladding and roof attachment.
- Deep root-like piles or rock anchors.

Snow and soil weigh less, but can accumulate deeply. Water loads depend on lower hydrostatic pressure but large volumes.

## 7. Fire safety

Open flame behavior at 0.16 g is not a minor variation of Earth fire. Reduced buoyancy changes flame shape, smoke rise and oxygen supply. A 24% oxygen fraction at 80 kPa has near-Earth oxygen partial pressure but may still change flame spread and ignition thresholds.

Cities use:

- Noncombustible urban cores.
- Pressure and oxygen isolation sectors.
- Wet and mineral firebreaks.
- Fast atmospheric sampling.
- Autonomous suppression aircraft.
- Smoke-safe refuge rooms with independent air.
- Landscape fuel management.

Global oxygenation waits for full-scale partial-gravity fire research.

## 8. Civil time

Human biology retains approximately 24-hour schedules. Terluna uses layered timekeeping:

- A standard civil day, likely 24 Earth hours for compatibility.
- A 29.53-day local solar month.
- Named environmental phases such as dawn, early light, high day, late light, dusk, early night and deep night.
- Local solar longitude for agriculture, weather and aviation.
- Earth-standard time for interplanetary coordination.

Communities may celebrate one sunrise and one sunset per month as major civic events.

## 9. Work and economy across the cycle

Economic rhythms follow environmental phases:

- Dawn: thaw management, planting, maintenance and migration.
- Day: solar-intensive industry, construction and high agricultural output.
- Dusk: harvest, storage, storm preparation and travel.
- Night: mining, indoor manufacturing, data work, culture, maintenance and nuclear-powered agriculture.

Continuous sectors use shifts. No worker remains awake for environmental “day” or “night” as a single human waking period.

## 10. Food systems

Food comes from a layered system:

- Open fields for robust staple crops.
- Controlled-environment farms preserving Earth photoperiods.
- Algae and microbial protein.
- Aquaculture in thermally stable reservoirs.
- Orchards and perennial storage crops.
- Night greenhouses using nuclear or stored energy.
- Strategic food reserves sized for multi-cycle failures.

Open agriculture is a supplement until the biosphere proves reliable.

## 11. Water and sanitation

Despite open seas, potable water remains a controlled utility. Systems separate ecological water from drinking, industrial and medical loops. Every settlement can isolate from the global hydrosphere during contamination events.

Waste nutrients are valuable but dangerous. Treatment plants recover nitrogen, phosphorus, carbon and water before release.

## 12. Transportation

Terluna supports an unusually broad transport mix:

- Walking, hopping and assisted mobility.
- Wheeled and tracked vehicles with low normal force and traction controls.
- Rail anchored for braking and wind stability.
- Rotorcraft and short-field aircraft.
- Airships and cargo balloons.
- High-altitude passenger aircraft.
- Suborbital rockets using controlled corridors.
- Tethers and mass drivers for space freight.

Low gravity reduces required lift but also tire traction and vehicle downforce. High-speed surface transport needs active adhesion or aerodynamic loading.

## 13. Medicine and emergency response

Hospitals must manage both ordinary disease and planetary-specific hazards:

- Decompression or atmospheric contamination.
- Low-gravity injury and vestibular disorders.
- Dust and metal exposure.
- Radiation during early phases or high-altitude operations.
- Fire and smoke with unusual buoyancy.
- Long-night mental health and sleep disruption.
- Imported-pathogen and engineered-organism incidents.

Medical transport connects every settlement to artificial-gravity treatment.

## 14. Psychology

The long environmental cycle can be beautiful and psychologically difficult. Design responses include:

- 24-hour indoor lighting.
- Access to artificial dawn and sunset.
- Night festivals and public winter gardens.
- Quiet dark-sky districts.
- Seasonal mental-health services.
- Routine travel between environmental phases where practical.

Nearside communities live under an almost fixed Earth, changing phase through the month; farside communities have no Earth in their sky. This difference will shape identity.

## 15. Nearside and farside cultures

### Nearside

Earth is a permanent geographic and emotional reference, with libration causing modest motion. Earthlight is brightest during local lunar night but is insufficient to replace sunlight for most ecosystems.

### Farside

The farside preserves direct isolation from Earth radio noise only if communications and atmospheric infrastructure are tightly managed. It may host protected astronomy, distinct political institutions and cultures oriented toward deep space.

## 16. Law and rights within the biosphere

Residents need rights concerning:

- Access to breathable air and artificial gravity.
- Climate and atmosphere transparency.
- Exposure to engineered organisms.
- Evacuation and refuge capacity.
- Participation in irreversible ecological decisions.
- Protection from unilateral pressure or composition changes.
- Cultural and scientific access to preserved lunar landscapes.

Atmospheric control operators hold power analogous to a combined water authority, grid operator, central bank and emergency government. Their authority requires unusually strong oversight.

## 17. Failure shelters

Every population center retains sealed refuges capable of surviving:

- Regional atmosphere contamination.
- Firestorms.
- Multi-week power disruption.
- Pathogen quarantine.
- Flood or freeze.
- Pressure decline.
- Solar-particle events during incomplete shielding phases.

Open-world confidence must never eliminate independent life support.

## 18. Human-settlement success criteria

A settlement system succeeds when residents can live ordinary lives outdoors while retaining:

- Proven partial-gravity mitigation.
- Independent emergency air and water.
- Fire-safe construction.
- Reliable food across multiple long cycles.
- Freedom of movement and civic participation.
- Health outcomes comparable with other human worlds.
- Cultural continuity without denying Terluna’s distinct environment.

<!-- END FILE: 07_HUMAN_SETTLEMENT_HEALTH_AND_CULTURE.md -->


---

<!-- BEGIN FILE: 08_INFRASTRUCTURE_INDUSTRY_AND_LOGISTICS.md -->

# Infrastructure, Industry and Logistics

## 1. Industrial premise

Terluna requires a preexisting lunar civilization capable of operating mines, refineries, power systems and transport networks at planetary scale. Terraforming is not the first lunar project; it is a late product of centuries of industrialization.

The core challenge is moving and processing roughly 10^18 kg of atmosphere plus 10^17–10^18 kg of water and nutrient material without catastrophic impacts or unsustainable losses.

## 2. Material hierarchy

### Locally abundant

- Oxygen bound in rock.
- Silicon, aluminum, iron, magnesium, calcium and titanium.
- Regolith for glass, ceramics, shielding and construction.
- Vacuum and sunlight during early industrial phases.

### Locally scarce or uncertain

- Nitrogen.
- Hydrogen and bulk water.
- Carbon.
- Phosphorus in concentrated accessible forms.
- Sulfur and other biological nutrients.
- Some catalysts and specialty elements.

The architecture uses lunar rock for most structure and imports volatile chemistry.

## 3. Nitrogen bottleneck

The reference atmosphere contains about **1.35 × 10^18 kg of nitrogen**. No demonstrated lunar deposit approaches this scale. Nitrogen can arrive as N2, ammonia, nitriles or nitrogen-bearing minerals.

Earth cannot be treated as a convenient source; removing a substantial fraction of its atmospheric nitrogen would be ecologically and politically unacceptable. Venus and Titan have large nitrogen inventories, but extraction from another planet or moon creates deep technical, scientific and ethical problems. The preferred conceptual path is a distributed outer-system and small-body volatile economy, contingent on actual resource surveys.

The nitrogen plan must include:

- Source inventory with uncertainty bounds.
- Ownership and preservation rules.
- Processing chemistry and contaminants.
- Packetized transport.
- Strategic reserves for millennia of losses.
- Alternatives if supply underperforms.

## 4. Oxygen industry

Oxygen extraction from regolith can use molten-salt electrolysis, molten-regolith electrolysis, hydrogen reduction, carbothermal reduction or other methods [S04]. Selection depends on energy, feedstock, metal coproducts and maintenance.

A coordinated industry produces:

- Atmospheric oxygen.
- Rocket oxidizer.
- Water oxygen.
- Iron, aluminum, titanium or mixed alloys.
- Glass and ceramics.
- Slag engineered for soil or construction.

The reference atmosphere requires processing about 1.5 × 10^18 kg of regolith at a 75% recovery assumption.

## 5. Water and hydrogen industry

Water systems combine:

- Protected extraction from selected polar deposits.
- Imported ice and hydrated materials.
- Hydrogen import and reaction with lunar oxygen.
- Industrial water recycling.
- Atmospheric water capture.
- Basin filling over centuries.

Hydrogen leakage is tracked isotopically. Every water-producing process includes a long-term replacement budget.

## 6. Carbon, phosphorus and sulfur

A biosphere needs far less carbon, phosphorus and sulfur than bulk atmosphere, but still orders of magnitude beyond early settlements. These materials arrive with volatile feedstocks or selected minerals and are recovered aggressively from wastes.

Carbon is stored in atmospheric CO2, soils, biomass, carbonates, fuels and strategic reserves. Phosphorus is especially important because it lacks a gaseous global reservoir; loss into inaccessible sediments can starve ecosystems.

## 7. Power architecture

A resilient system combines:

- Polar and equatorial solar arrays.
- Orbital solar-power stations.
- Fission for baseline and night power.
- Fusion if mature and economically real.
- Thermal storage in rock, salt and water.
- Chemical fuels and regenerative fuel cells.
- Superconducting or high-voltage continental grids.
- Local microgrids and black-start reactors.

Oxygen production alone can average hundreds of terawatts over a millennium in the reference screen. The complete program may demand petawatt-class industrial peaks or longer timelines.

## 8. Long-night energy storage

A 14.8-day night makes battery-only systems impractical for planetary industry. Storage uses multiple timescales:

- Seconds to hours: batteries, flywheels and grid controls.
- Hours to days: pumped water, compressed gas, thermal salts and chemical storage.
- Weeks: hydrogen, methane, ammonia, deep thermal reservoirs and nuclear generation.
- Years: strategic fuels and orbital power.

The atmosphere and hydrosphere themselves are thermal stores but cannot replace electrical generation for cities and control systems.

## 9. Mining geography

Industrial provinces are selected for:

- Favorable mineralogy.
- Proximity to power and transport.
- Low heritage and scientific value.
- Stable geology.
- Ability to contain dust and wastes.
- Access to future water and atmosphere systems.

Mining districts can be tens to hundreds of kilometers across. Restoration plans convert some pits into reservoirs, protected industrial museums, waste repositories or ecological basins.

## 10. Volatile transport chain

A safe chain is:

1. Survey and legally designate source bodies.
2. Extract and refine material near the source.
3. Divide it into independently controllable packets.
4. Move packets with high-efficiency propulsion or tethers.
5. Receive them in high Earth–Moon orbit.
6. Verify trajectory and composition.
7. Decelerate through propulsion, electromagnetic capture or controlled aerobraking.
8. Store in orbital or surface tanks.
9. Release only at the current pressure plateau.

No packet is large enough to create a civilization-ending impact if guidance fails.

## 11. Momentum and heat disposal

Terraforming logistics must close momentum, not only mass. Braking imported material produces heat and reaction momentum. Options include:

- Electric propulsion over long transfer times.
- Electromagnetic tethers exchanging momentum with reusable stages.
- Mass drivers ejecting inert counter-mass.
- Solar sails.
- Aerocapture after sufficient atmosphere exists.
- Export of manufactured products to balance flows.

The atmosphere itself cannot absorb arbitrary arrival energy without destructive heating.

## 12. Gas storage before release

A staged project stores vast gas inventories in:

- Subsurface caverns.
- Cryogenic tanks.
- Chemical forms such as ammonia, nitrates, oxides and carbonates.
- High-orbit depots.
- Enclosed regional atmospheres.

Chemical storage reduces pressure-vessel volume but adds processing energy. Strategic reserves remain isolated from the open atmosphere to recover from escape or contamination.

## 13. Atmospheric release system

Release nodes distribute gas gradually across latitudes and altitudes. They include:

- Composition metering.
- Isotope tracers.
- Emergency shutoff.
- Local weather constraints.
- Upper-atmosphere observation.
- Reversible capture where practical.

The program pauses at pressure plateaus rather than continuously racing toward the final atmosphere.

## 14. Planetary sensor network

Terluna requires continuous observations of:

- Surface pressure and composition.
- Vertical temperature and winds.
- Thermosphere and exosphere.
- Solar wind and ultraviolet flux.
- Ocean heat and salinity.
- Soil moisture and nutrient chemistry.
- Biomass and environmental DNA.
- Fire, methane and industrial emissions.
- Orbital drag and debris.

The data system is open, redundant and internationally audited.

## 15. Construction industry

Low gravity enables large spans, but atmospheric construction must resist wind. Major methods include:

- Sintered regolith foundations.
- Cast basalt and glass composites.
- Metal frames from oxygen-extraction coproducts.
- Buried pressure-safe districts.
- Large rotating artificial-gravity structures.
- Flexible dams and canal linings.
- Anchored towers and tether bases.

Early vacuum factories must either be moved, sealed or adapted before atmospheric release.

## 16. Space-access transition

A thick atmosphere ends the present concept of low lunar orbit. Orbital infrastructure moves outward before global pressure rises. The replacement network includes:

- High circular lunar orbits.
- Earth–Moon L1 and L2 logistics nodes.
- Other true Lagrange-point infrastructure where appropriate.
- Tethers extending above the dense atmosphere.
- Evacuated mountain mass drivers.
- High-altitude airports.
- Rocket corridors and spaceplanes.

Atmospheric bands are never named L1–L5.

## 17. Aviation logistics

Low weight enables very large aircraft and short takeoff, while the deep atmosphere enables routes at altitudes that would be near-space on Earth. Cargo systems can combine airship ports, railheads, mountain launch sites and sea transport.

Weather and icing may be more important limits than lift. The network needs alternate airports separated across environmental phases.

## 18. Communications and astronomy

The atmosphere, clouds and ionosphere alter radio and optical propagation. Infrastructure includes:

- Fiber and laser backbones.
- Stratospheric relays.
- Orbital communication constellations above drag.
- Farside radio-quiet governance.
- Adaptive-optics observatories.
- Orbital replacements for vacuum ultraviolet and X-ray astronomy.
- High-altitude farside astronomy platforms where atmospheric windows permit.

Terraforming sacrifices some unique surface-vacuum science, so replacements and archives must precede transformation.

## 19. Waste and circularity

At planetary scale, “waste” becomes geochemistry. The system tracks:

- Mining tailings.
- Salts from water processing.
- Nuclear materials.
- Persistent organic chemicals.
- Metals mobilized by soil weathering.
- Biological nutrients.
- Atmospheric pollutants.

Nutrient and volatile recovery is mandatory. Long-lived toxic compounds are restricted before they can enter the global water cycle.

## 20. Economic framing

A dollar cost is meaningless without specifying the future economy, energy price, automation, transport and project duration. More useful measures are:

- Kilograms imported per year.
- Joules per kilogram processed.
- Watts of continuous power.
- Machine-hours and replacement rates.
- Fraction of interplanetary industrial output.
- Atmospheric loss replacement fraction.
- Years of strategic volatile reserve.

The project is viable only when those flows are routine for its civilization.

## 21. Infrastructure success criteria

The industrial base passes when it can:

- Maintain multiple atmosphere plateaus without resource stress.
- Lose a major import convoy without catastrophic delay or impact.
- Black-start the planetary grid.
- Support sealed settlements independently of the open biosphere.
- Replace expected atmospheric losses for centuries.
- Preserve protected lunar regions.
- Operate under transparent multinational control.

<!-- END FILE: 08_INFRASTRUCTURE_INDUSTRY_AND_LOGISTICS.md -->


---

<!-- BEGIN FILE: 09_IMPLEMENTATION_ROADMAP.md -->

# Implementation Roadmap

## 1. Program principle

Terluna is not implemented by releasing a final atmosphere all at once. It advances through physically measured plateaus, each useful even if the next phase never happens. The architecture preserves a viable fallback to sealed cities and paraterraformed basins.

The durations below are scenarios, not forecasts. Automation, energy, transport and political legitimacy determine actual timing.

## 2. Phase 0 — Authority, archive and baseline

**Scenario duration:** 25–100 years.

### Objectives

- Establish internationally legitimate authority for irreversible lunar modification.
- Create a global geological, volatile, biological and heritage baseline.
- Sample and archive terrain that atmosphere, water and biology would alter.
- Define permanent vacuum, radio-quiet, heritage and polar-ice reserves.
- Resolve liability for cross-border climate and contamination harms.

### Deliverables

- High-resolution topography and subsurface maps.
- Global regolith and volatile sample archive stored on and off the Moon.
- Protected-zone treaty and enforcement system.
- Open environmental data platform.
- First-generation coupled climate and escape models.

### Stop gate

No large-scale atmosphere release or biological dispersal occurs without consent, protected areas, compensation mechanisms and an auditable baseline.

## 3. Phase 1 — Robotic industrial seed

**Scenario duration:** 50–200 years, overlapping later phases.

### Objectives

- Build power, mines, oxygen plants, metals production and transport.
- Demonstrate decade-scale autonomous maintenance.
- Create high-orbit volatile depots and safe packet capture.
- Construct sealed settlements and artificial-gravity research facilities.

### Deliverables

- Multi-terawatt then expanding power grid.
- Regolith oxygen plants with metal coproducts.
- Dust-controlled industrial districts.
- Evacuated mass-driver and tether demonstrators.
- Closed-cycle nutrient and water plants.

### Stop gate

Industry must survive multiple lunar cycles and major equipment failures without Earth-rescue dependence.

## 4. Phase 2 — Closed ecological cities

**Scenario duration:** 50–200 years of trials.

### Objectives

- Test soil, agriculture, animals and microbial consortia through many full lunar cycles.
- Determine artificial-gravity requirements.
- Study fire, smoke, rain and disease under lunar gravity.
- Operate communities at candidate pressures and gas mixtures.

### Required experiments

- Rotating and nonrotating multigenerational animal habitats.
- Large combustion chambers at 0.16 g.
- Kilometer-scale climate halls with crater and basin topography.
- Closed seas, aquifers and soil factories.
- Human residence with daily artificial-gravity protocols.

### Stop gate

No regional ecological release until food webs, nutrient cycles and health systems show predictable recovery from induced failures.

## 5. Phase 3 — Paraterraformed basins

**Scenario duration:** 100–400 years.

### Objectives

- Create roofed or membrane-contained regional landscapes.
- Operate open water, weather and aviation under containment.
- Test long-term ecosystem evolution at meaningful scale.
- Build social and legal experience governing shared atmospheres.

### Architecture

Candidate basins are capped by supported membranes, cable roofs, electrostatic dust barriers or hybrid structures. They can contain tens to thousands of square kilometers while using a small fraction of global gas mass.

### Stop gate

Regional climates must remain stable through infrastructure outages and deliberate perturbations. Failure cannot propagate globally.

## 6. Phase 4 — Thin global atmosphere

**Scenario duration:** 100–500 years.

**Pressure range:** roughly 0.5–5 kPa, initially unbreathable.

### Objectives

- Suppress some dust and moderate thermal extremes.
- Measure actual atmospheric escape, photochemistry and solar-wind interaction.
- Validate global wind and drag models.
- Begin the migration of orbital assets.
- Avoid biological release outside designated pilots.

### Instrumentation

The atmosphere is isotopically tagged and continuously observed from surface to exosphere. Pressure increases only after multiple solar cycles and extreme events are captured.

### Stop gate

Observed loss must agree with an affordable replenishment model. If loss is too high, global terraforming stops while enclosed basins continue.

## 7. Phase 5 — Pressure plateau ladder

**Scenario duration:** 200–1,000 years.

### Plateaus

- 10 kPa.
- 20 kPa.
- 40 kPa.
- 60 kPa.
- Optional 80 kPa reference target.

Each plateau lasts long enough to measure climate, escape, chemistry, orbital drag, erosion, dust, water stability and infrastructure effects. Gas is still dominated by inert buffer and controlled greenhouse species; oxygen remains limited until surfaces are passivated and fire systems mature.

### Stop gate at each plateau

- Escape and replacement budget closes.
- Climate has no unbounded trend.
- Orbital and launch systems are safe.
- Protected zones remain enforceable.
- Emergency gas capture or chemical storage retains meaningful capacity.

## 8. Phase 6 — Hydrosphere and soil expansion

**Scenario duration:** 200–1,500 years, overlapping pressure plateaus.

### Objectives

- Fill selected reservoirs and seas gradually.
- Establish rivers, wetlands, aquifers and controlled ice zones.
- Manufacture soil islands and connect them into corridors.
- Build carbon, nitrogen, phosphorus and sulfur reservoirs.

### Sequence

1. Industrial water reservoirs.
2. Sealed lakes.
3. Open basin lakes under nonbreathable atmosphere.
4. Regional seas.
5. Managed watersheds and rainfall.
6. Groundwater recharge.
7. Expanding soil and pioneer ecosystems.

### Stop gate

Water must not migrate irreversibly into cold traps, destabilize basins, contaminate protected ice or cause recurring global climate shocks.

## 9. Phase 7 — Controlled oxygenation

**Scenario duration:** 100–1,000 years.

### Objectives

- Raise oxygen in small steps.
- Form and monitor ozone.
- Passivate remaining reduced minerals.
- Introduce aerobic ecosystems.
- Validate combustion and corrosion controls.

### Oxygenation ladder

Oxygen partial pressure rises in increments rather than by fixed mole fraction. Each increment is held while soil sinks, wildfire behavior, atmospheric chemistry and biological productivity are measured.

### Stop gate

No further oxygen is added if regional fires become uncontrollable, oxygen sinks exceed production, ozone chemistry is dangerous, or ecosystems destabilize atmospheric carbon.

## 10. Phase 8 — Managed open biosphere

**Scenario duration:** centuries.

### Objectives

- Connect previously isolated soil and water ecosystems.
- Release pollinators, decomposers and selected animals.
- Establish migration and aerial-ecology corridors.
- Expand outdoor agriculture and settlement.
- Preserve sterile and low-biology reserves.

### Release discipline

Every organism has a documented genome, ecological function, geographic release boundary and recall or suppression plan. Global dispersal is assumed eventually possible, so only organisms acceptable at planetary scale are released.

### Stop gate

The biosphere must maintain bounded atmospheric and nutrient cycles through long observation periods without emergency planetary intervention.

## 11. Phase 9 — Permanent stewardship

**Duration:** indefinite.

Terraforming does not end when the air becomes breathable. Permanent functions include:

- Atmospheric replenishment and composition control.
- Climate and hydrological operations.
- Fire, pathogen and invasive-species response.
- Artificial-gravity health services.
- Protected-zone enforcement.
- Orbital-drag management.
- Continuous legitimacy review.

Terluna is a maintained world whose stewardship obligations last as long as the intervention’s consequences.

## 12. Cross-phase workstreams

### Atmospheric science

Runs from Phase 0 through permanent operations. Models are continuously corrected by pressure-plateau data.

### Human health

Begins before global atmosphere and may remain unresolved even after climate feasibility is proven.

### Governance

Precedes irreversible intervention and evolves as population and ecological complexity grow.

### Resource acquisition

Starts with small scientific and settlement flows, then scales only after the ethical and physical source inventory closes.

### Science preservation

Archiving and protected zones must occur before dust, weather, oxygen and water erase the original lunar environment.

## 13. Reversibility gradient

| Intervention | Approximate reversibility |
|---|---|
| Simulation and laboratory tests | Full |
| Sealed habitats | High |
| Enclosed basins | High to moderate |
| Thin atmosphere | Moderate, over long periods |
| Tens of kPa atmosphere | Low |
| Regional seas and global weathering | Low |
| Oxygenated open biosphere | Very low |
| Millennia of biological evolution | Effectively irreversible |

The burden of proof rises as reversibility falls.

## 14. Program no-go conditions

The global path stops if any of the following remains unresolved:

- Atmosphere loss consumes an unacceptable fraction of civilization’s resources.
- Nitrogen acquisition requires destructive appropriation.
- Climate models and pressure tests fail to maintain habitable refugia.
- 0.16 g is incompatible with healthy multigenerational life and artificial-gravity mitigation is impractical.
- Fire behavior cannot be bounded.
- Protected lunar science and heritage cannot be preserved adequately.
- International authority lacks legitimacy or enforcement.
- Regional pilots repeatedly collapse or contaminate reserves.

Stopping at paraterraforming is a successful fallback, not necessarily program failure.

## 15. Near-term research roadmap for a present-day creator

Even without planetary industry, the idea can advance through:

1. A canonical parameter registry.
2. A lunar topography-based basin and climate map.
3. A one-dimensional radiative-convective atmosphere model.
4. A first kinetic escape sensitivity study.
5. Low-gravity flight, rain and fire scaling experiments.
6. Long-photoperiod plant and ecosystem simulations.
7. A volatile-source mass-flow model.
8. Governance and protected-zone scenarios.
9. Visual maps and settlement concepts tied to the quantitative baseline.
10. Versioned decisions that clearly separate canon from provisional engineering.

<!-- END FILE: 09_IMPLEMENTATION_ROADMAP.md -->


---

<!-- BEGIN FILE: 10_RISKS_ETHICS_GOVERNANCE_AND_FAILURES.md -->

# Risks, Ethics, Governance and Failure Modes

## 1. Risk posture

Terluna combines planetary engineering with irreversible environmental transformation. Technical feasibility, moral legitimacy and political stability are separate requirements. Success in one cannot compensate automatically for failure in another.

## 2. Highest technical risks

### Atmospheric escape

A dense atmosphere may be lost slowly enough for routine replenishment or fast enough to make the project impossible. The uncertainty is decisive because the entire biosphere and economy depend on pressure.

### Nitrogen acquisition

The reference inventory requires about 1.35 × 10^18 kg nitrogen. A source can be physically accessible yet ethically unacceptable, scientifically protected or economically ruinous.

### Long-cycle climate instability

The atmosphere may fail to move heat effectively, or it may create excessive winds, clouds or greenhouse warming. Regional seas can moderate climate while also generating persistent storms.

### Partial-gravity biology

Humans, plants, animals and microbes may show multigenerational failures that short experiments miss. Artificial gravity can protect people but cannot rotate an entire outdoor biosphere.

### Fire and oxidation

Fresh regolith and infrastructure consume oxygen; later, dry biomass and unusual flame behavior can create large fires. Oxygenation is therefore a separate phase with its own gates.

## 3. Systemic ecological risks

- A pioneer organism becomes globally invasive.
- Nitrogen fixation or denitrification drives atmospheric imbalance.
- Decomposition consumes oxygen during late night.
- Methane or nitrous oxide causes unexpected warming.
- Pollinators or decomposers fail synchronously.
- Water ecosystems become anoxic under ice.
- Soil weathering mobilizes toxic metals.
- Pathogens cross between engineered species and humans.
- Evolution defeats genetic dependencies.

Ecological failures can be slow, nonlinear and difficult to reverse after global dispersal.

## 4. Human and settlement risks

- Chronic low-gravity disability.
- Developmental or reproductive harm.
- Wildfire and smoke stagnation.
- Basin pollution or cold-air pooling.
- Food failure during repeated long nights.
- Dependence on a centralized atmosphere authority.
- Inequality in access to 1 g facilities, safe climates and protected zones.
- Political coercion through air, water or climate control.

The right to breathable air and gravity-health mitigation must be treated as foundational civic infrastructure.

## 5. Space and Earth-system risks

- Imported packets become impactors.
- Atmosphere drag destroys satellites.
- Changed lunar albedo or thermal emission affects observations.
- Launch plumes contaminate Earth–Moon space.
- Large mass-flow systems create debris.
- Unmodeled changes affect precision navigation or astronomy.
- Earth-facing political conflict escalates over control of the Moon.

The added mass is too small to significantly change lunar gravity or orbit, but logistics and environment can still affect Earth–Moon operations.

## 6. Scientific loss

A global atmosphere would erase or compromise:

- Pristine vacuum exposure records.
- Surface solar-wind implantation.
- Some polar volatile archives.
- Unique dust and exosphere processes.
- Surface ultraviolet, X-ray and some radio astronomy advantages.
- Unweathered impact materials.
- Historical landing sites if not protected.

Archiving cannot preserve every spatial relationship or future research possibility. Some loss is irreducible.

## 7. Ethical question: who may transform the Moon?

The Moon is visible to and culturally significant for nearly every human society. Future lunar residents will bear direct consequences, while Earth populations, scientists and future generations hold legitimate interests.

A credible authority needs representation from:

- Lunar residents and workers.
- Earth states and publics.
- Scientific communities.
- Indigenous and cultural stakeholders with lunar traditions.
- Source-world or source-body stakeholders.
- Future-generation trustees.
- Independent ecological and safety institutions.

Technical ownership of equipment does not confer the moral right to alter a world globally.

## 8. International-law baseline

The Outer Space Treaty establishes state responsibility, due regard, consultation and avoidance of harmful contamination, but it was not written as a complete constitution for planetary terraforming [S19]. The Artemis Accords add current cooperative principles such as transparency, interoperability, emergency assistance, heritage protection and deconfliction, but they are likewise insufficient for global climate authority [S20].

Terluna would require new institutions governing atmosphere, biosphere, resources, liability, protected regions, migration and intergenerational obligations.

## 9. Required governance institutions

### Lunar Environmental Authority

Sets pressure, composition, release and protected-zone limits; publishes all data.

### Atmospheric Operations Service

Runs sensors and intervention systems under legally bounded authority.

### Planetary Court and Claims System

Handles transboundary damage, source-body disputes and compensation.

### Science and Heritage Trust

Controls protected sites and sample archives independently of industrial operators.

### Biosecurity Agency

Approves organisms, tracks genomes and commands quarantine or suppression.

### Future Generations Council

Can delay irreversible steps when evidence is incomplete.

### Resident Assembly

Gives inhabitants democratic control over daily environmental policy and emergency powers.

## 10. Consent model

Unanimity among all humanity is impossible, while simple majority rule is too weak for irreversible planetary change. A possible model requires concurrent supermajorities among:

- Lunar residents.
- Participating states or polities.
- Global population representation.
- Independent science and safety chambers.
- Protected-minority and future-generation institutions.

Each pressure plateau receives separate authorization. Consent to a thin atmosphere is not consent to an oxygenated biosphere.

## 11. Source-world ethics

Terraforming one world must not justify stripping another. Resource rules should prohibit:

- Destruction of potential ecosystems.
- Erasure of unique scientific records without compelling justification.
- Unilateral seizure of culturally significant bodies.
- Hazardous alteration of inhabited or potentially habitable worlds.
- Export that creates unacceptable environmental change at the source.

Distributed extraction from many low-value sources may reduce concentrated harm but increases monitoring complexity.

## 12. Irreversibility and burden of proof

The burden of proof increases with:

- Scale.
- Biological dispersal.
- Atmospheric pressure.
- Scientific loss.
- Dependence of human populations.
- Difficulty of restoring the previous state.

Early phases can use ordinary engineering risk tolerances. Global oxygenation demands a much stronger standard because rollback may take millennia and destroy the biosphere built around it.

## 13. Failure scenario: rapid atmospheric loss

### Sequence

Upper-atmosphere heating increases escape above projections; pressure declines faster than replenishment; water photolysis accelerates hydrogen loss; governments ration gas and close outdoor regions.

### Prepared response

- Halt releases and high-altitude leakage.
- Increase thermospheric cooling species if safe.
- Move people into sealed refuges.
- Protect seas with regional covers where possible.
- Convert atmospheric gases into chemical storage.
- Retreat to paraterraformed basins.

A strategic reserve must support this retreat.

## 14. Failure scenario: long-night ecological crash

### Sequence

Clouds reduce daytime productivity; plants enter night with inadequate stores; decomposition and animal respiration draw down local oxygen; food webs collapse late in night.

### Prepared response

- Emergency lighting and feeding in refugia.
- Regional oxygen injection.
- Reduced animal populations.
- More aquatic and fungal storage pathways.
- Revised canopy and cloud management.
- Preserve seed, embryo and microbiome banks.

## 15. Failure scenario: wildfire cascade

### Sequence

Drought, high oxygen fraction and strong winds create fires across connected biomes; low-gravity smoke behavior overwhelms evacuation routes.

### Prepared response

- Atmospheric and landscape compartments.
- Floodable firebreaks.
- Autonomous aircraft and orbital detection.
- Temporary regional oxygen dilution where technically possible.
- Sealed urban refuge networks.
- Post-fire erosion and nutrient recovery.

## 16. Failure scenario: imported-object accident

### Sequence

A volatile packet loses control and approaches an inhabited basin.

### Prepared response

- Packet mass capped below catastrophic thresholds.
- Multiple independent guidance and destruct systems.
- Off-plane receiving corridors.
- Ability to fragment, deflect or redirect to space.
- No direct surface delivery trajectories over inhabited regions.

## 17. Failure scenario: governance capture

### Sequence

One state, corporation or technical operator monopolizes atmosphere control and uses access to air, water or artificial gravity coercively.

### Prepared response

- Distributed physical control keys.
- Open-source operational models.
- Independent grid and gas reserves.
- Constitutional rights to life-support access.
- Transparent emergency actions with automatic review.
- Resident and international veto mechanisms.

## 18. Failure scenario: biosphere contamination of reserves

### Sequence

Spores, microbes or runoff cross a protected boundary and begin weathering pristine terrain.

### Prepared response

- Wide sterile buffer zones.
- Prevailing-wind and watershed-aware boundaries.
- Environmental DNA detection.
- Targeted sterilization or excavation where justified.
- Recognition that perfect recovery may be impossible.

## 19. Security risks

Atmosphere plants, climate controls, reservoirs, orbital depots and artificial-gravity hospitals are strategic targets. Security architecture must avoid creating a single planetary kill switch.

Resilience uses geographic distribution, offline controls, physical fail-safe states, independently governed reserves and strong prohibitions on environmental warfare.

## 20. Risk-ranking discipline

The CSV risk register separates impact and likelihood qualitatively. Future versions should add:

- Quantitative probability distributions.
- Exposure and consequence units.
- Leading indicators.
- Named owners.
- Mitigation readiness.
- Residual risk after controls.
- Cross-risk correlations.

Unknowns should not be scored as low likelihood merely because data are absent.

## 21. Ethical success criteria

Terluna is ethically defensible only if:

- It creates durable benefits that cannot be achieved more safely through contained habitats alone.
- Lunar residents freely participate in governance.
- Scientific and cultural losses are minimized and compensated.
- Source bodies are not destructively exploited.
- Future generations inherit functioning institutions and reserves.
- The program can stop when evidence turns unfavorable.
- No population is denied basic air, water or gravity-health protection for political reasons.

<!-- END FILE: 10_RISKS_ETHICS_GOVERNANCE_AND_FAILURES.md -->


---

<!-- BEGIN FILE: 11_RESEARCH_AND_VALIDATION_PROGRAM.md -->

# Research and Validation Program

## 1. Research objective

The program must turn Terluna from an internally consistent scenario into a sequence of testable models and demonstrations. Every major claim receives a verification path, an uncertainty range and a stop criterion.

## 2. Highest-priority question

**What is the loss rate and vertical structure of a 20–100 kPa lunar atmosphere across realistic compositions, solar activity and thermospheric temperatures?**

This question controls the economic feasibility of every open-world design. It requires a coupled model rather than extrapolation from the current exosphere or ancient millibar-to-kilopascal cases [S05][S06].

## 3. Atmosphere model hierarchy

### Level A — hydrostatic screening

- Pressure mass.
- Scale height.
- Lapse rate.
- Simple radiative balance.
- Condensation limits.

### Level B — one-dimensional radiative–convective chemistry

- O2/O3 photochemistry.
- CO2 and water profiles.
- Infrared cooling.
- Convective adjustment.
- Surface weathering sinks.

### Level C — global circulation model

- Lunar topography.
- 29.53-day moving solar forcing.
- Weak rotation.
- Clouds and precipitation.
- Regional seas, ice and soils.
- Dust and vegetation.

### Level D — upper-atmosphere and kinetic escape

- Thermosphere and exosphere.
- Solar ultraviolet and flares.
- Ionization, sputtering and pickup.
- Earth magnetotail passages.
- Multispecies diffusion.
- Hydrogen escape and water loss.

### Level E — coupled operations model

- Gas release and capture.
- Climate interventions.
- Orbital drag.
- Power and replenishment economics.
- Emergency scenarios.

Independent teams should implement competing models to expose structural uncertainty.

## 4. Atmospheric experiment ladder

1. Laboratory photochemistry and surface-reaction cells.
2. Lunar vacuum towers releasing trace tagged gases.
3. Sealed kilometer-scale atmosphere chambers at 0.16 g.
4. Large paraterraformed basins.
5. Sub-pascal global releases with isotope tracing.
6. Pascal-to-hundreds-of-pascals experiments.
7. Kilopascal pressure plateaus.
8. Multi-kilopascal climate demonstrations.

Each stage has predeclared recovery and termination procedures.

## 5. Climate questions

- Does the atmosphere develop superrotation?
- How strong are terminator fronts?
- Where do clouds persist during the long day?
- Which basins trap cold air or pollution?
- How much water freezes each night?
- Can deep seas and aquifers prevent climate collapse?
- What winds occur at 35–45 km?
- How sensitive is the climate to dust and vegetation?
- What is the full distribution of extremes, not just the mean?

## 6. Hydrology experiments

- Low-g drop, rain and snow microphysics.
- Open-channel flow and sediment transport at 0.16 g.
- Wave, tide and seiche behavior in crater basins.
- Ice growth and under-ice circulation through long nights.
- Gas exchange and oxygenation in deep low-g water.
- Fracture leakage and engineered aquifer sealing.
- Salt and nutrient concentration over repeated cycles.

Orbital centrifuge laboratories can provide continuous 0.16 g at scales unavailable in short parabolic flights.

## 7. Soil and geochemistry experiments

### Core measurements

- Reactive oxygen and metal species in fresh regolith.
- Oxidation and carbonate demand.
- Dust toxicity before and after weathering.
- Nutrient retention and leaching.
- Aggregate formation.
- Microbial weathering rates.
- Long-term root and fungal interactions.

### Scale progression

Gram-scale Apollo-like samples are insufficient. The program advances to tonnes, then landscapes, while preserving pristine control samples [S09][S16][S17].

## 8. Biosphere experiments

- Long-photoperiod plant growth for repeated 29.53-day cycles.
- Dual-clock gene expression.
- Night respiration and storage failure thresholds.
- Multispecies soil communities.
- Pollination and seed dispersal under low gravity.
- Aquatic ecosystems with long ice cover.
- Aerial organism flight and navigation.
- Fire ecology and post-fire recovery.
- Evolutionary escape from engineered dependencies.

The minimum useful ecological experiment spans many cycles and includes disturbances, not only stable operation.

## 9. Human research

### Partial gravity

- Bone and muscle dose-response across 0.16, 0.38, 0.5 and 1 g.
- Cardiovascular and ocular effects.
- Immune and microbiome changes.
- Pregnancy and developmental biology in animal models.
- Vestibular adaptation to mixed gravity.
- Artificial-gravity duration and rotation-rate tradeoffs.

### Environment

- Breathing and sleep at candidate pressures.
- Combustion products and smoke.
- Long-cycle lighting and mental health.
- Occupational exposure to treated soil and dust.
- Emergency decompression from an open atmosphere into shelters.

Human reproduction trials require exceptionally strong ethical safeguards and should follow extensive animal evidence.

## 10. Fire research

Fire is a critical independent gate. Required facilities include large 0.16 g centrifuge combustion laboratories capable of testing:

- Household and industrial materials.
- Forest and grass fuels.
- Flame spread at different O2 fractions and pressures.
- Smoke rise and stratification.
- Ember transport.
- Suppression droplets and foams.
- Basin-scale wind–fire interaction.

Small microgravity flames cannot validate landscape wildfire behavior.

## 11. Flight and weather research

- Wing and rotor performance in candidate atmospheres.
- Stall, gust and icing behavior.
- Airship structural design.
- Bird- and insect-like flight under 0.16 g.
- High-altitude turbulence and mountain waves.
- Emergency descent through a 40-km-deep weather column.

The 35–45 km platform concept is validated through long-duration balloons and aircraft at intermediate pressure plateaus.

## 12. Volatile resource program

The program creates a database of potential source materials with:

- Nitrogen, hydrogen, carbon, phosphorus and sulfur concentration.
- Total recoverable mass.
- Energy and delta-v.
- Scientific value.
- Contamination risk.
- Political status.
- Extraction and processing maturity.
- Delivery packet design.

The resource plan must show at least several times the final atmospheric inventory, including strategic reserve and losses.

## 13. Governance research

Technical models cannot answer legitimacy. Required work includes:

- Comparative constitutional designs for planetary utilities.
- Consent thresholds for irreversible phases.
- Rights of lunar residents versus Earth publics.
- Source-body environmental law.
- Liability for climate and biological damage.
- Protected-zone enforcement.
- Anti-monopoly design for atmosphere control.
- Future-generation representation.

Governance should be tested through simulations and real management of enclosed basins before global authority exists.

## 14. Verification metrics

### Atmosphere

- Pressure loss in Pa/year and fraction/year.
- Species-specific escape rates.
- Replenishment energy and mass fraction of economy.
- Composition variance by altitude and region.

### Climate

- Minimum and maximum inhabited-region temperatures.
- Ice fraction and water migration per cycle.
- Extreme wind and precipitation return periods.
- Ocean heat-content stability.

### Ecology

- Net primary production.
- Night storage margin.
- Oxygen and carbon budget closure.
- Nutrient leakage.
- Recovery time after disturbance.

### Human health

- Bone, muscle, cardiovascular, visual, immune and reproductive outcomes.
- Required artificial-gravity dose.
- Fire and air-quality mortality risk.

### Governance

- Auditability.
- Response time.
- Distribution of environmental benefits and harms.
- Protected-zone integrity.
- Ability to stop a phase legally and physically.

## 15. Model acceptance rules

- Publish inputs, source code and uncertainty.
- Compare at least two independent models.
- Calibrate against laboratory and pressure-plateau data.
- Use ensembles rather than one preferred run.
- Preserve failed models and explain rejection.
- Separate empirical values from design assumptions.
- Do not convert lack of evidence into a favorable default.

## 16. Initial computational projects

A present-day external workspace can begin with:

1. Python notebook reproducing mass, water and energy budgets.
2. One-dimensional radiative-convective model across 20–100 kPa.
3. Simplified longitudinal energy-balance model for the 29.53-day cycle.
4. Topographic basin-fill analysis using lunar elevation data.
5. Atmosphere loss sensitivity grid.
6. Volatile-source flow network.
7. Plant storage/respiration model across long light and dark periods.
8. Risk dependency graph linking hard gates.
9. Versioned canonical glossary and decision log.

## 17. Evidence grading

Every claim should carry one of these labels:

- **Recovered canon:** explicitly preserved from the project context.
- **Measured fact:** established lunar or physical measurement.
- **Published inference:** supported by a cited model or study.
- **Screening calculation:** simple calculation in this bundle.
- **Engineering hypothesis:** plausible design awaiting validation.
- **Speculative worldbuilding:** useful creative extension with low evidence.

This grading prevents attractive detail from being mistaken for established science.

<!-- END FILE: 11_RESEARCH_AND_VALIDATION_PROGRAM.md -->


---

<!-- BEGIN FILE: 12_QA_HANDOFF.md -->

# Terluna Q&A Handoff

Each question is self-contained. Answers are intentionally limited to one or two sentences for transfer into external tools.

## Q001. What is Terluna?

Terluna is a hypothetical transformation of Earth’s Moon into a maintained low-gravity world with breathable air, regional seas, biologically active soils, engineered ecosystems, permanent settlements and extensive atmospheric flight.

## Q002. Does the Terluna concept change the Moon’s orbit or gravity?

The reference concept keeps the Moon’s natural orbit, mass, synchronous rotation and 0.16 g surface gravity; it changes the surface environment rather than the Moon’s bulk mechanics.

## Q003. Why could a lunar atmosphere exist at all?

Hydrostatic gas can surround any body with gravity, and published work indicates ancient lunar volcanism probably produced a transient atmosphere near the kilopascal scale; the unresolved issue is how long a dense modern atmosphere would last [S05][S06].

## Q004. Is Terluna achievable with present-day technology?

No. It requires an interplanetary industrial economy able to move and process roughly 10^18 kg of volatiles, operate enormous power systems for centuries, and solve atmospheric escape, climate, ecology and partial-gravity health.

## Q005. What is the reference surface pressure for Terluna?

The bundle uses 80 kPa as a reconstructed reference case, with a broader 60–100 kPa open-world design range; the value is a modeling anchor rather than recovered canon.

## Q006. How much gas would an 80 kPa lunar atmosphere require?

A global 80 kPa atmosphere would contain about 1.87 × 10^18 kg of gas because lunar gravity is weak and requires a large mass column to produce pressure.

## Q007. What gases are in the reference Terluna atmosphere?

The dry reference mixture is 75% nitrogen, 24% oxygen, 0.9% argon, 800 ppm carbon dioxide and 200 ppm other controlled trace gases, with variable water vapor.

## Q008. Why does the Terluna reference atmosphere use 24% oxygen?

At 80 kPa, 24% oxygen provides a 19.2 kPa oxygen partial pressure near terrestrial physiological availability; combustion behavior still requires dedicated low-gravity testing.

## Q009. What is the largest atmospheric resource bottleneck?

Nitrogen is the dominant bottleneck because the reference atmosphere requires about 1.35 × 10^18 kg and the Moon has no demonstrated accessible reservoir near that scale.

## Q010. Could Terluna take its nitrogen from Earth?

Removing atmosphere-scale nitrogen from Earth would cause unacceptable environmental and political harm, so the concept requires non-Earth sources or a radically different low-pressure architecture.

## Q011. Where might Terluna obtain nitrogen?

Conceptual sources include ammonia- and nitrogen-bearing small bodies or outer-system resources; Venus and Titan contain large nitrogen inventories but raise severe transport, science and ethical issues.

## Q012. Could oxygen for Terluna come from lunar rock?

Yes in principle, because lunar regolith is rich in oxygen bound in oxides and laboratory extraction methods exist, although planetary-scale mining and energy requirements are extreme [S04].

## Q013. How much lunar regolith must be processed for the reference oxygen inventory?

At 45% contained oxygen and 75% recovery, producing about 4.93 × 10^17 kg of O2 requires processing roughly 1.46 × 10^18 kg of regolith.

## Q014. How much power would oxygen production require?

At an illustrative 20–50 MJ per kilogram of oxygen, oxygen production alone averages roughly 313–782 TW over 1,000 years, before most transport and construction loads.

## Q015. How long is a Terluna day?

The sunrise-to-sunrise solar cycle is about 29.53 Earth days, with an approximate 14.77-day illuminated half and 14.77-day dark half near the equator.

## Q016. How long is twilight on Terluna?

The canonical practical dawn or dusk transition is about five hours; it is an operational light interval inside the lunar solar cycle, while faint scattered light can extend beyond that threshold.

## Q017. Why was a previous 30-hour twilight estimate rejected?

The project canon was corrected to about five hours, so 30 hours should not be used unless describing a different very faint photometric threshold.

## Q018. How fast does the Terluna terminator move?

Near the equator, the sunrise or sunset line advances at roughly 4.3 m/s or 15 km/h, slowing by the cosine of latitude.

## Q019. Could organisms follow the Terluna terminator?

Some aircraft, flying animals or mobile communities could track favorable dawn or dusk conditions, but oceans, mountains, weather and land-use boundaries prevent a continuous simple route.

## Q020. Why might the long lunar night be survivable?

The reference atmosphere and hydrosphere have thermal capacities on the same order as a half-cycle’s absorbed solar energy, so heat storage and transport could prevent total freeze-out if clouds, circulation and water placement cooperate.

## Q021. Does the long night automatically freeze all Terluna water?

No. Deep seas, aquifers and insulated water beneath ice can remain liquid, although shallow water and exposed highlands may freeze and require climate management.

## Q022. What kind of hydrosphere does Terluna use?

The reference design uses regional seas, crater lakes, wetlands, rivers, aquifers and ice rather than an Earth-depth global ocean.

## Q023. How much water is in the reference Terluna hydrosphere?

The reconstructed reference inventory is 3 × 10^17 kg, equivalent to a global layer 7.9 m deep or about 79 m if concentrated over 10% of the surface.

## Q024. Does the Moon already contain enough water for Terluna?

The Moon contains scientifically important polar ice and other water signals, but known deposits are not a demonstrated source for a global hydrosphere [S03][S07][S08].

## Q025. How would Terluna receive imported water safely?

Water or hydrogen would be divided into controllable packets, refined and decelerated through orbital depots, propulsion, tethers or later aerocapture rather than delivered as giant impacts.

## Q026. Would Terluna have tides?

Earth’s tidal influence is mostly stationary on the nearside, while solar tides, libration and basin resonance can create slow multi-day or monthly water-level changes rather than familiar twice-daily tides.

## Q027. Would Terluna have rivers?

Yes. Water still flows downhill, but low gravity changes channel shape, sediment settling, flood propagation and erosion, so terrestrial river models require revalidation.

## Q028. Is lunar regolith usable as soil?

Raw lunar regolith is sharp, reactive, nutrient-poor mineral material rather than soil; plants have germinated in it but showed severe stress [S09][S16].

## Q029. How would Terluna manufacture soil?

Regolith would be screened, washed, weathered, chemically conditioned, amended with volatile nutrients and organic matter, inoculated with microbes and fungi, and cycled through pioneer crops.

## Q030. What happens to lunar dust after an atmosphere is added?

Electrostatic and ballistic dust behavior decreases, but ordinary windblown aerosols and dust storms become possible until rain, paving, soil crusts and vegetation stabilize the ground.

## Q031. How tall would the Terluna atmosphere be?

An Earthlike gas mixture at 285 K has a screening scale height near 50 km in lunar gravity, creating a weather and drag environment far deeper than Earth’s.

## Q032. Where is the top of Terluna’s main troposphere?

The project uses a provisional 35–45 km tropopause/lower-stratosphere range, pending a dedicated lunar general circulation model.

## Q033. What does “gravity waves” mean in Terluna weather?

It means atmospheric buoyancy waves generated by mountains, convection or terminator fronts, not gravitational waves in spacetime.

## Q034. Why should Terluna not label atmospheric layers L1 through L5?

L1–L5 conventionally name Lagrange points, so Terluna uses descriptive altitude bands to avoid confusing weather layers with orbital locations.

## Q035. What happens in the 35–45 km Terluna layer?

The concept envisions stratified research and flight zones with possible weak-shear windows, mountain and terminator gravity waves, sparse microbes, engineered algal films and long-duration science platforms.

## Q036. Would Terluna have strong winds despite low gravity?

Yes. Pressure gradients can drive strong winds, while buildings and trees weigh less and therefore need stronger anchoring relative to aerodynamic loads.

## Q037. Would Terluna have Earthlike jet streams?

Its slow rotation favors broader circulation and weaker familiar Coriolis organization, but high-altitude jets or atmospheric superrotation could still develop.

## Q038. Would rain fall slowly on Terluna?

For the same drop size, lower gravity reduces fall speed, while drops may grow larger before breakup; the net rain regime needs low-gravity microphysics experiments.

## Q039. Could Terluna have snow and ice storms?

Yes. Slow-falling snow, long-lived mixed-phase clouds and icing could be major hazards during dusk, night and dawn transitions.

## Q040. Would Terluna have an ozone layer?

An oxygenated atmosphere can form ozone, but its altitude, thickness and ultraviolet protection in a deep low-gravity atmosphere must be solved with photochemical models.

## Q041. How does atmosphere escape threaten Terluna?

Ultraviolet heating, sputtering, ion pickup, photochemistry and thermal escape can remove gases, especially light species; the rate determines whether replenishment is routine or prohibitive.

## Q042. Can Terluna reduce atmospheric escape?

Possible tools include thermospheric radiative cooling, composition optimization, magnetic or plasma interventions and continuous replenishment, but no unproven technology is assumed to solve the problem automatically.

## Q043. Would a dense Terluna atmosphere protect against radiation?

Its roughly 49,000 kg/m² mass column at 80 kPa should provide major shielding, but detailed particle-transport calculations are needed to account for secondary radiation [S10].

## Q044. Would the atmosphere destroy low lunar orbit?

A 50-km-scale-height atmosphere would create severe drag far above current low lunar orbits, so satellites and launch systems must move before pressure rises.

## Q045. How would spacecraft reach Terluna after terraforming?

The mature network uses high-altitude airports, evacuated mass drivers, tethers, controlled rocket corridors and higher lunar or Lagrange-point orbits.

## Q046. Does low gravity make aircraft easier to fly?

For the same aircraft and air density, stall speed scales with the square root of weight, giving an idealized lunar factor near 0.40 and enabling large payload or low-speed flight.

## Q047. Do airships gain six times the payload on Terluna?

No. Buoyant force and payload weight both scale with gravity, so displaced mass capacity is similar; airships mainly benefit from lower structural loads and a deep atmosphere.

## Q048. Why is aerial life important to Terluna?

Dense air and low weight create broad niches for soaring, gliding, ballooning, long-distance migration and high-altitude dispersal that are much harder on Earth.

## Q049. How would Terluna plants survive two weeks of daylight?

Engineered plants can use adjustable pigments, leaf folding, reflective surfaces, deep roots and internal rest cycles to prevent overheating and photodamage.

## Q050. How would Terluna plants survive two weeks of darkness?

They can store starches or oils, enter dormancy, shed leaves, use antifreeze chemistry and rely on bulbs, rhizomes or fungal networks until dawn.

## Q051. Would Terluna organisms use a month-long circadian rhythm?

Most would retain shorter internal rhythms while also using a second environmental clock for the 29.53-day growth, dormancy, migration and reproductive cycle.

## Q052. What animals are plausible in Terluna’s atmosphere?

Large soarers, cliff gliders, terminator migrants, ballooning organisms and aerial pollinators are plausible design directions, although ecology and physiology still limit size.

## Q053. Could Terluna support giant flying animals?

Low weight and dense air make larger fliers plausible, but muscle power, oxygen delivery, inertia, food supply, skeletal strength and safe launch still impose limits.

## Q054. What is the role of fungi in Terluna?

Fungi recycle stored biomass through the long night, build soil, move nutrients and support roots, but excessive late-night decomposition could consume dangerous amounts of oxygen.

## Q055. How is the Terluna biosphere introduced?

The sequence moves from sterile weathering to microbes, crusts, pioneer plants, wetlands and grasslands, complex aquatic systems, pollinators and only later large animals.

## Q056. Why must Terluna ecology remain managed?

Atmospheric oxygen, carbon, nitrogen, methane, water and nutrients can enter unstable feedbacks, so machines and institutions buffer the biosphere while it matures.

## Q057. Can engineered organisms be perfectly contained on Terluna?

No. Genetic dependencies reduce risk but mutation, horizontal transfer and physical dispersal prevent guaranteed permanent containment, so only planet-scale-acceptable organisms reach global release.

## Q058. Can humans live normally at 0.16 g?

The answer is unknown because no human has lived or reproduced for generations at lunar gravity; bone, muscle, cardiovascular, visual, immune and developmental risks remain critical [S14][S15].

## Q059. How does Terluna protect residents from low gravity?

Cities include rotating sleep habitats, clinics, maternity centers, schools and exercise facilities that provide regular partial or full Earth gravity.

## Q060. How large must a rotating 1 g habitat be?

A 250 m radius ring produces 1 g at about 1.9 rpm, while larger radii reduce rotation rate and gravity gradients.

## Q061. Why would Terluna buildings need strong anchors?

Their weight is one sixth of Earth weight while wind pressure can be substantial, so tall or broad structures can overturn unless tied deeply into rock.

## Q062. How would people keep time on Terluna?

Residents can use 24-hour civil days for sleep and work while tracking a 29.53-day environmental month divided into dawn, day, dusk and night phases.

## Q063. What would nearside Terluna look like at night?

Earth would remain nearly fixed in the sky and appear bright when the Sun is below the lunar horizon, creating culturally important Earthlight that is too weak to power most ecosystems.

## Q064. What makes farside Terluna distinct?

The farside has no Earth in its sky and can preserve valuable radio-quiet regions, making it culturally and scientifically different from the nearside.

## Q065. Would every Terluna settlement still need sealed shelters?

Yes. Sealed refuges protect against fire, contamination, pressure loss, pathogens, flood, freeze and grid failure even after the open atmosphere becomes reliable.

## Q066. What are paraterraformed basins?

They are large enclosed landscapes with their own atmospheres, water and ecosystems, providing most open-world benefits with far less gas and greater reversibility.

## Q067. Why is paraterraforming the main fallback?

It remains viable even if global atmospheric escape, nitrogen supply, governance or health makes an open lunar atmosphere unacceptable.

## Q068. What is the first irreversible Terluna step?

Large global gas release and biological dispersal begin the steep loss of reversibility, which is why archive, consent and protected-zone work must come first.

## Q069. What pressure plateaus does the roadmap use?

After a thin experimental atmosphere, the roadmap tests approximately 10, 20, 40 and 60 kPa before any optional 80 kPa final target.

## Q070. Why is oxygen added late?

Early oxygen would be consumed by fresh minerals, accelerate corrosion and create fire risk before soils, hydrology and combustion controls are ready.

## Q071. What could cause the Terluna project to stop permanently?

Unmanageable escape, unethical volatile sourcing, unstable climate, unsafe partial gravity, uncontrollable wildfire, ecological collapse or illegitimate governance can each end the global path.

## Q072. Would stopping at enclosed habitats mean Terluna failed?

No. A network of sealed cities and paraterraformed basins could deliver habitable lunar landscapes without accepting the irreversible costs of a global atmosphere.

## Q073. What scientific resources would global terraforming destroy?

Weathering and biology would compromise pristine vacuum exposure records, polar volatile archives, unaltered regolith, some astronomy and historical sites unless protected or archived.

## Q074. What law currently authorizes Terluna?

No current legal framework is sufficient for planetary terraforming; the Outer Space Treaty and Artemis Accords provide principles that would need major new institutions [S19][S20].

## Q075. Who should govern Terluna’s atmosphere?

A legitimate system should combine lunar resident control, international representation, independent science and biosecurity bodies, future-generation trustees and transparent distributed operations.

## Q076. Why is atmosphere control politically dangerous?

An operator that controls pressure, oxygen, water and artificial gravity can coerce populations, so no single state or corporation should hold a planetary kill switch.

## Q077. How should imported volatile accidents be limited?

Delivery packets must be small enough that one guidance failure cannot cause civilization-ending damage, with independent navigation, deflection and off-population receiving corridors.

## Q078. What is Terluna’s main climate validation test?

A coupled atmosphere–ocean–land model and staged pressure experiments must maintain habitable refugia through repeated lunar cycles and solar extremes without unbounded trends.

## Q079. What is Terluna’s main biological validation test?

Closed and regional ecosystems must survive repeated full lunar cycles, disturbances and late-night resource stress while closing oxygen, carbon, nitrogen and water budgets.

## Q080. What is Terluna’s main human-health validation test?

Multigenerational partial-gravity evidence must establish whether daily artificial gravity can prevent developmental, reproductive and chronic health harm.

## Q081. What is Terluna’s main governance validation test?

Institutions must demonstrate transparent control, equitable access to life support, enforceable protected zones and the legal ability to halt a pressure or biology phase.

## Q082. What does this bundle treat as recovered canon?

Recovered elements include the name Terluna, a long lunar light-dark cycle, about five hours of practical twilight, extensive aerial ecology, a 35–45 km upper weather layer and the terminology corrections on L1–L5 and gravity waves.

## Q083. Which major Terluna numbers are reconstructed rather than canon?

The 80 kPa pressure, 24% oxygen, 285 K mean target and 3 × 10^17 kg hydrosphere are quantitative baselines introduced for modeling and can be revised.

## Q084. What should an external workspace preserve when extending Terluna?

It should preserve the recovered canon and corrections, label new assumptions explicitly, keep the atmosphere-loss and nitrogen gates visible, and avoid presenting speculative biology or climate as established fact.

## Q085. What is the overall feasibility verdict for Terluna?

No known physical law rules out a maintained lunar atmosphere and biosphere, but the concept remains a far-future planetary engineering program whose critical resource, escape, climate, health and governance questions are unresolved.

<!-- END FILE: 12_QA_HANDOFF.md -->


---

<!-- BEGIN FILE: 13_EXTERNAL_WORKSPACE_CONTEXT.md -->

# External Workspace Context

## 1. Ready-to-paste project context

Copy the block below into an external workspace, AI project, wiki or research notebook.

--- BEGIN TERLUNA CONTEXT ---

**Project:** Terluna — terraforming Earth’s existing Moon.

**Core premise:** Terluna is a maintained open lunar world with breathable air, regional seas, manufactured soils, engineered ecosystems, permanent human settlements and extensive atmospheric flight. The Moon retains its existing orbit, 0.16 g surface gravity, synchronous rotation and 29.53-Earth-day sunrise-to-sunrise cycle. The project adapts biology and infrastructure to the natural long cycle rather than forcing a 24-hour planetary rotation.

**Recovered canon and corrections:**

1. The world is called **Terluna**.
2. Organisms adapt to prolonged daylight and darkness through resource storage, dormancy, migration and dual timing systems.
3. Aerial and gliding life is a major ecological feature.
4. The canonical practical dawn or dusk transition is about **five hours**, not thirty; those intervals are inside the 29.53-day cycle.
5. The main upper-troposphere/lower-stratosphere worldbuilding band is approximately **35–45 km**.
6. That high layer may contain weak-shear windows, atmospheric gravity waves launched by mountains or terminator fronts, sparse microbes, engineered radiation-hard algal films, and platforms for farside astronomy or nearside Earthlight photometry.
7. Do not label atmospheric layers L1–L5 because those labels mean Lagrange points.
8. “Gravity waves” means atmospheric buoyancy waves, not gravitational waves.

**Reference design values introduced for modeling, not fixed canon:**

- Surface pressure: 80 kPa, with a 60–100 kPa open-world range.
- Dry atmosphere: 75% N2, 24% O2, 0.9% Ar, 800 ppm CO2, 200 ppm other trace gases.
- Oxygen partial pressure: 19.2 kPa.
- Mean surface temperature target: 285 K.
- Water inventory: 3 × 10^17 kg, forming regional seas, lakes, rivers, wetlands, groundwater and ice.
- Main tropopause region: 35–45 km pending a dedicated climate model.

**Key quantitative consequences:**

- An 80 kPa lunar atmosphere requires about 1.87 × 10^18 kg of gas.
- The reference nitrogen mass is about 1.35 × 10^18 kg, making nitrogen the largest resource bottleneck.
- The reference oxygen mass is about 4.93 × 10^17 kg.
- Lunar regolith is roughly 45% oxygen by mass; at 75% recovery, the oxygen inventory requires processing about 1.46 × 10^18 kg of regolith.
- At 20–50 MJ/kg O2, oxygen production alone averages roughly 313–782 TW over 1,000 years.
- The atmospheric scale height is about 50 km at 285 K, and the dry adiabatic lapse rate is about 1.62 K/km.
- The 80 kPa column mass is about 49,300 kg/m², giving major shielding but also a very deep drag environment.
- The 3 × 10^17 kg hydrosphere is a 7.9 m global equivalent layer or about 79 m over 10% of the surface.
- The equatorial terminator moves about 4.3 m/s.

**Why the concept may be physically possible:**

- Gravity can support a collisional atmosphere, and ancient lunar volcanism probably produced a transient atmosphere near the kilopascal scale.
- A dense atmosphere is not necessarily lost instantly; lifetime depends strongly on composition, upper-atmosphere temperature and solar forcing.
- Terluna receives Earthlike solar flux.
- A dense atmosphere and regional seas provide heat transport and thermal storage across the long night.
- Lunar rock supplies oxygen and construction materials.
- Low gravity strongly favors aviation, gliding and large structures.

**Why the concept may fail:**

- Dense-atmosphere escape could demand unaffordable replenishment.
- Nitrogen and water imports are order-10^18 and order-10^17–10^18 kg problems.
- The long day-night climate may produce severe storms, freeze-out or overheating.
- Human health and reproduction at 0.16 g are unknown.
- Low-gravity fire behavior may make an oxygenated biosphere unsafe.
- Raw regolith is toxic, abrasive and infertile until extensively processed.
- A thick atmosphere eliminates current low lunar orbit and some unique lunar science.
- Global terraforming is politically and ecologically difficult to reverse.

**Biosphere design:**

- Use staged succession: physical weathering; microbes and biofilms; crusts and pioneer plants; fungi and detrital networks; grasslands, shrubs and wetlands; complex aquatic systems; pollinators and small animals; large animals only after long stability.
- Organisms retain shorter internal rhythms while treating the lunar environmental month as a growth/dormancy cycle.
- Plants store starches and oils, use bulbs and rhizomes, adjust pigments, fold leaves and enter night dormancy.
- Animals use fat, torpor, caches, burrows, migration and communal roosts.
- Deep water and under-ice habitats are major night refugia.
- Low gravity and dense air permit large soaring, gliding, ballooning and terminator-migrating organisms.

**Human design:**

- Outdoor air does not solve 0.16 g health.
- Cities include rotating sleep quarters, schools, clinics, maternity facilities and exercise habitats.
- Buildings need unusually strong anchoring because wind loads remain large while structural weight is low.
- Residents keep approximately 24-hour civil days within a 29.53-day environmental month.
- Nearside cultures live under a nearly fixed Earth; farside cultures have no Earth in their sky and may protect radio-quiet science.
- Every settlement retains sealed independent refuges.

**Infrastructure design:**

- Build lunar industry before global terraforming.
- Extract oxygen and metals from regolith; import most nitrogen, hydrogen, carbon and other volatile nutrients.
- Use packetized delivery through orbital depots; never deliver the inventory as giant impacts.
- Combine solar, nuclear, orbital power, thermal stores and chemical storage.
- Transition satellites out of low lunar orbit before pressure rises.
- Use high orbits, actual Earth–Moon Lagrange points, tethers, evacuated mass drivers, mountain launch sites and high-altitude aviation.

**Roadmap:**

0. Authority, archive and protected zones.
1. Robotic industrial seed.
2. Closed ecological cities and artificial-gravity research.
3. Paraterraformed enclosed basins.
4. Thin global atmosphere for escape and climate measurements.
5. Pressure plateaus around 10, 20, 40 and 60 kPa.
6. Hydrosphere and soil expansion under mostly inert atmosphere.
7. Controlled oxygenation.
8. Managed open biosphere.
9. Permanent stewardship.

**Hard no-go gates:**

- Dense-atmosphere escape and replenishment do not close.
- Nitrogen acquisition is destructive or infeasible.
- Long-cycle climate cannot maintain safe refugia.
- Human partial-gravity harm cannot be mitigated.
- Fire behavior is uncontrollable.
- Protected science and heritage cannot be preserved adequately.
- International and resident governance lacks legitimacy.
- Regional pilots repeatedly fail.

**Fallback:** A network of sealed cities and paraterraformed basins remains the default fallback and may be the rational endpoint even if global terraforming is technically possible.

**Evidence discipline:** Mark every addition as recovered canon, measured fact, published inference, screening calculation, engineering hypothesis or speculative worldbuilding. Do not silently convert the reference 80 kPa atmosphere or 3 × 10^17 kg hydrosphere into immutable canon.

--- END TERLUNA CONTEXT ---

## 2. Suggested system prompt for an AI workspace

Use the following instruction after the context block:

> Develop Terluna as rigorous far-future planetary engineering and internally consistent worldbuilding. Preserve recovered canon and corrections. Quantify material, energy, thermal, ecological and maintenance budgets. State uncertainty and identify hard gates. Treat the 80 kPa atmosphere and other reference numbers as revisable baselines. Never confuse atmospheric bands with Lagrange points, never call atmospheric gravity waves gravitational waves, and retain the canonical five-hour practical twilight. Prefer staged, reversible designs and preserve paraterraforming as a fallback.

## 3. Suggested workspace structure

- `Canon` — fixed names, corrected terms and explicit creator decisions.
- `Parameters` — values, ranges, units, status and confidence.
- `Physics` — atmosphere, climate, escape, hydrology and orbital effects.
- `Biosphere` — succession, traits, food webs, disease and fire.
- `Human systems` — health, artificial gravity, settlements and culture.
- `Infrastructure` — resources, power, transport and space access.
- `Roadmap` — phases, gates, owners and validation evidence.
- `Risks` — linked risk register and failure scenarios.
- `Sources` — primary literature and official references.
- `Decision log` — every change that affects canon or reference design.

## 4. Machine-ingestion files

Use `data/terluna_context.json` for the core context, `data/qa_handoff.jsonl` for retrieval or fine-grained import, and the CSV files for databases or spreadsheets. The Markdown files preserve the full reasoning and caveats.

<!-- END FILE: 13_EXTERNAL_WORKSPACE_CONTEXT.md -->


---

<!-- BEGIN FILE: 14_GLOSSARY.md -->

# Glossary

## Active stewardship

Permanent monitoring and intervention required to keep Terluna’s atmosphere, water and biosphere within safe bounds.

## Aerial ecology

The organisms and food webs that live, migrate, feed or reproduce substantially in the atmosphere.

## Aerocapture

Using atmospheric drag to reduce an arriving vehicle’s speed without relying entirely on propulsion.

## Agglutinate

Glassy welded lunar soil particle created by micrometeorite impacts; a contributor to regolith reactivity and abrasion.

## Artificial gravity

Apparent gravity produced by rotation or acceleration, used to mitigate chronic 0.16 g exposure.

## Atmospheric band

A descriptive Terluna altitude region such as lower troposphere or tropopause; never called L1–L5.

## Atmospheric gravity wave

A buoyancy-restored oscillation in a stratified atmosphere, often generated by mountains or convection.

## Atmospheric plateau

A deliberately held pressure level used to collect data before adding more gas.

## Atmospheric superrotation

Atmospheric circulation that moves around a world faster than the solid surface rotates.

## Biosecurity

Institutions and techniques that control pathogens, invasive species and engineered-organism release.

## Bond albedo

Fraction of total incident solar energy reflected by a world across all wavelengths and directions.

## Buffer gas

A mostly nonreactive gas that supplies pressure and moderates oxygen fraction; nitrogen is the reference choice.

## Canon

Project information explicitly fixed by the creator or recovered reliably from prior context.

## Carbon sink

A reservoir or process that removes carbon dioxide or other carbon from the atmosphere.

## Circadian rhythm

An internal biological cycle near one Earth day; Terluna organisms may retain it alongside a longer environmental clock.

## Closed ecological system

A materially constrained habitat that recycles air, water and nutrients without open exchange with a planetary biosphere.

## Cold trap

A persistently cold region where vapor condenses or freezes and can accumulate over time.

## Column mass

Mass of atmosphere above one square meter, equal to surface pressure divided by gravity.

## Coriolis effect

Apparent deflection of moving air or water in a rotating frame; weak on the slowly rotating Moon.

## Dry adiabatic lapse rate

Temperature change of an unsaturated rising or sinking parcel due to pressure change, approximately g/cp.

## Dual clock

Biological timing architecture combining a short internal rhythm with the 29.53-day environmental cycle.

## Earthlight

Sunlight reflected from Earth onto the lunar nearside, strongest during local lunar night.

## Engineering hypothesis

A plausible design claim that lacks direct validation and must not be presented as measured fact.

## Environmental month

One 29.53-Earth-day Terluna sunrise-to-sunrise cycle.

## Exobase

Altitude above which atmospheric particles can travel long distances without collisions; the transition to kinetic escape.

## Exosphere

The collision-poor outer atmosphere where particle trajectories are governed strongly by gravity and escape processes.

## Farside

The lunar hemisphere that normally faces away from Earth.

## Firebreak landscape

Water, mineral ground, low-flammability vegetation and access corridors designed to stop wildfire spread.

## Fixed nitrogen

Chemically reactive nitrogen such as ammonium or nitrate that organisms can use, unlike atmospheric N2 without fixation.

## General circulation model

A numerical model that simulates planetary atmosphere and often ocean, land, clouds and radiation in three dimensions.

## Global equivalent depth

Depth obtained by spreading a liquid inventory evenly across an entire world’s surface.

## Greenhouse gas

A gas that absorbs and emits infrared radiation, changing planetary heat loss.

## Hydrodynamic escape

Bulk atmospheric outflow driven by strong heating, capable of carrying multiple species.

## Hydrosphere

All water in seas, lakes, rivers, groundwater, ice and atmosphere.

## Isotope tracer

A distinctive isotopic composition used to identify gas sources, transport and loss.

## Jeans escape

Thermal loss of particles in the high-speed tail of a velocity distribution above the escape speed.

## Lagrange point

One of five positions in a two-body rotating system where gravitational and orbital effects permit special spacecraft behavior; labeled L1–L5.

## Libration

Apparent rocking of the Moon as seen from Earth, causing Earth to move slightly in the Terluna nearside sky.

## Mass driver

Electromagnetic launcher that accelerates cargo without onboard propellant, ideally inside an evacuated tube after atmosphere.

## Measured fact

A value or phenomenon established by observation or experiment rather than chosen for the design.

## Nearside

The lunar hemisphere that normally faces Earth.

## Nitrogen fixation

Conversion of atmospheric N2 into biologically usable compounds.

## No-go gate

A condition that must be satisfied before a program phase can proceed.

## Oxygen partial pressure

The fraction of total pressure contributed by oxygen, a major determinant of respiration and combustion.

## Paraterraforming

Creating large enclosed habitable landscapes without changing an entire world’s atmosphere.

## Partial gravity

Gravity between microgravity and Earth gravity; the Moon provides about 0.16 g.

## Photochemical escape

Atmospheric escape caused by chemical reactions that give particles enough energy to leave.

## Photometric twilight

Twilight defined by a measured light threshold rather than only solar angle.

## Pressure datum

Reference elevation at which a planet’s nominal surface pressure is specified.

## Published inference

A claim supported by a cited scientific model or interpretation but not necessarily directly measured.

## Radiative-convective model

A vertical climate model balancing radiation and convective heat transport.

## Recovered canon

Project content clearly present in accessible prior context.

## Reference design

A quantified baseline used for comparison, not an immutable final choice.

## Regolith

Loose fragmented rock, dust and glass covering the lunar surface.

## Resource storage organ

Biological structure such as a bulb, tuber or rhizome that stores energy and nutrients through the long night.

## Rossby deformation length

Scale at which rotation becomes important to atmospheric or oceanic flow; expected to be very large on Terluna.

## Scale height

Characteristic altitude over which pressure falls by a factor of e in an isothermal atmosphere.

## Screening calculation

A simple quantitative estimate used to identify scale and feasibility before detailed modeling.

## Seiche

Standing oscillation of water in an enclosed or partly enclosed basin.

## Sidereal rotation period

Rotation period relative to distant stars, about 27.3 Earth days for the Moon.

## Solar day

Time between successive local noons or sunrises, about 29.53 Earth days on the Moon.

## Solar-wind sputtering

Ejection of atmospheric or surface particles caused by energetic charged-particle impacts.

## Speculative worldbuilding

Creative extension that is physically motivated but has low evidentiary confidence.

## Synodic cycle

Cycle relative to the Sun; for the Moon it is about 29.53 Earth days.

## Terminator

Boundary between illuminated and dark terrain.

## Terminator front

Weather response associated with the moving dawn or dusk forcing zone.

## Thermosphere

Upper atmospheric region where absorption of high-energy radiation raises kinetic temperature.

## Thermal inertia

Resistance to temperature change due to heat capacity and heat transport.

## Tidal locking

Rotation state in which the same hemisphere generally faces an orbital partner.

## Tropopause

Boundary region between the convecting troposphere and the more stably stratified atmosphere above.

## Upper-atmosphere closure

Demonstration that thermospheric state and atmospheric loss yield an affordable replenishment burden.

## Volatile

Material such as water, nitrogen, carbon dioxide or ammonia that vaporizes readily and is scarce in bulk lunar rock.

## Weathering sink

Consumption of atmospheric gases through reactions with fresh rock and soil.

<!-- END FILE: 14_GLOSSARY.md -->


---

<!-- BEGIN FILE: 15_SOURCE_NOTES.md -->

# Source Notes

## Citation method

Narrative files cite sources by identifiers such as `[S05]`. Full URLs are also stored in `data/source_manifest.csv`. Sources establish lunar facts, published precedents and current legal principles; they do not validate the complete Terluna reference design.

## S01 — NASA Moon Facts

**Publisher:** NASA Science.  
**URL:** https://science.nasa.gov/moon/facts/

Used for lunar size, rotation, orbital and phase-cycle context, current near-vacuum exosphere, nearside/farside geometry and weak axial tilt. The Terluna solar cycle uses a more precise standard synodic value in the calculations.

## S02 — NASA Compare Earth and the Moon

**Publisher:** NASA Science.  
**URL:** https://science.nasa.gov/moon/by-the-numbers/

Used for surface gravity and escape-speed reference values. Small variation in published rounded constants does not affect order-of-magnitude conclusions.

## S03 — NASA Moon Water and Ices

**Publisher:** NASA Science.  
**URL:** https://science.nasa.gov/moon/moon-water-and-ices/

Summarizes evidence for polar ice and water detections from missions including LCROSS and Moon Mineralogy Mapper. It supports the existence of indigenous water, not a claim that enough recoverable ice exists for Terluna’s hydrosphere.

## S04 — ESA Turning Moon Dust into Oxygen

**Publisher:** European Space Agency.  
**URL:** https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/Turning_Moon_dust_into_oxygen

Used for the approximate 45% oxygen-by-mass character of lunar regolith and a molten-salt electrolysis demonstration. The 75% recovery and 20–50 MJ/kg system energy used here are screening assumptions, not values asserted by ESA for planetary production.

## S05 — Needham and Kring (2017)

**Title:** Lunar volcanism produced a transient atmosphere around the ancient Moon.  
**Journal:** Earth and Planetary Science Letters.  
**DOI:** https://doi.org/10.1016/j.epsl.2017.09.002

The paper estimates that volcanic outgassing could have produced a transient lunar atmosphere peaking around the kilopascal scale and persisting for tens of millions of years. It provides precedent for noninstantaneous lunar atmospheres, not a direct model of dense N2–O2 Terluna.

## S06 — Tucker et al. (2021)

**Title:** Lifetime of a transient atmosphere produced by lunar volcanism.  
**Journal:** Icarus.  
**DOI:** https://doi.org/10.1016/j.icarus.2021.114304

Used to show that atmospheric lifetime depends strongly on upper-atmosphere temperature and extreme-ultraviolet energy. The modeled ancient atmosphere and solar environment differ from Terluna, so extrapolation is deliberately limited.

## S07 — Li et al. (2018)

**Title:** Direct evidence of surface exposed water ice in the lunar polar regions.  
**Journal:** Proceedings of the National Academy of Sciences.  
**DOI:** https://doi.org/10.1073/pnas.1802345115

Supports direct spectral evidence for polar surface ice. Recoverability, continuity and preservation value remain separate questions.

## S08 — Honniball et al. (2021)

**Title:** Molecular water detected on the sunlit Moon by SOFIA.  
**Journal:** Nature Astronomy.  
**DOI:** https://doi.org/10.1038/s41550-020-01222-x

Supports molecular-water detections in observed sunlit high-latitude material at roughly hundreds of micrograms per gram. This is not a bulk water reserve.

## S09 — Paul et al. (2022)

**Title:** Plants grown in Apollo lunar regolith present stress-associated transcriptomes that inform prospects for lunar exploration.  
**Journal:** Communications Biology.  
**DOI:** https://doi.org/10.1038/s42003-022-03334-8

Plants germinated in real Apollo regolith but displayed strong stress and impaired growth. The study supports the claim that raw regolith interacts with plants but is not ready agricultural soil.

## S10 — Zhang et al. (2020)

**Title:** First measurements of the radiation dose on the lunar surface.  
**Journal:** Science Advances.  
**DOI:** https://doi.org/10.1126/sciadv.aaz1334

Provides direct radiation-dose context for the unmodified lunar surface. Terluna’s proposed atmospheric shielding requires separate transport calculations.

## S11 — Way and Del Genio (2020)

**Title:** Venusian habitable climate scenarios: Modeling Venus through time and applications to slowly rotating Venus-like exoplanets.  
**Journal:** Journal of Geophysical Research: Planets.  
**DOI:** https://doi.org/10.1029/2019JE006276

Provides a modern GCM example of slowly rotating terrestrial climate behavior. It does not model lunar gravity or Terluna’s moving 29.53-day forcing.

## S12 — Yang, Cowan and Abbot (2013)

**Title:** Stabilizing cloud feedback dramatically expands the habitable zone of tidally locked planets.  
**Journal:** Astrophysical Journal Letters.  
**DOI:** https://doi.org/10.1088/2041-8205/771/2/L45

Supports the mechanism by which persistent reflective clouds can limit day-side warming in slow or tidally locked climates. Terluna’s rotation, hydrology and gravity require a dedicated model.

## S13 — Colyer and Vallis (2018)

**Title:** Zonal-mean atmospheric dynamics of slowly rotating terrestrial planets.  
**Preprint:** https://arxiv.org/abs/1806.10494

Used for broad-overturning and weak-rotation dynamics. It offers theoretical context rather than a lunar forecast.

## S14 — NASA Human Research Program: Hazard — Gravity Fields

**Publisher:** NASA.  
**URL:** https://www.nasa.gov/hrp/hazard-gravity-fields/

Summarizes human risks from altered gravity. It supports caution but does not provide a known safe minimum gravity or artificial-gravity dose.

## S15 — Du et al. (2025)

**Title:** Lunar and Martian gravity alter immune cell interactions.  
**Journal:** npj Microgravity.  
**DOI:** https://doi.org/10.1038/s41526-024-00456-7

Provides partial-gravity cellular evidence relevant to immune uncertainty. It is not a multigenerational human outcome study.

## S16 — NASA What Is Lunar Regolith?

**Publisher:** NASA Science.  
**URL:** https://science.nasa.gov/biological-physical/what-is-lunar-regolith/

Used for the physical nature of fragmented sharp lunar material. Soil-development claims in the bundle remain engineering proposals.

## S17 — Linnarsson et al. (2012)

**Title:** Toxicity of lunar dust.  
**Journal:** Planetary and Space Science.  
**DOI:** https://doi.org/10.1016/j.pss.2012.05.013

Reviews health uncertainties and the need for realistic lunar-dust characterization. A future weathered Terluna soil would differ from pristine dust, requiring new toxicology.

## S18 — Crawford (2015)

**Title:** Lunar resources: A review.  
**Journal:** Progress in Physical Geography.  
**Preprint:** https://arxiv.org/abs/1410.6865

Provides general lunar-resource context and limitations. It does not identify a lunar nitrogen inventory sufficient for Terluna.

## S19 — Outer Space Treaty

**Publisher:** United Nations Office for Outer Space Affairs.  
**URL:** https://www.unoosa.org/oosa/en/ourwork/spacelaw/treaties/outerspacetreaty.html

Used for the existing international-law principles of responsibility, due regard, consultation and harmful contamination. The bundle does not claim the treaty alone resolves terraforming authority.

## S20 — Artemis Accords

**Publisher:** NASA.  
**URL:** https://www.nasa.gov/artemis-accords/

Used for current cooperative principles such as transparency, interoperability, emergency assistance, heritage protection and deconfliction. Terluna requires governance far beyond the Accords’ current scope.

## Additional evidence needed

The source set is intentionally compact and primary/official where possible. A research-grade next version should add:

- Multispecies atmospheric-escape literature.
- Lunar plasma and magnetotail models.
- Low-gravity combustion and fluid-mechanics data.
- Detailed lunar volatile and nutrient inventories.
- Cloud microphysics under partial gravity.
- Long-photoperiod ecology and controlled-environment agriculture.
- Artificial-gravity and multigenerational partial-gravity research.
- Lunar topography, crustal structure and basin-volume datasets.
- Comparative planetary-governance scholarship.

<!-- END FILE: 15_SOURCE_NOTES.md -->


---

<!-- BEGIN FILE: CHANGELOG.md -->

# Changelog

## Version 1.0 — 2026-07-12

- Reconstructed the Terluna concept from accessible Moon World context.
- Preserved the project name, five-hour practical twilight, aerial ecology, long light-dark biological adaptation, 35–45 km high atmospheric layer, and terminology corrections.
- Added an explicit evidence taxonomy separating canon, measurements, published inference, screening calculations and speculation.
- Introduced an 80 kPa quantitative reference atmosphere and alternative pressure range.
- Calculated atmosphere mass, gas inventory, scale height, lapse rate, water inventories, oxygen-mining scale and energy requirements.
- Added detailed climate, hydrology, soil, biosphere, human, infrastructure, governance, risk and research specifications.
- Added staged pressure plateaus and a paraterraforming fallback.
- Added machine-readable CSV, JSON and JSONL files, figures and reproducible Python calculations.
- Corrected an intermediate CSV error that had evaluated the text range `35–45 km` as subtraction.
- Corrected the illumination schematic so five-hour dawn and dusk intervals are counted inside the 29.53-day cycle.

<!-- END FILE: CHANGELOG.md -->
