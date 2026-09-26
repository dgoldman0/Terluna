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
