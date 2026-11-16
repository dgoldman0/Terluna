# The technology of the vacuum period, 2041–2241: a first roadmap

First draft of 16 November 2026, for the compact model of the vacuum economy; [technology.json](technology.json) holds its paths as data. Every economic study so far has held the founding design's engineering for 200 years, apart from a falling price of working parts and a slow growth of what a job does. This roadmap replaces that with paths: for each layer of the stack, a starting value, a slow, a central and a fast path to 2241, the step changes that could arrive and when, and what each changes in the engineering. Starting values marked [E§n] are the sourced values of the engineering lookup of 15 November (section n); [L§n] the local-content lookup; [T§n] the transport lookup. **Every path and arrival window here is the study's own projection, not yet sourced;** bounded lookups are to source the rates that move the results most.

The paths are written as values at 2041, 2070, 2100, 2150 and 2241; the model interpolates between them (geometrically where a quantity improves by a factor).

## 1. Food and engineered biology

Food is the settlement's largest load: in the founding design 14.5 kW of a resident's 25 kW average, almost all of it crop lighting (40 mol of light per square metre a day at 3.5 µmol/J, 1.25 kcal per mol of light, with losses of 1.3 and 1.15 and auxiliaries of 0.35). It is also where engineered biology moves fastest, and the project already treats its organisms as engineered.

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| LED efficacy, µmol of light per joule | 3.5 (design); 3.0 · 3.8 · 4.4 carried for the 2040s [F§3] | 3.8 · 4.0 · 4.2 | 4.1 · 4.4 · 4.9 | 4.5 · 4.9 · 5.2 |
| Crop food per mol of light, kcal | 1.25 (design), between wheat, 0.30 g of edible dry mass a mol, and potato, 0.75 g [F§2] | 1.3 · 1.4 · 1.5 | 1.45 · 1.7 · 2.4 | 1.7 · 2.2 · 3.2 |
| Share of the diet from electro-microbial food (hydrogen-oxidising and engineered microbes on CO₂ and electricity) | 0.03 (design: 18.7 g of protein a day) | 0.05 · 0.10 · 0.20 | 0.15 · 0.30 · 0.50 | 0.25 · 0.40 · 0.70 |
| Electricity to food energy, electro-microbial | 0.15 [F§1] | 0.15 · 0.16 · 0.18 | 0.18 · 0.20 · 0.24 | 0.20 · 0.23 · 0.26 |
| Food power per resident, kW (result of the rows above) | 14.5 | 12.6 · 10.6 · 8.5 | 9.5 · 6.3 · 3.0 | 6.6 · 3.9 · 1.5 |

[F§n] is the food lookup of 16 November (`vacuum-economy/research/technology_food/values.md`). It bounds these paths: the physical limit of a lamp is 5.27 µmol/J for 10% blue and 90% red light, and the best tested fixture gave 3.69 in 2022, with gains of about 10% a year over 2014–2022 projected to slow to 2–4%; crops convert at most 1.64–1.78 g of total biomass a mol, about 13% of the light's energy, or about 3.4 kcal of food a mol at a harvest index of one half, and engineered photosynthesis has gained 15–40% in field trials at ordinary CO₂; hydrogen-oxidising bacteria turn 14.6–16.5% of the electricity into whole dried biomass in two models, against about 1.3% for wheat and 4% for potato under lamps, with a theoretical limit of 25–27% for amino acids, and the only electricity-fed plant found makes 160 t a year (2024). The product can carry all of a diet's protein once its nucleic acids are cut to about 2%, which is 19–31% of the calories; shares above that need fats and sugars made from CO₂, for which the lookup found no energy balance, so the central and fast paths beyond about a third rest on that route.

Step changes: crops engineered for photosynthesis near the C4 limit; lighting tuned to the crop's action spectrum; food grown from electricity without light. What it changes: the settlement's power, the radiators, the farm floor (35 m² of a resident's 122 m²) and therefore the pressure hall built per resident.

## 2. Power: generation, storage and transmission

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Array specific power delivered to the surface, W/kg | 100 [E§2] | 150 · 200 · 400 | 300 · 1,000 · 3,000 | 1,000 · 3,000 · 10,000 |
| Cell efficiency | 0.32 (design) | 0.33 · 0.35 · 0.40 | 0.36 · 0.42 · 0.48 | 0.40 · 0.47 · 0.55 |
| Share of array mass made on the Moon (silicon, glass, aluminium from regolith) | 0 | 0.3 · 0.6 · 0.9 | 0.7 · 0.95 · 0.99 | 0.9 · 0.99 · 0.999 |
| Night storage, Wh/kg at system level | 550 [E§2] | 600 · 700 · 900 | 800 · 1,200 · 2,000 | 1,200 · 2,000 · 3,000 |
| Storage round trip | 0.40 [E§2] | 0.42 · 0.45 · 0.50 | 0.48 · 0.55 · 0.65 | 0.55 · 0.65 · 0.75 |
| Fission specific mass, kg/kWe | 150 [E§1] | 100 · 75 · 50 | 50 · 25 · 15 | 25 · 12 · 8 |
| Fission life, years | 10 [E§1] | 10 · 15 · 20 | 15 · 25 · 40 | 20 · 40 · 60 |
| Beamed power, electricity in orbit to the grid | 0.40 [E§3] | 0.42 · 0.45 · 0.50 | 0.45 · 0.55 · 0.65 | 0.55 · 0.65 · 0.75 |
| Beamed station, delivered W per kg in orbit | 100 [E§3] | 120 · 150 · 300 | 200 · 400 · 1,000 | 400 · 1,000 · 3,000 |

Step changes and windows: a lunar grid from the polar ridges, where the longest darkness is about 36 hours against 354 at the equator [E§2], by aluminium or superconducting cable (slow never, central 2080–2110, fast 2060–2080); a belt of arrays around the equator that needs no storage (fast only, after 2120); fusion, first on Earth and then shipped (slow after 2200, central 2090–2130 on the Moon, fast 2060–2080), with helium-3 from regolith as a later aneutronic fuel. What it changes: the mass and cost of firm power per resident, and whether beaming, local solar with storage, fission or fusion is the cheapest supply in each period; the model chooses among them by cost.

## 3. Launch from Earth

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Price to low Earth orbit, USD2026/kg | 600 [E§4]; the demand scenarios' 100 | 400 · 300 · 150 | 100 · 50 · 15 | 30 · 10 · 4 |

The floor of a fully reused chemical launcher is set by propellant and operations, a few dollars a kilogram; the fast path reaches it. Step change: launch that is not a rocket (fast only, after 2120). What it changes: every import, the passenger fare, the value of making things on the Moon, and the ceiling on the price of what the Moon sells in space.

## 4. Transport in space

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Tug and ferry exhaust speed, s of specific impulse | 450, H₂/O₂ | 450 · 450 · 900 | 900 · 900 · 2,000 | 900 · 2,000 · 5,000 |
| Lunar mass driver to orbit, kWh per kg launched | none | none · none · 3 | none · 2 (from about 2085) · 1.5 | 2 (from about 2060) · 1.5 · 1 |
| Active life of tugs and ferries, years | 10 [E§7] | 10 · 12 · 15 | 12 · 20 · 30 | 15 · 30 · 40 |

Nuclear thermal propulsion at about 900 s (central from about 2060), nuclear-electric and later fusion propulsion at thousands of seconds; electromagnetic launch from the lunar surface, which makes lunar material in orbit cheap and changes what the Moon can sell in space; aerobraking for cargo returning to Earth orbit [L§6]. What it changes: the propellant per tonne delivered, the reach to the sources (the velocity change a tug can afford), the size of the fleets, and the cost of moving lunar products.

## 5. Mining and processing

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Metal from regolith, kWh/kg | 35 [E§6] | 30 · 25 · 20 | 20 · 15 · 10 | 15 · 10 · 7 |
| Water extracted from asteroid rock, kWh/kg | 2 [E§5] | 2 · 1.8 · 1.5 | 1.5 · 1.2 · 1.0 | 1.2 · 1.0 · 0.8 |
| Water recovered from the rock heated | 0.7 (no source) | 0.7 · 0.75 · 0.8 | 0.8 · 0.85 · 0.9 | 0.85 · 0.9 · 0.95 |

What it changes: the power of the industry, the feed mined for each tonne of water and metal, and the plant at the sources.

## 6. Manufacturing, robotics and autonomy

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Imported share of new equipment by mass, once local production runs | 0.05 [L§1] | 0.05 · 0.03 · 0.01 | 0.02 · 0.005 · 0.001 | 0.005 · 0.001 · 0.0001 |
| Output per job, multiple of 2041 | 1 | 2 · 3 · 7 | 5 · 15 · 60 | 15 · 60 · 500 |
| Price of an imported kilogram of equipment, multiple of 2041 | 1 | 0.8 · 0.65 · 0.4 | 0.5 · 0.3 · 0.1 | 0.3 · 0.1 · 0.03 |
| Doubling time of the industry that builds industry, years | none | 20 | 8 | 3 |

Output per job carries robotics and artificial intelligence: the central path is about 5.5% a year to 2070, 3.8% to 2100 and 1.5% after, far above the 0.5–1.6% a year of the present scenarios. Step changes: chips made on the Moon (central 2090–2120), which removes the last large import; industry that copies its own machines from local material, with chips and a few elements from Earth. What it changes: the share of machinery bought from Earth, the residents' hands that the industry needs, and how fast the system can grow from its own surplus.

## 7. Habitat and life support

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Water recovered in life support and the farms | 0.98 [E§8] | 0.98 · 0.985 · 0.99 | 0.99 · 0.995 · 0.998 | 0.995 · 0.998 · 0.999 |
| Built mass per resident, multiple of the design's 335 t | 1 | 1 · 0.95 · 0.9 | 0.9 · 0.8 · 0.7 | 0.8 · 0.7 · 0.6 |
| Basic services power per resident, kW | 5 (design) | 5 · 4.5 · 4 | 4.5 · 3.5 · 2.5 | 4 · 3 · 2 |
| Life of ordinary equipment, years | 20 (design) | 20 · 22 · 25 | 25 · 30 · 40 | 30 · 40 · 60 |

Shielding and pressure set a floor on built mass: three metres of regolith overhead and a pressure envelope do not shrink with technology, so the reductions come from lighter halls, more floor per tonne and less farm floor. What it changes: the build per new resident, the upkeep, the water imported from the sources.

## 8. Computing

| Quantity | 2041 | Slow 2070 · 2100 · 2241 | Central 2070 · 2100 · 2241 | Fast 2070 · 2100 · 2241 |
|---|---|---|---|---|
| Computation per watt, multiple of 2041 | 1 | 4 · 10 · 50 | 20 · 100 · 1,000 | 100 · 1,000 · 30,000 |

Computing on Earth improves alike, so what the Moon can sell is set by its cost of power and cooling against Earth's, more than by the chips. What it changes: the computing line, and the autonomy that sets output per job.

## What the compact model does with this

Each quantity above is a path the model reads by year and by technology case (slow, central, fast), beside the demand scenario. The settlement's needs per resident follow food, basic services, built mass and water recovery; the industry follows processing energy, autonomy and local manufacture; the power supply is chosen in each period among the options by cost; transport follows launch price, exhaust speed and the mass driver. A run on the founding design's values held for 200 years remains as the static comparison.

## To source next

The paths that move the results most, by the founding design's loads: food per unit of light and the share of electro-microbial food; array specific power and its local manufacture; fission and fusion timelines; launch price after 2050; output per job under automation. Up to three bounded lookups, by layer, would source their rates and windows.
