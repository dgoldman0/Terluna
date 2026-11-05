# The Moon's energy supply: beamed power against its alternatives

5 November 2026 · `domain/provisioning`. Since 12 October the food and the sky towns take beamed power from the ring fleet, at 2.45 GHz to the towns' crowns and to closed receivers at the ground's loads ([power.md](power.md)). The question here is whether beaming through this deep, moist air makes sense: a drier, thinner world would suit it better, and producing the energy locally, local solar included, might still compete at a comparable scale. The air itself is a design condition: Earthlike, breathable, 1.0–1.2 atm, tall and dense. This study reviews what the project had analysed on the question, finds that analysis limited, and sets the supply options against one another on a common basis, the heat each leaves on the Moon included.

The screen, [energy_options.py](energy_options.py), writes [results/energy_options.json](results/energy_options.json) in about a second. It reads the domain synthesis's ledgers, the power screen, the heat-levers study, the first screen and the shared water scenario, and the design run's own sunlight at the ground, cloud and low winds, which [energy_extract.py](energy_extract.py) takes from the run's output into [results/energy_climate.json](results/energy_climate.json). Absorption runs through the power screen's ITU-R P.676-13 model. Numbers marked *product* come from a committed product, *computed* from this screen, *published* from a source listed in [energy_options_sources.json](energy_options_sources.json) with how far it was read, *recalled* from general knowledge not checked here, and *estimate* or *design guess* from a stated approximation. The product labels every outside value the same way.

## The answer

**The existing analysis was limited.** It weighed the beams against fusion on heat and against the first screen's storage and grid for the summit metropolis, and it sized the beams themselves in depth; it never set local solar under the design climate, wind, fission, raised receivers, other frequencies, optical links or loads made in orbit against the beams per delivered watt, on heat, land, storage, mass and failure together.

**Beamed power at 2.45 GHz should stand as the bulk supply for this world.** Every supply that is not the Moon's own sunlight, wind or water adds heat, and the beams add the least of those: 1.47–1.50 W for each watt delivered against 2.0–3.0 W for fission or fusion and 2.0–4.8 W for laser links. The dense air is a small part of that heat, 0.10–0.12 W, and it is the oxygen's: water vapour takes 3–5% of the absorption at 2.45 GHz, so a drier air would change almost nothing, and a 1.0 atm air would take about a sixth less. Local solar adds less heat, 0.19–0.79 W a delivered watt through its darker ground, but this air passes only 34% of the top-of-air sunlight to the ground, so each delivered terawatt takes 39,000–94,000 km² of panels and either a grid round the Moon or 300–330 TWh of storage for every lunar night. Each point of dimming that local solar saves costs 3.8–16 million km² of extra ground, 16–70% of the dry land, against the 0.55–1.38% of the land's rain and 0.07–0.30% of the growth that a point of dimming costs. Wind, rivers, tides, the Moon's heat flow and the seas' thermal gradients add no heat and reach about a terawatt together.

**The placement** (register, 5 November): computing runs in the array as far as latency allows, and the food the beams power is made on the Moon. Three-quarters of close computing can move, because the Moon's interactive services work within the 0.13 s round trip to the ring radius and only its millisecond control needs the ground ([computing](computing.md)). Moving it takes 0.06 points off the beams' heat, which falls to 28–152 TW on the Moon, 0.32–1.78 points, 3.0–22% of the resort's dimming (mix A2). The food would have been the larger lever: made at the ring radius it would take 0.12–0.97 points more, but its carbon, 1.5–6.2 times the managed import, and 2.3–3.4 Mt of phosphorus a year would arrive with it, and keeping both in their cycles would need freight back to the array or sinks on the Moon, while the food supply would depend on the array and its arrivals.

**What should change** (proposals): the ground's receivers should draw on the highest tiles, and better rectennas and tighter beams are worth as much heat as the air. Receivers raised to 20–40 km, another frequency or a thinner air each save at most 0.05 W a delivered watt; raising receivers earns its place only where it opens the closed ground, clears the reefs' layer or lets 5.8 GHz's smaller apertures work above the storms. Local supplies (rivers, roofs, the towns' own fuel cells) belong in the reserve for a fault across the ring, beside the week of stored hydrogen.

## What the project had analysed

| Source | What it covers on supply | Options set side by side | Measures | Evidence state |
|---|---|---|---|---|
| Register, 27 September | Fusion supplies most local power, supplemented by regional solar and wind | none | none | A planning assumption; fusion conditional ([engineering](../engineering/README.md)) |
| Register, 9 October | The planned climate kept, the added heat answered; industry and computing in orbit spare the Moon's heat; a terawatt used on the Moon adds 0.026 W/m² | heat placement | heat | Directions with first-screen arithmetic |
| Register and domain synthesis, 12 October | Array power for food and the sky towns, then beams on the recommended design; added dimming can be had by extracting more sunlight, the glow is the harder problem | array power against fusion | heat per delivered watt (1.38 against 2.0–3.0), land, exposure, reefs, a ring-wide fault | Chosen; screen arithmetic |
| First screen ([README](README.md), findings 1–4) | The night (71 TWh for the metropolis, 7.9 times Earth's pumped storage, 52 km³ through Korolev); an HVDC grid round the Moon losing 9–19%; 6,100 km² of level panels for the metropolis's mean; microwaves from the night half at 10 W/m²; rivers 0.87 TW gross; tides 25–50 GW; deuterium | storage, grid, beams for one city | energy, area, losses | Screen on an earlier climate's sunlight |
| [power.md](power.md) | Absorption (P.676-13), rain (P.838-3), tiles in view, transmitters, receivers' land, high latitudes, reserves | 2.45 against 5.8 GHz; lines against low beams near the poles | loss, heat per watt, land, tiles, faults | Computed on Earth-validated models and design guesses |
| [Heat levers](../research/studies/heat_levers/README.md) | Beams 0.39–1.84 points; food made in orbit saves 0.13–0.98 with its carbon 1.5–6.2 times the import; the top-up's fall; local sunlight, wind and rivers add no heat but cost land and storage | food in orbit against food on the Moon | points, carbon | Arithmetic on products; local supplies named, not sized |
| [Climate targets](../research/studies/climate_targets/README.md) | The heat beyond the sunlight, the beams 0.87–4.14 W/m² of it; the resort's 8.2–10.9 points | none | W/m², points, rain and growth a point | Arithmetic on the GCM |
| [Array industry](../research/studies/solar_shield_array/array_industry.md), [glow](../research/studies/solar_shield_array/glow.md) | Where heat is released; an absorbing dimmer returns 10–25% of its dimming as glow | orbit against Moon for heat | W/m² per TW | Screens |
| [Computing](computing.md) | Millisecond control on the Moon, the rest in orbit | placement | heat, latency | Screen |
| [Food systems](../research/studies/food_systems/README.md) | The blend's power-made and lamp-lit food, 42–56% of the beams' heat; the lamp-lit diet with array power or fusion | array against fusion for food | heat | Screen |
| [Engineering](../engineering/README.md) | Fission bootstrap, then fusion, then near-Sun collection with beaming (19 September reconstruction) | sequence only | none | Reconstruction |

The analysis is limited where it matters for the choice: the main alternatives were never set against the beams on a common basis. Local solar was never sized on the sunlight the design run now gives at the ground, under its cloud and deep air, and the night's storage was sized only for the metropolis. Wind was never read from the design winds. Fission, raised receivers, lower frequencies and optical links were not compared. The coupling between the beams' heat and the dimming that powers the tiles was implicit, and the absorbing dimmer's glow return was not carried onto the beams' points. Food made in orbit was counted for its heat and carbon, not its phosphorus or its descent.

## Findings

1. **The dense air takes 0.10–0.12 W of the 1.48–1.50 W of heat each delivered watt of ground beams leaves, and nearly all of it is the oxygen's** (*computed*).
   - **Where the heat goes.** Of each watt delivered to the ground at 2.45 GHz, the loads use 1.00 W, the rectennas lose 0.19 W, the spill round the receiver 0.14 W, the bus 0.05 W and the air 0.10–0.12 W, from the highest tiles to any usable tile. With no air at all the same link leaves 1.38 W. The towns' beams at 10 km leave 1.45–1.47 W, 0.07–0.09 W of it in the air.
   - **Why drier would change almost nothing.** Water vapour takes 3.2% of the overhead loss at 2.45 GHz through the 295 K column and 5.1% through the 301 K column; at 5.8 GHz it takes 15–23%. The absorber is oxygen's pressure-broadened relaxation spectrum, which is nearly flat between 1 and 10 GHz at this pressure, so 1 GHz loses 4.9% overhead against 2.45 GHz's 5.9%.
   - **Why thinner would help a little.** A 1.0 atm air at 301 K, set beside the design air for comparison only, loses 4.9% overhead against 5.8% for 1.2 atm at the same temperature. Carried to the ground's beams as used, the thinner air, a lower frequency or the warmer design climate each move the heat by about 0.02 W a delivered watt.
   - **Height helps more, and the line down takes most of it back.** From receivers at 20, 30 and 40 km the ground's beams lose 3.7–4.5%, 2.7–3.3% and 1.9–2.4%. A tether's conductor and converter take 2–3% (*estimate*), so the heat falls only to 1.43–1.49 W. The summit tower's crown at 35.6 km could carry such receivers for the metropolis; a sky town used as a relay gains nothing, since its beam down crosses the same lower air.
   - **What matters as much as the air.** Raising the rectennas' efficiency from 85% to 90% or the share each receiver captures from 90% to 95% would each cut the heat by about 0.08 W a delivered watt, two-thirds to four-fifths of the air's share; the 1977 record rectenna element reached 90.6% (*published*).

| Overhead loss | 1 GHz | 2.45 GHz | 5.8 GHz | 10 GHz |
|---|---:|---:|---:|---:|
| 1.2 atm, 295 K (the power screen's column) | 4.9% | 5.9% | 7.2% | 10.1% |
| 1.2 atm, 301 K (the design climate) | 4.8% | 5.8% | 7.5% | 11.8% |
| 1.0 atm, 301 K (comparison only) | 4.1% | 4.9% | 6.4% | 10.2% |
| Water vapour's share of it, 1.2 atm, 301 K | 1.0% | 5.1% | 23% | 47% |

| From a receiver at | 0 km | 4 km | 10 km | 20 km | 30 km | 40 km |
|---|---:|---:|---:|---:|---:|---:|
| 2.45 GHz overhead, 1.2 atm, 295 K | 5.9% | 5.2% | 4.4% | 3.2% | 2.3% | 1.7% |
| 5.8 GHz overhead | 7.2% | 6.2% | 5.0% | 3.5% | 2.5% | 1.8% |

2. **Every supply but the Moon's own sunlight, wind and water adds heat; the dimming that answers the beams' heat gives back about a fifth to a quarter of their power, and the resort's dimming can give all of it** (*computed* on the heat-levers and power products).
   - **A point.** A point of dimming is 1% of the hot world's absorbed sunlight, 2.25 W/m² or 85.5 TW over the Moon, as the climate targets and the heat-levers study count it. It intercepts 111 TW of light at the shield, which collectors at 20–30% turn into 22–33 TW, and it costs 0.55–1.38% of the land's rain and 0.07–0.30% of the growth through light (*product*).
   - **The beams per delivered TW.** The ground's beams need 0.0173–0.0176 points, 0.019–0.023 if the dimmer absorbs on the window rings and returns 10–25% of its dimming as glow (*product*, glow). Those points intercept 1.9 TW of light, which make 0.38–0.58 TW of electricity: 18–27% of the 2.1 TW the tiles must make for each delivered TW.
   - **The resort's dimming as a whole.** At 8.2–10.9 points it intercepts 910–1,210 TW of light, 182–242 TW of electricity at 20% conversion alone, more than the 47–224 TW the tiles need at both ends. If the dimmer is built as collectors, the beams need no light beyond what the climate already takes, and the night half's 851 PW stays a reserve. Their cost is the points: per delivered TW, 0.0095–0.024% of the land's rain and 0.0012–0.0052% of the growth.
3. **The deep, cloudy air passes a third of the sunlight to the ground, so local solar is land-hungry here, and the 354-hour night needs a grid round the Moon or a store refilled every lunar day** (*computed* on the design run).
   - **Sunlight at the ground.** Over years 10–19 of the design run the top of the air receives 300 W/m² on the cycle's mean and the ground 101.7 W/m², 34%; under a clear sky 142.9 W/m², 48%. The ground's loads, at C′'s latitudes, see 92.7 W/m²; cloud covers the land 36% of the time. The column's Rayleigh optical depth at 550 nm is 0.75 against Earth's 0.10 (*computed* on a *recalled* Earth value), so even an overhead clear-sky sun keeps under half its direct beam at 550 nm, which weighs on concentrating solar thermal more than on panels.
   - **Panels.** With today's 22% modules or 30% tandems (*published*, *estimate*), a 0.85 derate and 0–15% from tilting, a panel delivers 15.6–25.7 W/m² through a grid and 10.6–16.6 W/m² through hydrogen. A delivered TW takes 39,000–64,000 km² of panels with a grid, 60,000–94,000 km² with hydrogen, 12–21 TW of peak capacity, at a capacity factor of 5–9%, and 1.2–2.7 Gt of panels and mounts at Earth's material intensity. Panels repay the energy that made them in about a year in Europe (*published*); under this ground's dimmer sunlight they would take about one and a half to two years before any store (*estimate*, Europe's figure scaled by the sunlight).
   - **The night.** Without a grid a TW of load needs a store of about 300–330 TWh each lunar day: 18 Mt of hydrogen burned at half its heating value, 1.2 km³ of it at 20 MPa; or 230 km³ of water through Korolev's 3.74 km head; or 1.1–2.2 Gt of batteries. Heat stores, which carry Earth's solar thermal plants through hours (*recalled*), would have to carry 354. With a grid round the Moon, lines carry 48% of the energy from the day side over 4,000–5,500 km and lose 9–19% of it: 0.17–0.32 million km of 12 GW lines a delivered TW, 5–10 Mt of aluminium and 8,000–32,000 km² of corridors (*estimates* on the first screen's distances and a 1,100 kV line).
   - **Their heat.** Panels of albedo 0.05–0.10 on land the GCM gives 0.2, with 53–60% of the light reflected from the ground leaving the Moon (*estimate* from the heat-levers kernel), add 0.19–0.54 W a delivered watt through a grid and 0.30–0.79 W through hydrogen. On this world local solar adds heat too, an eighth to a half of the beams'. Floating panels on the seas, darker or brighter than open water by a few hundredths, add almost none, but they shade the seas' life and cut the evaporation that feeds the rain. Panels raised above the cloud cover a planet of albedo 0.40 and add 1.1–1.9 W a delivered watt, as much as the beams.
4. **Each point of dimming that local solar saves costs 3.8–16 million km² of extra ground** (*computed*).
   - **The exchange.** Against the ground's beams, local solar with a grid saves 0.011–0.015 points a delivered TW and takes 56,000–185,000 km² more ground: 3.8–16.4 million km² a point, 16–70% of the dry land. With hydrogen it takes 5.7–27.5 million km² a point. A point of dimming costs 0.55–1.38% of the land's rain and 0.07–0.30% of the growth.
   - **The sunlight life loses.** Per delivered watt the panels take 3.6–5.9 W of the ground's sunlight through a grid and 5.6–8.7 W through hydrogen, on their own land; the beams' dimming takes 0.65–0.66 W, thinly over the whole Moon.
   - **At C′'s loads.** The ground's 16–73 TW on local solar with a grid would take 1.0–14.0 million km², 9–120% of the 11.7 million km² outside the kept 30% and the parks, against the beams' 0.12–0.55 million km² of closed receivers. The low-cost end fits beside the blend's fields of 0.68–1.86 million km²; the high-cost end does not fit.
5. **Wind, rivers, tides, the heat flow and the seas' thermal gradients add no heat, and together they reach about a terawatt** (*computed* on the design run and *published* values).
   - **Wind.** The design run's lowest level, about 0.3 km up, blows 1.7 m/s over land in 3-day means and stays below a turbine's 3 m/s cut-in 89% of the time. Its power density, 7.8 W/m² of swept area, is 2–7% of the Global Wind Atlas's 323–410 W/m² even with the means' missing gusts counted as up to three times the cube (*estimate*). Scaled from the 0.5 W/m² Earth's wind farms realise (*published*), a TW would need 28–105 million km², more than the dry land. The surface stress dissipates 0.15 W/m² over land, 3.5 TW over all of it, the ceiling on everything the low winds carry. The cloud-resolving rings' instantaneous winds, every 3 hours over the earlier design case, agree: over land at 10 m the median is 0.9–1.2 m/s and the 99th percentile 3.8–5.2 m/s at the equator and 45°, and the polar rings are nearly calm. Fitted, they hold 6–12 W/m² at 10 m and 27–47 W/m² at 1 km, at most 14–27 and 57–121 W/m² with every slice of time at its slice's top speed (*computed* on the rings' products). On these winds the regional wind in the 27 September assumption has little to draw on.
   - **Rivers** carry 0.87 TW gross (the first screen, on an earlier climate), of which plants could take 0.22–0.44 TW, dispatchable through the lunar night.
   - **Tides** give 0.03–0.06 TW from the impounded seas (the first screen).
   - **The heat flow,** 14–21 mW/m² at the Apollo sites (*published*), is 0.53–0.80 TW over the Moon; tapping every square metre would give 0.05–0.12 TW of electricity, and mining the crust's stored heat faster than it flows would add 6.7–10 W of heat a delivered watt.
   - **The seas' thermal gradients** could give 0.24–0.50 TW if Earth's 7–14 TW of ocean thermal power scales by area (*published* on Earth), and only if cold deep water forms in the fresh seas, which the seas' circulation work has not yet shown.
6. **Fission and fusion leave 2.0–3.0 W a delivered watt, a third more to twice the beams'** (*recalled* efficiencies, *project* fusion figure).
   - **Fission** at 33–45% leaves 2.2–3.0 W and fissions 0.85–1.2 kt of heavy metal a delivered TW-year. The Procellarum KREEP Terrane holds 4–12 ppm thorium (*published*), so a TW-year needs 0.07–0.29 Gt of rock even in breeders, and the 28% sea floods most of that terrane (the first screen's Oceanus Procellarum).
   - **Fusion** leaves 2.0–3.0 W; deuterium–tritium burns 75–112 t of deuterium and 220–340 t of lithium-6 a TW-year, against the seas' 4.8×10¹⁴ kg of deuterium. The best laser shots reach a target gain of 4.13 and ITER plans deuterium–tritium operation from 2039 (*published*); no plant has made net electricity. It stays the conditional technology the engineering lane names.
   - **Where they win.** Neither wins on heat. They supply power independent of the ring: before the fleet transmits, and as a reserve.
7. **Optical links lose more of the beam, leave more heat, and stop in the cloud** (*estimate*). The deep column's Rayleigh optical depth at 850 nm is 0.13; with vapour and aerosol a beam at the ground's operating elevation keeps 61–81% of itself under a clear sky. With 50–70% lasers and 40–65% cells (*published*: 68.9% for the best cell, about 20% for the 2025 record link), a delivered watt leaves 2.0–4.8 W of heat and takes 3.2–9.6 W at the tiles. Cloud covers the land 36% of the time. Laser light scattered in the deep air also reaches the night sky the project keeps dark.
8. **Food made in orbit with the beams' power would be the largest lever on their heat** (*computed* on the heat-levers product).
   - **Heat.** The blend's food takes 7.4–56.3 TW delivered. Made at the ring radius, it removes 11.1–84.0 TW of heat from the Moon and brings back 0.63–0.94 TW as the food's own energy and 0.09–0.49 TW as its fall, 33,000–180,000 kg/s at 2.7 MJ/kg, shipped dry to fresh: 0.12–0.97 points saved, 0.016–0.017 points for each TW of food power moved. Moving three-quarters of close computing to orbit, leaving the millisecond control, saves 0.06 points on its own.
   - **What it costs.** Its carbon arrives at 0.47–0.71 Gt a year, 1.5–6.2 times the managed import that holds 400 ppm (*product*); its phosphorus at 2.3–3.4 Mt a year, two-thirds to all of what the whole population's sewage carries, as a net import that the towns' struvite and burial would have to keep out of the nutrient-poor waters; food security then rests on the array and on arrivals; the food plants join the array's industry.
   - **Status.** Not taken up: the food stays on the Moon and the blend as set, and three-quarters of close computing moves to the ring radius (register, 5 November; mix A2).
9. **Seven mixes against the loads** (*computed*; pairs span the low-cost end with each option's better values to the high-cost end with its worse).

| Mix | Heat on the Moon | Points | Share of the resort's dimming | Land's rain | Growth | Ground taken | Share of the land outside the kept 30% and the parks | Store |
|---|---|---|---|---|---|---|---|---|
| A. Beams as chosen, 2.45 GHz to towns and ground | 33–157 TW | 0.39–1.84 | 3.5–22% | 0.21–2.5% | 0.03–0.54% | 0.12–0.55 million km², closed | 1.0–4.7% | A week of critical loads, 425–1,270 TWh |
| A2. Beams as chosen, with the food on the Moon and three-quarters of close computing at the ring radius (the adopted placement) | 28–152 TW | 0.32–1.78 | 3.0–22% | 0.18–2.5% | 0.02–0.53% | 0.09–0.52 million km², closed | 0.8–4.4% | As A |
| B. Beams, with the food they power and three-quarters of close computing made at the ring radius | 17–69 TW | 0.20–0.81 | 1.9–9.9% | 0.11–1.1% | 0.01–0.24% | 0.04–0.10 million km² | 0.3–0.8% | As A |
| C. Beams to the towns; the ground on local solar with a grid round the Moon | 12–87 TW | 0.14–1.02 | 1.3–12% | 0.08–1.4% | 0.01–0.30% | 1.0–14.0 million km², with corridors | 9–120% | None; 2.7–23 million km of lines |
| D. Beams to the towns; the ground on local solar with hydrogen | 14–106 TW | 0.16–1.24 | 1.5–15% | 0.09–1.7% | 0.01–0.36% | 1.4–17.2 million km² | 12–147% | 4,800–21,900 TWh each lunar day, 0.29–1.31 Gt of hydrogen |
| E. Beams to the towns; the ground on fission or fusion | 41–269 TW | 0.48–3.15 | 4.4–38% | 0.26–4.3% | 0.03–0.93% | Small | Small | Fuel |
| F. Fusion for all, the 27 September assumption | 44–317 TW | 0.52–3.70 | 4.8–45% | 0.28–5.1% | 0.04–1.1% | Small | Small | Fuel |

10. **When each is best** (*interpretation* of the findings).
    - **Beams (A, A2, B)** are best while a point of dimming stays cheaper than the 3.8–16 million km² local solar would cover to save it, which holds while the shield can take more of the sunlight, as the direction of 12 October expects. Their weaknesses are the ring-wide fault, the closed ground under the receivers and the fleet's added points (a reflecting dimmer adds 0.022 to the window tiles' light-pressure coefficient a point, an absorbing one returns 10–25% as glow; *product*).
    - **A2 over A** throughout: the computing's results return as data, and only its millisecond control needs the ground.
    - **B over A2** only where the carbon and phosphorus can return: burial or kilns on the Moon for the carbon, struvite for the phosphorus, or freight back to orbit. A2, the adopted placement, keeps both in the Moon's own cycles.
    - **Local solar** is best for small, dispersed loads and for the reserve. The roofs and pavements of the settled land, 0.26–0.74 million km² (*product*), carrying panels on 30–50% of their area and storing in hydrogen, would deliver 0.8–6.2 TW, the size of the critical loads (2.5–7.5 TW), and would compete with the white roofs the heat-levers study ranks fourth. Carrying the ground's bulk loads (C, D) is best only if the shield's dimming becomes the binding limit, and then only at the low-cost end, where it takes land equal to the blend's fields.
    - **Fission or fusion (E, F)** are best only where the ring cannot serve: in the build before the fleet transmits, and as independent backup.
    - **Optical links, wind and panels above the cloud** are best nowhere on this world's ground.
11. **Fit with the world's design** (*interpretation*).
    - **Sky reefs and the towns.** The beams cross the reefs' layer at 10–20 km and step aside 0.2–14 times a day a receiver (power screen); receivers raised above 20 km would take the ground's beams out of that layer altogether, at the price of tethers in the flight bands. Local solar leaves the sky alone, but panels on the towns or above the cloud darken a bright planet.
    - **The shield.** Every option keeps the ring fleet's 4-lunar-radius shadow as built; the beams' 0.39–1.84 points are dimming taken within it.
    - **The cities.** The metropolis's receivers make a band 7–23 km wide round the city (power screen); local solar for the same load would cover 8–26 times the ground.
    - **The water.** Pumped storage for local solar would cycle 230 km³ of water through a 3.74 km head for each delivered TW every lunar day, thousands of km³ for the ground, against highland lakes the design keeps; the controlled water delivery's fall, 58–359 TW, stays unrecovered (heat levers).
    - **The dark night.** Microwaves and panels make no light; laser links scatter it.
    - **Conservation and the wild share.** The beams' receivers and local solar at the low-cost end fit in the land outside the kept 30% and the parks; local solar at the high-cost end would need land the kept share and the parks hold.
    - **Materials.** Panels and rectennas draw on what the far side's highland gives, aluminium, glass and silicates, with steel from distant dry mare or space (README, finding 8); the tiles come from the array's own industry on space resources. Hydrogen for stores comes from the seas' water or the Solar-System-wide operation's volatiles. Batteries at 1.1–2.2 Gt a TW would need lithium the Moon lacks, and iron-air would draw on iron the sea mostly covers. Fission's thorium lies mostly under the near-side sea; fusion's deuterium is in the seas and its lithium comes from the source industry (README, finding 4).
    - **Outside the supply.** The towns' hydrogen remade after leaks takes 1.7–17.8 TW delivered. At the high-cost end halving the leak would save about 13 TW of heat, 0.15 points, four times what raising every ground receiver to 40 km would save.

## The options on a common basis

Per delivered TW. Heat is what the option leaves on the Moon beyond the sunlight it would have absorbed anyway; a point is 85.5 TW of heat.

| Option | Heat, W a delivered W | Points a TW | Ground a TW | Night store a TW | Mass and materials a TW | How it fails | Maturity | Evidence |
|---|---|---|---|---|---|---|---|---|
| 2.45 GHz beams to the ground (chosen) | 1.48–1.50 | 0.0173–0.0176 | 7,500 km² closed, planting beneath | None; a week of critical loads, 12–19 TWh | Rectennas 0.03–0.12 Gt (*design guess*); about 10 Mt in orbit at the 1978 reference's 10 kg/kW | Ring-wide fault; storms cost 1–23% of a beam | Rectennas 82–91%; ground links of tens of kW; milliwatts in space; no system | *product*, *computed* |
| 2.45 GHz to the towns' crowns | 1.45–1.47 | 0.0170–0.0172 | On the towns | As above | On the towns | As above | As above | *product* |
| 5.8 GHz to the ground | 1.50–1.53 | 0.0175–0.0179 | As 2.45 GHz | As 2.45 GHz | Apertures 2.4 times smaller | Storms stop it | As above | *product*, *computed* |
| 1 GHz to the ground | 1.46–1.48 | 0.0171–0.0173 | As 2.45 GHz | As 2.45 GHz | Apertures 2.45 times wider | As 2.45 GHz | As above | *computed* at P.676's lower limit |
| Receivers raised to 30 km, lines down | 1.44–1.47 | 0.0169–0.0172 | Ground opened | As 2.45 GHz | Aerostats and tethers | A tether; lines in the flight bands | Tethered aerostats on Earth (*recalled*) | *computed* on *estimates* |
| The same at 5.8 GHz | 1.45–1.47 | 0.0169–0.0172 | Ground opened | As 2.45 GHz | Smaller apertures | Above the storms | As above | *computed* on *estimates* |
| Laser link | 2.0–4.8 | 0.023–0.056 | Small | Cloud outages | Lasers 3.2–9.6 W at the tiles | Cloud 36% of the time | 800 W over 8.6 km (2025) | *estimate* |
| Food made at the ring radius | 0.03–0.10 a W of food power moved | 0.0003–0.0011 | Opens the receivers' ground | Food stores | Food plants in the array | Supply through arrivals | Nothing at orbital scale | *computed* |
| Local solar with a grid round the Moon | 0.19–0.54 | 0.0022–0.0063 | 64,000–193,000 km² | None; 0.17–0.32 million km of lines | Panels 1.2–1.8 Gt, aluminium 5–10 Mt | Distributed; a cut line darkens the night side | PV and 1,100 kV HVDC commercial | *computed* on the design run |
| Local solar with hydrogen | 0.30–0.79 | 0.0035–0.0093 | 86,000–236,000 km² | 300 TWh, 18 Mt of hydrogen | Panels 1.8–2.7 Gt | The store carries a 354-hour night | Components commercial | *computed* |
| Local solar with pumped water | 0.20–0.56 | 0.0024–0.0065 | 59,000–166,000 km² | 326–329 TWh, 230 km³ cycled | Panels 1.2–1.9 Gt | Reservoirs swing every lunar day | Pumped hydro commercial | *computed* on *estimates* |
| Local solar with batteries | 0.19–0.52 | 0.0023–0.0060 | 56,000–153,000 km² | 332 TWh, 1.1–2.2 Gt | Panels 1.2–1.8 Gt | Distributed | Commercial at hours | *computed* on *estimates* |
| Panels above the cloud | 1.1–1.9 | 0.013–0.022 | Shades the ground | As local solar | Aerial platforms | As platforms | None | *estimate* |
| Floating panels on the seas | −0.02 to +0.01 | About 0 | Sea surface | As local solar | As local solar | As local solar | Commercial on reservoirs | *estimate* |
| Wind | 0 | 0 | 28–105 million km² | The winds are low by day and night | Turbines | Low winds | Commercial on Earth | *computed* on the design run |
| Rivers | 0 | 0 | Reservoirs | Dispatchable | Dams | Drought | Commercial | *product*; 0.22–0.44 TW in all |
| Seas' thermal gradients | 0 | 0 | Sea | Steady | Plants | Deep water's temperature | Pilots (*recalled*) | *estimate*; 0.24–0.50 TW in all |
| Tides | 0 | 0 | Barrages | Monthly | Barrages | A small monthly range | Commercial | *product*; 0.03–0.06 TW in all |
| Heat flow | 0 at the natural flow; 6.7–10 if mined | 0 to 0.12 | Wells | Steady | Wells | A weak flow, slow to recharge | Earth's hot spots only | *published*; 0.05–0.12 TW in all |
| Fission | 2.2–3.0 | 0.026–0.035 | Small | Fuel | 0.85–1.2 kt fissioned; 0.07–0.29 Gt of rock | Plant by plant | Commercial (*recalled*) | *recalled* efficiencies |
| Fusion | 2.0–3.0 | 0.023–0.035 | Small | Fuel | 75–112 t of deuterium | Plant by plant | Target gain 4.13; no plant | *product*, conditional |

## Method

- **Loads.** The synthesis's array case by component, over the Earth-service 1.376 W a delivered watt it was built on: residents 8.0–23.9 TW delivered, food 7.4–56.3, the towns' hydrogen remade 1.7–17.8, their control 0.28–2.8, trips 0.08, close computing 4.6; 22.1–105.5 TW in all, the power screen's total. The towns take their control, hydrogen, trips and the aloft half of the residents, 6.1–32.6 TW; the ground the rest, 16.0–72.9 TW.
- **Beams.** The power screen's chain, 70% at the transmitter, 90% captured, 85% at the rectenna and 95% at the bus, with the air's share as the power screen gives it for the ground's beams as used, from the highest tiles to any tile above 30°. Variants carry its load-weighted loss by the ratio of losses at the ground beams' mean operating elevation of 61°, computed through each column.
- **Columns.** The power screen's moon column rebuilt at 295 and 301 K, and a 1.0 atm column at 301 K with Earth's oxygen partial pressure, for comparison.
- **Local solar.** Sunlight at the ground, all-sky and clear, from the design run's surface shortwave over its surface albedo, by latitude band over land, weighted by where the ground's loads lie in C′. The lunar day's shape is binned by hour angle from the run's own Sun, smoothed over a tenth of the day. A constant load's store is solved for each band: the generation that balances the cycle at the store's round trip, and the longest run of deficit.
- **Heat of panels.** The albedo change times the share of light reflected from the ground that leaves the Moon, the heat-levers study's clear-sky kernel of 0.34 over the run's clear-sky transmission of 0.50, less a cloud's return.
- **Food in orbit.** The heat-levers product's food heat and carbon, the food's own energy at 2,500 kcal a person a day, the fall from the ring radius with the ring's orbital speed, and phosphorus at the 0.44 kg a person a year the metropolis's sewage carries (README, finding 11).

## Assumptions

The product lists each with its kind. The ones the conclusions lean on:
- the power screen's link chain and its 85% rectennas, at 5.8 GHz too;
- collectors turning dimmed light into electricity at 20–30%, as the climate targets take them;
- panels of 22–30% efficiency, a derate of 0.85, ground cover of 0.4–0.7 and albedo 0.05–0.10 against the GCM's land at 0.2 (*estimates*);
- 53–60% of light reflected from the ground leaving the Moon (*estimate*);
- hydrogen's 40% round trip, pumped water's 75–80%, batteries' 88%, iron-air's 40–45%;
- HVDC losses of 9–19% over 4,000–5,500 km, 12 GW lines and peak transfers of half to seven-tenths of the ground's load (*estimates*);
- a tether's line and converter at 2–3% (*estimate*);
- laser links at 50–70% wall-plug, 40–65% cells and 90% capture, with 0.05–0.30 of optical depth beyond Rayleigh's (*estimates* on *published* points);
- fission at 33–45% (*recalled*), fusion at 2.0–3.0 W a watt (*product*).

## Limits

- **The air.** P.676's line parameters are Earth's; the columns are clear and the clouds' liquid water is left out. The 1.0 atm column keeps the design air's other choices and is a comparison only; the air stays a design condition.
- **The ground's sunlight.** It is a T21 GCM's, about 170 km a cell, with a fixed land albedo of 0.2 and 3-day means. The day's shape is smoothed over a tenth of the day; the store per TW is a band's mean, and a site's darkest lunar day needs more; the summit's height and tilted or tracking panels' gain beyond 15% are left out. The grid case takes the loads' latitudes for the panels round every longitude; the near side is mostly sea, so its share of the panels would float or crowd its highlands.
- **The panels' heat.** It rests on a clear-sky kernel and an estimated cloud return. Panels' effects on evaporation, cloud and the rain are left out, and they could be regional and of either sign.
- **Wind.** The design run's winds are 3-day means of a T21 model, and the cloud-resolving rings that confirm them are two-dimensional, on flat ground, over the earlier design case; slope and valley winds over the real terrain, and airborne turbines high in the superrotating flow, are left out.
- **The seas' thermal power.** It is Earth's per area, and the fresh seas' deep temperature is open.
- **Food in orbit.** Its descent is counted only as heat; its route through the extended air and any limit on what descending freight may deposit high in it, in the spirit of S7, are left for the descent's own study, beside the water top-up's far larger flow.
- **Mass and maturity.** Masses are orders of magnitude from Earth's material intensities and the 1978 reference; maturity records today's evidence.
- **Literature.** Most sources were read as search summaries or abstracts, as recorded in the sources file.

## Next work

1. **The ground's sunlight at its worst.** A site's store through its darkest lunar day, and panels at height on the summit and the dry highlands, from the design run's 3-day records rather than band means.
2. **Wind over the terrain.** Slope, valley and coastal winds over the real relief, and the superrotating winds aloft for tethered airborne turbines, which the flat rings leave out.
3. **The reserve for a fault across the ring.** Rivers, roofs and the towns' fuel cells sized against the critical loads, in place of part of the week's hydrogen.
4. **The beams' larger losses.** Rectenna efficiency and the share each receiver captures, which outweigh the air.

## How to run it

```
climate/gcm/.venv/bin/python -m provisioning.energy_extract   # reads the design run (2.3 GB), writes results/energy_climate.json
python -m provisioning.energy_options                          # writes results/energy_options.json in about a second
python -m pytest provisioning/tests/test_energy_options.py
```

The products (schemas `terluna.provisioning.energy-climate/1` and `terluna.provisioning.energy-options/1`) bind their producers, inputs and sources by sha256, with their evidence statements, reading rules and assumptions. The tests regenerate the options product, check the inputs' hashes, the sums of each beam's heat, the coupling between the beams' heat and the dimming's electricity, the absorption's fall with height and its flatness in frequency, the store solver's energy balance and the mixes' sums; the extract's own regeneration runs when `TERLUNA_RERUN_EXTRACT` is set in the GCM's environment.
