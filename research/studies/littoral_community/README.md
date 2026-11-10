# The freshwater littoral community

9 November 2026, on domain/ecology. The study asks what a self-sustaining community of the Open Moon's fresh shallows needs: producers, decomposers, grazers and predators that feed one another and reproduce on their own through the 708.7-hour lunar day, on the shores of the held fresh seas and at the margins of the rain-fed lakes. It belongs to the aquatic programme's shelves, shallow reefs and lake margins ([plan](../aquatic_biogeochemistry/plan.md#habitats)).

The folder holds four pieces. The first is the shore water's physical setting through the design climate's lunar days, computed at five sites. The second carries that setting to every cell of the climate model as a map of the waters' warmth. The third is a review of founder species ([founders.md](founders.md)): which Earth species could found the community, what makes each easy to select, where each belongs and in what order they go in. The fourth is a community simulator, built and tested and then shelved: its mechanics work and its species traits are placeholders, so it produces no result about the Moon.

## The shore's water through the lunar month

[forcing.py](forcing.py) carries the corrected design climate (`A28R_301k_design_l15_snow_s938_turned_refit_precession`, years 10–19, the design years of [8 November](../../decisions.md)) to the water of five shore sites: four sea cells beside land, one in each of four latitude bands, each the band's coastal cell with the median ten-year surface temperature, and one equatorial land cell with the median land temperature for the lake margins. It writes [results/forcing.json](results/forcing.json) (schema `terluna.research.littoral-forcing/1`, with the hash of every code file and input, the shared constants and the design air's values it read) and hourly series on the research drive (`/media/projectspace/terluna-research/littoral-community/forcing/`), in about four minutes:

```sh
nice -n 10 env OPENBLAS_NUM_THREADS=1 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.forcing
```

Each site's water is computed two ways. A closed pan of 1, 2 or 4 m runs in the project's lake model ([lake.py](../../../biosphere/ecology/lake.py)), forced as the [cold-waters study](../cold_waters/README.md) forces its ponds, with the sunlight shaped within each 3-day window by the Sun's height so that the water warms from dawn and cools from dusk. A slab of the same depth takes the same surface fluxes and trades a tenth, three tenths or all of its water each day with the open water beside it: the run's sea surface at the sea shores, and a 20 m lake's surface at the lake margin. The closed pan and the fully exchanged slab bracket a real shore, whose exchange by waves, wind-driven currents and the monthly tide is still to be modelled.

### What it finds

**The shore's light and air.** The sea shores' open water, as the climate model carries it in a mixed layer 50 m deep, averages 28.9–30.4 °C and moves about a kelvin either way through the lunar day and night. Under the run's clouds, which pass 63–81% of the clear sky's sunlight, photosynthetic light just below the surface averages 27 mol m⁻² per Earth day at 8° and 9 at 64°, and it stays below 1 µmol m⁻² s⁻¹ for 261 hours of each lunar day at 8°, 244 at 30°, 210 at 47° and 109 at 64°. Oxygen crosses the surface at 1.2–1.8 m a day on the sea shores and 0.7 m a day at the sheltered lake margin, so a 2 m column exchanges its oxygen with the air over one to three days.

| Site | Open water, mean (range) | Dark hours a lunar day | Light below the surface | Cloud factor | Oxygen transfer velocity, mean (lowest tenth) |
|---|---|---|---|---|---|
| Sea shore, 8° N | 30.1 °C (29.1–31.0) | 261 h | 26.7 mol m⁻² d⁻¹ | 0.74 | 1.3 m/day (1.0) |
| Sea shore, 30° S | 30.4 °C (29.4–31.2) | 244 h | 23.8 mol m⁻² d⁻¹ | 0.80 | 1.8 m/day (1.2) |
| Sea shore, 47° N | 30.1 °C (29.5–30.6) | 210 h | 18.4 mol m⁻² d⁻¹ | 0.81 | 1.7 m/day (0.8) |
| Sea shore, 64° S | 28.9 °C (28.4–29.4) | 109 h | 9.0 mol m⁻² d⁻¹ | 0.63 | 1.2 m/day (0.8) |
| Lake margin, 14° S | 29.9 °C (25.7–39.3) | 261 h | 22.4 mol m⁻² d⁻¹ | 0.66 | 0.7 m/day (0.5) |

**Shallow water runs hot late in the lunar day, and exchange with the sea decides how hot.** A closed pan 1–2 m deep at the 8° sea shore swings from 23–24 °C at the end of the night to 38–41 °C late in the day, with 69–112 hours of each lunar day above 36 °C. Trading three tenths of its water a day with the sea's 50 m mixed layer holds a 2 m shore between 27 and 36 °C at 8° and 30°, and between 28 and 33 °C at 47°; trading all of it a day holds a 2 m sea shore below 34 °C at every site. These figures hold where waves, currents and the tide keep the sea's surface mixed into that layer; a calm sea's own surface warms through the lunar day, and a shore trading water with it cools only toward the surface's warmth ([the map below](#the-waters-warmth-across-the-moon)). Water 4 m deep stays below 34.5 °C at three tenths a day. At 64° the shore stays between 26 and 34 °C in every case, the closed 1 m pan included. The percentiles are over nearly ten years of hourly values of the water's mean temperature.

| Site | Water | 1st–99th percentile | Highest | Hours a lunar day above 34 °C; 36 °C; 38 °C |
|---|---|---|---|---|
| Sea shore, 8° N | 1 m, closed | 23.1–40.7 °C | 42.5 °C | 164; 112; 54 |
| Sea shore, 8° N | 2 m, closed | 23.6–38.4 °C | 39.9 °C | 144; 69; 12 |
| Sea shore, 8° N | 1 m, all its water changed a day | 27.4–34.9 °C | 36.0 °C | 34; 0; 0 |
| Sea shore, 8° N | 2 m, a tenth changed a day | 25.2–37.4 °C | 38.6 °C | 125; 41; 3 |
| Sea shore, 8° N | 2 m, three tenths changed a day | 26.7–35.7 °C | 36.5 °C | 71; 3; 0 |
| Sea shore, 8° N | 2 m, all changed a day | 28.3–33.0 °C | 33.6 °C | 0; 0; 0 |
| Sea shore, 8° N | 4 m, three tenths changed a day | 27.8–33.6 °C | 34.3 °C | 0; 0; 0 |
| Sea shore, 30° S | 1 m, closed | 24.1–38.7 °C | 42.2 °C | 161; 66; 17 |
| Sea shore, 30° S | 2 m, three tenths changed a day | 26.9–34.9 °C | 36.8 °C | 38; 1; 0 |
| Sea shore, 30° S | 2 m, all changed a day | 28.3–32.8 °C | 33.7 °C | 0; 0; 0 |
| Sea shore, 47° N | 1 m, closed | 26.1–36.1 °C | 38.2 °C | 42; 8; 0 |
| Sea shore, 47° N | 2 m, three tenths changed a day | 28.1–33.3 °C | 34.1 °C | 1; 0; 0 |
| Sea shore, 64° S | 1 m, closed | 26.4–32.3 °C | 33.9 °C | 0; 0; 0 |
| Sea shore, 64° S | 2 m, three tenths changed a day | 27.6–30.6 °C | 31.3 °C | 0; 0; 0 |
| Lake margin, 14° S | 1 m, closed | 23.3–41.5 °C | 44.0 °C | 218; 172; 109 |
| Lake margin, 14° S | 2 m, three tenths changed a day | 25.0–38.4 °C | 41.0 °C | 183; 95; 14 |
| Lake margin, 14° S | 4 m, three tenths changed a day | 25.5–36.9 °C | 38.9 °C | 148; 35; 0 |

**Lake margins run hotter than sea shores.** The lake site has no sea beside it: its open water is a 20 m lake whose own surface warms to 39 °C under the land cell's light wind, so exchange cools its margin little, and 2 m of sheltered margin spends 95 hours of each lunar day above 36 °C. Part of the contrast lies in the two open waters: the lake's is the lake model's surface, which holds the day's warm layer, and the sea's is the climate model's mixed layer, which holds none; the same column run at the 8° sea shore's cell puts its surface at 36.95 °C in its hottest hundredth of hours (the map below). The product's every case is in [results/forcing.json](results/forcing.json).

**Oxygen leaves a shallow floor as bubbles sooner than on Earth.** Water at 29.5 °C under the design air holds 7.7 mg/L of oxygen at saturation (8.3 at 25 °C, 6.9 at 36 °C). With the water's nitrogen at equilibrium, dissolved oxygen can exceed saturation only by the pressure of the water above it, and at lunar gravity 2 m of water adds 3.2 kPa to the air's 21.2 kPa of oxygen: the floor of a 2 m column bubbles at 115% of saturation and of a 4 m column at 130%, against 192% and 285% on Earth. A productive meadow's daytime oxygen therefore leaves as bubbles where Earth's shallows would carry it into the night.

**The cloud factor pairs the run's two kinds of output.** The run writes its all-sky fluxes as window means and its clear-sky diagnostics (`rsscs`, `rstcs`) as snapshots at each output time, as the climate lane records ([exoplasim_run.py](../../../climate/gcm/exoplasim_run.py), `CLEAR_SKY`): sound for global and lunar-day averages, and no partner for a single window's mean. The forcing's own check agrees: at every site the clear-sky surface field follows the Sun at the window's end (correlation 0.992–0.997) far more closely than the window's mean Sun (0.92), and it falls to twilight's 10–12 W/m² once the Sun has set by the window's end. So the study reads the run's clear sky as a function of the Sun's elevation from those snapshots and averages it over each window before dividing the window's all-sky sunlight by it.

### Limits

The setting is a model result at T21 on 3-day windows, a first physical setting for shore life. Weather within a window, waves and the shore's own circulation are absent; the slabs' exchange is a stated range, and at the sea shores it is exchange with the climate model's 50 m mixed layer, the cool end for the water a shore receives. On the sea shores the monthly tide moves the water level by 2.5–3.7 m ([ecology, finding 18](../../../biosphere/ecology/README.md#waters)), so the fixed depths stand for the shelf below low water. The lake margin takes its land cell's 10 m wind, a wind over land's roughness that a lake's fetch would raise, which makes its daytime heat an upper-side estimate. The climate model and CM1 differ on the air near the ground, a gap the project carries as a range ([register, 8–9 November](../../decisions.md)). Fourteen tests ([test_forcing.py](test_forcing.py)) check the product's binding and its own checks, the coast mask, the Sun's fit, the clear-sky reading and the slab against its limits and against the lake model, and drive the shelved simulator with one site's series.

## The waters' warmth across the Moon

[water_map.py](water_map.py) carries the forcing above from its five sites to every T21 cell of the corrected design climate at once, so that the sites' figures become a map. One Sun serves every cell, fitted on the row nearest the equator and carried along the longitude (the cells' own fitted noons lie within 0.002 hours of it); each cell reads its own clear sky from the run's end-of-window snapshots and takes its own cloud factors, air, humidity, sky infrared, pressure and 10 m wind. The water is a mixed layer of a stated depth, stepped hourly through the ten design years: still water 1, 2 and 4 m deep that trades none of its water; water mixed to 10 m, the body of a lake's upper water; and shore water 2 m deep that trades three tenths or all of its water a day with the open water's body. On the run's water cells, its seas and large lakes, that body is the run's own sea surface, a mixed layer 50 m deep; on its land cells it is the 10 m layer, and the waters stand at the cell's mean ground level. The run writes [results/water_map.json](results/water_map.json) (schema `terluna.research.littoral-water-map/1`, with the hash of every code file and input) in about four minutes:

```sh
nice -n 10 env OPENBLAS_NUM_THREADS=1 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.water_map
```

[render.py](../../../visualization/water-temperature-map/render.py) in the visualization lane draws the product as one sheet in the manner of the moisture map's: four equal-area maps in the turned Moon's frame, every cell a box of its class's colour, and the shares of the seas' area, the sea coast's length and the rain-fed lakes' area in each class of warmth. The sheet and its manifest stay on the research drive (`/media/projectspace/terluna-research/session-tools/water_temperature_map/`), out of Git; the manifest records the product's hash and that the classes as drawn hold the product's shares.

### A mixed layer and a calm surface

A mixed layer has one temperature through its depth, while calm deep water lies in layers: through the lunar day's 354 hours of sunlight its top few metres warm above the water below, as [the seas' oxygen study](../sea_oxygen/README.md) found in its columns, whose day's mixed layer is a median 3–6 m deep. The climate model's sea is one mixed layer, so its surface, and the 10 m layer on land, give the temperature of the water's body, which is the surface's too where waves, currents and the tide keep the surface mixed into it. To show what a calm surface does beside the mixed layers, the product runs the lake model's own column at the forcing's five sites, a still 20 m pan under the cell's unsheltered wind:

| Site | The body, hottest 1% of hours | Still layer 4 m deep | The lake model's calm column, surface | Still layer 2 m deep | The column's surface: hours a lunar day above 34 °C; 36 °C |
|---|---|---|---|---|---|
| Sea cell, 8° N | 30.91 °C | 35.44 °C | 36.95 °C | 38.47 °C | 87; 21 |
| Sea cell, 30° S | 31.17 °C | 35.15 °C | 35.69 °C | 37.34 °C | 47; 5 |
| Sea cell, 47° N | 30.59 °C | 33.19 °C | 33.38 °C | 34.89 °C | 3; 0 |
| Sea cell, 64° S | 29.39 °C | 30.51 °C | 30.64 °C | 31.52 °C | 0; 0 |
| Land cell, 14° S | 32.70 °C | 35.23 °C | 36.50 °C | 37.79 °C | 110; 18 |

At every site the calm column's surface in its hottest hundredth of hours, about seven hours of each lunar day, stands 1–6 K above the body and between the still layers 4 m and 2 m deep. So the map reads in two ways. The body's figures describe open water whose surface is kept mixed; the still layers 2 and 4 m deep describe enclosed shallows of those depths, and they bracket the surface of calm deep water at the five sites where the column ran.

### What the map shows

**On average every water is mild.** The seas and large lakes cover 14.6 million km², and their upper 50 m averages 29.8 °C, between 23.7 and 32.6 °C by cell. Every still layer averages within about 2 K of the body in its cell, mostly a little below it, and no water in any case averages 34 °C or more anywhere.

**The water's body stays below 36 °C through its hottest hours.** In its hottest hundredth of hours the seas' upper 50 m stays below 34 °C over all of the seas' area (33.2 °C in the warmest cell), and a lake's upper 10 m as one mixed layer stays below 36 °C over all 2.31 million km² of rain-fed lakes (35.8 °C in the warmest cell). Shore water 2 m deep that trades all of its water each day with the sea's body stays below 36 °C along all 44,580 km of sea coast; trading three tenths a day, it reaches 36 °C or more along 21% of the coast.

**Calm and enclosed water at low and middle latitudes runs hot late in the lunar day.** In their hottest hundredth of hours the still layers reach 36 °C or more over these shares:

| Water, weighted by | Still layer 4 m deep | Still layer 2 m deep | Still layer 2 m deep, 38 °C or more |
|---|---|---|---|
| The seas, by area | 18% | 50% | 26% |
| The sea coast, by length | 23% | 51% | 27% |
| The rain-fed lakes, by area | 38% | 89% | 58% |

The highest single hours of the ten years are 44–46 °C in still water 1–2 m deep near the equator.

**Latitude sorts the warmth.** The band means below weight each cell by its area; the first figure of each pair is for the run's water cells and the second for its land cells.

| Latitude | The body, hottest 1% | Still 4 m, hottest 1% | Still 2 m, hottest 1% | Still 2 m: hours a lunar day above 36 °C | Mean 10 m wind |
|---|---|---|---|---|---|
| 0–15° | 30.4; 33.4 °C | 34.3; 36.3 °C | 36.9; 39.1 °C | 50; 127 | 2.1; 0.6 m/s |
| 15–30° | 30.8; 32.7 °C | 34.4; 35.4 °C | 36.6; 38.2 °C | 41; 89 | 2.2; 0.7 m/s |
| 30–45° | 30.6; 31.9 °C | 33.7; 34.3 °C | 35.7; 36.8 °C | 22; 42 | 2.2; 0.8 m/s |
| 45–60° | 30.3; 30.4 °C | 32.4; 32.2 °C | 33.9; 34.2 °C | 2; 2 | 2.1; 0.8 m/s |
| 60–90° | 28.5; 28.2 °C | 29.5; 29.2 °C | 30.3; 30.2 °C | 0; 0 | 1.7; 0.9 m/s |

Still water 2 m deep stays below 38 °C in its hottest hundredth of hours in every cell from the row centred on 52.6° poleward, below 36 °C from the row at 58.1° and below 34 °C from the row at 63.7°. Over the land the same water runs 1–2 K hotter than over the seas within 45° of the equator, where the land cells' wind averages 0.6–0.8 m/s against the seas' 2.1–2.2 m/s.

### Limits

The map carries the forcing's limits to every cell: a model result at T21 on 3-day windows, with no weather within a window, no waves and no circulation along a shore. The still layers trade no water, lose no heat through their floors and take the climate model's 10 m wind, which holds no sea breeze, so they are the warm end for enclosed water of their depth; over the land that wind is a wind over land's roughness, which a lake's own fetch would raise. The calm column checks the still layers at five sites, four of them sea cells, and takes the forcing's pan, whose water attenuates visible light at 0.57 m⁻¹; clearer water spreads the sunlight deeper. A column for each of the 2,048 cells, about fifteen seconds apiece and some three hours on three threads, would turn that check into a map of the calm surface. The 10 m layer under the run's near-surface air averages 0.7 K below the run's own sea surface on the water cells (0.2–1.3 K by cell), which measures how closely the layers' exchange follows the climate model's. The coast's shares take each land cell's coast length and the rain-fed lakes' shares each land cell's lake area from the moisture map's settings, the coast reading the mean of its neighbouring water cells. The climate model and CM1 differ on the air near the ground, a gap the project carries as a range. Twelve tests ([test_water_map.py](test_water_map.py)) check the product's binding and shares, the array form of the lake model's bulk exchange against the lake model, the many-cell stepping against the forcing's one-site slab on a made-up climate, the agreement with the forcing's five sites, and the calm column's place between the still layers.

## The founders' review

[founders.md](founders.md) is a bounded literature review, 44 searches and two papers' pages read on 9 November, of the species that could found the littoral community; [sources.json](sources.json) gives each source with how far it was read. Its placements and its first community are proposals. In brief:

- **Heat sorts the candidates first.** Upper heat tolerance barely answers selection in a fish (0.04 °C a generation in zebrafish), while algae, water fleas, growth and timing answer within tens of generations, so the founders are animals that already stand the shore's heat.
- **The gentlest shores lie between 45° and 65°**, on shelves from 2 m, open to the sea: 27–34 °C and 110–210 dark hours. A first self-sustaining community is best tried there.
- **The hard, clear, phosphorus-poor design water favours stonewort meadows and a food web fed from the floor**, once a rich first decade has passed; fast colonisers fill that decade and recede.
- **Plants first, then crayfish, then their predators** holds against Earth's lakes on two conditions: meadows established and fenced before crayfish arrive, and predators that take crayfish small and moulting arriving with them.
- **The night's stores already exist** as turions, bulbils, tubers and resting eggs; their cue, moved to the lunar dusk, is the breeding target. No stonewort has been measured through 110–260 warm dark hours.
- **Founding stocks mixed from many sources** adapt and diversify, and aggressive founders bring spoiled catches, fluke-hosting snails, stripped small fishes and fouled intakes to design around.
- **Leads by guild:** warm-water stoneworts, greater duckweed and heat-tolerant azolla; water fleas and bladder snails; redclaw crayfish on open shores; spangled perch, jade perch and tilapia of the GIFT lines; stocked eels against the crayfish; mixed cichlid stocks for the crater lakes.

## The community simulator, shelved

[biosphere/ecology/littoral.py](../../../biosphere/ecology/littoral.py) steps a shallow community through the lunar month for an ensemble of settings at once. It was built on 9 November as this study's first approach and shelved the same day in favour of a review of founder species: a simulation's answers follow its species' rates, and those rates were placeholders. It is kept, with its tests, for the time when founders with measured rates exist.

**What it holds.** One mixed box of water over a lit floor, trading water with the open water beside it. Phytoplankton, a periphyton mat and a stonewort meadow with a carbohydrate reserve grow on light and phosphorus, with nitrogen fixed at a carbon cost where the water runs short; the meadow can idle its upkeep in the dark, lays calcite with its net photosynthesis and shelters juvenile animals. Detritus decays with or without oxygen, yielding methane that methanotrophs partly oxidise. Consumers eat by preference on a saturating response, build tissue at fixed carbon, phosphorus and nitrogen ratios, and respond to temperature through a thermal performance curve with heat and hypoxia mortality; each is one pool or juveniles and adults in De Roos and others' (2008) stage-structured form, where juveniles mature at a rate set by their growth and mortality and adults turn their surplus into offspring. Crayfish clip more meadow than they eat and stir the floor; oxygen exchanges with the air, with the open water, and leaves as bubbles at the floor's pressure.

**What its tests establish.** Sixteen test cases ([test_littoral.py](../../../biosphere/tests/test_littoral.py)) show that carbon, phosphorus, nitrogen and oxygen close to rounding through every exchange, in closed and open boxes and through anoxia, heat and starvation; that the maturation rate meets its analytic limits and reproduces a cohort's survival to maturity; that the bubble point follows gravity; that halving the step changes a lunar day's outcome by under 3%; and that each mechanism acts in its direction: light and darkness on oxygen, phosphorus on growth, idling on the meadow's night, heat and hypoxia on survival, clipping and stirring by crayfish, shelter and predation on juveniles, the breeding cue on births.

**What an example run does.** `example_community()` gives a zooplankter, a grazing snail, a crayfish and a predatory fish round placeholder traits, and `drive()` steps them through a plain lunar day (`lunar_forcing()`) or a site's series from the forcing above. `python3 -m biosphere.ecology.littoral` seeds them in sequence over three years in four boxes at a productive margin and prints each pool year by year, in about a minute and a half. In that run the grazers strip the periphyton, the crayfish strip the meadow within about a year of arriving and then decline, and the fish remove the crayfish from two of the four boxes and nearly from a third, the deepest keeping crayfish and fish under plankton-rich, turbid water. The sequence shows the pieces interacting; its timing and outcome belong to the placeholders.

**What it needs before use.** Traits for each founder from measurements, with their sources. A check against an Earth water whose history is known, such as a crayfish introduction's effect on a lake's plants. A propagule bank or immigration, since a pool grazed to nothing stays at nothing. The monthly tide's changing depth, plants' own periphyton, feeding by sight against the night's light, and more than one box where shore and open water differ. Until then its numbers describe its assumptions.
