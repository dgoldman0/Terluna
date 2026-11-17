# The atlas round

The atlas round builds the Open Moon's atlas down to its districts. It began on 16 November 2026 on its own branch, `domain/geography`, and the register's section of that day sets its frame ([the atlas round's frame](../../decisions.md#the-atlas-rounds-frame-16-november-2026)). Each step stands on the one before it, so the round works from the ground up:

1. **The surface redrawn under its water** (below): the floor and the seas' levels once the water's own gravity and the crust's sag and creep have settled, barring erosion.
2. **The drainage, the river network and the globe** on that floor, with the three upper seas filled and spilling and Moscoviense's water divided for the most diversity: the salt sea's share over its natural pass, and the rest through a cut in its northern ridge to Humboldtianum and on to Smythii–Marginis.
3. **The design climate's run on the final surface** at the Moon's settled turn under the final water, 6.97°, about 3–7 hours of machine time once the turned frame and surface are built. Whether the near-side sea may part at its pass is decided on it, and the globe is rebuilt on it.
4. **The south sea's ocean runs** on the new basin and climate, the paused oxygen run among them.
5. **The places**: districts defined by significance and wholeness, after a quick review of the heritage register's written evaluations; the young-Moon landforms to keep, chosen as a set; satoyama-like landscapes, and satoumi on the salt sea's shores; the riparian share of the land against Earth's, with the falls counted on finer heights; the lakes among fresh and soda; and the three-dimensional atlas on the present design climate.

## The surface redrawn under its water

The atlas fills LOLA's floor to one level of the dry Moon's geoid, so that the water below that level covers 28% of the surface. The water changes the floor and the level both: its weight bends the crust down, the mantle creeps away beneath the load, and the water's own gravity lifts the sea surface over deep water. Geography's crustal-loading study measured these for one common level on 4 November ([geography](../../../geography/README.md#the-crusts-response-to-the-water)). The round's frame adds what that study did not have: the three upper seas, Moscoviense, Orientale and Smythii–Marginis, filled to their rims above the common level, and the crust once its response has settled.

[redrawn_surface.py](redrawn_surface.py) solves that surface. Every cell below a common level is water, the level set so that those cells cover 28% of the surface as the atlas's do. Each upper sea stands 1 cm under its rim, the lowest pass from its water to another sea's, over all the land its water joins below the rim, as the [drainage network study](../drainage_network/README.md) defines both. The lakes keep the depths that the turned Moon's corrected drainage gives them and weigh on the crust as a fixed load wherever no sea covers them. The Moon answers through load Love numbers: the crustal-loading study's elastic response at every degree to 359, and the creep that the geology lane computes for the present Moon at degrees 2–50 ([geology](../../../geology/README.md)), its change from the water's arrival to the case's time. The sea-level equation is iterated as geography's study iterates it, and after each solution the catchments and passes are found again on the deformed heights, until no rim moves by more than 5 cm. The results go to [results/redrawn_surface.json](results/redrawn_surface.json) (schema `terluna.research.atlas-round-redrawn-surface/1`) and the grids to the research drive.

The round's surface stands at the crust's settled state, barring erosion (register, 16 November): LOLA's landforms as they are, under the sag and creep the water has finished causing. Under the present Moon's central temperatures the creep never quite stops. It grows by about the same amount every tenfold of time: the degree-2 load Love number h′ gains −0.09 by 1,000 years, −0.32 by 10 million and −0.53 by a billion. The geology lane holds today's temperatures, which it takes to describe the Moon for about 100 million years; after that the cooling interior stiffens the mantle and keeps the sag already made. The design case is therefore the central profile 100 million years after the water arrives, and the other cases bound it or trace the way there:

| Case | What it stands for |
|---|---|
| Atlas | The dry Moon's floor and geoid, deformed by nothing; it reproduces the drainage network study's rims and geography's sills |
| Elastic | The first years after the water arrives |
| Central, 1,000 years | The first millennium, the horizon of the build |
| **Central, 100 million years** | **The settled state, the design** |
| Central, 1 billion years | The most the present temperatures could add |
| Stiff, 1 billion years | The cool profile at its most |
| Soft, 10 million years | The warm mantle, carried as a stated risk; it settles within about 10,000 years |

All seven ran on 16 November in 25 minutes on one thread, 1.6 GB at most. Six settled within four rounds of finding the rims again; the warm mantle's did not (below).

| | Atlas | First years | First millennium | **Settled** | Billion years, central / stiff | Warm mantle |
|---|---|---|---|---|---|---|
| Common level against the atlas's | 0 | +36.8 m | +36.5 m | **+39.2 m** | +41.2 / +36.6 m | +61.1 m |
| Water, seas and lakes | 15.23 | 16.16 | 16.18 | **16.45** | 16.62 / 16.21 | 17.88 million km³ |
| Beyond the atlas's 13.37 million km³ | +1.86 | +2.79 | +2.81 | **+3.08** | +3.25 / +2.85 | +4.51 million km³ |
| Moscoviense's rim above the common level | 3,545 | 3,424 | 3,424 | **3,449** | 3,461 / 3,426 | 3,509 m |
| Orientale's rim | 2,342 | 2,351 | 2,352 | **2,341** | 2,335 / 2,352 | 2,283–2,337 m |
| Smythii–Marginis's rim | 767 | 726 | 726 | **720** | 720 / 724 | 735 m |
| South Pole–Aitken Sea, area | 2.31 | 2.40 | 2.40 | **2.41** | 2.41 / 2.41 | 2.41 million km² |
| its volume | 4.82 | 5.62 | 5.64 | **5.80** | 5.86 / 5.66 | 6.34 million km³ |
| its deepest water | 6,795 | 7,217 | 7,227 | **7,335** | 7,386 / 7,241 | 7,757 m |
| Near-side Sea, area | 5.19 | 5.08 | 5.08 | **5.08** | 5.08 / 5.09 | 5.12 million km² |
| Tenebrarum sill, under the sea | 75 | 313 | 319 | **368** | 385 / 326 | 460 m |
| Crisium passage | 155 | 203 | 205 | **220** | 230 / 207 | 338 m |
| Fecunditatis's western channels | 109 | 105 | 105 | **108** | 111 / 105 | 108 m |
| The near side's narrowest join, under the sea | 109 | 72 | 72 | **69** | 68 / 71 | 75 m |

The water's own gravity does most of the redrawing, within years of the water's arrival. Already in the first years the common level stands 36.8 m above the atlas's, the South Pole–Aitken Sea is 422 m deeper at its deepest and its sill to the Tenebrarum sub-basin lies 313 m down instead of 75. The creep to the settled state then adds less: 2.4 m to the common level, 118 m to the south sea's deepest water, 55 m to the Tenebrarum sill and 17 m to the Crisium passage, as the floor under the south sea's deepest water falls 189 m in all (24 m within years). The first millennium differs from the first years by under a metre in the level and by 21 thousand km³ of water, so the creep that matters runs over the following hundred thousand to hundred million years.

The water moves towards the south. Holding 28% at the common level keeps the area the same, so the shore moves both ways: at the settled state the common level floods 179,000 km² of land and gives up 178,000 km² (a read-only look, `atlas-round/redrawn_surface/look/shore_moves.py` on the research drive). The South Pole–Aitken Sea gains 76,000 km² of new sea along its shores and gives up none, as its own gravity draws water to it and its floor sinks most, while the Near-side Sea gives up 107,000 km² of its shallow margins and floods 24,000; the hollows held at the common level lose 21,000 km² net. Along today's shore the sea stands a median 2.8 m higher against the land, from 56 m lower (5th percentile) to 255 m higher (95th percentile). Nearly all of this comes within years: 69,000 km² gained by the south sea and 106,000 given up by the near side in the first-years case.

Measured from the redrawn common level, Moscoviense's rim stands 96 m lower than on the dry Moon and Smythii–Marginis's 47 m lower, while Orientale's holds within a metre. Part of each fall is the common level's own rise of 39 m; the rest is the ground's sag and the water's gravity about each rim. The three seas' deepest water is 180, 140 and 36 m deeper than on the dry Moon. Filled to their rims at the settled state they cover 844,000, 298,000 and 599,000 km² and hold 1.32, 0.38 and 0.36 million km³ beyond what the atlas holds in the same places at the common level, 2.06 million km³ together. Each joins several seas at its rim, so which way its overflow runs is the drainage's to find on the redrawn floor.

The settled surface takes 16.45 million km³ of water, 3.08 million km³ beyond the atlas's 13.37 million (10.72 million in the seas, 2.65 million in the lakes), a global layer of 81 m. The register's estimate of about 3.0 million km³ for the seas' 28% with the upper seas on top holds. The near side stays one sea: its narrowest join, now the neck across western Procellarum to Cognitum, Humorum and Nubium, lies 69 m under the settled sea, against the atlas's 109 m at Fecunditatis's western channels; whether the near-side sea's own balance holds it there is the design climate's run to show.

The bounds stay close to the settled state except the warm mantle. A billion years under the present temperatures would add 0.17 million km³ and 2 m to the level, and the stiff profile stays within 0.04 million km³ of the first millennium. The warm mantle, the least likely of the three profiles, takes 4.51 million km³ beyond the atlas with the common level 61 m up and the south sea 7.76 km deep, and there Orientale's rim does not settle: filled to one of two passes, its weight sinks the ground until the other is the lower, so its rim cycles between 2,283 and 2,337 m and the case is reported with that range.

What the round's next steps take from here: the deformed heights (the settled case's U and N on the research drive, with its water depths and the upper seas' extents), the common level and the three rims for the drainage on the redrawn floor; the water to deliver, 3.08 million km³ beyond the atlas; and the gaps between the first years and the settled state for holding the surface near its settled shape sooner (below).

## The drainage on the redrawn floor

The second step routes the design years' rain over the floor the first step redrew. [redrawn_drainage.py](redrawn_drainage.py) adds the settled case's change of the ground against the sea surface, U − N, to LOLA's heights at 16 pixels per degree, and takes the redrawn common level, 39.17 m above the atlas's, below which 28% of the surface lies. On those heights it runs the [drainage network study](../drainage_network/README.md)'s analysis unchanged: its geometry, rims, fills, chain, history and long profiles, on the design years' rain, on the same rain with the open water's evaporation raised 14.6% until the routing's budget closes, and at the ends of the rain range the climate lane carries. It then routes the state the round's frame sets and searches for the gravity routes by which a full upper sea's water could go elsewhere. The climate is still the design years' on the turned Moon at 8°, as the atlas's waters made it; the run on the final surface (step 3) replaces it. The record is [results/redrawn_drainage.json](results/redrawn_drainage.json) (schema `terluna.research.atlas-round-redrawn-drainage/1`). The design state's routing goes to the research drive for the satoyama map and the globe: every cell's receiver, discharge and drained area, with the lakes and the seas.

**The network keeps its shape on the redrawn floor** *(computed)*. The rims found on the redrawn heights reproduce the first step's to the centimetre, and every overflow goes where it went on the dry Moon: Moscoviense's to the south sea, Orientale's and Smythii–Marginis's to the near-side ocean. The floor has 274,986 catchments against the dry floor's 275,888. Heights are above each floor's own common level; years to fill are on the routing's budget, then the closed one:

| | Dry floor (network study) | Redrawn floor |
|---|---|---|
| Moscoviense: rim, area full, years to fill | 3,545 m, 837,000 km², 790 and 880 | 3,449 m, 844,000 km², 810 and 905 |
| Orientale | 2,342 m, 289,000 km², 795 and 954 | 2,341 m, 298,000 km², 819 and 988 |
| Smythii–Marginis | 767 m, 568,000 km², 453 and 537 | 720 m, 599,000 km², 474 and 563 |
| The south sea's rim; Humboldtianum's | 1,455 m; 938 m | 1,401 m; 942 m |
| Moscoviense's overflow | 694 km through 8 hollows, steep for 48 km, largest drop 654 m | 542 km through 7, steep for 44 km, 666 m |
| Orientale's; Smythii–Marginis's | 1,702 km through 19; 2,383 km through 8 | 1,806 km through 18; 2,419 km through 8 |
| The near-side ocean left to itself, closed budget; three quarters of the rain | settles 120 m down; 129 m | settles 69 m down; 133 m |

The near-side ocean stops higher on the closed budget because the floor the water redrew raised its narrowest join: the neck across western Procellarum lies 69 m under the settled sea (step 1), and there the ocean's loss to the air first outruns its supply, so left to itself it would part at that neck. Left to its rivers and Moscoviense's whole overflow, the south sea again fills to its rim and overflows into the near-side ocean with its salt, at year 4,665 on the routing's budget, against 4,565 on the dry floor.

**The design state** *(computed)*. With the three upper seas full and spilling and the other three held at the common level, each full sea's surplus is put in at the first cell past its pass and routed down with the rivers, filling the hollows on its way. Flows in m³/s, the routing's budget first and the closed one second:

| Sea | Level | Area | Its own rivers | From the full seas | Loss to the air beyond its rain | Balance |
|---|---|---|---|---|---|---|
| Moscoviense | +3,449 m | 844,000 km² | 54,400; 50,600 | | −1,300; 2,800 | passes on 55,800; 47,800 |
| Orientale | +2,341 m | 298,000 km² | 15,600; 13,900 | | 700; 2,100 | passes on 14,800; 11,800 |
| Smythii–Marginis | +720 m | 599,000 km² | 30,300; 27,400 | | −3,500; −600 | passes on 33,800; 28,000 |
| Near-side ocean | 0 | 5.08 million km² | 59,400; 51,300 | 48,500; 39,500 | 82,900; 106,700 | +25,100; −15,900 |
| South Pole–Aitken Sea | 0 | 2.41 million km² | 31,700; 28,100 | 55,800; 47,800 | 47,300; 56,400 | +40,200; +19,400 |
| Humboldtianum | 0 | 279,000 km² | 200; 0 | 0 | 6,500; 7,600 | −6,300; −7,600 |

Moscoviense's overflow brings the salt sea 1.7–3.6 times what it lacks. Held at the common level, the salt sea takes 15,600 m³/s of Moscoviense's 55,800 on the routing's budget and 28,400 of 47,800 on the closed one, so 40,200 or 19,400 m³/s must go another way, against 41,900 and 21,500 on the dry floor. Orientale's and Smythii–Marginis's overflows give the near-side ocean 25,100 m³/s more than it loses on the routing's budget and leave it 15,900 short on the closed one. Humboldtianum gets no overflow and needs 6,300–7,600 m³/s to stay at the common level. On less it settles lower, at the level where its loss to the air meets its supply:

| Water brought to Humboldtianum | Routing's budget | Closed budget |
|---|---|---|
| none | its deepest hollow, 2,409 m down, 4,600 km² | dry |
| 1,000 m³/s | 1,964 m down, 52,000 km² | 2,176 m down, 36,000 km² |
| 2,500 m³/s | 764 m down, 114,000 km² | 986 m down, 93,000 km² |
| 5,000 m³/s | 303 m down, 226,000 km² | 469 m down, 179,000 km² |
| 6,300 and 7,600 m³/s | the common level, 279,000 km² | the common level, 279,000 km² |

**Most of the hollows the atlas holds at the common level dry out on the design climate's balance** *(computed)*. Apart from the six seas, 2.23 million km² of the redrawn floor lies below the common level in hollows whose passes stand above it; the first step counted them full to the common level, 1.58 million km³. Routed by their own balance, as geography's drainage routes every hollow, with no water carried to them (the register's rule for the seas), hollows over 485,000 km², mostly within 30° of the turned equator, stay wet and fill above the common level to their own passes. Hollows over 1.75 million km², 1.48 million of it poleward of 30°, would dry, and with them 1.25 million km³ of the first step's water. The common level keeps its definition, the level with 28% of the surface below it, and the six seas' water covers 8.39 million km² of the floor below it, 22.1% of the surface; with the upper seas full and the lakes, seas and lakes cover 38.3%, against the atlas's 28% and its drainage's 9.8% of lakes. The climate run on the final surface takes the hollows at their own levels, and the water to deliver is counted again on it.

**The gravity routes out of the full seas** *(computed between 1.9 km cells)*. A full sea's water can go elsewhere only through a cut: a channel at the sea's level, cut down through the ground that stands above it. Two searches measure the cuts. The first takes each pass on the divide round the land draining to a full sea into land draining to another sea, and measures the deepest cut a channel needs from the shore over that pass and down the far side until the water runs on by itself. The second lets the water fall through hollows and cut again on the way, and finds the route to each sea whose deepest cut is shallowest, then among those the one whose cuts add up to the least. Places are on today's Moon:

| From | To | Deepest cut | Where | The way on |
|---|---|---|---|---|
| Moscoviense | the south sea | none (its pass) | 10.7° S, 175.1° E | 542 km through 7 lakes |
| Moscoviense | Humboldtianum | 438 m, one ridge cut 40 km long | 35.0° N, 139.7° E, 9 km from the shore | 1,979 km through 22 lakes, 305 km of it steep |
| Moscoviense | Smythii–Marginis, and through its overflow the near-side ocean | 438 m at the same ridge, then 89 m and 56 m from lakes beyond it | 35.0° N, 139.7° E; 34.2° N, 133.0° E; 31.2° N, 131.5° E | about 1,900 km to Smythii–Marginis, 4,300 km to the near-side ocean |
| Smythii–Marginis | Humboldtianum | 440 m, through Humboldtianum's own rim from water fallen to 502 m, with nine more cuts, ten adding up to 2,500 m | 45.5° N, 76.1° E for the deepest | about 1,050 km |
| Orientale | Smythii–Marginis, the south sea or Humboldtianum | 651 m | 33.1° S, 80.3° W | 3,400–10,000 km |
| Moscoviense; Orientale; Smythii–Marginis | Orientale; Moscoviense; the south sea | 1,056 m; 1,680 m; 707 m | | |

Moscoviense's northern ridge is the one modest way out of any full sea. Its cut is deepest where it crosses the divide, 438 m below a pass 3,887 m above the common level, and the passes beside it in the same ridge, 171, 56 and 16 m above the water, lie along the same trench. Past it the water falls 3,449 m to Humboldtianum through the dry northern far side, 39–57° of turned latitude, where the valleys get under 400 mm of rain a year. Two more cuts of 89 and 56 m, from lakes 3,028 and 2,627 m above the common level beyond the ridge, take it to Smythii–Marginis instead, and Smythii–Marginis's overflow on to the near-side ocean. From Smythii–Marginis, the nearest full sea to Humboldtianum, no route needs less than a 440 m cut, since the water falls to 502 m before it meets Humboldtianum's rim at 942.5 m; and its single shallowest cut from its own shore, 1,050 m deep, would be 1,755 km long. Orientale has no way out shallower than 651 m.

**The first satoyama map** stands on this drainage and the present design climate ([satoyama.md](satoyama.md#the-first-candidate-map), [satoyama_candidates.py](satoyama_candidates.py), [results/satoyama_candidates.json](results/satoyama_candidates.json)). Of 549,000 side valleys of 10–100 km², 50,600 have the order down the slope and the rain, 1.57 million km², a third of the hybrid fifth and nearly all within 30° of the turned equator. They could carry 204,000 km² of flooded fields, which give 0.89–1.51 Tg of methane a year: within the room under 1.87 ppm at the removal ledger's low-cost end and four times it at the high-cost end. The salt sea has 2,600 km² less than 10 m deep along 37,800 km of shore.

## The drainage with the works chosen

The drainage left three choices, made on 16 November (register, the atlas round's frame). Moscoviense's water is divided at two weirs: one at its natural pass sends the salt sea what holds it at the common level, and the northern ridge cut takes the rest. Humboldtianum takes its share from that cut, and the remainder goes on through the two further cuts to Smythii–Marginis and, with Smythii–Marginis's overflow, to the near-side ocean. A sill 24.7 m high at Moscoviense's pass and cuts 9.2 m and 6.5 m deep at Orientale's and Smythii–Marginis's hold each upper sea at its settled level from the start; at the settled state they change nothing in the routing. [governed_drainage.py](governed_drainage.py) routes the design rain over the redrawn floor with these works in place. Each branch's water is put in where it falls out of a cut into a hollow it then runs on from by itself, and taken out of the lake its next cut leaves from, so the rivers carry it where it runs. Each weir is set by routing again until the salt sea and Humboldtianum receive their need after the lakes on the way take their share; three routings settled them on each budget. The record is [results/governed_drainage.json](results/governed_drainage.json) (schema `terluna.research.atlas-round-governed-drainage/1`) and the routing, in the design product's form, `atlas-round/redrawn_drainage/governed_routing.npz` on the research drive. Flows in m³/s *(computed)*:

| | Routing's own budget | Closed budget |
|---|---|---|
| Moscoviense passes on | 55,800 | 47,800 |
| over its natural pass, to the salt sea | 15,600 | 28,400 |
| through the northern cut, for Humboldtianum (arriving) | 6,300 (6,300) | 7,700 (7,600) |
| through the northern cut, for Smythii–Marginis (arriving) | 33,800 (33,400) | 11,700 (11,100) |
| Smythii–Marginis passes on to the near-side ocean | 67,200 | 39,100 |
| The near-side ocean's balance at the common level | +58,500 | −4,800 |
| The level its balance would bring it to | its rim, 720 m up, where it joins Smythii–Marginis over 8.46 million km² | 62 m down, over 4.83 million km² |
| Lakes the works fill | 17,800 km² in 16 | 24,000 km² in 28 |

The northern river leaves the ridge cut at 35.2° N, 138.7° E with both shares, 19,500–40,200 m³/s, and within 53 km, the fifth lake on its way, enters a crater lake of about 5,500 km² standing 3,028 m above the common level and up to 2.1 km deep. Humboldtianum's share runs on from it, 1,950 km through 21 lakes in all, with 60 km of steep reaches in 21 places, and reaches Humboldtianum at 57.2° N, 111.7° E. Smythii–Marginis's share leaves the crater lake through the 89 m cut, runs 40 km down a side valley to a second crater lake of 2,500 km² at 2,627 m, leaves it through the 56 m cut and runs 2,400 km through 16 lakes to Smythii–Marginis at 25.8° N, 103.0° E. The lakes on the way take 50–190 m³/s of Humboldtianum's share and 450–590 m³/s of Smythii–Marginis's. Over Moscoviense's natural pass the overflow river to the salt sea keeps its course, 542 km through seven lakes, with 28–59% of its natural flow.

**On the closed budget the division keeps the near-side ocean whole** *(computed)*. Smythii–Marginis's share of Moscoviense's water brings the near-side ocean within 4,800 m³/s of its loss at the common level. It would then settle 62 m down with 95% of its area, 7 m above the neck across western Procellarum at which it would part without the division. On the routing's own budget it rises to its rim and joins Smythii–Marginis. The climate run on the final surface decides between them, and with it whether the near-side sea may part.

**The works fill lakes in dry country** *(computed)*. The new rivers fill 17,800–24,000 km² of lakes that the design state leaves dry or shrunken, 16–28 water bodies, mostly between 30° and 45° of turned latitude, where the valleys get under 400 mm of rain a year.

**The lakes fill more than the atlas's** *(computed)*. The governed routing's lakes hold 4.84 million km³ over 5.02 million km², against the 2.65 million km³ of the drainage the climate runs used. The wet hollows below the common level fill to their own rims once their own balance governs them, where the atlas held them at the common level, which only an outlet there could do. Mendeleev's basin, its floor 2.3 km below the common level, becomes a lake of 76,000 km² standing 3.8 km above it and holding 212,000 km³, and several far-side basins near the equator do the same. With the seas' 12.38 million km³ (the upper seas full), the final water comes to 17.2 million km³, 3.85 million beyond the atlas's 13.37 million, where the first step counted 3.08 million with the hollows at the common level and the lakes at their old depths.

## The settled turn on the final water

The design climate's run takes the Moon's settled turn (register, 16 November), put at 6.15° by the geology lane's orientation series, which loads the Moon with the atlas's water. [settled_turn.py](settled_turn.py) reproduces that series' settled turn exactly, 6.146°, from the series' own load and response, then loads the Moon with the final water instead: the redrawn seas with the upper seas full and the governed lakes, 1.72×10¹⁹ kg, block-averaged to the series' 4 px/deg grid. Through the same response (geology/orientation.py's evolve, the central profile's k′₂ and tidal k₂, the bulge re-forming about the moving axes) *(computed)*:

| | The series' load (atlas's water) | The final water |
|---|---|---|
| Settled turn, central at 100 million years | 6.15°, toward 166.7° W | 6.97°, toward 150.8° W |
| The point under Earth | 6.08° N, 4.98° E | 6.14° N, 0.84° E |
| The first years; the stiff profile at a billion years; the soft at ten million | 8.02°; 7.33°; 1.89° | 9.03°; 8.28°; 2.19° |

The final water turns the Moon 0.8° more at the settled state and leans its pole toward 150.8° W instead of 166.7° W; the final seas alone, without the lakes, would turn it 7.20°. The record is [results/settled_turn.json](results/settled_turn.json), with the settled frame for the turned surface. The design climate's run takes this turn, 6.97° (register, 16 November).

## The surface for the design climate's run

[turned_surface.py](turned_surface.py) builds the surface the run on the final surface stands on, as `climate/gcm/turned.py` built the atlas's for the runs so far, but from the redrawn floor and the governed routing. Every 16 px/deg cell carries its sea's level, its lake's or none, on the routing's own budget, the design climate's rain. A 4 px/deg cell then takes its share of wet cells as its water fraction, the mean of its sixteen heights as its height, and the means over its wet cells as its water's depth and surface, so the height of its dry part is exactly the mean of its dry cells. The surface is turned with settled_turn.py's rotation, its spin axis tilted 1.707° to the ecliptic in place of the elastic state's 1.613°, and averaged onto the GCM grids as before, with one change: the dry ground's height is averaged as it stands, because the redrawn floor has dry hollows below the common level, where every surface so far held such hollows as water. The runner's experiment A28S reads it, named by the governed routing's content *(computed)*:

| | The surface of the runs so far | The final surface |
|---|---|---|
| Seas; lakes; all water | 28.00%; 10.34%; 38.34% | 25.08%; 13.23%; 38.31% |
| PlaSim's T21 mask: land cells, land share | 1,381; 0.6160 | 1,419; 0.6169 |
| Land in the 10° bands north of 60° (turned) | 0.81–0.89 | 0.98–1.00 |
| Land in the 10° bands from 30° S to 10° N | 0.55–0.60 | 0.44–0.53 |
| Land in the 10° bands from 70° S to 30° S | 0.60–0.68 | 0.68–0.76 |

The same share of the Moon lies under water, but the water moves toward the equator. The far-side and polar hollows that the atlas held at the common level dry on the design rain, and the lakes of the equatorial rain belt fill above it, so the high latitudes keep almost none and the tropics gain. On the T21 grid, 21 land cells, 0.9% of the Moon, have mean ground below the water level, 398 m on average and 1,211 m at the deepest; the runner sets PlaSim's ground no lower than the water level and raises them to it. The record is [results/turned_surface_61afb0be43d71adf.json](results/turned_surface_61afb0be43d71adf.json).

The run starts from the corrected design run's year-49 climate carried onto the new surface. `climate/gcm/carry.py` now also turns the ground from one turned frame into another, here from the elastic 8.05° into the settled 6.97°, and records the state it carried, so the equinox keeps regressing from that run's count (`climate/gcm/precession.py`). On the T21 grid 199 cells turn from land to water and 237 from water to land, against 535 in all when the turned surface was first carried on 5 November, and each cell's state comes from 4.0° away on average.

## The design climate's runs on the final surface, and its drylands

The run at the corrected design run's sunlight, 0.938, warmed from 301.09 K to 301.67 K over its first nine model years and headed for about 302.0 K, so it was stopped after year 8 and run again from there at 0.930, the correction of 31 October from its years leaning cool. At 0.930 the climate cooled for some 35 years and settled at 300.31 K over years 35–49, its top of the air in balance; the surface turns out more sensitive to sunlight than the correction assumed, about 0.21 K per 0.001 of the sunlight scale, and slow to settle. Which sunlight the design climate takes is open (below). The runs' scripts and rules are on the research drive (`atlas-round/step3/`).

[drylands.py](drylands.py) reads the land as the moisture map reads it, rain over Priestley and Taylor's demand in UNEP's classes, each window's own land by area *(computed)*:

| | Final surface, 0.930, years 35–49 (300.3 K) | Final surface, 0.938, years 1–8 (301.5 K, still warming) | Atlas surface, corrected design, years 10–19 (301.0 K) |
|---|---|---|---|
| Hyper-arid; arid; semi-arid | 4.9%; 22.7%; 22.8% | 3.7%; 22.9%; 23.7% | 2.7%; 24.0%; 26.0% |
| Dry sub-humid; humid | 6.6%; 42.9% | 6.8%; 42.8% | 8.0%; 39.3% |
| All drylands | 57.1% | 57.2% | 60.7% |
| The waters' surface (all water cells) | 302.51 K | 303.69 K (years 6–8) | 302.90 K |

The final surface is less dry than the atlas's, 3.6 points less dryland, but has more hyper-arid land, the polar hollows that dry on the design rain. Over the 1.2 K between the two windows on the final surface the drylands hardly change: the cooler climate has 1.1 points less arid and semi-arid land and 1.2 points more hyper-arid. At 300.3 K the waters stand 0.4 K cooler than the design climate's on the atlas surface.

[river_losses.py](river_losses.py) asks what the rivers add, since PlaSim carries river water to the sea without losing any on the way. In the drylands, where open water loses more to the air than the rain brings (80% of the land on the routing's own budget), rivers of at least 100 m³/s cross only 2,570 km (2,990 on the closed budget); the large rivers' courses run mostly through lakes, which the surface already holds. Their channels, widened for lunar gravity, lose 7–45 m³/s, a belt of plants drinking their beds' seepage 17–94 m³/s, and the floodplains within 5 m of their river, 9,050 km² in 777 pieces of at most 704 km², about 130 m³/s; all of it is a third or less of the 32–46 km³ a year the ecology lane found shore greening to transpire. The young surface, barring erosion, has few flat plains beside its rivers. Carried onto the run's T21 grid, the floodplains add at most 4.1% of a cell's area in water, against the 47% a land cell needs to turn to water in PlaSim's mask, and the green land they and the corridors make covers 0.02–0.03% of the land, leaving the arid and semi-arid share at 45.5%.

## The drainage on the final surface's climate, and the near-side sea's neck

[final_drainage.py](final_drainage.py) routes the run at sunlight 0.930 (years 35–49, 300.3 K at the ground) over the redrawn floor with the works chosen on 16 November, as [governed_drainage.py](governed_drainage.py) routes the present design climate. The run's climatology is carried from its settled turned frame onto today's ground with the final surface's rotation, as `climate/gcm/turned_climatology.py` carries one (`climatology_A28S_301k_design_l15_snow_s930_settled_y35-49_today.npz` in `climate/gcm/products`, ignored by Git). The record is [results/final_drainage.json](results/final_drainage.json), and each budget's routing goes to `atlas-round/final_drainage/` on the research drive. The drainage on the present climate stays as it is beside it.

**The run's own surface brings the routing's budget nearer to closing** *(computed)*. The open water's evaporation needs raising by 7% to close it, against 14.6% on the present climate, now that the run's surface holds the upper seas full and the governed lakes. Flows in m³/s, the present design climate first and the run on the final surface second:

| | Routing's own budget | Closed budget |
|---|---|---|
| Moscoviense passes on | 55,800; 46,700 | 47,800; 42,500 |
| over its natural pass, to the salt sea | 15,600; 22,500 | 28,400; 28,300 |
| through the northern cut, for Humboldtianum | 6,300; 4,400 | 7,700; 5,300 |
| through the northern cut, for Smythii–Marginis | 33,800; 19,800 | 11,700; 8,800 |
| Orientale passes on | 14,800; 15,900 | 11,800; 14,200 |
| Smythii–Marginis passes on to the near-side sea | 67,200; 48,200 | 39,100; 34,400 |
| The near-side sea's own rivers | 59,400; 50,300 | 51,300; 46,000 |
| its loss to the air beyond its rain, at the common level | 82,900; 90,400 | 106,700; 101,800 |
| its balance there | +58,500; +24,000 | −4,800; −7,200 |
| The level its balance brings it to | its rim, 720 m up, on both | 62 m down; 68.7 m down, at its neck, where it parts |

The cuts' routes come out as before. The rain between 30° and 58° N and 100° and 140° E today, where the northern river runs, comes to 774 mm a year against 211 on the present climate, now that Moscoviense stands full in the run's surface, so the lakes the works filled on the present climate (17,800–24,000 km²) stand full on the new climate's own rain.

**On the closed budget the near-side sea parts at its neck** *(computed)*. [near_side_neck.py](near_side_neck.py) works out the two parts' balance and the neck between them. The neck is a sill 5.7 km long in western Procellarum at 8.5° N, 58.3° W today, its crest 68.7 m below the common level; its throat is 5.7 km wide with 2 m of water over the crest and 11.3 km wide with 7 m. Beyond it Cognitum, Humorum, Nubium and southern Procellarum, the south part, receive 6,500 m³/s from their own rivers and at the crest lose 19,800 to the air beyond their rain; the main part receives 88,200 and loses 77,100. The main part can spare 11,100 m³/s at the crest and the south part needs 13,300, so the sea parts. The main part stands 2.1 m over the crest, the neck pours 11,000 m³/s over it as a weir, and the south part falls about 100 m, to 166 m below the common level, keeping 829,000 of the 964,000 km² it covers at the crest. Each further 1,000 m³/s that the main part's rivers lose lowers the south part about 30 m more and the main part 0.1 m.

**The monthly tide** *(computed with a two-basin model)*. The equilibrium tide raises and lowers the south part's mean level against the main part's by up to 7.3 m in a month. With the neck wide open the south part would rise and fall 6.0 m a month, taking up to 11 million m³/s through the neck; on the atlas's surface, where the neck lay deep, the ocean models found monthly ranges of 6.8–8.7 m in Nubium, Humorum and Cognitum and currents of up to 1.6 m/s in Procellarum. Parted, the neck passes only the main part's spill, from nothing to 44,000 m³/s as the main part's tide lifts and lowers its water over the crest, at a mean speed through the throat of 1.3 m/s (2.4 at most), and the south part's own monthly range falls to 2 cm.

**Keeping one sea** *(computed)*. Two ways keep the parts joined:

| | Brings the near-side sea to | Takes |
|---|---|---|
| More water reaching the main part by gravity | 63.7 m down, the south part 0.5 m under the main one, the neck passing 13,500 m³/s with 5 m of water over its crest | 2,600 m³/s more than the main part now receives |
| A strait carrying the south part's need, 18 m of water over its bed, a drop of 0.5–1 m | 94.3 m down on both sides, 25.6 m under the crest: the main part 28 m lower than parted, the south part 72 m higher | a channel 1.4–2.0 km wide at its bed and 45 km long, 18.6 km of it through ground dry at that level, its cuts to 44 m: 1.1–1.6 km³ of digging |
| The same strait holding the level 10 or 20 m lower, for 600 or 1,400 m³/s of further losses | 104 or 114 m down | 2.7–3.8 or 3.9–5.4 km³ |

Through the strait the flow swings with the tide from 43,600 m³/s north to 38,600 south, at up to 1.1 m/s, and the south part's tide stays at 5 cm. A ship lane alone, 300 m wide and 10–20 m under the main part's level (10–33 million m³), leaves the parting as it is and runs at 3.7–4.2 m/s as part of the weir. These are first estimates: the bed drag is the ocean models' central 0.0025, the throat's section is taken over the whole sill's length, which overstates its friction, and the tide model leaves out the water's self-attraction (9% on the near side's seas) and the parts' own seiches.

## What is open

- **The near-side sea's parting** (register, 16 November) is due. On the closed budget of the run at 0.930 the sea parts at its neck, with a river pouring 11,000 m³/s over it and the south part 100 m lower; 2,600 m³/s more reaching the main part keeps it whole, and so does a strait of 1.1–1.6 km³ with the whole sea 94 m down. On the routing's own budget it rises to its rim. The design sunlight is open too, and a warmer climate moves these numbers.
- **The globe on the redrawn floor** is rebuilt once, after the design climate's run on the final surface (register, 16 November), and keeps its version 10 until then. Its builder (`visualization/atlas/render.py`) reads the atlas's water on the dry floor and geography's drainage; drawing the redrawn floor needs it to read the redrawn surface's seas and the governed routing's rivers and lakes. Its years control shows the ungoverned network's history on the dry floor, which the governed design no longer follows; what the rebuilt control shows is settled at the rebuild.
- **The design climate's sunlight on the final surface.** 0.930 settles at 300.3 K, and by the two runs 0.933 would give about 300.9 K and 0.934 about 301.1 K; whether the design temperature stays at 301 K or settles near 300 K, which cools the waters a little, is open, with the drylands above as one input. The drainage, the moisture map and the globe on the chosen climate follow.
- **Satoyama-like landscapes**: the criteria are in [satoyama.md](satoyama.md) and the first map on the redrawn floor's drainage is there too (above). Its candidates take a third of the hybrid fifth, nearly all of it in the equatorial rain belt, where they would fill 95% of that band's fifth; their flooded fields fit the methane room at its low-cost end and four times exceed it at its high-cost end. The salt sea has little shallow water for satoumi at 1.9 km cells.
- **The warm mantle's Orientale**: a bisection on Orientale's level would find the level at which its rim and its water agree, if that case comes to matter.

## Limits

- The water arrives at once in the response, where it fills over centuries; the settled state hardly depends on that, the first millennium's state a little more.
- The creep's Love numbers are incompressible and end at degree 50, beyond which the present Moon's thick lid barely flows; the elastic part is compressible and runs to degree 359.
- The lakes keep their depths as a load. The drainage rebuilt on the new floor gives them their own levels.
- The rotational and tidal bulge re-forming about the turned axes moves the sea surface by metres and is left out, as are the water the crust's pores take up and its uneven load, which in the crustal-loading study's relaxed case moves the straits by 64–78 m and the Tenebrarum sill by about 190 m, and the lateral changes in the crust and lithosphere.
- Every cell below the common level counts toward the 28%, as the atlas counts it, so hollows away from the seas are held at the common level until the drainage gives them their own.
- Nothing here is validated against a lunar observation; the Moon carries no load of this size.
- The drainage on the redrawn floor and its works ran on the present design climate on the turned Moon at 8°, whose rain fell on the atlas's waters. [final_drainage.py](final_drainage.py) routes the run on the final surface beside it; the products that read the drainage (the satoyama map, the rivers' losses, the globe) read the present climate's until the design sunlight is chosen. The design state holds the near-side ocean, the south sea and Humboldtianum at the common level whatever their balance, to measure what each lacks or gains; the governed state the round chooses replaces it.
- Cuts are measured between 1.9 km cells: a cut's depth is the cells' heights above the water and its length the path through them. Tunnels are not considered and no cut is priced.
- The first step loaded the crust with the hollows full to the common level; with most of them dry the load there is 1.25 million km³ lighter, a second-order change the next pass of the sea-level equation can take up.

## Running it

```
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.redrawn_surface    # about 25 minutes, 1.6 GB
```

It reads LOLA's heights and the geoid through `geography.fetch_inputs`, the geology lane's response table on the research drive (`geology/interior-response/interior_response_love.npz`) and the turned Moon's corrected drainage grid (`geography/turned_corrected/products/drainage_28pct_16ppd.npz`, ignored by Git, copied from main's checkout), and writes its grids to `atlas-round/redrawn_surface/` on the research drive, with a file of the cases finished so far after each case.

```
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.redrawn_drainage       # about 20 minutes, 4 GB
ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.satoyama_candidates     # under a minute
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.governed_drainage      # about 10 minutes
ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.settled_turn           # under a minute
ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.turned_surface         # under a minute
```

```
python -m climate.gcm.climatology A28S_301k_design_l15_snow_s930_settled:35-49   # then copied to ..._y35-49.npz
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.final_drainage         # about 10 minutes
ulimit -v 12000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.near_side_neck         # a few minutes
```

The turned surface writes its GCM files into `climate/gcm/products` (ignored by Git; in this worktree a link to main's checkout's) and its record into `results/`. The final surface's drainage carries the run's climatology onto today's ground there and reads it with the climate venv's Python, which has netCDF4. The run's scripts and rules are on the research drive in `atlas-round/step3/`: `prep.sh` (the base restart, one lunar day of a fresh start on the new surface, and the carry) and `run_design.sh` (the run under the design watch).

The drainage reads the first step's settled grids (checked against the hash its product records), the design climatology and geography's atlas, and keeps each piece on the research drive as it finishes (`atlas-round/redrawn_drainage/redrawn_drainage_so_far.json`); `--resume` takes the finished pieces from there and finds the geometry again, about 3 minutes. The design state's routing goes to `atlas-round/redrawn_drainage/design_routing.npz`, which the satoyama map reads and checks against the drainage's product.

## The checks

[test_redrawn_surface.py](test_redrawn_surface.py): the creep adds nothing when the water arrives and follows the response table between its times; it changes only the degrees the table holds; an upper sea stands at its rim over its own land while the rest of the water stands at the common level; and a rim is the lowest pass to any other sea, with every sea joined there named. On the product: the atlas case reproduces the drainage network study's rims to 0.1 m and geography's sills to 1 cm; every case holds 28% at its common level, the design case and the elastic one with their rims settled; and the elastic and design cases deepen the South Pole–Aitken Sea and its sill to the Tenebrarum sub-basin and need more water than the atlas.

[test_redrawn_drainage.py](test_redrawn_drainage.py), on a small walled basin above a sea with a hollow on the slope between: a full basin's surplus enters at the first cell past its pass, runs through the hollow and reaches the sea below, less only what the hollow's water loses to the air; a second, higher pass out of the basin is listed with the 10 m cut its level needs there, apart from the natural overflow; and on a graph of five catchments the route whose deepest cut is shallowest goes round through a 20 m cut where the direct way needs 30 m, its cuts adding up to the least. On the product: the redrawn floor reproduces the first step's rims, the full seas spill and what they send arrives below less the hollows' loss on the way, and the network's cases hold together. [test_governed_drainage.py](test_governed_drainage.py): along a route of cuts a branch is put in where it runs on by itself and taken out where a cut leaves a lake, and a cut onto a sea's water hands the sea its flow; a take-out at a lake the water passes leaves the sea below as it was and the river above it carries the flow; on a basin with a natural pass to one sea and a cut to another, the weir sends the first sea its need and the cut the rest. On the product: the held seas get their need, what the weirs release adds up to Moscoviense's surplus, the lakes on the way only take, and Smythii–Marginis passes on what reaches it. [test_settled_turn.py](test_settled_turn.py): the geology lane's own settled turn is reproduced from its load, the final water's masses add up and its frame is a rotation whose pole is the settled turn. [test_turned_surface.py](test_turned_surface.py): averaging 16 px/deg cells to 4 px/deg keeps the dry ground's height exactly as boundary.regrid derives it; on the product, the record carries the settled frame and its tilt, the water and the seas keep their shares through the turn, the files match their hashes, and the runner takes the record for experiment A28S. [test_satoyama_candidates.py](test_satoyama_candidates.py): a side valley ends where it reaches a large river, water or a closed hollow, and one too small is dropped; the paddy's methane rates are satoyama.md's; on the product, the candidates' flooded fields lie within their valley floors and 28% of their land, the bands add up and the methane follows from the area. [test_final_drainage.py](test_final_drainage.py): the run's climatology on today's ground is turned_climatology's carrying with the final surface's rotation, and the loader keeps the network study's runoff and open-water balance; on the product, the closing factor follows from the budget, the held seas get their need, the weirs release Moscoviense's surplus, and the near-side sea reaches its rim exactly when its balance is positive. [test_near_side_neck.py](test_near_side_neck.py): the level search agrees with the network study's tables, critical flow over a rectangular throat is B g^0.5 (2E/3)^1.5, friction alone sets a small head's flow, and the two-part balance conserves water; on the product, the throats widen with the level, the south part never stands over the main one, more loss lowers it, and the tide through a narrowed neck never exceeds the open neck's range.
