# Storm water on the drylands

1 November 2026, on `domain/hydrosphere`. PlaSim holds each land cell's water in one 0.5 m bucket that sheds water only when full, and on the design years the drylands keep that bucket 0.4–8% full, so the climate returns all of their rain to the air. The study works out what the drylands' storms would do instead. Phase 1 found how much of the rain runs off dry ground where it falls; phase 2a follows the wetted soil between storms and the runoff down the slopes and channels to the hollows, the sea and the soil. Both phases run offline on the design years `A28R_301k_design_l15_snow_s934`, years 10–19, following the plan of 1 November (`session-tools/agents/runoff/plan.md` on the research drive) and the register's direction of the same day ([the gap test's rules and the storm water](../../decisions.md#the-gap-tests-rules-and-the-storm-water-1-november-2026)): four designed surfaces side by side, young and mature crater hollows both worked out, no dam-building organisms, and the recharge reported beside the water table's 1% floor for dry land. [run.py](run.py) runs both phases in about 21 minutes and writes [results/storm_water.json](results/storm_water.json) (schema `terluna.research.storm-water-runoff/2`); no climate or cloud model runs.

## Phase 2a: the wetting between storms and the water delivered to the hollows

Phase 2a replaces phase 1's two treatments of the wetting between storms with a redistribution model at lunar gravity, routes each surface's runoff over geography's 1.9 km grid into its depression tree and the crater bowls below the grid, follows the water the slopes and channel beds take in below the root zone, and measures how much of the drylands' rain leaves the soil's evaporating store.

### What it finds

- **Between storms the wetted layer recovers most of the way.** With the layer spreading and drying between storms at random times within each mean, the drylands shed 6.2% of their rain from the sandy loam, 7.2% from the biological crust, 19.8% from the loam and 41.8% from the sealed loam under Walnut Gulch's storms *(model result)*:

  | Surface | Arid | Semi-arid | Dry sub-humid | Drylands (km³ a year) | Phase 1's bounds, semi-arid |
  |---|---|---|---|---|---|
  | Sandy loam | 4.3% (3.7–5.0) | 5.8% (4.8–7.1) | 7.3% (5.9–9.1) | 6.2% (5.1–7.6), 166 | 3.4–13.8% |
  | Biological crust on sandy loam | 5.1% (4.5–5.9) | 6.8% (5.8–8.1) | 8.4% (7.0–10.2) | 7.2% (6.1–8.6), 193 | 4.3–14.7% |
  | Loam | 14.7% (13.2–16.6) | 18.9% (16.6–21.9) | 22.7% (19.7–26.7) | 19.8% (17.3–23.0), 532 | 13.3–28.4% |
  | Sealed loam, seal formed | 34.0% (33.5–34.5) | 40.5% (39.9–41.4) | 46.2% (45.4–47.3) | 41.8% (41.1–42.7), 1,123 | 38.8–44.3% |
  | Loam sealing during each storm | 15.0% (13.6–16.9) | 19.3% (17.0–22.2) | 23.1% (20.1–27.0) | 20.1% (17.7–23.3), 542 | 13.8–28.7% |

  Each entry is the central estimate, with its range over the soil curves' spread in brackets. The central estimates lie in the lower third of phase 1's brackets: a storm's few millimetres spread downward within hours under capillary pull, which gravity leaves unchanged, and the surface dries, so the next storm meets soil that has made up much of its deficit. The semi-arid runoff rises by a factor of 1.7 over each storm on fresh soil for the sandy loam, 1.6 for the crust, 1.4 for the loam and 1.04 for the sealed loam, against phase 1's factors of 4.1, 3.4, 2.1 and 1.1 for one carried front. Under CM1's deeper storms the drylands shed 18.8%, 19.4%, 38.4% and 54.7%, with ranges of about a percentage point either way.

- **On cratered land most of the runoff ends in hollows.** Under Walnut Gulch's storms on young land, crater bowls below the grid take 48% of the sandy loam's runoff (19% falling inside them, 9% caught by their rims as the sheet runs down the slopes, 20% filling those the channels cross), the slopes soak back 42%, the grid's closed hollows receive 8%, the channel beds take 2.4% and under 0.1% reaches the sea or the humid land. The loam's larger storm runoff spends less of itself on the slopes and sends more on:

  | Surface (Walnut Gulch's storms) | Hollows: bowls below the grid | Hollows: the grid's depressions | Slopes' soil | Channel beds | Sea and humid land | To hollows, km³ a year |
  |---|---|---|---|---|---|---|
  | Sandy loam, young / mature | 48.0% / 38.0% | 8.0% / 12.9% | 41.6% / 45.7% | 2.4% / 2.7% | 0.1% / 0.6% | 93 / 84 |
  | Biological crust | 47.8% / 37.7% | 8.0% / 12.9% | 41.4% / 45.5% | 2.7% / 3.2% | 0.1% / 0.6% | 108 / 98 |
  | Loam | 59.7% / 49.1% | 11.4% / 18.2% | 25.7% / 28.3% | 3.0% / 3.5% | 0.2% / 0.9% | 378 / 358 |
  | Sealed loam | 69.9% / 58.4% | 14.6% / 23.3% | 10.1% / 11.1% | 5.1% / 6.0% | 0.2% / 1.1% | 949 / 917 |

  Shares are of each surface's runoff over all the drylands *(model result)*. Mature land, where bowls under 100 m have filled, holds fewer bowls (11% of the land inside them against 19%) and passes more to the grid's hollows, so the hollows receive almost the same water either way: 51–82% of the runoff on mature land against 56–85% on young. CM1's deeper storms soak far less into the slopes (2–9%) and send 87–96% of their runoff to hollows, 432–1,404 km³ a year. Of the drylands' 6.16 million 1.9 km cells, 1.59 million lie inside the grid's closed hollows already, which is why so little reaches a channel long enough to lose much to its bed: no flow path crosses more than 213 cells of slope before it ends. The water arrives at 96,793 of the grid's depressions, the median one receiving 71,000 m³ a year from the loam under Walnut Gulch's storms on young land (142,000 m³ on mature), with the largest hundred taking under 5% of the total.

- **Almost nothing drains below the root zone on the way.** The slopes' run-on soaks a few millimetres into soil already wetted by the storm, all of it within the top metre. The channel beds take 2–6% of the runoff, and only where a bed's yearly intake, a lunar day's share at a time, exceeds about 100 mm does any of it pass 1 m after a lunar day of redistribution in sand at lunar gravity (a tenth of 120 mm, a third of 200 mm) *(model result)*. The focused recharge from slopes and channels comes to under 0.001% of the rain on the sandy loam and the crust, 0.02–0.05% on the loam and 0.56–0.82% on the sealed loam under Walnut Gulch's storms, and 0.004–0.18% under CM1's, all of it below the water table's floor of 1% of the rain on dry land, which stays as it is. The water standing in the hollows recharges through their floors, which phase 2b works out.

- **After routing, the water leaving the soil's evaporating store passes the plan's thresholds on most surfaces.** The water reaching hollows, the sea and the humid land, with the recharge, equals 3.5% of the drylands' yearly evaporation from the sandy loam, 4.0% from the crust, 14.1% from the loam and 35.9% from the sealed loam on young land under Walnut Gulch's storms, and 8.6%, 10.1%, 36.6% and 98% of it in the wettest tenth, the third after noon *(model result)*; mature land differs by under a percentage point. Only the sandy loam stays below both of the plan's thresholds (5% of the yearly evaporation, 10% in any tenth), with the crust at the tenth's threshold; under CM1's storms every surface passes both, by 17–53% a year and 44–145% in the wettest tenth. This counts water moved from the slopes to the pools, not water lost to the air: the pools return much of it to the air in the days and weeks after the storms, later in the lunar day and from smaller areas *(estimate; phase 2b computes it)*. How much evaporation moves between tenths, and how much the pools' floors pass to the ground, is phase 2b's answer, and with it whether PlaSim's land scheme needs the change of the plan's option (b).

- **The chain at Earth's gravity matches Walnut Gulch's basin once its storms' runoff cores are about 5 km².** On a 129 km² basin of 1.9 km cells, Walnut Gulch's storms on Rawls' sandy loam and loam at the wilting point give 2.3% and 11.4% of a 305 mm year at the point, 1.3% and 7.9% off the slopes (the 0.18–5.42 ha watersheds measure 3–10%), and 1.2 and 11.8 mm a year at the outlet, a geometric mean of 3.8 mm against the measured 2 mm, inside the rule's factor of 3 *(model result against a measurement)*. The check sets the core area rather than testing it. The first version, which took each flow's duration where it formed, gave 9.6 mm; Lane's method takes the duration at the reach, which grows as the flood travels, and with that correction the chain gave 8.5 mm for cores of 50 km², failing the rule by a factor of 1.4. No measured core area was found, and the rule passes for cores of 3–9 km², so the value comes from the check.

### How it works

#### The wetting between storms

[wetting.py](wetting.py) carries the water each storm leaves in the topsoil as Smith, Corradini and Melone's (1993) rectangular profile, the form Ogden and Saghafian (1997) use between Green and Ampt storms: I mm of water above the soil's own content at a uniform saturation down to a sharp front. The front advances under gravity, at the layer's conductivity, and under the capillary pull of the drier soil below, Smith et al.'s profile factor p times the Brooks and Corey drive between the two contents over the front's depth. The surface dries at the smaller of the mean's Priestley and Taylor demand and the rate its desorptivity allows, S_d / (2√t), with S_d from the same curve (a sharp-front estimate). The next storm drives a new Green and Ampt front into the layer, with the deficit the layer has left and phase 1's front suction scaled by the drive from the layer's content to saturation; when it reaches the old front the two merge into one front carrying all the water into the dry soil. Each mean's storms fall at uniformly random times within its 2.95 days *(design choice)*, and the tables hold the expected runoff of exactly N storms for N up to 24 of Walnut Gulch's and 12 of CM1's, the Poisson mixture extrapolating beyond.

Gravity enters only through the conductivity (times g) and the bubbling pressure (over g). The capillary products, the diffusivity K dh/dθ, the drive, the desorptivity and with them the capillary spreading and the surface's drying are the same at any gravity; only the drainage under gravity falls to a sixth on the Moon *(calculation, tested)*. A wetted layer on the Moon therefore spreads almost as fast as on Earth while it is shallow and wet, and stays shallower once the capillary pull fades.

The curves are Rawls, Brakensiek and Saxton's (1982) class means as the GSSHA and SWMM5 tables reproduce them (read at that secondary): residual contents 0.041 and 0.027, bubbling pressures 146.6 and 111.5 mm and pore-size indices 0.378 and 0.252 (arithmetic means) for the sandy loam and the loam, and for the channels' sand 0.020, 72.6 mm and 0.694. Each curve's saturated conductivity is twice the Green and Ampt conductivity phase 1 uses, Rawls, Brakensiek and Miller's (1983) eq. 6 (read at the source), so that one soil underlies both models; the table's own 1.32 cm/h for a loam is twice that again. The range carries the bubbling pressure over one standard deviation of the class's suction either way (a factor of 4.1 for the sandy loam and 6.7 for the loam, Rawls et al. 1983, Table 2, read at the source), which sets both ends, and the profile factor at 2.0. Smith et al.'s factor of 1.7 is recalled, their paper not reached; with it the rectangular profile follows a numerical solution of Richards' equation for the same curves, at both gravities, to within 0.01–0.03 in surface saturation and a few millimetres in wetted depth over a day, with Richards' solution between the factors 1.7 and 2.0, while a factor of 1 misses by 0.04–0.05 *(calculation, tested)*. The pore-size index's spread was not found.

#### The routing

[routing.py](routing.py) and [phase2a.py](phase2a.py) route each surface's central runoff over geography's 1.9 km grid. `geography/drainage.py` now routes a supplied runoff field in place of the climatology's rain less evaporation (`compute(runoff=...)`, `--runoff`), keeping its default, with two new tests; the study supplies the loam's yearly runoff and takes the depression tree, the receivers and the drained areas back. Of the drylands' 6.16 million 1.9 km cells, 1.59 million lie below their depression's spill level, inside the grid's closed hollows; water arriving at any of them is delivered to that hollow. The other 4.57 million form a forest of slopes whose flow paths end in a hollow, the sea, or the humid land beyond the drylands.

Each storm's point runoff meets four things in turn.

- **The slope.** The excess runs 100 m down a slope of 0.05 *(design choices)* as a kinematic sheet under Manning's law with its gravity explicit, n = 0.05 on Earth *(an Earth value)* and 2.45 times that on the Moon: at the same discharge 1.71 times as deep and 0.58 times as fast *(calculation, tested)*. When the excess stops, the water still on the slope soaks into it at the soil's capacity while it drains over about half the equilibrium time; a slower, deeper sheet gives the slope longer to take it back.
- **The bowls below the grid.** The terrain study's saturation equilibrium, 1.54 D⁻² per km² with D in km at 0.69–3.9% of geometric saturation (its geometric middle, 1.6%), gives the bowls from 10 m to 1.9 km on young land and from 100 m on mature land, where the smaller ones have filled. They stand at half Pike's fresh depth of 0.196 D *(design choice: equilibrium craters are shallower)*. Water falling inside a bowl stays there; a sheet crossing one on its way down the slope fills it; and a channel crossing one loses up to the bowl's volume, the crossings Poisson along the path with a bowl's centre within half its diameter. Channels crossing bowls at random is the young land's picture; integrated drainage on mature land routes many channels around bowls, so the channels' capture is an upper bound there.
- **The channels.** Each reach runs to the steepest neighbour, 1.2 times the cells' distance *(design choice)*, in a bed w = 2.0 m (A / 1 km²)^0.5 wide *(design choice)*. Lane's (1983) method, NEH-4 chapter 19 (read at the source, eqs. 19-1 and 19-13 to 19-16), gives the reach's outflow a(x, w) + b(x, w) P above its threshold volume, with the unit channel's a = −0.00465 K D and k = −1.09 ln(1 − 0.00545 K D / P). The bed is Lane's class 2, clean sand and gravel at 3 in/h (2–5 in/h in his Table 19-1). The recast replaces K D, the depth a dry bed takes in during a flow, by K D times the ratio of the Green and Ampt depths the same dry bed takes in at the gravity in question and at Earth's (sand's front suction, 49.5 mm, Rawls et al. 1983): the capillary soaking stays Earth's and the drainage falls to a sixth, so a lunar bed takes in 0.72 of Earth's depth in a quarter-hour flow, 0.56 in an hour and 0.32 in ten hours *(calculation, tested)*. A flow's duration at a reach is its duration where it formed plus its travel from the farthest source at 1.5 m/s on Earth *(design choice)* times Manning's factor of 0.58 on the Moon.
- **The spread of a storm.** A storm makes its runoff over a core of 5 km²; where a reach drains more, one event brings it a core's worth and the reach sees proportionally more events. Syed et al. (2003, abstract) found that Walnut Gulch's storms mostly cover its 148 km² and that the core's coverage predicts the runoff; no core area was found, and the Walnut Gulch basin check below passes for cores of 3–9 km², so the check sets this value rather than testing it.

The storms enter as four classes of about equal runoff, drawn from each library's storms on the soil the design years give at storm time, and each cell's classes carry its central runoff in those shares. Every class falls on every dryland cell at once; the fates are summed per T21 cell where they happen and per depression where the water arrives, then weighted by each cell's events by tenth.

#### The recharge

The water the beds take in drains below the root zone, set at 1 m *(design choice)*, if the bed's yearly loss, spread over its area a lunar day's share at a time, reaches below that depth after a lunar day of redistribution and drying in sand at lunar gravity. Run-on soaks a few millimetres into the slopes, which the same calculation keeps within the root zone. The water standing in the hollows and seeping through their floors is phase 2b's.

### Evidence

| Value | Evidence state |
|---|---|
| Brooks and Corey's class curves, field capacity and wilting point | Measured on Earth (class means), read at a reliable secondary |
| Green and Ampt parameters, their suction's spread, Ks = 2 K | Measured on Earth (class means), read at the source |
| Smith et al.'s profile factor 1.7 | Recalled; checked against Richards' equation, which lies between 1.7 and 2.0 |
| The redistribution model's gravity scaling and its analytic limits | Calculation, tested |
| Lane's equations and bed classes | Regressions on Earth's measured channel losses, read at the source |
| Lane's recast for gravity | Calculation, tested |
| Walnut Gulch's per-flow channel losses (about 4 acre-ft a mile) | Measured (search summary of USDA-ARS reports) |
| Walnut Gulch's basin yield (about 2 mm a year at 149 km²) | Measured (Stone et al. 2008, abstract, via the plan) |
| The storms' runoff cores of 5 km² | Set by the Walnut Gulch basin check |
| Slope length and gradient, Manning's n, channel width, sinuosity, wave speed, root zone, bowl depth, storm timing within a mean | Design choices |
| The bowls' numbers | The terrain study's equilibrium law from measured lunar crater statistics (Gault 1970; Xiao and Werner 2015) |
| Runoff ratios, fates, recharge, leverage, hollows' inflow | Model result of this chain |

### Limits

- Every routing parameter is a design choice or an Earth value: the slopes' length, gradient and roughness, the channels' width, sinuosity and bed, the flood wave's speed and the root zone. The core area alone is constrained, by the Walnut Gulch basin, and Lane's coefficients are regressions on Earth's channels, recast for gravity by a calculation that no lunar measurement can test.
- The bowls catch channel flow wherever a channel crosses one at random. That is the young land's picture; on mature land integrated drainage routes many channels around bowls or through breached rims, so the channels' capture there (26–46% of the runoff under Walnut Gulch's storms) is an upper bound, and the grid's hollows would receive more of it. The bowls stand at half the fresh depth, and their capacity, 1.5–8.9 m of water spread over the whole land (3.7 m at the central count), never limits what one storm period brings. Their count is the terrain study's equilibrium share at its geometric middle; across its range of 0.69–3.9% of geometric saturation the land inside bowls spans 8.5–40% on young land and 4.8–24% on mature land, against 19% and 11% here, and the bowls' share of the runoff moves with it. The routing was not rerun at the range's ends.
- Each event class falls on every dryland cell at once, and each loss or delivery counts in the T21 cell where it happens. Flows crossing a T21 border put a little of one cell's runoff into its neighbour's account; over the drylands the fates sum to the runoff within 0.01%.
- The depression tree comes from `geography.drainage.compute` with its default climatology (run A) for the lakes' balance, where geography's committed drainage product uses the dimmer-shield climate. The heights, depressions, spill levels and receivers do not depend on it; the drained areas that set the channels' widths do, behind depressions that spill under that climate, and so does the lossless reference.
- Every hollow is terminal: water reaching one stays there in this phase. That holds while a hollow's inflow over a storm period stays far below its volume below the spill level, which phase 2b checks with its fill and spill.
- The wetting model keeps the soil of each mean from PlaSim's bucket and carries nothing from one mean to the next; the storms fall at random times within the mean's 2.95 days, against the demand averaged over the rainy means (3.1 mm a day). The large-scale rain falls on the mean's soil, as in phase 1's treatment of each storm on fresh soil.
- The recharge follows a rectangular profile in sand and in the slopes' soils, with roots drawing nothing below 1 m and the beds' water redistributing for one lunar day. Shallower roots or deeper, coarser alluvium would pass more; the floors of the hollows, where most of the water ends, are phase 2b's.

### What phase 2b needs

- **The hollows' inflow.** `research/runs/storm_water/hollows_2a.npz` (schema `terluna.research.storm-water-hollow-inflow/1`, rebuilt by run.py) holds, for each of the 96,794 depressions that receive water, its node in geography's depression tree, its floor cell, floor and spill heights and area below the spill level, and its inflow in m³ a year arriving in each tenth, per storm library, surface and hollow state. The product's cells give, per T21 cell and surface under Walnut Gulch's storms, the water reaching hollows (and the grid's share of it), the water lost to the soil and the part below the root zone, in mm a year by tenth, with CM1's storms in the class sums of `phase2a.routing.fates`.
- **The bowls' pools.** The bowls take most of the water. Their size distribution is in `phase2a.routing.bowls` and routing.py; each bowl's catchment is its own walls plus what its rim catches, a few times its pool's area, so most bowls hold water for days after a storm and act as seepage points, as the plan estimated.
- **The pools' budgets.** Penman's evaporation by day at water's albedo and the night's limit from the pool's heat and the GCM's air, from the radiation, air and wind by tenth that `extract.composite()` gives; seepage through clay-floored and sandy-banked floors at lunar gravity; fill and spill to the next depression; each hollow's area–volume curve from geography's grid.
- **The leverage, finished.** The evaporation the pools return by tenth, set against the slopes' evaporation the bucket gave, decides whether the drylands' evaporation moves beyond the plan's thresholds; the recharge through the hollows' floors joins the beds' recharge beside the 1% floor.
- **The biome maps' wet-patch layer**: the pools' wet-area fraction and persistence class per cell.

## Phase 1: the runoff from dry ground

Phase 1 works out how much of the rain would run off dry ground where it falls. [run.py](run.py) breaks each dryland cell's 3-day means into point storms, lets each storm soak into four designed surfaces by Green and Ampt's infiltration at lunar gravity, and writes the infiltration-excess runoff of every T21 dryland cell, by year and by tenth of the lunar day: a sandy loam, a loam, a sealed loam and a biological crust side by side over all three dryland classes. Its two treatments of the wetting between storms are the brackets phase 2a's redistribution model replaces; its other results stand.

### What it finds

- **The surface sets the runoff.** Under the Walnut Gulch storms, a family fitted to Earth's measured point storms, each storm falling on the soil the design years give at its 3-day mean, the drylands shed 3.5% of their rain from the sandy loam, 4.5% from the biological crust, 13.8% from the loam and 40% from the sealed loam *(model result)*. Over both storm libraries and both treatments of the wetting between storms the ranges are 3.5–21%, 4.5–22%, 14–41% and 40–56%. The arid land sheds a little less and the dry sub-humid land a little more, because its storms come more often and its soil holds more water:

  | Surface | Arid (4.47 million km², 70 mm a year) | Semi-arid (5.87 million km², 243 mm) | Dry sub-humid (2.20 million km², 434 mm) | Drylands (2,690 km³ a year) |
  |---|---|---|---|---|
  | Sandy loam | 2.8% (2.8–15.8) | 3.4% (3.4–20.4) | 4.0% (4.0–24.7) | 3.5% (3.5–21.4), 94–576 km³ |
  | Biological crust on sandy loam | 3.6% (3.6–16.3) | 4.3% (4.3–21.1) | 5.1% (5.1–25.4) | 4.5% (4.5–22.1), 122–594 km³ |
  | Loam | 11.1% (11.1–31.9) | 13.3% (13.3–39.4) | 15.5% (15.5–46.1) | 13.8% (13.8–40.9), 372–1,100 km³ |
  | Sealed loam, seal formed | 32.8% (32.8–45.4) | 38.8% (38.8–54.0) | 44.1% (44.1–61.5) | 40.0% (40.0–55.6), 1,075–1,497 km³ |
  | Loam sealing during each storm | 11.5% (11.5–32.8) | 13.8% (13.8–40.3) | 16.1% (16.1–47.1) | 14.3% (14.3–41.9), 386–1,126 km³ |

  Each entry is the yearly runoff over all the rain, area-weighted, first under Walnut Gulch's storms with each storm on the mean's soil (the configuration the Earth check tests), then the range over both libraries and both wetting treatments. Within the lunar day the runoff follows the convective rain: four tenths of the day, from shortly before noon to just after sunset, bring nine tenths of the semi-arid rain and of its runoff.

- **On the Moon a storm's depth counts as much as its peak.** Capillarity draws the first 20–30 mm of any storm into dry soil at nearly Earth's rate, because the conductivity falls to a sixth of Earth's while the suction head rises sixfold, and the sorptivity they set together stays the same. After that the capacity falls toward the conductivity itself, 1.8 mm/h in the sandy loam and 0.56 mm/h in the loam *(calculation)*. CM1's storms average 37.5 mm over 2.7 hours of rain, Walnut Gulch's 6.4 mm over 47 minutes, so the CM1 storms shed four times as much from the sandy loam as Walnut Gulch's do although their peaks are half as high. The 6 km window bounds the storms' intensities from below and their depths from above: its storms are the rain under each core's wettest column wherever that column moves, which a fixed point beneath the core receives only in part.

- **The wetting between storms is the largest open term.** Carrying one wetting front through all of a mean's storms multiplies the semi-arid runoff by 4.1 on the sandy loam, 3.4 on the biological crust, 2.1 on the loam and 1.1 on the sealed loam, against each storm meeting the mean's soil afresh *(model result)*. A semi-arid point takes about two of Walnut Gulch's storms in the four rainy tenths of each lunar day (0.36 of CM1's deeper ones), so its storms often follow each other within hours. Walnut Gulch's storms come days apart, and its runoff cannot separate the two treatments. How fast a wetted layer redistributes and dries between storms at lunar gravity decides where between the bounds the drylands lie; phase 2a's model puts them in the lower third.

- **Lunar gravity by itself raises the point runoff a little.** The same storms on the same soils at Earth's gravity shed 2.5% (sandy loam), 12.1% (loam), 38.4% (sealed loam) and 3.4% (biological crust) of the semi-arid rain, against 3.4%, 13.3%, 38.8% and 4.3% on the Moon *(model result)*. Short storms fall mostly within the capillary first hour, which gravity leaves nearly unchanged.

- **Lunar drops seal a bare loam slowly.** A raindrop's kinetic energy per millimetre of rain is 0.143 of Earth's in the drylands' air at 1.16 atm *(calculation)*. A loam whose seal forms afresh in each storm sheds 13.8% of the semi-arid rain, barely above the open loam's 13.3%; with Earth's drops it would shed 17.1%, and with the seal fully formed 38.8% *(model result)*. One e-folding of the seal's growth takes about 80 mm of lunar storm rain, against about 12 mm on Earth *(calculation from WEPP's stability factor)*, so an undisturbed lunar loam seals over roughly a year of semi-arid storms, while burrowing, roots and other disturbance can keep it open far more easily than on Earth.

- **At the point the runoff passes the plan's thresholds on most surfaces.** The runoff equals 3.5% (sandy loam) to 56% (sealed loam) of the drylands' yearly evaporation, and in the wettest tenth it reaches 9.5% to 150% of that tenth's evaporation. Only the sandy loam and the biological crust under Walnut Gulch's storms, each storm on fresh soil, stay near the plan's 5% a year and 10% in a tenth. This measures the water leaving the point. Run-on soaking into the slope below, the channels' beds and the pools return much of it to the soil and the air, and at Walnut Gulch the 149 km² basin yields under 1% of its rain against its small watersheds' 3–10% *(measured; abstract, via the plan)*. Whether the drylands' evaporation changes beyond the thresholds is therefore phase 2's answer, once the runoff is routed.

- **The chain passes its Earth check.** At Earth's gravity, Walnut Gulch's storms on Rawls' sandy loam and loam at the wilting point shed 2.3% and 11.4% of a 305 mm year whose winter third is too gentle to run off; the small watersheds measure 3–10%. The rule, set before the chain was run, asks the span of the two textures to overlap the measured range, which it does. The chain's geometric mean, 5.2%, lies within a factor of 1.06 of the measured 5.5%.

### How it works

#### The design years

[extract.py](extract.py) reads the design run `A28R_301k_design_l15_snow_s934`, years 10–19, the 1,200 3-day means of its year files on the research drive, and classes each land cell as `session-tools/design_run/design_years.py` does: its mean over the ten tenths of the lunar day, each tenth the mean of the 3-day means whose hour angle falls in it (from `climate.crm.gcm_site`'s subsolar longitude), and its rain over Priestley and Taylor's demand on that mean's net radiation, air and pressure, in UNEP's aridity classes. It finds the same shares: 19.1% of the land arid, 25.1% semi-arid and 9.4% dry sub-humid, 841 T21 cells together, with hyper-arid land at 4.2%. For each dryland cell it keeps every mean's rain and convective rain, evaporation, bucket water, net sunlight and longwave, air at 2 m, ground, surface pressure, 10 m wind and surface heat fluxes, with the tenth each mean falls in, in an ignored `research/runs/storm_water/design_years.npz` whose sidecar records the sha256 of each year file. `composite()` gives every field by tenth and over the day as `design_years.py` composites them; the pools of phase 2 read the radiation, air and wind from it.

#### The storms

[storms.py](storms.py) holds the storm libraries, each a set of point hyetographs at one-minute steps with weights, and the library is the input a later product replaces. A cell's 3-day mean of convective rain, D, brings a Poisson number of the library's storms to each point, with mean D over the library's mean storm depth; the expected runoff of the mean is the Poisson mixture of the runoff of that many storms.

- **Walnut Gulch** is a two-parameter family fitted to measured point storms. Each storm rises linearly to its peak over 5 minutes *(design choice)* and decays exponentially with a time constant τ, and storm depths are exponential with mean μ. Fitting the eight numbers of Osborn's (1983) Table 1, the seasonal counts and mean durations of point rain above 50, 100, 150 and 200 mm/h for a season of 203 mm of thunderstorm rain, gives μ = 6.4 mm and τ = 7.9 minutes. The fit gives 6.4, 1.3, 0.26 and 0.05 such events a season against the measured 4.6, 0.94, 0.29 and 0.06, lasting 6.7, 4.8, 4.1 and 3.6 minutes against 7.0, 4.3, 3.0 and 2.4. Held out of the fit, Mendez, Goodrich and Osborn's (2003) annual maxima test it: the family gives 115, 78, 47 and 24 mm/h over 5, 15, 30 and 60 minutes at the two-year return period against the measured 92, 63, 42 and 24, and 172, 116, 70 and 36 at ten years against 146, 95, 66 and 39. The family overstates the 5–30 minute intensities by 7–25% and matches the hour. It is discretised into 48 storms of equal probability.
- **CM1's 6 km window** reads any product of CM1's storm tracker (`terluna.climate.crm-storm-lives/1`) through `from_storm_lives()`. Each whole track of `box_0e_elec_corrected` (41 of its 68) becomes one storm: the rain under the core's wettest 6 km column at each ten-minute output. Its storms come from a box on wet ground at 0° N 0° E. When the semi-arid box at 44.3° S, 56.25° E runs and ten-minute windows are rerun from its restarts, `run.py --storms <its storm-lives product>` replaces this library.

The large-scale rain, a third of the semi-arid rain, falls after the mean's storms as one steady rain at 2 mm/h *(estimate)*. Rates of 1 and 4 mm/h move the semi-arid runoff by under 0.1 percentage points on every surface except the sealed loam, where 4 mm/h adds 2.8 points.

#### Infiltration on the four surfaces

[infiltration.py](infiltration.py) steps each storm in one-minute steps with Green and Ampt's capacity taken at each step's midpoint; rain first fills a 1 mm surface store *(design choice)*, the excess runs off, and the store soaks in after the rain. A crust is a thin layer of its own resistance over the soil, Green and Ampt's two-layer form, and the steady-rain test against Green and Ampt's exact solution holds within 0.5% at both gravities. Gravity enters through two terms. Every conductivity is k ρ g / μ and scales with g, every resistance with 1/g, and every suction head, a capillary pressure over ρ g, with 1/g; the product of conductivity and suction, and with it the sorptivity, stays the same at any gravity, while the time for gravity to take over from capillarity grows 36-fold *(calculation, each tested)*. The drops' energy per millimetre of rain scales with v², with g / ρ_air at a fixed drop size as the plan derives it, evaluated at the drylands' rain-weighted surface pressure; atmosphere/rain's solved fall speeds in the resort's air give 0.11–0.17 of Earth's energy for 2–4 mm drops, 0.137 at 3 mm against the scaling's 0.138 at the same pressure.

| Surface (design choice) | Soil and Earth parameters | On the Moon | Sources |
|---|---|---|---|
| Sandy loam | Conductivity 10.9 mm/h, front suction 110 mm, effective porosity 0.412 (porosity 0.453), field capacity 0.207, wilting point 0.095 | 1.81 mm/h, 665 mm | Rawls, Brakensiek and Miller (1983) class means, as Chow, Maidment and Mays (1988, Table 4.3.1) and HEC-HMS reproduce them; field capacity and wilting point from Rawls, Brakensiek and Saxton (1982), read in the GSSHA and SWMM5 tables that reproduce them |
| Loam | 3.4 mm/h, 88.9 mm, 0.434 (0.463), 0.270, 0.117 | 0.563 mm/h, 537 mm | As above |
| Sealed loam | The loam under a 5 mm seal of resistance 26.9 h (conductivity 0.19 mm/h) over a subsoil at 0.82 of its conductivity | 162 h (0.031 mm/h) over 0.46 mm/h | Rawls, Brakensiek, Simanton and Kohl (1990), Table 1 for a loam: a matric drop of 7 cm across the seal and a subcrust factor of 0.82, their eq. 4 for the seal's conductance; WEPP's 5 mm crust (Alberts et al. 1995) |
| Loam sealing during each storm | The seal's conductance relaxes from the open loam's to the formed seal's as exp(−C E) in the storm's cumulative rain energy E, with C = 0.0034 m²/J and Earth's rain carrying 29 (1 − 0.72 e^(−0.05 i)) J/m² per mm at i mm/h | E at 0.143 of Earth's | WEPP's crusting equation and its stability regression for a loam of 40% sand, 20% clay and 15 meq/100 g (Alberts et al. 1995, eqs. 7.9.3 and 7.9.8; fitted values 0.00006–0.0312 m²/J); Brown and Foster (1987), read in Nearing et al. (2017), which reproduces it |
| Biological crust | A 5 mm crust at half the sandy loam's conductivity over the sandy loam | 0.90 mm/h over 1.81 mm/h | Measured reductions of 39–77% in biocrusts' conductivity against bare soil (search summaries of recent field studies, among them J. Hydrol. 2019 and HESS 27:4247, 2023; the papers were not read); the ratio 0.5 is a design choice within them |

The soil at storm time reads the bucket's water as the water between the wilting point (empty) and field capacity (full) of each soil, since a Manabe bucket holds plant-available water. The semi-arid bucket's 13.6 mm puts the topsoil at 2.7% of that range, close to the wilting point. Reading it at PlaSim's wetness factor or leaving the topsoil air-dry moves the semi-arid runoff by a factor of 1.1–1.4.

#### The Earth checks and their sizing

The chain runs at Earth's gravity on Walnut Gulch-like inputs: two thirds of a 305 mm year as the Walnut Gulch storms on Rawls' sandy loam and loam at the wilting point, the winter's third at 2 mm/h, below both soils' conductivity. It gives 2.3% and 11.4% (air-dry, 2.0% and 9.7%; with the front carried through 3-day windows of an 80-day season, 3.8% and 15.4%), against the small watersheds' 3–10%. CM1's storms on the same soils at Earth's gravity give 6.6% and 26.9%.

The measured range spans a factor of 3.3, so the check resolves an error in the chain of about a factor of 1.8 either way. The lunar answer depends on each sub-model as follows (semi-arid, the factor between the lowest and highest runoff over each sub-model's alternatives, on the sandy loam, biological crust, loam and sealed loam):

| Sub-model | Alternatives | Factor | What tests it |
|---|---|---|---|
| Storm library | Walnut Gulch, CM1's 6 km window | 4.4, 3.6, 2.5, 1.4 | Osborn's counts (fitted) and Mendez's maxima (held out); nothing tests the lunar storms until the dry box runs |
| Wetting between storms | Each storm afresh, one front carried | 4.1, 3.4, 2.1, 1.1 | Nothing yet: Walnut Gulch's storms come days apart |
| Conductivity | Half and double | 6.8, 6.8, 2.6, 1.5 | The Walnut Gulch chain on two textures, within its resolution; the gravity factor analytically |
| Suction | Half and double | 6.3, 6.3, 2.5, 1.5 | As above |
| Surface store | 0 and 3 mm | 1.9, 1.9, 1.8, 1.7 | The Walnut Gulch chain carries the 1 mm store |
| Soil at storm time | Air-dry, wetness factor | 1.4, 1.3, 1.3, 1.1 | The Walnut Gulch chain at the wilting point and air-dry |
| Gravity | Earth's | 1.3, 1.3, 1.1, 1.0 | Analytic tests of each scaling |
| Large-scale rain | 1 and 4 mm/h | 1.0, 1.0, 1.0, 1.1 | Nothing (an estimate) |
| Seal formation | Formed, forming under lunar drops at 0.09–0.17 of Earth's energy, forming under Earth's | 2.8 on the sealed loam | The seal reproduces Rawls et al.'s crust factor exactly; the drops' energy against the rain solver |
| Crust conductivity | 0.35 and 0.61 of the soil's | 1.3 on the biological crust | The measured range of reductions |

The chain the Earth data test comes within a factor of 1.06 of the measurement, well inside the check's resolution. The two terms that carry the lunar answer by factors of 2–4, the storms and the wetting between them, are the two Walnut Gulch cannot test. The dry box's storms answer the first; a redistribution model at lunar gravity, or the pools' and soils' budgets of phase 2, answers the second. Mass conservation holds to a part in 10¹² on synthetic storms on every surface at both gravities.

### Evidence

| Value | Evidence state |
|---|---|
| The drylands' rain, convective rain, evaporation, bucket water and classes | Model result (PlaSim design years) |
| Walnut Gulch's storm counts, durations and annual maxima | Measured (Osborn 1983; Mendez et al. 2003) |
| The Walnut Gulch storm family | A fit to measured statistics |
| CM1's storms | Model result (CM1 on wet equatorial ground) |
| Soil texture parameters, the seal's matric drop and subcrust factor, WEPP's stability factor | Measured on Earth (class means, plot data and regressions) |
| Field capacity, wilting point, Earth's rain energy law | Measured on Earth; read at reliable secondaries in phase 2a, which confirmed the values phase 1 recalled |
| The biological crust's conductivity ratio | A design choice inside measured reductions read from search summaries |
| The four surfaces, the 1 mm store, the 5-minute rise | Design choices |
| The large-scale rain's 2 mm/h | Estimate |
| The gravity scalings and the drops' energy factor | Calculation, tested analytically and against the rain solver |
| Runoff ratios, storm counts, the share of storms that run off | Model result of this chain |
| Walnut Gulch's small-watershed runoff of 3–10% | Measured (Stone et al. 2008, from an abstract, via the plan) |

Values the build needed that it could not find in the time allowed: the conductivity of the smooth cyanobacterial crusts of hot deserts, which may seal more strongly than the crusts of the studies read; how densely a seal packs under drops with a seventh of Earth's energy, which here keeps Earth's formed seal; and Brooks and Corey's parameters with a redistribution model for the wetting between storms, which phase 2a supplies.

### Limits

- The runoff is the infiltration excess at a point on uniform ground. Run-on soaking into the slope below, the spread of storms over part of a cell and the channels' losses all lie in phase 2.
- The CM1 storms stand on wet ground at the equator and follow each core's wettest column; the Walnut Gulch storms are Earth's. Neither library is the drylands' own storms.
- The design years' numbers become final when the GCM and CM1 agree; `run.py --extract` rereads them.
- Water's viscosity and surface tension stay at the texture classes' reference temperature. Water at 30 °C, near the drylands' afternoon air, flows about a quarter more freely, which would raise every conductivity by that much on both worlds.
- The drops' energy factor holds at a fixed drop size. Lunar rain keeps larger drops, which survive breakup at lunar gravity, and at 4 mm the rain solver's speeds give 0.165 of Earth's energy per millimetre; the forming seal's sensitivity over 0.09–0.17 covers that range.

### What phase 1 left for phase 2

- **Routing.** A change to `geography/drainage.py` that lets it take a runoff field in place of the yearly rain less evaporation, with its own tests, then routing of each tenth's runoff on the 1.9 km grid into geography's depression tree, with overland flow at lunar gravity (1.7–1.8 times deeper and 0.55–0.58 times as fast as Earth's at the same discharge), seepage from runoff crossing the slope, and channel losses from Lane's method recast for a dry bed's capillary soaking and its drainage at a sixth of Earth's rate. Crater bowls below the grid act as local sinks from the terrain study's counts, for both young and mature hollows until the atlas round chooses its terrain, with no dam-building organisms yet.
- **The pools' budgets through the lunar night.** A water and energy balance for every closed depression by tenth: Penman's evaporation by day at water's albedo, the night's limit from the pool's heat and the GCM's air, seepage through clay-floored and sandy-banked floors, and fill and spill to the next depression; each hollow's share of lunar days wet and its depths.
- **Focused recharge** through channel beds and pool floors, reported beside the water table's 1% floor for dry land, which stays as it is.
- **The biome maps' wet-patch layer**: a sub-grid wet-area fraction and persistence class per cell.
- **What sustains or wears each surface**: its loss of soil to the storms' runoff against its formation (soil development, crust growth, the seal's re-forming over about 80 mm of lunar rain), and the salts left where pools dry.

## Running it and checks performed

```
/home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.storm_water.extract   # 20 s; needs netCDF4
/home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.storm_water.run       # both phases, about 21 minutes; the grid and one routing pass peak at 2.2 GB
python3 -m pytest research/studies/storm_water
```

[test_run.py](test_run.py) checks Green and Ampt's solution under steady rain at both gravities, each gravity scaling, the plan's first-hour uptake (30 mm on the Moon against 33–36 on Earth for a dry sandy loam, 15 against 16–17 for a loam), the seal against Rawls et al.'s crust factor, the drops' energy against the rain solver, mass conservation, the Walnut Gulch family against its measured statistics, the storm-tracker reader, the Poisson mixing, the reading of the bucket, the classes' agreement with the composites, and the product's binding to its code, inputs and constants. The literature was read for the values above, each a few minutes.

[test_phase2a.py](test_phase2a.py) checks that one storm on fresh soil in the wetting model is phase 1's storm exactly on every surface, mass conservation through storms and gaps at both gravities, the redistributed runoff between phase 1's two bounds, the constant-drive solution of the redistribution equation without gravity, the gravity-free limit and the capillary invariants, Green and Ampt's ponded intake against the sorptivity without gravity, the lunar layer staying shallower, the substeps' convergence, the rectangular profile against Richards' equation at both gravities (surface saturation within 0.03 and the depth holding nine tenths of the water within 10% over a day; the profile factor of 1.7 agrees to 0.01–0.03, a factor of 1 misses by 0.04–0.05), Manning's gravity factors, Lane's equations and their recast, a Walnut Gulch-like reach's loss of 2–5 acre-ft a mile, the bowls against the terrain study's law, the slopes' run-on bounds, the routing's water balance and topological order on a synthetic basin, the Walnut Gulch basin check, and the product's phase-2a fields; [geography/tests/test_drainage.py](../../../geography/tests/test_drainage.py) checks the supplied runoff field, and with the dimmer-shield climatology the changed `geography/drainage.py` reproduces geography's committed drainage product exactly (every summary field and grid array). Four products downstream of `geography/drainage.py` were rebuilt by the settle cascade (the lake kinds' census and product, lake methane and the resort study), each differing from before only in its hashes and, for three of them, the recorded run time. Two more pin the module's hash outside the cascade, the terrain study's routing and resources' sea phosphorus waters, and need their own reruns with geography's grid inputs: `make check` then passed 4,682 tests and failed those 2 binding tests. The literature for phase 2a was read for the values it needs, each a few minutes, and the sources say how; Smith, Corradini and Melone's paper and Ogden and Saghafian's were not reached, so their equation is used as recalled and checked against Richards' equation instead. No measured case of redistribution was found in the time allowed, and the scholarly reading of the sources remains to be done.
## Sources

- Alberts, E. E. et al. (1995). Chapter 7, Soil component. USDA Water Erosion Prediction Project, Hillslope Profile and Watershed Model Documentation, NSERL Report 10.
- Brown, L. C. and Foster, G. R. (1987). Storm erosivity using idealized intensity distributions. Transactions of the ASAE 30(2), 379–386 (read in Nearing et al. 2017).
- Celia, M. A., Bouloutas, E. T. and Zarba, R. L. (1990). A general mass-conservative numerical solution for the unsaturated flow equation. Water Resources Research 26(7), 1483–1496 (the method, recalled).
- Chow, V. T., Maidment, D. R. and Mays, L. W. (1988). Applied Hydrology, Table 4.3.1.
- GSSHA wiki, Table 9, as reproduced at swmm5.org (2015), Rawls and Brakensiek soil parameter estimates (read).
- Lane, L. J. (1983, reissued 2007). Chapter 19: Transmission losses. USDA NRCS National Engineering Handbook, Part 630 Hydrology (read: eqs. 19-1 to 19-16, Table 19-1).
- Mendez, A., Goodrich, D. C. and Osborn, H. B. (2003). Rainfall point intensities in an air mass thunderstorm environment: Walnut Gulch, Arizona. Journal of the American Water Resources Association 39(3), 611–621.
- Nearing, M. A. et al. (2017). Catena 157, 357–362 (read, eqs. 5–7, which reproduce Brown and Foster's equation and RUSLE2's change to it).
- Ogden, F. L. and Saghafian, B. (1997). Green and Ampt infiltration with redistribution. Journal of Irrigation and Drainage Engineering 123(5), 386–393 (not reached).
- Osborn, H. B. (1983). Timing and duration of high rainfall rates in the southwestern United States. Water Resources Research 19(4), 1036–1042.
- Rawls, W. J., Brakensiek, D. L. and Miller, N. (1983). Green-Ampt infiltration parameters from soils data. Journal of Hydraulic Engineering 109(1), 62–70 (Table 2 and eqs. 4–7 read in phase 2a).
- Rawls, W. J., Brakensiek, D. L. and Saxton, K. E. (1982). Estimation of soil water properties. Transactions of the ASAE 25(5), 1316–1320 (read in the GSSHA and SWMM5 tables, swmm5.org, 2015).
- Rawls, W. J., Brakensiek, D. L., Simanton, J. R. and Kohl, K. D. (1990). Development of a crust factor for a Green Ampt model. Transactions of the ASAE 33(4), 1224–1228.
- Smith, R. E., Corradini, C. and Melone, F. (1993). Modeling infiltration for multistorm runoff events. Water Resources Research 29(1), 133–144 (not reached; the redistribution equation and its profile factor recalled).
- Stone, J. J. et al. (2008). Long-term runoff database for a semiarid experimental watershed. Water Resources Research 44, W05S05 (abstract, via the plan).
- Syed, K. H. et al. (2003). Journal of Hydrology 271, 1–21, on the spatial structure of Walnut Gulch's thunderstorm rain and its runoff (abstract).
