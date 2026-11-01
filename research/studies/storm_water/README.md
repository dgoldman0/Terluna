# Storm water on the drylands, phase 1: the runoff from dry ground

1 November 2026, on `domain/hydrosphere`. PlaSim holds each land cell's water in one 0.5 m bucket that sheds water only when full, and on the design years the drylands keep that bucket 0.4–8% full, so the climate returns all of their rain to the air. This first phase works out how much of the rain would instead run off dry ground where it falls. [run.py](run.py) breaks each dryland cell's 3-day means into point storms, lets each storm soak into four designed surfaces by Green and Ampt's infiltration at lunar gravity, and writes the infiltration-excess runoff of every T21 dryland cell, by year and by tenth of the lunar day, to [results/storm_water.json](results/storm_water.json). The study follows the plan of 1 November (`session-tools/agents/runoff/plan.md` on the research drive) and the register's direction of the same day ([the gap test's rules and the storm water](../../decisions.md#the-gap-tests-rules-and-the-storm-water-1-november-2026)): a sandy loam, a loam, a sealed loam and a biological crust side by side over all three dryland classes, with the routing, the pools and the recharge left to phase 2. No climate or cloud model runs. The chain reruns in about four minutes when the design years or the storms change.

## What it finds

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

- **The wetting between storms is the largest open term.** Carrying one wetting front through all of a mean's storms multiplies the semi-arid runoff by 4.1 on the sandy loam, 3.4 on the biological crust, 2.1 on the loam and 1.1 on the sealed loam, against each storm meeting the mean's soil afresh *(model result)*. A semi-arid point takes about two of Walnut Gulch's storms in the four rainy tenths of each lunar day (0.36 of CM1's deeper ones), so its storms often follow each other within hours. Walnut Gulch's storms come days apart, and its runoff cannot separate the two treatments. How fast a wetted layer redistributes and dries between storms at lunar gravity decides where between the bounds the drylands lie.

- **Lunar gravity by itself raises the point runoff a little.** The same storms on the same soils at Earth's gravity shed 2.5% (sandy loam), 12.1% (loam), 38.4% (sealed loam) and 3.4% (biological crust) of the semi-arid rain, against 3.4%, 13.3%, 38.8% and 4.3% on the Moon *(model result)*. Short storms fall mostly within the capillary first hour, which gravity leaves nearly unchanged.

- **Lunar drops seal a bare loam slowly.** A raindrop's kinetic energy per millimetre of rain is 0.143 of Earth's in the drylands' air at 1.16 atm *(calculation)*. A loam whose seal forms afresh in each storm sheds 13.8% of the semi-arid rain, barely above the open loam's 13.3%; with Earth's drops it would shed 17.1%, and with the seal fully formed 38.8% *(model result)*. One e-folding of the seal's growth takes about 80 mm of lunar storm rain, against about 12 mm on Earth *(calculation from WEPP's stability factor)*, so an undisturbed lunar loam seals over roughly a year of semi-arid storms, while burrowing, roots and other disturbance can keep it open far more easily than on Earth.

- **At the point the runoff passes the plan's thresholds on most surfaces.** The runoff equals 3.5% (sandy loam) to 56% (sealed loam) of the drylands' yearly evaporation, and in the wettest tenth it reaches 9.5% to 150% of that tenth's evaporation. Only the sandy loam and the biological crust under Walnut Gulch's storms, each storm on fresh soil, stay near the plan's 5% a year and 10% in a tenth. This measures the water leaving the point. Run-on soaking into the slope below, the channels' beds and the pools return much of it to the soil and the air, and at Walnut Gulch the 149 km² basin yields under 1% of its rain against its small watersheds' 3–10% *(measured; abstract, via the plan)*. Whether the drylands' evaporation changes beyond the thresholds is therefore phase 2's answer, once the runoff is routed.

- **The chain passes its Earth check.** At Earth's gravity, Walnut Gulch's storms on Rawls' sandy loam and loam at the wilting point shed 2.3% and 11.4% of a 305 mm year whose winter third is too gentle to run off; the small watersheds measure 3–10%. The rule, set before the chain was run, asks the span of the two textures to overlap the measured range, which it does. The chain's geometric mean, 5.2%, lies within a factor of 1.06 of the measured 5.5%.

## How it works

### The design years

[extract.py](extract.py) reads the design run `A28R_301k_design_l15_snow_s934`, years 10–19, the 1,200 3-day means of its year files on the research drive, and classes each land cell as `session-tools/design_run/design_years.py` does: its mean over the ten tenths of the lunar day, each tenth the mean of the 3-day means whose hour angle falls in it (from `climate.crm.gcm_site`'s subsolar longitude), and its rain over Priestley and Taylor's demand on that mean's net radiation, air and pressure, in UNEP's aridity classes. It finds the same shares: 19.1% of the land arid, 25.1% semi-arid and 9.4% dry sub-humid, 841 T21 cells together, with hyper-arid land at 4.2%. For each dryland cell it keeps every mean's rain and convective rain, evaporation, bucket water, net sunlight and longwave, air at 2 m, ground, surface pressure, 10 m wind and surface heat fluxes, with the tenth each mean falls in, in an ignored `research/runs/storm_water/design_years.npz` whose sidecar records the sha256 of each year file. `composite()` gives every field by tenth and over the day as `design_years.py` composites them; the pools of phase 2 read the radiation, air and wind from it.

### The storms

[storms.py](storms.py) holds the storm libraries, each a set of point hyetographs at one-minute steps with weights, and the library is the input a later product replaces. A cell's 3-day mean of convective rain, D, brings a Poisson number of the library's storms to each point, with mean D over the library's mean storm depth; the expected runoff of the mean is the Poisson mixture of the runoff of that many storms.

- **Walnut Gulch** is a two-parameter family fitted to measured point storms. Each storm rises linearly to its peak over 5 minutes *(design choice)* and decays exponentially with a time constant τ, and storm depths are exponential with mean μ. Fitting the eight numbers of Osborn's (1983) Table 1, the seasonal counts and mean durations of point rain above 50, 100, 150 and 200 mm/h for a season of 203 mm of thunderstorm rain, gives μ = 6.4 mm and τ = 7.9 minutes. The fit gives 6.4, 1.3, 0.26 and 0.05 such events a season against the measured 4.6, 0.94, 0.29 and 0.06, lasting 6.7, 4.8, 4.1 and 3.6 minutes against 7.0, 4.3, 3.0 and 2.4. Held out of the fit, Mendez, Goodrich and Osborn's (2003) annual maxima test it: the family gives 115, 78, 47 and 24 mm/h over 5, 15, 30 and 60 minutes at the two-year return period against the measured 92, 63, 42 and 24, and 172, 116, 70 and 36 at ten years against 146, 95, 66 and 39. The family overstates the 5–30 minute intensities by 7–25% and matches the hour. It is discretised into 48 storms of equal probability.
- **CM1's 6 km window** reads any product of CM1's storm tracker (`terluna.climate.crm-storm-lives/1`) through `from_storm_lives()`. Each whole track of `box_0e_elec_corrected` (41 of its 68) becomes one storm: the rain under the core's wettest 6 km column at each ten-minute output. Its storms come from a box on wet ground at 0° N 0° E. When the semi-arid box at 44.3° S, 56.25° E runs and ten-minute windows are rerun from its restarts, `run.py --storms <its storm-lives product>` replaces this library.

The large-scale rain, a third of the semi-arid rain, falls after the mean's storms as one steady rain at 2 mm/h *(estimate)*. Rates of 1 and 4 mm/h move the semi-arid runoff by under 0.1 percentage points on every surface except the sealed loam, where 4 mm/h adds 2.8 points.

### Infiltration on the four surfaces

[infiltration.py](infiltration.py) steps each storm in one-minute steps with Green and Ampt's capacity taken at each step's midpoint; rain first fills a 1 mm surface store *(design choice)*, the excess runs off, and the store soaks in after the rain. A crust is a thin layer of its own resistance over the soil, Green and Ampt's two-layer form, and the steady-rain test against Green and Ampt's exact solution holds within 0.5% at both gravities. Gravity enters through two terms. Every conductivity is k ρ g / μ and scales with g, every resistance with 1/g, and every suction head, a capillary pressure over ρ g, with 1/g; the product of conductivity and suction, and with it the sorptivity, stays the same at any gravity, while the time for gravity to take over from capillarity grows 36-fold *(calculation, each tested)*. The drops' energy per millimetre of rain scales with v², with g / ρ_air at a fixed drop size as the plan derives it, evaluated at the drylands' rain-weighted surface pressure; atmosphere/rain's solved fall speeds in the resort's air give 0.11–0.17 of Earth's energy for 2–4 mm drops, 0.137 at 3 mm against the scaling's 0.138 at the same pressure.

| Surface (design choice) | Soil and Earth parameters | On the Moon | Sources |
|---|---|---|---|
| Sandy loam | Conductivity 10.9 mm/h, front suction 110 mm, effective porosity 0.412 (porosity 0.453), field capacity 0.207, wilting point 0.095 | 1.81 mm/h, 665 mm | Rawls, Brakensiek and Miller (1983) class means, as Chow, Maidment and Mays (1988, Table 4.3.1) and HEC-HMS reproduce them; field capacity and wilting point from Rawls, Brakensiek and Saxton (1982), recalled |
| Loam | 3.4 mm/h, 88.9 mm, 0.434 (0.463), 0.270, 0.117 | 0.563 mm/h, 537 mm | As above |
| Sealed loam | The loam under a 5 mm seal of resistance 26.9 h (conductivity 0.19 mm/h) over a subsoil at 0.82 of its conductivity | 162 h (0.031 mm/h) over 0.46 mm/h | Rawls, Brakensiek, Simanton and Kohl (1990), Table 1 for a loam: a matric drop of 7 cm across the seal and a subcrust factor of 0.82, their eq. 4 for the seal's conductance; WEPP's 5 mm crust (Alberts et al. 1995) |
| Loam sealing during each storm | The seal's conductance relaxes from the open loam's to the formed seal's as exp(−C E) in the storm's cumulative rain energy E, with C = 0.0034 m²/J and Earth's rain carrying 29 (1 − 0.72 e^(−0.05 i)) J/m² per mm at i mm/h | E at 0.143 of Earth's | WEPP's crusting equation and its stability regression for a loam of 40% sand, 20% clay and 15 meq/100 g (Alberts et al. 1995, eqs. 7.9.3 and 7.9.8; fitted values 0.00006–0.0312 m²/J); Brown and Foster (1987), recalled |
| Biological crust | A 5 mm crust at half the sandy loam's conductivity over the sandy loam | 0.90 mm/h over 1.81 mm/h | Measured reductions of 39–77% in biocrusts' conductivity against bare soil (search summaries of recent field studies, among them J. Hydrol. 2019 and HESS 27:4247, 2023; the papers were not read); the ratio 0.5 is a design choice within them |

The soil at storm time reads the bucket's water as the water between the wilting point (empty) and field capacity (full) of each soil, since a Manabe bucket holds plant-available water. The semi-arid bucket's 13.6 mm puts the topsoil at 2.7% of that range, close to the wilting point. Reading it at PlaSim's wetness factor or leaving the topsoil air-dry moves the semi-arid runoff by a factor of 1.1–1.4.

### The Earth checks and their sizing

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

## Evidence

| Value | Evidence state |
|---|---|
| The drylands' rain, convective rain, evaporation, bucket water and classes | Model result (PlaSim design years) |
| Walnut Gulch's storm counts, durations and annual maxima | Measured (Osborn 1983; Mendez et al. 2003) |
| The Walnut Gulch storm family | A fit to measured statistics |
| CM1's storms | Model result (CM1 on wet equatorial ground) |
| Soil texture parameters, the seal's matric drop and subcrust factor, WEPP's stability factor | Measured on Earth (class means, plot data and regressions) |
| Field capacity, wilting point, Earth's rain energy law | Measured, recalled without re-reading |
| The biological crust's conductivity ratio | A design choice inside measured reductions read from search summaries |
| The four surfaces, the 1 mm store, the 5-minute rise | Design choices |
| The large-scale rain's 2 mm/h | Estimate |
| The gravity scalings and the drops' energy factor | Calculation, tested analytically and against the rain solver |
| Runoff ratios, storm counts, the share of storms that run off | Model result of this chain |
| Walnut Gulch's small-watershed runoff of 3–10% | Measured (Stone et al. 2008, from an abstract, via the plan) |

Values the build needed that it could not find in the time allowed: the conductivity of the smooth cyanobacterial crusts of hot deserts, which may seal more strongly than the crusts of the studies read; how densely a seal packs under drops with a seventh of Earth's energy, which here keeps Earth's formed seal; and Brooks and Corey's parameters with a redistribution model for the wetting between storms.

## Limits

- The runoff is the infiltration excess at a point on uniform ground. Run-on soaking into the slope below, the spread of storms over part of a cell and the channels' losses all lie in phase 2.
- The CM1 storms stand on wet ground at the equator and follow each core's wettest column; the Walnut Gulch storms are Earth's. Neither library is the drylands' own storms.
- The design years' numbers become final when the GCM and CM1 agree; `run.py --extract` rereads them.
- Water's viscosity and surface tension stay at the texture classes' reference temperature. Water at 30 °C, near the drylands' afternoon air, flows about a quarter more freely, which would raise every conductivity by that much on both worlds.
- The drops' energy factor holds at a fixed drop size. Lunar rain keeps larger drops, which survive breakup at lunar gravity, and at 4 mm the rain solver's speeds give 0.165 of Earth's energy per millimetre; the forming seal's sensitivity over 0.09–0.17 covers that range.

## Phase 2

- **Routing.** A change to `geography/drainage.py` that lets it take a runoff field in place of the yearly rain less evaporation, with its own tests, then routing of each tenth's runoff on the 1.9 km grid into geography's depression tree, with overland flow at lunar gravity (1.7–1.8 times deeper and 0.55–0.58 times as fast as Earth's at the same discharge), seepage from runoff crossing the slope, and channel losses from Lane's method recast for a dry bed's capillary soaking and its drainage at a sixth of Earth's rate. Crater bowls below the grid act as local sinks from the terrain study's counts, for both young and mature hollows until the atlas round chooses its terrain, with no dam-building organisms yet.
- **The pools' budgets through the lunar night.** A water and energy balance for every closed depression by tenth: Penman's evaporation by day at water's albedo, the night's limit from the pool's heat and the GCM's air, seepage through clay-floored and sandy-banked floors, and fill and spill to the next depression; each hollow's share of lunar days wet and its depths.
- **Focused recharge** through channel beds and pool floors, reported beside the water table's 1% floor for dry land, which stays as it is.
- **The biome maps' wet-patch layer**: a sub-grid wet-area fraction and persistence class per cell.
- **What sustains or wears each surface**: its loss of soil to the storms' runoff against its formation (soil development, crust growth, the seal's re-forming over about 80 mm of lunar rain), and the salts left where pools dry.

## Running it and checks performed

```
/home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.storm_water.extract   # 20 s; needs netCDF4
/home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.storm_water.run       # about 4 minutes
python3 -m pytest research/studies/storm_water
```

[test_run.py](test_run.py) checks Green and Ampt's solution under steady rain at both gravities, each gravity scaling, the plan's first-hour uptake (30 mm on the Moon against 33–36 on Earth for a dry sandy loam, 15 against 16–17 for a loam), the seal against Rawls et al.'s crust factor, the drops' energy against the rain solver, mass conservation, the Walnut Gulch family against its measured statistics, the storm-tracker reader, the Poisson mixing, the reading of the bucket, the classes' agreement with the composites, and the product's binding to its code, inputs and constants. The literature was read for the values above, each a few minutes; the scholarly reading of the papers behind the recalled values remains to be done.

## Sources

- Alberts, E. E. et al. (1995). Chapter 7, Soil component. USDA Water Erosion Prediction Project, Hillslope Profile and Watershed Model Documentation, NSERL Report 10.
- Brown, L. C. and Foster, G. R. (1987). Storm erosivity using idealized intensity distributions. Transactions of the ASAE 30(2), 379–386 (recalled).
- Chow, V. T., Maidment, D. R. and Mays, L. W. (1988). Applied Hydrology, Table 4.3.1.
- Mendez, A., Goodrich, D. C. and Osborn, H. B. (2003). Rainfall point intensities in an air mass thunderstorm environment: Walnut Gulch, Arizona. Journal of the American Water Resources Association 39(3), 611–621.
- Osborn, H. B. (1983). Timing and duration of high rainfall rates in the southwestern United States. Water Resources Research 19(4), 1036–1042.
- Rawls, W. J., Brakensiek, D. L. and Miller, N. (1983). Green-Ampt infiltration parameters from soils data. Journal of Hydraulic Engineering 109(1), 62–70.
- Rawls, W. J., Brakensiek, D. L. and Saxton, K. E. (1982). Estimation of soil water properties. Transactions of the ASAE 25(5), 1316–1320 (recalled).
- Rawls, W. J., Brakensiek, D. L., Simanton, J. R. and Kohl, K. D. (1990). Development of a crust factor for a Green Ampt model. Transactions of the ASAE 33(4), 1224–1228.
- Stone, J. J. et al. (2008). Long-term runoff database for a semiarid experimental watershed. Water Resources Research 44, W05S05 (abstract, via the plan).
