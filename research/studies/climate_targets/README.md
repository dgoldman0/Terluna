# What each target climate takes

On 17 October the author opened the target climate as a question: whether the Moon is brought back to the design
temperature, with more dimming tested, or kept hotter and wetter provided it stays comfortable, with cooler and wintry
weather aloft, long coasts and perhaps a life chiefly of recreation
([register](../../decisions.md#the-target-climate-17-october-2026)). New runs wait on the choice. This study sets
three targets side by side across the domains, at the joint integration of the five branches
([record](../../integration/five-branches.md)): the hot world, where the corrected GCM settles; the lock, the design
temperature of 295 K; and between them the resort, a tropical resort Moon, which the study sets at 300 K. With them it
sets the heat the plan adds to the Moon beyond its sunlight, which every target answers, the winter each target leaves
on the summit tower, and the resort's best value near 300 K. [run.py](run.py) is arithmetic on the GCM's output and
the domains' committed products ([results/climate_targets.json](results/climate_targets.json)); no model runs.
[The three target climates compared](comparison.md) sets them out for the author's choice, with each one's case, plant
growth, the lakes and the levers that ease the shield.

## Where the corrected Moon settles

The design case (`A28_corrected_970`) holds 295.0 K at the design sunlight scale of 0.970 on PlaSim's earlier
longwave. With the longwave refitted by height and ice optics for the cold cloud, the 20-year test at the same sunlight
(`A28_corrected_970_lw_ice`) averages 306.5 K over years 10–19
([climate/gcm](../../../climate/gcm/README.md#where-the-corrected-moon-settles)). Its top of the air gains
1.65–2.16 W/m² less for each kelvin it warms (fitted over years 5–19 and 0–19), so it settles at 306.8–307.1 K, and that
restoring gain holds it. The radiative–convective model's Manabe–Wetherald columns, whose humidity falls off with
height as Earth's does, show no runaway up to 320 K
([radiative–convective column](../../../atmosphere/radiative_convective/README.md)). The hot world absorbs 225.4 W/m²
of sunlight, 236.2 W/m² under its clear sky. A fully saturated clear column there can send out at most 239.0 W/m², so
a saturated, cloudless hot world would only just balance; its unsaturated clear sky sends out 241.5 W/m². Its cloud
reflects 10.8 W/m² of sunlight and holds back 16.8 W/m² of infrared, 6.0 W/m² of net warming. Its balance rests on its
air staying unsaturated.

## The heat beyond the sunlight

| What the plan adds | As absorbed sunlight | Share of the absorbed sunlight | Basis |
|---|---|---|---|
| Methane and N2O at 1.9 and 0.1 ppm | 3.58 W/m² | 1.6% | The deep column's all-sky forcing of the held pair; the GCM held it with dimming of the same watts |
| The fleet's glow, both halves rolled | 3.89–4.24 W/m² | 1.7–1.9% | 5.34 W/m² of infrared at the GCM's efficacy of 0.645–0.71, and 0.61 W/m² of sunlight sent back |
| Beamed power used on the Moon | 0.87–4.14 W/m² | 0.4–1.8% | The power screen's 30–145 TW of beams arriving at the receivers and the towns' crowns, and the deep air's 2.6–12 TW of them |
| **Together** | **8.3–12.0 W/m²** | **3.7–5.3%** | |
| The water top-up's fall, the first thousand years | 1.5–9.5 W/m² | 0.7–4.2% | 1,500–4,100 km³ a year at the central values, falling 1.2–2.7 MJ/kg, all of it spent on the Moon |

The shares are of the hot world's absorbed sunlight, 225.4 W/m²; the glow study and the register give the glow's
dimming as 2.0% of the earlier design case's sunlight at the top of the air, 293.3 W/m² with 211.2 absorbed. The
gases, the glow and the power last as long as the plan does, and the plan answers them with dimming, so they come
on top of whatever a target needs. The top-up is the water the author keeps delivering to balance what the thawing
crust takes up ([geography](../../../geography/README.md#groundwater-how-fast-the-crust-takes-the-water)). Its fall
into the Moon's own well releases 1.8–4.0 W/m² through the years to 10,000 at the central values, and 12–28 W/m² in
the fast case's busiest centuries. Resources' first screen already asks the build's deliveries for energy recovery on
the way down: tethers, elevators or regenerative lowering
([resources](../../../resources/README.md)); the top-up needs the same for millennia after the build.

**The glow's efficacy runs from 0.645 to 0.71.** The 0.71 is the held branch of 13 October on the earlier radiation
([gases, glow and dimming](../gases_glow_dimming/README.md#what-it-finds)). The 0.645 is the glow pair of 17 October
on the corrected radiation, the design case's own before the longwave refit
([climate/gcm](../../../climate/gcm/README.md#the-gases-the-glow-and-the-dimming)), found the same way from the two
branches' yearly records, which the extract binds by their hashes. Over years 10–19 the branch with 13 October's
4.22% of added dimming stands at 293.45 K and the branch with its 3.52% at 294.41 K: 1.37 K for each point of dimming.
Carried to the design case's 295.02 K, the gases and the glow together take 3.07% of its sunlight, 7.03 W/m² of its
228.5 W/m² absorbed. The pair fixes only that total, since no gases-only branch flew on this radiation. The gases take
their forcing, 3.58 W/m², as the table books them and as 13 October's 1.69% was their forcing over that radiation's
211.2 W/m², and the glow takes the rest, 3.45 W/m², 0.645 of its 5.34 W/m² of infrared. The other split keeps the
gases at their 1.69% of this radiation's larger absorbed sunlight, 3.87 W/m², and leaves the glow 0.59, with the same
totals. The branches still drift over years 10–19, the one with more dimming cooling 0.11 K a decade, and the pair's
rerun on the final radiation waits on the target. At 0.645 the glow's term falls from 4.24 to 3.89 W/m², which takes
0.16 points off the dimming every target shares.

## The three targets

Each target's dimming is the test's gain per kelvin times the kelvins it sheds, as a share of the sunlight reaching
the Moon at the test's albedo. Its climate is interpolated between the two runs by the ground's temperature, and a
target beyond a run takes that run's own climate, so the 295 K column is the design case's.

| | The hot world: no further dimming | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Ground | 306.5 K, settling at 306.8–307.1 | 300 K | 295 K |
| Absorbed sunlight to shed | none | 11.8–14.7 W/m² | 20.0–25.5 W/m² |
| Its dimming | none | 5.2–6.5% | 8.9–11.3% |
| With the heat beyond the sunlight | 3.7–5.3% | 8.9–11.8% | 12.6–16.6% |
| Sunlight scale of a design case for the target, from 0.970 | 0.970 | 0.907–0.919 | 0.860–0.884 |
| Sunlight scale with the added heat answered too | 0.919–0.934 | 0.855–0.884 | 0.809–0.848 |
| Sunlight the target's own dimming takes; electricity at 20–30% | none | 618–773 TW; 124–232 TW | 1,052–1,340 TW; 210–402 TW |
| Growth lost through light to the target's own dimming | none | 0.37–1.9% | 0.62–3.3% |
| Rain on land, a day, every cell's land (the model's land cells) | 3.34 mm (2.94) | 2.72 mm (2.39); 2.48–2.64 with the dimming's own cut | 2.25 mm (1.98); 1.90–2.14 |
| Land under 0.5 mm a day, the same | 18.9% (20.1%) | 29% (33%) | 38% (42%) |
| Land wet-bulb at 2 m: median; 99th percentile | 26.9; 31.5 °C | 21.1; 26.0 °C | 16.6; 21.8 °C |
| Land wet-bulb 0.9 km up: median; 99th percentile | 26.3; 30.9 °C | 20.4; 25.3 °C | 15.9; 21.1 °C |
| Land air, median: at 2 m; 0.9 km up | 33.6; 32.1 °C | 26.1; 24.7 °C | 20.4; 19.0 °C |
| Seas' surface, median | 34.5 °C | 28.6 °C | 24.1 °C |
| Where the land's mean profile reaches 20 °C | 12.3 km up | 6.7 km up | below the lowest level, 0.9 km up |
| Where the land's mean profile freezes | 35.9 km up | 28.1 km up | 22.8 km up |
| Clear-sky absorbed sunlight, and its gap to a saturated column's 239.0 W/m² | 236.2; 2.8 W/m² | 220.8–223.8; 15.1–18.2 W/m² | 209.4–215.2; 23.8–29.5 W/m² |
| That gap with the added heat at full watts: answered; unanswered | 1.3–1.8; −10.7 to −7.4 W/m² | 13.6–17.3; 1.6–8.0 W/m² | 22.3–28.6; 10.3–19.3 W/m² |

The heat is in the 3-day means, each about a tenth of the lunar day, at 2 m and at the lowest model level, about
0.9 km up, on the model's own land: it runs 166 cells with a land fraction of 0.50–0.59 as water, and the land's heat
is weighted by land area and the seas' by water area. The output's air at 2 m is the lowest level's carried down the
dry adiabat to the ground's pressure, within 0.02–0.04 K on average over land, and the wet-bulb temperature at 2 m
carries the lowest level's specific humidity down with it, capped at saturation: together they read a well-mixed
lowest kilometre, as by day. The output's humidity is a snapshot: ExoPlaSim 3.4.2 computes it only at the steps where
mod(nstep, nafter) is zero (plasim.f90, line 3218) and adds that held value at every step of the mean (outmod.f90,
line 2556), so each record holds the humidity at its mean's first step, while the air, the pressure and the Sun's
cosine are true means. Each mean's humidity is therefore the mean of its start and its end, the next record's; each
run's last mean, with no next record, leaves the heat's statistics, which cover 1,199 of the 1,200. The wet-bulb
temperature is capped at the air's own temperature where the air is supersaturated, which with the humidity as
written held in 1.7% of the design case's land means at the lowest level and now holds in 0.24%. The land's mean
profile is the land-area mean of its columns' air at each height above the ground; the median of the columns' own
freezing heights is 36 km in the hot world and 23 km in the lock. In the hot world the wet-bulb temperature at 2 m
passes 28 °C for 31% of the land's time: 49% within 15° of the equator, 56% on land below 1 km and 4.6% on ground
above 4 km. It passes 31 °C for 2.1%. At 0.9 km the same shares are 22.3%, 38%, 44% and 1.7%, with 0.8% above 31 °C.
The first readout of 17 October, which took the humidity as written and counted every cell with a land fraction above
one half alike, put the lowest level's time above 28 °C at 17.6%.

The rain falls on a cell's land as well as on its water, so the land's rain and deserts count every cell's land
fraction and area with that cell's own rain; the model's own land cells alone give the dry bound, in brackets. The
18% of the boundary file's land that lies in cells the model runs as water is mostly coast and lakeshore, and those
cells rain 5.15 mm a day in the hot world. The T21 grid, about 170 km a cell, cannot resolve the coasts, and coast
and island hydrology, geography's groundwater work and the shore study to come, decides where between the two the
real coasts lie. The rain's own cut takes the GCM's 0.55% of the rain for each 1% of dimming at a held temperature
(the gases, glow and dimming study) and Earth's models' 1.38% (the light study). The costs in the table are each
target's own dimming. The 3.7–5.3% every target shares for the added heat costs a further 0.26–1.56% of the growth
through light and 2.0–7.3% of the land's rain, which brings the hot world's to 3.10–3.27 mm a day (2.72–2.88). The
gap to the saturated column counts the added
heat at full watts, the glow's 5.34 W/m² of infrared among them, where the dimming answers the glow at the GCM's
efficacy; its ends take the ends of the feedback, the glow's efficacy and the beamed power that leave the clear sky the
most load and the least.

## The dimming at every temperature

Each kelvin the Moon is held below where it settles takes 0.73–0.96% of the sunlight beyond the plan's additions: each
fit's own feedback, 1.65 and 2.16 W/m² a kelvin, over the test's 225.4 W/m² of absorbed sunlight, counted from that
fit's own balance. The sunlight scales are those a design case for the target would fly from 0.970, before and after
the heat beyond the sunlight is answered too; the product carries every half kelvin from 295 to 307 K.

| Ground | Its dimming | Sunlight scale of its design case | With the heat beyond the sunlight | Sunlight scale with it answered |
|---|---|---|---|---|
| 295 K | 8.9–11.3% | 0.860–0.884 | 12.6–16.6% | 0.809–0.848 |
| 296 K | 8.2–10.4% | 0.870–0.891 | 11.8–15.7% | 0.818–0.855 |
| 297 K | 7.4–9.4% | 0.879–0.898 | 11.1–14.7% | 0.827–0.862 |
| 298 K | 6.7–8.4% | 0.888–0.905 | 10.4–13.7% | 0.837–0.869 |
| 299 K | 6.0–7.5% | 0.897–0.912 | 9.6–12.8% | 0.846–0.876 |
| 300 K | 5.2–6.5% | 0.907–0.919 | 8.9–11.8% | 0.855–0.884 |
| 301 K | 4.5–5.6% | 0.916–0.927 | 8.2–10.9% | 0.865–0.891 |
| 302 K | 3.8–4.6% | 0.925–0.934 | 7.4–9.9% | 0.874–0.898 |
| 303 K | 3.0–3.7% | 0.935–0.941 | 6.7–9.0% | 0.883–0.905 |
| 304 K | 2.3–2.7% | 0.944–0.948 | 6.0–8.0% | 0.892–0.912 |
| 305 K | 1.6–1.7% | 0.953–0.955 | 5.2–7.0% | 0.902–0.919 |
| 306 K | 0.8% | 0.962 | 4.5–6.1% | 0.911–0.926 |
| 306.8–307.1 K, where it settles | none | 0.970 | 3.7–5.3% | 0.919–0.934 |

The two fits cross near 305.9 K, so the fit of years 5–19 gives the high end below it and the low end above it.

## Day and night

The output's czen, the cosine of the Sun's zenith angle averaged over each 3-day mean, tells the day's means from the
night's. A mean is the day's when the Sun is up for more than half of it. Each mean holds 143 of the model's half-hour
steps, and at the Sun's synodic rate, 29.53 days to a lunar day against the 29.79 days of the model's 1,430 steps,
those take the Sun through 36.3° of its hour angle. A mean that long centred on sunrise or sunset has a czen of 0.0786
times the cosine of the latitude, and the day's means are those above it. The cut is exact where the Sun's declination
is zero; the fixed orbit swings the declination through the tilt, 1.54° either way. The day and the night then take
half of the land's time each. In the land's 3-day means, the day's figure first and the night's second:

| Land: the day's means; the night's | The hot world | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Wet-bulb at 2 m: median | 27.3; 26.6 °C | 21.7; 20.6 °C | 17.4; 16.0 °C |
| Wet-bulb at 2 m: 99th percentile | 31.7; 30.9 °C | 26.3; 25.0 °C | 22.1; 20.5 °C |
| Air at 2 m: median | 34.7; 32.9 °C | 27.7; 25.3 °C | 22.4; 19.5 °C |
| Air at 2 m: 99th percentile | 39.7; 35.9 °C | 33.2; 28.8 °C | 28.3; 23.4 °C |
| Wet-bulb 0.9 km up: median | 26.7; 26.0 °C | 21.0; 19.9 °C | 16.6; 15.3 °C |
| Wet-bulb 0.9 km up: 99th percentile | 31.2; 30.4 °C | 25.7; 24.4 °C | 21.5; 19.8 °C |
| Air 0.9 km up: median | 33.2; 31.4 °C | 26.2; 23.9 °C | 20.9; 18.1 °C |
| Air 0.9 km up: 99th percentile | 37.9; 34.5 °C | 31.4; 27.4 °C | 26.5; 21.9 °C |

The night eases the air more than the humid heat. In the hot world the night's air at 2 m is 1.8 K cooler at the
median and 3.7 K at the 99th percentile, while its wet-bulb temperature falls 0.7 K at the median and from 31.7 to
30.9 °C at the 99th percentile. At a fixed humidity a kelvin of air moves the wet-bulb temperature 0.26–0.27 K at the
hot world's day medians and 0.37–0.39 K at the lock's, over the land's surface pressures. At 2 m the hot world's
wet-bulb temperature passes 28 °C for 38.7% of the day's land time and 23.8% of the night's, and 31 °C for 3.3% and
0.9%; at 0.9 km, 29.0% and 15.7%, and 1.4% and 0.3%. In the lock the night takes 1.4 K off the wet-bulb median at 2 m
and 2.9 K off the air's. The night's figures at 2 m read warm: the 2 m air is the lowest level's carried down, and in
calm night air the air near the ground cools below it, as the cloud-resolving ring's calm, saturated land nights show
([the CRM's account](../../../climate/crm/README.md#what-the-ring-shows)). In each place's coolest tenth of its means
by the air, 81% of them nights in the hot world, the 2 m air stands 1.1 K above the ground's own temperature, and
1.8 K in the lock. The resort is interpolated between the two runs, as in the targets' table.

## The weather shorter than three days

`A28_corrected_970_eddies` carried the design case on from its year-19 restart for two years, its years 0 and 1,
writing every 6.5 hours: eleven means to each 3-day mean, 1,320 a year. It saved the output every run saves, the air's
temperature, humidity and winds at every level, the surface pressure, the 2 m air, the rain, the fluxes and the Sun's
cosine. Each 6.5-hour mean takes its humidity from its start and its end, as the 3-day output does, and the last
record, with no next one, leaves the 239 whole 3-day means and the 2,629 6.5-hour means within them. Over those
means the land's heat in the 6.5-hour means stands against the 3-day means their elevens make, with each extreme's
excess and the ratio of its height above the median:

| Design case, land, its two 6.5-hour years | Median: 6.5 hours; 3-day means | 99th percentile: the same (excess; ratio) | Top: the same (excess; ratio) |
|---|---|---|---|
| Wet-bulb at 2 m | 16.61; 16.58 °C | 22.1; 21.9 °C (+0.29 K; 1.049) | 25.0; 24.1 °C (+0.97 K; 1.126) |
| Air at 2 m | 20.40; 20.39 °C | 28.0; 27.7 °C (+0.22 K; 1.029) | 30.8; 30.1 °C (+0.66 K; 1.067) |
| Wet-bulb 0.9 km up | 15.88; 15.85 °C | 21.5; 21.2 °C (+0.28 K; 1.048) | 24.3; 23.4 °C (+0.91 K; 1.117) |
| Air 0.9 km up | 18.97; 18.97 °C | 26.2; 26.0 °C (+0.21 K; 1.029) | 28.8; 28.3 °C (+0.57 K; 1.060) |

The 3-day means smooth the land's 99th percentile by 0.2–0.3 K and its top by 0.6–1.0 K. The 3-day means made from the
two years match the design case's own 3-day output over its years 10–19 within 0.1 K at the median and the 99th
percentile (the wet-bulb temperature at 2 m: 16.58 and 21.85 °C against 16.59 and 21.79 °C). In the 6.5-hour means
1.6% of the land's lowest-level means are supersaturated, against 0.21% of the 3-day means, and the wet-bulb
temperature takes the air's own there. The 6.5-hour means are themselves means, in cells about 170 km across. The
ratios are the design case's alone; stretching each target's 99th percentile above its median by them gives, as a
sensitivity, the 99th percentiles in 6.5-hour means, with the 3-day means' in brackets:

| Land, 99th percentile | The hot world | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Wet-bulb at 2 m | 31.7 (31.5) °C | 26.2 (26.0) °C | 22.0 (21.8) °C |
| Air at 2 m | 39.4 (39.3) °C | 33.0 (32.8) °C | 28.0 (27.8) °C |
| Wet-bulb 0.9 km up | 31.1 (30.9) °C | 25.6 (25.3) °C | 21.4 (21.1) °C |
| Air 0.9 km up | 37.7 (37.6) °C | 31.2 (31.0) °C | 26.3 (26.1) °C |

## Winter on the tower

The summit tower stands on the Moon's highest ground, at 5.4° N, 158.6° W, 11.6 km above sea level, and rises 24 km to
35.6 km ([summit tower](../summit_tower/README.md#the-central-port)). Its disks begin 12.3 km above the summit, the
second, at 13.8 km, holds the open winter parks, the blocks are sealed from the third, at 15.3 km, and the crown's
terminal stands at 23.4 km ([register](../../decisions.md#the-summit-port-and-metropolis)).

**How the tower meets the model.** T21 smooths the summit away. The cell holding it, centred at 2.8° N, 157.5° W, is
42% land, its land averaging 9.1 km above sea level in the boundary file, and the runner runs it as lake water at its
lakes' level. The model's own ground is 8.0 km there and 7.6 km at the summit's place, between the four cell
centres round it. The output carries no orography, so the extract takes the ground from the sea-level pressure, which
pyburn wrote by carrying the lowest level's air down from that ground at the standard atmosphere's 6.5 K/km; inverted
the same way, every mean gives the same ground to 0.01 m. Each column's levels rise hypsometrically from its ground
with the virtual temperature, heights take inverse-square gravity, and the tower is mapped onto the column by height
above sea level: its foot stands 4.0 km above the model's ground and its top 28.0 km. Each cell's air is interpolated
linearly in geopotential to heights every 250 m and weighted bilinearly between the four cell centres; the cell
holding the summit alone, or with its eight neighbours by area, moves the air at the summit and the top by 0.06 K at
most. The design case's levels over the summit stand at 8.5, 11.7, 17.6, 25.8, 36.2 and 48.7 km above sea level, so
the tower spans three of its layers and its air is interpolated across 8–10 km between levels.

The mean is over the 3-day means of years 10–19, through the lunar month, and the tenths are the warmest and coolest
tenth of them at each height. A target between the runs takes their air interpolated by the ground's temperature,
height by height, as the targets' climates do, so it mixes the two runs' physics:

| Ground | Mean air reaches 0 °C, above the summit | It reaches −5 °C | Tower in sub-freezing air | Air at the top, 24 km up: mean (cold; warm tenth) | Open winter disk, 13.8 km up |
|---|---|---|---|---|---|
| 295 K, the lock | 14.1 km | 19.4 km | 9.9 km | −9.3 (−9.9; −8.8) °C | 0.2 °C |
| 296 K | 15.1 km | 20.4 km | 8.9 km | −8.3 (−8.9; −7.7) °C | 1.1 °C |
| 297 K | 16.1 km | 21.6 km | 7.9 km | −7.2 (−7.8; −6.7) °C | 2.0 °C |
| 298 K | 17.1 km | 22.7 km | 6.9 km | −6.2 (−6.7; −5.6) °C | 2.9 °C |
| 299 K | 18.2 km | 23.9 km | 5.8 km | −5.1 (−5.7; −4.6) °C | 3.8 °C |
| 300 K, the resort | 19.3 km | 25.1 km | 4.7 km | −4.0 (−4.6; −3.5) °C | 4.7 °C |
| 301 K | 20.5 km | 26.2 km | 3.5 km | −3.0 (−3.6; −2.4) °C | 5.6 °C |
| 302 K | 21.7 km | 27.4 km | 2.3 km | −1.9 (−2.5; −1.4) °C | 6.5 °C |
| 303 K | 23.0 km | above 28 km | 1.0 km | −0.9 (−1.5; −0.3) °C | 7.4 °C |
| 304 K | 24.3 km | above 28 km | none | 0.2 (−0.4; 0.7) °C | 8.3 °C |
| 305 K | 25.6 km | above 28 km | none | 1.3 (0.7; 1.8) °C | 9.2 °C |
| 306 K | 26.9 km | above 28 km | none | 2.3 (1.7; 2.8) °C | 10.1 °C |
| 307 K, the hot world | 27.6 km | above 28 km | none | 2.9 (2.3; 3.4) °C | 10.6 °C |

Each kelvin of the target raises the freezing level 1.11 km and warms the top's air 1.06 K. The air over the summit
changes little through the lunar month: its warm and cold tenths lie 0.5–0.6 K either side of the mean at the top.

**The warmest target with a winter at the top.** The top's mean air reaches 0 °C at 303.8 K, and its warm tenth at
303.3 K. In the design case's 6.5-hour years the top's warm tenth stands 0.68 K above its mean, against 0.56 K in the
3-day means, and its warmest mean 1.25 K, against 0.82 K; carried to the targets, the warmest 6.5-hour mean at the top
reaches 0 °C at 302.6 K. With 1.5 km of the tower below freezing under the top, the depth of the design's winter parks
from 13.5 to 15 km, the warmest target is 302.65 K, and with 2 km it is 302.25 K. The margin taken here is the 2 km: it
holds a winter park's depth of sub-freezing air under the top, keeps the top's warmest 6.5-hour mean at about
−0.4 °C, and leaves about 0.4 K for the interpolation and the coarse layers. **The tower keeps a winter at its top up to
about 302 K.**

**The interpolation.** The test's own warming, year by year from 297.5 to 306.8 K, draws its own line of the top's air
against the ground, 1.10 K a kelvin over years 5–19 and 1.15 over years 0–19. On it the top is 0.3–0.6 K colder than
the interpolation at 300 K, 0.2–0.4 K colder at 302 K, and 0.5–1.0 K colder than the design case at 295 K: the
refitted longwave cools the upper air, so the interpolation, which keeps the earlier physics at the cool end, if
anything sets the freezing level high. The model's air is the free air of 3-day means in cells 170 km across; the
tower's own structure, its people's heat and the summit's sunlit ground warm the air next to it, and a design case
at the chosen target, with one physics, settles the figures.

**Against the metropolis brief.** The [brief](../../../habitation/summit_metropolis/README.md#parks-and-air-by-height)
puts the summit's air at about 12 °C, cooling 0.9 °C a kilometre, about 1 °C at the first disk, 12.3 km up, the
winter parks at 0 to −1 °C at 13.5–15 km, freezing at 13.9 km and −9 °C at the top. The design case gives 12.3 °C,
0.90 K a kilometre, 1.5 °C at the first disk, 0.5 to −0.8 °C at 13.5–15 km, freezing at 14.1 km and −9.3 °C at the
top. The brief's figures come from the earlier design case (`A28_dim5_moon`, years 20–29) through the site product
([site_winds.py](../../../climate/gcm/site_winds.py)), which gives 12.1 °C, 1.4 °C, 0.4 to −1.0 °C, 13.9 km and
−9.4 °C. That product took the model's ground as 6,946 m, from the sea-level pressure read as if the air below the
ground were isothermal; pyburn's own extrapolation returns 7,615 m. On the model's own ground the earlier design case
is 0.7 K warmer at the summit and 0.9 K at the top, freezing at 14.8 km, and the design case's air is 0.4–0.8 K
colder than that, its middle air cooled by the gas bands and the curved light. The two differences nearly cancel: the
design case lies within 0.3 K and 0.15 km of the brief's source at every height.

**At the resort the winter moves up the tower.** At 300 K the open winter parks at 13.8 km stand at 4.7 °C, and
winter begins at 19.3 km, in the top 4.7 km of the tower: the two highest disks, at 19.9 and 21.4 km and now sealed,
and the crown, with the freezing level's air at 0.66 atm. Open winter decks would move up to them.

## The resort's best value

The author lets the resort vary a little near 300 K, to whatever suits the other constraints best
([register](../../decisions.md#the-target-climate-17-october-2026)). The candidates against each constraint, the
dimming from the dimming curve and the climate interpolated as above:

| | 298 K | 299 K | 300 K | 301 K | 302 K | 303 K |
|---|---|---|---|---|---|---|
| Tower in sub-freezing air; the top's mean air | 6.9 km; −6.2 °C | 5.8 km; −5.1 °C | 4.7 km; −4.0 °C | 3.5 km; −3.0 °C | 2.3 km; −1.9 °C | 1.0 km; −0.9 °C |
| Land wet-bulb at 2 m: median; 99th percentile | 19.3; 24.3 °C | 20.2; 25.1 °C | 21.1; 26.0 °C | 22.0; 26.8 °C | 22.9; 27.7 °C | 23.8; 28.5 °C |
| Its own dimming; with the heat beyond the sunlight | 6.7–8.4%; 10.4–13.7% | 6.0–7.5%; 9.6–12.8% | 5.2–6.5%; 8.9–11.8% | 4.5–5.6%; 8.2–10.9% | 3.8–4.6%; 7.4–9.9% | 3.0–3.7%; 6.7–9.0% |
| Sunlight its own dimming takes | 792–1,000 TW | 705–887 TW | 618–773 TW | 531–660 TW | 444–547 TW | 357–433 TW |
| Rain on land, a day; land under 0.5 mm | 2.53 mm; 33% | 2.63 mm; 31% | 2.72 mm; 29% | 2.82 mm; 28% | 2.91 mm; 26% | 3.01 mm; 25% |
| Seas' oxygen at saturation | 8.05 mg/L | 7.92 mg/L | 7.79 mg/L | 7.68 mg/L | 7.57 mg/L | 7.46 mg/L |
| Eutrophic shallows' nights below 2 mg/L | 42–87% | 50–90% | 57–93% | 61–94% | 65–95% | 68–96% |
| Under a saturated column, the added heat answered | 17.1–21.8 W/m² | 15.4–19.5 W/m² | 13.6–17.3 W/m² | 11.9–15.0 W/m² | 10.2–12.7 W/m² | 8.4–10.5 W/m² |

The rain and deserts count every cell's land, before the dimming's own cut. The oxygen comes from ecology's waters'
oxygen in the three targets ([lake_oxygen](../lake_oxygen/README.md), on `domain/ecology` at bd0097c, bound by its
hash until the branches meet): saturation interpolated by the seas' median warmth, and the eutrophic shallows' nights,
across its three lakes, interpolated by the ground's temperature between its lock, resort and hot world, a guide to a
model run at three points only. Each kelvin warmer raises the freezing level over the summit 1.1 km and the ground's
wet-bulb 99th percentile 0.84 K, takes 0.73–0.96% of the sunlight off the dimming, 87–113 TW, and adds 0.1 mm of
rain a day on land with 1.6 points less desert, and costs the seas 0.12 mg/L of oxygen and the margin under a
saturated column 1.7–2.3 W/m².

**The recommendation is 301 K, within 300–302 K.** Only the tower's winter sets a line: with its 2 km margin it holds to
about 302 K, which makes 302 K the band's warm edge. The array's dimming, the rain and the deserts favour the warm end,
and comfort at the ground, the waters' oxygen and the margin under a saturated column the cool end, with no line among
them across the band. 301 K keeps 3.5 km of the tower in winter, the top at −3.0 °C and its warmest 6.5-hour mean at
−1.7 °C, a kelvin inside the tower's line; against 300 K it takes 0.7–1.0 points off the dimming, 87–113 TW of sunlight,
adds 3.5% to the land's rain and brings the ground's wet-bulb 99th percentile to 26.8 °C, 27.1 °C by day. 302 K takes as
much again off the dimming at the line itself; 300 K keeps the ground's humid heat 0.8 K lower. The author set the
resort at 301 K on 18 October.

## What it means

**The hot world needs the plan's dimming all the same.** Its 306.5 K carries none of the gases, the glow or the beamed
power. Left unanswered they add 8.3–12.0 W/m² as the GCM weighs them, about 4–7 K more at the test's feedback, and at
full watts they would ask its clear sky to send out 7–11 W/m² more than a saturated column can, leaving the balance
wholly to unsaturated air. Keeping the hot world still takes the 3.7–5.3% of dimming that answers them, a sunlight
scale of 0.92–0.93.

**Comfort at the ground follows the target closely.** In the hot world, at 307 K, the Moon's land is hot and humid at
every latitude, the poles included. At 2 m the wet-bulb temperature passes 28 °C for 31% of the land's time, 49%
within 15° of the equator and 56% on land below 1 km, and passes 31 °C, near where young, healthy adults can no longer
keep up light activity in humid heat (Vecellio et al. 2022), for 2% of the land's time; the seas are at a bath's
34.5 °C. Its comfort lives aloft, where the air is 20–25 °C from about 8 to 12 km up at 0.92–1.00 atm, and its
freezing air from 36 km up near 0.57 atm. In the resort, at 300 K, the ground is a warm tropical one, with seas near
29 °C, the 2 m wet-bulb's 99th percentile near 26 °C, 26.3 °C by day and 25.0 °C by night, 20 °C air 6.7 km up and
freezing air above 28 km. The lock, at 295 K, is the design case.

**Every kelvin below the hot world is dimming.** 300 K takes 5.2–6.5% of the sunlight beyond the plan's additions and
295 K 8.9–11.3%, each kelvin 0.73–0.96% of it. The dimmer would go from the corrected design case's 3% of the titania
stack's light to 11.6–14.5% for 300 K and 15–19% for 295 K, with the added heat answered. The light it takes is
power the array could use, 124–402 TW at 20–30%, but the fleet's tiles and photon authority already bind
([constraints review](../solar_shield_array/constraints_review.md)): a reflecting dimmer adds to the window tiles'
light pressure, and an absorbing one returns 10% of its dimming as glow with both halves rolled and 25% as built
([glow](../solar_shield_array/glow.md)).

**Dimming takes rain and light from the living Moon.** The hot world's land rains 3.34 mm a day against the lock's
2.25 and has half its desert, 18.9% against 37.6%; on the model's own land cells, 2.94 against 1.98 mm and 20.1%
against 42.3%. Each 1% of dimming takes 0.07–0.30% of the living Moon's growth through light, and its own 0.55–1.38%
of the rain at a held temperature beyond what the cooling does.

**The gap to a saturated column widens with dimming.** With the added heat answered, the clear sky's load sits
1.3–1.8 W/m² under the saturated column's ceiling in the hot world, 14–17 W/m² at 300 K and 22–29 W/m² at 295 K. The
hot world is stable in the GCM all the same, held by its restoring gain.

**The tower's winter caps the resort near 302 K.** Each kelvin raises the freezing level over the summit 1.1 km. The
lock keeps 9.9 km of the tower in winter, the resort at 300 K 4.7 km and 302 K 2.3 km, and from about 303.8 K the
tower's top is above freezing. With 2 km of winter kept under the top, the tower holds a winter up to about 302 K,
and the author set the resort at 301 K on 18 October, within 300–302 K. Above about 296 K the open winter parks of
13.5–15 km lose their winter, which moves up to the sealed disks and the crown.

**What the choice reopens.** The register's planned climate near 295 K (9 October), the mix of climates with the
hottest days pleasantly warm (25 September) and the 5% dimmer's 295 K; the design case for the chosen target, a 20-year
GCM run of about 1.5 hours on four ranks, and the glow pair from it, two more; then geography's lakes, the atlas globe,
and the studies that read the design climate: food, people, light, the gases and glow, the CM1 rings and the Jura
Island box.

## Assumptions and limits

- **Two physics.** The design case flew PlaSim's earlier longwave and liquid optics for all cloud, so the climate at
  300 K interpolates between two models; a cooler Moon under the refit would keep some of the refit's high cloud,
  whose amount at the cooler end is unknown. A run at the chosen target settles it.
- **The feedback** is the warming test's, applied to cooling; it is the GCM's for this state and Earth-tuned in its
  clouds and convection.
- **The saturated column's ceiling** is the radiative–convective model's clear, fully saturated column with a 200 K
  stratosphere.
- **The tower's air** is the free air of 3-day means over cells about 170 km across, on a ground that T21 puts 4.0 km
  below the summit, interpolated linearly across layers 8–10 km thick; the structure, its people and the summit's own
  ground warm the air beside it. The winter margin is a judgement: a winter park's depth under the top, the warmest
  6.5-hour mean below freezing and about 0.4 K for the interpolation.
- **The waters' oxygen** between ecology's three targets is interpolated; its lake model ran at those three alone.

## Running it

```sh
climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract   # about 80 s; needs netCDF4
python -m research.studies.climate_targets.run                              # under a second
```

The extract reads the two runs' output in `climate/gcm/runs`, the glow pair's yearly records, the design case's two
6.5-hour years, and the summit's cells in the test's twenty years and the earlier design case's years 20–29, about
3 GB at most, and writes `research/runs/climate_targets/extract.json` (ignored, on the research drive), with the
sha256 and size of every file it read; [test_run.py](test_run.py) checks the binding and the arithmetic. The run
reads ecology's waters' oxygen from Git at its pinned commit until the branches meet.
