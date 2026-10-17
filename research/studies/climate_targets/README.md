# What each target climate takes

On 17 October the author opened the target climate as a question: whether the Moon is brought back to the design
temperature, with more dimming tested, or kept hotter and wetter provided it stays comfortable, with cooler and wintry
weather aloft, long coasts and perhaps a life chiefly of recreation
([register](../../decisions.md#the-target-climate-17-october-2026)). New runs wait on the choice. This study sets
three targets side by side across the domains, at the joint integration of the five branches
([record](../../integration/five-branches.md)): the hot world, the warm Moon the corrected GCM settles on; the lock,
the design temperature of 295 K; and between them the resort, a tropical resort Moon, which the study sets at 300 K.
With them it sets the heat the plan adds to the Moon beyond its sunlight, which every target answers.
[run.py](run.py) is arithmetic on the GCM's output and the domains' committed products
([results/climate_targets.json](results/climate_targets.json)); no model runs.

## Where the corrected Moon settles

The design case (`A28_corrected_970`) holds 295.0 K at the design sunlight scale of 0.970 on PlaSim's earlier
longwave. With the longwave refitted by height and ice optics for the cold cloud, the 20-year test at the same sunlight
(`A28_corrected_970_lw_ice`) averages 306.5 K over years 10–19
([climate/gcm](../../../climate/gcm/README.md#where-the-corrected-moon-settles)). Its top of the air gains
1.65–2.16 W/m² less for each kelvin it warms (fitted over years 5–19 and 0–19), so it settles at 306.8–307.1 K, and that
restoring gain holds it. The radiative–convective model's Manabe–Wetherald columns, whose humidity falls off with
height as Earth's does, show no runaway up to 320 K
([radiative–convective column](../../../atmosphere/radiative_convective/README.md)). The warm Moon absorbs 225.4 W/m²
of sunlight, 236.2 W/m² under its clear sky. A fully saturated clear column there can send out at most 239.0 W/m², so
a saturated, cloudless warm Moon would only just balance; its unsaturated clear sky sends out 241.5 W/m². Its cloud
reflects 10.8 W/m² of sunlight and holds back 16.8 W/m² of infrared, 6.0 W/m² of net warming. Its balance rests on its
air staying unsaturated.

## The heat beyond the sunlight

| What the plan adds | As absorbed sunlight | Share of the absorbed sunlight | Basis |
|---|---|---|---|
| Methane and N2O at 1.9 and 0.1 ppm | 3.58 W/m² | 1.6% | The deep column's all-sky forcing of the held pair; the GCM held it with dimming of the same watts |
| The fleet's glow, both halves rolled | 3.59–4.24 W/m² | 1.6–1.9% | 5.34 W/m² of infrared at the GCM's efficacy of 0.59–0.71, and 0.61 W/m² of sunlight sent back |
| Beamed power used on the Moon | 0.87–4.14 W/m² | 0.4–1.8% | The power screen's 30–145 TW of beams arriving at the receivers and the towns' crowns, and the deep air's 2.6–12 TW of them |
| **Together** | **8.0–12.0 W/m²** | **3.6–5.3%** | |
| The water top-up's fall, the first thousand years | 1.5–9.5 W/m² | 0.7–4.2% | 1,500–4,100 km³ a year at the central values, falling 1.2–2.7 MJ/kg, all of it spent on the Moon |

The shares are of the warm Moon's absorbed sunlight, 225.4 W/m²; the glow study and the register give the glow's
dimming as 2.0% of the earlier design case's sunlight at the top of the air, 293.3 W/m² with 211.2 absorbed. The
gases, the glow and the power last as long as the plan does, and the plan answers them with dimming, so they come
on top of whatever a target needs. The top-up is the water the author keeps delivering to balance what the thawing
crust takes up ([geography](../../../geography/README.md#groundwater-how-fast-the-crust-takes-the-water)). Its fall
into the Moon's own well releases 1.8–4.0 W/m² through the years to 10,000 at the central values, and 12–28 W/m² in
the fast case's busiest centuries. Resources' first screen already asks the build's deliveries for energy recovery on
the way down: tethers, elevators or regenerative lowering
([resources](../../../resources/README.md)); the top-up needs the same for millennia after the build.

**The glow's efficacy runs from 0.59 to 0.71.** The 0.71 is the held branch of 13 October on the earlier radiation
([gases, glow and dimming](../gases_glow_dimming/README.md#what-it-finds)). The 0.59 is the glow pair of 17 October on
the corrected radiation, the design case's own before the longwave refit
([climate/gcm](../../../climate/gcm/README.md#the-gases-the-glow-and-the-dimming)), found the same way from the two
branches' yearly records, which the extract binds by their hashes. Over years 10–19 the branch matched watt for watt,
with 4.22% of added dimming, stands at 293.45 K and the held branch, with 3.52%, at 294.41 K: 1.37 K for each point of
dimming. Carried to the design case's 295.02 K, the gases and the glow together take 3.07% of its sunlight. The gases
keep their 1.69% of 13 October, since no gases-only branch flew on this radiation, which leaves the glow 1.38%:
3.15 W/m² of the design case's 228.5 W/m² of absorbed sunlight against its 5.34 W/m² of infrared. A point of dimming
takes about 8% more absorbed sunlight here than on the earlier radiation, the curved air's light among it, so the
same share counted in points, as 13 October counted it, is 0.55; with the gases answered by their forcing alone,
3.58 W/m², the glow's would be 0.65. At 0.59 the glow's term falls from 4.24 to 3.59 W/m², which takes 0.29 points off
the dimming every target shares.

## The three targets

Each target's dimming is the test's gain per kelvin times the kelvins it sheds, as a share of the sunlight reaching
the Moon at the test's albedo. Its climate is interpolated between the two runs by the ground's temperature, so the
295 K column is the design case's own climate.

| | The hot world: no further dimming | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Ground | 306.5 K, settling at 306.8–307.1 | 300 K | 295 K |
| Absorbed sunlight to shed | none | 11.8–14.7 W/m² | 20.0–25.5 W/m² |
| Its dimming | none | 5.2–6.5% | 8.9–11.3% |
| With the heat beyond the sunlight | 3.6–5.3% | 8.8–11.8% | 12.5–16.6% |
| Sunlight scale of a design case for the target, from 0.970 | 0.970 | 0.907–0.919 | 0.860–0.884 |
| Sunlight scale with the added heat answered too | 0.919–0.935 | 0.855–0.885 | 0.809–0.849 |
| Sunlight the target's own dimming takes; electricity at 20–30% | none | 618–773 TW; 124–232 TW | 1,052–1,340 TW; 210–402 TW |
| Growth lost through light to the target's own dimming | none | 0.37–1.9% | 0.62–3.3% |
| Rain on land, a day | 3.34 mm | 2.72 mm; 2.48–2.64 with the dimming's own cut | 2.25 mm; 1.90–2.14 |
| Land under 0.5 mm a day | 18.9% | 29% | 38% |
| Land wet-bulb at 2 m: median; 99th percentile | 26.9; 31.9 °C | 21.1; 26.2 °C | 16.7; 21.9 °C |
| Land wet-bulb 0.9 km up: median; 99th percentile | 26.3; 31.3 °C | 20.4; 25.6 °C | 15.9; 21.2 °C |
| Land air, median: at 2 m; 0.9 km up | 33.5; 32.1 °C | 26.1; 24.7 °C | 20.4; 19.0 °C |
| Seas' surface, median | 34.6 °C | 28.7 °C | 24.2 °C |
| Air at 20 °C over the land | 12.2 km up | 6.5 km up | below the lowest level, 0.9 km up |
| Freezing level over the land | 35.8 km | 28.0 km | 22.7 km |
| Clear-sky absorbed sunlight, and its gap to a saturated column's 239.0 W/m² | 236.2; 2.8 W/m² | 220.8–223.8; 15.1–18.2 W/m² | 209.4–215.2; 23.8–29.5 W/m² |
| That gap with the added heat at full watts: answered; unanswered | 1.0–1.8; −10.7 to −7.4 W/m² | 13.3–17.3; 1.6–8.0 W/m² | 22.0–28.6; 10.3–19.3 W/m² |

The heat is in the 3-day means, each about a tenth of the lunar day, at 2 m and at the lowest model level, about
0.9 km up. The land's figures are weighted by land area and the seas' by sea area. The output gives the air at 2 m and
no humidity there, so the wet-bulb temperature at 2 m takes the lowest level's specific humidity, capped at
saturation, as in air mixed through the lowest kilometre. In the hot world the wet-bulb temperature at 2 m passes
28 °C for 32% of the land's time: 46% within 15° of the equator, 57% on land below 1 km and 5.7% on ground above 4 km.
It passes 31 °C for 3.3%. At 0.9 km the same shares are 23.6%, 37%, 46% and 2.5%, with 1.6% above 31 °C. Until
17 October's correction the figures counted every land cell alike, so the small polar cells, a third of the land's
cells and a sixth of its area, weighed as much as the equatorial ones, and the lowest level's time above 28 °C read
17.6%; the coasts also took in polar cells across the pole. The rain's own cut takes the GCM's
0.55% of the rain for each 1% of dimming at a held temperature (the gases, glow and dimming study) and Earth's models'
1.38% (the light study). The costs in the table are each target's own dimming. The 3.6–5.3% every target shares for
the added heat costs a further 0.25–1.56% of the growth through light and 2.0–7.3% of the land's rain, which brings the
warm Moon's to 3.10–3.27 mm a day. The gap to the saturated column counts the added heat at full watts, the glow's
5.34 W/m² of infrared among them, where the dimming answers the glow at the GCM's efficacy; its ends take the ends of
the feedback, the glow's efficacy and the beamed power that leave the clear sky the most load and the least.

## The dimming at every temperature

Each kelvin the Moon is held below where it settles takes 0.73–0.96% of the sunlight beyond the plan's additions: each
fit's own feedback, 1.65 and 2.16 W/m² a kelvin, over the test's 225.4 W/m² of absorbed sunlight, counted from that
fit's own balance. The sunlight scales are those a design case for the target would fly from 0.970, before and after
the heat beyond the sunlight is answered too; the product carries every half kelvin from 295 to 307 K.

| Ground | Its dimming | Sunlight scale of its design case | With the heat beyond the sunlight | Sunlight scale with it answered |
|---|---|---|---|---|
| 295 K | 8.9–11.3% | 0.860–0.884 | 12.5–16.6% | 0.809–0.849 |
| 296 K | 8.2–10.4% | 0.870–0.891 | 11.7–15.7% | 0.818–0.856 |
| 297 K | 7.4–9.4% | 0.879–0.898 | 11.0–14.7% | 0.827–0.863 |
| 298 K | 6.7–8.4% | 0.888–0.905 | 10.3–13.7% | 0.837–0.871 |
| 299 K | 6.0–7.5% | 0.897–0.912 | 9.5–12.8% | 0.846–0.878 |
| 300 K | 5.2–6.5% | 0.907–0.919 | 8.8–11.8% | 0.855–0.885 |
| 301 K | 4.5–5.6% | 0.916–0.927 | 8.1–10.9% | 0.865–0.892 |
| 302 K | 3.8–4.6% | 0.925–0.934 | 7.3–9.9% | 0.874–0.899 |
| 303 K | 3.0–3.7% | 0.935–0.941 | 6.6–9.0% | 0.883–0.906 |
| 304 K | 2.3–2.7% | 0.944–0.948 | 5.9–8.0% | 0.892–0.913 |
| 305 K | 1.6–1.7% | 0.953–0.955 | 5.1–7.0% | 0.902–0.920 |
| 306 K | 0.8% | 0.962 | 4.4–6.1% | 0.911–0.928 |
| 306.8–307.1 K, where it settles | none | 0.970 | 3.6–5.3% | 0.919–0.935 |

The two fits cross near 305.9 K, so the fit of years 5–19 gives the high end below it and the low end above it.

## Day and night

The output's czen, the cosine of the Sun's zenith angle averaged over each 3-day mean, tells the day's means from the
night's. A mean is the day's when the Sun is up for more than half of it, which a czen above 0.078 times the cosine of
the latitude marks: the czen of a mean of a tenth of the lunar day centred on sunrise or sunset. The day and the night
then take half of the land's time each. In the land's 3-day means, the day's figure first and the night's second:

| Land: the day's means; the night's | The hot world | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Wet-bulb at 2 m: median | 27.1; 26.8 °C | 21.5; 20.8 °C | 17.3; 16.2 °C |
| Wet-bulb at 2 m: 99th percentile | 32.1; 31.6 °C | 26.5; 25.6 °C | 22.3; 21.0 °C |
| Air at 2 m: median | 34.6; 32.9 °C | 27.6; 25.3 °C | 22.2; 19.5 °C |
| Air at 2 m: 99th percentile | 39.6; 36.0 °C | 33.2; 28.9 °C | 28.2; 23.5 °C |
| Wet-bulb 0.9 km up: median | 26.5; 26.1 °C | 20.8; 20.1 °C | 16.5; 15.5 °C |
| Wet-bulb 0.9 km up: 99th percentile | 31.5; 31.1 °C | 25.9; 24.9 °C | 21.6; 20.3 °C |
| Air 0.9 km up: median | 33.1; 31.4 °C | 26.1; 23.9 °C | 20.8; 18.1 °C |
| Air 0.9 km up: 99th percentile | 37.9; 34.5 °C | 31.4; 27.4 °C | 26.5; 22.1 °C |

The night eases the air more than the humid heat. In the hot world the night's air at 2 m is 1.7 K cooler at the
median and 3.7 K at the 99th percentile, and its wet-bulb temperature 0.3 and 0.4 K: at a fixed humidity a kelvin of
air moves the hot world's wet-bulb temperature 0.28 K and the lock's 0.39 K. At 2 m the hot world's wet-bulb
temperature passes 28 °C for 35% of the day's land time and 29% of the night's, and 31 °C for 4.1% and 2.5%; at
0.9 km, 26.5% and 20.8%, and 2.0% and 1.1%. In the lock the night takes 1.1 K off the wet-bulb median at 2 m and
2.7 K off the air's. The resort is interpolated between the two runs, as in the targets' table.

## The weather shorter than three days

`A28_corrected_970_eddies` carried the design case on from its year 19 for two years, its years 0 and 1, writing
every 6.5 hours: eleven means to each 3-day mean, 1,320 a year. It saved the output every run saves, the air's
temperature, humidity and winds at every level, the surface pressure, the 2 m air, the rain, the fluxes and the Sun's
cosine. Over those two years the land's heat in the 6.5-hour means stands against the 3-day means their elevens make,
with each extreme's excess and the ratio of its height above the median:

| Design case, land, its two 6.5-hour years | Median: 6.5 hours; 3-day means | 99th percentile: the same (excess; ratio) | Top: the same (excess; ratio) |
|---|---|---|---|
| Wet-bulb at 2 m | 16.67; 16.64 °C | 22.2; 21.9 °C (+0.32 K; 1.058) | 25.1; 24.1 °C (+0.94 K; 1.122) |
| Air at 2 m | 20.37; 20.37 °C | 27.9; 27.7 °C (+0.22 K; 1.029) | 30.8; 30.1 °C (+0.66 K; 1.067) |
| Wet-bulb 0.9 km up | 15.93; 15.91 °C | 21.5; 21.2 °C (+0.32 K; 1.056) | 24.3; 23.5 °C (+0.81 K; 1.103) |
| Air 0.9 km up | 18.95; 18.94 °C | 26.2; 26.0 °C (+0.21 K; 1.029) | 28.8; 28.3 °C (+0.57 K; 1.060) |

The 3-day means smooth the land's 99th percentile by 0.2–0.3 K and its top by 0.6–0.9 K. The 3-day means made from the
two years match the design case's own 3-day output over its years 10–19 within 0.1 K at the median and the 99th
percentile (the wet-bulb temperature at 2 m: 16.64 and 21.86 °C against 16.69 and 21.94 °C). Stretching each target's
99th percentile above its median by these ratios gives, as a sensitivity from the design case alone, the 99th
percentiles in 6.5-hour means, with the 3-day means' in brackets:

| Land, 99th percentile | The hot world | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Wet-bulb at 2 m | 32.2 (31.9) °C | 26.5 (26.2) °C | 22.2 (21.9) °C |
| Air at 2 m | 39.4 (39.2) °C | 32.9 (32.7) °C | 28.0 (27.7) °C |
| Wet-bulb 0.9 km up | 31.6 (31.3) °C | 25.9 (25.6) °C | 21.5 (21.2) °C |
| Air 0.9 km up | 37.7 (37.5) °C | 31.2 (31.0) °C | 26.2 (26.0) °C |

## What it means

**The warm Moon needs the plan's dimming all the same.** Its 306.5 K carries none of the gases, the glow or the beamed
power. Left unanswered they add 8.0–12.0 W/m² as the GCM weighs them, about 4–7 K more at the test's feedback, and at
full watts they would ask its clear sky to send out 7–11 W/m² more than a saturated column can, leaving the balance
wholly to unsaturated air. Keeping the warm Moon still takes the 3.6–5.3% of dimming that answers them, a sunlight
scale of 0.92–0.94.

**Comfort at the ground follows the target closely.** In the hot world, at 307 K, the Moon's land is hot and humid at
every latitude, the poles included. At 2 m the wet-bulb temperature passes 28 °C for a third of the land's time,
nearly half of it near the equator and more than half in the lowlands, and passes 31 °C, near where young, healthy
adults can no longer keep up light activity in humid heat (Vecellio et al. 2022), for 3% of the land's time; the seas
are at a bath's 35 °C. Its nights cool the air at 2 m by 1.7 K at the median and its wet-bulb temperature by
0.3 K, which still passes 28 °C for 29% of the night's land time. Its comfort lives aloft, where the air is 20–25 °C
from about 8 to 12 km up at 0.92–1.00 atm, and its freezing air from 36 km up near 0.58 atm. In the resort, at 300 K,
the ground is a warm tropical one, with seas near 29 °C, the 2 m wet-bulb's 99th percentile near 26 °C, 26.5 °C by
day and 25.6 °C by night, 20 °C air 6.5 km up and freezing air above 28 km. The lock, at 295 K, is the design case.
Weather shorter than three days adds about 0.3 K to each target's wet-bulb 99th percentile, as the design case
measures it.

**Every kelvin below the warm Moon is dimming.** 300 K takes 5.2–6.5% of the sunlight beyond the plan's additions and
295 K 8.9–11.3%, each kelvin 0.73–0.96% of it. The dimmer would go from the corrected design case's 3% of the titania
stack's light to 11.5–14.5% for 300 K and 15–19% for 295 K, with the added heat answered. The light it takes is
power the array could use, 124–402 TW at 20–30%, but the fleet's tiles and photon authority already bind
([constraints review](../solar_shield_array/constraints_review.md)): a reflecting dimmer adds to the window tiles'
light pressure, and an absorbing one returns 10% of its dimming as glow with both halves rolled and 25% as built
([glow](../solar_shield_array/glow.md)).

**Dimming takes rain and light from the living Moon.** The warm Moon has 48% more rain on land than the design case
and half its desert. Each 1% of dimming takes 0.07–0.30% of the living Moon's growth through light, and its own
0.55–1.38% of the rain at a held temperature beyond what the cooling does.

**The gap to a saturated column widens with dimming.** With the added heat answered, the clear sky's load sits
1.0–1.8 W/m² under the saturated column's ceiling in the warm Moon, 13–17 W/m² at 300 K and 22–29 W/m² at 295 K. The
warm Moon is stable in the GCM all the same, held by its restoring gain; the narrow gap measures how far its balance
leans on its air staying unsaturated.

**What the choice reopens.** The register's planned climate near 295 K (9 October), the mix of climates with the
hottest days pleasantly warm (25 September) and the 5% dimmer's 295 K; the design case for the chosen target, a 20-year
GCM run of about 1.5 hours on four ranks, and the glow pair from it, two more; then geography's lakes, the atlas globe,
and the studies that read the design climate: food, people, light, the gases and glow, the CM1 rings and the Jura
Island box.

## Assumptions and limits

- **Two physics.** The design case flew PlaSim's earlier longwave and liquid optics for all cloud, so the climate at
  300 K and 295 K interpolates between two models; a cooler Moon under the refit would keep some of the refit's high
  cloud, whose amount at the cooler end is unknown. A run at the chosen target settles it.
- **The feedback** is the warming test's, applied to cooling; it is the GCM's for this state and Earth-tuned in its
  clouds and convection.
- **The heat at the ground** is in 3-day means, each about a tenth of the lunar day, so the means follow the day and
  the night and average away weather shorter than three days. The wet-bulb temperature at 2 m carries the lowest
  level's humidity down to the ground, as air mixed by the day's heating does.
- **Day and night** are told apart within those means: a day's mean spans 36° of the Sun's travel, and a mean that
  holds sunrise or sunset goes to the half that holds more of it.
- **The weather shorter than three days** is the design case's alone, two years of 6.5-hour means against their own
  3-day means. Carried to the hot world and the resort, it assumes their short weather stretches the tail as the
  design case's does. The 6.5-hour means are themselves means over T21's cells, about 170 km across.
- **The efficacies** are the GCM's: the gases' 1 from 13 October, and the glow's 0.59–0.71, from the glow pair on the
  corrected radiation of 17 October to the held branch of 13 October on the earlier radiation. The 0.59 keeps the
  gases' 1.69% of 13 October, since no gases-only branch flew on the corrected radiation, and its branches still drift
  over years 10–19, the matched one cooling 0.11 K a decade. The pair's rerun on the final radiation waits on the
  target.
- **The top-up's heat** spends all the fall energy on the Moon, the case without recovery.
- **The saturated column's ceiling** is the radiative–convective model's clear, fully saturated column with a 200 K
  stratosphere. It measures how far the balance leans on the air staying unsaturated; the same model's
  Manabe–Wetherald columns, with Earth's humidity falling off with height, mark no runaway below 320 K.

## Running it

```sh
climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract   # about a minute; needs netCDF4
python -m research.studies.climate_targets.run                              # under a second
```

The extract reads the two runs' output in `climate/gcm/runs`, the glow pair's yearly records and the design case's
two 6.5-hour years, about 2 GB at most, and writes `research/runs/climate_targets/extract.json` (ignored, on the
research drive), with the sha256 of every file it read; [test_run.py](test_run.py) checks the binding and the
arithmetic.
