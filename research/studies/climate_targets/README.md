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
| The fleet's glow, both halves rolled | 4.24 W/m² | 1.9% | 5.34 W/m² of infrared at the GCM's efficacy of 0.71, and 0.61 W/m² of sunlight sent back |
| Beamed power used on the Moon | 0.87–4.14 W/m² | 0.4–1.8% | The power screen's 30–145 TW of beams arriving at the receivers and the towns' crowns, and the deep air's 2.6–12 TW of them |
| **Together** | **8.7–12.0 W/m²** | **3.9–5.3%** | |
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

## The three targets

Each target's dimming is the test's gain per kelvin times the kelvins it sheds, as a share of the sunlight reaching
the Moon at the test's albedo. Its climate is interpolated between the two runs by the ground's temperature, so the
295 K column is the design case's own climate.

| | The hot world: no further dimming | The resort: 300 K | The lock: 295 K |
|---|---|---|---|
| Ground | 306.5 K, settling at 306.8–307.1 | 300 K | 295 K |
| Absorbed sunlight to shed | none | 11.8–14.7 W/m² | 20.0–25.5 W/m² |
| Its dimming | none | 5.2–6.5% | 8.9–11.3% |
| With the heat beyond the sunlight | 3.9–5.3% | 9.1–11.8% | 12.7–16.6% |
| Sunlight scale of a design case for the target, from 0.970 | 0.970 | 0.907–0.919 | 0.860–0.884 |
| Sunlight scale with the added heat answered too | 0.919–0.933 | 0.855–0.882 | 0.809–0.846 |
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
| That gap with the added heat at full watts: answered; unanswered | 1.7–1.8; −10.7 to −7.4 W/m² | 14.2–17.1; 1.6–8.0 W/m² | 22.8–28.4; 10.3–19.3 W/m² |

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
1.38% (the light study). The costs in the table are each target's own dimming. The 3.9–5.3% every target shares for
the added heat costs a further 0.27–1.56% of the growth through light and 2.1–7.3% of the land's rain, which brings the
warm Moon's to 3.10–3.27 mm a day. The gap to the saturated column counts the added heat at full watts, the glow's
5.34 W/m² of infrared among them, where the dimming answers the glow at the GCM's efficacy.

## What it means

**The warm Moon needs the plan's dimming all the same.** Its 306.5 K carries none of the gases, the glow or the beamed
power. Left unanswered they add 8.7–12.0 W/m² as the GCM weighs them, about 4–7 K more at the test's feedback, and at
full watts they would ask its clear sky to send out 7–11 W/m² more than a saturated column can, leaving the balance
wholly to unsaturated air. Keeping the warm Moon still takes the 3.9–5.3% of dimming that answers them, a sunlight
scale near 0.92–0.93.

**Comfort at the ground follows the target closely.** In the hot world, at 307 K, the Moon's land is hot and humid at
every latitude, the poles included. At 2 m the wet-bulb temperature passes 28 °C for a third of the land's time,
nearly half of it near the equator and more than half in the lowlands, and passes 31 °C, near where young, healthy
adults can no longer keep up light activity in humid heat (Vecellio et al. 2022), for 3% of the land's time; the seas
are at a bath's 35 °C. Its comfort lives aloft, where the air is 20–25 °C from about 8 to 12 km up at 0.92–1.00 atm,
and its freezing air from 36 km up near 0.58 atm. In the resort, at 300 K, the ground is a warm tropical one, with
seas near 29 °C, the 2 m wet-bulb's 99th percentile near 26 °C, 20 °C air 6.5 km up and freezing air above 28 km. The
lock, at 295 K, is the design case.

**Every kelvin below the warm Moon is dimming.** 300 K takes 5.2–6.5% of the sunlight beyond the plan's additions and
295 K 8.9–11.3%. The dimmer would go from the corrected design case's 3% of the titania stack's light to 12–15% for
300 K and 15–19% for 295 K, with the added heat answered. The light it takes is power the array could use, 124–402 TW
at 20–30%, but the fleet's tiles and photon authority already bind
([constraints review](../solar_shield_array/constraints_review.md)): a reflecting dimmer adds to the window tiles'
light pressure, and an absorbing one returns 10% of its dimming as glow with both halves rolled and 25% as built
([glow](../solar_shield_array/glow.md)).

**Dimming takes rain and light from the living Moon.** The warm Moon has 48% more rain on land than the design case
and half its desert. Each 1% of dimming takes 0.07–0.30% of the living Moon's growth through light, and its own
0.55–1.38% of the rain at a held temperature beyond what the cooling does.

**The gap to a saturated column widens with dimming.** With the added heat answered, the clear sky's load sits
1.7–1.8 W/m² under the saturated column's ceiling in the warm Moon, 14–17 W/m² at 300 K and 23–28 W/m² at 295 K. The
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
- **The efficacies** are the GCM's of 13 October on the earlier radiation; the glow pair on the corrected radiation of
  17 October gave about 0.59 for the glow, and its rerun on the final radiation waits on the target.
- **The top-up's heat** spends all the fall energy on the Moon, the case without recovery.
- **The saturated column's ceiling** is the radiative–convective model's clear, fully saturated column with a 200 K
  stratosphere. It measures how far the balance leans on the air staying unsaturated; the same model's
  Manabe–Wetherald columns, with Earth's humidity falling off with height, mark no runaway below 320 K.

## Running it

```sh
climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract   # about a minute; needs netCDF4
python -m research.studies.climate_targets.run                              # under a second
```

The extract reads the two runs' output in `climate/gcm/runs` and writes `research/runs/climate_targets/extract.json`
(ignored, on the research drive), with the sha256 of every file it read; [test_run.py](test_run.py) checks the binding
and the arithmetic.
