# What each target climate takes

On 17 October the author opened the target climate as a question: whether the Moon is brought back to the design
temperature, with more dimming tested, or kept hotter and wetter provided it stays comfortable, with cooler and wintry
weather aloft, long coasts and perhaps a life chiefly of recreation
([register](../../decisions.md#the-target-climate-17-october-2026)). New runs wait on the choice. This study sets
three targets side by side across the domains, at the joint integration of the five branches
([record](../../integration/five-branches.md)): the warm Moon the corrected GCM settles on, the design temperature
of 295 K, and 300 K between them. With them it sets the heat the plan adds to the Moon beyond its sunlight, which
every target answers. [run.py](run.py) is arithmetic on the GCM's output and the domains' committed products
([results/climate_targets.json](results/climate_targets.json)); no model runs.

## Where the corrected Moon settles

The design case (`A28_corrected_970`) holds 295.0 K at the design sunlight scale of 0.970 on PlaSim's earlier
longwave. With the longwave refitted by height and ice optics for the cold cloud, the 20-year test at the same sunlight
(`A28_corrected_970_lw_ice`) averages 306.5 K over years 10–19
([climate/gcm](../../../climate/gcm/README.md#where-the-corrected-moon-settles)). Its top of the air gains 1.65–2.16
W/m² less for each kelvin it warms (fitted over years 0–19 and 5–19), so it settles at 306.8–307.1 K. It absorbs
225.4 W/m² of sunlight. A fully saturated column can send out at most 239.0 W/m²
([radiative–convective column](../../../atmosphere/radiative_convective/README.md)); the test's clear sky already
sends out 241.5 W/m², so its balance rests on its dry air and its cloud. The same model's Manabe–Wetherald columns
show no runaway up to 320 K.

## The heat beyond the sunlight

| What the plan adds | As absorbed sunlight | Share of the absorbed sunlight | Basis |
|---|---|---|---|
| Methane and N2O at 1.9 and 0.1 ppm | 3.58 W/m² | 1.6% | The deep column's all-sky forcing of the held pair; the GCM held it with dimming of the same watts |
| The fleet's glow, both halves rolled | 4.24 W/m² | 1.9% | 5.34 W/m² of infrared at the GCM's efficacy of 0.71, and 0.61 W/m² of sunlight sent back |
| Beamed power used on the Moon | 0.87–4.14 W/m² | 0.4–1.8% | The power screen's 30–145 TW at the loads and the deep air's 2.6–12 TW of the beams |
| **Together** | **8.7–12.0 W/m²** | **3.9–5.3%** | |
| The water top-up's fall, the first thousand years | 1.5–9.5 W/m² | 0.7–4.2% | 1,500–4,100 km³ a year at the central values, falling 1.2–2.7 MJ/kg, all of it spent on the Moon |

The gases, the glow and the power last as long as the plan does, and the plan answers them with dimming, so they come
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

| | No further dimming | 300 K | 295 K |
|---|---|---|---|
| Ground | 306.5 K, settling at 306.8–307.1 | 300 K | 295 K |
| Absorbed sunlight to shed | none | 11.3–15.3 W/m² | 19.5–26.1 W/m² |
| Its dimming | none | 5.0–6.8% | 8.7–11.6% |
| With the heat beyond the sunlight | 3.9–5.3% | 8.9–12.1% | 12.5–16.9% |
| The GCM's sunlight scale for that, from 0.970 | 0.919–0.933 | 0.853–0.884 | 0.806–0.849 |
| Sunlight the target's dimming takes; electricity at 20–30% | none | 593–806 TW; 119–242 TW | 1,027–1,373 TW; 205–412 TW |
| Growth lost through light | none | 0.35–2.0% | 0.6–3.4% |
| Rain on land, a day | 3.34 mm | 2.72 mm; 2.47–2.65 with the dimming's own cut | 2.25 mm; 1.89–2.14 |
| Land under 0.5 mm a day | 18.9% | 29% | 38% |
| Land wet-bulb: median; 99th percentile | 25.9; 31.0 °C | 19.9; 25.3 °C | 15.4; 21.0 °C |
| Land air, median | 32.3 °C | 24.5 °C | 18.6 °C |
| Seas' surface, median | 34.0 °C | 27.7 °C | 23.0 °C |
| Air at 20 °C over the land | 12.1 km up | 6.4 km up | below the lowest level, 0.9 km up |
| Freezing level over the land | 35.7 km | 27.9 km | 22.6 km |
| Margin to the saturated column's 239.0 W/m², the added heat answered; unanswered | 13.6; 1.7–4.9 W/m² | 24.9–29.0; 12.9–20.3 W/m² | 33.2–39.7; 21.2–31.1 W/m² |

The heat at the ground is in the 3-day means at the lowest model level, about 0.9 km up, pooled over the land's cells.
In the warm Moon the wet-bulb temperature passes 28 °C for 17.6% of the land's time: 34% within 15° of the equator,
36% on land below 1 km and 2.3% on ground above 4 km. It passes 31 °C for 1.0%. The rain's own cut takes the GCM's
0.55% of the rain for each 1% of dimming at a held temperature (the gases, glow and dimming study) and Earth's models'
1.38% (the light study).

## What it means

**The warm Moon needs the plan's dimming all the same.** Its 306.5 K carries none of the gases, the glow or the beamed
power. Left unanswered they add 8.7–12.0 W/m², about 4–7 K more at the test's feedback, and bring the absorbed
sunlight with the added heat within 2–5 W/m² of what a saturated column can send out. Keeping the warm Moon still
takes the 3.9–5.3% of dimming that answers them, a sunlight scale near 0.92–0.93.

**Comfort at the ground follows the target closely.** At 307 K the Moon's land is hot and humid at every latitude,
the poles included. The wet-bulb temperature passes 28 °C for a third of the time near the equator and in the
lowlands, and 31 °C, near where young, healthy adults can no longer keep up light activity in humid heat (Vecellio
et al. 2022), for 1% of the land's time; the seas are at a bath's 34 °C. Its comfort lives aloft, where the air is
20–25 °C from about 8 to 12 km up at 0.92–1.00 atm, and its freezing air from 36 km up near 0.58 atm. At 300 K the
ground is a warm tropical one, with seas near 28 °C and the wet-bulb's 99th percentile near 25 °C, 20 °C air 6 km up
and freezing air above 28 km. At 295 K it is the design case.

**Every kelvin below the warm Moon is dimming.** 300 K takes 5.0–6.8% of the sunlight beyond the plan's additions and
295 K 8.7–11.6%. The dimmer would go from the design's 3% of the titania stack's light to 12–15% for 300 K and 15–19%
for 295 K, with the added heat answered. The light it takes is power the array could use, 119–412 TW at 20–30%, but
the fleet's tiles and photon authority already bind ([constraints review](../solar_shield_array/constraints_review.md)):
a reflecting dimmer adds to the window tiles' light pressure, and an absorbing one returns 14–25% of its dimming as
glow ([glow](../solar_shield_array/glow.md)).

**Dimming takes rain and light from the living Moon.** The warm Moon has 44% more rain on land than the design case
and half its desert. Each 1% of dimming takes 0.07–0.30% of the living Moon's growth through light, and its own
0.55–1.38% of the rain at a held temperature beyond what the cooling does.

**The margin against a runaway widens with dimming.** At 295 K the Moon absorbs 33–40 W/m² less than a saturated
column can send out, at 300 K 25–29 W/m², and the warm Moon 13.6 W/m², all with the added heat answered.

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
- **The heat at the ground** is in 3-day means; afternoons within them are hotter.
- **The efficacies** are the GCM's of 13 October on the earlier radiation; the glow pair on the corrected radiation of
  17 October gave about 0.59 for the glow, and its rerun on the final radiation waits on the target.
- **The top-up's heat** spends all the fall energy on the Moon, the case without recovery.

## Running it

```sh
climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract   # about a minute; needs netCDF4
python -m research.studies.climate_targets.run                              # under a second
```

The extract reads the two runs' output in `climate/gcm/runs` and writes `research/runs/climate_targets/extract.json`
(ignored, on the research drive), with the sha256 of every file it read; [test_run.py](test_run.py) checks the binding
and the arithmetic.
