# Where the fresh waters run coldest

4 November 2026. The seas' deep water settles at the temperature of the coldest water their surface makes, 25–30 °C at 301 K ([the seas' oxygen](../sea_oxygen/README.md)): with the Moon's axis 1.5° from the ecliptic there is no polar winter, and the lunar night cools the seas' surface by about a kelvin. This study looks for cooler water in the two places left, the lakes on the highest ground and shallow ponds cooling through the 354-hour night, as a first estimate on the design years, `A28R_301k_design_l15_snow_s934` years 10–19. [run.py](run.py) writes [results/cold_waters.json](results/cold_waters.json) (schema `terluna.research.cold-waters/1`, with the sha256 of every input) in about five minutes:

```sh
/media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.cold_waters.run
```

## What it finds

**The Open Moon's fresh water stays above about 15 °C everywhere.** The coolest standing water lies on the highest ground, 5–10 km above sea level, where deep lakes hold water near 21 °C and only shallow ponds go below 18 °C, in the lunar night's last days; the night's coldest hours, which the 3-day means leave out, could take the shallowest of them to about 15 °C.

**The air over the lakes is warm, and the lapse rate is small.** The design years' air warms going down by 1.04 K a kilometre between 0.9 and 5 km above the ground over land (0.93–1.14 across the land cells), against about 6.5 K a kilometre on Earth, so height cools a lake little. Geography's 27,865 rain-fed lakes cover 3.92 million km²; half their area lies above 4.0 km above sea level and the highest lake stands at 9.9 km. Carried to each lake's level from its cell's mean ground, the air over the lakes averages 24.9 °C, area-weighted, and 23.2 °C in each lake's coldest 3-day spell.

| Coldest 3-day spell of the air at the lake | Lakes | Their area |
|---|---|---|
| below 24 °C | 13,418 | 2.49 million km² |
| below 22 °C | 4,965 | 0.95 million km² |
| below 20 °C | 543 | 220,000 km² |
| below 18 °C | 15 | 583 km² |
| below 15 °C | none | |

The coldest air stands over lakes high above their cells' mean ground, mostly near the far side's equator: lakes at 5.5–5.7 km by 0.3° N, 139° E, 4 km above their cell's ground, take a coldest spell of 17.5 °C and a mean of 20.5 °C, and the highest lake, at 9.9 km by 2.5° N, 199° E, 17.7 and 19.1 °C. The lapse rates' spread across the land cells moves those coldest spells by about 0.4 K either way.

**A lake's water runs warmer than its air.** The project's lake model ([biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py)) run at the coldest lake's site, with the cell's means carried up to the lake's level, keeps a 30 m lake's surface at 26.7 °C on average and no colder than 21.9 °C in a 3-day mean: the lunar day's 354 hours of sunlight heat the water well above the air. Its floor averages 21.7 °C and never falls below 21.1 °C in the ten years, the coldest standing water a lake keeps. Deeper lakes' floors settle at their surface's coldest water as well, so the deep highland lakes hold water near 21 °C under a warmer surface.

**Shallow ponds on the high ground cool below 18 °C at the end of the night, and warm far above it by day.**

| Site | Water | Lowest | Lowest 3-day mean | Mean | Hours of each lunar day below 20 °C; below 18 °C |
|---|---|---|---|---|---|
| The coldest lake's, 5.7 km above sea level | 0.5 m pond | 16.9 °C | 17.1 °C | 25.8 °C | 244; 39 |
| | 1 m pond | 17.3 | 17.6 | 26.0 | 181; 4.5 |
| | 2 m pond | 18.9 | 19.1 | 26.5 | 24; 0 |
| The coldest land cell, 8.4 km above sea level | 0.5 m pond | 17.6 | 17.8 | 25.0 | 243; 1.1 |
| | 2 m pond | 19.4 | 19.7 | 25.7 | 2.8; 0 |
| A typical land cell, 41.5° N, 1.8 km above sea level | 0.5 m pond | 22.5 | 22.7 | 28.4 | 0; 0 |
| | 2 m pond | 23.3 | 23.5 | 28.7 | 0; 0 |

A lunar day lasts 709 hours. The ponds above are open to the 10 m wind; sheltered ponds, with 30% of it, run the same at night, where the model's 0.5 m/s floor and free convection carry the exchange, and up to 0.8 K warmer on average, from the day.

## Limits

- The estimate reads the design years before the GCM and CM1 agree on the rain and the air above it ([register, 1 November](../../decisions.md)); the night's air and the lapse rate are among what closing that gap could move.
- A lake far above its cell's mean ground takes the free air's lapse rate down from the cell, as the lake studies carry it, with the sky's infrared scaled by the fourth power of the air's temperature. The air over a high massif and the thinner, drier column above it are not modelled; the drier sky would radiate less and cool the high waters further, by an amount this estimate does not size.
- The forcing is the GCM's 3-day means made linear in time. The coldest hours within a mean run about 2 K below its mean over the coldest lakes (the run's minimum temperatures), which the ponds' lowest values do not include, and sunrise and sunset are smeared over a mean.
- The ponds are pans 1 hectare across with the lake model's water clarity and no shading of the water by day; shade would cool them by day and leave the night as it is.
