# The nearside sea through a lunar cycle

The nearside sea holds half of the Open Moon's standing water: 5.29 million km²
from Oceanus Procellarum and Mare Frigoris to Fecunditatis and Nectaris, mean
depth 648 m, along the Earth-facing coasts. This study follows its waves through
two lunar cycles of the recovered weather, with the physics of the
[Smythii–Marginis cycle](cycle.md), and measures the second cycle.
[nearside.py](nearside.py) runs the waves, [nearside_analysis.py](nearside_analysis.py)
measures them, and the [product](results/nearside.json) holds the statistics,
maps, reference points, swell and coastal power. [shore_month.py](shore_month.py)
follows four shores through the month with their tide
([product](results/shore_month.json)).

## The sea through the cycle

Through the second lunar cycle the nearside sea's significant wave height
averages 1.39 m over its area and time, against 0.88 m in Smythii–Marginis, and
its area mean stays between 1.02 and 1.84 m.

| Hs at least | Share of the nearside sea's area and time | At its widest | Smythii–Marginis, area and time |
|---|---:|---:|---:|
| 0.5 m | 94% | 98% | 81% |
| 1 m | 68% | 87% | 34% |
| 2 m | 17% | 40% | 1.3% |
| 3 m | 2.6% | 10% | 0.05% |
| 4 m | 0.18% | 1.9% | none |

The largest Hs is 4.91 m, in Mare Imbrium at 11.9° W, 39.9° N over 1,292 m of
water, 26.2 days into the cycle, with a 19.3 s mean period and a 22.2 s peak
period. Half the sea's nodes reach at least 2.9 m at some time in the cycle and
a tenth reach 3.9 m. The height exceeded 5% of the time is 2.4 m at the median
node, 1.7 times its mean height, and reaches 4.3 m in north-western Oceanus
Procellarum, by Harding. Where Hs is at least 0.5 m, mean periods run from 8.3 to
15.1 s (10th to 90th percentile, median 11.5 s), and metre-scale seas peak at
16.6 s at the median and 20.1 s at the 90th percentile. These waves are 34 to
130 m long and travel at 3 to 6 m/s. Weighted by wave energy, the mean period
runs from 10.3 to 14.9 s across the sea (10th to 90th percentile of nodes), longest
in western Oceanus Procellarum (14.5 s on average, 16.7 s at most, by von Braun)
and Mare Imbrium (14.1 s) and shortest in Serenitatis and Crisium (11.2 s).

North of 10° N the waves run steadily toward the south-west: their
energy-weighted mean direction points 247° counterclockwise from east, with a
resultant of 0.69. In the southern seas the direction varies through the cycle
(resultant 0.16), with a northward tendency. The seas travel toward the equator
from both sides, as an equatorward surface flow drives them.

| Reference point | Depth | Mean Hs | Largest Hs | Hours with Hs ≥ 1 m (longest stretch) | Hours with Hs ≥ 2 m | Swell's share of the energy |
|---|---:|---:|---:|---:|---:|---:|
| Oceanus Procellarum | 499 m | 1.75 m | 3.36 m | 565 (265) | 273 | 64% |
| Mare Imbrium | 1,124 m | 1.76 m | 4.32 m | 690 (582) | 201 | 58% |
| Mare Frigoris | 876 m | 1.17 m | 2.81 m | 410 (179) | 44 | 56% |
| Mare Serenitatis | 1,401 m | 1.48 m | 2.81 m | 545 (409) | 125 | 56% |
| Mare Tranquillitatis | 354 m | 1.21 m | 2.91 m | 368 (232) | 99 | 52% |
| Mare Fecunditatis | 532 m | 1.64 m | 3.73 m | 543 (248) | 204 | 34% |
| Mare Nubium | 767 m | 1.68 m | 3.19 m | 615 (425) | 220 | 56% |
| Mare Humorum | 930 m | 1.32 m | 2.45 m | 526 (217) | 47 | 57% |

The reporting cycle lasts 708.7 hours. In Mare Imbrium Hs stays at or above
1 m for 690 of them, 582 in one stretch.

## Swell and wind sea

The directional spectra at the reference points, every three hours, divide
into wind sea and swell by wave age: a component is wind sea where 1.2 × 28 u*
cos(θ − θ_wind) exceeds its deep-water phase speed, and swell otherwise. Swell
carries 34–64% of the wave energy through the cycle and dominates 26–67% of
the time. Its mean height is 0.75–1.32 m, 3.1 m at most in Imbrium, and its
median mean period 13–16 s. Fecunditatis, in the east, has the least swell
(34%).

## Coasts

The time-mean wave power travelling toward the coast, from SWAN's energy
transport at 1-degree sea nodes beside land, has a median of 87 W/m, a 90th
percentile of 319 W/m and a maximum of 731 W/m. The strongest coasts face the
equatorward seas:

| Coast | Land lies to the | Mean power toward it |
|---|---|---:|
| Western Oceanus Procellarum, by Russell | south-west | 731 W/m |
| North-western Oceanus Procellarum, by Harding | west | 654 W/m |
| Western Oceanus Procellarum, near 12° N | south-west | 595 W/m |
| South-western Mare Imbrium, by Delisle | south | 584 W/m |
| Southern Mare Imbrium, at the eastern end of Montes Carpatus | south | 583 W/m |
| Eastern Mare Fecunditatis, by Atwood | east | 575 W/m |
| The Tranquillitatis plateau's north-western shore, by Jansen | south-east | 520 W/m |
| The Aristarchus plateau's north-eastern shore, by Rupes Toscanelli | south-west | 506 W/m |

Each coast is named by the nearest gazetteer landmark: a sea, mountain, scarp or
valley, or a crater at least 20 km across. The landmarks lie 20–105 km from
their nodes; the coast near 12° N lies 143 km from Cardanus, the nearest.
The [coastal figure](../../visualization/waves/results/nearside_coasts.png)
draws every coastal node as a dot coloured and sized by its power. At 1-degree
spacing these are the sea-scale exposures of 30 km stretches of coast;
headlands and bays inside them need the nested grids of the
[shore study](shore.md).

## A month at four shores

Four coastal nodes show how the waves and the tide meet at a shore through the
month ([figure](../../visualization/waves/results/nearside_shores.png)): the
most exposed coast, in western Oceanus Procellarum by Russell; southern Mare
Imbrium at the eastern end of Montes Carpatus; the largest tide on the most
exposed tenth of the coast, in southern Mare Nubium by Pitatus; and a coast of
Mare Nectaris by Fracastorius near the median exposure.

| Shore | Depth | Mean power toward it | Largest Hs (day) | Hours with Hs ≥ 1 m | Tide range this month (typical) |
|---|---:|---:|---:|---:|---:|
| Western Oceanus Procellarum | 78 m | 731 W/m | 3.4 m (23.1) | 442 | 2.8 m (3.2 m) |
| Southern Mare Imbrium | 209 m | 583 W/m | 4.0 m (13.7) | 661 | 3.0 m (3.4 m) |
| Southern Mare Nubium | 101 m | 363 W/m | 3.5 m (13.7) | 403 | 7.7 m (6.3 m) |
| Mare Nectaris | 151 m | 111 W/m | 2.2 m (16.6) | 274 | 7.6 m (6.2 m) |

Storms reach the shores on different days. Southern Imbrium and southern Nubium
peak together on day 13.7, while western Procellarum stays mostly below 1 m
from day 7 to day 17.6 and then holds 2.7–3.4 m from day 20.5 to day 27, with
mean periods up to 20 s: swell about 100 m from crest to crest. The tide rises
and falls once through the month, by about 3 m on the northern shores and
7.6–7.7 m on the southern ones. Southern Nubium's high water, on day 14.6, comes
as western Procellarum's water stands lowest.

The GCM's month has no date, so the tide is a real one. The GCM's Sun turns
evenly once per synodic month on a clock that matches the GCM's own sunlight to
0.0001° over the year before the wave run. In each of the 234 months of
2026–2045 the real Sun stands once over the longitude where the GCM's Sun
stands at the middle of the wave month; the tide is the month of 7 February to
8 March 2038, within 23% of every shore's typical range, and its Sun stays
within 0.5° of the GCM's. No month is typical at every shore: across the
matched months the ranges at southern Imbrium and southern Nubium correlate at
−0.67. Air pressure moves the level by under 0.2 m and wind setup by about
1 cm.

## The calculation

The sea runs on the atlas's 1-degree nodes, 7,053 of them wet, on a 160 × 104
grid from 85.875° W to 73.125° E and 32.125° S to 70.875° N, in spherical
coordinates with the lunar radius. The grid crosses the 0-degree meridian, so
the GCM's longitudes extend periodically across it; across the meridian the
stress steps between neighbouring nodes match those elsewhere (99th
percentiles 0.0060 and 0.0083 Pa). 257 open-sea GCM cells feed the sea; ten
coastal nodes take their nearest open-sea cell, at most 234 km away. The
stress becomes SWAN's input wind through SWAN's own drag law every 15 minutes.

The physics and spectral grid are the Smythii cycle's: SWAN 41.51 in the
coupled-air OpenMP build (air 1.404 kg/m³), Komen growth and whitecapping, the
gravity-scaled AGROW knee, depth breaking with index 0.73, and 36 directions by
48 frequency intervals from 0.00497 to 0.497 Hz. The timestep is 300 s, the
author's choice on 2026-10-03 for a run of 3.3 hours on eight threads against
about 7.5 hours at Smythii's 150 s. Repeating the second cycle's stormiest
three days at 150 s from the main run's hotfile (hours 1,224–1,296 after a
day's adjustment, around the peak basin-mean stress at hour 1,242) changes
metre-scale seas by at most 3.6% in Hs and 9.0% in mean period, 0.12% and 0.25%
at the 95th percentile. Seas of a few decimetres change by up to 86% in Hs and
145% in mean period at single nodes, with 95th percentiles of 0.20% and 0.53%,
and the time above 1 m changes by at most 2.6 hours of the 72 at any node. One
node shallower than SWAN's 5 cm minimum depth stays dry and leaves the
statistics.

The record runs in 48-hour segments, each ending in a SWAN hotfile that starts
the next. A hotstart reproduces an unbroken run bit for bit: hours 4–6 of a
6-hour control match exactly after a restart at hour 3. The eight reference
points sit at the deep-water nodes (at least 300 m) nearest the IAU centres of
their maria, 5–18 km from them, except Tranquillitatis's, 120 km away where its
basin turns shallow. Of 1,888 three-hourly reference spectra with resolved Hs
of at least 0.1 m, 24 put more than 1% of their variance in the top two
frequency intervals; all 24 have resolved Hs under 0.49 m.

## Reproduction

The run uses the coupled-air OpenMP build on the research drive (see the
[synthesis](README.md#builds-and-storage)):

```sh
EXE=/media/projectspace/terluna-research/wave-runs/swan_builds/lunar-seas-swan/optimized_openmp/coupled_air/swan.exe
python -m climate.waves.nearside audit --executable $EXE
OPENBLAS_NUM_THREADS=1 nice -n 10 python -m climate.waves.nearside run --executable $EXE --step 300
python -m climate.waves.nearside_analysis assemble
python -m climate.waves.nearside_analysis check --executable $EXE
python -m climate.waves.nearside_analysis summarize
climate/gcm/.venv/bin/python -m climate.waves.shore_month check-sun   # netCDF4: the GCM's Sun clock
python -m climate.waves.shore_month
python -m visualization.waves.nearside
```

`run` resumes after any stop: completed segments are reused once their hashes
match, and an interrupted one is moved aside and run again. The runs and the
assembled hourly cube live under `research/runs/waves/nearside/`.
