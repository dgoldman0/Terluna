# The nearside sea through a lunar cycle

The nearside sea holds half of the Open Moon's standing water: 5.29 million km²
from Oceanus Procellarum and Mare Frigoris to Fecunditatis and Nectaris, mean
depth 648 m, along the Earth-facing coasts. This study follows its waves through
two lunar cycles of the recovered weather, with the physics of the
[Smythii–Marginis cycle](cycle.md), and measures the second cycle.
[nearside.py](nearside.py) runs the waves, [nearside_analysis.py](nearside_analysis.py)
measures them, and the [product](results/nearside.json) holds the statistics,
maps, reference points, swell and coastal power.

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
a tenth reach 3.9 m. Where Hs is at least 0.5 m, mean periods run from 8.3 to
15.1 s (10th to 90th percentile, median 11.5 s), and metre-scale seas peak at
16.6 s at the median and 20.1 s at the 90th percentile. These waves are 34 to
130 m long and travel at 3 to 6 m/s.

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
| Western Oceanus Procellarum, by Montes Agricola | south-west | 731 W/m |
| Western Oceanus Procellarum, by Mons Rümker | west | 654 W/m |
| Western Oceanus Procellarum, near 12° N | south-west | 595 W/m |
| South-western Mare Imbrium, by Mons Delisle | south | 584 W/m |
| Southern Mare Imbrium, along Montes Carpatus | south | 583 W/m |
| Eastern Mare Fecunditatis | east | 575 W/m |
| The Tranquillitatis plateau's north-western shore | south-east | 520 W/m |
| The Aristarchus plateau's north-eastern shore, by Mons Herodotus | south-west | 506 W/m |

The [coastal figure](../../visualization/waves/results/nearside_coasts.png)
colours the atlas shoreline by its nearest coastal node. At 1-degree spacing
these are the sea-scale exposures of 30 km stretches of coast; headlands and
bays inside them need the nested grids of the [shore study](shore.md).

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
about 7.5 hours at Smythii's 150 s. Repeating the second cycle's stormiest three days at 150 s from the main run's hotfile (hours 1,224–1,296 after a day's adjustment, around the peak basin-mean stress at hour 1,242) changes metre-scale seas by at most 3.6% in Hs and 9.0% in mean period, 0.12% and 0.25% at the 95th percentile. Seas of a few decimetres change by up to 86% in Hs and 145% in mean period at single nodes, with 95th percentiles of 0.20% and 0.53%, and the time above 1 m changes by at most 2.6 hours of the 72 at any node. One node shallower than
SWAN's 5 cm minimum depth stays dry and leaves the statistics.

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
python -m visualization.waves.nearside
```

`run` resumes after any stop: completed segments are reused once their hashes
match, and an interrupted one is moved aside and run again. The runs and the
assembled hourly cube live under `research/runs/waves/nearside/`.
