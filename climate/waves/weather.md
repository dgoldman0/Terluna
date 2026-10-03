# Waves through a recovered weather sequence

A fresh continuation of the saved GCM restart supplies **241 global snapshots,
three Earth hours apart, spanning 30 Earth days**. A selected seven-day Smythii–Marginis experiment
connects their surface stress to SWAN. This advances the sea-wide questions:
wave growth, propagation, exposure, persistence and changing direction through
the lunar day. Coastal breaking and run-up build on that wave climate.

In the run with four preceding days of weather, **significant height reaches
3.162 m**, with mean period **15.40 s** at that peak. It occurs at 84.125°E,
2.125°S, hour 68 of the selected week. At its widest extent, metre-scale seas
cover **86.0% of the resolved basin area**. The offshore reference point at
89.125°E, 1.125°S reaches 2.236 m; its mean period ranges from 8.78 to 15.62 s
over hours 24–168. These are conditional sea-wide results on the 1° grid.

The preceding weather matters to the early sea state. The two- and four-day
history extensions agree within 2.42% in height and 2.24% in mean period
after the first selected day. Starting the selected week from calm water
changes some early metre-scale heights by 23%. The peak of the main event
stays near 3.16 m in all four runs.

The [runner](weather.py), [checks](weather_checks.py) and
[product](results/weather.json) retain the episode selection, numerical
comparisons and input hashes. The [renderer](../../visualization/waves/weather.py)
displays the computed fields and histories.

## Atmospheric record

[wave_snapshots.py](../gcm/wave_snapshots.py) copies the settled
`A28_dim5_moon` year-29 restart, its executable, physical namelists, surface
maps and stellar spectrum into an isolated wave-study directory. It advances
1,446 half-hour model steps, saving instantaneous output every six steps.
The first saved step is 5; the final one is 1,445. The product retains this
step convention and measures its 30-day span from the first snapshot.

A control uses the same copied inputs with snapshots disabled. Its final
restart and ordinary atmospheric output are **byte-identical** to the snapshot
run. The source archive's input hashes also remain unchanged. Each atmospheric
run used eight MPI ranks for about 18 seconds. The serial decoder uses the
existing ExoPlaSim environment.

The ice and mixed-layer ocean restore their fields from the primary restart.
The optional `restart_dsnow` and `restart_xsnow` override files are absent
from this source archive. The runner carries them into its copy when present.

Recovered fields include eastward and northward winds at the lowest model
level, pressure, temperature, humidity, sea ice, roughness, friction velocity
cubed, and the two surface-stress components. Winds retain their model level,
sigma 0.9833. A separate streaming reader checks all 3,454,976 values of seven
surface fields directly against the binary snapshots. The decoded wind and
surface-stress directions agree to a cosine error below 2.4e-10.

The installed `fluxmod.f90` computes downward surface stress from the momentum
removed from the lowest atmospheric layer. Its diagnostic satisfies
`u_star^2 = |stress| R_d T_surface / p_surface`, with `R_d` supplied by the
planet namelist. The recovered diagnostic and this expression agree within
0.52% over the global record. The surface temperature advances within the
physics step; the snapshot records its later value. This density is the
model's **stress-closure density**. A moist 2 m density requires an additional
surface-layer diagnostic.

## Passing momentum to the sea

The coupling chooses a stress-equivalent SWAN input vector:

`rho_SWAN Cd_Wu(|Ueq|) |Ueq| Ueq = atmospheric surface-stress vector`.

The Wu coefficients come from the pinned SWAN 41.51 source and are recorded
in [shared/constants.json](../../shared/constants.json). SWAN converts its
wind input to friction velocity through that drag relation;
its Komen and AGROW options then determine wave growth.
The [SWAN physics manual](https://swanmodel.sourceforge.io/online_doc/swanuse/node28.html)
describes those options. This is a chosen one-way coupling of the two models.
An independently diagnosed 10 m wind and two-way wave-dependent roughness
remain further physical checks.

The existing SWAN build retains air density 1.40390 kg/m³. The atmospheric
stress-closure density over the mapped sea varies from 1.38853 to 1.40931 kg/m³,
with mean 1.39996 kg/m³ during the selected week. The inversion transfers the
dimensional stress using SWAN's fixed density. Spatial and temporal density
effects within the wave source terms remain a separate sensitivity.

Spatial interpolation uses open-water GCM corners, renormalizing their
bilinear weights. The whole saved record must keep a source cell ice-free.
An empty stencil would use the nearest open-water cell; every wet basin node
in this experiment has an ocean corner. Sixteen atmospheric cells contribute,
with the farthest contributing source about 231 km from its target. The
T21 circulation supplies regional weather; local coastal winds need finer
atmospheric grids.

Stress vectors are interpolated to 15-minute times before conversion to SWAN
input. At the intervening midpoints, SWAN's linear vector interpolation differs
from the intended stress by at most 3.7e-5 Pa; the 99th-percentile relative
error is 0.226% for stress of at least 1e-4 Pa. SWAN's saved wind components
and friction velocity provide a further check inside the wave solver.

## Episode and numerical checks

The selected week covers **days 23–30 after the first snapshot**, including
the month's largest area-weighted basin-mean stress on day 25.75. The
stress-equivalent input has median 2.83 m/s, 90th percentile 4.64 m/s and
maximum 7.22 m/s across sampled wet nodes and forcing times. These are
descriptions of this selected model episode.

The basin uses the existing flooded LOLA/GRAIL atlas, body 874, with a 1°
wave grid, 36 directions and 49 frequencies. Komen growth and whitecapping,
DIA interactions, the gravity-adjusted AGROW term and depth breaking retain
the previous pilot's settings. Ocean currents, changing water level and
bottom friction remain outside this calculation.

The timestep comparison uses 150 and 75 seconds with identical forcing.
A third run begins 48 hours earlier, at a 150-second timestep, to measure
dependence on the initial calm sea. Its initial comparisons changed some
metre-scale heights by over 20% during hours 24–30 of the selected week.
That finding prompted a fourth run with 96 hours of preceding weather.
Comparisons align the same weather times
and locations, require Hs of at least 0.1 m in both runs, and examine hours
24–168 of the selected week. Relative differences use the finer-timestep
or longer-history case as the denominator. A 5% bound applies to both Hs
and mean period. Additional windows expose how the history sensitivity
changes through the week; the original 24-hour comparison remains recorded.
The figure uses the run with 96 hours of preceding weather. The timestep
comparison retains its initially calm seven-day window.

Over hours 24–168 of the selected week:

| Comparison | Largest Hs change | Largest mean-period change | Both below 5% |
|---|---:|---:|---|
| 150 → 75 s timestep | 12.95% | 26.73% | Fails |
| Add 48 h of preceding weather to the calm start | 78.06% | 72.68% | Fails |
| Extend preceding weather from 48 → 96 h | 2.42% | 2.24% | Passes |

For the timestep pair, the 95th-percentile differences are 0.228% in height
and 0.464% in period. Where both runs have Hs of at least 1 m, the largest
changes are 0.446% and 1.760%. The full-field criterion fails in smaller seas:
the largest relative height change is 0.142 → 0.163 m in northern Marginis;
the largest mean-period change is 6.81 → 5.37 s in waves about 0.25 m high
on the southern margin. Timestep and spectral refinement of those conditions
remain further work.

The 48 → 96 h history comparison also retains the initial selected day:
including it gives maxima of 11.77% in height and 7.89% in period. From
selected hour 48 onward, those bounds fall to 0.752% and 1.032%. This measures
the decay of initialization effects for this weather sequence.

The wind files match over every common interval, within their written
precision, so the history and timestep comparisons keep their forcing fixed.
All four final offshore spectra pass the recorded edge checks. Spectral
coverage during weak-wave transients remains a separate check. The full
repository check passes 472 Python tests, 176 JavaScript tests, the layer
check and the provenance and ensemble validators; the
[verification record](weather_verification.json) retains the log and artifact hashes.

The [coastal refinement](coastal.md) has already established large local
geographic sensitivity. That remains relevant to this 1° basin calculation.
Precise coastal heights require further spatial refinement. Wave occurrence
across months and seasons requires representative longer histories and
checks of the atmospheric and wave closures at lunar gravity.

## Reproduction

Raw model files and full hourly wave fields are kept under ignored
`research/runs/waves/weather/`. Existing output directories require inspection
before replacement; completed SWAN cases are reused only after input and
output hashes match.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 climate/gcm/.venv/bin/python -m climate.gcm.wave_snapshots run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 climate/gcm/.venv/bin/python -m climate.gcm.wave_snapshots control
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 climate/gcm/.venv/bin/python -m climate.gcm.wave_snapshots export

python -m climate.waves.weather --executable /tmp/terluna-lunar-seas-swan/coupled_air/swan.exe --name basin_dt150 --step 150
python -m climate.waves.weather --executable /tmp/terluna-lunar-seas-swan/coupled_air/swan.exe --name basin_dt75 --step 75
python -m climate.waves.weather --executable /tmp/terluna-lunar-seas-swan/coupled_air/swan.exe --name basin_lead48_dt150 --step 150 --hours 216
python -m climate.waves.weather --executable /tmp/terluna-lunar-seas-swan/coupled_air/swan.exe --name basin_lead96_dt150 --step 150 --hours 264
python -m climate.waves.weather_checks
python -m visualization.waves.weather
```

The SWAN executable is the existing coupled build documented in the
[wave README](README.md). Each wave run uses one CPU thread. The
48-hour history extension follows the coarser-timestep run; the 96-hour
extension follows the finer-timestep run. Weather cases use at most two CPU
threads together. Repository checks also include short SWAN control runs.
