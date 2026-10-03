# Waves over a complete lunar cycle

This study follows Smythii–Marginis continuously through two solar cycles.
The first **29.53059 Earth days** initialize the wave field. Statistics use
only the second cycle, days **29.53059–59.06118** after the first atmospheric
snapshot. Across that reporting cycle, significant height reaches **3.375 m**,
with a **15.84-second mean period** at the peak. The offshore reference point
has Hs of at least 1 m for **309.6 hours**, including a continuous **98.1-hour**
episode. These are conditional results for one modeled cycle.

The [product](results/cycle.json) records height, exposure, persistence and
mean propagation direction from [cycle.py](cycle.py). The
[numerical audit](results/cycle_checks.json), produced by
[cycle_checks.py](cycle_checks.py), checks evolving spectra and matched
experiments. The
[renderer](../../visualization/waves/cycle.py) displays the second cycle's
heights and occurrence across the sea.

## The sea through this cycle

The largest height occurs at **83.125°E, 13.875°N**, in 345.5 m of water,
13.01 Earth days into the reporting cycle. Its peak period is 18.28 seconds.
The area and time mean of local Hs is **0.880 m**. At the greatest extent,
Hs of at least 1 m covers **84.7%** of the resolved basin area.

The reference point at 89.125°E, 1.125°S reaches Hs **2.119 m**. Its mean
period varies from **5.62 to 14.45 seconds**. The table averages basin
coverage over the full reporting interval and counts durations at that point:

| Hs threshold | Mean basin area at/above threshold | Reference hours | Longest reference stretch |
|---|---:|---:|---:|
| 0.5 m | 81.47% | 648.4 h | 482.7 h, continuing at the record's end |
| 1 m | 34.33% | 309.6 h | 98.1 h |
| 2 m | 1.31% | 7.1 h | 7.1 h |
| 3 m | 0.048% | 0 h | 0 h |

Seven stretches reach 1 m at the reference point; one is already underway
at the reporting boundary. The 98.1-hour stretch lies fully within the
record. Larger seas occupy smaller areas and shorter intervals, while
half-metre conditions are widespread through much of this modeled cycle.

## Atmospheric history and the two clocks

An isolated continuation of the corrected `A28_dim5_moon` year-29 restart
supplies **481 global snapshots**, three Earth hours apart, spanning **60
Earth days**. The atmosphere, ice and mixed-layer ocean carry their saved
states forward. The snapshot run and the control produce byte-identical
final restart and ordinary atmospheric-output files.

The new archive reproduces all **6,416,962 overlapping exported values**
from the earlier 30-day record exactly. An independent binary reader checks
**6,895,616 surface values** across the extended archive. These include the
surface stresses that drive the waves.

The first wave cycle begins at the first saved atmospheric snapshot, 2.5
hours after the restart under the GCM's step convention. The shared synodic
month defines the reporting boundaries. The solar period derived from the
GCM's rotation and orbit differs by 0.31 seconds across the two cycles.
The model's internal calendar day is about 29.8 Earth days; the solar clock
governs this analysis. Cycle phase here is elapsed time from the first
snapshot; local sunrise depends on longitude and the saved solar phase.

The CM1 histories retain their existing second-day selection, model days
29.5–59.0. This experiment uses the separate global GCM continuation and an
explicit first cycle of wave spin-up.

## Wave calculation and statistics

The calculation retains the [weather study's](weather.md) stress coupling,
lunar gravity, density, growth and dissipation settings. SWAN runs from calm
water for 1,419 hours on the 1° basin grid, with 36 directions and 49
frequencies. The 150-second timestep and hourly basin output carry the same
numerical definitions as the earlier week. Spectra at three wet nodes are
saved every three hours: the offshore reference and the two locations where
the earlier timestep test found its largest weak-wave differences.

The reporting interval lasts **708.73416 hours**. Linear interpolation places
its endpoints between hourly outputs. Threshold crossings are integrated on
the piecewise-linear height histories, giving hours above 0.5, 1, 2 and 3 m
and the longest continuous stretch within the reporting interval. Episodes
that reach either reporting boundary carry a censoring flag. A height that
only touches a threshold at one instant contributes zero duration.

Spatial fractions use cosine-latitude weights over resolved wet grid nodes.
The mean height is the time and area average of local significant wave
heights. A separate RMS height records the square-mean quantity associated
with wave energy. Propagation directions use an Hs-squared-weighted vector
mean, together with a resultant that measures the temporal consistency of
the hourly mean heading.
SWAN's Cartesian convention points toward propagation, counterclockwise from
east ([manual](https://swanmodel.sourceforge.io/online_doc/swanuse/node7.html)).
The product's forcing section retains the full input history, including
spin-up; its wave statistics and plotted interval use the second cycle.

One reporting cycle establishes occurrence within that simulation. Comparing
further cycles and seas will establish how these patterns vary. Coastal
geography, lunar calibration of wave growth and dissipation, currents,
changing water levels and beach run-up remain further work.

## Numerical comparisons

The strongest basin-mean stress during the reporting cycle occurs at archive
hour 1,011. Two additional runs cover hours 891–1,155, with 150- and 75-second
timesteps. Their comparison uses hours 1,011–1,155, after five days of wave
evolution in those shorter runs. The 150-second case also supplies a history
comparison against the continuous two-cycle run at exactly the same weather
times. The wind files are checked over every shared 15-minute input time.

| Comparison over archive hours 1,011–1,155 | Largest Hs difference | Largest mean-period difference | 5% criterion |
|---|---:|---:|---|
| 150 to 75 seconds, Hs ≥ 0.1 m in both | 16.29% | 29.37% | Fails |
| Same timestep comparison, Hs ≥ 1 m in both | 0.376% | 0.870% | Passes for this subset |
| Shorter history against continuous run, Hs ≥ 0.1 m in both | 0.440% | 0.913% | Passes |

Across the timestep comparison, the 95th-percentile differences are 0.133%
in height and 0.368% in mean period. Time above 1 m changes by at most
**0.316 hours (19 minutes)** at any node, and the time-averaged basin
coverage above 1 m changes by **0.0092 percentage points**. These checks
cover the selected six-day interval; complete-cycle timestep and spatial
refinements remain further calculations.

Three-hourly spectra are checked at the reference point and two weak-wave
locations. The existing spectral criterion requires less than 0.1% of
resolved variance in the lowest two frequency intervals, less than 1% in
the highest two, and the spectral peak inside the band. These checks apply
to spectra with resolved Hs of at least 0.1 m. The saved offshore endpoint
has its own check; the time series follows changing conditions through the
reporting cycle.

Of **708 spectra** at those three locations during the second cycle,
**129 fail the spectral criterion**. Their resolved Hs ranges from 0.118
to 0.507 m. The largest upper-band fraction is 6.59%. The reference point
accounts for nine failures; the two locations selected for weak-wave
sensitivity account for 70 and 50. Those counts describe the three probes.
The final reference spectrum passes its edge check. A wider frequency band
is the next numerical refinement for these smaller seas.

The optimized build uses GNU Fortran 13.3 with `-O3 -march=native` and
OpenMP. Over the earlier seven-day episode, its maximum changes against the
original executable are **0.512% in Hs** and **1.62% in mean period**, within
the existing 5% comparison bound for Hs of at least 0.1 m. Optimized serial
and two-thread outputs agree exactly for that week; two- and four-thread
outputs agree exactly across a separate 48-hour control. The long run uses
four threads and took 29.7 minutes. The shorter numerical cases use two
threads and took 8.8 and 17.2 minutes.

The atmospheric decoder retains the installed postprocessor's record parser
and transforms. Its local collector concatenates each variable once, reducing
repeated copying of large arrays. It matches the upstream collector across
90 variables and 2,809,946 values from eight complete snapshots; the full
export's overlap check covers the earlier 30-day record.

[cycle_verification.json](cycle_verification.json) records the actual checks,
source and artifact hashes, build comparisons and remaining work. Numerical
agreement and code checks establish their stated computational properties;
lunar physical calibration remains open.

## Reproduction

Raw runs live in ignored `research/runs/waves/cycle/`. Existing run folders
are retained and their hashes checked before reuse. Atmospheric source
archives, earlier wave runs and the original executable retain their bytes.

Capture the atmospheric record and its control:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.gcm.wave_snapshots run \
  --steps 2886 --output research/runs/waves/cycle
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.gcm.wave_snapshots control \
  --steps 2886 --output research/runs/waves/cycle
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 climate/gcm/.venv/bin/python \
  -m climate.gcm.wave_snapshots export --output research/runs/waves/cycle
```

Build in a fresh directory from the pinned archive and the existing density
record. `build_coupled(new_root, density, optimize=True, openmp=True)` in
[build.py](build.py) reads `new_root/swan4151.tar.gz`; its manifest records the
patches, compiler command and executable hash. The recorded build for this
run is under `/tmp/terluna-lunar-seas-swan/optimized_openmp/`.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 nice -n 10 python -m climate.waves.weather \
  --executable /tmp/terluna-lunar-seas-swan/optimized_openmp/coupled_air/swan.exe \
  --atmosphere research/runs/waves/cycle/atmosphere.npz \
  --run-root research/runs/waves/cycle --name basin_dt150 \
  --start-hour 0 --hours 1419 --step 150 --threads 4 --timeout 7200 \
  --diagnostic-spectra
```

For each refinement case, use `--start-hour 891 --hours 264 --threads 2`,
the names `lead96_dt150` and `lead96_dt75`, and their corresponding timesteps.
The analysis commands are:

```sh
python -m climate.waves.cycle
python -m climate.waves.cycle_checks
python -m visualization.waves.cycle
```
