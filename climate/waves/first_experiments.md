# The first wave experiments: strip, early growth and basin pilot

These experiments of 2 October 2026 set up SWAN at lunar gravity and measured
its numerical behaviour before any real weather drove it. The
[source audit](audit.md) holds the code review behind them; the
[synthesis](README.md) places them in the whole picture.

## Earth and lunar gravity over a 100 km strip

The first experiment compares Earth and lunar gravity over the same flat
100 km strip of water. The winds come from the corrected CM1 equatorial ring:
2.15, 4.00 and 5.84 m/s at 10 m, its water median, 90th and 99th percentiles.
Each case holds one wind steady in one direction for up to 48 hours from calm
water. The [results](results/waves.json) give significant wave height, peak
period and mean period at four fetches and six durations (a
[CSV](results/waves.csv) companion holds the same rows).

At 100 km fetch after 48 hours:

| Wind at 10 m | Earth Hs | Lunar Hs | Earth mean period | Lunar mean period |
|---|---:|---:|---:|---:|
| 2.15 m/s, water median | 0.099 m | 0.472 m | 1.24 s | 6.39 s |
| 4.00 m/s, water p90 | 0.296 m | 1.357 m | 2.12 s | 10.69 s |
| 5.84 m/s, water p99 | 0.573 m | 2.412 m | 2.92 s | 13.74 s |

Under the same wind, fetch and duration, lunar gravity gives waves 4.2–4.8
times as high, with periods about five times as long. Both worlds here use
SWAN's Wu drag law. At lunar gravity Charnock's law gives drag 0.85, 1.11 and
1.34 times Wu's at these three winds ([review](review.md)), so the lunar
heights at the two stronger winds sit 5–34% low and the median-wind case
8–15% high.

The numerical similarity test reduces gravity sixfold while enlarging length,
depth and duration sixfold and reducing frequencies sixfold, so the
dimensionless states should coincide. They agree within 0.303% in height and
0.326% in mean period with the modified source, against 5.44% and 3.91% for
stock SWAN. Earth outputs repeat exactly, and the modification changes the
lunar heights at the three main endpoints by at most 0.11%.

For the lunar 4.00 m/s case, endpoint differences are signed; the maxima are
absolute relative differences across hourly outputs from 1–48 hours and fetches
from 10–100 km:

| Change | Final Hs change | Final mean-period change | Largest earlier/local Hs change | Largest earlier/local mean-period change |
|---|---:|---:|---:|---:|
| Time step 300 → 150 s | 0.00% | 0.00% | 33.2% | 17.0% |
| Cell width 1000 → 500 m | −0.295% | −0.187% | 2.52% | 1.86% |
| Directions 36 → 72 | +0.147% | +0.094% | 0.479% | 0.372% |
| Low frequency halved, high frequency doubled | +0.074% | +0.094% | 0.906% | 5.61% |

The 48-hour endpoint is converged. Early growth depends on the timestep: at
100 km, halving it raises the 1-hour height from 0.132 to 0.176 m; the
difference falls to 4.3% after 6 hours and vanishes at the printed precision by
48 hours. The AGROW modification also changes early lunar heights. All 22
accepted runs have zero depth-breaking fraction after initialization, every
final spectrum passes the edge checks (the largest high-edge fraction is
0.683% of resolved variance), and the runs take about 7.5 minutes in all on
one thread.

### Design

| Choice | Value or treatment |
|---|---|
| Geometry | Cartesian, one-dimensional propagation strip; directional wave spectrum |
| Fetch samples | 10, 25, 50 and 100 km |
| Duration samples | 1, 3, 6, 12, 24 and 48 hours, starting from calm water |
| Depth | Uniform 1000 m, deep water |
| Gravity | Earth standard 9.80665 m/s²; Moon `GM/R²`, 1.62421887656 m/s², from shared constants |
| Densities | Air 1.28 kg/m³, SWAN's fixed default; water 1025 kg/m³ |
| Physics | Komen wind input and whitecapping; Wu drag; the discrete interaction approximation for four-wave interactions |
| Initial wind source | AGROW linear growth from `INIT ZERO`; stock and gravity-scaled frequency-knee variants |
| Grid | 1 km, 300 s; 36 directions; 49 logarithmically spaced frequencies |
| Frequency band | 0.03–3 Hz at Earth gravity, scaled in proportion to gravity |
| Propagation | First-order backward-space/backward-time scheme |
| Boundaries | No incoming swell; outgoing waves leave freely |
| Left out | Currents, varying bathymetry, bottom friction, ice, vegetation and capillary waves |

Depth breaking runs with `BREAKING CONSTANT 1.0 0.73`, and its diagnosed
fraction stays zero in deep water. With breaking switched off, SWAN 41.51 reads
an uninitialised variable; the [audit](audit.md) records the diagnosis.

## The source changes

The builder verifies the original SWAN 41.51 archive and makes independent
stock and modified serial builds. One change scales AGROW's 1 Hz frequency knee
with gravity, which restores dimensional similarity. A second replaces SWAN's
fixed air density, 1.28 kg/m³, with the mean moist air density over the
equatorial water columns between 75° E and 103° E in CM1's corrected ring,
1.40390 kg/m³ (range 1.38029–1.41668 kg/m³). [wave_forcing.py](../crm/wave_forcing.py)
recovers both 10 m wind components, surface pressure, 2 m temperature and vapour
mixing ratio from the ring's 237 equatorial snapshots, reading only the
required byte ranges and recording their hashes; its line product lives in
`research/runs/waves/inputs/`, described in [results/forcing.json](results/forcing.json).

## Early growth

[early_growth.py](early_growth.py) holds air density at 1.28 kg/m³ to isolate
the numerics. At 100 km fetch with the 4.001 m/s wind, the one-hour height
depends on the timestep:

| Timestep | One-hour Hs |
|---|---:|
| 300 s | 0.1320 m |
| 75 s | 0.2242 m |
| 15 s | 0.2800 m |
| 5 s | 0.2950 m |
| 2 s | 0.3037 m |
| 1 s | 0.3089 m |

From one to two hours, the 2 → 1 s change is at most 1.69% in height and 1.15%
in mean period over 10–100 km for heights of at least 1 cm. Over 20 minutes to
2 hours the height change reaches 4.63% and the mean-period change 12.23%: the
first minutes of growth need still smaller steps. At a 15 s timestep, doubling
the frequency intervals from 96 to 192 changes height by at most 0.451% and
mean period by 0.761% over 20 minutes to 6 hours. The
[full record](results/early_growth.json) keeps every comparison.

## Basin pilot and assumed slopes

[pilot.py](pilot.py) runs atlas body 874, Smythii–Marginis, at the 28% scenario's
water level, −1,654.2 m above the geoid, in spherical coordinates with the
lunar radius. Native atlas nodes every 1° and 0.5° (30.3 and 15.2 km at the
equator) keep their water depths; finite ground elevations keep dry nodes dry,
and other waters inside the rectangular crop are held dry.

Two episodes start from calm water: a 4.001 m/s wind toward the east for twelve
hours, and a wind that turns from east to north between hours 5 and 6 (SWAN's
vector interpolation briefly lowers the speed to 2.83 m/s). The basin maxima at
hour 12 are 1.222 m and 0.990 m. At hour 12, over matched wet nodes with Hs of
at least 0.1 m in both runs:

| Comparison | Largest relative Hs difference | Largest relative mean-period difference |
|---|---:|---:|
| 75 → 30 s timestep, 96 frequency intervals | 0.55% | 0.39% |
| 96 → 48 frequency intervals, coarse grid | 0.12% | 0.09% |
| 1° → 0.5° grid, 48 frequency intervals | **41.3%** | **22.7%** |

Grid refinement changes local heights by up to 0.35 m, while the offshore node
at 89.125° E, 1.125° S changes by less than 0.002%. The
[pilot product](results/pilot.json) records the episodes, maps, hashes and
differences; a fine-grid run stopped at its time cap is kept under
`research/runs/waves/pilot_interrupted/`.

The final offshore spectrum drives three stationary slopes from 50 m depth to
0.25 m, as `SPEC1D` frequency variance with each frequency's mean direction and
spread; bottom friction, triads, whitecapping, quadruplets and currents are
switched off to isolate shoaling and depth breaking. With breaking index 0.73,
at the first point where SWAN's breaking fraction `Qb` reaches 1%:

| Assumed slope (400 cells) | Depth | Significant height |
|---|---:|---:|
| 1:20 | 2.24 m | 1.100 m |
| 1:50 | 2.12 m | 1.041 m |
| 1:100 | 1.99 m | 0.958 m |

The 800-cell 1:50 case reaches that point at 2.18 m with Hs 1.057 m; at matched
depths the two grids differ by at most 1.51% in Hs. Breaking indices of 0.60
and 0.90 move the point to 2.49 and 1.74 m. [pilot_checks.py](pilot_checks.py)
writes the [numerical audit](results/pilot_checks.json): spectral edges,
stationary stopping, matched slope-grid errors and the basin resolution
differences. Tests compare metre and lunar-degree grids and check the
conservation of linear-wave energy flux before breaking. The
[SWASH study](runup.md) follows the waves through breaking to the shoreline.

## Reproduction

With the build described in the [synthesis](README.md#builds-and-storage):

```sh
python -m climate.waves.build --download --build-root <new build directory>
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.run --build-root <build directory>
TERLUNA_SWAN_EXECUTABLE=<build directory>/patched_agrow/swan.exe python -m pytest climate/tests/test_waves.py
python visualization/waves/plot.py
OPENBLAS_NUM_THREADS=1 python -m climate.crm.wave_forcing
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.early_growth --build-root <build directory>
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.pilot --build-root <build directory>
python -m climate.waves.pilot_checks
python visualization/waves/pilot.py
```

Builds and runs use one thread, with three-minute (strip), ten-minute (early
growth) and fifteen-minute (pilot) caps and a 2 GiB address-space limit.
Completed runs are reused after their input, executable and output hashes
match. The strip's actual repository checks are in [checks.json](checks.json)
and the pilot's in [pilot_verification.json](pilot_verification.json).
