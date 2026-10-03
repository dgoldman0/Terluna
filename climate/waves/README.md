# Wind waves at lunar gravity

The aim is the wave climate of the lunar seas: how winds build waves, how
energy travels across a basin, where terrain shelters or concentrates it,
and how conditions change through the lunar day. Coastal breaking, marine
travel and wave-driven mixing draw on that common picture.

SWAN (Simulating WAves Nearshore) runs here as a spectral wave model. The
experiments progress from prescribed winds to recovered atmospheric surface
stress. The first experiment compares Earth and lunar gravity over
the same flat, 100 km strip of water, with winds drawn from the corrected CM1
equatorial ring. It establishes a runnable wave calculation and measures its
numerical sensitivities.

The [results](results/waves.json) contain significant wave height, peak period
and mean period at four fetches and six forcing durations. Significant height
is the spectral quantity `Hs = 4 sqrt(m0)`; mean period is `Tm01 = m0/m1`.
The [CSV](results/waves.csv) is a tabular companion. Read it with the JSON's
evidence statement, units, source hashes and numerical checks. The
[renderer](../../visualization/waves/) plots the 100 km time series.

The follow-up below adds recovered equatorial wind vectors and surface density,
early-growth refinement, a two-dimensional Smythii–Marginis pilot and simple
coastal slopes. The original strip calculation and its outputs are retained.

The [coastal-resolution follow-up](coastal.md) nests eastern Smythii at
7.6, 3.8 and 1.9 km spacing and checks directional refinement. It also recovers
all three corrected CM1 ring histories and ten GCM years of regional wind
means, with each product's height and sampling interval retained.

The [weather-sequence study](weather.md) captures 241 global atmospheric
snapshots from a 30-Earth-day continuation, then drives a selected week of basin waves
with their surface stress. It checks the timestep and the influence of
48 and 96 additional hours of preceding weather.

The [continuous-cycle study](cycle.md) extends the global atmospheric record
to 60 Earth days and follows the waves through two solar cycles. It excludes
the first 29.53059 days and measures the second: Hs reaches 3.375 m, and the
offshore reference point has Hs of at least 1 m for 309.6 hours, with a
98.1-hour continuous episode. The study records occurrence, persistence,
period and direction, together with timestep, history and evolving spectral
checks. Small-wave frequency coverage and finer basin geometry remain open.

The [shore exposure study](shore.md) follows selected second-cycle spectra
through eastern Smythii's finer geography. It adds native 118 m terrain,
a regional propagation grid and local coastal refinements around the
shoreline spur. Full directional spectra measure incoming and outgoing
power at bathymetry-selected shore stations. Spatial, boundary and spectral
comparisons accompany the conditional exposure results.

## Findings from the first run

At 100 km fetch after 48 hours of the imposed wind, the modified build gives:

| Wind at 10 m | Earth Hs | Lunar Hs | Earth mean period | Lunar mean period |
|---|---:|---:|---:|---:|
| 2.15 m/s, water median | 0.099 m | 0.472 m | 1.24 s | 6.39 s |
| 4.00 m/s, water p90 | 0.296 m | 1.357 m | 2.12 s | 10.69 s |
| 5.84 m/s, water p99 | 0.573 m | 2.412 m | 2.92 s | 13.74 s |

Lower gravity produces heights 4.2–4.8 times as large and substantially longer
periods under these matched conditions. The Earth controls use the same
lunar-derived wind magnitudes, so this establishes a gravity effect within
the chosen model. Whether the actual lunar seas are rougher than Earth's
seas depends on their different wind histories, fetches and coupling physics.

The numerical similarity test agrees within 0.303% in height and 0.326% in
mean period with the modified source, against 5.44% and 3.91% for stock SWAN.
Earth outputs are exactly preserved. At the three main 48-hour endpoints,
the modification changes lunar height by at most 0.11%.

The lunar 4.00 m/s case has these sensitivity results. Endpoint differences
are signed; maxima are absolute relative differences across all hourly
outputs from 1–48 hours and fetches from 10–100 km.

| Change | Final Hs change | Final mean-period change | Largest earlier/local Hs change | Largest earlier/local mean-period change |
|---|---:|---:|---:|---:|
| Time step 300 → 150 s | 0.00% | 0.00% | 33.2% | 17.0% |
| Cell width 1000 → 500 m | −0.295% | −0.187% | 2.52% | 1.86% |
| Directions 36 → 72 | +0.147% | +0.094% | 0.479% | 0.372% |
| Low frequency halved, high frequency doubled | +0.074% | +0.094% | 0.906% | 5.61% |

**The early growth curves remain sensitive to numerical choices.** At
100 km, halving the time step raises the 1-hour height from 0.132 to 0.176 m;
the difference falls to 4.3% after 6 hours and disappears at the printed
precision after 48 hours. The AGROW modification also changes early lunar
heights substantially. These curves show conditional development, with the
strongest numerical support at the tested late endpoint. The follow-up adds
frequency-bin convergence tests.
Comparable refinement at the other wind strengths remains unperformed.

All 22 accepted runs have zero diagnosed depth-breaking fraction after
initialization. Every final spectrum passes the edge checks; the largest
high-edge fraction is 0.683% of resolved variance. The numerical runs take
about 7.5 minutes total on one thread in this environment, excluding builds
and the rejected disabled-breaking diagnostics.

## Experiment

The corrected equatorial ring supplies 10 m wind magnitudes of 2.15, 4.00 and
5.84 m/s: its water median, 90th percentile and 99th percentile. The source
holds 237 snapshots, three hours apart, across days 29.5–59. Each wave case
imposes one magnitude, a fixed direction and uniform wind for up to 48 hours.
Each wind percentile selects the magnitude of an imposed episode. Estimating
event durations and wave occurrence probabilities requires wind histories.

| Choice | Value or treatment |
|---|---|
| Geometry | Cartesian, one-dimensional propagation strip; directional wave spectrum |
| Fetch samples | 10, 25, 50 and 100 km |
| Duration samples | 1, 3, 6, 12, 24 and 48 hours, starting from calm water |
| Depth | Uniform 1000 m, an imposed deep-water reference |
| Gravity | Earth standard 9.80665 m/s²; Moon `GM/R²`, 1.62421887656 m/s², from shared constants |
| Densities | Air 1.28 kg/m³, SWAN's fixed default; water 1025 kg/m³, a seawater scenario |
| Physics | Komen wind input and whitecapping; Wu drag; default discrete interaction approximation for four-wave interactions |
| Initial wind source | AGROW linear growth from `INIT ZERO`; stock and gravity-scaled frequency-knee variants |
| Grid | 1 km, 300 s; 36 directions; 49 logarithmically spaced frequencies |
| Frequency band | 0.03–3 Hz at Earth gravity, scaled in proportion to gravity |
| Propagation | First-order backward-space/backward-time scheme |
| Boundaries | Incoming swell set to zero; outgoing waves can leave |
| Omitted processes | Currents, varying bathymetry, bottom friction, ice, vegetation and capillary-wave physics |

Depth breaking is explicitly enabled with `BREAKING CONSTANT 1.0 0.73` and
its diagnosed fraction must remain zero. Source review and a debug build
found an uninitialized-variable path when breaking was disabled. The supported
configuration avoids that path without adding a second source patch; the
[audit](audit.md) records the diagnosis and the checks.

The [input inventory](inputs.json) records the atmospheric and geographic
products. The first experiment uses a prescribed strip with a stated fetch and
depth. The follow-up adds the gridded atlas and full equatorial wind histories.

## Source change and evidence

The builder verifies the original SWAN 41.51 archive and creates independent
stock and modified serial builds. The modification scales AGROW's arbitrary
1 Hz frequency knee with gravity. The source audit explains why this restores
dimensional similarity and which other Earth-specific options are avoided.
The patch is a parameterization sensitivity; lunar air–water coupling and
growth/dissipation laws still need physical validation.

The two comparisons serve different purposes. The numerical similarity test
reduces gravity exactly sixfold while enlarging length, depth and duration
sixfold and reducing frequencies sixfold. It checks corresponding
dimensionless states. The main experiment keeps the physical basin and wind
duration fixed while changing gravity to the actual lunar value. Its wave
heights depend on the combined effects of gravity, fetch and duration.

Identical Earth controls must repeat exactly, and the modified build must
preserve stock Earth profiles. The modified similarity test requires height
and mean-period agreement within 1% wherever the reference height is at least
1 cm. Peak period remains a frequency-bin diagnostic on an approximately 10%
spaced grid.

For the lunar 4.00 m/s case, separate runs halve the time step, halve the cell
width, double the number of directions and widen the frequency band. The
wider band keeps approximately the same frequency spacing to test the band
limits. Acceptance requires changes below 5% in
the final 100 km height and mean period. Maximum and RMS differences across
the earlier times and inner fetches are retained too. These checks apply to
that case and endpoint; earlier outputs retain the reported timestep sensitivity.
Final spectra must have interior peaks and little resolved variance at the
band edges. The source audit explains how the spectral output floor is
handled and why integrated printed spectra differ from SWAN's height
diagnostic, which includes a high-frequency tail.

## Follow-up: recovered inputs, early growth and coastal waves

The local workspace holds the atlas grid and CM1 snapshots.
[wave_forcing.py](../crm/wave_forcing.py)
recovers both 10 m wind components, surface pressure, 2 m temperature and
vapour mixing ratio from the same 237 equatorial snapshots. It reads only the
required byte ranges and records **selected-byte hashes**, offsets and
source-file sizes. The resulting line product is kept in ignored
`research/runs/waves/inputs/`; its description and hash are
in [results/forcing.json](results/forcing.json).

Using CM1's moist ideal-gas closure, the air over the equatorial water columns
between 75°E and 103°E has mean density **1.40390 kg/m³**, ranging from 1.38029
to 1.41668 kg/m³. The coupled SWAN build uses that mean, with a separate,
recorded change to `PWIND(16)`. Density stays fixed through each run. CM1's
closure constants are recorded in `shared/constants.json` for interpreting
its thermodynamics. The wind record covers one equatorial line every three
hours; basin-wide winds and shorter gusts need additional forcing data.

### Early growth

[early_growth.py](early_growth.py) holds air density at the original 1.28 kg/m³
to isolate numerical refinement. At 100 km fetch, with the
4.001 m/s wind, the one-hour significant height changes as follows:

| Timestep | One-hour Hs |
|---|---:|
| 300 s | 0.1320 m |
| 75 s | 0.2242 m |
| 15 s | 0.2800 m |
| 5 s | 0.2950 m |
| 2 s | 0.3037 m |
| 1 s | 0.3089 m |

The [full record](results/early_growth.json) retains every comparison,
including failed criteria. From **one to two hours**, the 2→1 s change is at
most 1.69% in height and 1.15% in mean period over 10–100 km, excluding heights
below 1 cm. That window was examined after the broader check failed: over
20 minutes–2 hours, the largest height change is 4.63%, but the mean-period
change still reaches **12.23%**. The earliest transients are therefore still
unresolved. Predicting the response to a twenty-minute gust requires further
timestep refinement.

Frequency spacing is a separate check. At a 15 s timestep, increasing from
96 to 192 frequency intervals changes height by at most 0.451% and mean period
by 0.761% over the tested 20-minute–6-hour samples. These time and frequency
checks apply to this wind and numerical source choice; other wind strengths
and physical wind-wave onset remain open.

### Basin and assumed coastal slopes

[pilot.py](pilot.py) uses atlas body 874, Smythii–Marginis, at the 28% water
scenario's common level, −1,654.2 m above the geoid. It uses spherical
coordinates with the lunar radius explicitly supplied. Native atlas nodes
are sampled every 1° and 0.5° (about 30.3 and 15.2 km at the equator), retaining
their water depths. The finite ground elevations keep dry nodes dry; other
disconnected waters inside the rectangular crop are excluded. Small channels,
shoals and coastlines require finer grids. The bathymetry represents flooded
rock; beach and sediment profiles require separate assumptions.

Two imposed episodes start from calm water: a 4.001 m/s wind toward the east
for twelve hours, and a wind that turns from east at hour 5 to north at hour 6
and then holds north. SWAN interpolates the vector components during that
hour, so the speed briefly falls to 2.83 m/s. The turning episode therefore
changes both direction and speed.
The original coarse runs give basin maxima of 1.222 m and 0.990 m at hour 12.
Estimating the frequency of those wave heights requires wind histories.

The [pilot product](results/pilot.json) records the episodes, selected
offshore node, maps, source/build hashes and numerical differences. The
75→30 s timestep comparison uses 96 frequency intervals. A separate 96→48
interval comparison uses the coarse grid, and **both grids in the spatial
comparison use 48 intervals**. This separates spectral and spatial changes
while bounding run cost. The initially attempted fine-grid, 96-interval run
was interrupted after its projected cost exceeded the run cap; its partial
output is retained for diagnosis under ignored
`research/runs/waves/pilot_interrupted/`. The comparisons use the five completed
basin runs.

At hour 12, over matched wet nodes with Hs at least 0.1 m in both runs:

| Comparison | Largest relative Hs difference | Largest relative mean-period difference |
|---|---:|---:|
| 75 → 30 s timestep, 96 frequency intervals | 0.55% | 0.39% |
| 96 → 48 frequency intervals, coarse grid | 0.12% | 0.09% |
| 1° → 0.5° grid, 48 frequency intervals | **41.3%** | **22.7%** |

Differences use the second run as the denominator. The spatial comparison
changes local heights by up to 0.35 m and **fails the 5% bound**. Selecting surf
locations or ranking sheltered coasts requires finer coastal grids. At the selected
offshore node, final Hs and mean period each change by less than 0.002% under
that grid refinement, supporting its use in the slope experiment. The timestep
difference over the full hourly record reaches 16.5% in Hs. Early growth and
geographic directional refinement remain open.

The final spectrum at 89.125°E, 1.125°S from the 30 s run drives three
stationary slopes, from 50 m depth to 0.25 m. `SPEC1D` retains frequency
variance and each frequency's mean direction and spread. SWAN reconstructs
the directional distributions from those moments. Resolving multiple
directional peaks requires a full two-dimensional spectrum. Bottom friction, triads,
whitecapping, quadruplets and currents are omitted in these slope cases to
isolate propagation and depth breaking.

With the assumed breaking index 0.73, at the first sampled point where
SWAN's breaking fraction `Qb` reaches 1%:

| Assumed slope (400 cells) | Depth | Significant height |
|---|---:|---:|
| 1:20 | 2.24 m | 1.100 m |
| 1:50 | 2.12 m | 1.041 m |
| 1:100 | 1.99 m | 0.958 m |

The 800-cell 1:50 case reaches that diagnostic at 2.18 m with Hs 1.057 m.
At matched depths, the 400→800-cell comparison differs by at most 1.51% in
Hs (0.0047 m), below the stated 5% bound in this single refinement check.
All six stationary runs satisfy their
99%-of-wet-nodes stopping criterion, and all five offshore spectra pass the
stated spectral-edge checks.
Changing its breaking index to 0.60 or 0.90, on the 400-cell grid, moves the
diagnostic to 2.49 or 1.74 m respectively. The coefficients are assumed
sensitivity values; lunar calibration remains open. Here `Qb ≥ 0.01` marks
the chosen onset of breaking, and Hs describes the wave spectrum. The
calculations support metre-scale surf under the stated wind and slope assumptions.

[pilot_checks.py](pilot_checks.py) writes the [numerical audit](results/pilot_checks.json):
spectral edges, stationary solver stopping, matched slope-grid errors and
the separate basin resolution differences. Independent executable tests compare
metre and lunar-degree grids and check conservation of linear-wave energy flux
before breaking. Physical validation, run-up, wave-driven currents and individual
breaking waves remain future work. SWASH needs a source audit and pilot run.

## Reproduction

Use Python with NumPy, and `gfortran`, `gcc`, `make` and Perl on the executable
search path. The source download and builds go to an explicitly selected new
directory. SWAN source headers specify GNU GPL version 3 or later and credit
Delft University of Technology; keep those terms when distributing SWAN or a
modified executable. This repository records the download hash and patch,
without vendoring the upstream source or binary.

```sh
python -m climate.waves.build --download --build-root /tmp/terluna-swan
# Or replace --download with --archive /path/to/swan4151.tar.gz.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.run \
  --build-root /tmp/terluna-swan
TERLUNA_SWAN_EXECUTABLE=/tmp/terluna-swan/patched_agrow/swan.exe \
  OPENBLAS_NUM_THREADS=1 python -m pytest climate/tests/test_waves.py
python visualization/waves/plot.py
```

Builds and runs use one thread. Each run has a three-minute wall/CPU limit
and a 2 GiB address-space limit. Completed runs are reused only after input,
executable and output hashes match and the outputs pass the current parser.
An incomplete directory is retained for inspection and blocks replacement.
Raw runs live under ignored `research/runs/waves/`; selected results and
provenance live here. The actual repository checks are in [checks.json](checks.json).

On the author's machine, `research/runs/waves` links to
`/media/projectspace/terluna-research/wave-runs` on the large drive. Bulk inputs,
copied atmospheric continuations and raw wave outputs use that directory;
commands and recorded paths continue to use `research/runs/waves/`. Keep the
drive mounted for wave runs. Code and compact research reports stay in the
repository.

For the follow-up, the raw equatorial CM1 case and atlas grid must already be
available. The exporter stops with an error if either input is missing.
With the build above:

```sh
OPENBLAS_NUM_THREADS=1 python -m climate.crm.wave_forcing
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.early_growth --build-root /tmp/terluna-swan
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m climate.waves.pilot --build-root /tmp/terluna-swan
python -m climate.waves.pilot_checks
python visualization/waves/pilot.py
```

The pilot adds an isolated `coupled_air` build and leaves the stock and
AGROW-only executables intact. Its individual runs are serial, capped at
900 seconds and 2 GiB; the early-growth runs have a 600-second cap. The
follow-up's actual repository verification is in
[pilot_verification.json](pilot_verification.json), distinct from the
first experiment's check record and from the numerical audit.

## Next sea calculations

The [shore exposure study](shore.md) now supplies native 118 m terrain and
selected second-cycle coastal spectra. Extend those nests through evolving
weather, widen the frequency band for short waves, and refine the remaining
shoreline outliers. Explicit cross-shore profiles can then supply the
wave-resolving coastal calculations. The global three-hourly atmospheric
record supports continuous wave histories and comparisons between seas.
Extend the [second-cycle calculation](cycle.md) to further cycles and basins,
and widen the frequency band where the smaller waves reach its upper edge.
Short wind transients, parent directional resolution, and physical lunar wind
input and breaking remain open. Local CM1 boxes can address identified gaps
in gust fronts, coastal flow or storm winds. An individual-wave SWASH
calculation needs a source/gravity audit, resolved coastal geometry and
input-spectrum checks before interpreting run-up or breaking crests.

Temperature gradients, mixing and light penetration need their own ocean
calculations. Start with vertical heat/mixing and spectral underwater-light
columns, using the atmosphere's surface fluxes and illumination products.
Water absorption and scattering, dissolved material, particles and biological
optics supply the inputs to those calculations. The conventional 1%
photosynthetically active radiation depth provides an optical diagnostic;
biological compensation depth also depends on an organism's energy balance.
The abstract of
[Wu et al. (2021)](https://doi.org/10.1029/2020JC016874), read through the
[NOAA repository](https://repository.library.noaa.gov/view/noaa/55092), explains
the distinction between this optical threshold and biological compensation
depth; full-paper review remains pending. These ocean and optical calculations
remain unperformed.
