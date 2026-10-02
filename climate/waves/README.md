# Wind waves at lunar gravity

SWAN (Simulating WAves Nearshore) runs here as a spectral wave model driven by
prescribed winds. The first experiment compares Earth and lunar gravity over
the same flat, 100 km strip of water, with winds drawn from the corrected CM1
equatorial ring. It establishes a runnable wave calculation and measures its
numerical sensitivities. A new three-dimensional CM1 box is not required for
this comparison.

The [results](results/waves.json) contain significant wave height, peak period
and mean period at four fetches and six forcing durations. Significant height
is the spectral quantity `Hs = 4 sqrt(m0)`; mean period is `Tm01 = m0/m1`.
The [CSV](results/waves.csv) is a tabular companion. Read it with the JSON's
evidence statement, units, source hashes and numerical checks. The
[renderer](../../visualization/waves/) plots the 100 km time series.

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
strongest numerical support at the tested late endpoint. Frequency-bin
convergence and comparable refinement tests at the other wind strengths
remain unperformed.

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
These wind labels do not assign percentiles to the resulting waves, and the
ring summary does not establish that these winds persist for 48 hours.

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
| Boundaries | No prescribed incoming swell; outgoing waves can leave |
| Other processes | No currents, varying bathymetry, bottom friction, ice, vegetation or capillary-wave physics |

Depth breaking is explicitly enabled with `BREAKING CONSTANT 1.0 0.73` and
its diagnosed fraction must remain zero. Source review and a debug build
found an uninitialized-variable path when breaking was disabled. The supported
configuration avoids that path without adding a second source patch; the
[audit](audit.md) records the diagnosis and the checks.

The [input inventory](inputs.json) records the available atmospheric and
geographic products. The atlas summary gives sea areas and depth summaries,
but the local checkout lacks its gridded bathymetry and the rings' full wind
vector histories. The strip therefore represents a stated fetch and depth,
rather than the shoreline of a named lunar sea.

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
duration fixed while changing gravity to the actual lunar value. Finite
fetch and duration mean that its wave heights need not follow an exact
inverse-gravity ratio.

Identical Earth controls must repeat exactly, and the modified build must
preserve stock Earth profiles. The modified similarity test requires height
and mean-period agreement within 1% wherever the reference height is at least
1 cm. Peak period remains a frequency-bin diagnostic on an approximately 10%
spaced grid.

For the lunar 4.00 m/s case, separate runs halve the time step, halve the cell
width, double the number of directions and widen the frequency band. The
wider band keeps approximately the same frequency spacing; it does not test
frequency-resolution convergence. Acceptance requires changes below 5% in
the final 100 km height and mean period. Maximum and RMS differences across
the earlier times and inner fetches are retained too. These checks apply to
that case and endpoint; they do not establish convergence of every output.
Final spectra must have interior peaks and little resolved variance at the
band edges. The source audit explains how the spectral output floor is
handled and why integrated printed spectra differ from SWAN's height
diagnostic, which includes a high-frequency tail.

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

## Next sea calculations

Restore the existing gridded atlas and wind-vector histories to drive a
bounded two-dimensional basin case with changing wind direction and duration.
Audit the geographic radius and shallow-water physics before using lunar
coastlines or interpreting surf. Local CM1 boxes become useful if the basin
needs unresolved gust fronts, coastal flow or storm winds; they are an
additional forcing source, not a prerequisite for spectral wave modeling.

Temperature gradients, mixing and light penetration need their own ocean
calculations. Start with vertical heat/mixing and spectral underwater-light
columns, using the atmosphere's surface fluxes and illumination products.
Water absorption and scattering, dissolved material, particles and biological
optics must be specified before a light-depth result is meaningful. The
conventional 1% photosynthetically active radiation depth is a useful
euphotic-zone diagnostic; it does not define a universal boundary of all
biologically usable light. The abstract of
[Wu et al. (2021)](https://doi.org/10.1029/2020JC016874), read through the
[NOAA repository](https://repository.library.noaa.gov/view/noaa/55092), explains
the distinction between this optical threshold and biological compensation
depth; full-paper review remains pending. These ocean and optical calculations
remain unperformed.
