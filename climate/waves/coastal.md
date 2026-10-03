# Eastern Smythii: coastal resolution and recovered winds

This follow-up refines the eastern shore of the flooded Smythii–Marginis basin
and recovers the available atmospheric records over a wider area. The coast
uses the same prescribed twelve-hour, 4.001 m/s eastward wind as the basin
pilot. Wind recovery produces separate geographic datasets for subsequent
forcing experiments.

## Terrain and wave grids

[geography/coastal_grid.py](../../geography/coastal_grid.py) reads the
hash-checked LOLA 16-pixel/degree grid and evaluates the degree-200 GRAIL geoid,
static Earth tide and rotation directly on a small sector. The bounds are
91.03125–95.03125°E and 1.03125°S–2.96875°N. The grid has 65 × 65 nodes,
including 2,994 wet nodes connected to the western boundary. It keeps the
28% atlas scenario's water level, −1,654.2 m above the geoid. Isolated water
components are held dry. Heights represent flooded rock; sediment and beach
evolution require additional modelling.

Eleven of those wet nodes are shallower than 10 m; the shallowest is 0.199 m.
The kilometre grid samples the coastal approaches. Resolving the narrow
breaking zone requires finer cross-shore profiles and a stated beach geometry.

[coastal.py](coastal.py) supplies this same terrain to every coastal run.
The computational grids sample every fourth, second and native terrain node:
0.25°, 0.125° and 0.0625°, approximately 7.6, 3.8 and 1.9 km at the equator.
Holding the bottom input fixed isolates how the wave grid resolves that
terrain. The 0.5° basin parent retains the previous pilot's coarser atlas
bathymetry.

The parent writes full frequency–direction spectra at 32 points around the
coastal perimeter every 15 minutes. SWAN's
[NGRID](https://swanmodel.sourceforge.io/online_doc/swanuse/node31.html),
[NESTOUT](https://swanmodel.sourceforge.io/online_doc/swanuse/node32.html) and
[BOUNDNEST1](https://swanmodel.sourceforge.io/online_doc/swanuse/node27.html)
commands define the transfer. This retains directional structure through the
boundary exchange. The parent has 36 directions and 48 frequency intervals;
the child tests 36 and 72 directions at 3.8 km. Refining the child's directions
tests its propagation and source integration while keeping the parent's
directional resolution fixed.

All cases start from calm water, use a 75 s timestep, and retain the pilot's
CM1-derived air density of 1.40390 kg/m³. Wind input uses Komen/Wu and the
gravity-scaled AGROW knee. Depth breaking uses the assumed index 0.73.
Currents, bottom friction, triads, diffraction and wave-induced water-level changes remain
outside this experiment. The [original audit](audit.md) records the lunar
source changes and the remaining calibration questions.

The comparison uses final common wet nodes with Hs ≥ 0.1 m in both runs,
with a two-coarse-cell buffer inside the nesting boundary. Relative differences
use the finer run as denominator. Both maximum and 95th-percentile differences
are retained; the acceptance bound is 5% for the maximum height and mean-period
differences. This is a numerical check for the stated episode.

The parent basin table exactly matches the preceding 0.5° pilot. Every western
wet boundary sample has a parent spectrum. At 93.03125°E, 1.03125°S on the
southern boundary, the fine terrain exposes water where the parent supplies
`NODATA`. SWAN reads that sample as zero incident energy before boundary
interpolation (`swanmain.ftn`, the `SWN` boundary branch in the pinned 41.51
source). The output audit retains this mismatch for interpreting nearby waves.

## Results

The [complete product](results/coastal.json) and
[numerical audit](results/coastal_checks.json) retain all four coastal cases,
their final maps and each comparison. At hour 12:

| Comparison | Common interior nodes | Maximum Hs difference | 95th-percentile Hs difference | Maximum mean-period difference | 5% criterion |
|---|---:|---:|---:|---:|---|
| 7.6 → 3.8 km | 128 | 16.22% (0.131 m) | 2.88% | 9.49% | Fail |
| 3.8 → 1.9 km | 623 | 118.01% (0.623 m) | 7.50% | 67.34% | Fail |
| 36 → 72 child directions, 3.8 km | 623 | 0.57% (0.0057 m) | 0.45% | 0.31% | Pass |

The grid pairs use their own common nodes and two-coarse-cell boundary buffer.
The 1.9 km grid reveals a narrow shoreline spur. At 93.40625°E, 2.34375°N,
in 48.3 m of water, Hs changes from **1.150 m to 0.528 m**. Immediately west,
the native node at 93.34375°E stands 147.6 m above the water. The 3.8 km grid
samples wet nodes on either side and skips that dry node. This spur and its
sheltered water are a concrete target for the next local refinement.

The broad exposed-water maximum stays close to **1.223 m** across the three
36-direction grids. Individual coastal values remain sensitive to the grid,
including a 7.5% difference at the 95th percentile in the finest comparison.
The spatial check therefore remains open for selecting surf locations.
The directional check passes for the child at 3.8 km; parent directional
refinement and a coastal timestep check remain separate tasks.

At the finest resolution, one node reaches the chosen 1% breaking diagnostic:
94.34375°E, 0.78125°N, depth 0.199 m, Hs 0.114 m and breaking fraction 4.42%.
The sparse shallow nodes leave the metre-scale breaking profiles to finer
transects. The previous assumed-slope results remain conditional on their
offshore spectrum, slope and breaking coefficient.

The final boundary spectra pass the band-edge checks: the largest fraction
of resolved variance in the upper two frequency intervals is 0.0194%, with
interior spectral peaks. The independent flat-depth nesting test passes its
2% height-and-period bound. The [coastal figure](../../visualization/waves/results/coastal.png)
shows the changing shoreline and wave fields; the
[verification record](coastal_verification.json) separates numerical outcomes
from repository checks and the remaining physical work.

## Recovered wind coverage

| Product | Spatial coverage | Saved timing and level | Use |
|---|---|---|---|
| Equatorial CM1 ring | 1,816 columns around the equator; 70 water columns in the Smythii–Marginis rectangle | 237 instantaneous snapshots, three hours apart, days 29.5–59; 10 m winds | Wind histories along the sea's equatorial crossing |
| Two tilted CM1 rings | Great circles tilted 70°, ascending at 45°E and 135°E | 237 snapshots per ring, the same cadence and nominal model-day interval; 10 m winds | Higher-latitude coverage; both tracks pass outside this basin rectangle |
| Corrected GCM, years 20–29 | Native T21 subset, 8 latitude × 7 longitude nodes over 70–110°E, 15°S–30°N | 1,200 means, each over 143 half-hour steps (71.5 Earth hours); lowest model level at sigma 0.9833 | Regional circulation and geographic context |

[climate/crm/wave_coverage.py](../crm/wave_coverage.py) rotates CM1's along-ring
and left-of-track wind components into local east and north using the declared
great circle. It preserves each ring's independent dynamical realisation.
An independent Cartesian tangent calculation tests the rotation. Density uses
the same moist-gas closure as the earlier equatorial export. Selected-byte
hashes identify the six original surface fields at every snapshot.

[climate/gcm/wave_coverage.py](../gcm/wave_coverage.py) reads the archived
`A28_dim5_moon` means and hashes all ten source files. Its saved timestamps
retain the model's step convention. The lowest layer follows surface pressure
and atmospheric thickness. These are model-level averages; basin wave forcing
requires surface vectors sampled frequently enough to preserve wind episodes.
The inspected annual files omit `uas` and `vas`, and the run's snapshot
directory is empty.

The [ring record](results/ring_coverage.json) and
[GCM record](results/gcm_coverage.json) identify the ignored NPZ products, their
schemas, units, input hashes and reading rules. The
[coverage figure](../../visualization/waves/results/wind_coverage.png) shows
the tracks, the final GCM averaging window and one CM1 water-column history.

The next atmospheric input is a basin-wide history of 10 m eastward and
northward winds, hourly or three-hourly, with surface pressure, temperature
and moisture. A dedicated restart must save instantaneous fields and establish
a surface-wind diagnostic consistent with lunar gravity and the model's
surface-layer treatment. That history will support event durations and surf
occurrence calculations. The recovered ring tracks and GCM means provide the
available context for choosing and checking the run.

## Reproduction

The original pilot supplies the `coupled_air` SWAN build. The LOLA/GRAIL files,
28% atlas, three corrected CM1 runs and GCM years 20–29 must be restored first.
The GCM exporter uses the environment containing `netCDF4`.

```sh
OPENBLAS_NUM_THREADS=1 python -m geography.coastal_grid
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 python -m climate.waves.coastal \
  --build-root /tmp/terluna-lunar-seas-swan
python -m climate.waves.coastal_checks
OPENBLAS_NUM_THREADS=1 python -m climate.crm.wave_coverage
OPENBLAS_NUM_THREADS=1 climate/gcm/.venv/bin/python -m climate.gcm.wave_coverage
MPLCONFIGDIR=/tmp/terluna-lunar-seas-mpl-cache python visualization/waves/coastal.py
```

Runs are serial with a 2 GiB address-space bound. The native coastal grid has
a 1,500 s limit; the other cases have a 900 s limit. An initial native-grid
attempt was stopped when its progress projected beyond the original 900 s
allowance. Its partial files and interruption record remain under ignored
`research/runs/waves/coastal_interrupted/`; the comparison uses complete runs.
The cache verifies executable, input and output hashes before reuse. Raw runs,
terrain, recovered arrays and figures remain under ignored directories;
tracked results record their provenance. The independent nesting test compares
an embedded grid with its parent over uniform depth, requiring height and
mean-period agreement within 2% at interior nodes. Terrain tests check native
pixel orientation and exclusion of an isolated lake.
