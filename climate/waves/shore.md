# Wave energy at eastern Smythii's shore

This study follows wave energy around the shoreline spur near 93.4°E,
2.34°N. It connects the established second-cycle basin sea to native 118 m
terrain, then compares the two faces of the headland and a nearby southern
shore. The [numerical product](results/shore.json) records the individual
cases, spectra, energy fluxes and refinement checks.

The headland creates a strong contrast over about five kilometres. In the
east-northeast sea, the native 118 m wave grid gives the following values:

| Station | Significant height | Mean period | Directional power |
|---|---:|---:|---:|
| Western face | 1.842 m | 12.60 s | 600 W/m |
| Eastern face | 0.228 m | 2.89 s | 0.021 W/m |
| Southern shore | 1.334 m | 12.35 s | 210 W/m |
| Offshore reference | 1.997 m | 12.86 s | 639 W/m net transport magnitude |

The first three powers cross the local depth contour toward shallower
water. The offshore value gives the magnitude of the net transport vector.

The long-period incident sea reaches the western face and the southern
shore. Behind the headland, short local waves carry most of their energy
away from the eastern face. The incoming-power calculation reveals a much
stronger contrast than wave height alone. The tiny eastern incoming flux
belongs to this model's treatment of the shadow: diffraction and shoreline
reflection remain further physical calculations.

The short waves in that shadow retain numerical uncertainty. On the 474 m
grid, extending the upper frequency limit changes the eastern station's
mean period from 3.40 to 3.71 s; tightening solver accuracy gives 3.90 s.
The 118 m case has 8.88% of the eastern station's resolved variance in the
highest two frequency intervals, exceeding the 1% edge criterion. Heights
around a few tenths of a metre and the very small incoming power describe
the present sheltered-water result; its detailed short-wave spectrum needs
further convergence work.

The sheltered side changes as the incident sea turns. These comparisons
use the same 237 m wave grid and the same stations for all three phases:

| Day in the second cycle | Basin waves travel toward | Western Hs | Eastern Hs | Western incoming power | Eastern incoming power |
|---:|---|---:|---:|---:|---:|
| 5.97 | East-northeast | 1.787 m | 0.272 m | 570 W/m | 0.039 W/m |
| 12.97 | North-northeast | 1.025 m | 0.333 m | 187 W/m | 5.25 W/m |
| 14.97 | West-southwest | 0.250 m | 0.771 m | 8.35 W/m | 52.9 W/m |

The southern station's incoming power changes from 219 to 24.5 to
18.7 W/m across these phases. Geography and direction together determine
which shore receives the waves. The regional boundary and numerical
refinements were run for the first, stronger eastward case. The later
phases establish additional conditional examples; their own refinement
comparisons remain further work.

The [exposure map](../../visualization/waves/results/shore.png) shows native
terrain, wave height and energy transport. The
[three-phase map](../../visualization/waves/results/shore_phases.png) uses a
common power scale to show how the exposed side changes. The spectra's
independently integrated transport agrees with SWAN's vector output within
0.093% across all 49 reported station records.

## Terrain and water level

[shore_terrain.py](../../geography/shore_terrain.py) restores a pinned row
slice of the NASA LRO LOLA 256-pixel/degree height raster. The
[source record](../../geography/shore_inputs.json) identifies the archive
version, label, byte range, exact hashes and source credit. Its selected
0–6°N rows occupy 141.6 MB; three independent byte-range requests agree with
the corresponding downloaded bytes. The archive's
[label](https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_256_00n_90n_000_180.lbl)
describes the altimetry-derived raster and interpolation across coverage
gaps. The 118 m figure specifies raster spacing.

The water surface retains the atlas's 28% water scenario and its level of
−1,654.2 m above the geoid. The degree-200 GRAIL geoid is evaluated directly
at the terrain nodes. The inherited approximation between the LOLA
mean-Earth frame and the GRAIL principal-axis frame remains in use. Heights
represent flooded rock; sediment, erosion and beach formation remain further
geographic models.

The main square covers about 30.3 km on each side, from 92.751953125°E,
1.751953125°N. Its 257 × 257 native nodes include 54,136 wet nodes. A
four-connected water mask retains the sea connected to the western edge.
Every spatial refinement reads this same native bottom grid, so changing
wave-grid spacing isolates the propagation grid's response to the terrain.
An expanded 45.5 km square checks the placement of the coastal boundary.

## Weather phases and boundary spectra

The first lunar solar cycle supplies the basin spin-up. Each selected
boundary spectrum comes from the second cycle of the existing
[continuous wave history](cycle.md):

| Archive hour | Earth days into the second solar cycle | Selection |
|---:|---:|---|
| 852 | 5.969 | Strong eastward wave-energy proxy over the eastern basin |
| 1020 | 12.969 | Near the cycle's maximum basin significant height |
| 1068 | 14.969 | A later phase for a changed incident sea |

The first selection maximizes the regional mean of
`Hs² × Tm01 × max(cos(direction), 0)` over 90–94°E, 0–4°N in the second
cycle. The hourly maximum is at hour 853; hour 852 is the nearest saved
three-hourly boundary time. This selection is independent of the eventual
nearshore wave values.

The basin replay adds full frequency-direction spectra and transport
diagnostics. Its 3,103,551 retained wave-output values agree exactly with
the established calculation through hour 1080. A second replay retains the
full directional field at every resolved basin water node for the three
selected phases.

The 1° basin grid places land across some water visible in the finer
terrain. Direct nesting into the 30 km square leaves 25 wet boundary samples
without a source spectrum in the initial 64-point check. The study therefore
adds a wider eastern Smythii sector, about 182 × 174 km, before the small
coastal grids. It uses the same native terrain and the selected atmospheric
wind field.

[shore_region.py](shore_region.py) transfers spectra explicitly among wet
source nodes. Inverse-square weights use up to four neighbouring spectra
within two source-grid spacings. Each wet target must lie within 1.5 spacings
of an available source node; exceeding that distance stops the calculation.
Dry target points retain `NODATA`. The product records the distances and
weights' support. The regional outer boundary's largest nearest-source
distance is 1.008 basin cells, about 30.6 km. The local grids repeat this
procedure using the much finer regional spectra. This extension across
different resolved shorelines is a boundary assumption, with separate
regional-grid, boundary-position and interpolation controls.

Each coastal case holds its sampled wind and incident spectrum fixed while
the stationary wave solution settles. It measures spatial exposure to that
weather phase. Coastal timing, lag, duration and occurrence require an
evolving regional calculation with the preceding spectra.

## Energy reaching a shore

SWAN's default `TRANSP` output is a vector of wave-variance transport in
m³/s. Multiplying by the configured water density, 1025 kg/m³, and lunar
gravity, 1.62421887656 m/s², gives W/m: energy per second crossing a metre of
wavefront. The [SWAN output manual](https://swanmodel.sourceforge.io/online_doc/swanuse/node32.html)
documents its variance and energy conventions. Source inspection and an
executable control verify the conversion, the direction convention and the
deep-water group speed.

The station spectra provide a more specific coastal measure. For a unit
normal **n** toward shallower water, the incoming power is

`ρg ∫∫ Cg(f,h) E(f,θ) max(n · uθ, 0) dθ df`.

The outgoing contribution is integrated separately. Their difference is
the net shoreward flux. This separates wave energy arriving at a shore from
energy travelling along it or away from it. The normal comes from the local
gradient of the native flooded-rock depth. Frequency integration covers the
resolved spectrum; SWAN's height and mean-period diagnostics also include
its diagnostic high-frequency tail.

Comparison stations were chosen from bathymetry before running the coastal
waves. They share exact nodes on the 474, 237 and 118 m grids:

| Station | Longitude, latitude | Water depth | Purpose |
|---|---|---:|---|
| Western face | 93.251953°E, 2.361328°N | 22.38 m | Exposure on the western side of the spur |
| Eastern face | 93.408203°E, 2.330078°N | 49.04 m | Exposure behind the spur during eastward seas |
| Southern shore | 92.986328°E, 1.830078°N | 13.55 m | A nearby coast facing another direction |
| Offshore | 92.845703°E, 2.345703°N | 1270.77 m | Local incoming-sea reference |

The stations have different depths and local orientations; those quantities
accompany every reported flux. A coast's relative exposure can change as the
incident sea turns.

## Numerical scope

The inherited calculation uses Komen wind input and whitecapping, Wu drag,
the gravity-scaled AGROW modification, and depth breaking with index 0.73.
The atmospheric surface-stress conversion retains its fixed-density SWAN
build and spatial density correction. Currents, bottom friction,
diffraction, shoreline reflection, sediment evolution and swash are omitted.
Lunar growth and breaking coefficients still need physical calibration.

Stationary stopping requires the specified accuracy at 99% of wet nodes.
The output retains the final attained percentage even when a case reaches
its iteration cap. Spatial comparisons use shared wet nodes at least 2 km
inside the smaller domain and Hs ≥ 0.1 m in both cases. The recorded exposure
tolerances are 5% at the 95th percentile and 15% maximum for height; the
corresponding vector-flux tolerances are 10% and 30%. Depth-band results and
the largest local changes remain visible alongside the aggregate check.
Flux-vector errors use the refined vector magnitude with a 1 W/m floor;
the period diagnostic uses a 1 s floor. These floors bound relative errors
as the reference signal approaches zero.

The initial regional iteration remained oscillatory, satisfying the
criterion at 48.12% of wet nodes after 300 iterations. Frequency-dependent
under-relaxation with coefficient 0.01 reached 99.14% after 20 iterations.
This follows the [SWAN numerical method](https://swanmodel.sourceforge.io/online_doc/swantech/node48.html),
which reduces iterative updates while preserving the stationary equations.
The coastal comparisons use the converged regional product. The failed
attempt remains in the raw record. Tightening regional accuracy fivefold
passes the exposure tolerances; the same local accuracy check retains
outliers in small waves.

| Comparison | Height difference: 95th percentile / maximum | Transport-vector difference: 95th percentile / maximum | Exposure bounds |
|---|---:|---:|---|
| Wave spacing 948 → 474 m | 6.24% / 46.68% | 15.87% / 190.45% | Fail |
| Wave spacing 474 → 237 m | 3.57% / 35.08% | 10.44% / 115.43% | Fail |
| Wave spacing 237 → 118 m | 2.84% / 39.09% | 7.89% / 149.00% | Fail |
| Directions 36 → 72 at 474 m | 0.58% / 2.43% | 1.65% / 7.93% | Pass |
| Coastal domain 30 → 45 km | 0.96% / 4.87% | 2.30% / 10.29% | Pass |
| Regional spacing 1.9 → 0.95 km | 0.24% / 1.39% | 0.63% / 3.59% | Pass |
| Outer interpolation: four neighbours → nearest | 0.86% / 1.77% | 2.30% / 3.56% | Pass |
| Regional solver accuracy 0.005 → 0.001 | 2.55% / 6.32% | 6.22% / 13.46% | Pass |
| Local solver accuracy 0.005 → 0.001 | 4.82% / 18.21% | 6.47% / 45.24% | Fail |
| Upper frequency 0.497 → 0.973 Hz | 0.47% / 7.34% | 0.39% / 13.11% | Pass |

The broad exposed/shadowed pattern is stable under the tested boundary and
direction controls. Finer wave grids continue to change individual
nearshore values, especially in weak or sheltered water. Between 237 and
118 m, western incoming power changes from 570 to 600 W/m and southern
incoming power from 219 to 210 W/m. The maximum local height difference
across shared nodes still reaches 39%, so detailed shoreline rankings
retain a spatial convergence question.

The 118 m terrain resolves the headland and substantially more shallow
water: 234 active wave nodes are shallower than 10 m, and 45 reach the
chosen `Qb ≥ 0.01` breaking diagnostic in the eastward case. These node
counts describe the discretized calculation. The breaking and swash zones
on steep rock can still lie between
nodes. Beach-scale modelling needs explicit shore profiles and a gravity
audit of the chosen wave-resolving model.

## Reproduction and storage

All bulk files live under `research/runs/waves/shore/`, through the link to
`/media/projectspace/terluna-research/wave-runs`. Existing runs are reused
after input, executable and output hashes match. Incomplete runs remain
available for diagnosis. The basin replays use four threads; coastal runs
use two, with bounded wall time and address space.

Restore the terrain and produce the three geographic grids:

```sh
OPENBLAS_NUM_THREADS=1 python -m geography.shore_terrain --download
OPENBLAS_NUM_THREADS=1 python -m geography.shore_terrain \
  --bounds 92.501953125 1.501953125 1.5 1.5 \
  --output research/runs/waves/shore/expanded_terrain.npz
OPENBLAS_NUM_THREADS=1 python -m geography.shore_terrain \
  --bounds 90.751953125 0.001953125 6 5.75 \
  --output research/runs/waves/shore/regional_terrain.npz
```

Use the recorded optimized coupled-air executable from the
[continuous-cycle study](cycle.md#reproduction). Supply its path through
`--executable` to each simulation command:

```sh
python -m climate.waves.shore parent --executable /path/to/coupled_air/swan.exe
python -m climate.waves.shore fields --executable /path/to/coupled_air/swan.exe
python -m climate.waves.shore_region --hour 852 --stride 16 \
  --executable /path/to/coupled_air/swan.exe
python -m climate.waves.shore snapshot --hour 852 --stride 4 \
  --executable /path/to/coupled_air/swan.exe
python -m climate.waves.shore_checks
python -m visualization.waves.shore
```

The individual products retain the parameters of each completed refinement.
The [verification record](shore_verification.json) separates repository tests,
data-integrity checks, numerical comparisons and remaining physical work.
