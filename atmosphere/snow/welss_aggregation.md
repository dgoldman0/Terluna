# Author-model Earth aggregation comparison

The declared-source author-model rerun differs from archived growth-only totals by **0.8–2.1%**, while the aggregation case differs by **12.5–16.5%**. Aggregate number-density and mass-flux differences reach approximately **70%** and **68%**. This single-seed reconstruction therefore does not closely reproduce the archived aggregation population.

The calculation extends the Earth reference to the deposition, sedimentation and aggregation experiment in [Welss et al. (2024), Figure 21](https://doi.org/10.1029/2023MS003805). It executes the authors' McSnow code under declared settings and compares its profiles with their archived diagnostic-geometry case, `ob_agg0.nc`. A growth-only run of the same executable provides a control. These are author-solver reruns. The independently integrated growth characteristics remain a separate result in [welss_column.md](welss_column.md); an independently implemented coupled aggregation column remains outstanding.

## Declared source and experiment

The [build wrapper](welss_aggregation_build.py) checks 57 external files against their Git blob identities and SHA256 values in [the manifest](welss_aggregation_sources.json), then compiles revision `1e428d7e901cbc5b91f2040e575f736ec5963a74` with GNU Fortran 13.3, `-O3 -DNDEBUG -cpp`, without OpenMP. The dependency-ordered build uses the source set and release flags specified by the author's CMake configuration. The sole source patch replaces clock-based random seeding with a required integer argument; every entry in the compiler's seed vector receives 20261021. The product retains the exact patch, compiler and executable hashes, commands, namelist and input/output hashes. External source, executable and raw output remain outside Git.

The [run wrapper](welss_aggregation.py) follows the pinned [oblate run script](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/run/1d_obl) and explicitly records inherited process settings. The archive's experiment names distinguish `lwc0` for aggregation and `lwc2` for growth-only, so the runs retain those separate atmospheric water contents even though riming is disabled in both.

| Quantity | Declared value |
|---|---|
| Domain | 0–5,000 m; 250 cells; 20 m layers; 5 m² column area |
| Temperature | 273.15 K at surface; 0.0062 K/m lapse rate |
| Vapour | Source idealized atmosphere and interpolation; 5% over ice above 600 m, water saturation in the lower layer |
| Source | Uniform height 2,700–3,000 m; gamma mass distribution with shape 2 and mean 4.8 × 10⁻¹³ kg |
| Source normalization | 0.01 superparticles m⁻³ s⁻¹ × multiplicity 1,650 = 16.5 real particles m⁻³ s⁻¹; 75 superparticles per 5 s |
| Integration | 10 h; 5 s steps; 600 s profile outputs; seed 20261021 |
| Geometry | `habit0`, `agggeo1`, `snowgeo2`, diagnostic aspect ratio 1 |
| Collisions | `agg2`, geometric encounter area π(D₁+D₂)²/4, Böhm interception and Mitchell temperature-only sticking option 1; added velocity variance 0 |
| Growth control | `agg0`; liquid water 0.2 g/m³ below 600 m; aggregation run uses 0 |
| Other processes | Riming 0, fragmentation 0, multiplication 0, boundary source 0; original melt/shedding/drop-breakup defaults recorded explicitly |

For the selected aspect ratio 1 and `ckt0`, individual projected-area ratios affect drag and interception; the outer encounter area uses the diameters as shown above. Source interception efficiencies are retained even when they exceed 1. The prescribed vapour and temperature fields do not respond to ice growth. This is therefore a comparison of population calculations under a maintained water supply. These fixed-atmosphere results do not test the lunar column's water budget.

## Profiles and numerical scope

The source writer bins particles by `[z, z+20 m)` and reports number density, mass density, and their velocity-weighted moments. `Fn` and `Fm` are cell averages of `Nv` and `Nmv`, rather than face-crossing fluxes. Saved profiles contain total particles, pristine monomers and unrimed aggregates. The final particle CSV independently supplies multiplicity-weighted mass; its sum is checked against the rounded final binned mass profile.

Primary comparisons use the 31 snapshots from 18,000 through 36,000 s inclusive. Four final snapshots provide the alternative specified by public preprocessing. Both physical means, which retain zero populations, and the source preprocessing's positive-only means are saved. These alternatives matter particularly for rare aggregates near their source. The archived finite endpoint bins support unsmoothed comparison. Raw output is rounded to three significant digits; the wrapper averages the parsed values in float64 and does not claim byte-identical NCL processing.

The completed [growth control](results/welss_author_growth.json) and [aggregation run](results/welss_author_aggregation.json) give the following relative L1 differences, summed over populated reference bins below 3,000 m. These compare one declared-source stochastic realization with the archived author profiles.

| Profile | Growth-only | Aggregation |
|---|---:|---:|
| Number density | 2.076% | 12.521% |
| Mass density | 0.828% | 16.482% |
| Number flux | 0.817% | 14.651% |
| Mass flux | 0.862% | 15.251% |

Inside the nucleation belt, number-density differences are 8.004% for growth and 7.755% for aggregation; below the belt they are 0.781% and 13.713%. These regional measures prevent a small total-profile error from concealing a source-layer discrepancy.

For the Figure 21 categories, aggregation-run differences are:

| Category/profile | Physical raw mean | Positive-only mean |
|---|---:|---:|
| Monomer number density | 18.141% | 18.141% |
| Aggregate number density | 70.106% | 70.095% |
| Monomer mass flux | 32.206% | 32.206% |
| Aggregate mass flux | 67.642% | 67.642% |

Summed over the compared bins, aggregate number density is 70.010% higher and aggregate mass flux is 67.641% higher than the archive; monomer number density and mass flux are 14.920% and 32.206% lower. The main additional discrepancy develops below the nucleation belt.

Using four final snapshots instead changes total-profile number-density differences to 2.786% and 12.490%, respectively. The aggregation discrepancy persists under both averaging windows and positive-only processing. No parameter or processing choice was fitted to reduce it; these controls do not identify its physical or numerical cause.

The growth run completed in 761.6 s with 99,255 final superparticles. Its final binned mass differs from the full-precision particle sum by -0.01388%, consistent with the writer's rounding. Saved riming rates, liquid/rimed-ice/graupel category moments and final non-pristine mass components are zero.
The aggregation run completed in 1,088.4 s with 97,244 final superparticles. Its final binned mass differs from the full-precision particle sum by +0.00551%, consistent with the writer's rounding. Saved riming rates, liquid/rimed-ice/graupel category moments and final non-pristine mass components are zero.

All 12 author-build/configuration/profile/result tests pass. Independent review checked source/build identity and the archive-bin/category conventions. These checks establish the declared calculation and comparison, while sampling, discretization and archived-version differences remain unresolved.

The run reports all saved riming rates and liquid, rimed-pristine and graupel category moments, together with final pristine/rime/total particle masses. Those checks constrain unintended process activity at the saved times. They do not establish every internal stage or an exact whole-run external-water budget: the author's deposition diagnostic resets each 5 s, while saved profiles are 600 s apart. The pinned writer also assigns a scalar to the full `Fmono` vector, so that field is excluded from monomer-conservation claims.

A single seed and a single time step do not establish stochastic or discretization convergence. The common seed also does not enforce identical birth samples between cases, because collision sampling consumes additional random numbers. The present comparison preserves the recovered script and archive-label settings. Future source refinements must also preserve the integer births per update: the source floors its expected superparticle count. Changing multiplicity to 3,300 at 5 s, or halving the time step alone, would truncate 37.5 to 37 and alter the physical source. Those CLI combinations have not been admitted as controls. A coordinated change to multiplicity 825 and step 2.5 s retains 75 births per update and the same physical source rate.

## What the history variable means

For `habit0/agggeo1`, mass fixes diagnostic diameter and projected area, and aspect ratio remains 1. The selected sticking law depends only on temperature. Capacitance distinguishes monomers (`mm<1.01`) from aggregates: their unfilled coefficients are `1/π` and 0.25, with the source's void-volume branch retained. A boolean single-crystal/aggregate state therefore suffices for these selected closures. The author solver still transports and adds actual monomer counts, and other geometry or sticking options would need more history. The [independent component adapter](welss_aggregation_physics.py) is separate from this author execution.

## Reproduction

Run from the atmosphere worktree with GNU Fortran, Python, NumPy and SciPy available. The build fetches only missing manifest files and rejects mismatches. Use new scratch directories; both build and run commands reject existing destinations.

```sh
python -m atmosphere.snow.welss_aggregation_build --source /tmp/welss-reproduction/source --output /tmp/welss-reproduction/build --fetch
python -m atmosphere.snow.welss_aggregation run --build /tmp/welss-reproduction/build --output /tmp/welss-reproduction/aggregation --aggregation --seed 20261021
python -m atmosphere.snow.welss_aggregation run --build /tmp/welss-reproduction/build --output /tmp/welss-reproduction/growth --seed 20261021
python -m atmosphere.snow.welss_aggregation compare --run /tmp/welss-reproduction/aggregation --archive /tmp/lunar_snow_earth_population/ob_agg0.nc --output /tmp/welss-reproduction/aggregation.json
python -m atmosphere.snow.welss_aggregation compare --run /tmp/welss-reproduction/growth --archive /tmp/lunar_snow_earth_population/ob_hp0.nc --output /tmp/welss-reproduction/growth.json
```

The comparison requires separately recovered archive members with the pinned hashes; the build does not substitute data or download the complete 69 MB archive.

## Source and archive identity

The target is the diagnostic-geometry solid-line case selected by `f21.ncl` in the [author archive](https://doi.org/10.5281/zenodo.10077724). The growth control targets `ob_hp0.nc` selected by `f19.ncl`. Both 40,036-byte members were checked against extracted-member CRC32 and recorded SHA256 values. The full data ZIP checksum remains unverified; the separately retrieved complete plotting-script ZIP matched its published MD5. Archive credit is Welss, Siewert and Seifert (2023), CC-BY 4.0. Products retain derived model profiles and comparison metrics; external NetCDF arrays are not copied into Git.

Inspection covered paper sections 5.1–5.3 and Figure 21, the figure scripts and member attributes, and pinned source initialization, collision, growth, geometry, velocity, output and namelist paths. Exact source-file URLs and hashes are in the manifest. The pinned source dates to September 2024; archive processing dates to February 2023. The archived executable, compiler, seed and exact averaging window remain unidentified. Agreement or disagreement therefore measures this declared-source reconstruction, with those provenance limits attached. No laboratory or field validation is claimed.
