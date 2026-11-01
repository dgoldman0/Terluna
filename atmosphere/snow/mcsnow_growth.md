# McSnow crystal growth: starting shape and late growth failure

The subsequent [Earth growth assessment](earth_growth_assessment.md) completes the ventilation controls proposed here. It identifies the extra plate term as a major cause, tests separate humidity and seed-shape sensitivities, and compares transport with measured shape histories. This report preserves the preceding native baseline.

Native McSnow improves some Earth crystal-growth comparisons, but neither tested growth package provides a reliable replacement across the five temperatures. The most informative result is the revised package's dendritic crystal: from a 2 µm-radius spherical seed, its width is 1.59 mm at ten minutes against a fitted 1.36 mm, yet reaches 195.39 mm at thirty minutes against 4.17 mm. Its mass then exceeds the fitted value by a factor of 1777. These are failed empirical predictions of an isolated crystal, not evidence for real decimetre snowflakes.

The comparison follows checkpoint `6b5f3e6`, which contains the [SCALE growth controls](growth_benchmarks.md), [model review](growth_options.md) and recovered laboratory measurements. The [main product](results/mcsnow_growth_benchmarks.json) records 80 native executions; [five additional refinements](results/mcsnow_growth_refinements.json) investigate the sensitive cases. No lunar particle law has been replaced.

## Conditions and source identity

The [Python builder and runner](mcsnow_crystal_growth.py) links a [thin Fortran caller](mcsnow_crystal_driver.f90) to unchanged native deposition, habit, geometry and drag routines. It builds the pinned public 2024 revision and the versioned 2026 source archive separately. Both use habit option 1, Böhm velocity option 3 and enabled prolate interpolation. Original IGF1 and revised IGF2+ are tested as complete native packages; the revised option changes the growth table, small-crystal restriction and branching criterion together. The 2026 archive omits the growth tables, so both builds explicitly use the pinned 2024 tables. [Welss et al. (2024)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS003805), [2026 source archive](https://doi.org/10.5281/zenodo.17432967).

Each case holds temperature, pressure at 101000 Pa, and the same Murphy–Koop vapor partial pressure used by SCALE fixed. McSnow diagnoses 100.119–100.721% relative humidity over water under its own saturation convention; this difference is retained. Native Earth gravity and transport constants remain unchanged. One crystal grows without collisions, liquid collection, motion, vapor depletion or atmospheric heating. Terminal speed still supplies ventilation. This prescribed reservoir tests growth under maintained laboratory forcing; it supplies no finite atmospheric water budget.

Two initial histories answer different questions. Fitted restarts impose the empirical mass and full equatorial/polar dimensions at age three minutes, with spheroid volume πa²c/6. Spherical starts use native-density ice of radius 2 µm at age zero. The former shares SCALE's exact initial geometry; the latter tests the earlier growth history. The five admissible fitted states are used at −3.7, −5.3, −10.6, −14.4 and −20.1 °C. The other five fitted states still exceed bulk-ice density under the spheroid mapping and remain rejected. [Takahashi et al. (1991)](https://www.jstage.jst.go.jp/article/jmsj1965/69/1/69_1_15/_article/-char/en).

These are native process comparisons, not reproductions of the archived paper executable. Welss's ten-minute comparisons and author plotting scripts provide useful configuration evidence, but the source contains an explicit comment identifying an older implementation that omitted two mean-ventilation operations. The affected figure/executable history remains unresolved. The current calculations preserve those operations. The recovered observations also informed the candidate growth formulations, so agreement here would not constitute independent empirical validation.

## Growth depends strongly on the initial history

At −14.4 °C, the refined 2026 results are:

| Initial history and package | Width at 10 min, mm | Width at 30 min, mm | Mass/fitted mass at 30 min |
|---|---:|---:|---:|
| Laboratory fit | 1.361 | 4.175 | 1 |
| Spherical seed, original IGF1 | 1.041 | 6.306 | 3.01 |
| Spherical seed, revised IGF2+ | 1.587 | 195.389 | 1776.95 |
| Fitted three-minute state, original IGF1 | 5.142 | Failed before endpoint | — |
| Fitted three-minute state, revised IGF2+ | 5.132 | Failed before endpoint | — |

Sphere entries use 0.0125 s steps; fitted entries use 0.05 s. Both fitted trajectories already exceed fitted mass eightfold at ten minutes, subsequently grow uncontrollably and trigger floating-point exceptions. Every dendritic fitted restart in the original matrix fails across both releases, packages and tested steps. The last saved minute is age eighteen; its enormous dimensions are timestep-sensitive failure diagnostics and are not physical sizes or converged endpoints.

The original spherical seed is already too narrow and thick at age three: 0.390 times the fitted width and 2.744 times the fitted thickness. Its milder later growth therefore begins from a different shape. In the revised seed case, the 200 µm maximum-radius branching threshold is crossed at age 5.3875 minutes. The fitted restart begins at an equatorial radius of 198.8 µm and crosses after only 0.9 s. A restart alone nearly bypasses that package's prebranch history. The revised seed's fairly close ten-minute result followed by extreme thirty-minute growth shows why the shorter comparison is insufficient for the extended growth interval.

Other temperatures also require mass and both axes to be assessed together. With the 2026 revised package and a spherical seed at −5.3 °C, the thirty-minute column length is 1.926 mm versus a fitted 1.861 mm, but mass is 2.03 times and width 1.165 times the fit. Starting from the fitted three-minute state instead gives mass 1.083 times the fit while column length is only 0.379 times and width 1.556 times the fit. The original spherical package gives a 17.56 mm column after refinement. Thus a plausible maximum dimension or mass alone can conceal a failed shape history.

Compact-plate cases improve relative to SCALE in some configurations. At −10.6 °C, for example, the 2024 original fitted restart gives dimension and mass ratios of 1.294 and 1.423, compared with SCALE's 1.580 and 1.992. This local improvement does not resolve the dendritic failure or select one package across temperatures.

## What changes between versions, and what may drive the failure

The [source audit](results/mcsnow_growth_source_audit.json) identifies two 2026 deposition-density changes: newly deposited ice cannot be denser than the current apparent crystal density, and one capped-column branch uses 300 kg/m³. Both can affect porous fitted starts. Primary geometry and drag along this selected path are otherwise unchanged. Spherical-start trajectories agree between versions to roundoff; fitted starts show differences. This supports examining the newer density rules jointly. No separate intervention has isolated their individual contributions.

The native equations contain a specific coupled feedback. As a plate becomes thinner relative to its width, directional ventilation increasingly favors lateral growth. Meanwhile, the mean plate-ventilation correction contains `0.0028 * Xeq**1.5 / phi`, where `phi` is thickness/width. At fixed air properties this term scales as fall-speed^(3/4) × equatorial-radius^(3/2) × polar-radius^(−3/4). Together with thin-spheroid capacitance, it permits superlinear mass uptake if thickness, apparent density and fall speed vary slowly. Both vapor and heat-transfer resistance are present; missing latent-heat resistance is not the explanation found in the driver audit.

This is a mechanism supported by source analysis, not yet a unique causal diagnosis. Directional growth, deposition density, drag and capacitance evolve together. Both versions also use volume^(1/3) in the directional helper where the paper describes an equivalent-sphere radius; the helper does not restore the geometric conversion. That source/paper difference is preserved, and its effect has not been isolated. The eventual Böhm arithmetic overflow explains the program's termination, while the excessive growth precedes that overflow.

The next discriminating controls are to remove only the extra mean plate-ventilation term, and separately to set the directional factor ratio to one, each with unchanged-mode parity and timestep checks. Such controls would identify contributions; they would not constitute a calibrated repair. Source-history reconciliation and lower-pressure/humidity comparisons remain necessary before transferring a selected treatment to the lunar column. Aggregate fracture and storm optics remain unresolved.

## Numerical checks and reproduction

The main matrix has 72 completed executions and eight native failures. The targeted refinements add three completions and two failures, giving 85 recorded executions in total. All failures remain explicit, with raw logs and valid partial trajectories. Native scientific source files have no patches, and strict Fortran bounds and floating-point traps are enabled.

For completed fitted restarts, 1→0.2 s refinement changes endpoint dimension by at most 0.177% and mass by 0.163%. Initial sphere refinements of 0.2→0.05 s reached 6.64% dimension and 7.23% mass changes in the original column case. Additional 0.05→0.0125 s checks reduce its changes to 1.16% and 1.27%; the original dendrite changes by 0.731% and 1.21%, and the revised dendrite by 0.389% and 0.783%. These observed changes are much smaller than the large empirical errors. They do not establish an asymptotic convergence order.

For the failing 2026 fitted dendrites, 0.2→0.05 s changes ten-minute dimension by less than 0.194% and mass by less than 0.398%. Late trajectories remain strongly timestep-sensitive, so no final size or exact failure time is accepted. The completed and partial runs together have 24,795 raw snapshot values checked against their saved analyses, regenerated inputs, source/executable/output hashes, and consistent mass/axes/volume mappings. Deposition accounting closes within 2.22 × 10⁻¹⁶ relative error; this is bookkeeping within prescribed air, not a finite-reservoir water or energy conservation test. [Verification record](results/mcsnow_growth_checks.json).

The focused suite passes 51 tests covering the laboratory inputs, native mappings, source tampering, malformed or missing outputs, and partial-run classification. Full `make check` passed for the preceding committed checkpoint; the new continuation has separate focused, layer and provenance checks. External source archives and raw trajectories remain in the ignored cache.

```sh
python -B -m atmosphere.snow.mcsnow_crystal_growth build --version 2024 --public-source /path/to/pinned/mcsnow --output research/runs/snow_scale/mcsnow_growth/build_2024
python -B -m atmosphere.snow.mcsnow_crystal_growth build --version 2026 --public-source /path/to/pinned/mcsnow --archive research/runs/snow_scale/mcsnow_2026_source.tar.gz --output research/runs/snow_scale/mcsnow_growth/build_2026
python -B -m atmosphere.snow.mcsnow_growth_benchmarks --run
python -B -m atmosphere.snow.mcsnow_growth_benchmarks --refine
# Recheck existing raw files and regenerate the main numerical product:
python -B -m atmosphere.snow.mcsnow_growth_benchmarks
```
