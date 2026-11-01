# Reproduce the Earth benchmark, then continue snow-model development

## Handoff and scope

Prepared 23 October 2026 at the author's request to end the current session with a plan commit. The scientific checkpoint is `f18ee9b`. The interrupted continuation inspected records and source only; it ran no new scientific experiment and changed no model. This plan records the next work, not completed reproduction or an adopted repair.

Work exclusively in `/home/kir/Documents/Projects/terluna-atmosphere` on `domain/atmosphere`. The main checkout at `/home/kir/Documents/Projects/terluna` is not the work location. Start by checking the branch and working tree, then read the root working agreements, `research/README.md`, `research/status.json`, `atmosphere/README.md` and this plan.

The author has authorized reproducing the published Earth case, extending its duration and conditions individually, and using the findings to continue development toward the full lunar snow column. Continue beyond the first benchmark result to its next discriminating experiment. Keep the author's suspicion that the models cannot handle a full Earth storm separate from an established finding: the current single-crystal discrepancies do not demonstrate that broader claim.

## What has actually been used

- **Existing software:** compiled SCALE-SDM routines in a custom column adapter, and native McSnow deposition, habit, geometry and drag routines called through a custom Fortran driver. These are real source-code reuse, but not exact reproductions of an identified published executable.
- **New implementations:** the recent spatial diffusion, moving-interface and surface-transfer calculations implement published equations and hypotheses in our own code. They are not the original Libbrecht software. Their numerical checks do not establish autonomous Earth morphology or lunar validity.
- **Current limitation:** neither an accepted absolute lunar flake distribution nor a validated storm optical forecast exists. Preserve the earlier results and their evidence labels.

Read [McSnow growth](mcsnow_growth.md), [Earth growth assessment](earth_growth_assessment.md) and [model alternatives](growth_options.md) for the completed work. Reuse the existing runners and cached outputs before adding new infrastructure.

## 1. Establish the authors' actual Earth experiment

Start with [Welss et al. (2024)](https://doi.org/10.1029/2023MS003805), especially Section 4.2, Figure 10, the revised-growth comparison and the associated public code/data. The paper explicitly evaluates mass, both crystal axes and apparent density after ten minutes. These quantities are within its intended scope. Locate the author's configuration, driver, plotting scripts and raw output where available; distinguish them from a reconstruction based only on prose or plotted curves.

Build a short source-to-experiment record covering:

- Source revision/archive, source hashes, growth tables, compiler and options; identify which figure and model variant each build represents.
- Temperature grid, pressure, saturation convention, seed radius versus full diameter, seed mass/density, elapsed time, timestep, fall-speed option, ventilation and growth/branching settings.
- Axis definitions and units, geometry/volume convention, output averaging and plotted transformations. Compare against author calculations separately from comparison against laboratory data.

The paper specifies a 2 micrometre initial radius for its spherical-seed comparison. Verify all remaining settings from the author artifacts instead of inheriting our earlier choices. Our earlier driver imposed Murphy–Koop water saturation, which McSnow diagnoses as 100.119–100.721% relative humidity over water; an exact native-saturation control must distinguish that difference.

A decisive unresolved clue is in the pinned 2024 `src/mo_depo.f90`, around lines 393–394: the comment directs disabling the two multiplications of vapour diffusivity and heat conductivity by mean ventilation to reproduce an older paper version. Determine which experiments that comment concerns. Preserve the unmodified source and build any historical candidate separately with a recorded diff. Removing these two operations is not equivalent to removing only the extra thin-plate term. The existing control framework includes a unit-mean-ventilation option; inspect its scope before reusing it as a historical reconstruction.

The 2026 source archive lacks the growth tables, so the previous 2026 comparisons explicitly used the pinned 2024 tables. Do not silently call that combination the exact published package. Other source/figure provenance questions remain in the earlier archive records.

**Completion criterion:** reproduce the author output for the identified configuration with tolerances declared from output precision and numerical refinement. If only a figure is available, quantify digitization uncertainty and label the result accordingly. If an exact executable/configuration cannot be recovered, state the missing evidence and compare the documented reconstruction candidates; do not promote a candidate to an exact reproduction. An unexplained mismatch remains a reproduction issue before it becomes grounds for replacing the physics.

## 2. Extend a matched case and locate the first divergence

Keep the reproduced ten-minute case as a regression control. Extend only the duration to twenty and thirty minutes, then examine the other laboratory temperatures. Compare mass, both full axes, apparent density and fall speed throughout the growth history. Retain failures and partial trajectories. Use timestep refinement where a discrepancy or instability appears.

The existing comparison at −14.4 °C is a useful diagnostic, not the reproduction target: the revised spherical-seed package gives 1.587 mm at ten minutes against a fitted 1.361 mm, then 195.389 mm at thirty minutes against 4.175 mm. Removing only the extra mean plate-ventilation term yields 5.878 mm; also removing the directional ratio gives 4.676 mm. Those interventions do not jointly repair all temperatures and starting histories.

Trace the first departure through actual source terms: axis ratio, deposition density, branching state, capacitance, terminal speed, mean ventilation and directional growth. Separate the historical ventilation switch, native saturation and timestep effects before changing additional physics. Treat fitted three-minute restarts as a separate sensitivity, not an exact reconstruction of a tracked laboratory crystal. The laboratory curves summarize different crystals, and some mass/axis combinations depend on a hexagonal-prism versus spheroid mapping.

**Decision:** if the historical configuration reproduces the paper but fails later, delimit that model's demonstrated growth interval and identify the responsible extrapolation. If our driver or settings caused the mismatch, repair the interface and regenerate affected comparisons before changing physical laws. If several explanations remain, run the smallest control that distinguishes them.

## 3. Use that result to continue toward the storm and lunar column

After the reproduction assessment, carry out the next justified control rather than stopping at another model review. Prefer existing native code and a bounded, source-supported change over a new crystal solver. Preserve both native and altered versions so an apparent improvement can be checked against the reproduced benchmark and other habits.

First distinguish prescribed single-crystal growth from population growth with finite water, latent heating, sedimentation and collisions. Use the existing Earth column/aggregation controls and author column outputs to test the appropriate population processes; a good individual-crystal width alone does not validate a storm. Compare precipitation and suspended mass, number and size distributions, and the water/heat accounts appropriate to each experiment. Keep an imposed-temperature cooling term explicit.

Only then transfer the accepted components through controlled changes in local pressure, gravity, cloud depth, source flux and humidity history, using the existing matched SCALE lunar comparison where applicable. These are separate physical changes; the great altitude of the lunar snow zone is not itself a new microscopic growth law. If a component remains unvalidated, carry that uncertainty into the column result rather than reporting a calibrated lunar size distribution. Aggregate fracture, habit frequencies and visible storm optics remain separate unfinished requirements.

## Files and execution records

| Role | Existing material |
|---|---|
| Native McSnow builder and caller | `mcsnow_crystal_growth.py`, `mcsnow_crystal_driver.f90` |
| Growth matrix and source audit | `mcsnow_growth_benchmarks.py`, `results/mcsnow_growth_source_audit.json` |
| Ventilation interventions and parity controls | `mcsnow_ventilation_controls.py`, `mcsnow_ventilation_driver.f90`, `ventilation_benchmarks.py` |
| Laboratory fits and recovered measurements | `takahashi_growth.py`, `takahashi_measurements.py`, `results/takahashi_measurements.json` |
| Author column/aggregation comparison and provenance | `welss_column.md`, `welss_aggregation.md`, `welss_archive_sources.json`, `welss_aggregation_sources.json` |
| Existing full lunar comparison | `scale_lunar_methods.md`, `scale_lunar_results.md`, `scale_size_results.md` |
| Software source records | `solver_reuse_sources.json`, `growth_model_sources.json` |
| Ignored source/build/output caches | `research/runs/snow_scale/mcsnow_growth/`, `research/runs/snow_scale/mcsnow_ventilation/`, other existing `research/runs/` caches |

Keep new scientific code and records in `atmosphere/snow/`; external source bytes and raw runs remain ignored. Preserve pinned baselines. Update the current assessment and README when findings change their conclusions, identifying actual software reuse separately from new implementations.

Account for other machine activity. Use at most two numerical/build workers in total, with `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1`; use fewer when other sessions are active. Add only meaningful checks for any new interface or intervention. Run the repository-required `make check` before each commit and record actual results and evidence limitations. Commits at a meaningful reproduction finding or method change are authorized; pushes are not requested.

## Return to Earth growth and separate collection controls, 23 October 2026

The author directed the next work back to the remaining Earth growth error before a full lunar aggregation study. The [completed shape and collection diagnosis](mcsnow_shape_diagnosis.md) adds 48 trajectories, 148 fixed-state probes, 36 directional controls and 1,098 native pair evaluations. Equal-mass comparisons and fixed-state source rates isolate a shape-evolution failure. Unit directional ventilation improves the spherical dendrite but fails fitted restarts; prescribed fitted axis ratios expose remaining mass and initial-history errors and are not adopted as a growth law.

The fixed-property collection interface now passes independent implementation, symmetry, equal-particle and sticking controls. Its comparisons with fitted crystal properties are conditional kernel sensitivities, not observed aggregation or a completed population validation. Continue with a source-supported joint treatment of effective axis growth and deposition density, tested across habits, ages and starting histories. Keep the earlier numerical-archive reconstruction issue separate. The full finite-water Earth–Moon population comparison remains conditional on adequately constrained growth and aggregate geometry; the latest lunar humidity paths do not supply that validation.

## Automatic tuning checkpoint, 23 October 2026

The author requested an automatic search to locate the adjustments needed by the growth comparisons. The [completed tuner](mcsnow_tuning.md) runs fourteen bounded searches, including two shared starts, separate temperature/history fits and six neutral-control refits. The developed dendrite can predict the reserved thirty-minute mass and dimensions within 4.5%, while fitted prescriptions transferred to tiny seeds fail strongly. The common fit improves its training score by 25.5% but still leaves major axis errors. Controls compensate for each other, so the fitted numbers remain diagnostic.

Continue with a discriminating test of early crystal shape and growth allocation as a function of size/shape, constrained jointly by fixed-state rates and full histories. Test branching/hollowing and observation-to-envelope mapping alongside it. Preserve the training/withheld distinction and the native baseline; use the tuner to compare explicit candidate laws rather than promote a case-specific coefficient set. Independent early-state and humidity evidence, native drag/projected-area limits, and the original author configuration remain unresolved.
