# Automatic calibration of the remaining snow-growth errors

The bounded tuner finds close fits for individual crystal histories, but its best shared settings still fail the joint comparison across starting states. Fourteen searches and 13,218 cached parameter evaluations reduce the common training score from 0.18049 to 0.13446, a 25.5% reduction in squared logarithmic error. The corresponding RMS logarithmic error falls from 0.4248 to 0.3667. These are fitting scores, not percentage error bars or evidence that the tuned values are physical constants. No tuned replacement enters the lunar calculations.

The clearest success is the −14.4 °C three-minute restart: at the reserved thirty-minute age, its separately fitted settings predict mass, width and thickness within 4.5% of the laboratory fits. The tiny-seed fits remain less reliable at that later age. The tuned warm spherical-start case is 45.8% too massive, while the tuned spherical dendrite is 24.9% too thin. This contrast makes initial shape and its subsequent evolution a priority for the next model test.

## Which controls changed

The two independent common searches reach nearly identical minima. The selected common candidate is also compared with all six fixed-control refits on the same training objective; none of the withheld results is used to select it. Its controls are:

| Control | Original value | Best shared value |
|---|---:|---:|
| Transport multiplier | 1 | 0.766 |
| Axis-growth multiplier | 1 | 1.06 |
| Directional-ventilation exponent | 1 | ≈0 |
| Deposition-density multiplier | 1 | 0.863 |
| Branch-radius multiplier | 1 | 0.294 |
| Mean-ventilation enhancement | 1 | 1.17 |

The branch-radius factor corresponds to about 59 micrometres instead of 200 micrometres. The other multipliers are effective closure adjustments. In particular, a lower transport multiplier cannot be read as a measured correction to vapor diffusivity: it can compensate for an incorrect evolving crystal geometry. The earlier fixed-state diagnosis found only 3.7% excessive uptake when the late dendrite geometry was supplied directly.

Fixing each control at its original value and reoptimizing the other five tests whether its adjustment is indispensable within this search. The table compares increases in the squared-log fitting score against the best common score. The direct-reset column holds the other controls fixed.

| Control fixed at original value | Direct reset: score increase | Other controls refitted: score increase |
|---|---:|---:|
| Transport multiplier | 77.8% | 5.6% |
| Axis-growth multiplier | 5.1% | 4.0% |
| Directional-ventilation exponent | 5.6% | 3.7% |
| Deposition-density multiplier | 5.4% | 0.9% |
| Branch-radius multiplier | 55.8% | 11.2% |
| Mean-ventilation enhancement | 1.0% | 0.1% |

Earlier branching gives the largest advantage in these fixed-control comparisons, followed by the transport adjustment. The mean-airflow and deposition-density adjustments are particularly replaceable. Keeping the original directional correction raises the best score by only about 3.7%, despite the unconstrained optimum almost removing it. These results identify compensating controls rather than a uniquely determined physical repair. Several alternative fits reach a search boundary or the solid-density ceiling. The product also records local finite-difference sensitivities; these are not statistical parameter uncertainties.

## Different starting histories require different compromises

The following searches fit each history separately. The axis multiplier acts on the logarithmic polar-to-equatorial growth ratio. The final three columns are signed errors at thirty minutes, which was excluded from fitting. Polar dimension means thickness for the dendrite and length for the column.

| Temperature and start | Axis multiplier | Directional exponent | Training RMS log error | Mass at 30 min | Width at 30 min | Polar dimension at 30 min |
|---|---:|---:|---:|---:|---:|---:|
| -14.4 °C, tiny sphere | 0.805 | ≈0 | 0.0933 | +2.2% | +7.0% | -24.9% |
| -14.4 °C, three-minute restart | 1.836 | 0.0405 | 0.0302 | -0.7% | -4.5% | -3.9% |
| -5.3 °C, tiny sphere | 0.974 | 2 | 0.0747 | +45.8% | +9.8% | +33.6% |
| -5.3 °C, three-minute restart | 2.531 | ≈0 | 0.0580 | -3.5% | +2.0% | +14.5% |

The dendrite histories favor axis multipliers of 0.805 and 1.836; the column histories favor 0.974 and 2.531. Other controls also change, so these ratios alone do not establish a unique facet-growth correction. The branch-radius control is inactive in these column cases, which use the native immediate-hollowing rule; its separately fitted column value carries no information. Giving both starts at each temperature a common parameter set leaves RMS log errors of 0.1828 for the dendrite and 0.3544 for the column, above the corresponding individual-history fits.

Applying a fitted restart prescription to the tiny sphere is a stronger failure. The dendrite-restart prescription gives the spherical case only 23.1% of fitted mass and 26.4% of fitted width, with 2.72 times the fitted thickness at thirty minutes. The column-restart prescription gives 46.0 times the fitted mass and 178 times the fitted length from the spherical seed. These are numerically resolved consequences of the diagnostic controls, not credible predictions of real needles: halving the timestep changes the cross-history comparisons by at most 0.0291%. The same column prescription remains much closer for the volume-preserving three-minute restart, with mass −7.4%, width −2.4% and length −3.2%.

The next discriminating change should therefore test how the early crystal state and size/shape-dependent growth allocation interact with branching and hollowing. Fixed-state source rates should constrain that change alongside complete histories, so an uptake adjustment cannot merely hide an erroneous geometry. Unknown early laboratory shapes and the use of separate fitted populations remain competing explanations. The present searches do not identify a unique molecular attachment mechanism or establish that every possible constant-parameter family must fail.

## Predictions outside the fitting selection

Median absolute percentage errors across the declared comparison states are shown as original → shared-tuned. These are descriptive summaries with different state counts, not uncertainty intervals. Fall speed was excluded from fitting.

| Comparison group | States | Mass | Width | Polar dimension | Fall speed |
|---|---:|---:|---:|---:|---:|
| Fitting ages and starts | 12 | 29.7% → 17.1% | 19.8% → 12.3% | 21.2% → 29.1% | 40.9% → 19.5% |
| Reserved ages, same starts | 10 | 34.6% → 33.0% | 30.4% → 18.3% | 40.9% → 35.4% | 50.0% → 27.4% |
| Alternative volume mapping at the two fitted temperatures | 10 | 34.1% → 34.1% | 37.8% → 30.3% | 49.0% → 40.2% | 37.7% → 26.8% |
| Other eight temperatures, all admissible starts | 68 | 24.4% → 13.0% | 8.2% → 16.2% | 9.4% → 7.7% | 44.8% → 23.1% |

The shared adjustment improves several summaries, but worsens the training polar-dimension median and the width median at the other temperatures. At thirty minutes it still leaves the restarted dendrite 67.1% too wide and the restarted column 72.2% too short. The absence of a native execution failure across these reserved cases establishes runnable predictions, not joint physical agreement.

## Experiment and interpretation

The author requested an automatic search to identify which model controls can account for the remaining growth discrepancies. The [runner](mcsnow_tuning.py) patches a separate build of the already retained improved McSnow component: pinned 2024 revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`, with the extra oblate mean-ventilation term removed. The [batch adapter](mcsnow_tune_driver.f90) calls native growth, geometry and drag. The original author executable/configuration remains unidentified; these are documented custom-caller comparisons. Original source, growth tables and preceding numerical products remain unchanged. Every evaluated vector and its sampled outputs, including unsuccessful evaluations, are retained on the project-space drive.

| Control | Search interval | Operation |
|---|---:|---|
| Transport multiplier | 0.25–4 | Scale the capacitance used in both vapor and heat transfer |
| Axis-growth multiplier | 0.25–4 | Scale the effective logarithmic polar-to-equatorial growth ratio, leaving the intrinsic table used by deposition density unchanged |
| Directional-ventilation exponent | 0–2 | Raise the native directional correction to this power; 0 removes it, 1 preserves it |
| Deposition-density multiplier | 0.25–4 | Scale porous deposition density after native branching/hollowing, with an explicit upper limit at solid ice density; compact pre-branch growth remains solid |
| Branch-radius multiplier | 0.25–4 | Scale the native 200 micrometre branching radius, giving 50–800 micrometres |
| Mean-ventilation enhancement | 0–2 | Scale the enhancement above unity for both vapor and heat transport |

All six neutral values are one. The intervals are exploratory search limits, not measured uncertainty bounds. The transport multiplier can represent an error in geometric transport or other uptake assumptions; it does not uniquely identify a diffusion coefficient. Directional growth, effective deposition density and the intrinsic temperature table are deliberately separated for diagnosis. These controls are constant along a trajectory. Temperature, pressure, humidity and initial mass and axes are held at their declared values.

The common search fits four trajectories: −14.4 and −5.3 °C, each from the original 4 micrometre full-diameter sphere and from a fitted three-minute crystal with its mass and axes preserved. It minimizes the mean squared logarithm of model-to-fit ratios for mass, width and polar length at ages 5, 10 and 20 minutes. Each trajectory and each quantity receive equal weight. There are 36 residuals; fall speed is checked separately. Ages outside that fitting selection, alternative volume-preserving restarts, and the other eight temperatures do not enter the common fitting objective.

The inputs are the identified Takahashi et al. (1991) fitted histories already used in the [shape diagnosis](mcsnow_shape_diagnosis.md), with the same inadmissible-density exclusions and explicit prism-to-spheroid mappings. Their separately fitted segments summarize different crystals and have no pointwise error model; planar thickness is partly inferred from mass and basal area. Consequently, the loss is a diagnostic score, not a chi-squared statistic. Withheld ages share the same laboratory fits, and alternative mappings measure representation sensitivity. Agreement on either is weaker evidence than comparison with an independent experiment.

The search uses SciPy 1.16.2 [differential evolution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html) followed by bounded [least squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html). Multipliers use logarithmic search coordinates; the two airflow controls use linear coordinates. Two random seeds initiate the shared search. Separate fits then allow each temperature, and each temperature/initialization combination, its own six values. Six further searches fix one control at its neutral value and reoptimize the other five. A direct reset without reoptimization is reported separately, because other controls can compensate for the removed adjustment.

Each differential-evolution population has six members per free dimension. The generation budgets are 35 for each initial shared search, 25 for temperature and fixed-control searches, and 20 for individual histories. Local refinement allows 160 function evaluations, with a relative finite-difference step of 0.005. Its termination criteria and each global-search budget exit are retained. Finite budgets and multiple starts do not establish a global optimum. The result is the best identified compromise within this declared family, with local sensitivity and parameter-boundary diagnostics rather than confidence intervals.

Native search integrations use a 1/16-second step; selected solutions are repeated at 1/128 second. Every successful evaluation checks mass accounting, spheroid volume and the solid-density bound. A failed native trajectory receives a fixed penalty and remains in the journal; its last completed states are retained. Failure during a held-out case is isolated so that the remaining cases are still evaluated. A native drag-domain failure is an inadmissible numerical candidate, not an experimental rejection of a particular material property.

Reproduction requires the previously pinned McSnow 2024 baseline build and uses one numerical worker, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`:

```sh
python -B -m atmosphere.snow.mcsnow_tuning search
python -B -m atmosphere.snow.mcsnow_tuning analyze
```

Inputs, build identity, optimizer settings, sampled predictions, checks and cache hashes are recorded in [the numerical product](results/mcsnow_tuning.json). The resolved cache is under `/media/projectspace/terluna-research/atmosphere-snow/snow_scale/mcsnow_reproduction/growth_diagnosis/tuning/`.

## Verification

All twenty numerical checks pass. The unchanged adapter agrees with all 114 saved reference states to a maximum relative difference of 1.55e-15; the unit-directional control agrees exactly at its saved states. The largest selected-solution change between 1/16 and 1/128 second is 0.1836%. Across the saved outputs, the largest relative mass-ledger residual is 2.22e-16 and the largest spheroid-volume residual is 2.33e-15.

The journals retain 680 failed parameter-batch evaluations, principally native drag arithmetic failures for extreme proposed shapes. They receive penalties and do not qualify as fitted solutions. Twelve focused tests pass, covering held-out-data exclusion, equal objective weights, failure penalties, coordinate bounds, geometry conversion and the saved numerical controls. Exact raw reassembly matches the saved product. `make check` passes: 3,733 Python tests passed, 20 skipped, and all 191 JavaScript tests passed; provenance and ensemble integrity pass. A further layer check after staging the new files covers 1,323 files with zero violations. These activities are recorded in [the verification record](results/mcsnow_tuning_checks.json), separately from empirical validation.
