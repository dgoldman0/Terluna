# Diagnosing volume dependence in the lunar snow column

The completed [full-column calculation](event_column_results.md) conserved water while its arriving size distribution changed with the represented volume. This study tests the suspected coupling between sparse particle populations, local vapour allocation and aggregation. Its reader is the researcher deciding which numerical or physical change to investigate next. The experiments diagnose the prescribed equations; they do not establish lunar snowflake sizes or a new physical breakup law.

## Observations preserve the original calculation

The [diagnostic C++ driver](column_diagnosis.cpp) includes the frozen `event_column.cpp` and calls its particle-property, growth, collision, injection and settling routines. Only the orchestration and output conversion are repeated in the [diagnostic interface](column_diagnosis.py), allowing the original producer and its saved results to retain their identities. At a fixed seed, every original scientific field must match the original result exactly. The observer consumes no random numbers and changes no particles. The deliberately empty unit-test limit reports a zero relative water residual when both supplied and retained water are zero; the original interface leaves that 0/0 ratio undefined. Every production comparison has positive supplied water. Small-column tests cover the complete model and each process removal; the full baseline and time-control runs are also compared with their previously saved realizations.

Each growth substep has two Heun stages. The observer measures both stages using the same quadrature weights as the growth update. In a layer of volume V, the supply rate is c = qV and the sum of unrestricted particle demands is D. The admitted rate is min(D,c). An empty layer loses c of available capacity; an occupied layer with D < c loses c−D; a layer with D ≥ c uses its capacity. These categories therefore satisfy

`admitted + empty unused + occupied unused = available supply`.

The integrated admitted mass is checked against the actual deposition ledger. Occupied unused supply is a demand shortfall, not water deliberately discarded after uptake. The prescribed model has no stored vapour reservoir through which unused capacity can be recovered later. This classification accounts for the calculation's unused supply; it does not predict what the evolving population would do if that supply were redistributed.

The observer also records empty, single-particle and supply-limited time, the first two particle-number moments, the effective number of particles sharing growth demand, and the largest particle's share of that demand. The effective count is `(sum g)^2 / sum(g^2)` for individual unrestricted demands g; it measures concentration of growth demand and is bounded by the actual count. All time statistics are observations at the split scheme's growth stages, rather than exact continuous residence-time integrals. The half-step control tests their sensitivity to that sampling. Sticking, chipping and emitted-fragment counts are measured directly around each layer's collision substep.

The product stores a filling block and four gathering blocks, each resolved by layer. Layer zero borders the freezing level. Statistics from gathering are reported separately from startup, and whole-realization standard errors describe variation among four runs. The time blocks are not independent samples.

## Averaging demand before applying the cap

The [offline analysis](column_diagnosis_analysis.py) retains the observed particle histories and changes only the order of arithmetic. For each layer and each of the four gathering blocks it compares the integrated admitted demand with the smaller of integrated unrestricted demand and integrated supply. Because `min(D,c)` is a concave function of D, the latter is at least as large as actual admission. Their difference measures how much admission depends on demand varying across the two sides of the cap within that block.

For example, suppose a layer is empty for half an interval and demands 1 kg/s for the other half, with capacity 0.5 kg/s. Its average admission is 0.25 kg/s, whereas applying the same cap to its average demand gives 0.5 kg/s. This is an accounting example, not an additional snow simulation. The corresponding diagnostic on the saved columns uses their actual Heun-stage weights and water integrals.

This calculation smooths all within-block variability, including changes over time; it does not identify every fluctuation as a sampling error. It also leaves the particle masses, subsequent collisions and settling fixed. The result therefore quantifies the nonlinearity on those histories and supplies no feedback-corrected snowfall rate or size distribution. The analysis separately decomposes changes in actual deposition into changes in empty and occupied unused supply, and changes in mean arriving particle mass into mass-flux and number-flux ratios. Quantiles of the ensemble-mean mass CDF are retained as native-grid brackets, separately from means of realization quantiles.

## Controlled process removals

All primary runs retain the −15 °C cloud top, 30 layers, lunar gravity, inherited geometry, expected seed source, upwind settling generator and 132.589-hour duration. Baseline layer volumes are 0.05427, 0.21707 and 0.86829 m³. The three removals compare the first two volumes with four realizations at each setting:

| Experiment | Change | What persistence of the volume effect would establish |
|---|---|---|
| Complete column | All inherited processes retained | Reproduction of the original failure with passive diagnostics |
| No deposition | Disable vapour growth; retain encounters, sticking and chipping | Vapour allocation is not required for volume sensitivity under those altered particle histories |
| No aggregation | Disable all encounters, including collision-generated chipping; retain vapour growth | Collisions are not required under the resulting growth and settling histories |
| No chipping | Retain sticking; nonsticking encounters remove no mass and create no fragments | Chipping and fractional-fragment delay are not required |
| Half macrostep | Complete model at 30 s instead of 60 s, middle volume | Sensitivity of the diagnostic observations and distribution to splitting |

The primary matrix contains 40 full runs. An additional eight-run [sticking-only comparison](column_diagnosis_aggregation_only.py) disables both vapour growth and chipping, retaining sticking, the source and settling at the first two volumes with four realizations each. It addresses the residual upper-tail sensitivity in the growth-disabled comparison, where chipping remains active. The follow-up waits for the primary pool to finish, binds the completed primary product by hash, and reuses its unchanged solver and case runner. Together the studies contain 48 full runs.

Process removal changes the particle population and its subsequent rates. Differences among these experiments are consequently not additive contributions from independent processes. If the coupled model changes with volume while a reduced model does not, the finding identifies a dependence at the tested conditions; it does not prove a universal instability or identify an empirically correct repair.

The existing numerical screen is retained: maximum change of mean mass CDF at most 0.05, and changes in mean flux and mass-weighted 50th, 90th and 99th percentiles at most 5%. These thresholds are practical comparisons, not confidence bounds. CDF differences are evaluated at the exported diameter queries; the query grid is an output diagnostic and does not discretize particle masses in the solver. Four realizations may leave upper quantiles imprecise. A reduction of the tail is not itself an acceptance criterion, and no maximum flake diameter is introduced.

## Evidence and reproduction

The [runner](column_diagnosis_runs.py) records specifications, source and input hashes, compiler identity, actual completion or failure, individual outcomes and group uncertainties. It retains failed or missing cases explicitly and uses at most two numerical workers. Focused tests check original-output identity, independent water and event-count accounting, an empty-column limit and Poisson-source occupancy. These establish implementation properties separately from the full-column results.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.column_diagnosis_runs --run-all --workers 2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.column_diagnosis_aggregation_only --workers 2
python -m atmosphere.snow.column_diagnosis_analysis
python -m pytest -q atmosphere/tests/test_snow_column_diagnosis.py atmosphere/tests/test_snow_column_diagnosis_runs.py atmosphere/tests/test_snow_column_diagnosis_analysis.py
```

This extension uses the physical assumptions and source records of the frozen event column. The diagnostics above are local accounting identities and controlled numerical experiments. They add no external dataset, literature-derived closure or empirical validation.
