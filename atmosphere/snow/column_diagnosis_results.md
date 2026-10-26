# Why the prescribed lunar snow column changes with statistical volume

The controlled experiments identify a volume-dependent loss of vapour uptake when a sparse particle population fluctuates across the supply cap. They also retain an upper-tail refinement failure with sticking, source and settling alone. These are two distinct numerical findings; correcting water allocation alone has not been shown to make the size distribution converge. All numerical sizes in this report belong to the prescribed lunar column with inherited aggregate geometry. The study diagnoses its equations and numerical representation; empirical aggregate survival and predictive storm optics remain unresolved.

The study completes 40 primary runs and an eight-run sticking-only follow-up, each covering 132.589 hours with the last 66.295 hours used for arrival statistics. The prescribed cloud extends 16.281 km from −15 °C to freezing under the same lunar forcing as the [previous event-column study](event_column_results.md). Sixteen instrumented baseline and timestep realizations reproduce every original scientific output exactly; the observer changes no particle or random draw. The [primary product](results/column_diagnosis.json), [follow-up](results/column_diagnosis_aggregation_only.json) and [figure](../../visualization/snow/generated/column_diagnosis.png) retain the complete comparisons.

## Sparse layers change the water actually deposited

Each computational particle represents one particle in the chosen statistical volume. Collisions progressively reduce that population as snow moves toward the freezing level. In the complete model the lowest layer contains, on average, 0.414, 0.829 and 2.446 particles at successive fourfold volume increases. The smallest volume therefore samples a lower-cloud population through long empty periods and occasional particles with substantial growth demand.

The growth rule limits the sum of particle demands in each layer to that layer's prescribed supply. When the layer is empty, its current capacity goes unused. When it later contains a particle, that particle remains subject to the current cap; it cannot recover unused capacity from the earlier empty interval. The model keeps no vapour reservoir for that purpose. At the same concentration, changing statistical volume can consequently change uptake through the fluctuations of the represented population.

| Layer volume (m³) | Lowest-layer empty time (%) | Deposited vapour (mm/h) | Unused in empty layers (mm/h) | Unused in occupied layers (mm/h) |
|---:|---:|---:|---:|---:|
| 0.05427 | 64.27 ± 2.30 | 2.199 ± 0.016 | 0.571 ± 0.014 | 0.230 ± 0.003 |
| 0.21707 | 40.74 ± 2.44 | 2.373 ± 0.008 | 0.248 ± 0.008 | 0.380 ± 0.003 |
| 0.86829 | 7.53 ± 0.36 | 2.498 ± 0.004 | 0.026 ± 0.003 | 0.476 ± 0.004 |

Values are gathering-period means ± one standard error across four whole realizations. The three water terms account for the full column; the empty-time percentage describes its lowest layer.

The columns in the water account sum to 3 mm/h. Occupied unused supply rises with volume, partially offsetting the reduction in empty-layer losses. That occupied term measures insufficient instantaneous demand, so it must remain distinct from empty-layer loss. The accounted deposition increase is 0.17323 mm/h between the first two volumes and 0.12547 mm/h between the next two. Nothing in this attribution creates water beyond the prescribed allowance.

Applying the cap after averaging demand within each saved gathering block gives a second measure of the same nonlinearity. On the observed histories, the admitted rate under that arithmetic is 2.66657, 2.64210 and 2.63046 mm/h, whereas actual admission is 2.19927, 2.37250 and 2.49797 mm/h. The resulting gaps shrink from 0.46729 to 0.26959 to 0.13248 mm/h. Thus the complete model's increasing uptake accompanies a large reduction in the penalty associated with demand varying across the cap, even though the block-averaged-demand admission falls slightly.

The lowest layer indicates how sparse the population becomes, but its growth demand never reaches the cap at the recorded stages. Its averaging gap is therefore zero. The largest layerwise gaps occur higher up, around −3.75 to −4.25 °C, approximately 4.1–4.6 km above freezing, where demand repeatedly crosses the cap. At the largest volume the layer with the largest gap is almost always occupied, yet is supply-limited about 47% of sampled time. Fluctuating demand among occupied states therefore contributes as well as emptiness.

These averages retain the actual particle histories. They include temporal evolution as well as stochastic fluctuations within each 16.57-hour gathering block. Supplying the unused water to particles would change their later masses, collisions and settling. The calculation therefore quantifies the cap's effect on those histories and supplies no corrected size distribution. [The methods](column_diagnosis_methods.md) give the observation weights and accounting identity; [the offline product](results/column_diagnosis_analysis.json) preserves each realization and layer.

## Process removals separate the remaining questions

| Volume comparison | Maximum mass-CDF change (percentage points) | Flux change | Median change | 90th-percentile change | 99th-percentile change | Screen |
|---|---:|---:|---:|---:|---:|---|
| Complete, small → middle volume | 46.13 | +6.68% | +51.01% | +52.67% | +82.97% | Fail |
| Complete, middle → large volume | 28.99 | +5.08% | +26.59% | +69.97% | +104.24% | Fail |
| Vapour growth disabled | 3.55 | -0.26% | +1.64% | +10.57% | +33.35% | Fail |
| All collisions disabled | 0.13 | -0.01% | +0.05% | +0.09% | +0.19% | Pass |
| Chipping disabled | 43.52 | +6.50% | +51.45% | +43.95% | +77.55% | Fail |
| Both growth and chipping disabled | 2.68 | -0.21% | +3.15% | +1.03% | +22.50% | Fail |

Each volume comparison uses four independent complete realizations at both settings. The screen requires the maximum sampled mean mass-CDF change to be at most five percentage points, and mean flux and mass-weighted 50th, 90th and 99th percentiles to change by at most 5%. Quantiles in this table are means of realization quantiles. The thresholds are practical screens, not confidence bounds.

With collisions disabled, the lower layers remain well populated: the smallest-volume lowest layer holds about 1,738 particles on average. Every layer uses its full vapour capacity during gathering, and all five refinement metrics pass. The mean arriving 99th percentile is approximately 3 mm at both volumes. This establishes that deposition and settling can pass this volume comparison under the population produced without collisions; it does not establish their accuracy under the sparse population created by aggregation.

Disabling vapour growth leaves the seed ice, collisions, chipping and settling. Its mean mass median changes by 1.64%, but the 90th and 99th percentiles change by 10.57% and 33.35%. The mean 99th percentile rises from 14.18 ± 0.53 to 18.90 ± 0.48 mm, with one standard error across whole realizations. Quantiles of the ensemble-mean CDF give nonoverlapping native-grid brackets of 14.06–14.62 and 18.42–19.15 mm. Residual tail sensitivity therefore also appears when the CDF is averaged before taking its quantile. The complete model's vapour-allocation effect cannot be its sole explanation.

With chipping disabled and vapour growth retained, the mean median still changes by 51.45%, the 99th percentile by 77.55%, and the mean mass CDF by 43.52 percentage points. The median is 12.09 cm at the small volume and 18.31 cm at the middle volume. The main failure therefore does not require chipping or delayed fractional-fragment participation.

The eight-run follow-up disables both growth and chipping while retaining sticking, the seed source and settling. Its flux, median, 90th percentile and sampled mean CDF pass their screens, but the mean 99th percentile rises by 22.50%, from 15.95 ± 0.89 to 19.54 ± 1.48 mm. The ensemble-mean CDF gives 99th-percentile brackets of 16.41–17.05 and 19.90–20.68 mm. Thus the reduced sticking population also retains a failed upper-tail check without either removed process. Its four-realization uncertainties are appreciable; the experiment leaves the size of any systematic tail bias and its large-volume limit unresolved.

These removals change the evolving particle population and all its subsequent rates. Their effects cannot be added as independent fractions of the original failure. In particular, a pass after removing collisions diagnoses the dependence of this tested column on collisions; it establishes no universal instability or physical upper size.

## Finite sampling and the splitting control

The complete columns record only 102–114, 215–221 and 597–674 actual arriving particles per realization during the 66.29-hour gathering interval. The corresponding algebraic effective mass counts, calculated as `(sum m)^2/sum(m^2)`, are smaller still; these describe mass concentration without supplying independent samples. A mass-weighted 99th percentile can consequently be set by one heavy particle. Four realizations provide limited information about such rare outcomes, and particles sharing collisions and vapour are not independent uncertainty samples.

Finite quantile sampling cannot by itself account for the complete comparison: the mean mass CDF changes by 46.13 and 28.99 percentage points, and the median rises from 12.27 through 18.53 to 23.46 cm. The arriving number flux simultaneously falls from about 16,215 through 8,249 to 5,852 particles/m²/h while arriving water flux increases. The ratio of a larger mass flux to a smaller number flux raises the mean mass per arriving particle; additional uptake alone does not account for the size shift.

Halving the macrostep from 60 to 30 s at the middle volume changes deposition from 2.37250 to 2.36791 mm/h (−0.19%) and the lowest-layer empty fraction from 40.74% to 39.84%. The frozen-history allocation gap changes from 0.26959 to 0.27354 mm/h. These diagnostic averages are much less sensitive than the volume comparison. The size distribution still fails the existing timestep screen: its mean mass CDF changes by 6.45 percentage points and its 99th percentile by 5.63%, despite a 0.12% median change. Four realizations do not distinguish every small timestep difference from sampling variation; this failed control remains an outstanding numerical limit.

## What the diagnosis changes

The next numerical target is a representation whose collision population and local water uptake are stable together as statistical sampling improves. The controls justify prioritizing the sparse population created by sticking and the way its fluctuating demand meets the supply cap. They do not identify a universal upper diameter, establish mathematical gelation, or show that a missing fracture law is the cause of this numerical volume dependence.

A useful next diagnostic would separate the volume used to sample collision partners from the volume used to estimate capped vapour demand, with an explicit shared water ledger. That would test the uptake mechanism while allowing the subsequent population to respond, which the frozen-history calculation cannot do. It is a proposed diagnostic intervention: sharing supply changes the finite-volume model and must not silently become a physical water-transport assumption. A practical replacement for the small direct-particle population must also reproduce the existing event controls and pass larger-volume, time and distribution comparisons.

Aggregate geometry, contact failure and climate-supported moist histories remain additional physical questions after numerical convergence. The present study leaves all earlier snow populations and optical products unchanged. Its new products preserve the failed comparisons and impose no flake-size cap.

## Verification and records

All 48 full columns completed with at most two numerical workers. The largest absolute relative water residual is 1.57 × 10⁻¹⁵; the particle-number residual is zero. The observer reproduces all sixteen original baseline and timestep realizations exactly, and its integrated uptake and event counts match the solver ledgers. The earlier 122 scientific files retain their recorded hashes. The primary matrix and additional eight-case product identify their actual source files, forcing inputs, compiler and binary, settings and complete individual outcomes.

The focused numerical and reader suite passed all 39 tests. After inspection exposed a value label crossing the CDF threshold line, the renderer's label position was corrected and its nine reader tests passed again. Two final builds produced byte-identical PNG, SVG and manifest files. The final 2448 × 1908 PNG was inspected at original resolution for labels, water totals, scales, uncertainties and the distinction between the CDF screen and overall acceptance. The repository layer and provenance checks are recorded separately in the [verification record](results/column_diagnosis_checks.json). The full repository `make check` was not repeated for this continuation; its preceding result remains in the original event-column record.

Substantive review compared the observation weights and unchanged update order with the frozen driver, traced the water identities through the outputs, and read the process-removal and averaging arguments together. No new external dataset, physical closure, full scholarly source review or empirical validation was added. Work remains uncommitted on `domain/atmosphere`. Reproduction begins with the [methods](column_diagnosis_methods.md), followed by:

```sh
python -m atmosphere.snow.column_diagnosis_analysis
python -m visualization.snow.column_diagnosis
python -m pytest -q atmosphere/tests/test_snow_column_diagnosis.py atmosphere/tests/test_snow_column_diagnosis_runs.py atmosphere/tests/test_snow_column_diagnosis_analysis.py visualization/snow/tests/test_column_diagnosis.py
```

**Provenance, 25 October 2026.** [results/column_diagnosis.json](results/column_diagnosis.json),
[results/column_diagnosis_aggregation_only.json](results/column_diagnosis_aggregation_only.json) and
[results/column_diagnosis_checks.json](results/column_diagnosis_checks.json) were built with the rain solver as it stood
before 25 October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998) measured value. They keep
their values; their record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for byte as
[atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and their identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
