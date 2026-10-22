# Crystal growth across habits and a ventilation control

Native SCALE's discrepancy changes with crystal habit. In five Earth laboratory comparisons, it produces columns that are too short and plates that are too wide relative to the published fits, while all five endpoint masses are high. The −14.4 °C dendritic case remains the strongest discrepancy. Removing only the directional ventilation correction reduces its 30-minute width from 27.66 to 15.81 mm, against the fitted 4.17 mm. This identifies a contributing mechanism without establishing a corrected growth law.

The [new data and model review](growth_options.md) adds 169 recovered laboratory measurements, a recently published aggregate-geometry alternative and environmental controls for eventual lunar transfer. These findings leave the original full-column result intact as a numerical comparison: under matched forcing, SCALE produces larger lunar than Earth flakes. They further restrict interpreting either absolute size as a physical forecast.

## Matched laboratory comparisons

The [runner](scale_growth_benchmarks.py) extends the preceding isolated-crystal experiment at five temperatures in Takahashi et al.'s Table 2. Each starts at the three-minute fitted mass and full a/c dimensions, at 1010 hPa and water saturation, and ends at that row's last measured age. Eight identical particles remain suspended in fixed air. Native terminal-speed ventilation acts during growth; injection, collisions, liquid water, melting and transport are absent. Both 5 s and 1 s steps are run. These are prescribed crystal-growth experiments, not reproductions of a snowstorm. [Takahashi et al. (1991)](https://www.jstage.jst.go.jp/article/jmsj1965/69/1/69_1_15/_article/-char/en).

| Temperature, °C | Final age, min | Fitted maximum dimension, mm | Native dimension, mm | Native/fitted mass | Native/fitted fall speed |
|---:|---:|---:|---:|---:|---:|
| −3.7 | 20 | 0.1230 | 0.1864 | 1.94 | 1.59 |
| −5.3 | 30 | 1.8611 | 1.1397 | 1.41 | 2.11 |
| −10.6 | 12 | 0.2011 | 0.3176 | 1.99 | 1.23 |
| −14.4 | 30 | 4.1745 | 27.6620 | 19.10 | 0.80 |
| −20.1 | 12 | 0.2066 | 0.2758 | 2.02 | 1.83 |

Native entries use the 1 s step. Maximum dimension is max(a,c), preserving the long polar axis for columns; planar width and thickness remain separate in the [numerical product](results/scale_growth_benchmarks.json). Fall-speed fits come from Table 4. The initial output's speed precedes the imposed air reset and is excluded from that comparison. At −14.4 °C the endpoint speed is actually below the measured fit despite the excessive width, so an excessively fast endpoint fall speed alone does not explain the error.

Five other Table 2 fits cannot be initialized unchanged at three minutes: their mass and axes imply spheroid densities of 1020–1267 kg/m³, above SCALE's 916.8 kg/m³ bulk-ice value. The runner rejects those mappings explicitly. Independent fitted curves and real nonspheroidal crystals need not fit this representation; rejection does not identify erroneous observations. No density clipping or invented thickness was used to admit them.

The measured fits lack pointwise error estimates, and thin-crystal thickness was inferred from mass and basal area. Their piecewise joins are preserved rather than smoothed. The new [measurement product](results/takahashi_measurements.json) retains nearby individual observations separately from those fits and flags unresolved data issues. It neither supplies confidence intervals nor converts a comparison against fitted initial geometry into an exact replay of individual crystals.

## Isolating directional ventilation

A [separate diagnostic build](scale_growth_controls.py) records the native terms at every growth step. Mode 0 preserves the native calculation; mode 1 sets only the polar/equatorial ventilation-factor ratio to one. The mean ventilation expression, capacitance, deposition density, intrinsic growth function and terminal-speed law are unchanged. Their evaluated values can change as the crystal's geometry changes. Each mode is run at 5 s and 1 s.

| At age 30 min, −14.4 °C | Width, mm | Mass, mg | Polar dimension, μm | Integrated polar/equatorial logarithmic growth ratio |
|---|---:|---:|---:|---:|
| Published fit | 4.175 | 0.07286 | 25.80 | 0.5092 |
| Native, 1 s | 27.662 | 1.39186 | 12.55 | 0.1123 |
| Unit directional ventilation ratio, 1 s | 15.814 | 0.79774 | 22.00 | 0.2818 |

The native intrinsic ratio is 0.2818383. Its directional ventilation factor decreases from 0.732 to 0.137, reducing the effective growth ratio from 0.206 to 0.0387. Thus the native calculation increasingly favours lateral growth. Fixing that factor at one prevents this additional preference, cutting width by 43% and mass by 43%, while leaving width 3.79 times and mass 10.95 times the fitted values. The intrinsic ratio itself remains below the observed integrated ratio. This is direct evidence for one contribution within a coupled habit/geometry/uptake discrepancy, rather than a unique identification of a replacement formula.

Deposition density stays at 276.85 kg/m³ in both runs. The native mean ventilation factor rises from 1.17 to 5.96; the controlled run reaches 5.22. The law for mean ventilation therefore remains active. An arbitrary dimension ceiling would hide these interactions without repairing them, and none has been added. [Step traces and independent checks](results/scale_growth_verification.json) retain the actual terms, source identity and snapshot comparisons.

## Numerical and source checks

Across the five native temperatures, 5→1 s refinement changes endpoint dimension by at most 0.81% and mass by at most 1.48%; the maximum sampled trajectory changes are 1.18% and 2.29%. For the unit-ratio control, endpoint changes are 0.39% and 0.90%, with trajectory maxima 0.93% and 1.96%. These measured refinements are much smaller than the principal empirical discrepancies; an asymptotic convergence rate was not established.

Both new −14.4 °C baseline runs reproduce the corresponding earlier water-ledger bytes and all 473 native particle/history arrays exactly. Both instrumented mode-zero runs independently reproduce the new native baselines the same way. All 14 new executions have finite positive ice states, unchanged populations, reproduced physics inputs and identified executables. Four complete step traces agree with the saved particle mass, axes and density. Ice deposition closes within 7.7 × 10⁻¹⁵ relative error; no run uses more than 6.4 × 10⁻⁷ of the generous numerical deposition ceiling. The fixed Murphy–Koop water saturation appears as 100.044–100.071% on SCALE's convention, explicitly preserving that small difference.

The required full `make check` passed again for this growth checkpoint: 3271 Python tests, 20 skipped, 191 JavaScript tests, 1259 files with no layer violations, and provenance/ensemble integrity checks. The focused suite passes 80 tests. The [check record](results/growth_benchmark_checks.json) retains both the preceding `002d7b6` check and this completed check, alongside the native, data and source verification. External full papers remain in the ignored cache. Numerical/build concurrency was at most two workers, one thread each.

```sh
# Environment with NumPy and netCDF4; native builds/data must be present at recorded hashes:
python -B -m atmosphere.snow.scale_growth_controls build
python -B -m atmosphere.snow.scale_growth_benchmarks
python -B -m atmosphere.snow.scale_growth_benchmarks --controls
python -B -m atmosphere.snow.scale_growth_verify
# Recover the public observations from the identified cached workbook:
python -B -m atmosphere.snow.takahashi_measurements research/runs/snow_scale/growth_environment_sources/1991_takahashi.xlsx > atmosphere/snow/results/takahashi_measurements.json
```

The next repair comparison should use an existing habit-aware growth treatment, with identical initial geometry and air conditions, against measured mass and both axes across temperatures. McSnow's published growth alternatives supply such a candidate; its 2026 aggregate geometry can then address contact arrangements. Low-pressure, changing-humidity and natural-cloud evidence should test transfer beyond the saturated laboratory cases. A calibrated lunar size distribution, aggregate fracture law and storm optical forecast remain unresolved.
