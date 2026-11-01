# Moist exposure in archived lunar storm air

[moist_history.py](moist_history.py) follows **fixed-drift kinematic tracers** through an existing three-dimensional CM1 storm sequence. It calculates the temperature, ice-relative humidity and cloud-liquid exposure along each resolved path. These histories describe a selected weather episode under the older `A28_dim5` forcing; they do not supply a calibrated growing-flake population or the separate 301 K resort climate.

The distinction matters because a deep snow-bearing layer can contain dry air. [The earlier snapshot diagnostic](turbulence.md) established that difference spatially. The new calculation adds time-dependent advection, retaining the limits imposed by the saved fields. A snapshot's contiguous moist depth never becomes a residence time by itself.

## Available fields and sampling

The archive contains two corrected electrified-box storm windows, each with 216 complete three-dimensional snapshots at ten-minute intervals. The second window covers model days 40.50694–42.0. Its 64 × 64 horizontal grid has 6.011 km spacing; the selected winter air has approximately 2 km vertical spacing. Potential temperature, pressure, vapour, liquid, ice, snow and all three scalar-point wind components are available. The box is stationary and periodic horizontally.

The calculation uses its first 211 frames, a common **35-hour interval** ending at day 41.96528. This makes the 10-, 20- and 30-minute field sequences end at exactly the same time. The final five available frames are retained in the inventory but excluded from every control. Trajectories stop at the selected center-height boundaries, 20.409 and 59.794 km, or at the end of that time interval. A height exit is a numerical analysis boundary, not melting, arrival at the ground or a physical flake's death.

The full-cycle Morrison and electrified boxes supply three-hourly fields. Corrected two-dimensional rings also have three-hourly output. A separate 2 km horizontal box has fifteen-minute output, but its initialization, domain and configuration differ; it supplies no matched spatial refinement. The inspected ordinary GCM archive contains ten vertical layers and approximately three-day means. The instantaneous GCM wave products have three-hour cadence and export the lowest level, without the upper-air vertical velocity and cloud water needed here. [The source register](moist_history_sources.json) records the inspected headers, files and hashes; none of these other cases is silently substituted into the storm sequence.

## Conditional release and transport calculation

At 0, 6, 12, 18, 24 and 30 hours, the diagnostic samples at most 64 cells satisfying all three conditions: −20 to 0 °C, snow mixing ratio at least 10⁻⁵ kg/kg dry air, and ice-relative humidity at least one. A fixed random seed selects cells without replacement. Every selected position supplies four tracers, with prescribed downward drift speeds of 0, 0.3, 0.6 and 1 m/s. Their initial positions and release times remain identical across numerical controls.

The drift classes probe transport. No diameter, mass, habit or growth law is assigned to them. They do not sublimate in dry air, melt in warm air, collide or fragment. Consequently their continued presence cannot establish that a real snowflake would survive the same path. Zero drift samples the resolved wind; external nudging, imposed large-scale advection tendencies and unresolved mixing also prevent interpreting that path as a conserved material label.

Wind and thermodynamic fields use trilinear interpolation on scalar grid points and linear interpolation between saved times. Temperature and ice humidity are recomputed from interpolated potential temperature, pressure and vapour, using the existing CM1 thermodynamic constants and Murphy–Koop ice-saturation relation. An explicit midpoint step advances position. Midpoint quadrature accumulates cold, ice-saturated, near-saturated (`RHi ≥ 0.95`), dry (`RHi < 0.95`), snow-bearing and cloud-liquid exposure. All exposure masks include the −20 to 0 °C interval, which is wider than the earlier ideal −15 to 0 °C growth column.

The result also retains integrated humidity excess and deficit, integrated cloud liquid, sampled temperature/humidity extrema, and 24 illustrative temperature/humidity/cloud histories. First and longest uninterrupted moist episodes are diagnosed at the integration cadence. A separate flag marks a first episode still ongoing when the history is censored; another identifies an ongoing final episode. Crossing an analysis height within a step uses a linear segment approximation, tested by halving the integration step.

The six releases have **5–35 hours of available follow-up**. A release is a selected position, not an inferred flake birth or the onset of a moist cloud; all episode durations begin at release and recover no earlier exposure. Pooled exposure quantiles describe these selected, censored histories, with equal sampled-cell weights. They are neither lifetime quantiles nor particle-number, snowfall-mass or weather-frequency distributions. Per-tracer durations, releases, drift classes and censor flags are retained so readers can make matched comparisons without discarding incomplete histories.

## Numerical controls and evidence

The baseline uses 60-second steps and ten-minute fields. Three independent controls halve the step to 30 seconds, or retain 60-second steps while sampling fields every 20 or 30 minutes. These distinguish integration sensitivity from sensitivity to the available temporal sampling. Decimating an archive cannot recover fluctuations missing from the finest saved cadence, and no matched spatial-refinement experiment is available.

The result records the measured field-to-field displacement scale `|w| Δt / Δz` and its horizontal counterpart. A displacement exceeding a grid layer does not, by itself, invalidate interpolation; it identifies conditions where the saved fields poorly constrain intervening cloud structure. Consecutive moist samples at a fixed grid cell are retained separately as Eulerian spells. Their first-to-last sampled span is not a lower bound on continuous exposure: unresolved dry gaps can intervene. No positive continuous-duration lower bound follows without a temporal regularity assumption.

## Exposure in the selected storm window

The [saved result](results/moist_history.json) contains 354 release cells and 1,416 histories: 64 cells at five release times, and all 34 qualifying cells at hour 24. In the baseline, median accumulated ice-saturated exposure is 26–30.5 minutes across the four drift classes. The corresponding median cold exposure is 3.33–9.54 hours, so remaining within the selected winter temperatures does not imply continued access to saturated air.

| Prescribed downward drift | Median tracked time | Median first moist segment after release | Median / 90th-percentile cumulative moist exposure | Median dry winter exposure | 90th-percentile cloud-liquid exposure |
|---|---:|---:|---:|---:|---:|
| 0 m/s | 17.00 h | 19 min | 29.5 / 130.7 min | 4.31 h | 28.0 min |
| 0.3 m/s | 12.53 h | 20 min | 30.5 / 114.1 min | 6.50 h | 30.7 min |
| 0.6 m/s | 7.29 h | 20 min | 27.0 / 105.7 min | 3.25 h | 32.0 min |
| 1 m/s | 5.00 h | 19 min | 26.0 / 91.7 min | 2.22 h | 36.0 min |

Here moist means `RHi ≥ 1`; dry means `RHi < 0.95`; cloud liquid means at least 0.01 g/kg dry air. Each also requires −20 to 0 °C. The intermediate humidity band belongs to neither the moist nor the dry total. About 44–46% of histories encounter cloud liquid above the threshold, giving a zero median liquid-exposure duration in every class. These sampled-cell frequencies depend on the conditional releases and fixed drift.

The first segments are episodes in the interpolated resolved fields. Their 19–20-minute medians span only about two archived intervals. All first segments end before the history is censored in the baseline; their durations begin at release rather than at cloud entry or flake birth. A recorded zero segment means that the first integration midpoint already fails the moist mask, despite a qualifying release cell. Later re-entry accounts for the difference between first-segment and cumulative exposure. The full histories end at the time-window boundary in 561 cases, the lower analysis height in 705, and the upper height in 150. Only two are still moist at that endpoint. The zero-drift class has 298 time-window exits and 56 upper exits; neither count establishes a completed air-parcel lifetime. The pooled cumulative durations remain conditional on the unequal 5–35-hour follow-up and the chosen height domain.

The numerical controls show a larger effect from field sampling than from the integration step:

| Change from baseline | Paired cumulative-moist difference, median / 90th / maximum | Paired first-segment difference, 90th / maximum | Histories with a different terminal boundary category |
|---|---:|---:|---:|
| 60 to 30 s integration step | 0.5 / 0.5 / 8.5 min | 0.5 / 0.5 min | 0 / 1,416 |
| 10 to 20 min fields | 2.0 / 10.0 / 137 min | 6 / 100 min | 32 / 1,416 |
| 10 to 30 min fields | 4.0 / 23.0 / 192 min | 18 / 201 min | 70 / 1,416 |

These are absolute differences for the same release and drift, not differences between unpaired quantiles. The longest uninterrupted segment changes by at most one minute when the integration step halves. Coarser fields can alter vertical transport sufficiently to move a history from an early height exit to the end of the window; the largest resulting difference in tracked duration is about 34.4 hours. Thus the stable central exposure values do not establish accurate individual paths or the tail of residence time. The control addresses sensitivity to losing saved frames, while variability faster than ten minutes and finer than the grid remains unmeasured.

The sampling limitation is strongest in liquid-bearing winter snow. At those 9,345 sampled cell-times, the median `|w| Δt / Δz` is 1.31, its 90th percentile is 3.01, and 61.4% exceed one layer per ten-minute frame. The corresponding exceedance is 22.0% in moist snow and 3.18% in all sampled winter snow. These are resolved-air displacement scales, not a solver stability condition or a measure of unresolved turbulence. They explain why accurate integration through interpolated fields cannot by itself recover the cloud changes between saved times.

This calculation supplies a conditional set of temperature, humidity and cloud-exposure histories for testing growth assumptions. Converting them to flake sizes would additionally require evolving mass, area and fall speed, deposition and sublimation, liquid collection, and mechanically credible collision outcomes. Neither the 35-hour storm-window duration nor the longer idealized column stress tests can be assigned as uninterrupted growth time to these tracers.

## Reproduction and verification

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.moist_history
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest atmosphere/tests/test_snow_moist_history.py -q
```

The external archives must be available through `climate/crm/runs`. The producer copies only selected float32 fields into immutable memory, records their hashes, and verifies control-file hashes and snapshot file metadata afterward. Full climate fields remain external. All eight focused tests pass: seven analytical controls check periodic interpolation, thermodynamic conversion, constant-flow exposure, height and episode censoring, release identity and temporal interpolation; the saved-product check verifies source and constant identity, matched releases, finite arrays and exposure accounting. A second reader reviewed the transport and censoring logic. These checks establish implementation behavior; higher-resolution meteorology, evolving particle properties and empirical constraints remain separate requirements.
