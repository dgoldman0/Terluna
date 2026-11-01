# Crystal growth in archived lunar humidity

With the vapour, temperature, pressure and winds from an archived lunar storm, the improved McSnow control gives **0.724 mm median peak dimension** from 4 µm spherical seeds, compared with **5.986 mm** when water saturation is maintained along the evolving paths. Of 354 releases, **253 (71.5%) subsequently reach the tiny sublimation residual**, 74 leave the supported temperature interval toward colder air, one reaches freezing, four reach the declared diameter censor and 22 remain in the archived window when it ends. These are peaks during conditional individual histories; they do not describe snow arriving at the ground. Intermittent moist air and subsequent sublimation substantially reduce the growth sustained by continuous saturation.

The [numerical product](results/mcsnow_weather.json) records twenty executions, each following the same 354 release cells. Each crystal samples prescribed weather independently, without consuming its vapour or heating its air. This tests the improved component under the modelled lunar humidity rather than assigning every crystal a share of an assumed precipitation allowance. A population calculation still needs the incoming number and habit distributions, depletion, latent heating, collisions and breakup. The preceding [ideal-column comparison](mcsnow_lunar_growth.md) remains a separate experiment: its 3 mm/h was an imposed allowance, not a measured lunar water supply.

## What water is present

The source is CM1 case `box_0e_elec_corrected_storms_second`, forced by **`A28_dim5`**, over a 385 km periodic all-land square near 0° N, 0° E. The calculation uses 211 archived fields at ten-minute intervals, spanning 35 hours from model day 40.506944. It retains 21 vertical centres between 20.409 and 59.794 km, with roughly 2 km spacing around the snow region. This older storm archive resolves short moist episodes more usefully than the newer design-climate cases saved every three hours; it is not relabelled as the newer 301 K design climate. The existing [moist-history study](moist_history.md) supplies the source identities, interpolation and release cells. Release selection retains that study's saturation rule; the new growth and reported humidity diagnostics use native McSnow saturation.

Summing the archived vapour in cells whose centres lie between −20 and 0 °C gives a median **25.03 kg/m²**, equivalent to **25.03 mm of liquid water**, across all horizontal columns and saved times. The 10th–90th percentile range is **19.26–31.79 mm**. These are instantaneous stocks through the cold layers, not delivery rates. The corresponding median positive vapour excess above ice saturation is zero, with a 99th percentile of 0.721 mm; most columns therefore cannot sustain deposition throughout their cold depth without cooling or moisture transport.

Among cells already containing snow (`qs ≥ 10⁻⁵ kg/kg`), median ice-relative humidity is **89.47%**, with a 10th–90th percentile range of **74.58–104.34%**. Ice can persist temporarily in these cells while sublimating. The selected releases begin in moist snowy cells between −20 and 0 °C at six times, then move with resolved wind and their evolving native settling speeds. Counts across these selected cells are conditional diagnostic frequencies, not estimates of the storm's particle population.

## Humidity and starting-crystal comparisons

All table values use ten-second trajectory steps and ten-minute meteorological fields. The growth solver limits each native step to 0.5 seconds and 0.5% of the preceding crystal mass. Each row has 354 releases. “Peak” means the largest dimension attained before that path's exit; the 90th percentile is over all selected releases, retaining extreme outcomes.

| Seed diameter | Humidity treatment | Median peak, mm | 90th-percentile peak, mm | Tiny residual | 1 m censor |
|---|---|---:|---:|---:|---:|
| 4 µm | Archived vapour | 0.724 | 4.928 | 253/354 | 4/354 |
| 4 µm | Continuous water saturation | 5.986 | 90.387 | 0/354 | 20/354 |
| 4 µm | Archived vapour with water-saturation ceiling | 0.599 | 3.446 | 256/354 | 0/354 |
| 0.2 mm | Archived vapour | 0.360 | 1.194 | 260/354 | 0/354 |
| 0.2 mm | Continuous water saturation | 1.490 | 4.219 | 0/354 | 0/354 |
| 0.2 mm | Archived vapour with water-saturation ceiling | 0.355 | 0.816 | 260/354 | 0/354 |

The saturation control sets vapour pressure to the native water-saturation value at the local temperature throughout each path. The condensation-ceiling sensitivity instead uses the lesser of archived vapour pressure and water saturation, preserving dry periods while suppressing water supersaturation. Neither operation solves condensation or adds a particle population. Growth changes settling and hence later exposure, so the comparisons retain common release cells while allowing their subsequent paths to diverge.

Starting from compact **0.2 mm spheres** gives a median peak of **0.360 mm** in archived humidity, compared with **1.490 mm** in its saturation control. 260 of 354 reach the small residual, 2 reach freezing, 71 leave toward colder air, 21 reach the archive end and 0 reach the diameter censor. The starting mass and shape history therefore remain consequential even when both seed sizes sample the same release cells.

For the small seeds in archived humidity, median cumulative time at or above ice saturation is **41.99 minutes**, within a median tracked age of **1.056 hours**. The positive growth and subsequent mass loss are recorded separately. Reaching the residual ends a path with its small positive mass retained; it does not set the mass to zero or represent a simulated nucleation/recycling process.

## Extreme dimensions and weather resolution

Four small-seed paths reach the **1 m numerical censor** under archived vapour. Their final dimensions are approximately 1 m long by 0.53–1.32 mm wide, with aspect ratios of roughly 760–1,896. These highly elongated columns expose the continued extrapolation of the native habit and geometry law. They remain flagged model outcomes whose physical shape evolution and mechanical survival require validation. Removing the extra oblate mean-ventilation term in the preceding Earth work did not validate every habit under changing temperature and humidity.

Substantial supersaturation is present in the source fields before interpolation: **6.87%** of snowy −20 to 0 °C cells exceed water saturation, with a maximum water-relative humidity of **171.0%**. Interpolation and sparse saving can further alter the forcing seen by a trajectory; this audit establishes that interpolation alone did not create all the excess. The condensation ceiling removes all four small-seed metre censors and lowers the 99th-percentile peak from **839.6 to 9.67 mm**. Its largest remaining peak is still **157.5 mm**, so that sensitivity does not establish a physically valid extreme tail.

Halving trajectory steps from 20 to 10 seconds changes the median and 90th-percentile peaks by at most **0.357%** across the six principal seed/humidity combinations. All paired paths retain their termination category. Halving both the maximum native growth step and allowed fractional mass change alters the archived-humidity median and 90th percentile by at most **0.040%**. Reducing the mass-equivalent residual diameter from 0.2 to 0.1 µm preserves every termination category and leaves every peak unchanged. One extreme path changes from 322.8 to 231.4 mm under trajectory refinement, a 39.5% difference relative to the finer value; native growth refinement changes another peak from 696.8 to 726.4 mm, or 4.24%. Stable central quantiles therefore coexist with numerically sensitive extreme paths. These are measured numerical sensitivities, not empirical error bounds.

The meteorological cadence matters more than the central timestep change. Sampling the archive only every 20 and 30 minutes changes the small-seed median peak from **0.724 to 0.699 and 0.605 mm**; the 90th percentile changes from **4.928 to 4.759 and 4.198 mm**. These experiments omit saved information and cannot establish convergence to unresolved weather. The ten-minute result is conditional on the available spatial and temporal resolution.

## Native reuse, bookkeeping and reproduction

The [batch interface](mcsnow_weather_driver.f90) calls the pinned 2024 McSnow deposition, geometry and terminal-speed routines at revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`. The sole physical growth intervention remains removal of the extra oblate mean-ventilation term; directional growth, density and drag are retained, with lunar gravity **1.6242188766 m/s²**. The unused steady-column supply gate stays off. CM1 potential temperature and dry-air vapour mixing ratio are converted with the shared CM1 constants, while the native component retains its identified air constants and saturation equations. Pinned source archives and earlier scientific products remain intact.

The [trajectory runner](mcsnow_weather.py) uses midpoint motion with periodic horizontal interpolation, evolving terminal speed, adaptive native growth and explicit partial-step exits. It stops on reaching freezing or crossing below −30 °C, rather than continuing into melting or extending the habit table. It also stops at selected height boundaries, the archive end, the 1 m diameter censor or the mass of a dense 0.2 µm sphere. A 0.1 µm residual control tests the latter choice. The spheres' release histories are recorded with the humidity treatment actually selected.

The [analysis](mcsnow_weather_analysis.py) rereads raw trajectories and verifies release identity, elapsed time, axes, volume, density, aspect ratio, mass gained and lost, retained peaks and every exit. Native calls preserve signed latent heat, but air temperature is prescribed: that ledger reports heat exchanged with imposed air, not atmospheric energy closure. Two constant-air parity exercises reproduce the retained Earth and lunar histories exactly over 10 and 30 saved points. An independent constant-settling control recovers the analytic distance/speed travel time, and known pressure/mixing-ratio inputs recover the intended partial pressure and temperature.

The ignored cache is `research/runs/snow_scale/mcsnow_weather/`, on the project-space drive. It holds the selected weather bytes, per-frame hashes, build/compiler/table identities, each executed producer and request, raw arrays, run summaries and logs. The compact product contains the source and run identities. The `matrix` command requires an empty `runs/` directory; retain existing runs and their logs elsewhere before repeating it. The analysis command reassembles the product from retained outputs. From the repository root, with the external native sources and archived CM1 files installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_weather build
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_weather weather
# Run in a fresh run directory; completed runs are never overwritten.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_weather matrix
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_weather_analysis
```

All **20 executions / 7,080 conditional histories** pass the raw audit. The maximum particle-mass ledger residual, normalized by initial mass plus total mass gained and lost, is **2.88e-15**. The checkpoint passes `make check`: **3,697 Python tests passed, 20 skipped; 191 JavaScript tests passed; 1,317 tracked files have zero layer violations**. Historical provenance and ensemble integrity pass. The [check record](results/mcsnow_weather_checks.json) separates source identity, raw validation, numerical controls and substantive report reading from empirical validation.

These results support substantial growth during moist episodes followed by widespread sublimation in this archived storm. They narrow the ideal-column extrapolation without supplying a validated absolute snow-size distribution. A coupled population and better constrained crystal morphology remain necessary for snowfall rates, aggregate sizes and storm optics.
