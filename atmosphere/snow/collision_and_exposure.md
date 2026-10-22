# Collision spreading and moist growth histories

The collision diagnosis identifies a concrete numerical error. Under the same physical collision law, the fixed-pivot solver turns repeated capture of small particles into an excessively broad collector distribution. An independent calculation reproduces the mean growth while keeping the distribution much narrower. Separately, paths through the archived storm show brief moist exposure for many selected tracers. Together these results strengthen the case for revising both collision transport and the prescribed growth environment before interpreting the earlier 15–60 cm tail.

This work extends the [growth-transport repair](growth_repair.md). It preserves all earlier snow producers and saved products. Numerical controls, conditional weather histories and the published Earth comparison supply different kinds of evidence; none establishes a calibrated lunar flake-size distribution.

The [volume-dependence diagnosis](column_diagnosis_results.md) completes 48 full lunar columns. It measures reduced vapour uptake from a sparse population whose demand fluctuates across the local supply cap: deposition rises from 2.199 to 2.498 mm/h across the three statistical volumes under the same 3 mm/h allowance. Removing all collisions passes the tested volume screen. Removing growth, chipping, or both leaves failed upper-tail checks; the sticking-only 99th percentile changes by 22.5%, with appreciable four-realization uncertainty. This identifies a water-allocation contribution and a separate tail concern, while the coupled distribution and physical lunar sizes remain unresolved.

The [full-column event study](event_column_results.md) now completes 512 independent closed-box runs and 28 coupled Earth/lunar controls. The event solver passes the reference comparison and closes the water and number ledgers, while the lunar distribution fails volume convergence: mean mass medians are 12.3, 18.5, 23.5 cm and mean 99th percentiles are 19.8, 36.3, 74.1 cm across successive fourfold volume increases. These are numerical diagnostics under the inherited growth, geometry and chipping laws. They leave the physical 15–60 cm tail and predictive storm optics unresolved.

## The collision error is now measured independently

The [collision reference](collision_reference.md) compares the frozen snow kernel with exact-event, continuous-mass calculations. In its clearest control, an imposed 15.86 cm collector captures 0.209 mm particles from a maintained bath of 100 per litre. Deposition, settling and breakup are disabled. After 4.54 minutes, mean mass has increased by about 10.55% in both calculations, but their widths differ greatly:

| Collision representation | Mean mass / initial mass | Standard deviation / initial mass |
|---|---:|---:|
| Exact events, 8,192 independent histories | 1.105480 | 0.0629% |
| 96 mass pivots | 1.105511 | 19.395% |
| 191 mass pivots | 1.105492 | 13.197% |
| 381 mass pivots | 1.105483 | 9.158% |

At 96 pivots, the variance is about 95,125 times the exact-event sample variance while the means differ by only 0.00279%. Grid refinement reduces the spreading but leaves a large discrepancy. The grid calculation uses a matrix exponential here, eliminating time-step integration as an explanation. Its number and bath-water ledgers close, so accurate conservation also fails to detect this distribution error.

The mechanism is repeated mass allocation. A captured particle is far smaller than the distance to the next collector mass pivot. Allocating the merged mass between neighboring pivots preserves its mean but assigns a small fraction to a disproportionately large jump. Repeating that allocation broadens the distribution far beyond the randomness of actual captures. The [comparison figure](../../visualization/snow/generated/collision_comparison.png) shows the native number CDFs, sampling bands and grid dependence.

This control begins with a large collector. It measures the representation of its subsequent growth, not the formation or survival of a 15 cm flake. In ordinary closed boxes starting from small seeds, six-hour bulk moments approach the largest-volume exact-event result as the grid is refined. Rare tails remain uncertain. At twelve hours both the finite-particle reference and the mass-grid result remain sensitive to resolution. The direct event solver therefore supplies a diagnostic target, rather than a converged replacement for every full-column condition.

## The available weather does not establish prolonged moist residence

The [moist-history calculation](moist_history.md) follows 1,416 fixed-drift tracers through 35 hours of ten-minute CM1 output. Releases are selected snowy, ice-saturated cells between −20 and 0 °C at six times. Each position is followed with four prescribed downward speeds. The 354 positions and their release times remain identical across the numerical controls.

| Downward drift | Median cumulative moist exposure | 90th percentile | Median first moist segment after release |
|---|---:|---:|---:|
| 0 m/s | 29.5 min | 130.7 min | 19 min |
| 0.3 m/s | 30.5 min | 114.1 min | 20 min |
| 0.6 m/s | 27 min | 105.7 min | 20 min |
| 1 m/s | 26 min | 91.7 min | 19 min |

Moist exposure here requires both ice saturation and the stated temperature interval. These quantiles pool selected histories with equal sampled-cell weights and 5–35 hours of available follow-up. Releases are neither flake births nor inferred cloud entry; no pre-release history is recovered. Tracers retain their prescribed drift and survive dry air mathematically, so their trajectories do not establish that real flakes survive. The forcing is the older `A28_dim5` storm, distinct from the 301 K ideal snow column.

Halving the trajectory step from 60 to 30 seconds changes cumulative moist exposure by at most 8.5 minutes, with a paired 90th-percentile difference of 0.5 minute and no change in terminal censor category. Coarsening saved fields from ten to twenty or thirty minutes gives paired 90th-percentile differences of 10 and 23 minutes; 32 and 70 histories change their height-exit/time-limit category. Some tracking durations change by more than 34 hours. These are material meteorological sampling limits, and the first moist segment typically spans only about two saved intervals. Finer trajectory integration cannot recover missing cloud evolution.

The result consequently challenges the use of a continuously moist growth path as a representative default in this archived storm. It does not exclude rare prolonged exposure, establish lunar weather frequencies or determine a maximum aggregate size. A future growing-particle calculation must account for sublimation, changing settling speed and re-entry into moist regions, with adequately sampled weather.

## Published Earth aggregation: a remaining discrepancy

The [Earth aggregation component adapter](welss_aggregation_physics.md) matches 182 independently compiled cases from pinned author routines. It retains the monomer/aggregate capacitance distinction, rough-ice hydrodynamic interception, the source sticking intervals and geometric encounter area. Its source-level agreement does not establish a coupled population distribution.

The separate [author-model comparison](welss_aggregation.md) builds the public McSnow source and executes ten-hour deposition/settling and aggregation cases. The only source intervention makes the RNG seed explicit. The final comparisons use 31 raw profile snapshots from hours five through ten, inclusive. Integrated absolute differences against the archived author profiles are:

| Declared author rerun | Number density | Mass density | Number velocity moment | Mass velocity moment |
|---|---:|---:|---:|---:|
| Deposition and settling | 2.076% | 0.828% | 0.817% | 0.862% |
| With aggregation | 12.521% | 16.482% | 14.651% | 15.251% |

The aggregation case has a material 12.5–16.5% discrepancy, concentrated below the nucleation layer. Using the final four snapshots leaves total-profile differences of 12.5–16.8%, so that averaging alternative does not explain it. Separating particle categories reveals larger disagreements: the aggregate number-density profile differs by 70.11% and its mass velocity moment by 67.64%. Totals alone therefore understate the mismatch in the population composition. These differences are sums of absolute bin differences divided by the summed archive profile over its populated bins below 3 km. Regional, monomer/aggregate and alternative-averaging comparisons remain in the products and Earth report. The [profile figure](../../visualization/snow/generated/earth_comparison.png) shows the actual 20 m bins, with matching scales between cases.

This is an author-solver rerun under declared settings. The public source is from 2024, while the archived profiles were processed in 2023; their executable, seed and exact averaging times remain unidentified. One seed and one resolution do not establish stochastic or numerical convergence. The earlier independent deposition/settling characteristics remain a separate comparison, and an independent coupled aggregation column has not yet been tested against these profiles. Neither result is empirical validation of giant lunar aggregates, and none of these Earth closures has been adopted in the frozen lunar column.

## What a repair must demonstrate

The collision operator must represent many small captures without creating the rare large jumps measured above. A candidate continuous-mass or mixed transport/collision method can be checked first against the independent bath distribution, including its physical capture variance, then against the ordinary closed-box controls. Number changes, water conservation and agreement of means remain necessary, while distribution agreement supplies the additional test this diagnosis shows is essential. The subsequent [reconstruction controls](collision_repair.md) retain a failed beta candidate: its mean and variance are accurate at a cell crossing, while its number CDF differs from the exact analytic distribution by 18.4 percentage points. The entropy alternative reduces that crossing discrepancy to 0.185 percentage points in one completed control. Its original time-stepping refinements reach their runtime limits, but the subsequent [finite-time capture method](collision_jump_results.md) completes six analytic and six snow-kernel bath controls, with distribution and accounting criteria satisfied. It remains a replenished-bath component. The separate [event study](event_column_results.md) now passes the matched closed-box comparisons and completes the coupled column; its failed volume refinement prevents adoption of a predictive distribution.

Growth histories need a separate treatment. The archived paths provide conditional temperature, humidity and liquid-water histories with explicit censoring; their cadence and coarse spatial resolution prevent treating them as resolved cloud microphysics. Adding aggregate geometry and contact-network failure to a numerically spreading population would leave two explanations entangled. The present controls make the numerical part testable before selecting those physical laws.

## Verification and retained limits

All 241 snow tests passed after the completed products were installed. The layer
check covered 1,204 tracked and untracked source files without violations, and
`research/check.py --require-inputs` passed its retained-byte and numerical
checks. The 50 snow producers and JSON products present at the start of this
extension remain byte-identical. Work is confined to `domain/atmosphere`; the
main checkout was clean at that verification.

The [check record](results/collision_exposure_checks.json) binds the new products,
compiled component reference, source reviews and figure checks. The subsequent [checkpoint check](results/diagnostic_checkpoint_checks.json) records a full
repository `make check` pass before committing this extension, with 3059 Python
tests passed and 20 skipped. The staged layer check covers all 1,204 source files.
The earlier check records retain their original timing and scope. Software and
artifact integrity do not clear the failed lunar domain controls, late collision
nonconvergence, meteorological sampling limits or the Earth aggregation
comparison discrepancy. A calibrated lunar snow population and its predictive
optics remain outstanding.
