# Full prescribed lunar snow column: event solution and refinement

The full cloud-top-to-freezing column now runs with donor-depleting collision events, continuous particle masses, vapour-limited growth, chipping and settling. All **28 full-column controls** completed, including four independent realizations at each setting. The separate closed-box comparison passes its declared sampling checks. **The lunar size distribution fails volume convergence**, so this calculation supplies a complete numerical diagnostic rather than an adopted prediction of snowflake sizes or storm optics. The [column product](results/event_column.json), [closed-box product](results/event_closed_boxes.json) and [comparison figure](../../visualization/snow/generated/event_column_comparison.png) retain that distinction.

## Donor depletion passes the independent collision comparison

The new [event sampler](event_column.cpp) uses positive rate sums and rejected candidate encounters to sample the same finite-volume generator as the independently written direct-event reference. Every sticking event removes both parents, creates their actual summed mass and refreshes the rates. The [method note](event_column_methods.md) derives the sampler and distinguishes collision accuracy from the full column's volume, height, splitting and fragment-allocation approximations. The tested finite-time bath operators remain separate components; the full column extends the direct event reference instead of treating its donor population as replenished.

There are 128 realizations at each of the reference's four volumes, for 512 runs in total, with observations at 0, 2, 4, 6 and 12 hours. Every stored number and second-moment mean and mass-CDF query meets the five-combined-standard-error compatibility screen plus its numerical tolerance. At nonzero times, the largest standardized moment difference is 2.101; the maximum relative mass residual is 1.98e-15. At six hours and 2,560 initial particles, the second mass moment is (2.589 ± 0.112) × 10⁻¹² kg²/m³, compared with (2.939 ± 0.163) × 10⁻¹² in the direct reference. Uncertainties here are one standard error across realizations. At initialization, standardized ratios of tiny rounding differences to nearly zero sample errors are uninformative; the declared numerical tolerance handles those deterministic entries.

This establishes agreement with the reference at matched finite volume. The twelve-hour reference itself remains unconverged in volume. Its successful reproduction cannot certify an infinite-volume large-particle tail.

## Complete lunar runs retain the prescribed forcing

The lunar cloud extends 16.281 km from −15 °C to freezing. It receives the same expected source flux corresponding to 100 seed crystals/L, the same corrected Brown–Francis mass law and fixed representative dendrite area, and a maximum vapour supply of 3 mm/h. Its integration lasts **132.589 hours**, with the last 66.295 hours supplying the arrival statistics. That is the duration for which the cloud is maintained; individual particles have their own transit and growth histories. The upwind height-transport generator matches the earlier sectional column. The primary settings use 30 layers and a 60 s maximum split step.

Each computational particle represents one particle in its statistical volume. Enlarging that volume tests the finite-population approximation while retaining the imposed concentrations. The following values are diagnostic means of realization quantiles, with one standard error across four realizations; they are not physical forecasts or quantiles of a pooled sample.

| Volume per layer (m³) | Mass median (cm) | Mass 99th percentile (cm) | Arriving water (mm/h) | Mean particles in lowest layer |
|---:|---:|---:|---:|---:|
| 0.05427 | 12.27 ± 0.18 | 19.83 ± 1.23 | 2.27 ± 0.02 | 0.414 |
| 0.21707 | 18.53 ± 0.42 | 36.28 ± 4.23 | 2.42 ± 0.02 | 0.829 |
| 0.86829 | 23.46 ± 0.50 | 74.09 ± 5.18 | 2.54 ± 0.01 | 2.446 |

The largest setting therefore still has few particles near the freezing boundary despite many seed particles aloft. Finite-volume fluctuations change the population competing for local vapour as well as the sampling of its tail. The table's systematic shift exceeds the variation of the median among these realizations. It demonstrates unresolved volume dependence; it does not isolate a unique physical cause or prove gelation of the settling, fragmenting column. Large sizes arise with continuous masses under the inherited laws, so the old pivot-allocation error is not required for their appearance in this finite-volume model.

## The large particles respect the prescribed water supply

Across all Earth and lunar runs, the largest absolute relative water residual is **6.28e-16** and the number residual is **1.05e-16**. Whole-run water accounting, in kg/m², is:

| Lunar layer volume (m³) | Injected seed ice | Deposited vapour | Left through freezing level | Ice remaining |
|---:|---:|---:|---:|---:|
| 0.05427 | 7.281 | 287.293 | 283.938 | 10.636 |
| 0.21707 | 7.280 | 308.756 | 304.222 | 11.814 |
| 0.86829 | 7.288 | 325.294 | 319.704 | 12.879 |

A further retained fractional-fragment remainder closes each budget. At the largest lunar setting its total mass is bounded by 5.7e-06 kg/m². Its omitted instantaneous growth demand accumulated over the run is bounded by 0.00886 kg/m²; this bound excludes nonlinear effects of delayed participation. The source permits at most **397.768 kg/m² of vapour deposition** over the full run, plus seed ice. The reported deposition stays below that cap. Unused supply is not silently reassigned to snow.

The inherited [mass–size law](habits.py) assigns mass in kg as `0.0121 × D^1.9`, with maximum dimension D in metres. Under that extrapolation a 60 cm aggregate contains only 4.58 g of ice. Large dimensions can therefore coexist with a finite water budget; the unresolved question includes whether such sparse structures retain the assumed geometry and survive. This mass is a calculation from the imposed law, not a measurement of a lunar flake.

This answers the bookkeeping question within the prescribed scenario: the particles do not consume more water than the forcing supplies. The forcing continuously replenishes the cloud, so its accounting does not establish that a real lunar storm could sustain the assumed moist depth and supply. There is no numerical upper mass boundary or imposed maximum diameter. The largest possible particle in a finite realization is nevertheless bounded by that realization's total supplied water, a bound that moves when volume changes.

## Independent refinements and Earth transfer

The 30 s control halves only the macrostep. The 60-layer control doubles the represented horizontal area, holding each layer's statistical volume fixed while halving vertical spacing. Earth uses the same source concentration and geometry, with 1 mm/h vapour supply and its shorter cloud. The table reports differences in group means; its 5% quantile/flux and 0.05 CDF thresholds are practical numerical screens, not confidence bounds.

| Comparison | Maximum mean mass-CDF change (percentage points) | Flux change | Median change | 99th-percentile change | Screen |
|---|---:|---:|---:|---:|---|
| Lunar volume: 0.0543 → 0.2171 m³ | 46.13 | +6.68% | +51.01% | +82.97% | Fail |
| Lunar volume: 0.2171 → 0.8683 m³ | 28.99 | +5.08% | +26.59% | +104.24% | Fail |
| Maximum step: 60 → 30 s | 6.45 | +0.27% | +0.12% | +5.63% | Fail |
| Height: 30 → 60 layers, fixed layer volume | 3.03 | +2.02% | -1.38% | -1.25% | Pass |
| Earth volume: 0.0308 → 0.1231 m³ | 0.21 | +0.01% | +0.10% | -0.79% | Pass |

The Earth mass medians are 1.101 ± 0.004 and 1.102 ± 0.002 mm; their 99th percentiles are 4.703 ± 0.100 and 4.666 ± 0.066 mm. These are internal transfer controls. They do not add a new reproduction of the published Earth column simulations: the earlier [published-simulation comparisons](welss_column.md) and [author-model aggregation rerun](welss_aggregation.md) retain their separate evidence and limitations. Formula parity and agreement with a generic Earth scenario cannot validate the extrapolated lunar aggregate law.

Four realizations give limited precision for upper quantiles. A failed small refinement screen alone does not identify timestep or height bias independently of sampling variation. The strong lunar volume dependence already prevents admission of a size forecast. The full products retain standard errors, all four gathering blocks and complete area-weighted arrival moments so these issues can be examined without deleting the large particles.

## Interpretation and verification

The repair now reaches the complete prescribed column, with paired donor consumption and explicit water accounting. Its numerical controls leave the stable lunar population unresolved. Evolving aggregate geometry, interception, contact failure and climate-supported moist histories remain separate physical requirements. In particular the chipping law still acts primarily on the smaller collision partner and supplies no whole-aggregate failure mechanism. The calculated 15–60 cm range therefore remains a question about both the kinetic limit of the model and the applicability of its physics. No new snowstorm optical prediction has been adopted from this population.

The implementation and products are in `domain/atmosphere`. C++17 compilation uses `g++`, without fast-math flags, and caches the executable library by source, compiler and flags in the temporary directory. The 28 full runs used at most two numerical workers, with Python numerical libraries restricted to one thread. Their summed per-run elapsed times are 66.3 minutes; this is not the matrix's elapsed wall time or a controlled hardware benchmark. All earlier scientific producers and saved results remain intact.

Verification completed with 29 focused numerical and figure-reader tests passing; the seven figure-reader tests passed again after the final layout change. `make check` passed with 3,161 Python tests, 20 skips and 191 JavaScript tests, plus the repository provenance and ensemble integrity checks. A separate layer check included tracked and untracked sources: 1,226 files, zero violations. `research/check.py --require-inputs` passed, and 114 earlier scientific files retained their recorded hashes. Two builds produced byte-identical PNG, SVG and manifest files. The final 2448 × 1908 PNG was inspected at its original resolution for labels, units, uncertainty displays, numerical-failure notice and clipping. These checks establish software and artifact consistency; the failed scientific refinement screens remain visible.

The [verification record](results/event_column_checks.json) records the actual numerical, software, provenance, source-retention and figure checks. The [source ledger](event_column_sources.json) records the targeted method reading; no full scholarly review or empirical validation is implied. Reproduction commands are:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.event_column_benchmarks --boxes
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.event_column_benchmarks --run-all --workers 2
python -m visualization.snow.event_column
python -m pytest -q atmosphere/tests/test_snow_event_column.py visualization/snow/tests/test_event_column.py
```
