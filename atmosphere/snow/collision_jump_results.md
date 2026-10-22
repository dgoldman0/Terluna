# Finite-time capture: completed distribution controls

Finite-time capture now completes all six analytic cell-crossing controls and all six size-dependent snow-kernel controls. Each passes its declared checks of distribution shape, mean, variance, number and water accounting, and truncation bounds. This establishes a working maintained-bath component of the collision repair. The full lunar column still uses its preceding solver; closed aggregation, which consumes the donor population, has not yet been repaired or tested with this component.

## What changed

The preceding work identified errors in mass representation and a separate time-integration constraint, developed in the [reconstruction study](collision_repair.md). Fixed-pivot allocation gives small captures unrealistically large mass jumps, and the first bounded beta reconstruction retains accurate moments while distorting the distribution at a cell face. The entropy reconstruction improves that shape, but an explicit Euler stage can admit at most one expected capture for a point population. That restriction applies even when the point population has a negligible water mass, so small tails can force many short steps.

The new [constant-rate operator](collision_jumps.py) advances the actual Poisson distribution of capture counts before projecting its number, first mass moment and second mass moment into cells. The [size-dependent operator](collision_birth.py) evaluates rates along each actual-mass birth ladder and uses a positive uniformized transition. It therefore allows the capture rate to change as the collector grows within a macrostep. Both retain the bounded entropy reconstruction, add donor mass explicitly and stop particles at their first physical-domain exit. The [method note](collision_time_integration.md) derives the updates, truncation bounds and their relation to the literature. Earlier producers and saved products remain intact.

## Analytic crossings

The [analytic product](results/collision_jump_boundary.json) contains six completed controls. The baseline uses 191 cells, order-16 source quadrature, a 60 s maximum macrostep and a constant capture rate of 100/s. Separate controls halve the macrostep, double quadrature order, use 96 or 381 cells, or shift every grid pivot by a factor of 1.02. Each samples before, at and after its own first upper cell face, followed by a later observation. This checks crossing phase on each grid instead of keeping all grids at the 191-cell reference times. The baseline still matches the earlier 191-cell crossing at 207.396 s.

Across all saved comparisons the largest number-CDF error is 0.002584, or 0.258 percentage points, below the preset 0.005 threshold. Maximum absolute mean error is 1.6e-15 times the initial mass; maximum relative variance error is 2.5e-11. The analytic CDF is evaluated on both sides of each Poisson jump within twelve count standard deviations, with unqueried reference tails reported separately. Independent tests recompute the Poisson CDF through the incomplete-gamma identity. A continuous reconstructed density cannot reproduce arbitrary discrete low-count distributions; these high-count controls establish accuracy for their stated regime.

The matched baseline takes 8 macrosteps and 1.26 s, compared with 14,163 accepted SSPRK steps and 459.4 s in the earlier completed run. These are two recorded executions with different quadrature rules and system workloads, rather than a controlled hardware benchmark. One macrostep versus two half steps at the face changes the number CDF by 0.000214 and relative variance by 2.5e-12. The difference includes the intermediate reconstruction; the underlying exact constant-rate transition itself has the semigroup property.

## Size-dependent snow kernel

The [physical-bath product](results/collision_jump_bath.json) retains the frozen kernel, initial 15.8608 cm collector, 0.209337 mm donor, bath of 100 particles/L, and observations at 81.788 and 272.627 s. No deposition, settling or breakup acts in this control. The table lists the first and final nonzero observations in each paired entry.

| Mass cells | Quadrature order | Maximum step (s) | Ladder events | Number-CDF difference, first / final | Variance / reference sample variance, first / final | Macrosteps | Wall time (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 191 | 16 | 0.5 | 256 | 0.01175 / 0.00777 | 0.97294 / 1.00228 | 547 | 156.6 |
| 191 | 16 | 0.25 | 256 | 0.01175 / 0.00777 | 0.97294 / 1.00228 | 1096 | 192.8 |
| 191 | 32 | 0.5 | 256 | 0.01175 / 0.00777 | 0.97294 / 1.00228 | 695 | 452.8 |
| 96 | 16 | 0.5 | 256 | 0.01175 / 0.00777 | 0.97294 / 1.00228 | 547 | 18.8 |
| 381 | 16 | 0.5 | 256 | 0.01175 / 0.00778 | 0.97294 / 1.00228 | 547 | 201.3 |
| 191 | 16 | 0.5 | 384 | 0.01175 / 0.00777 | 0.97294 / 1.00228 | 547 | 240.4 |

The [initial matrix](results/collision_jump_bath_attempt.json) completed five controls and stopped the order-32 control at its 300 s wall-time limit. The [continuation runner](collision_birth_benchmarks.py) preserves those five completed runs unchanged and repeats only the unfinished control with a 600 s allowance. Its execution history binds the original attempt and keeps the numerical settings and acceptance gates unchanged. Initial controls used two single-thread workers; the continuation used one.

Every case meets the same mean and variance gates used for the earlier candidate: mean within five reference standard errors and variance within 5% of the independent sample variance. Every number-CDF discrepancy is below the reference's 0.015005 DKW half-width. That band has 95% coverage uniformly in mass at each saved time separately. It supplies compatibility with the 8,192-history sample, rather than a joint guarantee across cases or a proof of numerical convergence. The table reports variance relative to the reference's unbiased sample variance, while each candidate represents a population variance.

Number residuals are at most 2.4e-13 per initial collector; water residuals are at most 9.1e-14 of initial collector mass after including the prescribed donor input. These balances remain separate from the distribution checks. The bath is replenished: its added water is recorded as external supply. It does not impose a finite storm moisture inventory or demonstrate the formation of its initially large collector.

## Truncation, verification and remaining work

The saved controls separate omitted Poisson terms from reaching an artificial birth-ladder limit. The largest conservative centered-second-moment bound, divided by the reference variance, is 5.5e-24 in the analytic matrix and 4.6e-09 in the physical matrix, against a 10⁻⁴ admission threshold. These are transition-truncation bounds conditional on reconstructed sources. Entropy reconstruction, quadrature and floating-point errors are checked separately; neither the bounds nor the water ledger constitute a total solver-error certificate.

The [comparison figure](../../visualization/snow/generated/capture_comparison.png) shows analytic crossing errors, physical CDFs and sampling bands, each physical comparison error, and recorded runtimes. Its [renderer](../../visualization/snow/capture.py) verifies producer/input identities, native moments, stable variance and ledgers before plotting stored results. The [verification record](results/collision_jump_checks.json) records the actual software, provenance and artifact checks, including population scaling to 10⁻³⁰⁰, independent matrix-exponential and negative-binomial controls, and unchanged earlier scientific products.

All 321 snow tests pass, and the layer check covers 1,221 tracked and untracked source files without violations. The retained-byte and numerical provenance check passes with required inputs present. The full repository `make check` belongs to the preceding committed checkpoint and was not rerun for this uncommitted continuation.

A subsequent [full-column study](event_column_results.md) now completes donor-depleting closed aggregation and the coupled prescribed column through a separate continuous-mass event solver. It passes the matched finite-volume reference checks, while its lunar column remains sensitive to volume. This report retains the earlier maintained-bath component results and their own verification scope. The bath results do not resolve aggregate geometry, whole-aggregate fracture or moist growth histories, and supply no revised 15–60 cm size forecast or predictive snowstorm optics.
