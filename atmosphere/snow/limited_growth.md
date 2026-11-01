# Limited mass transport for vapour growth

The new limited transport substantially reduces the deposition spreading measured by the independent characteristic calculation. At 144 mass pivots, it reduces the excess relative mass variance over the once-projected reference by **93.1%** and reduces the mass above 5 mm from **14.47% to 0.001089%**. It preserves nonnegative populations and closes number and water ledgers. A remaining mass-resolution error is measurable, so this result supports testing the method in the coupled column rather than declaring that column converged.

[limited_growth.py](limited_growth.py) is a separate implementation. It preserves the checkpointed [sectional solver](sectional.py), the [first-order growth diagnostic](growth_transport.py), and their products. The [new product](results/limited_growth.json) binds its own source hashes and the frozen characteristic-comparison product. No physical rate coefficient, aggregate strength or size limit was fitted to obtain the improvement.

## Numerical formulation and water budget

The population remains a number concentration at each existing mass pivot. Arithmetic midpoints between adjacent pivots define the mass-cell faces, and boundary ghost masses continue the endpoint grid ratios. Cell-average number density is number concentration divided by cell width. A limited linear reconstruction uses the actual cell centroids, which differ from the pivot positions on a geometric grid. Its outgoing face value stays between the neighboring cell-average densities. Only reconstructed face values are limited; population entries are never clipped or renormalized.

For positive vapour growth, the outgoing number flux is the reconstructed density times the existing growth law evaluated at the face mass. If the uncapped number fluxes are `F_i`, their discrete deposition demand is

\[
P_{\mathrm{grid}}=\sum_i (m_{i+1}-m_i)F_i,
\]

including the final gap to the absorbing ghost. A common factor `min(1, P/P_grid)` caps all fluxes to the available vapour. Neighboring flux differences advance number concentration, while the same sum records deposited mass. A last-face transfer is recorded at the ghost mass in the upper-domain ledger and is never added to snowfall.

This explicitly changes the finite-grid growth quadrature from the old `sum n_i G(m_i)`. It is a numerical approximation to the same continuous growth equation. On refinement, the reconstruction and the discrete first-mass-moment quadrature must approach that equation; they are not assumed identical on a coarse grid. The unchanged positive SSPRK2 stepper uses an outgoing-flux CFL bound and rejects negative trial stages by reducing the step.

The reusable interface is `prepare(masses)`, followed by `LimitedGrowth(masses, face_growth, cfl=...)`. The preparation returns face masses for evaluating the physical growth law. `evaluate(number, production)` returns outgoing number transfers, the actual deposited-mass rate and its timestep limit. The object also exposes its upper ghost mass and pivot gaps for the surrounding column's ledgers.

## Characteristic comparison

The experiment uses the same fixed −10°C air, initial 0.1–0.4 mm population, vapour supply and closed-box duration as [growth_transport.md](growth_transport.md). Its 132.589-hour interval is a transport stress test based on the earlier simulation duration, not an inferred flake residence time. The continuously evolved physical population remains within 2.986–3.964 mm. Each grid is also compared with its own projected seeds transported continuously and projected once at the endpoint.

| Pivots | Once-projected reference variance | First-order variance | Limited variance | Limited mass above 5 mm | Limited mass 99th percentile |
|---:|---:|---:|---:|---:|---:|
| 96 | 0.03797 | 0.68189 | 0.09818 | 1.106% | 5.019 mm |
| 144 | 0.03083 | 0.41546 | 0.05751 | 0.001089% | 4.491 mm |
| 216 | 0.02769 | 0.26702 | 0.03830 | 0.000002229% | 4.047 mm |

Variance here means `Var(m)/mean(m)²`. The excess above the once-projected reference decreases by 90.6%, 93.1% and 95.6% relative to the first-order method. It remains positive in all three cases; the method did not produce an artificially narrower distribution than either reference. Errors against the continuous characteristics also decrease with mass-grid refinement. Quantiles remain discrete pivot quantiles.

At shorter times, the change is also substantial. With 144 pivots after 6 hours, relative variance is 0.13972, compared with 0.12715 for the continuous grid-matched reference and 0.39164 for first-order transport. Its mass 99th percentile is 1.034 mm, against the physical support maximum of 0.972 mm. At 24 hours the limited mass 99th percentile is 1.794 mm, while the physical support ends at 1.761 mm; the first-order value was 3.111 mm.

Reducing the 144-pivot maximum step from 30 to 15 to 7.5 seconds changes final relative variance by 3.39 × 10⁻⁷ and 1.13 × 10⁻⁷. The final change is about 0.00042% of the remaining difference from the once-projected reference. Timestep error is therefore small compared with the residual spatial error in this experiment. Both temporal controls and the predefined improvement gates passed. These gates require a measured reduction in error; they do not establish zero error or full-column convergence.

## Verification and limits

Smooth constant-growth translation of a compact distribution tests the reconstruction independently of the physical growth law. From 128 to 256 cells, its number-distribution L1 error converges at orders 1.94 on a uniform grid and 1.89 on a geometric grid. At 256 cells, errors in the known deposited mass are 8.56 × 10⁻⁷ and 6.37 × 10⁻⁶ respectively. Limiting reduces local order near sharp features, while these smooth controls approach second-order behavior.

Additional tests check linear-density reconstruction on geometric cells, nonnegative fluxes, empty donors, a sharp distribution with a tiny positive tail, dry conditions, the vapour cap, ghost exits, input validation and the saved product's identity. The mass and number identities were independently reviewed on 60 sparse irregular-grid cases: scaled residuals were at most 1.85 × 10⁻¹⁵ and 1.81 × 10⁻¹⁶. No negative transfer, negative allowed Euler stage or empty-donor flux was found in those checks.

The seven limited-transport tests passed, and all fourteen tests passed when run together with the characteristic comparison. The five physical-box runs took 17.93 seconds on one worker. Their largest relative water-ledger residual was 8.14 × 10⁻¹³. No committed baseline producer or result was changed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -m atmosphere.snow.limited_growth --budget-s 180
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -m pytest atmosphere/tests/test_snow_growth_transport.py atmosphere/tests/test_snow_limited_growth.py -q -p no:cacheprovider
```

This result repairs a demonstrated deposition-transport deficiency to a useful degree and supplies a tested alternative for the column. Collision discretization, the fixed aggregate geometry, whole-aggregate failure, hydrodynamic interception and realistic moist trajectories retain their separate uncertainties. Their effects cannot be inferred from a condensation-only box, and a reduced coupled tail would still require independent numerical and physical assessment.

**Provenance, 25 October 2026.** [results/limited_growth.json](results/limited_growth.json) was built with the rain
solver as it stood before 25 October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998)
measured value. It keeps its values; its record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for
byte as [atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and its identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
