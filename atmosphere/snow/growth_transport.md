# Mass-grid spreading during vapour growth

The existing first-order mass transport produces substantial numerical spreading in a controlled growth-only box. Continuous-mass trajectories finish within **2.986–3.964 mm**, whereas the 144-pivot calculation puts **14.47% of its mass above 5 mm** and a mass 99th percentile at **7.79 mm**. The difference persists when its timestep is reduced. This establishes a deposition-transport error independently of collisions; it does not establish how much of the earlier coupled column's runaway that error caused.

[The producer](growth_transport.py) and [product](results/growth_transport.json) leave the checkpointed sectional solver, comparison runner and their results unchanged. A continuous reference evolves particle masses with fixed number weights, while the comparison evaluates the same growth terms used by the sectional column. Neither includes collisions, settling, new seeds after initialization or a physical upper-size cutoff.

## Controlled problem and independent reference

The box holds the lunar air at −10°C, about 40.10 km and 525 hPa on the existing profile. Its initial population is 100,000 particles per cubic metre, log-uniform in diameter from 0.1 to 0.4 mm. The prescribed vapour source is 5.1186 × 10⁻⁸ kg m⁻³ s⁻¹, obtained from the heavy-column 3 mm/h supply divided by that column's 16.28 km depth. Geometry, fall speed, ventilation and diffusion are the existing laws. Initial diffusional demand is about 6.20 times the supplied vapour, so shared supply limits the growth.

The characteristic equations are

\[
\frac{dm_a}{dt}=G(m_a)\min\!\left(1,\frac{P}{\sum_b n_bG(m_b)}\right),
\]

where each number weight `n_a` stays fixed. An adaptive DOP853 integration follows the masses and the deposited-water ledger. Two zero-weight trajectories follow the endpoints of the physical input distribution, giving its support independently of quadrature-node spacing. There is no mass-grid projection during that evolution.

The longest diagnostic interval is the earlier calculation's complete 132.589-hour **simulation duration**, not an inferred flake residence time. This closed box accumulates 24.58 g/m³ of ice by then. It deliberately isolates prolonged numerical transport and is not a realistic cloud history. Intermediate outputs at 6, 24 and 72 hours show how the error develops.

Each mass grid has a second characteristic reference starting from exactly its own projected seed population. Those continuously grown masses are projected only once at each reported time. Comparing the fixed-grid population with that result separates repeated transport spreading from the unavoidable initial and endpoint interpolation between pivots.

## Measured spreading

At the final time, the physical characteristic distribution has relative mass variance `Var(m)/mean(m)² = 0.02525`. The following comparisons retain the grid-matched initial population and the single endpoint projection.

| Mass pivots | Relative variance after one endpoint projection | Relative variance after repeated grid transport | Mass above 5 mm | Mass 99th percentile |
|---:|---:|---:|---:|---:|
| 96 | 0.03797 | 0.68189 | 30.44% | 10.01 mm |
| 144 | 0.03083 | 0.41546 | 14.47% | 7.79 mm |
| 216 | 0.02769 | 0.26702 | 10.05% | 6.59 mm |

The declining variance and fixed-threshold tail demonstrate substantial mass-resolution sensitivity. Even the 216-pivot result remains much broader than its independently transported reference. The tabulated quantiles are discrete pivot quantiles; they do not warrant interpolation to more precise physical dimensions.

The error is already visible at shorter times. After 6 hours, the physical characteristic support is 0.519–0.972 mm. At 144 pivots, the reference after one endpoint projection has relative mass variance 0.13287, while repeated grid transport gives 0.39164 and a mass 99th percentile of 1.493 mm. After 24 hours, the physical support is 1.142–1.761 mm; the corresponding variances are 0.06597 and 0.38941, and the grid's mass 99th percentile is 3.111 mm.

The product also reports mass beyond each once-projected reference's highest occupied pivot. At the final time those fractions are 30.44%, 31.66% and 28.28% for 96, 144 and 216 pivots. Their thresholds differ—4.371, 4.097 and 4.047 mm—so their nonmonotonic sequence is not a convergence measure. The common 5 mm threshold and the variance comparisons provide the cleaner cross-grid comparisons.

## Numerical checks

Quadrature with 64 and 128 nodes and a tighter integration tolerance changed the final first two mass moments and maximum-support diameter by at most 2.31 × 10⁻¹⁰ relatively. Those moment/support controls pass the stated 10⁻⁶ comparison threshold. Discrete quadrature quantiles have their own sampling resolution; that threshold is not a claim about every reported quantile.

For 144 pivots, reducing the maximum step from 30 to 15 to 7.5 seconds changes the final relative mass variance by 1.01 × 10⁻⁶ and then 2.51 × 10⁻⁷. The final change is only 6.54 × 10⁻⁷ of the measured excess variance over the once-projected reference. The excess-support mass fraction has similarly small temporal changes. Both are well below the required 5% of the spatial discrepancy. The spreading is therefore resolved as a spatial-discretization effect in this box.

Analytical tests supply independent controls. For constant growth `G` on a uniform mass grid, first-order transport produces Poisson migration between pivots, with artificial mass variance `G Δm t`; the continuous monodisperse population has none. For `G(m)=a m` on a geometric mass grid of ratio `r`, the first moment is exact in the continuous-time discretization, but its relative variance is `exp[a(r−1)t]−1`. These tests show why decreasing the timestep cannot remove the mass-grid error. Further checks cover the supply cap, number conservation, endpoint tracers, one-time projection, physical seed quadrature and saved-product identity.

All seven focused tests passed. The complete numerical comparison took 7.90 seconds on one worker, within its 180-second limit. The largest relative water-ledger residual across the reported fixed-grid snapshots is below 10⁻¹¹; upper-domain mass loss is negligible in this box. Full producer and input hashes, including the prior comparison product used for the duration, are recorded in the new product. No empirical validation or new full-column convergence is claimed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -m atmosphere.snow.growth_transport --budget-s 180
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -m pytest atmosphere/tests/test_snow_growth_transport.py -q -p no:cacheprovider
```

## Consequence for the next repair

Mass conservation alone permits a badly broadened growth distribution. Repeatedly moving fixed pivots a short distance and immediately projecting them back would recover the same first-order spreading as the timestep tends to zero. A useful transport repair must reduce that spatial error and pass this characteristic comparison, alongside positivity and the vapour ledger, before its full-column populations are interpreted.

The retained collision law remains a separate issue. Fixed mass–size and area assumptions give superlinear large–small collection in their unbounded asymptote, and the installed chipping law does not erode a larger collector during encounters with smaller particles. This box measures neither that amplification nor its competition with settling and fragmentation. It establishes a significant numerical growth-transport deficiency that should be controlled before attributing the coupled runaway to physical closure alone.
