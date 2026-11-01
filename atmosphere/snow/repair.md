# Repairing the snow population calculation

This report records the numerical and literature work preserved in checkpoint `16cde54`. The subsequent [growth-transport diagnosis](growth_transport.md), [limited transport repair](limited_growth.md), [Böhm comparison](boehm.md) and [Earth benchmark specification](earth_benchmarks.md) extend those tests without changing the checkpoint producers or products. The latest [collision and moist-exposure synthesis](collision_and_exposure.md) records the subsequent findings and their limits.

The targeted literature pass and numerical repair address different parts of the failed snow calculation. The [geometry review](geometry_literature.md) and [strength review](strength_literature.md) identify physical mechanisms and their supporting measurements. The new [sectional solver](sectional.py) removes the inherited particle-weight inconsistencies and batched collision decisions while retaining the prescribed growth environment and physical rate laws. Its numerical distributions remain conditional on those laws. The [earlier audit](numerical_audit.md), stochastic populations and optical products remain reproducible records of the failure.

## Completed column controls

All three completed heavy-snow controls fail the finite mass-domain check. They retain the prescribed liquid-saturated lunar column, a −15°C cloud top, 100 initial crystals per litre, nominal vapour supply of 3 mm/h, 30 height layers and the same physical laws. Each runs for three initial-crystal crossing times (about 133 hours), with its second half used for the distribution averages. The percentages below account for all injected and actually deposited mass over the complete run; the resolved diameter boundary is numerical, not a proposed physical limit.

| Mass pivots | Largest resolved diameter (m) | Escaped size domain | Reached freezing level | Remaining in column |
|---:|---:|---:|---:|---:|
| 96 | 10.00 | 59.69% | 37.10% | 3.21% |
| 144 | 10.00 | 57.77% | 38.98% | 3.25% |
| 170 | 108.69 | 80.26% | 16.89% | 2.85% |

Extending the domain raises the escaped share from 57.77% to 80.26%. It also reduces the diagnostic arriving mass flux by 64.0% and the freezing-level projected-area concentration by 47.1%. The resulting optical estimate therefore depends strongly on where the numerical domain ends.

The relative mass-balance residual is at most 9.07e-13. Good mass accounting therefore coexists with a badly inadequate finite size domain. Extending the 144-pivot grid by 26 pivots retains its spacing and tests a larger domain, while the 96-to-144 comparison tests resolution sensitivity. Neither comparison supplies a converged physical population. Numerical spreading of a steep size distribution can promote large-particle growth, and the retained fixed geometry and chipping laws can amplify it; these controls do not uniquely divide the error between those mechanisms.

The mass-domain tolerance is 10⁻⁶. Every completed result exceeds it substantially, so the exported surviving-outflow quantiles and extinction are diagnostic outputs, not usable lunar snowfall or visibility predictions. The [product](results/repair.json) includes full mass ledgers and explicitly lists six unperformed controls: 216 mass pivots, a 15 s maximum step, half the adaptive CFL, 60 height layers, moderate lunar supply and an Earth reference. Those broader comparisons are deferred at this failed domain check; full mass-grid, time, spatial and domain convergence has not been established.

## What the literature changes

Approximately quadratic mass–size scaling has a physical basis. Geometric aggregation simulations and radar observations support a fractal dimension near two; the accompanying theory explains how unequal-size encounters can produce denser structures and thereby oppose preferential growth of the largest aggregates. This supports investigating large porous snow with a structured model. It supplies neither a fixed dimensional coefficient nor an assurance that the contacts survive. Direct observations retrieved in this pass establish aggregates around 4–6 cm, with measurement limits that prevent treating those sizes as a ceiling. The detailed sources and reading extents are in [geometry_sources.json](geometry_sources.json).

The fixed local combination `m ∝ D^1.9` and constant projected-area ratio omits that evolution. [kernel_diagnostics.py](kernel_diagnostics.py) calculates its encounter kernel against specified small partners. At −10°C against a fixed 0.2 mm particle, its local logarithmic slope with collector mass is about 1.06 at 10 cm and 1.04 around 60 cm. Values above one describe preferential collection in this particular calculation. The ideal aggregation theory identifies a related runaway mechanism, but these local slopes do not prove gelation in a finite column with settling, deposition and fragmentation.

An exponent-two comparison is matched to the original mass at 1 cm and labelled as a diagnostic. It changes the small partner's mass and speed as well as the collector's geometry. Its finite-size collision rates need not decrease, and its local slopes still exceed one over the displayed range. This comparison therefore does not establish a physical repair. The [kernel product](results/kernel_diagnostics.json) preserves the inputs and the derivative step-halving check.

The contact literature also supplies useful mechanisms. Demmenie et al. measured neck growth between approximately millimetre ice spheres touching without external compression at −3°C. Humidity affected that evolution. Other experiments measured pull-off forces after specified preloads, while direct snow-collision experiments measured habit-dependent fragment production. Their conditions differ from a freely suspended dendritic network. Indirect ice-aggregate simulations further show that strong sintered contacts can suppress energy-dissipating rolling, increasing rebound or brittle damage. The [strength review](strength_literature.md) consequently distinguishes contact age, neck dimensions, connectivity, deformation and largest-remnant mass. No universal decimetre failure diameter is inferred or installed.

The inherited chipping law also has a quantitative limitation: it creates at most 100 fragments of fixed 0.1 mm size per nonsticking encounter, a maximum chipped mass of 3.04 × 10⁻⁸ kg under the current mass law. Even when a 15 cm aggregate is the smaller collision partner, this removes at most 0.0092% of its mass; at 60 cm, at most 0.00066%. A larger collector meeting smaller particles loses none of its own mass in this treatment. These calculations describe the installed chipping rule, rather than measured strength. The kernel product records them separately from the geometry comparison.

## Independent numerical formulation

The corrected reference carries number concentration on a fixed grid of particle masses in each height layer. A collision contributes separate mean rates for sticking and chipping. A product between two mass pivots is represented by fractions at the adjacent pivots chosen to preserve its number and mass. This is numerical interpolation, not physical breakup. There are no packet counts to round, no absolute one-fragment admission threshold and no shared sticking draw for a batch of encounters.

The formulation follows the conservation principle of [Kumar and Ramkrishna's fixed-pivot method](https://doi.org/10.1016/0009-2509(96)88489-2). Their [moving-pivot paper](https://doi.org/10.1016/0009-2509(95)00355-X) explains a remaining limitation: coarse fixed grids can overpredict steep distribution tails. The current implementation uses explicit mass-grid comparisons to expose that error. Its analytical checks do not validate a coarse snow spectrum automatically.

The solver retains the corrected maximum-dimension geometry, representative dendrite area, lunar fall-speed and diffusion laws, temperature-dependent sticking, smaller-particle chipping and imposed water supply. All particles occupy their grid mass and therefore use the corresponding geometry and velocity throughout each evaluated stage. Positive adaptive second-order time integration advances collisions, conservative mass growth and conservative vertical transport together. A rejected stage reduces the timestep without changing the particle population.

The source number flux uses the velocity at the actual cloud-top state. The historical code used the uppermost cell midpoint, making the injected flux and initial-crystal crossing time depend slightly on vertical resolution. Sampling the boundary directly removes that dependence from the new refinement comparisons. It is an explicit numerical boundary correction alongside the deterministic treatment of the log-uniform injected sizes, whose total number and expected mass are preserved across mass grids.

## Finite domains and refinement

The last resolved diameter is a computational boundary, not a maximum allowed physical flake. An absorbing mass pivot immediately above the grid receives the fraction of a near-boundary product that would lie there; products beyond that pivot exit at their actual mass. Number and mass in those exits are recorded separately and never counted as snowfall. This conservative boundary treatment also avoids imposing an enormous encounter-frequency timestep merely because a practically empty top bin collides with tiny seeds.

The mass-domain check compares total escaped mass with all injected and deposited water, using a stated tolerance of 10⁻⁶. Passing it establishes only that mass-domain loss is small. The larger-domain comparison must also test area and extinction, and independent mass-grid, timestep and height-grid controls remain necessary. The upper-domain extension retains the same geometric spacing and adds larger pivots. The two construction formulas reproduce the shared pivots to floating-point precision, with relative differences up to 1.2 × 10⁻¹⁴.

The exported quantiles are discrete mass-pivot quantiles. Their changes include the grid's discrete steps. The raw maximum cumulative-distribution difference additionally treats roundoff-separated pivots as distinct points, which can create an artificial jump of a whole bin's mass in the domain-extension comparison. That metric is retained for transparency but excluded from the interpretation here; comparisons instead use the mass budget, mass flux, mass-weighted log diameter and optical area. Extinction is calculated from each pivot's arriving flux divided by its own freezing-level speed, with `Qext = 2`. The resulting 5% beam distance describes a homogeneous snow-only medium. It is not human visibility or a newly rendered storm.

## Verification and reproduction

The focused solver tests compare full constant- and additive-kernel coagulation distributions with analytical solutions, check second-order time refinement, and show reduction of known mass-grid error. Other controls check fractional-fragment retention, positivity, per-reaction number and mass accounting, domain exits, area-independent concentrations, grid-independent seed number and mass, vertical-source invariance and exact gathering intervals. Independent review also checked reaction conservation on irregular grids spanning twenty mass decades.

The completed check run passed all 80 focused snow tests. After assembling all three controls, the two repair-product tests passed again, including independent reconstruction of the mass budget and adequacy flag. The import-layer check covered 1,162 tracked and untracked source files with zero violations; the research integrity check passed for 110 retained files and its five protection-input hashes. The kernel diagnostic reproduced exactly after JSON normalization, and its derivative step-halving difference was 6.47 × 10⁻¹¹. [repair_checks.json](results/repair_checks.json) records these checks and the unperformed work. Before checkpoint `16cde54`, the full repository `make check` subsequently passed: 2,898 Python tests passed with 20 skipped, all 191 JavaScript tests passed, and the layer, research-integrity and ensemble checks passed. [checkpoint_checks.json](results/checkpoint_checks.json) records that later check. No new storm rendering was made.

The [comparison runner](repair.py) binds each case to its producer files, actual column input bytes and named shared constants. It rechecks that identity after a run and refuses stale caches during assembly. Results are cached in ignored `research/runs/snow_sectional` and exported to [results/repair.json](results/repair.json). The original population and audit products are not overwritten.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.repair --case moon_heavy_m144
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.repair --assemble
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.kernel_diagnostics
PYTHONDONTWRITEBYTECODE=1 python -m pytest atmosphere/tests/test_snow*.py -q -p no:cacheprovider
```

## Proposed physical comparisons

The literature supports several distinct experiments after the numerical errors are quantified. Each would retain the preceding calculation as its control, so that changing a tail can be attributed to the mechanism being tested.

1. **Interception before sticking.** Evaluate a published hydrodynamic collision efficiency using the local air properties and fall speeds, keeping the existing temperature-dependent sticking law separate. Compare the kernel for comparable and very unequal particles with the terrestrial cases covered by the parameterization. Its extrapolation to giant porous aggregates would remain explicit.
2. **Geometry that responds to aggregation.** Compare the fixed mass–size and area laws with a joint geometry model retaining monomer size, habit and aggregation history. The matched exponent-two diagnostic can identify sensitivity, but a reduced tail alone would not calibrate either geometry. Similar-mass particles also need a distribution of fall speeds if the model is to represent their collisions.
3. **Contacts along a moisture history.** Follow growth and sublimation along finite moist trajectories, retaining contact age and neck geometry separately from total mass. The uncompressed-contact experiments supply a relevant mechanism; the loaded pull-off measurements constrain different contact histories. Neither supplies a universal tensile coefficient for this lunar column.
4. **Whole-aggregate failure alongside chipping.** Resolve or parameterize the contact network's response to impact and unequal aerodynamic loading, reporting largest-remnant mass as well as fragment number. Where strengths and connectivity remain unmeasured, use declared sensitivity cases. A fitted maximum diameter would conceal that uncertainty.

These are proposed physical extensions, not changes already made to the snow equations. The completed implementation addresses numerical particle accounting and provides an independent reference; calibrated lunar snow sizes and storm optics remain unresolved.

**Provenance, 25 October 2026.** [results/repair.json](results/repair.json) and
[results/kernel_diagnostics.json](results/kernel_diagnostics.json) were built with the rain solver as it stood before 25
October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998) measured value. They keep their
values; their record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for byte as
[atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and their identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
