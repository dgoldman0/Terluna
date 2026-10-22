# What sets the snow sizes in the SCALE column

**The controlled comparison produces larger snow on the Moon than in its matched Earth column.** With the native particle laws, the lunar mass-weighted median is 1.59 mm; the most refined Earth counterparts give 0.744–0.746 mm. The corresponding 99th percentiles are 3.11 and 1.028–1.033 mm. The detailed Earth distribution still fails one timestep-refinement criterion, while all tested resolutions retain this size ordering. This experiment therefore finds no mechanism that makes lunar flakes smaller than the matched Earth flakes. It diagnoses the small absolute sizes within a prescribed model population; it does not establish the sizes of actual lunar snow or validate the earlier custom model's decimetre aggregates.

The [diagnostic runner](scale_size_diagnosis.py) and [data product](results/scale_size_diagnosis.json) retain the forcing factors, executable identities, population accounts and comparisons. The original [methods](scale_lunar_methods.md) explain the native adapter and its supplied cloud. This report addresses the physical interpretation of that configuration for the project's snow research.

## Gravity and cloud depth act in the same direction

The calculation separates gravity from depth while retaining temperature, pressure and vapour pressure at corresponding fractional cloud heights. Every gravity/depth control receives 10,232.43 seed crystals/m²/s, the same initial crystal mass and dimension distribution, and a column deposition ceiling of 3 mm/h water equivalent. Each lasts 132.589 hours; the final half supplies arrival statistics. The air remains prescribed, with no resolved wind, turbulence or liquid-droplet population.

| Gravity | Cloud depth (km) | Step (s) | Sample packets per minute per column | Snowfall (mm/h water) | Mass median (mm) | Mass 99th percentile (mm) |
|---|---:|---:|---:|---:|---:|---:|
| Earth | 2.696 | 60 | 4 | 1.722 | 0.714 | 0.991 |
| Earth | 2.696 | 15 | 4 | 1.761 | 0.731 | 1.008 |
| Earth | 2.696 | 15 | 12 | 1.822 | 0.742 | 1.025 |
| Earth | 2.696 | 15 | 24 | 1.837 | 0.744 | 1.028 |
| Earth | 2.696 | 5 | 12 | 1.832 | 0.746 | 1.033 |
| Moon | 2.696 | 60 | 4 | 2.523 | 0.932 | 1.542 |
| Earth | 16.281 | 60 | 4 | 2.784 | 1.077 | 1.880 |
| Moon | 16.281 | 60 | 4 | 2.762 | 1.590 | 3.114 |

Packets are numerical samples representing many physical crystals. Their multiplicities change with sampling rate so the physical source remains fixed. Entries are means of eight column summaries, not pooled percentiles or maximum flake sizes. The lunar baseline's median standard error is 0.0021 mm; the 5 s Earth counterpart's is 0.000083 mm. These sampling errors exclude physical model error and discretization error.

At either depth, reducing gravity increases the calculated dimensions. At either gravity, increasing depth does likewise. In the deep cloud, changing gravity hardly changes snowfall, yet lunar aggregation distributes that water among fewer, larger particles. The lunar run delivers about 1,317 particles/m²/s, averaging 7.76 original crystals per arriving particle; the Earth-gravity deep control delivers 3,782, averaging 2.71. The size difference therefore includes a change in particle population, beyond a change in assigned geometry.

The physical diagonal uses \(H_E=H_M g_M/g_E\). Writing fractional height as \(x=z/H\), the original dry-air hydrostatic relation becomes

\[
\frac{d\ln p}{dx}=-\frac{gH}{R_dT(x)}.
\]

Keeping \(gH\) and \(T(x)\) fixed preserves that pressure path. The diagonal retains the original profile's dry-air approximation and the small difference between its composition and SCALE's native thermodynamic constants. The two crossed cases deliberately impose nonhydrostatic air with dynamics disabled; they isolate processes rather than describe additional equilibrium atmospheres.

## The water allowance constrains the population, not a maximum diameter

The incoming flux was originally obtained from a nominal concentration of 100 crystals/L multiplied by the custom model's reference entry speed. SCALE receives that **flux**; it does not maintain a 100/L native concentration. Its evolving concentrations depend on settling, aggregation and residence time.

Let \(Q\) be the column deposition ceiling and \(F_0\) the incoming crystal number flux. Their ratio is

\[
\frac{Q}{F_0}
=\frac{3/3600\ \mathrm{kg\,m^{-2}\,s^{-1}}}
{10232.43\ \mathrm{m^{-2}\,s^{-1}}}
=8.144\times10^{-8}\ \mathrm{kg},
\]

or **81.44 μg of potential deposited ice per entering crystal**, plus its injected mass. The 3 mm/h ceiling and incoming crystal flux are selected scenario inputs, not measured or established lunar-storm limits. Finite water is a physical constraint; these particular values remain assumptions. Their ratio constrains the average water available to the population. An aggregate can collect many such crystals and their subsequent deposits. In the native lunar baseline, arrivals contain 75.04 μg per constituent crystal, including initial ice, and their mean particle mass is 0.582 mg. The imposed cap is not completely used.

For a steady, collision-free, monodisperse population using a uniform local allowance fully, \(n=F_0/v\) and the local supply is \(Q/H\). Let \(\ell\) measure distance fallen. Consequently,

\[
\frac{dm}{d\ell}=\frac{Q}{HF_0},
\qquad \Delta m=\frac{Q}{F_0}.
\]

A longer residence time and larger concentration compensate in this limit. With aggregation and unequal deposition demand, individual particles receive different shares; the integrated ceiling remains the same.

Likewise, slower differential settling alone does not imply less aggregation. The [native kernel](scale_lunar_methods.md#separating-particle-assumptions) has \(K=E(\sqrt{A_i}+\sqrt{A_j})^2|v_i-v_j|\). If all speeds scale by a common factor \(c\), then \(K\) scales by \(c\), concentration at fixed flux by \(1/c\), and the encounter exposure per falling distance, \(nK/v\), by \(1/c\). Lower speeds then increase exposure. Actual drag and evolving shapes prevent an exact common scaling, but this calculation identifies the concentration and residence-time terms needed to interpret the executed results.

## The particle laws still limit the interpretation of absolute size

The native lunar population has a mass-median apparent density of about 470 kg/m³. Its mostly compact geometry assigns comparatively small dimensions to each particle mass. Applying the old Brown–Francis dimension law to those same masses gives a pooled mass median of 6.07 mm, rather than approximately 1.59 mm. This arithmetic isolates the dimensional mapping; it does not recompute collisions, drag or growth with the alternative geometry.

The separately retained stronger-sticking run further demonstrates that the result depends on the particle laws. Replacing collection efficiency 0.1 with the original temperature-dependent law raises the lunar median to 3.56 mm and the 99th percentile to 7.47 mm. Its arrivals average 59.5 constituent crystals and 3.52 mg per particle. Its mass-median density is about 313 kg/m³. Thus both aggregation and mass-to-dimension conversion matter; neither comparison establishes a physical size ceiling. Even with stronger sticking, SCALE delivers about 172 particles/m²/s, whereas the custom ballistic/no-chipping controls deliver only 2.4–3.5. Applying the old size law to the stronger-sticking SCALE masses gives a median near 18 mm, still far below the custom column's 140–181 mm. The discrepancy therefore also involves the number of mergers and the feedback of geometry on further growth and collection. [Full comparison](scale_lunar_results.md)

SCALE's aggregate outcome assumes filling-in geometry and density-dependent compaction. Its authors identify shortcomings in this treatment, including excessive rounding and suppressed separation of constituent centres. They also report reasonable agreement between their terrestrial simulation's mass–dimension and velocity–dimension relationships and observational formulas. That evidence supports the framework's usefulness while leaving its applicability to this prescribed population open. [Shima et al., §§7.2 and 9.3.8](https://gmd.copernicus.org/articles/13/4107/2020/)

A specific lunar-transfer issue remains: the implemented compaction outcome depends on incoming masses, axes and densities, with a fixed limiting-density parameter, but does not use impact speed, collision energy or contact strength. At fixed incoming masses and shapes, reduced lunar collision energy does not alter that compaction rule; the evolving populations can still produce different outcomes. This is a candidate closure error; the present experiments have not measured its contribution. The 10 kg/m³ parameter itself cannot explain the much higher output densities without an additional sensitivity calculation. [Model equations 59–63](https://gmd.copernicus.org/articles/13/4107/2020/#section4)

## Changing aggregation and water per seed identifies the shared constraint

Disabling aggregation leaves every arriving particle with one original crystal. The lunar column then gives a mass median of **0.856 mm** and a number-mean mass of **0.08220 mg**. The Earth counterpart gives **0.708 mm** and **0.04885 mg** at the same 60 s step. The longer lunar exposure lets deposition nearly use the allowance per entering crystal; it does not generate an unlimited additional water supply. Lunar snowfall is 3.028 mm/h in this control, including separately accounted injected seed ice as well as deposition below its 3 mm/h ceiling.

Restoring aggregation changes the lunar number-mean mass from 0.08220 to 0.5823 mg and the mass median from 0.856 to 1.590 mm. Arrivals average 7.76 constituent crystals. The similar deposited mass per conserved crystal and the loss of particle number through merging account for the increased aggregate mass. The mass median and number-mean mass describe different weightings and need not have the same ratio.

Reducing only the entering crystal flux tenfold raises the allowance to **0.8144 mg per seed**. The resulting median is **3.316 mm**, with 99th percentile **4.676 mm**, while snowfall remains 2.519 mm/h. Arrivals average 3.30 constituent crystals, fewer than in the baseline, but their mean mass rises to 2.256 mg. Thus the larger dimensions in this sensitivity come despite less aggregation per arriving particle: each entering crystal receives substantially more water. Density and shape also evolve, so this is a full population response rather than a fixed-shape cube-root calculation.

Increasing only the water ceiling tenfold gives a **4.532 mm** median and **8.623 mm** 99th percentile in the refined 15 s run. Arrivals average 29.95 constituent crystals and 9.187 mg. Snowfall reaches 11.298 mm/h, well below the 30 mm/h allowance: the model does not capture every available share of the imposed source. Its shorter-step control removes the four clipped collision pairs seen initially and passes the full distribution screen, with a 4.40% change in the 99th percentile and a 1.35-percentage-point exact mass-CDF change. The original small-size range is therefore not a hard diameter cutoff.

## Verification and remaining scope

Thirteen new native runs and two retained references supply this diagnosis. All 45 tests for the affected SCALE adapters and diagnosis pass. The layer check includes tracked and untracked live sources; all 1,248 pass. Eight audited native/adapter source hashes match, and the research provenance check passes. Independent verification reproduces all 65 physics-input files for the new cases byte for byte and agrees with native precipitation histories to the declared relative tolerance of \(10^{-9}\). The new runs close their water ledgers to relative errors below \(6\times10^{-13}\). Apart from the coarse higher-water control's four clipped pairs, all new cases report zero omitted collision events; its refined replacement removes those omissions. Late constituent throughput is within 0.103% of injection. The [check record](results/scale_size_checks.json) retains individual results and late-time variability; the stronger-sticking reference's 99th percentile varies by 14.2% across gathering blocks, so steady bulk statistics do not establish a steady extreme tail.

The refinement screen requires changes no larger than 5% in snowfall and the 50th, 90th and 99th size percentiles, together with an exact mass-CDF distance no larger than 0.05. The CDF comparison gives each of the eight column-normalized distributions equal weight and evaluates every observed diameter. This exact calculation matters for the narrow Earth size distribution: the older sparse diameter grid can miss its largest difference.

| Refinement | Largest relative bulk/percentile change (%) | Exact mass-CDF difference (percentage points) | Screen |
|---|---:|---:|---|
| Lunar inactive capacity 128 → 192 per layer | 0.95 | 0.82 | Pass |
| Earth 60 → 15 s, 4 packets/min | 3.02 | 15.21 | Fail |
| Earth 4 → 12 packets/min, 15 s | 3.47 | 11.09 | Fail |
| Earth 12 → 24 packets/min, 15 s | 0.85 | 4.13 | Pass |
| Earth 15 → 5 s, 12 packets/min | 0.78 | 7.19 | Fail |
| Higher-water Moon 60 → 15 s, 4 packets/min | 4.40 | 1.35 | Pass |

For the last Earth timestep comparison, the sparse grid would report 5.00 percentage points and pass, while the exact maximum is 7.19 at 0.747 mm. The final Earth bulk quantities change by less than 0.8%, but the full distribution has not passed the declared timestep criterion. Further refinement is therefore needed for detailed Earth distribution claims. The lunar median remains more than twice the Earth median throughout these controls; this qualitative attribution does not depend on claiming complete distributional convergence.

A plausible source of the Earth sampling dependence lies in the prescribed-water adapter. It accepts the smaller of each cell's deposition allowance and its sampled positive native growth demand, then discards unused allowance. Fluctuations around that nonlinear minimum can reduce mean uptake. The shallow Earth cloud also holds fewer numerical particles per cell than the deep lunar cloud. Increasing Earth sampling raises uptake and mean mass while leaving constituent counts per arriving particle nearly unchanged, consistent with this mechanism. No cell-resolved demand-versus-cap trace was recorded, so its causal contribution remains unmeasured.

The experiments identify finite water per seed, the number of mergers and the assigned aggregate dimensions as the controls on the shared small-size regime. The active source has no millimetre diameter ceiling or implemented mechanical breakup that could impose a lunar-only limit. Establishing realistic absolute sizes still requires an observed Earth snowfall comparison for this configured population, weather-supported water and crystal sources, and tests of the aggregate geometry and compaction laws. In particular, comparison against an energy-dependent contact model at matched incoming masses and shapes would test whether ignoring lower lunar impact energies materially understates aggregate dimensions.
