# Earth deposition and settling against the Welss archive

[welss_column.py](welss_column.py) computes an independent, deterministic counterpart of the **oblate-labelled diagnostic-geometry control `ob_hp0`** in [Welss, Siewert and Seifert (2024)](https://doi.org/10.1029/2023MS003805), Figures 17 and 19. It follows the public McSnow revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`. The completed comparison gives integrated absolute differences of 0.55–1.77% across the four principal profiles; number density differs by 7.62% within the nucleation layer. These are comparisons between calculations under declared source equations. The experiment supplies neither empirical validation nor a replay of the unidentified executable that produced the archive.

## Prescribed source experiment

Ice nucleates uniformly between 2,700 and 3,000 m at 16.5 particles m⁻³ s⁻¹. Initial mass follows a gamma distribution with shape two and mean 4.8 × 10⁻¹³ kg. The source diagnostic geometry, monomer capacitance, Böhm terminal speed and ventilated depositional growth act on independent particles in a prescribed atmosphere. Aggregation and riming are disabled; vapour and temperature remain prescribed. The experiment lasts ten hours.

The atmosphere uses the pinned source's 250 layers over 0–5 km, its ideal initialization and subsequent linear interpolation with water-saturation adjustment. The source inputs are `h1=0`, `h2=600 m`, `ssat=0.05`, and `lwc0=0.0002 kg/m³`; `rt0` disables riming despite liquid water in the atmosphere. This preserves the source's pressure calculation, including its use of local water saturation pressure in the variable named `pv_srf`. Constants copied from the author model are identified as reference inputs, separate from Terluna's physical defaults.

[earth_benchmarks.md](earth_benchmarks.md) records the archive identities and recovered setup. The archived atmosphere supports 5% supersaturation over ice, and the source and run labels support enabled base ventilation. Exact run provenance remains unavailable: the profiles were processed in February 2023, public code history starts in March 2024, and atmospheric interpolation changed in July 2024. The available source averages four final snapshots, while the article specifies the last five hours. Arithmetic in the archived output supports a 31-snapshot mean, consistent with the latter interval at the stated ten-minute cadence. It does not identify the exact snapshot times or executable.

## Independent population calculation

For each nucleation height and initial mass, the solver integrates mass and age along descending height:

\[
\frac{dm}{dz}=-\frac{G(m,z)}{v(m,z)},\qquad
\frac{da}{dz}=-\frac{1}{v(m,z)}.
\]

It advances log mass to preserve positivity. Generalized Gauss–Laguerre quadrature integrates the gamma mass distribution over its infinite support; Gauss–Legendre quadrature integrates nucleation height within each 20 m source layer. The seed weights remain fixed. No particle distribution is repeatedly remapped onto a mass grid.

Continuous injection makes the arrival age determine whether a trajectory contributes at a given time. The five-hour-average weight is `clip((36000−age)/18000, 0, 1)`. The snapshot weight is the fraction of times 34,200, 34,800, 35,400 and 36,000 s reached by that trajectory. The final-time result is retained separately. This is the continuous-time expectation of the prescribed source rate; McSnow's five-second stochastic injection and finite-step process splitting are not replayed.

The source writer bins particles over `[z,z+20 m)` and calculates number density, mass density and their velocity-weighted moments. Accordingly, the characteristic calculation integrates `J W/v`, `J m W/v`, `J W` and `J m W` for `Nd`, `Md`, `Fn` and `Fm`, respectively, over each output bin and source distribution. `Fn` and `Fm` are these spatially averaged velocity moments, including inside the nucleation region. Deposition is integrated as `J G W/v`, in kg m⁻³ s⁻¹. For number density, `da/dz=−1/v` converts the bin integral into the difference of the arrival-weight antiderivative at the two endpoint ages. This avoids spatial quadrature error from the low speed of newborn crystals; the remaining moments retain independent bin-quadrature controls.

Raw model profiles and an exploratory five-point Lanczos pass are both retained in the frozen [physical product](results/welss_column.json). Inspection of the archive then established that its totals retain finite 0, 20, 2,960 and 2,980 m bins: a full `opt=0` pass would remove the first two, and the source zero-to-missing convention would also remove the last two. The separate [archive comparison](results/welss_archive_comparison.json) therefore compares raw profiles directly. Its metrics supersede the one-sided-filter diagnostics in the physical product. Figure 19 applies filtering at plotting time; reproducing that figure would require applying the same treatment to both curves.

All 150 finite archived number-flux values fit the arithmetic of a 31-snapshot mean of three-significant-digit inputs within two float32 units in the last place. Multiples such as 62 and 93 also fit; four snapshots fit only seven values. Together with the publication, this supports a five-hour averaging interpretation without recovering the exact processing configuration. No NCL list-aliasing bug is inferred from the presence of a filter block in later source code.

## Verification and reproducibility

[Independent Fortran fixtures](references/welss_microphysics_reference.json) exercise 44 source cases, including tiny seeds, the geometry transition and signed growth. Numerical controls check gamma normalization, analytic growth trajectories, constant-speed source-bin moments and finite-age averaging. The completed runs vary mass quadrature, source-height quadrature, bin integration and ODE accuracy independently. Compiled NCL fixtures check filter arithmetic and missing-value handling separately from whether the archive actually used that pass.

The runner checks the reference member's SHA256, records producer and source identities before integration, and checks them again afterward. External NetCDF bytes remain outside Git; result files contain model profiles, derived comparison errors, source credit and reference identity. Source reading extents and hashes are in [welss_column_sources.json](welss_column_sources.json).

```bash
python -m atmosphere.snow.welss_column \
  --reference /path/to/ob_hp0.nc \
  --output /path/to/welss_column.json
```

The five-case suite is reproduced with `python -m atmosphere.snow.welss_comparison --reference /path/to/ob_hp0.nc`; the separately bound raw comparison uses `python -m atmosphere.snow.welss_archive_comparison --reference /path/to/ob_hp0.nc`. Each completed case retains its raw profiles, source identities and settings. External author arrays remain outside Git.

## Results and numerical controls

The following relative L1 differences divide the sum of absolute model–archive differences by the archived profile sum. Uniform 20 m bins make this an integrated comparison. The full region contains 150 bins with lower edges 0–2,980 m; the nucleation belt contains the upper 15. These values describe the baseline calculation.

| Profile | Full populated region | Below nucleation belt | Within nucleation belt |
|---|---:|---:|---:|
| Number flux `Fn` | 0.6645% | 0.6683% | 0.5968% |
| Mass flux `Fm` | 0.8698% | 0.8701% | 0.4622% |
| Number density `Nd` | 1.7737% | 0.4958% | 7.6234% |
| Mass density `Md` | 0.5473% | 0.5472% | 0.5508% |

The spatially integrated agreement does not establish uniform local agreement. At the top populated bin, number density is 32.3% above the archive and mass flux differs by 17.8%. Number-density differences below the nucleation belt remain below 1.87% pointwise. The source-region excess therefore remains a specific unresolved discrepancy.

| Independent control | Change from baseline | Largest integrated raw-profile change |
|---|---|---:|
| Seed mass quadrature | 16 → 32 nodes | 0.0000684% |
| Nucleation height quadrature | 2 → 4 nodes per source layer | 0.005364% |
| Output-bin quadrature | 3 → 5 nodes | 0.001932%, for deposition |
| Trajectory integration | relative tolerance 10⁻⁷ → 10⁻⁹; maximum step 20 → 10 m | 0.0001617% |

The last column covers all five retained moments, including deposition. One refinement per variable measures sensitivity; it establishes neither a formal convergence order nor a rigorous error bound. Local sensitivity is larger near the upper source boundary: height refinement changes the top-bin number density by about 0.85%, while its excess over the archive remains 31.2%. The unresolved source-region discrepancy survives the tested numerical refinements.

Across the five cases, sampled ground-arrival ages range from 1.69 to 1.92 hours. Consequently the continuous last-five-hour, four-final-snapshot and final-time deterministic profiles coincide. The archived averaging window still affects stochastic sampling. The five physical integrations took 267 seconds in total, using one worker; their source and input identity checks passed.

## Compact crystal capacitance

One [separately bound sensitivity](results/welss_compact_capacitance.json) sets `C=D/2` only below the source's existing solid-sphere mass threshold, 1.70527 × 10⁻¹¹ kg. The tested source instead chooses `D/π` there because its calculated void volume is slightly positive. The intervention uses the sphere's geometric capacitance without fitting a parameter; diameter, area, speed, atmosphere, ventilation and nucleation remain unchanged. It retains baseline quadratures and uses relative ODE tolerance 10⁻⁹.

| Number-density comparison | Literal-source baseline | Compact-sphere capacitance |
|---|---:|---:|
| Full-profile relative L1 | 1.7737% | 0.9442% |
| Nucleation-belt relative L1 | 7.6234% | 3.1426% |
| Top-bin relative difference | 32.30% | 14.96% |
| Below-belt relative L1 | 0.4958% | 0.4640% |

The change materially reduces the source-region discrepancy, while a substantial local difference remains. It does not improve every quantity: full-profile mass-flux error changes from 0.8698% to 0.9196%, and mass-density error from 0.5473% to 0.5710%. This identifies a relevant source-geometry sensitivity without establishing which branch the archived executable used. The intervention has no independent quadrature-convergence study and is not adopted into the lunar model.

Two initial attempts encountered invalid internal log-mass trial stages and produced no admitted result. [The startup record](results/welss_capacitance_startup.json) retains their log hashes and checks a halved first step at the formerly failing source height. At relative tolerance 10⁻⁹, halving that step changes ground mass by at most 6.14 × 10⁻⁷ relative and age by 2.17 × 10⁻⁷, or 0.00141 s. These are measured sensitivities, not global error bounds.

The isolated driver bounds its initial trial step and returns a nonfinite derivative for unrepresentable trial states, allowing the RK error-estimation loop to reject them without clipping mass. SciPy's additional dense-output stages occur outside that rejection loop, so the wrapper is not a general safeguard. The admitted 104-second run encountered **zero invalid stages**, and all 3,765 saved profile values are finite. Future reuse should additionally check every dense-output coefficient.

This calculation tests deposition and settling under a published terrestrial configuration. Aggregation, lunar aggregate survival and empirical validation remain separate work. The frozen local lunar column has not adopted these reference closures.
