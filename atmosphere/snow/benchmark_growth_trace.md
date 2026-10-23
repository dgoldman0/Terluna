# Separating crystal growth rate from shape evolution

The new source trace identifies excessive vapour uptake as the principal amplification of the revised plate's ten-minute archive discrepancy. At −14.5 °C, removing only the extra oblate mean-ventilation term reduces the refined IGF2+ mass error from **49.10% to 9.21%** and width error from **25.91% to 6.35%**. At equal mass, the native and altered crystals have almost identical dimensions. The change slows their progress along essentially the same sequence of shapes.

This finding follows the [benchmark reconstruction](benchmark_reproduction.md) and [drag trace](benchmark_drag_trace.md). It supports investigating the uptake treatment while retaining the verified native drag. It does not reproduce the authors' original experiment: every tested candidate still differs from the numerical archive beyond its printed precision. The [product](results/benchmark_growth_trace.json) contains 60 completed executions, including archived-state restarts, actual-source traces, two explicit interventions and a thirty-minute laboratory comparison. No replacement growth law is adopted.

## Earliest available departure and the restart control

The earliest matching archived outputs are at three minutes; they cannot locate the first difference within those three minutes. At that age, the refined native IGF2+ reconstruction already has 5.90% excess mass and 2.74% excess width. Its thickness is 0.29% high. The crystal has not branched: the source first applies the porous deposition density at 332.203125 s. IGF1 also differs at three minutes, with 5.86% excess mass, 4.49% excess width and 3.05% deficient thickness.

Eight restart calculations supply the archived three-minute mass, envelope volume and aspect ratio, reconstruct the full dimensions, and continue for another 420 s. The fine native IGF2+ restart ends with 43.21% excess mass and 23.02% excess width, compared with 49.10% and 25.91% from the original sphere. Thus the early difference contributes to the endpoint error, but removing it leaves most of the later amplification. The corresponding IGF1 restart gives 11.65% excess mass and 4.91% excess width. The historical candidate remains too small in both variants. These restarts use rounded author **simulation outputs**, not an observed laboratory crystal followed through time.

## Equal mass distinguishes the growth rate from the shape path

The refined native IGF2+ crystal reaches the archived ten-minute mass at **535.677 s**. At that instant its width is only 1.21% high, thickness 2.09% low and volume 0.29% high. The much larger age-matched width error therefore mainly reflects faster mass accumulation. Those remaining shape errors still exceed archive rounding, so equal-mass agreement is a diagnostic rather than an age-matched reproduction.

There is a structural reason for the control result. For a pristine monomer at constant temperature, the source updates volume by deposited mass divided by deposition density. Its continuous-time limit gives

\[
\frac{dV}{dm}=\frac{1}{\rho_{\rm depo}},\qquad
\frac{d\ln\phi}{dm}=\frac{\Gamma_{\rm eff}-1}{\Gamma_{\rm eff}+2}
                       \frac{1}{\rho_{\rm depo}V}.
\]

For IGF2+, both deposition density and effective growth ratio are functions of the crystal and prescribed air state. Mean ventilation changes the uptake rate but does not enter these two relations: the directional calculation receives the Reynolds–Schmidt combination evaluated before mean ventilation modifies diffusivity. Consequently, removing the extra mean term changes elapsed time while preserving the limiting shape-versus-mass path. At the archived final mass, the fine native and altered widths differ by only 0.000050%, while their arrival times differ by 45.993 s. Their finite-step differences decrease under refinement. IGF1 has an additional dependency because its branching criterion uses the modified diffusivity; its width difference at equal mass is 0.0217%.

The mass-volume account supplies another conditional check. If a plate first deposits solid ice and then undergoes one lasting transition to density `rho_i Gamma`, its transition mass is

\[
m_b=\frac{m-\rho_i\Gamma V}{1-\Gamma}.
\]

Applying that assumption and the native constant-temperature growth factor to the archived IGF2+ endpoint gives `1.0949345 × 10⁻⁹ kg`, with a printed-rounding half-bound of `8.69 × 10⁻¹⁶ kg`. The native source trace branches at `1.0810928 × 10⁻⁹ kg`, only 1.26% lower. The same inference from the native final state recovers its traced transition mass to numerical precision. This constrains the branching explanation under the stated density law; it neither measures an archived branching event nor proves that the missing original driver used that law.

## Two distinct source interventions

The [runner](mcsnow_growth_trace.py) instruments the pinned 2024 deposition source without changing the native candidate's arithmetic. It separately compiles the historical two-operation deletion, removal of only the extra oblate mean term, and replacement of `V^(1/3)` by the physical envelope-volume sphere radius `(3V/4pi)^(1/3)` in directional ventilation. The latter follows the radius convention discussed in the [paper's directional-ventilation description](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS003805); it is a convention sensitivity, not an author-confirmed correction. It retains envelope volume rather than introducing a mass-equivalent volume.

At −14.5 °C and ten minutes, using a 1/64 s timestep, relative errors against the author archive are:

| Growth package and candidate | Mass | Full width | Full thickness |
|---|---:|---:|---:|
| IGF2+, native | +49.10% | +25.91% | +1.42% |
| IGF2+, historical, no mean ventilation | −20.93% | −11.74% | −4.33% |
| IGF2+, omit extra oblate mean term | +9.21% | +6.35% | −1.28% |
| IGF2+, physical directional radius | +52.07% | +28.67% | −0.61% |
| IGF1, native | +20.48% | +11.01% | −1.23% |
| IGF1, historical, no mean ventilation | −23.25% | −10.68% | −4.89% |
| IGF1, omit extra oblate mean term | +5.81% | +4.29% | −2.28% |
| IGF1, physical directional radius | +21.62% | +12.46% | −2.79% |

Sixteen fixed-state executions distinguish the direct effect of each intervention from the resulting trajectory feedback. At the archived ten-minute IGF2+ state, native mean ventilation multiplies vapour diffusivity by 1.70747 and heat conductivity by 1.75502; uptake is `2.64625 × 10⁻¹¹ kg/s`. Omitting the extra oblate term gives factors 1.31417 and 1.33473 and uptake `2.03117 × 10⁻¹¹ kg/s`, a 23.24% reduction. Capacitance, fall speed, deposition density and effective growth ratio remain identical at that supplied state. The historical deletion reduces uptake further to `1.54006 × 10⁻¹¹ kg/s`.

The physical-radius intervention instead preserves instantaneous uptake while changing the effective growth ratio from 0.166262 to 0.151727. The resulting plate becomes thinner and wider along its trajectory and the mass discrepancy grows. This separates an uptake control from a directional-growth control; neither is selected merely because one output approaches the reference.

## The next duration control

The extra-term removal was continued to thirty minutes at the laboratory comparison temperature of −14.4 °C, using the same 2 µm radius spherical seed, 101000 Pa and exact native water saturation. This is a comparison with the published fitted laboratory histories, distinct from the −14.5 °C numerical archive comparison above.

| Age | Full width | Width / laboratory fit | Mass / laboratory fit | Thickness / laboratory fit |
|---|---:|---:|---:|---:|
| 3 min | 0.23630 mm | 0.5944 | 0.8673 | 1.5895 |
| 10 min | 1.26432 mm | 0.9288 | 0.8436 | 1.2039 |
| 20 min | 3.35257 mm | 1.2145 | 1.1185 | 0.9639 |
| 30 min | 5.64284 mm | 1.3517 | 1.2600 | 0.8253 |

The retained native reconstruction reaches 113.515 mm at thirty minutes, against a fitted 4.175 mm. Removing the extra mean term suppresses that severe amplification but leaves a 35.17% width excess, 26.00% mass excess and 17.47% thickness deficit. Refining this control from 1/32 to 1/128 s changes the final width by 0.0073% and mass by 0.0144%. The remaining errors are therefore substantially larger than the measured numerical changes.

This establishes a useful boundary for subsequent work: the extra oblate mean term drives most of the late plate amplification, while a rate-only change leaves a residual shape-history problem and fails the joint laboratory comparison. The earlier multi-habit [ventilation controls](earth_growth_assessment.md#what-the-ventilation-controls-establish) retain their separate scope and evidence. A complete repair must satisfy mass and both dimensions across habits and ages before supplying an accepted population growth law. Finite-water, latent-heat, settling and collision controls remain separate requirements for the storm and lunar column; prescribed single-crystal growth does not establish their outcome.

## Reproduction and verification

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -B -m atmosphere.snow.mcsnow_growth_trace run
python -B -m atmosphere.snow.mcsnow_growth_trace assemble
```

The 60 executions comprise 32 seed trajectories, 16 fixed-state evaluations, eight archived-state restarts, two thirty-minute extensions and two output-cadence controls. All raw files reside under `research/runs/snow_scale/mcsnow_reproduction/growth_diagnosis/` on the relocated project-space drive. Source, compiler, executable, table, namelist and raw-file identities are retained. The original sources, earlier runners and preceding results remain intact.

Sixteen native/historical traced trajectories exactly reproduce all common saved states and air fields from their preceding executions. Each trace labels the pre-step evaluation time and the corresponding post-step state time. Its rate, deposition density and effective growth ratio are checked against the actual mass, volume and aspect-ratio increments. Seed trajectories retain the original 1 s timestep plus 1/16, 1/32 and 1/64 s controls; the largest final refinement change across compared fields and both ages is below 0.125%. Reducing output spacing from 0.25 to 0.0625 s changes selected IGF2+ matched-mass fields by at most `3.13 × 10⁻⁷` relative and arrival times by less than 0.000036 s, with exact equality at common output times.

The [check record](results/benchmark_growth_checks.json) distinguishes software checks, raw-output replay and substantive source/report reading from empirical validation. Forty focused tests pass. The repository-wide `make check` passes with 3,655 Python tests passed and 20 skipped, all 191 JavaScript tests passed, and provenance and ensemble integrity passed. A follow-up layer scan after staging includes both new Python files: 1,312 files checked, zero violations. The missing original generating configuration and the remaining archive and laboratory discrepancies are explicit outstanding evidence issues.
