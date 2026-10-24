# Warm-column deposition density and early growth constraints

The −5.3 °C benchmark now meets the requested ten-percent limit for mass, width and length from both the assumed 4 µm spherical seed and the fitted three-minute restart. The largest thirty-minute error falls from 23.36% to **9.33%**. A denser comparison at quarter-minute ages from three to thirty minutes, supplemented on both sides of the six-minute reference-fit join, gives a maximum of **9.73%**. All 18 checked cold states remain exactly equal to the preceding −14.4 °C result. This is a successful calibration of those two warm histories; the changed warm fall speeds and strong failures from other seed shapes remain material limits.

The [runner](mcsnow_column_density.py) and [product](results/mcsnow_column_density.json) preserve the [preceding growth-limit experiment](mcsnow_growth_limits.md), its native build and its checkpoint. Temperature is fixed at −5.3 °C, pressure at 101000 Pa and vapour at exact native water saturation. These single-crystal calculations have no population water depletion or atmospheric heat feedback. Mass and dimensions are compared with the [transcribed laboratory fitted histories](takahashi_growth.py), whose population and envelope assumptions carry into the calibration.

## Density and shape controls

The new hypothesis lets the density of newly deposited ice vary with the current enclosing-volume equivalent radius, independently of the existing shape transition:

$$w_\rho=\left[1+(R_\rho/r)^4\right]^{-1},\qquad
\rho_{dep}=\min\left[\rho_{ice},\rho_{native}\,d_{late}\,s_{early}^{1-w_\rho}\right].$$

This changes the envelope volume added per unit new mass. It does not compact existing ice. The rule acts only where the native inherent growth factor exceeds one; the cold plate branch retains its preceding treatment. The existing compact stage still uses solid deposition and a logarithmic polar-to-equatorial growth ratio of one. Its maximum full diameter is now fitted over 0–80 µm, rather than fixed at 40 µm. That diameter is an unmeasured stage limit, distinct from the assumed initial 4 µm seed.

The preferred diagnostic has these controls:

| Control | Value |
|---|---:|
| Common transport multiplier | 0.90715 |
| Early axis-growth multiplier | 2.9824 |
| Later axis-growth multiplier, fixed for this trial | 1.5 |
| Shape-transition equivalent-volume radius | 96.09 µm |
| Compact-stage maximum full diameter | 51.10 µm |
| Later native-density multiplier | 0.76151 |
| Early/later deposition-density ratio | 0.51949 |
| Density-transition equivalent-volume radius | 62.94 µm |
| Directional-ventilation exponent | 1.05 × 10⁻⁸ |
| Mean-ventilation enhancement multiplier | 0.21624 |

After the compact stage, the limiting deposition densities rise from about 0.396 to 0.762 times the native warm-column value. The shape multiplier decreases from its early to later limit. Both transitions use the declared fourth power, inherited for shape and newly tested for density. The warm branching-radius coordinate is fixed at its preceding, independently inactive value. Setting the new density ratio to one while holding the other fitted controls fixed raises the largest thirty-minute error to 20.93% and the largest sampled training error to 38.86%. That ablation establishes a contribution of the density transition within this parameter set; it does not identify a unique physical law.

## Calibration and diagnostic comparisons

Releasing density alone while retaining the 40 µm stage leaves a 14.23% late error. Releasing the stage limit gives a 5.15% late error but a 21.99% early error. Adding the three-minute seed endpoint and a gentler density change brings the initially checked ages within 6.67%, yet the interval check exposes a 13.33% length deficit immediately before six minutes. The separate fitted length segments drop by about 9.9% at that join. The reference describes fitted populations; the jump is not evidence of an individual crystal suddenly shortening.

The final objective therefore includes ages 3, 4, 5, 5.99895833, 10 and 20 minutes. The near-six-minute age is one search timestep before the join. Its weight is 1.25; other fitting ages have weight one. Each starting history has equal total weight, and mass and both dimensions have equal weight within each age. The three-minute restart is a supplied initial condition and contributes no fitting residual. An unweighted join trial gives a 10.33% interval maximum; weight 1.5 transfers the largest discrepancy to about nine minutes and gives 10.18%. Weight 1.25 balances these recorded discrepancies at 9.73%.

Fifteen and thirty minutes remain outside every numerical fitting objective. Ages 6, 8, 12 and 25 are also outside the final objective. All inspected comparisons have informed model design, selection or weighting, so none is presented as untouched independent validation. Scores from different fitting-age or weight families are kept separate. The final weighted RMS log error, 0.03945, is not directly comparable to the preceding score on fewer ages.

| Diagnostic | Largest 30-minute error | Largest reported wider-age error |
|---|---:|---:|
| Preceding warm fit | 23.36% | 27.56% at the original sampled ages |
| Evolving density, fixed 40 µm stage | 14.23% | 27.56% at the original sampled ages |
| Evolving density, fitted compact stage | 5.15% | 21.99% at the original sampled ages |
| Three-minute constraint, density ratio ≥0.5, native later shape limit | 6.10% | 13.33% in the interval check |
| Early constraints, later shape multiplier 1.5, unweighted join | 8.68% | 10.33% in the interval check |
| Same join with weight 1.5 | 9.84% | 10.18% in the interval check |
| Preferred join weight 1.25 | **9.33%** | **9.73%** in the interval check |

These maxima cover mass and both dimensions for the two fitting histories; the wider-age column identifies its sampling because the early trials were not all evaluated on the same dense grid. At thirty minutes the preferred spherical-start errors are +1.68% mass, +0.03% width and −1.15% length; the restart errors are −3.43%, +1.99% and −9.33%. The worst dense-grid discrepancy is restart length just before six minutes. The largest spherical-start discrepancy is +9.49% length near 9.5 minutes. This is sampled agreement, with a narrow margin below ten percent, rather than a bound over every time or environmental condition.

![Warm-column mass, width and length errors for the retained diagnostic comparisons](../../visualization/snow/generated/column_density.png)

The [figure source](../../visualization/snow/column_density.py) reads the product and checks plotted values against it. Green shading marks ±10%. Lines connect saved ages; the abrupt error change near six minutes follows the reference-fit join. The cold result is reported separately because its controls and temperature are separately prescribed. A [PDF](../../visualization/snow/generated/column_density.pdf) accompanies the PNG.

## What remains unresolved

The cold comparison retains its earlier thirty-minute maximum error of 4.75% and its earlier three-minute mass deficit of about 17%. Its exact parity follows from freezing the cold controls and leaving the native plate-density branch unchanged. There is no inferred interpolation between the two fitted temperatures.

Warm fall speed was excluded from the objective and remains a failure. Its thirty-minute error rises from about +3.2% in the preceding warm fit to +21.4% for the spherical start and +23.0% for the restart. Native drag is unchanged; the fitted crystal properties now produce different speeds. Thus improved mass and dimension histories do not jointly validate growth and settling.

Other geometry assumptions also fail. The volume-preserving restart, excluded from fitting, has a largest sampled mass-or-dimension error of 24.21%. At the same initial mass and volume, a 2:1 seed gives 7.20 times the fitted thirty-minute length; a 4:1 seed gives an unphysical runaway to about 356,000 times fitted length. The compact stage preserves an arbitrary seed's aspect ratio rather than forcing it toward one. The new fit therefore worsens a known seed-shape sensitivity while meeting the requested benchmark under its two stated starts.

Supplying the fitted crystal state directly still leaves instantaneous-rate errors: at 3.5 minutes the width rate is +41.7%; at 5.5 minutes the length rate is −23.5%; at thirty minutes uptake is +4.9%, width rate +20.5% and length rate −8.5%. The smooth transitions can compensate accumulated trajectory errors without matching every local growth rate. At 3.5 and 5.5 minutes, the reference-implied deposition density decreases from 334 to 247 kg/m³, whereas the model at those supplied states increases from 259 to 312 kg/m³. Those reference densities derive from fitted population histories and the chosen envelope mapping; they are not direct porosity measurements. Their opposite trend is a specific reason to retain the physical early-growth problem despite the closer integrated histories. Physical surface evolution, a credible compact-to-faceted transition, drag and independent humidity constraints remain necessary before transfer. The existing lunar calculations and aggregate predictions are unchanged.

## Reproduction and checks

The native calculation uses McSnow revision `1e428d7e901cbc5b91f2040e575f736ec5963a74` with the separately recorded diagnostic patches, driver, compiler and tables. The original authors' generating configuration remains unidentified, as recorded in the preceding source assessment. No new full-paper reading or mechanistic facet-spreading implementation is claimed here.

Seventeen fits retain the original local attempt, two 35-generation differential-evolution searches, compact-stage and density-bound sensitivities, early-age constraints and weighting controls. Each global population has six times the number of free coordinates and reaches its declared generation limit; bounded local refinements permit up to 400 evaluations. These finite searches establish recorded candidates, not a global optimum. The final source records each actual specification and all failed candidates.

The product records 7,769 parameter-batch evaluations, including 254 failed evaluations, and all 104 numerical checks pass. One retained interval attempt requested an off-grid output time and was rejected by the native driver; the corrected check uses aligned times. The largest measured refinement difference across the experiment is 0.3754%. For the preferred fit the dense interval and all-start comparisons change by at most 0.0546% between 1/32 and 1/128 s; the unobserved early history changes by at most 0.3086%. Mass-ledger and spheroid-volume residuals stay below 2.23 × 10⁻¹⁶ and 2.34 × 10⁻¹⁵. The two neutral-limit comparisons each reproduce 18 preceding states exactly, and every selected density family preserves all 18 cold states exactly.

Verification completed with 55 focused tests and a passing `make check`: 3,776 Python tests passed, 20 skipped, and 191 JavaScript tests passed. Thirty plotted curves were checked against the numerical product, and the PNG and PDF were visually inspected.

The [verification record](results/mcsnow_column_density_checks.json) identifies the exact product, independent objective recalculation, repository checks and figure inspection. Raw simulation journals, build records and a copied hash-verified checkpoint live on the project-space drive through the existing cache link. Earlier checkpoints and products retain their original bytes.

To reproduce from the pinned parent cache, use one numerical worker and run the actions in order:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
for action in build parity local global compact density_bounds seed_state late_control density_floor early_interval interval_profile joint_early joint_late joint_join join_weights analyze
do
  python -B -m atmosphere.snow.mcsnow_column_density "$action" || break
done
python -B -m visualization.snow.column_density
make check
```
