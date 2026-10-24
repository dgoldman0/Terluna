# Independent early and later crystal-growth limits

Allowing the early axis-growth multiplier to differ from its later value resolves much of the cold dendrite's conflict between the assumed spherical seed and the fitted three-minute restart. At −14.4 °C, training RMS log error falls from 0.1828 to 0.0267, a 97.9% reduction in squared-log loss. At the reserved thirty-minute age, mass and both dimensions are within 4.8% for both starts. The corresponding warm-column improvement is smaller: RMS log error falls from 0.0821 to 0.0741. The constrained refits below show that early shape alone does not reproduce these gains; deposition density, branching and transport remain consequential.

These are diagnostic fits to two temperature-specific laboratory histories. They identify a limitation of the preceding interpolation and quantify which additional adjustments this model uses. They do not identify measured attachment coefficients or establish a transferable lunar growth law. Unknown initial geometry and the fitted-history envelope convention remain material uncertainties.

## What changed

The [preceding state-growth experiment](mcsnow_state_growth.md) fixed the early multiplier at one. Its independent shape budget showed why that was restrictive: the assumed 4 µm sphere and the fitted three-minute cold state imply an equivalent logarithmic axis-growth ratio of 0.145–0.158, whereas the first measured-age segment implies 0.547. That geometric constraint does not depend on how quickly volume is supplied. The [new runner](mcsnow_growth_limits.py) retains McSnow revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`, the improved plate component, the six diagnostic controls and the fourth-power radius interpolation, while releasing its early limit:

$$q(r)=q_{early}^{1-w}q_{late}^{w},\qquad w=\left[1+(R/r)^4\right]^{-1}.$$

Here $r$ is the enclosing-volume equivalent radius. The multiplier acts on the native logarithmic axis-growth ratio, separately from the retained directional-ventilation exponent. It is not a molecular attachment probability or a ratio of linear facet velocities. Setting $q_{early}=1$ recovers the preceding state model; setting $R=0$ recovers the earlier constant adjustment. The warm tests retain the declared 40 µm maximum-full-diameter compact stage, during which deposition is solid and the logarithmic growth ratio is one. This preserves an initially spherical seed; it does not turn an arbitrary seed into a sphere or implement lateral facet spreading.

Both early and later multipliers initially range from 0.25 to 4, with logarithmic search coordinates; the transition radius ranges from 0 to 150 µm. The cold bound profiles separately lower the early bound to 0.0625 or raise the mean-ventilation enhancement bound from 2 to 4. Every search records its actual bounds, fixed coordinates and initialization. The native source and the earlier products remain unchanged in separate builds.

## Fitted controls and reserved histories

Two independent global search starts at each temperature return closely similar scores and active parameters. The warm branching-radius values differ because that control is inactive for these column trajectories, verified here by changing its multiplier to 0.25 and 4 with exactly unchanged sampled outputs for both training starts. The default-bound selections are:

| Control | −14.4 °C | −5.3 °C |
|---|---:|---:|
| Transport multiplier | 0.71884 | 0.92813 |
| Later axis-growth multiplier | 1.9577 | 2.6782 |
| Directional-ventilation exponent | 0.13653 | 4.0605e-06 |
| Porous deposition-density multiplier | 1.394 | 0.53824 |
| Branch-radius multiplier | 0.5968 | 0.86421 |
| Mean-ventilation enhancement multiplier | 2 | 0.082277 |
| Transition equivalent-volume radius (µm) | 24.697 | 20.845 |
| Compact-stage maximum diameter (µm) | 0 | 40 |
| Early axis-growth multiplier | 0.25 | 0.25 |

The cold fit asks for strongly flattened early growth followed by more thickness growth as the crystal develops. Its early multiplier reaches the lower bound and its mean-ventilation multiplier reaches the upper bound. The warm fit also chooses the lower early bound, an almost absent directional correction and little mean-ventilation enhancement. Initializing the warm searches with a larger early than later multiplier does not retain that proposed declining prescription at the identified minima. The warm fitted multiplier instead rises with size, leaving the early/late rate conflict described below.

| Temperature | Start | Mass at 30 min | Width at 30 min | Polar dimension at 30 min | Speed at 30 min |
|---|---|---:|---:|---:|---:|
| -14.4 °C | sphere | +1.0% | -4.4% | -4.5% | -5.5% |
| -14.4 °C | axes | +1.2% | -4.8% | -3.6% | -5.0% |
| -14.4 °C | volume | +0.4% | -12.9% | -7.7% | -3.8% |
| -5.3 °C | sphere | +0.5% | +2.1% | +23.4% | +3.3% |
| -5.3 °C | axes | -4.2% | +2.2% | +17.2% | +3.2% |
| -5.3 °C | volume | -7.7% | -2.6% | +0.4% | +5.2% |

The cold default fit retains a -17.0% mass error at three minutes, while its two dimensions are within 2.3%. The warm three-minute mass and length errors are +21.6% and +27.6%. Its thirty-minute length remains +23.4% from the sphere and +17.2% from the restart. Thus the smaller warm training score improves some reserved quantities and worsens others; it does not remove the warm trajectory discrepancy. The volume-preserving restart is an additional representation sensitivity, excluded from fitting.

## Which adjustments remain necessary in these searches

Three constrained families restore selected coefficients to their original values and refit the remaining controls. Each retains the temperature-specific compact-stage setting above. Their RMS log errors compare with the full default-bound adjustment:

| Free controls | −14.4 °C | −5.3 °C |
|---|---:|---:|
| Early/later axis multipliers and transition radius | 0.25093 | 0.21523 |
| Those shape controls, density multiplier and branching threshold | 0.08430 | 0.20823 |
| Those controls plus directional ventilation; original transport and mean ventilation | 0.08430 | 0.19409 |
| All eight controls | 0.02665 | 0.07412 |

Allowing density and branching removes much of the cold shape-only discrepancy. The warm case benefits much less while original transport and mean ventilation are retained: its best constrained fit overpredicts reserved thirty-minute mass by 128% from the sphere and 93% from the restart. At the supplied warm thirty-minute geometry, the same original coefficients give a mass rate 121% above the fitted rate, independently of the preceding trajectory. These joint restrictions identify compensation within the tested model; they do not isolate a unique error in any one physical coefficient. In particular, adding control freedom can expose a different optimum rather than measure the importance of that control by itself.

Finite-search behavior is checked explicitly. A less constrained search initially scored worse than an available nested-family solution. The original result is retained, and named local restarts begin from the better constrained candidates:

- `axis_m14p4_seed0_20261027_native_transport_mean_nested`: 0.084303 at the supplied nested solution → 0.084303 after refinement.

Each restart records its source fit and hash, preserves the fixed controls, and keeps the supplied candidate if refinement does not improve it. These controls prevent an optimizer exit from being mistaken for evidence against a permissible solution. They do not prove a global optimum.

## Bounds and supplied-state rates

| Cold search bounds | Training RMS log error | Early multiplier | Later multiplier | Transition radius (µm) | Mean-ventilation multiplier |
|---|---:|---:|---:|---:|---:|
| Default | 0.026654 | 0.25000 | 1.95773 | 24.697 | 2.00000 |
| Early lower bound 0.0625 | 0.026169 | 0.07222 | 1.90345 | 15.552 | 2.00000 |
| Mean-ventilation upper bound 4 | 0.021508 | 0.25000 | 1.88575 | 23.906 | 2.47846 |

Lowering the early bound changes the inferred early multiplier and transition radius substantially for little score improvement; the later multiplier stays relatively close. The overall early shape change is better constrained than its detailed path. Widening the ventilation range reduces training error further but increases reserved thirty-minute mass error to 7.1–7.5% for the two training starts. The bound profiles remain separate comparisons; reserved measurements are not used to select a replacement search range.

Resetting the early multiplier to one while holding each default fit's other controls fixed raises RMS log error from 0.0267 to 0.5150 in the cold case and from 0.0741 to 0.1459 in the warm case. Removing the warm compact stage while retaining its fitted controls raises RMS log error to 0.4819. These are fixed-parameter ablations, distinct from refitting either restricted family; the 40 µm threshold remains an assumed diagnostic value.

Supplied-state probes separate instantaneous growth from accumulated trajectory error. The following rates begin at the fitted thirty-minute state, preserving the fitted axes; width and polar entries compare logarithmic rates:

| Temperature | Mass rate error | Width rate error | Polar rate error | Deposition density: model / fitted (kg m⁻³) |
|---|---:|---:|---:|---:|
| -14.4 °C | +9.3% | -11.3% | -16.9% | 357.5 / 286.6 |
| -5.3 °C | -9.4% | +6.6% | +34.0% | 326.7 / 445.6 |

The cold thickness-rate error falls from 63.8% low under the preceding selected fit to 16.9% low here. At 3.5 minutes, however, its mass uptake is still 27.0% high at the supplied state. The warm effective logarithmic axis ratio rises from 3.33 at 3.5 minutes to 4.05 at thirty minutes, while the laboratory fitted segments imply 5.42 and 3.22. The warm model deposits ice at about 327 kg m⁻³ throughout these probes; the fitted implied deposition density is about 247 kg m⁻³ at 5.5 minutes and 446 kg m⁻³ at thirty minutes. These opposite tendencies identify a remaining growth-allocation and density constraint. The abrupt change between independently fitted segments is not an observed physical switch.

## Initial geometry and evidence limits

The same fixed-mass, fixed-volume seed-shape controls remain consequential. For the warm fit, initial aspect ratios 2 and 4 give 2.46 and 4.21 times the fitted thirty-minute length, improved from 4.67 and 80.98 under the preceding fit. The cold 0.5-aspect seed gives 1.53 times the fitted mass and 0.65 times the fitted thickness. All selected controls and any failed sensitivity runs remain in the product; these aspect ratios are sensitivity inputs, not measured seed distributions.

The laboratory comparison begins at three minutes, uses fitted populations of different crystals, and does not supply the assumed initial sphere. Some thickness estimates also depend on mass and projected area, so the three fitted quantities are not three fully independent measurements. Reserved ages come from those same fitted histories. The original author executable/configuration remains unidentified; the native source revision, driver, patches, compiler and tables define the actual reproducible experiment here. The [preceding source record](mcsnow_state_sources.json) retains the limited article reading and the recovered low-humidity facet-area data. Complete retrieval of the 2021 publisher method remains unresolved, and no facet-spreading implementation is claimed.

A physical successor needs joint constraints on early surface state, axis allocation and the density of added ice, along with credible uptake and ventilation under the observed geometry. The cold improvement supplies a much stronger diagnostic comparison for that work. The warm trajectory and supplied-state rates, unknown seed geometry, humidity response and transfer across temperature still prevent adoption of these fitted controls for full lunar aggregation.

![Fitted histories and the laboratory comparisons](../../visualization/snow/generated/growth_limits.png)

The figure compares each preceding temperature-specific fit, the independent early/later fit and the refit with original transport and mean ventilation. Separate laboratory-fit segments retain their joins. The plot begins at the three-minute comparison age; points mark sampled model ages and lines connect samples. [PDF figure](../../visualization/snow/generated/growth_limits.pdf).

## Reproduction and verification

Run one numerical worker with BLAS/OpenMP thread counts set to one:

```sh
python -B -m atmosphere.snow.mcsnow_growth_limits cold
python -B -m atmosphere.snow.mcsnow_growth_limits bounds
python -B -m atmosphere.snow.mcsnow_growth_limits warm_compact
python -B -m atmosphere.snow.mcsnow_growth_limits constraints
python -B -m atmosphere.snow.mcsnow_growth_limits nested
python -B -m atmosphere.snow.mcsnow_growth_limits analyze
python -B visualization/snow/growth_limits.py
```

The first four default searches use 35 differential-evolution generations; two bound profiles and six constrained searches use 25. Population size is six times the number of free coordinates, followed by bounded least-squares refinement with at most 240 function evaluations. All twelve global runs reach their declared generation limits; every local refinement reports a stopping criterion. Named nested followups use local refinement only. These finite searches do not establish global convergence. The objective gives equal weight to mass and both axes at ages 5, 10 and 20 minutes for both the spherical start and axes-preserving restart. Speed, other ages and the volume-preserving mapping are excluded from the score. Candidates must still execute through the requested thirty-minute horizon; failures receive a penalty. Other temperatures are used for exact neutral-limit checks, not for a new temperature-transfer fit.

The [numerical product](results/mcsnow_growth_limits.json) records 13 search/refinement results, 16,085 parameter-batch evaluations and 470 failed native candidates, retaining their raw journals and partial output. All 40 numerical checks pass. Four neutral/parent comparisons reproduce 114 sampled states each exactly. Selected histories and supplied-state rates are refined separately; the largest reported relative difference is 0.1437%. Maximum mass-ledger and spheroid-volume residuals are 2.22e-16 and 2.22e-15. The [check record](results/mcsnow_growth_limits_checks.json) records exact cached reassembly, 35 passing focused tests, and the passing full `make check`: 3,756 Python tests passed with 20 skipped, and all 191 JavaScript tests passed. The layer check reports 1,329 files and zero violations; provenance and ensemble integrity pass. Both PNG and PDF figures were inspected, with 54 plotted curves checked against their source values. A copied, hash-verified cache checkpoint preserves these journals and check logs before further searches append to the working cache. Simulation outputs remain on the project-space drive. Numerical consistency and improved fitted agreement are separate from empirical validation.
