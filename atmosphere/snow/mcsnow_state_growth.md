# Crystal growth through the early shape transition

An early compact stage followed by adjusted axial growth substantially reduces the warm-column discrepancy across both starting histories. A common adjustment across temperatures retains substantial mass and shape errors. The best identified state-dependent prescription reduces the squared-log score by 22.4% relative to the preceding best constant adjustment, from 0.134464 to 0.104290. At −5.3 °C, combining the later adjustment with a 40 µm compact-stage limit reduces the warm squared-log loss by 94.6% relative to its preceding constant-temperature fit, while keeping both starts in the objective. These diagnostic fits identify a useful early/later distinction; they establish neither a physical attachment law nor validated lunar crystal or aggregate predictions.

The geometric requirement is sharper than the optimization result. With the assumed 4 µm spherical seed, the fitted three-minute column requires an equivalent early logarithmic axis-growth ratio near 1.12; its first measured-age segment instead implies 5.42. For the dendrite, the corresponding ratios are 0.145–0.158 before three minutes and 0.547 after that age, with the range arising from the two envelope mappings. A single developed-crystal ratio cannot produce both stages from that spherical seed, regardless of how quickly total volume is added. Unknown laboratory seed geometry remains an alternative explanation.

## What the literature adds

The [2021 lateral-growth study](https://doi.org/10.1175/JAS-D-20-0228.1) treats facets spreading across a seed. Its [author archive](https://www.datacommons.psu.edu/commonswizard/MetadataDisplay.aspx?Dataset=6193) includes digitized measurements and model output. In four complete measured rows at −15 °C and 1.5% ice supersaturation, diameter increases from 23.81 to 26.51 µm while one pyramidal region's area falls by 95.8%. Thus a substantial surface transition accompanies only 11.4% growth in outer diameter. The separate author model output is retained as model output; it is not another experiment. These conditions differ from the present water-saturated benchmark and its unmeasured 4 µm seed.

[Thin-rim attachment experiments](https://doi.org/10.1016/j.jcrysgro.2013.04.037) and a [humidity-dependent morphology criterion](https://doi.org/10.1175/JAS-D-23-0131.1) motivate carrying surface structure and moisture conditions into a physical model. Radius is only a proxy for those quantities here. No rim width, facet-area fraction or molecular attachment coefficient is inferred from the present fit. The [source record](mcsnow_state_sources.json) distinguishes the retrieved archive, indexed article passages and unresolved full-article retrievals. The complete 2021 publisher HTML and PDF still return 403; the implemented equations below are independent diagnostic assumptions.

## A shape constraint independent of growth speed

The native shape update has

$$d\ln\phi=\frac{\Gamma-1}{\Gamma+2}\,d\ln V,$$

where $\phi=c/a$, $V$ is enclosing volume and $\Gamma=d\ln c/d\ln a$. Integrating this relation between the seed and the fitted three-minute state determines an equivalent constant ratio without specifying the time dependence of uptake. It is not the arithmetic mean of a varying ratio and is not a measured attachment coefficient. Both admissible envelope mappings give the early/later discrepancy above. The warm fitted ratio then decreases to 3.22 after six minutes. Its rise and subsequent decline are an additional constraint on a physical transition model; the independently fitted join does not locate an observed instantaneous switch.

The same calculation exhibits a degeneracy: holding the first measured-age ratio throughout growth would instead require an initially flattened seed. At the same initial envelope volume, the necessary aspect ratios are 0.113–0.117 for the dendrite and 0.0081–0.0092 for the later column. These are mathematical alternatives, not measured seed shapes. Takahashi et al. generated seeds by adiabatic cooling; the fitted histories begin at three minutes and describe different crystals. Consequently, neither an unknown germ nor an omitted surface transition can be excluded by these endpoint fits alone.

## Controlled changes to the native component

The [runner](mcsnow_state_growth.py) reuses McSnow commit `1e428d7e901cbc5b91f2040e575f736ec5963a74` and the previously retained removal of the extra oblate mean-ventilation term. The original figure executable and configuration remain unidentified; these comparisons use a declared custom driver and patches. It builds a separate executable and preserves the original model, previous products and reference fits. The generalized batch interface keeps existing default behavior and includes each alternate experiment in cache identity. Simulation artifacts remain on the project-space drive.

First, ten controls keep an initially spherical seed spherical, with solid deposition, below maximum full diameters of 0, 10, 20, 40 or 100 µm, using either the original controls or the preceding shared fit. Every positive threshold worsens the joint training score. Separately, a 40 µm threshold improves the warm spherical start at three minutes from mass/width/length ratios 1.825/0.826/2.180 to 1.289/1.018/0.906, but leaves length 71.6% too low at thirty minutes. The cold dendrite becomes much too thick. This opposite response motivates the temperature-specific combination below. This rejects that simple delay as a stand-alone repair; it does not test the mass redistribution involved in actual lateral facet spreading.

Next, both diagnostic families retain the six preceding controls and add one transition radius, $R$, between 0 and 150 µm. For the current enclosing-volume equivalent radius $r$, define

$$w=\left[1+(R/r)^4\right]^{-1}.$$

The `axis` family replaces the constant axis multiplier $q$ with $q^w$. The `axis_density` family additionally replaces the already bounded deposition-density fraction $f=\rho_{dep}/\rho_i$ with $f^w$, approaching compact deposition at small size. Both recover the previous six-control model exactly when $R=0$. The fourth-power transition, its radius bound and these interpolation choices are exploratory assumptions. They are not published kinetic coefficients. The original temperature table, native branching rule and shape-dependent directional ventilation remain active, with the six diagnostic controls applied as documented in [the preceding tuner](mcsnow_tuning.md).

Four searches use two random starts per family, one numerical worker, 35 differential-evolution generations with six population members per free parameter, and bounded least-squares refinement with a 200-evaluation limit. They fit the same four trajectories, three ages and three quantities as before: −14.4 and −5.3 °C; spherical and three-minute starts; ages 5, 10 and 20 minutes; mass and both dimensions. Reserved measurements, the volume mapping, other temperatures and speed do not enter the score. A candidate must nevertheless execute through the full requested horizon, including reserved ages; a native failure receives a penalty. Finite searches identify candidate minima, not a proven global optimum.

| Search | Training RMS log error | Transition radius (µm) |
|---|---:|---:|
| axis_20261024 | 0.3229 | 52.20 |
| axis_20261025 | 0.3229 | 52.47 |
| axis_density_20261024 | 0.3275 | 35.48 |
| axis_density_20261025 | 0.3275 | 34.94 |

The selected common family is `axis`. Its parameters are compared with the lowest-scoring warm trial, whose fixed compact-stage limit is 40 µm:

| Control | Common fit | Selected warm fit |
|---|---:|---:|
| Transport multiplier | 0.87742 | 0.9174 |
| Mature axis-growth multiplier | 1.9418 | 2.6372 |
| Directional-ventilation exponent | 0.88898 | 4.3748e-09 |
| Porous deposition-density multiplier | 0.89655 | 0.53248 |
| Branch-radius multiplier | 0.40768 | 1.0118 |
| Mean-ventilation enhancement | 0.58949 | 0.1091 |
| Transition equivalent-volume radius (µm) | 52.203 | 29.566 |

The transition radius above is an equivalent-volume radius. The branching control instead acts on the largest semi-axis and corresponds to 81.5 µm. These are different dimensions and operations. The warm branching-radius value is unidentifiable here: changing its multiplier from 0.25 to 4 leaves the selected warm trajectories unchanged, because native columns hollow immediately and only plates use the branching-size criterion. The other multipliers also remain compensating empirical adjustments; none is a measured correction to a material constant.

## Full histories and supplied-state rates

At thirty minutes, which was excluded from fitting, the common prescription gives these signed errors:

| Temperature and start | Mass | Width | Polar dimension |
|---|---:|---:|---:|
| -14.4 °C, 4 µm sphere | -31.0% | -14.7% | +27.3% |
| -14.4 °C, three-minute restart | +10.3% | +41.5% | -26.0% |
| -5.3 °C, 4 µm sphere | +81.3% | -14.2% | +89.6% |
| -5.3 °C, three-minute restart | -31.3% | +7.2% | -53.9% |

Median absolute percentage errors compare the preceding constant fit with the selected state-dependent fit. Different groups contain different state counts; these are descriptive summaries, not uncertainty intervals.

| Comparison | States | Mass | Width | Polar dimension | Speed |
|---|---:|---:|---:|---:|---:|
| training | 12 | 17.1% → 24.3% | 12.3% → 13.9% | 29.1% → 35.3% | 19.5% → 25.4% |
| withheld_age | 10 | 33.0% → 27.5% | 18.3% → 15.3% | 35.4% → 53.4% | 27.4% → 27.1% |
| geometry_sensitivity | 10 | 34.1% → 25.8% | 30.3% → 14.9% | 40.2% → 37.6% | 26.8% → 28.9% |
| other_temperature | 68 | 13.0% → 12.3% | 16.2% → 11.3% | 7.7% → 13.7% | 23.1% → 30.1% |

The common fit reduces the squared-log score without improving every part of the comparison. All four training medians above become worse, while the reduction in larger deviations improves the objective. Across the eight reserved temperatures, median polar-dimension error rises from 7.7% to 13.7% and speed error from 23.1% to 30.1%. Removing the radius transition while holding either family’s fitted controls fixed causes native execution failure; the preceding constant model remains viable after its own optimization.

Supplying the thirty-minute fitted geometry with the axes-preserving mapping directly tests the instantaneous rates separately from the simulated path. Width and polar entries below compare logarithmic rates. Density entries give model / fitted implied deposition density in kg m⁻³.

| Temperature | Mass rate | Width rate | Polar rate | Deposition density |
|---|---:|---:|---:|---:|
| -14.4 °C | -26.4% | +3.5% | -54.3% | 229.9 / 286.6 |
| -5.3 °C | +45.6% | -9.9% | +37.3% | 544.2 / 445.6 |

The common fit also changes uptake to compensate for shape. At the supplied cold thirty-minute state, the original improved component was only 3.7% high in mass uptake; the common state-dependent fit is 26.4% low. The trajectory-score improvement therefore does not establish a more accurate transport coefficient.

Two further searches allow temperature-dependent controls while preserving one set for both starts at each temperature. Two additional warm-column searches combine that adjustment with compact-stage maximum diameters of 20 or 40 µm. These declared trial values test the early/later tradeoff exposed by the first controls; they are not measured transition sizes. They use the selected common family and 25 generations, with the corresponding prior constant-temperature fit available as a fallback. Training RMS log errors are -14.4 °C, compact limit 0 µm: 0.1828 → 0.1828; -5.3 °C, compact limit 0 µm: 0.3544 → 0.3282; -5.3 °C, compact limit 20 µm: 0.3544 → 0.1577; -5.3 °C, compact limit 40 µm: 0.3544 → 0.0821. The cold search returns the preceding constant optimum with a zero transition radius. The smallest warm score comes from the 40 µm trial. Removing that stage with its other controls fixed raises RMS log error to 0.6403; reoptimizing the zero-cutoff model gives 0.3282. These temperature-specific adjustments are diagnostic comparisons at two temperatures; they define no interpolation or lunar transfer law.

| Temperature | Compact limit | Start | Mass at 30 min | Width at 30 min | Polar dimension at 30 min |
|---|---|---|---:|---:|---:|
| -14.4 °C | 0 µm | sphere | -12.6% | -4.4% | -24.9% |
| -14.4 °C | 0 µm | axes | +13.5% | +27.7% | -45.4% |
| -14.4 °C | 0 µm | volume | +10.6% | +14.6% | -46.7% |
| -5.3 °C | 0 µm | sphere | +18.6% | -13.9% | +61.7% |
| -5.3 °C | 0 µm | axes | -37.6% | +6.0% | -44.1% |
| -5.3 °C | 0 µm | volume | -37.6% | +0.1% | -49.3% |
| -5.3 °C | 20 µm | sphere | +30.9% | +0.9% | +56.4% |
| -5.3 °C | 20 µm | axes | -5.9% | +5.6% | +2.3% |
| -5.3 °C | 20 µm | volume | -8.1% | -0.1% | -9.7% |
| -5.3 °C | 40 µm | sphere | +4.1% | +2.6% | +27.9% |
| -5.3 °C | 40 µm | axes | -5.1% | +4.3% | +12.7% |
| -5.3 °C | 40 µm | volume | -8.3% | -0.8% | -2.9% |

The selected warm fit combines its early compact stage with a mature axis multiplier of 2.64, a porous deposition-density multiplier of 0.53, an almost zero directional-ventilation exponent and 0.11 times the mean-ventilation enhancement. Transport changes less, to 0.92 of its retained value. These are the values used by this fit, not independently identified necessities; other compensating parameter sets may exist. The branching-radius setting has the verified zero effect described above.

The reserved results show both progress and a specific remaining failure. At thirty minutes, the selected warm fit keeps mass and width within 5.1% for both trained starting histories, while length remains 27.9% high from the sphere and 12.7% high from the restart. At three minutes, which is also excluded from the objective, the spherical-start mass and length are still 19.8% and 31.8% high. With the axes-preserving mapping at the supplied warm thirty-minute geometry, mass uptake is 7.2% low, width growth 11.6% high and length growth 38.2% high; the implied deposition density is 323 versus the fitted 446 kg m⁻³. For the selected cold fit, uptake is only 1.4% low at that supplied state but thickness growth is 63.8% low, with deposition density 394 versus 287 kg m⁻³. These rate discrepancies direct the next physical test toward growth allocation and evolving density.

The cold result bounds this particular interpolation. At small radius, its axis multiplier approaches one, retaining the native early-growth limit; the fit cannot choose a separate early multiplier. A distinct early/later adjustment is therefore a remaining diagnostic control, motivated by the independent shape budget, before concluding that size dependence cannot help the dendrite. Additional fitted flexibility would still need independent surface-state and humidity evidence.

The product also retains fixed-volume, solid-seed aspect ratios of 0.25, 0.5, 1, 2 and 4 for both temperatures. These vary geometry without varying seed mass. They measure sensitivity to an unobserved input, not a distribution of naturally occurring seeds. Even the selected 40 µm warm fit gives 4.67 times the fitted thirty-minute length from an initially 2:1 seed, and 80.98 times from a 4:1 seed, at the same initial volume and mass. The original warm 4:1 seed, common-adjustment warm 2:1 and 4:1 seeds, and 20 µm warm trial with a 4:1 seed terminate numerically; these failures remain in the product. Thus the spherical-seed improvement leaves strong sensitivity to the initial shape. Withheld ages and restarts still use the same fitted laboratory populations; the two mappings test geometric representation rather than independent validation.

![Sampled histories and the laboratory fits](../../visualization/snow/generated/state_growth.png)

The plot includes the temperature-specific variant selected by its training score. It shows sampled model trajectories alongside the separately fitted laboratory segments; discontinuous fit joins remain separate. [PDF figure](../../visualization/snow/generated/state_growth.pdf).

## Reproduction and verification

Run one numerical worker with BLAS/OpenMP thread counts set to one:

```sh
python -B -m atmosphere.snow.mcsnow_state_growth scan
python -B -m atmosphere.snow.mcsnow_state_growth search
python -B -m atmosphere.snow.mcsnow_state_growth temperatures
python -B -m atmosphere.snow.mcsnow_state_growth seeds
python -B -m atmosphere.snow.mcsnow_state_growth analyze
python -B visualization/snow/state_growth.py
```

These commands require the earlier pinned McSnow build and the author archive in the recorded cache path. The archive’s download URL and SHA-256 are retained in the source record. Its first four complete observation rows are parsed explicitly; later rows contain blank series columns and are excluded instead of being shifted by whitespace parsing. The archive uses prepublication figure numbers: `fig3.output` corresponds to published Fig. 4, whose caption identifies the normal a-axis diameter. `fig5b.output` corresponds to published Fig. 6b; its metadata calls the dimensions semi-lengths, so the extraction uses only dimensionless ratios from that model-output file. Author files remain unchanged in ignored storage. The native build records compiler, source, table and patch identities. Every evaluated vector and unsuccessful execution remains in the raw journals.

The [numerical product](results/mcsnow_state_growth.json) records 11,536 parameter-batch evaluations, including 1,124 failed native candidates. All 16 numerical checks pass. Neutral controls reproduce the preceding original and constant-adjusted histories; selected trajectories and supplied-state rates are refined separately. The largest reported refinement change is 0.1720%. Maximum mass and spheroid-volume residuals are 2.22e-16 and 2.33e-15. The [check record](results/mcsnow_state_growth_checks.json) records the focused and repository-wide checks actually completed. Numerical checks, source reading and physical validation remain separate.
