# Remaining crystal-shape error and its effect on collection

The improved McSnow component still fails the joint mass-and-axis comparison with the Earth laboratory fits. At the supplied thirty-minute −14.4 °C laboratory-fit state, its instantaneous mass uptake is **3.7% high**, while its logarithmic width growth is **38.3% high** and thickness growth **71.8% low**. Matching mass along the model trajectory leaves substantial dimension errors. The next growth repair therefore has to address shape evolution as well as uptake. The fixed-property collection calculations execute consistently, but their results depend materially on the crystal properties supplied to them. No full lunar aggregate distribution is adopted.

This checkpoint contains **48 growth trajectories, 148 fixed-state probes, 36 directional-growth controls and 1,098 native pair evaluations**. The [growth diagnosis](results/mcsnow_shape_diagnosis.json), [directional controls](results/mcsnow_shape_controls.json) and [collection controls](results/mcsnow_pair_controls.json) retain the inputs, source identities, exclusions and numerical comparisons. The directional controls deliberately include ratios inferred from the comparison fits; their agreement is a causal diagnostic, not independent validation.

## Laboratory inputs and the model representation

The calculation uses the existing identified transcription of [Takahashi et al. (1991), Tables 2 and 4](https://doi.org/10.2151/jmsj1965.69.1_15): mass, two full dimensions and fall speed, at ten temperatures from −3.7 to −22 °C. Each comparison stays within its fit's measured age interval, beginning at three minutes. The fits summarize different crystals, have no pointwise uncertainty, and can jump where independently fitted segments meet. Their derivatives are evaluated only within individual segments; the product records every join. Thin planar thickness was inferred from mass and basal area, so mass and both dimensions are not three independent measurements.

Two mappings expose the representation uncertainty. The first preserves measured mass and both axes in a spheroid. The second preserves mass, the enclosing hexagonal-prism volume and aspect ratio, multiplying both model axes by `(9 sqrt(3)/(4 pi))^(1/3)`; only the reported axis comparison reverses this mapping. Native drag and collection use the actual model spheroid. Neither mapping supplies a measured projected-area history or establishes the adequacy of a spheroid for a branched crystal.

States whose mapped density exceeds native bulk ice density remain explicit exclusions. This removes twelve trajectory configurations and 52 fixed-state configurations, including repeated timestep or model controls; these counts are not numbers of rejected observations. Preserving the hexagonal volume makes four previously inadmissible three-minute states usable. The −16.5 °C fitted three-minute state remains inadmissible under both mappings. No mass or dimension is clipped to make a restart possible.

## Uptake and shape can be tested separately

All principal trajectories use the retained improved component: the pinned [Welss et al. (2024) McSnow](https://doi.org/10.1029/2023MS003805) revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`, with only the extra oblate mean-ventilation term removed. Native intrinsic growth, deposition density, directional ventilation, geometry and drag remain active. The atmosphere is held at each laboratory temperature, 101000 Pa and exact native water saturation. A 4 µm spherical seed supplies the earlier initialization control; admissible fitted three-minute states supply two separate restarts. These are documented custom-caller calculations, not a recovered author executable.

For each fitted smooth segment, write mass proportional to age raised to `b`, equatorial span to `a`, and polar length to `c`. The implied effective logarithmic axis-growth ratio is `c/a`; the implied deposition density is `(mass/envelope volume) b/(2a+c)`. These quantities characterize the fitted history under the selected envelope convention. They are not separately measured facet coefficients or new transferable growth laws. Some fitted segment rates imply deposition densities above bulk ice even when the instantaneous state is admissible; those inconsistencies remain flagged rather than becoming calibration targets.

The native routine is also evaluated at each supplied fitted state, before it can accumulate a trajectory error. At thirty minutes and −14.4 °C, preserving both measured axes gives an effective growth ratio of **0.10385**, against **0.50980** implied by the fitted history. Native deposition density is **256.43 kg/m³**, against an implied **286.58 kg/m³**. Bringing the width rate into agreement would require multiplying uptake by **0.723**, whereas bringing thickness growth into agreement would require **3.550**. A single uptake multiplier cannot satisfy both. With the volume-preserving mapping, the corresponding multipliers are **0.827** and **4.096**; the diagnosis survives that representation change.

The equal-mass comparison asks a different question: when does the computed crystal reach the mass fitted at a selected laboratory age, and what dimensions does it then have? It interpolates only within a computed, growing trajectory and records the bracketing times. For the −14.4 °C spherical-start case, matching the thirty-minute fitted mass leaves width **20.9% high** and thickness **18.4% low**. For the −5.3 °C spherical-start column, polar length remains **31.1% low** at matched mass. These comparisons shift model age; they diagnose the limits of a rate-only change rather than supplying successful age-matched predictions.

## Directional controls improve particular cases, without selecting a repair

Two additional controls retain the native uptake calculation, intrinsic growth table, deposition density and branching rule. The first sets directional ventilation to unity. The second supplies the effective axis-growth ratio derived from the fitted history, beginning only at three minutes and respecting the fitted phase change. The spherical case retains its native early growth before that age. This second control intentionally uses the comparison data and has no physical prescription for changing temperature or humidity.

The following values are final model-to-fit ratios at thirty minutes, using the finer timestep:

| Temperature and initialization | Directional treatment | Mass | Equatorial span | Polar length |
|---|---|---:|---:|---:|
| −14.4 °C, spherical seed | Retained improved model | 1.260 | 1.352 | 0.825 |
| −14.4 °C, spherical seed | Unit directional ventilation | 1.077 | 1.081 | 1.102 |
| −14.4 °C, three-minute axes restart | Retained improved model | 2.317 | 2.343 | 0.509 |
| −14.4 °C, three-minute axes restart | Unit directional ventilation | 1.978 | 1.857 | 0.692 |
| −14.4 °C, three-minute axes restart | Fitted effective ratio | 1.443 | 1.236 | 1.137 |
| −14.4 °C, three-minute volume restart | Fitted effective ratio | 1.426 | 1.130 | 1.085 |
| −5.3 °C, three-minute axes restart | Retained improved model | 0.900 | 1.406 | 0.315 |
| −5.3 °C, three-minute axes restart | Fitted effective ratio | 1.850 | 1.063 | 1.132 |

The spherical dendrite's smaller errors under unit directional ventilation do not persist across starting histories. Supplying the fitted ratio improves the restarted dendrite's dimensions and falling speed, while leaving its mass 44% high. For the restarted column, the same diagnostic brings dimensions closer while making mass 85% high. Applying the fitted ratios after three minutes to the spherical-start column produces 15.1 times the fitted mass and 7.05 times the fitted polar length. The phase-specific ratio is therefore rejected as an autonomous replacement. These results expose interactions among initial shape, subsequent shape evolution, deposition density and uptake; selecting a favorable endpoint would conceal those failures.

## Collection with prescribed crystal properties

The [pair driver](mcsnow_pair_driver.f90) calls unchanged native geometry, fall-speed and collection-efficiency routines. Its outer encounter area follows the native `ckt0` expression, including the prolate aspect correction. Each pair contains crystals from two selected ages at a common temperature; these are hypothetical coexisting crystals, not laboratory collision observations. Growth is disabled. The comparisons use Earth gravity, 101000 Pa, exact native water saturation, pure ice and no added turbulent speed spread. Temperature-only sticking option 1 and unit sticking remain separate sensitivities.

Thirty diagnostic-geometry pairs agree with the existing independent Python implementation to a maximum relative difference of **1.99 × 10⁻¹⁴** across the nonzero kernel, efficiency, encounter-area and relative-speed quantities. All 1,098 native evaluations are finite and nonnegative; reversing a pair preserves collection efficiency exactly, identical particles have zero differential-settling kernel, and changing sticking scales the kernel as specified. These establish implementation consistency for the tested conditions. They do not establish observed sticking efficiency or validate aggregate geometry.

Across the declared nonzero comparisons, replacing fitted crystal properties with the improved model's grown properties changes the collection kernel by factors **0.267–5.218**. This deterministic selection is not a population-weighted range or uncertainty interval. For the five- and ten-minute −14.4 °C spherical-start pair, the kernel is **0.267** times the value obtained from the fitted axes and masses. For the twenty- and thirty-minute −5.3 °C pair, it is **1.806** times that reference. Both references still use native drag and projected-area calculations.

Agreement in a total kernel can also conceal compensating errors. In the ten- and thirty-minute dendrite pair from the fitted-axes restart, encounter area is **5.00** times its fitted-state reference, relative falling speed **0.559** times and collection efficiency **0.458** times; their product gives a kernel ratio of **1.281**. An aggregation result alone would not reveal the incorrect crystal properties behind that closer number.

## Verification and next work

The growth trajectories use steps of 1/32 and 1/128 seconds, saving every quarter second. The largest sampled relative change in the growth diagnosis is **0.0434%**; the directional controls have the same maximum at saved comparison ages. Fixed-state probes take their source values from the initial evaluation. All twelve unchanged directional-control executions reproduce their preceding raw trajectory, air, source-trace, branch and step-diagnostic files byte for byte. Source-trace checks separately verify native deposition density and the prescribed effective ratio. All three products reassemble exactly from the retained raw outputs. Across 441,020 saved states, the largest particle-mass ledger residual relative to current mass is 2.22 × 10⁻¹⁶ and the largest spheroid-volume residual is 2.66 × 10⁻¹⁵. The repository-wide `make check` passes: 3,721 Python tests passed, 20 skipped; all 191 JavaScript tests passed; 1,321 source files have no layer violations; provenance and ensemble integrity pass. The [check record](results/mcsnow_shape_checks.json) separates these activities from scientific agreement.

The source and raw runs remain on the project-space drive, through the existing `research/runs/snow_scale/` links. Development builds and exploratory pre-three-minute ratio runs are retained separately and excluded from the final products. Reproduction uses:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_shape_diagnosis run
python -B -m atmosphere.snow.mcsnow_shape_diagnosis assemble
python -B -m atmosphere.snow.mcsnow_pair_controls run
python -B -m atmosphere.snow.mcsnow_pair_controls assemble
python -B -m atmosphere.snow.mcsnow_shape_controls run
python -B -m atmosphere.snow.mcsnow_shape_controls assemble
```

The next physical change must account jointly for effective axis growth and deposition density, and face the fitted restarts as well as the spherical seed. Humidity response, initial shape and independent mass/axis evidence must constrain that change before it is transferred to changing lunar weather. The checked pair interface can support this work with prescribed properties; evolving aggregate geometry, population water and heat feedback, breakup and the matched full Earth–Moon comparison remain open.

The subsequent [automatic tuning study](mcsnow_tuning.md) searches the coupled uptake, directional-growth, density and branching controls. Its separate-history successes and failed transfers quantify the remaining initialization dependence; fixed-control refits show why a fitted coefficient alone does not identify the missing physics.
