# Joint snow accuracy and non-regression requirement

**Replaced on 24 October 2026.** The author replaced this rule with a rule sized to the laboratory data. The model is compared with each crystal Takahashi and co-authors measured, at the crystal's own age and temperature, and passes on mass and on each measured axis when it lies within the scatter between repeat crystals; its fall speed is compared with the paper's fall-speed relation within the stated ±3 cm/s measurement accuracy ([snow README](README.md#the-data-sized-rule)). The 1,320 comparisons below set fitted population histories, growth rates differentiated from those fits, inferred deposition densities and unmeasured seed shapes against a 10% limit, while the fits carry no pointwise uncertainty and repeat crystals of one age differ by up to 2.2 times in mass at −5.3 °C and 2.4 times at −14.4 °C ([measurements](results/takahashi_measurements.json), crystals within 0.5 °C of each temperature). The audit and the candidates it tested remain below as the record of that rule.

## The audit as of 24 October 2026

No tested candidate meets the requested limit: every tracked physical error must be at most 10% and must not increase from the current model. The [acceptance product](results/mcsnow_joint_acceptance.json) freezes the preferred `joint_join_balanced` result in commit `ded734e`, product SHA-256 `304b90e52f4c0269cc0959d8157c91c9d0c2aee2def1d8f950bb0a9bce4db251`. The preceding [column-density experiment](mcsnow_column_density.md) remains a diagnostic calibration with known failures. Its mass-and-dimension result does not establish overall accuracy.

The audit covers 1,320 scalar comparisons at −5.3 and −14.4 °C. These comprise every recorded trajectory age and all four quantities (mass, equatorial dimension, polar dimension and fall speed), the alternate envelope mapping, five starting aspect ratios, and supplied-state growth rates, growth allocation and deposition density. Each candidate must preserve complete coverage and satisfy

\[
 |y_i/y_{i,\mathrm{reference}}-1|
 \leq \min(0.10, |y_{i,\mathrm{baseline}}/y_{i,\mathrm{reference}}-1|).
\]

An average score cannot override a failed comparison. The comparison allows only 32 double-precision machine epsilons of roundoff in relative error. Mass and volume checks retain their separate 10⁻¹⁰ numerical tolerances. Missing states, altered targets, failed native executions and nonfinite values prevent acceptance. The code recomputes errors from values and targets instead of trusting stored scores.

The targets retain their different evidence states. Supplied restart mass and axes are initial conditions; their agreement checks input preservation. Trajectories are compared with published fitted population histories that have already informed calibration. The different tiny seeds were not individually observed; their comparisons are robustness requirements against that same history. Local rates and deposition density are derived from the fitted histories and an envelope convention, rather than direct measurements of individual facet growth or porosity. The 10% requirement is a calibration criterion, not a reported observational uncertainty. This audit covers the two currently calibrated temperatures; it supplies no other-temperature or lunar validation.

## Baseline failures and structural constraints

The baseline exceeds 10% on 426 comparisons. Across the recorded trajectory ages and all three starts, its maximum absolute errors are:

| Temperature | Mass | Equatorial dimension | Polar dimension | Fall speed |
|---|---:|---:|---:|---:|
| −5.3 °C | 9.04% | 5.64% | 24.21% | 40.39% |
| −14.4 °C | 17.02% | 12.89% | 7.75% | 13.21% |

The warm 40.39% speed error occurs at three minutes from the spherical seed. The previously reported 21–23% concerns thirty minutes. The warm 24.21% length error belongs to the alternate volume mapping; the two principal starts retain the preceding 9.73% mass-and-dimension maximum. Starting-shape sensitivities and supplied-state rates introduce additional failures, including the previously recorded elongated-seed runaway.

Two necessary conditions explain why further tuning of the preceding density law cannot satisfy all these diagnostics. Its constant cold porous-deposition density faces targets from 231.02 to 395.58 kg/m³; even the best constant has a 26.26% maximum relative error. The warm axes-mapping targets require deposition density to decline from 334.49 to 247.10 kg/m³ between 3.5 and 5.5 minutes, then rise to 380.28 kg/m³ at eight minutes. Their 10% intervals exclude both a monotonically increasing and a monotonically decreasing density law. All three states exceed the permitted 80 µm compact-stage threshold. These are constraints on reproducing the derived diagnostic targets, with their fitted-history and geometry assumptions preserved.

## Drag, shape and density experiments

The [fixed-state drag product](results/mcsnow_joint_drag.json) compares ten drag choices for 354 supplied or predicted states. The independent Python translation matches the compiled native calculation within 1.12 × 10⁻¹⁵ relative difference. For columns, the native calculation blends a volume-equivalent cylinder and a prolate spheroid, using different characteristic widths. It also uses a different Best-number denominator from the expression identified as original in the pinned source comment. These are documented model choices; the original Böhm paper was not newly transcribed and verified here.

The promising drag sensitivity combines that original-documented denominator with the side projection of a horizontal hexagonal prism in the prolate branch. Averaging its projected width over azimuth gives 3/π times the corner-to-corner span, bounded by √3/2 and 1. The published literature motivates examining shape, orientation and characteristic length, while also distinguishing viscous-limit calculations from finite-Reynolds-number drag. The present geometry sensitivity remains unvalidated for evolving hollow columns. [Westbrook, small-crystal fall speeds](https://arxiv.org/abs/0801.1618); [McCorquodale and Westbrook, TRAIL part 2](https://doi.org/10.1002/qj.3936).

Reintegrating growth with this drag choice reduces the two warm thirty-minute speed errors from 21.37% and 22.95% to 5.59% and 6.89%. It changes ventilation as well: fitted-restart mass error grows from 3.43% to 5.37%, and length error grows from 9.33% to 10.28%. The comparison therefore fails the requested rule. Its maximum warm trajectory speed error is still 15.88%.

A separate compact-growth experiment replaces equal relative-axis growth by equal normal face advances. It preserves the spherical history and reduces the 4:1 seed's thirty-minute length ratio from about 355,799 to 1.428. Seven previously better comparisons nevertheless worsen, so it also fails acceptance. This result identifies persistence of the initial aspect ratio during the compact stage as a source of the runaway without establishing a complete early-shape law.

The third experiment lets added-ice density scale with the crystal's current bulk density, with separate early and late ratios, and adds a size-dependent transport control. Four local rate fits and three joint trajectory fits test these changes. The sharper transition brings maximum local diagnostic errors to 11.10% warm and 10.50% cold, but some previously smaller errors rise. Carrying those settings through complete histories exposes larger trajectory errors. Joint refits reduce their aggregate violation score but still fail individual constraints. The preserved-density refit also tests a diagnostic 0.85 interpolation between the two Best-number denominators. The warm branching-radius control is inactive in these column cases and cannot be inferred from their fits.

| Full comparison | Values above 10% | Values worse than baseline | Accepted |
|---|---:|---:|---|
| Neutral reproduction | 426 | 0 | No |
| Prism drag with existing growth controls | 288 | 409 | No |
| Equal normal growth during compact stage | 388 | 7 | No |
| State-dependent density, settings fitted to local rates | 685 | 818 | No |
| Joint warm fit retaining preceding density law | 205 | 615 | No |
| Joint cold fit with state-dependent density | 446 | 132 | No |
| Joint warm fit with state-dependent density | 501 | 596 | No |

Each row covers the complete contract; the single-temperature refits preserve the other temperature's baseline. The [joint experiment product](results/mcsnow_joint_growth.json) records 3,548 native parameter-batch evaluations, including four local fits, three trajectory fits, controls and refinement. Every native evaluation completed. The product selects no replacement model. These bounded local searches do not prove that every possible new growth law is infeasible.

The remaining task is a growth law that predicts shape and density consistently from both tiny seeds and supplied measured-age states while satisfying the individual limits. Drag and growth must be evaluated together because settling changes ventilation. A successful local-rate fit, a lower aggregate score or a better late-time speed cannot establish that result. Independent geometry and growth evidence remain necessary to distinguish a physical improvement from calibration to the fitted curves.

## Reproduction and verification

The new modules leave the earlier products and native builds intact:

```sh
python -m atmosphere.snow.mcsnow_joint_acceptance assemble
python -m atmosphere.snow.mcsnow_joint_drag run
python -m atmosphere.snow.mcsnow_joint_growth probe_fit
python -m atmosphere.snow.mcsnow_joint_growth screen
python -m atmosphere.snow.mcsnow_joint_search
python -m atmosphere.snow.mcsnow_joint_growth refine
python -m atmosphere.snow.mcsnow_joint_growth assemble
```

The cache is under `research/runs/snow_scale/mcsnow_reproduction/growth_diagnosis/joint_acceptance`, through the existing link to the project-space drive. Native source patches, compiler identities, executables, input configurations, every parameter evaluation and rejected fit are retained there. The first, broader shape-transition trial has its own preserved producer and build. Source retrieval records retain both downloaded PDF hashes and the failed Welss publisher retrieval. Reading covered the specified sections of the two drag papers; no full scholarly review is claimed.

Neutral reproduction agrees with all 1,320 baseline values within 1.12 × 10⁻¹⁵. The prism-drag calculation agrees with its independent fixed-state translation at both supplied three-minute states. Nine timestep refinements from 1/128 to 1/256 s change the tracked values by at most 0.191%, below the unchanged 1% refinement threshold. All twelve numerical checks pass, as do 44 focused tests. These numerical results are separate from the failed physical acceptance criterion. The full `make check` passed: 3,801 Python tests passed, 20 were skipped, and all 191 JavaScript tests passed; provenance and ensemble integrity checks passed. After the final provenance-label, control-name and refinement-coverage additions, the 44 focused tests passed again. The [checks record](results/mcsnow_joint_checks.json) retains the execution record, staged-file layer check and verified raw checkpoint.
