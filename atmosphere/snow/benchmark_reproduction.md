# Earth benchmark reconstruction and duration controls

The archived ten-minute McSnow results remain unreproduced with the identified public source and custom caller. The new [709-execution product](results/benchmark_reproduction.json) separates that source/configuration problem from saturation, timestep and late-growth effects. At exact native water saturation, the refined revised −14.4 °C thirty-minute width is 113.52 mm, still far above the fitted 4.17 mm. The earlier Murphy–Koop-saturation comparison gave 195.39 mm. Disabling the two mean-ventilation operations identified by a source comment gives 3.84 mm, but only 56.6% of fitted mass. Neither candidate is adopted as a repair.

This is the first execution checkpoint of the [handoff plan](benchmark_reproduction_plan.md), dated 23 October 2026. It preserves the earlier native builds and results. Its numerical archive comparison, laboratory-fit comparison and software checks establish different kinds of evidence; none establishes a complete Earth storm or an absolute lunar size distribution.

## Recovered experiment and missing provenance

The [source record](benchmark_reproduction_sources.json) identifies the paper, plotting scripts, numerical archive, public source revision, build settings and remaining gaps. The full 69,250,479-byte [Welss–Siewert–Seifert data archive](https://doi.org/10.5281/zenodo.10077724) matches its published MD5. External bytes stay in the ignored project-drive cache. The derived reference values retain the archive's CC-BY-4.0 credit.

The selected targets are the habit-ventilation curve in Figure 10 (IGF1) and the revised IGF2+ curve in Figure 13. Their filenames encode 1010 hPa, a one-second step, a 2 µm seed radius and 600 seconds of growth. Each file contains 121 temperatures from 273.15 to 243.15 K in 0.25 K steps; corresponding 180-second outputs provide an earlier comparison. The paper specifies water saturation. These are author numerical calculations, separate from the laboratory observations. [Welss et al. (2024), Sections 2.1, 4.1 and 4.2](https://doi.org/10.1029/2023MS003805).

The files have **nine columns**: elapsed time, temperature, two semiaxes, aspect ratio, mass, volume, speed and maximum dimension. The plotting scripts supply eight column names, leaving elapsed time as an implicit dataframe index. The new parser reads every field explicitly. It checks maximum dimension, volume and aspect ratio against semiaxis geometry within the printed rounding intervals, then doubles both semiaxes for comparison with the project's full-axis convention. Raw mass/volume density and the scripts' density clipped at 919 kg/m³ remain separate values. There is no averaging of individual-crystal endpoints.

The inspected public tree at revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`, its deposition-file history, the numerical archive and the plotting archive did not supply the generating single-crystal driver, complete namelist or executable. The public diagnostic source also did not identify that experiment. Consequently, filename settings and paper conventions cannot certify the remaining caller choices. Native moist-air density, the fall-speed selector, prolate interpolation, initial density and the mapping from paper variants to source switches remain explicit reconstruction choices. The HP output is used for fixed-state diagnostics; this checkpoint does not reconstruct every Figure 10/13 curve.

## Candidate builds and comparison criteria

The [runner](mcsnow_reproduction.py) reuses the pinned 2024 McSnow sources, tables, native dependency selector and compiler settings. Its custom caller adds a switch that sets vapour pressure to McSnow's own water-saturation value. With that switch disabled and the frozen initial conditions supplied, trajectory, air, branching and step-diagnostic files reproduce the preceding native baseline **byte for byte**. This checks the new interface against that baseline; it does not identify the authors' original caller.

Two main candidates keep the same directional growth treatment:

| Candidate | Mean transfer treatment |
|---|---|
| Native | Unmodified pinned deposition routine |
| Historical sensitivity | Only the two diffusivity/conductivity multiplications named by the source comment are disabled |

The comment does not identify a particular figure or configuration. Its suggested deletion is therefore tested as a historical candidate, not accepted as the published setting. Twelve additional cases replace only the thermal factor's `(L/(Rv T) − 1)` dependence with its leading `L/(Rv T)` term, motivated by the paper's Equation 4. Native unequal vapour/heat ventilation remains in the enabled-ventilation candidate, so this intervention is not a complete implementation of the printed equation. Every source/caller change, compiled source, table, command and executable is identified in a separate build record.

An archive match requires mass, both full axes, volume, aspect ratio, fall speed and maximum dimension each to lie inside its printed half-rounding interval. Timestep effects are reported separately: selected ten-minute cases use steps from 1 to 1/64 seconds, and the principal revised long cases extend to 1/512 seconds. The last-pair difference measures sensitivity; it is not a rigorous integration-error bound or proof of asymptotic convergence. Failed executions and partial trajectories are retained by the runner, although all 709 executions here completed.

## The short experiment still differs

For each candidate/package pair, **zero of 120 growing temperatures** matches all seven archived quantities at either three or ten minutes. The no-growth 273.15 K endpoint is excluded from these counts. At the archived one-second step, ten-minute median absolute relative errors are:

| Growth package | Candidate | Mass | Equatorial axis | Polar axis | Fall speed |
|---|---|---:|---:|---:|---:|
| IGF1 | Native | 8.36% | 2.64% | 1.76% | 4.03% |
| IGF1 | Historical sensitivity | 18.11% | 6.47% | 4.09% | 11.88% |
| IGF2+ | Native | 7.53% | 2.18% | 1.44% | 1.15% |
| IGF2+ | Historical sensitivity | 15.92% | 5.74% | 5.02% | 9.56% |

The saturation probes do not restore the archived values. In all twelve selected refinement sequences, every final residual exceeds the sum of its printed half-rounding interval and the last timestep change. The largest last-pair change over the seven quantities is 0.124%. For example, the native revised −14.5 °C mass remains 49.1% above the archive at 1/64 seconds. The thermal-factor probes also leave substantial discrepancies. These controls do not warrant treating the public source and custom caller as an exact published experiment.

To separate accumulated growth error from state-to-speed conversion, another control initializes the caller with archived mass, volume and aspect ratio and examines the **initial** fall speed, before integration. Full axes are reconstructed from volume and aspect ratio because all saved fields were rounded independently.

| Archived state at −14.5 °C | Archived speed, m/s | Reconstructed native speed, m/s |
|---|---:|---:|
| HP | 0.1839 | 0.177665 |
| IGF1 | 0.3564 | 0.183574 |
| IGF2+ | 0.2597 | 0.146190 |

Changing prolate interpolation has no effect on these plates. Sampling all eight corners of the mass/volume/aspect-ratio rounding box changes the IGF1 and IGF2+ speeds by less than 0.000000053 m/s across each box, far below the discrepancies. This sampled sensitivity is not a continuous uncertainty proof. At the selected −5.25 °C column state, the corresponding IGF1 and IGF2+ speeds with interpolation differ from the archive by less than 1%; the HP state differs by about 31%. The discrepancy therefore depends on state and variant. The archived output lacks projected area and the original configuration, so the result localizes an unresolved geometry/drag/caller question without identifying its cause.

## Longer exposure distinguishes an apparent size repair

Since no candidate matched the archive, these are **extensions of documented candidates**, rather than extensions of a reproduced paper configuration. Both growth packages run for thirty minutes at all ten available laboratory-fit temperatures, using two timesteps throughout and two further refinements for the revised −14.4 and −5.3 °C cases. The product compares mass, both axes and speed only within each fit's measured age interval; some intervals end before thirty minutes. Density is retained as model output, without turning separately fitted mass and dimensions into an independent density observation.

At native water saturation, the revised −14.4 °C width grows from 1.493 mm at ten minutes to 7.412 mm at twenty and 113.515 mm at thirty, against fitted widths of 1.361, 2.761 and 4.175 mm. The final mass is 43.168 mg, 592.5 times the fit. The last refinement, from 1/128 to 1/512 seconds, changes final width by 0.0434% and mass by 0.0877%. The late discrepancy persists under this numerical control.

The finest thirty-minute revised-package comparisons show why width alone cannot select the historical deletion:

| Temperature and candidate | Maximum dimension, mm | Dimension/fit | Mass/fit | Other axis/fit |
|---|---:|---:|---:|---:|
| −14.4 °C, native | 113.515 | 27.192 | 592.453 | 0.967 |
| −14.4 °C, historical sensitivity | 3.837 | 0.919 | 0.566 | 0.793 |
| −5.3 °C, native | 1.838 | 0.987 | 1.877 | 1.145 |
| −5.3 °C, historical sensitivity | 1.033 | 0.555 | 0.674 | 0.915 |

The other axis is plate thickness or column width. These are comparisons with [Takahashi et al. (1991), Tables 2 and 4](https://doi.org/10.2151/jmsj1965.69.1_15), through the existing identified transcription. The fits summarize different crystals and lack pointwise uncertainty; planar thickness was inferred from mass and area. The seed history is a model assumption. The earlier 240 source-term interventions retain their [separate interpretation](earth_growth_assessment.md#what-the-ventilation-controls-establish); disabling all mean ventilation here is a different intervention from deleting only its extra plate term.

## Consequence and reproducibility

This checkpoint rejects saturation mismatch, timestep refinement and the two-line historical switch as complete explanations of the archived comparison. It also retains the late plate-growth problem after correcting the saturation convention. It does **not** establish that McSnow cannot represent an Earth storm. Finite-water population growth, imposed cooling, latent heating, settling and collisions require their own coupled accounts and comparisons.

The next discriminating source test is the archived-state geometry/drag mapping, retaining mass, volume, aspect ratio, projected area and speed together. An original driver/configuration, or an identified source change that explains the short experiment, would make a later reproduction claim assessable. Further growth alterations must continue to face both archived short cases and joint laboratory mass/axis histories. The existing finite-water column controls remain available, but these candidates have not earned adoption as a lunar closure.

From the atmosphere checkout, with the identified baseline build and downloaded archives present:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m atmosphere.snow.mcsnow_reproduction_benchmarks run
python -m atmosphere.snow.mcsnow_reproduction_benchmarks assemble
```

The matrix runs serially. Existing runs are reused only after input, raw-file and build checks; assembly recomputes the saved analysis from raw CSV files. Generation-time producer hashes are retained in the original records, while the product identifies the current assembler. Raw files and builds use the lexical `research/runs/snow_scale/mcsnow_reproduction/` path, which follows the preserved link to the project drive. The [check record](results/benchmark_reproduction_checks.json) reports the software suite and raw-output audit separately from scientific agreement.
