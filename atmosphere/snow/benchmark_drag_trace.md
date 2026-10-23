# Archived crystal geometry and fall-speed trace

The identified native McSnow drag routine and the authors' archived Python routine agree when given the same crystal and air state. For the selected ten-minute IGF2+ plate at −14.5 °C, both give **0.146190 m/s**, while the numerical archive records **0.2597 m/s**. The [trace product](results/benchmark_drag_trace.json) records the geometry and intermediate drag quantities, together with separate mass and area controls. It locates the unresolved discrepancy in the archived speed's generating state, settings, output convention or source version; it does not identify which of those differs.

This follows the [benchmark reconstruction](benchmark_reproduction.md). The authors' `fig10.py` and `fig13.py` read the speed column but do not use it in their plots. Accordingly, the saved-speed discrepancy is assessed separately from reproduction of the plotted mass, axes, aspect ratio and apparent density. Neither this calculation nor the earlier reconstruction establishes a complete Earth storm or a lunar size distribution.

## One crystal through the native routines

The [diagnostic caller](mcsnow_drag_driver.f90) compiles the unchanged modules from pinned McSnow revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`. It supplies the archived mass, envelope volume and aspect ratio to one pristine monomer, with zero liquid, rime and frozen-rime mass. `m2d_particle` calls the native spheroid and projected-area routines. The caller then evaluates the native Reynolds routine, `boehmVel` and the production `vterm_snow` selector independently, without taking a growth step. Build, source, caller, executable and raw-file hashes are retained.

At the primary archived state, the trace is:

| Quantity | Value |
|---|---:|
| Ice mass | 4.206150 × 10⁻⁹ kg |
| Envelope volume | 1.3383342 × 10⁻¹¹ m³ |
| Polar/equatorial aspect ratio | 0.01470726 |
| Apparent density, mass/volume | 314.2825 kg/m³ |
| Reconstructed full width / thickness | 1.202296 mm / 17.68248 µm |
| Native projected area | 3.992423 × 10⁻⁷ m² |
| Area/enclosing-circle ratio, q | 0.3516607444 |
| Air density / dynamic viscosity | 1.359095964 kg/m³ / 1.643391549 × 10⁻⁵ Pa s |
| Best number / Reynolds number | 686.6428934 / 14.53575227 |
| Speed from each native call path | 0.1461897908 m/s |

The air uses 258.65 K, 101000 Pa and native water saturation. These remain explicit reconstruction choices where the original configuration is missing. Both full dimensions follow from volume and aspect ratio; the independently rounded archived semiaxes provide a separate geometry check. For this pristine oblate state, native projected area is `q π a²`, where `a` is the equatorial semiaxis and `q = (1 − phi) m/(rho_i V) + phi`. Native ice density is 919 kg/m³. The resulting q, diameter and mass enter the rough-particle Böhm relation directly. The three native call paths agree to rounding precision.

The [authors' plotting archive](https://doi.org/10.5281/zenodo.10077724) contains `utility_functions.py`, a second implementation of Böhm drag. The [runner](mcsnow_drag_trace.py) executes that inspected, hash-identified helper verbatim. Its oblate result agrees with the primary native trace to floating-point roundoff. Across all 14 selected states and four diagnostic controls, the largest relative difference is 9.98 × 10⁻⁸. The small residuals involve compact states: the helper uses double-precision decimal coefficients where two native coefficients have default Fortran precision, and its Best-number expression differs when q exceeds one. Independently rounded mass and volume place two selected compact states slightly above solid-ice density; their unmodified q consequently exceeds one by at most 1.73 × 10⁻⁷. Those input bytes are retained without clipping.

This agreement verifies the state-to-speed computation under the supplied geometry and air. The archive does not supply projected area, so it cannot establish that the generating experiment supplied the same drag inputs. The helper's prolate interpolation also differs from the pinned Fortran; this trace is restricted to oblate states.

## What the mass and area controls distinguish

At the same ten-minute temperature, the three variants give:

| Archived variant | Saved speed | Native geometry and mass | Solid-envelope mass and area |
|---|---:|---:|---:|
| HP | 0.183900 | 0.177665 | 0.338493 |
| IGF1 | 0.356400 | 0.183574 | 0.348891 |
| IGF2+ | 0.259700 | 0.146190 | 0.254349 |

All speeds are m/s. The last column preserves dimensions and aspect ratio while replacing mass by `919 V` and projected area by the full enclosing circle. It changes the archived mass by factors of roughly three. Its proximity to the IGF1 and IGF2+ speeds is therefore a diagnostic clue about a possible density/state difference, rather than a repair that preserves the reference crystal. The residuals are still 2.11% and 2.06%, far beyond the speed's printed ±0.00005 m/s interval. HP becomes 84.1% too fast.

The separate interventions identify what produces that proximity. For the primary IGF2+ state, changing area alone lowers speed to 0.121636 m/s; changing mass alone raises it to 0.302563 m/s. Both changes are needed for the 0.254349 m/s result. The neighboring −14.75 and −14.25 °C IGF1 and IGF2+ states retain approximately 2% residuals under the combined change, while the nearly dense three-minute states already differ much less under their archived mass: native speeds are 0.096274 versus 0.0965 m/s for IGF1 and 0.083224 versus 0.0836 m/s for IGF2+.

One source-supported air sensitivity uses the fixed density and viscosity in the authors' **column** plotting script, `fig16+18+20.py`: 1.28780 kg/m³ and 1.715772 × 10⁻⁵ Pa s. With the primary archived mass and native area, that helper gives 0.146012 m/s, leaving the large discrepancy. These column plotting values do not identify the single-crystal experiment's air. No air properties were fitted to the desired speed.

## Consequence for the growth comparison

Removing speed from the acceptance criterion does not recover the preceding growth calculations. Reanalysis of the existing product finds zero of 120 growing temperatures matching the printed precision of mass, both axes and aspect ratio under each candidate/package at either three or ten minutes. This subset is necessary for agreement with the plotted quantities; the separate all-fields-except-speed count also remains zero. The density plot's clipping is not used to excuse a mass or axis discrepancy.

The selected refined ten-minute native IGF2+ calculation at −14.5 °C has 49.10% excess mass and 25.91% excess width; its final timestep changes are 0.0333% and 0.0183%. For IGF1, the corresponding excesses are 20.48% and 11.01%, with last-step changes of 0.0731% and 0.0678%. Thus the growth discrepancy remains substantial independently of the archived speed field. These figures reuse the preceding identified runs; this checkpoint launches no new growth trajectories.

The immediate diagnosis supports keeping the native drag unchanged. It also prevents treating a saved, unplotted speed as proof that the authors' growth integrator used that value. Resolving the archive discrepancy requires the generating driver/configuration or an equivalent record of mass, projected area, characteristic length, air and speed at evaluation. The inspected public artifacts have not supplied that record. Subsequent growth work should trace the earliest mass-and-axis departure under the documented candidate, preserving the distinction between a source reconstruction and a reproduced original experiment. The solid-envelope substitution has not earned adoption in that work.

## Reproduction and checks

From the atmosphere checkout, with the previously identified baseline and archives present:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m atmosphere.snow.mcsnow_drag_trace run
python -m atmosphere.snow.mcsnow_drag_trace assemble
```

The runner compiles serially, performs 56 evaluations of 14 fixed states in one execution and verifies retained raw files before reuse or assembly. All raw files use `research/runs/snow_scale/mcsnow_reproduction/drag_trace/`, following the existing link to the project-space drive. The [check record](results/benchmark_drag_checks.json) separates interface tests, source comparison and repository checks from empirical validation. The report and the selected native and author routines were read in context; complete-paper review and the original experiment's provenance remain outstanding.

`make check` passes: 1,310 source files with zero layer violations, 3,642 Python tests passed and 20 skipped, all 191 JavaScript tests passed, and provenance and ensemble integrity passed. Reassembling the saved result from raw files and identified author archives gives exact equality. These checks establish software and record consistency; the scientific limitations above remain.
