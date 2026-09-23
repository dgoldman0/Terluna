# A3 — matched molecular atmosphere and frozen-cloud transport

Date: 22 September 2026, with execution records extending into 23 September UTC.
Scope: `immersion/light-cycle/shoreline`. This checkpoint is an additive numerical reference and diagnostic viewer. The production shoreline remains unchanged.

## Delivered result and qualification

A3 solves the shared A1/A2 molecular atmosphere and the fine A2 frozen cloud as one scattering problem. The reference follows molecular and cloud scattering throughout the full spherical domain. Photons may leave the cloud crop, scatter in the surrounding gas, and return. The gas optical state, cloud coefficients, solar spectrum, observer and phase functions are matched across the comparisons.

The calculations exposed a higher-order blue-band convergence gap in the saved A2 molecular reference. A3 retains that failed agreement check and measures it separately. The full coupled solver uses the shared molecular coefficients directly. A second Monte Carlo control evaluates the fixed-clear-molecular-source approximation without A2's moment grid, allowing molecular feedback to be distinguished from that inherited numerical discrepancy.

The implementation checks pass. Qualification of the saved A2 field and the production renderer remains open. Standard errors measure random sampling; they exclude atmospheric, cloud-grid, finite-crop and particle-model uncertainty.

## Inputs and identity

The serialized atmospheric profile is `data/a1/generated/moon-fair.profile.json`. It supplies the 1,737,400-m lunar radius, the 600,000-m optical top, pressure/temperature structure, upper continuation, and selected 300-DU ozone profile. A2's fine frozen field is `data/a2/generated/moon-fair-grid-fine.cloud.json`, a 97 × 65 × 97 vertex grid. The original 36.99-km horizontal crop and approximately 9.289–12.500-km Cartesian vertical extent are retained. Production sampling already includes its spherical-altitude calculation. Extinction outside the cloud box is zero; molecular transport continues to the atmospheric boundary.

This is the selected **fair-weather cloud field**. “Dense” in a ray name describes an optically thick part of that field; it does not designate the separate deep-convection sounding.

The saved A2 moment field is `data/a2/generated/moon-fair-molecular-standard-full.npz`. Profile identity is checked against both the molecular file and frozen cloud provenance. The main record contains hashes of all three inputs, the executed transport kernel, inherited numerical modules, and spectral constants. Profile SHA-256:

```
717deba2cb6f41afc291f1addb13b19c3ebbb6d98cc6dc21354ee02e567d86ca
```

Radiance is evaluated in all 48 supplied bins from 360 through 830 nm, spaced by 10 nm. A1's spectral solar irradiance, ozone cross sections, CIE 1931 functions and XYZ-to-linear-sRGB matrix are preserved with their source notices. No external spectral bytes are added by A3.

## Physical and boundary specification

Molecular density is evaluated from log-linear pressure and linear temperature using the same A1/A2 profile. Molecular scattering is the selected dry-air-equivalent scalar Rayleigh model applied to total number density. Ozone absorbs and has no scattering channel in this reference. The local extinction and scattering source are

```
chi = beta_R(lambda) * n(z)/n_ref + sigma_O3(lambda) * n_O3(z) + beta_cloud(x)
source = beta_R * n/n_ref * P_R[I] + omega_cloud * beta_cloud * P_cloud[I].
```

The cloud is conservative (`omega_cloud = 1`) and uses the exact A1/revision-04 two-component phase distribution: 85% Henyey–Greenstein with `g = 0.78 + 0.04 f_ice`, and 15% with `g = -0.25`. The frozen fair-weather field has zero ice fraction. Its trilinear interpolation and float32 coefficients are shared among the reference treatments.

The spherical ground is black and absorbing. Space is black except for the collimated solar source. Scattered radiance excludes the unscattered solar disk. This preserves the A2 reference boundary; the immersive scene's reflecting ground is a later boundary experiment. The phrase “production boundary” in the transport module's opening comment refers to the benchmark's execution boundary, not the unchanged shoreline's reflective surface.

Aerosols, water-vapor absorption, polarization, refraction, finite-disk illumination, airglow, thermal emission and cloud evolution are outside this experiment. The atmospheric and condensate profiles remain selected inputs, without a climate-equilibrium or precipitation prediction.

## Coupled estimator

`tools/a3/coupled_transport.py` implements backward null-collision path sampling in local production coordinates, with the planet centre at `(0, -R, 0)`. Every free flight is split at spherical majorant-shell crossings and cloud-box intersections. Molecular bounds use `p_max / T_min` over every enclosed profile node and interval endpoint. Ozone bounds include the triangular profile's interior maximum. The maximum frozen vertex extinction bounds trilinear cloud extinction.

Gas and cloud proposals are superposed Poisson processes. The algorithm evaluates the constituent selected by a proposal, reducing repeated profile lookup inside dense cloud without changing the sampling distribution. Accepted events are molecular scattering, cloud scattering or absorption. Null events continue the flight. The scalar Rayleigh angular CDF is inverted analytically; cloud directions use the inherited exact mixture sampler.

At every scattering event, next-event estimation evaluates sunlight through **both** media. Molecular and cloud transmittance are estimated using independent ratio-tracking processes and multiplied for the same physical segment. The product follows `T_gas+cloud = T_gas T_cloud`; spherical planetary occlusion is checked before estimating transmission. A ratio-tracking roulette protects efficiency at very low weights. Path roulette begins after twelve scattering events with survival probability 0.95 and compensating weight. Safety limits and majorant violations are counted as unresolved histories and fail the implementation gate.

The full coupled estimate contains all sampled scattering orders. It receives no added A2 sky background, avoiding double counting. The black-boundary conservative-volume and analytic absorption checks exercise normalization independently of the selected atmosphere.

## Treatments and controlled differences

**Clear molecular Monte Carlo.** The same spherical atmosphere, observer and solar input with cloud extinction set to zero.

**Cloud-only Monte Carlo.** The same frozen cloud, camera, Sun and spherical ground, with molecular scattering and ozone absorption removed. Because this cloud transfer is gray, one unit-irradiance estimate is multiplied by the fixed solar spectrum. Its spectral errors are perfectly correlated; color uncertainty preserves that correlation.

**Full gas + cloud Monte Carlo.** Both constituents participate along every free flight and solar visibility segment. Recorded scattering-history contributions are gas-only, cloud-only and mixed. Every category includes attenuation by both media. These categories are not separately rendered scenes.

**Saved A2 diffuse source + cloud.** Cloud transport is explicit. At the first backward molecular scattering event, the diffuse angular source is evaluated from A2's stored moments and the path terminates. The directly illuminated part of that event still uses gas-plus-cloud shadowing. This freezes the diffuse molecular source while allowing cloud attenuation and repeated cloud scattering.

**Grid-independent fixed-clear-source Monte Carlo.** The same frozen-source approximation is evaluated without the saved A2 moments. At the first molecular event, directly incident sunlight still crosses both media. The subsequent diffuse continuation traces the clear molecular atmosphere with the cloud disabled. Comparing this control with full coupling isolates the response of the diffuse molecular field to the cloud without mixing in A2 moment-grid error.

**Coupled first-order quadrature.** An independent deterministic integrator includes exactly one scattering event, of either constituent, with both media attenuating eye and Sun paths. Cloud line integrals split at every voxel plane; two-point Gauss integration is exact for the cubic trilinear interpolant along each segment. Molecular line integrals use A1's independently expressed shell-split Gauss integration.

**Matched production-expression replay.** The revision-04 diffuse, softened internal-light and compositing formulas are evaluated spectrally with the same cloud and molecular inputs. The constants `.24`, `.34`, `.05`, `.09` and `.16` remain unchanged. Ground reflection is set to zero to match the experiment. Solar and eye attenuation use the shared profile, including ozone, and the empirical diffuse falloff uses the profile-derived local scale height. A2 supplies the clear radiance and diffuse illumination. This is an explicitly adapted equation replay; no claim is made that these values are the unchanged RGB shader, the angular cache or a complete scene render. The A2 field error and physical feedback are reported separately from this comparison.

## Sampling, units and uncertainty

Eight matched observer rays are evaluated: six at a 45° solar elevation at the scene anchor and two at 6°. Observer origins are placed at radial altitude 2 m, retaining small local-angle differences across the spherical surface. Exact vectors are in the benchmark record. Cloud transmittance ranges from a clear gap to optical depth 5.093.

Each of four spectral Monte Carlo treatments uses 32,768 histories per wavelength and ray, split between two independent 16,384-history seed batches. The main calculation therefore contains **50,331,648 histories**. The eight gray cloud-only controls add 131,072 histories per ray, or 1,048,576 histories. All main and gray-control histories completed without unresolved safety limits or majorant violations.

Every wavelength/batch mean, standard error and seed is retained. Ordinary pooled variance includes within- and between-batch contributions. Spectral-to-XYZ uncertainties combine independent wavelength variances. Cloud-only errors instead use their known perfect spectral correlation. Difference error bars combine independent treatment variances. The 1,536 seed-batch comparisons have maximum discrepancy 3.171 estimated standard errors under the declared 6-SE diagnostic threshold. That check provides sampling consistency, not a physical error bound.

Raw spectral units are W m^-2 sr^-1 nm^-1. Y is photopic luminance in cd/m². Color conversion retains signed linear RGB. The display clamps negative channels and uses an exponential curve at a fixed scale of 18,000 with an exposure adjustment, followed by sRGB encoding. Exposure changes only small single-ray color swatches. Spectral curves and all exported evidence retain their original values.

## Coupling measurements

The feedback column compares full coupling to the grid-independent fixed-clear-source control. Error bars below are ±2 Monte Carlo standard errors in cd/m². They exclude model and discretization uncertainty.

| Ray | Cloud optical depth | Coupled Y ±2 SE | Feedback / fixed source | Feedback delta ±2 SE | Mixed-history fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| dense-sun45 | 5.093 | 9951.4 ± 62.3 | +3.00% | +290.2 ± 85.7 | 51.3% |
| moderate-sun45 | 1.729 | 9910.2 ± 49.5 | +3.62% | +346.6 ± 67.9 | 52.2% |
| clear-gap-sun45 | 0.000 | 7449.3 ± 33.5 | +2.60% | +188.5 ± 46.1 | 5.6% |
| edge-sun45 | 0.263 | 8002.1 ± 35.5 | +2.60% | +203.1 ± 48.6 | 23.0% |
| broken-sun45 | 1.307 | 10183.5 ± 46.2 | +2.28% | +226.9 ± 63.6 | 45.3% |
| oblique-sun45 | 0.031 | 9349.1 ± 40.0 | +2.74% | +249.0 ± 54.0 | 9.0% |
| dense-sun6 | 5.093 | 1660.5 ± 15.1 | -2.76% | -47.2 ± 21.3 | 87.2% |
| zenith-sun6 | 0.240 | 2113.7 ± 14.3 | +0.14% | +2.9 ± 20.1 | 23.3% |

At high Sun, the measured diffuse-field response raises the selected rays by 2.28–3.62%. The optically thick low-Sun ray decreases by 2.76%; the other low-Sun ray has no statistically resolved feedback at the present sample count. A single global brightness multiplier would erase these dependencies.

Mixed scattering supplies 51.34% of the dense high-Sun ray and 87.16% of the dense low-Sun ray. Those fractions describe the paths supplying radiance. Their magnitude differs from the incremental change caused by cloud feedback into the molecular source, because the fixed-source control already contains substantial gas/cloud interaction.

The matched production-expression replay is 39.58% below the full coupled dense high-Sun result and 39.28% below the dense low-Sun result. Relative to the Monte Carlo calculation using the **same saved A2 diffuse source**, those discrepancies are 37.29% and 38.66% below. This latter comparison removes the full-coupling/frozen-source difference from the contrast. Other rays differ substantially less, including an edge ray where the replay exceeds the same-A2-source calculation by about 4.50%. These are selected equation-replay comparisons, without a production-image accuracy bound.

First-order view-quadrature refinement changes the worst selected spectrum by 0.2658% in spectral L1. Production-expression replay refinement changes the worst selected spectrum by approximately 1.91e-8 relative spectral L1. These are integration sensitivities for those equations; they do not establish their physical adequacy.

## Independent A2 failure and refinement

Removing the cloud exposed a disagreement with A2's saved all-order molecular field, while independent first-order, pure-cloud, mixed-transmission and analytic tests passed. The clear check was expanded to two independent 1,048,576-history full-Monte-Carlo runs per wavelength, with roulette survival settings 0.95 and 0.99. A separate A2-source termination run reproduces A2 at the tested points, checking the field adapter.

At the 45° Sun / zenith diagnostic, the measurements are:

| Wavelength (nm) | Saved A2 | Independent MC ±2 SE | Dense deterministic refinement |
| --- | ---: | ---: | ---: |
| 400 | 0.120994 | 0.130682 ± 0.000304 | 0.129548 |
| 460 | 0.153182 | 0.159491 ± 0.000308 | 0.158759 |
| 580 | 0.087011 | 0.088035 ± 0.000190 | — |

Relative to the saved A2 values, independent MC is approximately 8.01% higher at 400 nm and 4.12% higher at 460 nm. The broader photopic luminance differences are smaller and ray-dependent; the eight main clear controls differ from A2 by roughly 1.4–3.0% when expressed relative to the MC values.

The same physical atmosphere was then evaluated with substantially denser molecular quadrature and tighter scattering-order convergence. The medium diagnostic uses 40 radial nodes, 24×24 angular quadrature, 128 shells and 384 solar samples; its two radiances are 0.128543 and 0.158247. The dense diagnostic uses 64 radial nodes, 32×32 angular quadrature, 192 shells and 512 solar samples. It reaches 66 orders and a final peak-moment increment of 9.04e-6.

The dense calculation moves closer to independent MC, leaving approximately 0.87% and 0.46% differences at 400 and 460 nm. Those residuals exceed Monte Carlo sampling error, so the dense diagnostic is a further convergence result, not an exact solution. Multiple numerical grids and the order threshold changed together; this experiment does not isolate one faulty interpolation expression. The independent evidence demonstrates that A2's earlier standard/fine comparisons were insufficient to bound its blue-band all-order error.

The failed inherited checks remain in `transport-tests.json` and `A3_VALIDATION.json`. No brightness gain, color correction or rescaling is applied to reconcile them. Existing A2 scientific records remain intact.

## Validation and display checks

Thirty analytic, phase, geometric, limiting-case, sampler and combined-transmission implementation checks pass. Four of six inherited-A2 clear agreement checks pass; the 400-nm and 460-nm checks fail. The aggregate implementation status separates these groups. The original script deliberately returns a nonzero status for the inherited disagreement. A3's validation summary preserves that failure instead of widening its tolerance.

The 159 inherited JavaScript tests pass on the unchanged shoreline runtime. Source/profile identity checks, exact scattering-history closure, Monte Carlo completion and deterministic integration-refinement checks pass. The maximum relative component-closure residual is 3.68e-14. Physical qualification of the saved A2 field and full production renderer remains explicitly false.

The exact generated A3 lab is exercised in Chromium at desktop DPR 1, desktop DPR 2 and 390-pixel mobile DPR 2. Tests cover all eight cases and three spectral comparisons per layout, actual curve and cloud-slice pixels, nonwhite default color swatches, backing resolution, exposure isolation, absence of layout overflow, and byte-identical JSON downloads. The corrected A2 lab is additionally exercised at all eight Sun states in both sky modes at desktop and mobile sizes. At its default +45° multiple-scattering state, the sky has about 30.93 display-code levels of luminance range and zero sampled clipped-white pixels. Together these are 104 state checks.

Local-file navigation is blocked by the browser's environment policy. That failed navigation record is retained. Successful checks insert the exact standalone HTML bytes into an in-memory document without changing policy settings. They establish rendering/control behavior in that environment, not direct local-file navigation or consumer-GPU performance. There are no JavaScript errors, console errors or external requests in the successful run. Final screenshots are inspected separately.

## Source preservation and previous display fix

The source baseline is reconstructed from the delivered A1 source/evidence archive followed by the delivered A2 archive. A3 adds new solver, test, methods and lab files. It also persists the previously delivered A2 display fix into `src/a2-accuracy-lab.js` and `index.a2-accuracy-lab.html`: scale 10,000 and backing-pixel sizing. The source matches the script/template of `Open_Moon_A2_Accuracy_Lab_fixed.html`.

The rebuilt A2 document retains the archive's historical validation metadata, which differs from the later metadata embedded in the standalone fixed HTML. Its numerical sky/fine/cloud/capture payloads and display code match the supplied fix, but the complete HTML hashes consequently differ. The reconciliation record identifies that metadata distinction. A3's independent browser record validates the repaired display; historical A2 validation is not rewritten.

Every production shoreline source, original atmosphere atlas, terrain, material, hydrology and renderer-lab file remains unchanged. The shoreline builder again produces SHA-256 `c7ead1fd83e0d55518d2951e08cd46d31f323280885125382a8013b5d74636d9`. The incremental patch contains ordinary source and compact evidence, with generated fields, full HTML documents, screenshots and per-job caches excluded from Git. The source/evidence archive includes the reusable generated inputs. No GitHub commit or push is performed by this checkpoint.

## Reproduction

Use the Python versions in `requirements-a1.txt`, plus Node.js and the Chromium executable used by the existing tests. Run from the bundled `immersion/light-cycle/shoreline` directory. The bundle includes A1/A2 input profiles, frozen fields, spectral constants and the molecular moment field.

```sh
python build_landscape.py
node --test tests/*.test.cjs > validation/a3/unit-tests.tap

# This intentionally exits nonzero while saved A2's blue agreement checks fail.
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python tools/a3/test_transport.py
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python tools/a3/test_extensions.py

# Main matched 48-band benchmark; individual job records allow same-input reuse.
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  python tools/a3/run_benchmark.py --samples 32768 --workers 3

# Independent clear checks and denser deterministic diagnostics.
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python tools/a3/diagnose_clear.py
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 \
  python tools/a3/refine_molecular.py --quality medium
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 \
  python tools/a3/refine_molecular.py --quality dense

python tools/a3/validate_a3.py
python build_a2_accuracy.py
python build_a3_accuracy.py
python tools/a3/check_lab.py
```

The final reference retains seeds and numerical tolerances. A different runtime may produce different random realizations; interpret agreement using reported sampling uncertainty. The diagnostic record includes an initial unoptimized kernel, its failed inherited checks and the interrupted high-sample development run. These are historical evidence, separate from the completed final benchmark and current-kernel test records.

## Remaining acceptance work

The next molecular work is independent convergence across the full spectrum, especially the blue end, before regenerating a production-ready clear field. The matched equation replay already reveals much larger cloud-lighting discrepancies in thick portions of the selected field. A replacement needs the same paired controls, reflecting-ground tests and a measured GPU implementation.

Coupled voxel and crop convergence remain separate from A2's earlier cloud-only convergence. The current 97×65×97 field is an exact specified benchmark input, not a continuum cloud solution. Additional camera positions, cloud realizations and solar elevations are also needed. Production GPU sampling, angular-cache interpolation, traversal/time updates, reflections and complete image comparisons remain open. Cloud thermodynamics, entrainment, mixed-phase energy, conserved condensate and precipitation remain independent work packages.

## Algorithm and numerical-source lineage

Pharr, Jakob and Humphreys, *Physically Based Rendering*, fourth edition, “Volume Scattering Integrators” and “Transmittance”, provide the null-collision, volumetric path sampling and transmission framework. The A3 implementation is project code; no PBRT source is copied.

- https://pbr-book.org/4ed/Light_Transport_II_Volume_Rendering/Volume_Scattering_Integrators
- https://pbr-book.org/4ed/Volume_Scattering/Transmittance

The A1 spectral and colorimetric constants and their BSD-3-Clause attribution remain in `A1_THIRD_PARTY_NOTICES.md` and `data/a1/BRUNETON_LICENSE.txt`, originating from Eric Bruneton's reference demonstration. The A2 deterministic solver is the existing independently adapted Terluna implementation.

- https://ebruneton.github.io/precomputed_atmospheric_scattering/

`A1_METHODS.md`, `A2_METHODS.md` and `CLOUD_METHODS.md` retain the earlier physical assumptions, methods and source notices.
