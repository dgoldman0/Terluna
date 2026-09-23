# A1 — shared atmospheric state and optical reference checkpoint

Date: 22 September 2026. Scope: `immersion/light-cycle/shoreline`.

## Delivered state

A1 adds a reproducible reference path, an offline inspection lab, and a frozen-cloud transport benchmark to the saved revision-04 source. The shoreline runtime, its original sky atlas, terrain, materials, hydrology, and renderer migration remain unchanged. Production adoption of the shared atmospheric state is an open integration gate.

The source baseline is the saved `Open_Moon_Clouds_Source.zip`, whose methods identify commit `032765c51b1aeb140c732e8f484e47efda15cf4d` as the preceding revision. The baseline builder reproduces the audited shoreline HTML exactly, SHA-256 `c7ead1fd83e0d55518d2951e08cd46d31f323280885125382a8013b5d74636d9`. All 33 payload entries present in the supplied source manifest were checked against their recorded hashes. A1's incremental patch is relative to those saved revision-04 bytes. It is intended to follow that revision's publication, with a normal patch check against the actual current tree.

A1 has three connected reference components. `src/atmospheric-profile.js` preserves the supplied thermodynamic sounding, adds an explicit upper extension, and serializes an immutable profile. `src/profile-optics.js` evaluates spherical optical columns from that profile. The spectral and cloud solvers under `tools/a1` produce inspectable numerical reference data. The separate `Open_Moon_A1_Accuracy_Lab.html` displays those data and exports the exact profile and raw linear-radiance records.

## Shared profile and domain

The profile preserves the original weather-column pressure, temperature, vapor mixing ratio and height at every sounding node. It derives density using the moist equation of state and records dry-gas, water-vapor and total molecular number densities. Between nodes, pressure is interpolated logarithmically; temperature and mixing ratio are interpolated linearly. Thermodynamic quantities are then derived together. Number density is evaluated as `p/(k_B T)`; the mass density continues to include the moist composition correction.

The finite sounding ends at approximately 103.58 km in the lunar fair-weather experiment, while pressure remains approximately 13.47 kPa. An explicit hydrostatic continuation therefore completes the finite optical domain. Default top heights are 600 km for the lunar cases and 120 km for the Earth controls. The extension uses spherical gravity, a temperature relaxation to the final sounding temperature, and a specified vapor-mixing-ratio decay. Lunar defaults use 1-km integration steps, an 80-km temperature-transition length, and a 45-km vapor-decay length. The Earth defaults use 200-m steps, 15-km and 8-km lengths respectively. These are selected upper-atmosphere experiments, with no thermosphere, escape, photochemistry, or climate-equilibrium claim.

The profile includes a selected 300-Dobson-unit ozone triangle between 10 and 40 km, peaking at 25 km. The zero-ozone comparison sets the amount to zero while retaining its atmospheric sounding. The aerosol extinction is zero. Molecular scattering uses a dry-air-equivalent cross section applied to total molecular number density. Water-vapor absorption and composition-dependent molecular polarizabilities remain outside this first reference. These assumptions are serialized with the profile and affect its SHA-256 identity.

The nine supplied states combine lunar, Earth, and zero-ozone lunar configurations with the fair-weather, deep-convection and inversion-fog soundings. The clear-sky calculations omit their condensate to isolate molecular and ozone transport. The original parcel-derived cloud boundaries are shown only as conditional diagnostics.

The reference manifest also evaluates top-height and upper-temperature sensitivity for every state. For the lunar fair-weather profile, the vertical molecular column is 64,831.03 reference-density metres, compared with 58,080 m for the inherited `1.2 × 48.4 km` exponential model: approximately 11.62% greater. Extending the selected finite domain from 600 to 900 km increases this particular integral by approximately 0.00220%. The 200, 250 and 300 K upper-temperature alternatives change the integral by approximately −0.01644%, +0.03507% and +0.08452%. These are selected-input sensitivities for a vertical integral; they provide neither a global atmospheric uncertainty interval nor a twilight error bound.

## Optical integration and spectral sky

The JavaScript optical integrator handles finite segments, top-of-atmosphere exits and solid-planet occlusion. Adaptive Simpson quadrature accumulates reference-density air metres and ozone Dobson units, splitting at important altitude shells and the tangent point. The result reports convergence and evaluation count. The transmission helper rejects unconverged integrals.

A separately expressed Python implementation uses Gauss–Legendre integration along spherical paths. A further SciPy adaptive integration splits at all profile knots and independently interpolates pressure and temperature. The validation records compare these three implementations at four ray geometries per state: upward from the surface, a near-surface tangent, an oblique ray from 10 km, and a downward-looking unblocked ray from 60 km. Analytic and limiting-case tests also cover the exact dry isothermal spherical hydrostatic profile, vacuum/zero-length paths, opacity positivity, ozone normalization, reciprocal segments, and segment additivity.

`tools/a1/reference_optics.py` computes spectral first-order scattering for 48 wavelengths from 360 through 830 nm. It uses the published solar and ozone tables described below, a Rayleigh phase function, wavelength-dependent Beer attenuation, and explicit planetary shadowing. The reference assumes a collimated solar beam and a black surface, with no molecular multiple scattering, refraction or aerosol scattering. The visual atlas omits the solar disk.

Each state contains eight solar elevations (−6, −2, 0, 2, 6, 15, 45 and 90 degrees), 33 Sun-relative azimuth samples and 25 elevation samples concentrated near the horizon. The view integration uses 96 samples; solar-path quadrature uses Gauss order 12. Six selected rays per state are compared with 384 view samples and order-24 solar quadrature. Across the nine states, the largest tested relative spectral-L1 difference is approximately 0.313%; for the lunar fair-weather case it is approximately 0.0241%. These are differences within the first-order approximation, with 54 selected comparisons rather than a full-atlas error bound. Spectral-bin refinement and angular-atlas interpolation remain separate unclosed tests.

Raw spectral radiance is retained in the generated `.spectra.npz` files. Spectral integration through CIE 1931 color matching functions produces linear sRGB radiance in photopic-weighted, cd/m²-equivalent units. Negative out-of-gamut channels are retained in the data. The separate direct-beam RGB array represents photopic-weighted normal illuminance, in lux-equivalent units. The lab's exponential display curve, exposure slider, clipping and sRGB encoding operate only on the displayed copy. The displayed colors are first-order diagnostics; the completed sky requires higher scattering orders.

## Frozen cloud and genuine multiple scattering

`tests/a1_frozen_probe.js` executes within the exact baseline shoreline HTML. It samples the actual production density GLSL into a 49 × 33 × 49 Cartesian vertex grid, using float-target GPU readback. The sample volume spans approximately 36.99 km horizontally in both directions and 9.289–12.500 km in Cartesian height. Production sampling retains its radial-altitude calculation. The volume is a finite crop; extinction outside the box is zero. The crop itself is a selected transport boundary.

`src/frozen-cloud-field.js` defines a solver-neutral float32, little-endian representation, x-fastest vertex ordering, a trilinear sampler and the equivalent GPU texture coordinates. An immutable copy prevents caller mutation from changing the sampled field or its conservative majorant. Extinction and ice-fraction arrays are shared by identity between the GPU test and offline reference. The phase function matches the selected revision-04 mixture: 85% Henyey–Greenstein with `g = 0.78 + 0.04 f_ice`, and 15% with `g = −0.25`. The supplied fair-weather volume contains liquid condensate; mixed-phase cloud validation remains open.

At 24 deterministic off-grid probes, CPU versus GPU sampling of the frozen field differs by at most `5.215 × 10⁻⁹ m⁻¹`. A separate comparison with the continuous production density gives a normalized summed absolute extinction discrepancy of approximately 10.86%, with a largest absolute discrepancy of `5.464 × 10⁻⁴ m⁻¹`. The first result establishes close agreement in sampling this frozen grid. The second exposes its coarse representation of the original field. The 24 probes provide no whole-field bound.

The GPU transport harness evaluates first-order cloud-only scattering along 12 upward-looking rays, at budgets of 128/32, 512/128 and 2048/512 view/shadow samples. Both solvers use unit incident normal solar irradiance at 45 degrees, a gray conservative cloud, zero molecular atmosphere, and black external boundaries. Neither exposure nor the production cache enters the comparison. This diagnostic GPU integrator intentionally isolates physical first-order transport; it does not reproduce the production renderer's empirical diffuse and softened multiple-scattering terms.

`tools/a1/cloud_reference.py` evaluates both first-order and multiple-scattering radiance using analog delta-tracked free flight, ratio-tracked solar visibility, next-event estimation, the two-component phase distribution and Russian roulette. Its majorant bounds the trilinearly interpolated grid. Safety caps are counted as unresolved paths and fail validation. Sample means and standard errors are reported per ray and scattering treatment.

The delivered reference uses 65,536 independent histories per ray per treatment. All 12 first-order estimates agree with the highest-budget GPU result within six estimated standard errors plus an absolute radiance allowance of `10⁻⁵ sr⁻¹`. Direct transmission agrees with a separate 16,384-step CPU integration within the selected 0.001 absolute tolerance. There were zero unresolved paths. Six analytic or sampling checks pass: vacuum, pure-absorption Beer transmission, homogeneous single scattering, a conservative volume with isotropic boundary illumination, phase normalization/moment integration, and the sampled phase moment/direction normalization.

Multiple scattering materially changes several rays in this controlled volume. That observation motivates calibration of the production approximation. The present benchmark supplies no multiplicative correction for the production image, whose molecular atmosphere, diffuse boundaries and heuristic illumination remain outside this comparison. Grid refinement, crop expansion and additional independent-seed/reference-solver checks remain open.

## Recorded checks and limits

The full JavaScript regression suite passes 159/159 tests, including 22 new profile/optics/field tests. The baseline's 137 tests are retained. An initial run before generating the required baseline HTML failed its missing-build fixture; after the builder ran, the full suite passed. That setup failure is retained separately from the final regression record.

The cloud capture used Chromium with the software SwiftShader backend. It recorded zero browser errors, external requests and context-loss events for the tested capture. This is no consumer-GPU performance measurement or sustained shoreline traversal result.

The new lab passes 18 state/control checks across all nine atmospheric states, profile-export hash verification, desktop and 390-pixel mobile layout checks, and zero browser errors/external requests. Desktop and mobile screenshots were inspected. The lab is a reference-data inspector; its angular sky atlas and extinction slice are explicitly labeled diagnostics.

`A1_VALIDATION.json` records source and output hashes, aggregate numerical differences, test counts and the open gates. Raw generated data remain outside the incremental source patch but are included in the source/evidence bundle. Ordinary source, tests, compact validation and methods remain individually editable.

## Reproduction

Run from this directory. The source/evidence bundle includes the baseline vendor and asset bytes needed for an offline rebuild. The Python packages are listed in `requirements-a1.txt`; the browser harness expects Chromium at `/usr/bin/chromium` and Xvfb on PATH. Node.js 22.16.0 was used for the JavaScript suite.

```sh
python build_landscape.py
node tools/a1/export_profiles.cjs
python tools/a1/reference_optics.py
python tools/a1/capture_field.py
python tools/a1/cloud_reference.py --samples 65536
python build_accuracy.py
python tools/a1/check_lab.py
mkdir -p validation/a1
node --test tests/*.test.cjs > validation/a1/unit-tests.tap
python tools/a1/validate_a1.py
```

The reference generator accepts `--ids` to run selected states in separate batches. Numerical validation is tolerant rather than tied to identical Monte Carlo realizations across runtimes. Profile serialization and the delivered lab's profile export are checked byte-for-byte. Production shader noise may differ on another backend; a fresh field capture receives a fresh identity, and all downstream comparisons must use that same file.

## Next acceptance gates

The next optical step is molecular multiple scattering using the same atmospheric profile and spectral inputs, with radiometric and energy-budget checks before its results feed the interactive sky, environmental lighting, eye paths and solar attenuation. Cloud work should first establish voxel and crop convergence, then compare the full production illumination approximation against the fuller reference under identical gas and boundary conditions. Cache interpolation, translation and time-update errors need their own image-level comparisons.

The cloud-physics work remains independently staged: entrainment, condensate loading, consistent mixed-phase energy, and explicit condensate/precipitation accounting. The current cloud boundaries continue to describe the original selected undiluted parcel experiments. Terrain, vegetation and renderer migration remain separate tasks.

## Equation, algorithm and numerical-data sources

Eric Bruneton, *Precomputed Atmospheric Scattering*, reference implementation and documentation: https://ebruneton.github.io/precomputed_atmospheric_scattering/ . The A1 solver is independently implemented. The documentation informed the reference structure and explicit spectral/color treatment.

The solar and ozone numerical tables were taken from the reference's documented demonstration constants: https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/demo/demo.cc.html . They are the binned ASTM G-173 extraterrestrial solar spectrum and the selected Bremen 233-K ozone cross sections used there. A1 preserves their 48-value grid; their use does not establish spectral-bin convergence or a predicted lunar ozone distribution.

CIE 1931 tabulations and the XYZ-to-linear-sRGB matrix were taken from the same reference's constants: https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/constants.h.html . The originating source's BSD-3-Clause notice is retained in `A1_THIRD_PARTY_NOTICES.md` and `data/a1/BRUNETON_LICENSE.txt`.

Pharr, Jakob and Humphreys, *Physically Based Rendering*, fourth edition, “Transmittance” and “Volume Scattering Integrators”: https://pbr-book.org/4ed/Volume_Scattering/Transmittance and https://pbr-book.org/4ed/Light_Transport_II_Volume_Rendering/Volume_Scattering_Integrators . These provide the algorithmic basis for the independently implemented null-collision reference. No PBRT source was copied.

Original parcel equation sources and their tests remain in `CLOUD_METHODS.md` and the unchanged revision-04 source.
