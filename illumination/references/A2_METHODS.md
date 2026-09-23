# A2 — molecular multiple scattering and cloud convergence checkpoint

Date: 22 September 2026. Scope: `immersion/light-cycle/shoreline`.

## Purpose and boundary

A2 extends the A1 reference path without changing the production shoreline renderer. A1 established one immutable thermodynamic/optical atmosphere, a first-order spectral sky, and a solver-neutral frozen cloud field. A2 answers two narrower questions before any visual integration:

1. How large is the molecular multiple-scattering correction when the A1 atmosphere is used consistently through the full visible spectrum?
2. How much of the frozen-cloud transport result is controlled by voxel resolution versus the finite horizontal crop?

The calculations remain conditional on the selected A1 lunar fair-weather profile. They are numerical reference experiments rather than a climate-equilibrium prediction or an empirical reconstruction of lunar weather. Cloud thermodynamics, morphology, condensate retention, and the A1 upper-atmosphere continuation retain their documented input status.

The generated A2 lab is an evidence viewer. Display exposure and sRGB encoding do not feed back into the stored radiance. The shoreline sky atlas, cloud-lighting approximation, cache, environment lighting, and tone mapping remain untouched by this checkpoint.

## Shared atmospheric state

A2 reads `data/a1/generated/moon-fair.profile.json` through A1's reference-profile loader. Pressure is interpolated logarithmically, temperature linearly, and the molecular density ratio relative to `101325 Pa / (k_B 288.15 K)` is evaluated as

```
n(z) / n_ref = p(z) / [k_B T(z) n_ref].
```

The finite sounding and its explicit A1 hydrostatic continuation supply the density at every path sample. The spherical lunar radius and the selected 600-km optical top come from the same serialized state. The selected ozone distribution is the A1 300-DU triangular profile; each path segment separately accumulates molecular scattering and ozone absorption.

The spectral inputs remain the 48 bins from 360 through 830 nm at 10-nm spacing. Molecular Rayleigh scattering uses

```
beta(lambda) = 1.24062e-6 (lambda / 1000 nm)^-4 m^-1
```

at the reference density, then scales locally with the shared molecular profile. Ozone absorption uses the same cross sections and solar spectral irradiance recorded by A1. Aerosol extinction, water-vapor absorption, refraction, Raman scattering, polarization, airglow, and a thermospheric/photochemical solution remain outside this reference.

## Successive molecular scattering orders

`tools/a2/molecular_multiple_scattering.py` generalizes the repository's existing spherical scalar-Rayleigh successive-order solver from an exponential gas profile to the A1 atmosphere. No source from the production renderer is used for the scattering solution.

For each scattering order, the angular source is represented by the exact scalar-Rayleigh zeroth- and second-moment reduction already used by the Terluna spectral transport prototype. The Rayleigh source factor is

```
3 / (16 pi) * [J + d_i K_ij d_j],
```

with axial symmetry about the Sun reducing the stored angular moments. Spherical rays are split at radial-shell crossings; segment attenuation uses the shared molecular density and ozone profile. The direct beam is evaluated on a separate, finer radial/solar-elevation table. The ground boundary is black in A2 so that the atmospheric energy and convergence checks isolate molecular transport. The Sun is collimated in this solver; the finite disk remains outside A2.

Scattering wavelengths are elastically uncoupled. The fine full-spectrum solve is therefore executed as eight independent six-band chunks and concatenated in wavelength order. This is algebraically equivalent to carrying all 48 uncoupled wavelengths in one array at the same angular and radial grids; the split avoids unnecessary memory and compilation scaling.

Three numerical grids are available:

| Quality | Radial nodes | Polar quadrature | Azimuth quadrature | Radial shells | Solar-path samples | Core solar/elevation spacing |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Draft | 16 | 8 | 8 | 32 | 128 | 3° |
| Standard | 22 | 12 | 12 | 48 | 192 | 2° |
| Fine | 28 | 16 | 16 | 64 | 256 | 1.5° |

Angular polar quadrature is split into four intervals around the horizon so near-horizontal transport receives additional resolution. Radial nodes and shell boundaries are concentrated toward the surface and include important sounding/ozone heights. Each order is accumulated into the total field. For the fine reference chunks, the solve continues until the largest tracked order increment is below 0.05% of the accumulated peak moment.

The required order depth is strongly wavelength dependent. The fine chunks required 56, 32, 20, 14, 11, 9, 8, and 7 orders for 360–410, 420–470, 480–530, 540–590, 600–650, 660–710, 720–770, and 780–830 nm respectively. Their largest final increment is `4.94725e-4`. The blue end therefore demands much deeper iteration than the red because the selected Rayleigh optical depth is much larger.

## Independent first-order cross-check

The multiple-scattering solver contains a different numerical path construction from A1, so its first scattering order can be compared against the previously delivered first-order reference without changing the physical closure. The full standard-grid atlas reproduces A1 with global spectral L1 difference `0.00104155` (0.104%). Among pixels whose A1 spectral sum exceeds 1% of the atlas maximum, the median relative spectral L1 is 0.0704% and the 95th percentile is 0.800%. This gate checks implementation consistency; it is not a comparison with empirical sky radiance.

Representative-wavelength development checks also showed the independent first-order mismatch declining with grid refinement: approximately 0.185% at draft, 0.107% at standard, and 0.0718% at fine on the selected test rays.

## Standard atlas and fine-grid convergence

The reusable standard numerical field is saved in `data/a2/generated/moon-fair-molecular-standard-full.npz`. `tools/a2/build_multiple_atlas.py` evaluates it on the same A1 angular atlas: eight solar elevations, 25 view elevations, 33 Sun-relative azimuths, and all 48 wavelengths. Both first-order and accumulated multiple-scattering spectra are retained in the `.npz` archive and the compact JSON lab payload.

Fine-grid validation is evaluated on seven deliberately different rays: zenith, toward-Sun mid-sky, low toward-Sun sky, low-Sun forward and side views, a horizon-Sun case, and a view away from the high Sun. After the eight fine spectral chunks are merged, standard-to-fine convergence is:

- overall spectral L1 across the seven rays: 0.3934%;
- largest individual-ray spectral L1: 0.5438%;
- largest individual-ray photopic-luminance relative difference: 0.8056%.

These figures quantify standard-versus-fine numerical sensitivity only. They are not total physical error bounds.

Molecular multiple scattering is a large physical correction in this selected atmosphere. The fine-grid multiple/first-order photopic-luminance ratios for the seven rays are approximately 2.090, 1.883, 2.105, 2.824, 5.359, 3.979, and 2.539. Ratios become especially large where first-order radiance is weak, so the raw spectra and absolute radiance remain the primary record. A ratio by itself must not be interpreted as a global brightness multiplier for the production image.

## High-Sun energy diagnostic

A2 numerically integrates the downward hemispheric multiple-scattered radiance at the surface and adds the attenuated direct beam. For solar elevations 15°, 45°, and 90°, the result is divided wavelength by wavelength by the local top-of-atmosphere horizontal spectral flux. The maximum ratios across the 48 wavelengths are 0.8129, 0.8853, and 0.9139 respectively. The validation gate requires all three high-Sun maxima to remain below 1.01.

This is a sanity check rather than a global conservation proof. Very low solar elevations are intentionally excluded because horizontal transport on a spherical planet can move energy laterally into a local column while the local projected top-boundary denominator approaches zero. A complete spherical energy closure would integrate over the planetary boundary rather than divide one surface point by its local horizontal incident flux.

## Cloud field resolution

The A1 frozen field was a 49 × 33 × 49 Cartesian vertex grid. `tests/a2_frozen_probe.js` executes the exact revision-04 production density GLSL and captures three grids over the same physical crop:

- 25 × 17 × 25;
- 49 × 33 × 49 (A1 resolution);
- 97 × 65 × 97.

Ninety-six common deterministic Halton probes evaluate the continuous production density and the trilinearly interpolated frozen field at identical off-grid locations. The normalized sum of absolute extinction differences divided by the continuous-field extinction sum decreases monotonically:

- 38.72% at 25 × 17 × 25;
- 16.00% at 49 × 33 × 49;
- 4.53% at 97 × 65 × 97.

The 96-probe A1 figure differs from A1's earlier 24-probe statistic because A2 intentionally uses a denser common probe set. Neither statistic is a rigorous whole-volume error bound.

A deterministic first-order cloud-only transport calculation then holds the ray set and sample budgets fixed while changing only the frozen field. Relative to the 97 × 65 × 97 field, summed first-order radiance L1 differences are 14.00% for the coarse grid and 5.77% for the A1 grid. Direct transmittance differences are 1.50% and 0.173% respectively. The fine grid is the comparison endpoint, not a proof of continuum convergence; another refinement remains an open gate.

## Finite cloud crop

Crop sensitivity is measured separately from voxel resolution. The horizontal crop sizes are 0.75×, 1.0×, 1.5×, and 2.0× the A1 extent. Horizontal dimensions are 37, 49, 73, and 97 vertices respectively, chosen so the horizontal voxel spacing remains exactly equal to A1. The vertical grid remains 33 vertices. A common set of central rays is used throughout.

The 0.75× crop strongly truncates first-order illumination paths: its radiance differs from the 2× reference by 57.82%. The 1.0× and 1.5× crops agree with the 2× result on these central deterministic rays to below `5e-7` relative L1. Direct transmittance is similarly stable. This establishes crop convergence for the tested central first-order paths, not for every observer position or cloud realization.

## Cloud multiple-scattering sensitivity

`tools/a2/cloud_convergence.py` reuses A1's null-collision Monte Carlo reference on four nonzero common rays with 32,768 histories per ray and configuration. It reports means and standard errors and records zero unresolved safety-cap paths.

Mean multiple-scattering sensitivity relative to the fine grid is 5.33% for the coarse grid and 2.30% for the A1 grid. Relative to the 2× crop, mean differences are 3.24% for 0.75×, 2.01% for A1, and 2.02% for 1.5×. These Monte Carlo differences are deliberately reported as sensitivity diagnostics rather than deterministic acceptance errors: several crop differences are comparable to their sampling uncertainty. Higher sample counts or an additional independent estimator are required before using them to set a production crop margin.

Every field capture records its profile identity, field hash, GPU backend state, and CPU/GPU trilinear sampler comparison. All six final captures report intact graphics contexts. The largest GPU/CPU frozen-field sampler extinction discrepancy remains below `1e-7 m^-1`.

## Validation gates

`A2_VALIDATION.json` is the machine-readable checkpoint. The numerical gates require:

- the complete JavaScript regression suite to pass;
- exact 360–830 nm / 10 nm coverage in the fine chunks;
- every fine spectral chunk to end below the 0.05% order-increment threshold;
- standard/fine spectral convergence below the selected 1% aggregate and 2% luminance gates;
- full-atlas first-order reproduction of A1;
- nonnegative higher-order added radiance;
- high-Sun downward spectral energy below the local incident-horizontal sanity bound;
- intact, hash-matched cloud captures with CPU/GPU sampler agreement;
- monotonic frozen-field improvement with grid resolution;
- the A1 grid to approach the fine first-order transport result;
- the A1 crop to agree with the 2× crop on the tested central first-order paths;
- zero unresolved cloud Monte Carlo histories.

The current suite passes 159/159 JavaScript tests. The standalone A2 lab is tested separately in Chromium. Its browser record verifies all eight solar elevations in both sky modes, all seven fine-grid rays, exposure isolation from the evidence arrays, validation export identity, desktop and 390-pixel mobile layout, and zero browser errors or external requests.

## Reproduction

Run from `immersion/light-cycle/shoreline` after the A1 evidence has been generated. A2 uses the same Python packages listed in `requirements-a1.txt`.

```sh
# Capture the six frozen cloud fields. The script may split heavy captures across
# fresh Chromium sessions on software renderers with limited long-lived resources.
python tools/a2/capture_cloud_convergence.py

# Deterministic and Monte Carlo cloud convergence.
python tools/a2/cloud_convergence.py --samples 32768

# Standard 48-band molecular field.
python tools/a2/molecular_multiple_scattering.py \
  --quality standard --all-wavelengths --orders 48 \
  --save-field data/a2/generated/moon-fair-molecular-standard-full.npz \
  --output data/a2/generated/moon-fair-molecular-standard-full.json

# Fine reference. Run each six-band interval with --quality fine and enough
# orders for the selected <0.05% increment gate, then merge them.
python tools/a2/merge_molecular_chunks.py

# Standard-grid full sky atlas and all numerical gates.
python tools/a2/build_multiple_atlas.py
python tools/a2/validate_a2.py

# Offline evidence viewer and browser checks.
python build_a2_accuracy.py
python tools/a2/check_lab.py
```

The delivered per-chunk JSON files record the exact quality, order count, wavelength interval, input profile hash, and convergence history. Reproduction should use those records instead of assuming one common order count is adequate across the spectrum.

## Interpretation and next gates

A2 establishes that the shared-profile molecular background needs substantial multiple scattering and that the standard A2 solution is close to the tested fine-grid solution. It also shows that the A1 frozen cloud grid remains visibly coarse as a representation of the continuous production density even though the A1 crop is already adequate for the tested central first-order paths.

The next transport milestone should combine the shared molecular atmosphere and a frozen cloud in one reference calculation. That coupled solver can then compare the production cloud illumination terms against a physically consistent gas-plus-cloud reference at matched Sun, observer, density field, particle phase function, and boundary conditions. Production cache/interpolation error should remain a separate measurement from physical-lighting error.

The next cloud-physics milestone remains independent: entrainment, condensate loading, mixed-phase thermodynamics, conserved condensate, and precipitation/virga accounting. A2 does not change the original conditional cloud-base/ceiling experiments.

## Source lineage

A2 inherits A1's recorded spectral data, colorimetry, shared atmospheric state, and external reference notices. `A1_METHODS.md` and `A1_THIRD_PARTY_NOTICES.md` retain those sources and licenses. The molecular successive-order formulation derives from the existing Terluna spectral transport prototype and has been independently adapted to the A1 profile. The cloud Monte Carlo algorithm and its source lineage remain documented in A1.

No external source bytes are added by A2. All new cloud fields, sky spectra, screenshots, and validation records are generated from the project source and the already pinned A1 inputs.
