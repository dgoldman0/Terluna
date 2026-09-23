# B1 — reference scene painter

## Scope and delivered interface

This checkpoint turns the A3 coupled transport primitives into full-frame light studies of an authored coastal scene. It is an offline numerical scene painter. The existing shoreline renderer and its data are preserved. The new viewer changes exposure only after radiance has been computed.

The first set contains a clear-air / broken-cloud Sun sweep at anchor elevations +60°, +30°, +12°, +6°, 0°, −2°, and −6°. Additional optical controls raise the same cloud by 20 km (at +6° and −2° Sun), multiply its extinction by 0.35, or multiply extinction by 2 (both at +6° Sun). Separate hemispheric skies and local incident-irradiance maps accompany these states. A black-surface control repeats the +12° fair-weather scene.

The local cloud field is a finite isolated crop. At low Sun, its ground shadow can fall tens or hundreds of kilometres away, outside the local map. An unshadowed local surface is consequently not evidence that cloud shadowing was omitted. Wider regional cloud fields remain a separate domain choice.

The Sun azimuth is fixed at 25° for this controlled sweep. These are sampled lighting states, not an ephemeris or a monthly weather animation. No formation/dissipation, precipitation, fluid evolution, trees, waves, buildings or eye adaptation are simulated.

## Inherited atmosphere and cloud

The serialized A1 fair-weather atmosphere supplies molecular density, planetary radius, ozone and upper continuation. A3's `next_collision`, `transmission`, scalar-Rayleigh sampling and cloud-phase routines are imported unchanged. Full path sampling uses molecular coefficients directly. A2's unconverged saved radiance field is not used as a source of scene illumination; its file is still read by the inherited input loader for provenance checks and the inherited spectral grid, solar irradiance and ozone cross sections. The angular radiance/moment values do not illuminate B1.

Straight rays retain A3's exclusions of refraction, polarization, Raman scattering, water-vapour absorption and aerosols. The only emitter is the Sun; Earthlight, airglow, thermal emission and artificial illumination remain outside this set.

The baseline condensate is A2's 97 × 65 × 97 frozen fair-weather field, with the same 36.99-km horizontal crop and approximately 9.289–12.500-km Cartesian vertical bounds. Gas continues outside that crop throughout the 600-km optical domain. The cloud has unit single-scattering albedo and the inherited 85% forward / 15% backward HG phase mixture.

The raised field is a Cartesian translation of those same coefficients by 20 km. Its optical phase, water content representation, horizontal extent and vertical optical-depth pattern remain fixed. It is a controlled height experiment, including an intentionally unchanged liquid-particle optical model. It is not a new thermodynamic solution for clouds at that temperature or altitude. The thin/thick controls scale extinction, not the cloud's geometric depth. Fog, a genuinely deeper layer, and deep convection remain future scene families.

## Geometry and surface boundaries

Coordinates match A3: local `(x,y,z)` in metres and a planet centre at `(0,-R,0)`. A 97 × 97 heightfield covers a 48-km square over a sea-level sphere. Fixed smooth ridges, an inlet, an island-like feature, and multiscale relief define an illustrative coast. The heightfield is triangulated; every eye, scattering and shadow path uses those same triangles. A 2-D grid traversal locates candidate triangles, with per-cell height bounds for empty-space rejection. Independent brute-force triangle tests validate the accelerator.

The mesh edge tapers below sea level. The sphere therefore closes the terrain body at the patch boundary. An initial preview exposed back-facing intersections when elevated terrain ended at an open vertical edge; the submerged boundary and its explicit check remove that geometry defect in the accepted runs. Development logs are retained separately.

The camera is placed 2 m along the geometric normal above the terrain at the selected coastal-overlook anchor near `x=-2200 m, z=-6500 m`. The anchor elevation is determined from the actual triangle intersection. The perspective camera has a 104° horizontal field of view, 10° yaw and 16° upward pitch. The camera is unchanged across the lighting experiments.

Pale shore, vegetation-coloured terrain and rock use explicit authored spectral Lambertian reflectances. They are proxies for bulk material colour, not measured terrain spectra or vegetation microgeometry. The far spherical boundary uses neutral reflectance 0.14. No empirical ambient-light term is added. Repeated surface/gas/cloud exchange follows the sampled light paths.

The smooth water sphere uses unpolarized dielectric Fresnel reflectance with a specified, nondispersive index 1.333. The reflected branch is traced with its Fresnel weight; transmitted energy goes into an absorbing deep-water reservoir. No seabed, underwater in-scattering, wave normals or shallow-water colour is supplied. This is a passive deep-water boundary with explicit energy partition, not a full ocean optical model.

## Transport and finite solar disk

Free flight uses the unchanged A3 superposed gas/cloud null-collision process. The nearest surface preempts any medium event beyond it. Absorption terminates the path. At a Lambertian surface, direct illumination is sampled explicitly and the continuation is cosine distributed with throughput multiplied by spectral reflectance. Water continues along its exact specular reflection with Fresnel throughput.

The Sun has a specified angular radius of 0.26667°. For unit normal incident irradiance its uniform disk radiance is `1 / [pi sin(radius)^2]`. Directions used in direct-light estimation are uniform in disk solid angle. The corresponding solid-angle integral factor is `2 / [1+cos(radius)]`. This keeps irradiance normalization explicit at horizon crossings and provides finite-disk shadowing.

Direct light is estimated once at every diffuse surface or volume event, with gas, cloud, planetary and terrain occlusion. Only camera or specular paths may encounter the emitted solar disk on escape. Those rules prevent double counting between explicitly sampled sunlight and diffuse continuations.

Unscattered camera-to-Sun and camera-to-water-to-Sun paths are evaluated separately with deterministic atmospheric/cloud line integrals and an 8 × 8 angular subpixel quadrature. The same terms are excluded from the stochastic camera calculation. This avoids rare, very large Bernoulli solar hits. These direct terms are also excluded from the presentation filter, preserving their angular footprint. Other paths involving a later water reflection, after gas/cloud or diffuse-surface scattering, remain stochastic and can have high variance.

Russian roulette begins after twelve interactions. Its survival probability is `min(0.95, max(0.05, throughput))`, with compensating throughput. Safety-limit, majorant and back-facing-intersection failures are counted. Accepted frames require zero incomplete histories. Completion establishes execution integrity; it does not establish image convergence.

## Spectral treatment and numerical uncertainty

The reference solar, ozone and CIE data still span all 48 supplied bins from 360–830 nm. The image calculation groups them into twelve consecutive four-bin groups. It integrates the supplied solar/CIE weights exactly within each group, evaluates Rayleigh extinction at that group's central wavelength, and interpolates the inherited ozone cross section there. Each group carries an independently evaluated transport coefficient. This is a **12-group spectral approximation**, not a 48-bin transport solve per pixel.

Random numbers are shared across groups for each pixel sample. This reduces some colour noise. Per-sample XYZ values are accumulated before estimating their variance, retaining the covariance introduced by shared random numbers. Camera rays are jittered within the pixel. Raw unit-transfer means, per-group sampling standard errors, XYZ values, XYZ sampling errors and surface-interacting contributions are retained in the frame NPZ files.

Tests first re-bin the measured A3 full spectra, then compare independent 48-bin and 12-group calculations on selected new surface-enabled B1 rays. These checks have limited ray/condition coverage. Their measured differences and Monte Carlo uncertainty are retained in the validation records. Spectral convergence across complete low-Sun images and other materials is still open.

Perspective frames use 160 × 96 traced pixels and 16 samples per spectral group (192 spectral histories per pixel). The default +12° clear/fair pair is additionally sampled at 64 histories per group; those refined results are delivered while the earlier 16-sample records are retained for comparison. Hemispheres use 80 × 40 pixels and 12 samples per group. These are intentionally coarse first scene studies. Rare specular paths and dim twilight have substantial sampling noise. The images are not certified radiometric reference images at a specified whole-image tolerance.

## Executed implementation and spectral checks

The 23 boundary/transport checks cover the triangle accelerator, sea intersections, cosine-direction sampling, passive reflectances, analytic Lambertian illumination, isotropic-sky water reflection, finite-Sun normalization, surface/molecule/cloud path completion, and the separation of direct solar terms. With terrain disabled and the original black spherical boundary restored, B1 reproduces A3's seeded path results to a maximum absolute difference of about 2.22e-16 in the selected test. The inherited 159 JavaScript regression tests also pass; they remain tests of the unchanged shoreline.

A separate check traces three surface-enabled scene rays at all 48 bins and at twelve groups, with 4,096 histories per wavelength/treatment. The relative Y differences and combined two-standard-error sampling allowances are approximately:

| Test ray | 12-group minus 48-bin Y | Combined 2 SE / 48-bin Y |
| --- | ---: | ---: |
| Fair cloud, +12° Sun, land-looking sample | −0.935% | 5.565% |
| Fair cloud, +12° Sun, upper-view sample | +1.591% | 2.981% |
| Raised cloud, +6° Sun | −1.449% | 3.593% |

These are reference-limited diagnostics. They do not establish a 1–2% spectral error bound or certify a full frame. Rebinning the earlier, much better sampled A3 spectra separately changes selected Y values by at most about 0.50%, but those rays do not include B1's new reflective boundaries. Raw results and the declared provisional gates are retained in `validation/b1/scene-spectral-check.json` and `boundary-tests.json`.

## Ground diagnostics

A 28 × 28 grid covers a 36-km square. Each query point is found by intersecting the same surface geometry; incident illumination is evaluated on its geometric normal. The direct term uses eight deterministic finite-disk directions and independent gas/cloud line integrals. The diffuse term uses cosine-distributed directions and the full scene path tracer, with the primary unscattered solar escape excluded. Multiplication by pi converts that directional average into incident irradiance.

The resulting data contain direct, diffuse and total incident XYZ, diffuse sampling errors, actual surface positions and normals. Y is in lux. The companion visibility map includes clouds, terrain and the horizon but excludes molecular attenuation; it therefore isolates direct-beam visibility, not the total illumination of a shadowed surface. Direct and diffuse terms are displayed separately. A small display-only marker identifies the observer and camera heading; map z increases downward. Full irradiance-map and finite-disk quadrature convergence remain open.

## Display and interpretation

Raw scene output is linear radiance. Y is photopic luminance in cd/m². Signed linear sRGB channels are preserved in the archives. Display clamps negative out-of-gamut channels, applies `1-exp(-RGB*2^EV/10000)` to the scene, then encodes sRGB. Ground irradiance maps use a common 60,000-lux-equivalent display scale. No white-balance correction or artistic colour grading is applied.

The default scene presentation averages nearby linear-light pixels within radius three, rejecting neighbours with different surface classes or substantially different surface depths. Cloud boundaries lack a geometric surface buffer and can soften. The analytically evaluated solar disk is added after filtering. This filtering is presentation only; it is neither additional transport sampling nor physical diffusion. The viewer exposes the raw samples and uses one exposure across the sequence. Spatial sampling, chromatic noise and the simple terrain remain visibly limiting.

The hemisphere and main frame are independent angular samples. The hemisphere is an equirectangular upper dome, with azimuth 0–360° horizontally and zenith-to-horizon vertically. It has no claim to angular image convergence. In particular, raised clouds can move mostly above the perspective frame; the dome is useful for seeing them.

## Reproduction and preservation

The complete schedule can be executed with `python tools/b1/run_checkpoint.py --workers 4`; `--viewer-only` rebuilds and checks the viewer from the supplied records.

Run from `immersion/light-cycle/shoreline`, with A3's generated input profile, cloud field and molecular input archive restored. `requirements-b1.txt` records the package versions used. The source/evidence distribution contains the necessary saved inputs; generating the images never requires external network access.

```
python tools/b1/test_b1.py
python tools/b1/check_spectral.py
python tools/b1/render_scene.py --scenarios clear fair --suns 60 30 12 6 0 -2 -6 --width 160 --height 96 --spp 16 --workers 4 --output data/b1/generated/final
python tools/b1/render_scene.py --scenarios high --suns 6 -2 --width 160 --height 96 --spp 16 --workers 4 --output data/b1/generated/final
python tools/b1/render_scene.py --scenarios thin thick --suns 6 --width 160 --height 96 --spp 16 --workers 4 --output data/b1/generated/final
# Repeat each scenario/Sun selection for the dome, with:
# --projection hemisphere --width 80 --height 40 --spp 12
python tools/b1/diagnostics.py --scenarios clear fair --suns 60 30 12 6 0 -2 -6 --size 28 --spp 12
python tools/b1/diagnostics.py --scenarios high --suns 6 -2 --size 28 --spp 12
python tools/b1/diagnostics.py --scenarios thin thick --suns 6 --size 28 --spp 12
python tools/b1/render_scene.py --scenarios fair --suns 12 --black --width 160 --height 96 --spp 16 --workers 4 --output data/b1/generated/controls
# Higher-sample default pair; replaces their earlier coarse records.
python tools/b1/render_scene.py --scenarios clear fair --suns 12 --width 160 --height 96 --spp 64 --workers 4 --output data/b1/generated/final
python build_b1_painter.py
python tools/b1/check_browser.py
python tools/b1/validate_b1.py
```

The package contains stable source files, compact evidence and the generated records. Generated radiance arrays, full viewer HTML and image files remain outside the incremental Git patch. Every frame pins the profile, baseline and transformed cloud, geometry, executed source and seed. A separate environment record pins the software versions. Baseline preservation checks compare inherited source bytes and the rebuilt shoreline HTML.

## Delivery checks

The delivered set contains 18 perspective lighting states, 18 independently sampled sky hemispheres, 18 incident-irradiance maps, and one black-surface perspective control. The accepted image/map records contain 84,054,528 spectral path histories, with zero unresolved histories. Another 737,280 histories support the independent selected-ray spectral comparison. Shared random numbers across spectral groups mean that the image histories are correlated across wavelengths.

All 19 implementation/delivery gates pass. They cover executed tests, exact input/source identities, finite arrays, correct dimensions/path counts, contribution accounting, complete scene sets, pixel-level browser behavior and preservation of the inherited source snapshot. The full-scene physical-accuracy qualification remains false. The standalone viewer passes 339 state/control checks across desktop, high-DPI desktop and mobile layouts, including data isolation from exposure and actual downloaded metadata bytes.

The local-file navigation attempt is recorded separately from in-memory rendering of the exact standalone HTML. No browser policy is changed to obtain a pass. This establishes the tested display/control path; it provides no claim about consumer-GPU rendering performance or sustained traversal in the production shoreline.

## Next gates

Increase image sampling, pixel resolution and the geometric detail near the observer. Validate full-image spectral grouping, solar subpixel integration, global energy closure, reflecting-surface and coupled cloud crop/grid sensitivity. Add finite-depth/rough water, then genuinely distinct fog, layered-cloud and convective fields with their stated physical assumptions. A real-time approximation can be tested against these full frames once the relevant reference errors are measured. The production shoreline and its caches have not been replaced by this checkpoint.

## Algorithm sources

Pharr, Jakob and Humphreys, *Physically Based Rendering*, fourth edition, provides the volumetric/surface transport framework, Lambertian reflection and dielectric Fresnel treatment. The B1 implementation is project code. No PBRT source was copied.

https://pbr-book.org/4ed/Light_Transport_II_Volume_Rendering/Volume_Scattering_Integrators

https://pbr-book.org/4ed/Reflection_Models/Diffuse_Reflection

https://pbr-book.org/4ed/Reflection_Models/Dielectric_BSDF

The inherited spectral inputs and their notices remain in `A1_THIRD_PARTY_NOTICES.md` and `data/a1/BRUNETON_LICENSE.txt`. `A1_METHODS.md`, `A2_METHODS.md` and `A3_METHODS.md` retain the earlier assumptions and unresolved convergence findings.
