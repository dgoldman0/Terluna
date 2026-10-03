# Optical comfort and surface vision

The committed light fluxes give dimmer horizontal matte ground and more strongly filled solar-only shadows on the design Moon at noon. Transferring the repo's computed angular sky patterns onto those fluxes gives 1.83 times Earth's ambient illumination for a level eye over dark ground, compared with 2.11 under a uniform sky. Looking down reverses the ordering. Finite pale surroundings can bring the two worlds' eye illumination close together, while physical roofs, side screens and shaded ground change both illumination and surface contrast. Optical comfort therefore needs a specified surface, gaze and shelter geometry.

**Evidence state:** reproducible screening from committed clear-sky products, with prescribed scene geometries. Sections 1–8 retain the initial flux-based calculation and imposed sky shapes. Section 9 adds computed angular patterns transferred from a different atmospheric model, and finite scenes with repeated matte reflection. The quantities describe light and selected contrast mechanisms. Human discomfort, visual performance in real landscapes and eye safety remain unvalidated. All albedos are spectrally neutral scenario values; a dark landscape represents open ground, rather than a modelled forest interior.

## Inputs and reproduction

The input is [surface_light.json](../../../illumination/surface_light/results/surface_light.json), checked at repository commit `88bb4f5e967804f6c0da01ba0ad2c89d95453516`. Its recorded producer hashes match the checked-out files. The main pair is the 1.2-atm design Moon behind the chosen titania film and 5% additional dimming, versus the unfiltered Earth control. Both controls are clear and aerosol-free. They establish a matched atmospheric comparison, rather than a distribution of typical Earth or lunar weather. The 1.0-atm design Moon supplies a pressure sensitivity. The atmosphere, shield and photopic spectral weighting are inherited unchanged.

The runner consumes the stored photopic band summaries, whose albedo feedback was evaluated on the original fine spectral grid. Reconstructing illuminance from the saved 1-nm spectra and reading rule agrees with all 364 upstream summaries within 0.177%; the largest difference is at 10° Sun and albedo 0.8 in the design case. This checks consistency between the stored forms, rather than physical validity. The runner also reads the committed flux comparisons with the spherical sky atlases. Neither external spectral downloads nor the full local simulation drive is required to reproduce this study.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.run
python -m pytest research/studies/optical_comfort -q
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/optical-comfort/plot.py
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.directional
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/optical-comfort/directional_scenes.py
```

The stable outputs are [optical_comfort.json](results/optical_comfort.json), schema `terluna.research.optical-comfort/1`, and [surface_scenes.csv](results/surface_scenes.csv). The JSON records the full input and implementation SHA-256 hashes, assumptions, units, reading rule and numerical convergence. It contains 273 combinations of world, solar elevation and landscape albedo, 80 gaze cases, plus local surfaces, angular masks, glare sensitivities, smooth-water reflection and low-Sun proxies. Its provenance depends on the input and implementation bytes, so a commit-message change leaves the product unchanged. The plotting code reads that product and writes ignored figures with a hash manifest.

| Evidence level | Quantities in this pass |
|---|---|
| Derived from the supplied horizontal fluxes | Horizontal matte luminance, ideal solar-only shadow fraction and contrast |
| Additional angular assumptions | Eye-plane and vertical-surface illumination, gaze dependence, angular masks, annular surface glare and ideal water highlights |
| Still unresolved | Actual scene radiance, adaptation, perceived discomfort, visual-task performance and exposure statistics |

## 1. Ground brightness and the light reaching the eye

Let B and D be direct and diffuse horizontal illuminance, and E = B + D. For a neutral Lambertian surface with reflectance rho,

\[
L_{\rm horizontal}=\rho E/\pi.
\]

For a uniform sky over a flat, uniformly illuminated landscape of albedo A, the vertical eye plane receives

\[
E_{v,\rm ambient}=(D+AE)/2,
\qquad
E_{v,\rm Sun}=B\cot h\max(\cos\Delta\phi,0),
\]

where h is solar elevation and delta-phi is the viewing azimuth relative to the Sun. The ambient term excludes the disk; it includes sky and ground. It is also the illumination on a vertical matte face oriented sideways to, or away from, the Sun. These are plane illuminances before eyelids, the nose, pupils and ocular transmission. They cannot be converted directly into retinal exposure or an adaptation state.

For the collimated Sun used in plane projections, its contribution to a vertical plane is zero when overhead. The finite solar disk would add a small edge contribution; the tables use the collimated approximation. The last two columns below additionally assume the uniform sky.

| Landscape albedo | Moon horizontal light, klux | Earth horizontal light, klux | Moon ground luminance, cd/m² | Earth ground luminance, cd/m² | Moon vertical ambient, klux | Earth vertical ambient, klux |
|---|---:|---:|---:|---:|---:|---:|
| 0.1 | 89.6 | 124.6 | 2,852 | 3,966 | 20.3 | 9.65 |
| 0.3 | 96.1 | 126.3 | 9,178 | 12,061 | 33.5 | 23.2 |
| 0.5 | 103.7 | 128.1 | 16,512 | 20,381 | 48.9 | 37.2 |
| 0.8 | 118.1 | 130.8 | 30,072 | 33,306 | 77.3 | 58.8 |

At albedo 0.1, noon ground luminance is 28% lower on the Moon, while the uniform-sky estimate of vertical ambient illumination is 2.11 times the Earth control. At albedo 0.8, those differences become 10% lower and 1.31 times higher. High-albedo ground is a bright extended source in both environments. Terluna adds stronger diffuse horizontal illumination and ground–atmosphere recycling.

The distinction between landscape albedo and a local object's reflectance matters. Painting a small wall changes that wall's reflected light; applying a new albedo to the whole atmospheric column would also change the sky and overstate the intervention. The runner keeps those operations separate.

For example, a matte surface of reflectance 0.5 in the albedo-0.1 noon landscape has these luminances:

| Orientation | Moon, cd/m² | Earth, cd/m² |
|---|---:|---:|
| Horizontal, fully exposed | 14,261 | 19,831 |
| Vertical, unobstructed surroundings | 3,235 | 1,537 |
| Horizontal, solar disk alone screened | 5,044 | 1,090 |

The brighter vertical surface is relevant to facades and the sides of objects. Actual trunks, leaves and building recesses also occlude the sky and ground, so their luminance needs their own geometry.

The 1.0-atm Moon gives 18.9 klux of vertical ambient light at noon over albedo 0.1, against 20.3 klux in the main design. This pressure change moderates the effect without removing it.

## 2. Looking toward, away from and below the Sun

At 30° solar elevation over albedo-0.1 ground:

| Quantity | Moon | Earth |
|---|---:|---:|
| Horizontal illuminance | 35.9 klux | 57.8 klux |
| Vertical eye plane, sideways or away | 12.1 klux | 5.86 klux |
| Vertical eye plane, toward the Sun | 38.6 klux | 95.6 klux |
| Direct contribution to the sunward eye plane | 26.5 klux | 89.8 klux |

Thus the uniform-sky scenario combines brighter ambient illumination with a substantially weaker direct source. The Sun remains compact and highly luminous. Its reduced flux does not establish that direct viewing is safe.

For gaze elevation e and viewing azimuth delta-phi, the same hemispheres give eye-plane illuminance

\[
E_{\rm eye}=\frac{D}{2}(1+\sin e)+\frac{AE}{2}(1-\sin e)
 +\frac{B}{\sin h}\max(\sin e\sin h+\cos e\cos h\cos\Delta\phi,0).
\]

Negative gaze elevations look down. With the Sun overhead, all the level and downward views in the following table exclude its collimated beam. The units are klux on the plane normal to the gaze; these are not vertical illuminances except in the level row.

| Gaze at noon | Moon, albedo 0.1 | Earth, albedo 0.1 | Moon, albedo 0.8 | Earth, albedo 0.8 |
|---|---:|---:|---:|---:|
| Level | 20.3 | 9.65 | 77.3 | 58.8 |
| 30° down | 14.6 | 11.1 | 85.9 | 81.7 |
| 60° down | 10.5 | 12.1 | 92.2 | 98.5 |
| Straight down | 8.96 | 12.5 | 94.5 | 104.6 |

The Moon–Earth ordering reverses in steep downward views in this model. On pale terrain, looking down increases the lunar eye-plane illuminance even as it makes that illumination lower than Earth's. A claim that the Moon exposes the eye to more light in every direction would therefore be wrong even under the original uniform-sky assumptions. Eye illuminance also remains distinct from discomfort: a change of gaze can change source positions and task contrast as well as the total flux.

The uniform-sky calculation is a transparent angular approximation. Horizontal diffuse illuminance alone does not specify vertical illumination. To expose that uncertainty, the study holds horizontal flux fixed and tests three axisymmetric shapes, with radiance proportional to 1, 1 + 2 sin(elevation), or 3 − 2 sin(elevation). Their vertical-to-horizontal diffuse ratios are 0.500, 0.396 and 0.645.

These shapes give 17.0–24.9 klux for a level gaze in the dark-ground lunar noon scene and 8.94–10.7 klux for Earth. This interval describes three imposed shapes, rather than a confidence interval or a bound on the real sky. Horizontal diffuse flux alone places no finite upper bound on vertical sky illumination if arbitrarily strong radiance near the horizon is allowed. The actual scattering solution must supply the angular distribution. The missing circumsolar structure, azimuth dependence and horizon brightness can change a particular view. The committed packed atlas supplies a reduced-resolution angular input in section 9. Its proxy profiles and unfiltered source still require a transfer assumption to use with the chosen design.

## 3. Shadows, shape cues and material contrast

Consider a small horizontal patch whose direct Sun is occulted while its sky exposure and surrounding landscape remain unchanged. Its illumination is D. The shadow-to-sunlit ratio is D/E, and the Michelson contrast between equally reflective sunlit and shadowed patches is

\[
C_{\rm shadow}=\frac{E-D}{E+D}=\frac{B}{B+2D}.
\]

| Scene | Moon shadow retains | Earth shadow retains | Moon shadow contrast | Earth shadow contrast |
|---|---:|---:|---:|---:|
| Noon, albedo 0.1 | 35.4% | 5.50% | 0.477 | 0.896 |
| Noon, albedo 0.8 | 51.0% | 9.97% | 0.325 | 0.819 |
| Sun 30°, albedo 0.1 | 57.3% | 10.3% | 0.272 | 0.814 |

This is a substantial reduction in shading cues and in the dark relief provided by a small solar blocker. It could soften the appearance of relief and small surface features. Filled shadows can also reveal detail that would be dark on Earth. The balance depends on the task.

For adjacent matte patches under the same illumination, ordinary reflectance contrast remains `(rho_target − rho_background) / rho_background`. Uniformly increasing diffuse light leaves that ratio unchanged. Predicting that all detail becomes washed out would conflate shadow contrast with material contrast. Veiling light in the eye, atmospheric path radiance, specular reflections and task geometry can each reduce contrast through a different mechanism.

## 4. What different forms of shade remove

A hat can remove the direct Sun when its geometry permits. Reducing ambient illumination additionally requires blocking the bright sky and the ground. An eye-level overhang that blocks only the upper sky leaves much of the vertical sky flux near the horizon.

For a uniform sky, the fraction of its vertical illumination originating below elevation alpha is

\[
f(\alpha)=\frac{2}{\pi}(\alpha+\sin\alpha\cos\alpha),
\]

with alpha in radians. Blocking sky above 45° therefore leaves 81.8% of the sky contribution to a vertical eye plane. Blocking above 15° leaves 32.6%.

The following are dark angular masks for a level gaze, with the distant ground remaining sunlit. All also screen the solar disk. Roof reflections and changes to nearby shaded ground are omitted; these are component sensitivities, not performance claims for a built awning.

| Lunar noon mask | Eye illumination, albedo 0.1 | Eye illumination, albedo 0.8 |
|---|---:|---:|
| Solar disk alone | 20.3 klux | 77.3 klux |
| Disk and sky above 45° | 17.4 klux | 71.9 klux |
| Disk and sky above 15° | 9.64 klux | 57.0 klux |
| Disk and all sky | 4.48 klux | 47.2 klux |
| Disk and 90% of each hemisphere's contribution | 2.03 klux | 7.73 klux |

On pale terrain, screening the overhead hemisphere alone leaves a large ground contribution. Recessed spaces, lateral screens and shade extending across nearby ground are consequently worth modelling. Their material reflectance and view geometry need to be included together. Vegetation can provide both shade and darker surfaces, although a forest floor requires an actual canopy-light calculation.

An ideal neutral filter scales scene luminance and illuminance together. It can lower the light level, but in the simple additive-veil model it scales both the target and the veil equally, leaving their contrast ratio unchanged. Selectively screening a glare source changes that ratio. Pupil changes, adaptation, spectral tint and real eyewear optics are outside this calculation. Polarization is also outside the current scalar sky model.

## 5. A restricted disability-glare calculation

The implemented CIE age-adjusted Stiles–Holladay relation is

\[
L_{\rm veil}=10\left[1+(\mathrm{age}/70)^4\right]
\frac{E_{\rm eye,source}}{\theta_{\rm deg}^2},\qquad 1^\circ<\theta<30^\circ.
\]

The official CIE description and Vos's author abstract provide the equation's conventions and validity range ([sources](sources.json)). The runner rejects angles outside this range. Ages 30, 50 and 70 are sensitivities; they do not describe all eyes of those ages. The solar disk is treated as a compact source at separations of 5°, 10° and 20°. No whole-sky glare score is constructed by feeding total diffuse lux into this point-source relation.

If an additive veil is approximately uniform across a small target and its background, the retained Weber contrast is `L_background / (L_background + L_veil)`. For an extended surround, this describes the change from the same eye and task under the reference surround; it does not model the baseline ocular transfer function. The kernel is achromatic and uses photopic units. Wavelength-dependent ocular scatter and the effects of the changed spectrum on colour discrimination remain outside this pass.

For a 50-year age parameter, a source separation of 10° and Sun elevation 30°, the direct Sun produces an estimated veil of 3,804 cd/m² on the Moon versus 12,869 cd/m² on Earth. An independently prescribed 1,000-cd/m² task background then retains 20.8% and 7.21% of its original contrast. This is a conditional near-Sun task geometry; it is not the horizontal landscape view in section 2, nor a safe-viewing calculation.

Surface glare is also screened through an imposed bright annulus between 5° and 25° from fixation. An equally illuminated matte surround with reflectance 0.8 around a task background of reflectance 0.1 gives an 8:1 luminance ratio. Integrating the **excess luminance above the task background** gives retained contrast of 82.3%, 79.3% and 70.7% for the three age parameters. The annular quadrature agrees at 64 and 128 nodes to about 2.4 × 10⁻¹⁶ in veil per unit excess luminance; a separate fine midpoint integration tests the implementation.

That normalized surface result is the same on Earth and the Moon when the luminance ratio and geometry are the same. An extended bright surface can impair a dark task, while the total illumination amplitude alone does not establish a larger fractional disability effect. The annulus is a sensitivity to one explicitly defined layout, rather than a claim about the complete field or visual acuity.

## 6. Water and other specular surfaces

The input ground model is Lambertian. It cannot produce a water highlight. The first pass therefore evaluates a separate limiting geometry: a horizontal water patch that resolves the reflected image of a uniform solar disk, viewed in the specular direction at a downward gaze equal to the Sun's elevation. The refractive index is 1.333, reflection follows the unpolarized Fresnel relation, and water-to-observer path attenuation is neglected. These values describe a small patch immersed in the albedo-0.1 illumination case; a sea's contribution to the full atmospheric boundary remains unsolved.

| Sun elevation | Fresnel reflectance | Lunar resolved highlight | Earth resolved highlight |
|---|---:|---:|---:|
| 90° | 2.04% | 17.4 million cd/m² | 35.3 million cd/m² |
| 30° | 5.97% | 26.9 million cd/m² | 91.1 million cd/m² |
| 10° | 34.8% | 39.9 million cd/m² | 340.5 million cd/m² |

A specular highlight can be orders of magnitude brighter than matte terrain even with the weaker direct Sun. In this matched geometry the lunar peak is lower than Earth's. The area, number, duration and position of highlights determine the practical glare problem, and require water-slope statistics and an observer. Wave height alone does not supply the slope distribution. Cox–Munk is a method lead; its terrestrial wind-to-slope relation has not been applied to lunar seas.

The same separation between diffuse reflectance and specular peaks applies to wet paving, wet leaves, glass, metals and ice. Their directional reflectance belongs in the next scene model. The table cannot establish whether marine glare exceeds the discomfort produced by the diffuse surroundings.

## 7. Low Sun and the long light cycle

The surface-light product's diffuse transport is plane-parallel. Its existing comparison with the spherical atlas shows a large low-Sun deficit. Applying the **direct and diffuse ratios separately** from the unfiltered comparison to the shielded design gives this sensitivity at landscape albedo 0.1:

| Sun elevation | Moon stored two-stream | Moon atlas-transferred proxy | Earth atlas comparison |
|---|---:|---:|---:|
| 10° | 9.37 klux | 12.7 klux | 15.7 klux |
| 5° | 4.30 klux | 8.44 klux | 6.49 klux |
| 1° | 0.806 klux | 5.83 klux | 1.67 klux |

The lunar proxy assumes that the component correction transfers across the shield spectrum and profile differences. It is neither a new shielded spherical solution nor a correction for arbitrary albedo. Its angular sky remains unspecified. The atlas's open global energy residual of up to 4.6%, along with its prescribed exponential profiles, still applies. The low-Sun values demonstrate a modelling dependency rather than precision predictions.

For zero solar declination and a level equatorial horizon, each 0–10° passage lasts 19.69 hours on the Moon and 40 minutes on Earth; each 0–30° passage lasts 59.06 hours and 2 hours respectively. There is one rising and one setting passage per solar day. At 60° latitude, the 0–10° interval lasts about 40 hours per passage. Terrain, lunar obliquity, refraction and calendar geometry can change those intervals.

These are durations of geometrical conditions. They are not continuous human exposure durations or accumulated physiological doses. The study does not extend its stored daylight samples below the horizon.

## 8. What can be concluded about comfort

The measurements needed for discomfort include luminance distribution, source size and position, task, adaptation and observer response. CIE 252 distinguishes contrast-driven and total-amount effects. It also emphasizes each model's applicable conditions. The original DGP experiment concerns indoor office conditions; its simple illuminance fit was expected to cover approximately 1,000–10,000 lux. Even the positive illuminance term and constant in the full DGP formula exceed 1 at the lunar dark-ground noon value. Extrapolating that model here, or clipping its result to 100%, would manufacture an outdoor discomfort probability.

The first-pass findings are consequently specific:

- Ordinary matte horizontal landscapes receive less noon light than the matched Earth control. There is no quantitative basis here for classifying everyday Terluna daylight as intolerable.
- Under the imposed angular skies, a level eye and upright surfaces receive more ambient illumination on the Moon. Steep downward gaze reverses the comparison under the uniform sky. Pale ground supplies a large part of the light in both worlds.
- Solar-only shadows provide less luminous relief and weaker shading cues. Diffuse fill also exposes details that would otherwise lie in dark shadows.
- Direct Sun and ideal matched water highlights are weaker than Earth's. Specular peaks remain a separate local-glare mechanism, and low-Sun geometry persists much longer.
- A scene's surface reflectance, orientation, view of bright surroundings and selective shading matter at least as much as a horizontal lux total. Their joint effect is the main design variable established by this pass.

The upstream low UV index describes its protected atmospheric scenario. It supplies no visible-light discomfort threshold or comprehensive eye-safety assessment. This study therefore carries no retinal-hazard, clinical or photobiological clearance.

## 9. Directional sky and finite scenes

The [directional product](results/directional_scenes.json), schema `terluna.research.optical-comfort-directional-scenes/1`, adds 630 angular views and 24 scene cases, each with six probes. The runner is [directional.py](directional.py), with angular recovery in [angular.py](angular.py) and surface transport in [scenes.py](scenes.py). A second [figure generator](../../../visualization/optical-comfort/directional_scenes.py) displays the stored results.

### 9.1 Recovering the sky already in the repository

The original atlas archives are ignored by Git, but [immersion/assets/sky/atmosphere.json](../../../immersion/assets/sky/atmosphere.json) retains a packed copy. Its bake preserves photopic luminance while mapping signed RGB into positive RGB, then subsamples and quantizes to half precision. Each of 174 Sun positions has 41 elevation samples across the full sphere and 33 azimuth samples over the symmetric half-plane; 21 elevation samples cover the sky. The original atlas had 81 × 65 angular samples. The diffuse panorama excludes the directly viewed solar disk, whose beam is stored separately.

We recover luminance as `0.2126 R + 0.7152 G + 0.0722 B`. The study interpolates bilinearly in sine of view elevation and azimuth; it normalizes by the exact horizontal integral of that interpolant. This is an explicit scientific reading rule, rather than the renderer's texture interpolation or display exposure. Across the 14 selected Moon/Earth frames, the recovered horizontal sky flux differs from the stored full-resolution flux by at most 0.326%. This checks integrated fidelity, without bounding fine angular features. The packed beam and diffuse tables also agree with the earlier surface-light atlas comparison to within 0.00044% at its stored angles. Archive hashes differ; that agreement concerns the recorded summaries, rather than proof of byte-identical original archives.

The Moon uses `moon_no_ozone`, the atlas case paired with the upstream lunar comparison; Earth uses `earth`. These atlases have prescribed exponential optical profiles, unfiltered sunlight and landscape albedo 0.1. For each world and exact solar angle, we scale its angular sky shape to the **current** diffuse illuminance and use the **current** direct beam separately. The resulting field is a hybrid sensitivity. It does not resolve the effects of the shield spectrum or solved vertical profiles on sky direction. All new scenes keep the atmospheric boundary at albedo 0.1; the local pale ground never resets the whole column's albedo.

At noon, the atlas sky alone gives vertical/horizontal diffuse ratios of 0.529 on the Moon and 0.783 on Earth, against 0.5 for a uniform sky. Earth's brighter horizon therefore moderates the Moon/Earth level-eye comparison. With the current dark-ground fluxes:

| Gaze at noon | Moon eye illumination, klux | Earth, klux |
|---|---:|---:|
| Level | 21.26 | 11.59 |
| 30° down | 15.36 | 12.66 |
| 45° down | 12.78 | 12.86 |
| 60° down | 10.73 | 12.83 |
| Straight down | 8.96 | 12.46 |

At 30° Sun, a level view's sky contribution also varies with azimuth. Its ratio to horizontal diffuse illumination ranges from about 0.52 to 0.63 on the Moon and 0.85 to 1.05 on Earth across the sampled headings. These are disk-excluded sky contributions; the direct beam must still be added when it lies in front of the eye plane. The downward-eye reversal and a brighter lunar ambient field in some views both survive using computed angular structure.

### 9.2 Geometry and transport

The six scenes share flat, dark ground with reflectance 0.1 extending outside the local structures. Coordinates are metres, with z up and the Sun along +x. The observer is at `(0, 0, 1.6)`, facing the opening toward −y. This level view excludes the directly viewed Sun at both tested elevations. A second eye plane points at the work surface, 38.7° downward. Scene dimensions and every reflectance are stored in the product.

| Scene | Prescribed local geometry |
|---|---|
| Open dark ground | Reflectance 0.1 throughout |
| Open pale patch | An 8 × 8 m ground patch, reflectance 0.8, centred on the observer |
| Pale courtyard | The same patch, with 3 m-high walls of reflectance 0.8 at x = ±4 m and y = +4 m; open toward −y |
| Canopy, dark floor | A 6 × 6 m opaque roof at z = 3 m, reflectance 0.3 on both faces |
| Canopy, pale floor | The same roof over the 8 × 8 m pale patch |
| Screened canopy, pale floor | The preceding canopy plus 3 m-high walls of reflectance 0.3 at x = ±3 m and y = +3 m |

Every surface obeys the diffuse relation `L = rho E / pi`. Irradiance includes the visible solar beam and the cosine-weighted integral of incoming radiance. Ray intersections determine sky exposure and solar shadows. At reflected surfaces the same transport is evaluated again, so a roof changes the illumination of nearby ground, and that altered ground changes the roof, walls and observer's view. The surface rendering equation and path estimator follow the primary reference recorded in [sources.json](sources.json).

Cosine-weighted paths use two fixed scrambled Sobol sequences of 131,072 samples per probe, with up to 18 surface reflections. The collimated Sun is evaluated explicitly at each surface and at the initial plane; diffuse paths cannot count it twice. The product separates unobstructed sky, reflected sky, reflected sunlight and directly incident sunlight, and also separates visible ground from structures. No volume scattering or extinction is added between local surfaces. Materials are neutral, opaque and matte; weather, specular reflection and vegetation transmission remain outside these scenes.

### 9.3 What reaches the eye and the work surface

Noon values below are klux. The workplane is an infinitesimal horizontal measurement patch at `(0, −1, 0.8)`; it adds no table geometry or scene-wide reflectance change.

| Scene | Moon level eye | Earth level eye | Moon workplane | Earth workplane |
|---|---:|---:|---:|---:|
| Open dark ground | 21.26 | 11.59 | 89.60 | 124.60 |
| Open pale patch | 39.47 | 36.93 | 89.60 | 124.60 |
| Pale courtyard | 41.07 | 40.84 | 89.48 | 127.56 |
| Canopy, dark floor | 11.62 | 7.64 | 12.40 | 4.44 |
| Canopy, pale floor | 19.66 | 15.07 | 15.72 | 7.81 |
| Screened canopy, pale floor | 12.87 | 10.23 | 6.75 | 3.27 |

Pale local surroundings raise reflected light substantially on both worlds. The level-eye comparison becomes almost equal in the pale courtyard, despite the dimmer lunar open horizontal surface. That is different from applying albedo 0.8 to the entire atmospheric column in section 1: the finite patch changes local reflection while the incoming sky remains the albedo-0.1 boundary.

A roof gives a larger reduction at the workplane than at a level eye looking out. On the Moon at noon, the dark-floor canopy reduces workplane illumination from 89.6 to 12.4 klux, while level-eye illumination falls from 21.3 to 11.6 klux. With a pale floor, the eye receives 19.7 klux under that same roof: 8.43 from unobstructed sky, 9.49 from visible ground and 1.73 from the roof. The ground term now includes its computed shade. The three screens reduce that total to 12.9 klux. Shade effectiveness therefore depends on the observer's view as well as the illumination of the task.

At 30° Sun, the beam enters below the bare roof's edge. At the central workplane, its ray reaches the roof's height 3.81 m to the side, beyond the roof's 3 m half-width. The roof blocks the noon beam but leaves this beam visible. The side screen intercepts it.

| Canopy scene, Sun 30° | Moon workplane, klux | Earth workplane, klux |
|---|---:|---:|
| Canopy, dark floor | 23.86 | 55.93 |
| Canopy, pale floor | 26.38 | 60.58 |
| Screened canopy, pale floor | 4.13 | 2.13 |

The open workplane at that Sun angle receives 35.9 klux on the Moon and 57.8 klux on Earth. With the dark-floor canopy, Earth retains about 97% of that illumination because its direct beam dominates; the Moon retains about 67%. Roof depth, height, heading and side screening therefore deserve joint treatment for the long lunar low-Sun intervals. These particular dimensions are test scenes, rather than optimized buildings.

### 9.4 Two visual tasks

The first task is a neutral mark of reflectance 0.1 on a work surface of reflectance 0.2. They share the same illumination, so their physical Weber contrast remains −0.5 in every scene; their absolute luminances change with the workplane values above. Ocular scatter, gloss and adaptation would require additional modelling.

The second task is a raised block 2 m wide, 1 m deep and 0.15 m high, with reflectance 0.2 on its tread and riser. It sits from y = −2 to −1 m. The two surface probes sample the centres of the tread and the viewer-facing riser, with the block itself present for ray visibility. Other scene probes omit the block. The signed shading contrast is `(L_tread − L_riser)/(L_tread + L_riser)`; a negative value means a brighter riser. At noon:

| Scene | Moon tread–riser contrast | Earth tread–riser contrast |
|---|---:|---:|
| Open dark ground | 0.620 | 0.830 |
| Open pale patch | 0.280 | 0.390 |
| Pale courtyard | 0.270 | 0.326 |
| Canopy, dark floor | 0.127 | -0.007 |
| Canopy, pale floor | 0.009 | -0.048 |
| Screened canopy, pale floor | 0.682 | 0.501 |

Open-ground shading is weaker on the Moon, consistent with the first pass. Inside a shelter, the angular light reaching each face becomes decisive. Under the dark-floor roof the Moon retains a shading difference while the Earth pair is nearly equal; under the pale-floor roof both pairs are weakly differentiated. The screened shelter directs enough light toward the tread relative to the riser to restore a strong difference. At 30° Sun the courtyard even reverses the Earth tread–riser ordering. Thus a general claim that lunar diffuse light always reduces shape contrast would fail in these controlled scenes.

These values measure one shading cue. They do not predict whether a person will detect a step: edge boundaries, texture, angular size, binocular cues and movement remain available, and adaptation and ocular scatter are absent. The task contrast here is also distinct from the two-horizontal-patch solar-shadow contrast in section 3.

### 9.5 Numerical review and remaining uncertainty

The focused tests include exact normalization of a bilinear sky, uniform-sky tilted-plane integrals, packed flux reconstruction, explicit ray intersections, analytic open-ground illuminance and the point-to-square view factor for a black canopy. Additional checks cover an unilluminated roof underside, linearity in source strength, positive repeated reflection and both component ledgers. These tests assess the specified calculation.

Doubling the directional quadrature changes sky illumination by at most 0.000137 of the horizontal diffuse illuminance. Across the 144 scene probes, increasing the nested samples from 32,768 to 131,072 per scramble changes illuminance by at most 0.407%; the two full scrambled estimates differ by at most 0.318%. A preliminary 65,536-path, 12-reflection run and the retained run differ in step contrast by at most 0.00170. A conservative radiance estimate on paths still live after 18 reflections gives a largest sampled omitted-tail contribution of 0.00091% of probe illumination. Those sampling diagnostics are not confidence intervals or a certified whole-scene error bound.

The dominant remaining uncertainty is the physical illumination boundary: the packed proxy sky needs a consistent spherical calculation with the chosen shield and solved columns. Cloud fields, glossy/wet surfaces and measured material spectra would change the scene inputs. The current results identify mechanisms and useful design comparisons without establishing human comfort thresholds.

## 10. Next data and experiments

| Next input | What it resolves |
|---|---|
| Spherical spectral source fields using the current shield and solved atmospheric profiles | A consistent directional sky and low-Sun flux, replacing the angular-shape transfer |
| Original atlas archives and source metadata | Full angular precision and provenance; the packed copy is already sufficient for the present screening |
| Specific site/building geometry, surface spectra and BRDFs | Replace the six prescribed scenes with locally supported environments, including wet surfaces and vegetation |
| Water directional spectra or slope distributions, with wind and viewing geometry | Highlight extent and variability; a lunar wave-to-glint calculation |
| Selected GCM/CRM cloud fields or representative cloud optical cases | Cloudy and foggy scenes, then the frequency and duration of particular lighting conditions |
| Defined visual tasks and an outdoor human-factors validation protocol | How the computed scenes affect perceived comfort and task performance |

The next high-value local calculation is a consistent angular clear-sky boundary, followed by a small selection of weather and surface cases chosen around the mechanisms above. A higher-resolution copy of the old atlas improves angular precision but retains its physical assumptions.

## Checks

There are 34 passing focused checks across the two runners. The first-pass checks cover independent solid-angle integration, tilted eye planes, gaze-dependent reversal, source-angle guards, annular integration, Fresnel limits, solar residence, separation of local reflectance from landscape feedback, and reproduction and provenance of the stored scene product. Numerical reproduction verifies the stated model; empirical visual comfort remains open. [checks.json](checks.json) records the actual test counts and repository-wide results. The existing Python failures reproduce on the untouched starting commit; the full gate also encounters missing inputs, dependencies and runtime restrictions.
