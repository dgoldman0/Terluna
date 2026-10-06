# Sea appearance

Figures of the [sea appearance study](../../research/studies/sea_appearance/README.md),
written as ignored files with provenance sidecars to `results/`.
`python visualization/sea-appearance/lighting.py` reads its lighting calendar and
map:

- `results/lighting_calendar.png`: the clear-sky light on the ground hour by hour
  through 7 February to 8 March 2038 at six coasts, from the sky lit by the Sun,
  from the Earth and from a 0.001-lux placeholder for stars and airglow, with
  their total, the modes of vision for an 18% grey surface, and the Sun's and the
  Earth's elevations.
- `results/earth_over_seas.png`: every water cell at 1-degree spacing coloured by
  the share of 2026–2045 the Earth's centre spends above its horizon, on the
  atlas's relief, with every water body outlined. The cividis ramp runs from dark,
  where the Earth never rises, to bright, where it always stands.

`python visualization/sea-appearance/slopes.py` reads the sea-slope product:

- `results/slope_spectra.png`: the curvature spectrum of the unified spectrum of
  Elfouhaily et al. (1997) at lunar and Earth gravity under the same 10 m winds,
  3 and 5 m/s, with the capillary scales of both worlds and the wavelength below
  which the study takes the waves from the wave runs.
- `results/sea_slopes.png`: the mean square slope through the month at the six
  coasts: the short waves hour by hour, all waves where the wave runs give
  spectra, Earth's clean sea under the same wind, and the hours too calm for short
  waves.

`python visualization/sea-appearance/waters.py` reads the water-colour product:

- `results/water_colours.png`: the remote-sensing reflectance of the biosphere's
  design-guess waters, and swatches of the light leaving each under the Moon's and
  Earth's clear daylight with the Sun 45° up, at one exposure on a display with
  Earth-daylight (D65) white, with their luminance and the depth where
  photosynthetic light falls to 1%.

`python visualization/sea-appearance/regimes.py` reads the regime product and its
panoramas on the research drive:

- `results/regimes_<coast>.png`, one per rendering coast: six strips through the
  month (the Sun at its highest, 4° up, 4°, 15° and 35° below the horizon, and
  the darkest hour), each the whole horizon 30° above and below from an eye 2 m
  above the water, with the sea in every direction. Colours are as calculated on
  a display with Earth-daylight (D65) white, unadapted; each strip has its own
  exposure, printed as the luminance shown white, and brighter highlights clip.

## Image-generation references

[image_guides.py](image_guides.py) prepares the Ingenii midday frame (IC-1) for
image generation from its stored scene product and LOLA coast grid. It projects
the terrain and interpolates the scene's display colours. Its terrain lighting
is a normal-based approximation, without cast shadows or interreflection. The
optional water guide calls the domain's wave realization on the scene's borrowed
Nectaris spectrum, projects its height and filtered slopes, and uses approximate
Fresnel contrast to show the wave structure. These are geometry, colour and
structure references; the spectral path tracer remains the radiance renderer.

```sh
OPENBLAS_NUM_THREADS=1 python visualization/sea-appearance/image_guides.py --water-structure
```

The guides and intermediate component images live in ignored
`results/illustrations/`. At the author's request, selected final images are also
preserved in Git under [illustrations/](illustrations/README.md). JSON sidecars
keep source hashes, prompts, generation passes and review notes. Images themselves contain no text or labels. Generated
terrain detail and wave texture are illustrative; review their framing, lighting,
shadows, atmospheric depth and wave scales against the inputs. Vegetation is
omitted for this optical study; the bare shore is a placeholder for land whose
ecology remains to be supplied by the domain research.

[wp6_guides.py](wp6_guides.py) prepares WP-6, the earthlit darkest hour at
western Procellarum. It projects the local LOLA coast, draws Earth's stored phase
and apparent size, and reconstructs the mean glitter path with the domain's
uniform-disk reflection model. Its optional wave guide uses the local SWAN
spectrum, without adding short wind waves. An instantaneous wave realization
need not show the mean glitter path; use the colour and structure guides together.
The guide does not solve full spectral transport, cast shadows or long-exposure
time integration. Its sidecars state these limits and compare the wave moments.

```sh
OPENBLAS_NUM_THREADS=1 python visualization/sea-appearance/wp6_guides.py --water-structure
```


## Nubium midnight and human vision

The supplemental [scene product](../../research/studies/sea_appearance/results/nubium-midnight.json)
puts nearly full Earth 53.25° above southern Nubium at the hourly minimum of the
Sun's elevation. A 65° by 87.4° portrait includes both Earth and the sea. Earth's
1.804° diameter becomes about 64 by 80 pixels in a 2048 by 3072 rectilinear frame;
the off-axis projection stretches its vertical dimension. A thumbnail compresses
this entire angular window. At matched scale Earth is about 3.67 times the Moon's
diameter at equal distance, regardless of the canvas resolution.

[nubium_reference.py](nubium_reference.py) supplies this supplemental scene to
the existing [sea reference tracer](../reference-renderer/seas/README.md). It uses
the full wave cascades, filtered facets, wave occlusion and interreflection,
terrain lighting and spherical-atmosphere transport. It saves XYZ radiance,
sampling errors, material masks and source hashes. The camera is offshore, 2 m
above the instantaneous long-wave surface; a standing platform is not modelled.

[nubium_human_view.py](nubium_human_view.py) projects the NASA maps with JPL
Horizons orientation and normalizes their integrated RGB-derived XYZ to the
calculated direct Earth illuminance. This replaces only the direct uniform disk;
water reflection still uses the tracer's uniform-Earth approximation. The
[external-input manifest](nubium-earth-inputs.json) records exact URLs, hashes,
credits, rights review and the saved Horizons query. Maps are historical texture
priors, not predicted March 2038 weather or spatially resolved spectral radiance.
Raw inputs stay in ignored `results/earth-inputs/`; mismatched hashes fail.

The common physical-radiance image is passed through the official Radiance
`pcond -h+` implementation: contrast sensitivity, mesopic colour loss, acuity and
veiling glare. The declared display is sRGB/D65, 100 cd/m² peak, 1000:1 contrast;
a 100:1 sensitivity view is also saved. A separate Earth-fixation sensitivity and
10° Earth-centred window show why one screen frame cannot represent looking
around and adapting. The main wide-view operator chooses a linear mapping and
clips the disk's bright detail. **That clipping does not establish that a person
would see a featureless Earth.** Display, gaze, observer and adaptation history
remain uncertain. The narrower view cannot include a horizon 53° below Earth.

The method is Ward Larson, Rushmeier and Piatko (1997),
[A Visibility Matching Tone Reproduction Operator for High Dynamic Range Scenes](https://graphics.cs.yale.edu/sites/default/files/1997tvcg_hdr_tonemapping.pdf),
implemented by [Radiance](https://github.com/LBNL-ETA/Radiance) at commit
`bcffc2b52d99b908adfdc2543e4cddf83268a125`; see the
[pcond manual](https://radsite.lbl.gov/radiance/man_html/pcond.1.html).
These are display-model assumptions, not a new domain model or empirical
validation of an Open Moon. Source reading and numerical checks do not replace
perceptual validation.

```sh
# Use an environment with numpy, scipy, Pillow and numba; restrict worker threads.
NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 <environment>/bin/python visualization/sea-appearance/nubium_reference.py --out visualization/sea-appearance/results/nubium-reference --width 2048 --spp-side 2
NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 <environment>/bin/python visualization/sea-appearance/nubium_human_view.py --render visualization/sea-appearance/results/nubium-reference/s_nubium__midnight.npz --build visualization/sea-appearance/results/nubium-reference/scene-build.npz --assets visualization/sea-appearance/results/earth-inputs --radiance-bin <Radiance-build>/bin --out visualization/sea-appearance/results/nubium-human
python visualization/sea-appearance/nubium_scale.py --out visualization/sea-appearance/results/nubium-human/nubium-scale.html
```

Build the official Radiance `pvalue`, `pfilt` and `pcond` executables from that
revision with a headless CMake build. Restore the three manifest inputs to the
named asset directory before running. The wrapper caches a completed scene build
and resumes rendering by row blocks. Keep a separate output directory for each
resolution/sample configuration, and never change sampling while partial row
blocks exist: the upstream tracer's partial-block reader does not enforce that
identity itself.

[nubium_scale.py](nubium_scale.py) writes a standalone interactive size comparison.
Its equal-field panels show geometric disks only, without labels inside them or
claims about the supplied Moon photograph's unknown lens, crop and glare.
[nubium_guides.py](nubium_guides.py) retains the earlier quick composition/wave
references and separately exposed Earth experiment. Those photographic guides
are preparatory work, not the accepted naked-eye treatment.

Selected images and their evidence/review records live in
[illustrations/](illustrations/README.md). The AI passes were reviewed for
composition, size, lighting, reflection distribution and waves. The generator
changed the geometry and returned 1024 by 1536 even when 2048 by 3072 was
requested; those outputs are retained as illustration candidates, not certified
physical views. The larger calculated rendering preserves the specified geometry.

The three focused view checks run with
`TERLUNA_RADIANCE_BIN=<Radiance-build>/bin python -m pytest visualization/sea-appearance/test_nubium_view.py`.
Without that environment variable the external Radiance conversion check skips;
geometry and midnight selection still run. The existing sea-kernel tests require
numba and were run separately from `make check`.


## Visibility-informed Nubium reference

The revised [wide view](illustrations/nubium-midnight-visibility-view.png) and
[Earth-centred view](illustrations/nubium-earth-visibility-view.png) retain
planetary structure supported by the [visibility screening](../../research/studies/sea_appearance/README.md#can-an-observer-resolve-earth).
The old `nubium-midnight-human-view.png` is retained for provenance and is
**superseded as an appearance reference**: its linear display mapping clipped
the bright disk and its coarse glare grid produced an angular halo.

[nubium_visibility_view.py](nubium_visibility_view.py) consumes the same
2048 × 3072 sea-tracer XYZ, replaces the direct Earth using observed EPIC
structure, and retains the solved beam's integrated XYZ. The historical Earth
is geographically reprojected using the saved Horizons orientation. Its fit to
EPIC geolocation has 0.47 native-pixel RMS and 0.90-pixel 99th-percentile error.
About 2.2% of target disk samples lie beyond the measured view and explicitly use
the nearest observed surface direction. A Lambert incidence ratio approximates
the phase change. Three narrow bands supply spatial RGB structure, normalized
to the solved integrated XYZ; this is **not spatial spectral colour validation**
or predicted weather. The sea still reflects the tracer's uniform Earth disk.

[nubium_eye_optics.m](nubium_eye_optics.m) calls the external model's CIE99
optical transfer function on direct Earth XYZ for age 24. It applies the same
achromatic kernel to each channel. The tiny FFT undershoot removed is recorded
(5.3e-14 of summed XYZ energy). [nubium_visibility_display.py](nubium_visibility_display.py)
projects that optical result into the wide camera by angular rays. It tapers the
halo from 4 to 5 degrees; the maximum computed halo there is 0.01191 cd/m²,
about 2.3% of the local sky. Background and water remain unblurred. This is an
explicit optical/display approximation, not a whole-retina simulation.

One monotone log-power luminance curve is used for both views and every pixel.
It maps the local 0.529 cd/m² sky to 3.5% of display white and uses a world-white
anchor 20% above the brightest source pixel. Gamut compression moves colours
toward neutral at fixed display luminance. **Those anchors are display choices,
not measured colour appearance.** There is no separate exposure for Earth.
The images contain no text or labels. The 10-degree square shows an angular
window toward Earth; the sea horizon is 53 degrees away and cannot appear in it.

[nubium_display_visibility.py](nubium_display_visibility.py) tests the revised
physical Earth and the displayed PNG against interior-texture removal. For the
PNG it assumes sRGB/D65, 100 cd/m² peak, 0.1 cd/m² black and a viewing distance
that preserves the stated angular field. A detectable difference establishes
retained detail, not equality of appearance, neural colour response or a match
to an arbitrary user's screen. The
[display record](illustrations/nubium-visibility-display.json) carries its
input hashes, constraints and check results.

```sh
# Reuse the existing physical sea build; restore the pinned EPIC input.
python visualization/sea-appearance/nubium_visibility_view.py --root "$PWD" --epic research/runs/sea_appearance/visibility-inputs/epic_1b_20151117002712_00.h5 --render visualization/sea-appearance/results/nubium-reference/s_nubium__midnight.npz --build visualization/sea-appearance/results/nubium-reference/scene-build.npz --assets visualization/sea-appearance/results/earth-inputs --out visualization/sea-appearance/results/nubium-visibility
# Supply optics-config.json with hdrvdp, image_package_list, compat, input,
# output and age=24. compat points to the study's octave_compat directory;
# input/output are nubium-earth-optics-input.mat / nubium-earth-optics.mat.
octave-cli --no-history --no-init-file --quiet visualization/sea-appearance/nubium_eye_optics.m optics-config.json
python visualization/sea-appearance/nubium_visibility_display.py --root "$PWD" --directory visualization/sea-appearance/results/nubium-visibility --build visualization/sea-appearance/results/nubium-reference/scene-build.npz
python visualization/sea-appearance/nubium_display_visibility.py --root "$PWD" --directory visualization/sea-appearance/results/nubium-visibility --backend-template research/runs/sea_appearance/visibility/backend.json --out visualization/sea-appearance/results/nubium-display-visibility
octave-cli --no-history --no-init-file --quiet research/studies/sea_appearance/visibility_backend.m visualization/sea-appearance/results/nubium-display-visibility/backend.json
```

The source and display checks run with
`python -m pytest research/studies/sea_appearance/test_visibility.py visualization/sea-appearance/test_visibility_view.py`.
External data and software are optional at test collection and required only for
the full study/display runs. The renderer still needs the existing solved sky,
numba, SciPy and Pillow, with h5py added for EPIC.

## Nobili lake

[nobili_reference.py](nobili_reference.py) supplies the admitted
[Nobili scene](../../research/studies/sea_appearance/results/nobili-night.json)
to the existing spectral terrain/water tracer. Its 50-degree, 4:3 view preserves
Earth above the measured western rim. The clear reference is 2560 by 1920 at
four samples/pixel; its lake realization has 0.177 m significant height. The
conditional cloud reference reflects the calculated angular cloud sky and
attenuated Earth, so sharp glints receive the beam attenuation separately from
diffuse light. Terrain diffuse band shapes, observer-based cloud attenuation
and clear-reference aerial perspective remain approximations. It is not a
fully coupled terrain/cloud transport solution.

[nobili_view.py](nobili_view.py) projects the existing NASA Blue Marble maps
with the new JPL Earth orientation. Historical cloud texture is illustrative.
Low-elevation Earth lies outside sRGB: the new positive-XYZ normalization
preserves its integrated physical beam without creating negative spatial
radiance. The existing RGB path remains available for earlier scenes. A global
monotone display curve and the CIE99 age-24 ocular kernel retain disk structure;
neural colour adaptation and exact naked-eye appearance are not validated.
The measured FFT negative-energy fraction is recorded and bounded below 1e-7;
this numerical undershoot is zeroed before display. No image contains labels.

```sh
# An environment with numba, scipy and Pillow; limit worker threads.
python visualization/sea-appearance/nobili_reference.py --out research/runs/sea_appearance/nobili/reference-full --width 2560 --spp-side 2
python visualization/sea-appearance/nobili_reference.py --out research/runs/sea_appearance/nobili/cloud-reference-full --width 2560 --spp-side 2 --cloud research/runs/sea_appearance/nobili/cloud-light.npz
python visualization/sea-appearance/nobili_view.py --reference research/runs/sea_appearance/nobili/reference-full --out research/runs/sea_appearance/nobili/clear-view
# Run nubium_eye_optics.m with the emitted optics-config.json, then repeat with --finish.
```

The optical pass needs the previously admitted HDR-VDP 3.0.7 environment.
For the cloud display pass, supply `--cloud` and use `--display-reference`
pointing to the clear display JSON to keep the same exposure convention.
Selected references, the first generated image, generation prompts and the
review record are preserved in [illustrations](illustrations/README.md).

Nobili uses the new `illumination/cloud_light/finite_source.py` source drivers
and `illumination/water_surface/full_spectrum.py` equilibrium spectrum. The
older solar-cloud driver and short-wave modules remain byte-for-byte intact
for existing published products. The Nobili run records retain their original
producer hashes from before this function-preserving packaging change.

## Additional cloud views

The [two cloud scenes](../../research/studies/sea_appearance/README.md#two-additional-cloud-views)
reuse the Nobili adapters with explicit scene, wave and ground-light products.
Fecunditatis uses a 50° portrait camera, Earth 42.645° high and a 1.937° disk.
Eastern Smythii uses a 90° landscape camera, with both Sun and Earth below the
horizon. Its legacy `__night` output suffix is a renderer identifier: the
selected scene is bright twilight. The offshore cameras have eyes 2 m above
water; no standing platform is established.

[nobili_reference.py](nobili_reference.py) reflects the conditional cloudy sky,
attenuates the direct Earth separately and uses the refined ground flux for
water-leaving light. The outer wave tile expands to at least eight peak
wavelengths, preventing the 33-m Smythii sea from losing its long-wave energy.
Its sampled 1.169-m significant height is close to the spectrum's 1.150-m value.
[coastal_cloud_guide.py](coastal_cloud_guide.py) traces the saved cloud geometry
independently; its grayscale is opacity, not brightness. Upper-cloud guides
exclude material below 20 km. They constrain the generated cloud footprint,
not resolved filaments or three-dimensional towers.

[cloud_display.py](cloud_display.py) filters the cloudy-minus-clear XYZ field on
a regular angular grid with a declared 4° Gaussian sigma. This suppresses photon
colour noise and sacrifices angular detail; it is a display approximation, not
an extra transport result. Raw groups remain in ignored research storage.
The independent flux checks use those raw groups, before smoothing.

The night image keeps the earlier Nobili global log-power display convention.
The twilight image uses one linear photographic exposure, mapping the clear
upper sky to 20% of display white, with a highlight shoulder beginning at 80%.
The earlier steep log-power daylight treatment was rejected. These operators
are not human colour-adaptation models and the two pictures cannot be used to
compare absolute brightness. The reference products carry the lux values.
The CIE99 age-24 kernel acts on the direct Earth only. An explicit horizon test
prevents painting a below-horizon Earth onto Smythii's sea; its zero direct
field bypasses the optics calculation.

Rebuild the scene/light/wave products using the study instructions, then:

```sh
STUDY=research/studies/sea_appearance
RUN=research/runs/sea_appearance/fecunditatis-earth-clouds
python visualization/sea-appearance/nobili_reference.py --scene "$STUDY/results/fecunditatis-earth-clouds.json" --waves "$STUDY/results/fecunditatis-cloud-waves.json" --ground "$STUDY/results/fecunditatis-cloud-ground-light.json" --cloud "$RUN/light/cloud-light.npz" --cloud-smoothing-deg 4 --width 2048 --spp-side 2 --out "$RUN/filtered-reference"
python visualization/sea-appearance/coastal_cloud_guide.py --scene "$STUDY/results/fecunditatis-earth-clouds.json" --width 768 --out "$RUN/structure"
python visualization/sea-appearance/nobili_view.py --scene "$STUDY/results/fecunditatis-earth-clouds.json" --prefix fecunditatis --reference "$RUN/filtered-reference" --cloud "$RUN/light/cloud-light.npz" --display-reference visualization/sea-appearance/illustrations/fecunditatis-earth-clouds-clear-display.json --out "$RUN/filtered-view"
# Run the existing HDR-VDP/Octave environment on the emitted optics-config.json:
octave-cli --no-history --no-init-file --quiet visualization/sea-appearance/nubium_eye_optics.m "$RUN/filtered-view/optics-config.json"
# Repeat the nobili_view command with --finish.
```

For Smythii, substitute its scene/wave/ground filenames, prefix `smythii`, width
2560 and the saved `smythii-twilight-clouds-clear-display.json`. That record
selects the linear exposure automatically. A first view pass writes its zero
optics product, so the Octave call is unnecessary. To recreate the clear
controls, omit cloud/ground arguments, render at width 1024 (Fecunditatis) or
1280 (Smythii), then use `--sky-reference` with Nobili's clear display for the
night view, or `--linear-display --sky-display .2` for twilight. External
NASA texture maps and the previously admitted optics environment are required.

Repeat the image-space checks from saved, committed references:

```sh
python visualization/sea-appearance/coastal_image_review.py --reference visualization/sea-appearance/illustrations/fecunditatis-earth-clouds-reference.png --display visualization/sea-appearance/illustrations/fecunditatis-earth-clouds-display.json --image visualization/sea-appearance/illustrations/fecunditatis-earth-clouds.png --out /tmp/fecunditatis-image-checks.json
```

Substitute the Smythii prefix for the other picture. These checks cover image
geometry, Earth detail ranks and broad displayed light levels, not radiometric
or perceptual validation. Horizon edges must also be inspected: bright wave
crests fooled the automatic detector in some early Earth passes. Selected
pictures, first passes, prompts and the remaining discrepancies are in the
[illustration record](illustrations/README.md#two-cloud-illustrations).
