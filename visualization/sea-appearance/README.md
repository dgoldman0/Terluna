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
