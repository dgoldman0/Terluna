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
