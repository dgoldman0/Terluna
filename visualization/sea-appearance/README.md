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
