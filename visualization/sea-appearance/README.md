# Sea appearance

Figures of the [sea appearance study](../../research/studies/sea_appearance/README.md).
`python visualization/sea-appearance/lighting.py` reads its lighting calendar and
map and writes ignored figures with provenance sidecars to `results/`:

- `results/lighting_calendar.png`: the clear-sky light on the ground hour by hour
  through 7 February to 8 March 2038 at six coasts, from the sky lit by the Sun,
  from the Earth and from a 0.001-lux placeholder for stars and airglow, with
  their total, the modes of vision for an 18% grey surface, and the Sun's and the
  Earth's elevations.
- `results/earth_over_seas.png`: every water cell at 1-degree spacing coloured by
  the share of 2026–2045 the Earth's centre spends above its horizon, on the
  atlas's relief, with every water body outlined. The cividis ramp runs from dark,
  where the Earth never rises, to bright, where it always stands.
