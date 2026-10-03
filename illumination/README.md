# Illumination and appearance

What light reaches the surface of an Open Moon and what its sky looks like: the
geometry of the Sun, Earth and stars, spectral sky radiance and surface irradiance,
and earthlight.

| Material | What it computes | Condition |
|---|---|---|
| [surface_light/](surface_light/) | Clear-sky direct and diffuse sunlight at the ground by wavelength (202–1000 nm): photosynthetic photons, colour, red:far-red, ultraviolet and ground brightness, for the design Moon, a 1.0-atm Moon and the Earth control | The atmosphere domain's radiation (1-nm ultraviolet-visible below 500 nm, line-by-line above) on its solved columns; two-stream, so light at Sun heights below about 20° is understated against the sky solver |
| [sky/](sky/) | Spherical, spectral, scalar multiple-scattering sky radiance and surface irradiance for Earth and two Open Moon optical profiles, Sun from −90° to +90° | Conditional on prescribed exponential optical profiles. Numerically checked (solver tests, atlas invariants, six noon Monte Carlo spot checks); a lunar global energy residual of up to 4.6% is open |
| [cloud_light/](cloud_light/) | Spherical molecular solar and observer paths to elevated features, using current column layers and an atlas-proxy control | Direct beam and atmospheric transmission only; cloud scattering and resolved twilight radiance remain open |
| [ephemeris.py](ephemeris.py) | Sun, Earth and star directions above a site; Earth's phase; earthlight as a fraction of sunlight | Mean-orbit geometry (synchronous rotation, lunar equator in the ecliptic, sinusoidal libration, Lambert-phase Earth of geometric albedo 0.367); good to a few degrees, not an ephemeris for dates |
| [stars/](stars/) | The Yale Bright Star Catalogue (9,096 stars): J2000 position, V magnitude, B−V | Catalogue data; the build pins the source file's hash |
| [geometry.py](geometry.py) | Angular sweep with a 29.53-day period, idealized equatorial horizon and a six-degree interval | Simple angular arithmetic |

`geometry.py` is a verbatim copy of the [planning diagnostic](../ensemble/planning/twilight_diagnostic.py);
`research/check.py` verifies the copy byte-for-byte. Roughly 11.8 hours through
six degrees is not a brightness curve or a universal sunset duration.

```sh
python illumination/geometry.py
```

## Products other work consumes

- The sky solver's atlases (`sky/data/*_atlas.npz`, generated, not committed). The
  [month-of-light viewer](../visualization/month-of-light/) displays them, and the
  immersion bakes them into its sky.
- The site-sky product (`ephemeris.product`, schema `terluna.illumination.site-sky/1`):
  the model's parameters for a site from [shared/scenarios/sites.json](../shared/scenarios/sites.json)
  plus golden samples. The immersion evaluates the same closed-form formulas for
  continuous time, and its tests check them against the samples.
- The bright-star catalogue ([stars/bright_stars.json](stars/bright_stars.json)).

Consumers read these products; they do not import the models.

The [optical-comfort study](../research/studies/optical_comfort/) consumes the
committed surface-light summaries to compare matte surface luminance, eye-plane
illumination, shadow contrast and ideal water glint with Earth. It keeps
landscape albedo separate from local surface reflectance and identifies the
remaining inputs for a resolved visual-comfort calculation. Its directional
extension recovers the packed atlas's angular patterns and transfers them onto
current light totals. Six finite matte scenes include cast shadows and repeated
reflection. A consistent shielded spherical sky remains open.

The [cloud-twilight study](../research/studies/cloud_twilight/) couples molecular
paths to corrected CM1 cloud heights and occurrence. It screens the long
illumination window, spectral filtering and sunsetward viewing geometry, with
the inherited sky atlas supplying separate ground-light context. Local cloud
microphysics and spatial fields are needed for a resolved appearance study.

For a near-side equatorial site 65° from the sub-Earth point, full Earth (earthlight
about 1.0 × 10⁻⁴ of sunlight above the atmosphere) comes near sunset, and the
midnight Earth is 72% lit. Through the reference Open Moon atmosphere that leaves
about 2 lux on the ground at midnight.

## Next work

- A date-accurate ephemeris, the lunar equator's 1.54° tilt, refraction, and
  horizon obstruction.
- Earthlight's own spectrum. Earthlit sky glow currently reuses the solar atlas at
  Earth's elevation, which treats earthlight as sunlight-coloured; the solver could
  compute it with Earth's reflectance spectrum as the source.
- Tie the optical profiles to the atmosphere domain's solved column instead of the
  exponential proxies.
- Benchmark low-Sun and night radiance. The libRadtran documentation
  (https://www.libradtran.org/doku.php?id=basic_usage) is a method lead, not an
  installed dependency.

Original NASA/USNO source admission remains as recorded in the immutable
[planning reference manifest](../ensemble/planning/planning_references/manifest.json).
See [status](../research/status.json).

Local follow-up on `study/optical-comfort` now checks native sky resolution, recovers second-cycle dusk microphysics and derives resolved coastal water slopes. [The study](../research/studies/optical_comfort/README.md#10-follow-up-with-local-data) records the remote baseline’s data access and the new local products.
