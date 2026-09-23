# Illumination and appearance

What light reaches the surface of an Open Moon and what its sky looks like: the
geometry of the Sun, Earth and stars, spectral sky radiance and surface irradiance,
and earthlight.

| Material | What it computes | Condition |
|---|---|---|
| [sky/](sky/) | Spherical, spectral, scalar multiple-scattering sky radiance and surface irradiance for Earth and two Open Moon optical profiles, Sun from −90° to +90° | Conditional on prescribed exponential optical profiles. Numerically checked (solver tests, atlas invariants, six noon Monte Carlo spot checks); a lunar global energy residual of up to 4.6% is open |
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
