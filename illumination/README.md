# Illumination and appearance

What light reaches the surface of an Open Moon and what its sky looks like: the
geometry of the Sun, Earth and stars, spectral sky radiance and surface irradiance,
and earthlight.

| Material | What it computes | Condition |
|---|---|---|
| [calendar/](calendar/) | Date-based topocentric Sun/Earth directions, a finite phase-dependent spectral Earth disk, and horizontal surface light through the full day and night | The static [light calendar](../visualization/light-calendar/) reads its versioned products; the date geometry is sampled against JPL across 2000–2500 |
| [surface_light/](surface_light/) | Clear-sky direct and diffuse sunlight at the ground by wavelength (202–1000 nm): photosynthetic photons, colour, red:far-red, ultraviolet and ground brightness, for the design Moon, a 1.0-atm Moon and the Earth control | The atmosphere domain's radiation (1-nm ultraviolet-visible below 500 nm, line-by-line above) on its solved columns; two-stream, so light at Sun heights below about 20° is understated against the sky solver |
| [sky/](sky/) | Spherical spectral sky radiance: established exponential-profile atlases, plus current solved-column skies for the shielded design Moon and Earth | The new product includes spectral and spatial checks, independent Monte Carlo and a global-energy ledger. [aerial.py](sky/aerial.py) runs the same integration from an observer's height and stops it at chosen distances: the light the air adds in front of a surface and the transmittance to it. Clouds, aerosols, polarization and refraction remain additional inputs |
| [cloud_light/](cloud_light/) | Regional evening cloud radiance, self-shadowing and multiple scattering on saved CM1 fields, with a fixed-observer history, dark viewing openings and faint-twilight context | Explicit 2-D cloud extrusions and particle optics; regional 3-D structure, terrain, Earth-source spectra and perceptual response remain further inputs |
| [water_surface/](water_surface/) | Reflection by rough water: the resolved waves' slope covariance from SWAN's spectra (printed spectra and restart files), the short waves from the unified spectrum of Elfouhaily et al. (1997), whose dispersion and capillary scales carry gravity and surface tension, and Fresnel facets with Gaussian slopes and Smith's shadowing for the glitter of the Sun and the Earth and the reflected sky; explicit sea surfaces drawn from those spectra on tiles from 2 km to 1.4 m, with Lagrangian (Gerstner) crests on the gravity waves and filtered-slope pyramids for rendering ([realization.py](water_surface/realization.py)) | Linear waves; the short-wave laws are Earth fits, checked at Earth's gravity against Cox and Munk's clean-sea slopes and evaluated at lunar gravity with the lunar seas' u*/c_m inside their fitted range |
| [water_column/](water_column/) | The colour of seawater: absorption and scattering of pure seawater, phytoplankton (Bricaud et al. 1998), dissolved organic matter and suspended regolith fines (Apollo soils through Hapke's model), the remote-sensing reflectance and the diffuse attenuation of daylight | Earth's ocean optics and its nadir reflectance model (Lee et al. 2002); fetched tables ([fetch_inputs.py](water_column/fetch_inputs.py)); the fines' optics rest on Hapke's isotropic and equivalent-slab approximations |
| [earthlight/](earthlight/) | The Earth's light at the Moon by wavelength and phase: Glenar et al.'s (2019) model spectrum of the whole Earth at the visual brightness Robinson et al. (2025) fit to observations (geometric albedo 0.242), the curve ephemeris.py and the immersion also follow | Spring-season model spectrum without the Earth's daily turn; observed visual phase curve from 5° to 144°, a Lambert sphere's shape beyond |
| [ephemeris.py](ephemeris.py) | Sun, Earth and star directions above a site; Earth's phase; earthlight as a fraction of sunlight | Mean-orbit geometry (synchronous rotation, lunar equator in the ecliptic, sinusoidal libration, the Earth's measured visual phase curve with sunlight's colour); good to a few degrees, not an ephemeris for dates |
| [stars/](stars/) | The Yale Bright Star Catalogue (9,096 stars): J2000 position, V magnitude, B−V; their places in the Moon's sky at a date and the naked-eye limit against a bright sky ([sky.py](stars/sky.py)) | Catalogue data; the build pins the source file's hash. Places by precession in longitude and the optical libration, without proper motion, parallax or aberration |
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
reflection. The [solved-column follow-up](sky/SOLVED_COLUMN.md) now supplies absolute angular radiance and flux from the same shielded spherical atmosphere, with updated gaze and finite-scene results.

The [cloud-twilight study](../research/studies/cloud_twilight/) now couples
236 second-cycle snapshots on each of three corrected CM1 rings and a highland
box to spectral cloud transport. It measures regional occurrence, selected
cloud colour and contrast, and a fixed observer through 8–80 hours after
sunset. A long-range case at 49 hours tests a dark foreground beneath a warm
horizon view. The products retain source hashes, photon-block uncertainty,
cloud-width and phase-function sensitivities, and a separate Earthlight proxy.

For a near-side equatorial site 65° from the sub-Earth point, full Earth (earthlight
about 6.3 × 10⁻⁵ of sunlight above the atmosphere) comes near sunset, and the
midnight Earth is 71% lit. Through the reference Open Moon atmosphere that leaves
about 0.7 lux on the ground at midnight.

## Next work

- Move remaining mean-orbit consumers to the [date calendar](calendar/) geometry;
  add refraction and horizon obstruction as distinct inputs.
- Earthlight's own spectrum in the older products. The earthlit sky glow of the
  atlas and ephemeris.py's earthlight have the Earth's measured brightness with
  sunlight's colour; [earthlight/](earthlight/) gives the spectrum.
- Extend the cloud scenes to regional profiles, three-dimensional cloud fields,
  terrain and Earth-source light; update the remaining consumers of the older illumination products.
- Benchmark low-Sun and night radiance. The libRadtran documentation
  (https://www.libradtran.org/doku.php?id=basic_usage) is a method lead, not an
  installed dependency.

Original NASA/USNO source admission remains as recorded in the immutable
[planning reference manifest](../ensemble/planning/planning_references/manifest.json).
See [status](../research/status.json).

Local follow-up on `study/optical-comfort` now checks native sky resolution, recovers second-cycle dusk microphysics and derives resolved coastal water slopes. [The study](../research/studies/optical_comfort/README.md#10-follow-up-with-local-data) records the remote baseline’s data access and the new local products.
