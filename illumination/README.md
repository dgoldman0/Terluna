# Illumination and appearance

What light reaches the surface of an Open Moon and what its sky looks like: the
geometry of the Sun and Earth, spectral sky radiance and surface irradiance, and
(still to come) earthlight, starlight and airglow.

| Material | What it computes | Condition |
|---|---|---|
| [sky/](sky/) | Spherical, spectral, scalar multiple-scattering sky radiance and surface irradiance for Earth and two Open Moon optical profiles, Sun from −90° to +90° | Conditional on prescribed exponential optical profiles. Numerically checked (solver tests, atlas invariants, six noon Monte Carlo spot checks); a lunar global energy residual of up to 4.6% is open |
| [geometry.py](geometry.py) | Angular sweep with a 29.53-day period, idealized equatorial horizon and a six-degree interval | Simple angular arithmetic |

`geometry.py` is a verbatim copy of the [planning diagnostic](../ensemble/planning/twilight_diagnostic.py);
`research/check.py` verifies the copy byte-for-byte. Roughly 11.8 hours through
six degrees is not a brightness curve or a universal sunset duration.

```sh
python illumination/geometry.py
```

## Products other work consumes

The sky solver's atlases (`sky/data/*_atlas.npz`, generated, not committed) are
this domain's main product. The [month-of-light viewer](../visualization/month-of-light/)
displays them, and the immersion bakes them into its sky. Consumers read the
atlases; they do not import the solver.

## Next work

- Site geometry: positions of the Sun, Earth (with phase and libration) and stars
  for a given selenographic site; horizon obstruction and shadows.
- Earthlight as a light source. On the near side a nearly full Earth lights the
  lunar night; the current solver has the Sun as its only source.
- Tie the optical profiles to the atmosphere domain's solved column instead of the
  exponential proxies.
- Benchmark low-Sun and night radiance. The libRadtran documentation
  (https://www.libradtran.org/doku.php?id=basic_usage) is a method lead, not an
  installed dependency.

Original NASA/USNO source admission remains as recorded in the immutable
[planning reference manifest](../ensemble/planning/planning_references/manifest.json).
See [status](../research/status.json).
