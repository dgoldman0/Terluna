# Illumination and appearance

[geometry.py](geometry.py) is a verbatim convenience copy of the existing [planning diagnostic](../ensemble/planning/twilight_diagnostic.py). It calculates an angular sweep using a 29.53-day period, an idealized equatorial horizon, a six-degree interval and a nominal solar diameter.

```sh
python illumination/geometry.py
```

**Condition:** simple angular arithmetic. Approximately 11.8 hours through six degrees is not a brightness curve or a universal sunset duration. A solar disk crossing is a different interval. No lunar twilight spectrum, sky color, multiple scattering, terrain horizon or ephemeris calculation has been implemented here.

## Next work

Compute location/altitude-dependent Sun-Earth geometry, horizon obstruction and shadows. Then apply a stated atmosphere and spectrum with appropriate spherical transmission/scattering and refraction. Benchmark low-Sun approximations before using them to describe the sky. Atmospheric and cultural interpretations belong to separate tests.

The libRadtran documentation at https://www.libradtran.org/doku.php?id=basic_usage is a method lead from the prior audit, not an installed dependency. Original NASA/USNO source admission remains as recorded in the immutable [planning reference manifest](../ensemble/planning/planning_references/manifest.json).

The stand-alone copy is checked against the planning original; deliberate later development should replace this duplication with an explicit reviewed implementation. See [status](../research/status.json).
