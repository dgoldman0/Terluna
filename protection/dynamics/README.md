# Ephemeris-based shield dynamics

The October 2026 [array study](../../research/studies/solar_shield_array/report.md)
adds a DE440s point-mass spacecraft model while preserving the historical
`protection/model.py` and its imported results.

`ephemeris.py` reads the pinned NASA/JPL kernel, supplies barycentric geometric
states and builds bounded Hermite interpolation for propagation. Epoch strings
are TDB calendar dates. It has no synthetic or circular fallback.
`model.py` computes inverse forces for actual-phase trajectories and independent
area nodes. `optical.py` supplies finite-Sun eclipse geometry, retarded ray
vectors, conservative target-motion coverage and the relaxed photon-momentum
control ball. Passive control has no sunward force component.

Restore the external input from NASA; its binary stays outside Git:

```sh
python -m pip install -r research/requirements.txt
python -m protection.dynamics.ephemeris --download
python -m pytest protection/dynamics research/studies/solar_shield_array
```

[ephemeris.json](ephemeris.json) records the source, attribution, size, published
MD5 and SHA-256. A mismatch stops retrieval. The short kernel covers 1849–2150;
use a separately pinned full DE440 input for later mission epochs.

The spacecraft force model uses shared constants and treats planetary systems
as point masses. It omits multipoles, detailed general relativity, solar-wind
forces, flexible structures and optical device physics. The array study records
the point-mass acceleration residual and numerical convergence. The ideal
optical ball is a lower bound on electrical thrust, with achievable materials,
exhaust paths and Earth-safe outgoing rays still to supply.
