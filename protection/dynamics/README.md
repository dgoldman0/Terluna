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
python -m pip install -r research/studies/solar_shield_array/requirements.txt
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

## Cycling membrane experiment

[cycling.py](cycling.py) adds an ideal specular sail, Earth/Moon eclipse unions,
retarded shadow-service geometry, circular-problem retrograde seeds and
service/feathering control. It is separate from the unchanged held-aperture
producer. The [cycling study](../../research/studies/solar_shield_array/cycling.md)
contains the collocation search and independent full-ephemeris propagation.
A 50-g/m² candidate remains bounded for three years; continuous fleet coverage
and practical attitude/formation control are not demonstrated.

## Simultaneous finite tiles

[fleet.py](fleet.py) adds independent common-epoch initial states, ideal safe
reflection, finite-square perspective shadows over the solar disk, seams and
swept collision/mutual-shadow guards. [fleet_exclusions.py](fleet_exclusions.py)
provides a vectorized equivalent of the reference pair scan for larger
populations. The [fleet study](../../research/studies/solar_shield_array/fleet.md)
records five 30-day populations, spatial maps, local minute-cadence coverage,
inventory exclusions and control bounds at 50 g/m². Continuous coverage remains
unachieved, and an independent replay differs by 898 m against a 50 m placement
allowance. Large-square stress probes expose significant finite-extent force
errors; neither practical formation control nor delivered power is established.

[active_formation.py](active_formation.py) adds bounded electric feedback for
separate 10 km tiles, finite-Sun mutual illumination and swept finite-square
separation. [active_retime.py](active_retime.py) integrates finite-extent gravity
and smooth thrust arcs for individual phase transfers. The
[active study](../../research/studies/solar_shield_array/active.md) records local
gap acquisition and independently replayed 1,400 km orbital retiming. Its local
coverage and control budgets leave the global fleet assignment unresolved.
