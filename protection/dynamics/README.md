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

[active_global.py](active_global.py) defines every member of a recurring
meridional-plane population at a common ephemeris date, with finite-square ray
intersections, return paths, constrained ideal beam normals and a radial
separation bound. Its inverse-dynamics force bracket permits any mutual
illumination fraction; it does not replace that missing optical solution with
full sunlight. The [global study](../../research/studies/solar_shield_array/active_global.md)
records its rejected high cost and the completed focused handovers; the planned
rigid-target feedback replay was not run. `natural_traffic.py` releases the
initial geometry into gravitational IVPs with ideal photon-force cancellation
and independently checked interpolation to integer members through three days.
The seven-day interpolation fails its placement allowance. The checked prefix
develops optical gaps and independently replayed tile intersections. Gap repair,
collision-avoidance control, continuous coverage and bounded-slew attitude
realizability remain separate gates.

`coverage_pilot.py` supplies common-epoch near-natural seeds, finite-square
cell capacities, fractional covering allocation and conditional propulsion
accounting. The [bounded pilot](../../research/studies/solar_shield_array/pilot.md)
varies population across orbital families and initial-phase bands, then refines
a coverage bottleneck. Its 17.29-million-tile capacity estimate retains
subcell, mutual-illumination and individual-placement gaps; a selected
attitude-rate witness fails the pilot cap. It supplies a cost threshold for
further research and accepts no operating fleet.

`natural_pattern.py` adds smooth transported roll, exact sampled-source
rectangle-union mutual illumination, an all-pairs guard on the neighbour
reduction, finite-area gravity and moving-window ray unions. The
[local-pattern study](../../research/studies/solar_shield_array/patterns.md)
tests initial velocity gradients and freely deforming arrays. Its six-hour
361-tile passage preserves sampled local coverage without electric translation;
fleet placement, interfaces, return closure and attitude engineering remain open.
