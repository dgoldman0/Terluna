# Simultaneous finite-tile fleet checkpoint

The Boolean service-time repair is committed and pushed as `2182374`.
This follow-up initializes every fleet member separately on 2026-10-04 TDB
and integrates all members through the same 30 days of DE440 forcing.
The populations use osculating lunar circular initial guesses in one or
several planes. No stored trajectory is shifted in time to manufacture a fleet.

The first product is [fleet_checkpoint.json](results/fleet_checkpoint.json).
It measures perspective shadow unions of finite 10 × 10 km squares across
16 directions on the finite solar disk every six hours. It evaluates the
four-lunar-radius protected disk and the solid lunar disk separately. The
intersection of those shadow unions measures area covered for every sampled
solar direction. Overlapping shadows count once. A 10 m seam plus a 50 m
placement allowance at each edge leaves a 9,890 m clear square.

| Initial population | Mean protected-region ray interception | Collision-guard conflicts | Retained after collision and mutual-shadow exclusions |
|---|---:|---:|---:|
| 576 members, one plane, 15,000 km | 0.005428% | 5,208 pairs | 61 |
| 576 members, 12 planes, 15,000 km | 0.002450% | 618 pairs | 48 |

Both initial populations remain between about 14,460 and 15,637 km over the
month. Neither provides full-Sun protection at any target area on the evaluated
solar quadrature. These small populations are method and family probes;
their physical area is already insufficient to cover the protected disk.
The planar concentration raises its average interception while leaving most
of the aperture empty. Plane diversity alone does not supply enough material
or coordinate the handovers.

The ray calculation uses finite square vertices, the finite Sun, perspective,
first-order interception-time tile positions, and Earth eclipses. The finite
Sun spreads an isolated 10 km tile's shadow over a roughly 140 km diameter at
15,000 km distance; its full-Sun umbra vanishes. Optical coverage therefore
requires coordinated neighbouring tiles over the entire solar disk.

Collision guards use circumscribed spheres around the squares, a 100 m
clearance and swept closest approaches between all 600-second outputs. They
include an acceleration allowance of 0.15 m/s², exceeding the measured maximum
of 0.024 m/s². A guard conflict is a conservative exclusion, not a claim that
the material surfaces necessarily collided. A second swept projected guard
excludes possible mutual solar shadows. Whole trajectories are removed by a
greedy graph selection. This is one feasible conservative subset selection;
it is not an optimum or a population limit.

Independent sail propagation assumes unshadowed incident light apart from
Earth and Moon eclipses. Consequently, a nominal population with mutual
shadows is a geometric probe whose interacting optical forces remain to be
solved. The strict subsets exclude those interactions under the stated bounds.
Their mean interceptions are only 0.000594% (61 members) and 0.000218%
(48 members). Their base optical inventories are 3.05e8 and 2.40e8 kg.

The ideal service/coast controller steers reflected beams away from Earth
and the protected Moon. Its guard includes the finite solar cone, tile corners,
100 km above Earth's surface, target motion during light travel and a 1 mrad
normal-pointing allowance. The sampled Earth margin is at least 0.0573° after
those allowances; the lunar margin exceeds 4°. The force uses the inherited
central-ray spectral actuator at 50 g/m². Finite-Sun force integration,
angular material response, flexible-tile dynamics and slew-limited attitude
control remain further work. The sampled plane commands reach about 0.150°/s
in the multi-plane pilot; these are lower bounds on actuator demand.

No electric correction is added. Collectors, habitat mass, torque hardware,
power transfer and delivered electricity have no demonstrated budget in these
runs. The architecture's 100 TW remains a separate design target.

## Reproduce

```sh
python -m pip install -r research/studies/solar_shield_array/requirements.txt
python -m protection.dynamics.ephemeris --download
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name planar_10km --planes 1 --phases 576 --arrangement planar
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name planes_10km --planes 12 --phases 48
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_analyze --names planar_10km planes_10km
python -m pytest protection/dynamics research/studies/solar_shield_array -q
```

The propagated raw histories are hash-checked and reproducible from the stored
configuration and deterministic initial-state generator. All initial phases
use the same time-dependent ephemeris. Numerical/source checks establish the
implemented calculation; continuous protection and an operational fleet remain
unachieved. Further families, denser geometry probes and quadrature refinement
are the next search steps after this checkpoint.
