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

`pattern_return.py` supplies gravity variational shooting maps, capped smooth
arcs, sparse endpoint/encounter programs and between-sample coplanarity
witnesses. `parallel_shadow.py` supports dynamic blocker lists and two-sided
finite-Sun momentum for arbitrary shared planes in uneclipsed local prefixes.
The [return checkpoint](../../research/studies/solar_shield_array/returns.md)
retains geometrically rejected endpoint proposals and independently repeats
an early clearance failure during the coast attitude transition. It does not
accept a complete handover, orbital return or recurring fleet budget.

`pattern_departure.py` and `departure_control.py` add temporary layer and
transverse preparation, smooth edge turns and convex impulse corrections in
a small reusable control family. The [departure study](../../research/studies/solar_shield_array/departure.md)
verifies a six-hour departure/coast from the actual service terminal states,
including coupled shadows, interval refinement and an independent replay.
`attitude_load.py` accounts for rigid-square Euler torque, explicit zero-net-
force electric couples and gravity/radiation disturbance envelopes. This
local ideal-actuator result leaves structural realization, hardware-mass
feedback, handover, return and the fleet budget unclosed.

`collection.py` adds per-tile first-intercept optical power, the redirected
spectral band and time-integrated energy. It shares the parallel-pattern
shadow geometry and checks its inferred photon momentum against the force
model. The [collection ledger](../../research/studies/solar_shield_array/collection.md)
applies it to the saved service and departure states and records the required
collection account for future population/control searches. The uneclipsed
prefix rejects body occultations; joint body/array visibility and specified
collector optics are required for broader trajectories and electrical yield.

`energy_coast.py` supplies smooth partial-turn schedules; `coast_control.py`
permits signed impulses in the existing three pattern modes. `tile_torque.py`
computes finite-area gravity torque and first moments of the same parallel
shadow unions used by the force model. The [energy continuation](../../research/studies/solar_shield_array/energy.md)
retains force-model failures, local control limits and the distinction between
a partial passage's actuation account and recurring fleet power.

`packing_control.py` varies the four starting depth levels and provides sparse
individual two-burn corrections, with measured force defects and explicit
nonlinear replay gates. `square_distance.py` checks global physical surface
separation and conditional interval bounds; `service_margin.py` evaluates
arbitrary solar-disk source positions against a complete receiver window.
The [packing pilot](../../research/studies/solar_shield_array/packing.md) tests
these geometry and energy freedoms while retaining collection accounting.

`cycle_control.py` continues the measured states into a smooth return-to-service
attitude and permits four independent acceleration arcs per tile. Its endpoint
LP provides a restricted-model lower bound and explicit failed encounter
branches. `mass_feedback.py` propagates additional uniformly distributed mass
with unchanged optical area and a per-member power allocation. The
[cycle-closure study](../../research/studies/solar_shield_array/closure.md)
tests return, mass and power constraints before expanding fleet coverage.
`compact_shadow.py` removes empty projected rectangles before evaluating the
same exact shadow union. All-pair geometry and original-force comparisons
check that this compaction preserves the optical and momentum calculation.

## Return correction and collection accounting

`return_escape.py` supplies smooth in-plane roll/recovery commands and exact partial-burn integrals. `compact_collection.py` preserves the original first-intercept and photon-momentum ledger with exact compact shadow unions, allowing denser grazing-source quadrature. The [encounter-directed study](../../research/studies/solar_shield_array/return_repair.md) retains transition failures, the migrated encounter and individual hardware overloads.

## Joint departure and arrival

`fast_parallel.py` and its compiled backend evaluate all-blocker rectangle
unions and first moments, cross-checked against the historical implementation.
`eclipse_parallel.py` adds joint body/array visibility on the same incoming
rays, with refined area quadrature in partial eclipse. `oriented_tiles.py`
provides independently oriented shadow cones, force/moment accounting and
finite-square contact, tested against parallel limits and independent ray
intersections. `joint_transfer.py` permits individual endpoint times and
power-capped electric control arcs in a reduced linear proposal model.

The [joint study](../../research/studies/solar_shield_array/joint_cycle.md)
distinguishes a full geometric return plus another service passage from its
failed actuator/supply and energy gates. Independent replay preserves all
361 physical tiles. A grouped-attitude proposal fails its reduced geometry
screen and receives no coupled feasibility credit. Per-member chronological
storage dispatch lives in `engineering/shield_storage.py`; added collectors
and storage are priced sensitivities, not accepted loaded trajectories.

## Service-led natural returns

`service_search.py` proposes natural service dates and charges their necessary
cadence inventory. `service_reduced.py` supplies cheap centre-force proposals,
explicitly omitting mutual shadows. `service_constraints.py` and
`service_arrival_controls.py` fit free terminal states to assigned optical
rays with early/arrival control arcs, individual ratings and energy ceilings;
the latter charges partial early burns exactly. `service_holes.py` adds
constraints from full receiver polygon differences. `oriented_clearance.py`
prunes all-pair contact tests with conservative projection bounds and agrees
with a brute-force oriented-square scan.

The [bounded continuation](../../research/studies/solar_shield_array/natural_return.md)
finds a 30.22% saving for the loaded twelve-hour passage with a slower turn,
but no operating return. Coupled dynamics reject the cheap sparse-ray proposal
and reproduce compression to the clearance floor in the slower-turn return.
Omitted shadow response cannot be treated as a small correction in these fits.

## Relative orbits

[relative_orbit.py](relative_orbit.py) describes tiles relative to a reference
orbit about the Moon. It gives osculating semi-major axes and periods, the
along-track drift per revolution of unequal orbits (about 3πΔa), and the burns
along the velocity that equalize energy. It also gives the quasi-nonsingular
relative orbital elements and the linear relative motion they define. The
closed-form minimum radial-normal separation of energy-matched relative orbits
follows D'Amico and Montenbruck (2006). The tests compare the linear solution
with Hill's equations and two-body propagation, and the separation with brute
force and the published form. The
[relative-orbit plan](../../research/studies/solar_shield_array/relative_orbits.md)
uses these tools to design formations whose relative orbits repeat and stay
separated, and verifies them in the coupled model.

## Ring bundles and keeping

[ring_bundle.py](ring_bundle.py) propagates nested rings of tiles that each
hold their own radial-facing attitude. Each tile's attitude is referenced to
its velocity and tilted about the orbit normal, so that neighbours on a ring
overlap like shingles. The sail force uses the compiled finite-Sun shadow and
eclipse geometry of `eclipse_parallel` and `fast_parallel` in the bundle
centroid's frame. Incidence and force direction use each tile's own normal.

`ContinuousRings` stands for complete rings with one propagated representative
per ring. Its neighbours on the ring and on the adjacent rings enter as copies
shifted by whole tile lags with an exact universal-variable two-body solution
(`kepler_shift`). Every tile is then shadowed as on a ring without ends. In
120 hours the tide moves real neighbours along their shared path by up to about
150 m per lag relative to those copies.

[ring_keeping.py](ring_keeping.py) holds rings to designed relative elements.
`geometric_elements` measures them exactly from vis-viva, eccentricity vectors
and unit normals. `plan` gives the near-circular impulsive two-burn and normal
correction. `ContinuousKeeping` applies low-thrust feedback that holds each
ring to its neighbour toward a free reference ring. The neighbour is compared at
the same phase after `gravity_shift` carries it along its orbit under the tide.
Its eccentricity and plane targets lie on axes fixed by an origin direction.
Under Earth's tide a whole stack of rings regresses together, so on those axes
every pair's plane difference turns with the stack, and holding it there works
against the common regression (see the solar shield array study's relative-orbit
plan).
