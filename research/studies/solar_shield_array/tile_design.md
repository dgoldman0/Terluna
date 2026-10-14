# Tile hardware, failures and replacement

14 October 2026 · `research/array-industry`

**A 10 km tile held by a rigid perimeter frame weighs about 200 times its allocation.** With edge stiffness,
three-axis momentum storage, electricity and the inertia they add, the reference annulus tile comes to
**2.684 kg/m²** and the window tile to **5.870 kg/m²**, against complete-tile allocations of 13.765 and 26 g/m². Its
100 m frame also misses the 150 m clearance between shingles. The result belongs to the rigid straight-edge
architecture under the assumptions below; supports with curved edges, smaller units and spinning membranes are compared
in [support_comparison.md](support_comparison.md).

The calculation is [tile_design.py](tile_design.py). Its versioned product,
[results/tile_design.json](results/tile_design.json), records the model and input hashes, the shared constants, the
assumptions, the numerical refinement and the software versions, and [tile_sources.json](tile_sources.json) holds the
primary sources it was checked against. The earlier products are unchanged.

## What was designed and what was tested

The engineering domain holds the [tile budget](../../../engineering/array_tile.py), the
[angular-momentum calculation](../../../engineering/array_momentum.py) and the
[replacement queues](../../../engineering/array_service.py); protection holds the
[local failed-tile dynamics](../../../protection/dynamics/tile_failure.py). This study couples them to the branch's
rolled attitude, optical films, radius layout and industrial inventory.

It is a first engineering assessment: explicit component sizing and rigid-body requirements, the membrane's shape and
frequency, local fault propagation and a count of the whole population's service. A flexible spacecraft following the
command, an allocated photon controller, removal and insertion trajectories and 28 million propagated tiles are outside
it.

The attitude test flies one circular lunar orbit at 20,000 km, with a period of about 70.5 h, under the branch's 80°
roll outside service and its cosine ramps: the largest angular speed is 1.066×10⁻⁴ rad/s and the largest angular
acceleration 1.538×10⁻⁸ rad/s². The attitude product the branch inherits flew the kept orbits near 15,000 km and stores
summary statistics, so the new hardware is tested on a time history of its own.

Two inherited inputs disagree with the hardware:

- The window's optical stack alone weighs **26.921 g/m²**, above the 26 g/m² complete-tile allocation: the match of one
  sail loading took the window tile at its stack's mass. The annulus film weighs 10.211 g/m².
- The attitude study's `store()` integrates torque components in the moving body axes and leaves out the transport
  term ω × h. Its 1.498×10¹⁰ N·m·s for the rolled law is the rotor sizing's starting allocation and floor here; the
  conserved calculation below adds the transport term and leaves out the radiation torque.

The annulus allocation is the glow product's 13.765 g/m² (the older industry estimate was 13.07 g/m²). Membrane
pressure takes that product's dynamic coefficients, 0.14612 for the annulus and 0.276 for the window, the dynamics'
convention that carries all of the window's turned-away light as reflected; the bare window stack's coefficient, about
0.139, is recorded beside it and gives a lighter load.

## The reference hardware

The tile is a 10 km square optical stack on a 1 g/m² load-bearing substrate, counted separately, under uniform biaxial
tension and held by a perimeter space truss. Each side has four chords, X braces, transverse posts and a 20% allowance
for joints. The reference truss is 100 m deep, with 50 m bays and chord tubes of 2 m radius. A further 20% covers mass
growth across the whole ledger, as a margin over the components designed here.

The substrate is taken at 100 MPa admissible stress and 1,400 kg/m³, and the frame at 100 GPa, 500 MPa admissible
stress and 1,600 kg/m³: assumed design values for a composite, set against the sources' fibre, laminate and typical
data-sheet values. Qualified allowables, load transfer and strain compatibility with the silica stack, and the added
substrate's spectrum and ageing are open.

For a square of side L under uniform normal pressure q and membrane tension T, the fixed-edge sine series gives a centre
sag of about 0.07367 qL²/T and a first membrane frequency of √(T/σ)/(√2 L). The pressure envelope is the solar
pressure plus σ(|α| + ω²)L/√2 with the local trim mass, a conservative bound that takes the rotational loading at its
largest magnitude and gives no credit for translational acceleration. On the rolled command the rotational part is 79%
of the annulus envelope. Tension is sized for 5 m of deflection with a factor of two for an ideal step response, which
covers a single step; repeated resonant forcing needs the flexible model. The reference first modes have periods of
about 533–542 s, the starting point for a coupled flexible controller.

**Edge stiffness sets the frame.** Each side carries the line load T, with a largest bending moment of TL²/8. The
chords are sized for the assumed strength and for Euler buckling and then for bending stiffness, and the braces for
shear stiffness, with half the allowed 5 m inward bow of the edge given to bending and half to shear. With four chords
the bending inertia is A_chord d², so the ideal bow is 5TL⁴/(384 E A_chord d²). Sized for member strength alone the
edges would bow far beyond the straight film edges the sag solution assumes. Joint and post compliance, local shell
buckling, whole-frame modes, deployment, fatigue and impact tolerance are open.

The momentum store has three independently reversible axes, two hoops to an axis attached to the perimeter, each axis
with the full allocated capacity. The reference hoops are 100 m in radius with a rim speed of 1,000 m/s: H = mrv,
energy mv²/2 and hoop stress ρv², so an assumed 1,600 kg/m³ rotor runs at 1.6 GPa against an assumed 2 GPa allowable.
A 50% allowance covers supports, bearings and motors; their detailed design and burst containment are open. The six
stations and the hoops' own inertia enter the full locked inertia tensor, its off-diagonal term included.

The power ledger holds motor electronics, generation, eclipse storage, 10 kV distribution conductors at 2% loss, heat
rejection, avionics and sensors. Storage covers the largest three-axis swing of spin energy in eclipse with no
regenerative credit: conservative power envelopes, ahead of an electrical schedule. Optical trim covers 25% of the
tile, with 2 g/m² of hardware on that area. Its target, a bistable device not yet qualified, is 5 J/m² an update, 20
updates an orbit and no holding power: 2.5 GJ per 10 km tile per orbit, about 9.85 kW on average. The driver and
storage also cover all of an orbit's updates falling in eclipse, an envelope in which the nominal update windows
overlap; the schedule of updates is a later design step. The NASA variable-reflectance source is a proposal for such a
device.

## Momentum and mass, solved together

The solver conserves the total angular momentum in inertial space:

`H_inertial(t) = H_inertial(0) + integral(R(t) torque_external_body(t) dt)`

`h_rotor_body(t) = R(t).T H_inertial(t) - I_locked omega_body(t)`

Motor torque changes the body components of the rotor momentum, and the bearings also carry the gyroscopic term ω × h.
Gravity gradient is the external torque. The calculation starts with no relative rotor momentum and applies no
radiation torque and no optical unloading: it measures what the attitude and gravity gradient alone ask of the store,
from a zero bias.

Rotor and motor allocations rise to 1.2 times the largest per-axis requirements; their masses change the inertia, and
the calculation repeats until the allocations agree within 0.1%, never falling below the inherited starting
allocation. This closes the rigid sizing loop. The rotor momentum ends the orbit away from where it began, at about
237 GN·m·s for the reference annulus and 519 GN·m·s for the window, so periodic operation needs external unloading.
The frame is sized for the membrane's support; the motor attachments' strength and the whole frame's out-of-plane slew
loads are open.

| Case | Annulus, g/m² | Window, g/m² | Meaning |
|---|---:|---:|---|
| Current complete-tile allocation | 13.765 | 26.000 | Inherited flight assumptions |
| Reference, 100 m frame depth | 2,683.5 | 5,870.1 | 5 m film sag and 5 m edge bow |
| 50 m frame depth | 10,625.9 | 23,231.6 | A shallower frame costs much more stiffness mass |
| 20 m permitted film sag | 684.6 | 1,500.5 | Lower tension; the optical shape to recheck |
| 50 m permitted edge bow | 284.8 | 626.6 | Relaxed straight edge |
| 400 m rotor radius | 1,973.5 | 4,317.4 | Lighter store; larger swept hardware |
| 3 km tile | 278.9 | 586.4 | Size sensitivity on the present formation |
| 1 km tile | 86.8 | 158.4 | Size sensitivity on the present formation |

The smaller tiles scale the frame geometry, sag and bow with the side, keep the housekeeping and avionics of a tile, and
start the storage and torque at L³ times the reference's, a radiation-torque-like allocation (inertia at fixed areal
mass goes as L⁴); the assembled-inertia calculation is then applied to each. At 1 km the rotors hold that starting
allocation, about eight times what the conserved calculation asks. Smaller tiles need their own count, coverage and
formation.

For the reference annulus tile the frame is about 150 million kg of the 268 million kg total, and the resized
allocation is about 7.44×10¹¹ N·m·s and 8.97×10⁷ N·m of motor torque on each axis; for the window tile it is about
1.63×10¹² N·m·s and 1.96×10⁸ N·m. At that mass the sail loading falls about 200-fold, and the orbit, steering and glow
would all have to be flown again.

Carried to the 28.233 million tiles of the present inventory, these tiles would weigh about **9,915 Gt**, against the
46 Gt of the optical inventory, and replacing whole tiles at the service policy below would move about
**760 Gt a year**, 24 million kg/s: the cost of this architecture, beside the film plant's own flows. Repairing
components and reusing frames is a different service model.

## Clearance and local failures

Parallel tiles at the 8.5 km pitch along the track and the 1.5° shingle stand 222.50 m apart along their normals before
any hardware thickness. A 50 m envelope leaves 172.50 m and a 100 m envelope 122.50 m, under the inherited 150 m gate:
the reference frame, whose side members lie over their neighbours' across the 1.5 km of shingle overlap, fails the gate
from the first epoch. Rotors 100 or 400 m in radius need their own swept geometry, which a whole-tile box overstates.

The fault calculation follows the failed tile's displacement from its healthy slot with the Hill–Clohessy–Wiltshire
equations of a circular lunar orbit. The neighbours hold their slots, and the residual is a constant difference of
acceleration in radial, along-track and orbit-normal coordinates. Finite parallel boxes and continuous event roots
locate the loss of clearance and the first contact, where the physical account ends.

| Envelope depth | Radial residual, m/s² | First loss of 150 m gate | First envelope contact |
|---|---:|---:|---:|
| 1 m | 10⁻⁵ | 63.01 min | 110.92 min |
| 50 m | 10⁻⁶ | 111.80 min | 311.14 min |
| 50 m | 10⁻⁵ | 35.35 min | 97.87 min |
| 50 m | 10⁻⁴ | 11.18 min | 30.95 min |
| 100 m | 10⁻⁵ | At the start | 82.47 min |

The accelerations are sensitivities spanning what a failed tile might be left with. Changing attitude, shadows, the
eccentric orbit, tides, debris and the neighbours' controllers are outside the experiment.

For the 50 m envelope and 10⁻⁵ m/s², an assumed independent backup acting both ways, capped at 10⁻⁴ m/s² and switched
on after 30 minutes, keeps the sampled clearance above 154.33 m and returns the tile to its slot; switched on after 60
minutes it avoids contact, after the safety gate has already been crossed. A tile back in its slot still has a failed
filter. The backup is an assumed actuator: as photon hardware, a tug or an independent thruster it would carry mass and
power of its own.

Lost attitude is faster. For the 50 m case, a conservative bound on surface displacement loses its 150 m guarantee after
about **182 s** at an angular error rate of 0.001°/s, or 18 s at 0.01°/s: the times at which separation stops being
guaranteed, which come before any contact.

## Replacement at the whole population

The population count follows all 28,233,094 slots for one year: uniformly staggered planned renewal every 15 years,
0.01 unexpected failures per working tile-year, dispatch capacity 20% above planned renewal, 27 days of planned transit
and one-day replacement of failed slots from local spares. Reliability, life and transit are declared scenarios.

That asks for about **5,153 planned exchanges and 773 failures a day** against 6,184 dispatches a day, 4.35% spare
capacity over all expected jobs. A planned exchange leaves the old tile working until its replacement arrives, while a
failure opens a slot at once. Jobs take finite transit, failed slots go first, and the queues conserve arrivals,
dispatches, completions and outstanding work. Failures leave already planned jobs in place, a conservative count.
Manufacturing and spare stocks are taken to keep up with dispatch.

| Scenario | Mean missing slots, fraction | Peak missing slots, fraction | Renewal backlog at year end |
|---|---:|---:|---:|
| Local spares, 1 day | 2.73×10⁻⁵ | 3.06×10⁻⁵ | 0 |
| Delivery from plant, 27 days | 7.11×10⁻⁴ | 7.48×10⁻⁴ | 0 |
| Capacity only 80% of planned renewal | 2.73×10⁻⁵ | 3.06×10⁻⁵ | 658,369 |
| Dispatch stopped for 14 days | 3.56×10⁻⁵ | 4.07×10⁻⁴ | 23,050 |
| One correlated loss of 10,000 slots | 2.91×10⁻⁵ | 3.84×10⁻⁴ | 0 |
| Failure rate 0.1 per tile-year | 9.83×10⁻³ | 1.91×10⁻² | 1,882,206 |

In the capacity-shortage case the missing slots stay few while the renewal debt grows, since overdue tiles keep the
constant failure rate. The local-spares case has about 140,000 replacements in transit at its peak, mostly planned
work, so local buffering carries an inventory and logistics of its own.

On an ideal periodic grid of 8,900 × 9,890 m optical rectangles at a pitch of 8,500 × 9,286.456 m, about 89.1% of the
ray area has a single covering tile. For independent outages of fraction q the added uncovered fraction is Σ w_m q^m
over the coverage multiplicities m: about 2.43×10⁻⁵ as a daily sampled mean with one-day spares, and 6.33×10⁻⁴ with
27-day delivery. These are the **ideal geometric gaps the outages add**. The seams, pinholes, oblique light and curved
geometry of the real stack, and the protection gates, belong to the shield's UV transmission count; the correlated
case, which breaks independence, has no such estimate.

Every service result takes the failed hardware as already contained. The contact windows above are minutes to hours,
so a one-day spare policy answers for the slot and containment for the failed tile. The counts and the physical
containment are kept apart: no failed tile is removed or replaced inside a trajectory.

## Design consequence and remaining work

The next hardware candidate relaxes the 10 km straight-edge frame: smaller independent units and tension-supported or
curved membranes with explicit load paths and optical shape ([support_comparison.md](support_comparison.md)). Wider
hoops lighten the store, leave the frame's stiffness where it is and enlarge the swept geometry. The window tile's
allocation needs its structure added to its 26.921 g/m² stack.

A candidate feeds an updated orbital layout and the net electricity gate once it has a coupled flexible-body and
controller calculation with the actual force and shadow history, the spectral properties of its complete stack, finite
component clearance through the rolled cycle, and independent fault containment. Manufacturing, deployed joints,
ageing, rotor containment, removal and insertion routes, and repair and reuse are open engineering tasks.

## Reproduction and checks

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.tile_design
python -m pytest engineering/test_array_tile.py engineering/test_array_momentum.py engineering/test_array_service.py protection/dynamics/test_tile_failure.py research/studies/solar_shield_array/test_tile_design.py
make check
```

The study takes about 6 CPU seconds and 90 MiB on one numerical thread; it was regenerated on 14 October 2026 after
the source check corrected the wording of [tile_sources.json](tile_sources.json), with every value unchanged (the first
run took about 7 CPU seconds). Its 50 focused tests pass:
the membrane series and stiffness, mass and energy accounting, conservation of angular momentum and frame invariance,
analytical local orbits and contact between samples, delayed recovery, queue conservation, optical multiplicity,
provenance and the study's integration. Doubling the attitude resolution from 1,441 to 2,881 samples changes the peak
momentum by 0.0101%, and cutting the local maximum step from 30 to 10 s moves the reported contact and clearance roots
by less than 0.1 s. These are numerical checks; empirical validation, physical qualification and a full reading of
the sources are separate steps.

`make check` passed with **2,053 Python tests passed and 16 skipped**, and all **191 JavaScript and calendar tests**
(129 immersion, 48 column and transport, 14 calendar). The layer check found no violations, and historical provenance,
protection verification and ensemble integrity passed. The full check ran on the final model, source, product and test
snapshot immediately before the commit; [tile_checks.json](tile_checks.json) records the commands, the counts and the
exact final model and test hashes.
