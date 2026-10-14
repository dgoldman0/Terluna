# Tile hardware, failures and replacement

14 October 2026 · `research/array-industry`

**The explicit 10 km perimeter-frame design does not meet the present fleet's
mass and clearance allowances.** Once edge stiffness, three-axis momentum
storage, electricity and their added inertia are included, the reference
annulus tile reaches **2.684 kg/m²** and the window tile **5.870 kg/m²**. Their
current complete-tile allocations are 13.765 and 26 g/m². This rejects the
specified rigid-edge architecture under the assumptions below; it is not a
lower bound on every possible sail or a finding that the shield is impossible.

The completed calculation is [tile_design.py](tile_design.py), its versioned
product is [results/tile_design.json](results/tile_design.json), and the focused
primary-source register is [tile_sources.json](tile_sources.json). The product
records its model and input hashes, shared constants, assumptions, numerical
refinement and software versions. No earlier product has been rewritten.

## What was designed and what was tested

The engineering domain owns the [tile budget](../../../engineering/array_tile.py),
[angular-momentum calculation](../../../engineering/array_momentum.py) and
[replacement queues](../../../engineering/array_service.py). Protection owns
the [local failed-tile dynamics](../../../protection/dynamics/tile_failure.py).
This study couples their results to the branch's rolled attitude, optical
films, radius layout and industrial inventory.

This is a preliminary engineering assessment: explicit component sizing and
rigid-body requirements, a membrane shape/frequency calculation, local fault
propagation and whole-population count simulation. It does **not** demonstrate
a flexible spacecraft following the command, an allocated photon controller,
a removal/insertion trajectory, or 28 million independently propagated tiles.
No candidate is recorded as a validated design.

The new attitude test is one **circular 20,000 km lunar orbit**, using the
branch's 80° roll outside service with its existing cosine ramps. The orbital
period is about 70.5 h. The maximum angular speed is 1.066×10⁻⁴ rad/s and
angular acceleration 1.538×10⁻⁸ rad/s². The inherited stored attitude product
used a different, perturbed orbit near 15,000 km. Its scalar summaries are not
a time history that can be replayed as the new hardware.

Two existing inconsistencies are now explicit:

- The window optical stack alone is **26.921 g/m²**, already above the 26 g/m²
  complete-tile allocation. The annulus film is 10.211 g/m².
- The earlier `attitude_schemes.store()` integrated torque components in moving
  body axes without the transport term. Its 1.498×10¹⁰ N·m·s rolled-law result
  is retained only as an initial hardware allocation. It is not a validated
  rotor requirement.

The annulus allocation uses glow's updated 13.765 g/m², rather than the older
13.07 g/m² industry estimate. Membrane pressure uses that product's dynamic
coefficients, 0.14612 and 0.276, to retain the redirector force convention.
Bare-film coefficients are also recorded; substituting their approximately
0.139 window value would understate the chosen dynamical load.

## The reference hardware

The tile is a 10 km square optical stack on a separately charged 1 g/m²
load-bearing substrate, with uniform biaxial tension and a perimeter space
truss. Each of the four sides has four chords, X braces, transverse posts and
a 20% joint mass allowance. Reference truss depth is 100 m, bays are 50 m and
chord tube radius is 2 m. A further 20% mass-growth allowance covers the complete
ledger; it is a margin, not a substitute for designing missing components.

The substrate target is 100 MPa admissible stress at 1,400 kg/m³. Frame targets
are 100 GPa modulus, 500 MPa admissible stress and 1,600 kg/m³. These are
assumed composite/substrate design values, **not qualified material allowables**.
The sources distinguish fiber properties, laminate properties and manufacturer
typical data. Load transfer and strain compatibility with the silica stack
remain unproved, as do the extra substrate's optical spectrum and aging.

For a square of side L under uniform normal pressure q and membrane tension T,
the fixed-edge sine-series solution gives centre sag approximately
0.07367 qL²/T and first membrane frequency √(T/σ)/(√2 L). The pressure envelope
includes solar pressure plus σ(|α|+ω²)L/√2, with local trim mass included.
This conservative bound discards the direction of rotational loading and any
cancellation by translational acceleration. Tension is sized for 5 m deflection
with a factor of two for an ideal step response. That factor does not bound
repeated resonance. Reference first-mode periods are approximately 533–542 s;
these frequencies are not a coupled flexible-control stability demonstration.

**Edge stiffness dominates.** Each side carries line load T, with maximum
bending moment TL²/8. Chord area meets assumed strength and Euler buckling,
then also bending stiffness; braces meet shear stiffness. Half the 5 m allowed
inward edge bow is assigned to bending and half to shear. With four chords,
the beam bending inertia is A_chord d² and its ideal bow is
5TL⁴/(384 E A_chord d²). Sizing only for member strength would conceal large
edge displacement and would not support the assumed straight film edges.
Joint/post compliance, local shell buckling, complete-frame modes, deployment,
fatigue and impact tolerance still need analysis.

The momentum store has three independently reversible axes, with two hoops
per axis attached to the perimeter. Each axis receives the full allocated
capacity. Reference hoop radius is 100 m and rim speed 1,000 m/s: H=mrv,
energy=mv²/2 and hoop stress=ρv². The assumed 1,600 kg/m³ rotor has 1.6 GPa
hoop stress against an assumed 2 GPa allowable. A 50% mechanical allowance
charges supports, bearings and motors, without claiming a detailed design or
burst containment. The six station positions and hoop self-inertia enter the
full locked inertia tensor, including its off-diagonal term.

The power ledger includes motor electronics, generation, eclipse storage,
10 kV distribution conductors at 2% loss, heat rejection, avionics and sensors.
The storage envelope allows the maximum three-axis spin energy swing in eclipse
without regenerative credit. These conservative power envelopes are not an
executed electrical schedule. Optical trim covers 25% of the tile, with a
2 g/m² hardware allowance on that area. Its unqualified bistable target is
5 J/m² per update, 20 updates per orbit and zero holding power: 2.5 GJ per
10 km tile per orbit, approximately 9.85 kW average. The driver and storage
also cover a hypothetical concentration of those updates into eclipse. Nominal
update windows overlap in that envelope; a realizable chronology is not claimed.
The NASA variable-reflectance source is a proposal, not proof of this device.

## Momentum and mass must be solved together

The new solver conserves total angular momentum in inertial space:

`H_inertial(t) = H_inertial(0) + integral(R(t) torque_external_body(t) dt)`

`h_rotor_body(t) = R(t).T H_inertial(t) - I_locked omega_body(t)`

Motor torque changes the body components of rotor momentum; bearings also
transmit the gyroscopic term ω×h. A gravity-gradient torque acts externally.
This diagnostic starts with zero relative rotor momentum and applies no
radiation torque or optical unloading. It is a declared no-unloading comparison,
not an upper bound over every environmental torque or an optimized initial bias.

Rotor and motor allocations are increased to 1.2 times the largest per-axis
requirements, their new masses change the inertia, and the calculation repeats
until allocations agree within 0.1%. Allocations never fall below the historical
starting envelope. This closes the specified rigid sizing loop only.
The final rotor momentum does not return to its initial value: its norm is
about 237 GN·m·s for the reference annulus and 519 GN·m·s for the window.
Periodic operation still needs external unloading. The frame calculation sizes
membrane support; motor-attachment strength and complete out-of-plane slew loads
are reported as unresolved, not supplied by the inertia iteration.

| Case | Annulus, g/m² | Window, g/m² | Meaning |
|---|---:|---:|---|
| Current complete-tile allocation | 13.765 | 26.000 | Inherited flight assumptions |
| Reference, 100 m frame depth | 2,683.5 | 5,870.1 | 5 m film sag and 5 m edge bow |
| 50 m frame depth | 10,625.9 | 23,231.6 | Shallower structure costs much more stiffness mass |
| 20 m permitted film sag | 684.6 | 1,500.5 | Reduced tension; optical shape not revalidated |
| 50 m permitted edge bow | 284.8 | 626.6 | Relaxed straight-edge approximation |
| 400 m rotor radius | 1,973.5 | 4,317.4 | Lower store mass; larger swept hardware |
| 3 km tile | 278.9 | 586.4 | Size sensitivity; formation not redesigned |
| 1 km tile | 86.8 | 158.4 | Size sensitivity; formation not redesigned |

The smaller-tile comparisons scale frame geometry, sag and bow with side,
retain housekeeping and avionics per tile, and initially scale storage/torque
as L³ to retain a radiation-torque-like allocation. Pure inertia at fixed
areal mass would scale as L⁴. The new assembled-inertia calculation is then
applied to each. These are not claims that the old rings can accept smaller
tiles with the same count or coverage.

For the reference annulus, the frame alone is approximately 150 million kg
of the 268 million kg total. The resized per-axis allocation is about
7.44×10¹¹ N·m·s and 8.97×10⁷ N·m motor torque. Window values are about
1.63×10¹² N·m·s and 1.96×10⁸ N·m. This mass changes sail loading so greatly
that the inherited orbit, steering and glow cannot simply be reused.

Copying these rejected reference tiles into the previous 28.233-million-tile
inventory would give about **9,915 Gt**, versus the historical 46 Gt optical
inventory. At the service policy below, whole-tile replacement would exchange
about **760 Gt/year**, or 24 million kg/s. These counterfactual figures expose
the cost; they are not a proposed fleet or film-production requirement.
Repairing components and reusing frames would require a different service model.

## Clearance and local failures

For parallel tiles at 8.5 km along-track pitch and 1.5° shingle angle, their
normal separation is 222.50 m before hardware thickness. A 50 m box envelope
leaves 172.50 m; a 100 m envelope leaves 122.50 m and already fails the inherited
150 m gate. The reference frame therefore fails this conservative local envelope
test at the initial epoch. Compact rotor stations require their actual swept
geometry: a whole-tile box containing 100 or 400 m radius rotors cannot prove
clearance, but its overlap alone does not prove rotor collision.

The fault calculation propagates the failed tile's displacement from its
healthy slot using lunar circular-orbit Hill–Clohessy–Wiltshire equations.
Neighbour slots are prescribed; the residual acceleration is a specified
constant difference in radial/along-track/orbit-normal coordinates. Finite
parallel box geometry and continuous event roots locate clearance loss and
contact. Physical interpretation stops at first contact.

| Envelope depth | Radial residual, m/s² | First loss of 150 m gate | First envelope contact |
|---|---:|---:|---:|
| 1 m | 10⁻⁵ | 63.01 min | 110.92 min |
| 50 m | 10⁻⁶ | 111.80 min | 311.14 min |
| 50 m | 10⁻⁵ | 35.35 min | 97.87 min |
| 50 m | 10⁻⁴ | 11.18 min | 30.95 min |
| 100 m | 10⁻⁵ | Already failed | 82.47 min |

These accelerations are sensitivities, not residuals measured for a failed
version of the newly sized tile. The fault experiment does not include changing
attitude, shadows, eccentric motion, tides, debris or neighbouring controllers.

For the 50 m envelope and 10⁻⁵ m/s² case, an **assumed independent bidirectional
backup** capped at 10⁻⁴ m/s² and activated after 30 minutes keeps sampled
clearance above 154.33 m and returns the tile to its slot. Activation after
60 minutes avoids contact in this run but cannot undo the earlier safety-gate
failure. Recovery does not repair the filter. The backup has not been supplied
as photon hardware, a tug, or a mass/power-accounted independent actuator.

Lost attitude can be more urgent. For the 50 m case, a conservative bound on
surface displacement loses its 150 m guarantee after about **182 s** at a
bounded angular error rate of 0.001°/s, or 18 s at 0.01°/s. Those are losses
of a geometric guarantee, not predicted collision times.

## Replacement at the whole population

The population simulation counts all 28,233,094 slots for one year. It assumes
uniformly staggered 15-year planned renewal, 0.01 unexpected failures per
working tile-year, dispatch capacity 20% above planned renewal, 27-day planned
transit and one-day failed-slot replacement from local spares. Reliability,
life and transit are declared scenarios, not measurements.

This requires approximately **5,153 planned exchanges/day plus 773 failures/day**,
against 6,184 dispatches/day: only 4.35% spare capacity over total expected
jobs. Planned replacement leaves the old tile working until arrival; a failure
opens a slot immediately. Jobs have finite transit, failed slots take priority,
and queues conserve arrivals, dispatches, completions and outstanding work.
Failures do not cancel already planned jobs, an explicit conservative count
approximation. Manufacturing and spare stocks are assumed to support dispatch.

| Scenario | Mean missing slots, fraction | Peak missing slots, fraction | Renewal backlog at year end |
|---|---:|---:|---:|
| Local spares, 1 day | 2.73×10⁻⁵ | 3.06×10⁻⁵ | 0 |
| Delivery from plant, 27 days | 7.11×10⁻⁴ | 7.48×10⁻⁴ | 0 |
| Capacity only 80% of planned renewal | 2.73×10⁻⁵ | 3.06×10⁻⁵ | 658,369 |
| Dispatch stopped for 14 days | 3.56×10⁻⁵ | 4.07×10⁻⁴ | 23,050 |
| One correlated loss of 10,000 slots | 2.91×10⁻⁵ | 3.84×10⁻⁴ | 0 |
| Failure rate 0.1 per tile-year | 9.83×10⁻³ | 1.91×10⁻² | 1,882,206 |

Low missing-slot counts in the capacity-shortage case hide a growing renewal
debt: the model does not invent an age-dependent failure rate for overdue tiles.
The local-spares case has about 140,000 replacements in transit at peak,
mostly planned work. Local buffering has an inventory and logistics cost.

In an ideal periodic 8,900×9,890 m optical rectangle grid at
8,500×9,286.456 m pitch, about 89.1% of ray area has only one covering tile.
For independent outages of fraction q, the additional uncovered fraction is
Σ w_m q^m over coverage multiplicities m. The one-day case's daily sampled
mean is approximately 2.43×10⁻⁵; 27-day delivery gives 6.33×10⁻⁴. These are
**additional ideal geometric gaps**, not total UV transmission or a protection
gate pass. Baseline seams, pinholes, angular illumination and real curved
geometry are absent. The correlated case is not assigned this independent
optical estimate.

Every service result assumes failed hardware has already been contained.
A one-day spare policy cannot address the minutes-scale envelope-contact
windows by itself. No failed object is silently deleted and replaced inside a
trajectory calculation; count accounting and physical containment are separate.

## Design consequence and remaining work

The next hardware candidate should relax the 10 km straight-edge frame:
compare smaller independent units and tension-supported or deliberately curved
membranes with explicit load paths and optical shape. Simply making the hoops
wider does not solve frame stiffness and enlarges their swept geometry.
The window stack/allocation inconsistency must also be resolved.

A candidate needs a coupled flexible-body/controller calculation with the
actual force and shadow history, spectral properties for its complete stack,
finite component clearance through the rolled cycle, and independent fault
containment before it can feed an updated orbital layout or a net electricity
gate. Manufacturing, deployed structural joints, aging, rotor containment,
removal and insertion routes, and repair/reuse remain open. These are engineering
tasks exposed by the calculation, not completed technology claims.

## Reproduction and checks

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.tile_design
python -m pytest engineering/test_array_tile.py engineering/test_array_momentum.py engineering/test_array_service.py protection/dynamics/test_tile_failure.py research/studies/solar_shield_array/test_tile_design.py
make check
```

The study uses about 7 CPU seconds and 90 MiB in the recorded run with one
numerical thread. Its 50 focused tests pass: membrane series and stiffness,
mass/energy accounting, angular-momentum conservation and frame invariance,
analytical local orbits and between-sample contact, delayed recovery, queue
conservation, optical multiplicity, provenance and study integration.
Doubling the attitude resolution from 1,441 to 2,881 samples changes peak
momentum by 0.0101%. Reducing local maximum steps from 30 to 10 s changes the
reported contact and clearance roots by less than 0.1 s. These are numerical
checks; empirical validation and full scholarly source review remain separate.

`make check` passed: **2,053 Python tests passed, 16 skipped**, and all
**191 JavaScript/calendar tests passed** (129 immersion, 48 column/transport,
14 calendar). Layer checks found zero violations; historical provenance,
protection verification and ensemble integrity passed. The full check was
rerun on the final model, source, product and test snapshot immediately before
commit. [tile_checks.json](tile_checks.json) records the commands, counts and
exact final model/test hashes.
No full-fleet simulation, physical qualification or manuscript review is claimed.
