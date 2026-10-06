# Solar shield and habitat array

The [relative-orbit plan](relative_orbits.md) of 6 October 2026 sets the
current direction. Gravity-only replays of the committed seed show that the
depth stack spreads the 361 tiles' semi-major axes over 59–65 km. That spread
drifts them up to 345 km apart along-track each revolution, which is the
expanded search's return error. One energy-matching burn of 0.42 m/s rms per
tile cuts that drift to 70 m. The hour-13.8 clearance failures come from the
pattern folding through its orbital plane a quarter orbit after mid-service.
The plan redesigns the formation in relative orbital elements, with equal-energy
coasts and e/i-vector separation, and verifies it in the coupled model against
the full protection gates. A zoned-aperture and collection study runs alongside.
Its decision gate chooses between the orbiting fleet and a held screen with a
zoned aperture.

The [expanded coupled-cycle investigation](expanded_cycle.md) adds
in-plane and arrival-acquisition controls and measures nine full-cycle response
columns. The local terminal matrix is full rank in those columns, with 34.2%
of the normalized residual outside its span. Exact replay exposes 27.6 km
terminal prediction error and fails clearance at hour 13.795 with ample
installed actuator headroom. A tenfold smaller step still mispredicts by
2.86 km. The 100 TJ next-service cutoff binds the local solve; relaxing it
leaves clearance unresolved. No cycle or physical resource bottleneck is
established. The [portable continuation package](local_continuation/README.md)
includes compact ephemeris inputs, saved optimizer state and verified replay
commands. The 36.186 TJ slower-turn baseline remains selected.

The preceding [coupled benchmark](coupled_optimization.md) and its eight-mode
failed search remain preserved. The current result calls for smaller signed
coupled stencils, localized separation freedom and exact trust-region replay.

The author directed this work on 4 October 2026: develop the combined solar
shield, power system and habitat fleet, and investigate the holding problem
with a full Sun–Earth–Moon ephemeris. The branch starts at main commit
`2f2e0a104a3a1ea962a6dead4179f522762f997e`.

[design.md](design.md) records the proposed architecture before dynamics work.
The [full-ephemeris study](report.md) compares the September held screen with
conserved solar momentum, variable sunward distance, transverse motion,
eclipses and propagated control. The selected monthly trajectory reduces the
50-g/m² array's mean holding power from about 309 to 238 TW under the reference
propulsion assumptions. It still consumes about 371,000 kg/s of propellant.

[Photogravitational cycling](cycling.md) now has propagated point-tile
candidates at the same 50-g/m² base mass. A 15,000-km retrograde seed with
filtering during service and feathering between passes stays bounded for
three years and averages about 15.7% useful projected area. A separate
optimized 8.9-day solar-sail arc closes within 17 m on independent replay.
These establish orbital motion and an inventory comparison.
[tacking.md](tacking.md) retains the initial literature review.

The [simultaneous fleet study](fleet.md) now tests five common-epoch populations
with multiple planes and radii, finite square shadows, solar-disk quadrature,
seams, swept collision and mutual-shadow exclusions, and Earth-safe ideal
reflections. Continuous coverage remains unachieved. A 3,072-member large-square
geometry probe averages 61.76% interception before exclusions; its selected
19-member subset averages 0.78% and reaches zero. Those 1,000 km squares need
extended-body dynamics; the 10 km replay also fails the 50 m phase-accuracy
gate. No sufficient inventory or collector/habitat power is demonstrated.
Four lunar radii remains primary, with three radii as a comparison.

The [active follow-up](active.md) returns to separate 10 km tiles. Bounded
electric control closes explicit gaps in a 625-tile local formation, and
independently replayed orbital transfers shift a tile's projected arrival by
about 1,400 km with 5.45–6.79 m/s of correction. Power and propellant are
reported. A collision-safe global assignment, recurring handovers and return
schedules remain to be solved; delivered power is still unevaluated.

The [global traffic follow-up](active_global.md) rejects a prescribed-circle
family whose large circulating inventory and continuous holding cost exceed
the held-screen benchmark. Its completed finite-square geometry tests and
interrupted checks are distinguished explicitly. A gravitational release
diagnostic measures three days of coverage loss and independently replayed
collisions when circle enforcement is removed; its seven-day interpolation
fails the position gate. Photon-force cancellation is only one cost component. A
minimum-intervention controller, its complete recurring budget, continuous
optical certification and bounded-slew attitudes remain unresolved.

The [coverage-led pilot](pilot.md) adds a bounded, common-epoch search over
54 near-natural orbital families with the same sail force and 10 km tiles.
Fractional population and initial-phase-band allocation, followed by a
coverage-bottleneck refinement, gives a 17.29-million-tile **area-capacity
estimate**. This supplies conditional fleet impulse budgets. Individual
placement, consistent mutual illumination, safe attitude histories and
recurring corrections remain unclosed; no operating fleet is accepted.

The [natural-pattern continuation](patterns.md) resolves one prescribed-roll
rate witness and propagates a 361-tile, naturally deforming local pattern
with mutual shadow force. It retains sampled finite-Sun coverage of a moving
20 km window for six hours with zero electric translational thrust. Its
kilometre-scale endpoint mismatch leaves acquisition, handovers, return and
recurring cost unresolved. This local result does not establish global
coverage or a sufficient fleet inventory.

The [bounded return search](returns.md) starts from those actual terminal
states. Its low-impulse endpoint proposals develop finite-square conflicts
during coast. Coupled propagation and a tighter replay reach the 100 m
clearance limit 47.69 minutes after service, as the common attitude changes
toward feathering. The full return and two-pattern handover were not executed;
transition timing, intermediate spacing and attitude choices remain coupled
design variables.

The [prepared-departure search](departure.md) resolves that early obstruction
in the ideal local model. A single smooth layer/row/column preparation burn
lets all 361 tiles turn and coast for three hours, with a 107.7 m conditional
swept clearance bound and 0.59 m independent replay disagreement. Translation
costs 1.41 m/s; an explicit electric torque-couple option, with disturbance
allowances, raises departure to 5.04 m/s equivalent impulse. Its provisional
power hardware changes the mass and still needs a dynamics replay. Handover,
full return, recurring cost and continuous global coverage remain open.

The [collection inventory](collection.md) tracks incident optical power and
integrated energy per tile on the saved service and departure trajectories.
It includes orientation, solar distance, both faces and finite-Sun mutual
shadows; the current prefix model rejects body eclipses. Subsequent inventory
and control searches must retain collection during every command stage,
including coast, with spectral allocation, electrical conversion, propulsion
and delivered power recorded separately.

The [energy-first continuation](energy.md) prioritizes energy at the existing
100 m clearance constraint. It compares longer Sun-facing passages and smaller
or slower turns with optimized preparation impulses. Coupled failures and a
measured-defect correction remain visible beside any accepted result.
Spatial shadow moments and finite-area gravity now supply an explicit rigid
attitude-load calculation for comparison with the earlier torque envelopes.

The [starting-arrangement pilot](packing.md) varies four-level depth spacing
and order, and offers two independent smooth burn vectors to every member.
It preserves the six-hour local service obligation before testing departure
cost, with coupled shadow forces, tighter replay, actual surface-distance
checks and a search for poorly covered solar-source/date combinations.
A 36 km total starting depth permits the full twelve-hour service/departure
with zero electric translation. Calculated control energy falls 82.13% from
the previous verified passage, with 1.76 km closest sampled surface separation.
Full return, handover and a recurring fleet cost remain unexecuted.

The [cycle-closure continuation](closure.md) uses those actual terminal states
to test four independent burn arcs, arrival timing and service reassignment.
It sizes proposed return power, propagates the resulting added hardware mass
from the initial epoch, and independently replays the decisive return event.
The twelve-hour passage survives a 25.4% additional-mass allocation with zero
electric translation and 1.72 km closest sampled surface separation. The
restricted return proposals require about 1.85 m/s of translation and fail
their clearance tests. Per-member propulsion and optical accounts distinguish
executed stages from proposed recurrence; no full cycle or fleet is accepted.

The [encounter-directed continuation](return_repair.md) tests early burn changes
and smooth in-plane rotations from the actual hardware-loaded states. A local
correction delays the 100 m stop by about 90 minutes, independently replayed
within 10.58 m, but reaches a new deep overlap and exceeds power allocations
on 160 tiles. Both roll transitions also fail. The original twelve-hour
passage is preserved. Comparable optical source accounts, explicit conversion
sensitivities and propellant/replacement trades retain the 23.835 TW target
beside 50 and 100 TW alternatives. No complete return or revised target is
accepted.

The [electromagnetic feasibility study](electromagnetic.md) uses the actual
mass-loaded group, time-resolved loads and extended square coils to compare
local generation/storage, microwave/laser and inductive transfer, magnetic
control and solar-wind scales. It favors local Sun-tracking generation and
short buffers, conditional microwave assistance and a narrow magnetic-trim
follow-up. A static return allocation improves about 10% after cooling, while
tested departure states do not improve. No new trajectory, plasma protection,
three-function architecture or relaxed operating target is accepted.
Its source register, compact results and `electromagnetic_checks.json` retain
the assumptions and executed checks separately from proposed follow-up work.

The [primary magnetic architecture comparison](magnetic_architecture.md)
reopens the four lunar stations. Finite regional loops of 500–1,000 km radius
greatly reduce the screened conductor/support mass at matched field requirements.
Stored atmospheric cases motivate testing a weaker moment, which also reduces
tile torque; storm compression and ion escape remain unresolved. Held upstream
sources retain power/mass failures near the array and substantial propellant
and wake-width penalties farther away. The bounded calculations are executed;
larger plasma and fleet campaigns remain paused for review.

The [joint departure and next-service investigation](joint_cycle.md) returns
to ordinary photogravitational/electric control, with the primary magnetic
shield set aside. Twelve bounded proposals test earlier velocity changes,
packing/assignment, staggered arrivals and grouped sail normals. A full
361-tile separated corridor returns at hour 45.77978 and completes another
six-hour optical passage on independent replay, with a 900.28 m conditional
clearance bound after numerical allowance. It is an ideal-force diagnostic:
145 tiles overload installed actuators, the executed sequence costs 1.53538 PJ,
and the early preparation loses the low-energy departure advantage. A cheaper
actuator-compliant execution stops at 100 m at hour 12.75076. Optimistic local
collectors with short buffers leave 27.067 TJ unserved in lunar eclipse.
No operating return, repeatable cycle, supply hardware, handover or fleet
budget is accepted. Source-bound results, exact initial states and the next
small local test are preserved; a larger fleet campaign remains paused.

The [service-led natural-return search](natural_return.md) broadens arrival
dates and frees endpoint states under energy and individual-power caps.
A slower departure turn preserves the loaded twelve-hour, zero-translation
prefix at **36.186 TJ**, **30.22% below** 51.859 TJ, with 1.722 km minimum
sampled clearance. Coupled replay still reaches 100 m at hour **13.799825**.
The cheap ray-fit candidate loses initial coverage and differs from coupled
propagation by 5.32 km; grouped sail pulses also fail their reduced screens.
No new return, next passage, delivered supply, handover or cycle is accepted.
The sixteen products preserve the failed branches, one duplicated assignment
variant, measured prefix accounts and the inventory penalty of longer returns.
The next small test targets early formation compression with coupled shadows.

The study couples protection, engineering, illumination and habitation. Its
runners, compact numerical results and interpretation belong together here;
dynamics components belong in the protection domain. Computation files and
external ephemeris kernels go in ignored `research/runs/`. The historical
protection implementation and its imported results remain byte-pinned.

## Reproduce the study

```sh
python -m pip install -r research/studies/solar_shield_array/requirements.txt
python -m protection.dynamics.ephemeris --download
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.magnetic_architecture_run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.validate
python -m research.studies.solar_shield_array.publish
python -m research.studies.solar_shield_array.tacking_screen
python -m research.studies.solar_shield_array.cycling_replay
python -m pytest protection/dynamics research/studies/solar_shield_array
python visualization/solar-shield-array/plot.py
python visualization/solar-shield-array/cycling_plot.py
```

The runner writes raw results to `research/runs/solar_shield_array/`.
`publish.py` checks source identity, relevant constants, time/area convergence
and the matching validation product before exporting compact
[holding](results/holding.json) and [validation](results/validation.json)
snapshots. The plot and its provenance are generated locally under
`visualization/solar-shield-array/results/`. `--quick` provides the centre-force
pilot without optimization, nodal-span propagation or the final exporter.

The final comparisons preserve a four-lunar-radius protected region and the
selected climate spectrum. Optical control is an optimistic lower bound;
10-g/m² films require a replacement optical stack. The 19-year test and ideal
feedback propagation establish numerical behavior for the tested families.
Final holding architecture, actual filter/PV forces, safe plumes, tile seams,
habitat capacity and long-term resource closure remain open.

## Checks

`make check` was run before this design-only commit: 545 Python tests passed,
55 skipped and 13 failed. The failures concern absent solar-spectrum inputs,
CM1 build configuration, the ExoPlaSim process lookup and ring-comfort tests.
No executable model was changed in that initial commit. Current targeted
and repository-wide check status, including the fleet phase-accuracy failure,
and the numerical
verification are recorded in [checks.json](checks.json) and
[validation.json](results/validation.json).

The cycling [report](cycling.md) gives the full run matrix, convergence limits,
source/constant provenance and reproduction commands. Its annual useful-area
statistic is numerically well resolved; orbit phase and perturbed tile
positions drift by kilometres, so metre-scale seam control is not established.
