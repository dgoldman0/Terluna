# Local handover and orbital return

## Bounded experiment

This continuation starts from remote
`4e8a7fd99039ffddcfeef6a52685af578edb1bb4` on
`research/solar-shield-habitat-array`. The remote head was checked before
starting. Main remains unchanged. The input is the source- and byte-verified
361-member terminal state from the tighter six-hour `patterns.json` replay.
The 50 g/m² optical allocation, physical 10 km squares, four-lunar-radius
target and existing climate spectrum remain fixed.

The objective is to determine whether that actual local pattern can hand off
its protected rays and reach another useful service passage at modest impulse.
Follow an evolving ephemeris reference through its next sunward passage;
do not prescribe a Sun-following circular orbit between services. Return date,
tile-to-lattice assignment, smooth maneuver coefficients and attitude roll are
decision variables. The cost is the sum of individual maneuver impulse,
converted separately to propulsion energy and propellant. Include both
departure and arrival preparation; retain unresolved recovery and repetition.

Start with a linearized shooting screen using actual terminal offsets and
velocities. Try at most 24 return-date/assignment/roll combinations. Rotating
the target assignment permits a tile to occupy another lattice location on
the next pass. Inspect finite-square encounters along the proposed arcs,
including near the relative-motion focusing times. Use offending pairs to
refine the trajectory/control proposal, up to five bounded correction steps.
Keep failed assignments and their witnesses. An infeasible restricted family
is not a general impossibility result.

Promising candidates proceed to individual nonlinear propagation with sail
force, mutual illumination and finite-Sun body eclipses. Use smooth burns,
an assumed translational acceleration cap of 0.001 m/s², a 0.1°/s attitude
rate cap, the existing reflected-beam exclusions and 100 m finite-square
clearance. A local handover uses a second independently propagated pattern
and tests the same receiver rays at simultaneous actual dates. Coverage,
encounters, endpoint defects and force-reduction assumptions must be checked
again on the propagated states. A tighter replay must agree within 50 m.

Allow **60 minutes of numerical wall time**, one calculation process with
one BLAS/OpenMP thread, 2 GiB address space and 1 GiB raw output. The budget
includes at most three nonlinear 361-tile return candidates, one 722-tile
short handover and one tighter replay of the strongest useful case. Check
resource/time limits between stages and save intermediate evidence. Repository
checks have a separate five-minute allowance. Stop at the first meaningful
feasibility or cost checkpoint, including an encounter or optical obstruction
that needs a wider control family.

Success means that a tested local handover/return preserves its assigned
finite-Sun rays, separation, safe beams and bounded attitude, and that the
measured all-member impulse leaves room beneath the conditional cost thresholds
in `patterns.md`. Sampled success and a continuous certificate are different
evidence states. A single return cannot establish long-term recurrence or a
complete fleet. The 238.35 TW held-screen benchmark retains its stated ideal
optics, mass closure and propulsion assumptions. Collection capacity and the
100 TW delivered-power ambition receive their own ledgers.

The plan was pushed before the search as
`ed1b05d231092972f62be13741a8b44392f8bac4`. The initial `make check` reached
the Python target but the refreshed execution environment lacked pytest.
The recorded pytest 9.1.1, Shapely 2.1.2 and jplephem 2.24 dependencies were
restored before numerical execution.

## Checkpoint result

**The tested return fails finite-tile clearance while changing from service
attitude to coast attitude.** The coupled propagation reaches the 100 m
clearance limit **47.69 minutes after the six-hour service passage**, with
tiles **338 and 339** as the limiting pair. A tighter propagation agrees
within **0.00990 m** and locates the same event within **0.000055 s**.
Integration stops at the clearance limit. This does not claim an actual
material impact in the coupled model.

The cheaper endpoint calculation finds a mean correction of about 0.4 m/s,
but it has extensive finite-square conflicts along the way. The selected
local encounter constraints are infeasible within its capped three-arc
parameterization. The result identifies a packing/attitude obstruction that
must be resolved before assigning a recurring operating cost. It leaves the
space of other schedules, lane separations and individual attitudes open.

The previously demonstrated six-hour service passage remains intact. No
two-pattern handover, complete orbital return or second service passage was
executed at this checkpoint. Those stages were withheld after the return
proposal failed its safety constraints, within the original stopping rule.

## Actual-date endpoint search

The input is every individual terminal position and velocity from the
16-source tighter replay in `patterns.json`, checked against the saved source,
constant and raw-array hashes. The natural reference starts at the mean of
those states at hour 6. Its point-mass gravity and ideal service/coast sail
force continue through the next sunward passage, near **hour 46.270** from
2026-10-04 TDB. No intermediate circular orbit is prescribed.

Three return dates lie 15 minutes before, at and 15 minutes after that
passage. Each uses four tile assignments, rotating tile identities through
0°, 90°, 180° or 270° of the same target lattice, and two smooth coast-roll
excursions, 0° or 45°. These are **24 geometry cases**, using twelve distinct
endpoint LPs because roll does not affect the gravity-only variational map.
The target velocities prepare the earlier local service pattern in the new
Sun frame. The next six-hour service dynamics have not been propagated.

Each tile has three non-overlapping sin² acceleration arcs. Their durations
are two hours, beginning 15 minutes after the saved service endpoint, centred
on the coast interval, and ending at the chosen return date. Each vector
coefficient has units of integrated velocity change. One shared 6×9 response
matrix and block-sparse constraints serve all 361 members. Endpoint position
and velocity are equality constraints. The LP minimizes the sum of absolute
inertial components; the reported impulse uses the Euclidean vector norms.
This objective is a surrogate for fuel minimization.

A conservative component box ensures that each vector acceleration stays
below 0.001 m/s². It is an inner approximation of the spherical acceleration
limit. The gravity variational map omits unequal shadow force and sail-force
gradients; it proposes controls for later nonlinear verification.

| Return date, hours | Identity-assignment mean proposal Δv | Other three assignments |
|---|---:|---|
| 46.020 | 0.43766 m/s | Endpoint LP infeasible under the component cap |
| 46.270 | 0.40028 m/s | Endpoint LP infeasible under the component cap |
| 46.520 | 0.39495 m/s | Endpoint LP infeasible under the component cap |

The two roll choices have the same endpoint impulse. Every solvable geometry
case has sampled finite-square clearance violations. The inherited normal
history also contains a later coast-rate excursion: coarse sampled maxima
range from 0.456 to 0.928°/s against the 0.1°/s cap. These are failures of this
chosen schedule. The attitude history is a smooth transported rotation spline
with the exact Sun-facing service planes at its boundaries; smoothness alone
does not enforce a rate limit.

## Encounter-directed correction

The selected diagnostic is case 1: the earlier date, identity assignment and
45° coast-roll excursion, which has the fewest coarse sampled clearance
violations. All available candidates already fail the attitude-rate screen;
this selection is for diagnosing their geometry, not for accepting a return.

Its 120 s refinement finds **166,185 pair/date clearance violations**. These
counts include repeated pairs and do not represent that many distinct impacts.
An adversary then roots the signed distance between each candidate pair's
parallel planes between snapshots and checks both finite-square edge
projections at each root. It finds **1,848 physical square-intersection events**
in the linear proposal, plus 21 other crossings within the clearance guard.
The retained worst event has less than 709 m separation in either in-plane
coordinate for two 10 km squares. Each intersection witness requires
coplanarity and overlapping finite edges.

The first correction adds 512 local constraints: 256 crossing-escape cuts and
256 static clearance cuts. At a coplanar crossing the cut chooses an in-plane
edge separation, so a small shift of crossing time alone cannot satisfy it.
The capped LP reports infeasibility. This rejects that chosen set of separating
half-spaces and three control arcs. Other axis choices, burn schedules,
attitudes and intermediate formations are not excluded. Further correction
iterations were not run after that failure.

## Coupled verification of the early obstruction

The diagnostic replay propagates all 361 actual states with the selected
control coefficients. A new parallel-plane illumination calculation permits
the common attitude to tilt throughout the transition. For each finite-Sun
source it constructs dynamic candidate blockers using conservative perspective
bounding disks, then evaluates exact rectangle unions on each receiving tile.
It uses no fixed lattice-neighbour assumption. Both illuminated faces are
handled with the signed normal momentum transfer. A 3×3 quadrature supplies
finite-area gravity. The ten-kilometre physical side is used for separation;
the 9.89 km clear side is used for optical force.

This implementation is tested against all-pair shadow unions at several
tilts, including both sides of grazing incidence, and against the previous
Sun-facing coupled model. It explicitly rejects partial Earth/Moon eclipses.
The entire executed prefix is uneclipsed. Extending it into a night-side
return would require the joint body/array illumination calculation promised
in the plan; that unexecuted extension is not supplied by this prefix model.

| Coupled diagnostic | Main prefix | Tighter prefix |
|---|---:|---:|
| Solar-force directions | 8 | 16 |
| Maximum integration step | 60 s | 30 s |
| Relative integration tolerance | 2×10⁻¹⁰ | 2×10⁻¹² |
| Time after service at the clearance event | 2,861.411041 s | 2,861.411096 s |
| Limiting pair | 338, 339 | 338, 339 |
| Clearance at stop | 100 m | 100 m |
| Minimum sampled reflected-beam margin | 99.416° | 99.416° |
| Maximum sampled attitude rate | 0.05404°/s | 0.05404°/s |

The two prefixes disagree by at most **0.00990 m** and
**7.18×10⁻⁶ m/s** at common comparison times. The 50 m replay gate passes.
The maximum executed translational command is **8.44×10⁻⁵ m/s²**. Thus the
verified early failure is the clearance constraint; the later attitude-rate
violation is a separate failure in the longer proposal.

At the stopping state, the limiting pair's in-plane projections are about
8,401 m and 28 m, with 100 m normal separation. Their centre distance is
8,402 m. The same saved positions with all squares restored to a Sun-facing
plane have **2,616 m** minimum clearance. This frozen-state diagnostic isolates
the geometric effect of the attitude transition; an alternative force history
was not propagated. If such a pair becomes coplanar while its centres remain
closer than 10 km, changing square roll alone cannot remove the overlap,
because each square contains a disk of radius 5 km.

The tighter prefix is also checked for coplanarity crossings between its
stored samples using a 30 s Hermite reconstruction and finite-square tests.
This audit finds no physical intersections before the clearance stop. It is
an additional search for missed events, not a continuous safety certificate.
The linear proposal differs from that coupled prefix by about **58.1 m**,
using a 120 s Hermite reconstruction of the proposal. That discrepancy already
exceeds the placement allowance, reinforcing the need for nonlinear force
verification before accepting an optimized endpoint trajectory.

## Impulse, energy and unresolved recurrence

The executed prefix starts the first smooth burn and stops before it ends.
Mean applied Δv is **0.02226 m/s**, with **0.06118 m/s** maximum. For the
1.805×10⁹ kg optical inventory and the existing 30 km/s exhaust, 70% efficiency
and 45° cant, that corresponds to **1.218×10¹² J** and **1,894 kg** of propellant.
Power-system mass, propellant feedback, attitude energy and complete hardware
are not included in that prefix conversion.

The selected endpoint-only proposal's 0.43766 m/s, divided by the elapsed
epoch-to-return time, would be 0.22825 m/s/day. Its geometry fails, so this
number cannot enter the fleet operating budget. The measured recurring impulse,
fleet propulsion demand, collection capacity and delivered power remain empty
in `return_budget.json`. No handover, complete return, recovery or recurring
maneuver frequency was established. Initial six-hour service still has the
previously measured zero electric translational impulse.

The inventory/impulse trade in `patterns.md` remains the comparison framework.
At 3 m/s/day, its 17.29, 25 and 35 million tile scenarios require 114.25,
165.24 and 231.34 TW respectively, with the stated added mass and thrust duty.
Those inventories are conditional. The 238.35 TW held-screen benchmark is the
annual ideal-optical, mass-closed reference at 50 g/m², with 300 W/kg power
hardware, 25% peak allowance and seven days of propellant, excluding storage
and full vehicle engineering. No collector capacity or credit toward the
100 TW delivered-power ambition follows from either calculation.

## Resource use and next decision

The screen took **48.1 s** and the two coupled prefixes **178.9 s**, about
**3.8 minutes** together. Peak recorded RSS was **252.7 MiB**. Execution used
one calculation process, one BLAS/OpenMP thread and one HiGHS thread. The new
compressed raw arrays occupy **21.75 MiB**. The
search stopped on the verified obstruction, leaving most of the 60-minute
allowance unused. One of the five allowed encounter-correction attempts
was made; its LP was infeasible. One nonlinear candidate prefix and one tighter
prefix were executed. The remaining return candidates, longer propagation
and 722-tile handover were not run.

Run the LP script in a fresh process as shown below. HiGHS retains its worker
pool globally and rejects a changed thread count after another solver has
initialized it. The analytic LP test uses the same isolated, one-thread launch
as this calculation. Mixing it after the older capacity LP in one Python
process is outside that launch contract.

The next useful optimization should include the timing of the attitude
transition and the pattern's intermediate depth/lane spacing, beginning early
enough to avoid this clearance loss. A second service pattern would permit the
outgoing pattern to change its optical shape during handover. Those changes
must be solved together, with the old and new patterns' mutual force and
finite-tile geometry checked at the same dates. More flexible maneuver timing
and per-tile attitude choices remain candidates; a particular repair has not
been established. The present result supports a concrete obstruction and a
reproducible test for it, while the lower-cost feasible region remains open.

The compact products are [return_screen.json](results/return_screen.json),
[return_prefix.json](results/return_prefix.json) and
[return_budget.json](results/return_budget.json). They preserve failed cases,
event witnesses, constraints, source/constant identities and raw-array hashes.
Prior pattern and capacity products remain unchanged. With their raw inputs
present, reproduce using:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.return_screen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.return_prefix
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.return_budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```

All **94 study/domain tests pass**, including an analytic double-integrator
shooting limit, smooth-arc normalization, between-sample impact detection,
dynamic versus all-pair shadow unions, the Sun-facing force limit, source
identity and the stopped-prefix ledger. `make check` reports **639 Python
passes, 55 skips and the same 13 recorded failures** involving missing spectral
inputs, CM1 configuration, process lookup and ring-comfort checks. The layer
check passes over 469 files; the later JavaScript, historical-provenance and
ensemble targets were not reached. `git diff --check` passes.

This is the requested feasibility checkpoint. `checks.json` records these
scopes. Continuous coverage, engineering realizability and full-fleet
acceptance remain separate gates.
