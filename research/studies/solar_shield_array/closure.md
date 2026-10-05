# Closing the low-energy passage

## Bounded continuation

The remote research branch was checked at
`85b9851fae79fb3304047481902f9086a8091b60` before this work. Main remains at
`2f2e0a104a3a1ea962a6dead4179f522762f997e`. The author asks to resolve the
return, handover, hardware-mass and full-coverage obligations left by
`packing.md`. Preserve that independently replayed local result and its
82.13% energy reduction.

Start from its actual twelve-hour terminal states. The objective is a
repeatable service cycle with a measured mass and impulse account, preserving
coverage and the 100 m physical clearance constraint. Use the existing 50 g/m²
base optics, 10 km physical squares, four-radius target and 238.35 TW held-screen
comparison. The 23.835 TW goal and intermediate inventory/cost scenarios remain
separate from optical collection and the 100 TW delivered-power ambition.

Begin with inexpensive return and attitude bounds. Follow a freely evolving
reference to the next actual-date sunward service passage. Vary arrival date,
return-to-service turn duration, member assignment and independent smooth
control vectors. Use up to twelve endpoint/attitude proposals and at most
four encounter-directed corrections per selected proposal. Permit four
non-overlapping acceleration arcs per tile, bounded at 0.001 m/s², with the
existing 0.1 degree/s attitude cap. Exact linear endpoints and selected
separating half-spaces propose a solution; coupled propagation must verify it.
Retain useful infeasible branches and distinguish them from global exclusion.

The next service packing is a decision constraint to test, not an unavoidable
requirement for every economical fleet. A return that cannot reacquire this
packing may still protect other rays under another assignment. Compare the
cost of this reusable pattern with the earlier optimistic inventory bound;
do not multiply a small fully protected receiver window into an accepted
global population. Overlapping finite-Sun rays and pattern interfaces need
joint geometric evaluation.

Account explicitly for the return-to-service attitude expenditure. Include
actual moments on executed trajectories, individual burn counts, preparation
and recovery obligations. Test the force response to additional hardware
mass, beginning with the previously stated 20% fixed-mass scenario and the
calculated propulsion-power allocation. State the assumed mass distribution;
a sensitivity run does not establish a complete vehicle design. Any eclipse
outside the available joint illumination model stops acceptance rather than
receiving an unmodelled force or collection credit.

If a return candidate passes the energy and geometry gates, propagate its
next service and a short two-pattern handover at simultaneous ephemeris dates.
Include cross-pattern shadows and finite-tile encounters. If a reusable
pattern cannot pass its own local cycle, preserve the limiting constraint
before instantiating a much larger fleet. Inventory, service assignments,
initial phases and other orbital families remain available to later searches.

The numerical ceiling is **4,800 s wall time**, one calculation process and
one numerical thread, **2 GiB address space** and **768 MiB new raw output**.
Allocate at most 600 s to cheap proposals/bounds, two coupled continuations,
one hardware-feedback passage, one short 722-tile handover and one tighter
replay of the decisive accepted or rejected trajectory. These are maxima,
not a promise to execute every stage. Reserve 600 s for load/collection and
geometric audits; check the remaining budget between stages. Repository
checks have a separate five-minute allowance. Publish this plan before
numerical search.

A positive cycle result must preserve assigned finite-Sun rays, safe reflected
beams, finite-tile separation, bounded attitude and thrust; it must reacquire
a service trajectory and agree with tighter replay within 50 m. Report added
mass and measured recurring obligations separately from unexecuted hardware,
recovery, recurrence and global coverage. Seek between-sample encounters and
missed rays. Stop at the first meaningful feasibility or cost checkpoint,
including an independently verified obstruction; no sampled result becomes a
continuous certificate.

## Execution decisions

The plan was pushed as `ffb8721696bcf93386871546b046e19a76edc551` after
`make check` reported 675 passes, 55 skips and the same 13 existing failures
in 29.97 s. The twelve endpoint proposals and three selected avoidance
corrections complete within the 600 s screen allowance. None passes the
restricted return geometry. The selected diagnostic is case 6, the next
reference passage with an eight-hour arrival turn and unchanged identities.

The first hardware sensitivity allocates power for the demonstrated twelve-hour
passage. The proposed return's stronger burns require another power allocation.
Before executing that return, add one more coupled twelve-hour hardware passage
with power sized from its proposed burns and attitude/disturbance envelopes.
This revises the hardware-passage count from one to two, within the unchanged
4,800 s wall, one-process, memory and output ceilings. Refit return controls to
those new actual terminal states. The one tighter replay will cover the entire
mass-loaded sequence from the initial epoch through the decisive return event.
Actual load accounting must check the installed power at every sampled member
and retain any insufficient margin as a separate constraint.

Before load accounting, set its disturbance and optical sampling to 120 s
with 16 solar sources and 240 s with 32 rotated sources. The dense return
geometry makes the original shadow-moment calculation expensive; this keeps
the account within its 600 s ceiling. Compare 120 s with 240 s at the same
16 sources, then compare source quadratures at matching 240 s dates.
Nominal attitude and smooth burn power retain two-second sampling, and this
choice does not change the trajectory or clearance checks. Treat any material
accounting discrepancy as uncertainty rather than an accepted saving.

## Return screen

The freely evolving feathered reference reaches the next sunward service
phase at **hour 46.27978** from 4 October 2026 TDB. Three arrival dates bracket
that passage by half an hour. Each uses a four- or eight-hour quintic turn
back to the service attitude, and either the existing tile identities or a
180-degree permutation of the target grid. Every member has four independent
two-hour acceleration arcs: immediately after hour 12, near mid-return,
six-to-four hours before arrival, and the last two hours before arrival.
There are **4,332 signed control coefficients**.

The proposal baseline propagates every actual member independently with
finite-area gravity and finite-Sun sail force, omitting mutual shadows.
The gravity response map adds control perturbations. The target is the same
36 km deep service packing, with a velocity gradient prepared for another
six-hour passage. Returning to that specific relative state is a restricted
service assignment. Its cost can change if another pattern or assignment is
allowed at the next passage.

| Arrival hour | Turn duration | Endpoint-only translation, m/s | Prior passage + translation + nominal arrival attitude, m/s |
|---|---:|---:|---:|
| 45.77978 | 4 h | 1.86973 | 3.31002 |
| 45.77978 | 8 h | 1.86969 | 2.96909 |
| 46.27978 | 4 h | 1.84741 | 3.28770 |
| 46.27978 | 8 h | **1.84758** | **2.94699** |
| 46.77978 | 4 h | 1.88465 | 3.32494 |
| 46.77978 | 8 h | 1.88500 | 2.98441 |

All six 180-degree assignments are infeasible under the conservative
component acceleration bounds. The two latest identity cases also let their
reference receiver leave the four-radius target during the proposed next
six-hour service. The remaining four identity cases meet the endpoint
equations but contain finite-square conflicts.

The eight-hour arrival costs about **0.34154 m/s** of nominal attitude
impulse. The selected combined-cost proposal's endpoint LP has component
objective 2.42986 m/s;
dividing by sqrt(3) gives a **1.40288 m/s** lower bound for the Euclidean
translation objective within those fixed linear endpoint equations. This
bound already exceeds the previous 0.53987 m/s remaining allowance. It does
not bound other service assignments, initial velocity fields or nonlinear
trajectories. The 2.94699 m/s total also leaves return disturbance torque,
encounter repair, handover, recovery and recurrence unpaid.

The three selected refinements use 120-second dates and search for actual
coplanarity between them. Cases 6, 2 and 4 have respectively **3,554, 3,625
and 3,554** predicted finite-square intersections, including repeated
pair/time events. Each adds 128 crossing-escape and 128 static separation
cuts. All three LPs report infeasibility for those selected half-spaces.
The predicted closest crossings place two centres only tens of metres apart
around hour 14.2. Beam margins remain positive and attitude rates remain
below the cap. These are specific geometry failures in the chosen return
family; an endpoint solution alone receives no operational acceptance.

## Hardware feedback before return

The first sensitivity adds 20% fixed hardware mass and the power equipment
sized for the previous twelve-hour passage. Its mean total-to-optical mass
ratio is about 1.212. Starting from the same initial states, the changed sail
acceleration moves individual terminal positions by up to **4.88 km** relative
to the optical-only result. The six-hour local coverage still passes, and
the departure's closest sampled surfaces remain **1.656 km** apart. This is
a distinct preliminary mass experiment, not another group of service tiles.

That power allocation does not cover the stronger proposed return burns.
The second allocation uses each member's larger prefix or proposed-return
peak, including nominal attitude and a conservative disturbance-couple
envelope, with 10% refit headroom and the scenario's 25% installation margin.
The power mass feeds back into the required power before propagation. The
resulting mass distribution is then carried from the original epoch, through
service and departure, into the return trial.

| Allocation for 361 tiles | Value |
|---|---:|
| Optical mass, retained at 50 g/m² | 1.805 billion kg |
| Assumed additional fixed mass | 361 million kg |
| Power hardware at 300 W/kg | 97.298 million kg |
| Total carried mass | **2.26330 billion kg** |
| Mean mass relative to optical allocation | **1.25390** |
| Installed propulsion capacity | **29.190 GW** |

The fixed and power masses are uniformly distributed over each rigid square;
gravity acceleration and specific rigid-body inertia therefore retain the
uniform-square model, while photon acceleration decreases with total mass.
Actual gravity and shadow torques are recomputed on the resulting trajectory.
The mass is constant over this diagnostic passage. Physical hardware layout,
flexibility, deployment, collector equipment, propellant depletion and a
resupply interval remain engineering obligations. The assumed fixed 20% is
a sensitivity allocation, not a demonstrated bill of materials.

The coarse second passage completes both six-hour stages with zero electric
translation, sampled local coverage and **1.720 km** minimum surface clearance.
Refitting the four return arcs to its actual hour-12 states costs **1.85801 m/s**
of proposed translation, with a 3.22909 m/s member maximum. Its 1,171 nonzero
scheduled burns have a maximum acceleration of **0.0005027 m/s²**. The target
keeps the selected relative packing and velocity pattern but follows the
new freely evolving mean. This refit is an endpoint proposal, not an accepted
return trajectory.

## Whole-sequence replay

The tighter calculation starts at the original epoch and passes each computed
terminal state directly into the next stage. It doubles force quadrature from
8 to 16 solar sources, halves the maximum integration step from 120 to 60 s,
and tightens relative tolerance from 2e-10 to 2e-12. It completes the hardware-
loaded six-hour service and six-hour departure with **zero translation**.

| Executed stage | Closest sampled physical surfaces | Conditional two-second interval bound | Maximum replay position difference, including reconstruction |
|---|---:|---:|---:|
| Service, hours 0–6 | 6,698.74 m | 6,669.95 m | 1.555 m |
| Departure, hours 6–12 | **1,719.63 m** | **1,683.84 m** | **6.994 m** |

Both stages retain positive reflected-beam margins. The largest sampled
attitude rate is **0.01172 degree/s**, below the 0.1 degree/s cap. All 49
service dates and 25 solar sources protect the complete assigned circular
receiver inside the target. Interval guards remain conditional on their
acceleration/rate envelopes; these tests do not certify continuous coverage
of the four-radius target.

A further **891 source/date probes** seek poorly covered directions between
those service samples: 533 grid probes followed by four bounded searches.
All cover the entire circular receiver. The smallest additional coverage
margin is **6,977.41 m**, including a 0.188 m circle-sagitta allowance, at
service time 745.03 s near the solar limb. This is a passed local adversarial
search, not an exhaustive source/time certificate.

The replay removes empty receiving-plane rectangles before taking the same
exact shadow union. This reduces cost where projected neighbours are numerous
but only a few actually shadow each tile. Comparisons with the original
backend at all three coarse stage endpoints differ by at most **3.47e-18 m/s²**
in total acceleration and **1.12e-16** in visibility. Separate tests compare
the compact union with unpruned all-pair geometry, including full occlusion
and near-edge-on orientations. The change accelerates the calculation without
introducing a different radiation-pressure model.

The coupled return stops at the **100 m surface-clearance floor** at
**hour 13.893535**, or **1 h 53 min 37 s after hour 12**. The tighter run finds
the same pair, tile indices **325 and 347** (zero-based), only **0.6523 s**
later than the coarse event. Maximum position disagreement, including raw
state reconstruction, is **9.04188 m**, inside the 50 m replay tolerance.
Their centres are about 12.87 km apart, but their body-frame displacement is
(-9,314.92, 8,885.43, 100.00) m: the parallel squares overlap in both edge
directions and approach along their common normal.

By this stop, all 361 first burns have begun and mean applied translation is
**0.97590 m/s**. Beam margin remains **44.68 degrees**, attitude rate remains
below the cap, and the thrust cap is not limiting. This is a clearance failure
of the proposed continuation, not a demonstrated physical collision: propagation
stops before impact. The all-pair two-second guard flags 50 intervals near
the event, with a 64.24 m conditional lower bound. That bound is not an
observed 64.24 m separation and does not certify the final intervals safe.
The earlier proposal intersections are retained separately as predictions.

The additional hardware mass therefore does not invalidate the demonstrated
low-translation service/departure, but the selected return neither closes
the geometry nor stays within the earlier 0.53987 m/s remaining allowance.
No second service or handover is propagated after this verified obstruction.

At the frozen terminal positions, changing the common Sun tilt to 0, 15,
30, 45 or 60 degrees produces minimum surface distances between **0.277 and
55.62 m**. Those orientations do not repair this snapshot. Rolling the
already feathered squares about their unchanged normal is more promising:
45 degrees gives only 104.05 m, but **60 and 75 degrees give 719.41 and
765.61 m**. The limiting pair changes. This is a geometric diagnostic;
the slew, its earlier swept clearance, changed mutual shadows, torque cost
and subsequent return remain unexecuted. A late instantaneous rotation is
not an available maneuver.

## Executed energy and optical inventory

The account integrates ideal edge-couple attitude propulsion and each
independent translation arc, with actual per-member mass, gravity torque and
reflected-band shadow moments. Service/departure contain no translation burns;
one four-hour departure turn is executed. Only the first 361 of 1,171 proposed
return burns start before the stop, and the eight-hour arrival turn is never
reached. A recurring maneuver frequency has not been measured.

| Executed stage | Propulsion energy, TJ | Mean first interception, bolometric equivalent, TW | Mean redirected band, optical TW |
|---|---:|---:|---:|
| Service, 6 h | 2.648 | 39.571 | 5.459 |
| Departure, 6 h | 49.211 | 17.057 | 2.353 |
| Return before clearance stop, 1.894 h | **67.479** | **0.0920–0.0950** | **0.01269–0.01310** |

The twelve-hour passage uses **51.859 TJ**, 25.10% more than the optical-only
41.455 TJ passage. Its mass-weighted equivalent impulse is 0.75609 m/s;
normalizing its actual energy to the original optical mass instead gives
**0.94807 m/s**. These normalizations differ because the carried mass differs.
Through the failed-return stop, measured energy totals **119.338 TJ**.
The corresponding ideal exhaust account is **185,637 kg**, or 0.00820% of
carried mass; depletion was not fed back into the constant-mass trajectory.

The largest 120-to-240-second energy change is 3.54e-5 relative. Comparing
16 and rotated 32 sources at matching dates changes energy by at most
1.83e-4 relative. Individual peak requirements, taking the larger value from
both evaluations, sum to **18.762 GW**. Every member fits its allocated power
with the full 25% installation margin; the smallest installed-capacity-to-peak
ratio is **1.3911**. This checks the **29.190 GW installed allocation** over
the executed stages, not the unexecuted rest of the return. It also does not
establish an energy supply, storage system or detailed thruster layout.

All 361 members retain separate time-resolved optical records, on both faces,
through every executed stage. The preliminary prefix-only hardware trial has
its own ledger and is not counted as additional inventory. Twelve-hour first
interception is **1.22315e18 J**, just 0.0165% below the optical-only passage.
Service/departure source refinement changes optical energy by at most 9.72e-5
relative. Near feathering, however, the return's 16/32-source discrepancy is
**3.14%**; the table shows both estimates, not a certified error bound. Finer
grazing-source quadrature remains a collection-account uncertainty.

These numbers are optical interception and redirection potential. The model
reflects its designated spectral band and transmits other light; installing
absorbers or capturing outgoing beams requires their momentum, conversion,
mass, storage and distribution to be modelled. No electrical collection
capacity or power delivered to the Moon is inferred. The **100 TW delivered-power
ambition** remains separate from propulsion consumption and installed
propulsion capacity.

## Conditional inventory trade

The following arithmetic repeats **only the measured expenditure through
the failed-return stop** once every two days, retaining the propagated mass
allocation. It charges no completed return, handover, recovery or global
coverage, and therefore does not describe an operational fleet.

| Assumed inventory | Measured-prefix expenditure per two days, mean TW | Installed propulsion capacity from this allocation, TW |
|---|---:|---:|
| 17.29 million | **33.07** | 1,397.67 |
| 25 million | **47.83** | 2,021.44 |
| 35 million | **66.96** | 2,830.01 |

Even at the optimistic capacity-relaxation inventory, this branch has already
spent more than the **23.835 TW** goal would allow under that recurrence.
Installed peak capacity is much larger than mean expenditure because the
burns are intermittent; it is not collection capacity or continuous demand.

A separate proposal-only comparison uses the earlier scalar 20% fixed mass,
10% thrust duty and power-hardware feedback scenario. The 2.94699 m/s nominal
case-6 proposal, repeated every 1.92832 days, gives **56.89, 82.28 and
115.20 TW** at the same three inventories: about **24%, 35% and 48%** of the
238.35 TW benchmark, before unpaid actual return loads, encounter repairs,
handover and recovery. These figures show the middle-ground cost scale of
this restricted proposal; they do not establish feasibility or bound another
family. The scalar `closed` field in these scenario rows refers only to its
mass/power algebra. All return and fleet acceptance fields remain false.

## Comparison and acceptance boundary

The **238.35 TW** benchmark is the moving-distance held screen in `report.md`:
50 g/m² base optics, independently held area quadrature, relaxed ideal optical
authority, finite-Sun coverage, 30 km/s exhaust, 70% electrical-to-jet
efficiency, 45-degree exhaust cant, 300 W/kg power equipment and a 25% peak
allowance. It includes feedback from power equipment and a seven-day
propellant buffer; storage is excluded. The new fixed 20% allocation and
constant carried mass are a different engineering sensitivity, not equivalent
hardware closure. Comparing scenario powers does not remove that difference.

Global coverage and handover are gated on a local return that survives its
own geometry and cost checks. The 17.29-million population remains an
area-capacity relaxation. A moving 20 km receiver protected by 361 tiles for
six hours cannot simply be copied across the four-radius target: actual
ephemeris phases, shared finite-Sun rays, pattern boundaries, cross-pattern
shadows and encounters must be evaluated together. No second service,
722-tile handover, acquisition or recovery is credited in this checkpoint.

The present return family restores one specified grid and velocity field,
keeps a common feathered orientation during early return and offers only two
identity assignments. Infeasible local separating cuts do not exclude other
escape directions. A useful subsequent search would change initial velocity,
coast attitude, return timing and the next service-ray assignment together,
allowing a different pattern on each passage. It should retain the measured
hardware-loaded passage as a comparator and test a recurring cost allowance
before expanding the simultaneous fleet calculation. That experiment has not
been executed here.

## Resources, checks and reproduction

This continuation stops at the verified local feasibility/cost checkpoint.
The seven completed producers take **4,389.34 s (73.16 min)**, within the
4,800 s ceiling. Maximum resident memory is **644.1 MiB** under the 2 GiB
address-space cap. New raw products occupy **202.85 MiB**, below 768 MiB.
All numerical producers run sequentially with one numerical thread; this
continuation has no numerical restarts. The additional hardware-passage
decision and accounting-grid choice were recorded before their execution.

The products preserve source, constant, parent-product and raw-data hashes.
Reproduction requires the existing source-bound packing and preceding raw
products described in their reports. With those inputs present, run these
stages sequentially with `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`:

```sh
python -m research.studies.solar_shield_array.closure_screen
python -m research.studies.solar_shield_array.closure_hardware
python -m research.studies.solar_shield_array.closure_sizing
python -m research.studies.solar_shield_array.closure_return
python -m research.studies.solar_shield_array.closure_validate
python -m research.studies.solar_shield_array.closure_budget
python -m research.studies.solar_shield_array.closure_audit
python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```

The targeted suite passes **141 tests in 12.32 s**, including analytic endpoint
control, individual mass/power feedback, unpruned shadow comparisons,
partial-burn energy accounting and all seven products' source/parent identities.
`make check` passes the layer check (**524 files, zero violations**) and
reports **686 Python passes, 55 skips and the same 13 existing failures in
30.60 s**. These concern missing pinned WHI spectral inputs, missing generated
climate configuration and the environment's process lookup. Later Makefile
targets are not reached. Checks have their own five-minute allowance;
`checks.json` records exact executed checks and unperformed work.

The research branch receives the tested code, compact products and this report.
Main is unchanged. No complete return, two-pattern handover, fleet operating
cost or continuous global coverage is accepted at this checkpoint.
