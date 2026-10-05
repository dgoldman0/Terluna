# Starting arrangements and independent maneuvers

## Bounded pilot

The remote head was checked as `8bb69c204ffb0ff39b9ec3c39da848a22fc13156`
before work on `research/solar-shield-habitat-array`. Main remains unchanged.
The author authorized changing the starting arrangement and allowing more
independent maneuvers after the partial-turn failures in `energy.md`.
Energy remains the objective; retain the existing 100 m clearance constraint.

The first question is whether a different service packing or independent
early corrections can make a substantially cheaper departure feasible.
Use the same actual-date reference passage and six-hour service obligation.
Keep all 361 physical 10 km squares, the 50 g/m² base optical allocation,
8.5 km projected pitch and the moving 20 km receiver window inside the
four-lunar-radius target. Compare four-level depth steps of 4, 8 and 12 km,
with both row/column parity assignments. The initial relative velocity field
comes from the existing six-hour shooting proposal at scale 0.95. Starting
states are design variables; their insertion and subsequent reacquisition
are explicit unpaid obligations until a complete cycle is demonstrated.

Test the four-hour 30-degree and four-hour 90-degree departure schedules
retained in `energy.md`. Select according to coverage and impulse. Allow two
non-overlapping two-hour smooth burns, starting at departure hours 0 and 3.
Every tile may choose its own three-component vector in each burn. The
resulting 2,166 control coefficients replace the three shared coefficients.
A sparse linear program minimizes the sum of absolute components as a fuel
surrogate; always report the sum of Euclidean vector norms and actual peak
acceleration. Conservative component bounds may approximate the spherical
0.001 m/s² acceleration cap from within. Retain the 0.1 degree/s attitude
cap and finite-Sun reflected-beam exclusions.

Begin with inexpensive trajectory/coverage proposals. Each changed packing
must then survive individual coupled service propagation, including finite-
area gravity and dynamic mutual-shadow sail force. Departure corrections use
gravity response maps around measured trajectories, with bounded updates and
encounter constraints. Solar-force changes remain in the nonlinear replay;
the response approximation cannot accept a candidate. Retain failed local
half-space choices and force defects. A solver failure is local evidence.

Search for encounters between output dates and refine their constraints.
Evaluate service rays at actual simultaneous ephemeris dates, including
rotated solar-limb and interior sources. Use the existing dynamic blocker
search rather than assuming the original neighbour stencil still applies to
a deeper or reordered packing. Reject body eclipses until joint body/array
visibility is implemented. Coverage assignments end when the service is
released: handover to a second pattern, complete return, recovery and
recurring cadence remain separate requirements. No released ray receives
coverage credit.

The numerical budget is **3,600 s wall time**, one calculation process and
one numerical thread, **2 GiB address space** and **512 MiB raw output**.
Allow six initial arrangements, twelve cheap arrangement/schedule screens,
at most five new six-hour coupled service passages (reuse the verified
baseline), four coupled departure trials or corrections in total, and one
independent tighter replay of the selected twelve-hour service/departure.
Allow at most four linear correction steps per departure proposal and
reserve 600 s for torque/collection accounting. Repository checks receive
a separate five-minute allowance. Check usage between stages; do not enlarge
these counts without recording the reason before execution.

A positive local checkpoint must preserve service coverage, complete the
six-hour departure, pass the sampled and swept geometry/beam/rate gates,
and agree with an independent tighter replay within 50 m. A calculated
service-plus-departure account below **2.11991 m/s equivalent impulse**,
half the previous 4.23982 m/s account, is the initial energy goal. The
optimistic 17.29-million-tile inventory allows only about 1.30 m/s for a
hypothetical complete two-day cycle at the 23.835 TW target. Preserve the
inventory trade and unpaid return/handover rather than assigning a recurrence.

Track each simulated tile's optical interception and energy through every
executed stage. Charge translation and calculated rigid gravity/radiation
torque, and record maneuver counts and provisional peak power-system mass.
Added hardware mass that has not been propagated remains a closure gap.
Keep optical collection, assumed electrical generation, propulsion demand
and the 100 TW delivered-power ambition separate. The 238.35 TW held-screen
comparison retains 30 km/s exhaust, 70% efficiency, 45-degree thrust cant,
300 W/kg power hardware, 25% peak allowance and seven days of propellant;
storage and complete vehicle engineering are excluded.

Stop at the first informative feasible local cost result or limiting
constraint. Continuous fleet coverage, operational recurrence and structural
realizability need their own evidence beyond this pilot.

Before the early checkpoint, `make check` passed the layer check (496 files,
zero violations) and reported 665 Python passes, 55 skips and the same 13
existing failures in 32.48 s. Later targets were not reached.
