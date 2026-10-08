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

## Executed screen and coupled candidate

The early plan was pushed as `2b54ebbfa56d379b4ccf35aa759653182496c6f8`.
The inexpensive screen executes all six arrangements and both retained turns.
Every arrangement passes the thirteen-date affine service-coverage proposal.
Changing row/column parity swaps which neighbouring tile occupies each depth
class; it preserves the projected starting grid and depth population.

The proposal model carries the measured force defects from the previous
failed departures, transports their remaining tails gravitationally and
adds the response to the changed initial state. This permits reuse of
expensive force evaluations for screening. Shadow forces and their changes
under a different arrangement are subsequently recomputed during propagation.

| Starting depth step | Depth order | Four-hour 30-degree turn | Four-hour 90-degree turn |
|---|---|---|---|
| 4 km | Original | Local correction fails | Local correction fails |
| 4 km | Swapped | Selected half-spaces infeasible | Selected half-spaces infeasible |
| 8 km | Original | Local correction fails | Local correction fails |
| 8 km | Swapped | Local correction fails | Local correction fails |
| 12 km | Original | Local correction fails | Proposal passes with zero translation |
| 12 km | Swapped | Local correction fails | Local correction fails |

The eleven failures describe the selected local half-space branches and
their limited refinement. Independent controls change more nearby pair
distances than the old three-mode controls. Several initial LP solutions
therefore create further encounters; the following constraint additions are
infeasible in those branches. These cases do not establish an energy floor
or exclude the starting arrangements under other schedules and branches.

The surviving 12 km arrangement retains the original parity order. Its four
depth levels span **36 km**, compared with the earlier 12 km total depth.
The projected pitch, clear width, tile count and base optical mass are
unchanged. Its approximate departure has 1.883 km minimum sampled separating
projection and no detected coplanar intersections. The LP selects zero for
all 2,166 translational control coefficients, with 0.68243 m/s nominal
attitude impulse for the four-hour turn.

The new service then passes coupled propagation with dynamic finite-source
shadow forces and finite-area gravity. All 49 tested service dates, each with
the solar centre, sixteen rotated limb directions and eight interior sources,
retain full coverage of the same moving receiver window. The window remains
inside the four-lunar-radius target. Electric translation remains zero.

From those actual service endpoints, another LP leaves every departure burn
at zero. Coupled propagation completes the six-hour departure, including the
turn and one hour of feathered coast. Its maximum departure position defect
from the revised proposal is **631.32 m**. Its two-second swept guard flags
no pair. Beam, attitude-rate and total acceleration checks also pass.

| Epoch-relative interval | Executed command |
|---|---|
| 0–6 h | Sun-facing natural service; zero translation |
| 6–7 h | Sun-facing natural departure; zero translation |
| 7–11 h | Smooth 90-degree turn about the first square edge; zero translation |
| 11–12 h | Feathered natural coast; zero translation |

All times are measured from 4 October 2026 TDB using the same DE440 ephemeris.
The service assignment ends at hour 6. Continued optical interception during
departure does not assign coverage to a replacement fleet member.

## Independent replay and geometry

The complete twelve-hour candidate is replayed with sixteen rather than eight
solar force directions, 60 s rather than 120 s maximum integration steps, and
relative tolerance 2e-12 rather than 2e-10. The departure begins at the tighter
service terminal states. The largest saved-interpolant position differences
are **1.995 m** during service and **8.988 m** during departure, including the
measured reconstruction allowances. Both meet the 50 m gate. The nonlinear
replay retains zero translational thrust throughout.

An independent distance calculation evaluates the Euclidean surface distance
between physical 10 km squares every two seconds. Its spatial-tree culling
uses a proven radius from an existing distance upper bound; it therefore
includes every pair that could improve the global minimum. An independent
polygon-distance test checks the square formula and culling. The interval
guard bounds relative centre motion, corner motion and interpolation error.

| Stage | Closest sampled surfaces | Conditional all-pair interval bound |
|---|---:|---:|
| Six-hour service | 6,698.57 m | 6,669.77 m |
| Six-hour departure | 1,760.43 m | 1,724.61 m |

The departure witness is pair 319/322 at hour 12, using zero-based member
indices. No interval falls below 100 m. The bounds assume acceleration no
greater than 0.15 m/s² and attitude rate no greater than 0.1 degree/s;
the largest evaluated values are 0.022035 m/s² and 0.011719 degree/s.
Angular acceleration peaks at 2.506e-6 degree/s². These generous bounds and
the source/integrator comparison support the local result, but do not prove
the continuous nonlinear dynamics or flexible-panel shape.

The earlier service guard's 100.25 m figure is a deliberately loose
circumsphere lower bound for an excluded pair, not a measured surface
approach. The new all-pair distance audit resolves that ambiguity. The
minimum sampled reflected-beam exclusion margins are 96.47 degrees during
service and 30.42 degrees during departure. They include the finite Sun and
Earth/Moon exclusion geometry; no near-exclusion beam is accepted.

The 49-date, 25-source service-coverage check passes again. A separate ray
adversary uses thirteen dates with 32 rotated limb directions, eight interior
directions and the solar centre, followed by four bounded local searches over
time and solar-disk coordinates. All **887** evaluations cover the entire
moving receiver disk. The worst found boundary margin is **6,977.42 m** at
744.41 s and a solar-limb source. This subtracts the 0.188 m circle-polygon
sagitta, so the positive margin covers the exact receiver circle for each
tested source and date. The search does not certify all intervening solar
directions and times. It also does not tile the four-radius target with
simultaneous patterns or verify their interfaces.

## Calculated control energy

The selected candidate's 2,166 independent burn coefficients remain zero.
The improvement comes from starting depth and the slower turn, not an
executed independent-burn solution. Each tile still needs continuous attitude
control against calculated gravity-gradient and reflected-band shadow torque.
The following comparison uses the same uniform rigid optical mass and ideal
electric edge-couple actuator as the previous actual-load account in
`energy.md`; the older 6.70269 m/s disturbance-envelope allowance is not the
comparison baseline.

| Twelve-hour account | Previous verified arrangement | Deeper arrangement |
|---|---:|---:|
| Service attitude, equivalent m/s | 0.03794 | 0.04090 |
| Departure translation, m/s | 1.40584 | 0 |
| Departure attitude, equivalent m/s | 2.79604 | 0.71696 |
| Total mean equivalent impulse, m/s | **4.23982** | **0.75786** |
| Propulsion energy for 361 tiles, TJ | 231.917 | 41.455 |

The local reduction is **82.13%**, exceeding the pilot's factor-of-two goal.
The turn's nominal attitude impulse falls from 2.72711 to 0.68243 m/s;
finite-area gravity and radiation moments give the actual actuator load.
The new service cost rises slightly because the changed shadows change torque.
Electric propulsion energy follows the stated 30 km/s exhaust, 70% efficiency
and 45-degree cant. Equivalent impulse includes opposing attitude-couple
forces; it is not centre-of-mass velocity change.

Loads and collection use sixteen solar directions at 60 s intervals, with
32 rotated directions at 120 s for refinement. Command torques and integrated
actuator costs use two-second dates. Time coarsening changes service/departure
attitude impulse by at most 8.65e-6 relative; the matched-time source check by
at most 2.33e-4. Independent optical accounting agrees within 3.54e-16 relative.
Three- versus five-node-per-axis gravity quadrature differs by at most
2.29e-14 N m/kg. No calculated disturbance exceeds the previous envelope.

Summing every member's individual peak demand gives **4.132 GW**, rather than
sizing all vehicles from the simultaneous fleet peak alone. At 300 W/kg plus
25% allowance this implies **17.22 million kg**, or **0.954%** of the pattern's
1.805-billion-kg optical mass. This is provisional power hardware only. Its
mass distribution, thrusters, conversion/collection equipment, storage,
propellant, structure and flexibility have not been added to the propagation.
The commanded rigid attitude is kinematically admissible; realization by an
engineered 10 km structure is unresolved.

Initial states are not free operational preparation. Relative to an initially
Sun-tracking lattice translating with the reference, their velocity difference
averages **0.998 m/s**, with a 1.745 m/s maximum. This is a state comparison,
not an executed acquisition maneuver. Establishing the deeper positions and
reacquiring these velocities on return remain unpaid. No recurrence or
maneuver frequency has been inferred from the twelve-hour prefix.

## Optical inventory and fleet comparison

Every one of the 361 physical members retains time-resolved optical power in
the saved inventory, including departure after its service assignment ends.
The same first-interception shadow geometry drives reflected-band force,
torque and collection accounting. Body eclipses are excluded throughout the
propagation. The ideal model transmits other wavelengths and includes no
additional absorption or multiple reflections. Installing collectors or
capturing outgoing beams would require their own momentum and mass account.

| Mean over the executed stage | Service, TW | Departure, TW |
|---|---:|---:|
| Projected aperture before mutual shadows | 48.053 | 26.888 |
| First interception, bolometric equivalent | 39.571 | 17.066 |
| Mutual-shadow loss, bolometric equivalent | 8.482 | 9.822 |
| Redirected spectral band, optical power | 5.459 | 2.354 |

The twelve-hour first-interception integral is 1.22335e18 J, about **6.00%**
above the previous verified twelve-hour schedule. The slower turn keeps more
projected aperture during departure; service interception changes by about
-0.010%. Source refinement changes the optical energy by at most 9.79e-5
relative. These are optical potentials of the simulated local inventory,
not electrical generation, available fleet collection capacity or power
delivered to the Moon. The **100 TW delivered-power ambition** remains separate.

The following is only a conditional trade. Assume this measured prefix
recurs once every two days, with the earlier scalar 20% additional fixed mass,
10% thrust duty and power-hardware feedback scenario. The remaining allowance
is for everything still unexecuted: return, handover, repacking and recovery.

| Assumed inventory | Prefix-only propulsion, TW | Entire two-day allowance at 23.835 TW, m/s | Allowance left after this prefix, m/s |
|---|---:|---:|---:|
| 17.29 million | 13.86 | 1.29773 | **0.53987** |
| 25 million | 20.05 | 0.89997 | **0.14211** |
| 35 million | 28.07 | 0.64407 | **-0.11379** |

The 17.29-million figure is a capacity relaxation, not a sufficient deployed
inventory. The 13.86 TW entry therefore cannot be claimed as an achieved
fleet operating cost. Even its installed propulsion power in the assumed
10% duty scenario is 173.29 TW, distinct from mean consumption and collection.
The held-screen benchmark remains **238.35 TW** under the assumptions stated
above. This local improvement leaves some room below its 10% target at the
optimistic inventory, but only a complete repeatable cycle can establish a
fleet saving. More inventory consumes that remaining allowance quickly.

## Resources, checks and reproduction

The search stops at this first positive local cost checkpoint. Twelve cheap
proposals, one changed service propagation, one coupled departure and one
complete tighter twelve-hour replay were executed. The five producers took
2,778.20 s; charging 120 s for an interrupted partial replay gives
**2,898.20 s**, within the 3,600 s numerical allowance. Maximum resident memory
was **685.1 MiB**, and new raw products occupy **136.03 MiB**.

The first tighter-replay attempt had an internal 1,200 s alarm that measured
progress showed would be insufficient. Its allowance was raised within the
existing total budget before the complete rerun. An initial process signal
failed because the managed process identifier was not visible to that shell;
the session interrupt then stopped it. The restart briefly overlapped the
partial run, conservatively charged as at most 30 s of two processes. Both
used one numerical thread. This was a departure from the one-process plan;
all subsequent numerical stages ran sequentially. The partial work is charged
at 120 s, larger than the provisional 87 s recorded in the replay product.

Products retain source, constant and parent hashes. Reproduction requires the
existing source-bound pattern, departure and energy raw products; their
generation is described in the preceding reports. With those inputs present,
run these producers sequentially with `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1`:

```sh
python -m research.studies.solar_shield_array.packing_screen
python -m research.studies.solar_shield_array.packing_run
python -m research.studies.solar_shield_array.packing_validate
python -m research.studies.solar_shield_array.packing_budget
python -m research.studies.solar_shield_array.packing_audit
python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```

The targeted dynamics/study suite passes **130 tests** in 8.70 s, including
independent response-map, control, exact-square-distance, ray-margin and
source-identity checks. `make check` passes the layer check (508 files, zero
violations) and reports **675 Python passes, 55 skips and the same 13 existing
failures** in 29.06 s. Those failures concern missing WHI spectral inputs,
missing generated climate configuration and the environment's process lookup.
Later Makefile targets were not reached. Details and the publication target
are recorded in `checks.json`.

No full return, two-pattern handover, recovery cycle, global placement,
hardware-mass feedback, additional optical physics or continuous fleet
certificate was executed. The next useful bounded experiment starts from
these actual twelve-hour terminal states and tries to close return and
handover within the remaining impulse allowance, preserving the new energy
gain. The unsuccessful local LP branches remain useful search evidence.
