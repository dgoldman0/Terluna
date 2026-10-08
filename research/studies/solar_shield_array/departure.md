# Opening a safe departure path

**A prepared local departure now passes the coupled numerical gates.** All
361 tiles complete a turn and three-hour free coast, with a 107.71 m
conditional swept clearance bound against 100 m and 0.59 m independent
replay disagreement. This resolves the tested early packing obstruction.
It does not establish a handover, full return or functioning fleet.

Translation costs 1.41 m/s per tile on average. An explicit ideal electric
torque-couple option raises the departure account to **5.04 m/s equivalent
impulse**, including a conservative disturbance allowance. Attitude actuation
and its added mass are now important cost variables; the translational result
alone cannot stand for the operating budget.

## Bounded experiment

This continuation starts from remote
`4b82a8386a1d4b207630f0181bdbbd83f6469bd7`, checked before work on
`research/solar-shield-habitat-array`. Main stays unchanged. All 361 actual
terminal states from the tighter six-hour local service replay remain the
initial conditions. The optical allocation is 50 g/m² and each physical tile
is a 10 km square with a 9.89 km clear side.

The previous return failed during its early attitude transition. At its
100 m clearance stop, the limiting centres were 8.4 km apart; restoring the
Sun-facing attitude at those same positions left 2.6 km of clearance. The
next objective is to find a dynamically reachable path out of that packing,
and measure its impulse before attempting another full return. A feasible
departure would remove one specific obstruction, without establishing a
handover, recurrence or a fleet solution.

Decision variables are the start and duration of relative-motion preparation,
temporary layer spacing, the direction and timing of a common attitude
transition, and smooth maneuver coefficients. Keep the reference centre on a
freely evolving ephemeris trajectory. Permit the local pattern to deform;
do not require circular Sun-following traffic. Begin with inexpensive finite-
square geometry and gravity-variational screens, then propagate promising
members together with sail force and dynamic mutual-shadow unions. Use the
offending pairs to refine intermediate geometry and control timing. Preserve
rejected cases and the constraint responsible.

The bounds are 100 m finite-square clearance, 0.001 m/s² translational
acceleration and 0.1°/s commanded attitude rate, with the established finite-
Sun reflected-beam exclusion. Physical tile edges set encounter constraints;
clear optical edges set intercepted light and sail force. Verify any shared-
attitude or low-dimensional control reduction against all 361 propagated
states. The available coupled model accepts uneclipsed prefixes only; reject
or explicitly implement and test any proposed extension that needs joint
body/array eclipses. No partial-eclipse force substitution is allowed.

Allow 60 minutes of numerical wall time, one calculation process, one
BLAS/OpenMP/solver thread, 2 GiB address space and 1 GiB raw output. Within
that budget, allow up to 96 trajectory/attitude proposals, six local control
refinements, four coupled prefixes of at most six hours each, and one tighter
replay. Static geometry scans may use up to 10,000 angle/spacing samples.
Repository checks have a separate five-minute allowance. Check resources
between stages and save intermediate products.

A useful positive checkpoint is a coupled departure that completes its
attitude transition, stays clear for a further coast interval, passes the
beam/rate/acceleration constraints and agrees with a tighter replay within
50 m. Search for missed encounters between snapshots and state the remaining
continuous-certification gap. A mean departure impulse below 3 m/s is a
screening target, leaving an explicitly unresolved allowance for return,
repacking, handover and recovery. It is not a recurring-cost acceptance limit.
If only more costly candidates work, retain the trade and its limiting
constraint. Stop at the first meaningful feasibility or cost checkpoint.

Coverage remains attached to actual receiver rays and dates. During any
continued Sun-facing service, evaluate the assigned local window with the
finite Sun; changing attitude releases that assignment and requires a second
pattern before it can become an operational handover. Do not count released
rays as protected or infer global coverage from a departure alone. A later
handover optimization must include both patterns' forces and encounters.

Keep propulsion energy, collection capacity and the 100 TW delivered-power
ambition separate. The comparison remains the 238.35 TW annual ideal-optical
held screen: 50 g/m², 30 km/s exhaust, 70% efficiency, 45° thrust cant,
300 W/kg power hardware, 25% peak allowance and seven days of propellant,
excluding storage and full vehicle engineering. The conditional inventory/
recurring-impulse table in `patterns.md` applies only after a repeating,
coverage-feasible schedule has actually been measured.

## Execution record

The plan was pushed as `cc7446118a8f3e1dc59d47274db128548daa2b3e` before
numerical search. Its `make check` run passes
the layer check and reports 639 Python passes, 55 skips and the same 13
previously recorded failures (missing spectral inputs, CM1 configuration,
process lookup, ring comfort and the dependent protection-architecture
check). Later targets are not reached.

## Geometry and control refinement

The static screen tests 20 spacing/axis/sign combinations at 91 angles each.
It uses the actual terminal positions, with a hypothetical displacement of
each of the four depth classes. Increasing the separation between successive
classes by 4 km opens a sampled 90° path about one square edge, with at least
731 m clearance for either turn direction. The other edge needs a larger
expansion. This establishes a geometric possibility, without pretending that
the displaced states have been dynamically reached.

The trajectory screen then tries 96 combinations: one- or two-hour expansion
burns, attitude delays of two, three or four hours, half- or one-hour turns,
four impulse coefficients and the two signs of the selected edge axis. All
use a freely evolving reference and gravity variational response from the
actual individual terminal states. The reference's unshadowed Sun-facing
sail is an explicit approximation in this proposal stage.

Twenty-four choices exceed the 0.001 m/s² acceleration cap. Every remaining
depth-only choice fails clearance; 36 also fail the reflected-beam exclusion.
The retained witnesses reveal a second obstruction: rows in the same depth
class focus together. In the diagnostic one-hour burn / four-hour delay case,
one such two-row spacing falls from about 14.0 km two hours after service to
10.19 km at four hours and 5.65 km at six hours. A depth-class impulse moves
both members alike and cannot repair this focusing.

The reserved refinement adds two transverse velocity modes. The three
parameters multiply the centred depth-class label, the normalized column
coordinate and the normalized row coordinate. Directions are the initial
solar frame; coefficients sum to zero over the 361 equal-mass members, so the
command does not prescribe a reference-centre acceleration. Tile positions
and attitudes are not held by formation feedback.

The selected schedule starts a two-hour sin² burn immediately at the saved
service endpoint. It then turns through 90° about one edge over an hour,
using 60 s orientation knots from a quintic relative-angle profile, and
coasts for another three hours. The rotation spline's actual rates and
accelerations are checked; its tiny interpolation residuals are included in
the attitude budget. The opposite sign directs some reflected
rays toward the excluded region and is rejected by the screen.

| Time after the six-hour service | Command |
|---|---|
| 0–2 h | Sun-facing; smooth layer, row and column preparation |
| 2–3 h | No translation thrust; one-hour turn toward feathering |
| 3–6 h | Feathered free coast; no translation thrust |

One convex correction selects 2,160 finite-square separating half-spaces
along the gravity proposal. The objective is the mean Euclidean norm of the
361 impulse vectors. The exact vector acceleration cap is a convex quadratic
constraint; no component-box approximation is used here. The three optimized
coefficients are **1.04549, 0.73873 and 1.29455 m/s**, in depth, column and row
order. They are mode coefficients, not three sequential burns.

The resulting mean translational impulse is **1.40584 m/s**. Peak commanded
acceleration is **0.000588 m/s²**. Thirty-second refinement finds at least
**389.0 m** separating clearance in the proposal, no detected coplanar
finite-square impacts and positive beam/rate margins. This is a local optimum
within the chosen half-spaces and three modes. The 96 failed cases, initial
guess, correction residuals and crossing witnesses remain in the products.

## Coupled propagation and the interval guard

The first coupled propagation completes all six hours with all 361 states,
finite-area gravity and dynamic finite-Sun mutual shadow forces. No electric
translation is applied after hour 2. It uses eight solar-force directions,
a 120 s maximum integration step and relative tolerance 2×10⁻¹⁰.

Its five-second geometry samples remain above 100 m, with a minimum of
133.56 m. The swept bound initially flags 27 possible pairs because its
0.1°/s attitude allowance consumes much of that narrow margin over five
seconds. This is a failed conservative exclusion, not 27 simulated impacts.
The original result is retained in `departure.json`.

Refining the same stored trajectory to two-second intervals removes those
flags: the minimum sampled clearance is 132.46 m and the conditional swept
lower bound is 107.73 m. This refinement uses a 30 s position/velocity Hermite
reconstruction. A separate 16-source, 60 s maximum-step propagation with
relative tolerance 2×10⁻¹² supplies the independent dynamics and interpolation
check in `departure_validation.json`.

That independent replay completes the whole prefix. The closest separating-
axis witness is tiles **186 and 187**, 8,968 s after the service endpoint,
during the turn. These quantities are conservative axis-based separation
bounds, not the exact minimum Euclidean distance between membranes.

| Validation quantity | Main run, refined geometry | Tighter replay |
|---|---:|---:|
| Solar-force directions | 8 | 16 |
| Maximum integration step | 120 s | 60 s |
| Geometry interval | 2 s | 2 s |
| Lowest sampled separating-axis clearance | 132.457 m | 132.440 m |
| Conditional swept lower bound | 107.728 m | 107.711 m |
| Possible pairs left by the swept exclusion | 0 | 0 |
| Detected finite-square intersections | 0 | 0 |

At two-second comparison dates, the largest position disagreement is
**0.5884 m**, and the largest velocity disagreement is 3.65×10⁻⁵ m/s.
The tighter trace's 30 s Hermite reconstruction differs from its dense
integration output by at most 1.75×10⁻⁶ m. The minimum sampled reflected-beam
margin is **20.310°**. The rate peaks at **0.046864°/s** and angular acceleration
at 4.01×10⁻⁵°/s². These are actual spline commands, including ephemeris motion.

The gravity proposal itself differs from the tighter dynamics by as much as
4.45 km in position, including 3.25 km of reference-centre displacement and
1.45 km of relative-pattern discrepancy. Its approximate force history is
unsuitable for accepting placement. The accepted local geometry belongs to
the coupled integration, with its independent replay, rather than to that
proposal.

The swept calculation uses the physical 10 km square, not the 9.89 km clear
aperture. Its centre-acceleration allowance is 0.15 m/s² and its orientation
allowance is 0.1°/s; the first coupled run's sampled maxima are 0.02246 m/s²
and 0.04686°/s respectively. These conditional derivative bounds and numerical
replays do not certify a perturbed, flexible structure or the full continuous
force model. The margin above 100 m is small enough to require further robust
design before treating this as an operational maneuver.

## Coverage assignment and actuation cost

The finite-Sun coverage check follows the same moving 20 km receiver disk.
At 0, 5, 10 and 15 minutes after the old service endpoint, that disk remains
inside the four-radius target and is fully covered for the tested source
centre, 16 limb directions and eight interior directions. By the 20-minute
sample its centre is 6,941.95 km from the target centre; the extra 10 km window
radius puts part of it outside the target. Later apparent coverage of that
moving disk cannot count as service of the target. No replacement pattern or
two-pattern handover has been propagated.

The translation figure alone omits the cost of rotating a 10 km membrane.
The additional budget treats each tile as a uniform rigid square with inertia
per unit mass `(side²/12, side²/12, side²/6)`. Its Euler torque includes both
angular acceleration and the gyroscopic term. An explicit ideal actuator uses
three opposite edge-midpoint force couples, each with a 10 km lever arm and
zero net translation force. For torque components τ, the total force charged
to those couples is `2 Σ|τ| / side`; both spin-up and braking are counted.

Gravity-gradient and reflected-band torque receive a separate pointwise norm
envelope. It permits arbitrary illumination asymmetry, using the tile's
circumsphere and the point-mass force-gradient bound. Its sampled time integral
is a conservative allowance for that model, not an exact shadow-torque history
or a continuous certificate. Absorption, thermal stresses, flexible dynamics,
plume routing, power storage and actuator installation remain outside it.

The rigid-body couples are one possible implementation. Reaction wheels or
distributed sail torque need their own mass, momentum-storage and force
models; neither has been credited as free control. Added power hardware also
changes sail acceleration and therefore requires a new coupled propagation.

| Six-hour departure account | Mean impulse or equivalent impulse |
|---|---:|
| Applied translation, one completed two-hour burn | 1.40584 m/s |
| Nominal rigid-body electric torque couples | 2.72711 m/s |
| Sampled gravity/radiation disturbance allowance | 0.90790 m/s |
| Sum for this actuator option | **5.04085 m/s** |

The one-hour turn reaches **29.17 MN m** nominal torque and **5,833 N** total
couple force per tile. Maximum angular momentum is 3.41×10¹⁰ N m s per tile.
The force-couple account includes spin-up and braking; it is not just the
13.94 MJ peak rigid rotation energy. A derivative-convention check compares
the angular velocity to `Rᵀ dR/dt` before applying Euler's equation.

For the 1.805×10⁹ kg optical inventory and the stated exhaust/efficiency/cant,
the departure option uses **2.757×10¹⁴ J** and about **429,000 kg** of propellant.
Its peak combined input with the disturbance allowance is **68.02 GW**.
The provisional power hardware alone is **283 million kg**, or **15.7%** of
the optical mass, using 300 W/kg and the 25% peak allowance. That mass has
**not** been fed back into the sail trajectories. Storage, thrusters, the
torque structure and full vehicle hardware are additional open terms.

The same rigid-body calculation on the preceding six-hour service trace adds
a **1.66184 m/s** sampled disturbance allowance; its nominal Sun-tracking
rotation impulse is negligible and electric translation remains zero. Thus
this actuator/envelope option totals **6.70269 m/s equivalent impulse** and
3.666×10¹⁴ J for the executed service-plus-departure geometry. The disturbance
envelopes are allowances, not measured actuator expenditure. Acquisition,
return, repacking, handovers, recovery and repetition remain unpaid.

The inventory comparison remains conditional: at 3 m/s/day the prior
17.29-, 25- and 35-million-tile scenarios require 114.25, 165.24 and 231.34 TW,
with their 20% extra fixed mass, 10% thrust duty, power-hardware and seven-day
propellant feedback. Their break-even recurring impulses against 238.35 TW
are about 5.97, 4.24 and 3.09 m/s/day. A single departure supplies neither a
repetition rate nor a full-cycle impulse. No number here is a collection
capacity or a contribution to the 100 TW delivered-power ambition.

## Resource use, checks and next decision

The recorded numerical stages took **1,013 s**, about **16.9 minutes**, plus
a small diagnostic reconstruction. Peak RSS was **724 MiB** and new compressed
raw output **32.28 MiB**. Execution used one calculation process and one
BLAS/OpenMP thread. This remains within the four-vCPU environment and the
60-minute / 2 GiB / 1 GiB budget. The search used 1,820 static geometry samples,
96 distinct proposals, one of six permitted local corrections, one coupled
candidate and one tighter replay. No enlarged fleet search was run.

The first meaningful positive checkpoint is a **nominal departure path**,
with a separate and substantially larger attitude-actuation account. Its
small remaining clearance margin needs robust optimization and hardware-mass
feedback. A slower turn, a partial-sail coast and other attitude actuators
are design variables; full feathering is not a requirement for a solution.
The next bounded search can vary those choices while solving handover and
return from the new actual terminal states. It must retain the beam/coverage
constraints and the cost of repacking rather than charging only another
low-impulse endpoint fit.

The source-bound products are `departure_screen.json`, `departure_refine.json`,
`departure.json`, `departure_validation.json` and `departure_budget.json` in
`results/`. The coarser interval-guard failure remains visible in the original
coupled product. The newer validation product records its resolution without
erasing it. Earlier service, return and global-fleet products are unchanged.
With the pinned service trace and ephemeris present, reproduce using:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.departure_screen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.departure_refine
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.departure_run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.departure_validate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.departure_budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```

All **103 study/domain tests pass**. New analytic tests check the symmetric
thrust pulse's displacement, zero-net-impulse preparation modes, the vector
acceleration cap, minimum-norm control, actual spline rates, rigid-body
spin-up/braking costs and two-sided disturbance envelopes. Product tests bind
sources and inputs, preserve failed reductions and reject invented recurring
costs. The final repository-wide check is recorded in `checks.json`.
`make check` passes the layer check over 480 files and reports **648 Python
passes, 55 skips and the same 13 existing failures** in 34.04 s. The failures
concern missing spectral inputs, CM1 configuration, process lookup,
ring-comfort checks and their dependent protection-architecture calculation.
The later JavaScript, provenance and ensemble targets are not reached.
`git diff --check` passes.

Full return, two-pattern handover, continuous fleet coverage, recurring
operating cost, hardware-mass propagation and structural realizability were
not executed or established. This checkpoint is ready for the decision about
expanding the search.
