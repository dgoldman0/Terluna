# Opening a safe departure path

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

The plan is checkpointed before numerical search. Its `make check` run passes
the layer check and reports 639 Python passes, 55 skips and the same 13
previously recorded failures (missing spectral inputs, CM1 configuration,
process lookup, ring comfort and the dependent protection-architecture
check). Later targets are not reached. Results and executed checks will be
added here; the numerical search and subsequent return/handover stages are
not yet completed.
