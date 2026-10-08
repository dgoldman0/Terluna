# Coupled constrained cycle optimization

Continue `d8de3e8438832bd51d916d32a9adef5c6d9d1dc4` on the existing research
branch. The author waived the preliminary remote-head check and authorized
an early checkpoint and results push. Main remains unchanged.

## Target and bounded experiment

The target is an executed repeatable cycle: assigned initial service,
departure, return, six-hour next service, and departure from its actual
endpoint. Short-horizon collision repairs initialize that problem. They must
report their terminal position/velocity debt to the next passage, rather than
treating a later collision as closure. The existing 1.535 PJ corridor remains
a geometric comparator with 145 actuator overloads.

Keep all 361 loaded physical 10 km squares, 9.89 km clear apertures, four
depth levels spanning 36 km, the unchanged spectral reflection/transmission
allocation, finite-Sun six-hour moving 20 km receiver within four lunar radii,
and the 100 m physical clearance floor plus numerical error. Preserve the
six-hour initial service. Retain the 36.186 TJ slower-turn and 51.859 TJ
previous twelve-hour comparators. No magnetic integration or new supply
hardware is introduced in this initial benchmark. Actual generation, storage,
delivery and unmet demand stay explicitly unknown until hardware is specified.

Allocate **1,800 CPU seconds**, **one numerical process/thread**, **2 GiB
address space** and **512 MiB new raw output**. Reserve at least 500 CPU seconds
for independent tightened propagation, clearance/optical/load accounting and
report products. Repository checks have the inherited separate five-minute
allowance per commit. Charge partial runs, compilation and numerical probes;
checkpoint coefficients, Jacobians, iterates and failed branches. Stop before
a larger optimization or training campaign. No machine learning is planned
until the physics evaluation benchmark establishes its value.

1. Reconstruct directly from the source-bound committed `joint_seed.json`,
   avoiding regeneration of the historical ignored raw chain. Reproduce the
   slower-turn service and hour-13.799825 limiting encounter. Measure coupled
   evaluation cost, deterministic repeatability, and directional sensitivity
   convergence before deciding the number of optimization variables/steps.
2. Use a small layer/row and selected-member translational basis, with smooth
   burns during departure. Test timing/departure freedoms only within the
   measured remaining budget. All effects and shadow forces include 361 tiles.
   Independent attitudes remain excluded from this first basis: the general
   oriented reference has force, torque and contact support but rejects body
   eclipses, and its cost needs a separate coupled benchmark.
3. Build local responses by perturbing controls in the coupled propagator.
   Use constrained trust-region steps and compare predicted geometry/state
   changes against nonlinear replay. Reject/shrink/rebuild inconsistent steps.
   Treat polygon/active-contact changes with finite perturbations, retained
   separating branches and exact-geometry replay, without smoothing the
   acceptance geometry. Separate any feasibility-restoration violations and
   slacks from accepted physical performance.
4. Impose per-member power, energy and finite-square separation constraints.
   Tie the short-horizon terminal conditions to the next service position and
   velocity and its subsequent departure. Clearly label any approximate tail
   transport as a terminal initialization model; it cannot certify service,
   recoverability or recurrence. Extend the coupled horizon when the local
   step and remaining computation permit.
5. Independently tighten Sun quadrature, step and tolerances for the selected
   candidate or obstruction. Refine close-approach intervals and deduct the
   measured reconstruction/replay error. Report actual executed energy and
   terminal state/attitude, interception chronology and per-member load, with
   off-service replacement required from hour 6. Complete next service and
   another departure remain separate acceptance gates.

Keep **238.35 TW** held-array comparison and **23.835 TW** target.
**17.29 million tiles** remains an unvalidated capacity relaxation. Longer
returns pay the conditional cadence/6-hour inventory multiplier, with spatial
coverage, assignment and handover still required. No failed-prefix fleet
extrapolation earns an operating-power claim.

## Implementation audit

`packing_control.response_maps` differentiates common centre gravity only.
Its response has no cross-tile blocks and no derivative of visibility,
reflected force or shadow moment. `measured_base` adds an observed replay
defect and transports its remaining tail with the gravity STM; it changes
the intercept, not the control Jacobian. Existing trust bounds in
`independent_impulses` bound coefficient changes but do not test a nonlinear
actual/predicted merit ratio.

`service_constraints.fit_free_controls` and
`service_arrival_controls.fit_service_arcs` already supply sparse HiGHS
inequalities, conservative individual acceleration/power caps, weighted L1
energy and exact partial-arc charging. They use prescribed smooth burn
windows and tile-specific vectors with the same gravity response map. Their
ray inequalities freeze a selected owner for sparse receiver/Sun rays;
they do not differentiate shadow-dependent dynamics. Free endpoints remove
the arrival lattice but do not remove that force-response limitation.

Reuse the tested DE440 environment, exact compiled parallel shadow unions,
joint eclipse geometry, finite-square distances, rigid torque accounting,
coverage polygons, smooth arcs and source-bound loaded states. The proposed
increment is a coupled control Jacobian plus an actual/predicted trust test,
not another rename of the existing defect-correction LP. The new low-dimensional
basis deliberately trades control freedom for affordable coupled derivatives;
failure constrains only its tested neighborhood and parameterization.

The early required `make check` passes the layer check and reports 745 Python
passes, 55 skips and the same 13 baseline failures in 32.13 s (missing spectral
and generated climate inputs, plus process lookup). Later Makefile targets
are not reached. The early checkpoint is `d55082e05d1622d82f1cb4121cf73a5f8f60854d`.

## Bounded result

**Coupled sensitivities substantially improve local prediction. The tested
control basis produces no accepted clearance repair, return, next service or
cycle.** A held-out 0.01 m/s perturbation is predicted within **0.04976 m** by
the coupled response, versus **3.65181 m** by the gravity-only response with
the same coupled baseline intercept. The comparison isolates the incremental
control response; it does not compare against the older wholly unshadowed
trajectory's kilometre-scale baseline error.

The final retained-contact trial fits installed actuator ratings and the
energy ceilings. Its tighter physical execution nevertheless reaches the
100 m floor at **hour 12.814874**, on **338/341**, earlier than the preserved
slower-turn encounter at 13.799825. Its first twelve hours cost **37.21538 TJ**,
about 2.84% above the preserved 36.18635 TJ passage. No improved physical
clearance/power/energy trade is established. The retained optimizer state
therefore remains the original zero-translation slower turn. Failed controls,
linear models and actual stopping states are preserved separately.

The primary function remains EUV/ionizing-radiation rejection. The spectral
products, apertures, seams and initial service assignment are unchanged.
Main magnetic protection remains outside this experiment.

## Measured cost and response quality

The committed epoch seed permits a fresh run without regenerating historical
ignored trajectories. The coarse original service costs 27.77 CPU s. Its
departure/restoration through hour 14.5 costs 47.61 CPU s initially and
40.05 CPU s on repeat. The two runs have **zero bitwise position or velocity
disagreement**. The first 100 m event is at 13.79961046 h, 0.773 s from the
previously recorded tighter baseline. Continuing this diagnostic past that
event creates optimization residuals and no physical service credit.

| Held-out response test, all 361 members | Coupled finite perturbation | Common gravity response with corrected baseline |
|---|---:|---:|
| Maximum position prediction error | 0.049761 m | 3.651807 m |
| RMS position prediction error | 0.010983 m | 0.444787 m |

The derivative uses a 0.02 m/s maximum modal perturbation; the independent
test uses 0.01 m/s. Halving the perturbation changes the position derivative
by at most 0.49761 m per dimensionless coefficient, while the full response
is 355.95 m. This supports a useful local response; it is not a derivative
convergence proof across every shadow transition. Subsequent replay tests
measure model error at the actual proposed step sizes.

The table covers the entire 6–14.5 h diagnostic window, including restoration
states after contact. Restricting the comparison to the first twelve hours
gives **0.02719 m versus 2.15843 m** maximum error and 0.00501 versus 0.23366 m
RMS error. The predictive improvement is present before any return contact.

The finite differences include the entire state of every tile in every force
evaluation, finite-area gravity, mutual-shadow force changes and the resulting
shadow/gravity torque account. A moved member can change another member's
force. No population reduction or independent-force replacement is used for
these responses. The unchanged compiled geometry takes exact rectangle unions
for each sampled Sun direction. The long diagnostic tail uses the available
joint body/array eclipse model with its area quadrature.

## Control space and whole-cycle target

Eight dimensionless parameters are tested. Six multiply inertial velocity
corrections along the common feathered normal at hour 13.8; a unit coefficient
has 0.1 m/s maximum modal amplitude. Each mode has zero mass-weighted mean.

| Freedom | What it permits | Restriction |
|---|---|---|
| Row parity, column parity, parity product | Three contrasts spanning the four depth classes | A single fixed thrust direction |
| Linear row gradient | A broad relative-velocity shear across all 19 rows | No independent column gradient or higher spatial modes |
| Opposite corrections on 325/347 and 338/341 | Two selected individual pair adjustments | Neighbouring members receive only the grouped modes |
| Common turn duration | Smooth 90-degree turn over 5.5–6.5 h, starting at hour 6 | Shared attitude; fixed 90-degree final tilt |
| Burn endpoint | A normalized smooth pulse from hour 6 to hours 11.5–12.5 | One pulse per departure; its sensitivity is zero at zero impulse |

The six-hour initial optical service is preserved. Initial states, hardware,
mass and assigned receiver are fixed. Independent sail normals, free tile
rolls, arbitrary per-member histories, in-plane thrust directions and arrival
arcs are absent from this first basis. The general oriented-square reference
supports shadow/force/torque/contact calculations but rejects body eclipses;
it is not substituted into a coupled return without further work.

The cycle schedule returns to Sun-facing attitude over eight hours before
the inherited **46.138697 h** next-service date, holds that service attitude
for six hours, and repeats the departure pulse and turn. The terminal dates
are next-service start, next-service end and hour **59.138697**, after the
subsequent departure. Targets transport the original service geometry and
relative velocity, and the actually changed first departure state, into the
corresponding later frames. Mean orbital motion is retained explicitly.

One complete **coupled restoration trajectory** supplies the baseline tail.
Changed early terminal states are initially transported along it with a
gravity state-transition map, including the repeated later pulse. The
normalized RMS of the three cycle residuals may not increase. Position and
velocity scales are 100 km and 5 m/s; their maximum physical debts are also
reported. Baseline maxima are approximately **319.42 km / 11.92 m/s** at
next-service entry and **344.51 km / 12.57 m/s** after the following departure.
These are substantial unresolved boundary defects. The tail transport is an
explicitly approximate initializer, not a coupled cycle Jacobian or a
recoverability certificate. Arrival timing, assignment and service packing
remain fixed in this bounded run; failure does not exclude alternatives.

## Trust steps, contact topology and retained failures

The constrained SQP subproblem combines coupled finite-difference responses,
exact thrust-vector norms, linearized actual attitude loads, per-member peak
power bounds, a 51.859 TJ first-twelve-hour ceiling and a 100 TJ estimate
through next service. The latter carries the baseline future attitude demand
but leaves uncomputed acquisition controls unpaid. It cannot certify a cycle
energy budget. There is no prolonged holding segment. The physical hardware
allocation is unchanged.

Clearance restoration has an explicit scalar slack, initially targeting a
400 m separating margin. It later targets 250 m in a shortened hour-14
initializer. Slack is reported in metres and cannot relax power or energy.
No virtual acceleration is inserted into the dynamics. Positive restoration
slack prevents acceptance. The solver's objective is compared with actual
coupled replay; poor agreement shrinks the trust region. Measured replay
secants update the coupled state/load response. Actual earlier encounters
reject a step even when its retained linear contact constraints improve.

No force/contact smoothing is used. Finite-difference stencils span changes
in shadow topology, while each subproblem freezes selected separating faces.
New contacts are detected in replay. A final rebuild unions the contact
branches exposed by all four rejected short-horizon trials, increasing retained
constraints from **19,239 to 21,331**. Source quadrature remains an explicit
numerical approximation, independently refined in the physical replay.

| Method stage | Replayed model agreement | Physical outcome |
|---|---|---|
| Initial hour-14.5 step | 17.73 m maximum state prediction error; merit ratio 0.991 | New earlier contacts; near-duplicate attitude knots also produce false power spikes; rejected |
| Second initial step | 30.59 m; merit ratio 0.989 | Fits sampled ratings; newly exposed contacts worsen the comparison; rejected |
| Four stable-timing hour-14 trials | 4.82, 0.992, 0.547 and 2.06 m errors; ratios 0.993–0.999 | All fit power/energy; first floor events near 12.815 h; all rejected |
| Retained-contact rebuild | 2.910 m; ratio 0.99187 | 286.149 m restoration slack; no feasible point; earlier encounter verified independently |

This is useful accuracy evidence and a failed local feasibility search.
Accurate derivatives do not make this eight-parameter contact problem feasible.
The remaining slack applies to the tested trust region and separating
branches, not to all trajectories or all controls. The actuator and energy
bounds have headroom in the final trial. Additional spatial/directional control
and cycle acquisition freedoms deserve testing before relaxing those bounds.

The failed timing implementation inserted turn knots within about 1e-10 s
of the regular attitude grid for coefficients extremely close to a bound.
Spline differentiation then created artificial torque spikes. The stable
implementation coalesces timestamps within **1 ms**, retaining the original
smooth angle law. Frame/rate/acceleration histories of the baseline and finite
derivative probes agree within 1e-12, allowing those expensive derivatives to
be reused explicitly. The failed trial remains preserved. A regression test
checks bounded angular acceleration at a nearly grid-aligned turn endpoint.

The first attempt also exposed a contact-merit issue: the minimum over a
newly selected, distance-limited contact set is not directly comparable to
the minimum over the preceding set. The continuation retains actual event
progression as a rejection gate, and the final rebuild retains the new rows.
This prevents a better score on one contact set being promoted into safety.

## Independent physical execution and complete prefix account

The decisive **rejected contact-rebuild trial** is independently propagated
from the original epoch with 16 instead of 8 Sun sources, 60 instead of 120 s
maximum steps and relative tolerance 2e-12 instead of 2e-10. Its unchanged
mass is **2,263,298,487 kg**. The physical run stops at the clearance floor.
The earlier restoration trace after that stop receives no execution credit.

| Gate or executed quantity | Result |
|---|---:|
| First service | 100% sampled full-receiver coverage; 49 dates × 25 sources; original seams and clear apertures |
| First twelve hours | 37.215376 TJ; 1,715.22 m sampled surface separation |
| Twelve-hour conditional clearance after numerical allowance | 365.278 m |
| Actual stopping hour | **12.814874263** |
| Limiting pair | **338 / 341**, zero-based |
| Terminal body-frame displacement | (-10,099.754, 5,429.061, 7.011) m |
| Coarse/fine event difference | 1.86036 s |
| Maximum position difference including reconstruction | 7.85889 m |
| Electricity through actual stop | **37.217094 TJ** |
| Attitude / translation energy | 36.174855 / 1.042239 TJ |
| Maximum individual power | 8.078903 MW |
| Minimum installed/required peak ratio | 6.19344; zero overloaded members |
| First-intercept bolometric-equivalent energy | 1.248760e18 J |
| Redirected spectral-band energy, already part of that interception | 1.722665e17 J |
| Ideal exhaust / modeled propulsion waste heat | 57,893.3 kg / 11.165128 TJ |

Adaptive all-pair square checks refine difficult intervals down to 0.5 s and
subtract twice the measured per-centre replay/reconstruction disagreement.
The interval ending 120 s before the stop retains a 138.835 m conditional
lower bound. The stop itself fails the numerical allowance, as expected
when the nominal distance is exactly 100 m. Bounds retain the inherited
0.15 m/s² centre-acceleration and 0.1 degree/s attitude-rate envelopes.
The tighter trajectory is a separate numerical integration of the same
physical model, not an independent physical theory.

At the final encounter date, pair 338/341 has 221.53 m normal separation in
the unmodified baseline and only 7.01 m in the changed trajectory; its in-plane
edge gap is about 100 m. This exposes the trade between normal compression
repair elsewhere and this neighbour's existing clearance. The basis has no
independent in-plane displacement mode to alter that edge encounter.
Coarse/fine executed control-energy and optical-interception integrals differ
by 0.00720% and 0.00515%, including their 1.86 s difference in stop time.

Actual terminal position, velocity, quaternion and angular velocity for all
361 members are in `results/coupled_validation.json`. Per-member power and
first-interception histories retain their complete chronology. Actual
generation, bus delivery, storage state, delivery loss and unmet electrical
demand remain **null** because no placed supply system has been specified.
Optical interception is not converted into free maneuver electricity.
Propellant depletion, flexible structure, thermal layout and detailed actuator
placement remain inherited model limits. No new hardware is priced as installed.

The first six-hour service remains the only credited passage. Replacement
coverage is needed from **hour 6**. There is no executed next service,
handover or recurrence. A 46.138697 h cadence would require at least
**7.68978 groups per six-hour service credit**, or eight undivided groups,
before spatial coverage and handover penalties. That is a conditional inventory
charge on an unachieved cycle. **238.35 TW**, **23.835 TW** and the unvalidated
**17.29-million-tile** capacity relaxation retain their existing meanings.

## Full-cycle diagnostic check

The final rejected controls are additionally propagated through the entire
restoration tail with coupled forces, including repeated departure controls.
This is explicitly post-contact diagnostic motion. Its 213.918 CPU s fit the
remaining allocation and test the cycle-tail approximation directly.

| Boundary | Maximum actual position debt | Maximum velocity debt | Terminal-proxy position error |
|---|---:|---:|---:|
| Next-service entry | 319.353 km | 11.9173 m/s | 88.757 m |
| Next-service endpoint | 343.309 km | 12.5599 m/s | 340.451 m |
| Following departure endpoint | 344.526 km | 12.5676 m/s | 1,100.231 m |

The proxy's small apparent residual improvement reverses under this coupled
check: normalized RMS is 1.1211418, compared with the baseline 1.1211312 and
the predicted 1.1211116. The terminal proxy therefore cannot accept a cycle,
even if the local clearance problem were solved. The actual coupled tail must
enter subsequent cycle-level sensitivities and constraints.

The tested centred receiver has only **6.996%** minimum finite-Sun coverage
over thirteen next-service dates and 25 sources per date. This excludes the
tested service assignment, not all translated receivers or changed packings.
The 101.560 TJ restoration-path control integral through the following
departure includes motion after failed contacts and has no operating-cycle
credit. The actual executed integral remains 37.217 TJ. Full return,
completed next service, recurrence, supply and handover remain unaccepted.

## Continuation and reproduction

`results/coupled_restart.json` stores exact float64 coefficients, basis,
starting and terminal arrays, trust radius, retained contact constraints,
contact and terminal Jacobians, and endpoint coupled responses in a checked,
compressed payload. `coupled_restart.restore()` reads it without the ignored
raw trajectories. A terminal state past contact is labelled restoration-only;
restart physical propagation from the saved actual hour-6 state and apply the
saved controls. Extending the horizon requires refreshed coupled responses.

The new method has eight parameters, about forty-second local coupled
evaluations, and demonstrated small-perturbation repeatability. This does not
justify a large training campaign. The next bounded experiment should add
early **in-plane** relative-motion control and arrival acquisition arcs,
retain exposed contacts from the outset, and differentiate the changed tail
with the coupled model. A small number of early directions followed through
the whole cycle, combined with late directions evaluated only after their
onset, is the smallest justified full-cycle sensitivity benchmark. Measure
its cost before increasing dimension. Keep the current exact service and
spectral requirements, per-member ratings, energy and inventory accounting.

Reproduction starts with the pinned DE440 kernel and committed loaded seed:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m protection.dynamics.ephemeris --download
python -m research.studies.solar_shield_array.coupled_benchmark
python -m research.studies.solar_shield_array.coupled_optimize
python -m research.studies.solar_shield_array.coupled_refine
python -m research.studies.solar_shield_array.coupled_contacts
python -m research.studies.solar_shield_array.coupled_finalize --fresh
python -m research.studies.solar_shield_array.coupled_restart --inspect
```

The attempt products retain the historical interruption and four-trial
continuation; a fresh run needs to preserve the preceding output as
`coupled_optimization_attempt.json` before `coupled_refine`, and the refinement
output as `coupled_refinement_attempt.json` before `coupled_contacts`.
The original interrupted attempt is not a required full six-iteration rerun.
The derivative artifacts can be regenerated individually using the stored
coefficients and `Experiment.evaluate`; running more iterations changes the
budget. Use the recorded controls for physical replay when reproducing the
reported case rather than assuming identical wall-time cutoffs across machines.

The original validation completed both IVPs and load accounts, then failed
JSON export on a NumPy boolean. `coupled_finalize` verified source/input/raw
hashes and recovered those completed trajectories, carrying their recorded
Hermite reconstruction errors into the geometric audit. A subsequent console
finalization error was corrected and its measured cost charged. `--fresh`
performs new integrations with the corrected typed writer. Neither recovery
changes controls, force histories or integrated loads. Scalar export and
near-grid timing have regression checks; partial products and their charges
remain recorded.

## Resources and checks

The bounded investigation charges **1,712.561 of 1,800 CPU seconds**:
1,517.852 s recorded by completed and partial producers, 90 s for the
interrupted optimization, 60 s for the failed validation export/audit,
24.708 s for the failed console finalization, and 20 s for ancillary
compilation and inspection. Peak resident memory is **472.54 MiB** under
the 2 GiB address-space cap; new raw output totals **290.11 MiB** under
512 MiB. Numerical evaluations use one process and one numerical thread.
The independent validation, recovery, whole-cycle diagnostic and final
audit are included. No larger optimization or learning campaign is started.

`coupled_checks.json` binds the nine result products and confirms unchanged
constants, scenario, spectral stack and UV transmission bytes. The portable
restart retains the unsuccessful neighborhood for continuation. This budget
establishes the value of coupled sensitivities, but does not establish a
feasible return or a complete repeatable cycle.

The final required `make check` passes the layer check on 615 files and
reports **754 passed, 55 skipped and 13 failed** Python tests in 30.46 s.
The failed test identities exactly match the early checkpoint: unavailable
spectral inputs, generated climate configuration and process lookup. The
failure stops the later JavaScript, provenance and ensemble Makefile targets;
those targets are not claimed as run. The two new test modules pass all
**nine** regressions independently in 0.77 s. They cover coupled/contact
responses, trust rejection, near-grid attitude timing, scalar export,
physical evidence gates, source bindings and portable restart integrity.
