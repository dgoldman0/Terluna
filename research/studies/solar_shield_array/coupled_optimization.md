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
are not reached. Results and reproduction details follow the bounded run.
