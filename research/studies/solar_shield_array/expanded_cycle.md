# Expanded coupled cycle investigation

Continue `762ce4d810e07ddc0788f315212c5a6f6bdca5fe` on
`research/solar-shield-habitat-array`. The author waived the preliminary
remote-head/concurrent-work check and authorized an early and final push.
Main remains unchanged. Primary magnetic integration remains set aside.

## Bounded decision experiment

Use the source-checked `coupled_restart.restore()` and `joint_seed.json`,
all 361 loaded 10 km physical squares, 9.89 km clear apertures and unchanged
spectral allocation. Keep the six-hour first service, 100 m surface clearance
plus measured numerical allowance, individual installed actuator allocations,
0.001 m/s² translation and 0.1 degree/s attitude limits. Preserve the 51.859 TJ
first-twelve-hour ceiling and 100 TJ through-next-service screening ceiling.
The selected zero-translation slower-turn comparator remains 36.186 TJ.

Allocate 1,800 numerical CPU seconds, one process and numerical thread,
2 GiB address space and 512 MiB new raw output. Reserve at least 400 CPU
seconds for independent replay, accounting and compact export; required
repository checks have the inherited separate five-minute allowance per
commit. Charge compilation, incomplete evaluations and diagnostics. Measure
one complete coupled evaluation and directional sensitivity first; use those
costs to limit differentiated directions and iterations explicitly.

Add independent in-plane separation/shear and arrival position/velocity
acquisition directions to the available early normal controls. Use smooth
finite burns and common finite-rate attitude changes. Arrival directions may
restart at their exact onset on the common coupled baseline, because their
prefix controls are identical. Every changed suffix still propagates all
361 interacting tiles through next service and the subsequent departure.
The full-cycle Jacobian must use coupled evaluations, never a gravity-only
transport of a changed prefix. Record any unmeasured columns separately.

Inspect singular values, effective rank, residual projection and finite-step
prediction error. Use constrained trust-region steps with energy in the
objective, explicit restoration slack and exact coupled trial replay.
Positive slack, missed optical service, resource violations, or contact
preclude acceptance. Continue past contact only for labelled optimization
diagnostics. Retain exposed contact choices and rejected trials. Diagnose a
small set of one-limit relaxations only after the actual-limits subproblem;
relaxed solutions never earn physical acceptance. Hardware changes require
new mass, photon forces and epoch replay; actuator ratings never imply supply.

## Repeatability gate

A safe service/departure/return/next-service/subsequent-departure sequence is
only a finite execution. Test transported relative position, velocity and
attitude at service entry, service exit and corresponding departure phase,
along with the receiver's actual-date position inside four lunar radii.
Use explicit finite-sequence closure tolerances of 50 m, 0.01 m/s and 1e-4 rad
as *search acceptance gates*, rather than established engineering tolerances.
The actual-date ephemeris is nonautonomous: small relative-state closure at
one return alone cannot establish indefinite recurrence. Recurrence evidence
requires further consecutive actual-date cycles from the computed endpoint,
without resetting state, preserving every service, safety and resource gate
and bounded defects. Report separately whether that extension was reached.

No generator, collector or storage system is placed in this experiment.
Generation, storage and unmet demand remain unknown. Track executed control
energy, individual loads and first-intercept/redirected light histories.
Preserve 238.35 TW held-array comparison, 23.835 TW operating target and the
unvalidated 17.29-million-tile estimate; do not extrapolate failed prefixes
into achieved fleet operation.

## Portable continuation

Preserve the expanded basis, coefficients, local responses, optimizer radius,
contacts, scenario and ephemeris identities, deterministic seeds, compact
ephemeris samples, dependencies and exact restart/replay commands. Provide a
CPU-budget-configurable runner. Check loading and representative replay with
the original ignored trajectory files absent. Prioritize larger local work
from the measured rank, model errors and active constraints. Stop after this
bounded investigation; do not launch the local campaign.

The early `make check` passes the layer check (615 files, zero violations)
and stops because this fresh environment lacks `pytest`. No Python test
result is claimed for that invocation. Restore the test dependency before
the final required check.
