# Energy-first continuation

The author prioritized reducing energy on 5 October 2026, with additional
clearance to be investigated after useful savings appear. Starting remote
head is `e65cc9c34a120926d3f5a7c6e590260e78275793` on
`research/solar-shield-habitat-array`; main remains unchanged.

## Objective and bounded first step

The operating target is 23.835 TW, ten percent of the 238.35 TW held-screen
comparison. Retain intermediate inventory/cost trade-offs. A recurring result
must include service, handover, return, repacking, recovery and maneuver
frequency, with collector and propulsion hardware mass in the dynamics.

Begin with the current departure's large attitude expenditure. Compare a
longer Sun-facing coast and smaller, slower turns, optimizing relative-motion
preparation for each. The first decision is whether this local schedule has a
substantially cheaper feasible region. Start from all 361 actual states after
the verified six-hour service, preserving that service and the previous
departure as comparisons. A local cost improvement can justify a subsequent
full-cycle search; it cannot establish recurring fleet power.

Decision variables are the common turn angle, axis and sign, its delay and
duration, and the three smooth depth/row/column impulse coefficients. A zero-
angle schedule is included. The reference remains freely evolving. Keep the
50 g/m² base optical allocation, 10 km squares and four-radius target.
Do not require a 90-degree turn or a completely feathered coast.

Use inexpensive gravity-variational proposals to minimize mean translational
impulse plus the explicit electric-couple attitude impulse. Seek finite-square
encounters and add local separating constraints. Rank proposal costs, then
verify promising candidates with coupled finite-Sun sail force and dynamic
mutual shadows. The proposal model never accepts placement on its own.

The actual clearance constraint stays at 100 m. Proposal padding is solely a
numerical allowance, and clearance is not an optimization reward. Keep the
0.001 m/s² translation cap, 0.1 degree/s attitude-rate cap, reflected-beam
exclusion and existing finite-Sun coverage checks. Released receiver windows
require a separate handover. The present coupled model rejects body eclipses;
an eclipse-containing return needs joint body/array illumination.

Compute rigid-square gravity-gradient torque by finite-area quadrature and
reflected-band torque from the first spatial moments of the same shadow unions
used for force. Compare the resulting electric couples with the earlier
conservative disturbance envelopes. Keep nominal torque, calculated model
disturbance and residual engineering uncertainty distinct. This is an ideal
rigid actuator account; hardware, flexibility and additional optical physics
retain their open status. Added power-system mass receives an explicit budget,
and an unpropagated mass increment cannot establish a hardware-closed result.

Track optical power and energy for every instantiated member through all
tested stages, using `collection.py`. No electrical generation or fleet
collection credit is inferred from a local patch. The held comparison uses
50 g/m², 30 km/s exhaust, 70% efficiency, 45-degree thrust cant, 300 W/kg power
hardware, 25% peak allowance and seven days of propellant; storage and complete
vehicle engineering are excluded.

## Resource budget and success criteria

Allow 45 minutes of numerical wall time, one calculation process and one
BLAS/OpenMP/solver thread, 2 GiB address space and 512 MiB new raw output.
Allow up to 48 distinct schedule proposals, at most four local corrections
per proposal, two coupled candidates of at most six hours and one independent
tighter replay. Torque and collection accounting may use up to 32 solar
directions, with temporal/source refinement on the selected trace. Repository
checks have a separate five-minute allowance. Preserve rejected proposals and
their limiting constraints; checkpoint between stages.

A useful positive result is a coupled local candidate that meets the existing
geometry, coverage-assignment, beam and rate gates and materially lowers the
same-model service-plus-departure energy. Require an independent replay within
50 m, refine between-sample encounters and expose the remaining continuous-
certification gap. A factor-of-two reduction in this partial account is a
screening goal. Compare it with the conditional full-cycle impulse allowance
at 17.29, 25 and 35 million tiles, leaving the unpaid return and handover visible.

Stop at the first meaningful energy or feasibility checkpoint. Failure in the
restricted starting pattern, controls or attitude schedules is local evidence.
The larger search still permits different initial states, orbital families,
phases, inventories, assignments and actuator architectures.

## Execution record

The plan precedes numerical search. Its `make check` run reports 655 Python
passes, 55 skips and the same 13 existing failures in 37.12 s; the layer check
passes and later targets are not reached. Numerical outcomes will be added
after the bounded experiment.
