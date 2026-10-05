# Encounter-directed return investigation

## Starting checkpoint and measured obstruction

On 5 October 2026 both the remote research head and the clean local branch
are `7e18c29a00342a5bcec143fa3e500147b462b333`. Main is
`2f2e0a104a3a1ea962a6dead4179f522762f997e`. Continue on
`research/solar-shield-habitat-array`, preserving the reviewed products and
their source identities. Read the research entry points, domain instructions,
`packing.md`, `closure.md`, `collection.md`, and the implemented propagation,
clearance, control, mass, torque and collection accounts before this plan.

Inspection of the actual tighter return states confirms pair 325/347 at
hour 13.893535. Its body-frame displacement is
(-9314.920, 8885.427, 100.000) m and normal separation is closing at
0.361015 m/s, including rotation of the frame. At hour 12 the displacement
was (-10616.375, 11904.536, 2537.498) m. The first return burns and inherited
relative velocities bring the squares into overlap along both edges as their
normal spacing closes. Pair 338/341 has an earlier 280.75 m sample at hour 13.
Repair must therefore consider the surrounding population and subsequent
approach, rather than only the stopping pair at its final geometry.

The preserved comparator carries all 361 physical 10 km squares, 9.89 km
clear sides, 50 g/m2 optical mass and the four depths spanning 36 km.
The propagated hardware allocation is 2.26330 billion kg for the pattern,
including 361 million kg fixed additional mass and 97.298 million kg power
hardware. Its initial epoch through hour 12 remains the accepted local
service/departure comparator. Return and recurrence remain open.

## Bounded computation plan

Numerical work is capped at **3600 seconds of wall time**, **one calculation
process**, **one numerical thread**, **2 GiB address space** and **512 MiB new
raw output**. The initial raw-state inspection receives a conservative
10-second charge. Repository checks have a separate five-minute allowance
per required pre-commit gate. Numerical producers run sequentially. Record
measured usage and partial attempts; preserve usable outputs at each stage.

1. Spend at most 600 seconds on twelve inexpensive proposals and encounter
   diagnostics. Compare the unmodified burn with removing or delaying its
   first arc, a few smooth in-plane turns, and encounter-local corrections.
   Examine whole-population physical-square distances through at least hour
   18, including intervening dates. Trace how each change alters the endpoint
   obligations. Refit a return target only when its prerequisite escape passes.
2. Allow at most three coupled trials, selecting their controls from the
   measured bottleneck and cost. Keep all mutual-shadow momentum and the
   existing per-member hardware mass. Any optical or mass change needs new
   coupled propagation from the first changed state. Continue each trial past
   its correction and inspect the next obstruction; never accept frozen
   endpoint geometry as an executed maneuver.
3. Reserve one tighter replay of the decisive candidate using 16 rather than
   8 solar force sources, maximum steps of 60 rather than 120 seconds and
   relative tolerance 2e-12 rather than 2e-10. Start from the corresponding
   already verified service/departure states, or replay earlier stages if
   their command or mass changes. Compare actual states, encounter dates and
   witness pairs, including reconstruction error.
4. Reserve 600 seconds for executed per-tile torque, propulsion and optical
   accounting, grazing-source refinement and the energy/material comparison.
   New raw products go under `research/runs/solar_shield_array/return_repair/`;
   source-bound compact results and the tested report stay with this study.

Acceptance retains **100 m physical surface clearance**, with a conditional
between-sample guard and subtraction of twice the per-centre replay and
reconstruction discrepancy. Require replay agreement below 50 m, thrust at
most 0.001 m/s2, attitude rate at most 0.1 degree/s, positive finite-Sun beam
exclusions and no unmodelled body eclipses. Keep finite-Sun coverage throughout
assigned service. A new feasible escape prefix establishes only that prefix.
A return requires arrival position, velocity and attitude plus the next
service trajectory; a repeatable cycle also pays recovery and repeat-epoch
changes. Attempt coupled handover only after those prerequisites pass.
Pause at the first meaningful verified return or energy-trade checkpoint.

## Comparable energy and material accounts

Keep **23.835 TW**, ten percent of the 238.35 TW held-screen reference,
visible. Evaluate intermediate operating allowances before any replacement
target is adopted. None is adopted by this plan. Separate first-intercept
sunlight, redirected band, conversion and capture assumptions, translational
and attitude propulsion, other operating and thermal loads, and delivered
electricity. Track all 361 members after service release, including feathered
coast. Integrate only actual executed trajectories; a full-cycle integral and
average require an actual full cycle.

Reconstruct the 586 TW extra-dimming account from 1235 W/m2 times 5% across
the lunar disk, and compare it with the held aperture and moving clear film
using identical spectral allocations. The local twelve-hour interception
fraction of about 59% cannot establish fleet yield. Report duty-cycle and
inter-pattern-shadow thresholds as explicit sensitivities, without claiming
that multiplying by the 17.29-million capacity relaxation supplies a fleet.
External capture of already redirected beams is an accounting case until
collector geometry, momentum, mass, conversion, heat and routing are solved.
Extra electricity does not pay for exhausted material: show propellant,
installed mass and replacement-lifetime sensitivities separately.

Before the early checkpoint, `make check` passed the layer check and reported
686 Python passes, 55 skips and the same 13 baseline failures in 42.81 s.
The failures concern missing pinned spectral inputs, absent generated climate
configuration and the environment's process lookup. Later Makefile targets
were not reached. No numerical implementation was changed for this checkpoint.
