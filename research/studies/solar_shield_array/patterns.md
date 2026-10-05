# Naturally deforming local patterns

## Bounded continuation plan

This continuation starts from remote commit
`19f0eee92ca5b5ed61dd26b3943d40bee5274cbb` on
`research/solar-shield-habitat-array`. Main remains unchanged.

The objective is to find a finite local pattern that preserves a moving
20 km diameter receiver window inside the four-lunar-radius target during a
six-hour service passage, with natural differential motion and little electric
correction. Optical allocation stays at 50 g/m² and physical tiles stay 10 km
square. This experiment tests a component needed by the global assignment;
it will not scale a moving local window into a global fleet claim.

First, replace the pilot witness's prescribed roll with smooth transported
roll and replay its force and reflected-beam history. The numerical attitude
rate limit is 0.1°/s. Report angular acceleration separately; torque authority,
flexibility and actuator mass still require engineering.

Next, use variational dynamics along an actual-date reference passage to
choose initial relative velocities. Initial pitch, depth staggering, population
and velocity gradients are design variables. The pattern can change shape
throughout the passage. Optical constraints use finite clear squares, solar
limb rays and a moving window that stays inside the protected disk. Select
candidate parameters according to the worst coverage sample, then propagate
all individual members with mutual interception changing their sail forces.
Check the assumptions behind any neighbour reduction against actual geometry.

The budget is one calculation process and one BLAS/OpenMP thread, at most
24 affine geometry trials, three six-hour coupled propagation cases with at
most 441 members, and one independent tighter pattern replay. Allow twenty
minutes of numerical wall time, at most 2 GiB RSS and 1 GiB of raw arrays.
Repository checks have their own five-minute allowance. Stop at the first
useful local coverage/control checkpoint and preserve failed candidates.

A successful local passage requires sampled full-Sun coverage, checked
finite-square separation between output samples, positive reflected-beam
margins, bounded smooth attitude commands, consistent mutual shadow force,
and numerical replay within 50 m. Search additional dates and source directions
for coverage failures before accepting the sampled result. Continuous spatial
and temporal certification remains a separate evidence state.

Report the measured correction impulse and propulsion energy, initial velocity
preparation and end-state differences. Acquisition, global handovers, return
formation, recovery and their recurrence remain distinct tasks unless explicitly
executed. Collection capacity and 100 TW delivered power remain separate from
propulsion. The 238.35 TW annual held-screen comparison retains its documented
ideal-optics and mass-closure assumptions in `pilot.md`.

This plan precedes the numerical work. Results and exact executed checks will
be appended at the checkpoint.
