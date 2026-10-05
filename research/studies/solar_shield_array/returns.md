# Local handover and orbital return

## Bounded experiment

This continuation starts from remote
`4e8a7fd99039ffddcfeef6a52685af578edb1bb4` on
`research/solar-shield-habitat-array`. The remote head was checked before
starting. Main remains unchanged. The input is the source- and byte-verified
361-member terminal state from the tighter six-hour `patterns.json` replay.
The 50 g/m² optical allocation, physical 10 km squares, four-lunar-radius
target and existing climate spectrum remain fixed.

The objective is to determine whether that actual local pattern can hand off
its protected rays and reach another useful service passage at modest impulse.
Follow an evolving ephemeris reference through its next sunward passage;
do not prescribe a Sun-following circular orbit between services. Return date,
tile-to-lattice assignment, smooth maneuver coefficients and attitude roll are
decision variables. The cost is the sum of individual maneuver impulse,
converted separately to propulsion energy and propellant. Include both
departure and arrival preparation; retain unresolved recovery and repetition.

Start with a linearized shooting screen using actual terminal offsets and
velocities. Try at most 24 return-date/assignment/roll combinations. Rotating
the target assignment permits a tile to occupy another lattice location on
the next pass. Inspect finite-square encounters along the proposed arcs,
including near the relative-motion focusing times. Use offending pairs to
refine the trajectory/control proposal, up to five bounded correction steps.
Keep failed assignments and their witnesses. An infeasible restricted family
is not a general impossibility result.

Promising candidates proceed to individual nonlinear propagation with sail
force, mutual illumination and finite-Sun body eclipses. Use smooth burns,
an assumed translational acceleration cap of 0.001 m/s², a 0.1°/s attitude
rate cap, the existing reflected-beam exclusions and 100 m finite-square
clearance. A local handover uses a second independently propagated pattern
and tests the same receiver rays at simultaneous actual dates. Coverage,
encounters, endpoint defects and force-reduction assumptions must be checked
again on the propagated states. A tighter replay must agree within 50 m.

Allow **60 minutes of numerical wall time**, one calculation process with
one BLAS/OpenMP thread, 2 GiB address space and 1 GiB raw output. The budget
includes at most three nonlinear 361-tile return candidates, one 722-tile
short handover and one tighter replay of the strongest useful case. Check
resource/time limits between stages and save intermediate evidence. Repository
checks have a separate five-minute allowance. Stop at the first meaningful
feasibility or cost checkpoint, including an encounter or optical obstruction
that needs a wider control family.

Success means that a tested local handover/return preserves its assigned
finite-Sun rays, separation, safe beams and bounded attitude, and that the
measured all-member impulse leaves room beneath the conditional cost thresholds
in `patterns.md`. Sampled success and a continuous certificate are different
evidence states. A single return cannot establish long-term recurrence or a
complete fleet. The 238.35 TW held-screen benchmark retains its stated ideal
optics, mass closure and propulsion assumptions. Collection capacity and the
100 TW delivered-power ambition receive their own ledgers.

This plan is recorded before executing the search. Results and actual checks
will be appended at the checkpoint.
