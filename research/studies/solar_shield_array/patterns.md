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

This plan was pushed before numerical work as
`4b222b227d95cd93aa357676b23091388ef39dab`. The remote head matched `19f0eee`;
there was no newer work to reconcile. Before this documentation checkpoint,
`make check` reported 624 passes, 55 skips and the same 13 existing failures.

## The local result

**A naturally deforming 361-tile pattern preserves the sampled 20 km diameter
window for six hours with zero electric translational thrust.** It includes
mutual shadow force and finite-tile gravity. The window moves across the lunar
receiver plane and remains inside the four-radius target. This is evidence for
a low-service-impulse local component; fleet-wide coverage and recurring cost
remain unresolved.

The calculation starts at the common DE440s date **2026-10-04 TDB**. No orbit
is copied to another date. The reference tile starts at 15,000 km lunar
distance, with a retrograde osculating circular velocity and phase −π/8 from
the Sun direction. Sail force contributes throughout the propagation. The
normal follows the actual Sun direction during this service passage; orbital
positions and velocities evolve freely.

The inexpensive size check is demanding even for this small receiver window.
At roughly 15,000 km sunward distance the solar limb displaces a shadow by
about 70 km. A 10 km radius receiver therefore requires a footprint about
160 km across at that distance. A single tile cannot serve this window over
the whole solar disk. The selected 19×19 pattern spans 162.89 km in clear
width initially, with 8.5 km pitch and four depth levels spaced by 4 km.
Each physical tile is 10 km square; its 9.89 km clear width reserves a 10 m
seam and 50 m placement allowance at each edge. The optical mass is
**1.805×10⁹ kg**, including every guard tile in the pattern.

## Variables and executed search

The reference integration also propagates the gravity state-transition
matrix, partitioned into position-from-position and position-from-velocity
blocks. For initial local offset \(d\), the proposal chooses

\[
 \delta q_0=R_0d,\qquad
 \delta v_0=s\,\Phi_{rv}(T)^{-1}
       [R_T-\Phi_{rr}(T)R_0]d.
\]

Here \(R\) is the actual-date Sun frame and \(T=6\) hours. This is an endpoint
shooting proposal, with gravity linearized around a sail-assisted reference.
The proposal omits the sail-force gradient and mutual shadows. Both the
full nonlinear dynamics and unequal illumination are included in the later
individual-member propagation; the affine path is never used as its force
target.

The 24 proposals combine 19 or 21 rows, pitches of 8.5, 8.75 or 9 km, depth
steps of 2 or 4 km, and gradient scales \(s=0.95\) or 1. Each has 13 dates
over six hours and central plus eight limb sources. All pass those optical
samples. The selector first requires coverage and finite-square separation,
then prefers fewer members and greater separation. It selects 19 rows,
8.5 km pitch, 4 km depth and \(s=0.95\). No globally optimal pitch, population,
orbit or velocity gradient is claimed.

Only the first selected candidate needs coupled propagation. The service
impulse objective reaches its zero lower bound for this initialized local
passage, so further electric-control optimization of this passage is not run.
Return and handover optimization would have different endpoint and coverage
constraints.

## Force and geometry consistency

Each member has its own state. A 3×3 Gauss surface quadrature supplies its
gravitational acceleration. The reflected spectral fraction is the existing
ideal filter-only model at 50 g/m²; no photon-force cancellation is applied.
For each solar-force source, upstream finite clear squares project to
rectangles on each receiving tile. Their exact union is removed once from
the illuminated area. The source-dependent incidence cosine and inverse-square
flux then give the normal photon force. Multiple reflections, absorption,
membrane flexibility and material spectra remain outside this ideal model.

The sparse eight-neighbour calculation is checked against all non-neighbour
pairs. The whole solar cone, finite source distance and a 50 m allowance widen
the excluded shadow cylinders. Relative motion, curvature and frame rotation
widen each 60 s interval further. The selected case retains over **6.4 km**
of non-stencil exclusion margin. A separate unit test compares reduced
illumination with all-tile shadow unions. Partial Earth/Moon eclipses reject
this local model; every propagated force evaluation remains fully uneclipsed.

The mean illuminated fraction is **0.82334**, and the least illuminated tile
sample is **0.47887**. Unequal photon force moves some tiles **2.46 km** from
the affine proposal. These are physical changes in the integrated trajectories,
not errors repaired by an uncounted controller.

Coverage uses unions of perspective-projected finite squares, including
first-order tile and lunar motion during light flight. Every source checks
the entire polygonal receiver disk, rather than a sparse set of receiver
points. The force calculation uses instantaneous inter-tile geometry;
moving-interceptor errors across the local pattern are small relative to its
50 m allowance and remain a model approximation.

The main optical checks use 25 dates 15 minutes apart and nine source points.
The adversarial pass adds 35 dates, including half-cadence dates and 60 s
samples near the worst sampled date, with 16 rotated limb directions, the
centre and eight interior directions. No missed receiver area is found
beyond floating-point roundoff. These are sampled source/time results;
continuous optical certification has not been executed.

The physical 10 km squares retain at least **2.685 km** separating-axis
clearance under the swept interval guard, against a required 100 m. Its
acceleration and attitude-rate bounds are 0.06 m/s² and 3×10⁻⁷ rad/s;
the sampled maxima are 0.0221 m/s² and 1.99×10⁻⁷ rad/s. This conditional
between-sample guard has no flagged pair. A rigorous enclosure of the assumed
derivative bounds remains distinct from these numerical checks. Sampled
reflected-beam margin exceeds **96°** to the protected Moon and atmospheric
Earth exclusion, including the finite Sun, tile corners and pointing allowance.

## Attitude witness

Pilot member 1546 failed the previous prescribed-roll rate cap. Minimum
normal-to-normal transport followed by a smooth rotation spline lowers its
peak rate from **0.12096 to 0.06467°/s**, below 0.1°/s. The maximum sampled
angular acceleration is **0.0001203°/s²**. A new finite-tile force integration
over its 1,800 s witness changes position by less than 2 μm and velocity by
about 1.1×10⁻⁹ m/s. The one-second beam checks retain **4.303°** margin; rate
and acceleration are checked at quarter-second cadence.

This resolves the particular numerical rate witness for a rigid two-sided
square. It does not demonstrate torque authority, power, stress, flexible
stability or reaction-momentum disposal. Roll changes projected square shapes;
the earlier fleet capacity coefficients must be recomputed before incorporating
this attitude into that allocation. Other members and their safe transitions
still require their own histories. The local Sun-facing pattern has a much
lower sampled attitude rate, about **0.0000114°/s**.

## What remains unpaid

The freely evolving initial velocity field differs from a lattice translating
with the reference centre by a mean **0.940 m/s**, at most **1.659 m/s**.
This is an insertion reference, not a performed acquisition maneuver. At the
end, relative to a Sun-frame lattice, the maximum position mismatch is
**2.976 km** and mean velocity mismatch is **1.028 m/s**, at most **1.801 m/s**.
Zero service thrust therefore does not establish a repeatable zero-cost cycle.

| Account | Executed status |
|---|---|
| One six-hour service passage | Zero electric translation impulse and propulsion energy |
| Initial acquisition | Velocity difference measured; acquisition trajectory and hardware uncomputed |
| Handovers to complementary service patterns | Unexecuted |
| Coast/return, next service preparation and recovery | Unexecuted |
| Recurring maneuver frequency and total fleet impulse | Unmeasured |
| Attitude actuation, flexible structure and full vehicle mass | Unclosed |
| Collection capacity and 100 TW delivered-power ambition | Unevaluated |

The 361 members cannot be scaled directly by local-window area into a fleet
inventory. Much of the patch is a finite-Sun guard border. Reusing borders
between adjacent windows would change interaction and placement constraints;
those interfaces have not been solved. The previous 17.29-million-tile
capacity relaxation remains conditional, with individual placement unresolved.

`patterns_budget.py` keeps these ledgers separate. On optical mass alone, the
insertion reference would cost **5.14×10¹³ J** and **79,948 kg** of propellant
at the reference exhaust velocity, efficiency and cant. This is a conditional
conversion of an initial velocity difference; the maneuver was not simulated.
Adding its mean velocity magnitude to the terminal mismatch gives a
**1.967 m/s velocity-only proxy**, which leaves the position mismatch unpaid.
It is neither a bound nor a measurement of recurring correction.

The existing inventory/impulse trade still determines the useful search region:

| Conditional inventory | 1 m/s/day | 3 m/s/day | 5 m/s/day |
|---|---:|---:|---:|
| 17.29 million tiles, earlier capacity relaxation | 36.93 TW | 114.25 TW | 196.57 TW |
| 25 million tiles, sensitivity only | 53.41 TW | 165.24 TW | 284.30 TW |
| 35 million tiles, sensitivity only | 74.77 TW | 231.34 TW | 398.03 TW |

These are mean **propulsion demands**, including the earlier scenario's 20%
additional fixed mass, 10% thrust duty, power-hardware mass feedback and
seven-day propellant buffer. Added mass changes sail response and requires
new coupled propagation. The comparison is **238.35 TW** for the annual
mass-closed ideal-optical held screen at 50 g/m², 30 km/s exhaust, 70%
efficiency, 45° cant, 300 W/kg power hardware and 25% peak margin. Neither
calculation includes complete vehicle engineering or storage. Installed
propulsion ratings in the data product are hardware ratings, not collector
capacity. No electrical collection or delivered-power credit is taken.

## Replay, checks and checkpoint

The coupled candidate uses DOP853 with at most 60 s steps, relative tolerance
2×10⁻¹⁰ and eight solar-force directions. The independent propagation starts
from the same individual initial states, halves the maximum step, tightens
the tolerance to 2×10⁻¹² and doubles the force directions to sixteen. Maximum
position disagreement is **0.661 m**, and velocity disagreement is
**6.48×10⁻⁵ m/s**. The 50 m replay gate passes. Both integrations retain the
same sampled coverage and separation results. This is a joint integrator and
solar-force quadrature sensitivity check; it does not isolate their separate
errors or independently validate the ephemeris/force model.

The full numerical run took **423.3 s**, peaked at **217.3 MiB RSS**, and wrote
**12.77 MiB** of compressed raw arrays. It used one process and one BLAS/OpenMP
thread. Of the predeclared budget, 24 affine trials, one coupled candidate and
one tighter replay were used. The remaining two coupled cases were not run
after reaching this checkpoint. The attitude replay was included in that
wall-time total. The cost ledger is an algebraic reduction of the saved product.

The first check of the new shadow/roll/geometry functions passed four tests:
independent polygon Boolean unions, all-tile versus reduced illumination,
transported-plane preservation, and a deliberately uncovered solar-limb ray.
All **85 study/domain tests** pass, including saved-product identity and
unperformed-return accounting. The repository check reports **630 Python
passes, 55 skips and the same 13 existing failures** (missing spectrum inputs,
CM1 configuration, process lookup and ring-comfort checks). The staged-source
layer check passes over 461 files. Make stops before the JavaScript,
historical-provenance and ensemble targets; they were not rerun at this
checkpoint. [checks.json](checks.json) retains the exact scopes.

Source and shared-constant hashes, raw-array hashes, individual initial states,
every affine trial and each optical check are retained in
[patterns.json](results/patterns.json). The cost ledger is
[patterns_budget.json](results/patterns_budget.json). Prior pilot products and
their failed witnesses remain unchanged. Reproduction requires the pinned
DE440s kernel and the pilot raw validation trace; regenerate the latter with
the commands in [pilot.md](pilot.md) if absent, then run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.patterns_run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.patterns_budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```

The experiment stops here. It demonstrates a sampled local region in which
natural differential motion, usable sail force and finite optical coverage
coexist without service thrust. It leaves a concrete terminal-state defect
for the next optimization. A bounded handover/return solve can now minimize
recurring impulse from those actual states while keeping simultaneous ray,
encounter and attitude constraints. Adjacent-pattern interfaces also need
their own joint geometry and illumination check before any local pattern is
reused in a global covering model. Neither extension, a continuous coverage
certificate nor a feasible fleet operating-cost result is claimed here.
