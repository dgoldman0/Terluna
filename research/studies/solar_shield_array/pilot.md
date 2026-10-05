# Coverage-led orbital pilot

## Scope agreed before computation

The objective is to locate a useful search region for continuous protection of
the four-lunar-radius target at lower recurring propulsion demand than the
238.35 TW held-aperture comparison. The first checkpoint is a necessary-condition
screen and a diagnosis of the constraints that prevent accepting its output.
The optical base remains 0.05 kg/m² and every physical tile is 10 km square.
The starting remote head was `fe5aac4cf3d141ee5aa99e65f71c4759b45d5a01`.

Starting radius, orbit-plane direction, phase and family population are design
variables. This pilot first propagates separate common-epoch initial states
with the filter's radiation pressure contributing to their acceleration.
It uses the existing ideal service/coast attitude law with Earth/Moon beam
guards. It does not impose Sun-following circles or cancel photon force.

The initial library has three radii (12,000, 15,000 and 20,000 km), six
meridional-plane directions and three plane offsets. A coarse 24-phase
library is followed, within budget, by 48 phases per family. All trajectories
start on 2026-10-04 TDB and run for three days. Different phases are separate
IVPs at actual ephemeris dates. The horizon tests approximately one or two
returns in the inner families; longer recurrence remains to establish.

At simultaneous dates, finite-square shadows are intersected with a partition
of the four-radius disk for several points on the finite Sun, including its
limb. A fractional covering program minimizes total optical inventory subject
to sufficient **summed shadow area in each cell at every sampled date/source**.
Each column is a uniformly populated set of independently propagated phases.
The program permits optimistic repacking within each cell and fractional
replication of the sampled set. Overlap therefore remains an explicit
relaxation. Its optimum is a restricted-library area-capacity bound, with phase
quadrature error; it is neither a realizable placement nor a universal lower
bound over all orbits. Population changes alone cannot establish ray coverage.

Binding cells and dates guide an additional plane family or a finer geometric
constraint if the first library leaves empty rows. The executed generation
steps and changes in bounds will be recorded. A sparse library that misses a
cell is a local search failure. It cannot rule out other starting states.

## Resource and acceptance gates

Use one calculation process, one BLAS/OpenMP thread, at most 4,096 initially
propagated members across the two libraries, at most 256 additional targeted
members, 1 GiB of raw arrays, 2 GiB RSS and 20 minutes of numerical wall time.
Allow at most two covering refinements, 10,000 covering rows and 256 columns.
An independent small replay may use at most 32 members and five minutes.
Record consumed resources and stop at the first informative bound or failed
gate. Repository checks have a separate five-minute allowance per checkpoint.

The screen must preserve the optical-area units, target area and full return
inventory. Refine temporal/source sampling and phase quadrature before
interpreting a numerical bound. Resolve positions within the inherited 50 m
allowance on an independent tighter replay; report any failure. Evaluate beam
margins and commanded attitude changes along the complete trajectories.
Conservative finite-tile swept encounter and mutual-shadow guards identify
incompatible simultaneous realizations. Shadow conflicts require coupled-force
repropagation before acceptance. A selected fractional multiplicity is never
treated as a physically colocated fleet.

Acceptance of an operating candidate requires individual tile identities,
finite-Sun ray unions, physical separation, safe reflected cones, realizable
attitude histories, consistent mutual illumination and force, and complete
service/handover/return/recovery costs. Search for missed rays and encounters
between samples and add the witnesses as constraints. Sampled success and a
continuous certificate remain separate evidence states. If these gates are
unfilled at this checkpoint, report the bound and the specific missing gate;
do not publish a feasible fleet or a recurring power estimate.

A subsequent control solve can vary initial states, assignments and smooth
correction arcs together, using multiple shooting or successive convex
approximations with replay. It will minimize fleet impulse and electrical
energy, subject to the acceptance constraints above. This stage is conditional
on the pilot and is not claimed executed by this plan.

## Cost comparison

Report inventory against several explicit fleet-average correction allowances
in m/s per day. These are conditional budgets, not measured maneuver
recurrence. Include a sensitivity for hardware mass and propulsion peak duty;
hardware changes would require a new coupled sail-dynamics solve. Account for
all circulating members, including the return hemisphere.

The 238.35 TW benchmark is the earlier annual, mass-closed, ideal-optical
held-aperture calculation at 50 g/m²: 30 km/s exhaust, 70% efficiency, 45° cant,
300 W/kg power hardware, 25% peak allowance and seven days of propellant.
It excludes storage and complete vehicle engineering. The same propulsion
conversion applies to this pilot's conditional impulse budgets. Propulsion
load, electrical collection capacity and the 100 TW delivered-power ambition
are separate quantities; collection and delivery remain unevaluated.

## Checkpoint record

This initial commit records the experiment before execution. Numerical results,
actual checks, useful failures and the decision about expansion will follow
here and in a compact, provenance-bearing result product. Main stays unchanged.

Before the plan commit, `make check` reported 617 Python passes, 55 skips and
the same 13 input/environment failures recorded at the reviewed head. The
layer check passed. Later Make targets were not reached. The pinned DE440s
kernel was downloaded and its recorded size and SHA-256 verified.
