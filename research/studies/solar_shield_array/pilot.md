# Coverage-led orbital pilot

**The pilot leaves a potentially useful propulsion budget, while continuous
protection remains unestablished.** The refined capacity allocation contains
17.29 million reference tiles, or 86.43 trillion kg of optical base mass.
At that inventory, a conditional 3 m/s per day of fleet-average correction
costs 114.25 TW with the additional-mass and peak-duty assumptions below.
The 238.35 TW comparison is reached at 5.97 m/s per day. Actual correction
recurrence has not been measured.

The first useful stop is the distinction between this capacity allocation and
an operating fleet. Refining one coverage bottleneck exposes empty subregions
in the preceding allocation; reweighting raises its inventory. An independently
replayed selected member also exceeds the assumed attitude-rate cap. Separation
and mutual illumination remain unresolved. The search is paused at these
specific constraints, with the numerical results preserved.

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

## Early checkpoint

The experiment plan was pushed as `89eca3edeaec8e7fa49a97c191f0e4b9eb97abf3`
before numerical execution. The current remote head matched the reviewed
commit; no newer remote work needed reconciliation. Main stays unchanged.

Before the plan commit, `make check` reported 617 Python passes, 55 skips and
the same 13 input/environment failures recorded at the reviewed head. The
layer check passed. Later Make targets were not reached. The pinned DE440s
kernel was downloaded and its recorded size and SHA-256 verified.


## Executed search and coverage refinement

The two libraries propagate **1,296 and 2,592 separate states** at their common
initial epoch, with 24 and 48 phases per family. The latter contains the former
phases as an exact nested subset. There is no time-shifted orbit replay or
interpolation between independently labelled phases. All members experience
lunar and external point-mass gravity, the filter-spectrum specular photon
force, natural eclipses and the ideal beam-guard attitude. A 3×3 surface
quadrature includes the finite tile's gravitational acceleration. Optical
force includes the clear-area fraction. Translational electric thrust is zero
in these isolated-member library propagations.

Over the three days, lunar distance stays between **11,857.7 and 21,356.2 km**.
The force model's largest sampled finite-extent correction is
6.19×10⁻⁹ m/s². This short interval does not establish the families' long-term
boundedness or their future service availability as the Sun direction changes.

The protected disk is divided into three equal-area annuli and eight angular
sectors. Each row of the covering program refers to a particular cell, actual
reception date and point on the solar disk. It asks that the **sum** of clear
square-shadow areas, after removing natural Earth eclipse, be at least the
available cell area. A square can contribute to neighbouring cells. The
perspective shadows include first-order tile and lunar motion during ray
flight. The target follows the earlier transverse receiver-disk convention;
the reflected-beam exclusion uses the full protected sphere. Central and limb sources impose separate constraints; they are not
used as an irradiance quadrature.

If a column's sampled phases have mean shadow intersection area
\(a_{cgts}\), its population \(N_g\) satisfies

\[
  \min \sum_g N_g, \qquad
  \sum_g N_g a_{cgts} \ge A_{cts}.
\]

This permits fractional replication and perfect rearrangement within each
cell. It loses the union geometry, individual assignment and interactions.
An optimum is a bound for the specified capacity columns. Finite phase
quadrature and changes in shadow-dependent force prevent treating it as a
universal physical lower bound. Increasing a weight does not create an
individual placement or distribute coincident copies over missing rays.

| Executed capacity model | Inventory, million tiles | Interpretation |
|---|---:|---|
| 24 phases, uniform families | Infeasible | Outer cell 22 is empty at day 1.5 for all five source points |
| 48 phases, six-hour dates | 16.486 | The added actual phases remove that particular empty-row failure |
| 48 phases, three-hour dates | 17.266 | More simultaneous dates increase the capacity requirement |
| Four independently populated phase bands per family | 14.700 | 216 population variables; initial-phase distribution can change |
| Phase bands passing sampled beam/rate guards | 16.074 | Reject columns exceeding the assumed 0.1°/s secant-rate limit |
| Add finer cells and solar-limb directions at day 2.5 | **17.286** | Reoptimized allocation after the strongest dual-priced date is refined |

The last step splits all eight sectors at that date and increases the solar
limb directions from four to eight, adding 432 constraints. Under the preceding
weights, **131 new constraints fail**, including a cell with zero capacity.
Reallocation raises inventory by **7.54%**. The final program has 3,432 rows and
216 columns. Its primal/dual agreement and column prices are checked. It passes
the tested area-capacity inequalities to numerical tolerance.

These changes do not demonstrate spatial or phase convergence. Other dates
retain the coarse partition, and every cell still admits subcell holes.
The coarse empty-cell failure is preserved as a phase-library limitation.
The finer phase library was sufficient for the capacity master, so the
optional targeted plane-generation branch was **not executed**. The completed
adaptive step used the coverage dual to choose additional constraints and
then changed the phase populations. No millions-member fleet was instantiated.

## Verification and rejected acceptance gates

The 48-phase and nested 24-phase trajectories differ by at most **12.52 m**.
An independent DOP853 replay of 12 selected/witness members reduces the maximum
step from 600 to 120 seconds and tightens relative tolerance from 2×10⁻¹¹ to
2×10⁻¹². Its maximum position difference is **5.08 m**, inside the 50 m placement
allowance, with a 0.263 mm/s maximum velocity difference. This is a numerical
check on those members and dates. It does not bound all members' errors.

All sampled ideal reflected cones clear the inflated Earth and protected-Moon
exclusions; the smallest library margin is **0.0573°**. These guards include
finite solar size, tile corners, body motion and the inherited pointing
allowance. Safe intermediate attitudes are still a separate constraint.

A ten-second replay around the steepest selected coarse attitude interval
finds **0.12096°/s** for member 1546, above the pilot's **0.1°/s** scenario cap.
This member belongs to an allocated phase band in the final refined solution.
The normal-plane rate in the same window is only **0.06467°/s**: the specified
square roll contributes to the failure. Free roll and a continuous attitude
schedule therefore deserve a local correction before increasing translational
thrust. Torque, angular acceleration, membrane dynamics and actuator mass were
not modelled, and a different cap requires an engineering basis.

The rate diagnostic minimizes over the eight equivalent orientations of an
unlabelled, two-sided ideal square. An initial reduction counted a 180° basis
flip as a slew; that reduction was corrected using the hash-verified existing
trajectories. The force code and trajectories were unchanged. Even the corrected
minimum-rotation diagnostic supplies only a necessary rate test. Real sided
coatings or labelled hardware may require a longer rotation.

The final allocation's **1,224 sparse phase representatives** receive a
600-second swept sphere guard, enlarged by an assumed 0.06 m/s² acceleration
bound and the 100 m clearance allowance. It flags **264 possible encounter
pairs**. Exact finite-square tests at saved instants find no intersections
among these representatives. This combination leaves encounters unresolved:
the conservative flags may be avoidable, and saved instants can miss a contact.
There is no inherited radial-shell separation certificate for these free paths.

The swept solar-shadow guard flags **38,812 possible interacting pairs**.
These are conservative possibilities, not a count of proved occultations.
Mutual shadows were not inserted into the isolated-member IVPs. The final
population weights therefore cannot be promoted to a simultaneous force and
coverage solution. No operating fleet passes the acceptance gates at this
checkpoint, and none of the estimates is a cost for an accepted fleet.

## Inventory and recurring impulse

The [conditional budget](results/pilot_budget.json) uses the same exhaust,
efficiency and cant as the held comparison. The main sensitivity adds fixed
hardware equal to 20% of optical mass, assumes thrust occupies 10% of each
member's time, and closes installed power hardware and seven days of propellant
against their own extra mass. The 20% allowance and duty factor are explicit
scenarios. This mass feedback has not been reinserted into the sail IVPs.

For mass \(M\) and fleet-average accumulated correction rate \(\dot v\),

\[
 P=\frac{M v_e\dot v}{2\eta\cos\theta}, \qquad
 \dot m=\frac{M\dot v}{v_e\cos\theta}.
\]

Every circulating member contributes to \(M\), including its off-duty return.
The assumed \(\dot v\) must eventually include service, acquisition, handovers,
avoidance, return adjustments, recovery and their recurrence. A count of cheap
individual maneuvers supplies none of that schedule by itself.

| Fleet-average correction | Propulsion at 17.286 million tiles | Propellant |
|---|---:|---:|
| 1 m/s per day | 36.93 TW | 57,445 kg/s |
| 3 m/s per day | 114.25 TW | 177,728 kg/s |
| 5 m/s per day | 196.57 TW | 305,782 kg/s |
| 10 m/s per day | 427.69 TW | 665,292 kg/s |

| Assumed inventory | Power at 3 m/s per day | Correction rate reaching 238.35 TW |
|---|---:|---:|
| 17.286 million tiles | 114.25 TW | 5.968 m/s per day |
| 25 million tiles | 165.24 TW | 4.243 m/s per day |
| 35 million tiles | 231.34 TW | 3.087 m/s per day |

The larger inventories are sensitivity scenarios. The first is the refined
capacity estimate with the limitations above. More inventory can make coverage
and repair easier, while reducing the allowable average correction per unit
mass. With no fixed 20% increment and continuous thrust, the first inventory's
scalar break-even is 7.75 m/s per day; with 10% thrust duty it is 7.04 m/s per
day. These values bound an engineering target under stated bookkeeping.

At 3 m/s per day in the main scenario, propulsion hardware is rated for
**1,428 TW installed peak capacity**, including margin. That rating sums each
member's hardware requirement; thrust need not occur simultaneously. The
calculation provides no collector design capable of supplying that load.
Collection capacity, transfer losses, storage, habitat loads and **100 TW
of delivered electricity remain unevaluated**. The model's optical momentum
budget also provides no collected electrical-power credit.

## Resource use, reproduction and stopping decision

The two library runs took 118.1 and 226.6 seconds. Correcting the attitude
reduction reused the first trace in 9.8 seconds. Four bounded validation
passes inspected the allocation/witness choices and checked the final optical
reduction. The final capacity reduction uses the same tile-local incident
vector for square roll as the finite-extent force calculation. Together, the
numerical work took about **ten minutes**, with one calculation process and
BLAS/OpenMP set to one thread. The highest recorded RSS was **600 MiB**. There
were 3,888 library IVPs across the nested grids and 17 distinct independently
replayed members across the four validation passes (48 replay instances).
No extra trajectory-generation batch was needed. Raw compressed products use
about 150 MiB under ignored `research/runs/solar_shield_array/pilot/`.

The products pin source hashes, named shared constants and raw bytes. The
force-model source chain and exact initial states are checked before a reduction
can reuse trajectories. Kernel size and SHA-256 were verified on loading.
The rate-reduction repair records its preceding producer and trace identity.
All new code and compact products belong to this branch; existing historical
models and results are preserved.

```sh
python -m protection.dynamics.ephemeris --download
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.pilot_run --phases 24
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.pilot_run --phases 48
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.pilot_validate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.pilot_budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest protection/dynamics research/studies/solar_shield_array -q
```

The useful result is an explicit inventory/impulse trade and a coverage-driven
allocation method that detects its own coarse-geometry failures. It supplies
**room to investigate a lower-power operating region**, conditional on keeping
the required recurring correction within that room. It does not establish that
such a feasible region exists.

Before expanding the search, the next experiment should give roll and normal
schedules bounded continuous dynamics, construct a small collision-compatible
pattern with coupled mutual illumination, and test exact rays through its
handover and return. Its missing-ray and encounter witnesses can generate new
initial states and control arcs. An energy-minimizing shooting/collocation or
successive-convexification solve belongs after that pattern is concretely
specified. **No constrained electric-control optimization, recurring recovery
schedule, continuous ray certificate or long-horizon fleet propagation was
executed in this pilot.** Those steps remain proposed for the next decision.

The final repository check reports **624 Python passes, 55 skips and the same
13 input/environment failures** seen at the starting head. All **79 study/domain
tests** pass, including area conservation, phase regrouping, LP dual closure,
square-orientation equivalence, budget feedback and saved-product provenance.
The layer check passes; Make stops before JavaScript, historical-provenance
and ensemble targets. [checks.json](checks.json) records those scopes and the
actual source/product hashes. `git diff --check` also passes.
