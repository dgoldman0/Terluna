# Joint departure and next-service investigation

Continue the reviewed `c47c96040fc386b1d980c82a036ec91830f47390` checkpoint on
`research/solar-shield-habitat-array`. The author waived the preliminary remote
head/concurrency check. Main stays unchanged. This study treats the moving
array as an optical EUV/ionizing-radiation filter and climate dimmer; primary
magnetic protection is outside this experiment.

**Executed result:** no candidate meets the combined safety, installed-power,
supply and low-energy requirements. A high-control separated corridor returns
all 361 tiles at hour **45.77978** and completes another six-hour shielding
passage in two independent integrations. Its uncertainty-adjusted conditional
surface-clearance bound is **900.28 m**. It nevertheless overloads **145**
installed actuators and costs **1.53538 PJ** through the next service endpoint.
Its first twelve hours cost **339.46 TJ**, versus the preserved **51.859 TJ**
loaded baseline, so it does not retain the earlier energy gain. The cheaper
executed alternative fits installed actuator power but stops at the 100 m
clearance floor at hour **12.75076**. No operating return, repeatable formation
cycle, handover or fleet budget is accepted. The search is paused at this
local result, before a larger fleet campaign.

## Budget and candidate families

Initial allocation: **1,800 CPU seconds**, one numerical process and thread,
**2 GiB address space**, **512 MiB new raw output**. Required repository checks
have a separate five-minute allowance before each commit. Inspection/code
preparation are separate from simulation time; numerical partial attempts are
charged. Review usage before each stage. Do not silently expand the search.
The historical full tighter regression alone took about 1,960 seconds; use
an explicitly cross-checked faster force implementation and checkpointed
regressions rather than rerunning every historical producer.

Preserve all 361 physical 10 km squares, 9.89 km clear apertures, 50 g/m²
optical allocation, four layers spanning 36 km, 100 m surface floor, 50 m
replay gate, 0.001 m/s² translational cap and 0.1 degree/s attitude cap.
Retain the four-lunar-radius optical target, six-hour moving 20 km receiver
service, stored transmission products and the 13.7950037% ideal redirected
spectral allocation. Coverage is geometric evidence conditional on that
unchanged spectral model; it does not validate angular thin-film response.

Three meaningfully distinct trajectory families, at most twelve inexpensive
full-population proposals and three coupled candidates:

1. Earlier velocity shaping, beginning after the six-hour assigned service,
   with the previous arrival packing as comparator. Optimize early escape and
   the next service endpoint together; do not treat changed departure states
   as free initial conditions.
2. Changed service packing and assignment: vary depth/shear/relative velocity,
   preserve ray service, and test staggered approach times. Include the motion
   and energy of every change, and continue into the next six-hour service.
3. Low-dimensional off-service steering: common tilts and depth/row-group
   independent normals. Independent normals require generalized finite-square
   shadow, photon force, moments and surface-contact tests before any reliance.
   Document restrictions and retain failed proposals, including power failures.

Reserve a tighter independent replay for the best complete candidate or a
limiting executed obstruction. A reduced trajectory or imposed kinematic path
can propose controls and costs; it cannot pass dynamics acceptance. Publish
this checkpoint before numerical search and results afterward.

## Baselines and accounting gates

Read `cycling.md`, `packing.md`, `closure.md`, `return_repair.md`,
`collection.md`, `electromagnetic.md`, their producers and result products.
Verify point-tile cycling separately from finite-formation closure. Reproduce
the hardware-loaded passage and original return witness (325/347, hour
13.893535). Preserve the strongest previous repair (151/186, hour 15.390986,
160 installed-power overloads) as an additional comparator.

The 82.13% reduction compares optical-mass prefixes. Subsequent comparisons
use the **2.26330 billion kg** loaded formation and **51.859 TJ** twelve-hour
account. Attitude torque and energy remain chargeable. No finite-formation
return, next service, handover or global operating budget was accepted at the
starting checkpoint.

Require all-member dynamics with actual-date DE440 forces, finite Sun,
mutual shadowing, finite-square clearance, adaptive close-approach/transition
checks, beam exclusions, attitude rates, torque, individual actuator ratings
and available electrical power. Include numerical disagreement in clearance.
Changed hardware must be carried from its first installation and replayed.
Track initial service, departure and re-entry coverage separately; unassigned
or uncovered times need replacement tiles and receive no fleet credit.

The supply candidate is local Sun-tracking generation and a short buffer,
using the explicit engineering assumptions in `electromagnetic.md`. Price
collectors, PMAD, buffers, conductors and heat rejection; propagate their mass
and photon momentum before acceptance. Count first-intercept light once,
separating redirected optical power, gross generation, bus electricity,
storage state/losses, propulsion/attitude energy and heat. No external beam
capture or transmitted-band collection is free. Preserve unknown accounts.

Report local executed results first, terminal state/attitude, actual next
service and its endpoint's support for repetition. One return is distinct
from recurrence. Keep **238.35 TW** held-array comparison and **23.835 TW**
operating target. **17.29 million tiles** remains a capacity relaxation. No
failed-prefix extrapolation is an achieved fleet demand. Pause before a
substantially larger fleet campaign.

The early required `make check` passed the layer check and reported 718 Python
passes, 55 skips and 13 existing failures in 29.03 s. Missing pinned spectral
inputs, generated climate configuration and process lookup remain the same
baseline blockers. Later Makefile targets were not reached.

## Recorded search decisions

The exact all-pair compiled rectangle union reproduces the original force,
first moments and collection without changing historical sources. The full
loaded regression costs 84.4 CPU seconds, leaving room within the original
1,800-second allowance.

Six initial endpoint proposals are followed by three corrected staggered-slot
proposals, two separated corridors and one grouped-attitude proposal, within the twelve-proposal
ceiling. The initial stagger screen disabled the final burn for earlier slots
and used a Taylor backcast; its failures apply to that additional restriction.
The corrected screen moves the last arc before the earliest slot and uses
independent dynamics to backcast each group's target. All original failures
remain preserved. The corridor family explicitly pays to capture departure
velocities and restore a next-service velocity field; geometry alone cannot
accept its inverse-dynamics controls.

The wider orbit enters lunar eclipse near hours 25–27. A joint body/array
area-ray model is added for that interval, with full/empty analytic limits
and spatial/source refinement. The existing uneclipsed model still supplies
an independent regression. A one-minute power buffer must be judged against
the measured eclipse load; earlier collection cannot be credited as available
maneuver power without a chronological storage and mass account.

The separated-corridor screen uses two variants; the twelfth proposal is a
four-depth-group early sail pulse. This replaces the unused third corridor
variant within the same twelve-proposal limit. The 1.08 expansion corridor
has positive sampled geometry and is selected as a deliberately high-control
comparator; its existing actuator ratings already fail the inverse screen.
Any subsequent propagation is an ideal-force diagnostic until those ratings
and actual supply pass. No resized hardware is accepted from an unpropagated
mass estimate.

For the selected independent replay, use 16 instead of 8 Sun sources, 32
instead of 16 area nodes per axis during partial eclipses, rtol 2e-12 instead
of 2e-10, and a 90 instead of 120 s maximum step. This spends the remaining
budget on full-cycle source/area refinement while retaining a distinct time
and tolerance check. The source-bound wrapper records those actual settings.

## What was reproduced and what causes the old failure

The source reports establish different evidence states:

| Source | Verified meaning |
|---|---|
| [cycling.md](cycling.md) | Three-year point-tile cycling and a separate 8.90548-day return, independently replayed within 16.39 m, concern point trajectories. Neither is a finite-tile formation cycle. |
| [packing.md](packing.md) | The 82.13% local control-energy reduction is the optical-mass comparison: about 41.455 versus 231.917 TJ. Attitude control remains nonzero. |
| [closure.md](closure.md) | With 2,263,298,487 kg total mass, about 25.4% above optics, six hours of service plus six hours of departure still need zero electric translation. The loaded twelve-hour cost is 51.859 TJ and the reviewed minimum sampled surface separation is 1,719.63 m. |
| [return_repair.md](return_repair.md) | The original encounter is near hour 13.893535; the strongest repair moves failure to hour 15.390986 and overloads 160 tiles. Neither completes a return. |
| [collection.md](collection.md), [electromagnetic.md](electromagnetic.md) | Intercepted and redirected light are source inventories. Generation, chronological bus delivery, storage, hardware mass and added photon forces require separate calculation. |

[joint_regression.json](results/joint_regression.json) freshly integrates the
loaded initial service, original departure and original return without state
resets. Minimum separations are 6,698.74 m through service, 1,720.33 m through
departure, then 100 m at hour 13.89335390 on tiles 325/347. The event differs
from the reviewed tighter run by 0.6523 s; maximum position disagreement is
9.0419 m. The original six-hour finite-Sun service check passes.

The new compiled parallel-square implementation uses the full blocker set and
the exact rectangle union, including first spatial moments. At six actual
saved states it agrees with the historical implementation within 4.34e-10 N
in force and 1.92e-6 N m in torque. The maximum optical difference normalized
by 1e11 W is 4.80e-12. This is an independent implementation cross-check, not
a new physical optical model or a per-tile relative error near darkness.

The decisive old encounter has body-frame centre displacement approximately
(-9,314.92, 8,885.43, 100) m. Both in-plane displacements are smaller than
the physical 10 km edge, while normal closing speed is about 0.361 m/s. The
repair's later encounter has projected centre distance about 7.84 km and
normal closing speed about 0.660 m/s. Two 10 km squares each contain a 5 km
radius disk; at that later projection, in-plane roll cannot eliminate the
overlap. These are physical-square constraints, not a point-distance or
optical-seam artifact.

At the initial epoch the minimum centre distance is 14,705.44 m: an
orientation-independent bounding-sphere gap of about 563 m. It falls to
13,983.38 m by hour 6 and 5,547.51 m by hour 12. The changing relative
velocities and depth crossings consume the initial clearance. A safe initial
packing alone cannot guarantee a safe turn or return. The restricted late
burns try to fix endpoint errors after this compression; local encounter cuts
can migrate the conflict to another pair. Individual power allocations then
restrict how rapidly the velocities can be changed. None of these witnesses
establishes impossibility for other trajectories or packings.

## Bounded family comparison

All proposals start from the actual loaded first-service endpoint. Service
before hour 6 is propagated and charged; no changed departure velocity or
arrival state is supplied for free. Endpoints include both position and
velocity, and commanded attitude transitions have finite duration. The search
tests twelve full-population proposals, rather than an unrestricted optimum.
Mean translation below is the proposal integral unless explicitly labelled
executed. The products retain solver status, attempted encounter cuts and
failed candidates.

| Family / product | Freedom tested | Result and restriction actually ruled out |
|---|---|---|
| Six initial proposals, [joint_screen](results/joint_screen.json) | Six smooth burn windows starting at hour 6; original packing comparator; 24/36 km layer spacings; identity or Hungarian assignment; arrival shifts of 0, -0.5 and -1 hour; endpoint boxes; common 75/90/105 degree off-service steering | Early original-packing proposal uses 1.61 m/s after one correction but retains conflicts; its next LP reaches the eight-second solver limit. Deep reassignment uses 2.17 m/s and selected cuts become infeasible. Deep identity reaches 3.28 m/s after six refinements and still conflicts. Initial staggered cases fail an additionally restricted last-burn/backcast model. These are local LP/iteration failures, not global exclusions. |
| Three corrected staggered proposals, [joint_slots](results/joint_slots.json) | Four arrival groups over 0.5 hour; final burn before the earliest arrival; independent backward propagation of group targets; 75/90/105 degree common steering | Endpoint proposals exist at 2.04/2.16/2.19 m/s, but all have finite-square conflicts. This corrects the first stagger screen's unnecessary restriction without making staggered arrival generally impossible. |
| Two separated corridors, [joint_bridge](results/joint_bridge.json) | Gradual capture of actual departure velocities; a spaced three-dimensional coast; explicit acquisition of next service position and velocity | 1.08 expansion, capture by hour 14, eight-hour approach: sampled gap 1.961 km, translation 21.18 m/s, 145 power overloads. 1.20 expansion, capture by hour 12, ten-hour approach: gap 3.311 km, 24.15 m/s, 222 overloads. These establish an expensive geometric route, not a low-cost optimum. |
| Four depth-group attitude pulse, [joint_grouped](results/joint_grouped.json) | Independent group offsets of -10, -3.333, +3.333 and +10 degrees, smooth between hours 6 and 18, combined with early electric controls | The reduced all-member trajectory first samples below 100 m at hour 12.6667, then samples intersection at hour 13.1667. It is rejected before costly coupled replay. This excludes this pulse/control combination only. |

The corrected staggered proposals are screened only; their staggered optical
handover is not executed or credited. New packings retain the 10 km physical
square and 9.89 km clear aperture and are tested against the same service
target. A changed depth does not relax spectral transmission or receiver
coverage. No useful solution emerged from the selected assignment/time
branches within this budget.

The cheaper early original-packing case was then propagated with coupled
shadow forces from the initial epoch. It reaches 100 m at hour 12.75076472 on
pair 263/300, with in-plane displacement (9,055.73, -7,617.85) m. Its complete
executed prefix costs **99.17197 TJ**, including the initial service, and
has zero installed-power overloads; the minimum installed/required peak ratio
is 1.8404. This satisfies the actuator-rating screen, not physical generation
availability. [joint_early_same_coarse](results/joint_early_same_coarse.json)
and [joint_failed_account](results/joint_failed_account.json) preserve the
stopped trajectory and its optical/energy account. A tighter replay was not
spent on this already-rejected branch; no fleet demand is extrapolated from it.

## Independent attitudes and eclipse geometry

[oriented_tiles.py](../../../protection/dynamics/oriented_tiles.py)
generalizes shadow area and first moments to independently oriented squares.
It clips each receiving aperture against the perspective cone of every
blocking square and the blocker's upstream half-space; this handles a tilted
blocker crossing the receiver plane. It unions shadows before integrating
photon force and torque, with each tile's own frame. Physical contact uses
vertex/face, edge/edge and edge/face intersection tests on the full 10 km square.
Tests cover both illuminated faces in the parallel limit, analytic parallel
clearance, transverse intersections, and tilted shadows against independent
direct ray/plane intersections, including first moments.

The four-group family uses this full geometry for clearance and selected
shadow/load snapshots. Its inexpensive trajectory uses unshadowed finite-Sun
forces. At the three checked states, the acceleration difference from the
oriented-shadow evaluation reaches 1.01e-5 m/s²; this comparison also changes
source quadrature and is not an isolated shadow-error estimate. The parallel
force limit agrees with the fast implementation within 5.73e-12 N. Nominal
attitude impulse rises from 0.6824 to 0.7348 m/s and maximum rate is
0.01245 degree/s; sampled radiation torque reaches 98,385 N m. Because it
already intersects, this proposal has no accepted coupled attitude/actuator
or hardware account. Arbitrary per-tile control, actuator placement and
flexible-film response remain outside this reduced family. No favorable
independent-attitude result is relied upon without replay.

The long corridor encounters lunar eclipse. [eclipse_parallel.py](../../../protection/dynamics/eclipse_parallel.py)
therefore tests body and tile occlusion on the same incoming rays, rather
than multiplying their separately averaged visibility fractions. Full and
empty body-shadow limits use analytic tests; partial shadows integrate tile
area. Analytic limits and an overlapping half-shadow case are tested. The
coarse/fine executions use 8/16 Sun sources and 16/32 area nodes per axis in
partial eclipses. Geometric ray overlap, force and intercepted energy share
this union, so light behind two occluders is not counted twice.

## Best complete geometric execution

The selected 1.08 corridor begins with the unchanged loaded initial states,
gradually captures the relative velocity by hour 14, holds expanded relative
positions, and approaches the next service over its final eight hours. It
turns toward feathering during hours 7–11 and returns smoothly to Sun-facing
service. Its inverse-model accelerations are saved, then applied open loop
to all 361 tiles in fresh nonlinear integrations. The actual trajectory is
never replaced by the prescribed path at a boundary. The fine replay starts
again at the initial epoch with tighter tolerances, altered time steps and
refined Sun/area quadrature; it does not refit controls.

| Executed stage | Time, hours | Minimum sampled separation | Maximum coarse/fine position difference including reconstruction | Conditional interval clearance after twice that disagreement |
|---|---:|---:|---:|---:|
| Initial service | 0–6 | 6,698.74 m | 1.56 m | 5,570.81 m |
| Departure and return | 6–45.77978 | 1,956.27 m | 20.79 m | 900.28 m |
| Actual next service | 45.77978–51.77978 | 6,693.69 m | 26.76 m | 5,517.17 m |

[Coarse](results/joint_corridor_108_coarse.json),
[fine](results/joint_corridor_108_fine.json) and
[audit](results/joint_audit.json) products distinguish these passed geometric
gates from failed operational acceptance. The all-pair interval bound uses a
0.15 m/s² acceleration envelope and the 0.1 degree/s attitude bound, with
recursive refinement down to 0.5 s if needed. It is a conditional numerical
bound, not an interval-arithmetic proof. Additional minimization around the
closest transfer sample finds 1,956.2673 m at hour 38.00568 on pair 2/3.
All three stages pass the existing 50 m position-disagreement gate.

The fine arrival differs from its prescribed target by at most **48.236 m**
and **0.0013064 m/s**. Maximum attitude discrepancy from Sun-facing is
6.73e-8 rad. Those actual states continue for six hours with zero electric
translation: service is not inferred from the endpoint alone. The historical
point-return velocity accuracy is not promoted to an unestablished
finite-formation terminal tolerance.

The next passage passes 49 dates with 25 finite-Sun sources each, the full
tile dimensions and seams, plus 553 additional source/date probes. The worst
found receiver-boundary margin is **5,322.08 m**, with zero missed area and a
0.1883 m circle-polygon sagitta allowance. The moving 20 km receiver remains
inside the four-lunar-radius target. These are numerical coverage tests,
conditional on the inherited spectral filter; continuous source/time and
material angular-response certification remain open.

Coverage is sampled through departure and re-entry as well as service.
The formation has **no assigned shielding credit from hour 6 to 45.77978**;
another population must protect its receiver throughout that interval.
Incidental interception while leaving, behind the Moon, or approaching does
not count as an executed handover. No replacement population was propagated.
At the second service endpoint, frame-normalized relative states differ from
the first endpoint by up to **237.16 m and 0.23793 m/s**. No reset is applied
and no following departure is executed. This is one geometric return and
second passage, not a demonstrated repeating cycle.

## Complete executed energy and the failed hardware gates

The following account uses the **loaded** mass, 30 km/s exhaust, 70% thrust
efficiency, the inherited 45-degree plume allowance and electric torque
couples. Attitude loads include rigid-body angular acceleration, finite-area
gravity and shadow moments. Detailed per-tile powers and optical integrals
are retained in [joint_account.json](results/joint_account.json) and hashed
raw arrays. Electrical integration is at 30 s; optical/disturbance comparisons
use 8 sources/120 s and 16 sources/240 s, with area refinement in eclipse.
Largest combined time/source energy discrepancy is 0.1083%, and the largest
optical discrepancy is 0.0709%; these are comparisons, not certified error
bounds or separately isolated quadrature tests.

| Stage | Propulsion and attitude electricity | Mean first-intercept bolometric-equivalent power | Mean ideal redirected-band power |
|---|---:|---:|---:|
| Initial six-hour service | 2.64803 TJ | 39.5703 TW | 5.45872 TW |
| Complete departure/return | 1,530.03811 TJ | 7.06454 TW | 0.974553 TW |
| Next six-hour service | 2.69236 TJ | 39.6956 TW | 5.47601 TW |
| Entire executed 51.77978 hours | **1,535.37850 TJ** | Integrated **2.723835 EJ** | Integrated **0.375753 EJ** |

Initial service plus one return costs 1,532.68614 TJ; the additional next
passage accounts for the difference. Sequence mean power is 8.23669 GW,
mean translation impulse 21.18365 m/s and attitude equivalent 1.18530 m/s.
Ideal exhaust is 2.38837 million kg and modeled propulsion waste heat is
0.460614 PJ. Propellant depletion, about 0.1055% of carried mass, is not fed
back into this constant-mass trajectory. Structural, thermal-deformation and
actuator-layout limits are not closed.

Maximum translational acceleration is 0.00050138 m/s², maximum command rate
0.011719 degree/s, and minimum sampled Earth beam clearance about 30.31 degrees.
These pass the inherited kinematic/beam limits. Maximum individual electrical
load is **112.511 MW** and maximum required torque **2.41239 MN m**. Aggregate
installed rating is 29.18955 GW, but 145 tiles overload: the worst installed
rating is only 47.58% of required peak. Summed individual required peaks are
24.33715 GW and simultaneous aggregate peak is 24.33047 GW. Aggregate capacity
cannot repair a particular tile's actuator limit, even with ideal bus routing.

Maintaining existing allocations and adding the 25% design margin would
require **19.6825 million kg** of additional actuator power equipment under
the inherited 300 W/kg assumption. This hardware is priced only; it is neither
installed nor replayed. The executed accelerations above ratings are explicitly
an ideal-force diagnostic, so operational return acceptance is false.

The initial deep packing and initial six-hour zero-translation service are
retained, but the new corridor starts costly preparation immediately afterward.
Its first twelve hours use **339.461 TJ**, **6.546 times** the loaded
51.859 TJ comparator. It consequently fails the requested preservation of
the low-energy departure gain, even apart from power and generation. The
geometric success must not be reported as meeting the user's full goal.

## Local generation and chronological storage sensitivity

The following is an **uninstalled optimistic supply sensitivity**, using the
generation/storage assumptions from electromagnetic.md: 30% PV conversion,
95% PMAD, 99% cable delivery, 300 W/kg array allocation, 200 Wh/kg and
1 kW/kg storage, 95% charge/discharge efficiency and 80% usable depth. Each
tile receives its own load and battery dispatch; generation has a 25% peak
margin and local Sun tracking. Collectors are assumed to receive fresh light
except for the calculated body eclipse. Their placement, mutual shadows and
combined optical geometry have not been designed. The artifact's
`local_supply_upper_bound` label denotes this optimistic model, not a rigorous
bound for all possible collector locations.

First-intercept optical energy in the table above belongs only to the
executed filter surfaces. The additional hypothetical collector light is
never added to it or to redirected filter energy. Actual generation,
delivered bus power, export and full operating-cycle cost remain **null** in
the result product. Redirected light is not treated as available electricity.

| Conditional hardware or availability | Result |
|---|---:|
| Collector area / generation hardware mass | 78.429 km² / 147.617 million kg |
| Additional conductors / battery heat rejection | 9.265 / 8.363 million kg |
| 60-second design buffer mass | 24.337 million kg |
| Energy unserved with that buffer | **27.067 TJ**, first at hour **25.5083** |
| Sampled partial/total eclipse interval | Hours **25.025–26.8917** |
| Eclipse-sized storage mass / nominal capacity | 73.802 million kg / 53.137 TJ |
| Assumed absorption photon force on added collectors | 320.447 N total |

The short-buffer mass is actually set by its power limit and supplies more
than a minute of energy; it still fails the eclipse load. The larger storage
dispatch starts and ends with 42.510 TJ usable energy and leaves only 0.009 J
numerical unmet demand in this optimistic supply model. The largest balance
residual is 0.164 J. This does not establish a realizable supply: the listed
mass, absorption momentum, reflected/thermal momentum, local shadows and
thermal geometry have not been propagated, and the actuators still overload.
No resized supply candidate is accepted.

Over the complete sensitivity sequence, available bus energy is 5.43420 PJ;
PV heat 11.55597 PJ, PMAD loss 0.288899 PJ, cable loss 0.054891 PJ and storage
loss 0.004363 PJ. Bus curtailment is 3.89446 PJ. Curtailment requires generation
regulation or an explicitly cooled dump and is not export. The apparent
positive net energy cannot deliver maneuver power through the eclipse with
the short buffer. Thermal and collector hardware are consequently measured
design sensitivities, not a completed power system.

## Scope, reproducibility and next bounded test

Keep **238.35 TW** as the held-array comparison and **23.835 TW** as the current
operating target. The **17.29-million-tile** result remains an optimistic
capacity relaxation. This investigation makes no fleet-demand extrapolation:
cadence, inventory, simultaneous full-target coverage, replacement assignment
and handover have not been demonstrated. It does not combine the optical
array with the primary magnetic shield or change the pinned spectral products.

The smallest useful next test is a **single local, power-constrained release
of the known safe 1.08 corridor**, optimizing relative velocity earlier and
reducing its long position-holding segment toward free flight. Start from
the saved loaded epoch and six-hour states, constrain every member's actual
rating, preserve the clearance margin and next optical passage, and suppress
eclipse maneuvers before considering larger batteries. Compare its first
twelve-hour cost against 51.859 TJ as well as its full sequence. A viable
proposal then needs one all-member replay with explicitly placed collectors,
their loaded mass/photon forces and chronological supply. Another departure
from the actual second-service endpoint is the next recurrence test; fleet
handover follows only after local hardware feasibility. This is a proposed
test, not a further campaign started here.

Source-bound producers and result products record CPU, memory and raw hashes.
The twelve producer invocations (including the exact-state export) use
**1,637.97 CPU seconds**, peak RSS **824.64 MiB**, and **478.49 MiB** of new
ignored raw output. An additional 20-second allowance is charged for brief
unmetered geometry/profiling checks and compilation, leaving about 142 seconds
of the 1,800-second numerical allocation. Repository checks are accounted
separately. There are twelve proposals, two new coupled candidate trajectories
and one tighter candidate replay, plus the regression. No hidden population
reduction is used for acceptance geometry.

[joint_seed.json](results/joint_seed.json) preserves exact float64 epoch,
first-service and original twelve-hour states, hardware ratios, the selected
corridor specification and hashes independently of ignored raw files.
The full historical raw chain remains necessary for the current producer
entry points; the seed and `joint_bridge.build` provide inputs for the next
local experiment without pretending the seed is a replay. Full regeneration
of missing historical runs needs a separately budgeted restoration.

With the existing historical raw results and pinned DE440 kernel available,
run sequentially with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`:

```sh
python -m research.studies.solar_shield_array.joint_regression
python -m research.studies.solar_shield_array.joint_screen
python -m research.studies.solar_shield_array.joint_slots
python -m research.studies.solar_shield_array.joint_bridge
python -m research.studies.solar_shield_array.joint_execute
python -m research.studies.solar_shield_array.joint_replay
python -m research.studies.solar_shield_array.joint_grouped
python -m research.studies.solar_shield_array.joint_account
python -m research.studies.solar_shield_array.joint_audit
python -m research.studies.solar_shield_array.joint_execute --case early_same
python -m research.studies.solar_shield_array.joint_failed_account
python -m research.studies.solar_shield_array.joint_seed
```

The result ledger includes twelve products from the twelve invocations above.
The exact count and independent resource entries are retained in
[joint_checks.json](joint_checks.json).
One post-run source maintenance change replaced the numerically identical
lunar-radius literal in `joint_screen.py` with the shared constant. The checks
file retains both source identities and every affected parent-product hash;
the old source is reconstructible and no numerical payload was recomputed or
changed. This is explicitly provenance maintenance, not a claimed rerun.

Targeted domain/study/storage tests passed (179 tests before the two final
provenance/export assertions). Final `make check` passes the layer check
(573 files, zero violations), then reports **734 passed, 55 skipped and the
same 13 baseline Python failures** in 30.38 s. All 16 added tests pass. Missing
pinned spectral inputs, generated climate configuration and process lookup
remain the baseline blockers; subsequent JS, provenance and ensemble Makefile
targets are not reached. The exact failures and log hash are retained in
joint_checks.json. This is not a passing whole-repository check.
