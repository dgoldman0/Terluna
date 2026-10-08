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

The plan preceded numerical search. Its `make check` run reported 655 Python
passes, 55 skips and the same 13 existing failures in 37.12 s; the layer check
passed and later targets were not reached. The following sections record the
executed experiment.

After 45 proposals and two coupled trials (about 500 s total), both nominal
leaders reached the 100 m stop. The cheaper four-hour turn lasted 5.49 h;
the two-hour turn lasted 3.57 h. The measured state defects now supply one
successive correction of the cheaper schedule. Before that retry, the count
budget is revised to three coupled candidates plus one tighter replay, keeping
the original 45-minute wall, process, memory and output limits. The added
trial tests a measured force-model defect with the existing three controls;
it does not enlarge the orbit or attitude search. Preserve both failed runs.

The measured-defect correction fails its chosen separating half-spaces even
with signed impulses in all three modes. The third coupled slot is therefore
allocated to the next surviving schedule already in the 45-case screen: a
four-hour 90-degree turn about the other square edge. This decision precedes
that propagation and retains the same total resource budget. The failed
correction and both cheap partial-turn failures remain separate products.

## Schedule search and the limiting constraint

The early plan was pushed as `374039df488132f67b1cb213499032100b5510e5`.
The search executed 45 distinct schedules. Zero turn appears once. The others
use 15-, 30-, 60- and 90-degree turns about either square edge and with either
sign. The 15–60 degree cases use one-, two- and four-hour turns; 90 degrees
uses one and four hours. The delay is two hours for the shorter turns and
one hour for the four-hour turns, leaving an hour of coast in the latter.
All use one two-hour preparation pulse and the three depth/column/row modes.
All states use the shared DE440 ephemeris at actual simultaneous dates,
with epoch 4 October 2026 TDB. Service occupies the first six hours; each
departure is tested from that same service endpoint.

Fifteen proposals pass the refined approximate geometry and beam/rate gates.
Twenty-five fail the selected local separating half-spaces or acceleration
constraint, and five fail the beam/rate screen. The zero-turn case reaches a
low-impulse first fit, then fails its refined encounter constraints. These
outcomes describe this control reduction and its chosen separating branches.

The leading proposal turns through 30 degrees over four hours. Its mean
translation is **0.68134 m/s**, with **0.22726 m/s** nominal electric-couple
impulse. Their **0.90860 m/s** sum is much smaller than the previous departure's
1.40584 + 2.72711 = **4.13295 m/s** nominal account. The two-hour 30-degree
variant has the same translation and a 1.13586 m/s nominal sum. These values
exclude disturbance loads and remain proposal costs.

Coupled propagation rejects both. The four-hour partial turn reaches the
100 m separating-axis clearance guard after **5.4866 h**; the two-hour partial turn reaches it after
**3.5708 h**. No physical intersection is claimed: the integration stops at
the chosen safety boundary. Two-second swept guards flag the ending intervals,
and coplanarity searches find no actual intersections in the retained prefixes.

The first failed trace differs from its gravity-variational proposal by as
much as **3.497 km**. A successive correction uses the measured state defect
through the stop and transports the endpoint defect gravitationally over the
remaining half hour. Its 2,160 selected separating half-spaces are infeasible
in the three-mode family, including when signed inward/outward impulses are
allowed. This is a local branch/family failure; it does not prove that the
30-degree schedule is physically impossible with other controls or packing.

The final coupled slot tests the next ranked untried family: a four-hour,
90-degree turn about the other square edge. Its nominal proposal account is
0.81075 m/s translation plus 0.68243 m/s attitude, or **1.49318 m/s**. It reaches
the 100 m stop after **4.4679 h**. Its independent tighter replay is consequently
not run.

| Coupled trial | Turn | Clearance-stop time after service | Limiting pair | Minimum sampled beam margin |
|---|---|---:|---|---:|
| 24 | 30 degrees in 4 h | 5.4866 h | 186, 187 | 23.01 degrees |
| 20 | 30 degrees in 2 h | 3.5708 h | 173, 174 | 42.59 degrees |
| 41 | 90 degrees in 4 h, other edge | 4.4679 h | 34, 35 | 30.42 degrees |

An all-pairs calculation on each saved endpoint confirms that its actual
parallel-square minimum surface distance is also 100 m, to 1e-9 m arithmetic
agreement, with the same limiting pair. This uses the Euclidean norm of the
nonnegative edge/normal separation excesses. It checks the stopped endpoint;
there is no propagation past the boundary and no collision assertion.

All three remain below 0.000360 m/s² commanded translation and 0.01172 degree/s
attitude rate. The reflected-beam, rate and thrust limits still pass at the
stops. **Finite-tile clearance is the rejection constraint.** The old verified
departure remains the feasible nominal comparison. No smaller-cost departure
has been accepted, and no larger clearance margin was optimized.

## Calculated torque and the cost distinction

The torque model integrates the gravity force over a uniform physical square
and the first spatial moments of the reflected-band shadow unions over its
clear aperture. It uses the same centre-ray irradiance/incidence approximation
as the force and collection models. The electric actuator supplies
`I alpha + omega cross (I omega) - gravity torque - radiation torque`.
For the stated three opposite edge force couples, its force account is
`2 sum(abs(torque components)) / side`.

This permits signed environmental torque to assist or oppose a command.
The earlier nominal-plus-envelope sum allowed the most adverse torque at
every time. Keep that conservative engineering allowance visible while
reporting the computed ideal-model history separately.

The already verified service passage has **0.03794 m/s** calculated attitude
impulse. Its isolated gravity-couple account is 0.03126 m/s and the isolated
radiation-couple account 0.01032 m/s; they do not add directly because the
actuator supplies their vector combination. The previous service allowance
was 1.66184 m/s. Nominal slow Sun tracking remains negligible.

On the old verified departure, calculated attitude control costs **2.79604 m/s**
including spin-up, braking and the signed disturbance history. Adding its
1.40584 m/s translation and the preceding service gives **4.23982 m/s** for
the same twelve-hour geometry. The previous nominal-plus-envelope account
was **6.70269 m/s**. The revised model account is about **36.7% lower**; the
maneuver itself has not improved. Its nominal turn still costs 2.72711 m/s.

The new trial's ledger ends at its clearance stop. Its shorter cost integral
cannot be compared as a saving over the complete old departure. The product
therefore leaves the energy-reduction and recurring-power entries unfilled.
Full return, handover, repacking, recovery, cadence and hardware-mass feedback
remain open. The roughly 1.30 m/s total allowance for a hypothetical two-day
cycle at the optimistic 17.29-million inventory and 10% power target remains
far below the verified prefix's revised 4.24 m/s. That allowance uses the
previous 20% fixed added mass and 10% thrust-duty scenario; the inventory is
a capacity relaxation and the two-day recurrence is an assumption.

The useful next trajectory search must change the inherited service/exit
packing or provide more independent control directions and times. Jointly
choose those with the return and handover, retaining energy as the objective
and 100 m as the clearance constraint. An alternative attitude actuator also
needs its own complete mass, momentum-storage, unloading and energy account.
The current failures do not establish a fleet-wide power floor.

For context, the same scalar inventory and hardware calculation gives these
complete-cycle allowances. The half-benchmark column keeps the intermediate
cost target visible. These are permitted average equivalent impulses, not
achieved operations or sufficient populations.

| Conditional inventory | At 23.835 TW, m/s per day | At 119.175 TW, m/s per day |
|---|---:|---:|
| 17.286 million | 0.6489 | 3.1232 |
| 25 million | 0.4500 | 2.1910 |
| 35 million | 0.3220 | 1.5798 |

The table solves `conditional_cost` with the propulsion settings in
`pilot_budget.json`, 20% fixed added mass and 10% thrust duty. A two-day cycle
gets twice the listed impulse allowance. Inventory can ease coverage or
handover scheduling while reducing the impulse available per member at a
fixed fleet power. Nothing here establishes which side of that trade wins.

## Torque refinement and optical inventory

The verified twelve-hour prefix costs **2.31917e14 J** in this rigid electric-
couple account. Its sampled peak propulsion demand is **64.01 GW** for the
361-tile patch. At 300 W/kg and 25% peak allowance, power hardware alone is
**2.6671e8 kg**, or **14.78%** of the patch's base optical mass. That increment
has not been propagated; storage, distribution, actuator structure and
flexibility remain outside this account.

Disturbance torques are evaluated at 60 s and the commanded rigid-body load
at 2 s. Coarsening the disturbance interpolation to 120 s changes integrated
attitude impulse by at most **0.000451%** across the three accounted traces.
Changing from 16 solar directions to 32 rotated directions, both at 120 s,
changes it by at most **0.00564%**. Three- versus five-point-per-axis gravity
quadrature differs by at most **2.05e-14 Nm/kg** at the checked hourly dates.
The torque helper and existing collection helper agree on total first
interception within **1.17e-17** relative error, and sampled torques remain
inside the earlier conservative envelopes. These checks refine the ideal
load calculation; they do not validate flexible hardware or certify the
failed trajectories.

Every rejected coupled trial retains a separate per-member optical ledger
through its own stop. The following means cover departure only. Each row is
an alternative experiment with an unequal duration; their powers cannot be
added into one fleet or compared as complete-cycle yields.

| Trial | Recorded departure, h | Mean first-intercept bolometric equivalent, TW | Mean redirected-band optical power, TW | Redirected-band energy, J |
|---|---:|---:|---:|---:|
| 24 | 5.4866 | 26.7386 | 3.6886 | 7.28562e16 |
| 20 | 3.5708 | 29.2473 | 4.0347 | 5.18649e16 |
| 41 | 4.4679 | 22.8561 | 3.1530 | 5.07143e16 |

These 16-source, 60 s accounts include mutual shadows and both faces.
Coarsening to 120 s changes each first-intercept energy by less than
0.000098%. Trial 41 also has the 32-source, 120 s load-account refinement;
its combined source/time change in interception is 0.00539%. Trials 24 and
20 have no independent source refinement. The older complete service and
departure optical comparison remains in `collection.md`.

The bolometric equivalent measures distinct intercepted rays. Only the
existing reflected-band fraction is credited as redirected optical power;
no additional absorbing collector or electrical efficiency is assumed.
Collection capacity, propulsion electricity and the **100 TW delivered-power
ambition** stay separate. A different collection design must feed its mass
and momentum change back into the trajectories.

## Resources and verification

The six completed producers record 1,118.69 s of numerical wall time. Including
a 4.68 s initial correction attempt and 193.18 s of accounting before fixing
the truncated-trace command construction, recorded numerical work totals
**1,316.55 s (21.94 minutes)** against the 2,700 s budget. Peak recorded RSS is
**843.43 MiB** and retained raw output **85.42 MiB**. Calculations ran one at a
time with one numerical thread and a 2 GiB address-space limit.

The accounting retry preserves each failed candidate's original commanded
turn duration and integrates only its executed prefix. It does not shorten
the intended turn to fit the stop. No tighter coupled replay, full return,
two-pattern handover, new fleet placement or hardware-mass propagation was
executed. Clearance stops are local rejection evidence, not continuous
coverage certification or a general impossibility result.

The final dynamics/study suite reports **120 passed in 7.34 s**. Its signed-
control check runs in a fresh process, matching the study runners: HiGHS
retains its first process-wide thread setting, which otherwise conflicts
with earlier tests using the default thread configuration. The executed
correction product reports infeasibility status 2, not that scheduler error.

The required `make check` passes the staged layer check (496 files, zero
violations), then reports **665 Python passes, 55 skips and 13 existing
failures in 28.00 s**. The failure set is unchanged: eight absorption/exosphere
input cases, CM1 terrain configuration, process discovery, two ring-comfort
cases and the protection-architecture design point. Later JavaScript,
provenance and ensemble targets are not reached. Source/parent-product hashes
and constants are checked by the passing study tests. No earlier pinned
model or numerical product was modified.

## Products and reproduction

`energy_screen.json`, `energy.json`, `energy_correction.json`,
`energy_replay.json`, `energy_budget.json` and `energy_collection.json` retain
the proposals, three clearance stops, failed signed correction and
torque/collection accounts.
The old source-bound products remain unchanged. Raw states, loads and per-tile
optical power remain under the ignored `research/runs/solar_shield_array/energy/`.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_screen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_correct
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_replay
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.energy_collection
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest protection/dynamics research/studies/solar_shield_array -q
make check
```
