# Encounter-directed return investigation

Three coupled corrections fail to complete the return. The strongest local
correction removes the original encounter and postpones the 100 m stop by
about 90 minutes, but develops a deeper overlap and exceeds installed power
on 160 tiles. Its tighter replay agrees within 10.58 m. The reviewed
twelve-hour passage remains intact. The energy comparison now distinguishes
the lunar-disk dimming estimate from the full held aperture, and keeps
23.835 TW alongside explicit 50/100 TW electrical and material sensitivities.
No new operating target or repeatable cycle is accepted.

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

## Decisions during the bounded search

The eleven timing/roll proposals took 439.54 s. The twelfth proposal's first
local LP used 0.54248 m/s mean translation, but adding later encounter cuts
produced 107,417 inequalities and reached its 15-second solver limit. This
is an unresolved solver result. Its 29.59-second attempt is retained. Continue
that same local proposal with at most four refinements, retaining the worst
currently violated date per pair/direction before adding cuts. Charge the
continuation to the remaining 130 seconds of the unchanged 600-second proposal
budget. Preserve both attempts and all tested failed branches.

## Executed proposal comparison

The first eleven proposals begin at the tighter hardware-loaded hour-12
terminal states. Each integrates independent finite-source forces through the
original hour-46.27978 arrival date. The first six return hours receive
20-second finite-square geometry and between-date coplanarity searches; the
whole proposed return receives 120-second surface checks. Mutual shadows are
omitted at this stage, so these are proposals. The subsequent coupled trials
recompute shadows. The unchanged proposal has 1,620 early near-crossing events;
removing its first burn raises that count to 3,555.

| Changed first-return command | Early mean translation, m/s | Nominal early attitude, m/s | Maximum endpoint position debt relative to unchanged proposal |
|---|---:|---:|---:|
| Omit the first burn | 0 | about 0.000001 | 316 km |
| Half the first burn | 0.48843 | about 0.000001 | 158 km |
| Three quarters of the first burn | 0.73264 | about 0.000001 | 78.7 km |
| Delay the first burn by 0.5 / 1 / 2 / 4 hours | 0.97686 | about 0.000001 | 50.8 / 101 / 199 / 375 km |
| Smooth +60-degree roll over 1.5 hours | 0.97686 | 2.42407 | Numerically unchanged |
| Smooth +75-degree roll over 1.5 hours | 0.97686 | 3.03009 | Numerically unchanged |
| Smooth -30-degree roll over 1.5 hours | 0.97686 | 1.21204 | Numerically unchanged |

Every proposal has subsequent square crossings. The roll proposals use a
quintic acceleration/deceleration history, followed by a hold and a proposed
four-hour recovery starting at hour 18. In-plane roll preserves the normal
and nearly preserves the independent model's force; its changed square edges
alter physical encounters and mutual shadows. Recovery is included in the
proposal's endpoint debt, but no coupled trial reaches that recovery. The
nominal attitude column excludes disturbance torque and is an actuator
account, not net spacecraft velocity change.

The twelfth proposal removes the endpoint constraint and permits two independent
early acceleration vectors on every tile, over hours 12–14 and 14–16. For each
approaching pair it initially preserves the strongest separating direction at
hour 12. These are local choices of collision-avoidance direction. The first
LP gives 0.54248 m/s, but creates other close approaches. The dense refinement
reaches its solver time limit. The compact continuation solves four LPs with
1,983, 5,756, 8,817 and 10,846 retained inequalities; their mean impulses rise
through 0.63402, 1.17518, 1.42416 and 1.60093 m/s. Each finds another conflict.
All records are retained. A failed branch or an iteration limit does not
supply a global lower bound.

Independent nonlinear propagation of the final local iterate still predicts
1.60093 m/s through hour 18 and a 566.06 km / 19.54 m/s maximum endpoint debt.
The later original burns are retained only to make that comparison explicit.
They have not been refitted to a different service assignment. No complete
arrival, next service, recovery or handover receives credit.

## Coupled trials and stopping states

Three corrections are propagated with the original per-member hardware mass,
the same optical reflection law, finite-area gravity, finite-Sun mutual
shadows and the 100 m stopping constraint. All 361 tiles remain in the force
and collection inventory after release from service. The unchanged six-hour
service and six-hour departure are reused from their verified trajectories;
there is no change before hour 12 and no additional assigned-service ray is
claimed during return.

The smooth 60-degree roll stops at hour **12.399292**, only **23 min 57 s**
after return begins. Pair 338/341 reaches 100 m during the transition. Its
body-frame displacement is (-10030.551, 7660.575, 95.219) m; the actual
Euclidean surface distance includes both the edge gap and normal gap. The
frozen final-angle clearance advantage cannot be reached along this tested
transition. Its later attitude recovery remains unexecuted.

The local-burn correction reaches hour **15.390970** in the coarse coupled
run, **89.85 minutes later** than the old stopping event. Pair 151/186 becomes
limiting at body displacement (7814.714, 652.301, 100.000) m. Mean translation
already applied is **1.57187 m/s**, with 518 started burns. Maximum commanded
acceleration is 0.0008165 m/s2 and the reflected-beam margin remains positive.
The independently tighter replay and the numerical-margin audit below decide
how accurately that obstruction is localized. Both trials stop at clearance;
neither is propagated into a physical impact.

## How to compare collection and operating allowances

The source reconstruction uses `scenario.json`, `design.md`, the held aperture
in `results/holding.json`, and the actual mass-loaded passage in
`results/closure_budget.json`. It places the lunar disk, the held aperture and
clear moving film on the same 1361 W/m2, face-on, unshadowed normalization.
Actual passage integrals use their simultaneous ephemeris dates. This avoids
mixing an annual holding demand, a twelve-hour collection fraction and a
hypothetical fleet cycle as though they were the same measurement.

The original extra dimming is 1235 W/m2 times 0.05 across the lunar disk.
The complete ideal redirected band is 1361 - 1235 times 0.95 = 187.75 W/m2,
or 13.7950037% of incident sunlight. Multiplying those quantities by the
larger held aperture gives a different optical account from the original
lunar-disk figure. The historical holding dynamics also allowed full optical
authority outside the climate window. Its potential collection therefore
needs an explicit outer-annulus allocation. Capturing the same redirected
band everywhere and capturing all permissible annular light are separate
sensitivities; neither establishes a collector consistent with the ideal
holding forces.

For the unchanged moving optical mode, write the hypothetical fleet mean as

`redirected optical power = tile count × clear area × solar irradiance × 0.137950037 × f`.

Here `f` is the full-cycle effective interception fraction: duty, orientation,
body occultations, mutual shadows and overlap between patterns all enter it.
Every tile contributes at every date, including tiles without a current
shielding assignment. The measured twelve-hour fraction fixes only that
passage. A two-day thought experiment with the remaining time given zero
collection divides it by four, before cross-pattern shading. It supplies a
threshold calculation and no operational duty cycle. A completed cycle does
not exist in the current evidence, so full-cycle energies and averages stay
null in the result products.

The electrical sensitivities capture only the existing redirected band with
external collectors. Captured optical power times the assumed conversion
efficiency gives gross bus electricity. Propulsion and explicitly assumed
other operations are subtracted at that bus, then the remaining power passes
through the assumed delivery link. Collector waste heat, missed optical
power and delivery loss are separate entries. No transmitted-band electricity
is counted. Device spectra, angular response, thermal design, collector and
routing mass, receiver geometry and radiation momentum remain unspecified.
Changing the shield's absorption, optical mode or loaded mass would require
another dynamics replay before accepting that design.

The third and final allowed coupled trial uses the already screened -30-degree
turn. A square has 90-degree rotational symmetry, so its final physical
orientation equals the +60-degree case while its nominal roll impulse is
halved. Its transition is different and requires its own propagation. This
uses the original three-trial limit and adds no new proposal family. The
single tighter replay remains assigned to the local-burn correction, which
has extended the original clearance-limited interval.

## Reconstructed optical resources and duty thresholds

The reconstructed source accounts are optical powers at 1 AU. They include no
conversion efficiency or capture credit.

| Area and allocation | Face-on first-intercept sunlight | Extra 5% dimming | Combined filter/dimmer band |
|---|---:|---:|---:|
| Lunar disk | 12.906 PW | **585.580 TW** | 1.78045 PW |
| Fixed 78,000 km held aperture | 229.320 PW | 10.4045 PW | 31.6347 PW |
| Selected held aperture, radius 7,454 km | 237.568 PW | 10.7787 PW | **32.7725 PW** |
| Same clear area as 361 physical tiles | 48.0571 TW | 2.18040 TW | 6.62948 TW |

For the actual six-hour service and six-hour departure, retain the original
16-source first-intercept means **39.571 / 17.057 TW** and redirected-band
means **5.459 / 2.353 TW**. The rotated 32-source check gives
39.5703 / 17.0551 TW and 5.45872 / 2.35275 TW. The reconstructed comparison
uses that latter check: 1.223108e18 J first interception and 1.687278e17 J
redirected optical energy over twelve hours, averaging 28.3127 and 3.90574 TW.
The same clear film continuously face-on and unshadowed at those actual dates
would intercept 2.076269e18 J. Thus the measured fraction is **58.9089%**.
The small change from the rounded prior figures comes from source/time
quadrature; the original ledger remains intact.

At the capacity-relaxation inventory of 17,285,560 tiles, the unshadowed
face-on normalization is 2.30109 EW first interception and 317.436 PW in the
unchanged redirected band. The unknown complete-fleet `f` must exceed
**10.3241%** to match the selected held aperture's same-band optical power.
A hypothetical twelve-hour passage every two days, with zero collection
elsewhere, supplies `f = 14.7272%` before cross-pattern overlap. It would have
to retain **70.1023%** of that already locally shadowed interception after
placing all patterns together to match the same-band held resource. Neither
the two-day cycle nor that overlap retention has been demonstrated.

The held model's broader annular allowance is approximately **218.982 PW**:
all incident light outside the enclosing climate window plus the redirected
band inside it. Matching that resource while keeping the moving array's
present uniform 13.795% optical mode requires **f = 68.9848%**. This is a
material distinction between source allocations. Both systems could use
different collection modes, but those modes need optical-momentum and mass
feedback. The old 586 TW figure supplies no basis for claiming that the
larger moving inventory beats the full held aperture.

## Electrical alternatives and recurring material

The sensitivity cases assume the following independent efficiencies. They
are illustrative inputs, with no selected spectral device or routing design.

| Case | Capture of redirected beams | Optical-to-electric conversion | Final delivery link |
|---|---:|---:|---:|
| Low | 25% | 20% | 70% |
| Middle | 50% | 40% | 85% |
| High | 80% | 60% | 95% |

The table below assigns the stated operating allowance entirely to propulsion,
uses the middle efficiency case, and sets other electric loads to zero.
Its final column gives the hypothetical effective fleet fraction sufficient
for **100 TW delivered after propulsion**. The result product also carries
0 and 10 TW other-load sensitivities, collector waste heat, missed beams,
link losses and explicit power deficits.

| Propulsion allowance | Share of held demand | Exhaust, kg/s | Redirected optical power required | Required fleet f |
|---|---:|---:|---:|---:|
| **23.835 TW, original target** | 10.0% | 37,077 | 707.410 TW | 0.22285% |
| **50 TW, alternative** | 21.0% | 77,778 | 838.235 TW | 0.26406% |
| **100 TW, alternative** | 42.0% | 155,556 | 1,088.235 TW | 0.34282% |
| 238.35 TW comparison | about 100% | 370,767 | 1,779.985 TW | 0.56074% |

At 50 TW, the low/high efficiency cases require respectively 1.2151% and
0.10190% effective interception to deliver 100 TW. At 100 TW they require
1.5301% and 0.13471%. These small fractions reflect the large assumed film
inventory; they do not establish its collection geometry or energy supply.

Relative to 23.835 TW, a 50 TW allowance needs **130.825 TW more redirected
optical power** under the middle capture/conversion case to preserve the same
net delivery. That is an additional **0.04121 percentage points** of effective
interception at the assumed inventory. Moving to 100 TW needs 380.825 TW more
redirected light, or **0.11997 percentage points**. The calculation supports
considering larger electrical allowances if actual collection is established.
It also increases exhaust by 40,701 or 118,479 kg/s respectively.

For comparison with the held array's own net output, the same middle
bookkeeping puts the 50 TW alternative's break-even `f` at **10.0275%** under
uniform-band static collection, or **68.6881%** under the broad static optical
allowance. These subtract the 238.347 TW held propulsion reference before
comparing net electricity. The compatibility of harvested beams with that
ideal holding controller remains an engineering input.

The larger inventory's propagated mass allocation scales to **1.08372e14 kg**,
including **8.64278e13 kg** of optical mass and **2.19444e13 kg** additional
modeled hardware. The held reference has 1.05628e13 kg of dry hardware after
removing its seven-day propellant buffer. A complete collector, storage and
routing system is absent from both mass comparisons.

| Assumed whole-hardware life | Moving replacement, kg/s | Held dry-hardware replacement, kg/s |
|---|---:|---:|
| 10 years | 343,411 | 33,471 |
| 30 years | 114,470 | 11,157 |
| 100 years | 34,341 | 3,347 |

These sensitivities assume complete replacement without recycling and identical
service lives. At 30 years, 50 TW propulsion plus moving-hardware replacement
uses about **192,248 kg/s**, compared with **381,920 kg/s** for the held case.
At 100 TW the moving total is **270,026 kg/s**. At a ten-year life, both those
moving alternatives exceed the held material flow. The equal-lifetime material
break-even is about **10.58 years** for 50 TW or **14.40 years** for 100 TW.
Actual component lives, recycling, manufacturing losses, collector mass and
sufficient inventory remain inputs to that comparison.

**The operating target remains 23.835 TW.** The 50 and 100 TW alternatives
have explicit electrical and material advantages and costs, but the present
results establish neither a complete safe return nor usable fleet generation.
They therefore supply comparison cases rather than an adopted replacement.

The first accounting producer reached its 540-second stage alarm after saving
complete local-correction and +60-degree accounts. Its unfinished third trial
is charged in full. Reallocate unused search/propagation time to permit a
660-second accounting continuation plus the existing 60-second audit, keeping
the **3600-second overall ceiling unchanged**. Archive the original producer
and partial product with hashes. Reuse only completed trial accounts whose
parent results, physics sources and raw data still match; recompute the
unfinished trial. This raises the accounting allocation within the overall
budget without adding trajectories or changing the physical model.

## Final clearance and hardware checkpoint

The -30-degree turn completes its attitude transition and reaches hour
**14.155954**, **15.75 minutes beyond the original stop**. Pair 305/343 then
reaches 100 m. Its in-plane centre separation is only **2.623 km**, and normal
spacing closes at about **0.705 m/s**. The +60-degree trial had stopped after
turning only **7.280 degrees**, demonstrating the significance of the complete
transition path despite identical final square geometry.

The local-burn correction's tighter replay stops at **hour 15.390986** on the
same pair 151/186. Event time differs by **0.05950 s**; maximum position
comparison, including reconstruction, is **10.5791 m** and maximum velocity
comparison is **0.0003573 m/s**. The 50 m comparison gate passes. At the old
failure time, the new global surface minimum is **468.07 m**, while original
pair 325/347 is **999.91 m** apart. This confirms removal of the original
obstruction in the propagated ideal model.

The uncertainty margin matters before the final stop. Subtracting both
centres' 10.5791 m replay allowance from the all-pair interval guard first
loses the 100 m requirement during **hours 14.958889–14.959444**. Thus a
conditional geometrical margin exists through about **2 h 57 min after
return starts**; the entire 3 h 23 min trajectory to the stopping event
cannot be labeled clearance-safe with that allowance. At the final stop,
the conservative coarse/fine bound after both allowances is only **39.71 m**.
The guard assumes acceleration at most 0.15 m/s2 and rigid attitude rate at
most 0.1 degree/s. The replay difference is measured numerical evidence,
rather than a proven error envelope for the physical structure.

The new local-burn witness has **7.840 km** in-plane centre separation and
normal closing speed **0.660 m/s**. Any 10 km square contains a 5 km-radius
inscribed disk. When two parallel squares' projected centres are less than
10 km apart, those disks overlap for every in-plane orientation, including
independent rolls. Their surface distance is then their normal spacing.
This geometrical obstruction applies to both migrated witnesses. Further
roll at those fixed centres cannot restore clearance. Useful next freedoms
are earlier changes in relative position/velocity, the initial velocity
field and packing, independent timing, or genuinely different plane normals.
The finite-Sun coverage and beam constraints must remain coupled to such a
change. Four local LP refinements do not exhaust those possibilities.

| Executed correction | Return time until 100 m floor | Return energy, TJ | Entire executed prefix, TJ | Ideal exhaust over prefix, kg | Hypothetical two-day scaling, TW |
|---|---:|---:|---:|---:|---:|
| Original reviewed return | 1.8935 h | 67.479 | 119.338 | 185,637 | 33.07 |
| +60-degree transition | 0.3993 h | 54.020 | 105.878 | 164,700 | 29.34 |
| -30-degree transition | 2.1560 h | 150.678 | 202.536 | 315,056 | 56.12 |
| Local-burn correction, tighter replay | 3.3910 h | 107.812 | **159.671** | **248,377** | **44.24** |

The three new rows use the rotated 32-source prefix account plus the
16-source/120-second return load integration; this changes the inherited
prefix cost by less than 0.001 TJ. Different stopping times prevent a
like-for-like saving claim from the incomplete totals. Scaling assumes the
17.29-million capacity-relaxation inventory and repeats just the measured
expenditure every two days. No completed return, recovery, acquisition,
handover or operational cycle is charged by that thought experiment.

All three corrections exceed individual installed power allocations:

| Correction | Tiles over installed capacity | Lowest installed/required peak ratio | Added power-equipment mass needed to keep the existing allocations and restore 25% margin |
|---|---:|---:|---:|
| +60-degree transition | 361 | 0.3229 | 120.48 million kg |
| -30-degree transition | 361 | 0.4807 | 74.77 million kg |
| Local-burn correction | 160 | 0.4037 | 34.23 million kg |

Those equipment additions are estimates on the existing trajectories. They
have **not** been installed or propagated. Redistribution of equipment could
reduce additive mass, but changes individual mass and sail acceleration too.
Even the local correction's sum of bare individual peaks, 27.09 GW, does not
ensure each tile fits the existing 29.19 GW aggregate allocation. These
physical power-allocation failures prevent accepting the ideal forced
trajectories as realized maneuvers. Added mass and changed momenta must be
replayed from their first change before accepting a revised vehicle.

## Executed collection and conditional electricity

All three correction traces retain every member's two-face optical histories,
including the entire feathered interval. Mean powers during the newly
executed return portions are:

| Correction | First interception, TW | Redirected optical band, TW |
|---|---:|---:|
| +60-degree transition | 0.093785 | 0.012938 |
| -30-degree transition | 0.094285 | 0.013007 |
| Local-burn correction | 0.094434 | 0.013027 |

These use 128 rotated solar sources. The 64/128-source integral differences
are **1.02–1.035%** in these grazing stages. This remains an uncertainty in
the small off-service collection term; the comparison is not a certified
error bound. Earlier service/departure source quadrature remains as recorded
above. Torque/propulsion energy is stable to better than 5.3e-9 in the
120/240-second coarsening check and 1.8e-10 in the matched-time source check
for these executed return loads. Nominal attitude and burn profiles are
integrated at two-second spacing.

For the local correction's entire executed 15.391-hour prefix, first
interception is **1.224261e18 J**, averaging **22.0956 TW**. Redirected optical
energy is **1.688868e17 J**, averaging **3.04808 TW**. Electric propulsion
consumes 159.671 TJ, averaging **2.88175 GW for the 361 tiles** over that
actual interval. These are partial-trajectory integrals and average powers;
full-cycle quantities remain unfilled.

The middle capture/conversion/delivery sensitivity gives **28,575 TJ** of
net delivered-energy potential over that same interval after its modeled
propulsion, averaging **515.7 GW**, with other electrical loads set to zero.
That assumes temporal balancing and successful external collectors. During
the return portion alone, generation in that sensitivity averages about
2.61 GW against about 8.83 GW propulsion, leaving a **76.0 TJ net energy
shortfall** to supply from earlier collection or transfer. Storage mass,
power routes, instantaneous capacity, thermal control and collector momentum
are unresolved. The low/high and both roll sensitivities retain the same
separate energy and heat accounts in `results/repair_synthesis.json`.

## Acceptance, reproduction and checks

**Preserved:** the reviewed deeper starting arrangement, original hardware
allocation, twelve-hour zero-translation passage and assigned finite-Sun
service evidence. **Verified obstruction:** the local correction clears the
old encounter but reaches a new deep overlap, with a tighter replay inside
the position gate. **Failed:** full return clearance, per-member installed
power for every correction, and every tested six-hour escape proposal.
**Unpaid:** revised hardware feedback, acquisition, remaining return, attitude
recovery, next service, handover, recurring duty, global placement and
continuous coverage, collection hardware, routing/storage and replacement
lifetimes. No handover was attempted because its prerequisites failed.
No new operating target, complete cycle or covering fleet is accepted.

The next bounded search should use the measured population and migrated
encounters to change approach velocities earlier, potentially during the
already released departure interval, while preserving the initial deep
packing as a comparator. Coupled independent normals require a more general
shadow/contact model. A changed return service assignment remains available;
its geometry and later consequences must be propagated. Keep the current
power allocation as an explicit per-member constraint or propagate a revised
allocation from the initial epoch. Compare 23.835, 50 and 100 TW alongside
propellant and replacement rates before selecting an operational allowance.

The completed producers, initial inspection charge and full interrupted
accounting charge total **about 2725 s (45.4 min)**, within the 3600 s ceiling.
Maximum recorded RSS is **567.13 MiB** and raw output is **212.89 MiB**, below
the 2 GiB address-space and 512 MiB output limits. Numerical producers run
sequentially with one thread. The accounting continuation reuses two completed
trial accounts only after matching parent hashes, all unchanged physics
sources, the archived original producer and each raw-data hash. Its final
trial is computed fully. The failed 540-second attempt remains charged.

Reproduction requires the existing raw products and ephemeris described in
`closure.md`. Run sequentially from the repository root with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`:

```sh
python -m research.studies.solar_shield_array.repair_screen
python -m research.studies.solar_shield_array.repair_local_screen
python -m research.studies.solar_shield_array.repair_refine
python -m research.studies.solar_shield_array.repair_run --case local_refined
python -m research.studies.solar_shield_array.repair_run --case roll_60
python -m research.studies.solar_shield_array.repair_run --case local_refined --replay
python -m research.studies.solar_shield_array.repair_run --case roll_minus30
python -m research.studies.solar_shield_array.repair_trade
python -m research.studies.solar_shield_array.repair_budget
python -m research.studies.solar_shield_array.repair_audit
python -m research.studies.solar_shield_array.repair_synthesis
python -m pytest protection/dynamics research/studies/solar_shield_array engineering/test_shield_power.py -q
make check
```

A fresh accounting run computes all trials with a 1200-second stage limit.
`--resume-manifest` is available for the source-bound interrupted account and
uses a 660-second continuation limit. Fresh computations need no cache.
The synthesis records the actual first-attempt charge when a cache was used.
The compact result products bind sources, constants, parents and raw files.
The targeted suite passes **150 tests in 11.70 s**, including conserved
optical/electrical partitions, partial-burn integration, smooth roll
boundaries, compact/full optical agreement, preserved failed branches and
all new product identities.

The final `make check` reports **695 passes, 55 skips and the same 13 baseline failures in 29.66 s**. Missing pinned WHI spectral data, generated climate configuration and the environment's process lookup account for those failures. Later Makefile targets are not reached. A four-test source-identity recheck passes after fixing fresh-run resource-summary handling; no dynamics or load result changed. After staging all new files, the layer check covers **539 files with zero violations**, and the staged whitespace check passes. `checks.json` records the gates and remaining work.
