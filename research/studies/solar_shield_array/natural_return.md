# Service-led natural-return search

**Result: a cheaper loaded twelve-hour passage, but no accepted return or
second service in this search.** Starting the 90-degree departure turn at
hour 6 and spreading it over six hours reduces measured electrical control
energy to **36.186 TJ**, from **51.859 TJ**: **30.22% less**. All 361 physical
tiles retain the initial six-hour optical service, zero electric translation
and about 1.72 km minimum sampled clearance through hour 12. The subsequent
free return reaches the 100 m floor at **hour 13.799825**. Longer natural
return dates and freer arrival geometry did not resolve that early encounter.

Continue `fc803222d34b29a242af471a5e36b4c1269586a3` on
`research/solar-shield-habitat-array`. Main remains unchanged. The primary
magnetic shield is outside this optical photogravitational/electric study.
The author authorized this bounded follow-up after reviewing the expensive
prescribed corridor. Its 1.53538 PJ over 51.77978 hours averages 8.23669 GW
for 361 tiles; the 238.35 TW held-array figure describes a different scope.
About 47% of that sequence's energy held relative positions for 23.8 hours.

## Bound and decision rules

Allocate 1,800 CPU seconds, one numerical process/thread, 2 GiB address space
and 512 MiB new raw output. Required checks have a separate five-minute
allowance per commit. Record unsuccessful proposals and partial executions.
Do not expand into a fleet campaign within this allocation.

Preserve all 361 physical 10 km squares, their 9.89 km clear apertures, loaded
mass and installed per-member ratings; the first six-hour service; the 100 m
surface floor plus numerical allowance; finite-Sun six-hour service of a
moving 20 km receiver inside four lunar radii; and the inherited spectral
rejection. Retain 238.35 TW held comparison, 23.835 TW target and the
17.29-million-tile capacity relaxation in their original evidence states.

Screen three families:

1. Freely deforming trajectories with smoothly commanded common attitudes,
   finding service opportunities over ten days without an arrival lattice.
2. Small early velocity changes, selected for separation and prospective
   service, with endpoint positions and velocities free. Pay for all changes
   from the actual loaded state and enforce individual actuator ratings.
3. Modest grouped off-service sail steering, using validated oriented-square
   geometry; any useful candidate requires coupled photon-force replay.

Use the loaded **51.859 TJ** first-twelve-hour account as comparator and a
**100 TJ** screening ceiling through the next six-hour passage. The latter
is an explicit search cutoff, not a derived fleet target. A reduced model may
only reject or propose. Acceptance needs actual all-member dynamics, adaptive
finite-square clearance, the executed service passage, optical/energy accounts,
individual power and a tighter replay. Added supply hardware must be placed,
charged and propagated before operational acceptance.

For a candidate with one six-hour passage every T hours, report the necessary
time-capacity inventory multiplier T/6, its corresponding tile/mass allocation,
and energy per credited service hour. These are conditional lower bounds;
spatial coverage, phases, replacement assignment and handover can require
more. Use energy/inventory Pareto comparison, with no presumed exchange rate
between electricity and hardware. A long coast never grants free inventory.

Start from the exact loaded epoch and service-terminal arrays in
`results/joint_seed.json`. No forced arrival reset, fixed two-day return,
prescribed held formation, or unexecuted attitude change earns acceptance.

## Baselines and the actual constraint

The source-bound seed and [preceding regression](joint_cycle.md) remain the
comparator. [Cycling](cycling.md) establishes separate point-tile cycling and
an independently replayed 8.905-day point-tile return. [Packing](packing.md)
establishes the optical-mass 82% local saving. [Closure](closure.md)
establishes the hardware-loaded twelve-hour passage: 361 squares, four depth
levels spanning 36 km, approximately 25.4% added mass, 1,719.63 m sampled
clearance and 51.859 TJ of attitude energy despite zero electric translation.
Its old return fails at hour 13.893535. [Return repair](return_repair.md)
moves the failure to hour 15.390986 and overloads 160 members. None of these
proves a finite-formation operating cycle. The loaded model here has
2.263298487 billion kg total mass and 29.189546 GW of installed actuator
power, with separate per-member allocations.

The new slower-turn failure is formation compression, not an insufficient
arrival-date search. At hour 13.799825, pair 325/347 has relative displacement
(-6,317.135, 8,952.981, 100.000) m in the common tile frame and relative
velocity (0.61148, -0.51907, -0.36875) m/s. Both in-plane projections overlap
the physical squares; normal separation is closing through the floor.
An affine fit from the initial positions has singular values
(2.6831, 0.3072, 0.06821), indicating strong compression along one direction.
Its 1.529 km RMS residual also shows that a single uniform deformation is
not the whole motion. These measurements do not prove that control or
independent roll cannot avoid the encounter. They identify where earlier
relative-velocity control must act.

## What was searched and what each failure means

`natural_fast_scout.py` searched ten days using centre gravity, a point-source
photon force and body visibility, omitting mutual shadows. All 361 positions
were retained. The initial finite-area, eight-source scout was stopped after
one complete family during the second because it consumed too much of the
allocation; its partial result and interruption record are preserved.

| Family or refinement | Executed search | Result and restriction tested |
|---|---|---|
| Free deformation | Sun-facing coast; 90-degree turns about either edge; four natural opportunities per family | First reduced-model clearance failures at 10.023, 12.748 and 10.023 h. Natural service dates occur near 46, 92, 139 and 185 h. Hypothetical Sun-facing coverage of the centred receiver at the first opportunity is only 6.6–9.0%; later centred probes are zero. This does not exclude every shifted receiver. |
| Two control arcs | Six first/third-return cases; four additional cases with control from hour 0 or 6 and a slower turn | Greedy ray ownership and the selected arc timings produce no feasible capped LP. Terminal positions and velocities are free. A coarse 51.948 TJ attitude estimate exceeds the 51.859 TJ prefix cutoff by only 0.17%; this is a screen result, not a physical exclusion. The slower turn resolves that cost issue. |
| Assignment changes | Four 8.5/9 km static ownership seeds with Hungarian assignment | Capped LPs fail, including energy-relaxed diagnostics retaining individual power limits. These failures constrain those ownership/timing branches, not all formations or an unrestricted minimum energy. |
| Early and arrival control | Four arcs at 0–6 h, 6–12 h, arrival minus 8–4 h and arrival minus 4–0 h | Static ownership seeds remain infeasible. Actual-deformation ownership makes the Sun-facing branch feasible in the sparse ray LP; selected encounter cuts reject the edge-facing first/third-return branches. |
| Full receiver holes and retained encounters | One refinement of the cheap Sun-facing branch | 11,739 violating pair/date samples and worst next-service source coverage 21.28%. Retaining 1,200 chronological collision cuts and 7,436 ray inequalities makes this local LP infeasible. |
| Four row-band attitudes | Smooth independent offsets spanning ±5 or ±15 degrees around the slower departure turn | Exact oriented-square checks on reduced trajectories find failures at 12.764 and 12.644 h. Three pre-contact oriented shadow/torque snapshots per case are retained. Neither receives coupled-trajectory acceptance. |
| Coupled physical executions | Cheap Sun-facing proposal; zero-translation slower turn, then tighter replay | Cheap proposal fails initial optical service and stops at 8.925 h. Slower turn passes the twelve-hour local prefix but stops at 13.800 h. No next service is executed. |

The first/third dates differ by about four days. The LPs constrain sampled
receiver-ray interception, not an arrival grid or prescribed terminal
velocity. They retain particular ray ownership and a small set of burn
windows. There is no new independent per-tile arrival-time search here;
the earlier staggered-slot failures remain in `joint_cycle.md`. Arrival
dates, receiver paths and assignments are proposals, not an exhaustive
continuous optimum. Every changed velocity is charged from the actual
loaded initial state; early-service changes must preserve that service too.

The dynamic-ray product has four rows but **three distinct variants**:
`face_1_8500` and `face_1_9000` use the same transported, actually deforming
8.5 km first-service seed. The wrapper ignores the static pitch label and
produces identical control arrays. Both attempts remain recorded and charged;
the second is not evidence for a separate 9 km dynamic packing.

The inexpensive fit reserved estimated attitude demand before imposing each
member's remaining power cap, weighted control-energy ceilings and the
twelve-hour comparator. Four-arc prefix charging integrates partial smooth
burns exactly, including burns during the first service. These caps cover
the reduced model only; coupled executions independently stop at an installed
rating violation. Generation and bus availability remain additional gates.

## Why the cheap proposal did not survive

The last sparse Sun-facing fit estimated **36.528 TJ through another service**,
including 8.753 TJ translation. It still had 7,876 violating pair/date samples
and only 75.34% minimum sampled next-service coverage. Expanding from sparse
receiver probes to whole receiver polygon differences at more Sun angles
exposed the much larger 21.28% worst coverage. The failed refinement retained
the previous controls for a diagnostic coupled execution; it did not create
a repaired candidate.

Coupled propagation from the exact epoch with those same controls covers at
worst **92.714%** of the assigned initial service receiver and reaches the
clearance floor at **8.925115 h**, pair 58/96. It differs from the linear
proposal by **5,318.37 m**. Nonlinear propagation with the same reduced force
differs from that proposal by only **20.91 m**, whereas coupled versus
nonlinear reduced propagation differs by **5,310.63 m**. This implicates the
omitted mutual-shadow force response. The comparison also changes finite-Sun
and finite-area forces, so it is not an exactly isolated shadow experiment.
For scale, the separate passive finite-area/source scout comparison differs
by 0.37 m at twelve hours, 3.47 m at 48 hours and 49.27 m at ten days.

Thus ray ownership was an avoidable restriction in earlier fits, but removing
it alone did not provide a usable trajectory. Sparse optical constraints and
uncoupled force response both generated false optimism. The estimated
36.528 TJ is **not** a measured cycle cost or a lower bound on one.

## Best executed local candidate and numerical limits

The slower common turn preserves the actual initial arrangement and six-hour
service, then turns over hours 6–12 instead of hours 7–11. DOP853 executions
use eight/16 solar samples, 120/60 s maximum integration steps and relative
tolerances 2e-10/2e-12; physical-square checks refine intervals to 0.5 s near
close approaches. Both start from the original loaded epoch. No endpoint is
reset. The tighter calculation is a separately propagated numerical replay,
not an independent physical model.

| Gate | Measured result |
|---|---|
| Initial service | Six hours; sampled finite-Sun coverage 100%; 300 extra limb/date probes have 8,334.88 m positive minimum margin, with 0.188 m circle approximation allowance |
| First twelve hours | 36.1863508 TJ, zero electric translation; 1,722.480 m minimum sampled physical clearance |
| Numerical allowance through hour 12 | 6.973 m position/reconstruction disagreement; conservative 60 s all-pair interval bound after allowance 374.010 m, above 100 m |
| Rates and loads to the stop | Maximum sampled rate 0.0078125 deg/s versus the 0.1 deg/s bound; maximum torque 1.11937 MN m; no individual actuator overload |
| Return failure | 100 m at 13.799825 h, pair 325/347; coarse/fine event difference 0.773 s and maximum position disagreement 8.918 m |
| Uncertainty at the failure | Conditional clearance after allowance falls to 73.510 m; the complete 13.8-hour prefix does not pass the 100 m floor with uncertainty |
| Terminal/next-service/repetition | Not reached, not executed, not assessed as a successful recurrence |

The 374 m interval bound is deliberately conservative at 60 s. It should not
be confused with the 1.72 km actual sampled separation or compared as a
geometric loss against the old, more finely subdivided interval bound.
Optical tests sample source angles and dates; no continuous optical
certificate is claimed. General oriented shadows, photon forces and contact
come from the tested `oriented_tiles.py` implementation; the new pruned
all-pair distance search agrees with an independent brute-force pair scan.
Grouped control remains only four row bands and one pulse shape, not 361
arbitrary attitude histories. Its estimated energy omits integrated shadow
torque and is explicitly incomplete.

## Complete accounts of the executed portions

These are local measured costs through each actual stop. A smaller number
from a shorter failed execution is not evidence of a cheaper operating cycle.

| Quantity | Slower-turn execution | Cheap Sun-facing execution |
|---|---:|---:|
| Executed duration | 13.799825 h | 8.925115 h |
| Electrical control energy | 36.189294 TJ | 10.127181 TJ |
| Mean electrical demand | 728.457 MW | 315.190 MW |
| First-intercept bolometric-equivalent energy | 1.249876e18 J | 1.155370e18 J |
| Redirected spectral-band energy, a subset of that interception | 1.724205e17 J | 1.593833e17 J |
| Mean electric translation impulse | 0 m/s | 0.0680678 m/s |
| Mean attitude equivalent impulse | 0.527611 m/s | 0.0791830 m/s |
| Ideal exhaust | 56,294.46 kg | 15,753.39 kg |
| Modeled propulsion waste heat | 10.856788 TJ | 3.038154 TJ |
| Maximum individual electrical demand | 7.82247 MW | 12.46046 MW |
| Minimum member installed/peak ratio | 6.21342 | 5.60903 |
| Members exceeding installed power | 0 | 0 |

The account integrates first-intercept light once, with mutual shadow unions
and the same physical apertures. Redirected energy is not added again as a
second collection source. Eight-source/120 s versus 16-source/60 s energy
integrals differ by at most 0.036%; optical integrals differ by at most
0.025% across these executed stages. The stored per-member chronology
supports peak checking as well as net energy. Propellant depletion is
reported but not propagated; rigid-body/ideal electric torque couples remain
the actuator model. Film absorption, radiator design, flexible structural
response and plume coupling are not closed by the waste-heat figure.

Actual generation, delivered bus power, storage state and full-cycle cost
are **null**, not zero. No collector or storage installation is added in
this round because both physical executions fail earlier geometric/optical
gates. [Collection](collection.md) and [electromagnetic](electromagnetic.md)
retain local Sun-tracking collection with short buffers as a candidate;
the prior corridor's 27.067 TJ eclipse deficit still demonstrates why cycle
net energy is insufficient. A surviving trajectory needs placed collectors,
conversion and transmission losses, storage dispatch, thermal loads, added
mass and changed photon forces replayed together before supply acceptance.

## Coverage and inventory consequence

Replacement shielding is required from **hour 6**; no continuation here earns
service credit during departure or return. There is no executed re-entry or
handover. The local spectral filter, seams, dimensions and finite-Sun receiver
requirement are preserved; the unprotected interval is not excused by the
formation's permission to leave service. No global exobase-protection claim
follows from this local calculation.

For illustration, the edge-a opportunities carry these necessary inventory
charges if one group could actually repeat a six-hour passage at that cadence:

| Start-to-start cadence | T/6 group lower bound | Undivided groups | Fractional tile lower bound | Loaded mass lower bound |
|---|---:|---:|---:|---:|
| 46.139 h | 7.690 | 8 | 2,776.0 | 17.404 billion kg |
| 92.417 h | 15.403 | 16 | 5,560.4 | 34.861 billion kg |
| 138.674 h | 23.112 | 24 | 8,343.6 | 52.310 billion kg |
| 184.911 h | 30.818 | 31 | 11,125.5 | 69.751 billion kg |

This table assumes interchangeable six-hour service credits. Spatial matching,
phase assignment, full-target coverage and handover are unproved, and all
listed free trajectories already failed before those dates. It is an inventory
penalty, not a validated fleet. With no complete cycle there is no measured
cycle energy per recurring service hour or achieved fleet demand. **238.35 TW**
remains the held-array comparison, **23.835 TW** the operating target and
**17.29 million tiles** an optimistic capacity relaxation.

## Reproduction, resources and next test

The sixteen explicit products in [natural_checks.json](natural_checks.json)
bind producers, constants, parent products and ignored raw arrays. Run the
`natural_*` producer modules individually in manifest order with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`; `natural_execute`
is invoked with `--case face_refined`, `--case free_slow`, then
`--case free_slow --fine`. Existing `joint_seed`/regression raw inputs must
first be reproduced and match their recorded hashes. `natural_scout` is an
intentionally interrupted attempt, not a required complete ten-day rerun.
The wrappers retain the exact reduced-model and ray-assignment adaptations.

The early checkpoint is `7bddd453a8f32da08ce66aa59be092eefbd7ccc4`.
The ledger charges the full 400 s allocation for the interrupted scout rather
than double-counting its 269.66 s partial product. Completed producers use **998.494 CPU s**; with the interrupted attempt and
20 s ancillary allowance, the charge is **1,418.494 / 1,800 CPU s**. Peak
recorded RSS is **663.734 MiB**, under the 2 GiB address-space limit; new raw
output is **138.943 / 512 MiB**. No larger simulation was launched.

All **11 new tests pass**. Final `make check` reports **745 passed, 55 skipped
and the same 13 baseline failures**, with 599 files passing the layer check.
Missing spectral/configuration inputs and the existing process-lookup failure
remain; later Makefile targets were not reached. An initial full-suite run
exposed HiGHS process-global thread-scheduler interference in two new tests.
Their setup now matches the fresh single-thread producer processes, and
infeasibility assertions require the actual infeasible status. All 33 saved LP
outcomes are optimal or infeasible, with no scheduler errors. Numerical
producers and their recorded products were unchanged by this test repair.

**Smallest remaining test:** one coupled, earlier layer/row-dependent velocity
phase adjustment on the slower-turn trajectory, targeting the measured
compression while retaining the initial full-receiver optical service.
Use the 15.673 TJ saved within the original twelve-hour budget as a ceiling,
not a promised allocation, and include mutual-shadow response in optimization
from the outset. Retain collision cuts and full receiver-hole constraints;
test the natural next passage only if that local maneuver survives. The
non-affine residual warns against assuming one velocity gradient will suffice.
Any survivor must then pass actual supply and the next-service terminal-state
and repetition tests. No larger fleet campaign is authorized by these results.
