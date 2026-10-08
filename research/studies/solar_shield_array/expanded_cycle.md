# Expanded coupled cycle investigation

Continue `762ce4d810e07ddc0788f315212c5a6f6bdca5fe` on
`research/solar-shield-habitat-array`. The author waived the preliminary
remote-head/concurrent-work check and authorized an early and final push.
Main remains unchanged. Primary magnetic integration remains set aside.

## Bounded decision experiment

Use the source-checked `coupled_restart.restore()` and `joint_seed.json`,
all 361 loaded 10 km physical squares, 9.89 km clear apertures and unchanged
spectral allocation. Keep the six-hour first service, 100 m surface clearance
plus measured numerical allowance, individual installed actuator allocations,
0.001 m/s² translation and 0.1 degree/s attitude limits. Preserve the 51.859 TJ
first-twelve-hour ceiling and 100 TJ through-next-service screening ceiling.
The selected zero-translation slower-turn comparator remains 36.186 TJ.

Allocate 1,800 numerical CPU seconds, one process and numerical thread,
2 GiB address space and 512 MiB new raw output. Reserve at least 400 CPU
seconds for independent replay, accounting and compact export; required
repository checks have the inherited separate five-minute allowance per
commit. Charge compilation, incomplete evaluations and diagnostics. Measure
one complete coupled evaluation and directional sensitivity first; use those
costs to limit differentiated directions and iterations explicitly.

Add independent in-plane separation/shear and arrival position/velocity
acquisition directions to the available early normal controls. Use smooth
finite burns and common finite-rate attitude changes. Arrival directions may
restart at their exact onset on the common coupled baseline, because their
prefix controls are identical. Every changed suffix still propagates all
361 interacting tiles through next service and the subsequent departure.
The full-cycle Jacobian must use coupled evaluations, never a gravity-only
transport of a changed prefix. Record any unmeasured columns separately.

Inspect singular values, effective rank, residual projection and finite-step
prediction error. Use constrained trust-region steps with energy in the
objective, explicit restoration slack and exact coupled trial replay.
Positive slack, missed optical service, resource violations, or contact
preclude acceptance. Continue past contact only for labelled optimization
diagnostics. Retain exposed contact choices and rejected trials. Diagnose a
small set of one-limit relaxations only after the actual-limits subproblem;
relaxed solutions never earn physical acceptance. Hardware changes require
new mass, photon forces and epoch replay; actuator ratings never imply supply.

## Repeatability gate

A safe service/departure/return/next-service/subsequent-departure sequence is
only a finite execution. Test transported relative position, velocity and
attitude at service entry, service exit and corresponding departure phase,
along with the receiver's actual-date position inside four lunar radii.
Use explicit finite-sequence closure tolerances of 50 m, 0.01 m/s and 1e-4 rad
as *search acceptance gates*, rather than established engineering tolerances.
The actual-date ephemeris is nonautonomous: small relative-state closure at
one return alone cannot establish indefinite recurrence. Recurrence evidence
requires further consecutive actual-date cycles from the computed endpoint,
without resetting state, preserving every service, safety and resource gate
and bounded defects. Report separately whether that extension was reached.

No generator, collector or storage system is placed in this experiment.
Generation, storage and unmet demand remain unknown. Track executed control
energy, individual loads and first-intercept/redirected light histories.
Preserve 238.35 TW held-array comparison, 23.835 TW operating target and the
unvalidated 17.29-million-tile estimate; do not extrapolate failed prefixes
into achieved fleet operation.

## Portable continuation

Preserve the expanded basis, coefficients, local responses, optimizer radius,
contacts, scenario and ephemeris identities, deterministic seeds, compact
ephemeris samples, dependencies and exact restart/replay commands. Provide a
CPU-budget-configurable runner. Check loading and representative replay with
the original ignored trajectory files absent. Prioritize larger local work
from the measured rank, model errors and active constraints. Stop after this
bounded investigation; do not launch the local campaign.

The early `make check` passes the layer check (615 files, zero violations)
and stops because this fresh environment lacks `pytest`. No Python test
result is claimed for that invocation. Restore the test dependency before
the final required check.

## Finding

**Unresolved search. No accepted return, next service, subsequent departure,
repeatable cycle or physical resource bottleneck is established.** The added
arrival controls reduce the diagnostic cycle error, but the formation still
fails clearance and the linear predictions degrade strongly during the next
service and departure. The measured finite-difference control space also
leaves substantial residual directions unresolved. Installed actuator
capacity has substantial headroom in both the actual-limits proposal and
its executed prefix. The 100 TJ next-service energy cutoff is active in the
local solve; lifting it does not resolve that solve's clearance slack.

The early checkpoint was pushed as
`f22a9f6f1382a0601a64c8d76b05a9236ea7dfb5`. The preserved physical selection
remains the original slower turn. The failed expanded controls and all local
diagnostics are retained separately. Main and the spectral products are
unchanged.

## Control directions and measured cost

The implementation offers 18 dimensionless parameters. All evaluations retain
361 coupled, hardware-loaded squares; there is no reduced tile population.

| Directions | Purpose | Evaluated columns in this run |
|---|---|---|
| Six saved normal modes | Change depth-class compression and selected pair separations | Row parity only; the other five remain available |
| Row gradient along in-plane u; column gradient along in-plane v | Change edge overlap independently of normal compression | Both |
| Opposite in-plane corrections on 325/347 and 338/341 | Localize an edge escape without forcing every row to move | Implemented, unmeasured |
| Three arrival position patterns | Correct each body-frame component of the baseline arrival-position error in an earlier acquisition arc | All three |
| Three arrival velocity patterns | Correct each component of the velocity error in a later arc, separating position and velocity authority | All three |
| Arrival date; arrival-turn duration | Change acquisition timing and common attitude history | Implemented, unmeasured |

Early modes have maximum amplitude 0.1 m/s per unit coefficient, act over
hours 6–12, and repeat at the following departure. Acquisition modes have
maximum amplitude 1 m/s per coefficient, with windows arrival minus 8–4 h
and arrival minus 4–0 h. Each is a mass-centred, per-tile residual pattern;
these six patterns do not supply arbitrary independent controls to every
tile. Commands use normalized sine-squared pulses. Timing modes retain the
specified nominal centre targets; general receiver-path and centre freedom
are further variables. Common attitude, the loaded initial packing and tile
identities remain fixed in the measured subset.

The initial full coupled continuation at 120-second maximum steps cost
**289.65 CPU s**, with 290.93 s including setup/export. The larger-step search
benchmark cost **127.24 CPU s**, using the same forces, 8 Sun samples,
16×16 area sampling in partial eclipse and relative tolerance 2e-10, with a
480-second maximum step. Its maximum disagreement with the 120-second
reference is **1.351 m / 5.53e-5 m/s** over the complete diagnostic sequence;
the first-event time differs by 0.0138 s. Its interpolation error is 0.0157 m.
This numerical comparison supports a cheaper search evaluation; it does not
establish derivative accuracy.

Three early full-cycle perturbations cost **112.04–125.67 CPU s each**.
The six acquisition perturbations cost **53.10–62.95 CPU s each**, restarting
from exact saved baseline states at their onset. Their unchanged prefixes
are reused. All changed suffixes include coupled next-service and subsequent-
departure dynamics. The common coefficient perturbation is 0.2, corresponding
to 0.02 m/s early and 0.2 m/s in acquisition. Nine of eighteen columns are
measured; no approximate tail Jacobian fills the missing columns.

## Rank, residual directions and trust failure

The terminal vector retains all 361 position/velocity residuals at next
service entry, exit and subsequent departure, normalized by 100 km and
5 m/s. The optimizer also includes finite-Sun whole-receiver coverage at
13 next-service dates × 25 sources, 34,121 retained contact half-spaces,
individual power/thrust constraints and both energy ceilings. Initial service
is unchanged by construction. Coverage and terminal slacks must vanish before
acceptance. No virtual accelerations are used.

| Measured matrix | Effective rank at relative threshold 1e-3 | Norm of chosen residual outside that retained span |
|---|---:|---:|
| Full-cycle terminal conditions | 9 of 9 measured columns | **34.246%** |
| Arrival-only terminal subspace | 6 of 6 | 34.527% |
| Three measured early controls, terminal conditions | 3 of 3 | 99.985% |
| Complete retained clearance residual | 8 of 9 | 99.634% |
| Early clearance through hour 15 | 2 of 3 early columns | 99.981% |
| Sampled next-service coverage deficits | 7 of 9 | 1.350% |

Terminal singular values span 41.49 to 0.724; its rank stays nine from
relative thresholds 1e-2 through 1e-6. Clearance and coverage ranks change
with threshold. Units and coefficient scales matter. The clearance projection
uses the negative part of the retained gap deficits: an unspanned deficit
vector is **not an infeasibility proof for inequalities**. These ranks concern
a finite-step, one-sided, locally frozen formulation, including diagnostic
post-contact states. The nonlinear tests below limit how literally its span
can be interpreted as infinitesimal control authority.

The actual-limits SLSQP terminates successfully in eight iterations **with
positive restoration slacks**: 1,157.5 m clearance slack and 2.85589 normalized
terminal slack. It predicts 100 TJ through next service, 134.429 TJ through
the following departure, and a maximum power/rating ratio of 0.2680. These
are failed-restoration predictions, not feasible-cycle costs.

| Diagnostic coupled replay | Original baseline | Constrained proposal | One tenth of proposal |
|---|---:|---:|---:|
| First nominal clearance-floor event, h | 13.799614 | 13.794785 | 13.799131 |
| Maximum position debt at next-service entry, km | 319.422 | 308.650 | 318.259 |
| At next-service exit, km | 343.293 | 320.463 | 340.997 |
| At subsequent-departure endpoint, km | 344.511 | 314.318 | 341.484 |
| Maximum velocity debt after subsequent departure, m/s | 12.5670 | 11.8554 | 12.4953 |
| Normalized terminal RMS | 1.12113 | 1.06354 | 1.11532 |

The complete constrained replay has a **27.571 km** maximum terminal
prediction error, growing from 585.7 m at arrival to 8.102 km at service exit.
It exposes **250 new contact rows**, and the retained-contact slack grows to
2,878.4 m instead of improving. Its actual/predicted merit ratio is **−4.59**.
The step is rejected and the trust radius is reduced. Its hypothetical next
service covers only **8.831%** at worst, despite the receiver remaining inside
four lunar radii; the prediction was 10.389%. No service credit is assigned.

The held-out tenth-step still has **2.861 km** maximum prediction error
(54.5 m at arrival, 834.9 m at service exit). Error falls by a factor of
9.64 when the step falls tenfold. That is consistent with substantial
finite-stencil bias or nonsmooth response, rather than an adequate local
derivative cured solely by shrinking the optimizer step. Two amplitudes do
not isolate stencil bias, force curvature and changing shadow topology.
The changing contact set is directly observed. The 1.35 m integration
benchmark disagreement is much smaller than these prediction failures.

The failures therefore include missing measured control directions, inaccurate
local full-cycle predictions and changing contact constraints. SLSQP solved
its slack-bearing subproblem; it did not stall or produce a feasible solution
that was subsequently mislabelled. A smaller, signed derivative stencil and
contact-specific control enrichment should precede a broad random search.

## Resource diagnostics and continuation

These are explicitly **local-model diagnostic solves**, made after rejecting
the actual-limits replay. Their large-step predictions are not promoted into
physical resource bounds. No extra relaxed-resource trajectory or hardware
installation is accepted.

| One limit changed | Predicted next-service energy | Clearance slack | Maximum power / original rating | Result |
|---|---:|---:|---:|---|
| Actual limits | 100.000 TJ | 1,157.50 m | 0.2680 | Failed restoration |
| Next-service energy ceiling 100 → 300 TJ; prefix ceiling unchanged | 143.871 TJ | 1,158.26 m | 0.4363 | Terminal RMS falls to 1.03748; clearance remains unresolved |
| Individual actuator ratings ×2; energy ceilings unchanged | 100.000 TJ | 1,157.50 m | 0.2680 | Identical coefficients and objective |
| Continue energy-relaxed start back to 100 TJ | 100.000 TJ | 1,157.50 m | 0.2680 | Returns within 2.3e-6 in coefficients of actual-limits result |

The 100 TJ cutoff affects the local arrival trade. Actuator capacity does not
bind these tested solves. This separates an active energy objective/constraint
from a demonstrated physical bottleneck: none is demonstrated here. More
power equipment cannot be inferred to cure the observed model and contact
failures. Doubling the installed 29.1895 GW allocation would add **97.2985
million kg** under the inherited 300 W/kg assumption, about 4.30% of the
loaded group mass. That mass and its changed photon acceleration have **not**
been propagated, so the rating relaxation is not a hardware-feasible case.
Available electrical supply remains unspecified in every case.

## Independent physical prefix and accounts

The actual-limits controls are independently replayed from the original epoch
using 16 Sun samples, a 32×32 partial-eclipse area grid, 60-second maximum
steps and relative tolerance 2e-12. The run stops at the first clearance event.
All later state, coverage and energy values above are diagnostic only.

| Executed gate or account | Result |
|---|---:|
| Initial service, 49 dates × 25 sources | 100% sampled whole-receiver coverage; inside four radii |
| Twelve-hour energy | **36.614024 TJ**, versus retained 36.186351 TJ |
| Twelve-hour minimum sampled square separation | 1,720.493 m |
| Twelve-hour adaptive bound after numerical allowance | 368.465 m |
| First nominal 100 m event | **13.79499558 h**, pair **325/347** |
| Coarse/fine event-time difference | 0.7580 s |
| Maximum position disagreement including reconstruction | **8.93656 m** |
| Conditional bound through 120 s before stop | 120.304 m |
| Final interval bound after numerical allowance | 73.383 m; fails 100 m gate |
| Electricity through stop | **36.616960 TJ** |
| Attitude / translation electricity | 36.189278 / 0.427681 TJ |
| Maximum individual demand / minimum installed-to-peak ratio | 7.90008 MW / **6.15369** |
| Overloaded members | 0 |
| Maximum executed translation / attitude rate | 1.0535e-6 m/s² / 0.00781248 degree/s |
| Maximum required torque | 1.11937 MN m |
| First-intercept light / redirected subset | 1.249874e18 / 1.724202e17 J |
| Ideal exhaust / modeled propulsion waste heat | 56,959.7 kg / 10.985088 TJ |
| Carried group mass | 2,263,298,487 kg; unchanged |

At the final encounter, the relative body-frame displacement is
(−6,314.451, 9,022.947, 100.000) m. Both in-plane projections still overlap.
The adaptive all-pair square audit refines close intervals to 0.5 s, subtracts
twice the measured positional allowance and retains the inherited conditional
0.15 m/s² acceleration and 0.1 degree/s rate envelopes. No safe extension of
the complete cycle is accepted. The source/date optical checks remain sampled
tests; the circle polygon has a 0.18825 m sagitta allowance.

Actual terminal positions, velocities, quaternion and angular velocity are
saved for every member. Commanded service-attitude error is zero by the
prescribed-attitude model; independent attitude tracking/axis allocation is
not solved. Torque is charged through the installed ideal electric-couple
capacity together with translation. Detailed torque-axis ratings, structure,+plume layout and propellant depletion remain inherited engineering limits.
The finite-sequence closure tests fail by kilometres and metres per second;
the next actual-date recurrence test is not reached. Zero further cycles are
credited. Generation, bus delivery, storage and unmet demand remain unknown.

## Resources and handoff

Numerical execution records **1,812.607 CPU seconds against 1,800**, an
**overrun of 12.607 s (0.70%)**. The final validation completed, but the
timer-only stage guard did not enforce its configured ceiling. Research
evaluations stopped on discovery. The final runner adds explicit process-clock
checks inside dynamics, load accounting, optimization and validation; it also
charges compilation children. The saved ledger retains this resource-control failure. Peak recorded RSS is **386.12 MiB** and new raw output is **290.34 MiB**;
the memory and output limits are satisfied. The two artifact exports take 8.06 CPU s
and required repository checks are recorded separately. `budget_guard.patch` preserves the exact source
change; reversing it reconstructs the pre-repair execution sources.

The compact [local continuation package](local_continuation/README.md) contains
the checked original restart/seed references, expanded basis, coefficient
vectors, contact and terminal/coverage Jacobians, trust state, six days of
pinned DE440 samples, dependencies, hashes, deterministic seed, exact replay
and validation commands, measured costs and all-member executed power,
torque and interception histories. The history is split into two hash-checked
byte parts to meet the upload limit; reassembly is verified byte-identical.
It can resume the local SQP directly or
rebuild coupled responses without the ignored Work trajectories. The preserved
physical selection is the zero-control slower turn; failed controls remain
explicitly labelled. No local campaign was launched.

Prioritize the larger local work as follows:

1. Benchmark signed 0.02 and 0.01 coefficient stencils through the complete
   sequence, especially acquisition-normal directions. Demand held-out error
   convergence at service exit and subsequent departure before trusting rank
   or resource conclusions. Cache each response separately and preserve an
   independent validation reserve.
2. Add the unmeasured localized in-plane pair directions and remaining normal
   modes, then enrich with modes aligned to actual deficient contact rows.
   Compare alternative separating faces with exact replay. Keep full-service
   coverage and next departure in every model.
3. Rebuild at any accepted restoration anchor. Extend the runner's one-anchor
   SQP workflow into a budgeted relinearization loop, with actual/predicted
   tests, contact retention and no state resets. Reuse the measured model only
   at its recorded anchor.
4. Test timing, centre/receiver path and assignment freedom after the local
   derivatives are reliable. The known high-cost corridor is a geometric
   initialization option with all of its costs retained. Compare only a few
   justified energy/actuator relaxations and continue survivors to real limits.
5. A safe finite sequence must then satisfy matched-phase error gates and
   further consecutive actual-date cycles. Supply hardware requires placement,
   chronological generation/storage accounting and original-epoch mass/photon-
   force replay. Fleet coverage and handover follow that local demonstration.

The **238.35 TW** held-array comparison, **23.835 TW** operating target and
unvalidated **17.29-million-tile** estimate retain their prior meanings.
There is no new fleet extrapolation or continuous cold-exobase protection
claim. Review this checkpoint before starting the larger local campaign.

## Repository verification

The portability `make check` passes the layer check (615 files, zero
violations), then reports **751 passed, 64 skipped and 13 failed** in 40.58 s.
The failed identities exactly match the previous checkpoint: missing spectral
inputs, generated climate inputs and process lookup. Later Makefile targets
are not reached. All six new tests pass. The full check and portable `inspect`
run with the original `research/runs` tree temporarily absent; the documented
360-second representative replay writes its outputs into a fresh directory.
The first check exposed a custom-output-directory bug, which is repaired
before this successful portability check. Extra skips reflect the deliberately
absent ignored data/kernel. After splitting the history for upload, the final
check with the original ignored data restored reports **760 passed, 55 skipped
and the same 13 failed** in 39.88 s; the layer check covers 620 files with zero
violations. The final package inspection verifies both history parts and the
saved optimizer model. All three check attempts fit within the inherited
five-minute check allowance. [expanded_checks.json](expanded_checks.json)
records exact checks, hashes, failures, resource limits and the overrun.
