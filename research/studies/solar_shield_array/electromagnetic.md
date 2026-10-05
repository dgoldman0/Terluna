# Electromagnetic infrastructure feasibility

## Checkpoint and bounded work

Continue the existing research branch from
`5dc49e11e9ba1ab782c3904793ce052f30b9e4b4`. The author explicitly waives the
preliminary remote-head/concurrent-work check. The clean local checkout is at
that commit. Main and all prior failed cases remain unchanged.

The decision is whether sharing power distribution, relative control and
solar-wind protection equipment earns its mass and operating complexity.
Compare the present electric-propulsion reference, tile-distributed equipment
and a hub-assisted group. Larger collectors and protection hubs are candidates.
Retain the 23.835 TW comparison target, the failed return and the distinction
between optical interception, generation, available bus power and delivery.

The environment exposes eight CPU equivalents and 8 GiB cgroup memory. New
numerical work is capped at **1,200 s wall time, one process, one numerical
thread, 2 GiB address space and 128 MiB new raw output**. Required repository
checks have a separate five-minute allowance per pre-commit gate. Charge
initial data inspection at 10 s and record subsequent numerical runs,
including failed attempts. No coupled plasma or new fleet propagation belongs
to this pass; pause for review before a substantially larger campaign.

Read the research/domain entry points, `return_repair.md`, `closure.md`,
`packing.md`, `collection.md`, their propagation and load implementations,
the protection architecture, and the September field/structure/plasma account.
Use source-bound saved hardware-loaded states and load histories, with the
strongest return correction's tighter replay. Analyze a central 3 by 3 tile
group (zero-based indices 160–162, 179–181, 198–200); retain all 361 members
when checking power availability, geometry and external field interactions.
Also inspect the original and migrated encounter pairs to avoid treating a
central group as an encounter bound.

1. Reconstruct time-dependent loads, actual separations and directions,
   individual power shortfalls and optimistic pooled-bus limits. Screen local
   collectors and batteries against microwave/laser links and near-field
   induction. Include receiving area, pointing, losses, heat, mass, bus limits
   and the source of energy. Existing propulsion capacity is not generation.
2. Integrate finite square-loop fields and forces at representative measured
   configurations. Check convergence and action/reaction, compare distant
   dipoles only where valid, include torque, inductance, support, conductor,
   cryogenics, fault energy, and unwanted interactions. Vary separation,
   force and group size. Internal forces receive no net-fleet momentum credit.
3. Screen solar-wind pressure, kinetic scales and particle rigidity. Compare
   tile coils and a larger hub with the existing lunar protection architecture.
   Pressure balance and single-particle bending are diagnostics; they establish
   neither a magnetopause nor atmosphere/dose protection. Keep EUV and energetic
   radiation requirements separate. Include external wind/photon momentum.
4. Rank candidates using these calculations and primary research. State
   demonstrated elements, engineering extrapolations and exploratory choices.
   Quantify replacement and propellant sensitivities. Specify the smallest
   coupled follow-up if supported, without claiming a changed trajectory until
   changed optical forces, added mass and field/plasma forces are propagated.

Compact, source-bound results and reproducible code will accompany this report.
Raw analysis products go in `research/runs/solar_shield_array/electromagnetic/`.
No full-cycle, global-coverage or accepted return result is implied by a
feasibility screen.

The early checkpoint's required `make check` reports 695 passes, 55 skips and
the same 13 baseline failures in 32.81 s. Missing spectral inputs, generated
climate configuration and process lookup account for those failures; later
Makefile targets are not reached. No numerical implementation has changed.
