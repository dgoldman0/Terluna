# Active assistance at 50 g/m²

**Active control can close local finite-tile gaps and produce roughly 1,400 km
of orbital retiming with several metres per second of correction. A complete
four-lunar-radius fleet schedule remains unresolved.** The author resumed the
investigation on 4 October 2026 to test whether modest active assistance could
improve the preceding 61.76% geometric interception result.

The active-control runs use separate **10 km × 10 km tiles**, retaining the
engineering reference scale. Every
trajectory is integrated at its actual date in the same DE440 ephemeris.
The base optical areal mass remains **0.05 kg/m²**, or **5 million kg per tile**.

The earlier [fleet probe](fleet.md) already changed sail attitudes during
service and coast. Its electric translational thrust was zero. The added
actuator here is bounded electric acceleration, with its energy and propellant
charged explicitly. The unchanged ideal spectral filter and reflected-beam
guards remain part of the force model.

## What the earlier gaps require

The [gap measurement](results/active_gap_scale.json) evaluates the actual
simultaneous footprints of the 3,072-square geometry probe at six dates, using
16 solar directions and a 129 × 129 receiver grid plus boundary witnesses.
At day 29.25, an identified ray is **1,406 km from the nearest existing
footprint**; the largest identified distance among these dates is **1,429 km**.
These are sampled witnesses, with no certification of the global maximum.
Their length scale belongs to that particular giant-square arrangement; a
differently populated 10 km fleet can have different gaps.

At day 29.25, even allowing every footprint an optimistic **10 km expansion**
raises its sampled ray reach only from **50.58% to 51.82%**. A 100 km expansion
raises it to **61.78%**. This diagnostic expands effective area and ignores
assignments and collisions, so it is a reach envelope rather than achieved
coverage. It shows why correcting seams alone will leave the large holes.

## A controlled formation of 625 finite tiles

The local formation has 25 rows and 25 columns at **9,750 m pitch**. Its four
checkerboard depth levels are separated by **1 km**. This permits neighbouring
optical footprints to overlap while their physical membranes remain separated.
Each 10 km square has a **9,890 m clear optical side**, retaining the existing
10 m seam and 50 m edge allowances.

A reference centre follows a natural 15,000 km lunar orbit through the sunward
service region. Each tile has its own position, velocity, gravitational force,
mutual illumination and electric command. The control law combines the force
needed to follow the assigned formation with damped position/velocity feedback;
its acceleration cap is **0.0005 m/s²**. This is a prescribed-geometry control
example, without an energy-optimality claim.

The moving receiver window is **100 km in diameter**, entirely within the
four-radius lunar aperture during this passage. Its area is about **0.0052%**
of the global target disk. The surrounding tiles supply
the finite-Sun margin. That window is a local coverage test, not a complete
shield aperture. No local-window area is multiplied into a claimed global fleet.

In the [six-hour maintenance run](results/active_formation.json), initial
position errors are 20 m. Every half-hour coverage check gives full interception
on the 32-direction solar grid. The swept finite-square guard is evaluated
every minute and flags **zero possible collision pairs**, with at least
**826 m conservative separating clearance** against a 100 m requirement.
Reflected cones retain positive Earth and Moon margins, including solar size,
tile corners, target motion and the assumed pointing allowance.

| Six-hour maintenance quantity | Result |
|---|---:|
| Tiles / base optical mass | 625 / 3.125e9 kg |
| Mean / peak tile acceleration | 141 / 281 µm/s² |
| Mean / largest tile integrated correction | 3.045 / 5.413 m/s |
| Mean / peak total propulsion power | 13.35 / 14.52 GW |
| Propellant over six hours | 448,566 kg |
| Peak-power hardware estimate | 60.49 million kg, 1.94% of base optical mass |

The lower-cap tests reach a 1 km displacement guard within an hour and stop
because the bounded optical-neighbour stencil no longer has its required
tracking bound. Their sampled receiver windows are still covered when they
stop. They establish failure to hold this assigned geometry, **not a proof
that those thrust caps cannot maintain useful coverage under a different
choreography**.

## Closing an explicit local gap

The [acquisition test](results/active_acquisition.json) begins with alternating
400 m and 300 m transverse offsets, placing every tile 500 m from its assigned
position and creating actual cracks between optical footprints. Initial mean
ray interception is **94.33%** on the 32-direction grid.

Bounded control reaches full sampled coverage by **50 minutes**; all later
ten-minute coverage checks remain complete through the three-hour test. All tile
position errors fall below 50 m by **61 minutes**. The swept guard flags zero
collision pairs, with at least **858 m** conservative clearance. Mean propulsion
power is **14.64 GW**, including the ongoing cost of holding the prescribed
formation. The mean integrated correction is **1.669 m/s per tile** over three
hours. Initial coverage during acquisition remains incomplete.

Mutual shadowing changes the force on downstream tiles. It is included here
through perspective unions of finite rectangles for each sampled solar
direction. The reference stencil was tested against all-pair polygon unions.
The force calculation assumes an ideal spectral mirror, transparent clear-band
transmission and no additional reflected-light interactions between tiles.

## Retiming a tile by about 1,400 km

The [orbital retiming runs](results/active_retime.json) address the larger
redistribution scale. A 10 km tile receives **four smooth six-hour thrust arcs**
over either three or seven days. The target is an advanced orbital phase
specified at the same final ephemeris date. Its actual controlled trajectory
is then integrated and independently replayed. A rotated target state supplies
a boundary condition; it is never counted as another fleet member.

The force calculation includes a 3 × 3 surface quadrature for finite-extent
gravity. The reference and controlled histories share their actual initial
epoch, with separate shooting runs for each lead time. These solutions target
arrival position and velocity; their thrust histories have not been minimized
for energy or propellant.

| Per 10 km tile | Three days | Seven days |
|---|---:|---:|
| Projected position shift at arrival | 1,398.1 km | 1,397.4 km |
| Total integrated correction | 5.450 m/s | 6.790 m/s |
| Independent arrival position error | 0.338 m | 0.020 m |
| Independent arrival velocity error | 0.0109 mm/s | 0.0048 mm/s |
| Peak electric acceleration | 149.85 µm/s² | 202.60 µm/s² |
| Electrical energy | 0.2294 GWh | 0.2858 GWh |
| Mean power over the full maneuver | 3.186 MW | 1.701 MW |
| Peak power | 22.71 MW | 30.70 MW |
| Propellant | 1,285 kg | 1,600 kg |
| Propellant / base optical mass | 0.0257% | 0.0320% |

The seven-day case has lower average power because the maneuver spans a longer
time, but it uses **more total energy and propellant** in this search. These two
cases do not establish a monotonic benefit from extra lead time. An initial
undamped two-arc shooting attempt missed its endpoint and was rejected; the
[search record](results/active_search.json) preserves that failure.

Both successful cases remain outside the protected sphere and retain positive
sampled Earth/Moon reflected-beam margins. They show individual reachability
at the identified gap's length scale. They do not allocate enough tiles to fill
a global gap, avoid other fleet members during transfer, or complete a
repeating return schedule.

## Energy, inventory and validation

All electrical figures use the existing propulsion assumptions: **30 km/s
exhaust speed, 70% efficiency and 45° cant**, with power hardware estimated at
300 W/kg and a 25% peak margin. Integrated impulse determines propellant and
electrical energy. Collector supply, habitat consumption and delivered power
are separate, unevaluated accounts. Power hardware and propellant are reported
without feeding their extra mass back into these 50 g/m² trajectories.

The [validation product](results/active_validation.json) confirms full sampled
coverage at 50 and 180 minutes with **128 solar directions**. The initial
94.33% mean-ray statistic retains its 32-direction resolution; the initially
perforated formation's full-Sun area is not claimed converged. Finite sampling
still leaves unsampled times and directions outside these checks.

The six-hour maintenance replay halves the maximum integration step, tightens
relative tolerance by 100-fold and doubles the force's solar quadrature.
Computed positions agree within **0.001 m**. The controller uses exact modelled
forces and states, so this primarily checks numerical execution of the
prescribed path. The omitted finite-extent gravity term for its 10 km squares
is at most **2.43e-9 m/s² at the checked dates**; the separate orbital-transfer
runs include that term explicitly.

Across the same first three hours, the initially gapped formation costs
**1.669 m/s** per tile, compared with **1.510 m/s** for the initially assembled
formation. Most of this example's thrust pays to maintain the prescribed
geometry throughout the passage.

Inventory makes small per-tile costs significant. As conditional bookkeeping,
applying the seven-day maneuver once to optical mass equal to one four-radius
**target disk's area** would average **2.58 TW** over that week. Applying it to
the entire **1.536e14 kg** inventory represented by the giant-square probe
would average **52.26 TW**. If a fraction *f* of that inventory needs one such
maneuver per week, the corresponding average is **52.26 f TW**. The required
*f*, sufficient fleet inventory, return costs and maneuver recurrence are
unsolved. These scaled masses are not coverage demonstrations or available
power supplies. At the study's specific-power assumption, the transfer's
peak-power hardware adds approximately **1.9%–2.6%** of a tile's base optical
mass before other systems are included.

The local controller assumes accurate state estimation and an actuator that
can deliver the commanded acceleration. Slew dynamics, flexible membranes,
thruster placement, exhaust paths, heat rejection and failure recovery remain
open. The refined numerical replay tests integration accuracy, not operational
navigation or hardware performance.

The results support investigating active assistance further. The next coupled
calculation must assign many finite tiles to changing receiver locations,
enforce collision-safe transfers and handovers, retain the shadow while other
tiles are moving, and count all return legs. The number of moves and their
recurrence determine the fleet's operating cost. The earlier 61.76% geometric
probe remains an unaccepted population with encounter conflicts.

## Reproduction

Install the study requirements and restore the pinned ephemeris as described
in [fleet.md](fleet.md). The gap diagnostic reads the earlier verified raw
`dense_double_1000km` trace. New raw trajectories remain in ignored
`research/runs/solar_shield_array/active/`.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_gap_scale
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_formation_run
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_formation_run --caps 0.0005 --suns 16 --max-step 30 --rtol 2e-12 --output active_formation_replay.json
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_acquisition
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_retime_run
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_validate
python -m pytest protection/dynamics research/studies/solar_shield_array -q
python visualization/solar-shield-array/active_plot.py
```

Products record producer sources, constants, raw trace hashes and their evidence
boundaries. [checks.json](checks.json) records the actual repository checks.
