# Simultaneous fleet coverage at 50 g/m²

**Continuous coverage remains unachieved.** The study now propagates distinct
three-dimensional populations at common dates, computes unions of finite-tile
shadows over the finite Sun, and tests seams, reflected beams and swept close
approaches. The expanded populations expose substantial spatial gaps and
conflicts. Numerical phase disagreement also exceeds the proposed placement budget.

The primary target remains **four lunar radii**, as in this branch's scenario.
The established three-radius protected region from [the September report](../../../protection/report.md)
is included as a comparison. The finite-Sun aperture margin is additional to
the protected radius. The author confirmed that four radii is an acceptable
primary target on 4 October 2026.

The investigation is **paused at the author's request on 4 October 2026**.
This report preserves the completed search, its failed gates and the work
needed before resuming. No expanded search is running.

## Repair and checkpoints

The Boolean trapezoidal-integration repair was pushed as **2182374** before
fleet expansion. NumPy added adjacent Boolean ordinates as Boolean values;
casting to floating point before integration fixes the plateau contribution.

| Raw service-time fraction | Before | Corrected |
|---|---:|---:|
| Annual inner cycling case | 8.00394% | 15.82345% |
| Three-year inner cycling case | 8.03610% | 15.88746% |

Seven affected diagnostics were regenerated from source-, configuration-,
constant- and SHA-256-verified traces. Useful projected-area fractions,
trajectories, return-arc closure and inventory comparisons stayed unchanged.
Regression tests cover constant indicators and mixed plateaus on irregular
sample times, as well as the raw-time bound on weighted service area.

The simultaneous finite-tile pilot was pushed as **93b8ef2**. All work remains
on `research/solar-shield-habitat-array`.

## Populations and achieved geometric coverage

Every member starts on **2026-10-04 TDB** and experiences the same ensuing
**30 days** of DE440 forcing. Initial guesses are osculating lunar circular
states; plane normals span the retrograde hemisphere in the initial lunar
frame. Position and velocity for each phase are integrated separately in the
coupled batch. No trajectory is replayed at shifted dates to construct a fleet.

The model retains lunar and planetary point-mass gravity and the ideal
filter-spectrum sail at **0.05 kg/m²**. The physical square affects service
entry, finite shadows and beam exclusion. Its centre supplies the translational
force approximation. The finite-extent check below quantifies that limit.

[The search manifest](fleet_search.json) records five completed populations
and one rejected population. The 15,000/20,000/30,000 km set crossed the radial
guard. Its exporter stopped before retaining a diagnostic trace, so the result
is attributed to that population; a particular failing plane, member or
crossing date remains unidentified. The tighter 15,000/16,000/17,000 km set
stays between **14,462 and 18,073 km** over the month.

Four-radius interception, before enforcing member exclusions:

| Population | Squares | Side (km) | Base optical mass (kg) | Mean ray interception | Lowest sample | Mean full-Sun area |
|---|---:|---:|---:|---:|---:|---:|
| 1 plane × 576 phases | 576 | 10 | 2.88e+09 | 0.00542766% | 0.00180539% | 0% |
| 12 planes × 48 phases | 576 | 10 | 2.88e+09 | 0.00244955% | 0.00148359% | 0% |
| 3 radii × 12 planes × 48 phases | 1,728 | 10 | 8.64e+09 | 0.00589811% | 0.00349135% | 0% |
| 24 planes × 64 phases | 1,536 | 1000 | 7.68e+13 | 44.1662% | 32.4546% | 37.0419% |
| 24 planes × 128 phases | 3,072 | 1000 | 1.536e+14 | 61.759% | 50.6083% | 55.2783% |

This table uses the common six-hour / 16-solar-direction analysis. All single-
radius populations start at 15,000 km; the three-radius family uses 15,000,
16,000 and 17,000 km. Optical mass is area times 0.05 kg/m², without payloads.

The 10 km cases test the reference tile geometry and orbital families. Their
inventories are deliberately small compared with a complete shield. The
1,000 km squares are finite-geometry stress probes with much larger optical
inventories. Their structure, flexibility and extended-body dynamics require
additional models. A 1,000 km square's area equals 10,000 reference cells;
the independent motions of those 10,000 cells have not been simulated here.

The nominal populations contain close-approach and mutual-shadow conflicts.
Their table entries describe geometric interception on the propagated centre
histories. They cannot be accepted as functioning fleets. A higher sum of
projected areas can coexist with large holes because the same rays encounter
several panels while other rays encounter none. In the 3,072-square probe, the mean sum is
**123.716% of the aperture**, but its union intercepts only **61.759%** of rays,
with a lowest sample of **50.608%**. This directly tests the failure of an
aggregate-area argument on simultaneous three-dimensional histories.

For the first 1,536-member large-square probe, refinement from six-hour/16-Sun
sampling to three-hour/32-Sun sampling changes mean interception from
**44.166% to 44.108%**. The lowest sampled value drops from **32.455% to
31.435%**. At its weakest refined date, day **28.125**, about **62.2% of
7,213 mapped positions** receive no filtering on the 64-direction solar grid.
Every tested population leaves the aperture incomplete at sampled dates
throughout the month. Those failures suffice to reject continuous coverage.

## Finite Sun, seams and temporal gaps

Each date uses actual three-dimensional tile positions and normals. Rays
from a quadrature over the solar disk are projected through finite square
corners onto the protected aperture. Tile motion during light travel is
included to first order in barycentric coordinates. The fixed reception-time
normal, solar limb/spectral structure, and the circular approximation to
Earth's projected eclipse silhouette remain model limitations.

For each solar direction, the calculation takes the **union** of square
shadows. Their intersection over solar directions measures the region covered
for every sampled direction. Thus ray interception and full-Sun-protected area
are distinct outputs. The saved spatial maps retain positions and their
fraction of intercepted solar directions.

A 10 m seam plus a 50 m placement allowance at each edge leaves a **9,890 m
clear square** for a 10 km tile: a **2.188%** area allowance. Tile overlap counts
once. At 15,000 km an isolated 10 km square has no full-Sun umbra; neighbouring
shadows must maintain the shared aperture across the solar disk.

The first large-square probe's weakest four-radius date gives **31.431%,
31.438% and 31.435% ray interception** with 16, 64 and 256 solar directions.
Full-Sun-protected area at that date changes from **25.569% to 25.110% to
24.936%**. The full-Sun statistic is more sensitive to solar angular resolution;
its engineering tolerance must ultimately follow the atmospheric leakage budget.

[The refinement product](results/fleet_refinement.json) contains both three-
and four-radius results, three-hour histories, solar-grid comparisons and
spatial maps. Both targets use the existing four-radius service controller;
no three-radius trajectory optimization is inferred. For the 1,536-square
probe, the three-radius mean interception is **42.130%**, with a lowest sample
of **26.657%**; four radii gives **44.108%** and **31.435%**. These fractions
refer to different aperture areas, so a smaller target need not have a larger
intercepted fraction under the unchanged controller.

The [local temporal audit](results/fleet_handover.json) recomputes days
**26.0–26.5** at **one-minute cadence with 32 solar directions**, interpolating
positions and velocities from the saved dynamics. Ray interception ranges from
**31.532% to 50.177%**; the aperture is incomplete at all 721 dates. This is a
local audit around the weakest original six-hour sample, not a whole-month
minute-grid convergence result. The finer month-wide run finds its lowest
value at a different date.

## Collision and mutual-shadow exclusions

Every 600-second propagation-output interval receives a swept close-approach
check. Circumscribed spheres around the physical squares enforce a 100 m
clearance. The linear-segment minimum is enlarged by an acceleration allowance
of **0.15 m/s²**; measured centre accelerations stay below about 0.024 m/s².
The 10 km pair guard is **14,242 m**. These are conservative exclusions;
material intersection has not been proved for every flagged pair.

A second swept guard rejects possible mutual solar shadows, including the
solar cone, relative depth and bounded Sun-direction change. Intercepted
spectral photons change a downstream panel's sail force. Nominal populations
with these interactions require coupled illumination/force propagation.

Whole member trajectories are removed with a deterministic greedy independent
set. Its objective favours low conflict count and can retain members with
little service; it supplies one conservative selection, with no minimum-
inventory or best-coverage claim. The geometry of the retained subsets is
then recomputed from their actual histories.

| Population | Possible collision / shadow pairs | Collision-only subset | Both guards: members | Both guards: mass (kg) | Mean ray interception |
|---|---:|---:|---:|---:|---:|
| 1 plane × 576 phases | 5,208 / 26,904 | 262 | 61 | 3.05e+08 | 0.000593625% |
| 12 planes × 48 phases | 618 / 28,768 | 357 | 48 | 2.4e+08 | 0.000217887% |
| 3 radii × 12 planes × 48 phases | 1,452 / 296,714 | 1,146 | 83 | 4.15e+08 | 8.49321e-05% |
| 24 planes × 64 phases | 204,493 / 487,104 | 49 | 17 | 8.5e+11 | 0.623114% |
| 24 planes × 128 phases | 820,038 / 1,949,873 | 48 | 19 | 9.5e+11 | 0.779212% |

Every subset passing both guards reaches **zero interception** at some sampled
dates. All 10 km subsets also have zero full-Sun-protected area throughout the
sampled history. No minimum sufficient inventory has been established.

These subsets pass the pair guards on the computed histories. The phase
accuracy failure below prevents promoting their geometric margins to a
physical formation-control demonstration. The vectorized encounter backend
is tested against the original pair-by-pair calculation and preserves its
pair minima; it makes the larger scans cheaper to reproduce.

## Earth-safe reflection and control requirements

The controller diverts the ideal reflected solar cone away from Earth and
the protected lunar sphere. Its exclusion includes tile corners, the finite
Sun, target motion during light flight, **100 km above Earth's surface**, a
**1 mrad normal-pointing allowance**, and a further angular margin. The
sampled Earth margin is at least **0.0573°** after those allowances. A slew-
limited actuator and its transition paths remain to be propagated.

The ideal filter diverts **13.795%** of bolometric sunlight. At 1 AU and
50 g/m², maximum face-on photon acceleration is **25.05 µm/s²** and maximum
transverse acceleration is **9.642 µm/s²**, at a normal cone of **35.264°**.
The [control-budget calculation](results/fleet_control_budget.json) separates
these analytical limits from executed maneuvers. An ideal bounded one-axis
rest-to-rest translation takes **1.27 hours for 50 m**, **1.79 hours for 100 m**,
and **17.89 hours for 10 km**. Gravitational coupling, solar availability,
coverage and safe-beam constraints belong in an actual avoidance solve.

At the 15,000 km circular speed of **572 m/s**, 50 m corresponds to **0.0875 s
of along-track phase error**. This is a state-estimation tolerance; the control
update interval needs its own error budget. The sampled plane commands in the
10 km multi-plane pilot reach **0.150°/s**. Sampling gives a lower bound on
actuator demand. Torques, structural flexure, hardware mass and electrical
attitude-control loads remain unbudgeted.

## Numerical and physical gates

An independent replay of eight distinct phases over the whole month, with
maximum step reduced from 1,200 to 300 seconds and relative tolerance tightened
from 2e-10 to 2e-12, differs by up to **898.3 m**. **The 50 m phase/placement gate
fails.** Reported sampled coverage can characterize these candidate geometries;
metre-scale seam closure and close-approach execution remain unresolved.

A 3 × 3 surface quadrature finds omitted finite-extent gravitational
accelerations up to **2.63e-9 m/s² for 10 km squares** and **2.66e-5 m/s² for
1,000 km squares** at the checked dates. The latter is comparable to the
entire available photon acceleration. The large-square trajectories therefore
serve as geometric probes, with a substantial extended-body force correction
still outstanding. [Validation results](results/fleet_validation.json) preserve
both the failed placement gate and the force residuals.

Further work must couple mutual illumination to translational force, integrate
finite-extent forces and real angular spectra, converge orbital phase to the
placement budget, and optimize initial phases and steering with collision and
coverage constraints inside the solve. Millions of independently moving 10 km
cells require a larger-scale numerical method. These runs establish neither a
sufficient fleet inventory nor a lower bound excluding better orbital families.

## Separate power and habitat accounts

No electric translational thrust is inserted into these runs. Collector and
habitat mass, attitude/control hardware beyond the assumed base allocation,
conversion losses, power transfer, replacement inventory and delivered
electricity have no demonstrated closure here. The architecture's **100 TW**
remains a separate design case. Additional mass changes the photon acceleration
and must be included in a subsequent coupled design.

## Reproduction

Install [requirements.txt](requirements.txt), restore the checksum-pinned
DE440s kernel, and use the configurations in [fleet_search.json](fleet_search.json)
with `fleet_run`. The ignored raw histories are checked against their source,
constant, configuration and trace hashes before analysis.

```sh
python -m pip install -r research/studies/solar_shield_array/requirements.txt
python -m protection.dynamics.ephemeris --download
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name planar_10km --planes 1 --phases 576 --arrangement planar
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name planes_10km --planes 12 --phases 48
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name families_inner_10km --radii-km 15000 16000 17000 --planes 12 --phases 48
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name dense_1000km --planes 24 --phases 64 --side-km 1000
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_run --name dense_double_1000km --planes 24 --phases 128 --side-km 1000
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_analyze --names planar_10km planes_10km --output fleet_checkpoint.json
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_analyze --names families_inner_10km --output fleet_families.json
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_analyze --names dense_1000km --output fleet_dense.json
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_analyze_fast --names dense_double_1000km --output fleet_dense_double.json
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_refine --products fleet_checkpoint.json fleet_families.json fleet_dense.json --names planes_10km families_inner_10km dense_1000km
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_handover --product fleet_dense.json --name dense_1000km
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.fleet_validate --names planes_10km dense_1000km --replay planes_10km
python -m research.studies.solar_shield_array.fleet_control_budget
python -m pytest protection/dynamics research/studies/solar_shield_array -q
python visualization/solar-shield-array/fleet_plot.py
```

The full diagnostics, sampling settings and exclusions remain in the JSON
products beside this report. [checks.json](checks.json) distinguishes passing
implementation checks, the failed formation-accuracy gate and the existing
repository input/environment failures.
