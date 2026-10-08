# Orbital tacking at the existing film mass

The first numerical cycling study supports pursuing the author's proposal
without reducing the **50 g/m² base optical allocation**. A tile starting on
a 15,000-km lunar retrograde seed remains between **13,972 and 16,003 km from
the Moon's centre over three years**, with gravity and solar radiation as its
only translational forces. It returns for **583 shadow-service passes**.

The successful control is simple: face the Sun while the tile can intercept
rays headed into the four-lunar-radius protected region, then feather the
membrane toward edge-on as it leaves. Gravity supplies the return motion.
The useful result is orbital cycling with scheduled optical force. Continuous
reflection and a simple rotating-energy feedback law fail in the tested
50,000-km cases. General solar-sail controllability has not been proved.

This establishes candidate **point-tile motion**, not a completed shield.
Population, orbital-plane diversity, spatial handovers and collision-free
UV coverage remain unsolved. The existing 238-TW held-aperture result and
these cycling results solve different geometrical tasks.

The numerical product is [results/cycling.json](results/cycling.json).
[cycling_plot.py](../../../visualization/solar-shield-array/cycling_plot.py)
renders its trajectory, range, service intervals and inventory comparison.

The follow-up corrects Boolean trapezoidal integration in the raw service-time
diagnostic. NumPy's Boolean addition counted adjacent active samples once
before division by two. Casting the indicator to floating point restores the
annual inner case from **8.00394% to 15.82345%** and its three-year case from
**8.03610% to 15.88746%**. Every affected diagnostic was regenerated from the
hash-verified trajectories. The useful projected-area fractions, trajectories,
return arc and inventory comparisons are unchanged. Regression checks include
all-active, all-inactive and mixed runs on an irregular time grid.

## What was computed

Every propagated case starts on 2026-10-04 TDB and uses the pinned DE440s
ephemeris, lunar gravity, the Sun, Earth and the other planetary systems.
The frame translates with the ephemeris Moon; its measured acceleration is
subtracted. Circular Earth–Moon dynamics supplies retrograde initial guesses
only. Propagation itself uses the full ephemeris point-mass model.

The actuator is an ideal spectral specular membrane with physical surface
areal mass σ. With incoming light direction **e**, surface normal **n** and
cone angle α between them,

\[
\mathbf a_{\rm sail}=\frac{2fS}{c\sigma}\cos^2\alpha\,\mathbf n,
\qquad \mathbf n\cdot\mathbf e=\cos\alpha\geq0.
\]

Its projected intercepting area is physical area times cos α. This is a
specific passive mirror law, narrower than the relaxed momentum ball used
in the held-aperture study. No electrical correction force is inserted into
the propagated trajectory.

Two spectral settings are compared. **Filter spectrum** reflects the
13.795% bolometric allocation permitted by the stored climate setup and
transmits the rest everywhere. **Full outside core** permits full reflection
when none of the incident solar ray bundle can reach the solid Moon; it uses
the same 13.795% limit within that central bundle. These are ideal spectral
actuators. The real TiO₂/PV/diffractive stack's angular transmission,
absorption and radiation force have not been fitted to this law.

The filter-spectrum experiment leaves much of the off-core light uncollected.
That is an optical-routing choice, not a reduction in film mass. Power
collectors and habitat hubs need their own orbital and momentum accounts;
this run does not supply the architecture's delivered 100 TW.

In the service-and-coast law, the sail stays face-on through the protected
ray corridor and rotates smoothly toward edge-on over another 1,000 km.
The transverse normal is chosen so reflected light turns away from the Moon.
The finite solar disk, both Earth and Moon occultations, retarded source and
blocker positions, and lunar motion during tile-to-Moon light travel are
included. Overlapping Earth/Moon silhouettes are counted as a union.
Reflected-ray clearance includes a target-motion margin. Earth-safe outgoing
beam directions have not been imposed.

The coverage statistic follows an infinitesimal membrane patch's incident
ray bundle. The fraction of the finite solar disk whose rays enter the
protected region is weighted by projected membrane area and available
sunlight, then averaged in time. It is a **useful projected-area fraction**,
not the fraction of the Moon protected by one tile or a demonstrated fleet.
Finite tile shapes, retarded moving-tile intersections and seams are absent.

## Propagated comparisons

Distances in this table are from the Moon's centre; the seed radius is not
held fixed. The 15,000-km family extends the previous held-screen search
inward while staying well outside the 6,950-km protected radius.

| Seed and optical setting | Duration | Lunar distance range | Useful projected-area fraction | Result |
|---|---:|---:|---:|---|
| 15,000 km; filter spectrum; service/coast | 1 year | 14,027–15,923 km | **15.687%** | Bounded throughout |
| Same case, longer propagation | 3 years | 13,972–16,003 km | **15.737%** | Bounded throughout |
| 30,000 km; filter spectrum; service/coast | 1 year | 27,427–32,292 km | 7.401% | Bounded throughout |
| 30,000 km; full outside core; service/coast | 1 year | 20,778–44,912 km | 7.791% | Bounded throughout, wider excursions |
| 50,000 km; full outside core; service/coast | About 211 days | Exceeds 300,000-km guard | — | Rejected |
| 50,000 km; continuous reflection with lunar beam guard | About 14.4 days | Crosses inner protection guard | — | Rejected |
| 50,000 km; rotating-energy feedback | About 28 days | Exceeds 300,000-km guard | — | Rejected |

The inner case makes **194 service visits in its first year**, about
**7.07 hours per visit**, with a maximum between-visit gap of **1.59 days**.
Its commanded membrane-plane rate during deployment/feathering reaches
about **0.067°/s** on the sampled history. That is a command requirement;
no attitude actuator, torque authority, flexible structure or attitude
hardware mass has been demonstrated. Eclipse force automatically falls with
available sunlight; no blackout electric-thrust substitute is used.

The three-year run supplies a finite-horizon boundedness test. It is not a
stability theorem or a demonstration over the 19-year nodal cycle, much less
the architecture's lifetime. Lunar/Earth multipoles, relativistic corrections,
real optical response, environmental forces and hardware loads remain outside
this point-mass experiment.

## A separate optimized return arc

A Hermite–Simpson collocation search varies free states and passive mirror
angles together. Its 50,000-km retrograde seed yields an **8.90548-day arc**
which returns to the same position and velocity in the instantaneous
Earth–Moon frames at its endpoints. The ephemeris is not made periodic.

Independent DOP853 propagation of the stored controls closes within
**16.39 m and 0.165 mm/s**, with a **2.031%** useful projected-area fraction.
The same initial state and elapsed time with photon acceleration removed
misses that return by **4,527 km and 16.33 m/s**. This gravity-only
counterfactual isolates the contribution of the sail schedule; it is not an
alternative filtering actuator.

Refining from 160 to 320 collocation intervals reduced the independently
propagated position discrepancy from about 171 m to 16 m. The optimizer hit
its evaluation budget. Acceptance rests on the separate 50-m / 1-mm/s replay
criteria and passive-force checks, not its success flag. The maximum scaled
dynamics defect is 1.89×10⁻⁹; the soft seed-position penalty is reported
separately. Initial state and every control knot are committed, so replay
does not require the ignored optimizer archive.

This optimized arc and the multi-year service-and-coast family are distinct
experiments. The former establishes an admissible finite return with active
steering. The latter demonstrates repeated service without an imposed
endpoint reset.

## Inventory implications

For a common comparison, divide the existing held aperture's physical base
area by each measured useful-area fraction. That reference aperture is
1.746×10¹⁴ m² and its unloaded optical mass is 8.73×10¹² kg.

| Annual candidate | Equivalent area multiple | Base optical inventory | Same mass in reference propellant years |
|---|---:|---:|---:|
| 15,000 km, filter spectrum | **6.37×** | **5.56×10¹³ kg** | **4.76** |
| 30,000 km, filter spectrum | 13.51× | 1.18×10¹⁴ kg | 10.08 |
| 30,000 km, full outside core | 12.84× | 1.12×10¹⁴ kg | 9.57 |

These are **reference-area inventory scales under perfect placement**.
The deployed physical area and spatial efficiency of a real cycling fleet
have not been solved. The reference aperture has a different solar distance
and finite-Sun margin, so this table is a comparison convention rather than
an exact minimum-mass fleet calculation. The last column compares mass with
the held case's 1.17×10¹³ kg/year of expended propellant. It is not an economic
payback or a forecast that the same fleet survives without replacements.
Attitude equipment, collectors, habitats, deployment, overlaps, reserves and
replacement mass are excluded.

The potential advantage is substantial: reusable orbital inventory can
replace the requirement that every tile continuously resist falling away
from a fixed sunward aperture. The calculation has demonstrated the tile
motion needed to investigate that trade at 50 g/m². It has not shown that
6.37 reference apertures are sufficient for continuous protection.

## Numerical limits and fleet coverage

The annual case was repeated with the maximum integration step reduced from
7,200 to 1,800 seconds, a tenfold tighter relative tolerance and twice as many
saved states. After evaluating both on the same 300-second ray grid, useful
area differs by **0.00113%**. Refining the ray grid from 300 to 120 seconds
changes that statistic by **0.0277%**. The photon force remains in the passive
hemisphere; the central spectral budget is preserved. The inner case's
reflected solar cone clears the protected lunar sphere by over four degrees.

Absolute orbital phase is less well resolved: the annual trajectories differ
by up to **5.70 km** between those integration settings. A separate initial
perturbation of 1 m in ICRF x and 1 mm/s in ICRF y grows to **47.7 km** of
position separation, while the orbit remains bounded and useful area changes
by only **0.0112%**. These are adequate checks of the reported area statistic
and large orbital-clearance margins. They do not establish metre-scale
formation control, predictable seam positions or collision-safe handovers.

Adding phase copies of one path cannot establish full coverage. In the
exactly planar limit, let **n** be the orbit-plane normal, **u** the sunward
direction, and **q = b + d u** a sunward tile position with d > 0. Then
**n·b = −d n·u**: projected tile centres from that plane occupy only one side
of a line through the target, or the line itself when the Sun lies in the
orbital plane. Finite-Sun broadening does not fill the four-radius disk.
Several orbital planes or other spatially distributed families are needed.

The next acceptance test must solve a three-dimensional population with
actual phases propagated in the time-dependent ephemeris, and require every
necessary UV ray to meet filtering material at every time. It must include
finite tiles, partial eclipse service, seams, collision avoidance and
Earth-safe reflected/routed beams. Simultaneous phase shifts cannot be
inferred by time-shifting one trajectory. Collector/habitat loads and
achievable optical/attitude control belong in that same closure.

## Reproduction and checks

The ephemeris binary stays outside Git and is verified against its recorded
size and SHA-256. The existing held-aperture products and historical optical
model were not changed.

```sh
python -m protection.dynamics.ephemeris --download
# Fast independent replay directly from the committed control schedule:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.cycling_replay
# All configured propagation cases; one worker by default, two at most:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.cycling_run --workers 2
# Recreate and refine the independently replayed return arc:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.cycling_search --radius-km 50000 --phase .5 --intervals 160 --evaluations 350
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.cycling_search --radius-km 50000 --phase .5 --intervals 320 --evaluations 300 --resume research/runs/solar_shield_array/cycling/candidate_50000_0_160_0.5.npz
python -m research.studies.solar_shield_array.cycling_analysis --arc research/runs/solar_shield_array/cycling/candidate_50000_0_320_0.5.npz
python -m pytest protection/dynamics research/studies/solar_shield_array -q
python visualization/solar-shield-array/cycling_plot.py
```

Raw products record configuration, dynamics source hashes, named constant
values and trace hashes. The exporter rejects changed sources, constants,
configuration or raw bytes. Model tests cover conserved reflected momentum,
Earth/Moon occultation unions, moving frames, circular-problem Jacobi
conservation, the retrograde return seed, shadow-service geometry and safe
feathering. Numerical trajectories are simulation evidence, not experimental
validation. Repository checks and their pre-existing blockers are recorded
in [checks.json](checks.json).
