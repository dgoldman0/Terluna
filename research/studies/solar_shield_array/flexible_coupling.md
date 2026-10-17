# Flexible array screening and coupled solver

14 October 2026 · `research/array-industry` · baseline `669ddd3`

The cheap screen is complete, and a reduced model couples membrane deformation, photon loading, attitude, orbital
motion and bounded control. Its numerical verification and the fixed-case solver benchmarks are recorded apart from
physical validation, and **the coupled design search comes next**. The rigid-tile flight, optical and support products
are unchanged.

## What the broader screen finds

The [screen product](results/geometry_screen.json) covers 2,700 mechanical designs and 16,200 spacing configurations.
The sides span 1–10 km, the edge retreats 0.5–10% of the side and the boom depths 0.3–1.5%, with three
compression-material scenarios, two tensile-cable scenarios and three prescribed load envelopes. Each mechanical design
is reused at shingles of 1.5°, 2° and 3° and two pitch factors: an algebraic enumeration, ahead of a propagated design
search.

Of the 2,700 mechanical cases, 2,650 converge; 30 diverge in sizing their own load and 20 reach the iteration limit,
and every failure is kept. Of the 16,200 spacing configurations, 2,404 pass the declared member-domain and 150 m
parallel-clearance screens. **None has even its sheet-plus-support subtotal within the inherited 13.765 g/reference m²
complete-tile allowance.**

| Prescribed normal load | Lightest passing subtotal, g/reference m² | Nominal clearance, m |
|---|---:|---:|
| Quarter of inherited inertial/gravity load | 21.472 | 159.156 |
| Half of inherited inertial/gravity load | 22.178 | 159.156 |
| Inherited inertial/gravity load | 23.590 | 159.156 |

All three use a 5 km side, a 125 m edge retreat, 37.5 m boom depth, the M55J compression proxy, the Spectra tension
proxy, a 3° shingle and the larger pitch. Their sufficient coverage construction needs about 122.395 million units,
4.335 times the earlier inventory. The subtotals leave out the control, power, rotor and service hardware; with the
inherited bus on every unit, housekeeping alone is 1.224 TW and the avionics inventory 0.612 Gt, and the alternative
area-scaled bus is a hypothesis.

The shingle angle sets the result. Within a ring the parallel shingles stand at the pitch times the sine of the
shingle, so the 150 m gate leaves the two supports between them only what is left over. By shingle angle, from the
screen's archived configurations (quarter to full inherited load):

| Shingle | Lightest passing subtotal, g/reference m² | Its tiles | Passing configurations |
|---|---:|---|---:|
| 1.5°, the flown design | 138.6–345.0 | 10 km, 2.5% retreat, 0.3% boom depth | 24 of 5,400 |
| 2° | 29.0–37.2 | 10 km, 7.5% retreat, 0.75% depth | 518 of 5,400 |
| 3° | 21.5–23.6 | 5 km, 2.5% retreat, 0.75% depth | 1,862 of 5,400 |

A 3° shingle with 5 km tiles keeps the 222 m between shingles of the present 10 km tiles at 1.5° with half the span.
It departs from the flown design: the attitude study's clearance envelope, shadow torques, momentum store and push are
all for 1.5°.

Without the clearance screen the smallest subtotal is 17.016 g/reference m², with a nominal clearance of only 6.973 m
and about 3.060 billion units, so lower component mass alone leaves the inherited formation open. The quarter- and
half-load cases are load hypotheses, ahead of an attitude law that achieves them, and the sunlight's pressure stays
unchanged in them. The screen is conservative, and its failures belong to the tested supports and attitudes.

The calculation reads the inherited film mass, dimensions, lattice, prestress assumptions and safety gate. It iterates
the acceleration envelope to contain finite tubes, braces and cables, with prescribed reserves for sag and member
deflection. An inscribed rectangle gives a sufficient nominal planar overlap on the undeformed geometry, with no
subtraction for the loaded anchors' retreat; a full-fleet optical and swept-clearance certificate needs more. Opaque
hardware, joints, neighbours out of parallel and service motion need geometry of their own. The material evidence keeps
the limits in [support_materials.json](support_materials.json) and [support_sources.json](support_sources.json).

The run took 310.34 s on one thread with a 163.16 MiB high-water mark; it was regenerated on 14 October 2026 after the
source check corrected the wording of its registers, with every record unchanged (the first run took 297.58 s). All
18,900 mechanical and spacing records, 34.78 MB, are kept in ignored JSONL archives under
`research/runs/solar_shield_array_geometry_screen/`, whose hashes, distributions, conditional Pareto sets and
representatives are in the committed product. Every case, failures included, is reproducible.

## Structural and dynamical model

[array_flexible.py](../../../engineering/array_flexible.py) builds a small-slope model of normal displacement on a
scalloped triangular membrane. Its energy includes the prescribed membrane prestress, tensioned edge cables, finite
compliance at the boom tips and an assumed twist restraint at the hub. Membrane and cable displacements are measured
from the hub's affine best-fit plane, which leaves three exact rigid normal modes, and the elastic modes are
mass-normalized and orthogonal to rigid translation and tilt.

The physics is narrower than a nonlinear finite-element sail model's: in-plane strain, wrinkling, thermal strain,
prestress loss, joint slip, local failure of the tube shells and deployment are outside it. The hub's restraint is an
uncalibrated effective stiffness. The sheet's mass matches the analytic curved area, while the support's mass is lumped
at member endpoints, and the product reports the inertia this gives against the earlier distributed-frame
calculation: a property of the model, which refining the membrane mesh leaves in place.

The structural result crosses into dynamics as a versioned, hash-checked `terluna.engineering.array-flexible/1`
product holding nodes, triangles, lumped masses, normal modes, modal frequencies, evidence and its reading rule, so the
dynamics read another lane's result through a product. The
[method and source notes](../../../engineering/flexible_sources.json) record how far each source was read; the
sources give the model its method context.

[flexible_fleet.py](../../../protection/dynamics/flexible_fleet.py) places the nodes at `x + R(r + Φq)`, propagating
the centre of mass, the quaternion, the body angular velocity and the modal coordinates. Its rotational and modal mass
matrix is

\[
\begin{bmatrix}I(q)+I_\mathrm{rotor}&G(q)\\G(q)^T&M_q\end{bmatrix},
\quad G_j=\sum_i m_i\,r_i\times\Phi_{ij}.
\]

For this product's normal-only modes on a flat reference G vanishes, since every r × Φ lies in the plane and the modes
are orthogonal to the rigid tilts; the rigid and elastic motions couple through the velocity terms and the forces, and
the random three-dimensional modes of the energy test exercise G itself. The acceleration solve keeps the centrifugal
and Coriolis terms `ω × (ω × r_i) + 2 ω × (Φ_i q̇)`. Node forces feed translation, torque and the generalized modal
forces from the same geometry. Lunar gravity acts at each node. Optional time-dependent point-mass perturbers subtract
the acceleration of the Moon origin, and callbacks for the source and for spherical occultation allow replays driven by
an ephemeris without downloading or synthesizing inputs. The fixed benchmarks state which environmental terms they use.

The diagnostics integrate external work and angular impulse in the chosen nonrotating Moon-centred frame, including
its origin-acceleration subtraction, which checks conservation with photons, damping, moving third bodies and motors.
Energy includes the supplied linear elastic potential. Agreement establishes the equations of this declared reduced
system; its structural realism is a separate question.

## Optical feedback and controls

[flexible_optics.py](../../../protection/dynamics/flexible_optics.py) recomputes face normals, areas, mutual shadowing
and multilayer transmission from the current deformed geometry. Absorption and specular reflection carry momentum and
straight transmission none. Each face force is shared among its three nodes, preserving its force, torque and virtual
work. Positive finite-Sun quadrature and the receiver-to-source rays use the same triangles.

The numerical approximations are explicit: illumination at face centroids, sampled solar directions, grey coefficients
and straight transmitted rays. Diffuse and secondary reflections, thermal recoil and unmeshed hardware surfaces are
outside the model. Chords along the inward-curved scallops lie outside the true curved boundary, so coarse meshes
overstate the edge area and can overstate coverage, a bias the exact normalization of the sheet's mass leaves in
place. Broadband attenuation and UV transmission are different quantities, the second needing the band's own
coefficients. The [optical source notes](flexible_optics_sources.json) give the primary sources and how far they were
read.

Four-sector reflection trim allocates a bounded change of reflectivity over a declared fraction of the film. Its
response is instantaneous; the actuators' bandwidth, switching electricity and thermal transients come later, and the
trim supplies a bounded set of wrenches. Partial boundary triangles take an area-homogenized trim fraction. The local
allocation holds the upstream transmission fixed, then the simultaneous commands trigger a full recomputation of the
illumination, with the residual force and torque and the allocation's status reported. Formation and attitude requests
come as controller callbacks, so the mechanics carry no imposed target motion.

Optional three-axis momentum wheels have finite motor torque and relative momentum capacity. Their canonical momentum
includes the body's rotation, and motor power comes from the relative wheel speed. A rotor's point mass belongs in the
structural product, and its own inertia is added explicitly. The wheels attach to the ideal mean reference frame, so
the local deformation of the hub attachment and the finite rotors' gravity gradients are outside the model. The tests
verify the exchange equations; a wheel design closed in mass and power, with periodic unloading, is a later step. The
fixed tile benchmarks use photon trim, with no added wheel hardware.

The integrator stops at the prescribed small-slope limit or at wheel capacity and keeps the terminal crossing state.
The sampled clearance is a lower bound on separation from axis-aligned bounding boxes, less the declared hardware
envelopes: a positive value establishes separation at that instant, and zero is inconclusive. Contact between outputs
needs continuous checks beyond successful integration and positive sampled coverage, and the model has no contact
response.

## Verification and benchmarking

The [check record](flexible_checks.json) separates executed tests, fixed benchmarks and the work still to do. The tests
cover independent square-membrane eigenfrequencies and static sag, mesh and modal refinement, affine rigid null modes,
mass and work identities, the deformed optical response, analytic reflection, multilayer rays, finite-source
convergence, free oscillation, the coupled free-body invariants, node gravity, motor exchange, saturation and the
terminal limits. An independent uncondensed kinetic-energy calculation recovers the coupled acceleration equations by
numerical differentiation, with general three-dimensional test modes and momentum wheels. Some tests are weaker than
their names: the internal-force check of `test_array_flexible.py` holds for any symmetric stiffness, the work ordering
holds by construction for nested modes, and the energy-minus-work closures of `test_flexible_fleet.py` check the
bookkeeping for any force law, the photon and tidal signs being tested in `test_flexible_optics.py` and
`test_flexible_environment.py`.

The final focused suite passes **119 tests**, and the required `make check` passed with **2,253 Python tests passed and
16 skipped, 191 JavaScript tests passed and 975 files through the layer check**, followed by successful provenance and
ensemble integrity checks; the Python suite took 452.62 s. These are numerical checks: scholarly source review,
hardware qualification and manuscript clearance are separate, and the exact logs and file hashes are in the check
record.

Benchmark timings and resolution changes are in [flexible_benchmark.json](results/flexible_benchmark.json). The
support case was fixed before the screen's results: the inherited 10 km M55J and Kevlar case with a 500 m edge retreat
and 100 m depth, kept as a diagnostic case. The mesh, modal, time-step, source and
integration-method comparisons keep the same declared physical inputs where they apply, and compare physical probes
since modal coordinates differ between bases. The product records CPU and wall time, the process's high-water memory,
solver calls, validity stops, the source, model and input hashes and the limits of the environment.

All 24 fixed cases completed. Building the ten structural products and running the cases took **125.25 s of wall time,
125.11 s of CPU and 125.82 MiB of peak process memory** on one numerical thread of an Intel i9-11980HK, with Python
3.13.7, NumPy 2.3.3 and SciPy 1.16.2. That is the run of 14 October 2026 regenerated after the source check, with a
climate model on four other cores; the first run, on a quieter machine, took 102.53 s of wall time and 102.47 s of CPU,
and every result other than the timings is unchanged. The benchmarks use lunar point gravity and a fixed solar
direction, and they leave the ephemeris, differential-gravity and eclipse helpers, which are tested separately, for
longer flights.

The modal check found a consequential limitation: **none of the 6, 12, 20, 40 or 80-mode approximations tested meets
the declared shape criterion** against the complete 170-mode normal basis on the 12-division mesh. Their largest
terminal probe differences are 18.806, 28.217, 18.212, 3.461 and 2.043 mm, against a criterion of 0.1 mm plus 1% of
the reference's largest terminal probe displacement, about 0.719 mm here: a numerical criterion chosen for the
benchmark. Part of the difference is the projection of the initial bump onto the
truncated basis: for six modes it is 1.82 mm at the first output, 0.81 mm at 20 s and 18.8 mm by 60 s. Mode
counts of 81–169 were not tested.

Complete discrete normal bases remove this truncation ambiguity. The 8-to-12-division mesh comparison changes the
terminal shape probes by 0.741 mm, just above its criterion, and the 12-to-16 comparison by **0.176 mm**, below its
0.720 mm criterion. Cutting the largest integration step from 10 to 0.5 s on the complete 16-division model changes the
terminal probes by about 1.2e-13 m. These results belong to the declared normal-displacement model and the fixed 60 s
case. The largest sampled slope grows by 25.35% between the complete 12- and 16-division meshes, from 6.46e-5 to
8.10e-5, so the displacement criterion is met while the surface normals and local slopes are still resolving. The
longer 621.74 s, trim, finite-source and integration-method cases keep six modes, so their shapes carry the truncation
above.

A separate read-only audit rebuilds the same nine probes from the saved states at every shared 5 s output. The complete
12-to-16 mesh difference peaks at **0.188 mm at 35 s**, below the 0.720 mm criterion, and time refinement stays
negligible; the 80-mode error peaks at 2.410 mm at 50 s. The audit's script, input hashes and results are recorded in
[flexible_checks.json](flexible_checks.json); it extends the comparison to every saved time, while unsampled times and
membrane locations stay outside it.

Integration costs in the regenerated run, without diagnostics and archive writing:

| Fixed case | Simulated time | Integration wall time | Actual RHS calls |
|---|---:|---:|---:|
| One tile, 6 modes, DOP853 | 30 s | 0.578 s | 167 |
| Same case, RK45 | 30 s | 0.242 s | 62 |
| Same case, Radau | 30 s | 0.584 s | 157 |
| Three tiles, 6 modes each, DOP853 | 30 s | 2.360 s | 167 |
| One tile, complete 12-division basis, DOP853 | 60 s | 3.223 s | 239 |
| One tile, complete 16-division basis, DOP853 | 60 s | 12.177 s | 317 |
| Same complete basis, maximum step 0.5 s | 60 s | 58.412 s | 1,529 |

RK45 and Radau agree with DOP853's terminal shape probes within 1.4e-11 m on that short six-mode case. Radau reports
81 function evaluations, and the instrumentation counts 157 including its numerical Jacobian. RK45 is fastest in this
comparison; a ranking of methods at complete modal resolution or over a service orbit is still to measure, and the
cheap six-mode timings understate the cost of a resolved structural calculation.

A first Radau run exposed a guard that rejected its Jacobian perturbations of inactive wheel coordinates. The fix
ignores those coordinates in the equations and diagnostics while still rejecting a nonzero initial momentum of an
inactive wheel; a regression compares implicit and explicit propagation, and the first failure and its source snapshot
stay in the ignored run archive. A startup check replaces one absolute tolerance across quantities of different units
with explicit tolerances for position, velocity, attitude, deformation, work and impulse, which cuts the actual RHS
calls from 491 to 95 over 2 s with a terminal position difference of 3.73 nm. The product records every tolerance and
the differences in work and impulse; the speedup belongs to this case.

The unforced six-mode case drifts in energy by at most 4.55e-13 J over 60 s. The separate analytic and independent
kinetic-energy tests matter here, since tiny relative errors against a large orbital energy can hide structural errors.
The bounded trim allocator completes and leaves sampled residuals of up to **3.213 N and 23.332 kN·m**, so its control
authority is short of the request, and full coverage at the 25 chosen receiver probes leaves the seams and complete
optical coverage to test.

The endpoint-lumped support inertia differs from the earlier distributed-frame reference by **24.83%, 24.19% and
23.97%** in relative Frobenius norm on the three meshes; on the finest mesh the principal inertias are about 24% too
large. Lumping a uniform member at its ends overstates its second moment, so a distributed support mass (for example
three-point lumping along each boom) or a stated sensitivity is needed before the control and torque predictions are
used quantitatively, along with the assumed hub and boom compliance. The conservation checks verify the declared lumped
model. The overestimate of the optical polygon's area falls from 0.243% to 0.061% across the meshes, and refining the
mesh leaves the support-inertia discrepancy where it is.

## Reproduce and continue later

Run from the repository root with one numerical thread:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.geometry_screen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.flexible_benchmark
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest engineering/test_array_flexible.py protection/dynamics/test_flexible_fleet.py protection/dynamics/test_flexible_fleet_energy.py protection/dynamics/test_flexible_optics.py protection/dynamics/test_flexible_environment.py research/studies/solar_shield_array/test_geometry_screen.py research/studies/solar_shield_array/test_flexible_benchmark.py
```

For the search, rebuild the structural product whenever the geometry, material or attached mass changes, then vary
spacing, initial orbit, attitude request and finite controls together. Keep the failed trajectories and budget each
resolution by the measured cost per right-hand-side evaluation. The three-tile timings scale with the shadow
candidates, the mesh and the source samples, so a dense full fleet costs more than a linear extrapolation of them.

Before a numerical optimum counts as feasible, the spatial and temporal sampling of shadows and clearance is refined,
the relevant ephemeris and attitude cycles are propagated, the hardware's optical surfaces and loads are counted, and
the support and attachment constitutive assumptions are replaced or bounded. Qualification, trials of the full fleet,
its faults and service, and the closure of manufacturing and resources lie beyond this numerical model. The work
supplies a tested reduced-model starting point and measured compute costs.
