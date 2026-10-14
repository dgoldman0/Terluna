# Flexible array screening and coupled solver

14 October 2026 · `research/array-industry` · baseline `669ddd3`

The cheap screen is complete. A new reduced model couples membrane deformation,
photon loading, attitude, orbital motion and bounded control; its numerical
verification and fixed-case solver benchmarks are recorded separately from
physical validation. **The coupled design search has not been run.** The old
rigid-tile flight, optical and support products remain unchanged.

## What the broader screen finds

The [screen product](results/geometry_screen.json) covers 2,700 mechanical
designs and 16,200 spacing configurations. Side lengths span 1–10 km, edge
retreats 0.5–10% of side, boom depths 0.3–1.5%, three compression-material
scenarios, two tensile-cable scenarios and three prescribed load envelopes.
Each mechanical design is reused at 1.5°, 2° and 3° shingle and two pitch
factors. This is an algebraic enumeration, not a propagated design search.

Of 2,700 mechanical cases, 2,650 converge; 30 diverge in own-load sizing and
20 reach its iteration limit. All failures are retained. Of 16,200 spacing
configurations, 2,404 pass the declared member-domain and 150 m parallel
clearance screens. **None has even its sheet-plus-support subtotal within the
inherited 13.765 g/reference m² complete-tile allowance.**

| Prescribed normal load | Lightest passing subtotal, g/reference m² | Nominal clearance, m |
|---|---:|---:|
| Quarter of inherited inertial/gravity load | 21.472 | 159.156 |
| Half of inherited inertial/gravity load | 22.178 | 159.156 |
| Inherited inertial/gravity load | 23.590 | 159.156 |

All three use a 5 km side, 125 m edge retreat, 37.5 m boom depth, M55J
compression proxy, Spectra tension proxy, 3° shingle and the larger pitch.
Their sufficient coverage construction requires about 122.395 million units,
4.335 times the earlier inventory. Those subtotals exclude complete control,
power, rotor and service hardware. At the inherited per-unit bus assumptions,
housekeeping alone is 1.224 TW and avionics inventory is 0.612 Gt; the
alternative area-scaled bus is expressly hypothetical.

Without the clearance screen, the smallest subtotal is 17.016 g/reference m²,
but its nominal clearance is only 6.973 m and inventory about 3.060 billion
units. Lower component mass therefore does not close the inherited formation.
The quarter/half-load cases do not demonstrate an attitude law that achieves
those loads; sunlight pressure remains unchanged. Failure of this conservative
screen does not prove every support architecture or attainable attitude fails.

The calculation reads the inherited film mass, dimensions, lattice, prestress
assumptions and safety gate. It iterates the acceleration envelope to contain
finite tubes, braces and cables, plus prescribed sag and member-deflection
reserves. An inscribed rectangle provides a sufficient nominal, undeformed planar
overlap construction; the optical core does not subtract loaded anchor retreat. This is not a full-fleet optical or swept-clearance certificate.
Opaque hardware, joints, nonparallel neighbors and service motion still need
their own geometry. Material evidence retains the limits in
[support_materials.json](support_materials.json) and
[support_sources.json](support_sources.json).

The one-thread run took 297.58 s with a 163.28 MiB process high-water mark.
All 18,900 mechanical and spacing records occupy 34.78 MB in ignored JSONL
archives under `research/runs/solar_shield_array_geometry_screen/`. Their
hashes, distributions, conditional Pareto sets and representatives are in the
committed product. Every case, including failure, is reproducible.

## Structural and dynamical model

[array_flexible.py](../../../engineering/array_flexible.py) constructs a
small-slope normal-displacement model on a scalloped triangular membrane.
Its energy includes prescribed membrane prestress, tensioned edge cables,
finite boom tip compliance and an assumed hub twist restraint. Membrane and
cable displacements are relative to the hub's affine best-fit plane, leaving
three exact normal rigid modes. Elastic eigenmodes are mass normalized and
orthogonal to rigid translation and tilt.

The physical assumptions are narrower than a nonlinear sail finite-element
model: no evolving in-plane strain, wrinkling, thermal strain, prestress loss,
joint slip, local tube-shell failure or deployment. Hub restraint is an
uncalibrated effective stiffness. Sheet mass matches the analytic curved area;
support mass is lumped at endpoints. The product reports the resulting inertia
discrepancy against the earlier distributed-frame calculation. This is a
model approximation, not numerical error removable by membrane mesh refinement.

The structural result crosses into dynamics as a versioned, hash-checked
`terluna.engineering.array-flexible/1` product, containing nodes, triangles,
lumped masses, normal modes, modal frequencies and evidence. It is not a
silent import of another lane's result. [Method and source notes](../../../engineering/flexible_sources.json)
record the specific reading extent and distinguish context from validation.

[flexible_fleet.py](../../../protection/dynamics/flexible_fleet.py) uses
world positions `x + R(r + Φq)`, with propagated centre of mass, quaternion,
body angular velocity and modal coordinates. Its rotational/modal mass matrix is

\[
\begin{bmatrix}I(q)+I_\mathrm{rotor}&G(q)\\G(q)^T&M_q\end{bmatrix},
\quad G_j=\sum_i m_i\,r_i\times\Phi_{ij}.
\]

The acceleration solve retains centrifugal and Coriolis terms
`ω × (ω × r_i) + 2 ω × (Φ_i q̇)`. Node forces feed translation, torque and
generalized modal forces from the same geometry. Lunar gravity acts at each
node. Optional time-dependent point-mass perturbers subtract the acceleration
of the Moon origin; source and spherical-occultation callbacks support
ephemeris-driven replay without silently downloading or synthesizing inputs.
The fixed benchmarks state which environmental terms they actually use.

Diagnostics integrate external work and angular impulse in the chosen
nonrotating Moon-centred frame, including its origin-acceleration subtraction. These provide checks
on conservation even with photons, damping, moving third bodies and motors.
Energy includes the supplied linear elastic potential. Agreement establishes
the equations of this declared reduced system, not its structural realism.

## Optical feedback and controls

[flexible_optics.py](../../../protection/dynamics/flexible_optics.py) recomputes
face normals, areas, mutual shadowing and multilayer transmission from the
current deformed geometry. Absorption and specular reflection supply momentum;
straight transmission supplies none. Each face force is distributed to its
three nodes, preserving its force, torque and virtual work. Positive finite-Sun
quadrature and receiver-to-source rays use the same triangle geometry.

The numerical approximations are explicit: face-centroid illumination,
sampled solar directions, grey coefficients and straight transmitted rays.
Diffuse/secondary reflections, thermal recoil and unmeshed hardware surfaces
are absent. Chords along inward-curved scallops lie outside the true curved
boundary, so coarse meshes overestimate edge area and can overestimate
coverage. Exact sheet-mass normalization does not remove this optical bias.
Broadband attenuation is not a UV transmission certificate.
[Optical source notes](flexible_optics_sources.json) state the primary sources
and the limited reading extent.

Four-sector reflection trim allocates a bounded change in reflectivity over
a declared fraction of the film. Response is instantaneous; actuator bandwidth,
switching electricity and thermal transients are not yet modeled. It cannot supply an arbitrary wrench.
Partial boundary triangles use an area-homogenized trim fraction. The local
allocation holds upstream transmission fixed, then the actual simultaneous
commands trigger a full illumination recomputation; residual force/torque and
allocation status remain visible. Formation and attitude requests are supplied
as controller callbacks, so no target motion is imposed on the mechanics.

Optional three-axis momentum wheels have finite motor torque and relative
momentum capacity. Their canonical momentum includes body rotation; motor
power is computed from relative wheel speed. Rotor point mass must already
be in the structural product. Additional rotor self inertia is explicit.
Wheels attach to the ideal mean reference frame: local hub attachment
deformation and finite-rotor gravity gradients are not resolved. These tests
verify exchange equations, not a mass-closed, power-closed wheel design or
periodic momentum unloading. The fixed tile benchmarks use photon trim and
do not invent wheel hardware to make their trajectories pass.

The integrator stops when the prescribed small-slope limit or wheel capacity
is reached and preserves the terminal crossing state. Sampled clearance is an
AABB separation lower bound minus declared hardware envelopes; a positive
value can establish separation at that instant, while zero is inconclusive.
Neither successful integration nor positive sampled coverage certifies no
contact between outputs. No contact response is implemented.

## Verification and benchmarking

The [check record](flexible_checks.json) distinguishes executed tests, fixed
benchmarks and work that remains unperformed. Tests cover independent square
membrane eigenfrequencies and static sag, mesh/modal refinement, affine rigid
null modes, mass/work identities, actual deformed optical response, analytic
reflection, multilayer rays, finite-source convergence, free oscillation,
coupled free-body invariants, node gravity, motor exchange, saturation and
terminal limits. An independent uncondensed kinetic-energy calculation
recovers the coupled acceleration equations by numerical differentiation,
including general three-dimensional test modes and momentum wheels.

The final focused suite passes **119 tests**. The required `make check` also
passes: **2,253 Python tests passed, 16 skipped; 191 JavaScript tests passed;
975 files passed the layer check**, followed by successful provenance and
ensemble integrity checks. The Python suite took 452.62 s. These checks do
not constitute scholarly source review, hardware qualification or manuscript
clearance; their exact logs and file hashes are in the check record.

Benchmark timings and resolution changes are in
[flexible_benchmark.json](results/flexible_benchmark.json). Their support case is
fixed before the screen's results: the inherited 10 km M55J/Kevlar case with
500 m edge retreat and 100 m depth. It is a diagnostic case, not a selected
feasible architecture. Mesh, modal, time-step, source and integration-method
comparisons retain the same declared physical inputs where applicable;
physical probes rather than incompatible modal coordinates are compared.
The product records CPU/wall time, process high-water memory, solver calls,
validity stops, source/model/input hashes and environment limitations.

All 24 fixed cases completed. Building the ten structural products and running
the cases took **102.53 s wall time, 102.47 s CPU and 126.20 MiB peak process
memory**, using one numerical thread on an Intel i9-11980HK. The environment
was Python 3.13.7, NumPy 2.3.3 and SciPy 1.16.2. These benchmarks use lunar
point gravity and a fixed solar direction; the ephemeris, differential-gravity
and eclipse helpers are separately tested but not exercised in these flights.

The modal check found a consequential limitation: **none of the tested
6, 12, 20, 40 or 80-mode approximations meets the declared shape criterion**
against the complete 170-mode normal basis on the 12-division mesh.
Their largest terminal probe differences are respectively 18.806, 28.217,
18.212, 3.461 and 2.043 mm. The comparison criterion is 0.1 mm plus 1% of
the reference's maximum terminal probe displacement, about 0.719 mm here.
This is a chosen numerical criterion, not a physical requirement or a
certified bound; mode counts 81–169 were not tested.

Using complete discrete normal bases removes this truncation ambiguity.
The 8-to-12-division mesh comparison changes the terminal shape probes by
0.741 mm, narrowly missing its criterion. The 12-to-16 comparison changes
them by **0.176 mm**, below its 0.720 mm criterion. Reducing the maximum
integration step from 10 to 0.5 s on the complete 16-division model changes
the terminal probes by about 1.2e-13 m. These results concern the declared
normal-displacement model and the fixed 60 s case; they are not a general
trajectory or physical-structure accuracy guarantee. In particular, the
maximum sampled slope increases by 25.35% between the complete 12- and
16-division meshes, from 6.46e-5 to 8.10e-5. Passing the displacement criterion
does not establish resolution of surface normals or local slopes. The longer 621.74 s,
trim, finite-source and integration-method cases retain six modes, so their
shapes have not been shown to be modally resolved.

A separate read-only audit reconstructs those same nine probes from the saved
states at every shared 5 s output. The complete 12-to-16 mesh difference peaks
at **0.188 mm at 35 s**, still below the 0.720 mm criterion; time refinement
remains negligible. The 80-mode error instead peaks at 2.410 mm at 50 s,
reinforcing its failure. The audit's script, input hashes and results are
recorded in [flexible_checks.json](flexible_checks.json). This extends the
comparison to all saved times, not to unsampled times or membrane locations.

Measured integration costs exclude diagnostics and archive writing:

| Fixed case | Simulated time | Integration wall time | Actual RHS calls |
|---|---:|---:|---:|
| One tile, 6 modes, DOP853 | 30 s | 0.630 s | 167 |
| Same case, RK45 | 30 s | 0.173 s | 62 |
| Same case, Radau | 30 s | 0.452 s | 157 |
| Three tiles, 6 modes each, DOP853 | 30 s | 2.331 s | 167 |
| One tile, complete 12-division basis, DOP853 | 60 s | 2.600 s | 239 |
| One tile, complete 16-division basis, DOP853 | 60 s | 9.801 s | 317 |
| Same complete basis, maximum step 0.5 s | 60 s | 46.424 s | 1,529 |

RK45 and Radau agree with DOP853's terminal shape probes within 1.4e-11 m
on that short six-mode case. Radau reports 81 function evaluations, but
instrumentation counts 157 including its numerical Jacobian. RK45 is fastest
for this particular comparison; no method ranking at complete modal
resolution or over a service orbit has been established. The cheap six-mode
timings cannot serve as the cost of a resolved structural calculation.

An initial Radau run exposed a guard that rejected its Jacobian perturbations
of inactive wheel coordinates. The fix ignores those inactive coordinates in
the equations and diagnostics while still rejecting nonzero initial inactive
wheel momentum. A regression compares implicit and explicit propagation;
the initial failure and source snapshot remain in the ignored run archive.
Another startup check replaces a single absolute tolerance across quantities
with different units by explicit position, velocity, attitude, deformation,
work and impulse tolerances. It reduces actual RHS calls from 491 to 95 over
2 s, with a 3.73 nm terminal position difference. The product records all
tolerances and work/impulse differences; this is not a universal speedup claim.

The unforced six-mode case has a maximum energy drift of 4.55e-13 J over
60 s. Separate analytic and independent kinetic-energy tests matter here:
tiny relative errors against a large orbital energy can conceal structural
errors. The bounded trim allocator completes, but still leaves sampled
residuals up to **3.213 N and 23.332 kN·m**; solver success does not establish
sufficient control authority. Likewise, 100% coverage at the 25 chosen
receiver probes does not test seams or certify complete optical coverage.

The endpoint-lumped support inertia differs from the earlier distributed-frame
reference by **24.83%, 24.19% and 23.97%** in relative Frobenius norm on the
three meshes. This approximation and the assumed hub/boom compliance must
be addressed before treating precise integration as precise attitude physics.
On the finest mesh the principal inertias are approximately 24% too large
relative to that reference; conservation checks verify the declared lumped
model. A distributed support-mass representation or explicit sensitivity
assessment is needed before using its control/torque predictions quantitatively.
The optical polygon area overestimate falls from 0.243% to 0.061% across
those meshes, but refining it does not remove the support-inertia discrepancy.

## Reproduce and continue later

Run from the repository root with one numerical thread:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.geometry_screen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.flexible_benchmark
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest engineering/test_array_flexible.py protection/dynamics/test_flexible_fleet.py protection/dynamics/test_flexible_fleet_energy.py protection/dynamics/test_flexible_optics.py protection/dynamics/test_flexible_environment.py research/studies/solar_shield_array/test_geometry_screen.py research/studies/solar_shield_array/test_flexible_benchmark.py
```

For the later search, rebuild the structural product whenever geometry,
material or attached mass changes; then vary spacing, initial orbit, attitude
request and finite controls together. Preserve failed trajectories and use
measured cost per right-hand-side evaluation to budget each resolution.
Do not extrapolate three-tile timings linearly to a dense full fleet:
intersection work depends on shadow candidates, mesh and source samples.

Before interpreting any future numerical optimum as feasible, refine both
spatial and temporal shadow/clearance sampling, propagate the relevant
ephemeris and attitude cycles, account for hardware optical surfaces and
loads, and replace or bound the support/attachment constitutive assumptions.
Qualification, full fleet/fault/service trials and manufacturing/resource
closure remain beyond this numerical model. The present work supplies a
tested reduced-model starting point and measured compute costs; it does not
claim a completed design search or a verified physical array.
