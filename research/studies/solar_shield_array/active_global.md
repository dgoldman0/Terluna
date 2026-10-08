# Rejected rigid traffic and the natural-orbit diagnostic

**The approximately 1,500 TW estimate came from an unsuitable control architecture.**
It continuously enforced Sun-following circular paths for a large circulating
inventory. The author challenged this comparison with the held screen, and
the extended validation was stopped. This estimate is not the cost of
occasional gap repair on naturally evolving trajectories.

This follows the author's request to count how many tiles must move, how often,
and what their transfers, simultaneous handovers and return paths cost. The
previous [local acquisition and retiming](active.md) did not determine that
global duty cycle. The earlier 61.76% large-square probe was not a usable
collision-safe passive fleet.

The new calculation defines an entire repeating population of separate 10 km
squares at 50 g/m². Every tile has an integer identity, a three-dimensional
path, an optical attitude and an actual-date inverse force. The finite-Sun
ray tests use those individual members. Returns are included in the inventory,
collision bounds and propulsion account.

**This is a costly candidate, not evidence that a little extra propulsion
closes the passive fleet's gaps.** It prescribes approximately circular,
Sun-aligned traffic and pays to maintain it. It is not an optimum over freely
evolving orbital families. Coverage remains a sampled result; bounded-slew
attitudes and a coupled, all-member optical/dynamical replay remain open.

## Rejected comparison

| Quantity | Prescribed circular family |
|---|---:|
| Separate 10 km tiles | 38,652,166 |
| Optical base mass at 50 g/m² | 193.26 trillion kg |
| Preliminary mean propulsion power | 1,480–1,520 TW |
| Preliminary propellant rate | 2.30–2.36 million kg/s |
| Conservative power after the assumed control-hardware/reserve mass closure | 1,641 TW |
| Cycle periods | 1.90–7.01 days |
| Front passages | 9.48 million per day |
| Nominal slot handover interval | 16.5–25.6 seconds |
| Minimum proven square clearance, including position allowances | 261.9 m |

The previously optimized held aperture uses approximately **238.35 TW** and
10.79 trillion kg including its control allowance. This candidate's optical
base mass alone is **17.9 times that total mass**. Its lower acceleration per
kilogram does not compensate for the inventory penalty. The held benchmark is
an annual, mass-closed ideal-optics calculation; the rejected budget is a
30-day, nine-date, coarse force quadrature with a separate mass allowance.
Those scope differences preclude a precision comparison, but do not support
continuing this roughly six-times-more-expensive architecture as a low-energy
solution.

The rough arithmetic is revealing: about 18 times the mass, at about one
third the propulsion input per kilogram, still costs about six times as much.
The rejected schedule needs 655–673 m/s of mean accumulated correction per
tile over 30 days. The earlier 5.45–6.79 m/s figures belonged to individual
arrival transfers; they never included this continuing global duty cycle.

The cost bracket permits any mutual-shadow fraction; it is not a full
illumination solution. A smaller 26,369,270-tile variant with 1 km radial
spacing is also rejected: 5,108 actual square intersections occur in the
tested return-pole neighbourhood at the initial date.

The interrupted global ray survey completed 26 dates through day 12.5,
testing 3,833,856 receiver/source/date combinations without a missed ray.
Five 120-second handover windows at 0.5-second cadence and one 32-second
window at 0.125-second cadence also completed without a missed ray. These
local windows extend through day 29.99 and are separate from the interrupted
global survey. They do not establish continuous global coverage.

The planned full-month force census, quadrature convergence and eight-day
rigid-target feedback replay **were not performed**. Their runnable producers
are retained; they must not be cited as completed checks. The preserved
negative result is [active_global_rejected.json](results/active_global_rejected.json).

## Natural-orbit release

The replacement diagnostic preserves the individual tile positions and gives
them lunar circular orbital velocities, then integrates the time-dependent
gravitational ephemeris. It does not enforce circular paths or Sun tracking
after initialization. Its only ideal translation control cancels the
reflected-band photon force. This is a force component bound, not a complete
fleet power budget: gap repair and collision avoidance are still absent.

Every sampled grid trajectory starts at the same epoch. The interpolated
family resolves actual integer members for the optical and intersection
tests. Interpolation is checked against separately integrated integer members.
The seven-day interpolation remains outside the 50 m placement allowance;
no geometry beyond the separately checked three-day prefix is admitted.
The current method interpolates vector coordinates in a rotating orbital frame,
while retaining initial phase as the trajectory identity. It does not shift
an ephemeris trajectory in time.

Initializing natural velocities is a new configuration, not a demonstrated
transfer of the rigid fleet into those states. The release calculation is
intended to locate the failures that a minimum-intervention controller must
address, before assigning that controller an operating cost.

The completed diagnostic uses 37,632 common-epoch grid IVPs, resolved to the
38,652,166 integer tile identities. Independent propagation of 143 members
checks 36 dates through day 3, plus three later dates. Maximum position
difference is 15.60 m through day 3, but 1,127 m by day 7. Accordingly the
following optical results stop at day 3. This is a sampled numerical check,
not a bound on every member's interpolation error.

| Elapsed days | Intercepted area/source rays | Receivers covering all 8 sampled Sun points | Intercepted rim/source rays |
|---|---:|---:|---:|
| 0 | 100.00% | 100.00% | 100.00% |
| 0.25 | 97.27% | 87.33% | 94.26% |
| 0.5 | 88.75% | 57.81% | 81.86% |
| 1 | 77.11% | 28.96% | 65.65% |
| 2 | 79.50% | 41.48% | 69.80% |
| 3 | 77.14% | 31.76% | 63.55% |

Each date tests 4,096 area receivers, 512 rim receivers and eight solar-disk
points. These coarse samples expose failures; they do not establish converged
full-Sun coverage fractions or resolve every handover. This is a different,
much larger inventory from the old 3,072-square probe and its 61.76% result.

A tested subset of 108,288 identities contains 64 intersecting square pairs at
day 1 and 16 at day 3. The subset is selected around the rejected schedule's
nominal return clock; natural motion may carry it away from the actual pole.
It is not a global encounter census. Independent, direct integration of eight
example pairs at each date confirms all 16 intersections, with advected
surface tangents and fresh reflected-beam attitudes. Thus these examples do
not depend on interpolating the trajectory grid. See
[natural_collision_replay.json](results/natural_collision_replay.json).

The mean upper bound for **photon-force cancellation alone** is 68.52 TW over
these three days, or 106,594 kg/s at the assumed propulsion settings. It omits
gap repair, collision avoidance, recovery after intervention and added hardware
mass. It is a component estimate for a fleet that fails its operating
constraints, **not a low-energy global operating budget**. No savings against
the held screen are demonstrated by that number. Full details and miss
witnesses are in [natural_traffic.json](results/natural_traffic.json).

## Geometry, identities and complete return paths

There are 2,304 meridional planes containing the instantaneous Sun–Moon axis.
Each plane has two alternating radial layers. The operational comparison uses
4.5 km between all distinct layers, and a nominal 9.5 km arc pitch between
successive tiles in a plane. Each individual square is 10 km on a side; its
optical clear side is 9.89 km, retaining the previous seam and placement
allowances. No 1,000 km physical panels are used.

Plane radii follow a folded ordering: adjacent plane azimuths have nearby
radii, including the azimuth wrap. The two radial parities share a clock. A
small, smooth phase correction registers their transverse positions during
the sunward passage and vanishes on return. Without that correction, depth
staggering creates a systematic transverse seam. All paths are explicit for
the entire orbit; no member disappears when it leaves the protected disk.

The frame rotates with the retarded Sun direction from DE440s. The clock
subtracts the component of frame rotation already in each orbital plane;
this avoids an unnecessary radial speed mismatch. Orbital senses are
retrograde relative to the ecliptic reference. Paths repeat in the changing
solar frame, with different periods in different planes. The external
ephemeris and required thrust do **not** repeat and are reevaluated at every
actual date. This is not a time-shifted replay of one propagated orbit.

The separate 1 km-layer comparison keeps a smaller inventory but is rejected:
finite-square separating-axis tests find actual intersecting return panels.
These are physical intersections at a tested instant, not merely overlapping
bounding spheres. The recorded count covers a return-pole neighbourhood, not
all possible fleet encounters.

## Coverage and handovers

Each test ray begins at a receiver point in the four-lunar-radius transverse
disk and ends at a finite solar-disk point. A spatial index locates nearby
curved orbital ribbons; final acceptance requires intersection with the
actual plane and clear boundaries of an identified square. Every candidate's
state and attitude are evaluated at the same reception epoch. First-order
tile and lunar motion during light travel are included. The index is tested
against an independent exhaustive finite-square calculation.

The main survey sampled the full receiver area and its exact outer rim through
day 12.5 before it was stopped. Additional short windows resolve actual member-to-member optical
handoffs, including solar-limb sources. These are spatial and temporal tests,
not a sum of projected areas. They also do not constitute a proof for every
receiver, every point on the Sun and every instant. A zero missed-ray count
must retain that sampling boundary.

The receiving disk is the same transverse receiver convention as the previous
fleet studies. The reflected-beam exclusion uses the full protected sphere.
Atmospheric scattering, spectral/angular coating errors and multiple
reflections between panels are not included.

## Separation and reflected beams

For different shells, every square is confined to a bounded radial interval.
The attitude law enforces a 17° maximum normal deviation from the radial
direction. A square corner can therefore depart radially by at most

\[
  h_r \le \frac{s}{\sqrt2}\sin17^\circ
       +\frac{s^2/2}{2(r_{\min}-s\sin17^\circ/\sqrt2)}.
\]

Twice this bound and two independent 50 m position errors are subtracted from
the 4.5 km layer spacing. Members on the same shell are two phase slots apart;
their enclosing spheres remain separated using an analytic lower bound on the
phase-warp derivative. These arguments cover the complete paths between
samples, including transfers into service and the return hemisphere. They
certify geometry conditional on the stated position and attitude envelopes.

The radial-shell proof belongs only to the prescribed schedule. It does not
apply after natural orbital release; the release model explicitly rejects
attempts to reuse that certificate.

Reflections avoid inflated Earth and protected-Moon cones, including the solar
angular radius, square corners, body motion and 1 mrad normal error. A smooth
out-of-plane detour replaces the discontinuous nearest-side choice around the
Moon on return. An in-plane fallback enforces the normal envelope when an
Earth-avoidance correction approaches grazing incidence. No reflected power is
credited as delivered electricity.

The Earth-avoidance law still selects ideal instantaneous normals. Its
history-dependent, bounded-slew realization is not established. In particular,
positive beam margins at the evaluated normals do not certify safe intermediate
attitudes through a branch change. Attitude mechanics, flexible membrane
dynamics and exhaust trajectories remain design gates.

## What the power account means

For a prescribed path, let \(b=\ddot q-g(q,t)\) and let \(s\) be the fully
illuminated ideal reflected-band photon acceleration. Unknown mutual
illumination is a fraction \(0\le\lambda\le1\). Required electric acceleration
is \(b-\lambda s\). Its lower and upper magnitudes are computed exactly on
that line segment:

\[
 a_{\min}=\min_{0\le\lambda\le1}|b-\lambda s|,\qquad
 a_{\max}=\max(|b|,|b-s|).
\]

Thus mutually shadowed panels are not all charged with unattenuated photon
force. The bracket encloses any such shadow pattern under the single-reflection
force model; it does not simulate the actual illumination field. Natural
Earth/Moon eclipses are included. A separate 3×3 quadrature checks finite-square
gravity corrections.

The preliminary force total uses spatial/phase quadrature. Its planned finer
spatial and temporal repeats and full integer-member census were stopped.
Propulsion uses the inherited 30 km/s exhaust velocity, 70% efficiency
and 45° thrust cant. These assumptions convert thrust into electrical input
and propellant. They do not establish a source for that electricity or safe
plume routing.

Distributed power hardware is sized from the sampled peak for each represented
plane and radial parity, with 25% margin. It is not sized solely from the
simultaneous aggregate peak. A separate conservative mass closure adds the
assumed 300 W/kg power-hardware allowance and seven days of propellant,
including their added thrust demand. This is not a complete vehicle mass
budget: collector, habitat, tank, thermal and detailed actuator requirements
have not been closed.

The unperformed feedback runner would propagate selected member identities
with bounded PD correction and an imposed photon disturbance. It is retained
as code only; it contributes no result to this report.

## What this establishes and leaves open

The study now supplies a concrete global inventory, cadence, complete return
paths, a finite-body separation argument and recurring propulsion costs for
one family. It also supplies an explicit rejection of the more compact
return geometry. The high cost belongs largely to maintaining prescribed
circles and plane alignment throughout the cycle. It is not a universal lower
bound on active orbital choreography.

A lower-energy search must let the members follow the actual perturbed
ephemeris more freely while retaining the spatial coverage, return ordering
and collision constraints. Cheap isolated arrival corrections cannot by
themselves establish that result. Continuous finite-Sun coverage certification,
safe bounded-slew attitudes, plume routing and a full hardware/resource account
remain unresolved. Collector/habitat demand and delivered electrical power
remain separate from every demonstrated result here.

## Reproduce and inspect

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_global_run --census
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_global_validate --mode geometry
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_global_validate --mode budget
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.active_global_validate --mode feedback
python -m research.studies.solar_shield_array.active_global_validate --mode assemble
python -m pytest protection/dynamics research/studies/solar_shield_array -q
```

The completed negative result and shared constants are hashed in
[active_global_rejected.json](results/active_global_rejected.json).
The longer reproduction commands above describe the interrupted study, not
completed validation. To reproduce the release diagnostic, run
`python -m research.studies.solar_shield_array.natural_traffic_run`, then
`python -m research.studies.solar_shield_array.natural_collision_replay`.
`python -m research.studies.solar_shield_array.reject_global` regenerates the
rejected comparison from the retained, identity-checked focused geometry and
interrupted global log; it does not rerun the abandoned full-month study.
Raw replay states and intermediate validation products remain under ignored
`research/runs/solar_shield_array/active/`. Repository checks are recorded in
[checks.json](checks.json); numerical implementation tests are distinct from
physical validation.
