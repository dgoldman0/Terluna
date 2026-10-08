# Full-ephemeris photogravitational shield study

Allowing the optical formation to move along the Sun–Moon line reduces the
calculated mean holding power from **309 TW to 238 TW**, about **23%**, under
the September propulsion assumptions and a 50-g/m² base optical allocation.
The selected trajectory moves between **67,700 and 105,500 km sunward** of the
Moon. Solar pressure provides a modest additional reduction. The remaining
**371,000 kg/s of propellant** keeps lifetime material use as a binding problem.

These results describe a four-lunar-radius protected region, the selected
climate spectrum and an optimistic optical controller. The material, power
equipment and propulsion parameters are design assumptions. The search covers
bounded trajectory families; it does not select a final array or establish a
passive solution. Habitats, optical device performance and safe plumes still
need their own models.

## The reference and the new calculation

The September model held a screen 78,000 km sunward of the Moon using a
circular phase sweep. This study reads NASA/JPL's **DE440s** ephemeris for
geometric positions and velocities of the Sun, Earth, Moon and planetary
systems. Its binary matches NASA's published MD5 and a recorded SHA-256.
The source and restoration command are in
[protection/dynamics](../../../protection/dynamics/README.md).

The primary interval starts on **4 October 2026 TDB** and lasts 365.25 days.
A 19-Julian-year run includes a lunar nodal cycle and changing distances,
inclinations and solar range. DE440s ends in 2150. A construction-era case near
2526 would use the longer DE440 kernel; this present-epoch study does not
predict geological stability or a billion-year supply history.

Point-mass gravity includes the Sun, Earth, Moon and the other planetary
systems. Earth's and the Moon's existing shared gravitational parameters are
preserved. New solar and planetary parameters come from DE440's kernel
comments. The difference between the point-mass external force at the Moon
and the acceleration derived from its ephemeris is about **1.5×10⁻⁹ m/s²**
at its largest sampled value, small beside holding accelerations near
10⁻³ m/s². Multipole gravity and the ephemeris's full relativistic force model
are outside the spacecraft model.

Let **q** be a tile's position relative to the Moon. Its required nongravitational
acceleration is

\[
\mathbf a_{\rm required}
=\ddot{\mathbf q}-\left(\mathbf g_{\rm tile}
-\ddot{\mathbf r}_{\rm Moon}\right).
\]

The model differentiates ephemeris velocities to obtain the reference
acceleration and uses actual ephemeris phase for the trial trajectories.
It calculates a vector residual after optical control. An independent
circular-limit test reproduces the earlier force expression.

## What sunlight can push

Incoming photons push away from the Sun. Consider a fraction *f* of the
incident optical energy available for redirection, with flux *S*, areal mass
*σ* and incoming propagation direction **e**. Momentum conservation gives

\[
b=\frac{fS}{c\sigma},\qquad
\mathbf a_{\rm optical}=b\left(\mathbf e-\mathbf v\right),
\]

where **v** describes an outgoing photon direction or a weighted mixture.
Allowing every passive outgoing direction produces a ball of radius *b*
centred at *b e*. The code projects the required force onto this ball.
Every acceleration in it has a nonnegative component along the incoming
photon direction. It cannot provide a sunward push.

This envelope relaxes real diffraction angles, mandatory absorption,
photovoltaic extraction, spectral rejection and Earth-exclusion directions.
Its residual thrust is therefore a **lower bound**. A real filter can require
more propulsion. In particular, a zero-momentum state inside this relaxed
description need not perform the required UV filtering.

The central climate window retains 1,173.25 W/m² of the nominal 1,361 W/m²,
so the control envelope can use at most **13.795%** of its incident spectrum.
The outer atmospheric-protection annulus may redirect the entire spectrum
because its rays miss the solid Moon. All tiles whose solar rays can reach
the surface retain the central constraint.

At the fixed reference distance, the required force has a sunward component
for about **78%** of the primary interval. For the selected moving trajectory
the share is about **95%**. Increasing reflective area alone cannot remove that
component. Even an unlimited passive optical actuator leaves the selected
trajectory's positive sunward component, averaging about **0.375 mm/s²**.
The complete required vector averages about 0.761 mm/s² there. Different
trajectories can change those values; this bound belongs to the tested path.

The fixed-distance centre calculation also compares fully absorbing and
fully reflecting optical states. They raise mean holding power to about
**322 and 336 TW**, respectively, while the ideal controller lowers it to
about 296 TW. Optical momentum must be scheduled with the dynamics.

## Aperture geometry and mass closure

The aperture covers the full angular solar disk across the protected region.
Its radius includes the finite-source penumbra, any transverse displacement,
ray-flight motion and a 50-m formation margin. The Moon moves about 8 km while
light travels from the fixed reference screen to it; the selected moving
screen needs a maximum first-order allowance of **10.7 km**. Source motion
also enters that allowance. Retarded solar and terrestrial positions are used
for photon momentum and eclipse geometry, with their Taylor evaluation checked
against exact kernel evaluations.

The final calculation integrates independently held area nodes across the
climate window and outer annulus. It recomputes gravity, optical authority,
eclipse visibility and power hardware at each node. A 96-node quadrature is
checked against 384 nodes for the selected trajectory. The scan and trajectory
fit use a cheaper centre-force approximation; the table below uses area
quadrature.

Power equipment and the seven-day propellant buffer add mass, which decreases
optical acceleration and increases required thrust. The calculation iterates
that coupling. It retains **30 km/s exhaust velocity, 70% electrical-to-jet
efficiency, 45° exhaust cant, 300 W/kg installed power equipment** and a 25%
peak thrust allowance. The cant is a budget for plume routing. Plume isolation
has not been demonstrated. Storage mass is excluded from this comparison.

| Case | Base optical mass | Mean holding power | Installed peak power | Propellant |
|---|---:|---:|---:|---:|
| Fixed 78,000 km, no optical thrust | 50 g/m² | 309.1 TW | 578.7 TW | 480,760 kg/s |
| Fixed 78,000 km, ideal optical control | 50 g/m² | 296.6 TW | 576.0 TW | 461,308 kg/s |
| Moving distance, ideal optical control, finer quadrature | 50 g/m² | **238.3 TW** | **550.5 TW** | **370,763 kg/s** |
| Same moving trajectory with a replacement filter | 10 g/m² | **40.4 TW** | **94.8 TW** | **62,768 kg/s** |

The moving aperture has a radius of about **7,454 km** and a held mass near
**1.079×10¹³ kg**. Its aperture is slightly larger than the fixed screen's,
so total held mass rises by about 1.3% even while mean power falls.
Peak installed power falls by about 5%. Daily energy and fuel benefit more
than installed hardware.

The 10-g/m² case applies a different mass to the same trajectory; it has not
been re-optimized for that mass. The stored filter already contains about
**25.9 g/m² of silica and titania**, before its antireflection layers, supports
or photovoltaics. That case therefore needs a new optical material system.
It still consumes about **2.0×10²¹ kg of propellant per billion years** if
the flow is continuously lost. The 50-g/m² moving case consumes about
**1.17×10²² kg** on the same extrapolation.

The 100-TW civilization design case is a useful delivered-power allocation.
It is separate from gross generation and the array's holding demand. The
much larger outer aperture can collect light whose rays miss the surface;
its electrical yield depends on the optical design. These holding results
specify a load for that power network, rather than a generation estimate.

## Trajectory search, eclipses and propagation

The distance search uses three harmonics of the actual lunar synodic longitude
difference, with bounds of 20,000–180,000 km. Two initial profiles converge to
similar local minima. A second family permits up to 2,000 km of transverse
motion and starts from the improved on-axis solution. It yields no meaningful
power reduction once the larger aperture is counted. This is a local result
for those basis functions and bounds.

The selected profile is closest to the Moon near full Moon and farther away
near the two quadratures. Applied through 19 years without refitting, its
centre-force closure averages **239.4 TW**, compared with **298.9 TW** for
the fixed reference with ideal optical control. The area-quadrature annual
comparison gives the final 238.3-TW figure.

The 19-year eclipse scan samples every ten minutes, then refines detected
contacts. It finds **38 passages**. The longest complete blackout lasts
**2.90 hours**, inside a **4.87-hour** penumbral passage near 7 August 2036.
There are times during that tile blackout when up to **4.5% of the Sun** is
still visible from the Moon. Protection cannot simply be suspended whenever
a shield tile loses local sunlight.

DOP853 propagation follows the chosen trajectory for 30 days using
interpolated DE440 states. With a 1-m initial position perturbation, a
1-mm/s velocity perturbation and ideal electrical feedback, the largest
position error is **1.53 m** and the full-Sun coverage margin stays positive.
The controller has unlimited electrical authority in this check; sensor
delay, flexible structures, collision avoidance and actual actuator limits
remain to test.

With the same feedforward thrust and no corrective feedback, the perturbation
grows to about **366 km** after 30 days and coverage is lost. With solar-only
control on the selected path, coverage is lost within the first sampled
ten-minute interval. Active formation control remains necessary in this
architecture.

A stress test suppresses electrical thrust through the worst penumbral passage.
The formation centre departs from its reference by about **60 km** and
recovers under feedback. The ideal complete aperture retains coverage in that
test. Independent tile seams, differential coasting, regional power links and
constrained recovery thrust are unmodelled. This result makes eclipse
coasting a candidate for a separate formation study; it does not set the
array's storage requirement.

## Checks and the next work

Changing the lunar acceleration difference interval from 30 seconds to 10
seconds changes acceleration by at most **5×10⁻¹² m/s²** in the sampled check;
90 seconds differs by about 4.2×10⁻¹¹ m/s². Ten-minute versus hourly sampling
changes annual centre-force holding power by less than 0.001%. The finer
aperture quadrature changes moving-array power by **0.025%**. Independent
unit-vector, phase, ray, eclipse, momentum and circular-limit tests cover the
physics components. Product tests bind the snapshot to its sources and the
shared constants it actually reads, recorded by name and value.

The next holding work should compare broader trajectory families and fleet
handover geometries under the same shadow constraint. Periodic Earth–Moon
orbits and more remote arrays need a demonstrated aperture history; a
libration orbit by itself supplies no continuous Sun–Moon alignment. External
momentum delivery deserves a supply model that tracks capture, return and
the atmosphere's plume limits. These routes can address lifetime propellant
more directly than further refinement of the present local minimum.

In parallel, the unified filter/PV and optical-routing designs need measured
or computed **T(λ), R(λ), absorption, generated electricity and complete
areal mass**. They can replace the relaxed control ball with achievable force
states. Habitat wheels should be sized on candidate economical trajectories
and compared with attached wheels using their full mass. The current result
supports keeping heavy settlements dynamically separate while those options
are developed.

## Sources and products

The numerical products are [holding.json](results/holding.json) and
[validation.json](results/validation.json). Their source and kernel hashes
are recorded. The original September protection model remains intact.

- Park, Folkner, Williams and Boggs (2021),
  [*The JPL Planetary and Lunar Ephemerides DE440 and DE441*](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf),
  *Astronomical Journal* 161, 105, DOI 10.3847/1538-3881/abd414.
  The abstract and coordinate discussion, and the kernel's technical comments,
  supply the ephemeris attribution and interpretation here.
- NASA/JPL NAIF [planetary kernels and checksums](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/).
- Simo and McInnes (2008),
  [*Solar sail trajectories at the Earth–Moon Lagrange points*](https://strathprints.strath.ac.uk/6990/).
  Its abstract is a lead for the next trajectory study; its communications
  orbits have not been implemented or treated as lunar shielding solutions.
