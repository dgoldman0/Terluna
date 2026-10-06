# Relative-orbit redesign of the shield formation

On 6 October 2026 the author reviewed the expanded coupled-cycle checkpoint
(`abbaf739`) and adopted a change of method for the formation work. This
document records the diagnosis behind the change, the plan, its decision gate
and the choices that remain open. The shield's protection governs every stage:
each candidate passes the same finite-Sun coverage, overlap and UV-transmission
tests as before, and a formation counts only when it protects.

The diagnosis numbers below come from gravity-only replays of the committed
seed states in [results/joint_seed.json](results/joint_seed.json). The Moon,
Earth, Sun and planets act as point masses, taken from the compact DE440 samples
in [local_continuation/ephemeris.npz](local_continuation/ephemeris.npz), with
the coupled model's own environment. [relative_diagnosis.py](relative_diagnosis.py)
reproduces them in about four CPU seconds and writes
[results/relative_diagnosis.json](results/relative_diagnosis.json). Sail force,
mutual shadows, finite-square contact and attitude are absent from these replays;
the coupled model that verifies every design includes them.

## Findings

### The return error is orbital-energy spread from the depth stack

The 361 loaded tiles have osculating semi-major axes about the Moon that span
59.0 km at the epoch, 58.1 km at the end of first service and 64.7 km at
hour 12, with standard deviations of 21–23 km. Depth along the Sun direction
explains the spread almost completely: a linear fit gives 1.57 km of
semi-major axis per kilometre of depth across the 36 km four-level stack, with
R² = 0.9999 at the epoch. The tiles' periods therefore differ by up to
18 minutes. Tiles whose semi-major axes differ by Δa separate along-track by
about 3πΔa per revolution, which predicts up to 290–323 km here.

Replaying the hour-12 states through one 46.48-hour revolution of the centroid
gives a maximum along-track repeat error of 345.1 km (rms 222.9 km). The coupled
replay in [coupled_optimization.md](coupled_optimization.md) reports 344.5 km
at the following departure endpoint, the same point in the cycle. The position
debt that the expanded search worked to remove is this secular drift.

One along-track burn per tile that equalizes the semi-major axes costs 0.42 m/s
rms and 0.64 m/s at most, sized by vis-viva and refined twice from the replayed
drift. It reduces the one-revolution along-track error from 345.1 km to 0.07 km.
The 3-D rms error falls to 2.58 km, mostly a radial residual of up to 5.3 km
from unequal eccentricities, which a second small burn addresses. The expanded
basis reached this direction only indirectly: its early patterns were probed at
0.02 m/s, and its acquisition arcs act in the last eight hours before arrival,
after the drift has grown to about 300 km.

Closing the cycle with energy management means equalizing energy after service,
restoring the service velocity field before the next passage and offsetting the
coast energy for the drift the service itself accumulates. At the present depth
that costs roughly 1.1 m/s rms per tile per cycle. The 23.835 TW operating
target allows 0.6489 m/s per day for the 17.286-million-tile inventory
([energy.md](energy.md)), about 1.26 m/s per 46.5-hour cycle. The burn scales
in proportion to stack depth.

### The hour-13.8 clearance failures are the pattern folding through its orbital plane

Each tile's cross-track offset oscillates once per orbit, so the service
pattern flattens about a quarter orbit after mid-service. Replayed from the end
of first service, the pattern's thickness normal to its orbit falls from
150.6 km at hour 6 to 1.4 km at hour 14.25. Over the same span the number of
tile-centre pairs closer than 10 km rises from zero to 3,707 at hour 14.17
before falling again. Replayed from the hour-12 state, the thickness reaches
1.7 km at hour 14.27. The coupled runs' first 100 m clearance events, at
13.795–13.800 h, come as the pattern closes, about half an hour before its
flattest point.

Relative to the centroid at hour 12, the tiles' relative eccentricity and
inclination vectors have median magnitudes of 45 and 47 km. The median phase
gap between them is 91°, the perpendicular arrangement in which radial and
cross-track separations vanish together. In the linear model, 1,170 of the
64,980 tile pairs have radial-normal minima below 100 m, and 4,947 below 1 km.
Energy matching leaves those counts about the same (1,231 and 5,866), so the
fold calls for a redesign of the e/i geometry itself.

In linear relative motion the in-plane and cross-track motions are independent.
Tiles that share depth and along-track position follow the same in-plane path
and meet at the fold for any optimizer setting. Formation flying keeps such
patterns apart with relative orbital elements. Equal semi-major axes make the
pattern repeat, and relative eccentricity vectors parallel to the relative
inclination vectors put the largest radial separation where the cross-track
separation vanishes. This is the e/i-vector separation of D'Amico and
Montenbruck (2006), flown on PRISMA and TanDEM-X. The minimum separations of
such relative orbits have closed forms, which turns the layout into a small
geometric design problem.

### Layout screen: finite patches meet, ring bundles stay clear

[relative_design.py](relative_design.py) screens layouts of parallel 10 km
squares and writes [results/relative_design.json](results/relative_design.json)
in about 40 CPU seconds. The screen runs in two steps: first linear relative
motion over one to three orbits, then ephemeris gravity over two orbits.

Every finite 2D patch with overlapping rows meets a neighbour somewhere in its
orbit. The geometry behind this is general. Tiles whose projections overlap
need a normal separation that keeps its sign around the whole orbit.

- **Radial-facing patches** (19 rows by 21 columns). For a common radial-facing
  attitude, the normal separation of two equal-energy relative orbits averages
  zero over an orbit. So it changes sign, and at that moment the pair's
  projections still overlap. Rows layered by relative eccentricity, at 2–10 km
  in depth with fold steps of 0.2–0.5 km, come within 0.3–2.4 m. Rows layered
  by semi-major axis come within 0 m, because rows of the same parity converge
  at the fold.
- **Sun-facing patches** keep a mean depth separation through the relative
  eccentricity. Neighbours along a row, though, become coplanar side by side
  near mid-service and meet at u = −2°.

A ring bundle avoids both problems. Each row is one orbit carrying a string of
tiles at equal energy, and the rows step outward in radius, so they slide past
one another and never reassemble. Two conditions keep it clear:

- **Along a ring**, a small shared tilt about the orbit normal makes
  neighbours overlap like shingles. Their normal separation is p sin α, at
  pitch p.
- **Between rings**, the radius step Δr has to beat the tilt across a full
  tile: Δr cos α − 10 km × sin α ≥ 150 m.

Bundles that meet both conditions stay clear for three linear orbits:

| Radius step, tilt | Minimum surface distance |
|---|---:|
| 0.40 km, 1.2° | 178 m |
| 0.45 km, 1.5° | 188 m |
| 0.60 km, 2.0° | 251 m |

During service these bundles overlap in projection by at least 0.60 km between
adjacent rings and 1.15 km along each ring. In projection the 19 by 21 segment
spans 163 km across-track at the service ends and 163–180 km along-track over
the service arc. Continuous rings therefore cover their band as long as they
keep passing.

In ephemeris gravity two further conditions appear:

- **Each ring has to start as time-shifted copies of one trajectory.** Started
  as two-body circular states, tiles on one ring take up different tidal
  responses by phase and diverge by up to 290 m in radius within 16 hours.
- **Each tile has to hold its own attitude**, radial-facing and referenced to
  its velocity. A shared frame loses the shingle offset to the perturbed
  orbit's flight-path angle and to the outer rings' slide away from the
  bundle centre.

With both, a bundle of 19 rings by 21 tiles stays clear for 93 hours of
ephemeris gravity. The minimum is 170.6 m at a 1.2° tilt and 215.4 m at 1.5°,
for radius steps of 0.6–1.0 km. The binding pair is the shingle along a ring.
Adjacent-ring radius steps vary by about ±150 m over that time as the rings
slide along their slightly eccentric common orbit. Sail force and mutual
shadows are absent from this replay. The coupled model adds both.

### Hardware headroom and attitude energy

In the expanded run's executed prefix the installed actuators carry 6.15 times
peak demand, doubling their ratings leaves the solution unchanged, and the
active limit is the 100 TJ next-service energy cap. Attitude used 36.19 TJ
against 0.43 TJ for translation. Two 90° turns per cycle at that cost scale to
roughly 20 TW across 17.29 million tiles, most of the 23.835 TW target.

### Plane tracking for full-disk coverage

The test patch orbits in the ecliptic: its mean orbit normal is antiparallel to
the ecliptic pole, so its plane always contains the Sun. Covering the full disk
takes inclined rings whose lines of nodes turn with the Sun at 0.9856° per day.
A simplified model places Earth and the Sun on circular coplanar orbits about
the Moon, with point gravity and no sail force, and propagates retrograde
circular orbits inclined 10–20° to that plane for 120 days. Earth's tide turns
their nodes as follows.

| Orbit radius | Node rate, inclination 10–20° | Node minus Sun after 120 days |
|---|---:|---:|
| 12,000 km | 0.494–0.470° per day | −61 to −64° |
| 15,000 km | 0.698–0.663° per day | −38 to −42° |
| 17,000 km | 0.848–0.805° per day | −19 to −25° |
| 19,000 km | 1.008–0.957° per day | +1 to −6° |
| 21,000 km | 1.177–1.116° per day | +22 to +15° |

Near 19,000 km the rings follow the Sun without steering. The useful service
fraction there is about 12%, against 15.4% at 15,000 km (arcsin(7,000 km/r)/π),
so about 28% more tiles. Steering inclined rings at 15,000 km would cost roughly
0.5–1.5 m/s per day. The pilot's three-day propagations were too short to see
the drift.

### Areal mass and collection

The September optical cell's 50 g/m² is 22 g/m² of silica substrate (10 µm),
3.9 g/m² of titania absorber and 24.1 g/m² of framing, coatings and metrology
([module catalogue](../../../protection/modules/catalogue.json)). The
[September report](../../../protection/report.md) calls the allocation a target
for a structural bill of materials still to be completed. Power and propulsion
sit in a separate holding pack of about 10.9 g/m² at 300 W/kg, with collectors
on about 1.09% of the aperture.

The annulus carries the same UV job as the climate window, because its rays
cross the upper air and exosphere above the limb out to four lunar radii. It
carries no climate job, since those rays leave the Moon. They continue to Earth
near new moons in eclipse season. A disk about 7,000 km in radius that blocked
visible light would cast a total eclipse some 10,000 km across on Earth for
hours, a few times a year. The annulus is therefore a thin UV absorber that
passes visible light. Titania on a thin support fits; the stack's 1 µm titania
layer is 3.9 g/m² of its 26.

The aperture intercepts about 240 PW. Collection sized to demand is small: the
100 TW civilization case at 300 W/m² needs about 0.2% of the aperture, and
holding about 1%. The 5% dimmer removes about 586 TW across the Moon's disk
([design.md](design.md)). Built as a band-selective semitransparent PV layer,
it would yield roughly 120–180 TW at 20–30% conversion, covering the
civilization case with light the climate design already removes.

### A year of ring planes under sail force

[ring_screen.py](ring_screen.py) follows single tiles on 70 rings for one year
from 4 October 2026 and writes
[results/ring_screen.json](results/ring_screen.json) in about four CPU minutes.
The rings are retrograde at radii of 15,000–21,000 km, tilted 0–20° from either
the ecliptic or the Moon's orbit plane. The model uses DE440s point-mass
gravity, the filter's sail force on radial-facing normals with the bundle
interior's 24% shadow loss, and finite-Sun eclipses. At every crossing of the
Sun meridian it records the ring's radius and its height off the Sun–Moon axis.

A radial-facing tile is lit on one face or the other around almost the whole
orbit, so the sail force keeps pumping the ring's eccentricity:

- **Crossing radius.** Over the year it ranges over 14,650–17,140 km for
  15,000 km rings and 18,546–23,263 km for 19,000 km rings, an eccentricity
  near 0.06–0.08.
- **Crossing height.** At every tilt it swings by 2,500–8,900 km at 15,000 km
  and by 3,500–5,200 km at 19,000–20,000 km. That is a common motion of the
  whole strip pattern against the Sun–Moon axis, together with the Moon
  orbit's 5.1° tilt to the ecliptic.
- **Steady drift.** Steering it away costs 0.06–0.84 m/s per day at 15,000 km,
  0.08–0.13 at 19,000 km, and 0.006–0.024 at 20,000 km for rings referenced to
  the ecliptic. With sail force, the radius where Earth's tide turns the planes
  with the Sun moves out from about 19,000 km to about 20,000 km.
- **Spread between rings.** Rings 5° apart in tilt separate by 224–2,528 km
  over the year. Adjacent rings in a bundle differ by about 0.035°, which scales
  this roughly to 2–18 km over the year. That is tens of metres per orbit,
  within the reach of continuous keeping.

A full-disk ring system therefore needs its common orbit placed on a frozen
solution, with eccentricity and node fixed relative to the Sun under tide and
sail force together. Alternatives are steering the whole pattern or oversizing
it by thousands of kilometres. These runs start from circular orbits, which are
not that frozen solution. Finding the frozen family is the next step of stage 3.

### A zoned aperture for the held screen

[zoned_aperture.py](zoned_aperture.py) repeats the published area quadrature
on the selected moving trajectory, with a base areal mass for each zone and a
chosen redirected fraction for annulus tiles. It writes
[results/zoned_aperture.json](results/zoned_aperture.json) in about ten CPU
seconds, and its 50 g/m² case reproduces the published 238.41 TW. The central
window covers 9.1% of the aperture (a radius of 2,246 km out of 7,454 km).
Even so, at 50 g/m² its tiles take 22.6 TW of holding power on their own.

A visible-passing annulus redirects about as much light as the climate stack,
13.8%, where the published optics let it redirect everything. Each annulus
tile then has less optical help against gravity.

| Central window / annulus | Mean holding power | Propellant | Held mass |
|---|---:|---:|---:|
| 50 / 50 g/m², published optics | 238.4 TW | 370,854 kg/s | 1.08×10¹³ kg |
| 50 / 50 g/m², visible-passing annulus | 249.7 TW | 388,388 kg/s | 1.09×10¹³ kg |
| 50 / 20 g/m² | 112.3 TW | 174,647 kg/s | 4.93×10¹² kg |
| 50 / 10 g/m² | 66.5 TW | 103,464 kg/s | 2.95×10¹² kg |
| 50 / 5 g/m² | 43.7 TW | 67,991 kg/s | 1.96×10¹² kg |
| 26 / 10 g/m² | 55.6 TW | 86,476 kg/s | 2.48×10¹² kg |
| 26 / 5 g/m² | 32.8 TW | 51,002 kg/s | 1.49×10¹² kg |

At the bare stack's 26 g/m² in the window and 5 g/m² in the annulus, the held
screen needs 32.8 TW, within 40% of the 23.835 TW target with no formation
flying, and about 1.6 Gt of propellant a year. Which of these areal masses can
be reached depends on two things: the UV transmission of the annulus film,
which comes from the protection domain, and a structural design for each zone.
The propulsion closure is the published one (30 km/s, 70%, a 45° cant,
300 W/kg and a seven-day buffer).

### The coverage requirement

The atmosphere's [loss response](../../../atmosphere/loss_response/README.md)
allows a UV transmission of a few tenths of a percent, averaged over time and
area, depending on the loss budget. Tiles must therefore overlap. Because the
budget is an average, it can also carry scheduled gaps such as handovers or a
failed tile awaiting cover.

## Plan

**Stage 1, tools and diagnosis (complete, 6 October).**
[protection/dynamics/relative_orbit.py](../../../protection/dynamics/relative_orbit.py)
supplies osculating semi-major axes, quasi-nonsingular relative orbital
elements, the linear relative motion, the closed-form minimum radial-normal
separation and energy-matching burns. Its tests check the linear solution
against Hill's equations and two-body propagation, and the separation against
brute force and the published form. [relative_diagnosis.py](relative_diagnosis.py)
reproduces the findings above as a data product, and
[test_relative_diagnosis.py](test_relative_diagnosis.py) binds that product to
its code, inputs and constants.

**Stage 2, the local redesign.** Redesign the 361-tile patch in relative-orbit
space with the same tiles: 10 km physical squares with 9.89 km clear apertures.
It keeps the six-hour finite-Sun coverage of the moving 20 km receiver window
inside four lunar radii and the 100 m clearance with numerical allowance. The
coast runs at equal energy, with burns that manage it each cycle. E/i-vector
separation carries the pattern through the fold. Attitude schemes are compared:
the present 90° turns, synchronized turns, and a steady rotation that keeps each
tile facing radially. The latter is face-on to the Sun at mid-service and
edge-on at the fold, and meets the Sun at up to about 25° at the ends of
service. Designs are built in the linear model and verified in the coupled
finite-square model with mutual shadows and sail force. Verification follows
the full sequence of service, departure, return, next service and subsequent
departure, with the existing recurrence gates of 50 m, 0.01 m/s and 1e-4 rad.

**Stage 2a, the layout screen (complete, 6 October).** The screen above shows
that the 361-tile patch, as a unit that reassembles each orbit, cannot stay
clear on equal-energy relative orbits. A bundle of ring segments can. On
6 October the author chose the ring bundle as stage 2's local unit, after
confirming that it differs from the rejected prescribed circular family of
[active_global.md](active_global.md). That family enforced Sun-following circles
with continuous thrust, about 22 m/s per tile per day and 1,500 TW in total.
The ring bundle flies natural orbits, and its 93-hour replay uses no thrust.

**Stage 2b, the coupled ring bundle.** The coupled verification runs a segment
of about 19 rings with per-tile radial-facing attitudes, sail force and mutual
shadows. It checks clearance and strip overlap through service and recurrence
over consecutive orbits. The radius steps are ordered as a staircase, so that
every ring has one shadowed edge and the shadow's push on the rings is alike.
[ring_bundle_run.py](ring_bundle_run.py) flies a finite segment, continuous rings
and a lone ring free. [ring_keeping_run.py](ring_keeping_run.py) (stage 2c)
flies continuous rings under low-thrust keeping. Their results are recorded with
the stage 2b and 2c findings.
Runs longer than the six days of compact samples use the pinned DE440s kernel,
which needs jplephem; the author approved installing both.

**Decision gate.** Suppose the cycle closes every coverage, clearance and
recurrence gate within about 1.3 m/s per tile per cycle, with attitude energy
inside the 23.835 TW target. Then the orbiting fleet proceeds to stage 3.
Otherwise the held screen with a zoned aperture from stage 4 becomes the
baseline.

**Stage 3, the global ring screen.** The year-long screen of single rings
(above, 6 October) shows the common annual motion of the strip pattern. Next
come the frozen family of the common orbit under tide and sail force, then
bundles on it. Inclined rings near the Sun-synchronous radius are screened with
the full DE440 ephemeris over a year. The screen
covers nested rings that never cross, the service fraction and inventory, the
recurring corrections, and handovers within the averaged UV budget.

**Stage 4, the zoned aperture and collection, in parallel.** The held
screen's holding power with a zoned aperture is computed above (6 October).
The remaining parts are these. The protection domain supplies the annulus
film's UV transmission. The shield's shadow on Earth is mapped at
eclipse-season new moons. Collection is sized to demand, including the dimmer
as semitransparent PV and small local collectors for holding.

The expanded-cycle continuation package stays preserved as the record of the
black-box approach, and this plan replaces its SQP campaign.

## Open choices

These choices are the author's. The recommendations come from the findings
above.

| Choice | Recommendation |
|---|---|
| Ring radius, about 15,000 km or near the Sun-synchronous radius | Screen 15,000–21,000 km with the full ephemeris in stage 3, leaning toward about 19,000 km |
| Attitude scheme | Radial-facing steady rotation for annulus tiles if their coverage holds; Sun-facing service for climate-window tiles, whose transmitted spectrum changes with incidence |
| Architecture: orbiting fleet or held screen with a zoned aperture | Set by the stage 2 gate |
| Areal-mass targets | Estimates of 10–25 g/m² for the climate window and 7–13 g/m² for the annulus, to be fixed by the annulus film's UV transmission and a structural design |
| Collection | The dimmer as semitransparent PV, if the climate model accepts a band-selective 5% in place of a flat one |
| Scheduled gaps | Spend part of the averaged UV-transmission budget on handovers once stage 3 sizes them |

## Sources

D'Amico, S. and Montenbruck, O. (2006). Proximity operations of formation-flying
spacecraft using an eccentricity/inclination vector separation. *Journal of
Guidance, Control, and Dynamics* 29(3), 554–563,
[doi:10.2514/1.15114](https://doi.org/10.2514/1.15114). Cited for the method;
stage 1 derives and tests the formulas it uses.
