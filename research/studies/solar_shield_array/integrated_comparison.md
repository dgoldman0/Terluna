# Integrated comparison of the shield candidates

On 6 October 2026 the author set the acceptance for the whole shield system
([research/decisions.md](../../decisions.md)). It has to meet the shield
requirements, give reasonable net-positive electricity for habitat, balance
its overall resource requirements, and be stable and safe. Every mass flow is
judged by its direction as well as its size: whatever a design expels must
stay away from the Moon's cold, dense exobase. This comparison replaces the
decision gate of the [relative-orbit plan](relative_orbits.md), which tested the
ring fleet's recurrence cost against a 23.835 TW target. The 23.835 TW figure
and the 100 TW civilization case remain as comparison references, with higher
demand as sensitivities.

A candidate has to pass every gate: the shield requirements, stability, safety
and net-positive electricity for habitat after its own loads. Candidates that
pass are weighed on resources: hardware, fresh feedstock, recycled throughput
and delivery. [integrated_ledger.py](integrated_ledger.py) gathers what the
products establish, adds the arithmetic that follows from stated assumptions,
marks what no model yet supplies and writes
[results/integrated_ledger.json](results/integrated_ledger.json).
[exhaust_isolation.py](exhaust_isolation.py) follows the held screen's exhaust
and writes [results/exhaust_isolation.json](results/exhaust_isolation.json).

On 7 October the author made the ring fleet the lead candidate and adopted S7
([decisions.md](../../decisions.md)). The held screen stays a documented
fallback, for which two literature checks continue. The ring fleet's form is
set:
- radial-facing tiles with centred filters, 25% trim and a momentum store;
- rings near 20,000 km, tilted from the Moon's orbital plane;
- one nested bundle, with annulus tiles near 13 g/m² that carry the 4 µm
  silica film.

## The candidates

**The held screen with a zoned aperture** holds independent tiles on the
selected moving trajectory, 67,700–105,500 km sunward of the Moon over each
month and 90,700 km on average
([the plan's zoned-aperture section](relative_orbits.md#a-zoned-aperture-for-the-held-screen)).
The window, 2,245 km in radius, carries the climate stack at 26 g/m²; the
annulus out to 7,454 km carries a 5 g/m² UV film that passes visible light.
Its area is about 1.75 million panel-equivalents of 10 km (1.78 million of
9.89 km clear aperture), with no packing layout yet.

**The ring fleet** flies nested ring bundles on natural orbits about the Moon,
each ring one orbit of shingled tiles, kept with photon forces. No full-fleet
run exists. Its inventory below follows from the bundle's geometry: one ring
for every 9.29 km height step across the aperture, and one tile for every
8.5 km of ring, with every tile of a ring that crosses the window carrying the
window's areal mass. The resources section adds the oversizing the planes'
motion needs and the annulus tiles' mass that nesting needs.

## The gates

Each state is an assessment of the cited evidence for the author: met,
conditional (met if a stated condition holds), open (no model yet decides it)
or not met. The ledger product carries the evidence behind each.

| Gate | Held screen, zoned | Ring fleet |
|---|---|---|
| UV transmission (O1, O8) | Open: a 0.1 µm titania film on silica stops the extreme and far ultraviolet from 3.5 g/m², but films of this kind pass solar X-rays: 3.7×10⁻⁵ of the sunlight below 175 nm at 5.8 g/m² and 1.0×10⁻⁵ at 10.2 g/m², against the window stack's 6.7×10⁻⁷. Traced along slant paths, the 4 µm film's X-rays add a tenth to a fifth of the sky glow's heat at solar maximum, and a gap heats the air 2.5–5 times the disk count, which O1 now replaces with the traced heat (below); continuous area with no seams | Open: with every tile of a 589-tile patch propagated, no receiver ray through its interior goes uncovered in four orbits (below); the annulus film as for the held screen; a full fleet's coverage, handovers and edges are unchecked |
| Protected radius (O2) | Met: the 7,454 km aperture covers four lunar radii with the finite Sun | Open: on a Sun-tracking eccentricity the sunward crossing radius holds within 800–4,100 km over a year. The planes shift the whole pattern about 1,600 km each way at 20,000 km, which 23% more rings cover; there the interior overlaps hold, while the edge strips recede by up to 1,600 km a year, beyond photon steering for the outermost rings. At 15,000 km the strips cross within weeks (below) |
| Window spectrum (E2, O4) | Conditional on the 26 g/m² window and the dimmer's form | The same, on the rings that cross the window |
| Earth's shadow and rejected light (O5, O7) | Open: the optical bound lets redirected light leave toward the Moon; the shadow on Earth is unmapped | Open: night-side tiles reflect back past the Moon, with light the screen has already filtered (to check); the shadow on Earth is unmapped |
| Stability | Conditional on continuous thrust (75.7 TW installed); an unpowered tile drifts about 2,700 km in its first day | Conditional: the kept bundle's edges first come within 150 m at day 8.1; its interior stays clear for 23 days with one tile per ring, and for four orbits with every tile of a patch propagated. Edge-on turns collide with the next ring; a filter centred on each tile cuts the steady shadow torque to a third, which 25% trim holds with a momentum store of about 2×10⁹ N·m·s per tile (below). Met if the centred filter, trim and store are built and keeping holds each tile within its overlap |
| Outflows (S7, proposed) | Not met as designed: without the magnets by the direct plumes and by the slow gas alone; with them by the slow gas and, for 1 kg/s at a 45° cant, the fast atoms (below) | Conditional: photon keeping releases nothing, and with 2° tilts and 25% trim it reaches 1.6 times the keeping demand over each orbit (below); met if a phase-scheduled law delivers it. Trim and a momentum store hold the attitude and release nothing. Electric keeping would release 0.4 t/s or more at or inside the planned magnetosphere |
| Holding cheaper than the atmosphere it saves (S6) | Not met: 51,000 kg/s against 300–80,000 kg/s of unshielded loss | Conditional on photon keeping |
| Net-positive electricity for habitat | Open: its own load is 32.8 TW mean | Open: its own loads are reflectivity trim and the momentum store, not yet sized for power |

## Where the held screen's exhaust goes

Holding a screen upstream of the Moon needs a sunward push that sunlight
cannot give, so all of its reaction mass leaves into the Moon's hemisphere. At
the zoned design point that is 51,000 kg/s carrying 23 TW of jet power. The
proposed requirement S7
([requirements.md](../protection_architecture/requirements.md)) keeps
outflows out of the escape region. The September report's energy diagnostic
([protection/report.md](../../../protection/report.md), section 8) lets a
tenth of the deposited energy do escape work against 0.94 MJ/kg, so each kg/s
of loss budget allows 9.4 MW deposited there: 4.1×10⁻⁷ of the jet power for a
1 kg/s budget and 4.1×10⁻⁶ for 10 kg/s.

**The direct path.** Each tile's propulsive need is split between two
thrusters canted from the resultant, in the plane that keeps both plumes off
the Moon's side. Sampled daily over the year at 96 nodes, and weighted by each
node's propellant flow, this is the share of fast ions that leaves on straight
lines into the protected sphere, whose half-angle from the screen is
3.8–5.9°:

| Cant | Propellant and power | Gridded ion (NEXT) | Hall (BPT-4000) |
|---|---:|---:|---:|
| 45° | 1 | 6.2×10⁻⁶ | 1.6×10⁻⁴ |
| 60° | 1.41 | 1.6×10⁻⁶ | 6.7×10⁻⁵ |
| 75° | 2.73 | 4.6×10⁻⁷ | 2.6×10⁻⁵ |

At 45° the worst node and date reach 7×10⁻⁵ for the gridded ion plume and
1.3×10⁻³ for the Hall plume. Canting each pair in the plane through the Moon
instead points one plume of a pair into the protected sphere at times and
sends 1–2% of the fast ions there. Hall plumes miss the allowance for 10 kg/s
by 6–40 times at every cant. A gridded ion plume meets it from about 60° and
comes to the 1 kg/s allowance at 75°, at 1.4 and 2.7 times the 45° design's
propellant and power. The plume profiles are read from published figures
(Young and others, 2019, for NEXT; Goebel and Katz, 2008, for the BPT-4000) to
about ±30%. These shares count fast particles on straight lines: the ions, and
the fast atoms that charge exchange makes from them on the same paths.

**The slow gas.** Flight-class Hall thrusters reach about 95% mass
utilization, so at least about 5% of the propellant leaves unionized, at
thermal speed; at 5–15% that is 2,550–7,650 kg/s here. Test particles left
every node every two days for a month, with each node's own velocity, a cosine
law about each thruster axis and the effusive speeds of oxygen molecules or
xenon atoms at 500 K, and flew for 60 days under DE440s point-mass gravity:

| Share of the unionized gas | Oxygen molecules | Xenon |
|---|---:|---:|
| Enters the protected sphere within 10 days | 1.2% | 4.6% |
| Enters it within 60 days | 2.2% | 8.1% |
| Passes within 10 Earth radii of Earth | 12.4% | 3.6% |
| Leaves beyond 3 million km | 30.4% | 11.5% |
| Still in flight at 60 days | 55.1% | 76.7% |
| Delivered, with a 30-day neutral lifetime | 1.5% | 5.5% |
| Delivered, kg/s | 38–114 | 140–420 |
| Kinetic power on arrival, MW | 30–91 | 97–291 |

The gas arrives at about 1.2 km/s at four lunar radii and gains speed as it
falls further in. Halving the integration step changed the fate of none of the
first launch's 2,304 particles. Neutral lifetimes of 10 to 100 days give
delivered shares of 1.0–1.9% for oxygen and 3.5–7.1% for xenon. The slow gas
alone deposits 3–30 times the 1 kg/s allowance, and at
the upper end three times the 10 kg/s allowance. The September report
proposed capturing the neutral gas locally; it would have to stop most of it
(below). A further 4–12% of the gas passes near Earth, which S7 also asks to
trace.

**The charged exhaust.** The other 85–95% of the flow is ions. 51 t/s is about
10⁵ times the solar wind's own mass flux through the aperture, so the screen
acts as a comet upstream of the Moon. Without a magnetosphere the wind loads,
slows and carries the ions downstream over the Moon, energising picked-up
oxygen to tens of keV, and the
[loss response](../../../atmosphere/loss_response/README.md) credits each
returning pickup ion with 1–10 sputtered molecules. No model here follows that
plasma, and it adds to both paths above.

**With the September magnets.** The charged-particle protection of the
requirements (C1), a lunar dipole of 1.5×10²¹ A·m², stands the solar wind off
at 10 lunar radii, 17,400 km. At its magnetopause, where about 71 nT balances
the wind ([protection/report.md](../../../protection/report.md), section 10),
exhaust ions gyrate on 140 km (oxygen molecules) to 575 km (xenon), and
picked-up oxygen on 1,900 km, so the magnetosphere holds them off. Xenon picked
up at the wind's speed gyrates on 7,700 km, a diameter close to the stand-off,
and is held poorly. What passes is neutral: the fast atoms that charge
exchange makes in each plume, which keep the beam's paths, and the slow gas.
Goebel and Katz's space-condition model of a 3 kW Hall thruster makes
charge-exchange ions at about a tenth of its beam; taking that tenth for both
thrusters, the capture of unionized gas each case needs is:

| Loss budget (allowance) | Thrusters | Cant | Fast atoms | Unionized gas to capture |
|---|---|---|---:|---:|
| 1 kg/s (9.4 MW) | Gridded ion | 45° | 14 MW | none suffices |
| | | 60° | 3.7 MW | 81–98% |
| | | 75° | 1.1 MW | 73–97% |
| | Hall | 45–75° | 58–359 MW | none suffices |
| 10 kg/s (94 MW) | Gridded ion | 45° | 14 MW | up to 73% |
| | | 60–75° | 1.1–3.7 MW | up to 68–69% |
| | Hall | 45–60° | 153–359 MW | none suffices |
| | | 75° | 58 MW | up to 88% |

The capture ranges span oxygen and xenon at 5–15% unionized, with a 30-day
neutral lifetime. Entry through the magnetosphere's cusps and by
reconnection, and the poorly held xenon, need a plasma model.

**What follows.** As designed, at a 45° cant and with no capture, the held
screen's exhaust breaks the proposed S7: without the magnets on both paths
computed here, before the plasma path is counted, and with them through the
fast atoms and the slow gas. With the magnets, meeting S7 takes gridded ion
thrusters, canted about 60° for a 1 kg/s budget or 45° for 10 kg/s, capture of
most of the unionized gas (81–98% for 1 kg/s, up to about 70% for 10 kg/s), and
a plasma estimate of what leaks past the magnetosphere. The held screen also
fails S6: its 51,000 kg/s of propellant exceeds the atmosphere's unshielded
loss of 300–80,000 kg/s everywhere but the top of that range.

## What holds the ring fleet: photon control

[photon_control.py](photon_control.py) replays the twelve kept orbits from
their checkpoints and writes [results/photon_control.json](results/photon_control.json)
in about two CPU minutes. The kept run's tiles carry 62.7 g/m² (the seed's
50 g/m² and its holding pack); lighter tiles gain sail reach and lose gravity
torque in proportion.

**Translation.** The filter's redirected band reflects specularly along each
tile's normal, so tilting the normal and trimming the band's reflectivity move
the sail force within a bounded set. Against the keeping law's demand on each
ring, with each tile's lit fraction behind the bundle's shadows:

| Tilt and trim | Interior: reach over demand per orbit | Interior: share met at the instant | Edges: reach over demand |
|---|---:|---:|---:|
| 1°, 10% | 0.67 | 44% | 0.38 |
| 2°, 25% | 1.58 | 63% | 0.89 |
| 5°, 25% | 2.13 | 70% | 1.13 |

The demand is the inflated one of the fixed keeping frame. A law that moves
its work to the lit phases of each orbit can deliver it with tilts of about 2°
for interior rings and 5° for the edges; 5 g/m² tiles have 12.5 times the
reach. Photon keeping holds on translation.

**Attitude.** It does not hold on attitude as built. Within a shingled ring each
tile lies under its neighbour along one edge, so the lit part of its redirected
band sits off centre: the centre of pressure lies 1.25 km from the tile's centre
(median; 3.1 km at the 95th percentile). That puts a radiation torque of about
52,000 N·m on every tile (median), whose pitch part averages 27,600 N·m over
each orbit, while the 1.5° shingle tilt holds each tile off
gravity gradient's equilibrium with a further 6,000 N·m (490 N·m at 5 g/m²).
Trimming reflectivity across a tile's halves would need more than the whole
band (128% at the median), and a momentum store would have to hold about
10¹⁰ N·m·s per tile. The ways of holding it are compared
[below](#what-holds-the-ring-fleet-the-attitude).

## What holds the ring fleet: the frozen orbit

A radial-facing tile is pushed toward the Moon on the day side and away from it
on the night side. Averaged over an orbit, that drives the eccentricity vector
at a steady rate across the Sun line, while Earth's tide turns the apse line
against the Sun; in the frame that follows the Sun the push balances at a forced
eccentricity along the Sun line, the balance that gives high area-to-mass
satellites their Sun-tracking (heliotropic) orbits (Colombo, Lücking and
McInnes, 2012). A ring started on a circle traces a circle about that
equilibrium, which is the swing the year-long screen found.
[frozen_rings.py](frozen_rings.py) finds the equilibrium by iteration over a
year of DE440s gravity with the sail force and shadows of the year-long screen,
starting each pass on the centre the previous pass found, and writes
[results/frozen_rings.json](results/frozen_rings.json) in about 15 CPU minutes.
Across ring radii of 15,000–20,000 km, tilts of 0–20° and both reference
planes:

| Tile areal mass | Forced eccentricity (apolune toward the Sun) | Sunward crossing radius spread over a year | Lowest perilune |
|---|---:|---:|---:|
| 62.7 g/m² | 0.04–0.07 | 600–2,400 km | 13,300 km |
| 26 g/m² | 0.08–0.15 | 800–4,100 km | 11,300 km |
| 5 g/m² | 0.38–0.48 | 1,800–10,000 km | 2,800 km |

Started on its circle, a 62.7 g/m² ring swings its crossing radius by
1,800–5,400 km; started on the forced eccentricity, by 600–2,400 km. What
remains of the eccentricity swing is the kind of correction reflectivity trim
supplies. Two things remain beyond it:
- **The planes.** Eccentricity cannot hold a ring's plane. At 19,000–20,000 km,
  near the radius where Earth's tide turns the planes with the Sun, the strip
  pattern still moves 2,400–3,400 km over the year; at 15,000 km tilted rings
  lag the Sun by 23–28 km a day. The pattern needs oversizing by roughly that
  much, or steering of its planes.
- **Light tiles.** At 5 g/m² the forced eccentricity reaches about 0.4, and
  the night-side perilune falls inside four lunar radii, as low as 2,800 km.
  Annulus tiles near the window's 26 g/m² keep the eccentricity at 0.08–0.15.
  By the averaged forcing, tiles turned edge-on outside a ±25° service arc
  would cut the push, and with it the forced eccentricity, about 3.4 times;
  the next section finds that the bundle's clearance rules those turns out.

## What holds the ring fleet: the attitude

[attitude_schemes.py](attitude_schemes.py) holds each tile's attitude in
several ways along the twelve kept orbits and writes
[results/attitude_schemes.json](results/attitude_schemes.json) in about four
CPU minutes. Its tiles as built reproduce photon control's torques.

**Clearance rules out edge-on turns.** Rings step outward by 1 km and slide
past one another, so each tile meets the tiles of the next ring at every
along-track offset. Tiles that share an attitude are parallel squares, and two
of them collide wherever they overlap in their own planes with less than
150 m between those planes. Checked against its own ring's neighbours and the
three nearest rings on each side, a tile may pitch 3° below and 6° above
radial-facing. Between 1° and 2° above, it lies flat against its own ring's
neighbours, because that pitch cancels the shingle. A tile turned edge-on in
pitch outside service meets the next ring's tiles over 82% of each orbit.
Rolling about the line through its ring's tile centres turns the ring's tiles
together and keeps the 222 m between shingles. The next ring then leaves one
side free on the day half of the orbit and the other on the night half, so a
rolled tile passes through radial-facing at the nodes.

**Rolling halves the push at most.** The edge-on turns would have cut the
drive on the eccentricity to 0.36 of the radial-facing tiles'. Rolls of 30°,
60° and 80° outside service, inside the clearance, cut it to 0.86, 0.57 and
0.48. They are costly. Turning a 10 km tile against gravity gradient and the
orbital rate takes 0.2–2 million N·m at 62.7 g/m² (median to 95th percentile),
and the momentum store grows to 2–6×10¹⁰ N·m·s.

**A centred filter holds the attitude.** The shadow torque comes from the
overlap: each tile's filter is shaded along the 1.4 km where it lies under its
neighbour. The filter can instead cover only the middle of each tile along the
ring, the 8.5 km pitch plus a narrow overlap, with the overhanging edges left
as bare film that passes the light. Each filter is then shaded only across
that overlap. The steady pitch torque falls from 27,500 N·m to:
- 6,300 N·m with a 240 m overlap;
- 9,200 N·m with 400 m;
- 12,200 N·m with 560 m.

With 25% reflectivity trim, the momentum store then stays near
1.9×10⁹ N·m·s per tile at every areal mass from 5 to 63 g/m². The tiles as
built need 1.5–3.0×10¹⁰ N·m·s, and more with every orbit. The sail force and
the push do not change. The overlap has to cover the neighbours'
along-track swing, ±160 m over five days in a free string (ring_bundle's lone
ring), and the parallax at the ends of service. A 400 m overlap keeps 240 m at
the far end of the swing.

**A moving mass also holds it.** On the tiles as built, the mass carries the
centre of mass toward the point where the torques cancel:
- 0.5 km of reach, with the same trim, leaves a store of 7.8×10⁸ N·m·s. With
  4.5 km of travel that mass is 12.5% of the tile's.
- 1 km of reach (29% ballast) leaves 5.2×10⁸ N·m·s.

Over the fleet's 17–23 Gt of tiles, that is 2–7 Gt of ballast.

**Steering by roll.** A 35° roll on the half orbit after service only, inside
the clearance, turns the ring's plane. At 62.7 g/m² and 15,000 km the sunward
crossing climbs or sinks 0.77 km a day, and the strip turns 0.002° a day. The
planes' section below sets this against their motion.

**One sail loading per nested bundle.** Rings 1 km apart cross wherever their
radii differ by more than that step. Nested rings must therefore share their
eccentricity vector to within the step over the radius, 7×10⁻⁵ at 15,000 km.
The forced eccentricity follows the sail force per unit mass, so every ring of
a nested bundle needs the same loading. With 5 g/m² annulus rings among
26 g/m² window rings, the two would be forced to eccentricities 0.3 apart and
cross by thousands of kilometres. Photon keeping reaches about a quarter of
the sail force and cannot hold a difference of that size. The annulus section
gives the tile mass that matches the window's loading.

So radial-facing tiles hold their attitude without added mass, given:
- a centred filter;
- 25% trim;
- a momentum store near 2×10⁹ N·m·s per tile, whose hardware is open.

Their push stays as the frozen orbit found it.

## The annulus film

[annulus_film.py](../../../protection/spectra/annulus_film.py) in the
protection domain evaluates films of titania on silica for the annulus, bare
and with the climate stack's anti-reflection and matching layers. It uses the
domain's thin-film code from 0.1 to 2,500 nm and writes
[annulus_film.json](../../../protection/spectra/annulus_film.json). Weighted by
the quiet-Sun spectrum (WHI 2008 below 202 nm, TSIS-1 above), three coated
films and the window stack compare as follows. The heat column is the upper
air's heat at solar maximum, from the escape model's film heat.

| Coated film | Mass | Worst, 10–175 nm | 2.5–175 nm | Below 175 nm | Heat | Visible passed | 2R + A |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.1 µm titania on 1 µm silica | 3.6 g/m² | 4×10⁻⁷ | 3×10⁻⁵ | 1.1×10⁻⁴ | 2.0×10⁻⁵ W/m² | 98.6% | 0.136 |
| 0.1 µm titania on 2 µm silica | 5.8 g/m² | 1.4×10⁻⁷ | 4×10⁻⁶ | 3.7×10⁻⁵ | 6.8×10⁻⁶ W/m² | 98.6% | 0.137 |
| 0.1 µm titania on 4 µm silica | 10.2 g/m² | 1.3×10⁻⁷ | 1.6×10⁻⁷ | 1.0×10⁻⁵ | 1.9×10⁻⁶ W/m² | 98.7% | 0.139 |
| Window stack, 1 µm titania on 10 µm silica | 26.9 g/m² | 2×10⁻⁴³ | 2×10⁻¹² | 6.7×10⁻⁷ | 1.2×10⁻⁷ W/m² | 97.8% | 0.139 |

- **A tenth of a micrometre of titania stops the extreme and far ultraviolet.**
  From 3.5 g/m² of coated film, no wavelength from 10 to 175 nm passes more
  than 10⁻⁶. The coating passes 97–99% of the visible light. Bare titania films
  reflect a quarter of the sunlight and pass 72–82% of the visible.
- **Soft X-rays set the silica.** Between 2.5 and 10 nm, the share of the
  sunlight below 175 nm that passes falls with the silica's thickness:
  3×10⁻⁵ through 1 µm, 4×10⁻⁶ through 2 µm and 1.6×10⁻⁷ through 4 µm.
- **Hard X-rays pass every light film.** Below 2.5 nm the quiet Sun sends
  2.4×10⁻⁵ W/m² of the band's 0.018 W/m². With them, a 5.8–10.2 g/m² film
  passes 1.0–3.7×10⁻⁵ of the sunlight below 175 nm. That is at or above the
  tight level of O1's swarm budget (3×10⁻⁵), against the window stack's
  6.7×10⁻⁷.
- **Traced, the films' X-rays add a fraction of the glow's heat.** The escape
  model's film heat counts every transmitted X-ray as heating the air above
  the base. By that count the 5.8 g/m² film adds 6.8×10⁻⁶ W/m² at solar
  maximum, against 8.4×10⁻⁷ W/m² for gaps at O1's standard level and
  1.2×10⁻⁷ W/m² for the window stack. Traced along their slant paths with the
  X-rays' measured solar cycle ([the X-ray check](#the-x-ray-check)), the
  4 µm film adds 2.5–6.2×10⁻⁷ W/m² at solar maximum and the 2 µm film
  0.9–2.3×10⁻⁶ W/m², a tenth to a fifth and a third to three fifths of the
  sky glow's heat.
- **Pressure.** Coated films have a pressure coefficient (2R + A) of 0.12–0.16,
  and bare ones about 0.5. The dynamics models carry the window's redirected
  13.8% as reflected, 0.276. An annulus tile therefore matches the window
  tiles' sail loading at 26 g/m² × 0.137/0.276 ≈ 13 g/m², which leaves room for
  the 10.2 g/m² film and about 3 g/m² of structure. The computed stack absorbs
  4.5% and reflects 4.7%; if the dimmer's 4.5% is reflected, the window's
  coefficient is 0.23 and the match is 15.5 g/m².

For the ring fleet as one nested bundle, the annulus tiles come to about
13–15.5 g/m². The fleet's mean areal mass then rises from 10.4 to 16–18 g/m².
For the held screen, a 5 g/m² annulus is the 2 µm film's mass, with the X-ray
share above. Its zoned design with a 10 g/m² annulus needs 55.6 TW and
86,000 kg/s ([zoned_aperture.json](results/zoned_aperture.json)), against
32.8 TW and 51,000 kg/s at 5 g/m².

## The ring planes

[plane_motion.py](plane_motion.py) flies nine rings across the aperture's
height at each radius from 15,000 to 22,000 km for a year, from both reference
planes. Each ring starts on its forced eccentricity, with 26 g/m² tiles. Each
sunward passage is measured at the Sun meridian and at 80% of the aperture's
half-width on either side. The rings' planes carry part of the strip pattern's
motion: the change in each crossing's elevation seen from the Moon, times its
distance. That part splits into a shift and a rotation common to the whole
pattern, and the differences between strips. The run writes
[results/plane_motion.json](results/plane_motion.json) in about 24 CPU minutes,
with a checkpoint after each pass.

**The common motion is a shift.** The tide holds the ring planes near the
Moon's orbit about Earth, while the Sun moves 5.1° above and below that plane
over a year. The whole pattern therefore moves up and down by about the radius
times sin 5.1°, and turns by up to 7–9° in the aperture plane. A pattern of
full-width strips covers the disk in any rotation, so oversizing has to cover
only the shift. Centred on its mean, the shift and the extra rings it needs
are:

| Radius | Shift each way | Extra rings |
|---|---:|---:|
| 15,000 km | 1,400 km | 20% |
| 19,000 km | 1,560 km | 22% |
| 20,000 km | 1,590 km | 23% |
| 21,000 km | 1,750 km | 25% |
| 22,000 km | 1,890 km | 27% |

These are for rings set from the Moon's orbit plane; set from the ecliptic,
the rings shift a little more.

**The differences come from tilt.** A ring tilted further from the Moon's orbit
plane precesses more slowly, so the outer strips lag the Sun. Three
measures follow at each radius:
- the change in overlap between neighbouring strips per 9.29 km height step;
- how far the edge strips move against the common motion;
- the photon steering each ring would need, against the 35° roll's reach at
  26 g/m².

At 20,000 km the photon steering needed is 0.2–0.4 of the reach for the
interior rings and 1.8–3.3 times it for the two outermost on each side. The
radii compare as follows:
- **15,000 km.** The outer strips drift 40 km a day against the middle ones,
  overlaps close by up to 24 km per step and the strips cross within weeks.
  Steering would need up to 21 times the photon reach.
- **19,000 km.** The edges close in at 15–26 m a day per step, and the edge
  strips move 2,100–3,500 km in the year. Steering would need up to 6.4 times
  the reach.
- **20,000 km.** The interior overlaps hold, opening by at most 615 m per step
  in the year, mostly the half-monthly breathing, against a 604 m overlap. The
  edge strips move up to 1,600 km.
- **21,000–22,000 km.** The overlaps open by 2.2–7.6 km per step in the year.

The edge strips' drift has a steady part: the top ring's plane turns 0.0095°
a day against the Sun, about 3.5° a year. Steering by roll needs the rolled schemes'
larger momentum store, about 2×10¹⁰ N·m·s per tile at full reach.

So the ring screen belongs near 20,000 km, with rings set from the Moon's
orbit plane. It then needs:
- about 23% more rings for the common shift;
- an overlap a little above 615 m;
- one of these for the outer rings: photon steering with a larger store, a
  radius profile that keeps each tilt turning with the Sun (the outer rings
  about 5% further out by the tide's cos *i* dependence, which needs the
  radius to grow toward both edges of the stack), or yearly oversizing of
  about another 23%.

The forced eccentricity also grows toward the stack's edges, from 0.11 in the
middle to 0.14 at the top at 20,000 km, about 7×10⁻⁵ per ring there. With
1 km radius steps, neighbouring rings near the edges would cross: at perilune
near the top and at apolune near the bottom. The radius step there has to
reach about 1.8 km.

## The bundle with every tile propagated

[bundle_validation.py](bundle_validation.py) flies a patch of the kept bundle,
31 tiles on each of its 19 rings at the 1 km radius step, for four orbits
(7.6 days). Every tile is propagated under the full forces:
- DE440s point masses, with finite-square quadrature;
- the filter's sail force on each tile's own normal;
- finite-Sun shadows among all 589 tiles;
- Moon and Earth eclipses.

Each ring's middle tile is kept as in the kept run. The other tiles are pulled
to their slots, the positions the kept run's copies assumed. The run writes
[results/bundle_validation.json](results/bundle_validation.json) in about
15 CPU minutes, with a checkpoint for each orbit.

- **Clearance.** No two tiles came within 150 m in the four orbits. The patch
  kept at least 153 m, at a corner, and its interior 211 m.
- **Receiver rays.** At every third sample of each service pass, 8,192 rays
  went from the plane through the Moon to points spread over the solar disk.
  Each crossed the bundle where its own ring and both neighbours still carry
  interior tiles; the rings slide about 9 km past each other per orbit, so the
  patch's ends shear apart. None of the 819,200 rays went uncovered. The
  least margin inside a tile's edge was 294 m, and a quarter of the rays
  crossed two or three tiles.
- **Overlaps.** Measured along the Sun's centre, they stayed at least 1,152 m
  along each ring and 579 m between rings.
- **Keeping.** Holding interior tiles to their slots took 4.1×10⁻⁷ m/s²
  (median) and 1.3×10⁻⁶ m/s² (95th percentile), 7% of the kept run's sail
  acceleration at normal incidence. Those tiles stayed within 121 m of their
  slots (40 m median). The patch's ends and edge rings, lit more than a
  complete ring's tiles would be, needed up to 9.6×10⁻⁶ m/s².

Within its four orbits, the patch confirms the kept run's two stand-ins: tiles
held near the copies' positions keep clear and cover every ray through the
interior. The validation has three limits:
- it runs four orbits, while the kept run's edges first came within 150 m at
  day 8.1;
- it is one patch at the middle of the aperture;
- the copies' slots are its target.

## Resources

Film replacement follows O8's 10–20 years. Worn film can be recovered and
remade; propellant is expelled for good. The held screen's columns are its
zoned designs with a 5 or a 10 g/m² annulus. The first ring-fleet column is the
bundle's geometry at 20,000 km with 5 g/m² annulus tiles. The second makes it
nestable: annulus tiles of about 13 g/m², which carry the window's sail
loading, and 23% more rings for the planes' common shift. The third adds
another 23% for the edge strips' recession over a year, in place of steering
the outer rings.

| | Held, 5 g/m² annulus | Held, 10 g/m² annulus | Ring fleet, geometry only | Nestable, outer rings steered | Nestable, edges oversized |
|---|---:|---:|---:|---:|---:|
| Tiles of 10 km | 1.75 million (equivalents) | 1.75 million | 22.4 million | 27.5 million | 32.7 million |
| Optical mass | 1.21 Gt | 2.00 Gt | 23.5 Gt | 45 Gt | 53 Gt |
| Propellant | 51,000 kg/s, 1.6 Gt a year | 86,000 kg/s, 2.7 Gt a year | none with photon keeping | none | none |
| Film replacement | 0.06–0.12 Gt a year | 0.10–0.20 Gt a year | 1.2–2.3 Gt a year | 2.2–4.5 Gt a year | 2.7–5.3 Gt a year |
| Thrust power | 32.8 TW mean, 75.7 TW installed | 55.6 TW mean, 128 TW installed | none | none | none |

At 15,000 and 19,000 km the geometry alone gives 17.4 and 22.2 Gt, but there
the planes' differences need more steering than photons supply (above).

The held screen burns 1.08 kg of propellant per held kilogram each year, about
2% of the Moon's mass over a billion years, and every kilogram added to it,
whether collector or habitat, burns the same. In 30 years the 5 g/m² design
expels 48 Gt, about what the nestable ring fleet holds in orbit. The ring
fleet has 16–19 times the held screen's tiles to make and replace, as recycled
throughput. Done by 30 km/s electric thrust, its keeping and steering would
release 0.4 t/s and 0.06–5.7 t/s of exhaust in lunar orbit, which is why
photon keeping was adopted.

## Electricity

The light is far larger than any demand case: the aperture intercepts about
240 PW. The dimmer's 586 TW of optical power would give 117–176 TW at 20–30%
conversion, and opaque collectors on 1% of the annulus 432–648 TW, at the cost
of about 1% dimming on Earth for a few hours at eclipse-season new moons unless
they turn edge-on then. What collection costs depends on where the collectors
ride.

| Demand | Held screen: generation with its own load | Collectors at 20–30% | Share of the annulus | Propellant if they ride the held screen |
|---:|---:|---:|---:|---:|
| 100 TW | 133 TW | 0.33–0.49 million km² | 0.2–0.3% | 11–17 kg/s per g/m², plus 70–104 kg/s for the light they absorb |
| 300 TW | 333 TW | 0.82–1.22 million km² | 0.5–0.8% | 28–42 kg/s per g/m², plus 174–262 kg/s |
| 1,000 TW | 1,033 TW | 2.5–3.8 million km² | 1.6–2.4% | 87–130 kg/s per g/m², plus 541–812 kg/s |

At 30 g/m² of collector that is 0.4–0.6 t/s more propellant at 100 TW and
3.1–4.7 t/s at 1,000 TW, all of it into the Moon's hemisphere. The ring fleet
needs the demand's collectors alone, 0.24–0.37 million km² for 100 TW, on its
tiles or flying free, with no holding cost. Its own loads are reflectivity
trim and the momentum store, not yet sized for power. Device yield, transfer losses, storage through eclipses and reserves
decide the power delivered to habitat and industry, and none is modelled yet.

## Inhabited infrastructure

Habitats need their own structure, radiation shielding and life support,
none budgeted here. They belong on natural orbits in either design: on the held
screen a habitat would burn about its own mass in propellant every year.

## What the comparison shows so far

The held screen meets the protected radius on its geometry, but as designed
its exhaust breaks the proposed outflow requirement, and its propellant fails
S6. With the September magnets it could meet the outflow requirement only
with gridded ion thrusters and capture of most of its unionized gas, pending a
plasma estimate of the magnetosphere's leak. Its 5 g/m² annulus, the 2 µm
film, adds a third to three fifths of the glow's heat at solar maximum; the
10 g/m² annulus of the 4 µm film adds a tenth to a fifth and raises its
propellant by 70%.

The ring fleet passes no gate yet and fails none, and the four checks of this
round narrow it to one form:
- **Radius.** Near 20,000 km, with rings set from the Moon's orbit plane. At
  15,000 km the planes tear the pattern apart within weeks.
- **Attitude.** Radial-facing tiles, because edge-on turns collide with the
  next ring. A filter centred on each tile, 25% trim and a momentum store of
  about 2×10⁹ N·m·s per tile hold it.
- **Nesting.** One sail loading across the nested bundle, which puts the
  annulus tiles near 13 g/m² and the forced eccentricity at 0.11–0.14.
- **Oversizing.** About 23% more rings for the common shift, and steering,
  a radius profile or more oversizing for the edge strips.

That is about 45 Gt of tiles. Every tile of a patch propagated with its own
forces kept clear and covered every receiver ray through the patch's interior
for four orbits. Its remaining items:
- the edge strips and the bundle's edge rings;
- the keeping in the regressing frame;
- a full fleet's coverage and handovers;
- the hardware of the momentum store.

Both candidates share O1's count of gaps, which the X-ray check tightens, and
the unmapped shadow on Earth.

## The X-ray check

The first of the author's next steps asked where solar X-rays that pass a light
annulus leave their energy, to settle whether the 4 µm film suffices and how O1
should count X-rays. [limb_heat.py](../../../atmosphere/middle_atmosphere/limb_heat.py)
in the atmosphere domain traces the sunlight that passes the ring fleet's
aperture through the Open Moon's tall air along slant paths, for the loss
response's six cases, and finds the state the heated air settles to
([its results](../../../atmosphere/middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07)).
The losses below are molecular escape with Earth's tide; the solar wind and the
exosphere step's losses are left out.

- **The 4 µm film's X-rays add a tenth to a fifth of the glow's heat.** The film
  passes X-rays near 1 nm, which the thermosphere takes only on rays tangent
  within about half a lunar radius of the limb. FISM2's daily spectra put that
  band's rise from quiet Sun to solar maximum at 13–21, where the escape model
  takes 100. At solar maximum the film then adds 2.5–6.2×10⁻⁷ W/m² behind the
  titania stack, against 3.1–3.8×10⁻⁶ W/m² from the sky's Lyman-α glow. With
  no gaps the cooler titania cases lose 1–2×10⁻⁵ kg/s (collisional upper air)
  and 0.003–0.005 kg/s (LTE); the window stack over the annulus would lower
  that to 6×10⁻⁶ and 0.002 kg/s. With all near-infrared heating the glow
  brings the loss to 0.7 kg/s, and the film raises it to 1.2–1.8 kg/s.
- **A gap heats the air 2.5–5 times what O1's count gives it.** Light through a
  gap reaches the air above the limb as well as the disk. Behind the titania
  stack its heat per unit transmission is 2.5–3.2 times the loss response's
  count around the quiescent air, and 3.1–3.9 at the heats where the cooler
  cases lose 1–100 kg/s; behind the 200-nm edge, 3.8–5.2. For the same
  allowance, gaps may pass 0.22–0.32 of what the count allows in the cooler
  titania cases and 0.12–0.27 with all near-infrared heating.
- **At the standard level the cooler titania cases stay under 1 kg/s, and the
  warmest passes 10 kg/s.** At O1's standard level of 2×10⁻⁴ and solar maximum
  the traced loss is 0.003–0.004 kg/s with collisional upper air and
  0.33–0.50 kg/s in LTE. With all near-infrared heating it is 30–36 kg/s, the
  exobase at 3.8–3.9 lunar radii, against the count's 2.3 kg/s; the tight level
  (3×10⁻⁵) gives 2.2–3.2 kg/s there.
- **The loss response's X-ray factors differ for gaps and films.** At solar
  maximum it multiplies a gap's X-rays by the ultraviolet's 2.5 and the titania
  film's by 100. FISM2 gives 4.6–6.9 for the whole band below 10 nm, which
  sets a gap's X-rays, and 11–24 for the films' bands.
- **Where the heat lies is the largest uncertainty left.** The traced light
  lands between the thermal column's 'low' and 'middle' heating shapes in
  depth. The losses above use 'middle', as the loss response does; with 'low'
  they fall 14–550 times. Both counts pass through the same column, so the
  shape moves their losses together.

The tracing uses the ring fleet's aperture at 20,000 km. The held screen's,
farther out, has wider zones and a solar smear of about 420 km at its edges,
and the same physics.

**The author's decisions.** The check leaves the 4 µm film in place. The window stack
over the annulus would lower the losses by a factor of 1.4–3. Behind the
titania stack, at the swarm's tight and standard levels, no case would change
its standing against a 1, 10 or 100 kg/s budget; the glow and the gaps set the
loss. Behind the 200-nm edge it would take the tight level at solar maximum
from 10–15 to about 6 kg/s. On 7 October the author accepted the two
recommendations that followed ([decisions.md](../../decisions.md)):
1. O1 counts gaps by the traced heat per unit transmission in place of the disk
   count, because it is the heat the air takes. The loss response takes it from
   `limb_heat.json`, and O1's levels and the design point were rederived on
   7 October ([below](#the-rederived-design-point)); the allowed transmissions
   fall to about a quarter to a third in the cooler titania cases.
2. The X-rays' solar-maximum factors are FISM2's measured ones (4.6–6.9 for
   the whole band below 10 nm, 11–24 for the films' bands) in place of the escape model's 2.5 for the
   band a filter passes and 100 for films. The atmosphere products that use the
   convention were rerun with them on 7 October.

Giving the thermal column the traced heating shape would settle the remaining
spread.

## The rederived design point

The loss response, its absorption and exosphere steps and the design point now
read the traced heat with FISM2's X-ray factors, at the ring fleet's protected
radius of 4 lunar radii
([requirements.md](../protection_architecture/requirements.md#retention),
[the loss response](../../../atmosphere/loss_response/README.md)). On the
author's order of 7 October two questions of how warm the upper air gets
followed the same day. The design point averages each state over solar cycles 23
and 24, year by year from FISM2's daily record measured band by band. The
thermal column takes the heat where the tracing puts it in height in place of
its fixed 'middle' shape. On 8 October the author had the column radiate in
the infrared, to see whether the runaway latches. The totals below are central
estimates with Earth's tide, behind the titania stack at O1's standard level and
averaged over the cycle unless stated; cycle times are the atmosphere's mass
over the loss.

- **The ultraviolet rises far less over the cycle than the escape model took.**
  At the maxima since 1957 the light from 10 to 175 nm is 1.25–1.64 times the
  quiet week, weighted by energy, where the escape model takes 2.5; the far
  ultraviolet from 122 to 175 nm, most of the energy, rises only 1.07–1.59.
  The glow's Lyman-α and the shortest light rise more than it took: 1.31–1.88
  against 1.5, and below 10 nm 4.6–9.6 against 6.0.
- **The measured solar maximum is milder than the escape model's, and cycle
  19's year harsher.** On 7 October the author had the solar-maximum columns
  take FISM2's measured maxima, the record taken back to 1947: the year around
  cycle 21's maximum of 1979–80, the strongest since FISM2 gained its MgII and
  Lyman-α proxies in 1978, with cycle 19's year, the strongest on record, and
  the escape model's 2.5 as stress cases. With the column's infrared cooling,
  at cycle 21's maximum the warmest titania case loses 1.5 kg/s at the standard
  level without magnets (66 billion years) and 0.090 kg/s with them, against
  1.7 and 0.13 kg/s in cycle 19's year and 1.8 and 0.19 kg/s under the 2.5; the
  cooler cases lose 0.62–0.82 kg/s. Before 1978 FISM2 rests on the 10.7 cm
  radio flux alone. No case runs away at the standard level at any maximum,
  behind either shield. R1's cycle mean, over cycles 23 and 24, is unchanged.
- **Placed where the tracing puts it, the heat warms the upper air far less.**
  The sky's glow, half or more of the heat at the swarm's levels, is absorbed
  within a few e-folds of pressure above the base, and the far ultraviolet
  through gaps in the lower thermosphere; heat laid low is radiated by CO2.
  The same heat leaves the exobase 17–42 K cooler than the 'middle' shape gives.
- **The runaway latches, and the infrared cooling puts it out of reach.** The
  column now radiates CO2's 15 µm band at the middle atmosphere's 400 ppm and
  the base's atomic oxygen and NO, each above what air at the base temperature
  emits ([infrared.py](../../../atmosphere/loss_response/infrared.py)). CO2
  takes most of the heat, so the air needs about three and a half times the
  heat to swell past the shadow's edge. Past a threshold, at 1.9–8.4×10⁻⁵ W/m²
  and an exobase of 4.5–5.8 lunar radii at quiet Sun, the unfiltered light
  beyond the aperture still outgrows the heat that swells the air, so the
  swollen state would outlast the maximum that caused it. At the standard level
  every maximum on record settles a tenth to three-fifths of the way to its
  threshold; at quiet Sun the gaps would have to pass 0.35–0.82% behind the
  titania stack to tip it. The oxygen atoms made above the base, and in the
  swollen state those the unfiltered light would split from the outer air's O2,
  are left out and could move the threshold either way
  ([the loss response](../../../atmosphere/loss_response/README.md)).
- **Every titania case stays within 2 kg/s without magnets.** With
  collisional upper air and in LTE the loss is 0.61–0.80 kg/s with no
  magnetosphere (163 to 124 billion years), almost all of it the solar wind's
  charge exchange and sputtering; with all near-infrared heating it is 1.4 kg/s
  (71 billion years), 1.5 kg/s in its worst year. At the maxima it loses
  1.5–1.8 kg/s. Without the infrared cooling that case lost 3.1–3.5 kg/s over
  the cycle and ran away in cycle 19's year.
- **For 1 kg/s the solar wind decides.** The September magnets bring every case
  far under 1 kg/s: 3×10⁻⁷–10⁻⁴ kg/s in the cooler cases and 0.044 kg/s with all
  near-infrared heating. A dipole of about 3×10¹⁹ A·m² holds the warmest case's
  exosphere to 1 kg/s at every maximum on record, and the collisional and LTE
  cases need none for that.
- **The ring fleet's screen may shelter the Moon from part of the solar
  wind.** The screen absorbs the wind that strikes it, and the wind closes in
  behind an absorbing screen over about eight of its radii. The ring fleet's
  screen, about 20,000 km out, lies under three screen radii from the Moon, so
  its wake may still have an empty core there, about 2.6 lunar radii in radius
  behind a 4-radius screen. Within it the wind's charge exchange falls away: the
  cooler cases lose 0.063–0.066 kg/s with no magnetosphere, and with all
  near-infrared heating 0.41 kg/s. The refill length is an assumption a plasma model has to test. The
  held screen, 68,000–106,000 km out, lies too far for its wake to reach the
  Moon.
- **Atomic oxygen adds little unless the upper air is warm and weakly
  mixed.** The glow and the light through gaps break enough O2 above the base to
  make 10–21 kg/s of oxygen atoms, and with the middle atmosphere's mixing scaled
  for the Moon almost all of it goes back down as odd oxygen: the atoms add
  0.002–0.003 kg/s in the cooler cases and 0.021 kg/s with all near-infrared
  heating (0.007 with the September magnets, which brings that case to about
  0.05 kg/s), about a sixth of what they added before the column's infrared
  cooling. The atoms made near the exobase leave hot and are about 0.001 kg/s
  of this. How strongly the Moon's upper air mixes decides the rest. An
  estimate from breaking gravity waves, calibrated on Earth's measured mixing,
  finds it about as strong as the Moon-scaled profile: the warm case's atoms add
  0.017 kg/s (0.006 with the magnets). Titan's measured mixing, carried to the
  Moon, gives 0.15 kg/s (0.044 with the magnets, for 0.088 kg/s in all, a
  trillion years); Titan's waves are driven by a hundredth of the sunlight, so
  that end may understate the mixing. Earth's own unscaled mixing, 0.46 kg/s,
  lies outside both. Hydrogen adds about 0.0025 kg/s.
- **The relaxed level runs away for the warmest air.** With a 4-radius shadow
  the air runs away once its exobase nears the shadow's edge, from a
  transmission of 0.13–0.82% behind the titania stack by case and activity and
  0.038–0.24% behind the 200-nm edge. The relaxed level's 0.25% holds the
  collisional and LTE titania cases at every measured maximum (1.0 and
  3.4 kg/s over the cycle without magnets), and runs away the warmest at every
  maximum, the LTE case under the escape model's 2.5 and the 200-nm edge in
  every case. The tables find those onsets in air that swells more than air
  heated where the tracing puts it, so they lean early.

The CO2 question this order left for third, how much of CO2's near-infrared
absorption heats the air, now moves the choices little: the cases it separates
all stay within 2 kg/s, and at 1 kg/s the solar wind sets the outcome. The outer
rings come next.

The [primary magnetic architecture comparison](magnetic_architecture.md) still
reads the design point from before the traced count. Its rerun needs the
closure trajectories in `research/runs/solar_shield_array/closure/`, which this
machine lacks, and its pinned-input test fails until then.

## Next steps

The author set this order on 7 October:

1. The X-ray question, in the atmosphere domain: where solar X-rays that pass
   a light annulus deposit their energy along slant paths above the limb. It
   tells whether the chosen 4 µm film suffices and how O1 should count
   X-rays. Done on 7 October ([the X-ray check](#the-x-ray-check)), and its two
   decisions carried through the loss response, O1's levels and the design
   point the same day ([the rederived design point](#the-rederived-design-point)).
   Before the outer rings, the author then had the upper air's warmth settled
   further: the solar-cycle mean and the heat placed where the tracing puts it
   were carried through the loss chain the same day, and they leave the near-
   infrared heating of CO2 little to decide (above). Atomic oxygen's escape, the
   measured solar maxima from 1947 and the upper air's mixing by gravity waves
   followed the same day, and on 8 October the thermal column's infrared
   cooling, which leaves the runaway a latch with its threshold beyond every
   maximum on record.
2. The outer rings: a radius profile that keeps each tilt turning with the
   Sun, with radius steps that grow toward the edges, flown for a year against
   steering by roll with the larger store.
3. The keeping in the regressing frame with the edge rings treated, on a
   bundle with centred filters, and the along-track keeping that holds each
   tile within its overlap.
4. The momentum store and the trim as hardware, and their power.

Noted by the author for later study (7 October): the ring fleet could extend its
protected radius by more active control as the solar weather requires,
spending more energy and propellant while it does; it need not hold its
maximum extent all the time, and its formation is already dynamic. A first
look at the same tiles finds little room. Each strip overlaps its neighbour by
about 615 m of the 9.29 km height step, most of it room for the stack's
half-monthly breathing, up to 5.6% of the step. Holding the breathing by
keeping, which the plan lists as the alternative, would take by estimate about
10⁻⁷ m/s² of out-of-plane sail force per tile, a few percent of the photon
reach; spreading the strips into the freed overlap, less a keeping margin,
would extend the shadow to about 4.1–4.2 lunar radii. The 23% of spare rings
covers the pattern's annual shift of about 1,590 km each way: near the two
times a year the shift passes zero the stack reaches about 4.9 lunar radii on
both sides, and otherwise its extra reach lies on one side. Holding that reach
on demand would mean turning every ring plane with the Sun's 5.1° yearly swing,
about 90 m/s per tile each half year, more than photon forces supply, and with
electric thrust a few hundred megatonnes of exhaust a year, which S7 counts.

For solar maxima the author keeps one dynamic option in view (8 October):
covering failed cells faster through the maximum years, the swarm's tight level,
with spare tiles and servicing ready across the fleet. At cycle 21's maximum it
brings the warmest titania case from 8.4 to 2.5 kg/s, and it keeps cycle 19's
year below the heat at which, in the screening model, that case's air passes the
shadow's edge and holds itself swollen after the Sun quiets
([the loss response](../../../atmosphere/loss_response/README.md)).

For the held screen, the literature can settle two conditions first: how close
thrusters come to ionizing all their propellant, neutralizer flow included,
and, from published comet, ion-release and mini-magnetosphere studies, a bound
on its plasma path with and without the magnets, which needs a model this
repository lacks.

## Running

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.exhaust_isolation
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.photon_control
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.frozen_rings
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.attitude_schemes
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.plane_motion
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.bundle_validation
python -m protection.spectra.annulus_film
python -m atmosphere.middle_atmosphere.fetch_limb_inputs --download
OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.limb_heat
python -m research.studies.solar_shield_array.integrated_ledger
python -m pytest research/studies/solar_shield_array/test_exhaust_isolation.py research/studies/solar_shield_array/test_photon_control.py research/studies/solar_shield_array/test_frozen_rings.py research/studies/solar_shield_array/test_attitude_schemes.py research/studies/solar_shield_array/test_plane_motion.py research/studies/solar_shield_array/test_bundle_validation.py research/studies/solar_shield_array/test_integrated_ledger.py protection/tests/test_annulus_film.py atmosphere/tests/test_limb_heat.py
```

The exhaust study takes about 27 CPU minutes, almost all of it the slow-gas
particles, and keeps one checkpoint per launch date in
`research/runs/solar_shield_array/exhaust_isolation/`. Photon control and the
attitude schemes read the kept run's checkpoints (about two and four CPU
minutes; without them they refly the twelve orbits, about an hour). The frozen
orbits take about 15 CPU minutes. The planes take about 24 CPU minutes, with a
checkpoint after each pass. The patch validation takes about 15 CPU minutes,
with a checkpoint for each orbit. The annulus films and the ledger take
seconds; the films need the protection inputs and the WHI spectrum restored.
The limb tracing takes about 50 CPU minutes and keeps each case's heat tables
in the ignored `atmosphere/middle_atmosphere/cache/limb_heat/`.

## Sources

- [protection/report.md](../../../protection/report.md), section 8: exhaust
  isolation and the energy diagnostic.
- [atmosphere/loss_response](../../../atmosphere/loss_response/README.md):
  returning pickup ions sputter 1–10 molecules each.
- J. A. Young, T. S. Matlock, M. Nakles and M. W. Crofton, *Far Field Plume
  Distribution and Divergence for NEXT: DART Mission*, AIAA SciTech Forum 2019
  ([NTRS 20190004942](https://ntrs.nasa.gov/citations/20190004942)), figure 5.
- D. M. Goebel and I. Katz, *Fundamentals of Electric Propulsion: Ion and Hall
  Thrusters*, JPL Space Science and Technology Series, 2008: chapter 7 (mass
  utilization of SPT thrusters) and chapter 8, figures 8-6 and 8-12 (plume
  current density).
- A. N. Heays, A. D. Bosman and E. F. van Dishoeck, *Photodissociation and
  photoionisation of atoms and molecules of astrophysical interest*, A&A 602,
  A105 (2017), tables 18–19.
- C. Colombo, C. Lücking and C. R. McInnes, *Orbital dynamics of high
  area-to-mass ratio spacecraft with J2 and solar radiation pressure for novel
  Earth observation and communication services*, Acta Astronautica 81, 137–150
  (2012): heliotropic orbits.
- T. N. Woods and others (2009), Geophys. Res. Lett. 36, L01101: the WHI 2008
  reference solar spectrum, quiet-Sun period, below 202 nm.
- O. Coddington and others (2021), *The TSIS-1 Hybrid Solar Reference
  Spectrum*, Geophys. Res. Lett. 48, e2020GL091709: the spectrum above 202 nm.
- P. C. Chamberlin and others (2020), *The Flare Irradiance Spectral
  Model-Version 2 (FISM2)*, Space Weather 18, e2020SW002588: daily spectra
  below 10 nm from LASP's LISIRD, for the X-rays' solar cycle.
- The films' optical constants as in
  [protection/README.md](../../../protection/README.md#short-wave-transmission-2026-09-25):
  CXRO atomic scattering factors (Henke, Gullikson and Davis 1993), fused
  silica (Franta and others 2016) and titania (Siefke and others 2016).
