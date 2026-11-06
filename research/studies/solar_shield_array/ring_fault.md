# The fault across the ring: causes, restart and separation

5 November 2026 · `research/array-industry`. Since 12 October the food and the sky towns take beamed power from the ring fleet, and the power screen ([power.md](../../../provisioning/power.md), finding 9 and "Reserves") calls a fault across the whole ring the design's one common point of failure: every beam stops, 2.5–7.5 TW of critical loads remain, and a week of them, 425–1,270 TWh, comes from stored hydrogen or water through Korolev's head. This study asks what can stop every beam at once, which of those faults can also reach the shade, how long a restart takes, and whether keeping the beams' control apart from the tiles' attitude and station-keeping stops a fault in the power from becoming a fault in the shade. The screen, [ring_fault.py](ring_fault.py), writes [results/ring_fault.json](results/ring_fault.json) in under a second from nine committed products; [test_ring_fault.py](test_ring_fault.py) checks it. Values are marked *product* (a committed product or its note), *computed* (this screen), *published* (a source in [ring_fault_sources.json](ring_fault_sources.json), with how far it was read), *recalled* (general knowledge not checked here) or *design guess*. Design rules and recommendations are proposals.

## The answer

**Ten classes of fault can reach the whole ring, and with the separation three of them can still open the shade: space weather across the fleet, shared navigation and timing, and a collision cascade at the rings' common line of nodes.** The film's UV filtering is passive: its titania absorbs whatever any controller does. The shade's position is active, because the ring fleet keeps station with photon forces, so a tile's attitude is its keeping. A film tile that loses attitude control touches its shingle neighbour within 45 minutes to 2.5 hours on the radial-facing law and 20–45 minutes on the rolled law the design flies (5–9 hours for the heavy reference tile), within 2–4 minutes if its momentum store dumps into it, and within seconds with the rolled law's larger store (*computed*). Without the separation every fault that reaches shared command or software (scheduling, a bad update, an intrusion, an operations error, shared navigation) reaches the tiles' attitude and keeping, and through them the shade. With it, two ways remain to defeat the separation itself: a flaw common to every tile's fallback controller, and an intrusion that breaks the separated command authority.

**The shade fails mechanically in hours; the upper air would answer in weeks to years; the beams' loads need power back in seconds to days.** Free rings come within 150 m of each other in about 2.5 hours at their edges and 6 hours inside at the design loading (*computed* from the bundle *product*). Losing tiles past the latch's thresholds takes 28,000–254,000 of the 28.2 million, which the service fleet replaces in 4.5–41 days at full dispatch, and the whole fleet in 12.5 years (*computed* on *products*). How fast the upper air swells is open: no product gives it, and a bracket on the column's heat capacity puts it at 15 days to 5 years (*computed*, order of magnitude). The climate under the deep air e-folds in 4–7 years (*product*). The latch is the catastrophe because it does not undo itself: the air stays swollen after the heat drops. Just short of the latch it already loses 140–190 kg/s, a cycle time of 520–700 million years, and the energy-limited bound for unshielded air at quiet Sun is 335–80,000 kg/s, cycle times of 1.2–290 million years, against the design's 1–100 kg/s and cycles of 1–100 billion years (*computed* on *product* values). The beams' clock is shorter at the loads: hospitals need power within 10 s, millisecond control and the towns' flight control within a second to a minute, and the food's cultures within minutes to an hour, so the first response belongs on site, as the power screen's first-hour batteries already have it.

**The separation holds against every fault of common control, at a small cost in mass and power and a large one in how slowly the attitude software may change.** With separate hardware, software, command authority and timing, staged updates and local autonomy, faults of the beams' reference, pilots, scheduling, updates, intrusion and operations reach the beams only. An independent controller, sensors and power on every tile at 10–100 kg and 10–100 W (*design guess*) come to 0.28–2.8 Mt and 0.28–2.8 GW over the fleet, at most 6×10⁻⁵ of its 46 Gt (*computed*). The attitude and keeping software may then change only in waves smaller than the correlated loss the shade can bear: 560–930 tiles without an independent fallback, which would take a fleet campaign 83–420 years, so that software is in effect frozen; with a dissimilar fallback that fails at most once in 10³–10⁴ tiles, a campaign takes 4–153 days (*computed*). Two physical couplings remain for the separation to answer. A transmitting tile's beams push it with 21–31 N in use, about a film tile's whole sail force, so a ring-wide beam stop is a simultaneous force step on every transmitter; and a whole-tile transmitter shares the rings' one sail loading only near 0.21 kg/m², 12–19 times lighter than the 1978 reference's (*computed*). The transmitters should fly as their own vehicles outside the shade's rings, or carry that step in keeping independent of the beams.

**With the separation, the week of stores covers the restart of every class but two.** The beams restart within a quarter of a day to three days after faults of their reference, pilots, scheduling, updates and operations, and a grid-wide black start adds up to a day (Iberia's grid took 12–24 hours in April 2025, *published*), so these classes need 1.25–4 days. An extreme solar proton event lasting six days needs 8 days, an intrusion up to 15, and a collision cascade outlasts any store (*design guesses* anchored on precedents). Local supplies carry part of the ground's critical loads indefinitely: rivers, tides and the heat flow 5–30%, and roofs with panels and hydrogen 20–100% (*computed* on the energy study's *product*). The towns carry their own critical loads on their reserve lift for 5–19 years (*computed* on the power screen). The recommendation is 8 days in the ground's store, 61–181 TWh more than the week, with local supplies and a faster, hardware-rooted restart from intrusion closing the rest.

## What the project had

| Source | What it gives this study | Evidence state |
|---|---|---|
| Register (main, read-only) | The ring fleet as lead candidate with photon keeping, one sail loading, the 4-lunar-radius shadow, the two-sided radius layout; formations of replaceable units; beams chosen 12 October; the whole system engineered for extreme uptime, millisecond control kept on the Moon (12 October); budgets with cycle times | Decisions and directions |
| [power.md](../../../provisioning/power.md), [power.py](../../../provisioning/power.py) | Critical loads by kind, the week's stores, the towns' hydrogen, single-tile losses (1–7 transmitting tiles a day, 10 minutes for a fresh beam), 3,541–25,270 whole-tile transmitters, their radio power | *Screen* on *design guesses*; the ring fault's causes and restart are its next work 5 |
| [Domain synthesis](../domain_synthesis/README.md) | The chosen distribution; retrodirective beams that focus only on a receiver's pilot | Chosen; screen arithmetic |
| [README](README.md), [constraints review](constraints_review.md) | The fleet as designed: 28.2 million 10 km tiles in 1,859 nested rings near 20,000 km; the annulus allocation 13.765 g/m²; photon authority near its limit at the edge rings | *Products* and estimates |
| [fleet.md](fleet.md), [pilot.md](pilot.md), [active.md](active.md) | Earlier architectures' coverage and conditional correction budgets in m/s a day | Superseded for the ring fleet; continuous coverage unresolved |
| [relative_orbits.md](relative_orbits.md), [ring_bundle.json](results/ring_bundle.json), [ring_keeping.json](results/ring_keeping.json) | Free rings first within 150 m at hours 5.4 (edge) and 12.3 (interior); keeping costs 0.117 m/s per ring per orbit on average and 0.254 at the top edge; the edge ring takes 15% more light | *Product*: 19 rings at 15,000 km, 62.7 g/m² tiles |
| [photon_control.py](photon_control.py) | Keeping demand 5.4×10⁻⁷ m/s² median, 5.7×10⁻⁶ peak; radial-facing tiles are gravity-gradient unstable | *Product*, kept run |
| [attitude_schemes.py](attitude_schemes.py) | Radiation torque on a balanced-filter tile 37,000 N·m median, 80,600 N·m 95th percentile; stores 1.9×10⁹ N·m·s (chosen) and 2.3×10¹⁰ (rolled); turn times that fix the inertia; pitch clearance of a few degrees | *Product*, kept run |
| [tile_design.md](tile_design.md) | A failed tile's times to the 150 m gate and contact; 89.1% of the aperture under one tile; 6,184 dispatches a day; the reference tile's inertia, 389 times a uniform plate's | *Product*; no validated tile |
| [ring_layout.py](ring_layout.py), [rolled layout](results/ring_layout_rolled.json) | The radius at which each height's ring turns with the Sun; rings far apart in the stack meet only at the common line of nodes; nesting needs eccentricity keeping | *Product* |
| [tacking.md](tacking.md), [electromagnetic.md](electromagnetic.md), [magnetic_architecture.md](magnetic_architecture.md), [return_repair.md](return_repair.md) | Literature and earlier control studies; photon recoil at 100 MW (0.46 N); no common-mode, timing or space-weather analysis | Context only |
| [array_industry.md](array_industry.md) | Trim devices over 25–89% of each tile change the redirected band's reflection by up to a quarter; the planned trim holds its state without power; rotors of 0.1–1.7 Gt | *Product* and plan |
| [protection/README.md](../../../protection/README.md), [transmission](../../../protection/transmission/README.md), [requirements](../protection_architecture/requirements.md) | The swarm's standard level, 1.9–2.2×10⁻⁴; a failed cell covered within a day to a week (O8) | *Product*; design levels are assumptions to test |
| [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | The latch: gaps passing 0.35–0.82% at quiet Sun behind the titania stack, 0.12–0.24% behind the 200-nm edge, tip the upper air; it stays swollen | *Product* (screening model); no response time |
| `provisioning/energy_options.md` on domain/provisioning (afe38eda), read-only | Rivers 0.22–0.44 TW, tides 0.03–0.06, the heat flow 0.05–0.12, roofs 0.8–6.2 TW; local supplies belong in the ring-fault reserve | *Product* on another branch |

The fleet's continuous coverage stays unresolved in these studies; this one takes the fleet as the products describe it.

## Findings

1. **The film filters and shades by itself; only the shade's position, and part of its dimming, depend on control** (*computed* on *product* text).
   - **Filter.** The UV transmission comes from uncovered aperture, since the intact film passes at most 1.7×10⁻⁴³ below 175 nm. The trim changes reflection in the redirected band and never the titania, so no trim state passes ultraviolet.
   - **Dimming.** If a trim device traded reflection for transmission, its full range would pass 0.29–1.29% of the sunlight to the Moon; if it traded reflection for absorption, the films would add 0.5–3.7 W/m² of glow. Devices that switch between specular and diffuse reflection put neither at stake. Both effects work on the climate, which answers over years.
   - **Position.** The tiles hold their slots with tilts and trim, so the shade's geometry is only as good as attitude and keeping control.

2. **A tile that loses attitude control reaches its neighbour in minutes to hours** (*computed* on the attitude and tile *products*). The shingle stands 222.5 m off its neighbour, so a corner meets the 150 m gate at 0.59° of pitch error (0.18° inside the tile study's 50 m envelope, which reproduces its 182 s at 0.001°/s) and the neighbour's surface at 1.8°.

| Tile | Inertia | Drift under the light's torque, to gate / to contact | Momentum store dumped into the tile (chosen / rolled store) | Wrong command at the 10⁶ N·m allocation |
|---|---:|---|---|---|
| Annulus allocation, 13.8 g/m² | 1.15×10¹³ kg·m² | 28–42 min / 49–73 min | 63 s / 5 s to the gate; 194 s / 16 s to contact | 8 min / 14 min |
| Mean tile, 16.4 g/m² | 1.37×10¹³ kg·m² | 31–46 min / 54–80 min | 75 s / 6 s; 231 s / 19 s | 9 min / 15 min |
| Kept run, 62.7 g/m² | 5.22×10¹³ kg·m² | 58–87 min / 1.7–2.5 h | 4.8 min / 23 s; 15 min / 72 s | 17 min / 30 min |
| Reference tile, 2.68 kg/m² with its frame | 4.46×10¹⁵ kg·m² | 3.3–5.3 h / 5.8–9.3 h | 6.8 h / 33 min; 21 h / 1.7 h | 2.7 h / 4.7 h |

   The torque is the balanced filter's radiation torque, 37,000–82,000 N·m (median to 95th percentile), with gravity's scaled by inertia to 20,000 km; the reference tile's gravity-gradient torque dominates. The edge-on turn times of the attitude study fix the light tiles' inertia as a uniform plate's to 10⁻⁴. A drive that brakes its rotor into the tile is the fastest path for light tiles, so the store's drives have to fail free. On the rolled law the design flies, the controller applies 1.0–4.9×10⁵ N·m to a light tile just to follow the roll (the attitude study's torque, scaled by inertia), so a rolled tile whose control stops departs from its commanded attitude with that torque missing and reaches the gate in 12–25 minutes and contact in 20–44, taking the shingle's angles as a stand-in for its clearance.

3. **Without keeping, free rings close on each other within hours and the stack's nesting comes apart over days** (*product*, scaled *computed*). In the bundle with sail force and mutual shadows, the top edge pair comes within 150 m at hour 5.4 and the interior at hour 12.3, at 62.7 g/m² and 15,000 km; at the design's 13.8 g/m² a steady differential force closes the same gap in 2.5 and 5.7 hours. The edge ring takes 15% more light than the interior, the asymmetry that drives it. A single tile among held neighbours, left with the keeping demand scaled to the design loading (2.5×10⁻⁶ m/s² median, 2.6×10⁻⁵ peak), crosses the gate in 25–82 minutes and touches in 70–227 minutes (the tile study's runs at 20,000 km, interpolated). Left on their forced eccentricities the rings cannot be nested at all (ring layout *product*). The kept run's cost agrees with the demand used here: 0.117 m/s per ring per orbit spread over the 45.8-hour orbit is 7.1×10⁻⁷ m/s², 1.3 times the photon study's median demand (*computed* on *products*).

4. **The rings' common line of nodes concentrates the collision hazard** (*computed* on the rolled layout *product*). The two-sided layout ranks the 1,859 rings by the radius at which each turns with the Sun, so 62% of radius neighbours come from opposite sides of the stack. At the nodes those neighbours pass 1 km apart at 160 m/s (median), 357 m/s (90th percentile) and up to 389 m/s, 60% of them faster than 100 m/s; strip neighbours pass at 0.23–0.25 m/s, and rings slide 0.037 m/s for each kilometre of radius step. A fragment whose orbit sweeps the stack's radii meets about 12 tiles an orbit near the nodes, two of them faster than 50 m/s (averaged over crossing phases). Whether fragments then multiply depends on how films, frames and rotors break at 10–400 m/s, which no product gives.

5. **Tens of thousands of lost tiles open the shade past the latch** (*computed* on the tile, swarm and loss *products*). The standard level already takes 1.9–2.2×10⁻⁴, 16–18% of the tightest latch transmission.

| Shield | Latch at quiet Sun | Tiles lost at once, contiguous / spread | Share of the fleet | Full dispatch to replace them |
|---|---|---|---:|---|
| Behind the 200-nm edge | 0.12–0.24% | 27,800–62,300 / 31,200–69,900 | 0.10–0.22% | 4.5–11.3 days |
| Behind the titania stack | 0.35–0.82% | 92,700–226,000 / 104,000–254,000 | 0.33–0.80% | 15–41 days |

   Contiguous losses uncover everything under them; spread losses uncover the 89.1% of each slot under one tile. Rebuilding the whole fleet at 6,184 tiles a day takes 12.5 years. The latch thresholds lean early, the loss product notes, because its tables were traced through air that swells more than the traced heating gives.

6. **No product gives the upper air's response time; a heat-capacity bracket puts it at 15 days to 5 years** (*computed*, order of magnitude). The column above the 0.3 Pa base holds 0.18 kg/m²; the air within three or four e-folds of pressure above the base stays near the base temperature, so the share that warms, e⁻⁴ to e⁻³, holds 3.5–9.6 J/m² per kelvin. Warming it 30–300 K at the threshold heats of 1.9–8.4×10⁻⁵ W/m² takes 15 days to 4.8 years. A time-dependent thermal column would replace this bracket.

7. **The beams push their transmitters, so a beam stop is a force step on the fleet** (*computed* on the power and attitude *products*). Radio power over the speed of light gives 21–31 N per whole-tile transmitter in use and 62–93 N at its rating, 0.73–1.09 times a film tile's median sail force of 28 N. Against the 444 N that sunlight puts on a whole tile absorbing it all, a ring-wide stop changes each transmitter's push by 4.7–7% in use and 14–21% at rating, the size of the edge ring's 15% light excess that closes free edges in hours. That collector shares the film tiles' acceleration only near 0.21 kg/m²; at the 1978 reference's 10 kg/kW it weighs 2.7–4.0 kg/m², 12–19 times more, and cannot hold a slot in rings of one sail loading by photons alone.

8. **Concentrated transmitters cannot open the shade even if all are lost; transmitters spread over a tenth of ten times the tiles can** (*computed*). The 3,541–25,270 whole-tile transmitters uncover 0.011–0.080% of the aperture if every one is lost, under the smallest latch margin of 0.098%; with two damaged neighbours apiece, 0.034–0.24%. The tenth-tile variant's 35,410–252,700 transmitters uncover 0.11–0.80% with no collateral, at or past the 200-nm edge's latch.

9. **The critical loads need power back in seconds to days, so the store serves the restart and on-site supplies serve the first minutes** (*published* where named, else *design guess*).

| Load | Critical share, TW | Power back within | Basis |
|---|---|---|---|
| Close computing, millisecond control | 1.16 | 0.01–1 s | *design guess*: its own uninterruptible supply |
| Health care (within residents) | part of 0.80–2.39 | 10 s | *published*: NFPA 110 Level 1, Type 10; Class 96 holds 96 hours of fuel |
| The towns' control | 0.06–0.55 | 1–60 s | *design guess*: on the towns' own batteries and fuel cells |
| Food's living cultures | 0.29–1.17 | 1–60 min | *design guess*: dense aerobic cultures use their oxygen in minutes (*recalled*) |
| Refrigeration of food and medicine (within residents) | part of 0.80–2.39 | 4–48 h | *published*: USDA, 4 h closed refrigerator, 48 h full freezer |
| Water and sewage pumps (within residents) | part of 0.80–2.39 | 2–24 h | *design guess*: tanks and reservoirs |
| Water lift | 0.001–0.028 | 1–3 days | *design guess*: reservoirs |
| Lamp-lit farms | 0.22–2.23 | 2–14 days | *design guess*: crops kept dark with air and water moving |

10. **A Halloween-class storm today puts satellites out of service at 2–130 times the rate per tile the shade can bear from one cause** (*computed* on *published* counts). The largest number of tiles one cause may leave uncontrolled at once, with two to four damaged neighbours each, is 5,600–9,300 behind the 200-nm edge and 18,500–30,900 behind the titania stack: a chance per tile and per event of 2.0–3.3×10⁻⁴ and 6.6–11×10⁻⁴, a tenth of that with a margin of ten. The storms of October–November 2003 put 2 of roughly 800–1,000 satellites out of service (0.20–0.25%, the denominator *recalled*), 1.8–2.3 times the loosest of these budgets and up to 130 times the tightest, and gave 28 anomalies (2.8–3.5%). The tiles must ride through such an event holding their slots, in a hardened safe mode or passively, far better than the satellites of 2003.

## The fault classes

How each function fails is given without the separation: it *stops*, *holds* its last state or *drifts* under a wrong reference. Restart is the beams' restart in days, without and with the separation, before the grid's black start. Classes, restarts and likelihoods are *design guesses* anchored on the precedents named (ring_fault_sources.json).

| Class | Reaches, without → with the separation | Spread and detection | Beams / attitude and keeping / film | Restart, days | Likelihood and precedents |
|---|---|---|---|---|---|
| F1 The beams' phase and timing reference | beams → beams | Seconds; pilot lock and spot efficiency fall | stops / – / – | 0.1–7 → 0.01–0.25 | Occasional; Galileo's timing facility, 2019 (about six days); GPS time, 2016 (12 h) |
| F2 The receivers' pilots | beams → beams | A lost pilot defocuses its beam at once; receivers report | stops / – / – | 0.1–3 → 0.1–3 | Occasional; a shared design across 0.7 million crowns and up to 48,000 patches |
| F3 Scheduling and command software | beams, attitude → beams | Beams hold their assignments and fade as tiles set, 2–20 h; watchdogs | holds / holds / – | 0.01–1 → 0.01–0.5 | Expected; Starlink's control plane, 2025 (2.5 h) |
| F4 A bad update or latent timer | beams, attitude → beams | As fast as the rollout, or everywhere at once; canaries, a dissimilar stack | stops / drifts / – | 0.5–10 → 0.1–3 | Expected unless no update reaches the fleet at once; CrowdStrike, 2024 (reverted in 78 min, days to repair 8.5 million machines); GLONASS ephemerides, 2014 (11 h) |
| F5 An intrusion into command | beams, attitude → beams | At the attacker's pace once trusted access is held; authentication, rate limits | stops / drifts / – | 3–60 → 1–14 | Occasional to expected; Viasat KA-SAT, 2022 (about 30,000 modems replaced) |
| F6 Space weather across the fleet | beams, attitude → beams, attitude | Simultaneous: solar protons bring upsets and latch-ups within hours of a flare and stay days; charging in Earth's magnetotail, which the Moon crosses a few days each month (*recalled*); radiation monitors | stops / holds / – | 0.5–7 → 0.5–7 | Expected: Oct-1989-class events a few times a solar cycle; TDRS-1's 243 attitude-memory upsets (1989); Anik E1/E2's wheels (1994); Galaxy 15 (2010) |
| F7 Shared navigation and timing for attitude and keeping | beams, attitude → attitude | As fast as the broadcast; cross-checks against neighbours' ranging, stars and the Sun | stops / drifts / – | 0.5–7 → 0 | Occasional: GLONASS 2014, GPS 2016, Galileo 2019 in six years among four systems |
| F8 A collision cascade | all three → all three | Fragments return to the nodes every half orbit, 35 h; the tiles' ranging | stops / stops / stops | 180–3,650 | No precedent at this scale; the concern of crowded low Earth orbit (Kessler and Cour-Palais 1978, *recalled*) |
| F9 A meteoroid stream | film → film | Hours, fleet-wide; forecast | – / – / drifts | 0 | Expected yearly; a day of storm holes 6×10⁻⁸–6×10⁻⁶ of the aperture, at most 0.65% of the smallest latch margin (*computed*) |
| F10 An operations error | beams, attitude → beams | As fast as a command is accepted; two-person and two-system authority | stops / holds / – | 0.01–2 → 0.01–1.5 | Expected; Iberia's black start, 2025 (12–24 h); Starlink, 2025 |

Hardware that fails tile by tile stays with the power screen's single-tile losses; a common design flaw that wears out with age arrives over months, where trend monitoring, staggered production and dissimilar suppliers answer it.

## The shade's clock

Once attitude and keeping control is lost across the ring, its sequence runs in four stages (*computed* on *products*):
- **Minutes.** Stores that dump into their tiles bring film tiles to contact in 12 s to 4 minutes; wrong commands at the torque allocation in 12–30 minutes.
- **Hours.** Drift under the light's torque brings film tiles to contact in 45 minutes to 2.5 hours, and in 20–44 minutes on the rolled law (5–9 hours for the heavy reference tile); free edge rings come within 150 m in about 2.5 hours, the interior in about 6.
- **Days.** Contacts spread through the stack; radius neighbours meet at the nodes at up to 390 m/s, and the nesting comes apart as the rings relax toward their forced eccentricities.
- **Weeks to years.** Damage beyond 28,000–254,000 tiles carries the UV transmission past the latch; the upper air would swell over 15 days to 5 years by the bracket; the climate warms over 4–7 years as the dimming goes; the swollen air loses more than the 140–190 kg/s found just short of the latch, up to the unshielded energy-limited bound of 80,000 kg/s, cycle times of 1.2–700 million years.

A ring-wide loss of attitude and keeping that lasts beyond the first hour or two therefore brings tiles into contact across the fleet faster than any service can follow, and rebuilding the fleet takes 12.5 years, longer than every bracket on the upper air. A correlated loss of up to the latch's count is different. Spares kept in the formation and the service fleet's 6,184 dispatches a day replace the 200-nm edge's count in 4.5–11 days, inside the bracket's lower end of 15 days, and the titania stack's in 15–41 days, at that end. That sizes a recovery requirement: spares in the formation of order the edge's latch count, 28,000–62,000 tiles (the tile study's local-spares case already holds 140,000 replacements in transit at its peak, *product*).

All of this stays inside the 4-lunar-radius design. Recovering from a swollen state is not shown within it and stays open.

## The beams' clock and the stores

The power screen's stores run for 7 days, after the US NRC's practice for safety-related emergency power on Earth (Regulatory Guide 1.137). Each class below needs the upper end of its restart plus up to a day for the grid's black start (*computed*; restarts are *design guesses*).

| Class | Days without / with the separation | Energy with the separation, TWh | Ground store after local flows, TWh | Within the week |
|---|---|---|---|---|
| F1 Reference | 8 / 1.25 | 76–226 | 44–164 | yes |
| F2 Pilots | 4 / 4 | 243–723 | 139–526 | yes |
| F3 Scheduling | 2 / 1.5 | 91–271 | 52–197 | yes |
| F4 Bad update | 11 / 4 | 243–723 | 139–526 | yes |
| F5 Intrusion | 61 / 15 | 910–2,710 | 522–1,970 | no |
| F6 Space weather | 8 / 8 | 485–1,450 | 279–1,050 | no, by a day |
| F7 Shared navigation | 8 / 0 | 0 | 0 | yes |
| F8 Collision cascade | 3,651 / 3,651 | 221,000–660,000 | | no |
| F9 Meteoroid stream | 0 / 0 | 0 | 0 | yes |
| F10 Operations error | 3 / 2.5 | 152–452 | 87–329 | yes |

Energies run from the low-cost end's 2.53 TW to the high-cost end's 7.53 TW of critical load. The ground's column subtracts rivers, tides and the heat flow (0.30–0.62 TW), pairing the least load with the most local supply and the most load with the least. Those flows carry 5–30% of the ground's critical loads (2.07–5.78 TW); roofs and pavements carrying panels stored in hydrogen raise that to 20–100% (*computed* on the energy study's *product*). The towns carry theirs on 0.10–0.39% of their reserve lift a week, enough for 253–974 weeks.

## The separation

These design rules are proposals.

1. **Hardware.** Each tile carries an attitude and keeping controller of its own, with its own sensors, a few square metres of photovoltaics and a battery, on a bus no beam electronics share. The beams' collectors and arrays never power attitude.
2. **The transmitters.** Transmitters fly as their own vehicles outside the shade's nested rings, so that their recoil, their mass and any fault of theirs never touch a shade tile. Where a transmitter stays in a ring, it matches the rings' sail loading (about 0.2 kg/m² for a whole-tile collector) and its keeping carries the recoil step with the beams on or off. Whole-tile transmitters stay preferred over the tenth-tile variant (finding 8).
3. **Software.** The attitude and keeping software is small, verified and dissimilar: two stacks of different design, processors and suppliers, with a simple hardware arbiter that hands control to the untouched one within seconds.
4. **Command authority.** The beams' command system has no write path to attitude or keeping. Fleet-wide attitude moves need several independent authorities, and every tile enforces in hardware a rate limit and a veto on any command that brings it toward its neighbours.
5. **Timing and navigation.** The beams take their phase from the receivers' pilots and their own oscillators; attitude and keeping take their geometry from ranging and sighting of neighbours (precedents: Starling's swarm navigation at 0.6–1.3% of range, CanX-4/5's sub-metre control, *published*), with absolute references only as slow corrections that neighbours cross-check. No shared clock is trusted unchecked.
6. **Updates.** Nothing is updated across the fleet at once. Waves stay below the correlated loss the shade can bear: 560–930 tiles without an independent fallback, so in practice that software is frozen; with the dissimilar fallback failing at most once in 10³–10⁴ tiles, 0.56–9.3 million tiles a wave and a campaign of 4–153 days. Only one stack changes at a time, after canaries and a soak of one to three days.
7. **Local autonomy.** Each tile holds its slot from its own relative sensing with no fleet-wide input for longer than any shared service's restart, at least a week, and through a six-day solar proton event in a hardened safe mode that keeps attitude.
8. **Actuators that fail safe.** Momentum-store drives fail free and cannot return stored momentum faster than a set rate; trim devices hold their state without power, as the planned branch already intends; trim trades specular for diffuse reflection, keeping the dimming independent of trim state (finding 1).
9. **Beams that fail safe.** Retrodirective beams focus only on a pilot and spread when it stops (US DOE and NASA 1978); a beam douses when its path is obstructed (Dickinson's 1978 NASA patent); a transmitter that loses command stops radiating, unlike Galaxy 15, which relayed for nine months while uncommandable (*published*).

**The cost** (*computed* on *design guesses*): 0.28–2.8 Mt and 0.28–2.8 GW over the fleet, 6×10⁻⁶–6×10⁻⁵ of a mean tile's 1,640 t and of the 46 Gt inventory, beside the trim's 0.006–1 TW. The larger costs are complexity and pace: two software stacks verified to independence, several command authorities, attitude software that changes over months at best, and transmitters with vehicles and orbits of their own, which this screen leaves unsized.

**What remains after it:** space weather, which reaches every tile's electronics together; shared navigation and timing, bounded by relative sensing but present wherever tiles trust a common reference; the fallback's and arbiter's own design, common to every tile; and debris at the nodes.

## What remains catastrophic

The faults that can still open the shade across the ring are an extreme solar particle or charging event that leaves more than 5,600–30,900 tiles uncontrolled at once; a latent flaw common to the fallback controllers; a corruption of shared navigation or timing that relative sensing fails to bound; an intrusion that defeats the separated command authority; and a collision cascade at the nodes.

For a chance of 1% over the horizon (*design guess*), such faults may occur at most 10⁻⁵ a year over a thousand years, 10⁻⁶ over ten thousand, 10⁻⁸ over a million and 10⁻¹¹ over a billion, the protection requirements' life (S2). Aviation holds a catastrophic failure condition to 10⁻⁹ per flight hour, about 9×10⁻⁶ a year of continuous operation, and nuclear power its large releases to 10⁻⁶ per reactor-year (*published*). Beyond ten thousand years the target falls below the most demanding targets engineering sets itself, so the design has to exclude these faults by physical margin: tiles that hold their slots through any common event; spares in the formation that cover a latch-sized loss within the upper air's response; dissimilar control; and a stack whose nodes stay clear.

The ring-wide classes the design can exclude from the shade are F1–F5 and F10, the faults of the beams' common control, which the separation confines to the beams, and F9, a meteoroid storm, which cannot approach the latch. F7 can only be bounded, by relative sensing.

## Assumptions and Earth values

- **Earth values leaned on.** Restart times of Earth's software services, satellite navigation systems, a satellite network under intrusion and a national grid, used as precedents for the Moon's systems; the Halloween 2003 satellite counts, from orbits inside or near Earth's magnetosphere, while the ring lies outside it most of each month and sees solar protons unshielded; NFPA 110 and USDA cold holding; the US NRC's 7-day fuel practice behind the power screen's week; the FAA and NRC probability targets, used as comparisons. No lunar value here is bounded by an Earth limit.
- **Torques and clocks.** The attitude study's torques for radial-facing balanced-filter tiles and for the rolled law on the kept run at 15,000 km; the radiation torque does not change with areal mass, and the gravity-gradient and attitude torques go with inertia (gravity's also with the cube of the mean motion). Free-ring times scale from a 19-ring bundle with 0.6 km steps as the inverse square root of the acceleration. The rolled tile's clearance takes the shingle's angles as a stand-in.
- **Coverage.** The tile study's ideal lattice for the single-cover share; contiguous losses as the conservative case.
- **Nodes.** Each height's ring at the radius where it turns with the Sun, interpolated from 17 points; the matched layout shifts and stretches these radii and keeps their order, which the test confirms on the 13 flown rings.
- **Design guesses.** Restart times, tolerances, the separation's 10–100 kg and 10–100 W a tile, two to four damaged neighbours per uncontrolled tile, a margin of ten, a soak of one to three days, fallback failure of 10⁻³–10⁻⁴, an accepted chance of 1%, a whole-tile transmitter at 0.2 kg/m², and the meteoroid storm's tenfold to thousandfold flux for a day (*recalled*).

## Limits

- The screen runs no dynamics; contact, tearing and fragmentation are not modelled, so the time from first contacts to a latch-sized loss is bracketed only by its stages.
- The upper air's response is an order-of-magnitude bracket on an assumed warming; the products' thermal column is steady.
- Fragment encounters average over crossing phases, while a real fragment's crossings lock to its orbit for many orbits.
- The local supplies come from domain/provisioning (afe38eda), outside this tree; the screen carries their values with the product's hash.
- Most sources were read as search summaries, as the sources file records.

## Open questions, with recommendations

1. **The upper air's response time.** Recommendation: a time-dependent thermal column in the atmosphere domain, run on the latch's thresholds; it decides whether spares can outrun the latch.
2. **How the stack breaks at 10–400 m/s.** Recommendation: a fragmentation screen of films, frames and rotors, and a node-crossing model of the fragments.
3. **The transmitters' vehicles and orbits.** Recommendation: fold into the power screen's next work 2, with the recoil and sail loading above.
4. **The tiles' safe mode through an extreme proton event.** Recommendation: size the hardened electronics and a passive attitude hold against the October 1989 and Halloween 2003 environments seen outside Earth's magnetosphere.
5. **The days of storage.** Recommendation: 8 days in the ground's store, with local supplies and a hardware-rooted restart path for intrusion covering up to 15.
6. **Spares in the formation.** Recommendation: of order 28,000–62,000 tiles, the 200-nm edge's latch count, placed for dispatch within days.

## How to run it

```
python -m research.studies.solar_shield_array.ring_fault      # writes results/ring_fault.json in under a second
python -m pytest research/studies/solar_shield_array/test_ring_fault.py
```

The product (schema `terluna.research.ring-fault/1`) binds its producer and its ten inputs by sha256, with its constants, evidence statement, reading rule and design. The tests check the product's currency, the inertia against the attitude study's turn times, the gate angle against the tile study's lost guarantee, the scaling of drift and free-ring times, the rolled law's faster departure, the kept run's cost against the photon study's demand, the latch counts against margin, cover and dispatch, the node speeds and the matched layout's ring order, the recoil, the stores against the power screen's week, the fault table's completeness, the upper-air bracket, the meteoroid storm, the swollen state's cycle times, the transmitters' matched loading, the Halloween comparison, and that every precedent the fault table cites is in the sources file with its reading.
