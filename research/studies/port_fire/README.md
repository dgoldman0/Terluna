# Fire safety in the central port

The [central port](../summit_tower/README.md#the-central-port) at 0.16 g: how smoke and heat move through its spaces,
when detectors and sprinklers answer, how long people have and how long they need, and what the tower's shafts,
water, sealed zone and firefighting from the air ask. The author asked for four pieces of CPU-light work on
2026-09-27:
1. scaling Earth's fire-engineering formulas to lunar gravity and each zone's air;
2. design fires and a rule for materials;
3. getting out;
4. the tower's systems.

```
OPENBLAS_NUM_THREADS=1 python -m research.studies.port_fire.run   # results/port_fire.json, a few seconds
```

The models are in [engineering/fire](../../../engineering/fire/README.md), the runner is [run.py](run.py), the numbers
are in [results/port_fire.json](results/port_fire.json) and the sources in [sources.json](sources.json).

**Since 28 September.** The [infrastructure review](../infrastructure_review/README.md) of 8 October reruns this
study with main's corrected climate, whose winds aloft are about twice as strong. The wind-fed fires in the upper
zones grow: a stores compartment burning out in the sealed zone at the median wind releases 818 MW for 1.2 h in place
of 484 MW for 2.1 h. Smoke, sprinklers and getting out stay as below.

**Evidence.** This is first-order fire-engineering arithmetic.
- **What is tested.** The correlations are Earth's, and the model reproduces them in Earth air (its tests).
- **What is extrapolated.** Froude scaling carries them to lunar gravity. That scaling is sound for buoyant flow, and
  untested at partial gravity beyond small samples.
- **What is assumed.** The lunar design fires are a bracket. Crowd flows, the time people take to start moving, and the
  aerial evacuation figures are assumptions.

Nothing here is a fire strategy or a certified design.

## What carries over from Earth

Fire engineering's formulas for plumes, ceiling jets, smoke filling, flashover and vent flows follow Froude scaling:
spaces of the same size behave alike when their fires' dimensionless heat release is the same. So in a given space, a
fire at lunar gravity behaves as an Earth fire several times larger, played about 2.5 times more slowly.

| | Earth | Open floors, 0–3 km | Enclosed, 3–10 km | Enclosed, 10–15 km | Sealed, 15–24 km |
|---|---|---|---|---|---|
| Air | 1 atm, 20 °C, 20.95% O₂ | 0.92 atm, 10 °C (outside air) | 0.83 atm, 20 °C | 0.74 atm, 20 °C | 0.95 atm, 20 °C (pressurised) |
| A fire acts as an Earth fire of | 1× | 2.7× | 3.0× | 3.3× | 2.6× |
| Slowed by | 1 | 2.46 | 2.46 | 2.46 | 2.46 |
| Flame height, same fire | 1 | 1.9 | 2.0 | 2.1 | 1.9 |
| Air drawn into the plume at 2 m, same fire | 1 | 0.69 | 0.65 | 0.62 | 0.68 |
| Heat release that flashes a room over | 1 | 0.63 | 0.58 | 0.55 | 0.62 |
| Largest fire a room's openings can feed | 1 | 0.32 | 0.28 | 0.25 | 0.32 |

The oxygen and the pressure add to what gravity does:
- **Flames about twice as tall.** Air with 17.5% oxygen gives a fifth less heat per kilogram burned, so a flame
  draws in more of it, and the flames lengthen.
- **Hotter, slower smoke.** The plume takes in about two-thirds of Earth's air, so its smoke is hotter and less
  diluted, and it moves more slowly.
- **Starved room fires.** A room's door and window breathe about a third as much, since buoyant flow through an
  opening goes as the square root of gravity. Its hot layer holds its heat and flashes over at a smaller fire, and
  a fully developed fire is starved of air.

## 1. Smoke, heat and sensors in the port's spaces

The spaces:
- **A room.** A hotel or short-stay room of 30 m² under a 3.0 m ceiling.
- **A compartment.** A floor compartment of 2,000 m² under a 3.5 m ceiling, the size of Earth's usual smoke reservoir.
- **A hall.** A hall of 5,000 m² through three of a band's four storeys, under a 12 m ceiling.

The sensors are smoke detectors, judged by the usual 13 K rise of the ceiling jet, and quick-response sprinklers (RTI
50 (m s)^½, 68 °C). The table compares Earth with the enclosed floors at 3–10 km for EN 1991-1-2's design fire, with
the slower lunar estimate (section 2) in brackets. The other zones differ by less than a tenth.

| Space and fire | | Smoke detector | Sprinkler opens | Smoke at head height, no sprinklers | Smoke at head height, sprinklers | Flame height at sprinkler |
|---|---|---|---|---|---|---|
| Room, hotel fire (medium growth) | Earth | 54 s | 168 s, 315 kW | 59 s | 59 s | 1.0 m |
| | Moon | 46 s (75 s) | 138 s, 210 kW (220 s, 135 kW) | 76 s (104 s) | 76 s (104 s) | 2.4 m (1.8 m) |
| Compartment, offices (medium) | Earth | 120 s | 204 s, 462 kW | 9.2 min | 17.0 min | 1.2 m |
| | Moon | 96 s (161 s) | 164 s, 299 kW (266 s, 196 kW) | 8.9 min (14.7) | 27.0 min (35.4) | 2.7 m (2.0 m) |
| Compartment, shops (fast) | Earth | 68 s | 123 s, 670 kW | 6.2 min | 13.4 min | 1.3 m |
| | Moon | 60 s (96 s) | 104 s, 485 kW (164 s, 299 kW) | 4.9 min (8.9) | 19.5 min (27.0) | 3.2 m (2.4 m) |
| Hall, events (fast, 500 kW/m²) | Earth | 131 s | 242 s, 2.6 MW | 15.1 min | 43 min | 2.8 m |
| | Moon | 105 s (176 s) | 183 s, 1.5 MW (309 s, 1.1 MW) | 14.7 min (23.9) | 74 min (beyond 90) | 5.5 m (4.4 m) |
| Hall, terminal (slow) | Earth | 439 s | 837 s, 1.9 MW | 33 min | 55 min | 1.6 m |
| | Moon | 307 s (553 s) | 550 s, 0.84 MW (1,018 s, 0.72 MW) | 36 min (54) | beyond 90 min | 3.8 m (3.1 m) |

What this shows:
- **Sprinklers open sooner, on smaller fires.** The ceiling jet is slower but hotter, and a sprinkler bulb heats
  faster in a hotter jet. They open on fires 0.4–0.7 times the size they would on Earth. Detectors judged by
  temperature also answer sooner.
- **Point detectors lag.** At detection the ceiling jet moves at 0.2–0.8 m/s, against 0.5–2 m/s on Earth, so smoke
  takes longer to enter them. Aspirating detection draws air in through pipes by fan, so it does not wait on the
  jet. The July 2026 notes' "fast atmospheric sampling" is the same idea.
- **With sprinklers, smoke takes longer to reach people.** In compartments and halls it takes about 1.5 times as long
  to come down to head height as on Earth, 20–27 minutes in compartments and over an hour in halls.
- **Without sprinklers, a hot layer governs.** With Earth's fast fires, a compartment becomes untenable a little sooner
  on the Moon (4.9 against 6.2 minutes for shops). The less diluted smoke reaches 200 °C overhead before it reaches
  heads. With the slower lunar estimate it takes longer than on Earth (8.9 minutes).
- **Flames reach the ceiling.** The fires the sprinklers hold stand 2.7–3.6 m tall in 3.5 m compartments, against
  1.2–1.8 m on Earth, and 2.4 m in a room with a 3 m ceiling, against 1.0 m. Ceilings need to be noncombustible, and
  sprinklers must reach flames that touch the ceiling.
- **Smoke control is smaller or the same.** Holding the smoke above 2 m against a sprinkler-held fire takes 1.5–2.2
  m³/s of exhaust per compartment, against 2.4–3.1 on Earth. Natural vents need about the same area, 1.6–2.1 m²
  against 1.2–1.4, because the hotter layer drives its own venting harder.

**The room a fire starts in.** A room flashes over at 1.7 MW on the Moon, against 2.9 MW on Earth.
- **The same fire flashes it over sooner:** at 392 seconds against 514. The slower lunar estimate takes 784 seconds.
- **Sprinklers hold it far below that,** at 8–14% of the flashover size.
- **Its occupants have a minute or a little more** before smoke reaches their heads, as on Earth. They depend on the
  room's own alarm and sprinkler.
- **A fully developed room fire is starved:** 2.1–2.7 MW against Earth's 7.5 MW, burning 69–89 minutes against 25.

## 2. Design fires and materials

EN 1991-1-2 Annex E gives each occupancy a growth rate, a peak heat release per square metre and a fire load. The
lunar estimate takes one growth class slower, twice the time to reach 1 MW, and scales the peak by 0.65. That factor
is Peatross and Beyler's measured fall in burning rate at 17.5% oxygen.

Why a slower estimate:
- **Measured slowing.** In partial-gravity flights, upward flame spread over thin fuels fell roughly in proportion
  to gravity (Feier et al. 2002).
- **Less heat to the fuel.** A fire's heat feedback to its fuel falls with weaker buoyancy.

The design keeps Earth's fires as the other bracket.

| Port use | Stands in | Earth: time to 1 MW, peak, fire load | Lunar estimate |
|---|---|---|---|
| Travel: terminals, docks, lounges | transport (public space) | 600 s, 250 kW/m², 122 MJ/m² | 1,200 s, 162 kW/m² |
| Hotels and short stays | hotel (room) | 300 s, 250 kW/m², 377 MJ/m² | 600 s, 162 kW/m² |
| Commerce and trade | office | 300 s, 250 kW/m², 511 MJ/m² | 600 s, 162 kW/m² |
| Shops, markets, food | shopping centre | 150 s, 250 kW/m², 730 MJ/m² | 300 s, 162 kW/m² |
| Events, conferences | theatre (cinema) | 150 s, 500 kW/m², 365 MJ/m² | 300 s, 325 kW/m² |
| Services and stores | library (the Annex's heaviest load) | 150 s, 500 kW/m², 1,824 MJ/m² | 300 s, 325 kW/m² |

**Which materials burn.** Lunar gravity lets materials burn in less oxygen. The measured shifts are from drop-tower
centrifuge tests (Ferkul and Olson 2011), tabled by Miller et al.:

| Material | Limit at 1 g (upward test) | At lunar gravity | Shift |
|---|---|---|---|
| Mylar G film (70 kPa) | 20.0–21.2% O₂ | 14.1–15.6% | 5.6–5.9 points |
| Ultem 1000 (70 kPa) | 23.0–23.5% | 19.9–21.0% | 2.5–3.1 points |
| Nomex HT90-40 (101 kPa) | 22.1–23.5% | 19.9–21.0% | 2.2–2.5 points |

- **Rule for critical materials.** Materials in escape routes, concealed spaces, cable insulation and the sealed
  zone stop burning in NASA-STD-6001's upward test at 23.5% oxygen. That is the design's 17.5% plus the largest
  measured shift. Testing at 1 atm adds margin, since the port's air is at 0.74–0.95 atm and a material's oxygen
  limit rises as pressure falls (NASA White Sands tests).
- **Everywhere else.** Rooms, shops and halls follow Earth practice for linings and furnishings, with sprinklers.
- **Top up by pressure.** Topping up the middle zones' air should raise its pressure and keep the 17.5%. Enriching
  it instead to match the summit's breathing would take 21.3% oxygen at the top of the 3–10 km zone, 23.6% at 10–15
  km and 28.3% at 24 km, and 27–36% to match Earth's sea level. Those fractions reach or pass the 23.5% at which
  critical materials must stop burning, and NASA's planned lunar habitat air (34% at 56.5 kPa) is the condition the
  fire literature flags.

## 3. Getting out

**The floor plan.**
- **Plates.** Floor plates are strips or rings about 40 m deep with open air on both faces, so no point is more than
  20 m from a rim.
- **Compartments.** Each is 2,000 m², a 50 m length of strip, with two double doors into its neighbours. A walk to a
  door is at most about 45 m.
- **Moving people.** People first move sideways into the next compartment, which serves as a refuge as it does in
  Earth's hospitals. The rim is the way out when a whole band must be cleared.

**Walking at 0.16 g.** People walk as fast as on Earth: at lunar gravity the change from walking to running comes at
1.42 m/s (De Witt et al. 2014). Their grip for braking and turning is a sixth of Earth's, so crowd flow through doors
is taken as 0.6–1.0 of Earth's. The time people take to start moving after the alarm is assumed at 1–3 minutes where
they know the place and 2–5 where they are strangers.

The table gives the people in the compartment or hall where a fire starts, when it is full. It compares the time
they need (RSET: detection, alarm, starting, moving) with the time until smoke is untenable, on Earth and at 3–10 km,
for Earth's design fire.

| Use and space | People | RSET, Earth / Moon | Untenable with sprinklers | Margin with sprinklers | Margin without |
|---|---|---|---|---|---|
| Terminal hall | 167 | 11–14 / 8–12 min | 55 min / beyond 90 | +41 min / — | +20 / +24 min |
| Commerce compartment | 167 | 4–6 / 4–6 min | 17 / 27 min | +11 / +21 min | +3.1 / +2.8 min |
| Shops compartment | 500 | 5–8 / 5–9 min | 13 / 20 min | +5.2 / +10.5 min | −2.0 / −4.2 min |
| Events hall | 1,667 | 6–9 / 6–10 min | 43 / 74 min | +33 / +64 min | +5.8 / +4.6 min |
| Stores compartment | 40 | 3–5 / 3–6 min | 13 / 20 min | +8.1 / +14 min | +0.9 / −0.7 min |

- **Sprinklers are the difference.** With sprinklers, every space leaves a clear margin, larger on the Moon than on
  Earth. Without them, fast fires in shops and stores beat the people on both worlds, as they do in Earth codes that
  require sprinklers there.
- **Hotels and short stays.** Occupants of the room a fire starts in depend on its own alarm and sprinkler (section
  1). Everyone else is kept apart from the fire by its walls and door.

**Clearing a whole band by air.** This is needed only when a band itself is at risk. The planning figures are a berth
every 60 m of rim, ships of 100 people, and a ship at each berth every 5 minutes.
- **A full-size band below 15 km.** It has 400,000 m² of floor, 10 km of rim and 166 berths, and lifts off 55 people a
  second.
- **Clearing times.** A full band clears in 4 minutes for hotels, 10 for commerce or short stays, 30 for shops and 40
  for events packed to capacity.
- **The sealed zone.** Its bands are smaller: 74,000 m² at 19.5 km, with 30 berths. They clear in the same times.
- **Canopies and wings.** These are the backup for people able to use them.

## 4. The tower's systems

**Shafts.** With 20 °C inside, the stack effect builds 6–12 Pa per 100 m of shaft, against 134 Pa per 100 m in an
Earth tower in winter.
- **Stairs and lift lobbies.** They are held at least 12.5 Pa above the fire floor (NFPA 92), and a person can open a
  door against at most about 85 Pa (NFPA 101's 133 N).
- **How far a shaft can run.** An unbroken shaft can run 620–1,200 m before the stack pressure alone fills that
  window, against 54 m on Earth. Stairs and lift shafts can therefore run through about three bands between breaks.

**Water.** At lunar gravity a water column's pressure grows 16 kPa per 100 m, against 98 kPa on Earth.
- **Pressure zones.** Standard sprinkler pipework (1.21 MPa) serves 745 m of height per pressure zone, against 123 m
  on Earth, so 33 zones cover the tower.
- **Tanks.** Each band's tank holds 186 m³, NFPA 13's ordinary-hazard sprinkler and hose demand for 90 minutes.
- **Lifting cost.** Pumping water to the top costs 10.8 kWh per m³.

**The sealed zone.** It is pressurised with ordinary air to the summit's 0.95 atm, and its envelope carries 25–37 kPa.
- **If the envelope fails,** the air inside matches Earth at 4,290–5,670 m. The FAA gives 20–30 minutes of effective
  performance at 5,500 m (18,000 ft), and half that after a rapid decompression.
- **So the zone needs** pressurised refuges or drop-down masks, as airliners have.
- **Its materials** follow the critical-materials rule and airliner cabin practice (FAA heat release limits).

**Firefighting from the air.** With stations every 3 km up the tower, crews reach any band within 135 seconds,
including a minute to launch. NFPA 1710 allows 240 seconds of travel for the first engine on Earth. Crews enter at the
rim, and the water comes from the band's own tank.

**Fully developed fires and the frame.** A floor compartment can burn out behind failed windows along both faces.
Here, at 3–10 km, EN 1991-1-2's fire load burns as follows:

| Compartment | Fire load | Earth, still air | Moon, still air | Moon, median wind (2.5 m/s) | Moon, 90th-percentile wind (4.3 m/s) |
|---|---|---|---|---|---|
| Offices | 1.0 TJ | 198 MW, 1.4 h | 56 MW, 5.1 h | 263 MW, 1.1 h | 460 MW, 0.6 h |
| Shops | 1.5 TJ | 198 MW, 2.0 h | 56 MW, 7.3 h | 263 MW, 1.5 h | 460 MW, 0.9 h |
| Stores | 3.6 TJ | 198 MW, 5.1 h | 56 MW, 18 h | 263 MW, 3.8 h | 460 MW, 2.2 h |

- **Still air.** The fire burns a third as hard for three to four times as long.
- **Wind.** It feeds a fire at any gravity. The port's modest winds, with medians of 2–4 m/s in the ring model, drive
  through a floor whose faces have failed and give it more air than Earth's buoyancy would. This is the wind-driven
  fire of Earth's tall buildings. In NIST's tests, a 9–11 m/s wind drove gases above 400 °C through the corridor and
  stairs of a seven-storey building.
- **So the frame's steel within a band** needs protection rated for both cases: the large wind-fed fire, and the long
  fire in calm air.
- **Windows stay whole** in the fire compartment as long as possible, and there are cross-walls that stop the wind
  running through a floor. Together these keep a compartment fire in the calm case.

## What is left out

- **Fire at partial gravity at room scale.** No fire of room size has burned at partial gravity. Samples small enough
  for 20–30 s parabolic flights and 5 s drop towers are all there is, and NASA's FM2 is planned as the first fire on
  the Moon. The design fires stay a bracket until full-scale tests at 0.16 g exist. The July 2026 notes already had
  oxygenation wait for them.
- **The scaling's own limits.** Froude scaling carries buoyant flow. Radiation, soot, heat lost to walls and laminar
  flow at low Grashof numbers do not scale with it. The zone model keeps 70% of the convective heat in the layer, so
  it runs hot for low ceilings.
- **Sprinkler sprays at 0.16 g.** Drops fall more slowly and spread wider, and how a spray reaches and cools a fire
  has not been tested.
- **People.** Crowd flows at 0.16 g and the time people take to start moving are assumptions.
- **Smoke outside the port.** Wind-driven smoke through failed facades and plumes rising to the band above are not
  modelled. Fire carried between bands by brands, and fires on ships at berth, are not modelled either.
- **Next.** CFD runs (FDS) of a band section with its wind would test the scaling and the wind cases. The same work
  then extends to the megacity.
