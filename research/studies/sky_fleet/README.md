# The Open Moon's flyers

The second step toward the distribution of the Open Moon's flyers, after the
[sky-ship study](../sky_ships/README.md) found how large a ship the air carries. This study answers, with light
calculations:
- what flies at each size under a sixth of Earth's gravity, from canopies, micro gliders and wings flown on a
  person's own power to sky boats, sky ferries, regional and long-haul ships, freighters and high platforms;
- how many of each the summit metropolis and its port need, at the busy hour and over a day;
- how they share the air by height;
- which long-haul ship suits the port's crown, with drawings of the choices to scale.

```
OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_fleet.run   # results/sky_fleet.json, about a second
```

The runner is [run.py](run.py), the numbers are in [results/sky_fleet.json](results/sky_fleet.json) and the sources in
[sources.json](sources.json). The models are in [engineering/flight](../../../engineering/flight/README.md): the
buoyant hull fitted to four airships' weight statements, and, new here, a flyer's drag polar (power, sink and glide
against speed), a canopy's descent, and the energy and mass of an electric flyer that takes off vertically. The study
reads the sky-ship study's product (the fitted hull, the flight band's winds, the crown's long-haul traffic), the
summit tower study's (the port's traffic, its canopies, its site, and the frame, terminal and disks of its chosen
form), the design run's air over the whole Moon
([climate/gcm/global_winds.py](../../../climate/gcm/global_winds.py)) and the cloud-resolving equatorial ring's daytime
mixed layer and storms ([climate/crm](../../../climate/crm/README.md)). The design run is main's corrected
`A28_dim5_moon`, the ring is main's corrected equatorial ring (`ring_ring_equator.json`) and the crown takes the
chosen form's frame; the study was rerun on these inputs on 8 October 2026 in the joint integration, and the
[infrastructure review](../infrastructure_review/README.md) records the September figures. The drawings are made in
[visualization/sky-fleet](../../../visualization/sky-fleet/README.md).

**Status.** The author adopted these findings on 2026-09-28 as the planning basis for the Open Moon's flyers
([decisions register](../../decisions.md#sky-ships-and-flyers)); the long-haul class for the crown, whether winged
liners join it, and the sky boats' share of trips, ownership and traffic system stay open.

**Since 28 September.** Against the September figures:
- the high platforms near 70 km meet winds of 25–29 m/s and need about seven times the power to hold station, and
  the 200 m platform's envelope needs 1.2 times today's laminates' strength;
- the day's thermals reach about 10 km, with cloud base near 7 km;
- liners in the band meet a steady westerly of about 11 m/s;
- the crown's long-haul traffic is about a fifth higher in the summit tower study's port programme, so the 750 m
  liner needs nine berths (eight before) and the crown holds 2,300–4,100 t of hydrogen at berth (about 2,100 t
  before). That traffic follows the square-lattice port's programme by height, whose wider corrected base gives
  long-haul travel 61% of the travel floor (53% before); the chosen form's programme, and the crown's berths with
  it, are still to reconcile.

The infrastructure review finds lightning starting at the flight band's height, which bears on the crown's berths
and on the hydrogen ferries in storms.

**Evidence.** First-order sizing of the kind used to compare concepts.
- **What is tested.** The drag polar reproduces a search for least power and least drag, keeps the glide ratio
  under any gravity and scales speed, sink and power with it as theory says; the canopy model gives the summit tower
  study's canopies; the hover model gives the published hover estimate of the Joby S4's 2015 design at a figure of
  merit of 0.58; the long-haul ships here are the sky-ship study's; the berths hold their swing.
- **What is carried over.** Earth's published flyers set the reference performance: gliders' polars, people's
  sustained power, the Joby S4's mass and battery, small airships' payload shares and drag. Lunar gravity enters
  through the physics.
- **What is assumed.** The lunar designs' sizes (a micro glider, a pedalled wing, the electric classes' payloads and
  ranges), the ships' passenger masses and speeds, and the metropolis's travel, a planning case built on Earth's
  transit cities. No flyer has been designed, and the berth layout at the crown is a first arrangement for comparing
  ship sizes.

## The answer

- **Wings and rotors get cheap under a sixth of the gravity; buoyancy stays as it is.** Hovering takes a fifteenth
  of Earth's power, a kilometre on a wing a sixth of the energy, and a glider glides as far 2.46 times more slowly.
  Buoyant lift per cubic metre is Earth's. So the small flyers are winged and rotor-borne, and from ferries up the
  flyers are buoyant ships, which hold still at their berths for free.
- **People fly on their own power.** A pedalled wing of 11 m span and 22 kg flies level on 63 W, within what an
  untrained adult holds for three hours; Earth's human-powered aircraft needed 29–35 m of span and 200–250 W. A
  foot-launched micro glider of 9 m and 4.5 m² glides 15 to 1 at 8 m/s: launched from the port's lowest ring it
  reaches anywhere in the metropolis. A canopy 4.8 m across lands a person at 4 m/s.
- **Sky boats used as cars are light and frugal.** A four-seat electric sky boat weighs 684 kg, hovers on 11 kW,
  cruises at 40 m/s on 4 kW and flies 10 km on 0.74 kWh, a twenty-fourth of what the same trip takes on Earth. A
  buoyant sky boat would be a hull 33–39 m long using three to seven times the energy.
- **Large sky ferries carry a metro's load.** A 250 m hydrogen ferry carries 1,210 people at 25 m/s on 11.4 Wh a
  passenger-km; a line every two to three minutes carries 24,000–36,000 people an hour each way at 42 km/h with its
  stops.
- **The metropolis's sky at the busy hour.** With a Tokyo-like split and sky boats in cars' place, about 290,000 sky
  boats, 33,000 air taxis, 280,000 personal wings and gliders and 890 ferries are aloft over the 3,200 km². All the
  city's flying uses 1.0 GW, half a percent of its power. Stacked in lanes at the spacing NASA's urban air mobility
  analysis uses, the air from 300 m to 1.5 km holds sky boats for 7–8% of trips; a car-like 10% needs closer,
  automated spacing. Owned at Tokyo's rate for cars, sky boats would need berths on three-fifths of the city's land;
  shared, 1.1 million serve the same trips from 65 km².
- **The Moon's tall air stacks the classes over 80 km:** wings and gliders below 300 m, sky boats to 1.5 km, ferries
  to 2.5 km above the ground; gliders soaring to about 10 km by day; regional ships at 22–27 km; long-haul liners and
  freighters in the flight band at 35–45 km; winged liners near 55 km; high platforms near 70 km.
- **The long-haul ship for the crown.** Four liner sizes, 500 m to 1.5 km, all sound, are drawn at the crown to
  scale. The 750 m liner is the recommendation: nine berths on two levels above the terminal, arms of 669 m, a
  ship every 13–22 minutes, 6,990 passengers each (about as many as the largest cruise ship), 19 Wh a passenger-km
  and 257 t of hydrogen a ship. The 500 m liner serves three and a half times the destinations at the cost of 32
  berths and arms of 936 m.
  Winged liners that take off on tilting rotors cross to a mean destination in 5.2 hours for 14 Wh a passenger-km,
  about the largest airships' energy at five times their speed, and are worth adding. Berths in the tower's lee sit in
  its wake, which needs a wind study before the layout is fixed.

## Flight at every size under a sixth of Earth's gravity

Gravity enters each kind of lift differently. For the same flyer in the same air (the models' scaling, checked by
their tests):

| Lift | What lunar gravity does | Moon ÷ Earth |
|---|---|---|
| Gliding on a wing | the glide ratio stays; speed and sink go as the square root of gravity | glide × 1; speed and sink × 0.41 |
| Powered flight on a wing | power at a lift coefficient goes as gravity to 1.5; energy per kilometre at a lift-to-drag ratio as gravity | power × 0.067; energy per km × 0.17 |
| Hovering on rotors | power goes as gravity to 1.5 | × 0.067 |
| Climbing | energy goes as gravity | × 0.17 |
| Descending under a canopy | speed goes as the square root of gravity | × 0.41 |
| Buoyancy | lift per cubic metre stays; the structure sized by weight shrinks | lift × 1 |

- **Wings and rotors get cheap, buoyancy stays as it is.** A buoyant ship needs the same volume per tonne as on
  Earth and pays the same drag at a given speed; only its weight-driven structure shrinks, which counts at the large
  sizes the sky-ship study found. Everything flying on wings or rotors needs a fifteenth of Earth's power to hover
  and a sixth of the energy per kilometre.
- **So the small flyers are winged and the large ones buoyant.** A person flies on their own power; a family's sky
  boat hovers on 11 kW; a buoyant craft with the same payload is a hull 33–39 m long.
  From ferries up, buoyant ships hold still for free at their berths and carry hundreds to tens of thousands of
  people.
- **Earth's flyers, six times larger or six times smaller in area.** The same glider under lunar gravity glides as
  far, 2.46 times more slowly. A glider of the same shape with a sixth of the wing and drag area flies at Earth's
  speed and sink: the Moon's micro gliders.

## The classes

The distribution, from the smallest flyer to the largest. Masses and powers are the models' for the Moon; heights
are above the metropolis's ground unless given above sea level.

| Class | Lift | Size | Mass | Carries | Speed | Power and energy | Where it flies |
|---|---|---|---|---|---|---|---|
| Round canopy | drag | 4.8–5.8 m across laid flat | a few kg | a person | lands at 4 m/s | none | from any rim of the port, any tower |
| Micro glider | wing | 9 m span, 4.5 m² | 18 kg | a person | 7.9 m/s; glides 15 to 1 | none: height and thermals | 30–300 m, thermals to about 10 km by day |
| Personal wing | wing, pedalled | 11 m span, 7 m² | 22 kg | a person | 6.7 m/s | 63 W at the pedals | 30–300 m |
| Personal flyer | rotors and wing, electric | 1 m rotors | 164 kg | one | 30 m/s | 2.6 kW hover, 0.7 kW cruise | 0.3–1.5 km |
| Sky boat | rotors and wing, electric | six 1.7 m rotors | 684 kg | four people | 40 m/s | 11 kW hover, 4 kW cruise; 0.74 kWh for 10 km | 0.3–1.5 km |
| Air taxi or sky van | rotors and wing, electric | eight 2.6 m rotors | 2.1 t | 12 people or 1.2 t | 45 m/s | 34 kW hover, 14 kW cruise | 0.3–1.5 km |
| Parcel and cargo drones | rotors | 0.2–0.8 m rotors | 8–162 kg | 5–100 kg | 20–25 m/s | 0.1–2.6 kW hover | low, in their own lanes |
| Sky ferry | buoyant, hydrogen | 250 m, 42 m across | 210 t lift | 1,210 people | 25 m/s | 1.1 MW; 11.4 Wh a passenger-km | 1.5–2.5 km |
| Regional ship | buoyant, hydrogen | 350 m, 58 m across | 440 t lift | 1,770 seated | 35 m/s | 23 Wh a passenger-km | 20–30 km above sea level; docks 10–15 km above the summit |
| Long-haul liner | buoyant, hydrogen | 500 m – 1.5 km | 1,000–27,000 t lift | 2,000–56,500 in cabins | 30 m/s | 12–27 Wh a passenger-km | flight band, 35–45 km; berths at the crown |
| Winged liner | wing, rotors tilting | 25 m span | 120 t | 300 seated | 145 m/s | 14 Wh a passenger-km; 5.6 MW to hover | cruise near 55 km |
| Freighter | buoyant, hydrogen | 500 m – 1.5 km | 1,000–27,000 t lift | 840–23,000 t | 25 m/s | 12–38 Wh a tonne-km | flight band |
| High platform | pressure hull | 100–200 m | 12–94 t lift | 4–28 t | holds station | 115–750 kW against the wind | near 70 km |

### People on wings: canopies, micro gliders and personal wings

Earth's gliders and human-powered aircraft set the reference. Each is fitted with a parabolic polar that holds its
published best glide and its speed (or its least sink), and the model's other figures are set beside the published
ones:

| Flyer | Span, area, mass | Earth: least sink (published) | Earth: power (published) | The same flyer on the Moon | With a sixth of the wing |
|---|---|---|---|---|---|
| Paraglider, EN-B (MAC PARA Eden 7) | 9.1 m, 19.3 m², 68 kg | 0.95 m/s (1.05) | – | glides at 4.4 m/s, sinks 0.39 m/s | 3.7 m, 3.2 m²: Earth's 10.8 m/s and 0.97 m/s |
| Sailplane, 18 m (HpH 304S Shark) | 18 m, 11.7 m², 600 kg | 0.52 m/s at 445 kg (0.47) | – | 14 m/s, sinks 0.25 m/s | 7.3 m, 1.9 m² |
| Musculair 2 | 19.5 m, 11.7 m², 78 kg | 0.29 m/s (0.27) | 251 W at 10 m/s (250) | 18 W at the pedals, at 4.4 m/s | 7.9 m, 1.9 m²: 37 W |
| MIT Light Eagle | 34.7 m, 31 m², 110 kg | 0.21 m/s | 275 W at 7.8 m/s (about 200) | 19 W at the pedals, at 3.2 m/s | 14 m, 5.1 m²: 39 W |

Earth's aircraft flown on a person's power needed 29–35 m of span and 200–250 W, a trained athlete's output, to cross
the English Channel (Gossamer Albatross, 1979) or to fly 115 km from Crete to Santorini (Daedalus 88, 1988). Under
lunar gravity the same craft needs under 20 W, but flies at 3–4 m/s, slower than the wind. The lunar designs keep the
speed and spend the gain on a smaller wing:

- **The personal wing:** 11 m span, 7 m², 22 kg, pedalled through a propeller. It flies level at 6.7 m/s on 63 W at
  the pedals and glides 22 to 1. Sedentary and recreational adults hold about 70 W for three hours, enough for
  7.2 m/s; an untrained man's critical power (155 W) gives 10.4 m/s or a climb of 0.45 m/s, and a fit young man's hour
  (215 W) 11.7 m/s or 0.73 m/s. On Earth the same wing would need 912 W.
- **The micro glider:** 9 m span, 4.5 m², 18 kg, foot-launched. It glides 15 to 1 at 7.9 m/s and sinks 0.53 m/s; on
  Earth it would need to fly at 19 m/s. A glider launched from the port's rings, 3 km up, reaches 45 km: the whole
  metropolis.
- **Canopies.** A flat circular canopy's drag coefficient is 0.75–0.80 on its cloth area, and it inflates to about
  0.7 of its flat diameter (Knacke). A person of 80 kg lands at 4 m/s under a canopy 4.8 m across laid flat
  (3.3 m inflated) in the air at the tower's foot and 5.8 m (4.0 m) at its top; on Earth it takes 11.4 m. The summit
  tower study's canopies, 3.7–4.5 m across with a coefficient of 1.3, are about a tenth larger than they need be.
  Steerable ram-air canopies glide 3 to 5 (Knacke; Airborne Systems), descend 2.46 times more slowly than on Earth
  at the same loading, and can be a sixth of the area for Earth's speeds.
- **Thermals.** Soaring is a daytime activity. On the ring's land the mixed layer, which thermals fill, deepens from
  0.3 km at night to 9.9 km in the late morning, with its cloud base near 7 km; it is deeper than 2 km for
  335 hours of the lunar day, 138 of them dry, before the afternoon's storms. The climate work takes the Moon's
  convection as Earth's stretched six times in size with the same speeds, so its thermals are six times wider, and
  a lunar glider's turns are six times wider too: 52 m at 7 m/s banked 30°.

### Sky boats, air taxis and drones

Electric flyers that lift off on rotors and cruise on a wing, sized for their trips with the technology of the
Joby S4: an airframe of 36% of the take-off mass (fitted to its certified 4,800 lb with 1,000 lb of payload over its
100-mile target), 235 Wh per kg of battery pack, a lift-to-drag ratio of 14, a figure of merit of 0.6, and motors,
controllers and rotors at 2 kW per kg.

| On a 10 km trip (90 s of hover, a climb to the lane) | Moon: mass | Hover | Cruise | Trip | Earth, the same mission: mass | Hover | Trip |
|---|---|---|---|---|---|---|---|
| Personal flyer, one seat, 60 km range | 164 kg | 2.6 kW | 0.7 kW | 0.14 kWh | 257 kg | 60 kW | 1.9 kWh |
| Sky boat, four seats, 150 km range | 684 kg | 11 kW | 4.0 kW | 0.74 kWh | 1.8 t | 420 kW | 18 kWh |
| Air taxi or sky van, 12 seats or 1.2 t, 200 km | 2.1 t | 34 kW | 14 kW | 2.7 kWh | 9.4 t | 2.2 MW | 100 kWh |
| Parcel drone, 5 kg | 8 kg | 0.1 kW | – | 0.01 kWh | 11 kg | 2.5 kW | 0.07 kWh |
| Cargo drone, 100 kg | 162 kg | 2.6 kW | 0.6 kW | 0.14 kWh | 237 kg | 55 kW | 1.7 kWh |

- **Sky boats used as cars are light and frugal.** The four-seat sky boat carries 32 kg of battery for 150 km. The
  same mission on Earth needs 1.8 t, 554 kg of it battery, and 18 kWh for the same 10 km trip; the Joby S4 aims to
  fly 100 miles with five aboard at 4,800 lb.
- **A buoyant sky boat costs more.** Earth's small airships keep 24–39% of what they lift for payload (Zeppelin NT,
  Skyship 600), so four people in helium need 1,050–1,740 m³: a hull 33–39 m long and 8–10 m across. With the
  drag of Earth's non-rigid ships (NACA Report 397) it uses 210–300 Wh per km at 15 m/s and 380–530 at 20 m/s,
  against 74 Wh per km, hover and climb included, for the electric sky boat at 40 m/s.
- **Slow rotors.** Under a sixth of the weight, the same rotor carries its craft at the same thrust coefficient with
  2.46 times less tip speed. Earth's air taxis already turn their rotors slowly to be quiet: the Joby S4's 2015
  design ran its tips at 370 ft/s against a Robinson R44's 705.

### Sky ferries

Rigid hydrogen ships for the metropolis's backbone, beside metro and rail: the sky-ship study's modern hull in the
city's air (1.16 kg/m³, a kilometre above the summit), cells full at 3 km above it, 150 kg for each passenger with a
bag and a share of the seats and decks, two hours of fuel between refills, and stops every 4 km with two minutes at
each.

| Length | Across | Lift | Passengers | Energy a passenger-km | Line, every 2 / 3 / 5 minutes | Turning radius |
|---|---|---|---|---|---|---|
| 150 m | 25 m | 45 t | 250 | 19 Wh | 7,500 / 5,000 / 3,000 an hour each way | 690 m |
| 200 m | 33 m | 108 t | 610 | 14 Wh | 18,300 / 12,200 / 7,300 | 920 m |
| 250 m | 42 m | 210 t | 1,210 | 11.4 Wh | 36,300 / 24,200 / 14,500 | 1.2 km |
| 300 m | 50 m | 363 t | 2,110 | 9.6 Wh | 63,300 / 42,200 / 25,300 | 1.4 km |

- **Metro capacity at metro speed.** With its stops a ferry averages 42 km/h, about Moscow's metro (41 km/h) and
  above London's Underground (33 km/h). A line of 250 m ferries every two or three minutes carries 24,000–36,000
  people an hour each way; Hong Kong's and Tokyo's busiest lines carry 45,000–71,000.
- **Holding still is free.** A buoyant ferry waits at its stop without power, which makes stations quiet.

### Regional ships

The port's regional docks, 10–15 km above the summit, handle 17,600–29,200 passengers an hour (the summit tower
study's travel floor split by trip). Rigid ships seated at 200 kg a passenger for trips of a few hours, cruising at
35 m/s on routes of up to 1,500 km, with a planning mean of 500 km (four hours):

| Length | Passengers | Energy a passenger-km | Departures an hour | Destinations with a ship every 4 hours / daily | Ships cycling through the port |
|---|---|---|---|---|---|
| 245 m | 580 | 34 Wh | 15–25 | 61–101 / 364–605 | 150–250 |
| 350 m | 1,770 | 23 Wh | 5.0–8.3 | 20–33 / 119–198 | 50–80 |
| 500 m | 5,320 | 16 Wh | 1.7–2.7 | 7–11 / 40–66 | 20–30 |

The 350 m class gives each of 20–33 destinations a ship every four hours, from five to nine docks busy at once.

### Freighters and high platforms

- **Freighters.** Rigid hydrogen ships at 25 m/s on the long-haul route carry 840 t of cargo at 500 m (38 Wh a
  tonne-km), 6,900 t at 1 km (18 Wh) and 23,000 t at 1.5 km (12 Wh). They can set the tower's frame segments
  anywhere up it.
- **High platforms.** Pressure hulls near 70 km, in air of 0.44 kg/m³, for the astronomy and photometry platforms
  of the decisions register: a 100 m hull lifts 12 t and its envelope uses three-fifths of today's laminates'
  strength; a 200 m hull lifts 94 t and its envelope needs 1.2 times that strength. Holding station takes
  115–190 kW and 460–750 kW against the 99th-percentile and highest winds at that height (25 and 29 m/s).

## How many: the metropolis at the busy hour

The planning case: 100 million people within 32 km of the tower. Trips follow Earth's transit cities: 2.8 trips a
person a day (2.0 in Tokyo's region and among London's residents, 3.7 in the Paris region, walking counted), with
the busiest hour carrying 10% of them (8–13%: urban roads' design hour, Hong Kong's morning peak). The split follows
Tokyo's 23 wards (rail 51%, car 8%, bicycle 13%, walking 24%), with sky boats in cars' place:

| Mode | Share of trips | Trips a day | Flight | Aloft at the busy hour | Per km² | Fleet | Energy a day |
|---|---|---|---|---|---|---|---|
| Walking and cycling | 30% | 84 million | – | – | – | – | – |
| Metro and rail | 37.5% | 105 million | – | – | – | – | – |
| Sky ferries | 12.5% | 35 million | 9 km, 13 min | 890 ferries | 0.3 | about 1,000 in service | 5.1 GWh |
| Sky boats | 10% | 28 million | 8 km, 8 min; 1.3 people a flight | 293,000 | 91 | 1.1 million shared, or 25 million owned | 14.8 GWh |
| Air taxis | 2% | 5.6 million | 10 km, 11 min; 3 a flight | 33,000 | 10 | 62,000 | 4.9 GWh |
| Personal wings and gliders | 5% | 14 million | 4 km, 12 min | 279,000 | 87 | 10 million owned | on people's own power |
| Other | 3% | 8.4 million | – | – | – | – | – |

- **Energy.** The city's flying uses about 25 GWh a day, 1.0 GW on average: half a percent of the metropolis's
  200 GW at 2 kW a person. Across the low and high cases, 0.7–1.4 GW.
- **Sky boats are shared.** Owned at Tokyo's rate for cars (222 private cars per 1,000 people; Paris 260), the
  metropolis would keep 25 million sky boats, whose berths at 60 m² each would cover 1,500 km², three-fifths of its
  land. Shared, each flying 20 trips a day, 1.1 million serve the same trips from 65 km² of berths, on roofs and
  façades.
- **The air holds sky boats for 7–8% of trips at today's standards.** Sky boats and air taxis put 101 flyers over
  each square kilometre at the busy hour, 122 m apart on average in six layers between 300 m and 1.5 km. NASA's
  analysis of urban air mobility lays routes 457 m apart and needed 549 m between aircraft, with 61 m vertically, to
  keep the unmitigated collision risk under 10%. Lanes spaced that way, stacked every 61 m from 300 m to 1.5 km, hold
  76 flyers a square kilometre: sky boats for 7.5% of trips. Automated separation at a near miss's distance (152 m)
  holds eleven times more. A car-like share needs closer spacing than today's analyses allow, or lanes stacked
  higher.
- **Wings are dense but slow.** Personal wings and gliders put 87 flyers over each square kilometre, 93 m apart in
  three layers below 300 m, meeting at 7 m/s: the spacing of a busy cycle street, 13 seconds apart.
- **Ferries are few and large.** About 890 are aloft at the busy hour, one to every 3.6 km².

## How the flyers share the air

[visualization/sky-fleet/figures.py](../../../visualization/sky-fleet/figures.py) draws a section to scale through
the metropolis (`build/heights.png` there): the ground, the tower, each class's band, the daytime mixed layer and the
storm tops. The Moon's tall air stacks the classes over 80 km:

| Height | What flies there | Why there |
|---|---|---|
| 30–300 m above the ground | Personal wings, micro gliders, canopies, drones in their own lanes | Short climbs on a person's power; slow flyers kept below the powered traffic; glider and parkour zones designed in |
| 0.3–1.5 km above the ground | Sky boats, air taxis | Climbs of a minute or two; many layers |
| 1.5–2.5 km above the ground | Sky ferries | Their own corridors, descending to stations |
| Up to about 10 km above the ground, by day | Soaring gliders | The cloud-resolving ring's daytime mixed layer over land reaches 9.9 km, with cloud base near 7 km; it is deeper than 2 km for 335 hours of the lunar day, 138 of them before the afternoon's storms; at night it is 0.3 km deep |
| 10–15 km above the summit (22–27 km above sea level) | Regional ships at the port's regional docks | The regional docks the author placed |
| 35–45 km above sea level | Long-haul liners and freighters; the crown's berths | The flight band: above the median storm top (22 km), in air like Earth's at 4.3–5.8 km |
| Near 55 km | Winged liners in cruise | Faster in air of 0.57 kg/m³ |
| Near 70 km | High platforms | Above most storms: the tallest tenth reach 38 km |

- **Storms.** On the ring, storms form over land through about 217 hours of each lunar day's afternoon and evening,
  any place is under rain 5.7% of the time, and a rain spell lasts three hours at the median and twelve at the
  90th percentile. The summit, the rainiest ground in the climate run, sees more. Small flyers land for storms and
  the backbone carries on.
- **Turning.** A lunar glider at 7 m/s banked 30° turns on 52 m, against 18 m for an Earth glider at 10 m/s; the
  Moon's thermals are wider in the same proportion.

## The long-haul ship for the crown

The crown's long-haul traffic is 37,700–62,600 passengers an hour (the sky-ship study, from the port's travellers).
Four sizes of rigid hydrogen liner, all within the sizes the air carries reasonably, cruising at 30 m/s in the flight
band on trips of up to 48 hours with a cabin, services and a share of the crew for each passenger (400 kg). The mean
trip between random points on the Moon is a quarter of its circumference, 2,730 km: 25 hours at 30 m/s.

| | 500 m | 750 m | 1 km | 1.5 km |
|---|---|---|---|---|
| Across | 83 m | 125 m | 167 m | 250 m |
| Passengers | 2,010 | 6,990 | 16,740 | 56,500 |
| Departures an hour | 9.4–15.6 | 2.7–4.5 | 1.1–1.9 | 0.3–0.6 |
| Destinations with a ship every 4 hours / daily | 38–62 / 225–374 | 11–18 / 65–107 | 5–7 / 27–45 | 1–2 / 8–13 |
| Berths at the crown, at the busy end | 32 on 3 levels | 9 on 2 levels | 4 on 2 levels | 2 on 1 level |
| Longest arm | 936 m | 669 m | 238 m | 358 m |
| Hydrogen in one ship | 76 t | 257 t | 610 t | 2,059 t |
| Hydrogen berthed at the crown | 2,440 t | 2,320 t | 2,440 t | 4,120 t |
| Energy a passenger-km, with the hotel load | 27 Wh | 19 Wh | 16 Wh | 12 Wh |
| Load at the nose, weathervaned in the design wind (38.6 m/s) | 0.19 MN | 0.43 MN | 0.76 MN | 1.7 MN |
| Ships cycling through the port | 510–850 | 150–240 | 60–100 | 20–30 |

**The drawings.** [visualization/sky-fleet](../../../visualization/sky-fleet/README.md) draws each class berthed at
the crown from one camera, its berth plan to scale, and all four in elevation beside the Hindenburg
(`build/long-haul-classes.png` and `build/long-haul-classes.html` there), on the tower as the summit tower study sized
its form: 48 diagrid members leaning 2–3° from vertical near the top, the sealed disks as spoked wheels, the 120 m
core, and the eight-storey terminal at the flight band's lower edge (23.4 km above the summit). The berths are a first
arrangement for comparing sizes: ships lie nose-in on the lee side of arms running north and south from the frame
above the terminal, noses into the west wind, which blows there with a steadiness of 0.97–1.00 (the mean vector over
the mean speed, 20–30 km above the summit, in the design run's 3-day means). Side by side they are spaced so that,
swinging together with the wind, their fins meet only past 40°; levels are stacked a hull, its fins and 40 m apart,
from just above the terminal's roof up to the frame's top at 24 km, and lifts would continue up the frame from the
terminal to them. The arms are deep trusses, a seventh of their length deep; they replace the eight arms of 0.9–1.3 km
first drawn, which the metropolis brief found too slender.

**What separates the classes.**
- **Frequency.** Airlines answer growing demand more with frequency than with larger aircraft, and gain more market
  share from frequency than from size (Givoni and Rietveld 2009; Wei and Hansen 2005); airliners of more than 100
  seats averaged 185 seats a flight in 2024. A 500 m liner leaves every four to six minutes and serves a destination
  every four hours for 38–62 of them; a 1.5 km liner leaves every two to three hours.
- **Hydrogen in one ship.** The crown holds 2,300–2,400 t of hydrogen at berth with any class up to 1 km and 4,100 t
  with two 1.5 km liners; one ship holds 76 t at 500 m and 2,059 t at 1.5 km, 4.6 and 123 times the Hindenburg's.
  The transport rules the author set for space traffic cap each packet's energy so that each failure stays bounded;
  the same reasoning favours smaller ships.
- **Energy.** Larger ships use less per passenger: 27 Wh a passenger-km at 500 m, 12 Wh at 1.5 km.
- **The crown.** Larger ships need fewer and shorter arms: 32 berths on arms of up to 936 m at 500 m, 9 on arms of
  669 m at 750 m.
- **Structure.** The hull's share for payload barely changes over these sizes (0.85–0.86); the sky-ship study's best
  length is 750–900 m.

**Speed.** Faster liners cut the mean trip and cost energy and passengers:

| Cruise | Mean trip | 500 m: passengers, energy | 750 m | 1 km | 1.5 km |
|---|---|---|---|---|---|
| 30 m/s | 25 hours | 2,010, 27 Wh | 6,990, 19 Wh | 16,740, 16 Wh | 56,500, 12 Wh |
| 40 m/s | 19 hours | 1,840, 48 Wh | 6,580, 31 Wh | 15,960, 24 Wh | 54,550, 17 Wh |
| 50 m/s | 15 hours | 1,600, 82 Wh | 6,020, 50 Wh | 14,900, 37 Wh | 51,950, 25 Wh |

**Winged liners.** Under a sixth of the gravity a winged liner that takes off and lands on tilting rotors is cheap to
fly. With 300 seats at 400 kg each (the liners' mass per passenger), a lift-to-drag ratio of 17 and a wing loading of
3,000 Pa, it weighs 120 t on a 25 m span, cruises at 145 m/s near 55 km and crosses to a mean destination in
5.2 hours (10.4 at most), for 14 Wh a passenger-km with the same hotel load: about the energy of the 1–1.5 km
airships at five times their speed. It hovers on 5.6 MW for the minute or so of each landing. Carrying half the
crown's long-haul passengers would take 31–52 landing pads with an hour's turn; at 70 m square each, about 0.25 km²
of deck.

**The wake.** Berths in the lee sit in the tower's wake. Above the terminal the frame is an open lattice about
350 m across, its 48 members of 6 m covering about a quarter of its circumference, and the terminal below it is a
solid ring 400 m across; together they shed eddies of about their own width, every few minutes in a steady wind. At
the Empire State Building's mast in 1931 a Navy airship could not moor for the building's eddies, and the terminal was
abandoned; a ship moored to a high mast must be flown all the time it is there (Walker 1975). Before the layout is
fixed, the crown needs a study of its wind, in a computational model or a tunnel; the arms can move the berths out
of the wake.

**Recommendation.** All four classes are sound structure, well within the sizes the air carries. The 750 m liner is
the one the drawings favour: in the September drawings its eight ships on two levels at the crown read as a
composition, at the hull's best size, each carrying about as many people as the largest cruise ship (Icon of the Seas:
7,600 at most, 364 m long); the corrected traffic needs nine at the busy end. It leaves every 13–22 minutes and serves
65–107 destinations daily from arms of 669 m. What it costs against the 500 m liner: under a third of the departures
and of the destinations served every four hours, and 257 t of hydrogen in each ship instead of 76 t. The 500 m liner
is the choice if the long-haul network must reach many cities often; 32 of them berth on arms of up to 936 m, a busy
harbour in the sky. Winged liners are worth adding for fast trips whichever class the arms serve; they would take most
of the travellers in a hurry, and leave the airships the overnight and leisure trips they serve best.

## What this leaves out, and what comes next

- **The author's choices.** The long-haul class for the crown, whether winged liners join it, and how much of the
  metropolis's travel sky boats and personal wings take; each moves the numbers above.
- **The crown's wind.** A computational or tunnel study of the tower's wake at the crown, before the berths are fixed;
  then sizing the arms and the terminal in the tower model.
- **Traffic.** An automated traffic system for sky boats at the densities above: lanes, layers, how craft keep apart,
  and what happens in storms and at night.
- **Designs.** None of these flyers is designed. The next steps would be a sky boat, a personal wing and a ferry at
  the level of Earth's concept studies, and the liners' cabins and hydrogen safety at the size chosen.
- **Moon-wide travel.** The other cities, their ports and the routes between them, which set the long-haul network
  and the freighters' work.
- **Night.** The 354-hour night: no thermals for gliders, lit lanes, and how the metropolis travels through it.
