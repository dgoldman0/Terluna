# The summit metropolis

A beautiful, sustainable mega-metropolis centred on the Moon's largest sky port, on the far side's summit
highland. This brief collects what the author decided on 2026-09-27 (entered in the
[decisions register](../../research/decisions.md#the-summit-port-and-metropolis)), the first estimates behind
those choices, and the design notes still open. The port itself is sized in the
[summit tower study](../../research/studies/summit_tower/README.md) and its fire safety in the
[port fire study](../../research/studies/port_fire/README.md).

The density, ring and heat figures below are first estimates worked by hand from those studies' products. No
model of the metropolis exists yet.

**Since 28 September.** The figures below rest on main's corrected design run, A28_dim5_moon, and on the summit
tower, port fire, sky-ship and sky-fleet studies as rerun on it on 8 October 2026 in the joint integration. The
[infrastructure review](../../research/studies/infrastructure_review/README.md) of 8 October records the September
figures and reads this brief against later work:
- the electrified storms and their heights;
- storm protection, with lightning that is rare and very large, and flashes that start at the crown's height;
- the night, with about 82 hours of practical dusk after sunset and about 190 dark hours;
- the coasts, and the infrastructure off the Moon.

## The site

- **The summit.** The atlas's highest cell, at 5.4° N, 158.6° W, stands 11.6 km above sea level on a gently
  rounded rise several hundred kilometres across. The ground falls about 3 km within 100 km.
- **Water.** The drainage estimate fills the craters and basins with rain-fed lakes. They cover 24% of the ground
  within 25 km of the summit, 41% within 100 km and 57% within 300 km, against 10.3% of the whole Moon.
  The Korolev lake (141,000 km², up to 6.7 km deep) has its shore 84 km away and 3.7 km below the summit.
- **Rain.** The summit's highland is the rainiest ground in the corrected design run: the four cells round the
  summit rain 9–17 mm a day, 3–6 m a year, and the rainiest cells nearby 19–20, against 4.8 mm a day for
  equatorial land. Rain follows ground height across the equatorial belt (correlation about 0.7). The magnitude is
  uncertain: the run's cells are 170 km across, and the CM1 ring along the equator, which is flat and at sea level,
  rains about four-fifths as much as the GCM over equatorial land. Afternoon storms reach 22 km at the median.
- **Air.** At the summit the free air averages 12 °C and 0.95 atm (1.16 kg/m³). The air temperatures in this brief
  are the GCM's. CM1 and the GCM disagree on how warm and humid the air over land is near the ground, and the
  climate work is paused, with comfort carried as the range between them
  ([decisions register](../../research/decisions.md#research-practice)).
- **Heritage.** No site in the conservation register lies within 1,400 km; the nearest is Chang'e 6.
- **Sky.** Earth never rises on the far side, and the 354-hour night has no earthshine.

## The port and its tower

The port is commercial (trade, travel, hospitality and short stays, with no homes and no industry) and holds
about 1.5 million people on a regular basis. Long-haul sky ships dock in the flight band, so the tower rises
24 km, to 35.6 km above sea level, on a base 8.2 km wide. Regional docks sit at 10–15 km, and every block of
floors has berths.

- **Frame.** A round diagrid gathers into six splayed legs, as on the Eiffel Tower, through a deep transfer ring a
  kilometre or two up. The six lift spines run down the legs. The flared legs carry the wind's overturning, as
  Eiffel's profile does, and the ground between them stays open for the city. In the hand-sized form (below) each
  of the six feet stands on about 40,000 m² at 1 MPa.
- **Floors.** Rings run round the frame in the lower tower, and disks round a central core from about 12 km up, as
  in the author's picture of the port (chosen on 2026-09-28, replacing full rings). Each is a block of a few
  storeys with a park and open-air shops on its roof. On the first estimate's levels, the eight below 12 km
  (3.1–11.5 km) are rings and the twelve from 12.3 to 22.2 km are disks, 7.0 km² at the lowest and 0.2 km² at the
  top. A disk holds far more floor than a ring, so fewer, more widely spaced disks may hold the programme.
- **Wide lower rings.** The lower rings, where the open air is mildest, are wide, with fewer storeys under the
  park. At 150 m instead of the first estimate's 40 m, ring 0's park grows from 0.83 km², a strip 21 km round, to
  about 3 km², close to Central Park's size. Their decks cantilever from the frame. In lunar gravity a 150 m
  cantilever carrying 3 t/m² (a park, two storeys and structure) bends like a 61 m one on Earth, about the 67 m
  the Marina Bay Sands SkyPark cantilevers; its trusses taper, deep at the frame and thin at the inner edge. The
  bending at the root, 21 km round ring 0, goes into the frame, whose members need stiffening for a node level
  above and below: 1.5 Mt of steel for the six rings in the hand-sized form, under 2% of its 98 Mt. Cables from
  the frame above would need about a tenth of that.
- **Crown.** The sealed long-haul terminal sits in the flight band; the summit tower study's sized form makes it
  eight storeys round the frame's top at 23.4 km. Its docking arms, as first drawn, cantilever 0.9–1.3 km at about
  90 times their depth, which is not credible structure. The
  [sky-fleet study](../../research/studies/sky_fleet/README.md#the-long-haul-ship-for-the-crown) sizes the arms for
  each long-haul class instead: deep trusses running north and south from the frame above the
  terminal, with ships nose-in on their lee side. For the crown's 37,700–62,600 long-haul passengers an hour it
  compares 32 berths of 500 m liners on three levels (arms up to 936 m), 9 of 750 m on two (669 m), 4 of 1 km
  (238 m) and two of 1.5 km, drawn to scale in visualization/sky-fleet. It recommends the 750 m liner, with winged
  liners added for fast trips; the choice is open. Berths in the lee sit in the tower's wake, which needs a wind
  study before they are fixed.

The summit tower study sizes this form by hand in [form.py](../../research/studies/summit_tower/form.py)
(commit aac8c64): the diagrid's members, the legs and arches, the transfer ring and a buried tie between the feet,
six rings of cantilevered trusses from 3.1 to 10.7 km, seven disks from 12.3 to 21.4 km as spoked wheels round the
core, and the terminal. At the corrected design gust of 32.4 m/s it takes 98 Mt of steel, with a first period of
99 s and 37 m of sway at the top in the service gust. It is concept sizing by hand calculation, and a structural
analysis has yet to check it. The study's square lattice, with a band every 200 m, stays as a comparison: 46 Mt of
frame steel on a 9.3 km base.

## Parks and air by height

With the design's oxygen (17.5%), the air at the top of each zone compares with Earth's at these heights:

| Above the summit | Like on Earth | Mean air temperature | Open air |
|---|---|---|---|
| 0–3 km | 1,900–2,400 m (Mexico City) | 12 to 9 °C | Parks and open-air shops for everyone |
| 3–10 km | up to 3,500 m (La Paz) | down to 3 °C | Parks for residents and acclimatised visitors: cool alpine gardens |
| 10–15 km | up to 4,300 m (El Alto) | down to −1 °C | A few high gardens for short visits, each beside a pressurised refuge; winter parks on the disks from about 13.5 km |
| 15–24 km | up to 5,700 m | down to −9 °C | None: the blocks are sealed |

The summit is the Moon's highest ground, and the air cools about 0.9 °C per kilometre of height, a seventh of
Earth's rate, since that rate scales with gravity. So the tower's foot averages about 12 °C against the lowlands'
21–24 °C, and the metropolis runs from there to about 16 °C on its lowest ground, 4 km lower. The disks begin at
12.3 km, where the air averages about 1 °C, so none of their open decks is temperate. The mild open-air parks are
the ground and the lower rings: about 9 °C on ring 0 and 6–8 °C on rings 1 and 2.

Every interior is kept at a comfortable pressure and temperature, so enclosed floors take topped-up air from
about 3 km up. High parks close for the afternoon storms and their lightning.

**Winter.** The corrected design run keeps ground-level air at about 19–22 °C from the equator to the poles, so
the tower's upper levels hold about the only winter on the Moon, and it lasts all year, since the Moon has no
seasons. Each height also keeps its temperature through the lunar day: in the run the air at the tower's heights
changes by only 0.1–0.4 °C between its warmest and coolest hours (run A28_dim5_moon, model years 20–29, read at the
site; the site product does not yet carry temperature by hour angle). The mean air falls through freezing at about 13.9 km above
the summit, so the open winter parks sit on the disks from about 13.5 to 15 km, at about 0 to −1 °C, and the high
gardens below them, from 10 km, are raw and near freezing. Above 15 km, in the sealed blocks, the air is −1 to
−9 °C. The storms bring wet snow, sleet and graupel, and night cloud leaves rime.
Snowmaking and chilled rinks give reliable snow and ice. At a sixth of Earth's gravity snowflakes fall at under
half Earth's speed and snow lies fluffier; skiing works on Earth's slope angles, more slowly.

## The metropolis

**Size.** About 100 million people at about 40,000 per km², within about 32 km of the tower: a rough planning
scale, chosen on 2026-09-27. At a sixth
of Earth's gravity height stops limiting density: the same steel reaches 2.5–3 times as high, the winds are light,
and sky boats reach upper floors, which shrinks the lift cores that cap Earth's towers. A green district whose
towers cover a fifth of the ground and average 30 storeys, at about 70 m² of floor per person for home, work and
services, houses about 85,000 people per km². Averaged with parks, lakesides, schools and lower districts, a
metropolitan density of 30,000–50,000 per km² is viable. The ground each needs, with the lakes left open:

| Population | 20,000 per km² | 30,000 per km² | 40,000 per km² | 50,000 per km² |
|---|---|---|---|---|
| 50 million | out to 32 km | 26 km | 23 km | 20 km |
| 100 million | out to 48 km | 38 km | 32 km | 29 km |

At 40,000 per km², 100 million fit on the summit's crown and upper slopes, leaving the lower rise and Korolev's
shore for farms, foraging and wild land. Density is capped by daylight between towers (the Sun takes a week to
cross the sky, so shadows linger, though the bright diffuse sky helps), by moving people, and by the city's own
heat: at 2 kW per person, 40,000 per km² releases about 80 W/m², a large addition to the sunlight the ground
absorbs, which would noticeably warm the summit's cool air.

**Character.** The metropolis uses the third dimension that low gravity and dense air open up. Small sky boats
serve much as cars do, and gliders fly throughout. Glider and parkour zones are designed in, so the rest of the
architecture does not invite either.

**Transit.** A high-capacity backbone of metro, rail and large sky ferries carries the metropolis, with sky boats
and gliders alongside, as cars run alongside metros on Earth. The port alone moves 60,000–100,000 travellers an
hour.

**Power.** Fusion supplies most of the Moon's local power, following the recovered industrial architecture's
sequence (fission to start, then deuterium–tritium, then deuterium–deuterium or advanced fusion). It is a
conditional technology here: no reactor has been designed or validated. Its steady output suits the 354-hour
night. Solar supplements it by day: level panels give about 33 W/m² averaged over the lunar cycle, so local
sunlight alone could carry about 10–20 million people at 2 kW each. The summit's wind is weak near the ground
(16 W/m²) and stronger aloft (about 160 W/m² 20 km up): screened slow rotors over the port's whole face would give
about 575 MW, 140% of the port's own use and about 0.3% of the metropolis's 200 GW.

**Food.** Regional agriculture and foraging feed the metropolis, drawing on the biosphere work: the canopy model's
crop stand fixes about 12.5 g of carbon per m² per day under the Moon's sky at the equator. The land each person needs has not
been computed.

**Design notes, proposed and not yet decided.** Architecture for a very wet, stormy climate: arcades, deep eaves,
drained roofs, canals and rain gardens, with drainage as a main job of the city's form. Districts set round the
crater lakes. The city flows in under the tower between its legs.

## Flyers

The author adopted the flyer findings of the [sky-ship](../../research/studies/sky_ships/README.md) and
[sky-fleet](../../research/studies/sky_fleet/README.md) studies on 2026-09-28
([decisions register](../../research/decisions.md#sky-ships-and-flyers)). What they mean for the summit:

- **People on wings.** A pedalled wing of about 11 m flies level on 63 W, within what an untrained adult holds for
  three hours. A foot-launched micro glider of 9 m and 4.5 m² glides 15 to 1, so one launched from ring 0, 3 km up,
  reaches anywhere in the metropolis. Canopies about 4.8 m across laid flat land a person at 4 m/s from any rim of
  the port; the summit tower study's canopies, 3.7–4.5 m across with a coefficient of 1.3, are about a tenth larger
  than Knacke's coefficients need. Thermals rise through a daytime mixed layer up to about 10 km deep on the
  equatorial ring's land, with cloud base near 7 km, and there are none at night.
- **Sky boats.** Four-seat electric sky boats of about 700 kg hover on 11 kW and fly 10 km on 0.74 kWh. In a
  planning case after Tokyo's 23 wards, with sky boats taking cars' tenth of 280 million trips a day, about 290,000
  are aloft at the busy hour, 91 over each km². NASA's urban air mobility spacing, in lanes stacked every 61 m from
  300 m to 1.5 km, holds sky boats for 7–8% of trips; a car-like share needs dense automated traffic management.
  Owned at Tokyo's rate for cars, their berths would cover three-fifths of the land; shared, 1.1 million serve the
  same trips from 65 km² of berths on roofs and façades.
- **Sky ferries.** Rigid hydrogen ferries of about 250 m carry 1,200 people at 25 m/s. Lines every two or three
  minutes carry 24,000–36,000 people an hour each way at 42 km/h with their stops, about half what the busiest metro
  lines of Hong Kong and Tokyo carry. About 890 are aloft at the busy hour.
- **Energy.** All the metropolis's flying uses about 1.0 GW, half a percent of its 200 GW.
- **Heights.** Wings, gliders and canopies fly below about 300 m above the ground, sky boats and air taxis to
  1.5 km and ferries to 2.5 km; regional ships use the docks 10–15 km above the summit, long-haul liners and
  freighters the flight band, and winged liners cruise near 55 km.
- **Storms.** Small flyers land for the afternoon storms, and the backbone carries on.

## The explorable map

An Unreal map of the metropolis and its landscape, with the tower at the centre, is being built in
`/media/projectspace/tmp` on a copy of the game project (github.com/dgoldman0/terluna-game), outside Git until
it reaches milestones:

- A few representative districts in full detail, for example the precinct under the tower, a dense core, a
  lakeside neighbourhood, terraces on a crater wall and an outer district down the rise. The rest of the
  metropolis is generated by the same rules as low-detail massing, labelled as such.
- Terrain from measured topography (LOLA with Kaguya stereo, SLDEM2015, about 59 m), with generated fine detail and
  the drainage estimate's lakes. Erosion is not modelled.
- The day first; the night comes later.
- Labelled placeholders: vegetation (the engineered plants are undecided), ships (no ship class is designed) and
  weather. The sky is the game's research sky without Earth.

## Open

- A structural analysis of the hand-sized form, taking up what the hand sizing leaves out (joints and nodes,
  fatigue, construction stages, gust dynamics, ice, the ground under the feet and the tie), and the crown's docking
  arms for the chosen long-haul class.
- The number and spacing of the disks, and the width of the wide lower rings (150 m is an example).
- The long-haul ship class: four sizes compared and drawn on the sized form (sky-fleet study, visualization/sky-fleet).
  The study recommends the 750 m liner, with winged liners on landing pads for fast trips. A wind study of the
  crown's wake comes before the berths are fixed.
- The sky boats' share of trips, their ownership and their traffic system.
- The land per person for food, and the power and heat plan for the chosen population.
- The metropolis through the 354-hour night: about 82 hours of usable twilight at each end and about 190 hours that
  need lamps, by the sea-appearance branch's definition of dusk.
- From the [infrastructure review](../../research/studies/infrastructure_review/README.md):
  - storm protection: the tower as a lightning conductor reaching to the height where flashes start, the crown's
    hydrogen berths, lightning protection for charges beyond Earth's highest level, the backbone's ferries in storms,
    and a storm warning service;
  - the frame and the crown sized with a design gust that varies with height: over the port's upper third the
    corrected winds would give 36–39 m/s, against the single 32.4 m/s;
  - wind devices in the upper frame, which give the port more than its own use in the corrected winds;
  - the high platforms' height or drift in the stronger winds near 70 km;
  - the link between the summit and orbit, and what launches release into the upper air.
