# Decisions

The author's decisions about the Open Moon (Terluna), with where each is recorded
and whether it stands. A later decision replaces an earlier one, and the
replaced entry names what replaced it. Modelling choices, reference values and
results appear here only where the author adopted them.

Sources:

- **July**: the recovered canon of the July 2026 knowledge bundle
  ([archive/july_2026_knowledge_bundle/](../archive/july_2026_knowledge_bundle/README.md);
  its file `01_RECOVERED_CANON_AND_CORRECTIONS.md` and `data/decision_log.csv`).
- **Sep 4**: the user requirements of the September feasibility report
  ([baselines/feasibility/report.md](baselines/feasibility/report.md), section 1).
- **Sep 9**: the September protection report
  ([protection/report.md](../protection/report.md), section 3).
- **Sep 19**: the 19 September 2026 industry discussion, known through a
  reconstruction written on 2026-09-26
  ([engineering/reference/industrial_architecture/](../engineering/reference/industrial_architecture/README.md),
  `RECOVERED_INDUSTRIAL_ARCHITECTURE.md`, section 8). The author confirmed its
  six traffic-safety rules on 2026-09-26.
- **Author**: decisions the author stated in working sessions, with the date.
- The conservation study keeps its own decisions, D1–D11, in
  [its README](studies/conservation/README.md#4-decisions); they are not repeated here.

## The world

| Decision | Source | Status |
|---|---|---|
| The world is Earth's own Moon, transformed in place, and the project is called Terluna. | July | Stands. The paper ensemble's working title is *Constructing and Sustaining an Open Moon*. |
| The Moon keeps its orbit, gravity, synchronous rotation and 29.53-day solar cycle. Life and culture adapt to the long day and night, storing resources across them. | July | Stands |
| Construction within 500 years of 2026, and operation for at least 10⁹ years. | Sep 4 | Stands |
| An open atmosphere with no pressure dome over the Moon. | Sep 4 | Stands. Local enclosed places, such as the undersea Tranquility community (conservation D6), are separate. |

## Light and time

| Decision | Source | Status |
|---|---|---|
| Practical dawn and dusk each last about five hours, inside the 29.53-day cycle. This replaced an earlier 30-hour figure. At the equator five hours is about 2.54° of solar elevation. | July | Stands; settled on 2026-09-26 against the ensemble plan's nominal 12.5-hour sunset, a design case. It still needs a photometric definition for the actual atmosphere. |

## Atmosphere

| Decision | Source | Status |
|---|---|---|
| Surface pressure 1.2 atm (121,590 Pa). | Sep 4 | Stands; the design pressure of the climate runs. The July bundle's 80 kPa was its own reference value. |
| An exobase near 250 K is preferred, with 260 K a candidate limit. | Sep 4 | Stands. The solved upper air behind the titania film gives 142–193 K before the sky's Lyman-alpha glow is counted ([middle atmosphere](../atmosphere/middle_atmosphere/README.md)). |
| Water, nitrogen and oxygen come from the Solar-System-wide resource operation. | Conservation D4 | Stands; replaced the September baseline's oxygen from lunar rock. |
| Atmospheric layers take descriptive names, never L1–L5, which stay reserved for Lagrange points. | July | Stands |
| "Gravity wave" means an atmospheric buoyancy wave; write "atmospheric gravity wave" on first use. | July | Stands |
| A band at 35–45 km was placed as the top of the main troposphere, provisionally, pending a climate model. | July | Settled on 2026-09-26: the band stays as a named flight band, and its physics comes from the models. The 1-D column puts the tropopause near 86 km (288 K surface); in the 5% climate the band holds air at −10 to −20 °C and 0.58–0.48 atm ([climate/gcm](../climate/gcm/README.md)). |

## Water and climate

| Decision | Source | Status |
|---|---|---|
| Standing water covers 28% of the surface. | [shared/scenarios/water.json](../shared/scenarios/water.json); conservation D10 | Stands; selected on 2026-09-25 within the author's 25–35% range. |
| The shield passes 5% less sunlight than the titania stack at every wavelength, and the Moon settles near 294 K. | [climate/gcm](../climate/gcm/README.md) | Stands; chosen on 2026-09-26. |

## Protection

| Decision | Source | Status |
|---|---|---|
| The protection architecture: a solar-filter complex with an industrial hub; magnetic protection; protection that grows around occupied destinations and Earth–Moon traffic corridors, secondary to protecting the lunar atmosphere and surface; Earth–Sun L1/L2 hubs for freight, power and industry. | Sep 9, carrying the original 2025 architecture | Stands |
| No superconducting planetary ring. | Sep 4 | Stands |
| The protection is designed across a range of total loss budgets, 1–100 kg/s, until the designs show which rate works best. Each budget is stated with its atmospheric cycle time, the time the loss and the resupply that matches it take to replace the whole atmosphere: about 100 billion years at 1 kg/s, 10 billion at 10 kg/s and 1 billion at 100 kg/s. | Author, 2026-09-26 | Stands |
| A surface radiation dose of at most 0.027 mSv/day. | Sep 4 | Stands |
| Protection hardware is a formation of replaceable units, never a single megastructure. Ageing and damage are handled by replacing units, supplied over the long term. | Author, 2026-09-26 | Stands |
| The protection design covers the Moon and its atmosphere. Destinations and corridors come later, and the older Earth–Moon–L1–L2 megastructure concept stays outside it. | Author, 2026-09-26 | Stands |

## Transport and safety

The author confirmed these rules on 2026-09-26 and regards them as set, to be reopened only for good reasons.

| Decision | Source | Status |
|---|---|---|
| Passive failures miss planets: nominal transfer paths and their failure dispersions avoid inhabited bodies unless capture has been verified. | Sep 19 | Stands |
| Earth never lies behind the lunar catcher, so a missed lunar capture cannot become an Earth-impact trajectory. | Sep 19 | Stands |
| Major braking happens remotely: high-energy packets approach the inhabited Moon only after capture and deceleration. | Sep 19 | Stands |
| Packet energy is capped. Mature dense traffic uses packets of roughly 10⁶–10⁸ kg in preference to 10¹²–10¹⁵ kg units, so each failure's consequences stay bounded; the exact cap awaits risk optimization. | Sep 19 | Stands |
| Planetary-scale surveillance, tracking, interception and tug capability, and assigned arrival corridors are part of the infrastructure. | Sep 19 | Stands |
| No civilization-threatening kinetic payload enters a planetary intercept corridor before verified capture. | Sep 19 | Stands |

## Life and people

| Decision | Source | Status |
|---|---|---|
| Aerial and gliding life is a major feature of the world. | July | Stands |
| High-altitude life is sparse: microbes or engineered radiation-hard films. | July | Stands |
| Long-endurance high-altitude platforms serve farside astronomy and nearside earthshine photometry. | July | Stands |

## The summit port and metropolis

The author set these on 2026-09-27. The [summit tower study](studies/summit_tower/README.md#the-central-port)
sizes the port, and the [metropolis brief](../habitation/summit_metropolis/README.md) collects the tower's form,
the metropolis and the first estimates behind them.

| Decision | Source | Status |
|---|---|---|
| The Moon's largest sky port stands on the far side's summit (5.4° N, 158.6° W). It is commercial: trade, travel, hospitality and short stays, with no homes and no industry. It holds about 1.5 million people on a regular basis. | Author, 2026-09-27 | Stands |
| Long-haul sky ships dock in the flight band, so the tower rises 24 km, to 35.6 km above sea level. Regional docks sit at 10–15 km, and every block of floors has berths. | Author, 2026-09-27 | Stands |
| The frame is a round diagrid that gathers into six splayed legs, as on the Eiffel Tower, through a deep transfer ring a kilometre or two up. The ground between the legs stays open. | Author, 2026-09-27 | Stands; the tower model has yet to size the legs and the ring |
| The floors are rings in the lower tower and disks round a central core from about 12 km up, as in the author's picture of the port, each a block of a few storeys with a park and open-air shops on its roof. They replace the full rings of 2026-09-27. The lower rings, where the open air is mildest, are wide (at 150 m, ring 0's park would be about 3 km²), with fewer storeys under the park and no terraces. Their decks cantilever from the frame, chosen for their look; cables from the frame above would be lighter. | Author, 2026-09-28 | Stands; the tower model has yet to size the rings, the disks, the core and the cantilevers |
| Open-air parks follow the air: for everyone up to 3 km above the summit, for residents and acclimatised visitors up to 10 km, and a few short-visit high gardens, winter parks among them, from 10 to 15 km. Above 15 km the blocks are sealed. Every interior is kept at a comfortable pressure and temperature. | Author, 2026-09-27 | Stands |
| The disks' open decks are all cold, about 1 °C at 12 km and colder above, and give the Moon its winter: the open winter parks sit just below freezing at about 13.5–15 km, and it is colder still above 15 km, where the blocks are sealed. The mild open-air parks are the ground and the lower rings. | Author, 2026-09-28 | Stands |
| A beautiful, sustainable mega-metropolis surrounds the port and centres on it, planned for about 100 million people at about 40,000 per km², within about 32 km of the tower. | Author, 2026-09-27 | Stands, as a rough planning scale |
| The metropolis uses the third dimension: small sky boats serve much as cars do, gliders fly throughout, and glider and parkour zones are designed in. | Author, 2026-09-27 | Stands |
| A high-capacity backbone of metro, rail and large sky ferries carries the metropolis, with sky boats and gliders alongside. | Author, 2026-09-27 | Stands |
| Fusion supplies most of the Moon's local power, supplemented by regional solar and wind. | Author, 2026-09-27 | Planning assumption; fusion is a conditional technology here ([engineering](../engineering/README.md)) |
| Regional agriculture and foraging feed the metropolis. | Author, 2026-09-27 | Planning assumption; the land per person has not been computed |

## Sky ships and flyers

The author adopted these findings on 2026-09-28 as the planning basis for the Open Moon's flyers. The
[sky-ship study](studies/sky_ships/README.md) and the [sky-fleet study](studies/sky_fleet/README.md) work them out
at first order, and the [metropolis brief](../habitation/summit_metropolis/README.md#flyers) collects what they mean
for the summit. No flyer has been designed. The long-haul class for the crown, whether winged liners join it, and
the sky boats' share of trips, ownership and traffic system stay open.

| Decision | Source | Status |
|---|---|---|
| The rule of similarity: under a sixth of Earth's gravity, a flyer six times larger in every length, flying at the same speed in air of the same density, carries the same stresses from its weight and from the same gusts. From 9 km up the Moon's air has the density of Earth's at a sixth of the height above 9 km, so the flight band (35–45 km) holds the air of Earth's 4.3–5.8 km. | Author, 2026-09-28 ([sky-ship study](studies/sky_ships/README.md)) | Stands |
| Rigid hydrogen ships are reasonable up to about 3 km in the flight band and at their best near 750–900 m; their structure would take all their lift near 7 km with 1930s materials and near 30 km with modern ones. Pressure hulls reach 280–350 m there, and wings carry their weight with the An-225's share at about 710 m of span. | Author, 2026-09-28 ([sky-ship study](studies/sky_ships/README.md)) | Stands; the hull model is fitted to four airships' weight statements |
| Wings and rotors are cheap under lunar gravity: hovering takes a fifteenth of Earth's power and a kilometre on a wing a sixth of the energy, while buoyant lift per cubic metre is Earth's. The small flyers are winged or rotor-borne; from ferries up, the flyers are buoyant ships, which hold still at their berths for free. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands |
| People fly on their own power. A pedalled wing of about 11 m span flies level on about 63 W, within what an untrained adult holds for three hours; a foot-launched micro glider of about 9 m and 4.5 m² glides 15 to 1 at 8 m/s; a round canopy about 4.8 m across laid flat lands a person at 4 m/s. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands |
| A four-seat sky boat is an electric flyer that lifts off on rotors and cruises on a wing: about 700 kg, 11 kW to hover and about 0.74 kWh for a 10 km trip, a twenty-fourth of the same trip on Earth. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands; sized with the Joby S4's technology |
| Large sky ferries are rigid hydrogen ships of about 250 m carrying about 1,200 people at 25 m/s; a line every two or three minutes carries 24,000–36,000 people an hour each way. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands |
| At a car-like tenth of the metropolis's trips, sky boats fly more densely than NASA's urban air mobility spacing allows, which holds 7–8%; owned at the rate cars are owned in Tokyo, their berths would cover three-fifths of the metropolis's land, while a shared fleet of about 1.1 million fits in 65 km². | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands as a finding; the share of trips, ownership and the traffic system are open |
| The flyers stack by height: wings, gliders and canopies below about 300 m above the ground, sky boats and air taxis to 1.5 km, sky ferries to 2.5 km; gliders soar to about 10 km by day; regional ships at the regional docks, 22–27 km above sea level; long-haul liners and freighters in the flight band; winged liners near 55 km; high platforms near 70 km. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands, as a first arrangement |
| Winged liners that take off and land on tilting rotors cross the Moon fast for little energy: a 300-seat liner reaches a mean destination in about 5 hours for about the energy per passenger-km of the largest airships. | Author, 2026-09-28 ([sky-fleet study](studies/sky_fleet/README.md)) | Stands as a finding; whether they serve the crown is open |
