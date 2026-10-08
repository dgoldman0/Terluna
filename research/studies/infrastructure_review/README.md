# The infrastructure work beside the other branches

This branch (`habitation/infrastructure`) worked out four things on 27–28 September 2026, from main as it stood on
27 September (70eed7c): the summit sky port and the metropolis round it, the port's fire safety, and the Open
Moon's sky ships and flyers. Four lines of work have moved on since:
- **Main** corrected its climate model and paused the climate work. It also added the waves and tides of the seas
  and two optical studies.
- **The atmospheric electricity branch** electrified the storms.
- **The sea-appearance branch** gave dawn and dusk a photometric meaning and found that the twilight wraps round
  the Moon.
- **The solar-shield branch** designed the shield's ring fleet, its power and industry, and the habitat array off
  the Moon.

This review reads each of them against the infrastructure work:
- what holds and what changes;
- what they add: storm protection, which the branch had not covered, the night, the coasts, and the
  infrastructure off the Moon;
- the questions they open.

The author asked for it on 8 October as the branch's closing work: "We'll cap this branch with that investigation."
The branch's merge into main waits for the joint integration of the open branches.

**How it was made.** The review reads each branch's committed documents and results at the commits listed under
[Sources](#sources). It adds three measurements, each stated where it is used:
- the branch's own wind exporters run on main's corrected climate run, committed here as two climate products;
- the four studies rerun unchanged on a temporary copy of this branch, with main's corrected inputs in place of
  the September ones ([How the numbers were made](#how-the-numbers-were-made)). The studies' committed results keep
  the September inputs until the integration reruns them;
- light arithmetic on the other branches' products: strike counts, twilight hours, fields, hydrogen inventories.

Every lunar number taken from another branch is that branch's model result or estimate, with the caveats it gives.

## In brief

- **The September climate inputs came from a faulty run.** The design run the studies read, `A28_dim5`, had
  Earth's 23.7° tilt and a Sun clock that gave the summit's model cells 7–8% too much sunlight. The corrected run,
  `A28_dim5_moon`, keeps the air's temperature within about half a degree below 50 km (a degree or so warmer at
  60–70 km) and its density within half a percent. Its westerly winds aloft are about twice as strong: in the
  flight band the median 3-day wind is 9–13 m/s (4–5 before) and the highest 20.5 m/s (13.4).
- **With main's corrected inputs the studies' numbers move:**
  - the tower's design gust, from 27.3 to 32.4 m/s;
  - the port's chosen form, from 90 to 98 Mt of steel;
  - screened slow rotors over the port's whole face, from 74% to 140% of the port's own use;
  - the moored ships' design wind, from 34.6 to 38.6 m/s, and their service wind from 13.6 to 24.4 m/s;
  - the power the high platforms near 70 km need to hold station, about seven times larger.

  The panels (32.8 W/m² against 33.0), the lakes round the summit and the Korolev lake stay as they were. The
  summit's model cells rain 9–17 mm a day (11–22 before).
- **Lunar lightning is rare and very large.** The electrified storms flash 0.022 times per km² a year, a hundredth
  of Earth's average. A fifth of the flashes strike the ground, lowering a median 149 C and up to 1,025 C. Earth's
  highest lightning-protection level is sized for 300 C.
- **The port's crown stands where the flashes start.** Flashes start at a median 35.8 km above sea level, and the
  crown stands at 35.0–35.6 km. A grounded steel frame 24 km tall, reaching into a storm's charging layer, would
  very likely start lightning whenever a charged storm comes near. This is new ground for every branch, and the
  ships berthed at the crown hold about 2,100 t of hydrogen.
- **The flight band is above most storms and inside the strongest.** The storms that flash have updrafts of
  12–36 m/s and tops at 66–98 km; the sky-ship study's design gust was 14.3 m/s. Ships keep out of flashing
  storms, which give about an hour's warning between a storm's core forming and its first flash.
- **The night has a long twilight.** At the summit, practical dusk (while ordinary outdoor activity needs no lamps)
  lasts about 82 hours after sunset and again before sunrise. The dark part of the 354-hour night is about 190
  hours, and below a tenth of a lux about 116. The crown keeps the Sun about a day longer on each side of the night.
- **The coasts have monthly tides and gentle seas.** Typical monthly tidal ranges are 3.7 m over the nearside sea
  and 5.5–6.0 m in its southern maria. For the same wind, waves grow higher and slower than Earth's, but the winds
  are light, and wave power runs about two orders below Earth's exposed coasts.
- **Off the Moon, the shield branch plans a fleet of about 28 million tiles and a habitat array, and puts
  energy-intensive industry and computing in orbit.**
  - Its power links run from orbit to orbit, and transport between the surface and orbit is still open.
  - Its outflow rule S7 covers the protection systems.
  - Its one large work on the surface is the regional magnets: four circuits of 500–1,000 km radius near the poles,
    with 9–15 GW of refrigeration running through the night and fields that set a corridor along their routes.
- **The fleet's hydrogen is small beside the air's own.** The summit's fleet holds tens of thousands of tonnes of
  hydrogen. The air holds about 105 million tonnes at the 0.53 ppm the escape model takes.

## 1. Main's corrected climate

### What changed in the inputs

On 29 September main found two faults in the GCM's sunlight (climate/gcm README):
- **The tilt.** ExoPlaSim's default orbit gave the runs Earth's 23.7° tilt in place of the Moon's 1.54°.
- **The Sun clock.** PlaSim's Sun clock swept 392.5° a model day and then skipped back 32.5°. The far side around
  164° E took up to 17% more sunlight than it should, and longitudes 148–180° E saw noon twice a lunar day.

Along the row at 2.8° N, the model cells round the summit received 399–405 W/m², against 373 W/m² once corrected
(`climate/results/gcm/sun_check_A28_dim5.json` and `_moon.json` on main). The corrected design run,
`A28_dim5_moon` (years 20–29), settles at 294.9 K against 294.3 K. The author kept the 5% dimmer shield and
accepted the corrected climate's drier land and smaller lakes (decisions register, 2026-09-29).

What the studies read, against main:

| Input | This branch | Main | What changes |
|---|---|---|---|
| GCM winds and air at the summit and over the Moon | `site_winds_A28_dim5_summit.json`, `global_winds_A28_dim5.json` | The corrected run's winds have no product on main. This review exports both from `A28_dim5_moon`, years 20–29, with the branch's own scripts | Temperature within about half a degree below 50 km (a degree or so warmer at 60–70 km) and density within half a percent; the winds aloft about twice as strong (table below) |
| The equatorial CM1 ring | `ring_ring.json`, forced by the faulty run | The same values for everything the studies read (it gained fog and latitude fields). Its successor, `ring_ring_equator.json`, is forced by the corrected run and also takes the GCM's rising air and land wetness column by column | Storms more numerous and shallower at the median, stronger winds aloft (below) |
| Lakes | `drainage.json` before the correction | Lakes cover 10.3% of the Moon (11.6% before) | None round the summit (below the table) |
| Sunlight at the ground | `surface_light.json` | Rerun at 294.9 K | Level panels 32.8 W/m² (33.0) |
| Rain at the summit | 20–24 mm a day in the rainiest cells nearby | 9–17 mm a day in the four cells round the summit (11–22 before), 19–20 in the rainiest nearby (22–24); equatorial land 4.8 (4.6) | The summit stays the rainiest ground, with about three-quarters of the rain |

Round the summit the lakes are the same in both grids:
- the Korolev lake keeps 140,652 km², its surface at 7,849 m and depths to 6,670 m;
- lakes cover 24%, 41% and 57% of the ground within 25, 100 and 300 km of the summit.

The winds over the whole Moon, as 3-day means over ten model years (`global_winds_A28_dim5.json` and
`global_winds_A28_dim5_moon.json`, here):

| Height above sea level | Median, September → corrected | Highest |
|---|---|---|
| 0 km | 3.8 → 3.5 m/s | 11.5 → 9.8 m/s |
| 10 km | 2.6 → 2.6 | 9.1 → 8.4 |
| 20 km | 2.9 → 3.9 | 9.6 → 10.3 |
| 30 km | 3.5 → 6.8 | 10.8 → 14.0 |
| 35 km | 4.1 → 8.7 | 10.2 → 15.3 |
| 40 km | 4.7 → 10.8 | 11.8 → 18.4 |
| 45 km | 5.3 → 12.8 | 13.2 → 20.5 |
| 55 km | 6.5 → 16.9 | 16.2 → 25.3 |
| 70 km | 7.8 → 20.9 | 15.9 → 29.3 |

At the summit, at the crown's height (35.4 km above sea level), the median goes from 4.1 to 7.8 m/s and the highest
from 8.3 to 11.2 m/s, blowing from the west with a steadiness of 0.99 (0.95 before). Near the summit's ground the
winds stay as they were or grow lighter.

The equatorial ring forced by the corrected run, against the ring the studies read:

| | The studies' ring | `ring_equator` |
|---|---|---|
| Storm cloud tops: median, 90th percentile, highest | 28, 68, 82 km | 22, 38, 84 km |
| Storms raining at once along the ring | 2.9 | 6.6 |
| Rain on land; share of the time a land place is under rain | 2.1 mm a day; 3.3% | 3.9 mm a day; 5.7% |
| Highest updraft; at 40 km | 12.4; 14.3 m/s | 15.6; 15.5 m/s |
| Highest 3-hourly wind at 40 km; its 99th percentile | 20.6; 8.9 m/s | 23.0; 15.0 m/s |
| Storm cloud in the flight band | 0.38% of land columns | 0.69% |
| Daytime mixed layer, deepest; cloud base | 11.7 km; 9.9 km | 9.9 km; 7.3 km |

The two rings differ in setup as well as forcing, and both are flat, two-dimensional and at sea level. The storm
heights the four studies and the metropolis brief quote (28, 68 and 82 km) are the first ring's.

### What it changes in the studies

The four studies rerun, unchanged, with main's corrected inputs (the GCM products exported here, `ring_equator`,
main's drainage and sunlight):

| Result | September inputs | Corrected inputs |
|---|---|---|
| Highest modelled wind over a 10 km tower | 13.6 m/s | 16.1 m/s |
| Design gust (× 1.4 × 1.2 × 1.2); service gust | 27.3; 19.0 m/s | 32.4; 22.5 m/s |
| The port as a square lattice: frame steel; base | 32.8 Mt; 8.2 km | 46.0 Mt; 9.3 km |
| The port's chosen form: steel; sway in the service gust | 90.0 Mt; 30 m | 98.0 Mt; 37 m |
| The reference 10 km tower: frame steel | 2.9 Mt | 3.5 Mt |
| Wind energy 20 km and 30 km above the summit | 39 and 64 W/m² | 158 and 469 W/m² |
| Screened slow rotors over the port's whole face | 295 MW, 74% of its use | 575 MW, 140% |
| The square lattice's steel with those rotors stopped for storms | 75.6 Mt | 111.3 Mt |
| Screened slow rotors over half the face | 147 MW, 37% | 287 MW, 70% |
| Level panels | 33.0 W/m² | 32.8 W/m² |
| A stores compartment burning out in the sealed zone at the median wind | 484 MW for 2.1 h | 818 MW for 1.2 h |
| Moored ships in the band: design wind; service wind | 34.6; 13.6 m/s | 38.6; 24.4 m/s |
| A 750 m liner: head-on mooring load; useful share | 0.34 MN; 86.1% | 0.43 MN; 85.8% |
| Design updraft in the band (the 2-D ring) | 14.3 m/s | 15.5 m/s |
| Largest pressure hull in the band at 25 m/s | 350 m | 330 m |
| A high platform near 70 km: 99th-percentile and highest wind | 12.7 and 15.9 m/s | 24.9 and 29.3 m/s |
| Power to hold a 100 m and a 200 m platform at the 99th percentile | 15 and 62 kW | 115 and 460 kW |
| The 200 m platform's envelope, as a share of today's laminates' strength | 0.39 | 1.2 |
| Storm hours over equatorial land in a lunar day | 157 | 217 |

Notes on these numbers:
- **Gravity sizes most of the chosen form.** Its steel grows by 9% where the square lattice's grows by 40%. The port's
  programme figures (passengers an hour, berths) move by up to a fifth only because the wider square lattice carries
  36.0 km² of floor in place of 34.9; that is the sizing method's bookkeeping, and the author's 1.5 million people
  stand.
- **The design gust is one value for the whole height.** It comes from the highest wind over the first 10 km above
  the summit. Over the port's upper third the corrected ring's highest 3-hourly winds reach 17.7–19.5 m/s
  (14.0–15.3 before), which would give 36–39 m/s with the same factors. A design gust that varies with height is the
  next refinement for the frame.
- **The tower as a power plant.** On 27 September the author asked whether wind devices in the frame could at least
  power the tower. In the corrected winds they can:
  - screened slow rotors over the whole open face give 140% of the port's use, for 65 Mt more steel than the open
    frame (in the square lattice's sizing);
  - over half the face they give 70%, for 28 Mt more.

  The air at the crown's height carries 195 W/m² in 3-day means, against 34 before. Low on the tower, where most of
  the face is, the wind stays weak, so the devices' yield comes mostly from the upper frame.
- **A steady westerly in the band.** At the band's median of about 11 m/s, a liner cruising at 30 m/s gains about a
  third of its ground speed eastbound and loses it westbound. An east–west round trip then takes about 15% longer
  than in still air, and the mean trip of 25 hours grows by several hours on such routes. Faster cruise, flying
  lower (the median is 5 m/s at 25 km), or routing that rides the wind eastbound answer it.
- **High platforms.** At 70 km the median wind is now 21 m/s. Platforms there hold station on several hundred
  kilowatts, or drift round the Moon with the westerly, or move to a height where the wind slackens. The GCM's wind
  peaks near 75–80 km and weakens above.
- **Gliders.** By the corrected ring the day's thermals reach about 10 km, with cloud base near 7 km, and the land
  sees more storm hours.

### What holds

- **The air by height.** The freezing level at 25.4 km above sea level (13.6 km above the summit) and the oxygen
  each height gives hold, and with them the parks by height, the winter on the disks and the fire study's air.
- **The water.** The summit's site, its high lakes and the pumped-storage head of 3.74 km to the Korolev lake hold.
- **The panels.** The clear-sky yield holds within 1%. The optical studies add that the solved sky gives several
  times the two-stream model's light with the Sun within a few degrees of the horizon, which lifts the level panel's
  mean by a few percent and east and west walls by more. Cloud cover of 17–25% brings the all-sky mean to about
  29 W/m². The summit's own cloud and its thinner air above are still to compute.
- **The sky ships.** Rigid ships of up to about 3 km in the band, and the rule of similarity, hold: the hull sizes
  barely move with the corrected winds.
- **The fire study.** Smoke, sprinklers, egress and water hold; only the wind-fed fire case grows, in the upper
  zones.

**Comfort.** Main paused the climate work on 30 September because CM1 and the GCM disagree about the air near the
ground over land (decisions register, Research practice). On the rings' land CM1 gives 42 comfortable hours a lunar
day against the GCM's 387 (18–26 °C with a dewpoint at most 15 °C), and both put the most comfortable land poleward
of 30°. The summit's free air averages about 12 °C, below the band's lower edge, and the summit's own comfort is
still to read. A dry highland box at 44.7° S gives 183 comfortable hours a lunar day, rising with height to 262 at
5–6 km. The metropolis's climate carries the range between the models until an independent GCM or new evidence
decides it.

## 2. Storm protection

This branch treated storms as weather to wait out:
- small flyers land for the afternoon storms;
- high parks close for the storms and their lightning;
- the tower "will be struck by lightning and needs the protection Earth's tall towers have";
- long-haul ships dock in the flight band "so they never descend through the storms".

The atmospheric electricity branch has since electrified the storms. Its main run, `box_0e_elec_corrected`, is a
box 385 km square at 0° N, 0° E: all land, flat at sea level, on 6-km columns, through two lunar days. It gives the
storms an infrastructure meets, and its heights are heights above sea level.

### The storms and their lightning

| | Main run (10th–90th percentile where given) |
|---|---|
| Flashes | 530 in two lunar days, all in the afternoon storms and none at night: 0.022 per km² a year, about a hundredth of Earth's average |
| Ground strikes | 106, a fifth of the flashes, all negative: 0.0044 per km² a year |
| A ground strike | Lowers a median 149 C (63–434 C, largest 1,025 C) and releases 66 GJ (23–162 GJ, largest 325 GJ) |
| A flash in cloud | Moves 108 C (48–218 C, largest 703 C) and releases 86 GJ (42–185 GJ, largest 758 GJ). Its channels span a median 29.8–41.8 km in height over 144 km², about 12 km across; the largest spreads over 1,156 km² |
| Where flashes start | A median 35.8 km above sea level (33.8–39.8 km), in a field of 152–188 kV/m, where the storm's potential averages 250–460 MV (up to 1.2 GV) |
| Charge | The main negative charge sits at 40–42 km (−13 to −14 °C), carried mostly by snow. Positive charge lies beneath it at about 32 km (−4 °C) and above it at 58–64 km. The −7 °C level is at 33.8 km |
| Storms that flash | About a tenth of storm cores, the strongest: peak updrafts of 12–36 m/s (medians 19–21), tops at 66–98 km, and 3–49 million tonnes of graupel and hail between 31 and 53 km. The other cores peak at a median 6 m/s with tops at 54–56 km |
| Warning and duration | The first flash comes 0.7–1.5 h after a storm's core forms. A storm flashes for a median half hour, at most 0.8 times a minute, and at most one or two storms flash at once in the box |
| Rain | 17–87 mm/h at the peak beneath flashing cores, about 0.8 h after their busiest flashing |
| Hail and graupel aloft | Hail 12–13 mm across, falling at 4.8–7.7 m/s; graupel 3.7–4.4 mm; both in the −5 to −30 °C layer |
| Thunder | At a ground strike 93–97 dBA, a 53 Hz rumble lasting about 27 s, with 103–107 dB at 31.5–63 Hz (about 4.5 Pa). Beneath a flash in cloud 60–68 dBA. Above a daytime background to roughly 40–50 km |
| The ground's field | Within six hours of a flash, a median 4.9 kV/m at the strongest point and 61.5 kV/m at most (with point discharge off) |

Under Takahashi's charging law the storms invert their charge, and their largest flashes grow: a ground strike to
1,489 C and 2.4 TJ, a flash in cloud to 1,964 GJ.

What the evidence is:
- one cloud model, on 6-km columns, at one equatorial lowland site;
- calculated lunar numbers, with Earth's benchmark storm checking the code and its schemes;
- a 2-km box has yet to flash, and coarser grids make fewer, larger flashes;
- the electrified runs cover lowland only. The summit, the rainiest ground in the GCM, may see more storms than the
  lowland box, or fewer by CM1's drier count.

### Lightning is rare and very large

Earth's lightning-protection standard (IEC 62305-1) sizes its highest protection level for flashes that lower up to
300 C. Earth's ordinary ground strikes lower 5–30 C. The lunar ground strikes' median is half the standard's
ceiling, their top tenth exceeds 434 C, and the largest lowers 1,025 C (1,489 C under Takahashi's law).

The charge a strike lowers sets how deep its arc melts into metal where it attaches. So roofs, hulls, cables,
joints and earthing are sized for charges beyond Earth's highest level, against events any one place meets seldom.
At the box's rate the metropolis's 3,200 km² takes about 14 ground strikes a year, with about 70 flashes overhead.
The lightning scheme gives each flash's charge and energy but not its peak current, so the currents and the heating
they bring are still to estimate.

### The tower's crown stands where the flashes start

The crown's terminal stands 35.0 km above sea level and the frame's top 35.6 km. In the electrified box:
- the flashes start at a median 35.8 km;
- the −7 °C level is at 33.8 km;
- the main negative charge sits at 40–42 km.

The summit's free air keeps about the same heights for its temperatures. The freezing level is 25.4 km above sea
level, and −9 °C stands at the crown in the corrected site product. So under a passing strong storm the tower's top
quarter would stand in the storm's charging layer, with the main negative charge 4–6 km above the crown. That is a
first reading, since no run has stormed over ground near 11.6 km.

A grounded conductor 24 km tall holds its top near the ground's potential while the air round it sits at the
storm's. Two Earth findings bear on what follows:
- tall structures start an upward leader when the potential of the air at their tops reaches about 1.6 MV (Rizk
  1990);
- structures taller than a few hundred metres start most of the flashes that strike them (Rakov and Uman 2003).

A storm that brought the air at the crown even to a few megavolts, a hundredth of the potential its own flashes
start at, would have the tower start flashes, probably well before the storm's field reached breakdown. The tower
would very likely start lightning whenever a charged storm came near. It is a lightning conductor reaching into the
charge.

Tall conductors are new ground for the atmospheric electricity branch. The questions they open:
- **How often the tower starts flashes, and how much charge it drains.** Early, smaller flashes and continuing
  currents through the frame may take the place of the storm's own large ones nearby.
- **What the frame must carry.** The current, charge and heating that the frame, its joints, the legs and arches, the
  six feet, the buried tie and the ground beneath must take.
- **Corona.** The upper frame gives off corona under storms; in the point-discharge test, vegetation alone gave off a
  median 3.4 A over the 148,000 km² box. The corona brings radio noise and changes the field round the tower.
- **The crown's ships.** They are berthed at the top of the conductor, with about 2,100 t of hydrogen.
- **The people.** The sealed blocks lie inside the steel frame's cage, and the open decks below 15 km close for
  storms already.

The branch's charge fields can answer the first of these. An electrostatic calculation that sets a grounded 24 km
conductor into the main run's charge structure gives the field at its top through a storm's life. From that comes
when the tower would launch leaders and how much charge it could carry.

### Ships in the flight band

The flight band sits above the median storm and inside the strongest tenth:
- Storms that flash reach 66–98 km, with peak updrafts of 12–36 m/s. The sky-ship study took 14.3 m/s, the strongest
  updraft at 40 km in the two-dimensional ring, as its design gust.
- Their graupel and their 12–13 mm hail fill the band, their flashes start in it, and the field there reaches
  150–240 kV/m before a flash.
- A conducting ship 750 m long spans tens to a hundred megavolts in such a field, depending on how it lies in it.
  On Earth, most lightning that strikes aircraft is started by the aircraft (Mazur 1989).

The branch leaves open how craft keep clear of storms and whether they would start flashes. It also leaves open
what dose the runaway electrons of the storms' strongest fields give in the band.

What follows for the flyers:
- **Ships keep out of flashing storms.** The storms give warning:
  - a storm flashes 0.7–1.5 h after its core forms, time for a liner at 30 m/s to cover 76–160 km;
  - at most one or two storms flash at once over 148,000 km².

  Routing round storm cores, as Earth's aircraft do at cruise, keeps the sky-ship study's sizing, which holds for
  ships that meet the two-dimensional ring's gusts. A ship built to ride out a strong storm's core needs its gust
  structure sized for up to 36 m/s, a rerun of the sky-ship model that takes seconds.
- **The crown's berths.** Ships berthed at the crown take longer to leave than a cruising ship takes to turn away,
  and they sit at the top of the tower's conductor. The answers, alone or together:
  - procedures, with ships casting off when a core forms within range;
  - air terminals standing above the berths;
  - hulls whose conducting outer cover keeps a flash outside the gas cells;
  - helium in the long-haul liners.

  The sky-fleet study already left hydrogen safety at scale and a helium supply open; the lightning raises their
  priority.
- **Charge lingers.**
  - Near the ground the air lets a charge relax in about an hour (the aerosol study), where near Earth's ground it
    takes minutes.
  - Cloud holds separated charge for about a day.

  Charge on hulls, mooring lines, refuelling hoses and people therefore stays longer than on Earth, and bonding and
  earthing during hydrogen handling at berths matter more.
- **The backbone in storms.** The metropolis brief has the backbone carry on through storms while small flyers land.
  Metro and rail carry on. The sky ferries are hydrogen ships of 250 m flying 1.5–2.5 km above the ground, under
  storms whose ground strikes pass through that air, so the ferries need the same rules as the liners: hold or
  reroute round flashing cores, or be built to take a strike.
- **Small flyers land for storms,** as the brief has them. The hour between a core forming and its first flash gives
  them time, and gliders soaring the afternoon thermals keep clear of growing cores.
- **High platforms and winged liners.** Flashing storms' tops at 66–98 km reach the platforms' 70 km and the winged
  liners' 55 km, and a tenth of the flashes in cloud reach above 51 km. Jets above storm tops are a question for the
  electricity work's stage 3.

### The metropolis on the ground

- **Strikes.** The metropolis takes about 14 ground strikes a year at the box's rate (0.003–0.006 per km² a year by
  the aerosol study's count), each of a charge Earth's highest protection level does not cover.
  - On Earth, structures start upward flashes of their own from heights of about 100 m (Rakov and Uman 2003).
    Under storms the ground's field reached tens of kilovolts per metre in the main run, so the metropolis's tallest
    towers would take most of its strikes.
  - Air terminals, down conductors, earthing and surge protection are sized for the charge, and the grid for rare
    large surges.
- **Thunder.** A ground strike within a few kilometres is a loud, deep rumble: 93–97 dBA, about 4.5 Pa at 63 Hz. Every
  flash within 40–50 km is heard over the day's background, and light glazing and panels feel the low frequencies.
- **People outdoors.** Flashing storms give about an hour's warning after their cores form, and the high parks and
  open decks close for storms, as the brief has them. A storm warning service reading the storm cores is part of the
  city's infrastructure.
- **Nitrogen.** Lightning fixes 0.02–4.6 kg of nitrogen per km² a year, 0.2–47% of what Earth's lightning fixes.
  That is a soil and ecology question more than an air-quality one.

### Rain, hail and wind

- **Rain.** Beneath flashing cores the rain peaks at 17–87 mm/h, heavy to violent by Earth's usual classes (violent
  above 50 mm/h). The architecture the metropolis brief proposes for a wet, stormy climate (arcades, deep eaves,
  drained roofs, canals and rain gardens) answers it. The summit's total stays open between the GCM's 9–17 mm a day
  and CM1's drier count.
- **Hail and ice aloft.** Hail and graupel stay above the freezing level, 25.4 km above sea level. They reach the
  tower's upper decks and the ships in the band, and melt before they reach the summit's ground.
  - A ship at 30 m/s meets 13 mm hail at about that speed: covers that take it, or routing round cores.
  - Supercooled cloud water (0.03–0.17 g/m³ in the storms) ices ships and the frame, beside the rime that night
    cloud leaves.
- **Wind.** The storms' surface gusts in three dimensions have yet to be read for structures. The tower's design gust
  keeps its three factors until a run that resolves the summit's terrain gives its own winds.

### Space weather

The four regional magnets give the Moon a magnetosphere that stands off the solar wind near 10 lunar radii. On
Earth, solar storms drive currents in the upper air, and their fields induce currents in long conductors on the
ground: power lines, pipelines, rails. Whether the Open Moon's grids meet the same is unmodelled. The plasma above
149 km, the magnets' field lines and the solar wind's potentials are open in the electricity branch's stage 3 and in
the shield branch's plasma items.

## 3. The night

The branch's brief says the far side's 354-hour night has no earthshine, and it left the night itself open ("the
metropolis through the 354-hour night"). The sea-appearance branch and main's optical studies now give the night
its light (clear sky, a sea-level column, no refraction):
- **The evening is long.** The Sun sinks 0.51° an hour and takes about an hour to set. Practical dusk, while ordinary
  outdoor activity needs no lamps, ends at 2.98 lux, the light at which Earth's civil twilight ends in the same
  solver. The Moon reaches it with the Sun 41.4° down, 82 hours after sunset at the equator. Day vision holds for the
  first 46 hours. The author adopted this definition on 4 October; it replaces the July canon's five hours.
- **The summit gets dark.** At latitude φ the Sun sinks at most 90° − φ below the horizon, so most of the Moon keeps
  a twilit night. Near the far side's equator, where the summit stands (5.4° N), the night grows dark round midnight.
  By the twilight tables:
  - about 190 hours lie beyond practical dusk;
  - about 116 hours lie below a tenth of a lux;
  - night vision begins about five days after sunset;
  - the darkest light is about 0.001 lux.

  The far side has no Earth in its sky. Starlight and airglow are a 0.001-lux placeholder.
- **Height brings the Sun back.** A point 40 km up keeps the Sun's centre in view for about a day after sunset at the
  ground, its light turning deep red. By the same geometry the summit's ground keeps the Sun up to about 13 hours
  longer at each end of the night (less where the land round it rises), and the crown about 22–23 hours. The
  tower's upper floors and the panels on its upper frame have a shorter night than the city below.
- **What it means for the city.** Lamps are needed for about 190 hours of each lunar month. Outdoor life, flying and
  the parks have about 82 hours of usable twilight at each end of the night, and 46 hours of day vision. The pumped
  storage still carries the full night for power, since twilight gives panels next to nothing. The light calendar
  (sea-appearance branch) computes the light at any point and date from 2000 to 2500, with named light phases, which
  answers part of the habitation lane's open item on civil time.
- **What the night sky holds.** The shield's tiles are lit round almost the whole of their orbits, and the shield
  branch notes that night-side tiles reflect light back past the Moon, to be checked. How bright the lit ring stack
  makes the night sky is still to compute, and it matters for the far side's dark nights and for the astronomy
  platforms.

## 4. Seas, coasts and lakes

The branch built inland, on the Moon's highest ground, so coasts, harbours and sea transport are new ground for it.
Main's wave and tide work (climate/waves and geography's tide product, forced by `A28_dim5_moon` through one measured lunar
cycle) gives what coastal infrastructure would be designed against:
- **Tides come once a month.** The Earth's bulge is built into the atlas's sea level. What remains is a monthly swell
  and rock from the orbit's eccentricity and the librations, with 7 cm from the Sun. Typical monthly ranges:
  - 3.7 m over the nearside sea (6.0 m at the 95th percentile of its area, 9.8 m the largest in 19 years);
  - 5.5–6.0 m in Fecunditatis, Nubium and Humorum (8.7–8.8 m in 19 years);
  - 2.1 m in the South Pole–Aitken sea;
  - 0.23–0.55 m in Smythii and Smythii–Marginis.

  High water comes about two weeks after low water. On beaches of 1:50 to 1:100 a 7.6 m range moves the waterline
  380–760 m.
- **Waves grow higher and slower than Earth's for the same wind,** 4.2–4.8 times the height and about five times the
  period for the same wind, fetch and duration. But the winds over the seas are light: 2.6 m/s at the median and
  5.8 m/s at the 99th percentile. Over the nearside sea:
  - Hs averages 1.39 m, is at least 1 m 68% of the time, and reached 4.9 m at most in the measured cycle;
  - waves are 34–130 m long and travel at 3–6 m/s;
  - swell carries a third to two-thirds of the energy.
- **Shores.** Wave power toward the coast averages 87 W/m at the median coastal node and 731 W/m at the most exposed
  (western Procellarum), with hourly storm peaks of 2.6–3.6 kW/m. Wave power runs two orders below Earth's exposed
  coasts and gives little energy.
  - At the Smythii headland, run-up reaches 2.07 m on rock at the 98th percentile and about 1 m on beaches.
  - Under Froude scaling a wave of a given height loads a structure about a sixth as hard as on Earth while running
    up as high. So water levels (tide, run-up, and surge of up to 0.18 m from air pressure) govern coastal works more
    than wave forces.
- **Loose coastal material.** The same waves move grains about six times larger than on Earth, and fine grains
  settle about six times more slowly. Shores therefore rework toward gentle profiles and coastal water runs turbid.
  Sediment transport is not computed.
- **Open.** These remain:
  - a storm and extreme climate: one measured cycle gives no return periods and no design wave;
  - the far-side and limb seas beyond Smythii–Marginis (the South Pole–Aitken sea holds 22% of the water);
  - tide and run-up together;
  - the nearside shores, where storm waves run about twice the Smythii forcing and tides 4–10 times larger.

**The lakes.** The wave and tide work covers the seas at sea level. The rain-fed highland lakes, the Korolev lake
among them, still have their waves, tides and seiches to compute. The summit tower study's pumped storage moves about
105 million m³ through each lunar cycle to carry the port through the night (0.53 km³ for a city of a million), from
a lake that holds 400,000 km³. The lake's own behaviour is new ground: its seiches, wind waves on a 141,000 km² basin
7.8 km above sea level, and its shoreline works. Smythii–Marginis, about twice the Korolev lake's area, reached Hs
3.4 m in its strongest weather and is the nearest computed analogue.

## 5. Off the Moon

The infrastructure branch stayed on the Moon. The solar-shield branch (head eeaecbc) designs what lies off it, and
it reaches the surface in two places: the regional magnets and the upper air's loss budgets.

### What the shield branch builds

- **The shield.** A ring fleet leads since 7 October:
  - about 1,860 nested rings, each one orbit of shingled 10 km tiles, at 19,600–21,900 km, tilted up to about ±26°
    from the Moon's orbit plane;
  - about 28 million tiles and 46 Gt, kept in place by photon forces with no propellant;
  - films replaced at 2.2–4.5 Gt a year by a film plant drawing 0.15–2.9 TW near 2,000 K, with 4,000–8,000 tile
    exchanges a day.

  Failed cells are covered within a day to a week from spares kept in the formation. The author keeps faster covering
  through solar maxima in view.
- **The habitat array.** Toroidal wheels round an optical disk, a research direction since 4 October:
  - a 5 km wheel turns at 0.42 rpm for 1 g, a 50 km wheel at 0.13 rpm;
  - habitats belong on natural orbits, apart from the heavy shield;
  - their structure, shielding, life support, size, population and orbits are unbudgeted.
- **Power.** The aperture intercepts about 240 PW, and 0.24–0.37 million km² of collectors on the ring fleet would
  give 100 TW. Every link the branch sizes runs from orbit to orbit, for example 100 MW at 5.8 GHz over 100 km at 51%
  from bus to bus. "Net-positive electricity for habitat" is an acceptance gate still open for both candidates, and
  the branch's receivers are all in orbit.
- **Industry and computing.** The heat screen finds that energy-intensive industry and computing belong in orbit:
  - each terawatt of computing needs 830 km² of radiator;
  - heat used on the Moon adds 0.026 W/m² per terawatt to its budget, against 0.19% of that from orbit.

  The planned branch `research/array-industry` decides which industries go to orbit and which to the Moon, under S7
  and the six traffic-safety rules.

### What it means for the Moon's infrastructure

- **Power plan.** The metropolis's plan keeps fusion with solar and wind, since the shield branch's power stays in
  orbit so far. Power from the shield to the ground is an open option:
  - it needs ground receivers, beam limits and a path through the tall air and its storms;
  - it adds its heat to the Moon's budget, so 100 TW used on the Moon would add 2.6 W/m², about the forcing that
    human activity puts on Earth's climate today.

  The metropolis's 200 GW adds 0.005 W/m² to the Moon's mean. Its local heat, about 80 W/m² over the city, is the
  question the brief already names.
- **Industry placement.** With energy-intensive industry and computing in orbit, the lunar cities' own demand per
  person falls and the freight link between the surface and orbit becomes central.
- **Transport between the surface and orbit.** What exists is the July bundle's sketch in main's status record: low
  lunar orbit ends with the atmosphere; L1/L2 logistics nodes, high orbits, tethers above the dense air, evacuated
  mass drivers, high-altitude airports and rocket corridors. The infrastructure branch's flyers
  already reach high: winged liners near 55 km, platforms near 70 km, and a port crown at 35.6 km, the highest built
  point.
  - Whether a high port, a tether foot or a launch site belongs on the summit's tower or beside it is a new
    question.
  - So is what launches through the tall air release into the upper air. S7 traces every outflow of the protection
    systems; launches from the surface have no rule yet, though main's July sketch already names rocket corridors.
- **The magnets on the surface.** The regional circuits are a screening design, not rerun since 7 October:
  - four loops of 500 or 1,000 km radius, each closing on its own cap near a pole;
  - 514 or 156 MA-turns in bundles of 10 or 3 m radius, about 4 or 1 Gt for the set;
  - 15 or 9 GW of refrigeration running through the night;
  - forces of 1.3–4.4×10¹⁰ N between circuits;
  - routes of 12,600–25,100 km in all, crossing latitudes from about 52° (or 15°) to 85°.

  The summit lies outside every tested route. The polar lands the routes cross are where both climate models put
  the most comfortable land, so routes and settlements meet there. The branch leaves "relief, heritage areas,
  settlements and sea crossings" unassigned.

  A straight-wire estimate of the bundle's field gives about 0.1 T at 1 km and 1 mT at 100 km for the 500 km
  circuits. The 0.5 mT line that MRI suites keep people with pacemakers behind would then lie 60–200 km from the
  routes, a corridor for land use, instruments and flyers. The magnets' night load (14.8 GW through half a lunar cycle
  is 18.9 PJ) asks for a regional grid or firm generation; the magnetic study says as much.
- **Hydrogen and the upper air.**
  - The escape model holds hydrogen at 0.53 ppm at the ground and loses about 0.0025 kg/s of it, 79 t a year.
  - At 0.53 ppm the air holds about 105 million tonnes of hydrogen.
  - The summit's fleet holds tens of thousands of tonnes: 31,000–54,000 t in the 750 m liners cycling through the
    port (2,100 t of it berthed at the crown at any time), and about 17,000 t in the sky ferries (about 17 t in each
    250 m ferry's 200,000 m³).

  Hydrogen lost from the fleet's cells therefore adds a small share to the air's own, and the loss budget of
  1–100 kg/s is far above the hydrogen's part in it. `atmosphere/loss_response/oxygen.py` can settle it by raising
  the ground's hydrogen.
- **Light at the ground.** The panels' spectrum follows the shield: the titania stack's window times the 5% dimmer,
  1,173 W/m² at the top of the air. Two choices still open would change the 22% panel estimate: the dimmer's form
  (even, or taken from the near infrared, or adjustable) and the ultraviolet edge (O3). Annulus X-rays deposit above
  the 0.3 Pa base, so flight heights receive none of them.
- **Failure.**
  - If gaps passed 0.35–0.82% of the Sun's light at quiet Sun, the upper air would swell into a state that holds
    itself, losing 140–190 kg/s. Every solar maximum on record settles a tenth to three-fifths of the way to that
    threshold at the standard level.
  - Where failed tiles and fragments end up is part of S7 and not yet traced. Tiles of 13–26 g/m² would slow high
    in the tall air and drift down.

### What the Moon's infrastructure asks of the shield

- Power delivery to the surface, if the author wants it, with its receivers and its heat placed on the Moon.
- A rule for launch exhaust and for anything lifted to orbit, in the spirit of S7, together with the six
  traffic-safety rules for craft leaving the Moon.
- The night sky's brightness under the lit stack.
- The magnets' routes laid on real terrain beside the settlements, with field corridors.
- Spare tiles and service craft based somewhere, and the freight to keep the film plant fed (2.3–4.6 Gt a year),
  which the transport link between the surface and orbit would carry if the feedstock comes from the Moon.

## 6. Decisions the other branches bear on

| Decision | What the other branches show | Standing |
|---|---|---|
| Long-haul sky ships dock in the flight band, so the tower rises 24 km (2026-09-27) | The band is above the median storm and inside the strongest tenth: flashes start at 35.8 km, at the crown's height. The band's winds are about twice as strong in the corrected climate | Stands. The crown's berths need storm procedures and lightning protection beyond Earth's highest level; for the author: whether that changes how ships berth at the crown |
| The round diagrid on six legs, and the floors as rings and disks (2026-09-27/28) | The corrected design gust gives 98 Mt of steel against 90 | Stands; the form is resized with the corrected winds and a design gust that varies with height |
| Open-air parks follow the air (2026-09-27) | The brief closes high parks for storms. Flashing storms give about an hour's warning, and storm hours over equatorial land rise in the corrected ring | Stands |
| Fusion supplies most of the Moon's power, with regional solar and wind (2026-09-27, planning assumption) | The shield's power stays in orbit so far. Wind devices in the port's upper frame give 140% of its use in the corrected winds. The polar magnets add 9–15 GW of steady load | Stands |
| A high-capacity backbone of metro, rail and sky ferries (2026-09-27) | Ferries are hydrogen ships under storms whose ground strikes pass through their air | Stands; for the author: the ferries hold or reroute round flashing cores while metro and rail carry on |
| Rigid hydrogen ships up to about 3 km in the band; ferries as rigid hydrogen ships (2026-09-28, findings) | Hulls barely change with the corrected winds; mooring loads rise by a quarter and the service wind nearly doubles. Lightning starts in the band | Stand as findings; hydrogen safety and a helium supply rise in priority |
| The flyers stack by height (2026-09-28, a first arrangement) | Thermals reach about 10 km with cloud base near 7 km in the corrected ring. Winds near 70 km now average 21 m/s, and the 200 m platform's envelope exceeds today's laminates | Stands as a first arrangement; for the author: the platforms' height, or drifting platforms |
| Practical dawn and dusk last about five hours (July canon, on this branch) | The sea-appearance branch's photometric definition (2026-10-04): 82 hours at the equator, 98 at 30° | Replaced at the integration by the sea-appearance branch's row |
| The protection architecture and the six traffic-safety rules (September) | The shield branch adopted S7 for the protection systems' outflows | Stand; launches from the surface are not yet under a rule |

## 7. Open questions, by lane

**Atmosphere (electricity).** Each depends on the electrified storms:
- the field at the top of a grounded 24 km conductor in the main run's charge structure, and when the tower would
  start flashes;
- storms over high ground, with an electrified box over the summit massif;
- the flashes' peak currents;
- whether craft start flashes;
- the runaway electrons' dose in the flight band;
- how solar storms reach long conductors on the ground (stage 3 and the plasma above).

**Climate.**
- **The summit's own weather in three dimensions.** Slope winds, rain on the windward slopes and the storms over the
  massif, which the speed-up factor of 1.2 stands in for. The highland box ran at 44.7° S, far from the summit.
- **Storm statistics from the corrected rings and boxes.** The studies' storm heights come from the first ring.
- **Comfort.** It stays a range between the models until the climate work reopens.

**Engineering.**
- **Lightning protection beyond IEC 62305's highest level:** charges to about 1,000–1,500 C, earthing in the
  highland's soils, the tower as a conductor, and hydrogen ships hardened against strikes.
- **The port's frame** resized with the corrected winds and a design gust that varies with height.
- **Wind devices** in the upper frame with the corrected winds.
- **The high platforms'** height or drift.
- **Hail and icing** on ships and the upper frame.
- **The crown's wind and wake**, already open, now in a stronger and steadier westerly.

**Habitation.**
- the city through roughly 190 dark hours a month, with lighting, activity and rhythm;
- a storm warning service;
- the backbone's rules in storms;
- civil time with the light calendar's phases;
- coastal settlements with monthly tides;
- daylight between towers, where the optical studies show that shade keeps about a third of the open light at noon
  (against a twentieth on Earth) but no street or block geometry is computed.

**Geography.** The highland lakes' waves and seiches, the Korolev lake first.

**Protection and the shield.** The planned `research/array-industry` branch takes up:
- power delivered to the surface;
- launches, exhaust and the transport link between the surface and orbit;
- the lit stack in the night sky;
- where failed tiles go;
- the magnets' routes, corridors and night load beside settlements;
- hydrogen from the fleet in the escape model.

## 8. Joining main

- **Climate products.** This review commits the corrected run's summit and whole-Moon winds:
  `climate/results/gcm/site_winds_A28_dim5_moon_summit.json` and `global_winds_A28_dim5_moon.json`. At the
  integration the four studies switch to them, to `ring_ring_equator.json`, and to main's drainage and sunlight
  products. Then they rerun (seconds each) and their READMEs take the new numbers. The table in section 1 shows what
  to expect.
- **Storm figures.** The studies and the metropolis brief quote the first ring's storm heights (28, 68 and 82 km)
  and the summit's September rain (20–24 mm a day). Both change: the corrected ring, the electrified box and the
  corrected climatology give the figures above.
- **Constants.** The products should record the shared constants they read by name and value (main's practice since
  4 October). The sky-fleet study reads two, the Moon's radius and the synodic month.
- **Files changed on more than one line of work.**
  - `research/decisions.md`, `research/status.json` and `research/README.md` on every branch;
  - `climate/gcm/README.md`, `habitation/README.md` and `visualization/README.md` on main;
  - `engineering/README.md` on the shield branch.

  Each holds additions from both sides.

## How the numbers were made

**Winds of the corrected run.** From this worktree, with main's GCM environment:

```sh
/path/to/main/climate/gcm/.venv/bin/python -m climate.gcm.site_winds A28_dim5_moon:20-29 summit 5.375 201.375
/path/to/main/climate/gcm/.venv/bin/python -m climate.gcm.global_winds A28_dim5_moon:20-29
```

They read `climate/gcm/runs/A28_dim5_moon/model/MOST.000{20..29}.nc` and take under a minute together.

**The studies with corrected inputs.** The studies were rerun on 8 October in two temporary copies of this branch
made with `git archive`:
- **The baseline copy** kept the September inputs and reproduced every committed result exactly, outside the
  producer hashes.
- **The corrected copy** took main's inputs under the names the code reads:
  - `ring_ring_equator.json` as `ring_ring.json`;
  - the corrected-run wind products as the `A28_dim5` ones;
  - main's `drainage.json` and `surface_light.json`;
  - main's atlas and drainage grids in `geography/products/`.

Each copy ran:

```sh
OPENBLAS_NUM_THREADS=1 python3 -m research.studies.summit_tower.run
OPENBLAS_NUM_THREADS=1 python3 -m research.studies.summit_tower.form
OPENBLAS_NUM_THREADS=1 python3 -m research.studies.port_fire.run
OPENBLAS_NUM_THREADS=1 python3 -m research.studies.sky_ships.run
OPENBLAS_NUM_THREADS=1 python3 -m research.studies.sky_fleet.run
```

The copies were measurement only; this branch's studies and results are unchanged.

**The rest.**
- The summit's rain comes from main's `climate/gcm/products/climatology_A28_dim5{,_moon}.npz`, at the four T21 cells
  round 5.4° N, 201.4° E.
- The lake shares come from the two drainage grids, weighted by area within great-circle distances of the summit.
- The strike counts take the box's rates over the 3,200 km² within 32 km of the tower.
- The twilight hours take the solar depressions of the sea-appearance branch's definitions at 5.4° N, with the Sun
  on the Moon's equator.
- The Sun at height takes the horizon's dip, arccos(R/(R + h)), without refraction.
- The hydrogen in the air takes the column of 1.2 atm under lunar gravity over the whole Moon at 0.53 ppm.
- The magnets' field is that of a straight wire, μ₀NI/2πr, near the bundle.

## Sources

Branches and commits read on 8 October 2026:
- main at 2f2e0a1: climate/gcm, climate/crm, climate/waves, geography (lakes and the monthly tide), illumination,
  the optical comfort and cloud twilight studies, and the decisions register.
- `research/atmospheric-electricity-plan` at 24d236b: the atmospheric electricity study, climate/crm's electrified
  storms, the aerosol study, and `climate/results/crm/elec_box_0e_elec_corrected.json` and its storm windows.
- `study/sea-appearance` at 46681f4: the sea-appearance study, the lighting calendar, earthlight and the decisions
  register.
- `research/solar-shield-habitat-array` at eeaecbc: the solar-shield study (its integrated comparison,
  `array_industry.md`, the magnetic architecture, the electromagnetic study and the design brief), the protection
  requirements, the loss response and the decisions register.

Earth references used here for lightning:
- IEC 62305-1:2010, *Protection against lightning. Part 1: General principles*: the highest protection level's flash
  charge of 300 C.
- F. A. M. Rizk (1990), "Modeling of transmission line exposure to direct lightning strokes", *IEEE Transactions on
  Power Delivery* 5(4): the potential at which a tall structure starts an upward leader.
- V. A. Rakov and M. A. Uman (2003), *Lightning: Physics and Effects*, Cambridge University Press: upward flashes
  from tall structures.
- V. Mazur (1989), "Triggered lightning strikes to aircraft and natural intracloud discharges", *Journal of
  Geophysical Research* 94(D3): aircraft starting the flashes that strike them.
