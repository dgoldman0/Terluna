# A sky tower on the summit

The author proposes the Moon's largest skyport as a very tall tower near its highest point, on the far side's
highlands: many layers for the city's life, docking for sky ships on the way up, a frame that lets the wind
through, and perhaps wind devices built into the frame (slow-turning blades, bladeless masts) to power the tower
and the city around it while easing the wind's load. This study answers three questions with light calculations:

- how the summit's wind, as the climate models now give it, bears on how tall the tower can be;
- how a frame that lets the wind through stands up to it;
- whether devices in the frame would give significant energy, and whether they ease the load.

The last section, [the central port](#the-central-port), sets out the design the author chose from these results.
The runner is [run.py](run.py), the numbers are in [results/summit_tower.json](results/summit_tower.json), and the
models are in [engineering/towers](../../../engineering/towers/README.md). Everything is closed-form sizing and
integrals over wind distributions; the whole study runs in about a second. Sources are in
[sources.json](sources.json).

## The site

- **The summit.** The atlas's highest cell, at 5.375°N, 158.625°W, stands 11,589 m above sea level (9,935 m above the
  geoid) on its 7.6-km grid. The cell averages 7.6 km of ground. LOLA's highest point, 10,786 m above the mean lunar
  radius near Engel'gardt crater (LROC 2010), rises above that average.
- **A broad rise of highland.** The summit tops a gently rounded rise of land several hundred kilometres across.
  Within 100 km of the summit the median ground is 8.6 km above sea level, and within 300 km, 6.5 km. All of it
  is land.
- **High lakes.** The drainage estimate's rain-fed lakes fill the basins around the rise. The largest nearby fills
  the Korolev basin: 141,000 km², up to 6.7 km deep, its surface 7,849 m above sea level. Its shore is 84 km from
  the summit and 3.7 km below it. A lake of 32,000 km² lies 156 km to the north, 6.1 km below the summit, and a
  small crater lake lies 9 km away.
- **The air.** The design run's free air at the summit's height averages 12 °C and 0.95 atm, with a density of
  1.16 kg/m³, close to Earth's at sea level. Ten kilometres higher it averages 3 °C, 0.78 atm and 0.99 kg/m³.
  With the design's oxygen (17.5%), a person at the summit breathes what they would at about 1,900 m on Earth,
  and at the top of a 10 km tower what they would at about 3,500 m.
- **Dense to the top.** The air thins slowly with height, so the whole tower stands in dense air.

## The wind at the summit

Two models give the wind, matched by height above sea level:

- **The design GCM.** Run A28_dim5, model years 15–24, from 3-day means of both wind components at T21. Its ground
  at the site is the 170-km cell's mean, 6.9 km above sea level.
  [climate/gcm/site_winds.py](../../../climate/gcm/site_winds.py) exports it as a product.
- **The cloud-resolving ring.** CM1 in two dimensions along the equator, flat at sea level, over two lunar days,
  sampled every 3 hours. It holds storms, but only the east–west wind.

The ring's distribution is fitted with a Weibull law, and its scale is raised for the north–south wind it lacks.
The factor is the cube root of the GCM's ratio of the speed's mean cube to the east–west component's, 1.03–1.22 at
these heights. The fitted law puts the 99th percentile about 10% below the ring's own, so its wind power is on
the low side.

| Above the summit | Above sea level | Air (kg/m³) | GCM 3-day: median / highest / 50-year | Ring 3-hourly: median / 99th pct / highest | Wind power (W/m²) |
|---|---|---|---|---|---|
| 0 km | 11.6 km | 1.16 | 2.0 / 5.5 / 6.0 m/s | 1.9 / 6.1 / 12.8 m/s | 16 |
| 5 km | 16.6 km | 1.07 | 1.7 / 4.4 / 4.7 m/s | 2.3 / 6.3 / 13.2 m/s | 27 |
| 10 km | 21.6 km | 0.99 | 1.9 / 5.3 / 5.7 m/s | 2.8 / 6.5 / 13.6 m/s | 27 |
| 20 km | 31.6 km | 0.83 | 3.4 / 7.5 / 8.1 m/s | 4.0 / 7.5 / 15.3 m/s | 39 |
| 30 km | 41.6 km | 0.70 | 4.7 / 8.7 / 9.3 m/s | 5.2 / 9.0 / 20.2 m/s | 64 |

**What the models show.** The wind is light: medians of 2–3 m/s, a 99th percentile of 6–7.5 m/s, and the highest
wind in either model 13–15 m/s over the first 20 km. Aloft the day and night differ little. Low down the wind
turns often: the GCM's mean vector is 40–60% of its mean speed below 5 km, rising to 90% at 20 km, where the
eastward flow dominates.

**The design gust.** The highest modelled wind over a 10 km tower is 13.6 m/s. The study raises it by three
factors for what the models cannot show:

- 1.4 for a three-second gust against a 3-hourly snapshot of 6-km columns;
- 1.2 for the flow speeding up over the summit's rise;
- 1.2 for going from two lunar days to a return period of decades.

That gives a design gust of 27.3 m/s. The service gust for sway is 19.0 m/s. A 50 m/s gust, an Earth-like design
value, is also run.

## How the wind bears on height

**The frame.** The frame is a square lattice of high-strength steel (S690, allowable stress 323 MPa). Its members
fill 12% of each face. Its drag follows the lattice-tower standard TIA-222-G for round members, with the wind
along a diagonal (force coefficient 2.06 on the member area).

**The floors.** Four storeys of 50,000 m² every 200 m give 1 km² of floor per km of height, at 500 kg per m² of
floor. The frame is at least 320 m wide at the top, so every layer holds its floors within half the frame's plan.

**The checks.** For each height and design gust the model tries 30 base widths and 4 flares. It keeps the lightest
frame that:
- sways less than 1/500 of its height in the service gust;
- has a whole-frame buckling factor of at least 3;
- has a base no wider than half its height.

Frame mass in megatonnes (the floors add 0.5 Mt per km of height):

| Height above the summit | No wind (sway check only) | 13.6 m/s | 27.3 m/s (design) | 50 m/s | Base at 27.3 m/s | Carbon composite at 27.3 m/s |
|---|---|---|---|---|---|---|
| 1 km | 0.014 | 0.014 | 0.016 | 0.030 | 0.44 km | 0.002 |
| 5 km | 0.37 | 0.42 | 0.56 | 1.2 | 1.3 km | 0.076 |
| 10 km | 1.6 | 1.9 | 2.9 | 6.2 | 2.3 km | 0.46 |
| 15 km | 4.1 | 5.0 | 7.9 | 17 | 4.0 km | 1.6 |
| 20 km | 8.3 | 10 | 17 | 37 | 6.0 km | 3.4 |
| 30 km | 24 | 30 | 52 | 114 | 9.0 km | none within the checks |

The same floors on Earth, with the U.S. Standard Atmosphere (NOAA/NASA 1976) and the same 19 m/s service gust for
sway:

| Height | No wind | 50 m/s |
|---|---|---|
| 1 km | 0.090 Mt | 0.099 Mt |
| 5 km | 3.7 Mt | 4.6 Mt |
| 10 km | 33 Mt | 40 Mt |

**The wind raises the frame's mass but sets no height limit up to 30 km.** At the design gust it adds 76% to a
10 km frame, 107% to a 20 km frame and 120% to a 30 km frame, over gravity alone. At 50 m/s the frames of 10–30 km
weigh 3.8–4.9 times the gravity-only ones. The effect grows with height because the air stays dense aloft and the
wind's moment grows with the square of the height.

**Wind sets the base width.** At the base of the 10 km tower, the design gust makes up 35% of each leg's force;
gravity makes up the rest. The service gust's sway is what decides how wide the base must be. Without wind, Earth's
frames of 7.5 km and more sit at the buckling minimum (factors of 3.4–3.5), where self-weight governs; the Moon's
gravity-only frames keep factors of 21–111.

**The Moon allows about 2.5 to 3 times the height for the same steel.** The 20 km frame at the design gust (17 Mt)
weighs about what Earth's 7.5 km frame does at 50 m/s (15 Mt). Earth's 10 km frame (40 Mt) weighs between the
Moon's 20 km and 30 km frames. Steel's self-weight length, the height of a column that its allowable stress can
hold up, is 25 km on the Moon and 4 km on Earth.

**Carbon-fibre composite frames** are five to eight times lighter. At 30 km, though, the lighter frame would need
a base wider than half its height to keep its sway within the limit.

## How to build it for the wind

**The reference tower** rises 10 km above the summit, with its top 21.6 km above sea level:
- base 2.3 km wide, flaring toward the ground (width 320 m + 2.0 km × (1 − z/H)²);
- 2.9 Mt of steel in the frame, and 5 Mt of floors carrying 10 km² of floor space in 50 layers;
- a first natural period of 76 s, and 17.7 m of sway in the service gust (1/565 of its height);
- a peak acceleration at the top of roughly 0.04 m/s² in that gust, about 4 thousandths of Earth's gravity;
- a whole-frame buckling factor of 21.5;
- footings of about 4,900 m² under each leg at 1 MPa.

For comparison, the usual guidance on tall buildings takes 5 thousandths of Earth's gravity as the threshold of
perception, and limits once-a-year peaks to 5–7 thousandths in homes and 9–12 in offices (Burton, Kwok and
Abdelrazaq 2015). The flared outline follows the same logic as Eiffel's: his tower's profile answers the wind's
moment, and its idealised form is exponential (Weidman and Pinelis 2004).

**An open frame is the largest saving.** The open lattice with its floor bands has a drag of 0.36 per unit of
outline area. A closed tower of the same outline has 1.3: 3.6 times the overturning moment. The lightest closed
tower of this height needs 8.2 Mt of steel, 2.8 times the open frame's.

**Vortex shedding.** The floor bands would shed vortices in step with the frame's first mode only at about 60 m/s,
4.4 times the highest modelled wind and well above the design gust. Setbacks that vary the cross-section with
height "confuse the vortices" a tower sheds (Irwin 2008), as on the Burj Khalifa; varied and perforated band edges
would do the same here.

**Stays and dampers.** Stays anchored in the highland around the base could narrow it; the model leaves them out. The tallest
guyed masts, such as the Warsaw radio mast (646 m, 1974–1991; Wikipedia), carried their lateral loads on stays. Dampers such
as Taipei 101's 660-tonne pendulum (the tower's operator) add to the frame's own damping.

**Cold, ice and lightning.** The top's air averages 3 °C, so rime will form in cloud at night. Storm cloud on the
ring tops out at 28 km at the median and 82 km at most, so the tower will be struck by lightning and needs the
protection Earth's tall towers have.

## Wind devices in the frame

**The resource.** The wind carries 16–39 W/m² over the first 20 km above the summit. The U.S. wind atlas counts
300–400 W/m² at 50 m as class 3, which it calls suitable for most turbine applications, and its lowest class ends
at 200 W/m² (Elliott et al. 1986).

**The layout.** For the reference tower, devices fill half of the frame's open face, 4.6 km², facing the wind.
Their output is reduced by two factors:
- 0.8 for the wind that goes around a porous tower;
- 0.7 for the mesh that keeps fliers out of the slow rotors.

**The devices.**
- **Fast rotors.** These follow NREL's 5-MW reference turbine: a power coefficient of 0.45 near its peak of 0.48,
  cutting in at 3 m/s (Jonkman et al. 2009).
- **Slow rotors.** These take the top of the published range for American many-blade rotors, 0.15–0.30 at a
  tip-speed ratio of 1–1.5 (Fraenkel 1986). Other authors put them nearer 0.15 (Ragheb and Ragheb 2011), so their
  output here is an upper estimate.
- **Savonius rotors.** These take FAO's 0.15.
- **Bladeless masts.** These take 0.05 of the wind's power on their frontal area, just below the 6% peak of a
  wake-oscillator model (Breen, Mallik and Adhikari 2025). No independent measurement of such devices exists in
  what was found; industrial units so far deliver 1–100 W.

| Device | Output per m² | Mean output | Share of the tower's use (114 MW) | Share of a city of a million at 2 kW each (2 GW) |
|---|---|---|---|---|
| Fast three-blade rotors (not safe for fliers) | 6.7 W | 31 MW | 27% | 1.5% |
| Slow many-blade rotors, screened | 3.1 W | 14 MW | 12% | 0.7% |
| Savonius rotors on vertical axes, screened | 1.6 W | 7.2 MW | 6% | 0.4% |
| Bladeless oscillating masts | 0.7 W | 3.4 MW | 3% | 0.2% |

**The tower's own use.** This assumes 100 kWh per m² of floor a year. That is 60% of the median US office's
167 kWh (ENERGY STAR 2024).

**Turning time.** Every device turns 45% of the time; the rest of the time the wind is below 3 m/s.

**Earth's towers with turbines.**
- **The Bahrain World Trade Center's** three 29-m turbines were forecast to supply 11–15% of its two towers'
  power and to run about half the time (Norwin and Atkins 2008). No measured output was found.
- **London's Strata SE1** built its turbines in, and they have barely turned (NBS 2011; Wikipedia, Strata SE1). The Pearl River Tower's
  design expected 1% of its energy from its turbines (Epstein 2008).
- **A UK trial of 26 building-mounted turbines** found output 15–17 times below what wind maps and power curves
  predicted. Some high-rise units were switched off after noise complaints (Encraft 2009).

**The devices add load.** A device that takes energy from the wind takes momentum from it too, and pushes that
load onto whatever holds it:
- **Stopped for the storm.** Slow rotors raise the lightest frame's steel by 45% (4.2 Mt against 2.9). Bladeless
  masts nearly double it (5.6 Mt).
- **Turning in the service gust.** The masts oscillate across the wind, which raises their drag coefficient from
  about 1.2 to about 2 (oscillating cylinders reach 2.5 at their largest swings; Griffin 1985). They then put 1.7
  times the open frame's design-gust moment on the tower.

The load reduction the author looked for comes from the open frame itself.

**Regenerative dampers.** One device both eases the wind's load and makes power: a tuned mass damper whose damping
is a generator. In a simulated 76-storey building it harvested hundreds of watts to kilowatts in mean winds of
8–25 m/s, while damping the building as well as a passive damper (Shen et al. 2018). That helps the structure, but
yields only a trickle of power.

**Fliers.** Band's collision model (Band 2012) is applied in a 4 m/s wind, to fliers crossing downwind at their
lunar airspeeds from the ecology register:

| Rotor | Moth, 1.5 m/s | Bird of 100 g, 4.7 m/s | Person on wings (assumed 6 m span, 6 m/s) | Blade tip speed |
|---|---|---|---|---|
| Slow, 12 blades, 10 m radius | 44% per crossing | 50% | 100% | 4 m/s |
| Fast, 3 blades, 40 m radius | 7% | 7% | 20% | 28 m/s |

**Slow rotors are no safer.** Slow many-blade rotors strike slow lunar fliers on most crossings, gently. Fast rotors
strike less often and lethally. Band's model shows why: part of the chance of a strike depends only on how much of
the disc the blades fill, whatever their speed. Lunar fliers are slow, so they spend long in the disc.

**What the numbers leave out.** Earth's birds avoid most rotors (98% or more for many species; Band 2012), which
these per-crossing figures do not count. Painting one blade black cut fatalities by over 70% at a Norwegian wind
farm (May et al. 2020). Screens make either kind of rotor safe.

## What could power it

**Panels.** Averaged over the lunar cycle, the design's clear-sky sunlight at sea level (202–1000 nm, from the
illumination product) gives 33 W per m² of level panel and 15 W per m² of east- or west-facing wall, at 22%
efficiency. The summit, with less air above it, gets more. On a square metre, then, panels give 5–11 times what
the best flier-safe wind device gives. Panels on half of the reference tower's east and west faces (9.9 km²)
would give about 150 MW, more than the tower's own use, all of it by day. A city of a million would need about
61 km² of level panels on the highland around the tower.

**Storage through the night.** The night lasts 354 hours. The summit stands 3.74 km above the Korolev lake, so
water pumped up by day and run back down at night delivers 1.35 kWh per m³ at 80% round trip. A plant on Earth
would need a 620 m head to deliver that much. The tower's night needs 30 million m³ lifted, and a city of a
million's night 0.53 km³; the Korolev lake holds 400,000 km³.

## What this leaves out, and what comes next

- **Winds over the real summit.** A terrain-resolving run over the summit's rise and the Korolev basin (CM1 in three
  dimensions, or a large-eddy model) would replace the speed-up and gust factors with the summit's own winds.
  Longer runs would give its extremes.
- **The structure.** Member and joint design, the dynamic gust response with dampers, stays, the order of
  construction and which materials can be made on the Moon, and footings in the porous highland crust. Large
  members in these winds are in supercritical flow, where the standard's round-member factor is lower (0.40
  against 0.57); the model keeps the higher one.
- **Docking.** Sky-ship loads, approach paths and the airflow around the layers where ships dock.
- **People.** Comfort in sway at 0.16 g: the same sideways acceleration tilts the apparent vertical six times more
  than on Earth. Also fire and evacuation from great height.
- **Power.** A power plan for the summit city: panels, pumped storage with the Korolev lake, and supply from the
  wider industrial network.

## The central port

The author's choices, on 2026-09-27:
- **Use.** The port is commercial: trade, travel, hospitality and short stays, with no homes and no industry.
- **Long-haul docks.** Long-haul sky ships dock in the flight band, so they never descend through the storms.
- **People.** The floors hold about 1.5 million people on a regular basis, in roughly the mix below.

`port()` in run.py sizes it.

**Height.** The named flight band begins 35 km above sea level
([decisions](../../decisions.md)), 23.4 km above the summit ground. The tower therefore rises 24 km, and its top
stands 35.6 km above sea level.

**Frame.** The floors follow the frame: up to 100,000 m² a storey, and 2% of the frame's plan where that is less, in
four storeys every 200 m. At the design gust the lightest steel frame for them has:
- a base 8.2 km wide, narrowing as (1 − z/H)^1.5 to 320 m at the crown;
- 33 Mt of steel in the frame and 17 Mt in the floors;
- a first period of 45 s, and 44 m of sway at the top in the service gust (1/550 of its height);
- footings of about 28,600 m² under each of its four legs at 1 MPa.

The industrial architecture moves bulk material at about 2.7×10⁸ kg a second, so the frame is two minutes of that
stream.

| Height above the summit | Width | Floor | Air at the top of the band | Oxygen like Earth at | Use |
|---|---|---|---|---|---|
| 0–3 km | 8.2–6.8 km | 6.0 km² | 0.89 atm, 9 °C | 2,430 m | Open public floors and terraces: markets, food, big event halls, the ground interchange, some hotels |
| 3–10 km | 6.8–3.8 km | 13.9 km² | 0.78 atm, 3 °C | 3,520 m | Enclosed floors with topped-up air: most hotels (pressurised above about 6 km), commerce, conferences, short stays |
| 10–15 km | 3.8–2.1 km | 10.0 km² | 0.70 atm, −1 °C | 4,290 m | Enclosed floors: regional docks and terminals, commerce, short stays |
| 15–24 km | 2.1–0.32 km | 5.0 km² | 0.59 atm, −10 °C | 5,670 m | Sealed, pressurised long-haul terminals, customs, lounges and transit rooms, under the crown's docking arms |
| Total | | 34.9 km² | | | |

**Programme and occupancy.** These are the areas per person assumed for each use:

| Use | Floor | Floor per person present | Present at a busy hour |
|---|---|---|---|
| Travel: terminals, docks, customs, lounges, ship servicing | 7.5 km² | 30–50 m² | 120,000–200,000 |
| Hotels | 7.0 km² | 30–45 m² | 124,000–186,000 |
| Short stay: compact and transit rooms | 3.0 km² | 12–20 m² | 120,000–200,000 |
| Commerce and trade | 7.0 km² | 12–20 m² | 245,000–408,000 |
| Shops, markets, food and drink | 4.0 km² | 4–8 m² | 250,000–499,000 |
| Events, conferences, leisure, observation | 3.5 km² | 3–6 m² | 233,000–466,000 |
| Services and back of house | 3.0 km² | 50–100 m² | 24,000–48,000 |

The busy hour counts each use at an assumed share of its capacity: 80% for travel, hotels and services, 70% for
commerce, 50% for shops and food, 40% for events. That gives 1.1–2.0 million people present. The author set the
design figure at about 1.5 million on a regular basis. With every space full, the tower holds 1.9–3.5 million, and
its exits must be designed for that.
- **Overnight.** The hotels hold about 116,000 rooms (186,000 guests), and short stays add 200,000 beds.
- **Travellers.** 120,000–200,000 are in the port at a busy hour. At two hours each, that is 60,000–100,000
  passengers an hour.

**Programme by zone.** The uses are placed by where people arrive and what air they need, then scaled so that each
use and each zone adds up (`fit_programme`, iterative proportional fitting). The floor is in km²:

| Use | 0–3 km | 3–10 km | 10–15 km | 15–24 km | Total |
|---|---|---|---|---|---|
| Travel | 0.6 | 0.4 | 2.5 | 4.0 | 7.5 |
| Hotels | 1.0 | 5.2 | 0.8 | – | 7.0 |
| Short stay | – | 1.0 | 1.6 | 0.4 | 3.0 |
| Commerce and trade | 0.5 | 3.5 | 3.0 | – | 7.0 |
| Shops, markets, food and drink | 2.0 | 1.4 | 0.4 | 0.2 | 4.0 |
| Events, conferences, leisure | 1.4 | 1.4 | 0.7 | – | 3.5 |
| Services and back of house | 0.5 | 1.0 | 1.1 | 0.4 | 3.0 |
| Floor | 6.0 | 13.9 | 10.0 | 5.0 | 34.9 |
| Present at a busy hour | 267,000–516,000 | 451,000–799,000 | 302,000–527,000 | 95,000–164,000 | 1.1–2.0 million |

- **Travel by trip.** Long-haul terminals (4.0 km²) fill the sealed zone under the crown's docking arms in the flight
  band. Regional docks and terminals (2.5 km²) sit at 10–15 km. Local interchange (1.0 km²) is at the base and in the
  hub bands.
- **Sleep.** Hotels sit low, where the air is best: like 1,950–2,430 m on Earth at the base, and about 3,000 m at
  6 km. Hotel floors above about 6 km are pressurised.
- **Crowds.** Markets, food and the big event halls fill the open floors at the base, next to the ground interchange.

The areas per person are rough planning figures from Earth practice, still to be checked against published
standards. Terminals take the most, because their docks, baggage and servicing are counted in.

**Power.** The floors use about 400 MW at 100 kWh per m² a year. The port's face is large: its open face is 77 km²,
and the wind carries about 1.95 GW through it. The Betz limit on that is 1.16 GW.

| Devices | Half the open face | The whole open face | The whole face, with the summit speed-up of 1.2 | Frame steel, devices stopped for storms |
|---|---|---|---|---|
| Screened slow rotors | 147 MW (37% of the tower's use) | 295 MW (74%) | 527 MW (132%) | 55.5 Mt (half face), 75.6 Mt (whole face) |
| Fast rotors (not safe for fliers) | 317 MW (79%) | 634 MW (159%) | 1,158 MW (291%) | |
| Bladeless masts | 35 MW (9%) | 70 MW (18%) | 125 MW (31%) | |

- **Wind.** Screened slow rotors over the whole face cover three-quarters of the tower's use in the modelled winds,
  and all of it if the summit's winds run a fifth stronger. The price is 43 Mt more steel, 2.3 times the open
  frame.
- **Panels on the tower.** Louvers over half the face that turn edge-on for storms add 14% to the steel (37.5 Mt).
  Solid panels there would add 160% (85.5 Mt).
- **Panels on the ground.** Level panels on the highland around the base need 12 km² for the tower's average use.
- **Storage.** Carrying the tower through the 354-hour night takes 105 million m³ of water pumped between the
  summit and the Korolev lake.

**A sky-ship world.** On the Open Moon lift is cheap, and the third dimension is ordinary access.
- **Building.** Sky ships can set frame segments anywhere up the tower.
- **Leaving.** Every layer's rim is an exit, by ship, by air taxi, on wings or under a canopy. A falling person
  reaches 18–22 m/s in the tower's air, against 43 m/s on Earth. A round canopy 3.7–4.5 m across lands them at
  4 m/s; on Earth that takes 8.9 m. The ride down takes 12 minutes from 3 km, 42 from 10 km and 100 from the top.
- **Fire.** Flight does not settle fire safety at 0.16 g. The [port fire study](../port_fire/README.md) works out
  how flames and smoke behave under weak buoyancy in the port's spaces, the compartments that give people time to
  leave, and what the tower's shafts, water, sealed zone and firefighting from the air ask.

**Form, chosen on 2026-09-27 and 2026-09-28.** The author chose a round diagrid that gathers into six splayed legs.
Its floors are rings in the lower tower and disks round a central core from about 12 km up, each a block of a few
storeys with a park and open-air shops on its roof, and the wide lower rings cantilever from the frame
([decisions](../../decisions.md#the-summit-port-and-metropolis),
[metropolis brief](../../../habitation/summit_metropolis/README.md)). The frame and floor figures above still
assume the square lattice and a band every 200 m; the tower model has yet to size the new form.
