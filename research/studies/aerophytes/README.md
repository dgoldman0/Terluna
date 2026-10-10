# Photosynthetic aerophytes: the real limits and how a giant grows

The author asked on 9 October 2026 for giant photosynthetic aerophytes that grow over centuries as far as the physics
allows, with fittings, seams, gusts, damage, water and storms worked out, and with their water drawn from the air.
Aerophytes were called floaters until the author named them on 10 October 2026. The findings, all *screen* values
from [run.py](run.py) unless marked *derived*:

- **Size.** A single round aerophyte survives from about 50 m to 700 m across at 10 km with cellulose-class tendons, to
  1 km with 300 MPa fibre, and to 1.5 km and 3 km with those fibres in canopy-like light; collagen tendons, which
  creep at a tenth of their strength, hold it to 100–200 m. A colony of 60–150 m modules has no structural size limit,
  because each module's pressure head is set by its own height and the ties that hold the colony weigh 0.4–0.9 kg/m²
  at any width. Every hull carries its storm gusts, 100–300 Pa of spare pressure for trim and a wet fire-barrier skin,
  and every round body is divided into enough slack gas cells to survive losing one. Taking a burning cell's heat at
  one cell radius, the near end of its design range, moves the smallest colony module to 100 m.
  [Structure](structure.md).
- **Age.** A giant's age is the time its own light takes to make its gas. A 200-m body takes 12–304 years to grow, a
  1-km body 29–435 years as a colony and 73–1,186 years as a round body, and a 2-km body 36–544 years as a colony and
  225–702 years as a round giant, counted from a 40-m juvenile or one founding 60-m module. Round giants grow almost
  linearly in radius and slow as their tendons grow; colonies double every 3.6–55 years. [Growth](growth.md).
- **Water from the air.** Mean rain meets the reference need of 1.25 kg/m² a day equatorward of 35–38°, half the
  Moon; cool leaves and a CO₂-concentrating mechanism cut the need fivefold and widen that to 51°, and humid layers to
  65–69°. Rain comes with the afternoon storms every 15–30 days. A 2-km body's 300 Pa of spare pressure banks
  4.8 kg/m² of it, enough for the concentrating trait's 2–4 kg/m² between storms, while small and flat bodies must
  drink as they lose, from sorbent skins that take 0.2–2.4 kg/m² a day from the 64% air. [Water](water.md).
- **Storms and fire.** A soft aerophyte rides updrafts it cannot fight, so storms are escaped by route, hour and
  latitude; they rain in the tropical afternoon and flash only by day. Hydrogen burns at 4–71% in the lunar air, the
  98% lifting gas cannot burn, and a hole's flammable zone reaches at most 240–380 hole diameters, less for tears
  whose jets become buoyant plumes. A 200-m body meets 0.02–1.4 lightning strikes in 300 years and a 2-km body
  0.9–23, so giants are made of small gas cells or modules in wet skins that a strike burns one at a time.
  [Storms](storms.md).
- **Sailing and the Sun.** The mean wind blows east above 2 km, so no layer carries an aerophyte west with the Sun. A
  wing of a tenth of the body's area, hung on a weighted tether 5–9 km below it, sails it 0.4–0.8 m/s across the
  shear, 30–70° of latitude in a lunar day, for 0.03–0.7 kg/m² of hanging tether and ballast. Riding the eastward
  wind shortens the solar cycle to 21 days at 10 km and 15 days at 20 km at the equator. In the CM1 rings the low
  evening flow lingers near sunset for three to six weeks; a westward branch 5–10% weaker would cut that to 10–18
  days. [Sailing](sailing.md).
- **The canopy nursery.** Megaforest crowns can raise juveniles to 50–60 m, the smallest round bodies that float with
  their reserve, in one to seven years from fermented host sugar and their own light, with water and phosphorus drawn
  from the host's sap as mistletoes draw theirs. Beyond that the gas outgrows any host (a 200-m body's trimmed gas is
  26–243 hectare-years of forest growth), the juvenile's drag starts to threaten its tree, and the adult grows on its
  own light. [Growth](growth.md).
- **Gas biology.** Fish secrete oxygen into swim bladders at over 100 atmospheres, guanine plates make the bladder
  wall a hundred times tighter, and the Portuguese man-of-war makes its float gas in a gas gland; plated barriers
  could cut an aerophyte's permeation 10–50 times. [Gas biology](gas_biology.md).
- **People.** Aerophytes give membrane for gas cells, fibre, shade and roosts before food; their tethers, storm
  excursions and slow falls are the hazards. [Resources](resources.md).

Independent read-only reviews of the new components on 9 October 2026 found errors that changed headline numbers:
gas-cell walls counted as hydrogen-loss area, the spare superpressure and the fire-barrier water left out of the hull
and the mass, colony modules without the round bodies' seams, tethers sized without their own weight, venting and
water stores overstated, and periods, closures and lightning ranges read too coarsely. This version corrects them, and
the [checks](checks.json) list each; two further reviews checked the new sources against what could be read.

## What a giant looks like

The giant the physics allows is wide and flat, a sky reef. A colony of wet gas modules 60–150 m across (80–150 m at
20 km) lies one layer deep on a hexagonal raft at 10–20 km, a kilometre or two wide, its upper faces green with
photosynthetic tissue and its undersides furred with sorbent and fibre fringes that drink from the air. Weighted wings
hang on tethers kilometres below to sail it across the wind, and fringes comb cloud. Modules are specialised as a
siphonophore's zooids are: floats, collectors, gas-makers, sails and the reproductive modules that bud the next
generation. Storms and lightning take modules, and the colony regrows them as a tree regrows branches, so the colony
can grow for centuries. Round giants are the other form: single bodies of hundreds of metres to a few kilometres,
divided into eight to over a hundred slack gas cells inside a wet skin, slower and older, banking storm rain in their
tall gas columns. Migrating flyers roost on both and bring the land's and seas' phosphorus up to them.

## How the study is built

Run `python -m research.studies.aerophytes.run`; check with `python -m pytest research/studies/aerophytes`. The
[product](results/aerophytes.json) (schema `terluna.research.aerophytes/1`, which continues
`terluna.research.floater-viability/2`) binds 15 producer files, 13 source registers and the inherited products by
SHA-256 and records the shared constants its code reads. [Checks](checks.json) record what was run. The renderer is
[visualization/aerophytes/plot.py](../../../visualization/aerophytes/plot.py); its images stay outside Git.

| Component | Holds |
|---|---|
| [mechanics](mechanics.md) | Hydrostatic head, sphere and lobe screens, neutral trim, punctures |
| [envelope](envelope.md) | Hydrogen permeation, barrier measurements, inflation work |
| [biology](biology.md) | Upkeep, net production, night stores, hydrogen pathways, buds |
| [environment](environment.md) | Air, transpiration, droplets, leaf heat |
| [photoperiod](photoperiod.md) | Twilight, Sun-following geometry, light measures |
| [navigation](navigation.md) | Airspeed, drag work, wind tolerance |
| [trim cycle](trim_cycle.md) | Signed mass ledger, venting and compression |
| [structure](structure.md) | Joints, sustained load, gusts, damage, creep, colonies and the size limits |
| [growth](growth.md) | Ages, the canopy nursery, phosphorus, lifespan |
| [water](water.md) | Rain, vapour, dew, droplets, the trim rule and the water map |
| [storms](storms.md) | Updrafts, rain surges, flammability, lightning and fire |
| [sailing](sailing.md) | Shear, tethered wings, sea drogues, day length, sunset lingering |
| [gas biology](gas_biology.md) | Swim bladders, guanine barriers, the man-of-war's gas gland, siphonophores |
| [resources](resources.md) | Membrane, food, fibre, shade, hazards |

## Air and units

The study reads the sky-ship product's air; heights are above sea level. Per-area quantities use projected area
πR² (or a raft's hexagonal cell); membrane and gas-contact areas are separate; diameter is 2R.

| Height | Air density | Mean temperature | Lift per m³ of lifting mixture |
|---|---:|---:|---:|
| 10 km | 1.204 kg/m³ | 14.2 °C | 1.098 kg/m³ |
| 20 km | 1.024 kg/m³ | 5.1 °C | 0.933 kg/m³ |
| 40 km | 0.734 kg/m³ | −13.9 °C | 0.669 kg/m³ |

The air is dense and tall: at 10 km its pressure scale height is 52 km and its density scale height 62 km, and the
mean air reaches 0 °C at 25.7 km. The
lifting mixture is 98% hydrogen and 2% ambient air by mole; the inherited rounded density puts the hydrogen
inventory 1.05% below ideal xpV/RT at 10 km, and every ledger uses the density convention. The atmosphere holds about
400 ppm CO₂.

## The reference organism

**The reference is a light organism on a long-lived inert envelope.** Its envelope has a 5 MPa assembled skin (the
structure note derives it from 20 MPa coupons, 0.7 joints and a decade's sustained load), a 25 µm structural floor,
a separate 100 µm barrier at 0.173 Barrer, 1,500 kg/m³ materials, one gas–air partition, 20 Pa of base pressure and a
10 m/s gust (the coupled giants use the storm gusts of 14.4–17.0 m/s). Per projected square metre it carries 1 kg of
living dry tissue, 0.1 kg of community dry matter, 90% water in both, 5 kg of free water and its night reserve, and it
uses at most 70% of its full lift. It fixes 2 kg C/m² a year gross, keeps up living tissue at 0.003 kg glucose per kg
a day (a quarter of that in the dark), feeds 0.1 kg C/m² a year to its consumers, renews 10% of its inert structure and
half its living tissue a year at 1.39 kg glucose per kg, and makes hydrogen at 1% of a 200 W/m² mean light by day and
by fermentation at 2.1 mol/mol by night.

**At 10 km the reference sphere passes its mass, light and carbon gates at 40, 60, 100 and 150 m.** It supports
18.70 kg/m² (23.5 t) at 40 m and 24.59 kg/m² (193 t) at 100 m, mostly gas by volume:

| Barrier | Diameter | Share of full lift used | Construction carbon left | Parent-funded daughter time |
|---|---:|---:|---:|---:|
| 0.173 Barrer, 100 µm | 40 m | 0.639 | 1.029 kg C/m²/yr | 7.04 years |
| 0.173 Barrer, 100 µm | 100 m | 0.336 | 0.735 kg C/m²/yr | 16.14 years |
| 0.039 Barrer, 25 µm | 40 m | 0.619 | 1.076 kg C/m²/yr | 6.32 years |
| 10 Barrer, 100 µm | 40 m | 0.747 | −9.668 kg C/m²/yr | none |
| 10 Barrer, 1,000 µm | 100 m | 0.432 | −0.503 kg C/m²/yr | none |

With 10 kg of living dry tissue per square metre no size passes. The daughter time pays construction, reserves and the
gas photons from the parent's fixed area. The 40-m body's 1.75 t of hydrogen permeates at 0.062 kg a day, a 77-year
cycle time for its gas.

## The limits, in brief

**Size.** The low edge, about 45 m, is where a sphere's gas column, two thirds of its diameter, can no longer lift the
living tissue, the water store and the wet skin with a 30% reserve; even a massless envelope fails below 31 m, and
below 36 m with its wet skin. Round bodies' tendons take a share of the lift that grows in proportion to radius over
the fibre's breaking length, (5π/4)ρgR/σ, and renewing them every century eats the surplus at 0.3–3 km, depending on
the fibre and the light. Gas cells answer damage: a body needs supported mass ÷ dumpable water of them, never fewer
than eight; 22 at 1 km with 300 MPa fibre (52 with cellulose-class fibre) and 64 at 2 km. An arrested 4-m tear empties
one in minutes in a small body and in hours in a giant. Hydrogen replacement costs a steady 0.2 kg C/m² a year at
every size, because walls between cells lose no gas. Fatigue matters most for small bodies, whose stress swings most
with the gusts; a giant's stress is 97% steady, and creep at its renewal interval sets its allowable.

**Gas and age.** Photolytic hydrogen costs 3.8–9.5 kg of carbon per kilogram, fermented 8.9–21.9 kg. A colony's new
square metre costs 10.7–34.8 kg C and earns 0.89 kg C a year at the reference light, so it doubles in 8.4–27.3 years
(3.6–6.3 in canopy-like light).

**Water.** The reference needs 458 kg/m² a year, 37 kg per lunar cycle against 5 kg stored. Rain supplies the tropics;
vapour and humid layers supply the poles; venting gas to follow a 1 kg/m² daily swing would cost 104 W/m², so the
swing lives within the hull's spare pressure: 0.1 kg/m² in a 40-m body, 0.5 in a 200-m body and 4.8 in a 2-km body at
300 Pa.

**Storms and fire.** A body at 64% gas fill rises with an updraft until full at 37–41 km; holding against the storms'
90th-percentile updraft takes 3.5–9 kg/m² of extra weight. Venting to stay full costs 1.6% of the gas per kilometre of
forced rise. Rain must be shed at up to 70 kg/m² an hour. Oxygen permeating in reaches 5% of an unscrubbed module's
gas in 70–350 years, and a lining of hydrogen-oxidising microbes keeps it out for 4–18% more hydrogen.

## Twilight and the Sun

**Poleward of about 60° the clear ground never falls below 1 µmol/m²/s of PAR, and yet at 60° the reference ground
plant runs a carbon deficit for about 308 hours a cycle.** The PAR threshold is reached at 60.4°, 61.9° allowing the
declination range; at 60° the plant spends 34 hours below it. The diffuse twilight carries 5.44 W/m² of PAR with the
Sun 10° down and 1.11 W/m² at 20°. [Photoperiod](photoperiod.md).

**Holding the Sun takes 4.30 cos(latitude) m/s westward at 10 km, free on a matching wind and 90 W/m² for a sphere in
calm equatorial air.** A 0.5 m/s mismatch costs 0.141 W/m², 1 m/s 1.13 W/m² and 2 m/s 9.05 W/m² at 25% efficiency, and
bursts cost 1/d² more. With fixed light and production, the continuous-light case leaves 0.984 kg C/m² a year of
construction carbon against 1.029 for the half-cycle night. The project's winds offer no westward mean flow, so
aerophytes ride the eastward wind to shorten their nights, linger in the evening return flow near the ground, or sail
north and south between the rain and the twilight. [Navigation](navigation.md), [sailing](sailing.md).

## What would decide it next

1. **A joint wet membrane test.** Hydrogen, oxygen and nitrogen permeation, sustained-load strength and creep at the
   renewal interval, seams and joints, tear arrest and healing time, on one assembled living-compatible film, with and
   without guanine-like plates.
2. **The latitude and solar-hour wind product.** It places the westward branches, the shear and the rain sector along
   real paths, and replaces the global-mean profile and the five rings used here.
3. **A trajectory model.** Water, carbon, gas and trim along a sailing path through storms, twilight and humid
   layers, with the trim rule enforced hour by hour.
4. **A population model.** Recruitment, juvenile survival, fragmentation, lightning and storm losses, and the
   phosphorus cycle through roosts, litter and aerophyte falls.
5. **Storm records over centuries.** Return levels of gusts and updrafts beyond one lunar day of one ring.
