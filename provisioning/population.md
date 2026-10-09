# Population and living across the horizons

The author's working scale of 9 October 2026: about 20–30 billion people across the Solar System by the frontier of
thousands to tens of thousands of years, affluent by modern standards and that the norm, most of them on Earth, on
the Moon and in the array, planned over the first thousand years and the several thousand after
([register](../research/decisions.md)). Numbers marked *screen* come from [population.py](population.py) and
[results/population.json](results/population.json); *derived* marks arithmetic in the text;
[population_sources.json](population_sources.json) lists every source and how far it was read. How people live
together, work and decide belongs to the human paper, *Living on an Open Moon*; this note gives the physical and
provisioning situation they live in.

## Interpretation corrected on 9 October 2026

The original A/B/C placements below are historical illustrative inputs. They did not evaluate aerial residential
capacity and provide no basis for preferring 16–20 billion residents on Earth or limiting the Moon to 2–5 billion.
The author reopened the distribution. New joint accounting uses
[shared population comparisons](../shared/scenarios/population.json) and the
[inhabited-volume study](../research/studies/inhabited_volume/README.md), with explicit surface and aerial residents.
No population allocation has been adopted. `population.py` and its saved product remain unchanged for reproduction.

Runoff-per-person thresholds below are terrestrial water-stress indicators, not physical population ceilings.
They do not model irrigation, reservoirs, fresh seas, reuse, regional redistribution or ecological allocation.
Likewise, heat-tolerance and crop-area rows are conditional screens, not demonstrated carrying capacities.

The [ways-of-living study](../research/studies/ways_of_living/README.md) takes up the reopened distribution with
eighteen settings built from the project's designs. Each is placed on the 28% atlas with its air, warmth, rain, light
and hazards, its costs per billion residents and the shared cases' totals by mix. There, water beyond the runoff
comes from the fresh seas for 3.2–9.5 GW a billion, and heat is answered by added dimming at 0.009–0.014% a
terawatt, so the runoff and heat rows below read as costs. How people live in each setting is in the
[daily-life study](../research/studies/daily_life/README.md).

## Findings

1. **20–30 billion takes growth at the pre-industrial world's pace or slower.** The UN puts today's 8.2 billion at a
   peak near 10.3 billion in the mid-2080s and 10.2 billion in 2100 (UN DESA 2024); IHME's reference scenario peaks
   at 9.73 billion in 2064 (Vollset et al. 2020). From 2100, 20–30 billion by the year 3000 takes 0.075–0.12% a year,
   by 5000 0.023–0.037% and by 12,026 0.007–0.011% (*screen*). The world grew 0.055% a year from year 0 to 1750 and
   0.18% from 1500 to 1750, against 2.04% in the late 1960s (UN 1999). Over a 30-year generation 0.12% a year is a
   net reproduction rate of 1.037 (*screen*), fertility just above the replacement 2.1 (UN DESA 2024).
2. **The Moon fills after the build.** From 10 million people on the Moon in 2526 (assumption), growth of 1–2% a
   year by births and arrivals reaches a billion in 2756–2987 (*screen*). Through the build the surface is an
   industrial site and protection flies before the air is delivered
   ([industrial architecture](../engineering/reference/industrial_architecture/RECOVERED_INDUSTRIAL_ARCHITECTURE.md),
   section 9), so the array's habitats and enclosed bases house the build's workers.
3. **Affluence is 2–6 kW of primary energy per person, centred on Europe's 3.7.** In 2025 the world used 2.3 kW a
   person, the EU 3.7, high-income countries 6.0 and the United States 8.6 (OWID 2026); the EU's final energy is
   2.6 kW (Eurostat). Decent living takes 0.48 kW of final energy, 0.41–0.58 by country (Millward-Hopkins et al.
   2020). All needs that Vogel et al. (2021) assess are met at 1.9 kW of final energy, and a doubling beyond it adds
   under 5%; at Europe's two-thirds of primary energy reaching its users, that plateau is 2.9 kW of primary energy
   (*derived*). The band runs from the 2000-watt society's 2 kW (Schulz et al. 2008) to the high-income average,
   with America's 8.6 kW as a high case; it gives the metropolis brief's 2 kW a basis.
   Computing at 0.1–1 kW a person in orbit (assumption) adds 2–30 TW, 40–640 times today's data centres
   (*screen*; IEA 2025).
4. **The chosen heat tolerance constrains uncompensated local energy use.** A terawatt used on the
   Moon puts 0.026 W/m² on it, 13.5 times what it puts on Earth (*screen*). At Europe's energy with industry in
   orbit a person releases 2.8 kW on the Moon, and a warming of 0.1 K holds 1.4–2.2 billion people, 0.5 K
   6.8–11 billion (*screen*). Each terawatt is answered by 0.009–0.014% more dimming; 1% more answers 72–115 TW,
   the heat of 26–41 billion such people (*screen*). The scenarios' 2–5 billion lunar residents need 0.05–0.19%,
   beside the 1.9–3.4% the fleet's glow needs ([first screen](README.md#heat-and-energy)).
5. **Where energy is used matters more than how it arrives.** Two-thirds of Europe's primary energy reaches its
   users as final energy (Eurostat), and all of it ends as heat where it is used. Moving industry to orbit, a quarter of
   final energy, cuts a person's heat on the Moon by a quarter. Beaming the rest down replaces the plants' losses
   with the receivers': 1.38 W of heat per watt delivered (*screen*, from the shield study's link efficiencies),
   a further tenth off (*screen*); for electricity alone the gain is larger, since a thermal plant converting 40%
   (assumption) releases 2.5 W per watt. Sunlight converted on the ground adds heat through the panels' albedo
   alone; fusion and any supply that adds to the sunlight absorbed count in full (Flanner 2009).
6. **Dense cities reach their limit before the planet does.** At the metropolis's 40,000 people per km² the band
   releases 60–240 W/m² over the city, by placement (*screen*); Earth's hottest city cells release 100–577 W/m²
   (Allen et al. 2011). The metropolis's 100 million over a 170 km GCM cell give 5–21 W/m² (*screen*), above the
   3 W/m² at which Flanner's model cells warm 0.15–0.24 K (Flanner 2009). The shield's dimming is global; a city's
   own heat is answered in the city.
7. **The runoff stress index crosses its historical threshold near 4 billion.** The Moon's rivers carry 7,240 km³ a year to the seas (*screen*,
   from the [drainage product](../geography/results/drainage.json)), 16% of Earth's 45,500 km³ (Oki & Kanae 2006),
   close to its dry land's 18% of Earth's 130 million km² of ice-free land (IPCC 2019). Falkenmark's stress line of
   1,700 m³ a person a year then falls at 4.3 billion and the scarcity line of 1,000 m³ at 7.2 billion (*screen*;
   Falkenmark et al. 1989); Earth reaches them at 26.8 and 45.5 billion. For scale, the average consumer's water
   footprint is 1,385 m³ a year and an American's 2,842 (Hoekstra & Mekonnen 2012). Rain-fed lakes cover 24% of the
   ground within 15° of the equator and 0.1% at 45–60° ([geography](../geography/README.md)). The seas and lakes
   start fresh and freshen slowly, by the resources domain's first screen (branch `domain/resources`), so while
   they stay fresh, water drawn from them and returned supplies beyond the runoff; lifting it to the mean land
   height of 2,345 m takes 1.06 kWh/m³ (*screen*), before any desalination.
8. **Food land is ample on the Moon by area, and its rain sets the rain-fed share; on Earth, food within planetary
   boundaries binds first.** The Moon's
   23.4 million km² of dry land feeds 8.4–22 billion at the first screen's 131–341 m² a person if crops take
   Earth's 12% share of land (*screen*); ruminants and flooded paddies are designed out, since methane lives for
   centuries in its air ([README](README.md#food)). Earth's 4.8 billion ha of farmland (FAO 2023) feeds 34–37
   billion on plant-rich diets of 0.13–0.14 ha a person and 4.4 billion on the American diet of 1.08 ha (*screen*;
   Peters et al. 2016; Poore & Nemecek 2018). A transformed food system feeds 10.2 billion within four planetary
   boundaries (Gerten et al. 2020); 16–20 billion on Earth are 1.6–2.0 times that, so part of Earth's food grows off
   the land, for instance microbial protein on electricity at over ten times any crop's protein per area (Leger et
   al. 2021), on energy the array can send.
9. **Earth's own heat is small on sunlight and reaches Flanner's 2100 scenario on thermal plants.** 16–20 billion
   at 3.7 kW on thermal plants put 0.12–0.14 W/m² on Earth, 0.09–0.11 K at AR6's 3 K for 3.93 W/m², and at 6 kW
   0.19–0.24 W/m² (*screen*; IPCC AR6). Flanner's 2100 scenario of 0.19 W/m² warmed continental regions 0.4–0.9 K
   (Flanner 2009); human forcing today totals 2.72 W/m² (IPCC AR6). Thermal plants reach that level at 26 billion
   people at Europe's energy and 16 billion at the high-income average (*screen*).
10. **The array's light is no limit; its habitats are its largest mass.** The whole population's energy, 40–260 TW
    across scenarios and band, is 0.005–0.03% of the 851 PW the fleet's night half intercepts and the Moon never
    needs ([joint ledger](../research/studies/joint_synthesis/results/ledgers.json)), on 0.1–0.9 million km² of
    collectors (*screen*). The Stanford torus shields each person with 990 t of lunar soil (Johnson & Holbrow
    1977), so 100 million residents take 25–99 Gt and a billion 248–990 Gt, the low ends for a denser design with a
    quarter of the torus's shield per person (assumption): 2–21 times the fleet's 46 Gt, or 4–42 days of the build's matter
    stream of 2.7×10⁸ kg/s (*screen*).
11. **By rocket, commuting weekly to the array takes five to twelve times a person's energy, and a rotation each
    lunar cycle one to three times.** The ring radius is 14.1 hours away on a minimum-energy transfer, 7.3 on a
    faster ellipse and 5.9 on a parabola, for 2.58–3.10 km/s from the surface in vacuum (*screen*). The orbital
    energy alone is 0.75 kWh/kg, 300–750 kWh for a traveller with 400 kg of craft and cabin (the sky fleet's
    long-haul share) to 1 t (assumption). On the minimum-energy transfer, hydrogen rockets with 0.5 km/s of losses
    in the air (assumption) take 2.95–7.7 MWh of electricity a round trip (*screen*; IEA 2023; L3Harris). Weekly
    commuting then uses 18–46 kW per commuter, a rotation each lunar cycle 4.2–11 kW and every other cycle
    2.1–5.4 kW (*screen*), against Europe's 3.7 kW. Lifting the traveller's 400 kg to the high platforms at 70 km is
    12.6 kWh (*screen*): lift is cheap, and orbital speed is the cost.
12. **Affluence lasts the horizons as a plateau.** At 2% a year from today's 18.9 TW, energy use reaches Flanner's
    2100 level on Earth in 82 years, 1 W/m² on Earth in 165 and all the light the fleet intercepts in 576–588
    (*screen*). Sustained affluence holds energy per person near the needs plateau, with population growth under
    about 0.1% a year.

## Historical constraint indicators

At Europe's 3.7 kW, with industry in orbit and thermal plants for the rest (*screen*):

| On the Moon | Billions | On Earth | Billions |
|---|---:|---|---:|
| Heat, unanswered, at 0.1 K | 1.4–2.2 | American diets on today's farmland | 4.4 |
| Runoff at the stress line, 1,700 m³ a person | 4.3 | Food within four planetary boundaries | 10.2 |
| Runoff at the scarcity line, 1,000 m³ a person | 7.2 | Thermal plants at Flanner's 2100 level | 26.3 |
| Crops on Earth's 12% share of land | 8.4–22 | Runoff at the stress line, 1,700 m³ a person | 26.8 |
| Heat answered by 1% more dimming | 26–41 | Plant-rich diets on today's farmland | 34–37 |

The Moon's population at a warming of the planned climate, in billions, with industry in orbit (*screen*):

| Energy a person | Released on the Moon | At 0.1 K | At 0.5 K | With 1% more dimming | City at 40,000 per km² |
|---|---:|---:|---:|---:|---:|
| 2.0 kW, the 2000-watt society | 1.5 kW | 2.5–4.0 | 12.6–20 | 48–76 | 60 W/m² |
| 3.7 kW, the EU | 2.8 kW | 1.4–2.2 | 6.8–11 | 26–41 | 111 W/m² |
| 6.0 kW, high-income | 4.5 kW | 0.8–1.4 | 4.2–6.7 | 16–26 | 181 W/m² |
| 8.6 kW, the United States | 6.5 kW | 0.6–0.9 | 2.9–4.7 | 11–18 | 259 W/m² |

## Scenarios

Historical placements of the author's 20–30 billion, in billions. These were selected without a carrying-capacity comparison; use the shared cases for new accounting:

| | Earth | Moon | Array | Elsewhere | Total |
|---|---:|---:|---:|---:|---:|
| A, Earth-centred | 16 | 2 | 0.1 | 1.9 | 20 |
| B, Earth and Moon | 18 | 4 | 0.5 | 2.5 | 25 |
| C, spread | 20 | 5 | 1 | 4 | 30 |

What each takes at Europe's 3.7 kW (*screen*):

| | A | B | C |
|---|---|---|---|
| Whole population's energy | 74 TW | 92 TW | 110 TW |
| Heat on the Moon, industry in orbit | 5.6 TW, 0.15 W/m² | 11 TW, 0.29 W/m² | 14 TW, 0.37 W/m² |
| Added dimming to answer it | 0.05–0.08% | 0.10–0.15% | 0.12–0.19% |
| Moon's runoff per person | 3,620 m³ | 1,810 m³ | 1,450 m³ |
| Moon's dry land in crops | 1.1–2.9% | 2.2–5.8% | 2.8–7.3% |
| Heat on Earth, thermal plants | 0.12 W/m², 0.09 K | 0.13 W/m², 0.10 K | 0.14 W/m², 0.11 K |
| Earth against food within boundaries | 1.6 times | 1.8 times | 2.0 times |
| Earth's runoff per person | 2,840 m³ | 2,530 m³ | 2,280 m³ |
| Stocks at 335 t a person, Earth and Moon | 6,000 Gt | 7,400 Gt | 8,400 Gt |
| Shield for the array's residents | 25–99 Gt | 124–495 Gt | 248–990 Gt |

The 335 t a person are developed countries' stocks in use (Vélez-Henao & Pauliuk 2023, from Krausmann et al.
2017), so the scenarios hold eight to eleven times the world's 792 Gt of 2010; decent living needs 42 t. Scenario C
puts the Moon's runoff below the stress line; A and B stay above it.

## Living

- **On the Moon.** The [ways-of-living study](../research/studies/ways_of_living/README.md) sets out the settings
  with their costs: the summit metropolis, coastal and lake towns, wetland-edge and canopy towns, the rain belt's
  farming country, highland towns and districts held up by towers or terrain, nearside and limb towns, the far
  side's twilit-night places, high-latitude towns, the undersea domes, and roaming, moored and streamlined sky towns.
  How many live in each is the author's choice. The metropolis brief's 70 m² of floor a person for home, work and services is three times the decent-living 24 m²
  (Vélez-Henao & Pauliuk 2023). A person using Europe's final energy, less industry, draws 650 kWh through each
  night, 590 times Earth's pumped storage per person (*screen*), so fusion, a day-side grid or power from orbit
  carries the night ([README](README.md#heat-and-energy)).
- **In the array.** Habitats fly natural orbits clear of the rings
  ([integrated comparison](../research/studies/solar_shield_array/integrated_comparison.md#inhabited-infrastructure)).
  Most people adapt to 4 rpm (Globus & Hall 2017): Earth's gravity at 56 m of radius at 4 rpm and 224 m at
  2 rpm, the Moon's at 9 and 37 m (*screen*), so one wheel can offer both. The torus shields to 2.5 mSv a year,
  5 with its factor-two margin, against the Moon's surface limit of 9.9 mSv (the register's 0.027 mSv a day,
  *derived*); it feeds each person from 20 m² of farm under continuous sunlight and gives 3 kW of electricity
  (Johnson & Holbrow 1977): 20,000 km² of farm and 3 TW for a billion residents (*screen*). The ring radius lies outside
  Earth's magnetosphere most of each month
  ([array industry](../research/studies/solar_shield_array/array_industry.md)).
- **Between them.** People spend about 1.1 hours a day travelling (Schafer & Victor 2000) and the array is half a
  day away, so its workers commute in rotations and the tinkerers who work at its industry most of the time live
  there. An ascent at the orbital minimum, by tether or elevator from altitude, brings a weekly commute to
  1.8–4.5 kW (*derived*).

- **Through the cycle.** The [daily-life study](../research/studies/daily_life/README.md) gives each setting's light
  through the 708.7-hour cycle. People keep a 24-hour day indoors, with dark rooms in the 354-hour day and
  daylight-level light through the night (11–45 kWh a person each night). The array's 0.13 s round trip falls inside
  conversation's 0.21 s gap and Earth's 2.56 s does not; the resource operation's crews work 14 minutes to 5.7 hours of
  light away.

## Earth service

- **Power.** From the Moon's distance a 10.5 km transmitter, about one tile, reaches a rectenna of the 1978
  reference size (10 × 13 km) at 2.45 GHz, and 4.4 km does at 5.8 GHz, against 1 km from geostationary orbit
  (*screen*; US DOE & NASA 1978; Rodenbeck et al. 2021). The link delivers 51% bus to bus
  ([shield study](../research/studies/solar_shield_array/results/electromagnetic.json)), on 4,800–7,200 km² of
  collectors per TW delivered (*screen*). Receiving land limits it: at the 1978 rectenna's 49 W/m² each TW takes
  20,400 km², 64 TW per 1% of Earth's ice-free land (*screen*), and a rectenna sees the Moon about half of each
  day. Criswell (2002) proposed 20 TW for 10 billion people from lunar power bases on the same principle.
- **Computing.** Light takes 1.28 s each way (*screen*), which suits training, batch work, simulation and
  archives; interactive services stay on Earth. A laser link sent 622 Mbit/s from lunar orbit in 2013 (Boroson &
  Robinson 2014).
- **Goods.** Freight down to Earth costs little energy and falls under the six traffic-safety rules
  ([register](../research/decisions.md#transport-and-safety)).
- **Earth's own arrays.** Earth needs no protection shield. A sunshade at Sun–Earth L1 blocking 1.8% of sunlight is
  about 20 Mt (Angel 2006), 1.6–3.2 days of the film plant's output (*screen*), and power arrays for Earth in
  geostationary orbit need transmitters a tenth of the Moon-distance size. The array can serve Earth with computing
  and supplementary power from its first centuries, and its film plant, photon-kept formations, collectors and links
  are the industry that would build Earth's own.

## What the plan needs

1. The author's choice of a distribution among the ways of living, weighed on their costs in the shared cases
   ([ways of living](../research/studies/ways_of_living/README.md)), with the energy band's sensitivities kept.
2. An answer to human heat on the Moon beside the glow's, counted by the heat each supply adds.
3. Regional domestic and agricultural water balances for every case, with recovery and supply alternatives; the 4-billion stress index is not a capacity limit.
4. A way up from the surface for people, and habitats in the array with their own orbits, shielding and food.
5. Food on Earth partly off the land if Earth holds more than about 10 billion.
6. Energy per person held at a plateau, and population growth near zero over the thousands of years.

## Open questions, in order of what they settle

1. **The distribution among the ways of living.** Each setting's costs and each case's mixes are in the
   [ways-of-living study](../research/studies/ways_of_living/README.md), and how people live in each is in the
   [daily-life study](../research/studies/daily_life/README.md); distributions for the author to choose among come next.
2. **How the Moon answers human heat.** The added dimming the shield can give beside the glow's; panel albedo for
   sunlight converted on the ground; a day-side grid's 9–19% losses; fusion's share. It constrains a chosen case's local heat and
   power system ([questions](questions.md), 1 and 3).
3. **The Moon's water beyond its runoff.** How long the seas and lakes stay fresh and how much can be drawn from
   them and returned, desalination's energy if they salt, farms against cities, groundwater, and settlement away from
   the tropics; it replaces the 4–7-billion stress-index shorthand with actual water budgets, with geography, resources and ecology.
4. **Heat in dense districts.** What a district of 40,000 people per km² may release and how it sheds it; it sets
   the metropolis's density and energy.
5. **People between the surface and orbit.** The ascent study, tethers or an elevator from altitude, and the
   rotations they allow ([questions](questions.md), 5); it sets commuters against residents.
6. **The array's habitats.** Orbits, shield per person by design, spin, food and life support; it sets the array's
   population and mass.
7. **Earth's service and Earth's own arrays.** Rectenna land, round-the-clock supply, computing bandwidth, and an L1
   shade or geostationary arrays for Earth, outside the shield's scope.
8. **For the human paper.** How rotations between the Moon and the array shape families, work and belonging;
   whether the array's residents form communities and nations of their own, as the Moon's nations form from
   Earth's; who decides the shared budgets of heat, dimming and water; how people across Earth, the Moon and the
   array meet and trade.

`python -m provisioning.population` rewrites the product; `python -m pytest provisioning/tests/test_population.py`
checks it.
