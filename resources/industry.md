# The array in industry, and the flows of the near horizons

On 9 October 2026 the author asked how the array connects to industry in general and to the resource flows of the
continuing Solar-System-wide extraction, and to plan on nearer horizons: the first thousand years with the build, and
the several thousand after. The population frame is 20–30 billion people across the Solar System, affluent by modern
standards, mostly on Earth, the Moon and in the array. This note places industry, sizes the flows after the build and
sets those horizons against the billion-year frame.

Numbers marked *screen* come from [industry.py](industry.py) and its product
[results/industry.json](results/industry.json); *derived* marks arithmetic on cited values; literature is cited as
(Author year), with every source and how far it was read in [industry_sources.json](industry_sources.json). Every flow
that drains an inventory carries its cycle time, the time the flow and the matching resupply take to replace it.

## Findings

1. **Concentrated industry flies on platforms of its own; the tiles carry distributed hardware.** The bundle's one sail
   loading, 94 g/m² per unit pressure coefficient, makes window tiles 26.0 g/m² and annulus tiles 13.07 g/m². An
   annulus tile holds its 10.21 g/m² film and 2.86 g/m² for structure, trim, store and wiring, 5.9 Gt over the fleet,
   of which the momentum store's rotors take 2–29% (*screen*). Photon keeping reaches a quarter of the sail force
   (integrated comparison), so a tile could carry at most 410 t beyond its match, paid from keeping of which
   eccentricity holding already uses up to 19–26% at the outer rings. Hardware spread over every tile moves the one
   loading together: each percent of the fleet's mass, 0.46 Gt, adds 0.023–0.046 Gt a year to the film plant, and each
   percent of the fleet's sunlight it absorbs puts 1.66–2.91 W/m² of glow on the Moon (provisioning branch, finding 1),
   the heat of 63–110 TW used there (*screen*).
2. **Where heat is released ranks the places across more than a factor of a million.** A terawatt used on the Moon adds
   0.026 W/m² to its mean; released evenly at the ring radius, 5.0×10⁻⁵ W/m², 530 times less, and 8.7–12×10⁻⁵ W/m²
   from the tiles' Moon-facing faces; at the Sun–Earth L1/L2 hubs, 1.1–1.9 million km from the Moon,
   0.6–1.6×10⁻⁸ W/m², another 3,100–8,900 times less (*screen*). The Moon's count applies to energy brought from outside
   its own sunlight, wind and rivers: fusion, fission and power beamed from orbit. Panels, turbines and dams on the Moon
   move heat its sunlight already brings, apart from the panels' darker albedo.
3. **The plant's tiles travel on their own sails, and plants spread over the stack's tilts keep few in transit.** The
   film plant renews whole tiles, 2.3–4.6 Gt a year, 3,850–7,730 exchanges a day; the rings tilt −23.8° to +24.1° and
   span 19,600–21,400 km. From one plant in the middle plane a leg takes 110 m/s at the median ring and 216 m/s at the
   edges; at a quarter of the tiles' 2.4×10⁻⁵ m/s² push the median leg takes 211 days, keeping 1.6–3.3 million tiles,
   6–12% of the fleet, in transit. Eight plants spread over the tilts cut that to 14 m/s, 27 days and 0.7–1.5%.
   Electric tugs at 30 km/s would release 535–1,070 kg/s of exhaust inside the protected region for one plant and
   69–138 kg/s for eight, outflows S7 governs (*screen*). The planned array-industry branch names one film plant; the
   eight plants are this screen's variant.
4. **Computing on the integrated ledger's scale is a gigatonne industry.** Today's servers weigh 12.8 kg per kW
   (130 kg for 10.2 kW; NVIDIA 2024), and their radiators 2.5–11.6 Mt per TW (array heat screen). 100 TW of computing
   is 1.5–2.4 Gt of hardware and 1,000 TW is 15–24 Gt, a third to a half of the fleet; on a 4–6 year server life (The
   Register 2022) they need 0.21–0.32 and 2.1–3.2 Gt a year, the film plant's size at the high case. With its
   collectors' heat each terawatt used in orbit releases 2.9–4.4 TW, so 1,000 TW near the fleet puts 0.15–0.22 W/m² on
   the Moon if it radiates evenly, and 5–7×10⁻⁵ W/m² from the hubs (*screen*). Control and interactive services fit
   beside the fleet, 0.13 s from the Moon and back; bulk computing fits the hubs, 7.5–12.6 s.
5. **The film's recovery fraction sets the shield's share of the operation.** At 99% recovery the plant takes
   732–1,470 kg/s of fresh feed and replaces the fleet's whole mass every 1,000–2,000 years; at 99.9%, 73–147 kg/s and
   every 10,000–20,000 years (*screen*). Earth recycled 4 of the 62 Gt it processed in 2005 (Haas et al. 2015); the
   plant's loop is one material stream. A thousand years of 99% feed, 23–46 Gt, lifted with the fleet from the airless
   Moon at 2.7 MJ/kg takes 20–40 GW over a century beside the fleet's own 40 GW (*screen*).
6. **On the near horizons the air needs no resupply except its CO2.** Escape takes 0.044–1.4 kg/s, cycle times of
   70 billion to 2 trillion years. Rock takes 165–5,270 kg/s of O2 (cycle 3.7–118 million years), 0.008–0.24% of the
   air's oxygen in years 1,000–10,000. Organic burial at lake rates would take 373–4,710 kg/s of nitrogen (cycle
   17–215 million years; Mendonça et al. 2017 with the C/N of Meyers 1994), 0.004–0.05% of it in those years. The CO2
   turns over in centuries: weathering at the corrected runoff takes 0.30–2.96 Gt a year and lake-rate burial
   0.74–1.87 Gt, so the air's 1.9×10¹⁵ kg lasts 400–1,860 years, 650–6,500 by weathering alone. Kilns returning the
   weathered carbon take 38–380 GW on the Moon (0.001–0.010 W/m²); the buried carbon needs a return path or
   6,400–16,200 kg/s of imported carbon, and its burial adds 0.54–1.36 Gt of O2 a year to the air (*screen*).
7. **The crust decides whether the water import ends with the build.** Escaping hydrogen takes 0.022–0.024 kg/s of
   water, a cycle time of 2×10¹³ years. If the crust's pores fill after the seas, 5.0×10¹⁸–5.9×10¹⁹ kg over
   1,000–10,000 years is 1.6×10⁷–1.9×10⁹ kg/s, 2–210% of the build's water stream (*screen*).
8. **The Moon's people build from lunar rock, import copper, and meet the heat budget through their energy.** At 43 t
   a person for decent living (Vélez-Henao & Pauliuk 2023) to 335 t in industrial countries (Krausmann et al. 2017),
   1–6 billion people hold 43–2,010 Gt, 1–57 mm of the dry land's soil, turned over at 0.53–31 Gt a year with building
   lives of 65–81 years (Berglund-Brown et al. 2025). Their copper at 140–300 kg a person (UNEP 2010) is 0.14–1.8 Gt,
   7–86% of Earth's identified resources (USGS 2025): 0.9–11 m of all the dry land's soil at the highlands' 4.6 ppm
   (Meyer 2010), or 1,100–14,200 Gt of CI-like asteroids at 127 ppm (Lodders 2003). Copper lasts 25–40 years in a
   building, so at 50% end-of-life recycling (UNEP 2011) the top-up is 1.8–36 Mt a year, up to 1.6 times Earth's 2024
   mine output; aluminium, the highlands' metal, substitutes in power cable (USGS 2025). At 119–270 GJ of primary
   energy a person (Germany to the United States; Energy Institute 2026) they use 3.8–51 TW: 0.10–1.35 W/m² and
   0.06–1.4 K if fusion or beams supply it (*screen*).
9. **The array's people set the largest continuing material stream.** SP-413's torus shields 10,000 people with 9.9 Mt
   at 4.5 t/m², 990 t a person, with 15 t of structure (Johnson & Holbrow 1977). For 2–14 billion people that is
   2,000–14,000 Gt of habitat, 45–406 fleets with their goods, and grown over 1,000–5,000 years, 12,700–446,000 kg/s.
   Water shields at least as well per tonne (Durante & Cucinotta 2011), and as water the whole habitat mass is
   0.003–0.1% of what the build delivers; lifted from the airless Moon it would take 1.7–12 TW for a century
   (*screen*).
10. **After the build the operation shrinks by a factor of 2,400–380,000, and the billion-year frame turns on return
    paths.** The build delivers 1.07–4.83×10⁹ kg/s, 11–48 packets of 10⁸ kg a second, braked at the hubs at
    19–70 MJ/kg: 20–340 PW, 4–67% of the 19 September civilization-scale envelope of 501 PW. Afterwards the air's
    top-up, the film's feed, the array's habitats and the Moon's copper come to 1.3×10⁴–4.5×10⁵ kg/s, 11–390 such
    packets a day and 0.25–31 TW of braking; the crust's water and imported carbon would hold it larger (*screen*).
    Over 10⁹ years the rock takes 6.5–206 times the air's oxygen, lake-rate burial 4.7–59 times its nitrogen and the
    film's fresh feed 0.7–15 air masses (*screen*; first screen).

## Places

| Place | To leave or reach | Sunlight | Heat on the Moon per TW | Light time to the Moon | At hand | Suits |
|---|---|---|---|---|---|---|
| The Moon's surface | 2.82 MJ/kg to escape, 2.70 to the ring orbit; after the air, 75 t of air per m² | 293 W/m² at the top of the air | 0.026 W/m² from outside energy | — | Highland Al, Ca, Si, O; mare Fe, Ti before flooding; water, air | Food, buildings, the Moon's goods and stocks, kilns |
| Platforms near the fleet, 20,000 km | 0.12 MJ/kg to lunar escape | 1,361 W/m²; the night half intercepts 851 PW | 5.0×10⁻⁵ W/m² even; 8.7–12×10⁻⁵ from tile faces | 0.067 s | Worn tiles, 2.3–4.6 Gt a year | Film plant, spares, collectors, control computing, habitats |
| Sun–Earth L1/L2 hubs | 14–45 m/s of phasing to a lunar encounter (Chen et al. 2016); about 1 m/s a year to hold (Roberts 2011) | 1,361 W/m² | 0.6–1.6×10⁻⁸ W/m² | 3.7–6.3 s | What arrives | Capture and braking, packets, bulk computing, heavy industry for Earth and Moon |
| Earth | 62.6 MJ/kg to escape, 33.1 to low orbit | — | — | 1.28 s | Everything, at the deepest well | People, knowledge, precision parts |
| Main belt, 2.7 AU | 4.8 km/s out; arrives at 6.2 km/s, 19 MJ/kg | 187 W/m² | — | 14–31 min | Silicates, metals, carbonaceous bodies | Metals, glass feed |
| Jupiter's trojans, 5.2 AU | 5.6 km/s out; 8.8 km/s, 39 MJ/kg | 50 W/m² | — | 35–52 min | A candidate volatile population (question 1) | Volatiles |
| Kuiper belt, 40 AU | 3.7 km/s out; 11.8 km/s, 70 MJ/kg, 46 years | 0.85 W/m² | — | 5.4–5.7 h | Ices, bound nitrogen | Volatiles refined at the source, hydrogen kept there |

Transfers are coplanar Hohmann orbits to 1 AU without gravity assists (*screen*).

## Industries by place

| Industry | Place | Why |
|---|---|---|
| Film plant and spare tiles | Platforms near the fleet, several across the tilts | Its feed is the fleet; heat stays in orbit; vacuum for evaporation and deposition |
| Collectors and radiators | Platforms near the fleet and at the hubs | The night half's light near the fleet; radiators edge-on to the Moon |
| Computing and its hardware | The hubs for bulk work, the fleet's platforms for control | Heat 3,100–8,900 times lower at the hubs; 0.13 s round trip near the fleet |
| Habitats | Natural orbits near the fleet | A massive habitat keeps an orbit of its own (design brief) |
| Refining of imported volatiles | The source bodies | The hydrogen stays at the source (the domain's finding 2); heat is rejected there |
| Metals and glass for the Moon | The Moon | Highland aluminium and calcium; mare iron and titanium mined before flooding |
| Metals and glass in orbit | The airless Moon, then asteroids through the hubs | After the air, launches cross 75 t/m² |
| Capture, braking and repackaging | The hubs | Major braking happens remotely (traffic rules); the recorded architecture's freight hubs (S5) |

## Flows after the build

| Flow | Rate | Cycle time | Years 500–1,000 | Years 1,000–10,000 | Over 10⁹ years |
|---|---|---|---|---|---|
| Air to space | 0.044–1.4 kg/s; designed 1–100 | 70 billion–2 trillion years | 0.7–22 Mt | 13–400 Mt | 0.04–1.4% of the air |
| O2 to rock | 165–5,270 kg/s | 3.7–118 million years | 0.0004–0.014% of the O2 | 0.008–0.24% | 6.5–206 air-oxygens |
| N burial, lake rates | 373–4,710 kg/s | 17–215 million years | under 0.003% of the N2 | 0.004–0.05% | 4.7–59 air-nitrogens |
| CO2 to weathering and burial | 1.0–4.8 Gt a year | 400–1,860 years | 27% to all of the air's CO2 without a return | — | — |
| Water as escaping hydrogen | 0.022–0.024 kg/s | 2×10¹³ years | 0.4 Mt | 7 Mt | 0.005% of the water |
| Crust's water, if late | 1.6×10⁷–1.9×10⁹ kg/s | ends when full | — | — | — |
| Film fresh feed, 99–99.9% | 73–1,470 kg/s | fleet replaced in 1,000–20,000 years | 1.2–23 Gt | 21–420 Gt | 0.7–15 air masses |

## People's stocks

Scenarios of this screen within the author's frame; with Earth at its projected peak of 10.3 billion (UN 2024) they
total 13, 20 and 30 billion. The provisioning domain's scenarios (branch `domain/provisioning`) place more people on Earth,
16–20 billion, and 0.1–1 billion in the array, whose habitats then come to 99–990 Gt; the author's placement settles
which holds.

| Scenario (Moon, array) | Moon's stocks | Turnover | Moon's energy, and warming if brought in | Copper top-up | Array's habitats | Grown over 1,000–5,000 years |
|---|---|---|---|---|---|---|
| 1 and 2 billion | 43–335 Gt | 0.53–5.2 Gt/yr | 3.8–8.6 TW, 0.06–0.23 K | 1.8–6 Mt/yr | 2,010 Gt | 12,700–63,700 kg/s |
| 3 and 7 billion | 129–1,000 Gt | 1.6–15.5 Gt/yr | 11–26 TW, 0.19–0.68 K | 5.3–18 Mt/yr | 7,040 Gt | 44,600–223,000 kg/s |
| 6 and 14 billion | 258–2,010 Gt | 3.2–31 Gt/yr | 23–51 TW, 0.37–1.4 K | 10.5–36 Mt/yr | 14,100 Gt | 89,200–446,000 kg/s |

## The horizons

| Flow | Kind | Years 0–500 | Years 500–1,000 | Years 1,000–10,000 | 10⁹ years |
|---|---|---|---|---|---|
| Air and seas | Ends | 1.1–4.8×10⁹ kg/s; 20–340 PW braked | — | — | — |
| Crust's water | Ends, timing open | in the build, if fast | 1.6×10⁷–1.9×10⁹ kg/s, if slow | — | — |
| Fleet | Built, then renewed | 46 Gt lifted, years 150–350 | 2.3–4.6 Gt a year; 73–1,470 kg/s fresh | same | 0.7–15 air masses fresh |
| Carbon | Continues | the biosphere's 1.5–1.6×10¹⁵ kg of CO2 | kilns of 38–380 GW, or 0.9–4.2×10⁴ kg/s imported | same | a return path |
| Air's other losses | Continue | — | under 0.014% of any gas | under 0.24% | 6.5–206 oxygens, 4.7–59 nitrogens |
| Moon's stocks | Grow, then turn over | — | toward 43–2,010 Gt | 0.53–31 Gt a year | turnover |
| Array's habitats | Grow with population | — | 1.3×10⁴–4.5×10⁵ kg/s | same | small renewal |
| Computing | Grows with demand | — | 0.21–3.2 Gt a year | same | renewal |
| Copper | Grows with population | — | 1.8–36 Mt a year | same | recycling |

The airless window closes as the air builds, years 250–480 in the protection report, so whatever orbit needs from
lunar rock is cheapest before then.

## Open questions, in order of what they settle

1. **The array's population and its habitats' shielding**: passive mass at 4.5 t/m², magnetic or water shields, and
   the area each person has. This sets the largest stream after the build.
2. **When the crust takes its water** (the domain's question 10). Late filling keeps the operation near build scale for
   millennia.
3. **The carbon's return paths** (question 6): kilns for the weathered carbon, a return for buried carbon, or an import.
   This is the near horizons' one fast drain on the air.
4. **The film plant's recovery fraction, number and siting** (question 5, with the planned `research/array-industry`
   branch): how tiles travel between rings and plants, and the share of S7 their traffic may use.
5. **Computing demand and the mass per terawatt of hardware built for orbit**, which decide whether computing becomes
   the array's largest industry.
6. **The Moon's energy sources.** The share of fusion and beamed power sets how much of the Moon's own use the planned
   climate takes.
7. **Element-by-element imports scaled by population** (question 11): copper and zinc first, with aluminium's reach as
   a substitute.
8. **The hubs' catchers and energy recovery for 20–340 PW of braking**, and S7 extended to industrial platforms near the
   fleet.
9. **Pre-air stockpiles**: feed reserves and habitat material mined from the future seafloor before it floods, within
   the conservation programme's gates.

## Files and running

`python -m resources.industry` writes [results/industry.json](results/industry.json) (schema
`terluna.resources.industry/1`) in under a second from the first screen, the array heat screen, the integrated ledger
and ring layout, the joint ledgers, the geography atlas and drainage, the CO2 study and the 19 September calculations;
`python -m pytest resources/tests/test_industry.py` checks that it binds them and that its arithmetic holds. The CO2
study's weathering is rescaled to the corrected runoff (factor 0.822, as the first screen rescales the oxygen sink), so
it runs a fifth below the 0.36–3.6 Gt a year of the domain's finding 5.
