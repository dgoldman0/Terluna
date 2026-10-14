# The seas' make-up through time

13 October 2026. Resources' question 9 asks how the seas' composition changes through time, with the delivered water's
salts a design lever. The author asked for the model with fresh, salted and phosphorus-added water, and for it to look
past the first thousands of years to the stable or dynamic state the seas could hold over the world's billion-year
life. This box model follows the common-level waters and the rain-fed lakes that feed them from the water's delivery
to a billion years.

`python -m research.studies.sea_chemistry.run` writes [results/sea_chemistry.json](results/sea_chemistry.json)
(schema `terluna.research.sea-chemistry/1`, with its producer, the sha256 of its eight inputs, its evidence statement
and reading rule) in about fifteen seconds: 54 runs, each a delivery, a set of sinks, a fate for the sediment and one
of the carbon cycle's three weathering corners. [sources.json](sources.json) records each published value and how far
it was read.

## What it finds

- **The first thousands of years follow the carbon cycle.** Fresh water starts at pH 5.6 and takes weathering's
  carbonate: pH 8.1 by year 1,000 and 8.7 by year 10,000, at 0.13 and 0.43 g/kg. The seas reach three times calcite
  saturation in year 1,320–2,740 (central 1,750), where the
  [carbon cycle](../carbon_cycle/README.md) finds 1,306–2,725 (central 1,742). Brackish water, Earth's seawater diluted
  tenfold, moves from pH 7.4 to 8.5 over the same time, and salted water stays near 35.1 g/kg and pH 8.0 while its
  calcium falls a sixth as calcite settles.
- **Left alone, the seas turn alkaline.** The Moon's rock gives the rivers calcium, magnesium and a little sodium and
  potassium, each with the carbonate that balances it, and almost no chloride or sulphur. Calcite takes the calcium,
  and the rest of the alkalinity stays in seas that have no outlet. Fresh seas reach pH 9 in 13,000–34,000 years.
  Magnesium clays forming at pH 9 hold them near it for about a million years, and the sodium's soda then carries them
  on: 1 g/kg after 0.5–0.6 million years, Earth's salinity as a soda water after 24 million and a saturated soda brine,
  like Earth's soda lakes at about 300 g/kg, after 200 million. Without the clays, the magnesium keeps its carbonate in
  solution and the same path runs about fifteen times faster.
- **Delivered salt sets only the start.** Chloride has no sink, so a salted sea keeps Earth's salinity, but
  weathering's alkalinity gathers in it all the same: pH 9 after 0.21–0.30 million years and a tenth saltier after
  0.17–0.24 million. Brackish seas join the fresh seas' path within ten million years.
- **Sediment fills the waters first.** The lakes catch nearly all the land's sediment and are full in 0.2–6.8 million
  years. The seas take the particles that pass and their own deposits of calcite, opal and magnesium minerals,
  0.7–2.3 km³ a year, and are full in 0.7–12 million years. A Moon without plate tectonics has nothing that renews its
  basins. The [terrain study](../terrain_through_time/README.md) routes the rock lake by lake: at first the lakes take
  38% of it, dry hollows 41% and the seas 21%, so the lakes' capacity is 99% filled in 2–25 million years and the
  seas in 0.8–11 million.
- **Every lasting state is a held one.** Four flows at weathering's catchment rate keep the seas wherever the
  delivery sets them, for the billion years:
  - calcite, 296 Mt a year, which calcifiers build as they do on Earth;
  - magnesium, 126 Mt a year of magnesite made where the rivers enter the lakes, the carbon cycle's managed route, or
    242 Mt a year of magnesium clays at pH 9;
  - soda and sodium sulphate, 20 Mt a year from soda works on the rivers;
  - silica, 476 Mt a year of opal, which silica-building life holds at 0.1 mmol/kg.

  Held this way, fresh seas stay hard fresh water at 0.14 g/kg and pH 8.3, brackish seas at 3.65 g/kg and pH 8.3, and
  salted seas Earth's seawater at 35.1 g/kg and pH 8.0. Open water for the billion years also needs the land's
  sediment, 0.7–8.3 Gt a year, returned to it, which the conservation study's gates price at 0.11–1.1 GW of ideal
  lifting power for 0.71–7.1 Gt, and the seas' own deposits, 0.7–2.3 km³ a year, taken out.
- **Phosphorus follows the lakes' trap.** While the lakes catch it, and dredging their deltas for people's phosphorus
  keeps them catching, the mixed seas hold 0.4–23 mg/m³ fresh and 0.9–116 salted, against Earth's ocean mean of 70;
  salt water frees part of the particles' phosphorus. Once the lakes are full the seas get the land's whole loss:
  7–350 mg/m³ fresh and 20–2,600 salted. Earth's mean delivered with the water falls to the rivers' level by year
  10,000 at a 1,000-year residence and by year 100,000 at 20,000 years. An alkaline sea with its calcium spent can
  hold far more, since apatite no longer takes the phosphate (Toner & Catling 2020).
- **The carbon.** Weathering takes 0.109 Gt C a year from the air at the catchment rate. Calcite buries 0.036 Gt C a
  year, and the magnesite route 0.018 more; the rest returns to the air. The air's CO2 is held at 400 ppm here, so
  the resupply the carbon cycle sizes for its first 10,000 years runs on at that rate.

## The seas over a billion years

Central weathering, seas kept open. Salinity in g/kg and pH; a dash where the brine is past the constants' range.

| Delivery and sinks | 1,000 years | 10,000 | 100,000 | 1 million | 10 million | 100 million | 1 billion |
|---|---|---|---|---|---|---|---|
| Fresh, natural (magnesium clays) | 0.13, 8.1 | 0.43, 8.7 | 0.81, 9.0 | 1.6, 9.1 | 15, 9.8 | 153, – | 291, – |
| Fresh, magnesium kept | 0.13, 8.1 | 0.43, 8.7 | 2.5, 9.4 | 23, 10.0 | 221, – | 298, – | 298, – |
| Fresh, held | 0.06, 7.9 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 |
| Brackish, natural | 3.6, 8.1 | 3.9, 8.5 | 4.4, 9.0 | 4.8, 9.0 | 18, 9.8 | 156, – | 294, – |
| Brackish, held | 3.6, 8.0 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 |
| Salted, natural | 35.1, 8.0 | 35.2, 8.0 | 36, 8.6 | 40, 9.0 | 49, – | 187, – | 293, – |
| Salted, held | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 |

| Milestone, years | Fresh, natural | Fresh, magnesium kept | Salted, natural |
|---|---|---|---|
| Calcite settles in the seas | 1,320–2,740 | 1,320–2,740 | from the start |
| pH 9 | 13,000–34,000 | 13,000–34,000 | 211,000–298,000 |
| 1 g/kg | 501,000–631,000 | 17,000–40,000 | |
| Earth's salinity, or a tenth above it | 24 million | 1.5–1.7 million | 168,000–237,000 |
| A saturated soda brine (300 g/kg) | 200 million | 13–14 million | 178 million |

Ranges run across the weathering corners. The natural seas' composition after a few million years is a soda water:
sodium and potassium with carbonate, calcium below 0.01 mmol/kg and magnesium in clays. Without the clays it is a
magnesium-soda water.

## Sediment and deposits

| | Earliest | Latest |
|---|---|---|
| The land's erosion | 8.3 Gt a year (the soil loss the phosphorus particles imply) | 0.71 Gt a year (10 m/Myr, the conservation study) |
| Dry density of what settles | 500 kg/m³ | 1,500 kg/m³ |
| The seas' deposits | 2.3 km³ a year (natural seas, 80% porous) | 0.68 km³ a year (held seas, 50% porous) |
| Lakes full | 196,000 years | 6.8 million |
| Seas full | 0.73 million years | 12 million |

Without their own deposits the seas would be full in 0.83–29 million years. When the seas fill, their water spreads
over the low ground the sediment builds: the volume stays, and the deep basins go.

## How it is set up

- **The waters.** The common-level waters, 10.7 million km³, and the rain-fed lakes, 3.1 million km³, as well-mixed
  boxes at the design case's sea surface (297.0 K) and land air (294.3 K), each in CO2 equilibrium with the air at
  400 ppm and 1.216 bar. The lakes take 96.5% of the land's runoff and pass it on through 445 years of residence.
- **Weathering**, imported from the carbon cycle: Dessert's basalt law over the land's discharge, 9.1×10¹² eq a year
  at the catchment rate, with the developed soils starting at 1, √10 or 10 times it and ageing toward it on White and
  Brantley's law; erosion keeps the catchment rate renewed. The highland soil sets the cations: calcium 65% of the
  alkalinity, magnesium 33%, sodium and potassium 2%. Its export gives 0.87 mol of silica and 0.0025 of sulphate for
  each equivalent.
- **The chemistry**, imported from the carbon cycle: Plummer and Busenberg's constants with Davies activities in fresh
  water, Earth's seawater constants in the salted sea, and its carbonate solver. Calcite settles where a water stands
  above three times saturation, or a salted sea's delivered 4.37. Magnesium clays form at pH 9 (*design guess*;
  Tosca and Masterson 2014 find them forming from alkaline, silica-rich water at 25 °C) with one and a half silicons
  for each magnesium, and their carbon returns to the air. Amorphous silica holds the silica without
  silica-building life.
- **The deliveries.** Earth's reference seawater (Millero et al. 2008) at full strength or a tenth of it, or pure
  water, with Earth's ocean mean of phosphorus as a variant.
- **Past the constants' range** (an ionic strength of 0.5 mol/kg in fresh water, 43 g/kg in a salted sea) the model
  keeps the masses alone: calcium, magnesium (with clays) and silica stay where the minerals held them, and at
  300 g/kg soda, sodium sulphate and magnesium carbonates deposit what arrives (*design guess* on Earth's saturated
  soda lakes).
- **Phosphorus**, the sea phosphorus study's mixed box: the dissolved input that passes the lakes, the land's whole
  dissolved loss once they are full, and in a salted sea 10–30% of the particles' phosphorus as well, with burial over
  1,000 or 20,000 years.

## What it leaves out

- **Ion pairs and brines.** Davies activities hold to an ionic strength of about 0.5; the soda brines past that are
  given by mass balance, and their minerals' solubilities are not computed.
- **The minerals' kinetics.** Calcite at three times saturation stands for calcifiers; the clays' threshold, the
  magnesite route and the opal of silica-building life are design choices, and dolomite and zeolites are not followed.
- **The water.** Its volume is held fixed: escape to space, the water clays bind and the spreading of filled seas over
  low ground are outside the model.
- **Places.** One mixed sea stands for the common-level waters; the sea phosphorus study shows the crater waters and
  Orientale taking most of what the rivers bring, and the South Pole–Aitken Sea little.
- **The land's own long future.** Erosion holds its rate while the seas stay open; how the land's relief and its
  weathering wear down over hundreds of millions of years without uplift belongs to geography.

## Files

| File | Holds |
|---|---|
| [run.py](run.py) | The box model: deliveries, weathering, the lakes and seas, minerals, sediment, phosphorus and the works |
| [test_run.py](test_run.py) | Its binding, regeneration, balances and the states it should hold |
| [sources.json](sources.json) | The published values and how far each was read |
| [results/sea_chemistry.json](results/sea_chemistry.json) | The product |

```
python -m research.studies.sea_chemistry.run
python -m pytest research/studies/sea_chemistry
```
