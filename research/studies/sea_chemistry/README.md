# The seas' make-up through time

13 October 2026. Resources' question 9 asks how the seas' composition changes through time, with the delivered water's
salts a design lever. The author asked for the model with fresh, salted and phosphorus-added water, and for it to look
past the first thousands of years to the stable or dynamic state the seas could hold over the world's billion-year
life. This box model follows the common-level waters and the rain-fed lakes that feed them from the water's delivery
to a billion years.

`python -m research.studies.sea_chemistry.run` writes [results/sea_chemistry.json](results/sea_chemistry.json)
(schema `terluna.research.sea-chemistry/1`, with its producer, the sha256 of its ten inputs, its evidence statement
and reading rule) in about fifteen seconds: 54 runs, each a delivery, a set of sinks, a fate for the sediment and one
of the carbon cycle's three weathering corners. [sources.json](sources.json) records each published value and how far
it was read.

## What it finds

- **The first thousands of years follow the carbon cycle.** Fresh water starts at pH 5.6 and takes weathering's
  carbonate: pH 7.6 by year 1,000 and 8.5 by year 10,000, at 0.04 and 0.31 g/kg. The seas reach three times calcite
  saturation in year 4,000–6,600 (central 6,600), where the [carbon cycle](../carbon_cycle/README.md) finds
  3,975–6,555 (central 6,555): only a share of the weathering reaches the waters ([springs](../springs/README.md)), and
  while the crust fills its floors draw off the waters' load with the water they take, which delivered water of the
  seas' own make-up replaces. Brackish water, Earth's seawater diluted tenfold, moves from pH 7.4 to 8.4 over the same
  time, and salted water stays near 35.1 g/kg and pH 8.0 while its calcium falls 7% as calcite settles.
- **Left alone, the seas turn alkaline.** The Moon's rock gives the rivers calcium, magnesium and a little sodium and
  potassium, each with the carbonate that balances it, and almost no chloride or sulphur. Calcite takes the calcium, and
  the rest of the alkalinity stays in seas that have no outlet. Fresh seas reach pH 9 in 30,000–60,000 years. Magnesium
  clays forming at pH 9 hold them near it for about a million years, and the sodium's soda then carries them on: 1 g/kg
  after 0.7–1.3 million years, Earth's salinity as a soda water after 27–38 million and a saturated soda brine, like
  Earth's soda lakes at about 300 g/kg, after 220–300 million. Without the clays, the magnesium keeps its carbonate in
  solution and the same path runs about fifteen times faster.
- **Delivered salt sets only the start.** Chloride has no sink but the floors' water, which delivered seawater replaces,
  so a salted sea keeps Earth's salinity, but weathering's alkalinity gathers in it all the same: pH 9 after 0.30–0.53
  million years and a tenth saltier after 0.25–0.42 million. Brackish seas join the fresh seas' path within ten million
  years.
- **Sediment fills the waters first.** The lakes catch nearly all the land's sediment and are full in 0.2–6.8 million
  years. The seas take the particles that pass and their own deposits of calcite, opal and magnesium minerals,
  0.4–1.4 km³ a year, and are full in 0.8–16 million years. A Moon without plate tectonics has nothing that renews its
  basins. The [terrain study](../terrain_through_time/README.md) routes the rock lake by lake: at first the lakes take
  38% of it, dry hollows 41% and the seas 21%, so the lakes' capacity is 99% filled in 2–25 million years and the
  seas in 0.9–15 million.
- **Every lasting state is a held one.** Four flows at the mature Moon's weathering, the 59% of the catchment rate
  that reaches the waters, keep the seas wherever the delivery sets them, for the billion years:
  - calcite, 150 Mt a year, which calcifiers build as they do on Earth;
  - magnesium, 95 Mt a year of magnesite made where the rivers enter the lakes, the carbon cycle's managed route, or
    182 Mt a year of magnesium clays at pH 9;
  - soda and sodium sulphate, 15 Mt a year from soda works on the rivers;
  - silica, 282 Mt a year of opal, which silica-building life holds at 0.1 mmol/kg.

  Held this way, fresh seas stay hard fresh water at 0.14 g/kg and pH 8.3, brackish seas at 3.65 g/kg and pH 8.3, and
  salted seas Earth's seawater at 35.1 g/kg and pH 8.0. Open water for the billion years also needs the land's
  sediment, 0.7–8.3 Gt a year, returned to it, which the conservation study's gates price at 0.11–1.1 GW of ideal
  lifting power for 0.71–7.1 Gt, and the seas' own deposits, 0.4–1.4 km³ a year, taken out.
- **Phosphorus follows the lakes' trap.** While the lakes catch it, and dredging their deltas for people's phosphorus
  keeps them catching, the mixed seas hold 0.4–23 mg/m³ fresh and 0.9–116 salted, against Earth's ocean mean of 70;
  salt water frees part of the particles' phosphorus. Once the lakes are full the seas get the land's whole loss:
  7–350 mg/m³ fresh and 20–2,600 salted. Earth's mean delivered with the water falls to the rivers' level by year
  10,000 at a 1,000-year residence and by year 100,000 at 20,000 years. An alkaline sea with its calcium spent can
  hold far more, since apatite no longer takes the phosphate (Toner & Catling 2020).
- **The carbon.** The weathering that reaches the waters takes 0.065 Gt C a year from the air in the mature Moon,
  59% of the catchment rate's 0.109. Calcite buries 0.018 Gt C a year, and the magnesite route 0.014 more; the rest
  returns to the air. While the crust fills, its floors take 129 Gt C of the waters' carbon down by year 10,000 and
  1,118 by a million years. The air's CO2 is held at 400 ppm here, so the resupply the carbon cycle sizes for its
  first 10,000 years runs on at that rate.
- **Closed lakes turn to soda, and ocean-like water needs its chlorine brought** (15 October;
  [below](#closed-lakes-soda-lakes-and-salt-seas)). A lake with no outlet in the dry belts becomes a soda lake at a pH
  near 10, reaching the ocean's salinity in about 100,000 years at 10 m deep. Seawater's main salt has to come from off
  the Moon: 36 Mt for each km³ of sea.

## The seas over a billion years

Central weathering, seas kept open. Salinity in g/kg and pH; a dash where the brine is past the constants' range.

| Delivery and sinks | 1,000 years | 10,000 | 100,000 | 1 million | 10 million | 100 million | 1 billion |
|---|---|---|---|---|---|---|---|
| Fresh, natural (magnesium clays) | 0.04, 7.6 | 0.31, 8.5 | 0.81, 9.0 | 0.96, 9.0 | 11, 9.7 | 116, – | 294, – |
| Fresh, magnesium kept | 0.04, 7.6 | 0.31, 8.5 | 1.3, 9.2 | 14, 9.8 | 163, – | 294, – | 294, – |
| Fresh, held | 0.02, 7.3 | 0.12, 8.2 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 | 0.14, 8.3 |
| Brackish, natural | 3.6, 7.8 | 3.8, 8.4 | 4.4, 9.0 | 4.6, 9.0 | 14, 9.7 | 119, – | 296, – |
| Brackish, held | 3.5, 7.6 | 3.63, 8.2 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 | 3.65, 8.3 |
| Salted, natural | 35.1, 8.0 | 35.2, 8.0 | 35.3, 8.2 | 40, 9.0 | 44, – | 148, – | 294, – |
| Salted, held | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 | 35.1, 8.0 |

| Milestone, years | Fresh, natural | Fresh, magnesium kept | Salted, natural |
|---|---|---|---|
| Calcite settles in the seas | 4,000–6,600 | 4,000–6,600 | from the start |
| pH 9 | 30,000–60,000 | 30,000–60,000 | 299,000–531,000 |
| 1 g/kg | 708,000–1,259,000 | 38,000–75,000 | |
| Earth's salinity, or a tenth above it | 27–38 million | 1.8–3.0 million | 251,000–422,000 |
| A saturated soda brine (300 g/kg) | 224–300 million | 16–24 million | 200–282 million |

Ranges run across the weathering corners. The natural seas' composition after a few million years is a soda water:
sodium and potassium with carbonate, calcium below 0.01 mmol/kg and magnesium in clays. Without the clays it is a
magnesium-soda water.

## Sediment and deposits

| | Earliest | Latest |
|---|---|---|
| The land's erosion | 8.3 Gt a year (the soil loss the phosphorus particles imply) | 0.71 Gt a year (10 m/Myr, the conservation study) |
| Dry density of what settles | 500 kg/m³ | 1,500 kg/m³ |
| The seas' deposits | 1.4 km³ a year (natural seas, 80% porous) | 0.39 km³ a year (held seas, 50% porous) |
| Lakes full | 196,000 years | 6.8 million |
| Seas full | 0.77 million years | 16 million |

Without their own deposits the seas would be full in 0.83–29 million years. When the seas fill, their water spreads
over the low ground the sediment builds: the volume stays, and the deep basins go.

## Closed lakes, soda lakes and salt seas

15 October 2026. The author asked whether the load that gathers where plants by the fresh seas draw groundwater could
make salt lakes, how such lakes would compare with Earth's ocean, and so whether ocean-like seas for saltwater life
would come cheaply. This section answers from the study's results and published values (*screen*); the shores
themselves are [geography's](../../../geography/README.md#shores-by-the-fresh-seas).

- **A closed lake gathers what its rivers bring, and the Moon's water makes soda.** Weathering gives the waters
  calcium, magnesium and 2% sodium and potassium, each with the carbonate that balances it, and almost no chloride or
  sulphur. Its alkalinity exceeds its calcium, so on the chemical divides of evaporating waters (Hardie and Eugster
  1970) a closed lake loses its calcium to calcite, its magnesium to clays and its silica to opal, and keeps its sodium
  and potassium with carbonate. Rain from fresh seas brings no chloride, where on Earth sea spray salts even inland
  rain. The Okavango's overflow ends in the Makgadikgadi Pans, whose brine is worked for both soda ash and common salt;
  the Moon's pans would give soda alone.
- **How fast.** The common-level waters and the lakes turn over once in about 1,900 years (13.8 million km³ against the
  rivers' 7,240 km³ a year). A closed lake 10 m deep that evaporates 1.5 m a year turns over in 6.7, so it runs the
  left-alone path roughly 290 times faster: 1 g/kg in 2,500–4,400 years, Earth's salinity as soda in 94,000–133,000
  and a saturated soda brine in 0.8–1.0 million. A shallower lake goes faster in proportion to its depth. The scaling
  leaves out the lake's own floor, its exchange with groundwater and the rain on it.

  | | Earth's ocean | A closed lake on the Moon, in time | The held seas |
  |---|---|---|---|
  | Salinity | 35 g/kg | past 35 g/kg towards a saturated brine near 300 | 0.14 g/kg |
  | Mainly | sodium chloride, magnesium sulphate | sodium carbonate and bicarbonate | calcium, magnesium bicarbonate |
  | pH | 8.1 | about 10 | 8.3 |

- **Soda lakes hold a rich, specialist life.** On Earth they are among the most productive waters: Arthrospira in
  Lakes Bogoria and Natron feeds the lesser flamingos, Mono Lake's brine shrimp and alkali flies feed millions of
  migrating birds, a tilapia lives in Lake Magadi at pH 10 and the pearl mullet in Lake Van's soda water (*Earth
  evidence, recalled*). Corals, shellfish, sea urchins and nearly all sea fish need seawater's chloride, calcium and
  magnesium at a pH near 8, which a soda lake lacks. The soda works can make one in years: a cubic kilometre at 30 g/kg
  holds 30 Mt of soda, two years of the works' 15 Mt a year at the mature Moon's weathering.
- **Ocean-like water needs its chlorine brought.** Seawater at 35.1 g/kg holds 36 Mt of salt per km³, 55% of it
  chloride, which lunar rock cannot give. The separate seas, at the atlas's volumes:

  | Sea | Volume | Salt for Earth's salinity |
  |---|---|---|
  | Moscoviense | 71,000 km³ | 2,600 Gt |
  | Orientale | 75,000 km³ | 2,700 Gt |
  | Humboldtianum | 331,000 km³ | 11,900 Gt |
  | Smythii–Marginis | 351,000 km³ | 12,600 Gt |
  | Near-side | 3.42 million km³ | 123,000 Gt |
  | South Pole–Aitken | 4.84 million km³ | 174,000 Gt |

  Chlorine is scarce off the Moon too. CI chondrites carry about 700 ppm (Lodders 2003), with several times as much
  sodium, and comet 67P's ice carried hydrogen chloride at about 10⁻⁴ of its water (Dhooghe et al. 2017), about 1%
  of seawater's chloride (*recalled*). At a few hundred ppm, the 10²⁰–10²¹ kg of small bodies the nitrogen plan
  processes ([register](../../decisions.md), 9 October) hold 3×10¹⁶–7×10¹⁷ kg of chlorine, two to forty times what
  the four smaller seas need (1.6×10¹⁶ kg), if it is kept. Aquariums need little: a million cubic metres of seawater
  holds 36,000 t of salt.
- **Once in, the salt stays.** Chloride has no sink but the floors' water, which delivered water of the sea's own
  make-up replaces, and the four flows hold a salted sea at 35.1 g/kg and pH 8.0. Each separate sea is its own basin,
  so salting one leaves the others fresh, and salting is simplest while a sea fills. A salt sea's own balance of
  rivers and evaporation then sets whether its rivers freshen it, and the sea phosphorus study finds Orientale and the
  crater waters taking much of the rivers' flow. Its shores turn Earth-like, with a salt wedge under them and salt in
  the spray.
- **The author's direction.** The seas stay fresh or near fresh, the author leaning to aquariums for saltwater life,
  and the soda lakes stay as their own local ecosystems
  ([register](../../decisions.md#salt-water-and-the-fresh-seas-15-october-2026)).
- **Open: the delivered water's own chlorine.** This study delivers pure water. Ice like 67P's carries about
  0.28 g/kg of chloride. Delivered as it comes, it would hold the seas near that, about 0.5–0.6 g/kg in all with the
  cations that balance it, at the line between fresh and brackish water, and the shores' soils would gather it where
  plants draw ([resources' question 17](../../../resources/questions.md)). Cleaning it from the water keeps the seas on
  this study's path and yields chlorine for salt water. The land's closed lakes stay soda lakes either way, since
  evaporation leaves the seas' salts behind and only spray carries them inland.

## How it is set up

- **The waters.** The common-level waters, 10.7 million km³, and the rain-fed lakes, 3.1 million km³, as well-mixed
  boxes at the design case's sea surface (297.0 K) and land air (294.3 K), each in CO2 equilibrium with the air at
  400 ppm and 1.216 bar. The lakes take 96.5% of the land's runoff and pass it on through 445 years of residence.
- **Weathering**, imported from the carbon cycle: Dessert's basalt law over the land's discharge, 9.1×10¹² eq a year
  at the catchment rate, with the developed soils starting at 1, √10 or 10 times it and ageing toward it on White and
  Brantley's law; erosion keeps the catchment rate renewed. The springs study's share of it reaches the waters, 36% at
  1,000 years, 42% at 10,000 and 59% in the mature Moon at the centre (each corner takes stage 1's crust of the same
  end), with its make-up: the highland soil's calcium, 65% of the alkalinity, magnesium 33% and sodium and potassium
  2%, less the calcium the groundwater lays in the rock, so 45–56% calcium arrives. Its export gives 0.87 mol of silica
  and 0.0025 of sulphate for each equivalent.
- **The floors.** Stage 1's uptake through the lakes' and seas' floors, 2,400 km³ a year at first and nothing after ten
  million years, draws off their water with its dissolved load. The water delivered to replace the seas' carries the
  delivery's own make-up, fresh, brackish or salted; the lakes' comes from their rivers.
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
