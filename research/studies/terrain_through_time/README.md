# The terrain through time

13 October 2026. The author asked how much of the young Moon's terrain, its small pools, lakes and landforms, is planned
kept and how much wears away, and what could reasonably rebuild it. For the built land and the gardens they proposed
blocks, buildings and terraces made from the returned sediment, and asked what would serve the wild. This study works
out how long the lakes, hollows, small craters and relief last without upkeep, what keeping or renewing them costs,
and where travertine and dam builders could keep pools going on their own. The set of landforms to keep is the author's
choice in the coming atlas round ([register](../../decisions.md#the-terrain-through-time-13-october-2026)).

Two runners write two products:

- `python -m research.studies.terrain_through_time.routing` reruns geography's drainage routing, checks it against
  [geography/results/drainage.json](../../../geography/results/drainage.json) and writes
  [results/terrain_routing.json](results/terrain_routing.json) (schema `terluna.research.terrain-routing/1`) in about
  a minute and a half with 5 GB. Its grids, each lake's and hollow's fill times, the deltas at every mouth, the
  channels and the local relief, go to `research/runs/resources/terrain/terrain_16ppd.npz` (ignored by Git) for the
  atlas round.
- `python -m research.studies.terrain_through_time.run` reads that summary with the sea chemistry and the conservation
  gates and writes [results/terrain_through_time.json](results/terrain_through_time.json) (schema
  `terluna.research.terrain-through-time/1`) in seconds.

[sources.json](sources.json) records each published value and how far it was read.

## What it finds

- **Only conservation's heritage sites are kept so far.** Its register's places and a vacuum enclosure a few hundred
  metres across, holding a young crater, are the planned terrain; its zones and 31 science targets are sampled before
  the water comes. The register's direction that machines maintain the wild names no landforms. Everything else
  follows erosion, at the conservation study's 10–100 m of rock a million years.
- **Left alone, the small waters go first and the large lakes last.** The 1.9 km grid holds 28,480 rain-fed lakes over
  3.9 million km². At the central rate an eighth of them are filled within 1,000 years, 43% within 10,000 and 79%
  within 100,000. Most of these are the small ones, so 98% of the lake area is still open water at 10,000 years, 86%
  at 100,000 and 40% at a million. Lakes under 10 km² last a median 4,600 years (1,300–28,000 across the cases), the 50
  largest 1.5 million (0.5–6.9 million), and 99% of their capacity is filled within 2–25 million years.
- **The pools below the grid fill fastest.** In the saturation equilibrium of the Moon's small craters the young land
  holds 106–600 craters of 10 m and more in every km², and 1–6 of 100 m and more. Soil creep and what the pools grow
  and settle halve a 10 m crater's depth in 57–750 years, a 30 m crater's in 460–5,800, a 100 m crater's in
  3,700–44,000 and a 1 km crater's in 84,000–870,000. The airless Moon takes 44 million years to halve a 100 m crater;
  rain, creep and life work a thousand to ten thousand times faster.
- **The relief wears down over millions of years.** The land lowers 1.2–12 mm a century, fastest in the wet tropics.
  Its local relief, about 1.1 km within 10 km, halves in 8 million years within 15° of the equator and 34 million
  poleward of 60° at the central rate (4–170 million across the cases).
- **Deltas move the shorelines first.** Where rivers meet standing water the rock builds out as deltas, filling the
  nearest water to the surface. At the central rate delta plains cover 54,000 km² of the lakes within 1,000 years and
  256,000 km² within 10,000, 6.5% of their area (1.3–18% across the cases), against the 87,000 km² the lakes taken as
  whole boxes lose by then. The common-level waters, the seas with the small waters at their level, gain 16,000 and
  104,000 km² of new land. The largest deltas grow up to 38 km² in the first thousand years and 130–280 km² in ten
  thousand, reaching 9–13 km from their mouths; some fill their whole lake or crater water, as at Einthoven and
  Schlüter. The big rivers' deltas start small, since lakes upstream catch their rock until they fill.
- **The seas fill in 0.9–15 million years.** A fifth of the land's rock reaches them past every lake, and their own
  deposits add more; the [sea chemistry](../sea_chemistry/README.md)'s box gives 0.8–16 million.
- **Keeping all of it is one loop and costs a few gigawatts.** At the start 38% of the land's rock settles in the lakes,
  41% in the 119,740 dry hollows of the drier middle latitudes and poles, and 21% in the seas. Dredging the lakes'
  deltas and carrying the rock back uphill keeps the lakes open and the relief standing. At the central 3.55 Gt a year
  that is 0.19 GW of lift to where the rock came from, or 0.31 GW to each catchment's highest third; 33 billion
  tonne-km of haulage, 0.11 GW at rail freight's energy; and 0.9–2.7 km³ a year of dredging at 0.5–4.3 GW of fuel
  energy. Over the whole erosion range the loop takes 0.2–11 GW and moves 1.4–21% of what the world now quarries as
  sand and gravel. 1,281 lakes catch half of what the lakes catch.
- **The wild can take the rock back as a fall of dust.** Spread where it came from, the returned rock lands at
  1–14 mm a century; spread on each catchment's highest third, 4–42 mm a century. Loess fell on China's Loess Plateau
  at about 10 mm a century in the Holocene and 21 in the last glacial maximum. The rivers, slopes and lakes downstream
  then shape the land as they would under slow uplift.
- **The seas' stone builds the built land.** The seas must shed 0.23 km³ a year of calcite, magnesite, opal and salts
  to stay fresh, enough for terraces 10 m high over 23,000 km² every thousand years, beside the land's own
  0.27–2.7 km³ a year of rock.
- **Travertine builds while the soils are young.** At the catchment rate the rivers carry 1.26 meq/L of alkalinity and,
  at the air's CO2, stand at 0.76 times calcite saturation: they build tufa on 0.2% of their length, where warm ground
  or evaporating lakes harden them. Tufa needs about twice that load (2.0 times for Earth's streams' threshold, 1.6
  times for Plitvice's rule). The young soils supply it. In the carbon cycle's central case, 94% of the river length
  builds tufa at first, 72% after 100 years, 24% after 300 and 0.5% after 1,000. With the fastest soils, 99% builds for
  the first 1,000 years, 90% at 3,000 and 24% at 10,000; the slowest soils build none. Builders designed to calcify from
  twice saturation keep 43% of the rivers building at 1,000 years in the central case and 30,000 with the fastest soils.
  Once the soils reach the catchment rate, lasting barriers need springs carrying twice the rivers' load, or builders
  that calcify from undersaturated water as freshwater shellfish do. The springs carry 0.9–1.3 times the soil water's
  load and come to saturation or below it in air ([springs](../springs/README.md)), so the barriers rest on builders.
  Plitvice's barriers rise about 13.5 mm a year, 17 times faster than its lake floors fill, so a lake behind a living
  barrier deepens.
- **Dam builders keep the small waters the grid cannot see.** The grid resolves only channels carrying 1 m³/s and
  more, 266,000 km of them. Under lunar floods, which move a stream bed's material 1.8 times as easily as Earth's at
  the same discharge and slope, dams persist on 3.4–7.8% of that length, 9,000–21,000 km with 36,000–310,000 dams
  and 125–1,100 km² of ponds. The small streams below the grid are the builders' main ground. On Earth, beaver held
  13% of a 294 km² boreal landscape in ponds behind 835 dams in 1986, and North America held 60–400 million beaver.

## The lakes through time

Rain-fed lakes still open water after each span, central erosion with the earliest and latest cases in brackets.

| | 1,000 years | 10,000 | 100,000 | 1 million | 10 million |
|---|---|---|---|---|---|
| Lakes, by count | 88% (73–97) | 57% (35–82) | 21% (9.7–48) | 3.3% (0.3–15) | 0 (0–1.4) |
| Lake area | 99.7% (99.1–100) | 98% (94–99.5) | 86% (72–97) | 40% (11–79) | 0 (0–24) |

| Lake size | Lakes | Area | Median life, central (earliest–latest) |
|---|---|---|---|
| Under 10 km² | 15,346 | 57,000 km² | 4,600 years (1,300–28,000) |
| 10–100 km² | 9,976 | 308,000 km² | 34,000 years (9,700–198,000) |
| 100–1,000 km² | 2,528 | 733,000 km² | 203,000 years (60,000–1.1 million) |
| 1,000–10,000 km² | 580 | 1.65 million km² | 677,000 years (206,000–3.5 million) |
| Over 10,000 km² | 50 | 1.18 million km² | 1.5 million years (0.49–6.9 million) |

A lake counts as lost when it is filled to marsh, under 2 m of open water or half its water. Without the lakes' own
accretion the central median life is 27,000 years against 15,000; with erosion following runoff in place of rain,
which sends 84% of the rock to the lakes, it is 9,600. The dry hollows fill as well: 19% of them within 10,000 years
and half within 100,000 at the central rate.

## Deltas

Delta plain at the water's surface, central erosion with the earliest and latest cases in brackets.

| | 100 years | 1,000 | 10,000 |
|---|---|---|---|
| In the lakes, km² | 9,400 (1,000–34,000) | 54,000 (8,500–167,000) | 256,000 (50,000–694,000) |
| Share of the lake area | 0.24% | 1.4% | 6.5% (1.3–18) |
| Lake area the boxes lose, central | 830 km² | 11,000 km² | 87,000 km² |
| At the common-level waters, km² | 2,200 (210–9,600) | 16,000 (1,900–63,000) | 104,000 (15,000–357,000) |
| Mouths with over 1 km² of delta | 2,200 | 16,000 | 76,000 |

The largest after 10,000 years, central case:

| Where | Into | Rock at the start | Depth at the mouth | Delta at 1,000 years | At 10,000 | Reach |
|---|---|---|---|---|---|---|
| Orientale's ring (Montes Cordillera) | a lake of 1,450 km² | 2.6 Mt a year | 63 m | 32 km² | 277 km² | 13 km |
| Einthoven | a lake of 240 km² | 0.47 | 15 m | 31 | 236, the whole lake | 12 |
| Schlüter | a common-level water of 240 km² | 0.10 | 2 m | 15 | 235, nearly all | 12 |
| Neander E | a lake of 580 km² | 0.39 | 9 m | 23 | 218 | 12 |
| Mendeleev | a lake of 1,640 km² | 0.75 | 36 m | 16 | 193 | 11 |
| Mendeleev | a lake of 670 km² | 0.62 | 2 m | 14 | 183 | 11 |
| Orientale's ring | a lake of 1,120 km² | 0.02 | 36 m | 0.4 | 180 | 11 |
| Walker | a common-level water of 450 km² | 0.01 | 9 m | 0.6 | 180 | 11 |
| Mendeleev | a lake of 180 km² | 0.46 | 41 m | 23 | 179, the whole lake | 11 |
| Orientale's ring | a lake of 420 km² | 0.07 | 59 m | 1.6 | 174 | 11 |

Mouths that start with little rock grow once lakes and hollows upstream fill and pass theirs on. Shallow water near the
mouth, on the flooded maria's edges and in the lakes' shallow arms, lets a delta spread fastest. Every mouth's delta
is in the grids.

## Small crater pools

| Crater | Craters this size and larger per km² | Half the depth gone: creep | Lakes' own accretion | Together |
|---|---|---|---|---|
| 10 m | 106–600 | 60–810 years | 980–9,800 | 57–750 |
| 30 m | 12–67 | 540–7,300 years | 2,900–29,000 | 460–5,800 |
| 100 m | 1.1–6.0 | 6,000–81,000 years | 9,800–98,000 | 3,700–44,000 |
| 300 m | 0.12–0.67 | 54,000–730,000 years | 29,000–294,000 | 19,000–209,000 |
| 1 km | 0.01–0.06 | 0.6–8.1 million years | 98,000–980,000 | 84,000–870,000 |

Fresh craters (Pike 1977) diffused on a radial grid lose half their depth at κt = 0.024 D², a quarter at 0.055 D² and
nine tenths at 0.11 D²; Minton et al.'s (2019) fit to diffused craters gives half at about 0.016 D². Equilibrium
craters are already shallower than fresh ones and go sooner. Which of them hold water depends on the water table:
groundwater reaches the ground across three quarters of the mature tropics and two fifths of 15-30°, and little
poleward of 45° ([geography](../../../geography/README.md#groundwater-the-water-table-on-land)).

## The relief

| Band | Land | Local relief | Lowering, central | Relief halves, central (range) |
|---|---|---|---|---|
| 0–15° | 5.44 million km² | 940 m | 79 m a million years | 8.2 million years (4.1–41) |
| 15–30° | 5.58 million km² | 1,110 m | 77 | 10 (5.0–50) |
| 30–45° | 5.05 million km² | 1,130 m | 57 | 14 (6.9–69) |
| 45–60° | 3.69 million km² | 1,140 m | 35 | 23 (11–114) |
| 60–90° | 3.63 million km² | 1,200 m | 25 | 34 (17–169) |

Ahnert's (1970) proportion for Earth's mid-latitude basins halves relief in about 4.5 million years.

## Upkeep

| The land's rock | 0.71 Gt a year (10 m/Myr) | 3.55 (50 m/Myr) | 7.1 (100 m/Myr) | 8.3 (phosphorus' soil loss) |
|---|---|---|---|---|
| Lift back where it came from | 0.04 GW | 0.19 | 0.38 | 0.44 |
| Lift onto each catchment's highest third | 0.06 GW | 0.31 | 0.61 | 0.71 |
| Haulage over 9.3 km, at rail freight's energy | 6.6 billion t-km, 0.02 GW | 33, 0.11 | 66, 0.23 | 77, 0.26 |
| Dredging every lake | 0.18–0.54 km³, 0.10–0.85 GW | 0.90–2.7, 0.51–4.3 | 1.8–5.4, 1.0–8.5 | 2.1–6.3, 1.2–9.9 |
| Spread where it came from | 1.2 mm a century | 5.8 | 12 | 14 |
| Spread on the highest third | 3.6 mm a century | 18 | 36 | 42 |
| Share of Earth's sand and gravel | 1.4–1.8% | 7.1–8.9% | 14–18% | 17–21% |

The lift is ideal, at lunar gravity, over the mean 1.0 km from each eroding cell to where its rock first settles, or
1.7 km to its catchment's highest third. Dredging takes the lakes' share at 500–1,500 kg/m³ and the sea phosphorus
study's hopper dredging, 0.41–1.15 kg of fuel a cubic metre.

## Travertine

The share of river length whose degassed water passes Earth's streams' tufa threshold, and in brackets the share where
builders designed to work from twice saturation build, by years since the soils were first wetted.

| Soils, against the catchment rate | 0 | 100 | 300 | 1,000 | 3,000 | 10,000 | 30,000 |
|---|---|---|---|---|---|---|---|
| Fastest, 10 times the rate | 100% (100) | 100% (100) | 100% (100) | 99% (100) | 90% (99) | 24% (90) | 0.5% (43) |
| Central, 3.2 times | 94% (99) | 72% (97) | 24% (90) | 0.5% (43) | 0.3% (15) | 0.2% (1.6) | 0.2% (0.8) |
| Slowest, at the rate | 0.2% (0.5) | 0.2% (0.5) | 0.2% (0.5) | 0.2% (0.5) | 0.2% (0.5) | 0.2% (0.5) | 0.2% (0.5) |

| Load against the rivers' mean | 1 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|
| Alkalinity, meq/kg | 1.26 | 1.89 | 2.51 | 3.77 | 6.29 |
| Calcite saturation, degassed | 0.76 | 2.4 | 5.3 | 16 | 62 |

Springs build tufa where they carry about twice the rivers' mean load: at 1.5 times hardly any land's springs do, at 2
times half (the warmer half), at 3 times nearly all. How much more than the rivers the Moon's springs carry needs a
groundwater model. Stonewort meadows in hard-water lakes carry 440–690 g of calcium carbonate a square metre (Pełechaty
et al. 2013), the kind of builder a designed one would follow.

## Dam builders

Length of the 266,000 km of resolved channels where dams are built at baseflow and stand through the two-year flood.

| | Lunar floods (equal Shields number) | Earth's discharge times slope |
|---|---|---|
| Dams persist, steady flows (baseflow 0.3, flood 4 times the mean) | 20,800 km (7.8%) | 45,700 km (17%) |
| Dams persist, flashy flows (0.2 and 10 times) | 8,900 km (3.4%) | 20,800 km (7.8%) |
| Built with breaches and blowouts, steady flows | 37,700 km (14%) | 78,300 km (29%) |

At 4–15 dams a km and Earth's typical 0.35 ha a pond, the steady lunar case holds 83,000–310,000 dams and
290–1,100 km² of ponds, nearly all within 30° of the equator.

## How it is set up

- **The routing.** Geography's drainage at 16 pixels a degree (1.9 km) with the design climate's runoff, reproduced to
  the cubic metre a second, the lakes' area and volume and the count of 27,865 lake bodies.
- **The units.** Resources' 25,610 spilling lakes, merged where their water joins, and each of the 122,610 closed
  depressions, 2,870 with a closed lake and 119,740 dry: 28,480 lakes, against the drainage product's 27,865 connected
  bodies, since closed lakes inside connected water count apart. Each has its capacity (a spilling lake's water, a
  closed depression's whole hollow), its lake area, its inflow and the unit downstream, through its outlet or, once
  full, its pass.
- **Erosion.** The conservation study's rock, 0.71, 3.55 and 7.1 Gt a year, spread over the land in proportion to local
  relief, the highest less the lowest ground within about 5 km (Ahnert 1970; Montgomery & Brandon 2002), and to rain
  raised to 0.31, the BQART model's exponent on water discharge (Syvitski & Milliman 2007). Rock is 2,600 kg/m³.
- **The fill.** In steps of a twentieth of a decade to a billion years, each unit keeps what reaches it on Siyam's form
  of Brune's median curve at lunar settling, as the [sea phosphorus](../../../resources/sea_phosphorus.md) study keeps
  particles, or all of it in a closed depression, settled at 500–1,500 kg/m³, with 0.1–1 mm a year of the lakes' own
  accretion, from a deep oligotrophic lake's to Earth's small lakes' (*design guess* on Earth's rates). A full unit
  passes everything on. The earliest case pairs fast erosion with loose sediment and fast accretion, the latest slow
  with dense and slow.
- **Small craters.** Pike's fresh profile with its ejecta thinning as the cube of distance, diffused linearly; soil
  creep's diffusivity 0.003–0.04 m² a year, from Oregon's humid forest to California's semi-arid grassland (Roering et
  al. 2007); equilibrium at 0.69–3.9% of geometric saturation (Xiao & Werner 2015).
- **Travertine.** The carbon cycle's carbonate solver at the land's air temperature and the air's CO2, with the
  rivers' load from Dessert's law spread by the climate run's surface temperature, diluted by rain on the lakes and
  concentrated by their evaporation; tufa at a saturation index of 0.75 (Merz-Preiß & Riding 1999), Plitvice's three
  times saturation, and designed builders at twice (*design guess*); the young soils' factor from the carbon cycle.
- **Dam builders.** The beaver-dam capacity model's thresholds (Macfarlane et al. 2017) as discharge times slope,
  carried to the Moon at equal Shields numbers, which allows 0.41 of Earth's value; baseflow and flood ratios to the
  mean flow are *design guesses*.
- **Deltas.** At every cell where the land's rock first meets a lake or sea, and where full lakes and hollows
  upstream pass theirs on as they fill, a delta fills the nearest water of the same lake or sea to its surface, cell
  by cell outward, at each cell's depth (at least 1 m). It keeps half to all of what settles (*design guess*: fines
  spread over the floor or pass on, and lunar settling is six times slower), paired from the latest case to the
  earliest. The deltas in one lake, or in one connected sea or small water at the common level, together stop at its
  area; reach is the radius of a half-disc of the delta's area.
- **Upkeep.** Each eroding cell's lift to where its rock first settles, or to its catchment's highest third, and its
  straight distance there; lunar gravity; rail freight at 107 kJ a tonne-km (IEA 2019).

## What it leaves out

- **What the grid cannot see.** Channels under 1 m³/s, ponds, springs and craters under 1.9 km enter only through the
  crater statistics and Earth's references.
- **The deltas' shapes and the lakes' layers.** Deltas grow as fans from fixed mouths; switching channels, levees and
  lobes are not followed, and the lakes' fill times take each lake as one box. How deep lakes layer is for ecology's
  lake model (question 10).
- **Uplift and the crust.** The land has none, and the crust's response to unloading and to the water is outside it.
- **Wind and landslides.** Dust carried by wind, and the landslides that lunar gravity's taller cohesive slopes delay,
  are not followed; the dry land's storms enter through rain alone.
- **Weathering of the returned rock.** The carbon cycle keeps weathering at the catchment rate on the ground erosion
  renews; rock spread back on the heights renews it the same way.
- **The kept set.** Which landforms to keep, at what upkeep, is for the atlas round.

## Files

| File | Holds |
|---|---|
| [routing.py](routing.py) | The routing, erosion, units, fill through time, channels and the rivers' load |
| [run.py](run.py) | Lakes, small craters, seas, relief, upkeep, travertine and dam builders |
| [test_run.py](test_run.py) | Both products' bindings, balances, regeneration (the routing's on request) and physics |
| [sources.json](sources.json) | The published values and how far each was read |
| [results/terrain_routing.json](results/terrain_routing.json) | The routing summary |
| [results/terrain_through_time.json](results/terrain_through_time.json) | The product |

```
python -m research.studies.terrain_through_time.routing    # 1-2 minutes, 5 GB, with geography's inputs restored
python -m research.studies.terrain_through_time.run
python -m pytest research/studies/terrain_through_time     # TERLUNA_TERRAIN_REGENERATE=1 reruns the routing
```
