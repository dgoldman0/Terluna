# Jura Island

15 October 2026, on `research/jura-island`. The author chose Jura Island, the arc of Montes Jura round the flooded
Sinus Iridum, for a close look at its climate. The island is laid out mostly hybrid: its bay-facing slopes and shores
one of the Moon's great gardens, the high arc, the outer coasts and the islets wild, and one harbour town on the inner
bay ([register](../../decisions.md#jura-island-15-october-2026)). The climate box is CM1's `box_jura`
([the CRM's account](../../../climate/crm/README.md#jura-island)). This study holds what the island's shores and heights
add to it.

## The island at 16 px/deg

`python -m research.studies.jura_island.shores` writes [results/shores.json](results/shores.json) (schema
`terluna.research.jura-shores/1`, with its producer, the hashes of its inputs, its evidence statement and reading rule)
in about a second. Heights come from LOLA at its native 16 px/deg, about 1.9 by 1.3 km at 45° N, less the degree-200
GRAIL geoid, against the atlas's sea level of −1,654 m.

- **Size and height.** 57,824 km² (58,357 in the atlas's 4 px/deg cells) at 38.3–51.7° N and 25.2–47.0° W, with a mean
  height of 713 m. The summit stands 2,717 m above the sea at 44.41° N, 39.28° W: 2,514 m in the atlas's cells and
  2,135 m at the box's 6 km. The frame round it holds 93 islets, the largest 519 km².
- **A steep island.** Little of it lies near sea level, and its first 1–2 km inland already stand a median 112 m up,
  a tenth of them below 19 m and a tenth above 355 m.

  | Height above the sea | 0–10 m | 10–25 m | 25–50 m | 50–100 m | 100–250 m | 250–500 m | 0.5–1 km | Above 1 km |
  |---|---|---|---|---|---|---|---|---|
  | Share of the island | 0.8% | 1.3% | 2.3% | 5.2% | 17.1% | 22.1% | 23.6% | 27.6% |

  | Within | 2 km of the sea | 5 km | 10 km | 20 km |
  |---|---|---|---|---|
  | Share of the island | 10.7% | 29.6% | 54.5% | 83.1% |
  | Median height | 83 m | 206 m | 346 m | 471 m |

- **A long coast.** About 4,200 km at this resolution, four-fifths of its 5,379 km of pixel edges.

## Its shores

The fresh sea under every shore feeds its water table, and at lunar gravity plants draw on it from 15–25 m above the
sea (geography's shore screen, "Shores by the fresh seas", on domain/resources). The product's belt screen takes that
reach, and the widest belt the ground can feed on a flat shore (0.56 km at the groundwater model's central
transmissivity, 2.5 km in sandy ground), with each coastal pixel's mean height over half its width as its shore's
slope (*screen*):

- **A green rim.** The belts cover 1,084–2,627 km², 1.9–4.5% of the island. On its steep shores most are 100–200 m wide,
  and only the low flats reach the ground's supply limit; 715–1,230 km² of the island lies below the reach.
- **A monthly marsh.** The Near-side Sea's monthly tide, 3.7 m typical, floods and drains a band of every shore once a
  month, a freshwater tidal marsh.
- **Gardens beyond the rim.** The gardens on the bay-facing slopes can draw on the bay itself: fresh water needs no
  desalting, and lifting a cubic metre 250 m takes 0.11 kWh at lunar gravity before the pumps' losses.

## The heights

Rising air cools by 1.6 K a kilometre at lunar gravity against Earth's 9.8, so for the same humidity its cloud base
stands six times higher, and the 2.7 km arc lifts air as a 450 m hill would on Earth (*derived*). The arc catches cloud
only from the most humid sea air, and the island's green comes mostly from the sea beneath it. The bare box tests
whether the arc catches any.

## What follows

Open; the author has not chosen these. The recommendation is to take them in this order:

1. **The bare box runs to its end.** It gives the sky's share of the island's water (rain, cloud, dew and the shore
   nights' humidity) and the baseline a planted run is compared with.
2. **The shores' water table at tens of metres,** resources' question 16 on domain/resources, with Jura Island as its
   worked case.
3. **The island's cover:** the belts from 2, the heights' cover from the box's rain, cloud and dew, and the gardens
   placed with the author.
4. **A planted box run,** if the cover's share of the box's shore columns is large enough to matter, against the bare
   run. CM1's land categories already take a class for each column.

## Files

| File | Holds |
|---|---|
| [shores.py](shores.py) | The island's heights, its nearness to the sea, its coast and the belt screen |
| [test_shores.py](test_shores.py) | A cone island measured as a cone, and the product's binding and regeneration |
| [results/shores.json](results/shores.json) | The product |

```
python -m research.studies.jura_island.shores
python -m pytest research/studies/jura_island
```
