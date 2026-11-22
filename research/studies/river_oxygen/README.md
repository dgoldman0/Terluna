# The oxygen the rivers carry down into the fresh seas

22 November 2026. The [sea oxygen study](../sea_oxygen/README.md) runs the fresh seas' deep water as columns with no rivers in them, and the [river deltas study](../river_deltas/README.md) found the rivers cooler than the seas and plunging along their floors. Cool river water that sinks takes its oxygen down with it, and no study had counted that supply. This one counts it on the design surface as it stands since 18 November ([atlas round](../atlas_round/README.md#the-design-on-the-corrected-rule-on-the-climate-of-its-own-surface)), for its six fresh seas; the South Pole–Aitken Sea is salt and is treated in the [south sea study](../south_sea_ocean/README.md#the-rivers-temperature-and-the-storms-winds).

[mouths.py](mouths.py) measures the seas and their large rivers on the design routing ([results/mouths.json](results/mouths.json), schema `terluna.research.river-oxygen-mouths/1`). [waters.py](waters.py) runs the project's lake model for the lake each river leaves and for the sea beside its mouth ([results/waters.json](results/waters.json), schema `terluna.research.river-oxygen-waters/1`). [run.py](run.py) sets each river's water beside its sea and writes [results/river_oxygen.json](results/river_oxygen.json) (schema `terluna.research.river-oxygen/1`). [sources.json](sources.json) records what was read from outside the project.

Every figure is a calculation on the project's products and model runs; nothing is measured. The headline bound needs only the rivers' flows and the seas' sizes. Where the water goes rests on temperatures that differ by about a kelvin between two runs of a one-dimensional model, and is given with its range.

## What it finds

**The rivers bring each fresh sea about a thousandth of its volume a year** *(computed from the design routing and the redrawn floor)*. The six fresh seas hold 6.24 million km³ under 6.83 million km², and their rivers bring 192,000 m³/s, 6,060 km³ a year, 0.097% of their volume.

| Sea | Area, km² | Volume, km³ | Rivers, m³/s | km³ a year | Of its volume a year |
|---|---|---|---|---|---|
| Near-side sea, main part (Uamoana) | 3,839,000 | 2,784,000 | 77,600 | 2,448 | 0.088% |
| Near-side sea, southern part (Bultengis) | 973,000 | 369,000 | 6,000 | 190 | 0.051% |
| Humboldtianum (Yuldustengis) | 279,000 | 333,000 | 6,200 | 197 | 0.059% |
| Moscoviense (Hanakosa) | 844,000 | 1,440,000 | 52,800 | 1,665 | 0.116% |
| Orientale (Bahari Masariki) | 298,000 | 469,000 | 18,600 | 586 | 0.125% |
| Smythii–Marginis (Dagat na Kaseh) | 599,000 | 846,000 | 30,900 | 975 | 0.115% |

The southern part of the near-side sea also takes 9,400 m³/s of the main part's surface water through the neck. Earth's two great tropical rift lakes receive as much river water for their size or more, 0.35% of its volume a year for Lake Malawi and 0.074% for Lake Tanganyika, and both hold no oxygen below about 200 m *(sources.json; Tanganyika's depth from search summaries)*.

**If every river sank with all its oxygen, it would meet 6% of the deep water's demand in nutrient-poor seas** *(computed; the bound)*. River water saturated at its temperature under the design air holds 7.7 g of oxygen in each m³. All of it, carried below the nights' mixing, is 47 Mt of oxygen a year, against a demand of 774 Mt a year for the floor and deep water of the seas' 6.62 million km² deeper than 50 m in the sea oxygen study's oligotrophic state (117 g for each m² a year), 1,451 Mt in its mesotrophic state and 2,176 Mt in its eutrophic state: 6.0%, 3.2% and 2.1%. By sea the oligotrophic share is 1.3% in the near-side sea's southern part, 4.4% in its main part, 4.8% in Humboldtianum and 11–13% in Moscoviense, Orientale and Smythii–Marginis. If each volume of river water also took two volumes of saturated surface water down with it, the shares would be three times these, 18% over the six seas. So the rivers cannot be the deep water's main supply of oxygen under any arrangement of their temperatures; the wind's mixing, which the sea oxygen study rests on, brings the deepest tenth of each basin it ran to the surface every 0.3–3.8 years where it holds.

**Every large river reaches its sea out of a lake, and carries that lake's water** *(computed)*. The 73 rivers followed one by one carry 180,000 of the 192,000 m³/s, and the 14 of 1,000 m³/s or more carry 161,000. Each leaves a lake large enough to bring the flow through it to the lake's own temperature between 4 and 309 km above its mouth, 6 to 2,240 m above the sea, and reaches the sea within about two days. A river as deep as these, 38 m for the 30,400 m³/s of the largest, loses little of that water's temperature on so short a course: the rivers keep 29–100% of their source's departure from the local water, the great ones 55–100%.

**The river water is as warm as the sea's surface or warmer, so it stays in the upper layers** *(model result: lake-model columns on the earlier design years)*. Taken from its lake's surface down to its own channel's depth, the river's water averages 0.6–2.8 K warmer than the water at the foot of the sea's upper 40 m, by sea and weighted by flow, and lies within 0.6 K below to 1.8 K above the sea's top 10 m. Thirteen of the fourteen great rivers are warmer on average than that foot water. The lakes stand in the land's air, which the long lunar day heats more than the sea's; a lake's height cools its water by 1.04 K for each km, and the river's fall gives back 0.39 K of that as the heat of its work. In the central case 62% of the river water settles within 5 m of the surface, 34% deeper in the upper 40 m, which the nights mix, and 4.4% is denser than the whole sea column and goes below the nights' mixing; 0.2% is colder than any water the sea column ever holds.

| Sea | At the surface | In the upper layers | Below the nights' mixing | Across the cases |
|---|---|---|---|---|
| Near-side sea, main part | 49% | 42% | 8.9% | 6.2–27% |
| Near-side sea, southern part | 93% | 5% | 1.6% | 1.4–2.3% |
| Humboldtianum | 66% | 34% | none | none |
| Moscoviense | 75% | 23% | 1.2% | 1.0–6.7% |
| Orientale | 97% | 2% | 1.0% | 0.8–3.8% |
| Smythii–Marginis | 44% | 54% | 2.3% | 1.8–7.4% |
| The six seas | 62% | 34% | 4.4% | 3.2–14.5% |

One river does most of the sinking: the 17,600 m³/s that enters the near-side sea at 15° N, 288° E from a lake of 11,300 km² eight kilometres above its mouth. Its water averages 0.24 K colder than the sea's foot water, and 31% of it goes below. The rivers out of lakes shallower than 10 m, whose water cools most through the night, send up to 37% of theirs below, and they carry 1,300 m³/s between them.

**What goes below supplies about half a percent of the deep demand** *(computed from the same columns)*. In the central case 4.3 Mt of oxygen a year goes below the nights' mixing in the six seas, with the surface water the sinking rivers take in: 0.55% of the oligotrophic demand, 0.29% of the mesotrophic and 0.20% of the eutrophic. Across the cases the oligotrophic share runs from 0.24% to 0.90%. By sea it is 0.79% in the near-side sea's main part, 0.52% in Smythii–Marginis, 0.33% in Moscoviense, 0.28% in Orientale, 0.04% in the near-side sea's southern part and nothing in Humboldtianum. The supply equals the oligotrophic demand of 36,500 km² of floor, 29,500 km² of it in the near-side sea's main part, of the 6.62 million km² deeper than 50 m. The sinking river water's own decaying organic matter could use up to 0.7 Mt of it a year, at the world's rivers' mean load; water that has stood for years in a lake has spent most of that already. At the central rate the water going below would take 6,000 years to replace the deep water of the near-side sea's main part and 18,000–54,000 years in the other seas.

**The great rivers would not plunge at their mouths** *(computed with the plunge criterion at lunar gravity)*. An inflow plunges where the water is deep enough for its densimetric Froude number to fall to about a half. For the twelve great rivers whose water is at times denser than the sea's surface, that depth is 260–630 m, since the density difference is small and the Moon's gravity gives it a sixth of the force; on Earth the same river and contrast would plunge at 0.55 of that depth. These rivers therefore leave their mouths as jets and mix into the upper layers, where the classes above put them.

**The same holds for the lakes** *(a first look, computed on the research drive in `river-oxygen/first_look/lakes_flushing.py`, not a product)*. The design surface's 31,375 lakes off the seas hold 5.95 million km³ at a mean depth of 907 m. Half of that volume lies in lakes whose through-flow would take more than 890 years to replace it, and 11% in lakes it replaces within a century, so the rivers do not renew the deep lakes' water either.

## What this changes in the earlier studies

The river deltas study gave each river the land air's temperature at its mouth, 1–4 K cooler than the seas' surface, and found 69–92% of the mature rivers' water plunging through the lunar day. That holds for a river that runs long enough over land to take the air's temperature. On the design surface the large rivers leave lakes within a few days' flow of their mouths and carry lake water, which this study finds as warm as the sea's surface or warmer. Their plumes float or spread in the upper layers, which moves the deltas they build towards that study's floating-plume cases; the screen is not run again here.

The sea oxygen study's columns hold no inflows, and this study finds that leaving them out costs its deep water under 1% of its demand.

## How it is set up

- **The seas and rivers** are the design routing's (`design_routing_16ppd.npz` on the research drive, the design on the corrected lake rule and the estimated design climate), on the redrawn floor at 16 pixels a degree. A mouth is a cell off the seas whose water goes straight into a sea. The rivers kept one by one are the largest that carry 95% of a sea's river water, at most 20 in a sea; the rest, 12,100 m³/s over the six seas, take the kept rivers' shares.
- **A river's source** is the first lake up its course that holds at least 170,000 m² for each m³/s through it, the area over which a water surface exchanging 24 W/m² for each kelvin brings the flow to the lake's own temperature. Smaller lakes on the way count towards what the river loses of its source's temperature.
- **The temperatures** come from [biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py) run as the [cold waters study](../cold_waters/README.md) runs its lake, on the design years `A28R_301k_design_l15_snow_s934`, years 10–19, with two years' spin-up: for each river a column of its source lake at the lake's level, to the lake's mean depth and at most 40 m, under the land cell nearest the lake; and a column of 40 m of the sea at the sea's level under the water cell nearest the mouth. There are 125 distinct columns for the 73 rivers, each kept on the research drive (`river-oxygen/columns/`) with its temperature profile every six hours. Eight of the 73 sea columns stand on land cells, in the upper seas that the earlier run holds as land.
- **The river's temperature at its mouth** is the mean of its lake's water from the surface to the river's own depth, a travel time earlier, moved towards the lower ground's warmth by the share it has lost and warmed by its fall. The channel is the river deltas study's mature channel at the mean flow, depth 0.153 Q^0.535 m and width 3.875 Q^0.465 m at 1.69 m/s.
- **Where it goes**: the river takes in one volume of the sea's top 10 m for each of its own, and the mixture settles at the first depth of the sea column whose water is as dense. Water denser than the whole column goes below the nights' mixing. The foot of a 40 m column holds the coldest water the nights make there, which is what the deep water of these seas is made of (sea oxygen study).
- **The cases** change one choice at a time: no surface water taken in, a fifth of a volume (a mouth confined in a drowned valley) and two volumes (the Rhône in Lake Geneva); the exchange with the air at 15 and 40 W/m² for each kelvin; the channel half and twice as deep; the lake's top metre or its whole column drawn; and no surface water with each of those two draws.
- **The oxygen** is saturation at each water's temperature under the air at the sea, with the design air's oxygen share, by the lake model's own function. **The demand** is the sea oxygen study's areal rate for the floor and the deep water together, 0.32, 0.60 and 0.90 g for each m² a day in its three nutrient states, over each sea's area deeper than 50 m.

## Assumptions and limits

- **Two climates.** The flows are the design routing's on the design climate as estimated on its own surface. The temperatures are from the earlier design years, 0.7 K warmer in the mean, on today's Moon with the earlier atlas's waters, because the design climate has no record written every three days. The contrast between a lake and the sea beside it is read on that earlier climate.
- **A kelvin decides.** The lake and the sea are two columns of a one-dimensional model under one climate model's weather at T21. Whether river water slips over or under the sea's own water turns on differences of a kelvin or less between them, within what such models can be wrong by. The classes are therefore a reading with a range, and the bound of 6% does not depend on them.
- **Land air and sea air.** A lake takes its cell's air as the land has it, although a lake of thousands of km² makes its own; with the sea's air it would be nearer the sea's temperature and settle in the upper layers all the same.
- **No flow is computed.** The classes say where water of a given temperature would lie in a sea layered as the sea column is. The jet at the mouth, the underflow on the slope and its mixing there are not modelled; the mixing where an inflow plunges is a bracket from Earth's lakes and reservoirs. The slope's mixing does not weaken under the Moon's gravity (the General Lake Model's inflow equations hold gravity only in the underflow's first thickness), so water that did run down a slope would take in more of the sea's water, and what it takes in below the upper layers carries no new oxygen.
- **No sediment.** A load of 0.13 kg/m³ would make river water as dense as water 0.3 K colder, and some of Earth's rivers reach their lakes' floors by their load alone (Lake Brienz). The rivers here leave lakes, which keep their sediment, and floods are outside the mean flows; the young Moon's bare catchments would differ.
- **The demand** is Earth lakes' areal rate applied to all the floor deeper than 50 m, as the sea oxygen study applies it.
- **The small rivers**, 6% of the water, take the kept rivers' shares. Rivers out of shallow lakes send more of their water below, so the small rivers' share may be higher; if all of their water sank it would add 0.4% of the oligotrophic demand.
- **Mean flows.** The routing gives the mean year; the atlas round's year-by-year rivers are not read.

## The checks

[test_run.py](test_run.py): the channel's fits close on the flow; a river keeps less of its source's temperature the farther and shallower it runs and through small lakes; a fall of a kilometre warms water by 0.39 K on the Moon and 2.35 K on Earth; the river draws its lake's water to its own depth; colder water settles deeper and mixing only lifts it; the plunge lies (g_Earth/g_Moon)^(1/3) deeper on the Moon; the demand is the sea oxygen study's; a column looks for a cell of its own kind; each sea receives what the design routing's record says and measures its area; the kept rivers carry at least nine tenths of each sea's water and each leaves a lake of the stated size; the shares add to one and the supply stays under its bound in every case; the products bind their inputs; and both products regenerate (the calculation where the research drive holds the columns, the mouths with `TERLUNA_RIVER_OXYGEN_REGENERATE=1`).

Checks performed on 22 November: mouths.py ran in 17 s; the 125 columns ran in 31 minutes on one thread; run.py ran in 33 s; the tests passed (14, and the mouths' regeneration run once by itself). The literature lookup was bounded at 30 searches and 25 documents opened; its table, with what it could not verify, is on the research drive in `session-tools/agents/river_oxygen/values.md`. No independent review has read the study.

## Running it

```sh
ulimit -v 14000000
OPENBLAS_NUM_THREADS=1 nice -n 10 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.river_oxygen.mouths   # about 20 s
OPENBLAS_NUM_THREADS=1 nice -n 10 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.river_oxygen.waters   # about 35 minutes; skips columns already on the drive
OPENBLAS_NUM_THREADS=1 nice -n 10 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.river_oxygen.run      # about 35 s
```

The first reads the design routing and LOLA's heights; the second reads the design years' files on the research drive (it needs `netCDF4`, which the ocean environment holds) and writes the columns to `river-oxygen/columns/` there.
