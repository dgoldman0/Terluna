# The seas' phosphorus budget on the design surface

22 November 2026, on `domain/ecology`. The register holds the open water of the lakes and seas nutrient-poor, with productive habitats where nutrients gather ([the waters' trophic design](../../decisions.md#the-design-climate-and-the-waters-19-october-2026)), and the living system is meant to take up what enters. No study had said how much enters each sea, how much the living system must take, or what the supply does to the deep water's oxygen. [The phosphorus study](../../../resources/sea_phosphorus.md) of 11 October followed the land's phosphorus through the lakes of the atlas's drainage; the oxygen studies take a nutrient state as given. This study closes the budget on the design surface: [loads.py](loads.py) routes the land's phosphorus through the design surface's lakes to each of its seven seas ([results/loads.json](results/loads.json)), and [run.py](run.py) sets that supply, with the dust and people's losses, against what each sea can take ([results/sea_nutrients.json](results/sea_nutrients.json), schema `terluna.research.sea-nutrients/1`). It counts phosphorus alone; nitrogen is taken to follow it, as nitrogen fixers make up a shortfall under an air of nitrogen *(recalled; not budgeted)*.

## What it finds

**The lakes keep the land's phosphorus on the design surface as they did on the atlas** *(computed)*. The land loses 70–187 kt of phosphorus a year in solution and 1.4–4.0 Mt on particles (the phosphorus study). Of the land's runoff, 90% reaches a sea and 10% ends in lakes that pass nothing on; only 3.5% reaches a sea without crossing a lake. With each lake keeping phosphorus as Earth's lakes do, 5.6% of the dissolved loss reaches the seas and 3.6–4.3% of the particles (6.5% and 3.6–4.2% on the atlas): 4.0–10.7 kt a year in solution and 51–169 kt on particles, which settle at the river mouths. The three upper seas pass on 3% of what they receive and the main Near-side Sea 1%, so the seas feed one another almost nothing.

**The rivers bring each sea a thirtieth or less of what the oxygen studies' nutrient-poor state would need** *(computed)*. Per m² of sea a year:

| Sea | Rivers, in solution | Rivers' particles, at the mouths | Rivers and dust | Needed to sustain 34 g C/m² a year |
|---|---|---|---|---|
| Near-side Sea, main part (Uamoana) | 0.4–1.2 mg | 5–18 mg | 1.1–29 mg | 150 mg |
| Near-side Sea, southern part (Bultengis) | 0.4–0.9 | 3–12 | 1.1–29 | 117 |
| South Pole–Aitken, the salt sea (H'Yamo Melaha) | 0.04–0.10, and 0.02–0.2 freed from particles | 0.2–0.7 | 0.8–28 | 168 |
| Humboldtianum (Yuldustengis) | 0.02–0.05 | under 0.02 | 0.7–28 | 159 |
| Moscoviense (Hanakosa) | 1.1–3.0 | 17–52 | 1.8–31 | 160 |
| Orientale (Bahari Masariki) | 2.0–5.4 | 29–94 | 2.7–33 | 155 |
| Smythii–Marginis (Dagat na Kaseh) | 0.6–1.5 | 9–28 | 1.3–30 | 157 |

The dust is the phosphorus study's inherited range, 0.7–28 mg per m² a year, the same over every sea; Earth's oceans receive 0.9–10 from the air. It is the largest input to open water and the least known: over the seas it is 6–260 kt a year against the rivers' 4–11 kt in solution. The salt sea and Humboldtianum receive next to nothing from rivers, since their water comes through Moscoviense and long chains of lakes.

**In a layered sea without an outlet, the supply sets the export and the oxygen spent below** *(the budget's rule; reasoning on the column)*. Phosphorus entering the lit layer is taken up, sinks as matter and is set free where the matter decays. Above about 250 m the mixing returns it within a year or so; below, over centuries or longer. The lit layer therefore strips itself, within decades, until what its sinking matter loses below 250 m each year equals what enters. From then on the export is the supply over the share of the export that passes 250 m, 14–20% in these seas with the Moon's slower sinking, and the water below 250 m spends 106 mol of oxygen for each mol of phosphorus supplied, 3.4 mmol for each mg, whatever the export. So:

| Sea | Export on the rivers alone | On rivers and dust | Oxygen spent below 250 m, rivers and dust | With the canopy colony, the export is multiplied by |
|---|---|---|---|---|
| Near-side Sea, main part | 0.1–0.3 g C/m² a year | 0.3–6.7 | 0.004–0.10 mol/m² a year | 20 |
| Near-side Sea, southern part | 0.1–0.3 | 0.3–8.4 | 0.004–0.10 | 85 |
| The salt sea | 0.01–0.06 | 0.2–5.8 | 0.003–0.10 | 9 |
| Humboldtianum | under 0.01 | 0.2–6.1 | 0.003–0.10 | 12 |
| Moscoviense | 0.2–0.7 | 0.4–6.7 | 0.006–0.11 | 10 |
| Orientale | 0.5–1.2 | 0.6–7.4 | 0.009–0.11 | 11 |
| Smythii–Marginis | 0.1–0.3 | 0.3–6.5 | 0.004–0.10 | 11 |

Earth's poorest ocean gyres export 11–12 g of carbon per m² a year. On their rivers these seas export a thousandth to a tenth of that, and at the dust's top a half to three quarters. The oxygen studies' nutrient-poor state, 34 g, with 0.57 mol of oxygen demand passing 250 m, is five times richer than the supply sustains at the dust's top, and hundreds to thousands of times richer than the rivers alone sustain. The seas are held nutrient-poor by the lakes and by their own layering, before any living system takes anything; the open question is the opposite one, how poor they may be allowed to be.

**On its rivers every sea keeps its deep oxygen for more than three thousand years unaided, and at the dust's top for about a thousand** *(computed: the largest export, and the supply behind it, with which the unaided column keeps 2 mg/L at every level)*:

| Supply a sea can take, mg P/m² a year | For 300 years | For 1,000 years | For 3,000 years |
|---|---|---|---|
| Near-side Sea, main part | 82 | 35 | 17 |
| Near-side Sea, southern part | 49 | 18 | 7 |
| The salt sea | 96 | 45 | 24 |
| Humboldtianum, Moscoviense, Orientale, Smythii–Marginis | 101–109 | 47–51 | 24–26 |

The rivers' 0.02–5 mg lie under every limit. Rivers and dust at the dust's top, 28–33 mg, pass the thousand-year limit in the southern Near-side Sea, the shallowest, and lie under it elsewhere; they pass every three-thousand-year limit. The 1,000-year limits correspond to an export of 5–11 g of carbon. This is the budget's target for the living system: what enters the open water of each sea, by every route together, stays under these lines, or the sea needs the canopy colony's trunk and capsules.

**Three inputs decide, and none of them is the rivers** *(computed, on the ranges named)*:

- **Dust.** At 0.7 mg per m² a year it is nothing; at 28 it uses most of the thousand-year room. It falls on the whole surface, so no habitat at a mouth can take it. The open water takes it, and it is what the open sea's life would mostly live on.
- **People.** Beside rivers and dust the seas have room for 90–350 kt a year in all for a thousand-year hold, 11–44 g for each of 8 billion people. The food design loses 44–180 g a person a year in every form. Sent through the rivers and lakes, 5.6% of the dissolved part arrives, which fits. Sent straight into a sea from its shore, each gram a person adds 0.87 mg per m² at 8 billion, and a shore city's whole loss does not fit. People's phosphorus is the one input large enough to break the budget, and it enters at points, where a living system can take it.
- **The water the seas are filled with.** The phosphorus in the top 250 m at the start is stripped into the water below within decades and takes its oxygen once. Each mg/m³ in the delivered water takes 0.7–0.9 mol of oxygen per m². At 10 mg/m³ that is a tenth of the oxygen above 2 mg/L that the water at 250–900 m holds (a quarter in the southern Near-side Sea); at Earth's ocean mean, 70 mg/m³, it is 60–83% of it, and more than all of it in the southern Near-side Sea. The seas must be filled with water poor in phosphorus, under about 10 mg/m³, a requirement no study had set.

**What the living system has to take, and where** *(computed loads against the seas' own shallow floor; the uptake rates are recalled)*. The rivers' dissolved phosphorus, spread over each sea's floor shallower than 20 m, is 0.01–0.6 g per m² of that floor a year, and all of the rivers' phosphorus with the particles 0.1–11 g, highest in Moscoviense and Orientale, whose shallows are small (5,000 and 2,700 km²). Marshes and meadows on Earth lay away roughly half a gram to a gram per m² a year *(recalled)*, so the mouths' habitats can take the dissolved part on a tenth or more of the shallows, and the particles stay in the deltas, where the phosphorus study already places people's dredging. The budget assigns the living system three jobs:

1. **At the mouths and shores, take people's phosphorus before it reaches open water.** This is where taking up what enters is needed and sized: outfalls and farm drains led through reed beds, shellfish beds and meadows, with their harvest carried back to the land.
2. **In the open sea, recycle.** The supply is fixed by dust and rivers; how much life it carries depends on how many times each atom of phosphorus is used before it sinks below 250 m. The [canopy colony](../canopy_colony/README.md)'s feeders and trunk return what would be lost, and multiply the export 9–20 times in the deep seas for the same supply. With the colony the deep seas on rivers and dust would export 1–7 g of carbon at the dust's low end and 50 g or more at its top, and the oxygen spent below stays what the supply sets, with the larger turnover in the water the trunk flushes.
3. **Bury.** In the end the sea loses phosphorus only to its floor, with the calcifiers' carbonate and the mud. Over tens of thousands of years the water below gathers what the supply brings and begins to send it back up, and burial is what holds the sea where it is.

**What this changes in the canopy colony study.** Its central case, 34 g of carbon, stands for a sea fertilised to the lakes' nutrient-poor class, which the supply reaches only if people or the delivered water bring it there. On rivers and dust the colony is not needed for oxygen for one to three thousand years; it is needed for life, as the recycler that makes a very poor sea productive, and its trunk then carries the oxygen its own recycling costs. Its sizes fall with the export nearly in proportion.

## How it is set up

- **The land's loss** is the phosphorus study's: 0.063–0.178 g per m² of dry land a year, 4.5–4.8% of it in solution, riding the land's runoff at one concentration.
- **The routing** is the design routing's (the atlas round's product of 18–21 November). What a dry cell adds to the flow through it is its runoff; the six cells where an upper sea's water enters are left out of it. A lake is the cells at one level; its residence is its volume on the redrawn floor over the flow leaving it. The design surface's lakes number 31,375; half the water leaving them leaves lakes that hold less than two months of their outflow, and a tenth leaves lakes that hold 39 years or more.
- **The lakes' keeping** follows the phosphorus study: a lake passes on 1/(1 + √τ) of the dissolved phosphorus it receives, τ its residence in years (Vollenweider 1976), and 1 − exp(−β I/C) of the particles (Siyam's form of Brune's curves), with settling slowed sixfold by lunar gravity in the central and high cases. A lake that passes no water on keeps everything.
- **The seas that pass water on** are the three upper seas and the main Near-side Sea, with the flows of the design routing's record; each passes on what a water of its residence would, 890–9,400 years.
- **The column** is the canopy colony study's for each sea: its own floor, the salt sea's mixing, the Moon's slower sinking. The share of the export passing 250 m is the flux there over the sea's open water; the supply a sea can take is found by halving on the export until the unaided column just keeps 2 mg/L everywhere through the horizon.
- **Dust** is the phosphorus study's inherited range, all of it counted as usable; people's losses are that study's design target.

## Assumptions and limits

- **Earth's lake relations on lunar lakes.** The keeping of dissolved phosphorus is an empirical relation for Earth's lakes at a steady flow; lake beds that turn anoxic give phosphorus back (the lake studies), which this routing does not carry. The lakes' volumes are the routing's levels over the redrawn floor.
- **One concentration in all runoff.** Wet and dry lands, soils and farms differ; people's losses are not placed on the map.
- **The rule's horizon.** No return from below 250 m is counted. That holds for centuries to a few millennia at these supplies; over tens of thousands of years the deep water gathers enough to send phosphorus back, and the export then rises towards what burial allows. The 250 m line is the south sea run's; the fresh seas borrow it, and a fresh sea whose cold end overturns it returns its phosphorus as well as its oxygen.
- **Dust is a range inherited from an early scenario**, counted as wholly usable where only part of dust's phosphorus dissolves *(recalled)*. It decides the upper half of every table here.
- **Not counted:** phosphorus carried by flying and swimming animals between land and water (the phosphorus study's carriers), the fall of aerial life onto the seas, harvests, and nitrogen.
- **The oxygen limits are the canopy colony study's column's**, with its limits: one circulation run's mixing over three model years, fitted, on the atlas's basin.
- Nothing here is validated against an observation.

## The checks

[test_run.py](test_run.py): water is carried down a made-up drainage in order and a lake keeps its share; a lake is the cells at one level and its water leaves where a receiver lies outside it; the supply sets the export and the oxygen spent below, by hand; the colony's shares add up to what passes 250 m; the safe export is the edge of holding; salt water frees part of the particles and fresh water none. On the products: the lakes keep most of the land's phosphorus and the runoff is all accounted for; the inputs are bound by hash; the rivers bring far less than the nutrient-poor state takes and every sea holds for three thousand years on its rivers; and the budget regenerates. The loads regenerate with `TERLUNA_SEA_NUTRIENTS_REGENERATE=1`, which reads the research drive.

## Running it

```sh
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.loads    # about a minute
OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.run                           # about 15 seconds
```
