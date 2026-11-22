# The seas' phosphorus budget on the design surface

22 November 2026, on `domain/ecology`. The register holds the open water of the lakes and seas nutrient-poor, with productive habitats where nutrients gather ([the waters' trophic design](../../decisions.md#the-design-climate-and-the-waters-19-october-2026)), and the living system is meant to take up what enters. No study had said how much enters each sea, how much the living system must take, or what the supply does to the deep water's oxygen. [The phosphorus study](../../../resources/sea_phosphorus.md) of 11 October followed the land's phosphorus through the lakes of the atlas's drainage; the oxygen studies take a nutrient state as given. This study closes the budget on the design surface: [loads.py](loads.py) routes the land's phosphorus through the design surface's lakes to each of its seven seas ([results/loads.json](results/loads.json)), [dust.py](dust.py) reckons the phosphorus that wind-blown dust brings from the Open Moon's own ground and winds ([results/dust.json](results/dust.json)), and [run.py](run.py) sets that supply, with people's losses, against what each sea can take ([results/sea_nutrients.json](results/sea_nutrients.json), schema `terluna.research.sea-nutrients/1`). Phosphorus is the budget's measure; nitrogen is set beside it within bounds.

## What it finds

**The lakes keep the land's phosphorus on the design surface as they did on the atlas** *(computed)*. The land loses 70–187 kt of phosphorus a year in solution and 1.4–4.0 Mt on particles (the phosphorus study). Of the land's runoff, 90% reaches a sea and 10% ends in lakes that pass nothing on; only 3.5% reaches a sea without crossing a lake. With each lake keeping phosphorus as Earth's lakes do, 5.6% of the dissolved loss reaches the seas and 3.6–4.3% of the particles (6.5% and 3.6–4.2% on the atlas): 4.0–10.7 kt a year in solution and 51–169 kt on particles, which settle at the river mouths. The three upper seas pass on 3% of what they receive and the main Near-side Sea 1%, so the seas feed one another almost nothing.

**The rivers bring each sea a thirtieth or less of what the oxygen studies' nutrient-poor state would need** *(computed)*. Per m² of sea a year:

| Sea | Rivers, in solution | Rivers' particles, at the mouths | Rivers and dust | Needed to sustain 34 g C/m² a year |
|---|---|---|---|---|
| Near-side Sea, main part (Uamoana) | 0.4–1.2 mg | 5–18 mg | 0.5–3.8 mg | 150 mg |
| Near-side Sea, southern part (Bultengis) | 0.4–0.9 | 3–12 | 0.4–3.6 | 117 |
| South Pole–Aitken, the salt sea (H'Yamo Melaha) | 0.04–0.10, and 0.02–0.2 freed from particles | 0.2–0.7 | 0.1–3.0 | 168 |
| Humboldtianum (Yuldustengis) | 0.02–0.05 | under 0.02 | 0.05–2.7 | 159 |
| Moscoviense (Hanakosa) | 1.1–3.0 | 17–52 | 1.2–5.7 | 160 |
| Orientale (Bahari Masariki) | 2.0–5.4 | 29–94 | 2.1–8.1 | 155 |
| Smythii–Marginis (Dagat na Kaseh) | 0.6–1.5 | 9–28 | 0.6–4.2 | 157 |

The dust is this study's reckoning, 0.03–2.7 mg per m² a year, the same over every sea (below). The salt sea and Humboldtianum receive next to nothing from rivers, since their water comes through Moscoviense and long chains of lakes, so the dust is most of what they get.

**Dust brings 0.03–0.4 mg of phosphorus per m² a year, a twentieth or less of the range the project had inherited** *(a first reckoning from the sand study and the aerosol study's dust note)*. The phosphorus study carried 0.7–28 mg, Earth's deposition, as a sensitivity input. The Open Moon's own figure follows from two things worked out since. [The atlas round's sand study](../atlas_round/sand.md) finds that plants and the kept crust hold the dry land, that 26,000 km² of ground is bare and moving (10,000–48,000 across its cases), and that a wind over such ground carries 35 m³ of sand past each metre of width a year (5 on the climate run's means alone, 140 at the gusts' high end). The aerosol study's dust note gives the dust that saltating sand raises as a share of the sand carried, with Earth's share scaled by gravity as its central choice. Bare moving ground then raises 0.1–0.9 kg of dust per m² a year, the rate of Earth's active deserts, 2–24 Mt in all, and its phosphorus, spread evenly over the Moon, is 0.03–0.36 mg per m² a year. With the gusts, the dust share and the moving ground all at their high ends it is 2.7 mg, and all at their low ends 0.002. Two things could raise it. Ground laid bare by people adds 0.1–1.4 mg for each percent of the arid and semi-arid land left bare and loose, so one percent of that land tilled dry or trampled outweighs all the natural dust. And the note's high scheme, Earth's dust share unscaled by gravity, multiplies every figure by six. Particles of living origin are not counted.

**In a layered sea without an outlet, the supply sets the export and the oxygen spent below** *(the budget's rule; reasoning on the column)*. Phosphorus entering the lit layer is taken up, sinks as matter and is set free where the matter decays. Above about 250 m the mixing returns it within a year or so; below, over centuries or longer. The lit layer therefore strips itself, within decades, until what its sinking matter loses below 250 m each year equals what enters. From then on the export is the supply over the share of the export that passes 250 m, 14–20% in these seas with the Moon's slower sinking, and the water below 250 m spends 106 mol of oxygen for each mol of phosphorus supplied, 3.4 mmol for each mg, whatever the export. So:

| Sea | Export on the rivers alone | On rivers and dust | Oxygen spent below 250 m, rivers and dust | With the canopy colony, the export is multiplied by |
|---|---|---|---|---|
| Near-side Sea, main part | 0.1–0.3 g C/m² a year | 0.1–0.9 | 0.002–0.013 mol/m² a year | 20 |
| Near-side Sea, southern part | 0.1–0.3 | 0.1–1.1 | 0.001–0.012 | 85 |
| The salt sea | 0.01–0.06 | 0.02–0.6 | 0.0003–0.010 | 9 |
| Humboldtianum | under 0.01 | 0.01–0.6 | 0.0002–0.009 | 12 |
| Moscoviense | 0.2–0.7 | 0.3–1.2 | 0.004–0.020 | 10 |
| Orientale | 0.5–1.2 | 0.5–1.8 | 0.007–0.028 | 11 |
| Smythii–Marginis | 0.1–0.3 | 0.1–0.9 | 0.002–0.014 | 11 |

Earth's poorest ocean gyres export 11–12 g of carbon per m² a year. On rivers and dust these seas export a thousandth to a sixth of that. The oxygen studies' nutrient-poor state, 34 g, with 0.57 mol of oxygen demand passing 250 m, is twenty to thousands of times richer than the supply sustains. The seas are held nutrient-poor by the lakes and by their own layering, before any living system takes anything; the open question is the opposite one, how poor they may be allowed to be.

**On rivers and dust every sea keeps its deep oxygen for more than three thousand years unaided** *(computed: the largest export, and the supply behind it, with which the unaided column keeps 2 mg/L at every level)*:

| Supply a sea can take, mg P/m² a year | For 300 years | For 1,000 years | For 3,000 years |
|---|---|---|---|
| Near-side Sea, main part | 82 | 35 | 17 |
| Near-side Sea, southern part | 49 | 18 | 7 |
| The salt sea | 96 | 45 | 24 |
| Humboldtianum, Moscoviense, Orientale, Smythii–Marginis | 101–109 | 47–51 | 24–26 |

Rivers and dust, 0.05–8 mg, lie under every limit, the tightest being the 7 mg that the southern Near-side Sea, the shallowest, can take for three thousand years against its 0.4–3.6. The 1,000-year limits correspond to an export of 5–11 g of carbon. On the inherited dust's top, 28 mg, the seas would have held for about a thousand years. This is the budget's target for the living system: what enters the open water of each sea, by every route together, stays under these lines, or the sea needs the canopy colony's trunk and capsules.

**What decides is what people add and what the water starts with** *(computed, on the ranges named)*:

- **People.** Beside rivers and dust the seas have room for 320–360 kt a year in all for a thousand-year hold, 41–44 g for each of 8 billion people. The food design loses 44–180 g a person a year in every form. Sent through the rivers and lakes, 5.6% of the dissolved part arrives, which fits many times over. Sent straight into the seas, each gram a person adds 0.87 mg per m² at 8 billion: the design loss's low end just fits and its high end is four times too much, before any city's loss is gathered into the one sea it stands on. People's phosphorus is the one input large enough to break the budget, and it enters at points, where a living system can take it; the register now keeps it from entering a sea directly.
- **People's ground.** Dust from land laid bare is the second way people reach the open water, and no shore habitat can take it: each percent of the dry land left bare and loose adds 0.1–1.4 mg per m² of sea a year, as much as the rivers bring.
- **The water the seas are filled with.** The phosphorus in the top 250 m at the start is stripped into the water below within decades and takes its oxygen once. Each mg/m³ in the delivered water takes 0.7–0.9 mol of oxygen per m². At 10 mg/m³ that is a tenth of the oxygen above 2 mg/L that the water at 250–900 m holds (a quarter in the southern Near-side Sea); at Earth's ocean mean, 70 mg/m³, it is 60–83% of it, and more than all of it in the southern Near-side Sea. The seas must be filled with water poor in phosphorus, under about 10 mg/m³, a requirement no study had set.

**What the living system has to take, and where** *(computed loads against the seas' own shallow floor; the uptake rates are recalled)*. The rivers' dissolved phosphorus, spread over each sea's floor shallower than 20 m, is 0.01–0.6 g per m² of that floor a year, and all of the rivers' phosphorus with the particles 0.1–11 g, highest in Moscoviense and Orientale, whose shallows are small (5,000 and 2,700 km²). Marshes and meadows on Earth lay away roughly half a gram to a gram per m² a year *(recalled)*, so the mouths' habitats can take the dissolved part on a tenth or more of the shallows, and the particles stay in the deltas, where the phosphorus study already places people's dredging. The budget assigns the living system three jobs:

1. **At the mouths and shores, take people's phosphorus before it reaches open water.** This is where taking up what enters is needed and sized: outfalls and farm drains led through reed beds, shellfish beds and meadows, with their harvest carried back to the land.
2. **In the open sea, recycle.** The supply is fixed by rivers and dust; how much life it carries depends on how many times each atom of phosphorus is used before it sinks below 250 m. The [canopy colony](../canopy_colony/README.md)'s feeders and trunk return what would be lost, and multiply the export 9–20 times in the deep seas for the same supply. With the colony the deep seas on rivers and dust would export 0.1–20 g of carbon, up to the range of Earth's poorer open ocean, and the oxygen spent below stays what the supply sets, with the larger turnover in the water the trunk flushes.
3. **Bury.** In the end the sea loses phosphorus only to its floor, with the calcifiers' carbonate and the mud. Over tens of thousands of years the water below gathers what the supply brings and begins to send it back up, and burial is what holds the sea where it is.

**What this changes in the canopy colony study.** Its central case, 34 g of carbon, stands for a sea fertilised to the lakes' nutrient-poor class, which the supply reaches only if people or the delivered water bring it there. On rivers and dust the colony is not needed for oxygen for three thousand years or more; it is needed for life, as the recycler that makes a very poor sea productive, and its trunk then carries the oxygen its own recycling costs. Its sizes fall with the export nearly in proportion.

**Nitrogen arrives in excess or is made up by a small fixing in the sea** *(bounds; the lake relation for nitrogen is not carried)*. Life takes up 7.2 mg of nitrogen with each mg of phosphorus, so the supply above calls for 0.3–59 mg of nitrogen per m² a year across the seas. Undisturbed land on Earth loses 0.1–0.5 g of nitrogen per m² a year to its rivers (Lewis and colleagues 1999 find 5.1 kg per hectare for tropical watersheds, falling with the runoff). How much of that the lakes keep is not known here: if they kept all the nitrogen of every river that crosses one, the rivers would bring the seas 0.01–45 mg per m² a year at the low loss, and if they kept none, 20–4,300 at the high. Lightning adds 0.07–4.6 mg at the lightning study's site, by its two ways of counting. In the worst pairing, the least arriving with the most phosphorus, the sea's own nitrogen fixers would have to supply 14–21 mg per m² a year, 1–1.5 mmol, which Earth's open ocean is reported to exceed many times over *(recalled, not confirmed)*. So phosphorus is the budget's measure, on one condition the design must meet: that nitrogen fixers live in the seas. Where the deep water keeps its oxygen, the loss of nitrogen to the air from oxygen-free water does not arise.

## How it is set up

- **The land's loss** is the phosphorus study's: 0.063–0.178 g per m² of dry land a year, 4.5–4.8% of it in solution, riding the land's runoff at one concentration.
- **The routing** is the design routing's (the atlas round's product of 18–21 November). What a dry cell adds to the flow through it is its runoff; the six cells where an upper sea's water enters are left out of it. A lake is the cells at one level; its residence is its volume on the redrawn floor over the flow leaving it. The design surface's lakes number 31,375; half the water leaving them leaves lakes that hold less than two months of their outflow, and a tenth leaves lakes that hold 39 years or more.
- **The lakes' keeping** follows the phosphorus study: a lake passes on 1/(1 + √τ) of the dissolved phosphorus it receives, τ its residence in years (Vollenweider 1976), and 1 − exp(−β I/C) of the particles (Siyam's form of Brune's curves), with settling slowed sixfold by lunar gravity in the central and high cases. A lake that passes no water on keeps everything.
- **The seas that pass water on** are the three upper seas and the main Near-side Sea, with the flows of the design routing's record; each passes on what a water of its residence would, 890–9,400 years.
- **The column** is the canopy colony study's for each sea: its own floor, the salt sea's mixing, the Moon's slower sinking. The share of the export passing 250 m is the flux there over the sea's open water; the supply a sea can take is found by halving on the export until the unaided column just keeps 2 mg/L everywhere through the horizon.
- **Dust** is alpha times the drift potential times the sand's bulk density over the ground the sand study finds bare and moving, with alpha 0.165 of Earth's 10⁻⁵–10⁻⁴ per metre, the soils' phosphorus grade of 0.048–0.057%, and an even fall over the Moon. The budget takes it from the central case's low end to every term at its high end. People's losses are the phosphorus study's design target.
- **Nitrogen** is 16 mol for each mol of phosphorus taken up; the rivers' part is bounded by the runoff that crossed no lake and by all the runoff; lightning's is the lightning study's product.

## Assumptions and limits

- **Earth's lake relations on lunar lakes.** The keeping of dissolved phosphorus is an empirical relation for Earth's lakes at a steady flow; lake beds that turn anoxic give phosphorus back (the lake studies), which this routing does not carry. The lakes' volumes are the routing's levels over the redrawn floor.
- **One concentration in all runoff.** Wet and dry lands, soils and farms differ; people's losses are not placed on the map.
- **The rule's horizon.** No return from below 250 m is counted. That holds for centuries to a few millennia at these supplies; over tens of thousands of years the deep water gathers enough to send phosphorus back, and the export then rises towards what burial allows. The 250 m line is the south sea run's; the fresh seas borrow it, and a fresh sea whose cold end overturns it returns its phosphorus as well as its oxygen.
- **Dust is a first reckoning.** Its three factors are each known within a factor of five to thirty: the gusts inside a model cell, the share of saltating sand that becomes dust at lunar gravity, which no experiment has measured, and the ground left bare. It is counted as wholly usable where only part of dust's phosphorus dissolves *(recalled)*, and as falling evenly, where much falls near its source on land. Pollen, spores, smoke and the fall of aerial life are not counted; on Earth such particles carry about a tenth of the air's phosphorus *(recalled)*, and above the Open Moon's forests and its aerial life they may be the larger part.
- **Not counted:** phosphorus carried by flying and swimming animals between land and water (the phosphorus study's carriers), the fall of aerial life onto the seas, and harvests. Nitrogen is bounded, not budgeted: the land's loss is Earth's, and no rule for the nitrogen a lake keeps is carried.
- **The oxygen limits are the canopy colony study's column's**, with its limits: one circulation run's mixing over three model years, fitted, on the atlas's basin.
- Nothing here is validated against an observation.

## The checks

[test_run.py](test_run.py): water is carried down a made-up drainage in order and a lake keeps its share; a lake is the cells at one level and its water leaves where a receiver lies outside it; the supply sets the export and the oxygen spent below, by hand; the colony's shares add up to what passes 250 m; the safe export is the edge of holding; salt water frees part of the particles and fresh water none; the dust raised is the share of the sand carried, and its phosphorus is spread over the Moon. On the products: the lakes keep most of the land's phosphorus and the runoff is all accounted for; the inputs are bound by hash; the rivers bring far less than the nutrient-poor state takes and every sea holds for three thousand years on its rivers; and the budget regenerates. The loads regenerate with `TERLUNA_SEA_NUTRIENTS_REGENERATE=1`, which reads the research drive.

## Running it

```sh
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.loads    # about a minute
python -m research.studies.sea_nutrients.dust                                                             # seconds
OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.run                           # about 15 seconds
```
