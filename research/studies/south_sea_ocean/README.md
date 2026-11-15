# The south sea as an ocean: its salt, its water and its rivers

15 November 2026, on `domain/hydrosphere`. The author asked to test the viability and the cost and benefit of making the south sea, the South Pole–Aitken Sea, an Earth-like ocean, the more so if it connects to brackish and fresh waterways, which would also give eels and other fish that hatch in salt water their sea. The request followed a conversation on why the Moon's waters turn to hard water and soda lakes, kept on the research drive (`conversation-exports/chatgpt-salt-seas-2026-11-15/`). The fresh seas and the leaning toward aquariums for saltwater life are the register's rows of [15 October](../../decisions.md#salt-water-and-the-fresh-seas-15-october-2026); the [seawater lakes](../seawater_lakes/README.md) of 6 November are the smaller alternative this study is weighed against.

This first stage answers what the project's products already settle: how much salt the sea holds, whether its water balance lets it keep that salt, which rivers reach it, what its floor takes, and how salt changes its mixing. [inflows.py](inflows.py) repeats geography's river routing under three rules (3.5 minutes, 5 GB) and writes [results/inflows.json](results/inflows.json) (schema `terluna.research.south-sea-inflows/1`). [run.py](run.py) reads it with the design years' climatology, the atlas and geography's crust uptake, in seconds, and writes [results/south_sea_ocean.json](results/south_sea_ocean.json) (schema `terluna.research.south-sea-ocean/1`). Where the salt comes from is worked out on `domain/resources` (`resources/sea_salt.md` there) and summarised below. No ocean or climate model ran; the sea's circulation, its brackish reaches, its spray and its living value are named under [what is open](#what-is-open) with what each needs.

## What it finds

**The sea holds 174,000 Gt of salt at seawater's strength, 55% of it chloride** *(computed)*. The sea covers 2.32 million km², 6.1% of the Moon, from 22° S to 87° S today and from 14° S to 82° S on the turned Moon. It is 2.1 km deep on average and 6.7 km at its deepest, and its 4.84 million km³ are 53% of the six large seas' water. At seawater's make-up of 35.165 g/kg and 29.5 °C:

| | Chloride | Sodium | Sulphate | Magnesium | Calcium | Potassium | Bicarbonate | Minor ions | All salt |
|---|---|---|---|---|---|---|---|---|---|
| g per kg of sea water | 19.35 | 10.78 | 2.71 | 1.28 | 0.41 | 0.40 | 0.14 | 0.08 | 35.165 |
| In the sea, Gt | 95,650 | 53,300 | 13,410 | 6,350 | 2,040 | 1,970 | 690 | 400 | 173,800 |

A sea at half that strength, 17.6 g/kg like the Black Sea's surface *(recalled)*, holds 85,800 Gt, and one at a fifth, 7.0 g/kg like the central Baltic *(recalled)*, 34,100 Gt. Magnesium, calcium and sulphur are plentiful in lunar rock; chloride is the ion that has to be brought, and sodium comes with it (below).

**The sea loses 0.71 m of water a year to the air beyond its rain** *(computed on products)*. Over the model's water cells under the sea, which hold 96% of its area, the design years give 0.51 mm of rain a day against 2.45 mm of evaporation, at a mean surface of 30.2 °C. The loss is 1,640 km³ a year, 51,900 m³/s, and the whole sea's volume in 2,950 years. The climate lane carries the model's rain as uncertain by a factor of 0.75 to 2 against its cloud-resolving boxes; applied to the rain on the sea, the loss runs from 1,200 to 1,750 km³ a year. The sea is driest where it is coldest:

| Latitude on the turned Moon | 19° S | 30° S | 42° S | 53° S | 64° S | 75° S | 80° S |
|---|---|---|---|---|---|---|---|
| Rain, mm a day | 1.23 | 1.26 | 0.51 | 0.17 | 0.07 | 0.06 | 0.08 |
| Evaporation, mm a day | 2.86 | 2.94 | 2.81 | 2.42 | 1.66 | 1.21 | 0.97 |
| Surface, °C | 30.7 | 31.4 | 30.8 | 30.0 | 29.1 | 27.7 | 27.1 |

A sea that loses more water than its rain returns keeps its level only on what flows in, and keeps every gram of salt that flows in with it. The seawater lakes study found the easiest salt water to hold is a terminal lake in a dry belt, where every flow of water and salt runs inward and only evaporation leaves. The south sea is that setting at more than a thousand times the area of the largest lake that study sited.

**As the atlas draws the waters, almost no river reaches the sea; once the craters on the way overflow, rivers bring 61% of its loss, and by way of Moscoviense 176%** *(computed on products)*. geography's drainage counts every hollow below the common sea level as sea, so a river that ends in a crater below that level never arrives. [inflows.py](inflows.py) repeats the routing on the same design years under three rules:

| Rule for the waters below the common level | Rivers into the south sea | Share of its loss (over the rain range) | Its level, left to itself |
|---|---|---|---|
| Every hollow takes its rivers in at the common level, as geography's drainage has it | 26 km³ a year | 2% (1–2%) | falls 0.69 m a year |
| The six large seas take rivers in; other hollows fill by their own balance and spill | 998 km³ a year | 61% (57–83%) | falls 0.28 m a year |
| Only the Near-side and South Pole–Aitken seas take rivers in; the four smaller seas fill and spill too | 2,884 km³ a year | 176% (165–240%) | rises 0.54 m a year |

The first rule cannot hold for the hollows that matter. The four that take 5,500–32,900 m³/s are craters of 3–280 km² at the common level, and the water between that level and their rims, 800–2,000 km³, is two to six years of their rivers. Full, they stand 2.5–4.5 km above the common level and send their rivers on. Under the second rule three rivers then enter the sea on its north and east shores: 16,000 m³/s at 15° S, 174° W on the turned Moon, about the Mississippi's flow *(recalled)*; 6,300 m³/s at 41° S, 140° W; and 4,000 m³/s at 16° S, 171° E.

The third rule adds the far side's largest river system. It drains 2.1 million km² of the far side's wet tropics into the Moscoviense Sea at 59,000 m³/s, 83 times what that sea evaporates. Moscoviense spills 3,545 m above the common level, and the 272,000 km³ between the two levels are 145 years of its surplus. Its overflow runs south and enters the south sea at the third river's mouth, 16° S, 171° E, with 63,800 m³/s in all, a third of the Amazon's flow *(recalled)*.

**So the sea's level can be held by rivers alone, by choosing how much of the Moscoviense system's water goes south** *(reading of the computed values)*. Between 1,000 and 2,900 km³ a year can reach the sea against a loss of 1,640 (1,200–1,750). The balance point is about 640 km³ a year beyond the second rule's rivers, a third of what the Moscoviense route carries. Sending it means either letting Moscoviense fill to its rim, which makes it a lake 3.5 km higher and a quarter larger, or cutting its outflow at a chosen level and place; neither is designed here. Whatever the sea receives beyond its loss raises it, and at 1.35 km above the common level it would join the nearside ocean ([geography](../../../geography/README.md)), so the inflow has to be governed and its surplus sent elsewhere.

**The fresh water enters at the north, which is where brackish water would lie** *(reading; not modelled)*. All the large mouths fall on the sea's northern half, its warm and rainier end, between 15° S and 41° S on the turned Moon, and the far south receives no river. A terminal sea fed at one end grades from fresh at the mouths to its saltiest water at the far end: the Caspian runs from about 1 g/kg off the Volga to 13 g/kg in the south, and the Mediterranean from 36 to 39 g/kg eastward *(recalled)*. The rivers therefore give the sea its fresh and brackish waterways without a fishway: estuaries at three or four mouths, a river system of 2.5 million km² above the largest, and Moscoviense and the far-side lakes at its head. How large the brackish reach is, and whether the far south turns saltier than sea water, depends on the sea's circulation and waits for the ocean models set up on this branch ([climate/ocean](../../../climate/ocean/README.md)).

**Evaporation would drive the sea's deep water, where cooling drives the fresh seas'** *(computed; the consequence is a reading)*. At 29.5 °C a gram of salt per kilogram weighs as much as 2.2 K of cooling. The sea's net evaporation salts a surface layer 100 m thick by 0.25 g/kg a year if nothing mixes it, the weight of 0.55 K of cooling each year, where the [sea oxygen study](../sea_oxygen/README.md) found the fresh sea's deep water renewed by cold spells of 0.5–2.9 K. The cold end at 80° S, 27.1 °C against 31.4 °C at 30° S, is 1.45 kg/m³ denser at equal salt, which 1.9 g/kg of extra salt at the warm end would offset. Surface water therefore grows dense over most of the sea and sinks, as in the Mediterranean and the Red Sea *(recalled)*, and the deeps should stay ventilated. Under the river plumes the opposite holds: river water is 26 kg/m³ lighter than the sea and caps it. Sea water also holds about a sixth less oxygen than fresh water at the same warmth (6.3 against 7.7 mg/L, the sea oxygen study's table). Whether the deep basins keep their oxygen is a result only once the lake model's columns run with salt.

**The crust takes salt with the water it draws from the floor, a fifth of the sea's salt within 10,000 years at the central case** *(computed on products)*. geography's stage 1 follows the water the dry crust takes from under each sea floor. Under this sea, with what it carries at seawater's strength:

| Crust | By year 1,000 | By year 10,000 | By year 100,000 | Salt a year, years 1,000–10,000 |
|---|---|---|---|---|
| Slow | 24,000 km³, 0.5% of the sea | 195,000 km³, 4% | 1.12 million km³, 23% | 0.7 Gt |
| Central | 64,000 km³, 1.3% | 1.03 million km³, 21% | 2.53 million km³, 52% | 3.8 Gt |
| Fast | 503,000 km³, 10% | 2.48 million km³, 51% | 3.47 million km³, 72% | 7.9 Gt |

The salt goes down with the water and does not come back. Holding the sea at seawater's strength therefore takes the fill and then the crust's share on top, 37,000 Gt more by year 10,000 at the central case; if the water alone is replaced, the sea freshens by the same shares. Stage 1's columns are fresh water, and what salt does in the hot rock at depth is outside them.

**The rivers' own load is small.** At the held fresh sea's 0.14 g/kg, 1,640 km³ of river water a year bring 0.23 Gt of dissolved load, a millionth of the sea's salt, most of it calcium and magnesium carbonate that calcifiers lay down. The slow gathering of alkalinity in a salted sea is the [sea chemistry](../sea_chemistry/README.md)'s: pH 9 after 0.30–0.53 million years if left alone.

**For eels the sea offers what the lakes could not, with its warmth still to check** *(reading of published and recalled values)*. European eel eggs are activated and reared at 30–40 g/kg in water with sea salt's full mix of ions, and fertilisation failed in plain sodium chloride solutions; larvae from the fourth day tolerate half-strength water (Sørensen et al. 2016; Benini et al. 2022; Okamura et al. 2009; summaries, in the lookup on the drive). A sea at half strength would therefore not serve the spawning, while its brackish north serves the later larvae. The pressure of the eels' spawning depth, 150–200 m on Earth, lies 0.9–1.2 km down at lunar gravity (the seawater lakes study), which a sea 2.1 km deep on average has and no sited lake had. The open check is temperature: Japanese eel eggs hatch at 19–28 °C and none at 31 °C (the register's row of 6 November), and the sea's surface stands at 27–31 °C, so the spawning layer has to be the cooler water at depth or in the south.

**The salt has one source, the feed the project already plans to process, and whether it suffices rests on a grade nobody has measured** *(resources' calculation on `domain/resources`, summarised)*. Lunar ground cannot give the chloride: at the 7–70 ppm inferred for soils it would take 2–19% of the Moon's mass. Near-Earth asteroids cannot either: the whole population holds the chloride of at most 3,500 km³ of sea water, 0.07% of the sea. The chloride rides in the outer-system feed for the air and the water. If the air's nitrogen came from rock like CI chondrites, at 717 ppm of chlorine, that feed would carry 7.6 times the sea's chloride (1.2 times at the grade of one disputed 2017 study); if from ice like comet 67P's, whose gas alone has been measured, 13–93% of it, with the water's own chloride adding 4%. The rock that holds the sea's chloride, 1.3×10²⁰ kg, also holds 13 times its sodium. The salt is 1.3% of the standing water's mass and one to five years of the delivery stream.

**The other seas' balances do not close as drawn either** *(computed on products; outside the question asked)*. On the same reading:

| Sea | Net evaporation | Rivers in, six-seas rule | Share of its loss | If left to itself |
|---|---|---|---|---|
| Near-side | 0.69 m, 3,660 km³ a year | 1,890 km³ a year | 52% (95% with the smaller seas spilling) | falls 0.33 m a year |
| South Pole–Aitken | 0.71 m, 1,640 km³ | 998 km³ | 61% | falls 0.28 m a year |
| Humboldtianum | 0.73 m, 209 km³ | 8 km³ | 4% | falls 0.71 m a year |
| Smythii–Marginis | 0.38 m, 108 km³ | 1,071 km³ | 10 times | fills to its rim, 767 m up, in 230 years |
| Orientale | 0.35 m, 31 km³ | 520 km³ | 17 times | fills to its rim, 2,342 m up, in 470 years |
| Moscoviense | 0.35 m, 23 km³ | 1,900 km³ | 83 times | fills to its rim, 3,545 m up, in 145 years |

The atlas's one common level is a storage geometry, as geography states. Under the design years' rain three of the six seas would rise to their rims within centuries and the other three would fall unless water is moved between them. That belongs to the atlas round; it is recorded here because the south sea's supply depends on it.

## What is open

- **The circulation, the brackish reaches and the deep oxygen.** A run of the sea with its rivers, evaporation and tide in one of the ocean models set up on this branch, and the lake model's columns with salt in the density. Neither is sized yet.
- **Spray and salt in the rain.** A salt sea sends salt inland where the fresh seas send almost none. Earth's one measured transect in hand, Australia's, gives 101 kg of chloride per hectare a year at the shore, 13 at 100 km and 2.8 at 500 km (Wilkins et al. 2022, in the lookup), and the sea's own catchment returns what falls in it. The lunar case differs in wave breaking, bubble bursting and the slow settling of particles in tall, dense air, and has to be worked out before the cost to soils, to the fresh lakes downwind and to aerial craft can be stated.
- **The shores.** The fresh seas feed green belts from the groundwater under their shores ([coast water](../coast_water/README.md)). Around a salt sea that groundwater is salt and the belts become salt shores; their length and area have not been taken off the products.
- **The living value.** What an ocean of this size adds in marine life, food and culture beside the seawater lakes and aquariums, and what the conservation programmes would do with it, is not assessed here.
- **A weaker sea.** Half strength halves the salt and loses the eels' spawning; the reefs' and the open sea's limits (corals grow at 25–42 g/kg, summaries in the lookup) set how far a cheaper sea can go.
- **The rain.** Every balance here rests on the model's rain at T21, which the climate lane carries as a range; the rivers were not rescaled with it.

## How it is set up

- **The sea** is the atlas's body 1417 at the 28% water share, on the 4 pixel/degree grid; its areas by latitude band and its place on the turned Moon follow the frame record's rotation.
- **Rain, evaporation and surface temperature** are the design years, 10–19 of the corrected run on the turned Moon (`A28R_301k_design_l15_snow_s938_turned_refit_precession`), read at the model's water cells under the sea, each atlas cell given to the nearest model cell. The model's sea is a 50 m slab of fresh water.
- **The rivers** are geography's routing of the same years' runoff over LOLA's 16 pixel/degree heights (`geography.drainage.route`), run three times with different outlet waters. The first run reproduces geography's own product, 199,973 m³/s to the waters below the common level. A basin's spill level is the level of the depressions under its cells, which is one value across each basin.
- **The salt** uses the sea chemistry's seawater (Millero et al. 2008) and the seawater lakes study's density (EOS-80).
- **The floor** multiplies the sea's area in each latitude band by stage 1's water held under a sea floor in that band.

## Assumptions and limits

- The climate, the lakes' open-water balance and the sea level stay the design years' in all three routings. A hollow that fills or a sea that spills would change the climate that made its rivers, so the rules bracket where the rivers end and predict no level.
- Salt lowers water's vapour pressure by about 2%, which the evaporation here leaves out; so are the sea's own gravity and the crust's sag, which move its shore by tens of metres.
- The times to fill leave out the growing lakes' evaporation, which is a few km³ a year for the craters and tens for the smaller seas.
- Earth's seas named for comparison are recalled values and serve as reference points.
- Nothing here is validated against a lunar observation.

## The checks

[test_run.py](test_run.py) and [test_inflows.py](test_inflows.py): seawater's ions sum to its salinity and its charges balance; the salt agrees with the sea chemistry's 174,000 Gt and the seawater lakes' 35.94 Mt/km³; the sea keeps its area on the model's grid and turns with the frame; the balance reads only the model's water cells; the floor's take grows with time and carries its salt; waters join across the grid's seam; the six seas are named by the atlas's centroids; the routing reproduces geography's discharge; each rule passes more water to the two great seas than the one before; and both products regenerate (the routing's only with `TERLUNA_SOUTH_SEA_REGENERATE=1`).

## Running it

```
ulimit -v 9000000; python -m research.studies.south_sea_ocean.inflows   # 3.5 minutes, 5 GB
python -m research.studies.south_sea_ocean.run                          # seconds
```

Each run's river mouths go to `research/runs/south_sea_ocean/`, ignored by Git.
