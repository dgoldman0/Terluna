# The atlas round

The atlas round builds the Open Moon's atlas down to its districts. It began on 16 November 2026 on its own branch, `domain/geography`, and the register's section of that day sets its frame ([the atlas round's frame](../../decisions.md#the-atlas-rounds-frame-16-november-2026)). Each step stands on the one before it, so the round works from the ground up:

1. **The surface redrawn under its water** (below): the floor and the seas' levels once the water's own gravity and the crust's sag and creep have settled, barring erosion.
2. **The drainage, the river network and the globe** on that floor, with the three upper seas filled and spilling, Moscoviense's water divided for the most diversity and a channel from Smythii–Marginis's surplus to Humboldtianum.
3. **The design climate's run on the final surface** at the Moon's settled turn, 6.15°, about 3–7 hours of machine time once the turned frame and surface are built. Whether the near-side sea may part at its pass is decided on it.
4. **The south sea's ocean runs** on the new basin and climate, the paused oxygen run among them.
5. **The places**: districts defined by significance and wholeness, after a quick review of the heritage register's written evaluations; the young-Moon landforms to keep, chosen as a set; satoyama-like landscapes, and satoumi on the salt sea's shores; the riparian share of the land against Earth's, with the falls counted on finer heights; the lakes among fresh and soda; and the three-dimensional atlas on the present design climate.

## The surface redrawn under its water

The atlas fills LOLA's floor to one level of the dry Moon's geoid, so that the water below that level covers 28% of the surface. The water changes the floor and the level both: its weight bends the crust down, the mantle creeps away beneath the load, and the water's own gravity lifts the sea surface over deep water. Geography's crustal-loading study measured these for one common level on 4 November ([geography](../../../geography/README.md#the-crusts-response-to-the-water)). The round's frame adds what that study did not have: the three upper seas, Moscoviense, Orientale and Smythii–Marginis, filled to their rims above the common level, and the crust once its response has settled.

[redrawn_surface.py](redrawn_surface.py) solves that surface. Every cell below a common level is water, the level set so that those cells cover 28% of the surface as the atlas's do. Each upper sea stands 1 cm under its rim, the lowest pass from its water to another sea's, over all the land its water joins below the rim, as the [drainage network study](../drainage_network/README.md) defines both. The lakes keep the depths that the turned Moon's corrected drainage gives them and weigh on the crust as a fixed load wherever no sea covers them. The Moon answers through load Love numbers: the crustal-loading study's elastic response at every degree to 359, and the creep that the geology lane computes for the present Moon at degrees 2–50 ([geology](../../../geology/README.md)), its change from the water's arrival to the case's time. The sea-level equation is iterated as geography's study iterates it, and after each solution the catchments and passes are found again on the deformed heights, until no rim moves by more than 5 cm. The results go to [results/redrawn_surface.json](results/redrawn_surface.json) (schema `terluna.research.atlas-round-redrawn-surface/1`) and the grids to the research drive.

The round's surface stands at the crust's settled state, barring erosion (register, 16 November): LOLA's landforms as they are, under the sag and creep the water has finished causing. Under the present Moon's central temperatures the creep never quite stops. It grows by about the same amount every tenfold of time: the degree-2 load Love number h′ gains −0.09 by 1,000 years, −0.32 by 10 million and −0.53 by a billion. The geology lane holds today's temperatures, which it takes to describe the Moon for about 100 million years; after that the cooling interior stiffens the mantle and keeps the sag already made. The design case is therefore the central profile 100 million years after the water arrives, and the other cases bound it or trace the way there:

| Case | What it stands for |
|---|---|
| Atlas | The dry Moon's floor and geoid, deformed by nothing; it reproduces the drainage network study's rims and geography's sills |
| Elastic | The first years after the water arrives |
| Central, 1,000 years | The first millennium, the horizon of the build |
| **Central, 100 million years** | **The settled state, the design** |
| Central, 1 billion years | The most the present temperatures could add |
| Stiff, 1 billion years | The cool profile at its most |
| Soft, 10 million years | The warm mantle, carried as a stated risk; it settles within about 10,000 years |

All seven ran on 16 November in 25 minutes on one thread, 1.6 GB at most. Six settled within four rounds of finding the rims again; the warm mantle's did not (below).

| | Atlas | First years | First millennium | **Settled** | Billion years, central / stiff | Warm mantle |
|---|---|---|---|---|---|---|
| Common level against the atlas's | 0 | +36.8 m | +36.5 m | **+39.2 m** | +41.2 / +36.6 m | +61.1 m |
| Water, seas and lakes | 15.23 | 16.16 | 16.18 | **16.45** | 16.62 / 16.21 | 17.88 million km³ |
| Beyond the atlas's 13.37 million km³ | +1.86 | +2.79 | +2.81 | **+3.08** | +3.25 / +2.85 | +4.51 million km³ |
| Moscoviense's rim above the common level | 3,545 | 3,424 | 3,424 | **3,449** | 3,461 / 3,426 | 3,509 m |
| Orientale's rim | 2,342 | 2,351 | 2,352 | **2,341** | 2,335 / 2,352 | 2,283–2,337 m |
| Smythii–Marginis's rim | 767 | 726 | 726 | **720** | 720 / 724 | 735 m |
| South Pole–Aitken Sea, area | 2.31 | 2.40 | 2.40 | **2.41** | 2.41 / 2.41 | 2.41 million km² |
| its volume | 4.82 | 5.62 | 5.64 | **5.80** | 5.86 / 5.66 | 6.34 million km³ |
| its deepest water | 6,795 | 7,217 | 7,227 | **7,335** | 7,386 / 7,241 | 7,757 m |
| Near-side Sea, area | 5.19 | 5.08 | 5.08 | **5.08** | 5.08 / 5.09 | 5.12 million km² |
| Tenebrarum sill, under the sea | 75 | 313 | 319 | **368** | 385 / 326 | 460 m |
| Crisium passage | 155 | 203 | 205 | **220** | 230 / 207 | 338 m |
| Fecunditatis's western channels | 109 | 105 | 105 | **108** | 111 / 105 | 108 m |
| The near side's narrowest join, under the sea | 109 | 72 | 72 | **69** | 68 / 71 | 75 m |

The water's own gravity does most of the redrawing, within years of the water's arrival. Already in the first years the common level stands 36.8 m above the atlas's, the South Pole–Aitken Sea is 422 m deeper at its deepest and its sill to the Tenebrarum sub-basin lies 313 m down instead of 75. The creep to the settled state then adds less: 2.4 m to the common level, 118 m to the south sea's deepest water, 55 m to the Tenebrarum sill and 17 m to the Crisium passage, as the floor under the south sea's deepest water falls 189 m in all (24 m within years). The first millennium differs from the first years by under a metre in the level and by 21 thousand km³ of water, so the creep that matters runs over the following hundred thousand to hundred million years.

The water moves towards the south. Holding 28% at the common level keeps the area the same, so the shore moves both ways: at the settled state the common level floods 179,000 km² of land and gives up 178,000 km² (a read-only look, `atlas-round/redrawn_surface/look/shore_moves.py` on the research drive). The South Pole–Aitken Sea gains 76,000 km² of new sea along its shores and gives up none, as its own gravity draws water to it and its floor sinks most, while the Near-side Sea gives up 107,000 km² of its shallow margins and floods 24,000; the hollows held at the common level lose 21,000 km² net. Along today's shore the sea stands a median 2.8 m higher against the land, from 56 m lower (5th percentile) to 255 m higher (95th percentile). Nearly all of this comes within years: 69,000 km² gained by the south sea and 106,000 given up by the near side in the first-years case.

Measured from the redrawn common level, Moscoviense's rim stands 96 m lower than on the dry Moon and Smythii–Marginis's 47 m lower, while Orientale's holds within a metre. Part of each fall is the common level's own rise of 39 m; the rest is the ground's sag and the water's gravity about each rim. The three seas' deepest water is 180, 140 and 36 m deeper than on the dry Moon. Filled to their rims at the settled state they cover 844,000, 298,000 and 599,000 km² and hold 1.32, 0.38 and 0.36 million km³ beyond what the atlas holds in the same places at the common level, 2.06 million km³ together. Each joins several seas at its rim, so which way its overflow runs is the drainage's to find on the redrawn floor.

The settled surface takes 16.45 million km³ of water, 3.08 million km³ beyond the atlas's 13.37 million (10.72 million in the seas, 2.65 million in the lakes), a global layer of 81 m. The register's estimate of about 3.0 million km³ for the seas' 28% with the upper seas on top holds. The near side stays one sea: its narrowest join, now the neck across western Procellarum to Cognitum, Humorum and Nubium, lies 69 m under the settled sea, against the atlas's 109 m at Fecunditatis's western channels; whether the near-side sea's own balance holds it there is the design climate's run to show.

The bounds stay close to the settled state except the warm mantle. A billion years under the present temperatures would add 0.17 million km³ and 2 m to the level, and the stiff profile stays within 0.04 million km³ of the first millennium. The warm mantle, the least likely of the three profiles, takes 4.51 million km³ beyond the atlas with the common level 61 m up and the south sea 7.76 km deep, and there Orientale's rim does not settle: filled to one of two passes, its weight sinks the ground until the other is the lower, so its rim cycles between 2,283 and 2,337 m and the case is reported with that range.

What the round's next steps take from here: the deformed heights (the settled case's U and N on the research drive, with its water depths and the upper seas' extents), the common level and the three rims for the drainage on the redrawn floor; the water to deliver, 3.08 million km³ beyond the atlas; and the gaps between the first years and the settled state for holding the surface near its settled shape sooner (below).

## What is open

- **Holding the settled shape sooner** (16 November, to consider soon). The crust's own creep cannot be hastened, but the visible surface could be built near the settled state and held there while the crust creeps beneath: the water delivered to match, each upper sea's level set by an outlet or a raised sill, the straits cut or built to their settled depths. From the first years to the settled state the common level rises 2.4 m, Moscoviense's rim 25 m, Orientale's and Smythii–Marginis's fall 9 and 6.5 m, the water grows by 0.29 million km³, the Tenebrarum sill deepens 55 m and the Crisium passage 17 m. The drainage step sizes the outlets.
- **Satoyama-like landscapes** (16 November): their criteria, what a satoyama-like mosaic needs on the Moon's land, water and climate, are worked out beside the drainage step, and their candidates mapped with it and refreshed on the final climate. The methods chosen limit methane, as the register holds methane at 1.9 ppm by biological means with technology making up the difference.
- **The warm mantle's Orientale**: a bisection on Orientale's level would find the level at which its rim and its water agree, if that case comes to matter.

## Limits

- The water arrives at once in the response, where it fills over centuries; the settled state hardly depends on that, the first millennium's state a little more.
- The creep's Love numbers are incompressible and end at degree 50, beyond which the present Moon's thick lid barely flows; the elastic part is compressible and runs to degree 359.
- The lakes keep their depths as a load. The drainage rebuilt on the new floor gives them their own levels.
- The rotational and tidal bulge re-forming about the turned axes moves the sea surface by metres and is left out, as are the water the crust's pores take up and its uneven load, which in the crustal-loading study's relaxed case moves the straits by 64–78 m and the Tenebrarum sill by about 190 m, and the lateral changes in the crust and lithosphere.
- Every cell below the common level counts toward the 28%, as the atlas counts it, so hollows away from the seas are held at the common level until the drainage gives them their own.
- Nothing here is validated against a lunar observation; the Moon carries no load of this size.

## Running it

```
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.redrawn_surface    # about 25 minutes, 1.6 GB
```

It reads LOLA's heights and the geoid through `geography.fetch_inputs`, the geology lane's response table on the research drive (`geology/interior-response/interior_response_love.npz`) and the turned Moon's corrected drainage grid (`geography/turned_corrected/products/drainage_28pct_16ppd.npz`, ignored by Git, copied from main's checkout), and writes its grids to `atlas-round/redrawn_surface/` on the research drive, with a file of the cases finished so far after each case.

## The checks

[test_redrawn_surface.py](test_redrawn_surface.py): the creep adds nothing when the water arrives and follows the response table between its times; it changes only the degrees the table holds; an upper sea stands at its rim over its own land while the rest of the water stands at the common level; and a rim is the lowest pass to any other sea, with every sea joined there named. On the product: the atlas case reproduces the drainage network study's rims to 0.1 m and geography's sills to 1 cm; every case holds 28% at its common level, the design case and the elastic one with their rims settled; and the elastic and design cases deepen the South Pole–Aitken Sea and its sill to the Tenebrarum sub-basin and need more water than the atlas.
