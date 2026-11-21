# Made sand seas on the Open Moon: a siting study

20 November 2026. The register's row of this date has one to three sand seas made in the hyper-arid land, fed with the sand in the sediment that the terrain's upkeep dredges from the lakes and returns uphill, and leaves the wild land's natural wind-blown sand open. This paper sites them. It asks where in the hyper-arid land a made sea can lie, how much sand it takes and from how far, what the winds make of it in ten thousand years, and what it covers. [sand_seas.py](sand_seas.py) computes it on the fields of the sand study ([sand.md](sand.md), [sand.py](sand.py)) and writes [results/sand_seas.json](results/sand_seas.json) (schema `terluna.research.atlas-round-sand-seas/1`) and the grid product `/media/projectspace/terluna-research/atlas-round/sand_seas/sand_seas_16ppd.npz`; [test_sand_seas.py](test_sand_seas.py) checks it; pictures of the seas as laid are in `/media/projectspace/terluna-research/atlas-round/sand_seas/maps/`.

Labels: *project result* is an output of a committed or branch study, named with its file; *calculation* is arithmetic done here, in sand_seas.py unless stated; *proposal* is a rule or a size chosen here for the author to weigh; *inference* is reasoning from these without a calculation. Figures without a label are calculations of sand_seas.py, recorded in results/sand_seas.json. Everything rests on the sand study's fields and carries that study's evidence state: nothing here is observed.

## The answer

**Two seas fit the hyper-arid land: 21,556 km² of star dunes on the level basalt plain round 22.8° N, 42.4° W, and 4,774 km² of transverse ridges on the floor of the walled plain round 70.3° N, 126.7° E. Together they hold 491 km³ of sand, and bringing it from the ends of drainage where the upkeep already gathers sediment takes about 0.6 GW of haulage for five thousand years.**

The hyper-arid land is three patches, 362,457 km² in all. The basalt plain is the one made for a sand sea: 81% of it is level, its winds carry three times the sand of the other patches' (a drift potential of 101 m³ per metre a year against 32 and 28) and blow from many sides, and their net drift converges on it. The far side's polar margin is four fifths of the hyper-arid land and nearly all rough highland; its one level floor of any size is that of a walled plain about 210 km across (Schwarzschild on today's maps), which holds a sea as a bowl holds water. The third patch, in the north-west, is small and sits on the hyper-arid divide.

The sand is the costly part. Dry land sheds little, so the sources near either site are thin: the ends of drainage within 500 km of the basalt plain gather 0.006 km³ of durable sand a year, and within 500 km of the walled plain 0.005 km³. The sand study's first figure, a sea of 10,000 km² twenty metres deep in about 1,200 years with 0.45 GW of haulage, took the lakes' sand of the whole Moon as lying 500 km away. With the sources where they lie, the basalt plain's 390 km³ arrive in five thousand years from ends up to 1,782 km off, 1,250 km on average, at 0.53 GW; the walled plain's 101 km³ from within 924 km at 0.07 GW. The two together draw 27% of the durable sand the loop gathers on the whole Moon through those years, and their haulage is 5.4 times the loop's own *(calculation: 177 billion tonne-kilometres a year against the terrain study's 33 billion)*. Slower is far cheaper: at ten thousand years the same two seas take 0.22 GW. The sand kept from the uphill return is 302 km³ of rock, 14 mm off the land's mean height *(calculation: 491 km³ at 1,600 kg/m³ as rock of 2,600 kg/m³ over 22.13 million km²)*.

Once laid, the seas stay. The basalt plain's net drift is 23 m³ per metre a year towards the south, 0.28 of the total, and converges; its downwind margin sheds 2.1 km³ in ten thousand years, 0.5% of the sea. The walled plain's drift is steadier, 29 m³ per metre a year towards the south, and its margin sheds 0.8 km³, 0.8%. Neither sends sand straight into water. By the sand study's rules the plain's sand stands as star dunes 1.5 km apart and 75 m high at its mean thickness, and the walled plain's as transverse ridges 1.8 km apart and 90 m high; in both the sand at hand sets the size, which the drift reaches in 5,700 to 5,800 years *(calculation: the spacing squared over forty times the net drift that moves)*.

The sand's colour follows its sources. The basalt plain's own neighbourhood sheds dark basalt sand, 84% of what lies within 250 km, and the far sources that a sea of this size needs are highland feldspar and glass: its sand is 21% basalt overall, pale on dark ground, and darker on its south-east side (43%) than on its north (3 to 4%). The walled plain's sand is highland sand throughout.

Two findings reach back to the sand study. Its rule for the rivers' sand credits each end of drainage with up to 5,000 km² of the ground upstream, so chains of lakes count the same ground more than once: the ends would gather 0.88 km³ a year where the land sheds 0.39 by the same yield. Followed to the first standing water, the wind's offer of river sand falls from 19.4 to 17.2 km³ in the central case and from 31.0 to 27.1 km³ with the sand left open, and the open case's map from 249,033 to 247,018 km² under sand; its conclusions stand, and its cases wait to be rerun on the corrected rule. And the loop's lake deltas gather 0.172 km³ of durable sand a year by this count, beside the 0.17 km³ that study drew from the terrain study's shares; the closed hollows gather 0.190 km³ more.

## A. The hyper-arid ground

**Hyper-arid land, where rain is under a twentieth of the demand, is 362,457 km², 1.64% of the land, in three patches.** The dryness is the moisture map's rain over demand, carried to the 1.9 km grid as the sand study carries it; a patch is joined ground of 3,000 km² or more.

| Patch | At, on today's Moon | Area, km² | Dryness | Drift potential, m³/m/yr | Net drift: share of the total, towards | Convergence, m per thousand years | Basalt | Level share | Heights within 10 km spread by (median) |
|---|---|---|---|---|---|---|---|---|---|
| The far side's polar margin | 69.2° N, 145.4° E; 59.0° to 74.2° N | 284,904 | 0.035 | 32.3 | 0.82, 192° | 0.006 | 0% | 3% | 314 m |
| The basalt plain | 24.2° N, 40.2° W; 19.6° to 28.7° N | 60,829 | 0.041 | 100.7 | 0.29, 194° | 0.184 | 95% | 81% | 26 m |
| The north-west | 65.8° N, 55.9° W; 63.3° to 68.1° N | 16,405 | 0.049 | 27.8 | 0.37, 232° | 0.045 | 0% | 15% | 207 m |

Ground counts as level where the heights within about 10 km spread by less than 60 m. The convergence is the sand a flux at the wind's capacity would lay down, from the net drift smoothed over 9.5 km. The basalt plain is the dry floor of Oceanus Procellarum east of the Aristarchus massif, with the near-side sea to its north and its south-west; lakes lie in its east and rivers cross its south-east. The polar margin's level ground is the floor of the walled plain at 70° N, 121° E, in two pieces of 3,588 and 2,169 km² either side of an inner crater.

## B. The sand at hand

**The terrain's upkeep gathers 0.362 km³ of durable sand a year over the whole Moon: 0.172 km³ in the deltas of lakes that spill and 0.190 km³ on the floors of closed hollows. The seas' shores take 0.025 km³ more, which the loop leaves.** The yield is the sand study's: the terrain study's lowering by latitude, spread by rain to the power 0.31, four tenths of it sand-sized and half of that lasting, 0.387 km³ a year off all the land. Here each land cell's sand goes to the first standing water or dry hollow floor on its way down the design routing, so no ground is counted twice and the ends together hold what the land sheds. Basalt is 9% of the lake deltas' sand and 6% of the closed hollows', by the rock of the ground each cell's sand came from.

For a sea the nearest ends are taken first. Haulage is priced as the terrain study prices its loop's, 107 kJ a tonne-kilometre, over great-circle distances, with sand of 1,600 kg/m³ in bulk.

| Within | The basalt plain: sand, km³ a year | basalt | The walled plain: sand, km³ a year | basalt |
|---|---|---|---|---|
| 250 km | 0.0020 | 84% | 0.0014 | 0% |
| 500 km | 0.0062 | 81% | 0.0053 | 0% |
| 1,000 km | 0.0181 | 65% | 0.0246 | 0% |
| 2,000 km | 0.1008 | 17% | 0.0917 | 1% |

## C. The seas

**Each sea's outline follows a stated rule round a seed, and its sand lies as a sheet smoother than the ground: deep in the hollows, thin over the swells.** The rules are proposals. A sea keeps 7.6 km from standing water and 5.7 km from a channel carrying 20 m³/s or more, and lies on hyper-arid ground only. Its surface follows the ground smoothed over 15 km, so a hollow under it holds up to 2.5 times the mean thickness and a swell that stands a mean thickness above its surroundings carries none. It thins to its edge over 2.7 km on its upwind side and 17 km downwind, where sand laps out over the ground.

| | The basalt plain | The walled plain |
|---|---|---|
| Rule of the outline | where the net drift converges fastest, less for swells and towards the hyper-arid land's edge | the lowest ground of the floor, up to where it would spill |
| Centre, on today's Moon | 22.77° N, 317.64° E | 70.30° N, 126.67° E |
| Area, km² | 21,556 | 4,774 |
| Mean thickness inside the taper, m | 25 | 30 |
| Sand, km³ | 390.0 | 100.8 |
| Reaches from its centre, km | 128 | 49 |
| Ground, m above the dry Moon's geoid (10th, 50th, 90th percentile) | −1,610, −1,549, −1,397 | −3,100, −2,921, −2,801 |
| Dryness; plant cover by the globe's curve | 0.038; 4.1% | 0.034; 3.6% |
| Drift potential, m³/m/yr | 87.2 | 32.6 |
| Net drift: m³/m/yr, towards, share of the total | 23.0, 188°, 0.28 | 29.0, 191°, 0.89 |
| Convergence, m per thousand years | 0.324 | −0.018 |
| Sheltered share | 0.01 | 0.09 |
| Level share of its ground | 95% | 60% |
| Relief within 5.7 km: under the sand, ridden over, stands through | 64%, 33%, 3% | 22%, 59%, 20% |
| Dunes at the mean thickness | star dunes, 1,500 m apart, 75 m high | transverse ridges, 1,800 m apart, 90 m high |
| Dunes as mapped: forms; crest spacing (10th, 50th, 90th percentile) | star 13,770 km², linear 7,786 km²; 624, 1,079, 1,820 m | transverse 4,774 km²; 772, 1,368, 1,939 m |
| Shed by the downwind margin in ten thousand years | 2.08 km³, 0.5% | 0.77 km³, 0.8% |
| Sand outside the outline (fringe and plume): area, volume, reach from the centre | 6,119 km², 22.7 km³, 142 km | 1,882 km², 5.6 km³, 124 km |
| Region and districts, by the redrawn atlas's codes, which the sand study uses | L15: L15.2, a corner of L15.1 | L06: L06.1 |
| Districts of the page's draft of 18 November | L15.2 (16,158 km² of the sea), L15.1 (5,398 km²) | L07.1 (4,162 km²), L03.4 (611 km²) |
| Nearest heritage site | Lunokhod 1 and Luna 17, 506 km | Hakuto-R Mission 1 crash site, 1,322 km |

Relief counts as under the sand where the ground within a square 5.7 km a side rises by less than the mean thickness, as ridden over up to three times that, and as standing through beyond. A dune's form follows the net drift's share of the total, by the sand study's classes; its spacing is the smaller of what the drift gathers in the years and three dune heights of sand, and its height a twentieth of its spacing. On the basalt plain the share lies on the divide between star and linear dunes, 0.3, so the map parts the sea between them.

**The haulage, for each span of delivery.** The plan drawn is five thousand years.

| Delivered in | The basalt plain (390 km³): from within; mean haul; power; basalt; share of the loop's sand | The walled plain (101 km³): the same |
|---|---|---|
| 1,000 years | more than all the Moon's ends give | 2,114 km; 1,366 km; 0.75 GW; 1%; 28% |
| 2,000 years | 2,778 km; 1,868 km; 1.98 GW; 12%; 54% | 1,413 km; 961 km; 0.26 GW; 1%; 14% |
| 5,000 years | 1,782 km; 1,250 km; 0.53 GW; 21%; 22% | 924 km; 638 km; 0.07 GW; 0%; 6% |
| 10,000 years | 1,335 km; 924 km; 0.20 GW; 38%; 11% | 689 km; 459 km; 0.03 GW; 0%; 3% |

**Other sizes on the basalt plain, at the same thickness.** The power is the haulage while it lasts; the basalt's share is the sand's.

| Area, km² | Sand, km³ | In 2,000 years | In 5,000 years | In 10,000 years |
|---|---|---|---|---|
| 5,000 | 90 | 0.24 GW, 34% basalt | 0.06 GW, 65% | 0.02 GW, 77% |
| 10,000 | 181 | 0.65 GW, 19% | 0.18 GW, 41% | 0.06 GW, 65% |
| 20,000 | 362 | 1.77 GW, 12% | 0.48 GW, 23% | 0.18 GW, 41% |
| 40,000 | 724 | 5.35 GW, 7% | 1.28 GW, 14% | 0.48 GW, 23% |

A small sea built slowly is a dark one, of the plain's own basalt sand; a large one built fast is pale. On the walled plain a sea of 10,000 km² would take 0.20 GW over five thousand years and one of 20,000 km² 0.57 GW, but its floor holds 4,774 km² before the sand would spill.

**The sand's make-up by side.** Each source's sand is set down on the side of the sea nearest it. Clockwise from north in eight directions, the basalt plain's sand is 4, 14, 24, 43, 32, 13, 10 and 3% basalt, from 0.002 km³ a year in the north to 0.018 km³ in the west; the map carries the whole sea's share at its centre and each side's share at its edge.

## D. The map

**With the seas set down, 275,799 km² lie under sand, 1.25% of the land, against 247,018 km² without them; dunes stand on 56,946 km² against 19,566 km², and 36,561 km² of cells are more than half sand where none was.** Both are the sand study's open-sand case (the crust's upkeep absent from sand, the central winds and supply, ten thousand years) with each end of drainage given the sand of its own ground. Each sea's body lies where it was put, and what its downwind margin sheds is carried as the land's own sand is. The seas send no sand into water: the loss to water is 3.9 km³ with them and without.

The product `sand_seas_16ppd.npz` holds, on the design surface's grid, the sand product's seven fields for that case (`sand_share`, `active_share`, `form`, `crest_azimuth_deg`, `crest_spacing_m`, `drift_azimuth_deg`, `sand_thickness_m`), so the globe's look reads it in place of `sand_16ppd.npz`; `made_sea` (each sea's number inside its outline); `dune_height_m`; `made_sand_part` and `made_sand_basalt` (which of a cell's sand a made sea brought, and the basalt's share of it); and `sand_share__open` and `sand_thickness_m__open`, the same case without the seas. How the globe draws a sea is in visualization/atlas/README.md.

## E. Decisions for the author and open inputs

**1. The seas' sizes and the span of their building.** It is open because the haulage grows fast with both: the sizes drawn (21,556 and 4,774 km²) and the span (five thousand years) are this study's proposals. The options are in part C's tables: smaller or slower seas for a fraction of the power, a larger plain sea for more. Recommendation: keep the sizes drawn and take the span from the date the atlas shows; five thousand years at 0.6 GW sits inside the upkeep's 0.2 to 11 GW, and ten thousand at 0.22 GW costs a third as much. It decides the outline and the sand's colour the globe draws. A rerun with other sizes takes three minutes and a recut of the globe's close view two minutes for the three packs that hold the seas.

**2. A third sea.** It is open because the register allows up to three. The options: none; the north-west patch at 65.8° N, 55.9° W, where a floor of 2,413 km² would carry linear dunes on the hyper-arid divide; a second sea in the polar margin's rough ground, dunes among craters. Recommendation: none for now. The two drawn give the three forms between them, and the north-west's dryness of 0.049 leaves its sand half held by plants.

**3. The seas' names.** They carry none here, and the page labels neither.

**4. Whether the seas enter the design surface.** It is open because the globe draws them while the surface, its routing and the districts' fact sheets do not know them. A sea 25 to 30 m deep changes no river's course at the 1.9 km grid by itself, but the districts' texts and the land's use would name it. It can wait for the regions' redraw.

**5. The sand study's cases on the corrected rule.** Its six cases stand on the capped catchment. A rerun on each end's own ground takes six minutes of machine time and changes the river sand by about an eighth; its paper's tables would follow.

**6. Dune shape at lunar gravity.** The forms, the height as a twentieth of the spacing and the spacing as three dune heights of sand are Earth's rules. A pass through the literature on dune shape in low gravity would firm up the height and the lee face.

**Inputs this study lacked**, in the order of what each would change: the terrain study rerun on the design routing, in place of its lowering by latitude; a carrier and routes for the haulage; how the sand is sorted from the dredged sediment and what that costs; the engineered plants' answer to a deep body of loose sand; the wind over open ground below the run's scales, which sets the dunes' forms.

## Limits

- Everything the sand study lists as a limit holds here: one run's winds at 170 km and 3 days, the display cover curve, Earth's dune rules.
- A sea is set down whole and stands ten thousand years. The years of its building are not followed, and its dunes reach the size given 5,700 to 5,800 years after its sand lies.
- A sea's body stays where it is put. Its own slow march downwind, the moving net drift over its thickness, is 4 to 5 km in ten thousand years on either site and is not drawn.
- The margin sheds along the net drift only. Winds from many sides fringe a sea with a sheet a few kilometres wide, which enters as the map's last smoothing over 4.7 km.
- Haulage is priced over great-circle distances at one energy per tonne-kilometre. No route, grade or carrier is designed.
- The outline rules are simple and could be drawn otherwise; a sea's volume and haulage depend on its area and thickness and hardly on its exact outline.

## Files

- [sand_seas.py](sand_seas.py): the reckoning. [test_sand_seas.py](test_sand_seas.py): seven tests (the following of drainage, an outline's rule, the taper, the haulage's arithmetic, what a margin sheds, the dunes' rule, the written product).
- [results/sand_seas.json](results/sand_seas.json): the record, with the patches, the supply, each sea's figures, haulage tables and other sizes, and the two maps' sums.
- `/media/projectspace/terluna-research/atlas-round/sand_seas/sand_seas_16ppd.npz` (127 MB): the grid product. `maps/sea_plain.png`, `maps/sea_walled.png`: each sea's sand and thickness over the relief. `explore/`: the first look at the hyper-arid land.
- [sand.py](sand.py) was restructured for this study into shared fields and one-case functions, with sand set down by hand and the rivers' sand as inputs a case may take; rebuilt in memory, its six cases reproduce the product of 19 November array for array.
