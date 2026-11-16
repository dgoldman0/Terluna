# The drainage network as the basins fill

16 November 2026, on `domain/hydrosphere`. The register's row of 15 November on the drainage network asks for the network that forms as the Moon's basins fill to be worked out as a study of its own, ahead of the atlas round: the lakes' levels and the rivers' courses together, from the south sea up to Moscoviense and through the other basins that fill and spill, with the globe page regenerated on its result. The atlas round, which settles the common sea level, receives the network's natural outcome and the choices it leaves.

[network.py](network.py) computes it on the project's products in about 19 minutes of machine time and writes [results/network.json](results/network.json) (schema `terluna.research.drainage-network/1`). It follows each of the six large seas of the atlas above and below the common level. For each it finds the rim it spills over, the basins its water joins on the way up, the water and the years it takes to fill, and where its overflow goes. It then follows all six together from the atlas's state for 10,000 years and traces the long profiles of the largest rivers and of every overflow river. The south sea study's routing ([research/studies/south_sea_ocean](../south_sea_ocean/README.md)) gave candidate connections under three fixed rules; this study replaces the rules with each sea's own level.

## What it finds

**The three smaller seas fill and spill whatever the rain; the two large seas' natural fate is not settled by the routing** *(computed on products; the second half is a reading of the budget below)*. Moscoviense, Orientale and Smythii–Marginis receive five to more than a hundred times what they lose to the air at the common level, so they fill to their rims and overflow in every case, within about 110 to 2,900 years. On the design rain the near-side ocean and the south sea receive in rivers 47–71% of what they lose, the range between the routing's own budget and the closed one, a margin inside the routing's error and the rain's range. In an ungoverned network they could fall by a hundred metres or more, or rise for millennia, and the atlas round has to choose their levels.

**The routing's water budget does not close, and the large seas' balances lie inside its error** *(computed)*. The climate model's own budget closes: over the whole Moon its rain and evaporation agree within 0.04%, and its land's rain less evaporation, 142,000 m³/s, equals the loss over its water. The routing counts runoff only where a land cell's rain exceeds its evaporation. On the 16 px/deg grid, with all six seas at the common level and the hollows below it as lakes, it routes 109,000 m³/s from land whose rain less evaporation is 46,000 m³/s. The cells it leaves out are land where the climate model's evaporation exceeds its rain; most of its coarse cells near water mix land and water, so their land carries the water's evaporation *(a reading; the cells were not mapped)*. With the lakes adding 96,000 m³/s of rain beyond their evaporation, the six seas receive 202,000 m³/s against a loss of 139,000 m³/s at the common level, 46% more. Raising the open water's evaporation by 14.6% everywhere closes the budget to about 2%. That case is given beside the routing's own and the rain range's ends, and it is the better guide to the two large seas *(reading)*.

**Each smaller sea's water floods neighbouring basins below its rim, and Moscoviense's floods a lowland nearly twelve times its own area** *(computed)*. A sea that rises above the common level spreads over every basin its water reaches below its lowest pass. Those basins' own hollows stand in the atlas as separate waters at the common level:

| Sea | In the atlas | Rim above the common level | Overflows to | At the rim | Water to fill | Years to fill: routing's budget, closed budget |
|---|---|---|---|---|---|---|
| Moscoviense | 64,800 km² | 3,545 m | the south sea | 837,000 km² | 1.31 million km³ | 790, 880 |
| Orientale | 84,000 km² | 2,342 m | the near-side ocean | 289,000 km² | 373,000 km³ | 800, 950 |
| Smythii–Marginis | 268,000 km² | 767 m | the near-side ocean | 568,000 km² | 443,000 km³ | 450, 540 |
| Humboldtianum | 277,000 km² | 938 m | the near-side ocean | 547,000 km² | 512,000 km³ | does not fill |
| South Pole–Aitken | 2.31 million km² | 1,455 m | the near-side ocean | 2.98 million km² | 4.11 million km³ | only with Moscoviense's overflow |
| Near-side | 5.19 million km² | 767 m | Smythii–Marginis, the same pass | 8.72 million km² | 5.82 million km³ | does not fill on the closed budget |

Moscoviense's water spreads little as it rises, to 138,000 km² at 3,250 m above the common level. Between about 3,270 and 3,320 m it joins the lowland to its south around the Freundlich–Sharonov basin *(the basin's name recalled; its centre near 18° N, 175° E on today's Moon matches the largest hollow joined)*, whose own hollows lie 2.3–2.7 km below the common level, and from 3,321 m the two fill as one water. The south sea study counted Moscoviense's own 272,000 km³ and 145 years. The route from the south sea up to Moscoviense therefore opens after about 800–900 years, through a lake of 837,000 km² standing 3.5 km above the south sea. The south sea's rim of 1,455 m lies in the 100 m step above the −0.3 km at which geography's water inventory joins it to the near-side ocean.

**The order and dates of the connections** *(computed: the six seas followed together from the atlas's state)*:

| Case | Smythii–Marginis full | Moscoviense full | Orientale full | The south sea | The near-side ocean |
|---|---|---|---|---|---|
| The routing's budget | 430 years | 770 | 780 | falls 76 m by year 500, then rises; full and overflowing at year 4,565 | rises to 741 m by year 10,000 |
| Budget closed | 515 | 870 | 940 | falls 143 m by year 1,000, then rises to 1,400 m by year 10,000 | falls 120 m within 500 years and stays there |
| Three quarters of the rain | 1,265 | 1,615 | 2,935 | falls to −817 m by year 10,000 | falls 129 m and stays there |
| Twice the rain | 110 | 225 | 175 | full at year 615 | joins Smythii–Marginis at year 435 |

The four large craters on the rivers into the south sea fill in two to six years, as the south sea study found, before any of these.

**Humboldtianum dries towards its deepest hollow** *(computed)*. Its rivers bring 249 m³/s against a loss of 6,450 m³/s at the common level, so on its own rivers it falls to its deepest hollow, 2.3 km down on the routing's budget and to its floor on the closed one. Only twice the rain fills it, over 10,600 years.

**The near-side ocean stands on its maria's floors, so a fall beyond 130 m empties it** *(computed)*. On the closed budget its rivers and the overflows of the seas above it hold it 120 m below the common level, where it keeps 89% of its area; on three quarters of the rain they hold it 129 m below. Its maria's floors lie between about 130 and 150 m below the common level: at −129 m it keeps 70% of its area, joined to its deepest point through channels the 1.9 km cells resolve, and at −150 m only 9%.

**The rivers fall in staircases of lakes with steep reaches between them** *(computed at 1.9 km cells)*. The largest river into the south sea, 16,000 m³/s at 15° S, 174° W on the turned Moon, traced 2,360 km upstream to where its flow falls to a hundredth of the mouth's, rises 7.7 km through 28 lakes. Its steep reaches total 99 km in 22 places, and its largest drop between neighbouring cells is 489 m. The south sea's next river, 6,300 m³/s, has 100 km of steep reaches in 32 places; Moscoviense's great river, 59,000 m³/s, has 112 km in 21 places through 25 lakes. Moscoviense's overflow, once full, runs 694 km from 3,545 m down to the south sea through nine hollows it fills on the way. Their surfaces step down at about 2,770, 1,810, 750 and 110 m, and its steep reaches total 48 km in 10 places, the largest drop 654 m between neighbouring cells. Orientale's overflow runs 1,700 km through 20 hollows to the near-side ocean, and Smythii–Marginis's 2,380 km through nine. The south sea's own overflow would follow the lowest passes 10,500 km through 16 hollows to the near-side ocean, a long way round. Where fish or boats could pass is therefore decided at these steep reaches. At 1.9 km cells they mark where falls and rapids lie, and a fall's own height within a cell is not resolved.

**For the salt sea, the natural network is the opposite of the design** *(reading)*. Left to its own rivers, the south sea falls by 1.1 km on the routing's budget and 1.6 km on the closed one. The Moscoviense overflow that would hold it arrives only after 800–900 years, then overfills it for millennia until it spills over its rim at 1,455 m into the near-side ocean, carrying its salt with the overflow. The salt sea's design holds its level with a governed share of the Moscoviense system's water. That share can come early only if it is carried past the Freundlich–Sharonov lowland, which otherwise has to fill with Moscoviense first, 1.3 million km³ in all.

## What the atlas round receives

- **The common level is a storage geometry, and the seas would not keep it.** Each sea's natural level differs. The near-side ocean drops 120–130 m on the closed budget and the dry end of the rain, to the edge of its maria's floors. The south sea falls, then rises over millennia once Moscoviense overflows. Humboldtianum dries, and the three smaller seas rise to their rims.
- **Choices it leaves** *(not designed here)*: hold each upstream sea below its rim by an outlet cut at a chosen level, or let it fill; govern how much of the Moscoviense system goes south and by which way; feed Humboldtianum or let it shrink; hold the near-side ocean by moving water to it.
- **The hollows the atlas holds at the common level** become a staircase of lakes as the seas above them fill, among them the Freundlich–Sharonov lowland at 3.5 km above the south sea.

## How it is computed

- **The heights** are LOLA's at 16 pixels per degree above the GRAIL geoid, as geography's drainage takes them, against the atlas's common level of −1,654 m. The water's own gravity and the crust's sag under it, which the atlas round folds in, are left out, as in every routing so far.
- **The seas** are the six largest waters below the common level, named as the south sea study names them; each one's deepest cell is its seed.
- **Where a sea's water reaches.** Every cell drains to its steepest neighbour, every pit gathers a catchment, and geography's routines find the lowest pass between each pair of catchments that touch: 275,888 catchments with no outlet at all. Geography's flood, run from one sea's catchment alone, gives every catchment the lowest chain of passes that links it to that sea. A cell joins the sea's water at the higher of its own height and that chain's highest pass. Sorting the cells by that level gives the sea's area, volume and loss to the air at every level, with every basin its water floods on the way up.
- **The rim and the overflow.** A sea's rim is its lowest pass to another sea's water. Where several seas lie beyond one pass, the overflow goes to the one its spill path reaches. That path is found by following the water downhill through a routing in which this sea is a hollow like any other and the others take the rivers in.
- **The rivers.** Geography's routing with one change, in [routing.py](routing.py): each sea takes the rivers in at its own level, where geography's takes every outlet as one water at one height. Every hollow and lake that is not a sea fills by its own balance and spills, as geography's drainage does. The change lives in this study because six products elsewhere bind the hash of `geography/drainage.py`; the atlas round, which rebuilds geography's drainage and those products, is the place to fold it in. The rivers a sea receives at the common level come from a routing with all six there. Those it receives at its rim come from one with it at its rim and the others at the common level.
- **Filling** runs from the common level to the rim in 400 steps, each step's volume over the rivers less the loss to the air at that level. A sea whose loss meets its rivers first stops there.
- **The chain and the history.** The seas are taken from the highest rim down. A full sea's surplus at its rim goes to the sea its overflow reaches, and two seas full at one pass join. The history follows all six together from the atlas's state in steps of five years. A sea's own rivers rise from those at the common level to those at its rim as it rises, and an overflow arrives five years after it leaves.
- **The long profiles.** The three largest rivers into each sea are traced upstream from the mouth along the larger tributary at each confluence, until the flow falls under a hundredth of the mouth's, the surface standing at each lake's level. Each overflow river is traced from the rim down its spill path to the sea below, its surface at the level each hollow on the way fills to. A reach counts as steep where the surface drops more than 20 m a kilometre over one 1.9 km cell.
- **The climate** is the design years' (`A28R_301k_design_l15_snow_s938_turned_refit_precession`, years 10–19). Land runoff is rain less evaporation where positive, and open water evaporates at the climate model's zonal mean over its sea cells, as in geography's drainage. The rain range scales the rain on land and water alike, with evaporation held; the closed budget scales the open water's evaporation alone.

## Assumptions and limits

- The climate stays the design years' as the atlas's waters made it. A sea that shrinks or grows changes the rain and evaporation that feed it, and that response is outside the study.
- The routing's budget is open by 46% at the common level. The closed case raises the open water's evaporation uniformly, a closure chosen here to balance it; which waters and latitudes carry the extra evaporation is beyond these products.
- The routing is a steady annual balance. Groundwater, infiltration losses and travel times are left out, as in geography's drainage, and so is the crust's uptake under new lakes, which geography's stage 1 computes for the seas.
- The water each state needs is not checked against the delivered inventory. The near-side ocean rising 741 m on the routing's budget would take nearly the 5.8 million km³ its rim holds above the common level, which the atlas's inventory may not hold.
- The history's rivers are interpolated between the common level and the rim, the seas above a full sea reach it instantly, and two seas that join are not followed further. Where two seas fill to one pass, as the near-side ocean and Smythii–Marginis do on twice the rain, the chain lists the first it reaches as overflowing into the other; the history has them join.
- Steep reaches are found at 1.9 km cells: they mark where falls and rapids lie, and a fall's form within a cell is left unresolved.

## The checks

- The study's routing routes exactly as geography's does when every sea stands at one level, on a slope with a hollow that spills and one too dry to; with two seas at different levels each sets the spill of its own hollows and every litre reaches a sea ([test_network.py](test_network.py)). The first version changed `geography/drainage.py` itself and reproduced the south sea study's three rules identically, but six products elsewhere bind that file's hash, so the change moved here.
- The study's own tests: a cell joins a water at its height or the pass on the way; a sea's table gives area, loss and volume and the level its loss meets a supply; a basin of one area fills in its volume over its surplus; a long profile measures a river's fall, steep reaches and lakes; an upper sea fills and overflows into the one below at the time its volume and surplus give; the product's budget, chain and history hold together.

## Running it

```
ulimit -v 14000000; nice python -m research.studies.drainage_network.network     # about 19 minutes
```

Its routing reads the climatology and geography's products from the repository.

## What is open

- **The globe page.** The register asks for the globe to be regenerated on the network. No single settled state comes out of it: the upstream seas fill over centuries, and the two large seas' levels are the atlas round's choice. Which state the globe should draw waits on the author.
- **The water inventory.** Whether the delivered water holds each state's seas: the upstream seas alone take 2.1 million km³ to fill.
- **The climate's response** to seas that shrink or grow, which a design climate run on a new surface would give.
- **The passages.** The steep reaches on the rivers between the south sea and the lakes above it, read at finer heights, decide where eels or boats would need passes built.
