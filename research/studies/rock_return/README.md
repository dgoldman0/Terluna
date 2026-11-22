# The living return of rock

22 November 2026. Rain moves the land's rock downhill into hollows, lakes and seas, and a Moon without uplift has nothing that lifts it back. The [terrain study](../terrain_through_time/README.md) priced the return for machines: the deltas dredged and the rock carried uphill. This study asks whether animals could do that work in the wild, by joining two habits that are each common on Earth, carrying earth for nests and travelling far as part of a life, and whether such a loop would hold itself. The register's row of 22 November sets the aim it is judged against: the hybrid and used land may have designed and engineered upkeep, and the wild should if possible keep a dynamically stable equilibrium on its own.

**On the design surface the land loses 3.2 Gt of rock a year (0.6–6.5 across the terrain study's range), and flying a kilogram of it the 9.4 km back and the 1.0 km up costs a carrier 18–73 kJ at lunar gravity, against 280–1,900 kJ for the same way on Earth. Carriers weighing about 5 Mt alive, eating 0.05–1.2% of the land's plant growth, would return all of it, and a loop in which they multiply where sediment gathers returns exactly what arrives and recovers from losing half its carriers in 40–550 years, while the lakes take 11,000 years to lose a hundredth of their area. Energy and food do not limit the living return. Reach and placement do: flyers and waders can get at 46–66% of the rock, divers at the delta fronts raise that to 61–97%, and the rest lies on deep floors; and the rock has to be set down within 30–400 m of the slope it left, which asks for carriers that are many and spread, or that scatter their loads.**

**Followed on the map for a hundred million years, the return as designed keeps the ground's forms for about ten million years and then turns on two things the design had left open. If worn ground can be given more than it lost, the fastest-wearing ground draws the rock to itself: the walls of craters and the steps above lakes rise into ramparts, 1% of the ground gains 0.7–0.9 km in thirty million years, single spots grow without limit, and the plains between are drained. If ground takes no more than it lost, nothing grows and the craters stay round, and country where hollows lie close keeps its map for the whole hundred million years; high country farther from any hollow than the carriers' range still wastes away. Neither reading raises the map's own summits, which stand farthest from where rock comes to rest and lose height in every run.**

Two runners write two products. [rock.py](rock.py) follows the rock on the design surface and writes [results/rock.json](results/rock.json) (schema `terluna.research.rock-return-rock/1`) in about four minutes. [run.py](run.py) sizes the carriers from it with the functions in [carriers.py](carriers.py) and writes [results/rock_return.json](results/rock_return.json) (schema `terluna.research.rock-return/1`) in seconds. [sources.json](sources.json) records each source and how far it was read. Nothing here is validated: no animal carries sediment kilometres uphill, and flight with a load at lunar gravity has not been measured.

## What it finds

**The job** *(computed: [rock.py](rock.py))*. The terrain study's central lowering, 79 m a million years at the equator falling to 25 near the poles, spread by local relief over the design surface's 22.1 million km² of dry land, takes 3.23 Gt of rock a year, 146 g from each m², a mean lowering of 56 m a million years. It carries 1.6–1.8 Mt of phosphorus, the whole of what the land loses bound to particles. Where each cell's rock first comes to rest:

| First settles | Share of the rock | Mean way back | Mean height back |
|---|---|---|---|
| On dry ground, in a hollow where a course ends | 41% | 10.2 km | 1,070 m |
| In a lake that spills | 39% | 8.3 km | 850 m |
| In a lake that passes nothing on | 14% | 10.6 km | 1,280 m |
| In a sea | 6% | 8.1 km | 870 m |
| All | 100% | 9.4 km (a tenth under 2.4, a tenth over 19) | 1,000 m (a tenth under 94, a tenth over 2,330) |

The mean way back agrees with the terrain study's on the earlier drainage, 9.3 km and 1,030 m, though the surfaces and routings differ. The seas' share is smaller here, 6% against 21%, because the design surface's lakes stand in front of them. With the lakes passing on what Brune's curves pass, the lakes that spill keep 37% of the rock, the waters and hollows that pass nothing on 57%, and the seas receive 6%. Of the 25,452 lakes that spill, 1,075 keep half of what the lakes keep.

**Whose land it leaves** *(computed, on the land-use screen)*. With 30% of the land wild, the wild land sheds 1.03 Gt a year, the hybrid land 0.63, the land people fill 0.65 and the open land beyond their need 0.92. With half wild the figures are 1.71, 0.62, 0.50 and 0.41. The wild's rock settles much as all the rock does, 43% of it on dry ground. The lakes are shared: 28% of what they keep comes off wild land (50% with half wild), and only 6% of the lake area (17%) lies in lakes fed nine tenths or more by wild land.

**Who can reach it** *(computed from the [river deltas screen](../river_deltas/README.md)'s shares and the floors' depths; the reach of a dive is a design guess)*. Rock on dry ground is all within a flyer's reach. Of a river's load into standing water, the delta plain keeps 8% where the soils are clay and 42% where they are sand, the plain and the delta front together 34% and 95%, and the rest spreads over the floor. Of the floors, 28% of the spilling lakes' area, 41% of the closed lakes' and 11% of the seas' lies shallower than 200 m, taken here as a diving bird's reach by breath.

| Share of all the rock within reach of | Deltas of clay | Deltas of sand |
|---|---|---|
| Flyers and waders: dry hollows and delta plains | 46% | 66% |
| With divers at the delta fronts | 61% | 97% |
| With divers on the floors shallower than 200 m | 73% | 98% |
| Left on the deeper floors | 27% | 2% |

Pressure is no obstacle to going deeper. Water presses a sixth as hard here: 200 m down the pressure is 4.2 atmospheres, and the pressure a seal bears at 1,500 m in Earth's sea is met here at 9 km, deeper than any of the seas. Breath is the limit, so the deep floors ask for divers of a seal's or a whale's endurance. Raising a kilogram of grains through 1,000 m of water costs 4–5 kJ, less than flying it.

**Why Earth has no such animal, and why the Moon could** *(computed with the [sky ecology](../sky_ecology/README.md) study's flight rule, glide ratio 15–25 and muscle efficiency 0.2–0.25; the loads are design guesses)*. A flyer pays its weight over its glide ratio for every metre flown and its weight for every metre climbed. From the height it has gained, a glide covers the way home for 66–81% of the rock.

| Metabolism for each kilogram delivered | Moon | Earth |
|---|---|---|
| The load alone | 9–13 kJ | 54–80 kJ |
| With the carrier, load a twentieth of its body | | 1,180–1,880 kJ |
| With the carrier, load a quarter of its body | 46–73 kJ | 280–440 kJ |
| With the carrier, load 0.6 of its body | 25–38 kJ | |
| With the carrier, load equal to its body | 18–28 kJ | |

A kilogram of delta sediment holds 0.7–1.5% organic carbon, 14–180 kJ that an animal can use if it digests 5–30% of it *(design guess)*. On Earth that is 0.01–0.6 of what the flight costs, so carrying sediment home never pays; on the Moon it is 0.2–10 times the cost, so a carrier that feeds on what the mud holds can come out ahead. The same wing bears six times the mass here, which is why loads near the carrier's own weight are entertained.

**How many carriers, and what they eat** *(computed; trips and loads are design guesses)*. A thousand loaded trips a year take 680 hours at 15 m/s, 8% of the year. At that pace and a load 0.6 of the body, the carriers for all the land's rock weigh 5.4 Mt alive whatever their size, and size trades their number and spacing against their food:

| Carrier, by weight | Number | Apart | Hauling adds to upkeep | Dry food a year | Of the land's plant growth |
|---|---|---|---|---|---|
| 60 g, a hornero's | 90 billion | 16 m | 1.4–2.2% | 379–381 Mt | 0.9–1.2% |
| 1 kg, a crow's | 5.4 billion | 64 m | 3.5–5.4% | 157–160 Mt | 0.4–0.5% |
| 10 kg, a goose's | 540 million | 200 m | 7–11% | 78–81 Mt | 0.2–0.25% |
| 1.72 t, the sky ecology study's albatross scale | 3.1 million | 2.7 km | 37–58% | 19–22 Mt | 0.05–0.07% |

The sky ecology study already counts on 1.8–4.8 million flyers of the albatross scale fed by the sky and on 2.0–10.7 million to return the land's dissolved phosphorus, so the rock's carriers are of the size of guilds the design has. The land's plants grow 1,440–1,900 g of dry matter per m² a year, and the carriers bring back the phosphorus that left with the rock, so plants that feed them fruit are paid in what they lack. The slow and fast erosion cases scale every number here by 0.2 and 2. At 300 trips and a quarter of the body the albatross-scale carriers number 25 million; at 3,000 trips and a load equal to the body, 0.6 million.

**Staying put** *(computed from the terrain study's range of soil creep; the slope at which loose soil slides, about 30°, is a design guess)*. Returned where it left, the rock is a layer of 5.6 mm a century, 11 mm as loose soil, about what fell as loess on China's plateau. It cannot be piled on the main crests. Soil creep passes downslope a flux of its diffusivity times the slope, so a line of returned rock that feeds the ground for a distance to each side needs a slope of the loss times that distance over the diffusivity, and at 30° it can feed only so far:

| Soil creep | Slow erosion | Central | Fast |
|---|---|---|---|
| 0.003 m² a year, humid forest | 154 m | 31 m | 15 m |
| 0.01 m² a year | 514 m | 103 m | 51 m |
| 0.04 m² a year, dry grassland | 2,055 m | 411 m | 205 m |

Feeding a kilometre to each side would need slopes of 1.4–19, which no soil holds. So the rock has to be set down every 60–800 m at the central erosion, on every spur and slope: the smaller carriers of the table above stand that close, and the great ones would have to scatter a thousand loads a year over the 7 km² each serves. What is set down more thickly has to stand as rock, and that is the work of cement. The rock that leaves the land in solution is the stock for it: the seas would otherwise shed 150 Mt of calcite, 95 Mt of magnesite and 280 Mt of opal a year, lime at 4–38% of the particles' mass (8% at the central erosion) and opal at 4–44%. The carriers' own nitrogen is no source of cement. Bacteria that split urea lay one calcite for two nitrogen atoms, and the 0.2–0.7 Mt of nitrogen a year in the albatross-scale carriers' food could lay 0.7–2.4 Mt of calcite, 0.02–0.07% of the rock they return. What can lay it is sized under "The cement" below.

**The loop holds itself** *(a two-variable model: [carriers.py](carriers.py); its rates are design guesses)*. Let sediment arrive at a steady rate, let carriers take it at a rate that saturates with the deposit standing in reach, and let them multiply in proportion to what they find and die at a steady rate. The loop then rests where the carriers return exactly what arrives, whatever the erosion: doubling the erosion doubles the carriers and leaves the standing deposit unchanged, since the deposit at rest is set by the carriers' own rates alone. Every setting with growth above death gives a stable rest, approached through damped swings.

| Carriers' yearly increase and loss | Deposit standing at rest | Settles with a time of | After losing half the carriers, return within 5% in | After the erosion doubles |
|---|---|---|---|---|
| 5% and 2% (slow) | 67 years of supply | 222 years | 545 years | 238 years |
| 10% and 3% | 13 years | 37 years | 108 years | 74 years |
| 30% and 10% (fast) | 5 years | 15 years | 38 years | 27 years |

At rest the carriers find 30–40% of their fill, so the loop keeps 13–18 Mt of them where full work would need 5.4. Against these decades and centuries the waters are slow: with nothing returned, a hundredth of the lakes' area takes 11,000 years to fill, a tenth 460,000 years and half 10 million. A loop that fails for a thousand years costs the lakes a thousandth of their area, and the sediment is still in the deltas when the carriers return.

**What it does for the waters** *(computed, each lake a box that does not change as it fills, sediment at 1,000 kg/m³)*. The share of the lakes' area still open:

| | 100,000 years | 1 million | 10 million | 100 million |
|---|---|---|---|---|
| Nothing returned | 96% | 84% | 51% | 10% |
| Flyers and waders (8–42% of the lakes' rock returned) | 96–97% | 85–88% | 53–64% | 11–15% |
| With divers at the delta fronts (34–95%) | 97–99% | 87–97% | 61–89% | 14–66% |
| With divers on floors shallower than 200 m (53–96%) | 97–100% | 90–98% | 68–91% | 17–71% |

Over hundreds of millions of years every percent left on the floors counts: lakes whose rivers bring sand are held, and lakes whose rivers bring clay are held only with divers that work the deep floors. The seas hold 12.0 million km³. The land's rock alone would fill them in 17 million years (Orientale) to 206 million (the South Pole–Aitken Sea), and their own stone, 0.23 km³ a year, in 52 million; together 28 million. So the seas last only if the stone is laid within reach, upstream by living barriers and in the shallows by reef builders, where carriers can take it.

**The wild and the managed.** At a thousand trips and a load 0.6 of the body, the wild land's rock takes 1.0 million carriers of the albatross scale with 30% of the land wild and 1.65 million with half, or 1.7–2.8 billion of a crow's weight. The hybrid and used land's 1.1–1.3 Gt a year would cost machines 0.25–1.6 GW at the terrain study's rates, and the same carriers would serve it where people prefer. The open land beyond people's need, 0.4–0.9 Gt, belongs to neither aim yet. Because most lakes are fed from land of several kinds, the wild's carriers and the managed land's upkeep meet at the same deltas.

## The return on the map

**One moment, the whole surface** *(computed: [placement.py](placement.py); the rules and ranges are design guesses, and a range here is a square on the grid)*. Every cell's rock is gathered where it first settles, and carriers set it down on the land in range by a rule. The share of the rock set down on the wrong cell, which is half the summed difference between what cells get back and what they lost:

| Rule of placement | Range 10 km | 20 km | 50 km |
|---|---|---|---|
| On worn ground, in proportion to its wear | 30% | 12% | 5% |
| On ground worn by district (wear judged over 10 km) | 32% | 16% | 12% |
| Evenly over the land in range | 33% | 22% | 21% |
| On crests | 64% | 57% | 56% |

The worn-ground rule gives each cell back what it lost once the carriers' range covers the way from the hollows; with 20 km half the land is left unbalanced by under 7 m a million years and a tenth by over 33. It does not return a catchment's own rock to it, and need not: at 20 km only 14% of the rock is set down on land draining to the place it was taken from, since worn ground takes rock from whichever hollow is near. The misplacing is local. Taken in blocks of 240 km, under 1.5% of the rock is in the wrong block under any rule.

**Through time, two regions** *(a landscape model: [evolve.py](evolve.py) and [process.py](process.py); not validated)*. Two squares of the design surface, 730 km from north to south at 1.9 km a cell, are stepped forward: water cuts the ground it crosses by the stream-power law, scaled so that each region's first mean lowering is the rock module's; what is cut is laid in the first hollow on its way from the bottom up; carriers take the rock where it comes to rest and set it down round that place, 86% of it within 20 km. Korolev's country (4° S, 203° E) is the Moon's highest ground, lowering 96 m a million years, and stands mostly under lakes in the design, so it is run with a lake in every hollow that spills and cuts its outlet. The southern highlands (45° S, 350° E) are dry, lowering 48 m a million years, crowded with craters, and are run with dry hollows that keep what enters them.

The design says that carriers set rock down on bare, freshly worn ground. It does not say whether ground can be given more than it lost, and that decides the shapes, so the return is run in two readings: **shared by wear**, each place's rock divided over the ground round it in proportion to the wear; and **no more than it lost**, the rock offered in the same way but each cell taking only what it lost, what it turns away being offered again to ground still short and otherwise left lying. With all the rock reached:

| | 10 million years | 30 million | 100 million |
|---|---|---|---|
| **Southern highlands, shared by wear:** forms kept | 0.99 | 0.93 | 0.17 |
| ground moved (root mean square) | 87 m | 224 m | 1,150 m |
| highest ground, from 2,190 m | 2,040 m | 3,210 m | no limit |
| 1% of the ground has risen over | 290 m | 740 m | 2,860 m |
| the thirty highest summits (median) | −76 m | −225 m | −1,020 m |
| **Southern highlands, no more than it lost:** forms kept | 0.99 | 0.95 | 0.84 |
| ground moved | 97 m | 312 m | 1,130 m |
| highest ground | 2,045 m | 1,970 m | 1,840 m |
| the thirty highest summits | −26 m | −88 m | −718 m |
| largest hollows: depth kept, roundness (0.21 at the start) | 0.99, 0.21 | 0.97, 0.21 | 0.90, 0.23 |
| **Korolev, shared by wear:** forms kept | 0.97 | 0.89 | 0.01 |
| highest ground, from 10,030 m | 9,340 m | 13,330 m | no limit |
| 1% of the ground has risen over | 540 m | 910 m | 1,430 m |
| the thirty highest summits | −237 m | −554 m | −2,240 m |
| **Korolev, no more than it lost:** forms kept | 0.95 | 0.60 | 0.11 |
| highest ground | 9,045 m | 7,930 m | 4,020 m |
| the thirty highest summits | −90 m | −1,860 m | −7,310 m |

"Forms kept" is the correlation between the starting ground and the later ground once each is taken against its own mean over about 30 km: 1 when the map's forms are unchanged. With nothing returned it is 0.3–0.5 at ten million years and 0.1–0.2 at thirty, the hollows filled and their outlines broken (roundness 0.4, where a disc is 0).

- **In plan.** For ten million years the ground moves by 90–350 m and the map is its own under either reading. After that, shared by wear, the rock gathers where wear is fastest: on the walls and rims of craters in the dry region, which become ring ramparts round every crater, and on the steps where rivers drop to lakes in the wet one, which become walls along the shores. The plains and uplands between them are drained. Taken no more than it lost, the crowded southern highlands keep their craters, rims and plains for a hundred million years, changed mainly along the model's open edge; Korolev's wide uplands, tens of kilometres from the nearest lake, are not reached by rock set down within 20 km of where it lies: the region's ground lowers by 0.8 km on average in thirty million years, its crests by 1.4 km.
- **In height.** No reading grows the map's own peaks. The highest summits lose height in every run, by 26–240 m in ten million years, because they stand farthest from the hollows. What grows, shared by wear, is new: ground on steep walls and shore steps that starts low or in the middle of the range and climbs past the old summits, by 1–3 km within thirty million years in the model's worst cells. The spread of heights narrows from the top and gains a thin tail of spikes. The model sets no limit to a slope, so its spikes stand for ground that in fact would steepen to what its cement holds and then shed. Taken no more than it lost, the spread of heights keeps its shape and loses a little from the top.
- **Craters.** In the southern highlands the largest hollows keep their depth and roundness for thirty million years under both readings, and for a hundred million taken no more than it lost. In Korolev they keep them shared by wear (depth 0.99 of its start at thirty million years) and lose a third of their depth taken no more than it lost, as the country round them wastes. Shared by wear the craters grow sharper than they began; with nothing returned they are gone.
- **Reach decides as much as placement.** In Korolev's lake country, with the rock of the delta fronts reached (95%), forms kept is 0.88 at ten million years and 0.45 at thirty; with flyers and waders alone (42%), 0.57 and 0.20.
- **Range.** Shared by wear with 86% of the rock set down within 50 km in place of 20, forms kept at thirty million years is 0.94–0.98 in both regions under either water.

**Other rules, for comparison** *(evolve.py)*. Rock set down on crests runs away within ten million years (relief 3–57 times its start by thirty). Rock set down evenly, or by wear judged over a district of 10 km, blurs the map: forms kept 0.2–0.6 at thirty million years, the hollows' depth 0.4–0.9 of its start. Pictures of every case are on the research drive, `terrain-return/evolution/`: the files ending `_as_designed_shared_by_wear.png` and `_as_designed_no_more_than_it_lost.png` show each region's map at the start and after 30 and 100 million years, where the ground rose and fell, and how it is spread over height; `korolev_wet.png` and the three like it set the rules side by side.

**What the model cannot see.** A cell is 1.9 km, so where on a slope rock is set down, and soil creep, are below it. No slope fails, rain is even, the rock has one strength, and returned rock wears as the ground does. The region's edge is open: rock that crosses it is offered only to ground near the edge, which wastes the margins of the long runs. Straight valleys aligned with the grid in the runs with nothing returned are the grid's.

## The cement

*(computed: [cement.py](cement.py), the arithmetic of a delegated sizing, read and run before it was adopted; its outside sources were read mostly as search summaries, as [sources.json](sources.json) records)*

**Thinly spread rock needs a living cover and no mineral cement.** A year's return is a layer 0.11 mm thick as loose soil. Bare on dry land it would be stripped in 2–10 months at the storm-water study's upper rates, and its loose grains blow at a wind of 3.3 m/s where storms gust to 5–9. A living cover with 1.5–7.3 g of organic binder per m² a year, 0.08–0.5% of the plants' growth, holds it, if the load is set down on that cover and not on bare ground.

**Rim rock needs little strength at lunar gravity.** A layer steeper than its friction angle stands by cohesion in proportion to its weight, a sixth of Earth's, while a cement is as strong on either world. Each metre of thickness needs 0.2–0.8 kPa dry and 0.9–1.2 kPa saturated at slopes of 35–60°, against 1.1–5.1 kPa dry on Earth. A cohesion of 10 kPa holds 16 m at 45°, against 2.7 m on Earth, and the roots of a natural forest give 26–94 kPa. On two tentative anchors from sand cemented by bacteria, 12 kPa takes a cement of 0.2–4% of the mass and 120 kPa 2–8% *(design guess)*.

**The stock is there; water decides whether carbonate stays.** Weathering frees calcium worth 10.8 g of calcite per m² a year, 7.4% of the returned rock, and the seas' calcite carried up with the mud would add 4.6%; highland particles carry iron worth 6.3% of their mass as goethite. A cement of 3% needs 4.4 g per m² a year and 0.5 g of fixed carbon. Water leaving the dry classes of land can carry off 0.02–18 g of calcite per m² a year, so a 3% cement lasts 390–21,000 years in semi-arid ground. Water leaving the humid land, 42% of all the land, can carry off 37–169 g, so there it lasts 23–105 years.

**The rates.** Earth's soils lay carbonate at 1–5 g per m² a year, which would cement a century of return in 88–438 years, and iron crusts thicken 3–140 times more slowly than the rock returns. The fast calcifiers, stoneworts at 440–690 g per m² and tufa at millimetres a year, work under water held above saturation, which the uplands lack. One fed trial cemented a metre of crushed iron crust in 16 months.

So the dry 58% of the land can rest on carbonate laid 1–15 times faster than Earth's soils lay it, and the humid 42% needs an iron cement laid far faster than Earth's, or roots alone. No organism turns upland soil to rock within centuries, and none is known to lay silica as cement, though opal is the largest stock.

## What the design asks of the organisms

The sizing turns the aim into requirements, each a target no organism meets today:

1. **A carrier that takes sediment where water leaves it and sets it down on the ground it came from**, the dry hollows and delta plains first, within 30–400 m of every slope: many small carriers, or great ones that scatter their loads as birds that cache seeds scatter theirs.
2. **A choice of bare, freshly worn ground for its nests or caches**, so that the return goes where the loss is. This is what would balance the loop place by place; the model here balances only the whole.
3. **Numbers that follow the sediment**: carriers that multiply where sediment stands, as the loop assumes.
4. **Upland plants that feed it**, 0.05–1.2% of their growth by the carrier's size, in exchange for the phosphorus returned.
5. **Divers for the delta fronts and floors**, of a bird's breath to 200 m and of a seal's or whale's beyond.
6. **Builders that lay the dissolved rock within reach**, in barriers upstream and reefs in the shallows.
7. **Placement that stops when the ground is made good.** A carrier that goes on setting rock down where wear is fastest builds ramparts and drains the plains; one that leaves ground once it is covered keeps the map.
8. **A range that covers the way from the hollows to the highest ground**: 20 km serves country crowded with hollows, and wide uplands need 50 km or more.
9. **A cover that takes the fall and binds it**, with 1.5–7.3 g of binder per m² a year, and for rim rock a carbonate cement on dry land and an iron cement, or roots, on wet land.

## How far Earth's organisms go toward each part

| Part | On Earth | Source |
|---|---|---|
| Lifting sediment in bulk | Gray whales excavate at least 172 Mt a year on the Bering shelf, about three times the Yukon River's load | Nelson and colleagues 1987 |
| Building with carried earth | A hornero of about 60 g builds a mud nest of up to 80 times its weight each season, gathered nearby | a popular account of a study |
| Carrying loads home over long ways | Seabirds fly meals to a fixed nest over hundreds of kilometres | *recalled* |
| Bound earth that lasts, evenly spread | About 200 million termite mounds over 230,000 km² of Brazil, 690–3,820 years old | Martin and colleagues 2018 |
| Sand turned to stone by bacteria | Treated sand reaches about 1.3–19.6 MPa in engineering trials | reviews of the method |
| Living cover that slows the loss | Sri Lanka's steep, wet, forested highlands lose 2–11 m a million years; the Great Smoky Mountains 25–30 | von Blanckenburg and colleagues 2004; Matmon and colleagues 2003 |

The last row bears on the size of the job: the central case here lowers the land 56 m a million years, so cover as good as Earth's best would cut the rock to return several times.

## How it is set up

- **The rock.** Each dry cell's lowering is the terrain study's central rate at its latitude times the cell's local relief over the mean relief of the dry land in its band of latitude, the relief being the highest less the lowest surface within 5 km. The terrain study spreads its rock by relief and by rain to the power 0.31; this surface's rain is not carried here, so the bands' rates stand for the rain.
- **Where it settles.** Down the design routing's receivers to the first sea, lake or lowest point that passes nothing on; the way back is the great-circle distance and the height from that place's water or ground to the cell.
- **What the lakes keep.** Siyam's form of Brune's curves at lunar settling, as the [sea nutrients](../sea_nutrients/README.md) study carries particles; a lake's years to fill are its volume at 500–1,500 kg/m³ over what it keeps.
- **Reach.** The river deltas screen's central shares for the mature Moon's soils; sediment beyond the delta front is taken as spread evenly over the floor.
- **Flight.** The sky ecology study's rule and ranges, with the carrier's body flown out and up with its load and flown back level only where a glide from the height gained falls short. Upkeep is that study's allometry of birds' field metabolism by weight.
- **Food.** Usable energy of 15 kJ per gram of dry food and 1–3% nitrogen *(design guesses)*; the land's plant growth is the biomass comparison's.
- **The loop.** Sediment S and carriers B: dS/dt = supply − duty·B·S/(S + H), dB/dt = B·(growth·S/(S + H) − death), with the duty 600 body weights a year and H the supply of 10–100 years.

## Assumptions and limits

- **No erosion is modelled on this ground.** The lowering is a range for Earth carried to the Moon by the conservation study and placed by relief. Living cover, the reach of lunar floods and the soils' own strength could each move it.
- **The combination has no precedent.** Each habit exists; no animal joins them, and none carries mineral loads near its own weight. Loaded flight at lunar gravity is unmeasured.
- **Lakes are boxes.** They do not shrink, shoal or pass more on as they fill, and their own deposits are not added, so the years to fill run long beside the terrain study's.
- **The loop is a model of the whole.** It shows that the return can follow the supply; it does not show that the rock goes back to the slope it left. Where it is set down depends on the carrier's choice of ground.
- **Dust.** Sediment set down dry on open slopes could blow, and the seas' phosphorus budget rests on little dust. The returned rock must be bound or wet when it is placed.
- **The landscape model is coarse.** See "What the model cannot see" above. It shows the kinds of change each reading makes and how large; it predicts no place.
- **The cement's sources were read thinly**, and its link from cement to strength is a guess on two anchors. The deep divers are named and not sized.
- **No independent review.**

## Files, running and checks

| File | What it is |
|---|---|
| [rock.py](rock.py) | The rock on the design surface: where it settles, the way back, whose land, the lakes and floors |
| [carriers.py](carriers.py) | Flight cost, reach, creep's limit, the urea route, the loop and its recovery |
| [run.py](run.py) | The sizing from the products |
| [placement.py](placement.py) | The rules of placement on the whole surface at one moment, and the regions cut for the runs through time; writes [results/placement.json](results/placement.json) |
| [evolve.py](evolve.py) | The landscape model and the rules side by side on two regions for 30 million years; writes [results/evolution.json](results/evolution.json) |
| [process.py](process.py) | The return as designed, in two readings, for 100 million years; writes [results/process.json](results/process.json) |
| [cement.py](cement.py) | The cement's arithmetic, with its own checks |
| [test_run.py](test_run.py), [test_map.py](test_map.py), [test_cement.py](test_cement.py) | 59 tests, and the regeneration of rock.json from the research drive when `TERLUNA_ROCK_RETURN_REGENERATE=1` |
| [sources.json](sources.json) | Each source and how far it was read |

`ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.rock_return.rock` reads the design routing and the land-use grid on the research drive, as does `placement`; `python -m research.studies.rock_return.run` reads only committed products. `evolve` and `process` read the regions `placement` writes to the drive and want numba for speed (the research drive's `venvs/optics`); each takes about twelve minutes on six cores, and their snapshots and pictures stay on the drive. The tests check the arithmetic of each function, that the loop rests where the return equals the supply and is stable over its settings, that the stored sizing is what the code gives, that the way back agrees with the terrain study's within a tenth, and the figures this page states.
