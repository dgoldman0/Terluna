# The living return of rock

22 November 2026. Rain moves the land's rock downhill into hollows, lakes and seas, and a Moon without uplift has nothing that lifts it back. The [terrain study](../terrain_through_time/README.md) priced the return for machines: the deltas dredged and the rock carried uphill. This study asks whether animals could do that work in the wild, by joining two habits that are each common on Earth, carrying earth for nests and travelling far as part of a life, and whether such a loop would hold itself. The register's row of 22 November sets the aim it is judged against: the hybrid and used land may have designed and engineered upkeep, and the wild should if possible keep a dynamically stable equilibrium on its own.

**On the design surface the land loses 3.2 Gt of rock a year (0.6–6.5 across the terrain study's range), and flying a kilogram of it the 9.4 km back and the 1.0 km up costs a carrier 18–73 kJ at lunar gravity, against 280–1,900 kJ for the same way on Earth. Carriers weighing about 5 Mt alive, eating 0.05–1.2% of the land's plant growth, would return all of it, and a loop in which they multiply where sediment gathers returns exactly what arrives and recovers from losing half its carriers in 40–550 years, while the lakes take 11,000 years to lose a hundredth of their area. Energy and food do not limit the living return. Reach and placement do: flyers and waders can get at 46–66% of the rock, divers at the delta fronts raise that to 61–97%, and the rest lies on deep floors; and the rock has to be set down within 30–400 m of the slope it left, which asks for carriers that are many and spread, or that scatter their loads.**

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

Feeding a kilometre to each side would need slopes of 1.4–19, which no soil holds. So the rock has to be set down every 60–800 m at the central erosion, on every spur and slope: the smaller carriers of the table above stand that close, and the great ones would have to scatter a thousand loads a year over the 7 km² each serves. What is set down more thickly has to stand as rock, and that is the work of cement. The rock that leaves the land in solution is the stock for it: the seas would otherwise shed 150 Mt of calcite, 95 Mt of magnesite and 280 Mt of opal a year, lime at 4–38% of the particles' mass (8% at the central erosion) and opal at 4–44%. The carriers' own nitrogen is no source of cement. Bacteria that split urea lay one calcite for two nitrogen atoms, and the 0.2–0.7 Mt of nitrogen a year in the albatross-scale carriers' food could lay 0.7–2.4 Mt of calcite, 0.02–0.07% of the rock they return. The cement has to be laid by organisms that calcify with light, as the terrain study's stonewort meadows carry 440–690 g of calcium carbonate per m² and Plitvice's barriers rise 13.5 mm a year, and by the living crusts and roots that hold ordinary ground. That step is not sized here.

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

## What the design asks of the organisms

The sizing turns the aim into requirements, each a target no organism meets today:

1. **A carrier that takes sediment where water leaves it and sets it down on the ground it came from**, the dry hollows and delta plains first, within 30–400 m of every slope: many small carriers, or great ones that scatter their loads as birds that cache seeds scatter theirs.
2. **A choice of bare, freshly worn ground for its nests or caches**, so that the return goes where the loss is. This is what would balance the loop place by place; the model here balances only the whole.
3. **Numbers that follow the sediment**: carriers that multiply where sediment stands, as the loop assumes.
4. **Upland plants that feed it**, 0.05–1.2% of their growth by the carrier's size, in exchange for the phosphorus returned.
5. **Divers for the delta fronts and floors**, of a bird's breath to 200 m and of a seal's or whale's beyond.
6. **Builders that lay the dissolved rock within reach**, in barriers upstream and reefs in the shallows, and calcifiers that cement what is set down.

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
- **Cement and the deep divers are named, not sized.**
- **No independent review.**

## Files, running and checks

| File | What it is |
|---|---|
| [rock.py](rock.py) | The rock on the design surface: where it settles, the way back, whose land, the lakes and floors |
| [carriers.py](carriers.py) | Flight cost, reach, creep's limit, the urea route, the loop and its recovery |
| [run.py](run.py) | The sizing from the products |
| [test_run.py](test_run.py) | 24 tests, and the regeneration of rock.json from the research drive when `TERLUNA_ROCK_RETURN_REGENERATE=1` |
| [sources.json](sources.json) | Each source and how far it was read |

`ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.rock_return.rock` reads the design routing and the land-use grid on the research drive; `python -m research.studies.rock_return.run` reads only committed products. The tests check the arithmetic of each function, that the loop rests where the return equals the supply and is stable over its settings, that the stored sizing is what the code gives, that the way back agrees with the terrain study's within a tenth, and the figures this page states.
