# The three domains together

The resources, ecology and provisioning branches merged on 12 October 2026, and four refreshes of main's studies
followed on the merged tree ([integration record](../../integration/domains.md)). This study fits the branches together
at the author's choices: distribution C′ (8 billion on the Moon, 4 billion on the ground and 4 billion aloft, under the
compromise rule of 30% of the land kept wild and another 20% hybrid), fed by the food blend (power-heavy 40%, power-made
50%, fully lamp-lit 10% of the diet), with heritage meat from the homesteads. Its [ledgers](ledgers.py)
([results/ledgers.json](results/ledgers.json)) recombine the committed products in five places where the domains meet:
heat, land, methane and N2O, phosphorus, and the sky's layers, with power from the array divided between the towns and
the ground. They are arithmetic on screen products and take under a second; every figure keeps its source's evidence
state. The distributions study placed C′ with the power-heavy mix, so here the blend's food takes its place in the
totals; placing the blend's fields band by band is that study's next run (next work, item 2).

## What each branch brought

- **Resources:** where every material comes from and how its cycle closes: the imports from the Solar-System-wide
  operation, lunar rock, the inventories and flows, nitrogen extraction, the industry's kilns and weathering, and
  phosphorus from the seas, with the lake deltas dredged as the people's distribution system beside a natural cycle.
- **Ecology:** the living sky and the land: the sky's bands and the people's share of them, the aerophytes and sky
  reefs with their guards, symbioses and living sails, the sky's food web and phosphorus, the taxonomic system (on
  hold), living light, and the machines in the living world with the starting constraints on what may live wild.
- **Provisioning:** people on the Moon: the heat and computing ledgers, the inhabited volume and the drifting towns,
  ways of living and daily life, eight food systems in seven mixes and the author's blend, cities on the land and in
  the lava tubes, the distributions of 25–30 billion with the author's C′, the wild share and the megaparks, and the
  network.

## The author's Moon in five ledgers

At 8 billion on the Moon, pairs from the low-cost to the high-cost end:

| Ledger | C′ with the blend | C′ as placed with the power-heavy mix |
|---|---|---|
| Heat on the Moon, with fusion and close computing | 40–284 TW, 1.1–7.5 W/m² | 28–133 TW |
| Added dimming that answers it | 0.35–3.9% | 0.24–1.9% |
| The same with power from the array | 30–145 TW; 0.26–2.0% | 22–77 TW; 0.19–1.1% |
| Dry land used, with 1 kg of heritage meat a person | 1.24–2.89 million km², 5.3–12% | 1.95–7.71 million km², 8.4–33% |
| Wild share of the land outside the kept 30% and the parks | 75–89% | 34–83% |
| Methane at the soils' central sink, steady state | 22–29 ppm, people 7–30% of the sources | |
| N2O, rise a year | 0.36–0.55 ppb, people 12–23% of the sources | |
| Phosphorus dredged back, share of what the lakes receive | 20–51% and 8–27% at the two land losses | the same |

## Heat

| With fusion, TW | Low-cost end | High-cost end |
|---|---|---|
| Residents, from efficient to high-income use | 12.1 | 36.2 |
| Food: the blend | 14.8 | 171 |
| Sky towns: control | 0.55 | 8.4 |
| Sky towns: lifting hydrogen remade | 3.4 | 53.7 |
| Trips down, food lifted and water | 0.16 | 0.32 |
| **Total** | **31.0** | **270** |
| Close computing, a hundredfold today's data centres a person | 9.2 | 14.0 |
| **With close computing** | **40.2, 1.06 W/m², dimming 0.35%** | **284, 7.5 W/m², dimming 3.9%** |
| The ring fleet's infrared glow | 5.7 W/m², 1.9% of the design sunlight | 9.9 W/m², 3.4% |

The blend buys its small fields with energy: C′'s heat with fusion rises 1.7–2.3 times over the power-heavy mix's,
and food becomes 48–63% of it, the largest part at both ends. Power from the array leaves its conversion losses in
space, so the same Moon gives 30–145 TW with close computing, about half the high-cost end's heat; with the array's
power the residents' 11.0 TW stay just above food's 10.2 TW at the low-cost end. The planned climate keeps the
investigated design case near 295 K, so every terawatt is answered: people's heat at 0.0087–0.014% of the sunlight a
terawatt, provisioning's factor, and the fleet's glow as provisioning gives it, its W/m² over the design sunlight's
293 W/m². The glow is the larger load. People and the fleet together would need 2.3–7.3% of the sunlight taken away
if dimming answered both, 2.2–5.4% with the array's power. The fleet's answer is still open: films that absorb less, a
heat mirror on its Moon-facing faces, redirecting the glow, or added dimming ([plan](../../plan.md), next work
item 3). The parks' automation and amenities, which the register leaves to be sized, are not in the ledger.

## Power from the array: a distribution recommendation

The author chose power from the array for food and the sky towns, given a good recommendation for distributing it
(register, 12 October). Everything received ends as heat on the Moon, so the intake is the heat ledger's array parts
(*screen*, [ledgers](results/ledgers.json), `array_power`):
- **The sky towns take theirs on their own tops.** Their control, hydrogen remade, trips and the aloft half of the
  residents come to 8.3–45 TW. The 696,000 towns' tops, a kilometre across, make 547,000 km², so the towns take it at
  15–82 W/m², without lines or night storage. From the ring, 18,300–19,900 km away, a transmitter 2.2–2.4 km on a side
  puts 95% of a 5.8 GHz beam on one top (5.1–5.6 km at 2.45 GHz), a fraction of one 10 km tile.
- **The ground takes the rest at its loads.** Food, the ground's half of the residents and close computing come to
  22–100 TW. At the 1978 reference rectenna's 49 W/m² mean that needs 0.45–2.05 million km² of receivers, at its
  230 W/m² peak 96,000–436,000 km², 0.4–1.9% of the dry land; strung over the blend's fields the same power arrives at
  33–54 W/m², near that mean, beneath a mesh that shades the crops. A beam from one tile lands in 220–240 m at 5.8 GHz,
  so receivers can sit where the load is: at the power-made food plants, the lamp-lit farms and the cities' industry,
  joined by a grid, with several tiles behind each receiver for the uptime the author requires.
- **What it costs.** All of these intensities exceed the 10 W/m² that ICNIRP (2020) sets for the public, so the
  receivers and the towns' tops are closed spaces, with people beneath the rectenna's shield and flyers kept out of
  the beam columns; what aerophytes drifting through a column take is to be checked. The ring stands above the
  horizon everywhere except within about 5° of its plane's poles, where the grid serves. The night stays dark, and
  the Moon takes only the receivers' losses and the air's absorption as heat, which the array's factor already counts.

The recommendation, for the author: beams straight to the towns' tops; ground receivers at the loads at up to the 1978
peak, joined by a grid; and a choice between closed receiving land, 0.10–0.44 million km² at that peak, and meshes
over the fields at the cost of their shade. A screen of receivers, beams, the grid and reserves through the lunar
cycle firms it up (provisioning, question 3).

## Land

The dry land is 23.4 million km². The compromise rule keeps 7.02 million km² wild and 4.68 million km² as megaparks
and great gardens, leaving 11.7 million km² for everything else.

| Use, million km² | Low-cost end | High-cost end |
|---|---|---|
| Settlements of C′: the metropolis, tower districts, highland towns, scattered communities | 0.26 | 0.74 |
| The blend's fields: rain belt | 0.50 | 1.43 |
| The blend's fields: highlands | 0.18 | 0.44 |
| Heritage pasture at 0.5, 1 and 2 kg a person a year | 0.15, 0.30, 0.60 | the same |
| Used, at 1 kg | 1.24 (5.3% of the dry land) | 2.89 (12%) |
| Wild outside the kept land and the parks, at 0.5–2 kg | 87–91% | 73–77% |
| Wild in fact, the kept 30% included | 73–75% | 66–68% |

The blend's fields take 27–40% of the power-heavy mix's land, so the wild share outside the kept land rises from
34–83% to 73–91% across the heritage cases. Land wild in fact is the kept 30% and the land no settlement, field or
pasture uses; the unused part stays open to later use. The blend grows no staple crops in the middle latitudes' wet
lowlands. The homestead belt there, 1.39 million km², holds the homesteads among their own heritage and hobby farms:
11%, 21% or 43% of it at 0.5, 1 or 2 kg of meat a person a year. With power-heavy, 16 of C′'s 32 corners overfilled
a band at the high-cost end; the blend's smaller fields leave room in the totals, and whether every corner then fits
band by band is the distributions study's next run.

## Methane and N2O

| | Low | Central, low-cost end | Central, high-cost end | High |
|---|---|---|---|---|
| Methane sources, Tg a year: the wild's lakes and wetlands | 12.8 | 58.9 | 58.9 | 218 |
| Farms and waste, the lunar design's 0.3–3 kg a person | 2.4 | 2.4 | 24.0 | 24.0 |
| Heritage cattle (0.5, 1, 1 and 2 kg a person) | 0.81 | 1.71 | 1.71 | 3.4 |
| Soils' uptake, Tg a year per ppm | 4.77 | 2.92 | 2.92 | 1.07 |
| Steady state | 3.4 ppm | 21.6 ppm | 29.0 ppm | 230 ppm |
| Removal that holds Earth's 1.85 ppm, Tg a year | 7.2 | 57.6 | 79.2 | 244 |
| Share of the atmosphere a full-capture removal from open air passes a year | 0.22% | 1.8% | 2.4% | 7.5% |
| Least power for complete separation | 0.49 GW | 3.9 GW | 5.4 GW | 16.7 GW |

| N2O, Tg N a year | Low-cost end | High-cost end |
|---|---|---|
| Sources: the wild's soils and waters / the blend's farms | 0.96 / 0.13 | 1.30 / 0.38 |
| Rise a year, and the time to reach Earth's present 0.34 ppm | 0.36 ppb, 950 years | 0.55 ppb, 620 years |
| Removal that holds any level | 1.09 | 1.68 |
| Left if complete denitrification cuts every source by 50–95% | 0.055–0.55 | 0.084–0.84 |
| Share of the atmosphere a full-capture removal passes a year at 0.34 ppm; least power | 0.11%; 48 MW | 0.16%; 74 MW |

The wet Moon's own waters and soils make most of both gases. In the bare column under the stack the air holds no OH
(the joint synthesis), so methane leaves only into soils, here every kind of dry land, fields, towns and parks
included; a living column may make some OH (ecology, finding 9). N2O lasts 45,000 years to 3.9 million years on the
light through the shield's gaps, so on any human timescale it has no sink (finding 7). People add 7–30% of the
methane at central values and 12–23% of the N2O. The heritage cattle are counted by their enteric methane; their
manure and pastures are not.

Holding Earth's present levels would take 58–79 Tg of methane a year at central values and all of the N2O. The
power that separation needs is small; the burden is the air: a removal from open air would have to pass 1.8–2.4% of
the whole atmosphere through its contactors each year for methane at full capture, and 0.11–0.16% for N2O. The
author chose biological means coupled with technologies that make up the difference; the research direction recorded
with that choice puts biology on the land and seas first, treatment where the gases are concentrated next, and
removal from open air last, and these figures show the burden open air would carry. The denitrification bracket is
one N2O-respiring bacterium's cut on Earth's farm soils (Hiis et al. 2024, Nature 630, 421–428), applied here to every
source, the wild waters included. How far each gas may rise is set by its warming in the deep column, still
uncomputed (next work, item 1); the levels in these tables are Earth's, for reference.

## Phosphorus

| Tg P a year | Wild land's loss 1.47 | Wild land's loss 4.16 |
|---|---|---|
| People's losses, the blend's two ends | 0.35 / 1.44 | 0.35 / 1.44 |
| Into the waters | 1.82 / 2.91 | 4.51 / 5.60 |
| Kept by the lakes | 1.75 / 2.80 | 4.32 / 5.36 |
| Dredged back to the farms, share of what the lakes receive | 20% / 51% | 8% / 27% |
| Carriers returning the wild land's dissolved loss | 0.021 | 0.094 |

The blend does not change people's phosphorus: the diet and its losses are the food study's, eaten 3.5–4.8 Tg a
year with 0.35–1.44 Tg lost, so resources' steady-state ledger stands for C′ with the blend. The wild land's loss
takes Earth's river concentrations at the lunar-cycle register's former runoff; at the rerun register's runoff it
would be about a fifth lower, and the dredged share 10–57%. The towns' sewage holds 2.1–4.3 Tg a year to recover, and
the blend's power-made food cycles 1.1–3.5 Tg through microbial biomass in closed plants. The sky's own need at 1%
aerophyte cover (sky ecology) is 0.0003–0.050 Tg a year for a mature sky with designed pools, 0.004–0.061 Tg with an
Earth forest's floor and 0.062–0.42 Tg if nothing were caught aloft, and 0.009–0.53 Tg while the sky grows: small
against the land's loss unless the aerophytes recycle poorly. The developed fines let the productive coasts' light
reach 12.4 m to 1% (9.5 m before) and the river mouths' 2.3 m (1.6 m), which deepens the seas' lit floor beside the
deltas that are dredged.

## The sky's layers

| Where | What is there |
|---|---|
| 0–10 km above sea level, the low sky | 13% of the air; cloud on the ring 0.2–3.4% of the day |
| Up to 2.5 km above the ground | Local flight: wings, gliders and canopies to 0.3 km, sky boats and air taxis to 1.5 km, ferries to 2.5 km; over the summit metropolis, whose ground is 11.6 km above sea level, that is 11.6–14.1 km |
| Up to about 10 km above the ground, by day | Soaring gliders in the daytime mixed layer (sky fleet) |
| About 10 km above sea level | The sky towns' reference height: 696,000 towns of about 5,750 people (50 t of district a resident) in C′, 5.2 km apart if spread evenly, 4.3 km in the busiest band |
| 10–22 km above sea level, the storm layer | The sky reefs at 10–20 km; storms top out at 22 km at the median |
| 20–30 km above sea level | Regional ships, with their docks 10–15 km above the summit (sky fleet) |
| 22–35 km above sea level, cool air | Aerial life sparse above 26 km; the network's drifting mesh and fixed posts near 30 km, where the ring puts cloud 1.9% of the day and 0.9% of the night; one storm in ten tops 38 km |
| 35–45 km above sea level, the flight band | Liners and freighters, storm cloud 0.6–0.7% of the time; the summit port's crown at 35.6 km |
| 17,900–20,200 km above the surface | The ring fleet, at radii of 19,600–21,900 km, 0.13 s away by round trip |

The sky is shared between 10 and 30 km. The summit metropolis and the air traffic over it, 11.6–14.1 km above sea
level, the sky towns near 10 km, the sky reefs at 10–20 km, gliders and the storms all use that layer, and regional
ships fly at 20–30 km. The towns share the reefs' wind tracks, and the author's rule has towns give way to reefs. By
the distributions study's encounter rule, with reefs a kilometre or two across making up half the aerophyte cover, a
town passes 0.04–1.2 reefs a day at 1% cover, drifting 0.1–1 m/s apart; the distributions study's 13–132 passes a day
count 100 m bodies making up the whole cover. Sizing the towns' steering reserves for those meetings is next work
(item 5). The network's 30 km sits above the median storm and below the tallest: a drifting node spends about 1–2%
of its time in cloud, one storm in ten reaches 38 km, and the largest storm on the ring tops 84 km. The heights come
from the products named in the [ledgers](results/ledgers.json), the reefs' band from sky ecology, the network's and
the ring fleet's from [provisioning/network.md](../../../provisioning/network.md), and the regional ships' and gliders'
from the [sky fleet](../sky_fleet/README.md).

## What it means for the whole

The author's choices fit together on the land. C′ with the blend leaves 66–75% of the dry land wild in fact, the kept
30% and the land nothing uses, beside the megaparks' fifth, at every corner of the totals. Heat is the cost of that
land, and it is the largest tension between the domains: the blend's power-made and lamp-lit food make food the
largest part of people's heat with fusion, and at the high-cost end people and the ring fleet together would need
several per cent of the sunlight taken away, which the plants then lose. Power from the array halves people's share
at the high-cost end.

The gases are first an ecological design problem. The wild's lakes, wetlands and soils make most of the methane and
N2O, and a wet Moon under a shield with no ultraviolet window keeps both for centuries or longer, so the biomes' own
budgets (oxidizers in the waters, denitrification that ends in N2) decide the air before any machine does. The
warming of both gases in the deep column sets how much of that work is needed.

Phosphorus closes in steady state: people's losses go to the lakes and come back by dredging, 8–51% of what the lakes
receive, while the natural cycle keeps its carriers and the lakes bury the rest. Between 10 and 30 km the sky is
shared by the summit metropolis, the towns, the reefs, the gliders, the regional ships and the storms, and the
network flies above most of it.

## Next work, in order

Machine time is for this laptop, one core, alongside the author's other work.
1. **The gases' warming in the deep column** (ecology, question 3): the line-by-line radiative–convective column
   with methane and N2O lines added (HITRAN, fetched and hash-recorded like the column's other inputs) and swept over
   the ranges above. A few hours. It sets the levels the gases may reach, and so how much biology and technology
   must remove.
2. **C′ placed with the blend:** the distributions study takes the blend as its food and places it band by band,
   with the homestead belt holding the heritage farms. Minutes; the distributions README follows.
3. **The fleet's night and its glow** (plan, item 3; the planned branch research/array-industry): the night-side
   attitude, films and a heat mirror sized against the glow, with this heat ledger as the people's part of the same
   budget. About an hour of CPU for each set of kept orbits.
4. **The wet Moon's methane and nitrogen:** methane from the lakes by depth with oxidation in their waters and
   sediments, and the biomes' nitrogen cycle ending in N2. Screens of minutes, with the literature.
5. **Traffic and the network:** the sky shared between 10 and 30 km, with the towns giving way to reefs (passes,
   moves and their energy at 0.1–2% cover), the metropolis's traffic and the regional ships; and a constellation model
   of the mesh with the ring's storms and the link budgets. Under an hour.
6. **The aerophytes** (ecology's next): one assembled wet membrane, then trajectories through the design winds
   (hours on the saved winds); one specified living sail; a population model with phosphorus return; measured
   tissue phosphorus and the seas' food web.
7. **The night:** the dark-night requirement's value, lamps for the dark belt and a beetle-class glowing plant.
8. **Water and energy for the parks and the high-latitude towns:** the parks' water, automation energy and night
   light, and the high-latitude towns' water.
9. **The signal vocabulary** of the machines in the wild, with its channels and regional dialects.
10. **The sea-appearance illustrations** re-rendered from the developed-fines waters (the reference renders in
    minutes to hours; the generated images wait for the author's look).
11. **The editorial pass** on the Codex-written prose against the author's standards.
12. **The taxonomy**, on hold at the author's request.

## Decisions

The author decided the list on 12 October
([register](../../decisions.md#the-domain-synthesiss-decisions-12-october-2026)):
1. The levels at which methane and N2O are held wait for their warming in the deep column, which is the next work.
2. Food and the sky towns run on power from the array, given a good distribution recommendation; the one above waits
   for the author.
3. Added dimming can be had by extracting more of the sunlight; the fleet's glow is the harder problem to work out.
4. Heritage meat is 1–2 kg a person a year, 21–43% of the homestead belt before high-efficiency feed production cuts
   it; wild culls and fishing can supply part of it if parasite risk is kept low.
5. Sky towns fly in the reefs' layer. If sky reefs become part of the useful balancing system they get more
   preference; otherwise sky towns get a bit more precedence. A somewhat higher reef density is acceptable.
6. The network's stability decides whether its nodes ride through the storm tops that reach 30 km or climb round
   them.
7. The taxonomy stays on hold.

## Running it

```
python -m research.studies.domain_synthesis.ledgers      # writes results/ledgers.json
python -m pytest research/studies/domain_synthesis
```

The product binds its producer and its thirteen inputs by sha256: the distributions, the food blend, provisioning's
population, ecology's first screen and people products, resources' sea phosphorus, sky ecology, the inhabited volume,
the sky ships, the sky fleet, the summit tower's form, and the shield's ring layout and integrated ledger. The tests
regenerate it exactly and check the sums it rests on.
