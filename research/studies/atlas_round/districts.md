# Regions and districts: the ruleset and the first draft

The register defines a district by significance and wholeness: a region that reads as a whole and meaningful place and stands apart from the others (16 November 2026). On 17 November the author asked for that definition to be formalized into rules so that a first draft of the atlas could be drawn, and chose its scale, its naming, how it shows, and that the coastal seas are districts of their own linked to their land. This file holds the rules, how the draft applies them ([districts.py](districts.py)), the tongues its names come from ([district_names.py](district_names.py)) and what the draft is not yet.

## The rules

The rules draw on how geographers already draw such places: Kevin Lynch's districts, areas known from inside by a common character; landscape character areas, whose boundaries follow natural lines rather than administrative ones; and the way ecoregions are drawn, with a boundary where several layers change together.

**Wholeness.**
1. *One centre.* A district gathers round one thing it is about: a crater or basin floor, a mare plain, a range or massif, a valley and its river, a lake, a bay or strait, an island.
2. *Natural edges.* Its boundary follows lines the land draws, preferring a shoreline, a drainage divide, a rim crest, a scarp, then the edge of a climate and cover zone; where nothing marks the edge, the narrowest neck.
3. *One character.* Most of its land shares one landform type and one moisture zone, or one gradient its centre organizes.
4. *One piece,* with its lakes and rivers.

**Significance.**
5. *A reason to be named,* in a sentence: its landform, its water, its living character, its works, its heritage. A district whose sentence would repeat a neighbour's merges with it.
6. *A name in the Moon's own tongues* (the author's choice, 17 November): no Latin, English or mainland European names. Each region takes a tongue, loosely as Earth's peoples sit in like landscapes, and its districts are named in it.
7. *Significance is judged.* The data propose the units and their facts; written judgment of each place follows, the draft's own descriptions being short and computed.

**Standing apart.**
8. Neighbours differ in their centre and in at least one of: the waters they drain to, landform type, moisture zone, rock.
9. Their shared edge lies where those differences happen.

**Scale and nesting.**
10. Two levels: regions, the large wholes, and districts within them, never crossing a region's edge.
11. Districts of about 5,000–150,000 km² (the author's choice for the first draft); smaller wholes are features inside a district, larger ones regions; small islands may stand alone.
12. *The coastal seas are districts of their own* (the author's choice), the water within 60 km of the shore, each linked to the land districts that face it and those whose rivers reach it, so that a coast and its land can be planned together, as satoumi needs (the hybrid the author chose).

The page shows regions when zoomed out and districts with their names when zoomed in (the author's choice), with short descriptions.

## How the draft applies them

[districts.py](districts.py) draws on the design surface, today's grid at 8 pixels per degree (3.8 km), on the design climate's moisture map, the LROC mare boundaries, the design routing, the cuts, the falls and the heritage register.

- **Centres:** the IAU gazetteer's named craters of 60 km or more on land, the lakes of 1,000 km² or more, fill-in points about 200 km apart where neither stands, and every island. The gazetteer only places the centres; its names are not used.
- **Edges:** districts grow from their centres over the ground by a watershed, so that they meet on ridges and rim crests. Crossing the edge of a moisture zone costs as much as climbing 600 m, the edge of the mare basalt 400 m. The moisture zones are the map's classes on its log aridity, smoothed over 0.7 of the climate run's grid as the globe smooths its climate fields.
- **Standing apart:** neighbours sharing their moisture zone, their rock and the waters they drain to, whose boundary stands less than 200 m above the higher of their mean ground, merge while the whole stays under 150,000 km².
- **Size:** districts under 5,000 km² join the neighbour they share the most edge with, among those the join keeps within the maximum. Districts over 150,000 km² split round centres spread by k-means, each cell going to its nearest centre over the ground, where every 100 m climbed between cells counts as another cell's length, so that the parts meet on the steep lines. Islands under 3,000 km² join the nearest district across the water.
- **Regions:** the same watershed from the basins of 300 km or more and fill-in points about 750 km apart. Each district goes to the region holding most of it. Regions under 250,000 km² join a neighbour, and regions over 1.2 million km² split by grouping whole districts, until neither remains.
- **Coastal seas:** each cell of the sea within 60 km of the shore goes with the land district whose shore is nearest. The pieces join within one sea up to 40,000 km². A coastal sea is a bay where land within 120 km lies along ten of sixteen directions from its label point, and a strait where two land masses of 20,000 km² or more lie within 60 km of it. Each is linked to the land districts it faces and to those that send it at least 5% of their area's rivers.
- **Descriptions:** computed from the facts. They give the moisture zone and landform, the lakes, the largest river, the steep reaches, the cut and the heritage sites, and the water each district drains to, by name.

## The tongues

Each tongue has a head for each kind of centre (lake, crater, river, plain, coast...) and a modifier for each distinguishing trait (moisture zone, rock, height, lakes, falls, the near side's Earthlight...). A name joins them as the tongue joins words. [district_names.py](district_names.py) lists every word.

| Tongue | Roots | Where | Pattern |
|---|---|---|---|
| Terlunan | the author's conlang v0.6: Tahitian, Māori, Hawaiian, Hakka, Sumerian, Akkadian | shores and islands, and the seven seas | modifier-head, fused: Mitifanga, salt bay |
| Proto-Olivolingvo | the author's Esperanto-based tongue with Hebrew loans | the polar lands and the semi-arid uplands | h' + noun + adjective: H'Naharo Keraha, the frozen river |
| Sahran | Arabic and Tamazight | the deserts | head + modifier: Wadi Nur, Adrar Abyad |
| Tagalog-Hebrew | Tagalog heads, Hebrew modifiers, the author's idea | the lowland rainforests | head, -ng or na, modifier: Lawang Or, Bundok na Lavan |
| Pwani | Swahili, without noun-class agreement | the humid woodlands and savannas | head + modifier: Ziwa Nyota |
| Dala | Kazakh, Kyrgyz, Mongolian | the steppe plains | modifier-head, fused: Karakol |
| Runa | Quechua | the high mountains and their cloud forests | modifier-head, fused: Yanaqucha |
| Satogo | Japanese | the humid hill country | modifier-head, fused, voiced: Kurogawa |

The seas are named in Terlunan: Honuamoana, the Earth Sea (the near-side sea); Itimoana, the little seas (its southern basins); Mitimoana, the salt sea (South Pole–Aitken); Fetiamoana, the star sea (Humboldtianum); Teiteimoana, the high sea (Moscoviense); Pi'omoana, the arch sea (Orientale); Tahamoana, the side sea (Smythii–Marginis).

Terlunan's v0.6 dictionary lacks some words the atlas needs, so these are coined here in its manner, from its donor languages, as proposals: *loko* lake (Hawaiian, since *roto* is taken by inside), *siyenloko* thousand lakes, *awa* river, *awaawa* valley, *mania* plain, *hiwi* ridge (Māori), *tahatai* shore (*taha* side + *tai* tide), *koawa* channel, *matie* grassland (Tahitian), *rere* waterfall, *puna* spring, *teitei* high, *hohonu* deep, *maroke* dry, *maaku* wet, *kohu* mist, *'ute'ute* red (beside *'ere'ere* black and *'uo'uo* white), *galpoka* great crater. The other six tongues are built for the atlas from everyday landscape words of their roots, without sacred or protected terms. Some names match places on Earth, as names from real roots will. All are working labels, for the people who come to live there to keep or replace.

## The first draft (17 November)

On the design climate (the run at sunlight 0.930, years 35–49) the draft holds 540 districts on land in 42 regions, their median 38,400 km² (10th–90th percentile 6,900–120,500), and 148 coastal seas in the 7 seas: 115 open coasts, 23 bays and 10 straits. The tongues take 7 regions each for Terlunan and Dala, 6 each for Proto-Olivolingvo and Satogo, 5 each for Sahran and Tagalog-Hebrew, 4 for Runa and 2 for Pwani, with the seas in Terlunan. The product is `atlas-round/districts/districts_8ppd.npz` on the research drive, recorded in [results/districts.json](results/districts.json). The globe shows it with `--districts` (visualization/atlas/README.md).

**What the draft is not yet.** Its significance is computed, not yet written. A district's centre is often a landform worth keeping, so the draft also lists candidates for the kept-landforms register, which is chosen in this round. Settlement and population follow capacity, and the heritage treatments follow the author's review of the register. The boundaries carry the 8 px/deg grid's steps. Five islands stand under 5,000 km² and three districts up to 4% over 150,000 km². The coastal seas' 60 km is a working width. The draft is redrawn on each new design climate: the districts' moisture zones and the rivers' ends come from it.
