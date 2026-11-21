# The land's use on the Open Moon: what is fixed, what the record offers, and a first distribution

The register knows three kinds of land: wild land, hybrid land (heavily manicured megaparks and great gardens with their amenities) and the land where people live and work densely. On 19 November 2026 it opened three questions for the atlas ([register](../../decisions.md#the-atlas-rounds-frame-16-november-2026), line 1084): how the hybrid land works as gardens and resorts, whether the managed biomes lie with habitation or in the hybrid land, and how the three kinds are distributed over the Moon. This paper sets out what the register fixes (part A), what the record holds on garden and resort country (part B), where each managed living system falls (part C), a first distribution computed on the design surface in three variants (part D), and the decisions that remain (part E).

The short answer is that the shares are fixed and the map is open. With the food the register has chosen, the 4 billion people on the ground and everything that feeds and powers them take 5.4–14.6% of the land, and the hybrid fifth another 20%, so that about two thirds of the land carries no use under either wild share the register allows. A sorting by stated rules meets every share and every need at the compromise the register allows, 30% of the land wild beside the hybrid fifth; at the preferred half wild it meets them all except the room for the rain belt's fields at the food blend's high-cost end. The record gives gardens a developed design and gives resorts little beyond the climate's name.

Every figure carries its state: *(register)* a rule of the register, with its line in `research/decisions.md`; *(study)* a study's computed result or stated scenario, with its file; *(arithmetic)* arithmetic done here on products, kept in [results/land_use.json](results/land_use.json) by [land_use.py](land_use.py); *(proposal)* a proposal of this paper. The record was read at commit `615becb8` of `domain/geography`. Latitudes are today's unless called turned; the climate's bands follow the turned frame of the design surface, whose poles stand 6.7° from today's.

## A. What is fixed

**The register fixes the kinds, their shares, the direction of settlement and three places; it holds no map and one rule for a managed living system.** Hybrid land is defined by what it is for, leisure in great gardens kept by machines. The wild is at least half of the land and of the dark sky by preference, or 30% beside another 20% of hybrid land. People on the ground live mostly in the drier highlands, in the summit metropolis above all, with homestead communities among the farms of the middle latitudes. The places are Jura Island, laid out mostly hybrid; the summit metropolis; and the crater garden towns, whose form and number are adopted and whose sites are unchosen. The one rule for a managed living system puts much of the satoyama-like activity in the hybrid land.

### The rows

The standing is the register's own word for each row. Rows marked as working rules, directions or recorded are the register's rows in force; two rows are open and one is provisional.

| Date | The rule, in one sentence | Standing | Line |
|---|---|---|---|
| 25 Sep | The Moon has a mix of climates and biomes, and its hottest days are pleasantly warm. | Stands as an objective | 64 |
| 27 Sep | The Moon's largest sky port stands on the far side's summit (5.4° N, 158.6° W), commercial, with no homes and no industry. | Stands | 168 |
| 27 Sep | The tower's open-air parks follow the air: for everyone to 3 km above the summit, for residents and acclimatised visitors to 10 km, a few high gardens and winter parks to 15 km. | Stands | 172 |
| 27 Sep | A metropolis of about 100 million people at about 40,000 per km² surrounds the port within about 32 km. | Stands as a rough planning scale | 174 |
| 27 Sep | The metropolis uses the third dimension, with glider and parkour zones designed in. | Stands | 175 |
| 27 Sep | Regional agriculture and foraging feed the metropolis. | Planning assumption | 178 |
| 8 Oct | The ring fleet avoids lighting the Moon's night; light over cities may be acceptable as fine spotlights. | Stands as a direction | 130 |
| 9 Oct | The population's distribution is reopened and the inhabited volume evaluated, aerial settlement included, before an allocation is chosen. | Research direction | 151, 220 |
| 9 Oct | Much floating habitation roams. | Working direction | 152 |
| 9 Oct | The aerial biomes are considered almost first. | Research direction | 154 |
| 11 Oct | The Moon carries 4–8 billion people with plenty of it kept wild, the aerial life above all. | Working range | 277 |
| 11 Oct | 2–4 billion of them live aloft. | Working case | 278 |
| 11 Oct | On land people live mostly in the drier highlands, the summit metropolis above all, with many communities scattered elsewhere, and cities can be taller. | Direction | 279 |
| 11 Oct | The food system counts present and emerging agriculture, hydroponics, aquaponics and the fog deserts. | Research direction, studied | 281 |
| 11 Oct | The large lava tubes are evaluated as first habitats and lasting cities; the undersea cities are connected with them. | Research direction | 282, 284 |
| 11 Oct | The division of labour and the synergies among places are evaluated, so that a region's weakness becomes a strength. | Research direction, studied | 283 |
| 11 Oct | The distribution lies closest to C: 8 billion on the Moon, 4 on the ground and 4 aloft; homestead communities live among the farms of the middle latitudes, about 30–49°, and the rest of that belt's people go where they are better sustained. | Built as C′ | 285 |
| 11 Oct | At least half of the land and of the dark sky kept wild is preferred; it can come down to 30% wild with another 20% hybrid. | Working rule | 286 |
| 11 Oct | Hybrid land is heavily manicured megaparks and great gardens, food gardens among them, for leisure and entertainment, with small shops, restrooms, boutiques and places to rest; residents on the ground and aloft and visitors use them, and heavy automation keeps them. | Working rule | 287 |
| 11 Oct | Every kind of land is heavily managed, the wild for inspection, evaluation and preservation; wild land need not serve regular human use. | Direction | 288 |
| 11 Oct | Organisms differ by region across a connected world, and the parks and gardens stand apart from the wild by their intensive manicuring. | Direction for ecology | 296 |
| 12 Oct | Parks and gardens are built as far as possible round the lower water: reservoirs and water features in the fog deserts, dry-adapted planting in the highlands above all. | Direction | 289 |
| 12 Oct | Gardens are designed for richness within dry land: water placed sparingly where it shows, dry-adapted and silver planting, dryland food gardens, the lunar cycle as the show, shade fitted to each slope. | Direction | 290 |
| 12 Oct | The megaparks and great gardens are almost purely automated, their amenities and food-service shops included. | Direction | 302 |
| 12 Oct | The network's fixed posts watch the edges of the wild; fibre and small towers serve settled ground and the parks. | Design direction | 305 |
| 12 Oct | Food is a blend of three mixes (40% power-heavy, 50% power-made, 10% lamp-lit), with a little heritage meat from the homesteads and hobby farms. | The chosen mix | 307 |
| 12 Oct | The comparison run keeps the parks reading of hybrid land, in place of a reading as working countryside open to grazing and the homesteads' homes. | Stands | 308 |
| 12 Oct | The wild is monitored and maintained by machines, which also keep wildlife and habitation apart. | Direction | 314 |
| 12 Oct | The shepherds' signals stay out of the parks' and the night's show. | Direction | 317 |
| 12 Oct | What lives wild follows general constraints on an open world: the biomes are shaped by climate and the parks by manicuring. | Answers a question in general terms | 319 |
| 12 Oct | Beamed power reaches closed receivers at the ground's loads, with planting beneath the mesh. | Chosen | 330 |
| 12 Oct | Heritage meat is 1–2 kg a person a year. | Working range | 332 |
| 13 Oct | How long the terrain lasts and what renews it is worked out: returned sediment as blocks and terraces for built land, sediment carried back uphill for the wild, dredged deltas and living builders; only the heritage sites are kept so far. | Done as a study | 358 |
| 15 Oct | Jura Island is laid out mostly hybrid: its bay-facing slopes and shores a great garden, its high arc, outer coasts and islets wild, one harbour town on the inner bay. | Direction for the island's design | 378 |
| 17 Oct | The target climate is weighed as a world for recreation, with relief aloft and along its long coasts; the middle target, the resort, is in keeping with the hybrid fifth. | Answered and named | 425, 426 |
| 17 Oct | The low-gravity architecture suits the semi-arid parks, gardens and their amenities better than a city. | Provisional | 448 |
| 19 Oct | The resort at 301 K is the design climate (since 18 November estimated at 300.31 K on the design surface, line 1071). | Set | 492 |
| 19 Oct | Open water is nutrient-poor, with productive habitats at deltas, shelves, lake margins and in farms. | Set | 493 |
| 31 Oct | The hybrid regions carry considerable weight in the economy, with travel among the three worlds. | Done as a study | 598 |
| 31 Oct | The activities that suit the Open Moon alone or best are identified. | Done as a study | 605 |
| 3 Nov | The parks design's premise of 93% humid desert nights is superseded; its revision waits. | Done as a note | 693 |
| 12 Nov | Crater garden towns are 500–1,000 m rim to rim, terraced round a lake, of 900–15,000 residents; the 750 m example's estate is 1.8 km². | Adopted | 883, 884 |
| 12 Nov | The eventual buildout takes the best 2.5% of all catalogue craters 200–1,000 m across, 0.1% of them in the vacuum period, whose communities later stand open. | Adopted | 892–895 |
| 12 Nov | Four lava-tube cities are studied with local intensive agriculture and aquaculture. | Study programme | 902, 910 |
| 16–17 Nov | A district is a whole and meaningful place; a coastal sea is planned with the land that faces it, as satoumi needs. | Recorded | 1020 |
| 16 Nov | The heritage register's written evaluations are reviewed before districts are sited near its places. | Recorded | 1022 |
| 16 Nov | The atlas looks for satoyama-like mosaics of managed woodland, fields, ponds, channels, wetlands and settlements, and for satoumi on the salt sea's shores. | Recorded | 1028, 1031 |
| 16 Nov | The atlas estimates the riparian share of the land and maps the falls where they meet the satoyama candidates. | Recorded | 1025, 1032 |
| 16 Nov | The mosaics' criteria are set now and their candidates mapped on the design surface. | Done as a first map | 1041 |
| 16 Nov | Their methods limit methane, and methane-eating microbes are boosted at the sources. | Recorded | 1042, 1043 |
| 16 Nov | Much of the satoyama-like activity lies in the hybrid land, managed to limit methane and to draw lasting output. | Recorded | 1044 |
| 16 Nov | The hybrid fifth may grow a little where that serves sustainability, management and people. | Open | 1045 |
| 18 Nov | Each region holds its coastal seas within 60 km as its waters. | Recorded | 1074 |
| 19 Nov | The regions and districts are drawn again, their fact sheets widened with living cover and heritage. | Recorded | 1083 |
| 19 Nov | How hybrid land works as gardens and resorts, where the managed biomes lie, and the general distribution are worked out. | Open | 1084 |

### What the studies made of the rules

The food and distributions studies applied the wild share on 11 and 12 October, on a land of 23.4 million km², and their figures are the record's working numbers *(study)*.

- **The hybrid fifth** is 4.68 million km² of parkland, 585 m² for each of 8 billion residents, closed to farms, grazing, the homesteads' homes and towns ([distributions](../distributions/README.md), "Under the author's wild share"; [food systems](../food_systems/README.md)).
- **The people of C′**: the metropolis's 0.1 billion on 2,464 km²; tower districts of 0.77–1.45 billion on 8,950–16,900 km²; highland towns of 1.15–1.93 billion on 115,000–193,000 km²; scattered communities of 1.0–1.6 billion on 100,000–571,000 km²; in all 264,000–736,000 km² settled. Of the middle latitudes' lowland people 0.24–0.40 billion move to the high-latitude lowland, and 0.6–12 million stay as homesteaders ([distributions](../distributions/README.md), "People and land", "With the author's food, the blend").
- **Their food**: the blend's fields take 0.68–1.86 million km² at 8 billion under a kept half or the compromise, and 0.68–2.17 million km² under half wild with a hybrid fifth; heritage meat takes 298,000–595,000 km² of pasture at 1–2 kg, 21–43% of the homestead belt; the food gardens' harvest is outside the accounts ([food systems](../food_systems/README.md), "The author's blend").
- **All together**, farms, pasture, settlement and the homesteads' homes take 5.3–14% of the land, and 73–89% of the land outside a kept 30% and the parks stays wild ([distributions](../distributions/README.md)).
- **The dark sky**: towns, cities and sky towns put 0.49–17.2% of the land under a light-polluted sky, so that at least 62.8% keeps a dark sky even with every park lit ([distributions](../distributions/README.md), "Light").
- **The high ground**: 3.90 million km² stands 4 km or more above the sea, 1.35 million km² of it poleward of 30°, where both climate models find the most comfort ([land cities](../land_cities/README.md)).
- **The parks' use and upkeep**: 52 visits a resident a year, 416 billion visits, 243 per km² a day; 2.86 billion visitor trips a year to a hybrid region and 206 million visitors present on the Moon on average; operating and renewal budgets of 0.70–21.1 trillion U a year across three assumed levels of provision ([amenities](../world_economy/amenities.md), all scenarios). At Earth's staffing the parks would need 50 million to 1.5 billion gardeners, which is why machines keep them ([habitation](../../../habitation/README.md#parks-and-gardens)).
- **The mosaics**: 56,693 candidate side valleys over 1.77 million km², with 890,000 km² of woodland ground, 527,000 km² of valley floor and room for 215,000 km² of flooded fields, whose methane of 0.94–1.59 Tg a year stands against a room of 0.35–2.5 Tg ([satoyama.md](satoyama.md); [results/satoyama_final.json](results/satoyama_final.json)).

**On the present land** *(arithmetic)*. The design surface holds 22.13 million km² of land at 16 px/deg, 5.4% less than October's, because its lakes hold more. A fifth of it is 4.43 million km², 0.25 million km² less than the 4.68 million the register's rows and the studies carry, and 553 m² a resident at 8 billion. The people's and the food's areas are carried in km² as the studies give them, so their share of the land rises from 5.3–13.6% to 5.6–14.4%.

### Where the record disagrees with itself

1. **Hybrid land and working land.** The October studies keep the hybrid fifth clear of farms, grazing, homes and towns, and the row of 12 October (line 308) sets aside a reading of hybrid land as working countryside. The row of 16 November (line 1044) puts much of the satoyama-like activity, whose mosaics hold fields, ponds and settlements, in the hybrid land. The later row governs the mosaics; it leaves their villages and the ordinary fields unstated.
2. **Dry gardens and wet mosaics.** The garden rows build the parks round the lower water in dry land (lines 289, 290). The mosaics need 800–3,000 mm of rain a year. Both stand, so the hybrid land holds both, and inside one fifth the mosaics' candidates take 40% of it *(arithmetic)*.
3. **Zone by zone.** The studies hold the wild and hybrid shares in every zone. Between the turned equator and 15° the candidate valleys cover 105% of that band's fifth *(arithmetic)*, so the hybrid share cannot be held band by band with the mosaics inside it.
4. **The drier highlands and the summit.** The settlement row places people in the drier highlands, the summit metropolis above all (line 279). The summit's massif is the wettest high ground: 13.1 mm of rain a day within 32 km of the tower in the distributions study, and on the design climate a wet-class cell with 5,140 mm a year *(arithmetic on the moisture map)*.
5. **How much high ground is dry.** The studies call the 1.35 million km² of high ground poleward of 30° the dry high ground. On the design surface that ground measures 1.43 million km², of which 0.64 million km² lies in the moisture map's dry classes *(arithmetic)*.
6. **The fog desert.** The parks design and the studies place 1.65 million km² of parks in a fog desert with 93% humid nights. The design climate's desert nights run 57–60% humid and fog forms over 0.4% of the desert land (line 693); the zone's name outlived its premise.
7. **The shared cases.** [shared/scenarios/population.json](../../../shared/scenarios/population.json) holds five comparison cases of 9 October and no C′; the distribution the register built on 11 October lives in the distributions study alone.

### The summit metropolis and the crater garden towns

The metropolis is habitation by the register's own words: about 100 million people at about 40,000 per km² within about 32 km of the port *(register, line 174)*. On the design surface that circle holds 2,464 km² of land, the rest of its 3,217 km² being crater lakes *(arithmetic)*, which matches the area the studies give the metropolis. Its parks, the tower's among them, are the parks a dense city needs and belong to it. What surrounds it is unassigned: the metropolis brief leaves the lower rise and Korolev's shore for farms, foraging and wild land as a first estimate, and feeding it from the region is a planning assumption *(line 178)*.

A crater garden town is habitation by what it is: homes, work, shops and community rooms for 900–15,000 people on terraces round a lake, with outer gardens that bring the 750 m example's estate to 1.8 km² *(register, lines 883, 884)*. The record does not say which kind of land the towns stand in. Their brief names watered uplands and the transition between wetter and drier regions as the first settings to examine, which is also the mosaics' margin. Their number is a fraction of a catalogue count that lies on another branch (`domain/provisioning`, `research/studies/crater_gardens/`), so no area can be given here; the superseded scenario of 20,000–100,000 towns would put 36,000–180,000 km² under estates *(arithmetic on line 885's figures)*.

## B. What hybrid land does

**The record supports five kinds of garden country in some depth and treats resorts mostly by name.** The gardens have a design direction, a worked island, a candidate map and an account of their use. The resort is the name the register gave the climate; resort country as a place has the opening question of 19 November behind it and scattered physical facts about the coasts and the water.

| Kind | What people do there | What it needs | What on the design surface suits it | Source and depth |
|---|---|---|---|---|
| Dry gardens of the highlands | Walk, rest, meet, eat at the shops, come back to see the garden change through the lunar cycle | Little water: stone and gravel rivers, a few small shaded features, mulches, the garden's own treated water; shade structures; machines for upkeep | Dry classes on high ground, in clear, drier air, beside the highland cities; dry-class ground 4 km or more above the common level is 0.87 million km² *(arithmetic)* | Register lines 289, 290; [habitation](../../../habitation/README.md#parks-and-gardens). A design direction; its water, energy and night light are unsized |
| Water gardens of the dry lowlands | The same, round still pools and rills | Reservoirs; water led or lifted. A fully watered garden in the desert takes 0.71 million m³ a year per km², 32–66 MW per 1,000 km² when lifted from the groundwater *(study)* | Arid and semi-arid lowland with fresh water near; 2.24 million km² of dry land lies within 25 km of a sea or great lake *(arithmetic)* | Register lines 289, 290, 693; [resort at 301 K](../resort_301k/README.md). The premise of humid nights is corrected and the revision waits |
| Garden shores and islands, after Jura | Arrive by the harbour town or by air; gardens on the slopes and shores that face sheltered water, under a wild crest | Water drawn from the fresh sea, 0.11 kWh for each m³ lifted 250 m before pump losses; shore belts where plants reach the water table 15–25 m above the sea; a gateway town | Dry coasts and islands with a bay and high ground behind; the near side's coasts keep Earth in the sky; 0.67 million km² of dry land lies within 5 km of a sea or great lake *(arithmetic)* | Register lines 378, 379; [the island study](../jura_island/README.md). One worked layout; its climate box found the island dry, with hot afternoons and about a third of the cycle comfortable |
| Satoyama-like mosaics in garden form | Tend and harvest with the machines, walk the valley, eat what grows there | Rain of 800–3,000 mm a year; a side valley of 10–100 km² with gentle floor and steep or high ground; water from its own runoff, stored and reused down the slope; drained paddies, composted straw, clean ponds and boosted methane-eating microbes | Humid and wet valleys within 30° of the turned equator: 1.77 million km² of candidates, 0.98 million of it within 15° *(study, arithmetic)* | Register lines 1028, 1041–1044; [satoyama.md](satoyama.md). Criteria and a candidate map; what grows and how it is run belong to the ecology lane |
| Satoumi shores of the salt sea | The record does not say | Shallows for eelgrass and shellfish; the coast and its land planned together | The South Pole–Aitken Sea's shore, where 2,609 km² lies less than 10 m deep along 37,800 km *(study)* | Register line 1031; [districts.md](districts.md), rule 12. Given in passing |
| Food gardens | Grow and eat within the parks | Unsized; their harvest is outside the food accounts | Unsaid | Register lines 287, 290; [amenities](../world_economy/amenities.md). Given in passing |

**What all of them share.** Every garden is a monthly performance: bloom at each dawn and dusk, silver foliage, living light and night scent through the 354-hour night *(register, line 290)*. Small shops, restrooms, boutiques, food service and places to rest stand among them *(lines 287, 302)*, and machines do the keeping, much of it at night by infrared so the living light keeps the dark *(habitation)*. On the present fifth the amenities study's 416 billion resident visits a year come to 257 a km² a day *(arithmetic)*. Water is the largest unsized need: if the 2.62 million km² of dry gardens in the first variant below were watered as the resort study waters a desert garden, they would take 1,860 km³ a year, 24% of the 7,640 km³ the rivers bring the seas *(arithmetic)*, which is the measure of why the register builds them round the lower water.

**The night as part of the place.** The moisture map sorts the land by its darkest hour: 10.7% of it has dark nights, below 0.01 lux, on the far side's low latitudes; 31.2% dim nights; 38.8% nights lit by the Earth; and 19.3% nights that keep twilight *(arithmetic on the map)*. The record uses the night three ways: living light and night scent in every garden (register), the dark belt's dark nights as a strength of its region ([distributions](../distributions/README.md), synergies), and night gardens, performances and observatories as activities that suit the Moon ([distinctive activities](../world_economy/distinctive_activities.md), an idea).

### Resorts

The register names the resort in three rows, all about the climate: a warm Moon whose relief lies aloft, in the summer heat some people seek and along its long coasts, perhaps chiefly a world for recreation, in keeping with the hybrid fifth *(lines 425, 426, 492)*. The fourth is the question of 19 November itself. What the record holds beside them:

- **Visitors.** In the amenities study's connected case 3.58 billion arrivals a year stay 21 nights, 80% of them using a hybrid region, and 206 million visitors are on the Moon on average *(study, scenario)*. The flow is unsized against transport and arrival ports.
- **The coasts.** The land cells of the design surface carry 56,400 km of coast: 15,100 km with a monthly tide under 2 m, 30,000 km with 2–4 m and 11,300 km above 4 m ([results/moisture_map_final.json](results/moisture_map_final.json)). The tide rises and falls once a month, and on beaches of 1 in 50 to 1 in 100 it moves the waterline 140–280 m ([ways of living](../ways_of_living/README.md)). Waves travel at 0.41 times the speed of Earth's at the same wavelength *(register, line 71)*.
- **The water's warmth.** On the corrected design climate the seas' and large lakes' upper 50 m averages 29.8 °C, 23.7–32.6 °C by cell, and still water 2 m deep reaches 36 °C or more along 54% of the sea coast in the hottest hundredth of hours *(register, lines 866, 869)*. The water is warm everywhere and at times too warm in the shallows; no fresh water falls below about 15 °C.
- **The air.** The classes' mean air runs from 25.0 °C on wet land to 27.5 °C on hyper-arid land, and the warmest three-day mean is 34–35 °C in every class *(arithmetic on the map)*. Both climate models put the most comfortable land poleward of 30° and on high ground.
- **Things to do in and by the water.** Deep viewing rooms and tall aquaria, slower waves and wave-riding made for them, tall fountains ([distinctive activities](../world_economy/distinctive_activities.md), "Water as an activity space"; ideas of a study).
- **Somewhere to stay.** The register's amenities end at places to rest. Visitor accommodation appears once, in a table cell of the distinctive-activities study.

The record holds no beach as a place: the sea-appearance study's coast is a placeholder without beaches or surf, the coast-water study leaves beaches and sediment out, and nothing treats bathing water, lodging or the siting of a resort. The distribution in part D shows the consequence. Its gardens follow the register into the drylands, so the humid and wet half of the land, with 1.63 million km² within 5 km of a sea or great lake against the drylands' 0.67 million km² *(arithmetic)*, holds hybrid land only where a candidate valley lies. Whether a warm, wet coast can be resort country is decision 7 of part E.

**Given only in passing**, and used here as ideas: canopy walks, flying clearings and water terraces between garden levels; aerial gardens on aerophyte hosts; long-light and long-night festivals; large testing campuses ([distinctive activities](../world_economy/distinctive_activities.md)); the falls and river corridors of the filled seas' overflows, 850 steep reaches on the design surface ([results/falls.json](results/falls.json)), which the register asks to be mapped where they meet the mosaics *(lines 1015, 1032)*; and the settings the conservation study lists for communities, among them islands, near-side coasts under the Earth and far-side waters.

## C. The managed biomes

**The register's definitions settle the gardens, the grazing and the homesteads, half-settle the mosaics, and leave the fields unstated.** The test the rows share is purpose: hybrid land is kept for leisure and what large public gardens are for, with food gardens inside it; land worked to feed or house people is outside it.

| Managed living system | Area the studies give it | Kind by the register's own definitions | State |
|---|---|---|---|
| Megaparks and great gardens, with their pools, rills, reservoirs and food gardens | The hybrid fifth: 4.68 million km² on October's land, 4.43 million km² now *(arithmetic)* | Hybrid (lines 287, 289, 290) | Settled |
| Satoyama-like mosaics: woodland, gardens, rain-fed fields, flooded fields, ponds, channels, wetlands | 1.77 million km² of candidate valleys, with room for 215,000 km² of flooded fields, 106,000 km² at Sri Lanka's share *(study)* | Much of it hybrid (line 1044); the rest unstated | Half settled |
| The mosaics' settlements | Unsized | Hybrid land holds amenities; homes lie outside it in the studies' reading and by line 308 | Open |
| Satoumi on the salt sea | 2,609 km² of water under 10 m deep *(study)* | Water, planned with the land that faces it; the land's kind unstated | Open |
| The blend's fields | 0.68–1.86 million km² at 8 billion, 0.68–2.17 million under half wild *(study)* | The studies keep fields out of the wild and the hybrid land; no row names their kind | Open |
| Vertical farms and power-made food | Buildings of 1.9–22 Gt, no land figure *(study)* | With habitation, by what they are | Unstated, uncontested |
| Heritage pasture and farms | 298,000–595,000 km² at 1–2 kg *(study)* | Outside hybrid land (line 308), in the homestead belt (lines 285, 307) | Settled |
| Homesteads, hobby farms and their homes | 0.6–12 million people among the heritage farms *(study)* | Outside hybrid land (line 308) | Settled |
| Grazing on dry land | An ecology allowance of 0–1,000 m² a person; none in the blend *(study)* | Outside hybrid land (line 308) | An idea outside the chosen food |
| Receivers' ground with planting beneath | 36–59% of the settled land, 95,000–434,000 km² *(arithmetic)* | At the ground's loads, with habitation (line 330) | Settled as to place |
| Productive waters and aquatic farms | Water; the sea-supported mix lies outside the blend *(study)* | No land kind | Outside the land's kinds |
| Kept lakes and landforms; dryland rivers | Upkeep of 0.2–11 GW; 25 gravity routes, unchosen *(study)* | Management of every kind of land (lines 288, 358); the rivers unstated | Open for the rivers |
| A crater garden town's terraces and outer gardens | 1.8 km² an estate; the number rests on a count elsewhere | The town is habitation; its setting unstated | Open |
| The metropolis's parks | Ring 0's park about 3 km²; the tower's parks to 15 km | The city's own (line 172) | Settled |

### The readings that remain, and what each does to the shares

Three readings are possible for the mosaics and the fields *(arithmetic)*.

1. **The mosaics in the hybrid land, the working estates with habitation.** The candidate valleys take 1.77 million km² of the fifth and gardens the other 2.62 million km²; fields, pasture and homesteads lie in the open land. Wild 30% or 50%, hybrid 20%, open 50% or 30%. This is variant A or B of part D.
2. **Every managed biome with habitation.** The hybrid fifth is gardens alone, 4.39 million km², all of it in the drylands since the garden rows point there, and the mosaics are the first of the rain belt's farm country. The shares are those of reading 1; the humid and wet 44% of the land holds no hybrid land at all. This is variant C. It sets aside line 1044.
3. **Every managed biome in the hybrid land.** Mosaics, fields and pasture together are 2.75–4.23 million km², 12.4–19.1% of the land. Inside one fifth they leave gardens 0.2–1.7 million km²; beside a full fifth of gardens they raise the hybrid land to 32–39% and leave habitation its towns and receivers, 1.6–5.3% of the land. This reading sets aside line 308 for the pasture.

Between readings 1 and 3 lies the growth the register left open: the mosaics beside a full fifth of gardens, hybrid land at 28.0%, 6.20 million km².

**Recommendation** *(proposal)*: reading 1, sorted by purpose and access. A managed living system belongs to the hybrid land where people walk in it as a garden and machines keep it: the food gardens, and the mosaics in their garden form with their woodland, ponds and channels. It belongs with habitation where it is an estate worked to feed or house people: the blend's fields and food plants, the heritage pasture, the homesteads. The mosaics' villages then stand as habitation at the mosaics' edge, as Jura's harbour town stands at its garden's. The reason is that this reading keeps every row in force, the garden rows, line 308 and line 1044 together, and that the alternative readings each set one aside.

## D. A first general distribution

**Under either wild share the land has room to spare; what the variants decide is which unused land is promised to the wild, and where the mosaics count.** Variant A, 30% wild beside the hybrid fifth, meets the compromise's shares and every need. Variant B, half wild beside the hybrid fifth, meets them all except the rain belt's room at the food blend's high-cost end, where it falls 0.96 million km² short. Variant C, the mosaics placed with habitation, meets them all and leaves the humid and wet land without hybrid land.

### The rules

| | Rule | State | Source |
|---|---|---|---|
| 1 | The land within 32 km of the summit port is the metropolis: habitation. | Register | Line 174 |
| 2 | Jura Island takes its layout: the high arc, the outer coasts and the islets wild, the rest hybrid, the harbour town too small to draw. | Register; the high arc as land more than 1 km above its sea and the outer coast as land within 5 km of the sea outside the bay are proposals | Line 378 |
| 3 | The satoyama-like candidate valleys are hybrid land. | Register | Line 1044 |
| 4 | Gardens fill the rest of the hybrid share from the drylands, the land nearest a sea or great lake first. | Register for the dry land; proposal for the order | Lines 287–290; the order follows Jura's garden shores and the cost of lifting water |
| 5 | In every moisture class, and within it in every night regime, the wild keeps its share, taken from the land ranked last for dense living. | Register for the share of the land and of the dark sky; study for holding it zone by zone; proposal for the night regimes and the order | Line 286; the food and distributions studies |
| 6 | The land left is open to habitation and its farms, ranked: dry highlands poleward of 30°, the high-latitude lowland poleward of 48.6°, the homestead belt at 30–48.6°, the rain belt within 30°, dry highlands within 30°, other land. | Register for the direction; study for the bands | Lines 279, 285; C′ and the land-cities study |
| 7 | Of the open land, C′ at 8 billion with the blend fills the best-ranked land of each band up to the band's need; the rest is open land beyond that need. | Study figures; arithmetic for the needs by band | Distributions; food systems |

The rules run on the design surface's 10.8 million land cells of about 1.9 km. Each cell reads its moisture class, night regime and rain from the climate model's cell that holds it in the turned frame; 15.7% of the land lies in cells the model holds as water and reads the nearest land cell, as the moisture map assigns it. The classes' shares of the land agree with the map's own count to 0.04 points. High ground stands 4 km or more above the common level; the drylands are the map's four dry classes; a great lake is one of 1,000 km² or more.

**Reading the maps.** On the cell maps the shores, the garden rims and the valleys are the land's own: the rims follow real shores at a set distance and the candidate valleys are drainage units of 10–100 km². The straight edges and rectangular steps between wild and open land are the climate model's cells, about 170 km across and slightly oblique to today's parallels because they lie in the turned frame, together with the limits of turned latitude that rule 6 uses. The fine mixture of kinds in the lake country is the candidate valleys set among cells whose rank is the same across a whole climate cell. A second set of maps therefore draws each of the redrawn atlas's land districts, whose edges follow divides, rims and shores, in the kind that holds most of its land, striped with its second kind where none holds more than half.

### What the surface offers

| The design surface's land *(arithmetic)* | km² | Share |
|---|---:|---:|
| All land | 22.13 million | 100% |
| Hyper-arid, arid, semi-arid, dry sub-humid | 0.64, 4.68, 5.30, 1.81 million | 2.9%, 21.1%, 24.0%, 8.2% |
| Humid, wet | 3.23, 6.47 million | 14.6%, 29.2% |
| Dry highlands poleward of 30° (dry classes, 4 km or more up) | 0.64 million | 2.9% |
| High-latitude lowland poleward of 48.6° | 6.93 million | 31.3% |
| Homestead belt: 30–48.6°, dry sub-humid to wet, below 4 km | 2.40 million | 10.9% |
| Rain belt: within 30°, humid and wet | 6.75 million | 30.5% |
| Dry highlands within 30° | 0.24 million | 1.1% |
| Other land | 5.17 million | 23.4% |
| Satoyama-like candidate valleys | 1.77 million | 8.0%, 40% of a fifth |
| Dry land within 5, 25, 50 km of a sea or great lake | 0.67, 2.24, 3.39 million | 3.0%, 10.1%, 15.3% |

### The variants

The three differ only in what the register leaves open. **A** takes the compromise, 30% wild beside the hybrid fifth, with the mosaics in the hybrid land. **B** takes the preferred half wild beside the same hybrid fifth, as the register's comparison run does. **C** is A with the managed biomes placed with habitation. The record also carries **A+**, A with a full fifth of gardens beside the mosaics, as the measure of the fifth grown.

| Of the land *(arithmetic)* | A | B | C | A+ |
|---|---:|---:|---:|---:|
| Wild | 30.0% | 50.0% | 30.0% | 30.0% |
| Hybrid | 20.0% | 20.0% | 20.0% | 28.0% |
| Open to habitation and its farms | 50.0% | 30.0% | 50.0% | 42.0% |
| of which C′ fills at its high end | 14.6% | 11.7% | 14.6% | 14.6% |
| Hybrid land in candidate valleys; in gardens, million km² | 1.77; 2.62 | 1.77; 2.62 | 0; 4.39 | 1.77; 4.39 |
| How far inland the gardens reach | 37 km | 39 km | 101 km | 101 km |
| Hybrid land a resident at 8 billion | 553 m² | 553 m² | 553 m² | 775 m² |
| Hybrid share of the hyper-arid, arid, semi-arid, dry sub-humid classes | 14, 16, 25, 37% | 15, 16, 26, 34% | 32, 30, 38, 45% | 32, 30, 39, 50% |
| Hybrid share of the humid, wet classes | 8, 21% | 8, 21% | 0, 0% | 8, 21% |

By turned latitude the kinds lie as follows, as wild / hybrid / open and, in brackets, the part of the open land that C′ fills *(arithmetic)*:

| Turned latitude | A | B | C |
|---|---|---|---|
| 0–15° | 28 / 33 / 39 (31) | 47 / 32 / 20 (20) | 25 / 18 / 57 (27) |
| 15–30° | 35 / 28 / 37 (16) | 59 / 28 / 13 (13) | 31 / 19 / 49 (17) |
| 30–45° | 40 / 20 / 40 (12) | 49 / 20 / 30 (13) | 39 / 25 / 37 (16) |
| 45–60° | 24 / 13 / 63 (9) | 42 / 14 / 44 (8) | 22 / 26 / 52 (8) |
| 60–90° | 23 / 5 / 72 (4) | 52 / 5 / 43 (4) | 33 / 12 / 56 (5) |

In A the equatorial band is a third garden and mosaic country, and the high latitudes are mostly open land that nobody needs: the polar lowland is C′'s first destination for those who move, so the rules rank it high and the wild takes its share elsewhere.

### The checks

| What is fixed | A | B | C |
|---|---|---|---|
| The wild share of the land, in every moisture class *(register, line 286)* | 30% in each: the compromise, below the preferred half | 50% in each | 30% in each |
| The wild share of the dark sky, as the land under an unlit sky with every park lit *(the distributions study's reading)* | 61.8–79.5% | 61.8–79.5% | 61.8–79.5% |
| The wild share of the land with the darkest nights *(proposal for the map's reading)* | 30%; fields and settlement fill 44% of that land | 50%; 31% | 30%; 25% |
| The hybrid share *(register, lines 286, 287)* | 20.0% | 20.0% | 20.0%, none of it in the humid and wet classes |
| Room for C′'s people and food, 5.4–14.6% of the land, 5.9–16.6% with the receivers | Fits in every band | Fits at the low end; the rain belt is 0.96 million km² short at the high end | Fits in every band |
| Jura's layout *(register, line 378)* | Kept | Kept | Kept |

The room by band, as the open land in the band against C′'s need there, in thousand km² *(arithmetic; the needs are bounds built from the studies' ranges)*:

| Band | Need, low to high | Open in A | Open in B | Open in C |
|---|---:|---:|---:|---:|
| Dry highlands poleward of 30°: tower districts and highland towns | 124–210 | 537 | 476 | 398 |
| High-latitude lowland: its communities and those who leave the middle latitudes | 53–300 | 5,140 | 3,259 | 4,045 |
| Homestead belt: heritage pasture | 298–595 | 1,783 | 1,380 | 1,669 |
| Rain belt: the blend's fields and its communities | 728–2,131 (2,441 under half wild) | 3,150 | 1,477 | 4,923 |

B's shortfall has a plain cause. Half of the wet class is wild, a fifth of it is mosaic, and the blend's high-cost end under half wild asks 2.17 million km² of fields where the compromise asks 1.86 million. The low-cost end needs 0.73 million km² and fits twice over. Three things would close the gap, each a decision in part E: counting the mosaics' harvest towards the diet, since their valleys could carry 215,000 km² of flooded fields and 527,000 km² of valley floor; letting fields use the 0.79 million km² the homestead belt keeps open beyond its pasture; or holding the wild share over the Moon in place of class by class.

Four further readings of the result *(arithmetic)*:

- **The unused land.** After the hybrid fifth and C′'s needs with the receivers, 63–74% of the land carries no use in A and 62–74% in B. A promises 30 points of it to the wild and leaves 33–44% of the land open and unneeded; B promises 50 and leaves 12–24%.
- **The shores.** Rule 4 takes 97.7% of the drylands' land within 5 km of a sea or great lake for gardens, so the only wild dry shore left is Jura's outer coast. The humid and wet shores are 32% wild, 6% hybrid and 62% open in A.
- **Reach.** Half of A's hybrid land lies within a sky boat's 66.7 km trip of land that C′ fills, 78% within the sky ferry's 180 km route, and 11% farther than 300 km.
- **By whole districts.** If each district took the one kind that holds most of it, A's land would be 31% wild, 15% hybrid and 54% open, in 153, 100 and 296 districts, with no kind above half in 73 of them. The hybrid share falls because the garden rims are strips that seldom fill a district; 68% of the land lies in its district's kind.

### The regions

Shares of each region's land as wild / hybrid / open, with the part of the open land that C′ fills in brackets, in percent *(arithmetic)*. Regions carry the redrawn atlas's codes and anchors; their names are a draft and are left out. The label point is today's and the last figure is its turned latitude.

| Region | Anchor | Label point today; turned latitude | Land, km² | Main classes | A | B | C |
|---|---|---|---:|---|---|---|---|
| L01 | crater Peary | 90° N, 85° E; 83° N | 656,313 | semi-arid 55%, arid 38% | 0 / 6 / 94 (2) | 20 / 6 / 74 (2) | 0 / 11 / 89 (2) |
| L02 | crater Poncelet | 60° N, 157° W; 67° N | 821,393 | arid 50%, dry sub-humid 25% | 26 / 6 / 69 (18) | 50 / 6 / 44 (18) | 27 / 16 / 57 (19) |
| L03 | crater Poncelet | 66° N, 76° W; 67° N | 715,873 | arid 92%, hyper-arid 8% | 6 / 14 / 81 (6) | 24 / 14 / 62 (6) | 14 / 28 / 58 (6) |
| L04 | crater Baillaud | 64° N, 50° E; 57° N | 556,104 | arid 100% | 11 / 32 / 58 (9) | 10 / 33 / 57 (9) | 2 / 58 / 40 (9) |
| L05 | crater Poncelet | 59° N, 122° W; 65° N | 540,799 | semi-arid 57%, arid 43% | 59 / 1 / 39 (1) | 75 / 1 / 24 (1) | 66 / 4 / 30 (2) |
| L06 | lake at 53° N, 114° E | 57° N, 131° E; 57° N | 1,052,465 | semi-arid 51%, arid 24% | 32 / 15 / 54 (6) | 52 / 15 / 33 (6) | 31 / 30 / 39 (6) |
| L07 | crater Gerard | 38° N, 119° W; 43° N | 710,140 | humid 62%, dry sub-humid 32% | 54 / 2 / 44 (5) | 68 / 2 / 30 (4) | 57 / 2 / 40 (6) |
| L08 | ground at 30° N, 83° W | 33° N, 87° W; 37° N | 373,069 | humid 100% | 34 / 0 / 66 (8) | 35 / 0 / 65 (9) | 34 / 0 / 66 (9) |
| L09 (island group) | Montes Jura | 44° N, 40° W; 42° N | 67,678 | dry sub-humid 72%, semi-arid 28% | 48 / 52 / 0 (0) | 48 / 52 / 0 (0) | 48 / 52 / 0 (0) |
| L10 | Mare Tranquillitatis | 33° N, 76° E; 28° N | 776,965 | humid 44%, wet 43% | 19 / 40 / 41 (21) | 44 / 40 / 16 (11) | 17 / 16 / 67 (36) |
| L11 | lake at 16° N, 135° E | 32° N, 127° E; 32° N | 457,985 | wet 81%, humid 15% | 7 / 18 / 75 (58) | 21 / 18 / 61 (56) | 6 / 9 / 84 (38) |
| L12 | lake at 38° N, 170° E | 42° N, 155° E; 45° N | 191,821 | humid 84%, wet 16% | 29 / 33 / 38 (26) | 34 / 34 / 32 (23) | 24 / 40 / 37 (30) |
| L13 | lake at 14° N, 160° W | 28° N, 165° W; 34° N | 506,615 | wet 94%, humid 6% | 43 / 17 / 40 (34) | 49 / 17 / 34 (33) | 43 / 4 / 53 (36) |
| L14 | lake at 16° N, 90° W | 19° N, 102° W; 24° N | 471,047 | humid 79%, wet 14% | 67 / 3 / 30 (0) | 81 / 3 / 17 (0) | 67 / 4 / 29 (1) |
| L15 | ground at 9° N, 47° W | 21° N, 40° W; 19° N | 500,994 | semi-arid 55%, arid 28% | 45 / 45 / 10 (1) | 52 / 47 / 1 (0) | 17 / 75 / 7 (7) |
| L16 | lake at 1° N, 128° W | 8° N, 143° W; 14° N | 365,160 | wet 100% | 26 / 18 / 56 (45) | 37 / 18 / 45 (45) | 26 / 6 / 67 (33) |
| L17 | lake at 2° S, 84° W | 6° N, 78° W; 9° N | 436,177 | semi-arid 40%, humid 27% | 51 / 30 / 19 (0) | 70 / 28 / 2 (0) | 44 / 44 / 12 (5) |
| L18 | Sinus Medii | 12° N, 30° W; 9° N | 351,174 | wet 74%, humid 19% | 31 / 46 / 22 (21) | 47 / 47 / 6 (6) | 31 / 5 / 64 (43) |
| L19 | Sinus Amoris | 13° N, 43° E; 6° N | 236,272 | wet 91%, humid 9% | 38 / 40 / 22 (22) | 57 / 40 / 3 (3) | 37 / 2 / 61 (39) |
| L20 | lake at 16° N, 120° E | 10° N, 112° E; 8° N | 474,058 | wet 82%, humid 18% | 18 / 29 / 54 (35) | 50 / 27 / 23 (23) | 18 / 6 / 76 (23) |
| L21 | lake at 8° S, 177° W | 12° S, 170° W; 6° S | 403,517 | wet 89%, humid 11% | 21 / 17 / 62 (62) | 42 / 16 / 42 (42) | 21 / 4 / 74 (47) |
| L22 | lake at 2° S, 110° W | 1° S, 107° W; 4° N | 456,658 | wet 65%, humid 28% | 38 / 26 / 36 (16) | 63 / 22 / 15 (15) | 38 / 18 / 44 (11) |
| L23 | Mare Vaporum | 3° S, 12° E; 9° S | 615,594 | wet 72%, humid 27% | 24 / 44 / 32 (32) | 50 / 44 / 6 (6) | 22 / 10 / 68 (36) |
| L24 | lake at 6° S, 152° E | 9° S, 149° E; 7° S | 471,897 | wet 98%, humid 2% | 19 / 15 / 66 (62) | 41 / 15 / 45 (45) | 19 / 2 / 80 (25) |
| L25 | lake at 28° S, 141° W | 19° S, 150° W; 13° S | 350,343 | wet 77%, humid 19% | 22 / 31 / 47 (42) | 45 / 31 / 24 (24) | 22 / 17 / 62 (15) |
| L26 | Lacus Aestatis | 21° S, 67° W; 20° S | 837,478 | humid 33%, dry sub-humid 32% | 38 / 24 / 38 (10) | 64 / 24 / 11 (9) | 38 / 25 / 37 (14) |
| L27 | lake at 26° S, 31° E | 17° S, 0° E; 23° S | 1,044,000 | semi-arid 64%, humid 21% | 47 / 33 / 21 (4) | 54 / 34 / 12 (4) | 36 / 39 / 25 (11) |
| L28 | lake at 22° S, 52° E | 22° S, 55° E; 29° S | 411,047 | humid 52%, wet 34% | 6 / 39 / 56 (27) | 26 / 39 / 36 (34) | 6 / 19 / 75 (33) |
| L29 | lake at 17° S, 68° E | 19° S, 84° E; 24° S | 375,621 | humid 52%, wet 46% | 20 / 34 / 46 (21) | 36 / 33 / 31 (31) | 21 / 11 / 68 (30) |
| L30 | Lacus Solitudinis | 17° S, 117° E; 18° S | 303,277 | wet 71%, humid 26% | 28 / 31 / 41 (19) | 63 / 30 / 6 (6) | 28 / 14 / 58 (20) |
| L31 | lake at 35° S, 98° W | 33° S, 119° W; 27° S | 573,453 | humid 59%, wet 31% | 30 / 31 / 39 (19) | 54 / 30 / 16 (13) | 30 / 23 / 47 (21) |
| L32 | Lacus Timoris | 39° S, 44° W; 40° S | 565,955 | dry sub-humid 46%, semi-arid 40% | 55 / 7 / 37 (0) | 76 / 8 / 16 (0) | 50 / 17 / 33 (0) |
| L33 | lake at 34° S, 96° E | 40° S, 61° E; 46° S | 680,443 | semi-arid 94%, dry sub-humid 6% | 48 / 3 / 48 (4) | 56 / 4 / 40 (4) | 43 / 13 / 44 (5) |
| L34 | lake at 34° S, 96° E | 40° S, 122° E; 40° S | 394,118 | semi-arid 54%, dry sub-humid 46% | 19 / 34 / 47 (7) | 39 / 35 / 26 (6) | 12 / 50 / 39 (8) |
| L35 | lake at 20° S, 162° E | 35° S, 148° E; 32° S | 447,143 | humid 52%, semi-arid 46% | 15 / 50 / 35 (12) | 34 / 49 / 17 (14) | 13 / 58 / 28 (11) |
| L36 (island group) | island at 41° S, 173° E | 42° S, 173° E; 36° S | 3,953 | arid 100% | 0 / 100 / 0 (0) | 0 / 100 / 0 (0) | 0 / 100 / 0 (0) |
| L37 | crater Schiller H | 53° S, 89° W; 49° S | 1,021,768 | semi-arid 100% | 26 / 5 / 69 (4) | 58 / 5 / 36 (4) | 28 / 11 / 61 (5) |
| L38 | ground at 42° S, 27° E | 48° S, 20° E; 54° S | 830,861 | arid 39%, semi-arid 35% | 28 / 0 / 72 (0) | 64 / 0 / 35 (0) | 30 / 0 / 70 (0) |
| L39 | ground at 60° S, 64° W | 67° S, 0° W; 72° S | 538,570 | arid 77%, semi-arid 23% | 42 / 0 / 58 (0) | 95 / 0 / 5 (0) | 79 / 0 / 21 (0) |
| L40 | crater Chamberlin | 66° S, 77° E; 70° S | 879,140 | arid 80%, semi-arid 20% | 21 / 15 / 64 (8) | 50 / 16 / 34 (8) | 30 / 29 / 41 (5) |
| L41 | crater Cabeus | 90° S, 26° E; 83° S | 661,403 | arid 70%, wet 18% | 19 / 18 / 63 (12) | 26 / 19 / 55 (12) | 23 / 33 / 44 (11) |

In A the regions richest in hybrid land are the dry plain above the near-side sea (L15, 45%), the wet near-side lowlands round Sinus Medii and Mare Vaporum (L18, L23, 44–46%), the southern lake country at L35 (50%) and Jura (L09, 52%). The regions C′ fills most are the far side's wet lake countries (L21, L24, 62%; L11, 58%), where the blend's rule of the wettest land first puts its fields; they are also the land with the darkest nights. The record holds the same shares for 549 of the redraw's 554 land districts; the other 5 are lake districts, whose cells hold no land at this resolution.

### The fixed places against the general rules

**Jura** *(arithmetic)*. The redrawn atlas's island group at Montes Jura holds 67,700 km² of land: an island of 61,800 km² rising 2,750 m above its sea, and 119 islets. The land that drains into the bay, inside the gazetteer's circle of Sinus Iridum, is 15% of the island (19% for a circle 1.15 times as wide); 27% of the island stands more than 1 km up and 16% lies in the outer coast's belt. The garden the row names, bay-facing land below the high arc, is therefore 7% of the island, and half of the island, the back country that drains away from the bay, is named by neither part of the row. Counted as hybrid on the row's word that the island is mostly hybrid, the layout comes to 52% hybrid and 48% wild. The general rules alone would make the island 86% garden and 13% open with almost no wild land, since nearly all of it is dry land within 37 km of the sea: they reproduce the garden island and lose the wild crest and coast, which is what a fixed layout adds.

**The summit** *(arithmetic)*. The metropolis's 2,464 km² of land stand 9.4 km above the common level at 12° N of turned latitude, in a wet cell. Without rule 1 the general rules would call the ground open land of the rain belt, so the rule and the rules agree on its kind and differ on its reason.

**The heritage sites** *(arithmetic, each site's own 1.9 km cell)*. Of the conservation register's 27 sites, 18 lie under water on the design surface and 9 on land. Against the 28% atlas their computed fields were made on, Apollo 16 and SLIM now lie in lake cells, and Lunar Prospector's, LCROSS's and Hakuto-R's sites lie on land. In A the nine land sites fall in wild land (Apollo 14, IM-1), hybrid land (Apollo 12, Hakuto-R) and open land (Luna 9, Chandrayaan-3, Lunar Prospector, LCROSS, Luna 25). No rule here treats a site; the register reviews their evaluations before districts are sited near them.

### The maps

The maps are on the research drive in `/media/projectspace/terluna-research/session-tools/agents/atlas_land_use/maps/`, equal-area, in today's frame with the near side in the middle and the regions' borders and codes: `land_use_A.png`, `land_use_B.png`, `land_use_C.png` and `land_use_A_plus.png` cell by cell; the same four with `_districts` by district; and `land_use_variants.png` and `land_use_variants_districts.png` with A, B and C together. The grid behind them is `land_use_16ppd.npz` in the folder above.

### Limits

The climate behind the classes is an estimate, read on cells of about 170 km, and the nearest-cell assignment errs dry in the lake districts. The needs by band are bounds from the studies' printed ranges, worked on October's land and October's zones; the studies have not been run on the design surface. Distances to water are measured in 5° bands and capped at 300 km. A town, a crater garden town and a harbour are smaller than a cell, so the open land is where such places may stand, and no site is chosen. The order that places gardens and the order that places the wild are proposals, and decision 7 shows what the first of them costs.

## E. Decisions for the author

1. **The wild share the atlas is drawn on.** The register prefers half and allows 30% beside the hybrid fifth; both stand. *Options:* A, B, or B with one of the three ways of closing its gap. *Recommendation:* half wild, with the rain belt's shortfall closed by decision 3. With the blend, two thirds of the land has no use under either share, so the choice only decides whether the unneeded land is promised to the wild, and the register's preference costs no room that the low-cost end needs. *Cost:* C′'s rain-belt placement at the high-cost end must be worked again on the design surface; nothing in the atlas waits on it.

2. **Where the managed biomes belong.** Line 1044 and line 308 point different ways for the mosaics' fields and villages, and no row names the ordinary fields' kind. *Options:* the three readings of part C. *Recommendation:* reading 1, by purpose and access, because it keeps every row in force. *Cost:* the mosaics' villages become small habitation places inside garden country, which the atlas must be allowed to draw; reading 2 would take hybrid land out of the humid and wet 44% of the land.

3. **Whether the mosaics' harvest counts towards the diet.** The food accounts leave garden harvest out, and line 1044 asks lasting output of the mosaics. *Options:* leave it out; count it once the ecology lane has yields under the methane methods. *Recommendation:* count it. The candidates' flooded fields and valley floors are of the size of B's shortfall, and food grown where people walk is what the row describes. *Cost:* the methane room binds first: at the removal ledger's high-cost end it holds 22% of the flooded fields the valleys could carry.

4. **Whether the hybrid fifth grows.** Line 1045 leaves it open. *Options:* the mosaics inside one fifth (gardens 2.62 million km², 11.8% of the land); the mosaics beside a full fifth of gardens (28.0%, 6.20 million km², A+); growth only where a band's mosaics outrun its fifth. *Recommendation:* hold one fifth over the Moon with the mosaics inside it. The gardens that remain already line every dry shore to 37 km inland, and A+ turns half of the dry sub-humid land to garden. *Cost:* none to room; the equatorial band is then a third hybrid and the polar bands a twentieth.

5. **Whether the shares hold everywhere or over the Moon.** The register states the shares for the land; the studies hold them zone by zone; the mosaics cannot fit a band's fifth at the equator. *Options:* both shares in every zone; both over the Moon; the wild in every moisture class and night regime with the hybrid share over the Moon. *Recommendation:* the third, as the rules here have it, so that each climate's living cover and the darkest nights keep their wild while gardens go where they suit. *Cost:* B's shortfall is a consequence of holding half of the wet class wild; holding the half over the Moon would remove it and leave the rain belt's wild land below half.

6. **Kinds of land by whole districts or by cells.** Cells give exact shares and edges that belong to the climate grid; districts give natural edges and coarser shares, with hybrid land falling from 20% to 15% in A because garden rims seldom fill a district. *Options:* by cells; by districts, each in its leading kind; by districts as the atlas's unit, each carrying its three shares, with garden shores and mosaic valleys drawn inside districts as features. *Recommendation:* the third. A district is the atlas's whole place and its edges are the land's; the shares stay exact in its facts; and the two forms of hybrid land that are narrow by nature keep their shape. *Cost:* about an hour to rank districts for the hybrid share and redraw, one minute of machine time, once the redraw's districts are final.

7. **The order that places gardens, and resort country on wet coasts.** Nearest water first gives gardens 98% of the dry shore and leaves the wild none; and every garden lies in the drylands, so a warm wet coast is hybrid only in its mosaic valleys. *Options:* keep the order; turn gardens to sheltered water after Jura, bays and great lakes, and leave open coasts wild, using the coastal seas' records of which waters are enclosed; in addition, name resort shores on humid and wet coasts as hybrid land. *Recommendation:* the second for the drylands. For wet coasts the record offers no row to build on, so the paper lists the question and fills nothing in: 1.63 million km² of humid and wet land lies within 5 km of a sea or great lake. *Cost:* an hour of work for the second option; the third needs a word on what a resort shore is, and the missing studies of beaches and bathing water below.

8. **What the open land beyond the need is.** A leaves 33–44% of the land open and unneeded, B 12–24%. *Options:* draw all open land as habitation's; draw the filled land as habitation and the rest as open land held for growth; draw the rest as wild in fact. *Recommendation:* the second, as the maps here do, since it shows where people are and what is promised to no one. *Cost:* a fourth tone on the atlas.

9. **Jura's back country.** The row names 7% of the island as garden and 43% as wild, and half is unnamed. *Options:* the back country as garden country, which keeps the island mostly hybrid at 52%; as wild, which makes the island nine tenths wild. *Recommendation:* garden country, on the row's own word, with the high arc at 1 km and the coast belt at 5 km confirmed or replaced. *Cost:* none beyond the island.

10. **Where visitors sleep and where garden towns may stand.** The register gives hybrid land places to rest and no lodging, and gives crater garden towns a form and no setting. *Options:* lodging inside hybrid land as an amenity; in habitation places at the gardens' gates, as Jura's harbour town is the gateway, and aloft; crater garden towns as such gates inside garden and mosaic country, counted as habitation. *Recommendation:* gates and sky, with garden towns allowed as enclaves outside the kept wild. It keeps hybrid land free of homes, as line 308 has it, and uses forms the register has adopted. *Cost:* 206 million visitors at the scattered communities' 100–357 m² a person would take 20,600–73,500 km², under 0.4% of the land *(arithmetic)*.

11. **The dark sky on the map.** The share is met as land under an unlit sky, and the rules also hold it within the land of darkest nights. In A the blend's fields still fill 44% of that land, the far side's wet low latitudes. *Options:* leave it; rank the darkest land last for fields. *Recommendation:* leave it until the fields are placed on the design surface, since the studies count lit sky by homes, with none for fields, and the register's lighting direction (line 130) governs what is lit. *Cost:* none now.

**Inputs the record lacks**, none of them filled in here: the crater catalogue's count on this branch; C′ and the food blend run on the design surface and its moisture map in place of October's land and zones; the parks' water, machine energy and night light, unsized since October; any study of beaches, bathing water and lodging; yields and methane for the mosaics from the ecology lane; the review of the heritage evaluations; and a climate for the summit's massif that says how wet the metropolis's ground is.

## Running it

```sh
python3 -m research.studies.atlas_round.land_use                    # about 80 s on two processors; results/land_use.json and the grid on the drive
python3 -m pytest research/studies/atlas_round/test_land_use.py -q  # 13 tests, under a second
```

The record binds the script's hash, the hashes of the moisture map, the satoyama map, the surface's record and the heritage register, and the hashes of the five products it reads from the research drive.
