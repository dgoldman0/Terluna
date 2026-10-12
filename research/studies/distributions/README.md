# Distributions of people: Earth, the Moon, the array and the colonies

**11 October 2026.** Step six of the population work: coherent distributions of 20–30 billion people across Earth, the
Moon, the solar-shield array and the Solar System's budding colonies, for the author to choose among. They build on
the author's decisions of the day ([register](../../decisions.md#people-on-the-moon-and-beyond-11-october-2026)): 4–8
billion on the Moon with plenty kept wild, the aerial life above all; 2–4 billion of them aloft in sky towns that move
under power with tracking and air traffic control; the ground's people mostly in the drier highlands, the summit
metropolis round its tower above all, with many communities scattered elsewhere and taller cities; the array and the
colonies carrying a decent share, with Earth at 10–16 billion a bit much. Every account here is given outside a kept
30% and outside a kept half. The author then chose C′: C with the middle latitudes' lowland people moved where they are
better sustained and homestead communities among the farms there. It stands beside A–D below. The author also set the
wild share: at least half of the land and of the dark sky wild is preferred, and it can come down to 30% wild with
another 20% hybrid, a compromise between keeping 30% and keeping half, the hybrid land being heavily manicured
megaparks and great gardens. The study applies each rule to every zone and to the sky. C′ and A–D are also given under
the compromise, with the earlier reading, half wild with another 20% hybrid, beside them for comparison at the
author's request.

[run.py](run.py) writes [results/distributions.json](results/distributions.json) (schema
`terluna.research.distributions/1`) in under a second from committed products on this branch and seven files on other
branches, bound by commit and SHA-256. Numbers marked *screen* are computed by this code; *derived* marks arithmetic on
cited values; other numbers come from the product named beside them. Every range runs over the 16 corners of its case:
the highland share, the tower share and the scattered communities' density each at both ends, with every other input at
its favourable and at its unfavourable end; C′ adds the homestead density at both ends, 32 corners. Two primary sources
were read for the colonies ([sources.json](sources.json)); the checks are in [checks.json](checks.json).

## Findings

**The author's distribution, C′, is C with the middle latitudes' lowland people moved, at C's costs** (*screen*).
Ranked on spare land, water energy, comfort under both climate models, rain and lightning, and night lighting, the
high-latitude lowland homes poleward of 48.6° sustain people best, with a mean rank of 1.9 of 5. They take everyone
who leaves the lowland homes at 30–48.6°: 0.25–0.40 billion outside a kept 30%, and 0.22–0.40 billion outside a kept
half, where 0–28 million stay as homesteaders. Homesteads stay only where the middle latitudes hold farms, at 8 billion
outside a kept half at the high-cost end: 2.8–28 million people among 1.39 million km² of fields. C′ keeps C's heat,
energy, water lift and wild land, and puts 0.49–14.7% of the dry land under a light-polluted sky against C's
0.50–15.4%. Every corner fits outside a kept 30%. Outside a kept half the high-cost end fits with the remedy's poleward
and twilit placements once the homesteads' homes count within the fields they work, 0.02–0.71% of their land. Under
the author's compromise of 30% wild with another 20% hybrid, the hybrid fifth holds 4.68 million km² of megaparks and
great gardens, farms and towns keep to the same half, and C′ fits as it does outside a kept half.

1. **The Moon carries 4–8 billion people with 61–94% of its dry land wild outside a kept 30%, and 31–92% outside a
   kept half** (*screen*). With half its people aloft, on the food study's power-heavy mix (vegetables and fruit in
   lamp-lit vertical farms, the animal part as microbial protein, the rest on rain-fed fields), 4 billion leave 84–94%
   of the dry land wild outside a kept 30% and 77–92% outside half. Six billion leave 76–91% and 61–88%, and 8 billion
   64–88% and 34–84%. With 6 billion on the ground and 2 aloft, 8 billion leave 61–87% and 31–82%. The fields and the
   towns fit every band of the land at every corner outside a kept 30%, and at 4 and 6 billion outside half. At 8
   billion outside half every high-cost corner overfills: the rain-fed fields fill the middle-wet zone, which the
   scattered communities' land also covers, and at 7 or 8 of those corners the 30–45° highlands or the 0–30° land
   overfill as well. The favourable end fits, leaving 80–84% wild with 4 aloft and 76–82% with 2. Every corner fits,
   for twice the water lift, when the dry highland towns go poleward of 45° first, the scattered communities gather on
   the twilit and Earthlit homes and the middle band's communities keep to the fog-desert side of their homes. The mix
   draws no water from the seas outside a kept 30%.

2. **The ground's 2–6 billion live in the summit metropolis, the dry highland cities and towns scattered across the
   Moon, on 0.13–1.1 million km², 0.56–4.7% of the dry land** (*screen*). The summit metropolis is a city of a hundred
   million round the great sky port, its districts among crater lakes on the far side's equatorial massif, on 2,464 km²
   at 40,600 a km². The dry highland cities stand in cool, bright air poleward of 30°: tower districts of thirty storeys
   on the ridges at 85,700 a km², and highland towns on crater walls above high lakes at 10,000 a km². With 60–75% of
   the ground's people in the highlands and a third to half of the dry highlands' people in towers, these cities hold
   1.1–1.4 billion at 4 billion, 1.7–2.15 at 6 and 2.3–2.9 at 8. They settle 6.5–10%, 10–16% and 14–22% of the dry high
   ground outside a kept 30%, and 9.1–15%, 14–22% and 19–30% outside half. With 6 billion on the ground, 3.5–4.4
   billion settle 21–33% and 29–46%. The scattered communities, 0.50–2.4 billion, spread over the other homes in
   proportion to their land at 100–357 m² a person, 23–112 people a km² of home. At 8 billion outside a kept half the
   high-cost end's fields reach about 326,000 km² of the 30–45° highlands, and the dry high ground is then up to 78%
   used with 4 aloft and 94% with 2.

3. **The Moon uses 11–130 TW of primary energy; with fusion it takes 9.2–120 TW of heat, answered by 0.080–1.7% more
   dimming, and with power from the array 7.7–71 TW and 0.066–0.98%** (*screen*). Primary energy counts the residents'
   2–6 kW a person, with industry in orbit, and the electricity of the food, the sky towns and the water at fusion's
   input of 2.0–3.0 W a watt. At the favourable end the residents give 66–74% of the heat. At the high-cost end with
   half aloft, remaking the sky towns' hydrogen after 10% of it is lost a year gives the largest share, 45% of the heat
   with fusion (27–54 TW), and the residents 30%. With 2 of 8 billion aloft the residents give the largest share at
   both ends, 41–74%. The fleet's glow needs 1.9–3.4% more dimming (population study). Close computing at a hundred
   times today's data centres a person adds 2.3–4.6 TW and 0.040–0.19%.

4. **Water closes everywhere for gigawatts** (*screen*). Domestic water is 73–438 km³ a year, and lifting it takes
   0.36–42 GW. The dry high ground's own runoff, 80.6 km³ a year and nearly all of it at 30–45°, equals 105–400% of its
   people's gross domestic use at 4 billion and 33–126% with 6 billion on the ground. Raised from the nearest water,
   their supply takes 1.0–3.0 GW a billion at 30–45° and 6.1–19 GW a billion poleward of 45°, a tenth of that at 90%
   recovery. The chosen mix grows on 839–4,610 km³ of rain a year.

5. **Phosphorus comes from the lake deltas and back from the towns, with the natural cycle kept beside them**
   (*screen*, on resources' routes). The food loses 176–1,440 kt a year at the lunar design's recycling, 44–180 kt a
   billion. Dredging it from the deltas moves 0.079–3.0 Gt of sediment a year for 0.029–9.4 GW of fuel and takes
   4.6–107% of the particles the deltas catch from the land; the share can pass 100% because the people's own losses
   reach the same lakes. Struvite carries the dry highlands' and the sky's share, 1.1–9.8 Mt a year, up for 0.097–1.6
   MW. Beside it, carriers return the wild land's dissolved loss, 21–94 kt a year, through 0.60–5.3 million
   albatross-scale roosting flyers, and harvests take at most 30% of any water's bioavailable input.

6. **The dry highland people light almost none of the sky, and with half aloft the sky towns light the most of what
   remains at the high end** (*screen*, on living light's areas). A person on the dry high ground at 30–45° lights 0–36
   m² of polluted sky, against 28–1,260 m² across the highland homes as a whole, and poleward of 45° none. With half
   aloft, 4, 6 and 8 billion put 0.35–7.9%, 0.43–11.6% and 0.50–15.4% of the dry land under a light-polluted sky with
   the metropolis lit like a town of 10,000, and 0.83–10.4%, 0.91–14.2% and 0.98–17.9% with it lit at Earth's area a
   person, against 0.97–15%, 1.1–21% and 1.3–28% for living light's highlands mix with half aloft. At the high end the
   residents aloft light 1.1–2.2 million km²; at the favourable end the metropolis lit like a town lights 47,900 km²,
   against 25,500–51,000 km² from the residents aloft. Living light for 30–70% of outdoor lighting brings the shares
   with the metropolis as a town to 0.14–6.9%, 0.15–10.2% and 0.17–13.5%.

7. **Two to four billion aloft live in 348,000–696,000 sky towns, shading 0.72–1.44% of the Moon beside the
   aerophytes' central 1%** (*screen*). At 50 t a resident they stand 6.2–8.7 km apart over the sky outside a kept 30%.
   In the busiest band, the equatorial, where the 10 km winds gather 1.47 times its share, they stand 5.1–7.2 km apart.
   There a town comes within three diameters of another 1.0–20 times a day as neighbours drift 0.1–1 m/s apart, and
   0.50–10 times from its own powered moves. Across the fleet, at the density each town meets on average, that is
   0.066–5.3 million conflicts a day for tracking and automated traffic control to resolve. A town that shares its
   height with aerophytes at 1% cover passes 13–130 of them a day.

8. **Six of the thirteen undersea communities lie within 250 km of a mapped submerged tube, and joined they hold
   1.07–28.6 million people** (*screen*). Tranquility Base lies 170 km from rille M, Taurus-Littrow 152 km from rille O,
   Lunokhod 1 139 km from rilles D and E1, Lunokhod 2 103 km from rille R, Chang'e 3 205 km from the Sinus Iridum
   gravity low and Chang'e 5 173 km from rille E2; Hadley-Apennine lies 18 km from Rima Hadley's centre. The
   land-cities inventory gives rilles D and E1 one position, so both join Lunokhod 1, as every tube within 0.5 km of
   the nearest does. The six domes hold 22,400 people and their tubes 1.04–1.60 million at a third of their size and
   19.2–28.5 million at full size. All the submerged and coastal tubes hold 1.7–37 million (3.7–94 million with the
   gravity candidates), pumped dry for 2.45–245 MW per 50 km of tube under 500 m of sea (land cities).

9. **The undersea and tube cities' mineral wealth is the sea floor's sediment, and their seepage is fresh sea water**
   (*screen*). A 50 km tube under 500 m of sea seeps 0.058–5.8 km³ a year, the water of watered fields for 0.050–22
   million people. At the Near-side Sea's concentration after 1,000 years it carries 0.0029–1.2 t of phosphorus a year,
   against the 114–841 t its 2.6–4.7 million residents' food loses. The Near-side Sea's own sediment receives 3.5–10 kt
   of particulate phosphorus a year from its rivers, enough for 19.5–226 million people at the lunar design's losses:
   0.52–130 times the people the submerged and coastal tubes hold, and 0.21–62 times with the gravity candidates.

10. **The fog desert below the dry highland cities can feed them from glass that drinks piped sea water** (*screen*,
    on the food study's systems). A billion highland residents fed so take 0.25–0.97 million km² of vented glass,
    189–671 km³ of sea water a year piped a mean 69 km and 0.24–1.2 TW, with 62–155 kt N of N2O a year against 20–65 kt
    N a billion on the power-heavy mix. At 6 billion the dry highlands' 1.7–2.15 billion would take 8.7–43% of the fog
    desert, keeping their share of the rain belt wild. Their food climbs 3.5 km to the cities for 0.12 GW a billion.

11. **The dark belt keeps the Moon's dark nights, and the twilit and Earthlit lands need few lamps** (*screen*). The
    far side's low latitudes hold 3.10 million km² of dry land whose nights fall below 0.01 lux, 13% of it, where
    luminous signalling and dark-sky astronomy work. The metropolis lies in it and lights 1.5–3.2% of it like a town of
    10,000, 5.1–22% at Earth's area a person. Walking lamps burn 165 hours a cycle over the summit, 16 in the nearside
    towns and none in the high latitudes, and scattered communities four times as dense on the twilit and Earthlit
    homes light 44–45% of the sky the spread puts under their lamps.

12. **Earth, the array and the colonies carry their own costs** (*screen*, *derived*). Earth at 8–10 billion eats within
    its boundaries on 10–14 million km² of farmland, 22–29% of today's. Each billion beyond 10.2 needs 0.695–2.79 TW of
    power-made food, which the array sends down to 14,000–57,000 km² of rectenna while it spares 1.3–1.4 million km² of
    farmland. Each billion in the array needs 248–990 Gt of shield, 5.4–21 times the fleet, and 1,005 Gt of habitat,
    grown at 6,400–32,000 kg/s over 1,000–5,000 years. Each billion in the colonies uses 3.4–14 TW of primary energy
    with power-made food where sunlight is weak, 60 TW of power on Kalpana One's lamp-lit farms, and 2.0–5.5 kWe a
    person in a Mars city with 150 m² of greenhouse each; its free-space habitats take 1,005–2,330 Gt of local material.

## The four distributions and the author's C′

Earth holds 9 billion in each, within the 10.2 billion that eat inside four planetary boundaries. Each Moon row reads
its case outside a kept 30% unless it names the half (*screen*).

| | A: the Moon at 4 billion | B: the Moon at 6 billion | C: the Moon at 8 billion | D: 30 billion, the Moon on the land | C′: the author's distribution |
|---|---|---|---|---|---|
| Earth, billion | 9 | 9 | 9 | 9 | 9 |
| The Moon: ground + aloft | 2 + 2 | 3 + 3 | 4 + 4 | 6 + 2 | 4 + 4 |
| The array | 7 | 6 | 5 | 8 | 5 |
| The colonies | 5 | 4 | 3 | 5 | 3 |
| Total | 25 | 25 | 25 | 30 | 25 |
| Primary energy used by everyone, TW | 60.1–233 | 60.3–246 | 60.4–258 | 71.2–274 | 60.4–258 |
| The Moon's primary energy, TW | 11.1–65.4 | 16.7–98.1 | 22.3–131 | 20.2–99.7 | 22.3–131 |
| Heat and dimming, fusion | 9.18–59.6 TW, 0.0797–0.828% | 13.8–89.3 TW, 0.120–1.24% | 18.4–119 TW, 0.159–1.66% | 16.3–88.0 TW, 0.141–1.22% | 18.4–119 TW, 0.159–1.66% |
| Heat and dimming, power from the array | 7.65–35.3 TW, 0.0664–0.491% | 11.5–53.0 TW, 0.0996–0.736% | 15.3–70.6 TW, 0.133–0.982% | 13.9–56.5 TW, 0.120–0.785% | 15.3–70.6 TW, 0.133–0.982% |
| Dry land wild outside a kept 30% / half | 84–94% / 77–92% | 76–91% / 61–88% | 64–88% / 34–84% | 61–87% / 31–82% | 64–88% / 34–84% |
| Corners that fit outside half; their wild land | 16 of 16; 77–92% | 16 of 16; 61–88% | 8 of 16; 80–84% | 8 of 16; 76–82% | 16 of 32; 80–84% |
| Land settled | 0.13–0.37 million km² | 0.20–0.55 | 0.26–0.74 | 0.40–1.1 | 0.26–0.74 |
| Dry high ground settled, kept 30% / half | 6.5–10% / 9.1–15% | 10–16% / 14–22% | 14–22% / 19–30% | 21–33% / 29–46% | 14–22% / 19–30% |
| Dry land under a light-polluted sky, the metropolis as a town / at Earth's area | 0.35–7.9% / 0.83–10.4% | 0.43–11.6% / 0.91–14.2% | 0.50–15.4% / 0.98–17.9% | 0.43–13.6% / 0.91–16.1% | 0.49–14.7% / 0.97–17.2% |
| Sky towns; their shade on the Moon | 348,000; 0.72% | 522,000; 1.08% | 696,000; 1.44% | 348,000; 0.72% | 696,000; 1.44% |
| Conflicts a town a day in the equatorial band, drifting / powered | 1.0–10 / 0.50–5.0 | 1.5–15 / 0.76–7.6 | 2.0–20 / 1.0–10 | 1.0–10 / 0.50–5.0 | 2.0–20 / 1.0–10 |
| Conflicts a day across the fleet | 0.066–1.3 million | 0.15–3.0 million | 0.27–5.3 million | 0.066–1.3 million | 0.27–5.3 million |
| Hydrogen held in the air | 1.4–14 ppm | 2.1–21 ppm | 2.8–28 ppm | 1.4–14 ppm | 2.8–28 ppm |
| Phosphorus need; share of the deltas' catch | 176–720 kt; 4.6–54% | 264–1,080 kt; 6.9–81% | 352–1,440 kt; 9.2–107% | 352–1,440 kt; 9.2–107% | 352–1,440 kt; 9.2–107% |
| The array's shield | 1,740–6,930 Gt | 1,490–5,940 Gt | 1,240–4,950 Gt | 1,980–7,920 Gt | 1,240–4,950 Gt |
| Its habitats against the fleet's mass | 152 times | 130 times | 109 times | 174 times | 109 times |
| Days of the build's matter stream, shield | 74–297 | 64–255 | 53–212 | 85–339 | 53–212 |
| The colonies' primary energy, power-made food | 17.0–72.2 TW | 13.6–57.8 TW | 10.2–43.3 TW | 17.0–72.2 TW | 10.2–43.3 TW |
| The colonies' free-space habitats | 5,030–11,700 Gt | 4,020–9,330 Gt | 3,020–7,000 Gt | 5,030–11,700 Gt | 3,020–7,000 Gt |
| Earth's farmland, within its boundaries | 11.7–12.6 million km² | the same | the same | the same | the same |

Everyone's energy is primary energy: the residents' 2–6 kW a person on Earth, the Moon and in the array, with the
Moon's and the colonies' electricity and Earth's food off the land at fusion's input of 2.0–3.0 W a watt. At its top
it is 0.032% of the 851 PW the fleet's night half intercepts.

Earth's own range moves each total (*screen*). At 8 billion the totals fall to 24 and 29 billion and everyone's energy
by 2.0–6.0 TW. At 11 billion they rise to 27 and 32 billion; the 0.8 billion beyond 10.2 eat 0.556–2.23 TW of
power-made food from the array, sent down to 11,400–45,600 km² of rectenna, which heats Earth by 0.77–3.1 TW
(0.0011–0.0046 K, or 1.4–5.6 TW and 0.0021–0.0083 K from thermal plants), and everyone's energy rises by 5.1–18.7 TW.

**A, the Moon at 4 billion, the array and the colonies carrying more.** The Moon keeps 2 billion on the ground and 2
billion aloft, the least the author set for the sky, and so the most wild land and sky, the least heat and the least
lit sky of the four; the array holds 7 billion and the colonies 5. It rests on the array's 1,740–6,930 Gt of shield,
the largest continuing material stream after the build, and on colonies whose figures come from designs never built.

**B, the Moon at 6 billion.** The Moon holds 3 billion on the ground and 3 aloft, the food study's own split; the array
holds 6 billion and the colonies 4. It sits between A and C on every account, and every band of the land fits at every
corner under both kept shares.

**C, the Moon at 8 billion.** The Moon holds 4 billion on the ground and 4 aloft, the most the author set for the sky;
the array holds 5 billion and the colonies 3. Its 696,000 towns shade more of the Moon than the aerophytes' central
cover, and its heat and dimming are the largest of the four. Outside a kept half its high-cost end overfills the
middle-wet zone, the 30–45° highlands or the 0–30° land at 8 of the 16 corners, and the favourable end's 8 leave 80–84%
of the dry land wild. All 16 fit when the dry highland towns go poleward of 45° first, the scattered communities
gather on the twilit and Earthlit homes and the middle band's communities keep to the fog-desert side of their homes.

**D, 30 billion, the Moon at 8 billion mostly on the land.** The Moon holds 6 billion on the ground and 2 aloft,
keeping the sky at its least occupied; the array holds 8 billion and the colonies 5. It asks the most of the dry high
ground, the land and the array: 3.5–4.4 billion settle 21–33% of the dry high ground outside a kept 30%, the scattered
communities reach 2.4 billion, and the array's 8 billion need 1,980–7,920 Gt of shield. Outside a kept half only its
favourable end fits as spread (8 of 16 corners, leaving 76–82% wild); placed in C's three ways, all 16 fit for
2.23–82.7 GW of water lift. C puts more heat on the Moon, 18.4–119 TW against D's 16.3–88.0.

Every distribution shares these open questions: whether people grow up healthy at a sixth of Earth's gravity; the crop
model's water use, which sets the rain-fed land within a factor of two; comfort over the land, where CM1 and the GCM
disagree and both put the most comfort on high ground poleward of 30°; the sky towns' hulls, gas barriers and traffic
rules; and the colonies' capacity. The product holds each distribution's assumptions and uncertainties as lists.

**C′, the author's distribution.** C's 25 billion, with the Moon at 8 (4 on the ground and 4 aloft), the array at 5 and
the colonies at 3. The lowland towns of the middle latitudes move to the high-latitude lowland homes, which rank best on
the study's measures, and homesteads stay among the farms wherever the middle latitudes hold any. Its heat, energy,
water lift and wild land are C's, and it lights a little less of the sky. Every corner fits outside a kept 30%; outside
a kept half the favourable end fits, 16 of the 32 corners with the homestead density at both ends, and the high-cost
end fits with the remedy's poleward and twilit placements once the homesteads' homes count within the fields they
work. The next section gives it in full.

## The author's distribution, C′

On 11 October the author chose C with one change
([register](../../decisions.md#people-on-the-moon-and-beyond-11-october-2026)): "let's drop the population there and
distribute it where it's better sustained. Other than that closest to I guess option C", then "I didn't say all towns
necessarily, though homestead communities make the most sense there." C′ keeps C's totals: Earth 9 billion, the Moon 8
with 4 on the ground and 4 aloft and the metropolis at its 100 million, the array 5 and the colonies 3. In the middle
latitudes, the lowland homes of the 30–48.6° band over the middle-wet zone and the fog desert, homestead communities
live among the farms and the rest of the belt's people move (*screen*).

### Homesteads among the farms

Homesteads live at the density of Earth's family farms, a household of four on each 20–200 ha of field, 2–20 people a
km². The power-heavy mix puts fields in the middle latitudes only at 8 billion outside a kept half at the high-cost end,
where they fill the whole middle-wet zone outside the kept half, 1.39 million km². There 2.8–28 million homesteaders
live among the fields, their homes and gardens on 278–9,920 km² at the scattered communities' 100–357 m² a person. Each
household works its own 20–200 ha, 5–50 ha a resident. Planting, harvest and repair crews come from neighbouring towns.
Those corners fit only with the dry highland towns poleward of 45°, 206 km from the zone's mean latitude of 38.2°:
beyond a four-seat sky boat's 66.7 km with its 6 km climb and the sky ferry's 180 km route, and within the regional
ships' 1,500 km, on which a crew member's trip takes 4.80 kWh each way at 23.3 Wh a passenger-km. Crew towns on the
84,200 km² of highland that stays free at 30–45° at that end would bring the crews within a sky boat's reach: a 50 km
trip with the climb takes 4.16 kWh a boat, 1.04 kWh a seat (land cities). Everywhere else the fields lie
equatorward of 30°, the middle latitudes hold no farms, and all of the belt's lowland people move.

### Where they go

The other homes, ranked on the study's own measures (*screen*; 1 is best, and ties share their mean rank):

| Home | Spare land outside a kept 30% / half, the least over C's corners | Water with nothing recovered, W a person | Comfortable hours a lunar day, CM1 / GCM | Rain, mm a day; lightning, flashes a km² a year | Walking-lamp hours a cycle; lit sky, m² a person | Mean rank, kept 30% / half |
|---|---:|---:|---:|---:|---:|---:|
| High-latitude lowland homes, poleward of 48.6° | 3.93 / 2.76 million km² | 3.15–9.46 | 50 / 613 | 0.125; 0.0001 | 0; 0 | 1.9 / 1.9 |
| Dry highland cities, 60–90° | 118,000 / 75,100 km² | 6.25–18.8 | 249 / 611 | 0.125; 0.0001 | 0; 0 | 3.0 / 2.6 |
| Dry highland cities, 45–60° | 175,000 / 111,000 km² | 6.06–18.2 | 141 / 550 | 0.125–1.06; 0.0001–0.0019 | 0; 0 | 3.2 / 2.8 |
| Dry highland cities, 30–45° | 450,000 / 0 km² | 1.01–3.03 | 141 / 550 | 1.06; 0.0019 | 58; 0–36 | 3.0 / 3.3 |
| Lowland homes equatorward of 30° | 2.21 million / 0 km² | 3.15–9.46 | 0–2.6 / 160–233 | 2.7–5.0; 0.0034–0.0067 | 92–104; 23–1,730 | 3.9 / 4.4 |
| The lowland homes at 30–48.6° they leave, for comparison | | 3.15–9.46 | 38–100 / 501–516 | 1.06; 0.0023–0.0056 | 16–102; 0–1,730 | |

Spare land comes from the band check, with the belt's lowland land moved out. Water comes from land cities' lifts for
the dry high ground's bands and from ways of living's lift from the seas for lowland homes, which bounds the homes near
lakes from above. Comfort comes from land cities' height-corrected rings at each home's latitudes and mean height, rain
from ways of living's latitude bands, and lightning from main's representative places: the highland box at 44.7° S,
the polar plain at 75° N, southern Mare Nubium and the eastern Smythii headland for the low latitudes, and Russell and
the South Pole-Aitken coast for the belt itself. Walking lamps and lit sky come from living light, with the dry high
ground's rule. Rain and lightning count against a home, as in the author's choice of the drier highlands, and the water
measure carries supply. A home's rank on comfort is the mean of its CM1 and GCM ranks, on weather the mean of its rain
and lightning ranks, and on night lighting the mean of its lamp and lit-sky ranks; ties in the mean go to the dry
highlands.

Each person who moves goes to the best-ranked home that still has spare land outside the kept share at that corner, at
that home's own density. The high-latitude lowland homes rank first at both kept shares and hold everyone who moves:
0.25–0.40 billion outside a kept 30%, and 0.22–0.40 billion outside a kept half, where 0–28 million stay as
homesteads. They join the high-latitude towns, where long low light wraps the whole night and the polar highs are bases
for the polar sampling programme (ways of living). The ground's people stay 60–75% in the highlands, as in C.

### C′ against C as spread

| Outside a kept 30% / half | C as spread | C′ |
|---|---:|---:|
| People in the lowland homes at 30–48.6° | 0.25–0.40 billion | none / 0–28 million homesteaders |
| Moved to the high-latitude lowland homes | – | 0.25–0.40 / 0.22–0.40 billion |
| The high-latitude lowland homes' own land used | 0.668–3.81% / 0.935–5.34% | 1.29–7.34% / 1.70–10.3% |
| Dry land wild | 64–88% / 34–84% | the same |
| Land settled | 264,000–736,000 km² | the same |
| Heat with fusion; added dimming | 18.4–119 TW; 0.159–1.66% | the same |
| Water lifted | 0.753–27.5 GW | the same |
| Dry land under a light-polluted sky, the metropolis as a town / at Earth's area | 0.503–15.4% / 0.982–17.9% | 0.493–14.7% / 0.972–17.2% |
| The same with living light | 0.166–13.5% / 0.490–16.0% | 0.164–12.9% / 0.488–15.4% |
| Corners where every band fits | 16 of 16 / 8 of 16 | 32 of 32 / 16 of 32 |

The moved people keep their land a person, their energy and their water, so the Moon's totals stay C's. They gain
comfort, 50 and 613 hours a lunar day under CM1 and the GCM against 38–100 and 501–516, nights that need no walking
lamps, rain of 0.125 mm a day against 1.06, and 0.0001 flashes a km² a year against 0.0023–0.0056. Their own homes'
rain is 0.0972 mm a day, so their water comes from lakes and seas a median 32.8 km away (ways of living). Were they to
join the dry highland cities instead, the best-ranked dry band with room (30–45° outside a kept 30%, 60–90° outside a
kept half), the lift would rise to 0.782–28.1 GW and 0.993–32.0 GW and the dry high ground would be 16–23% and 22–33%
settled.

### The corners

Outside a kept 30% all 32 corners fit. Outside a kept half (*screen*):

| C′ outside a kept half | Corners that fit, of 32 | Overfill at the high-cost end | Dry land wild | Heat with fusion | Water lifted | Light-polluted sky, the metropolis as a town / at Earth's area |
|---|---:|---|---:|---:|---:|---:|
| As placed | 16 | the middle-wet zone by 278–9,920 km² at 16 corners, the 30–45° highlands by 13,900–39,700 km² at 12, the 0–30° land by 12,900 km² at 4 | 34–84%; 80–84% where every band fits | 18.4–119 TW | 0.753–27.5 GW | 0.493–14.7% / 0.972–17.2% |
| With the dry highland towns poleward of 45° first and the scattered communities gathered on the twilit homes | 16 | the middle-wet zone by 278–9,920 km² at 16 corners | 34–84% | 18.4–119 TW | 1.47–54.6 GW | 0.454–11.8% / 0.932–14.3% |
| With the homesteads' homes counted within the fields they work as well | 32 | none | 34–84% | 18.4–119 TW | 1.47–54.6 GW | 0.454–11.8% / 0.932–14.3% |

At 8 billion outside a kept half the high-cost corners fit without a fog-desert placement once the homesteads' homes
count within the fields they work. Beside the fields those homes overfill the full middle-wet zone by 278–9,920 km²,
0.02–0.71% of it, against C's 13,000–74,200 km² of towns; within the fields they take that share of the fields' land,
which the crop model's factor-of-two uncertainty in the rain-fed land covers. With the remedy's poleward and twilit
placements every corner then fits, leaving 34–84% of the dry land wild, for 18.4–119 TW of heat with fusion, 1.47–54.6
GW of water lift and 0.454–11.8% of the dry land under a light-polluted sky with the metropolis lit as a town
(0.932–14.3% at Earth's area).

### Under the author's wild share

The author prefers at least half of the land and of the dark sky wild, or 30% wild with another 20% hybrid as a
compromise ([register](../../decisions.md#people-on-the-moon-and-beyond-11-october-2026)), and the study applies each
rule to every zone and to the sky. The hybrid land is heavily manicured megaparks and great gardens, food gardens among
them, for leisure and entertainment, with the small shops, restrooms and places to rest that so much land needs;
residents on the ground and aloft use them, and visitors too, and heavy automation keeps them. Farms, grazing, the
homesteads' homes, cities, towns and sky towns keep outside a rule's wild and hybrid shares. At the author's request
the earlier reading, half wild with another 20% hybrid, which leaves farms and towns 30% of each zone and of the sky,
runs beside them for comparison with the same parks (*screen*; product keys `kept_50`, `wild_30_hybrid_20` and
`wild_50_hybrid_20`):

| C′, 32 corners | Half wild, preferred | 30% wild and 20% parks, the compromise | Half wild and 20% parks, for comparison |
|---|---|---|---|
| Corners that fit as placed | 16 | 16 | 16 |
| With the remedy's poleward and twilit placements | 16 | 16 | 16 |
| With the homesteads' homes within their fields as well | 32 | 32 | 32 |
| Farms, settlement and the homesteads' homes, share of the dry land | 8.2–33% | 8.2–33% | 8.5–20% |
| Megaparks and great gardens | none | 4.68 million km² | 4.68 million km² |
| Water from the seas for food, km³ a year | 0–19.6 | 0–19.6 | 0–3,250 |
| Heat with fusion; added dimming | 18.4–119 TW; 0.159–1.66% | 18.4–119 TW; 0.159–1.66% | 18.4–121 TW; 0.159–1.68% |
| Sky towns' share of the sky open to them; spacing in the busiest band | 33%; 4.29 km | 33%; 4.29 km | 55%; 3.32 km |
| Conflicts a day across the fleet, drifting | 0.74–7.4 million | 0.74–7.4 million | 1.2–12 million |
| Dry land under a light-polluted sky with the twilit placements, the metropolis as a town / at Earth's area | 0.454–11.8% / 0.932–14.3% | 0.454–11.8% / 0.932–14.3% | 0.454–11.8% / 0.932–14.3% |

Under the compromise the parks take the hybrid fifth, 4.68 million km², about 585 m² for each of the Moon's 8 billion
residents before visitors. Farms and towns keep to the same half as under a kept half, so C′'s land, energy, heat,
water and sky are the half's, and every corner fits with the remedy's poleward and twilit placements once the
homesteads' homes count within their fields. As placed, the high-cost end overfills the middle-wet zone by
278–9,920 km² at 16 corners, the 30–45° highlands by 13,900–39,700 km² at 12 and the 0–30° land by 12,900 km² at 4.
Towns, cities and sky towns put at most 14.7% of the dry land under a light-polluted sky as placed with the metropolis
lit as a town, and 17.2% with it at Earth's area a person. Even with every park lit at least 62.8% keeps a dark sky,
more than the preferred half and the compromise's 30%; the parks' own night light is still to size.
A–D fare under the compromise as under a kept half: A and B fit as placed, and C and D fit with the remedy's three
placements.

With half wild and another 20% hybrid, farms and towns keep to 30% of each zone. The power-heavy mix's high-cost end
then grows 1.06 million km² of its crops' equivalent by watering its fields, with 3,250 km³ a year from the seas and
0.53 TW of lift (food study), and the sky towns, kept to 30% of the sky, have 1.7 times as many conflicts a day. C′
fits every corner the same way as under the compromise. B then needs the remedy's three placements, and D fits at 10
of its 16 corners with them.

## The Moon at 4, 6 and 8 billion

### The places

The summit metropolis is a city of a hundred million round the Moon's great sky port, its districts set among crater
lakes on the far side's equatorial massif and flowing in under the tower between its legs. The dry highland cities
stand in cool, bright air poleward of 30°: tower districts on ridges, thirty storeys of homes, work and gardens on a
fifth of their ground, the streets between them open to at least a quarter of the sky, and highland towns on crater
walls above high lakes, under long twilit nights. The scattered communities are harbour towns on the seas and great
lakes, wetland and canopy towns, market villages among the rain belt's farms, nearside towns under the fixed Earth,
limb towns where it rises and sets, and the far side's and the high latitudes' twilit towns. Beneath the seas, domed
communities over the heritage sites join the tube cities in the submerged lava tubes of the maria. Above them all the
sky towns ride the winds at 10 km and move under their own power, their parks breathing the air of Earth's 1,610 m
(daily life).

### People and land

Per case, outside a kept 30% unless named (*screen*). Each range runs over every combination of the highland share,
60–75% of the ground's people, the tower share, a third to half of the dry highlands' people, and the scattered
communities' 100–357 m² a person.

| | 4 billion, 2 aloft | 6 billion, 3 aloft | 8 billion, 4 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|---:|
| Summit metropolis | 0.1 billion, 2,464 km² | 0.1, 2,464 | 0.1, 2,464 | 0.1, 2,464 |
| Tower districts, 85,700 a km² | 0.367–0.700 billion, 4,280–8,170 km² | 0.567–1.08, 6,610–12,500 | 0.767–1.45, 8,950–16,900 | 1.17–2.20, 13,600–25,700 |
| Highland towns, 10,000 a km² | 0.550–0.933, 55,000–93,300 | 0.850–1.43, 85,000–143,000 | 1.15–1.93, 115,000–193,000 | 1.75–2.93, 175,000–293,000 |
| Scattered communities, 100–357 m² a person | 0.500–0.800, 50,000–286,000 | 0.750–1.20, 75,000–428,000 | 1.00–1.60, 100,000–571,000 | 1.50–2.40, 150,000–857,000 |
| Undersea domes | 48,600 people | 48,600 | 48,600 | 48,600 |
| Land settled | 131,000–366,000 km² | 198,000–551,000 | 264,000–736,000 | 398,000–1,110,000 |
| Food land, power-heavy | 0.832–2.27 million km² | 1.24–3.40 | 1.66–5.23 | 1.66–5.23 |
| Share of the dry land used | 4.12–11.3% | 6.17–16.9% | 8.23–25.5% | 8.80–27.1% |
| Dry land wild outside a kept 30% / half | 84–94% / 77–92% | 76–91% / 61–88% | 64–88% / 34–84% | 61–87% / 31–82% |
| The same outside half, the corners that fit | 77–92% | 61–88% | 80–84% | 76–82% |
| Highlands used, kept 30% / half | 11–26% / 15–35% | 16–38% / 22–53% | 21–50% / 30–79% | 24–54% / 33–85% |
| Dry high ground settled, kept 30% / half | 6.5–10% / 9.1–15% | 10–16% / 14–22% | 14–22% / 19–30% | 21–33% / 29–46% |
| Dry high ground used with the fields, kept half | 9.1–15% | 14–22% | 19–78% | 29–94% |

The food study's town allowance, 25–100 m² a person on the highlands, is taken out of its land and this study's
settlement put in. The high ground is the 3.90 million km² of dry land 4 km and more above the sea; the dry high ground
is its 1.35 million km² poleward of 30°, 80% of it on the far side, under 1.36 mm of rain a day and 0.0047 flashes per
km² a year, against 6.16 mm and 0.022 flashes equatorward (land cities). The summit's massif is the Moon's wettest and
stormiest ground: 13.1 mm a day within 32 km of the tower and 0.046 flashes per km² a year at the port (main's places).

### The land bands

Each band of the dry land outside the kept share is checked with the chosen mix's fields and this study's settlement
together (*screen*). The fields follow the food study's fill, wettest band first. In the 30–45° band they take the
highland share of what is free there once the food study's own towns are spread over the highland bands by land, and
at 0–30° they stay within the food study's own fill. Rain-fed fields in the middle-wet zone cannot move to the fog
desert, so the zone is a band of its own. The scattered communities' land lies 47.5% at 0–30° and 25.2% in the middle
band, 30–48.6°, by their homes' latitudes. The products give no overlay of the homes on the food study's zones, so
their middle-band land covers the middle-wet zone by area, 51.5% of it on the zone. The dry highland towns spread over
the dry high ground's bands by land.

| Share of each band used, kept 30% / half | 4 billion, 2 aloft | 6 billion, 3 aloft | 8 billion, 4 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|---:|
| All the dry land at 0–30°, 7.71 / 5.51 million km² | 11–31% / 16–44% | 17–47% / 23–77% | 22–71% / 31.1–100.2% | 22–73% / 31.6–103% |
| Highlands at 30–45°, 574,000 / 410,000 km² | 6.5–10% / 9.1–15% | 10–16% / 14–22% | 14–22% / 19.0–110% | 21–33% / 28.9–125% |
| Highlands poleward of 45°, 373,000 / 267,000 km² | 6.5–10% / 9.1–15% | 10–16% / 14–22% | 14–22% / 19–30% | 21–33% / 29–46% |
| Middle-wet zone, 1.95 / 1.39 million km² | 0.33–1.9% / 0.47–2.7% | 0.50–2.9% / 0.70–4.0% | 0.67–3.8% / 0.93–105% | 1.0–5.7% / 1.4–108% |
| Fog desert and polar lowland, 5.77 / 4.12 million km² | 0.34–2.0% / 0.48–2.7% | 0.51–2.9% / 0.72–4.1% | 0.68–3.9% / 0.96–5.5% | 1.0–5.9% / 1.4–8.2% |
| Corners where every band fits, of 16 | 16 / 16 | 16 / 16 | 16 / 8 | 16 / 8 |

At 8 billion outside a kept half, the high-cost end fills the food study's wet zones to their limit. Its rain-fed
fields take the whole middle-wet zone outside the kept half, 1.39 million km², and reach about 326,000 km² of the
30–45° highlands by the food study's own rule. With the scattered communities' land on it as well, the middle-wet zone
overfills at all 8 high-cost corners, by 13,000–74,200 km² with 4 aloft, 16,200–92,800 km² with 3 and
19,500–111,000 km² with 2. With the dry highland towns there, the 30–45° highlands overfill by 13,900–39,700 km² at 6
of those corners with 4 aloft, by 13,900–71,800 km² at 8 with 3 aloft and by 34,200–104,000 km² at 8 with 2 aloft.
Where 60% of the ground's people are in the highlands and the scattered communities spread at 357 m² a person, the
0–30° land overfills too, by 12,900, 80,700 and 149,000 km². Those corners cannot be built as placed. Only the
favourable end fits, and it gives the feasible range beside the full one: 80–84% of the dry land wild with 4 aloft,
78–83% with 3 and 76–82% with 2, for 18.4, 17.3 and 16.3 TW of heat with fusion.

Three placements together make all 16 corners fit at 8 billion outside half. The dry highland towns go poleward of
45° first, which answers the 30–45° highlands; the scattered communities gather four times as densely on the twilit
and Earthlit homes, which answers the 0–30° land; and where the fields fill the middle-wet zone, the middle band's
scattered communities settle the fog-desert side of their homes. The poleward move alone, or with the twilit
gathering, fits the favourable end's 8 corners. With all three (*screen*):

| Outside a kept half, so placed | 8 billion, 4 aloft | 8 billion, 3 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|
| Corners that fit, spread by land / so placed | 8 / 16 | 8 / 16 | 8 / 16 |
| Highlands poleward of 45° used | 48–77% | 61–97% | 73–100% |
| Dry highland people left at 30–45° | none | none | 0–0.62 billion |
| All the dry land at 0–30° used | 31–97% | 31–98% | 31–99% |
| Middle-wet zone used; the fog-desert side of the middle band's homes used | 1.3–100%; 1.3–16% | 1.7–100%; 1.7–20% | 2.0–100%; 2.0–24% |
| Water lifted | 1.47–54.6 GW, against 0.753–27.5 spread by land | 1.85–68.7, against 0.948–34.6 | 2.23–82.7, against 1.14–41.7 |
| Dry land under a light-polluted sky, the metropolis as a town | 0.458–12.1%, against 0.503–15.4% | 0.413–10.4%, against 0.469–14.5% | 0.367–8.69%, against 0.434–13.6% |
| The same, the metropolis at Earth's area | 0.937–14.7%, against 0.982–17.9% | 0.891–12.9%, against 0.947–17.0% | 0.845–11.2%, against 0.912–16.1% |

The wild share stays as it is, since the same land is used, and the heat rises by at most 0.12 TW. The lift doubles
because poleward of 45° the land yields 0.04 km³ of runoff a year, and the water climbs a median 4.5–4.7 km from lakes
64–76 km away, against 0.75 km from 33 km at 30–45° (land cities). The 45–60° band then holds 1.73, 2.18 and 2.62
billion dry highland people, who take no lit sky here; lit at 36 m² a person they would add at most 62,000–94,000 km²,
0.27–0.40% of the dry land (*derived*).

### Food

The chosen mix is the food study's power-heavy mix at its full range. By the food study's zones it is the one mix with
fields that fits both ends of every input at 4, 6 and 8 billion under both kept shares (*screen*, from the food
product), and the band check above holds it at every corner except the 8 high-cost corners at 8 billion outside
half. The food menu below gives every other mix.

| | 4 billion | 6 billion | 8 billion |
|---|---:|---:|---:|
| Food land, the towns taken out | 0.832–2.27 million km² | 1.24–3.40 (to 4.03 outside a kept half) | 1.66–5.23 (to 6.97 outside a kept half) |
| Power | 0.539–3.39 TW | 0.808–5.09 | 1.08–6.78 |
| Rain used; water from the seas | 839–2,300 km³; none | 1,260–3,460; none | 1,680–4,610; up to 19.6 outside a kept half |
| N2O | 0.0258–0.0841 ppb a year | 0.0388–0.126 | 0.0517–0.168 |
| Phosphorus lost | 0.176–0.720 Mt a year | 0.264–1.08 | 0.352–1.44 |
| Phosphorus in microbial biomass, cycled | 0.365–1.18 Mt a year | 0.548–1.77 | 0.731–2.36 |

### Energy, heat and dimming

Primary energy counts the residents' 2–6 kW a person and every electrical load at fusion's input of 2.0–3.0 W a watt.
With fusion each watt of electricity heats the Moon by 2.0–3.0 W and the residents' energy by 1.51–4.52 TW a billion;
with power from the array by 1.38 W and 0.686 of the primary energy. Each terawatt of heat is answered by
0.0087–0.014% more dimming (food study, population study). *Screen*:

| | 4 billion, 2 aloft | 6 billion, 3 aloft | 8 billion, 4 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|---:|
| Residents' primary energy, 2–6 kW a person | 8.00–24.0 TW | 12.0–35.9 | 16.0–47.9 | 16.0–47.9 |
| Electricity: food | 0.539–3.39 TW | 0.808–5.09 | 1.08–6.78 | 1.08–6.78 |
| Electricity: sky towns' control, 1–10% of the time at 5 m/s | 0.138–1.38 TW | 0.207–2.07 | 0.276–2.76 | 0.138–1.38 |
| Electricity: hydrogen remade, 1–10% lost a year | 0.852–8.86 TW | 1.28–13.3 | 1.70–17.7 | 0.852–8.86 |
| Electricity: trips down and back, food lifted | 0.0393 TW | 0.0590 | 0.0786 | 0.0393 |
| Electricity: water lifted | 0.365–13.4 GW | 0.559–20.5 | 0.753–27.5 | 1.14–41.7 |
| Electricity in all | 1.57–13.7 TW | 2.35–20.5 | 3.14–27.4 | 2.11–17.1 |
| Primary energy | 11.1–65.4 TW | 16.7–98.1 | 22.3–131 | 20.2–99.7 |
| Heat, fusion | 9.18–59.6 TW, 0.242–1.57 W/m² | 13.8–89.3, 0.364–2.36 | 18.4–119, 0.485–3.15 | 16.3–88.0, 0.430–2.32 |
| of it the residents | 6.04–18.1 TW, 30–66% | 9.06–27.1, 30–66% | 12.1–36.2, 30–66% | 12.1–36.2, 41–74% |
| of it the food | 1.08–10.3 TW, 12–17% | 1.62–15.4, 12–17% | 2.16–20.6, 12–17% | 2.16–20.6, 13–23% |
| of it the sky towns' control | 0.276–4.18 TW, 3.0–7.0% | 0.414–6.27, 3.0–7.0% | 0.552–8.36, 3.0–7.0% | 0.276–4.18, 1.7–4.8% |
| of it the hydrogen remade | 1.70–26.8 TW, 19–45% | 2.56–40.3, 19–45% | 3.41–53.7, 19–45% | 1.70–26.8, 10–31% |
| of it trips down and food lifted | 0.0786–0.119 TW, 0.20–0.86% | 0.118–0.179, 0.20–0.86% | 0.157–0.238, 0.20–0.86% | 0.0786–0.119, 0.14–0.48% |
| of it water lifted | 0.000729–0.0406 TW | 0.00112–0.0620 | 0.00151–0.0834 | 0.00228–0.126 |
| Added dimming, fusion | 0.0797–0.828% | 0.120–1.24% | 0.159–1.66% | 0.141–1.22% |
| Heat, power from the array | 7.65–35.3 TW | 11.5–53.0 | 15.3–70.6 | 13.9–56.5 |
| Added dimming, power from the array | 0.0664–0.491% | 0.0996–0.736% | 0.133–0.982% | 0.120–0.785% |
| Close computing, 0.577 TW a billion | 2.31 TW, 0.0401–0.0972% | 3.46, 0.0601–0.146% | 4.62, 0.0801–0.194% | 4.62, 0.0801–0.194% |

Each share's range spans the two ends: with half aloft the residents give 66% of the heat at the favourable end and 30%
at the high-cost end, the hydrogen remade 19% and 45%. The hydrogen lost each year also stays in the air while the soils
take it up over about two years: 1.4–28 ppm at 2–4 billion aloft (ways of living), against the 0.53 ppm the air holds
now. Over the dense districts the heat is answered in the city: 60–181 W/m² over the metropolis and 129–388 W/m² over
the tower districts (ways of living).

### Water

Domestic water is 50–150 litres a person a day, 18–55 km³ a year a billion (*screen*):

| | 4 billion, 2 aloft | 6 billion, 3 aloft | 8 billion, 4 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|---:|
| Domestic water | 73.2–219 km³ a year | 110–329 | 146–438 | 146–438 |
| The dry high ground's people's use | 20.1–76.7 km³ a year | 31.1–118 | 42.1–159 | 64.1–241 |
| Its own runoff, 80.6 km³ a year, against that use | 105–400% | 68.4–259% | 50.7–191% | 33.4–126% |
| Power to lift the water | 0.365–13.4 GW | 0.559–20.5 | 0.753–27.5 | 1.14–41.7 |
| of it the dry high ground's | 0.333–12.7 GW | 0.515–19.6 | 0.697–26.4 | 1.06–40.1 |

The lift raises all the water from the nearest water. The dry high ground lifts it a median 0.75 km at 30–45° and
4.5–4.7 km poleward of 45°, for 1.0–3.0 GW a billion residents at 30–45° and 6.1–19 GW a billion poleward of 45° with
nothing recovered (land cities' lifts by band); spread by land its people need 3.0–9.1 GW a billion. The metropolis
takes the equatorial high ground's median lift, 0.46 km from water 7.6 km away, 0.61–1.8 GW a billion, and its own
ground's runoff, 11.7 km³ a year, is 2.1–6.4 times its people's use. Sky towns catch 137 litres a resident for each
millimetre of rain on their tops, so 0.0036–0.11 mm a day covers the 0.50–15 litres of makeup after 90–99% recovery.

### Phosphorus

At the lunar design's 0.044–0.18 kg lost a person a year, from resources' sea-phosphorus routes (*screen*):

| | 4 billion | 6 billion | 8 billion |
|---|---:|---:|---:|
| Eaten | 1,760–2,400 kt a year | 2,640–3,600 | 3,520–4,800 |
| Lost, the new phosphorus needed | 176–720 kt a year | 264–1,080 | 352–1,440 |
| Sewage recovery, 60–90% of excreta | 1,060–2,160 kt a year | 1,580–3,240 | 2,110–4,320 |
| Sediment dredged from the lake deltas | 0.0791–1.50 Gt a year | 0.119–2.25 | 0.158–3.00 |
| Its fuel | 0.0295–4.70 GW | 0.0442–7.05 | 0.0589–9.40 |
| Share of the particles the deltas catch from the land | 4.60–53.7% | 6.89–80.6% | 9.19–107% |
| Struvite to the dry highlands and the sky | 1.08–4.85 Mt a year | 1.64–7.34 | 2.20–9.84; 1.92–9.13 with 2 aloft |
| Its lift to the dry highlands | 0.0969–0.505 MW | 0.150–0.775 | 0.203–1.05; 0.308–1.59 with 2 aloft |

The deltas catch 1.34–3.83 Tg of particulate phosphorus a year from the land and keep 66–175 kt more from solution;
half the land lies within 50 km of standing water and 72% within 100 km, so sediment suits the farms near the lakes and
struvite the far land and the sky (resources). The design target needs the towns to return at least 92.5% of their
phosphorus at its low end. Beside these the natural cycle keeps its own: carriers return the wild land's dissolved
loss, 21–94 kt a year outside a kept 30% to half, as 0.60–5.3 million albatross-scale flyers roosting over it or lake
insects at 11–104% of Lake Mývatn's highest emergence, and fish, seaweed and sorption take at most 30% of any water's
bioavailable input, 21–56 kt a year in all.

### Light

Land under a light-polluted sky at the night's darkest moment with shielded lamps, from living light's areas a person
by home (*screen*):

| | 4 billion, 2 aloft | 6 billion, 3 aloft | 8 billion, 4 aloft | 8 billion, 2 aloft |
|---|---:|---:|---:|---:|
| Metropolis lit like a town of 10,000 / at Earth's area a person | 47,900–98,600 / 160,000–693,000 km² | the same | the same | the same |
| Dry high ground at 30–45°, 0–36 m² a person | 0–30,500 km² | 0–46,800 | 0–63,200 | 0–95,900 |
| Scattered communities | 9,370–638,000 km² | 14,000–956,000 | 18,700–1,280,000 | 28,100–1,910,000 |
| Residents aloft, 13–540 m² a person | 25,500–1,080,000 km² | 38,300–1,620,000 | 51,000–2,170,000 | 25,500–1,080,000 |
| Share of the dry land, the metropolis as a town | 0.354–7.88% | 0.429–11.6% | 0.503–15.4% | 0.434–13.6% |
| Share of the dry land, the metropolis at Earth's area | 0.832–10.4% | 0.907–14.2% | 0.982–17.9% | 0.912–16.1% |
| With living light for 30–70% of outdoor lighting, the metropolis as a town / at Earth's area | 0.141–6.94% / 0.465–9.40% | 0.153–10.2% / 0.478–12.7% | 0.166–13.5% / 0.490–16.0% | 0.154–11.9% / 0.478–14.3% |

Each home's lit land is capped at its land, and each lighting design takes its own figures. The dry highland people at
30–45° take the far side's twilit-night places' figure, the home at those latitudes and side, and poleward of 45° they
take none. Living light lights the metropolis over 27,300–91,900 km² as a town and 103,000–666,000 km² at Earth's area.

### The sky

Sky towns at 10 km with 50 t a resident (*screen*); kept 30% unless named:

| | 2 billion aloft | 3 billion | 4 billion |
|---|---:|---:|---:|
| Towns; at 20 t and 100 t a resident | 348,000; 139,000 and 696,000 | 522,000; 209,000 and 1.04 million | 696,000; 278,000 and 1.39 million |
| Shade on the Moon; within a three-diameter grid | 0.72%; 8.25% | 1.08%; 12.4% | 1.44%; 16.5% |
| Sky left wild outside a kept 30%, by shade / grid | 99.0% / 88.2% | 98.5% / 82.3% | 97.9% / 76.4% |
| Sky left wild outside a kept half, by shade / grid | 98.6% / 83.5% | 97.8% / 75.2% | 97.1% / 67.0% |
| Hydrogen; all the lifting gas | 14.9 Gt; 19.3 Gt | 22.4; 29.0 | 29.8; 38.6 |
| Spacing over the sky outside a kept 30% / in the equatorial band | 8.7 / 7.2 km | 7.1 / 5.9 | 6.2 / 5.1 |
| Spacing outside a kept half / in its equatorial band | 7.4 / 6.1 km | 6.0 / 5.0 | 5.2 / 4.3 |
| Towns under a median storm, the equatorial band | 8.79 | 13.2 | 17.6 |
| Conflicts a town a day in the equatorial band: drifting / powered / uncorrelated winds | 1.0–10 / 0.50–5.0 / 30 | 1.5–15 / 0.76–7.6 / 44 | 2.0–20 / 1.0–10 / 59 |
| Conflicts a town a day at the density towns meet, drifting / powered | 0.76–7.6 / 0.38–3.8 | 1.1–11 / 0.57–5.7 | 1.5–15 / 0.76–7.6 |
| Conflicts a day across the fleet, drifting / powered | 0.13–1.3 / 0.066–0.66 million | 0.30–3.0 / 0.15–1.5 | 0.53–5.3 / 0.27–2.7 |
| The same outside a kept half | 0.19–1.9 / 0.093–0.93 million | 0.42–4.2 / 0.21–2.1 | 0.74–7.4 / 0.37–3.7 |
| Aerophytes passed a town a day at 0.1%, 1% and 2% cover | 1.3–13, 13–130, 26–260 | the same | the same |

A conflict is a pass with the centres within three envelope diameters, 3 km, the inhabited-volume grid's spacing. The
rate is density × relative speed × swept width, the sky-ecology study's rule for aerophytes, which reproduces that
study's 3.3–33 passes a day for a 100 m body at 1% cover. Neighbours a few kilometres apart ride nearly the same wind,
and the zonal-wind product resolves 170 km, so their relative drift is taken as 0.1–1 m/s, the sky-ecology study's
figure for bodies on the same tracks; the 3-day winds' spread at 10 km, 1.5 m/s east and 1.4 m/s north, would give
two uncorrelated towns 2.9 m/s. Powered moves are the inhabited-volume roaming case's 5 m/s for 1–10% of the time.

The per-town rates are the equatorial band's, the busiest: at 10 km the winds gather 1.47 times its share of towns
there, and 66–87% of passive routes share a 171 km model cell with another (sky ecology), which the model cannot part
further. The fleet totals take the density each town meets, the bands weighted by their occupancy, 1.12 times the
even density, and count each conflict once.

### How many aloft

The same Moon with 2, 3 or 4 billion aloft (*screen*, kept 30%):

| | 6 billion: 2 / 3 / 4 aloft | 8 billion: 2 / 3 / 4 aloft |
|---|---|---|
| Dry land wild | 75–91% / 76–91% / 77–92% | 61–87% / 62–88% / 64–88% |
| Dry high ground settled | 13.6–21.6% / 10.0–16.0% / 6.48–10.4% | 20.6–32.8% / 17.1–27.2% / 13.6–21.6% |
| Heat, fusion; added dimming | 12.7–73.7 TW, 0.111–1.03% / 13.8–89.3, 0.120–1.24% / 14.8–105, 0.128–1.46% | 16.3–88.0, 0.141–1.22% / 17.3–104, 0.150–1.44% / 18.4–119, 0.159–1.66% |
| Dry land under a light-polluted sky, the metropolis as a town | 0.394–10.7% / 0.429–11.6% / 0.463–12.5% | 0.434–13.6% / 0.469–14.5% / 0.503–15.4% |
| Sky towns | 348,000 / 522,000 / 696,000 | 348,000 / 522,000 / 696,000 |

Each billion moved aloft frees 3.5–5.6 points of the dry high ground outside a kept 30%, adds 1.0–16 TW of heat with
fusion, most of it the hydrogen remade, and lights more of the sky than it would on the dry high ground.

## The food menu

Every mix of the food study at half aloft, outside a kept 30%, with this study's settlement in place of the food
study's towns (*screen*, from the food product). "To x%" is how far towards its high-cost end a mix fits before its
first zone fills.

| Mix | 4 billion | 6 billion | 8 billion | Power | N2O, ppb a year | Wild dry land |
|---|---|---|---|---|---|---|
| Rain-fed first | to 84.3% | to 56.2% | to 42.1% | 0.0250–0.0775 TW | 0.0225–0.0975 | 54–93% |
| Greenhouse-heavy | both ends | both ends | to 96.5% | 0.505–3.71 TW | 0.147–2.07 | 33–93% |
| Fog-desert glass | to 57.2% | to 32.2% | to 16.0% | 1.42–4.07 TW | 0.251–0.792 | 75–92% |
| Power-heavy | both ends | both ends | both ends | 0.539–6.78 TW | 0.0258–0.168 | 64–94% |
| Fully lamp-lit | both ends | both ends | both ends | 20.8–425 TW | 0.0813–0.403 | 96–99% |
| Power-made | both ends | both ends | both ends | 2.79–22.4 TW | 0.00663–0.0329 | 96–99% |
| Sea-supported | neither end | neither end | neither end | – | – | – |

Ranges span the three numbers over each mix's feasible range. Rain-fed first grazes up to 3.37 million km² of the fog
desert at its feasible end; greenhouse-heavy sends 1.13–35.7 Mt of phosphorus a year through fish sludge and draws
38.9–2,710 km³ of sea water; fog-desert glass draws 980–2,630 km³; the sea-supported mix asks more phosphorus of the
waters than their life can spare. Outside a kept half, fog-desert glass cannot be built at 8 billion, and
greenhouse-heavy fits to 91.5% of the way at 6 billion and to 63.9% at 8. N2O has no sink in the design air, so every
rate accumulates: 25.8–168 ppb in a thousand years on the power-heavy mix (food study).

## The tube and undersea cities joined

The thirteen undersea communities over the submerged heritage sites, each with the nearest lava tube to which the
land-cities inventory gives a capacity (*screen*; people at 70 m² a person, from a third of the rille's or collapses'
size to full size):

| Undersea community | Depth | Nearest tube | Its setting | Distance | Tube holds |
|---|---:|---|---|---:|---:|
| Apollo 11, Tranquility Base | 467 m | Rille M, Mare Tranquillitatis | under the sea | 170 km | 0.125–3.1 million |
| Apollo 15, Hadley-Apennine | 453 m | Rille N, Mare Tranquillitatis | coastal | 656 km | 0.311–7.68 million |
| Apollo 17, Taurus-Littrow | 1,037 m | Rille O, southern Mare Serenitatis | under the sea | 152 km | 0.0117–0.305 million |
| Surveyor 1 | 378 m | Rille I, Oceanus Procellarum | dry | 429 km | 0.103–2.58 million |
| Lunokhod 1 and Luna 17 | 865 m | Rilles D and E (E1), Oceanus Procellarum, at one position | under the sea | 139 km | 0.124–3.04 million |
| Lunokhod 2 and Luna 21 | 1,141 m | Rille R, eastern Mare Serenitatis | under the sea | 103 km | 0.0639–1.6 million |
| Luna 16 | 807 m | Rille L, Mare Tranquillitatis | dry | 382 km | 0.308–7.38 million |
| Luna 24 | 2,149 m | Rille L | dry | 636 km | as above |
| Chang'e 3 and Yutu | 1,211 m | the gravity low by the Sinus Iridum pit | under the sea | 205 km | 0.699–20.0 million |
| Chang'e 4 and Yutu-2 | 4,027 m | Rille Q, southern Oceanus Procellarum | under the sea | 3,370 km | 5,210–135,000 |
| Chang'e 5 | 747 m | Rille E (E2–E4), northern Oceanus Procellarum | under the sea | 173 km | 0.019–0.499 million |
| Chang'e 6 | 3,312 m | Rille Q | under the sea | 2,980 km | as above |
| Blue Ghost 1 | 2,233 m | Rille L | dry | 730 km | as above |

Hadley-Apennine lies 18 km from the centre of Rima Hadley, a rille the literature ties to tubes, which the land-cities
inventory gives no capacity, and Chang'e 5 lies 96 km from Rima Sharp's. Each dome holds 3,740 people on half its
floors, 48,600 across the thirteen (ways of living). The six within 250 km of a tube hold 22,400 in their domes and
1.04–1.60 million in their tubes at a third of their size, 19.2–28.5 million at full size, 14.3–20.0 million of it in
the Sinus Iridum gravity low if it is a tube. A city of this kind joins its dome's harbour and access tower to the
tube's pumped vault, under the near-side sea and beneath the Earth.

## Synergies and division of labour

Quantified where the products allow (*screen*); each proposal is labelled once.

**The tubes' seepage and the undersea cities' membranes** (the author's example). The 50 km reference tube under 500
m of sea seeps 0.058–5.8 km³ a year through its grouted ring, and under 2 km of sea 0.22–22 km³, all of it pumped up for
2.45–245 MW and 33–3,300 MW. That stream equals the water of watered fields feeding 0.050–22 million and 0.19–82
million people, which the coast's farms can take beside the tube city's 2.6–4.7 million people. The seas are fresh and
poor in phosphorus for thousands of years: at the Near-side Sea's 0.049–0.20 mg/m³ after 1,000 years the stream
carries 0.0029–1.2 t of phosphorus a year under 500 m of sea, against the 114–841 t its residents' food loses. The
minerals that serve the farms lie in the sea floor's sediment: the Near-side Sea's rivers bring it 3.5–10 kt of
particulate phosphorus a year, enough for 19.5–226 million people at the lunar design's losses, and the mouths of all
the seas receive 50.7–168 kt, enough for 0.281–3.81 billion. Proposal: the undersea and tube cities of the near-side
sea dredge and work that sediment for the coast's farms and send struvite inland, and their vaults hold the lamp-lit
farms of the power-heavy mix beside the people they feed, out of the weather and at one temperature.

**The lake deltas for the highland and lowland farms.** The deltas catch 1.34–3.83 Tg of particulate phosphorus a year,
and 4–8 billion people need 0.176–1.44 Tg of it back; dredging moves 0.079–3.0 Gt of sediment for 0.029–9.4 GW of fuel
(the phosphorus tables above). Lowland farms within 50 km of a lake take sediment, half the land; the dry highland and
sky towns take struvite, 7.9 kg a kilogram of phosphorus, lifted for under 1.6 MW.

**The fog desert's glass drinking piped sea water.** The fog desert, 4.82 million km² at 30–60° below 4 km with 0.18 mm
of rain a day, lies below the dry highland cities and a mean 69 km from the seas. Vented glass there grows a billion
residents' plant foods and fish feed on 0.246–0.974 million km² with 189–671 km³ of sea water a year, 0.019–0.0882 TW
of it for the pipes and 0.243–1.2 TW in all, with 62.4–155 kt N of N2O a year (food study). For the dry highlands'
people:

| | 4 billion | 6 billion | 8 billion |
|---|---:|---:|---:|
| Dry highland people | 1.1–1.4 billion | 1.7–2.15 | 2.3–2.9 |
| Glass; share of the fog desert | 0.271–1.36 million km²; 5.6–28% | 0.418–2.09; 8.7–43% | 0.566–2.82; 12–59% |
| Sea water | 208–939 km³ a year | 321–1,440 | 435–1,950 |
| Power | 0.267–1.68 TW | 0.413–2.58 | 0.559–3.48 |
| N2O | 68.6–217 kt N a year | 106–333 | 144–450 |
| Food lifted 3.5 km to the cities | 0.136–0.173 GW | 0.210–0.265 | 0.284–0.358 |

Proposal: the dry highland cities draw part of their food from the fog desert's glass below them, which leaves their
share of the rain belt wild and puts the desert's sun to use, at 62.4–155 kt N of N2O a billion against the
power-heavy mix's 19.8–64.6.

**The dark belt and the twilit high latitudes.** The dark belt is the far side's low latitudes, 3.10 million km² of dry
land whose nights fall below 0.01 lux, 0.0014 lux at the summit at their darkest; the highland box at 44.7° S keeps
1.44 lux and the polar plain 485 lux (main's places). Walking lamps burn 165 hours a cycle over the summit, 84 in the
highland homes as a whole, 92–104 in the low-latitude homes, 16 in the nearside towns and none in the high latitudes,
where every night keeps walking light (living light). The metropolis lights 47,900–98,600 km² of the dark belt like a
town of 10,000 and 160,000–693,000 km² at Earth's area a person. Placing the scattered communities four times as
densely on the twilit and Earthlit homes lights 44–45% of the sky the spread lights: at 6 billion, 6,180–426,000 km²
against 14,000–956,000. Proposal: the dark belt outside the metropolis's lit circle holds the Moon's dark-sky
reserves, its far-side radio and optical astronomy and its luminous ecology, the metropolis lit like its small towns;
the twilit high latitudes and the Earthlit near side take the scattered communities that want bright nights.

**The array's power and computing behind Earth's food.** Each billion people on Earth beyond 10.2 billion eats
0.695–2.79 TW of power-made food. From the array that is 14,200–57,000 km² of rectenna at the 1978 reference's 49 W/m²
and 3,340–20,100 km² of collectors, 0.96–3.8 TW of heat on Earth against 1.7–7.0 TW from thermal plants, and
0.000082–0.00033% of the night half's light; it spares 1.3–1.4 million km² of farmland (population study). The array's
computing serves Earth's training, science and batch work over the 2.56 s round trip (computing note).

**More weaknesses that serve.**
- *The summit's rain.* The wettest ground on the Moon, 13.1 mm a day, fills the metropolis's lakes, a quarter of its
  ground, and its own runoff is 2.1–6.4 times its people's water; the wet equatorial high ground's 2.54 million km² at
  6.16 mm a day is the Moon's water tower.
- *The sky towns' rain.* A town's top catches 137 litres a resident for each millimetre, so 0.0036–0.11 mm a day closes
  its water after recovery; proposal: the towns gather rain for their own farms' water.
- *The array's shield as water.* Built of water, the array's whole habitat mass is 0.003–0.1% of what the build
  delivers (resources' industry note); proposal: the habitats' shield doubles as the array's water store.
- *The tubes' stable rock.* A tube vault holds one temperature through the 354-hour day and night (land cities);
  proposal: tube vaults hold the towns' cycle of food from the dusk harvest and the lamp-lit farms.

## Earth, the array and the colonies

**Earth** (*screen*, from the population and ways-of-living products):

| Earth, billion | 8 | 9 | 10 | 11 | 12 | 16 |
|---|---:|---:|---:|---:|---:|---:|
| Within four planetary boundaries | yes | yes | yes | to 10.2 | to 10.2 | to 10.2 |
| Farmland, plant-rich diet | 10.4–11.2 million km² | 11.7–12.6 | 13.0–14.0 | 13.3–14.3 | 13.3–14.3 | 13.3–14.3 |
| Share of today's farmland | 22–23% | 24–26% | 27–29% | 28–30% | 28–30% | 28–30% |
| Food off the land, power-made | – | – | – | 0.556–2.23 TW | 1.25–5.02 TW | 4.03–16.2 TW |
| Its rectenna from the array | – | – | – | 11,400–45,600 km² | 25,500–103,000 km² | 82,300–330,000 km² |
| Primary energy, 2–6 kW a person and the food off the land at fusion's input | 16.0–47.9 TW | 18.0–53.9 | 20.0–59.9 | 23.1–72.7 | 26.5–87.1 | 40.1–145 |
| Warming from thermal plants, the residents' energy | 0.024–0.072 K | 0.027–0.081 | 0.030–0.090 | 0.033–0.099 | 0.036–0.11 | 0.048–0.14 |
| Warming from the food off the land, array / thermal plants | – | – | – | 0.0011–0.0046 / 0.0021–0.0083 K | 0.0026–0.010 / 0.0047–0.019 | 0.0083–0.033 / 0.015–0.061 |

**The array**, per billion residents (*screen*, *derived* from the population product and resources' industry note):

| Account | Per billion |
|---|---|
| Shield, a quarter of the Stanford torus's to all of it | 248–990 Gt, 5.4–21 times the fleet's 46 Gt |
| Habitat with the torus's 15 t of structure a person | 1,005 Gt |
| Grown over 1,000–5,000 years | 6,400–32,000 kg/s |
| Days of the build's matter stream, the shield | 10.6–42.4 |
| Farm under continuous sunlight; the torus's electricity | 20,000 km²; 3 TW |
| Energy at 2–6 kW a person; its heat on the Moon | 2.0–6.0 TW; 0.00010–0.00030 W/m² |
| Collectors for that energy | 4,900–22,000 km² |
| Commuting by rocket, per million commuters | 17.6–45.7 GW weekly; 4.2–10.8 GW on rotations each cycle |

**The colonies**, per billion residents (*derived* from the project's figures, Globus et al. 2007 and Donofrio & Kirk
2021):

| Account | Per billion |
|---|---|
| Living at 2–6 kW a person | 2.0–6.0 TW |
| With power-made food where sunlight is weak, its electricity at fusion's input | 3.39–14.4 TW |
| On lamp-lit farms, Kalpana One's 60 kW of continuous power a resident (50 for agriculture) | 60 TW |
| A Mars city's electricity, 2 kWe (the Mars Society's minimum) to 5.46 kWe a person | 2.0–5.5 TW |
| Its greenhouses, 150 m² a person at 1,600 kcal a day; 234 m² at the lunar diet's 2,500 | 150,000–234,000 km² |
| Free-space habitats of local material, the Stanford torus to Kalpana One's 2,330 t a person | 1,005–2,330 Gt |
| Mirrors beyond the torus's 20 m² of sunlit farm a person: main belt, 2.7 AU / Jupiter's distance, 5.2 AU / 40 AU | 126 / 521 / 32,000 m² a person |
| Light time from the Moon: main belt / Jupiter's distance / the Kuiper belt | 14–31 min / 35–52 min / 5.4–5.7 h |

With the farm's own 20 m² the collecting area is 146, 541 and 32,000 m² a person. Kalpana One houses 3,000 people in a
cylinder 250 m in radius and 325 m long at 2 rpm, behind up to 10 t/m² of lunar or asteroid shielding, perhaps seven
million tons in all, with 170 m² of living area at 1 g a resident and 50 m² of lamp-lit farm to feed each (Globus et
al. 2007). City NPC plans a million people in Mars's lava tubes and craters at 5 GWe from forty 141 MWe reactors, with
15,000 ha of greenhouse and lava-tube districts 25 km long and 100 m across for 25,000 people each (Donofrio & Kirk
2021, a design-competition plan).

## Assumptions

Stated once; the product's `assumptions` block holds them all.

- **Design guesses**, each a range: 60–75% of the ground's people in the highlands, the metropolis included (the
  author's "mostly"); a third to half of the dry highlands' people in tower districts of 30 storeys; relative drift of
  0.1–1 m/s between neighbouring sky towns (the sky-ecology study's figure for aerophytes on the same tracks); and in
  C′, 2–20 homestead residents a km² of the middle latitudes' fields, a household of four on each 20–200 ha, the span of
  Earth's family farms from mixed-farming country to the mechanised grain belts. The first two vary independently over
  the corners, and the homestead density over C′'s 32.
- **Placement.** The metropolis keeps its decided 100 million on its 2,464 km²; the other highland people live on the
  dry high ground poleward of 30°, spread over its bands by land, or poleward of 45° first where the bands overfill;
  the undersea domes hold 48,600; the scattered communities spread over the other homes in proportion to their land at
  100 m² a person (the lunar design's compact town) to 357 m² (the affluent spread), or four times as densely on the
  twilit and Earthlit homes. Their middle-band land covers the middle-wet zone by area, or keeps to the fog-desert side
  of the homes where the fields fill the zone. Tube cities count among the scattered communities and cover no land. In
  C′ the scattered communities of the 30–48.6° lowland homes go to the best-ranked home with spare land, except
  homesteads among the farms there, whose homes take the scattered communities' land a person beside the fields, or
  within the fields they work in the variant that fits the high-cost corners.
- **Food.** The power-heavy mix at its full range, its town allowance replaced by this study's settlement and checked
  band by band.
- **Energy.** The provisioning band, 2–6 kW a person of primary energy with industry in orbit, and every electrical
  load at fusion's input; fusion supplies the Moon's power, with power from the array beside it; sky towns at 10 km,
  50 t a resident, powered 1–10% of the time at 5 m/s, losing 1–10% of their hydrogen a year; trips down for a fifth
  of the residents a day (ways of living).
- **Water.** 50–150 litres a person a day, raised from the nearest water; 0% and 90% recovery bound the lift.
- **Land rules.** 30% and half of every zone and of the sky kept wild, the Kunming–Montreal and Half-Earth levels in
  ecology's account; the author's compromise of 30% wild with another 20% hybrid; and half wild with another 20%
  hybrid, kept for comparison. The hybrid part is the author's megaparks and great gardens, and farms, grazing, the
  homesteads' homes, cities, towns and sky towns keep outside a rule's wild and hybrid shares.
- **Elsewhere.** Earth's people within its boundaries eat a plant-rich diet on 1.3–1.4 million km² a billion; the
  array's residents need a quarter to all of the Stanford torus's shield; the colonies use the provisioning band with
  power-made food, and each distribution sets their people, with no capacity computed for them, as in the shared
  cases.

## Departures from the known designs

- **The middle-wet fields.** The band check holds the middle-wet zone's fields to the room the food study's fill gives
  them, the zone's land outside the closed share. The fill reaches that room exactly when it fills the zone, and the
  food product's three significant figures can put the field area up to 0.5% above it: 835,000 km² against 834,000
  with half wild and another 20% hybrid. Held to the room, C′'s corners fit once the homesteads' homes count within
  their fields (16 of 32 otherwise), and with the remedy B, C and D fit 16, 16 and 10 of their corners (8 each
  otherwise).
- **The highland people's latitudes.** The food study spread its towns over every highland band at 25–100 m² a person,
  and land cities set half to all of the ground's people in the highlands at one density. Here 60–75% live in the
  highlands, all but the metropolis's 100 million poleward of 30°, in towers and towns.
- **The band check.** The food study fitted its mixes zone by zone with its own towns; here each band is checked with
  this study's settlement in place, the middle-wet zone as a band of its own, and the corners that overfill are marked
  unbuildable.
- **The summit's weather.** The author places the metropolis among the drier highlands; its massif is the Moon's
  wettest and stormiest ground, and the dry highlands lie poleward of 30°, so the metropolis keeps its site and its
  scale and the dry high ground takes the rest of the highland people.
- **The scattered communities' land**, 100–357 m² a person, where land cities took 10,000 a km².
- **The lit sky of the dry high ground** takes the far side's twilit-night places' 0–36 m² a person at 30–45° and none
  poleward of 45°, the band edge nearest living light's 48.6°. Lit at 36 m² a person, the whole 45–60° band would add
  at most 37,000 km², 0.16% of the dry land, at 8 billion with 2 aloft as spread, and 62,000–94,000 km², 0.27–0.40%,
  with the dry highland towns poleward of 45° first (*derived*).
- **The sky's conflicts** use the sky-ecology study's encounter rule and its relative drift, with the inhabited-volume
  grid's three diameters as the conflict distance.
- **Lunokhod 1's tube.** Rilles D and E1 share one position in the land-cities inventory, and both join it.
- **C′'s measures** take representative values: lightning at main's places, rain by latitude band, comfort from the
  land-cities rings by latitude and height, and the lowland homes' water from the seas, which bounds homes near lakes
  from above.
- **C′'s homesteads' homes**, counted within the fields they work where those fill the middle-wet zone, take 0.02–0.71%
  of the land the food study's fill gives to crops.
- **Earth** holds 9 billion in each distribution; the shared cases held 8–14.

## Inputs and bindings

Committed on this branch, hashed in `producer.inputs` and each feeding the basis: the food-systems, land-cities,
ways-of-living, daily-life, inhabited-volume and conservation heritage-register products, the provisioning population
and computing products, the zonal-wind product and the shared population cases. Bound on other branches and read
through Git:

| Input | Branch and commit | sha256 |
|---|---|---|
| `resources/results/sea_phosphorus.json` | domain/resources, 98d5ae9 | `9cf6927649942ed6…` |
| `resources/industry.md` | domain/resources, 98d5ae9 | `3da7e7b7d7dfa144…` |
| `resources/results/phosphorus.json` | domain/resources, 98d5ae9 | `159462a220b6661e…` |
| `biosphere/ecology/results/people.json` | domain/ecology, 4932a68 | `f8eb2892b733fdd9…` |
| `research/studies/living_light/results/living_light.json` | domain/ecology, 4932a68 | `5c0e0ae098f13fea…` |
| `research/studies/sky_ecology/results/sky_ecology.json` | domain/ecology, 4932a68 | `d21e6c16c3ba7ded…` |
| `research/studies/joint_synthesis/results/places.json` | main, 80497dd | `9707e4d1c7332ba3…` |

The runner refuses bytes whose sha256 differs; the full hashes and Git blobs are in `run.py` and in the product's
`producer.bound_inputs`. The values taken from the industry note are transcribed in `TRANSCRIBED` with the phrase each
comes from, and the tests find every phrase in the bound bytes.

## Reproduce and check

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n 10 python3 -m research.studies.distributions.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n 10 \
    python3 -m pytest -q research/studies/distributions
```

The tests check that the product regenerates; the producer, input and bound hashes and the shared constants; that the
transcribed phrases occur in the bound note, the Trojan and Kuiper-belt rows among them, and every published figure
names a source; that every committed input feeds the basis; that people are conserved at every combination of the
guesses, settled areas follow their densities and the guesses vary independently; that the food land equals the food
study's fields with its towns taken out, and the chosen mix fits the food study's zones; that the land bands partition
the dry land and flag overfill, the middle-wet zone is its own band, the fields reach the dry high ground by the food
study's rule and the three placements together fit; the energy bases and the heat by component; that the phosphorus
follows the need; that the lit-sky arithmetic keys the metropolis by design, reproduces living light's own spread shares
and caps each home; that the spacing reproduces ways of living's 5.51 km and the encounter rule the sky-ecology study's
passes; that the fleet totals take the density towns meet; the great-circle distance and the tube pairs and ties;
Earth's boundary rule and its sensitivity line; the colonies' mirrors beyond the farm; that every distribution keeps the
author's ranges and its energy sums; that C′ keeps C's totals, keeps homesteads at the guessed density among the
middle-wet fields, places everyone who moves within spare land in the ranked order and conserves people, that the
ranking follows its measures, that the high-cost corners outside a kept half overfill by the homesteads' homes alone
once the remedy's placements apply and fit with those homes within the fields, and that the crews' distance to the 45°
line lies beyond a sky boat's range and within a regional ship's; that the land rules match the food study's and keep
a kept half apart from the compromise; that under the compromise C′'s results match the kept half's beside the parks
and its corners fit as the README gives them under every rule; that the middle-wet fields are held to the fill's room
and lose at most the food product's rounding; that the hybrid part is parkland, kept from grazing and the homesteads'
homes; that the wild shares count against the land and sky outside the closed share; that A–D fit under the
compromise and the comparison as stated; and the README's words and lines. They skip cleanly
when an input is absent. [checks.json](checks.json) records the test run and the three rounds of the independent review.

## What remains

- The dry high ground's own weather over real terrain beyond the one CM1 box at 44.7° S, and the wet massif's storms
  round the metropolis.
- The water and the land of the dry highland cities poleward of 45°, where 8 billion outside a kept half would put them.
- In C′, the high-latitude towns' water from lakes and seas a median 32.8 km away, the work of the middle latitudes'
  farms beyond the homesteads, and the homesteads' own design.
- The parks under the author's compromise: their water on dry land, in the dry highlands and in the fog desert and
  polar lowland, which hold 1.65 million km² of them; the energy of the automation that keeps them; their amenities;
  and their night light within the hybrid fifth of the dark sky.
- An overlay of the scattered homes on the food study's zones, and the middle band's towns beside the middle-wet fields.
- The heights of the sky towns, the aerophytes and the sky reefs, and the traffic rules on the shared tracks; the
  relative drift of neighbours at kilometre scale, below the 170 km winds' resolution.
- How fast the dredging of the lake deltas can run without harming their life, and the struvite and concentrate routes.
- The joined undersea and tube cities' passages, pumps and harbours, and the fresh seas' minerals over thousands of
  years.
- The array's habitats as designs with orbits, shields and food, and the share of commuters to residents.
- The colonies' capacity and their habitats' mass beyond the two designs read.
- A lifetime's health at a sixth of Earth's gravity, and the warming of the N2O and methane every mix adds.
