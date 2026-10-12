# Cities on the land: the lava tubes and the highlands

**11 October 2026.** On 11 October the author set step six's working range: 4–8 billion people on the Moon, 2–4
billion of them aloft, and on land most of them in the drier highlands, the summit metropolis round its tower above
all, with many communities scattered elsewhere and cities taller at lunar gravity. The author also asked for the large
lava tubes as the first phase's habitats and as lasting cities, protected from water percolation
([register](../../decisions.md#people-on-the-moon-and-beyond-11-october-2026)). This study answers both on the 28%
[atlas](../../../geography/README.md#atlas-at-the-selected-water-share) and checks two statements made to the author:
that the known skylights and tube candidates lie in mare basalt in the low basins, so that many would end up under
water or below the water table, and that a tube 1 km wide and 50 km long with ten floors holds about 7 million people
at 70 m² of floor each.

[run.py](run.py) writes [results/land_cities.json](results/land_cities.json) (schema `terluna.research.land-cities/1`)
in about 25 seconds, with [tubes.py](tubes.py), [groundwater.py](groundwater.py), [highlands.py](highlands.py) and
[grid.py](grid.py). Numbers marked *screen* are computed by this code; *derived* marks arithmetic on cited values; other
numbers come from the product or source named beside them. Sources are in [sources.json](sources.json), with how far
each was read, and the checks in [checks.json](checks.json).

## Findings

### Lava tubes

1. **The known tube candidates lie in the low basins' mare basalt, and the 28% sea covers most of them: the statement
   holds.** Of the 17 pit sites in and beside the maria (the LROC atlas's 16 mare pits, two of them nested, and the two
   pits at the maria's edge that Wagner & Robinson 2014 counted with them), 13 lie under the sea, 54 m to 2.25 km down,
   and the Mare Insularum pit under a rain-fed lake 5 m deep (*screen*). Three stand on land: the Marius Hills pit 60 m
   above the sea and 15 km from its shore, the Schlüter pit 157 m up beside a lake, and the Mare Tranquillitatis pit 784
   m up, 48 km inland. Of the 19 places where Coombs & Hawke (1992) measured 65 roofed rille segments, 12 lie under the
   sea, 1 on the coast and 6 on dry land; the three collapse chains Sauro et al. (2020) measured, two in the Marius
   Hills and one by Gruithuisen, lie 35–481 m under the Near-side Sea; two of the four GRAIL gravity candidates lie
   under it too (*screen*).
2. **On land the rain lifts the water table over every pit floor within years to four centuries, unless the basalt
   drains as freely as Hawaii's young flows.** Solved in two dimensions on the atlas, with the seas and lakes at their
   levels and the ground holding it, the central case, Oahu's 36% of the rain recharging basalt of the Columbia
   Plateau's median conductivity (a sixth of it at lunar gravity), brings the water table to the ground at 4 of the 12
   land places and to 6–330 m below it at the other 8 (*screen*). It stands over all three pit floors and both radar
   floors, the Tranquillitatis floor by 40–90 m, over every floor estimate of rilles A, B, I and N, and over the deeper
   estimates of rilles G and K (*screen*). Rising from the sea's level when the rain begins, it reaches the Marius Hills
   floor in 2–5 years and the Tranquillitatis floor in 165–412, for specific yields of 0.04–0.10; dry basalt beneath
   adds up to 420 years (*screen*). In permeable basalt with the GCM's coastal cells evaporating more than they rain, it
   stays within 90 m of the sea's level: the Marius Hills floor then stands 0.5–15 m above it and the Tranquillitatis
   floor 520–570 m, while the deeper rille floors still reach below. In tight, thin basalt it stands at the ground at
   all 12 places; rilles G, K and L look 10–197 m dry only because their points sit on rises above their cells' mean
   ground.
3. **Water presses on a lunar tube with a sixth of Earth's weight, so Earth's sealing practice holds tubes far deeper.**
   100 m of water presses 162 kPa, as 16.6 m does on Earth, and 2 km of sea 3.25 MPa, as 331 m does, deeper than the
   Seikan Tunnel's 240 m (*screen*). Earth's tunnel guidelines call for waterproofing sheets where water stands more
   than 15 m above the crown (ITA 1988); the same pressure stands 91 m deep on the Moon (*derived*).
4. **Keeping a large tube dry takes megawatts.** The author's 1 km by 50 km tube under a 25 m roof, with the water table
   at the ground, takes 0.17–17 m³/s of seepage through a grouted ring of Earth practice, pumped for 0.03–2.9 MW; under
   500 m of sea 1.8–184 m³/s and 2.5–245 MW; under 2 km of sea 33–3,300 MW (*screen*). Its 2.6–4.7 million residents
   would use 5.2–28 GW (*derived*). Above the water table only the percolating rain arrives, 0.08–1.6 m³/s pumped for
   0.01–0.28 MW. The tightest ring passes 0.14 L/m² a day, above the 0.04–0.08 L/m² a day Earth allows into tunnels and
   underground public spaces (NASEM 2024), so the living space sits inside membranes, with the seepage drained behind
   them.
5. **The estimate's arithmetic holds; as a vault the same tube holds 2.6–4.7 million.** 1 km by 50 km by ten floors is
   500 km², 7.1 million people at 70 m² (*derived*). A tube 50 m high is a vault: at 4–5 m storeys it holds 9–12 levels
   narrowing upward, 363–467 km² of floor, and with half to 70% of that for homes, work and services, 2.6–4.7 million
   people (*screen*). A kilometre's span is stable under a roof of 2 m in Blair et al.'s (2017) models and of 100 m in
   Theinat et al.'s (2020). The mapped candidates are narrower and shorter: roofed rille segments 0.03–1.3 km wide and
   up to 17.6 km long (Coombs & Hawke 1992), tubes of 170–860 m inferred from collapse chains whose intact lengths run
   4–27 km (Sauro et al. 2020), and a radar-sounded tube about 30 km long in the Marius Hills (Kaku et al. 2017, as
   Naito et al. 2020 report it).
6. **The mapped tubes hold 2–52 million people, at most 1.3% of the Moon's 4–8 billion, and those on dry land 0.6–14
   million.** The rille segments, collapse chains and radar tube, each tube counted once, hold 2.3–4.1 million at a
   third of their rilles' size and 29–52 million at full size, at 70 m² a person (*screen*): 0.58–1.0 and 8.1–14 million
   on dry land, 0.81–1.4 and 9.0–16 million under the sea, and 0.93–1.7 and 12–22 million in places that straddle the
   coast, the Marius Hills tube among them. If the four GRAIL candidates are tubes of the size gravity can see, several
   kilometres wide and tens of kilometres long, each holds 0.70–20 million on ten levels, and all tubes together 5–129
   million, 1.3–34 million of them on dry land (*screen*).
7. **Before the air, a tube roof is the Moon's best shelter.** The airless surface takes 160–416 mSv a year of cosmic
   rays and up to 2,190 mSv in one solar storm (Naito et al. 2020). Under 6 m of rock no cosmic-ray dose is observable
   (De Angelis et al. 2002), and inside a tube 50–75 m from its skylight the dose falls below 1 mSv a year (Naito et al.
   2020), against the Open Moon's limit of 0.027 mSv a day. A 5 m roof weighs 1.5 times Earth's atmosphere per square
   metre and a 20 m roof six times (*screen*). Holding 1 atm of air by its weight takes 19–21 m of basalt and 0.5 atm
   9–10 m (*screen*); thinner roofs hold their air in anchored liners or in habitats set inside. The rock stays near −20
   °C deep in a tube (Hörz 1985) and near 17 °C in the shade beyond a skylight (Horvath et al. 2022).
8. **The air ends that role early, and the seas end it for most tubes.** Once the surface pressure reaches 16.8 kPa,
   about a seventh of the design's 121.6 kPa, the air weighs as much per square metre as Earth's whole atmosphere
   (*screen*), which holds Earth's cosmic-ray dose at sea level to 0.39 mSv a year. The tubes matter through the
   construction order's first steps, power and processing, protection and extraction while the surface is an industrial
   site, until the gases accumulate; then the seas fill the basins that hold most of them.
9. **A lasting tube city brings its daylight down, pushes its air through and never lets its pumps stop.** Heliostats
   feeding light guides need mirrors over 1.4–4.5% of the lit floor at the equator and 2.1–7.1% at 45° for 250–500 lux
   through the 354-hour day (*screen*); at night the tube lights itself as every home does, 11–45 kWh a person (daily
   life). A city of 7 million moves 35,000–70,000 m³/s of outdoor air for 18–117 MW of fans (*screen*), since the stack
   effect is a fifth of Earth's: 3.9 Pa across a 100 m shaft with a 5 K difference, against 21 Pa. If the pumps stop,
   the lowest level fills in 3.8 years under waterlogged ground, 4 months under 500 m of sea and 34 days under 2 km
   (*screen*). About five shallow moonquakes a year shake the whole Moon, the largest above magnitude 5, and a site sees
   one above 4.5 within 100 km about once in 400 years (Oberst & Nakamura 1992).

### Highland and tall cities

10. **3.90 million km² of dry ground stands 4 km or more above the sea, 91% of it on the far side; 1.35 million km² of
    it lies poleward of 30°, where both climate models put the most comfort.** It falls from 1.77 million km² at 4–5
    km to 2,800 km² above 10 km (*screen*). The ways-of-living study gives 1.99 million km² of it to highland towns and
    terrain-held districts once the metropolis, the shores and the wet margins have taken theirs.
11. **The comfortable highlands are the dry ones, and the far side's equatorial massif round the summit is the
    wettest and stormiest ground on the Moon.** Rain over the high ground falls from 7.6 mm a day within 15° of the
    equator to 2.1 mm at 30–45° and 0.34 mm at 45–60°, and rises with height to 12 mm above 9 km; lightning runs from
    0.044 flashes per km² a year above 10 km to 0.0073 at 30–45° and 0.001 at 45–60° (*screen*, from the design
    climatology and the lightning map). Within 100 km of the summit the dry land takes 13.3 mm a day and 0.048
    flashes; the wettest and stormiest cell, 19.8 mm and 0.072 flashes, lies 7.5 km up at 2.6° S, 208.6° E, 327 km
    from the summit (*screen*). The summit port takes 0.046 (0.011–0.083), the highland box at 44.7° S 0.0019 (main's
    places).
12. **Both models find the high ground poleward of 30° the most comfortable, and CM1 finds comfort rising with
    height.** On ground 4–9 km up poleward of 30°, CM1's rings corrected for height give 141–351 comfortable hours of
    the 709-hour lunar day and the GCM 550–695; equatorward of 30°, 16–54 against 229–425 (strict to loose; climate
    products). The 3-D box over the highlands finds 82 strict hours below 1 km rising to 262 at 5–6 km. The GCM's land
    air at the highland bands' mean heights is 18.8–19.9 °C, and 14.8–17.6 °C with the register's 1 °C a kilometre.
13. **For the same steel per floor, towers rise three to four times as high at a sixth of the gravity.** In the
    engineering lattice model sized for towers of 50 m floor plates, a 30-storey frame takes 2.2 kg of steel per m² of
    floor on the highlands against 10.7 on Earth, a 100-storey 9.0 against 40.6 and a 400-storey 59 against 222, so the
    frame of Earth's 120–800 m towers reaches 2.6–3.8 times as high (*screen*); the summit tower study found 2–2.5 times
    at 10–30 km. The design gust presses 0.44 times Earth's design wind against a weight of 0.17, so the wind takes
    26–56% of a leeward leg's force (11–25% on Earth) and towers of 400 storeys anchor their windward legs (*screen*).
14. **Heat and daylight cap a district's density.** Towers on a fifth of the ground at 30 storeys house 85,700 people
    per km² at 70 m² of floor, the metropolis brief's tower district, and release 129–387 W/m²; at 100 storeys they
    house 286,000 and release 431–1,290 W/m², beyond Earth's hottest city cells at 100–577 W/m² (*screen*;
    provisioning). Between towers 50 m wide on that grid, the street sees at least 25% of the sky at 30 storeys and 8%
    at 100 (*screen*).
15. **2–6 billion on the ground take a small share of the high ground at city densities.** All 6 billion in the
    highlands at 10,000 per km² settle 600,000 km², 15% of the high ground and 22% of it outside a kept 30%; at the
    metropolis's 40,000 per km² 150,000 km² (3.8%), at the tower district's 70,000 km² (1.8%) (*screen*). Within 100 km
    of the summit, 18,700 km² of dry land would hold 0.75 billion at the metropolis's density, and within 200 km 2.6
    billion; lakes cover 41–48% of that ground (*screen*).
16. **Water reaches the highlands for 0.9–2.7 GW per billion people from the nearest lakes and 7–21 GW from the sea;
    recovering 90% cuts each tenfold.** The median highland home stands 0.66 km above its nearest standing water and 5.1
    km above the sea (*screen*). A cubic metre lifted 5 km takes 3.2 kWh, against the 2.4 kWh California's State Water
    Project spends bringing water to Southern California over the Tehachapis' 600 m lift (*derived*). Within 30° of the
    equator the high ground's own runoff equals the domestic water of 89–268 billion people, at 30–45° of 1.5–4.4
    billion, and poleward of 45° almost none, so the comfortable dry highlands there lift their water a median 4.5–4.7
    km from a median 64–76 km away, and a tenth of them from beyond 136–146 km (*screen*).
17. **Sky boats and ferries make the height ordinary.** A four-seat sky boat climbs 6 km for 2.3 kWh, as much as 83
    km of cruise, so with that climb its 150 km range becomes 67 km; a 250 m ferry line carries 24,000–36,000 people an
    hour each way at 11.4 Wh a passenger-km; a regional ship covers 300 km in 2.4 hours for 7 kWh a passenger
    (*screen*, from the sky-fleet product).
18. **Half to three quarters of the ground's people in the highlands leaves 0.5–3 billion for the scattered
    communities, a light scatter over 21.4 million km².** Spread over the other settings' homes in proportion to land,
    they live at 23–140 people per km² of home, and at 10,000 per km² their towns cover 0.2–1.4% of it (*screen*).

## The lava tubes

### The places

A tube city on dry ground is a vaulted avenue under tens to hundreds of metres of stone: a long nave lit by mirrored
shafts and skylight gardens, its levels stepping back as the vault narrows, quiet and at one temperature through the
lunar day and night. Before the air it is the Moon's safest home; after it, a city of a few million among the Moon's
billions, with the sea or the water table never far below.

### The candidates

The [inventory](tubes.py) gathers five kinds of evidence (sources.json gives each source and how far it was read):

| Kind | What it rests on | Candidates | Sizes |
|---|---|---|---|
| Mare pits | LROC NAC imaging, shadows and stereo (LROC pit atlas; Wagner & Robinson 2014; Ashley et al. 2011; Haruyama et al. 2009, 2010) | 16 mare pits (two nested) and two at the maria's edge | mouths 12–175 m, inner depths 5–105 m; voids seen 7–20 m beyond the walls, the Marius Hills void 17 m high and the Tranquillitatis opening about 60 m |
| Roofed rille segments | Lunar Orbiter and Apollo photographs (Coombs & Hawke 1992) | 65 segments on 18 rilles (19 places), 9 rated prime | 0.13–17.6 km long, 0.03–1.3 km wide, minimum roofs 15–525 m |
| Collapse chains | Kaguya-LOLA terrain at 59 m (Sauro et al. 2020) | Marius Hills South and North, Gruithuisen | collapses 527–858 m wide and 87–190 m deep; 4.1–27.2 km intact; tubes up to 800 m wide and 190 m high, possibly a half to a third of that |
| Gravity | GRAIL gradiometry and cross-correlation (Chappaz et al. 2014, 2017; Sood et al. 2016; Blair et al. 2015) | Marius Hills, south of Rima Sharp, Sinus Iridum, Vallis Schröteri | detectable only at widths of several kilometres, heights of hundreds of metres and lengths of tens of kilometres; 1–4 km wide at Schröteri and Sharp |
| Radar | SELENE sounder (Kaku et al. 2017, 2018); LRO Mini-RF (Carrer et al. 2024) | Marius Hills rille A, west of Rima Mairan, the Tranquillitatis conduit | Marius Hills about 30 km long at 75–150 m depth (as cited by Naito and Sauro); Tranquillitatis tens of metres long |

The named rilles the literature ties to tubes are placed at their IAU centres from the gazetteer archive the geography
domain holds: Vallis Schröteri, Rima Sharp, Rima Marius, Rima Hadley, Rimae Prinz and Rimae Sulpicius Gallus.

### Roofs

Blair et al. (2017) find tubes stable to about 5 km wide in rock laid down in one thick layer and 3.5 km in thin layers
under 50–100 m roofs, and a kilometre's span stable under a 2 m roof; Theinat et al. (2020), with a fractured rock
mass, need a 100 m roof for a kilometre and 1 km for 4 km (both read from Williams 2021's reproductions of their
maps). All 54 rille segments with a measured roof lie within Blair et al.'s stable spans; 47 meet Theinat et al.'s
thinnest roof, and seven fall short, among them four of the nine prime segments south of Gruithuisen, 720–970 m wide
under 30–46 m (*screen*).

### Against the atlas

Each candidate's ground is LOLA's 16 px/deg height less the degree-200 geoid, set against the atlas's sea level, −1,654
m above the geoid, and the drainage estimate's rain-fed lakes. The atlas fills every closed hollow below that level, so
a candidate whose own ground lies below it is under water: under the sea where its atlas cell or one of the eight
around it is sea or an atlas lake. Rille H lies so, 23 m below the level just behind the atlas's coast. Under a lake
means a rain-fed lake covers it; coastal means its atlas cell or one of the eight around it, 7.6 km across, holds
standing water; dry means farther. The Coombs & Hawke positions are given to the half degree, so the product adds the
sea's share of a box a degree across: their dry and coastal places lie in boxes 0–88% under the sea (*screen*).

| Kind | Under the sea | Under a lake | Coastal | Dry |
|---|---:|---:|---:|---:|
| Mare and edge pits | 13 | 1 | 1 | 2 |
| Rille places | 12 | 0 | 1 | 6 |
| Collapse chains | 3 | 0 | 0 | 0 |
| GRAIL candidates | 2 | 0 | 0 | 2 |
| Radar voids | 1 | 0 | 0 | 2 |
| Named rilles at their centres | 4 | 0 | 0 | 2 |

Over the submerged candidates stand 54 m of water at the Southwest Fecunditatis pit and up to 2.25 km at the Runge pit
in the Smythii–Marginis Sea; the Mare Ingenii pit lies under 1.9 km of the South Pole–Aitken Sea, Lacus Mortis under
372 m of the Near-side Sea, and rille P under 623 m of a crater lake at the common level. A region whose candidates
lie both on land and under water, or whose box lies 10–90% under the sea, straddles the coast: the Marius Hills, its
pit 60 m up while its chains and rille J lie 35–169 m down and its gravity and radar boxes 36% and 51% under the sea,
and rilles D, E1, H, L and N (*screen*).

### The water table

The seas are fresh, so no salt wedge holds the groundwater down: the water table meets each sea or lake at its level and
rises inland under the recharge. [groundwater.py](groundwater.py) solves it in two dimensions on the atlas within 400 km
of each land candidate: unconfined Dupuit flow over a flat base b below the window's lowest water, K∇²(H²/2) + R = 0 for
the saturated thickness H (Woessner & Poeter 2020), with the seas and lakes held at their levels, no flow across the
window's edge, and the ground holding the water table wherever the mound would rise above it, as springs and wet ground.
The solve reproduces the strip's analytic mound, (h + b)² = b² + (R/K)·x(2L − x), within 10⁻⁶ m, and the water table at
every land place moves by less than 0.1 m between windows reaching 200 and 600 km out (*screen*). Hydraulic conductivity
is K = kρg/μ, so the same rock passes a sixth of Earth's water at lunar gravity (Woessner & Poeter 2020; Clifford 1993):
a mound small against the saturated thickness beneath it stands 6.0 times as high as Earth's for the same rain and one
that dwarfs it 2.5 times; at the land places the rain's mound stands 4.2–6.0 times as high (*screen*). Three cases:

| Case | Recharge | Earth conductivity | Basalt below the lowest water | Specific yield |
|---|---|---|---|---|
| Low | a tenth of the rain, at most the GCM cell's rain less evaporation | 1.8×10⁻² m/s, Hawaii's dike-free flows at the top of their usual range (Hunt 1996) | 1,620 m | 0.04–0.10 |
| Central | 36% of the rain, Oahu's island-wide share (Engott et al. 2017) | 2.5×10⁻⁴ m/s, the Columbia Plateau's median for its basalt units (Kahle et al. 2011) | 740 m | 0.04–0.10 |
| High | half the rain, the top of Hawaii's 10–50% (Gingerich & Oki 2000) | 3×10⁻⁷ m/s, the Columbia River basalt's bulk low end (Burns et al. 2016) | 100 m | 0.04–0.10 |

The basalt thicknesses are the nearside mean and its one-sigma bounds from GRAIL (Gong et al. 2016), and the specific
yields Hawaii's (Hunt 1996). Recharge never exceeds the rock's conductivity. The rain is the GCM's, carried onto the
atlas; its cells, 170 km across, mix land and sea, and at most mare candidates they evaporate more than they rain, so
the low case gives them no recharge. Gong et al. measure the basalt down from the mare surface, while the solve sets
the base that far below the lowest water; under high ground such as the Tranquillitatis pit's the basalt may end above
that base, and the saturated column then runs on into the crust. The water table is the atlas cell's, 7.6 km across,
set against each candidate's own ground at 16 px/deg; a candidate on a rise within its cell stands farther above the
water table, one in a valley nearer to it.

The rise starts when the rain begins, from the water table at the sea's level under the land, by the Boussinesq equation
S∂h/∂t = ∇·(KH∇h) + R stepped implicitly with the ground holding it in every step. In a strip it keeps within 4.3% of a
fully nonlinear one-dimensional run up to 95% of its rise, and it settles on the steady solve (*screen*). Its implicit
steps grow to 25–29 years over the centuries, which runs the rise times 3–5% long: the independent review gets 274 years
to the Tranquillitatis floor where the product gives 289. Its times scale with the specific yield; the high case's tight
rock barely drains sideways, so its times are the recharge's first rate, the recharge over the specific yield. From dry
basalt the water first fills the rock below the sea's level, in the central case at most 40–420 years more (*screen*):

| Place | Ground above the sea | Floor below the ground | Water table below the ground: low, central, high | Central: floor under water, 90% of the rise |
|---|---:|---:|---|---|
| Marius Hills pit | 60 m | 45–60 m | 60 m, at the ground, at the ground | 2–5 years, 17–42 years |
| Mare Tranquillitatis pit | 784 m | 125–175 m | 694 m, 85 m, at the ground | 165–412 years, 149–373 years |
| Schlüter pit (coastal) | 157 m | 57–75 m | 86 m, 51 m, at the ground | 7–17 years, 6–14 years |
| Rille A | 199 m | 126–231 m | 126 m, 6 m, at the ground | 10–26 years, 40–99 years |
| Rille B | 121 m | 58–397 m | 71 m, at the ground, at the ground | 9–23 years, 19–46 years |
| Rille G | 295 m | 112–650 m | 294 m, 175 m, 13 m | the deeper floors only; 49–124 years |
| Rille I | 80 m | 138–359 m | 59 m, at the ground, at the ground | from the start; 10–24 years |
| Rille K | 534 m | 72–752 m | 533 m, 91 m, 10 m | the deeper floors only; 80–200 years |
| Rille L | 423 m | – | 389 m, 329 m, 197 m | –; 11–28 years |
| Rille N (coastal) | 19 m | 204–737 m | 19 m, 12 m, at the ground | from the start; 2–5 years |
| Vallis Schröteri (centre) | 95 m | – | 95 m, at the ground, at the ground | –; 74–184 years |
| Rimae Prinz (centre) | 260 m | – | 260 m, 145 m, at the ground | –; 51–128 years |

The Marius Hills radar tube, 75–150 m down, and the Tranquillitatis radar conduit, 130–175 m down, share their pits'
water tables: the first lies below the sea's level from the start, the second goes under in 162–404 years. The gravity
candidates at the Marius Hills and Vallis Schröteri share those places' water tables. "From the start" means the floor
lies below the sea's level, under water as soon as the basalt beneath is full. At the Schlüter pit the central water
table stands only 6–24 m over the floor, and with the neighbouring rain-fed lake treated as land the shallowest floor
sits 0.9 m above it. In the high case a rille floor goes under at the first rate in 6–99 years and the Tranquillitatis
floor in 107–268 (*screen*). The crust beneath the seas fills from below as well, which the geography domain's
groundwater estimate leaves to time and permeability
([geography](../../../geography/README.md#groundwater-a-first-estimate)); while it fills, the water table near the shore
can stand below the sea's level.

### Keeping a tube dry

Water at depth presses ρgh at lunar gravity, 1.62 kPa a metre (*screen*):

| Head | Pressure | The Earth head of the same pressure |
|---|---:|---:|
| 15 m | 24 kPa | 2.5 m |
| 100 m | 162 kPa | 16.6 m |
| 500 m | 812 kPa | 83 m |
| 2,000 m | 3.25 MPa | 331 m |

Earth's practice carries over at those pressures. Norway's subsea tunnels, down to 287 m below the sea, are drained,
pre-grouted against 2–3 MPa and pumped, with a sump for at least 24 hours of the allowed 300 litres a minute per km
(NFF 2009); grouted Oslo tunnels pass about 2 litres a minute per 100 m, ungrouted ones 10–40 (NFF 2011). The Seikan
Tunnel's three pump stations drain about 10 million m³ a year, with 84 hours' margin if every pump stops (Takai 2012).
A watertight lining carries the full head, and the waterproofing sheets Earth uses above 15 m of water stand at 91 m
on the Moon. A tunnel-waterproofing guide serves 5–30 m of head with membrane systems, 30–60 m with membranes and
watertight concrete and more with double, actively controlled membranes (Sika); on the Moon the same pressures stand
at 30–181 m, 181–362 m and beyond 362 m (*screen*). In Norway's subsea road tunnels lighting and ventilation
dominate the running costs and pumping costs little (NFF 2001).

The study gives the author's tube, 1 km wide, 50 m high and 50 km long under a 25 m roof, a grouted ring 5 m thick of
Earth practice's 10⁻⁹–10⁻⁷ m/s (Katibeh & Aalianvari 2012 put the practical minimum near 10⁻⁸), passing a sixth of
that at lunar gravity, and pumps of the inhabited-volume product's 70% (*screen*):

| Setting | Mean head | Seepage | Through the envelope | Pumping |
|---|---:|---:|---:|---:|
| Above the water table | – | 0.08–1.6 m³/s of percolating rain | – | 0.01–0.28 MW |
| The water table at the ground | 50 m | 0.17–17 m³/s | 0.14–14 L/m² a day | 0.03–2.9 MW |
| Under 500 m of sea | 550 m | 1.8–184 m³/s | 1.6–157 L/m² a day | 2.5–245 MW |
| Under 2,000 m of sea | 2,050 m | 6.8–684 m³/s | 5.9–585 L/m² a day | 33–3,300 MW |

Even the tightest ring passes more than the 0.04–0.08 L/m² a day allowed into underground public spaces and tunnels
(NASEM 2024), so the city lives inside membranes and linings, with the seepage gathered behind them and pumped.

### What a tube holds

A level counts where a semi-elliptical vault over a flat floor leaves a full storey of headroom. With storeys of 4–5 m
and half to 70% of the floor for homes, work and services (the rest streets, courts and light wells, a design guess),
the author's tube holds 363–467 km² of floor and 2.6–4.7 million people at the metropolis brief's 70 m² a person,
1.3–2.3 million at twice that floor and 5.2–9.3 million at half of it (*screen*). Its levels add up to the estimate's
ten full-width floors once the vault stands 54–67 m high (*screen*).

| Candidate | Small tube: people at 70 m² | Large tube |
|---|---:|---:|
| The 65 roofed rille segments, 168 km | 1.7–3.0 million | 24–42 million |
| Gruithuisen chain, 27 km intact | 0.53–0.93 million | 4.9–8.6 million |
| Marius Hills chains, 10.5 km intact | 0.06–0.11 million | 0.57–1.0 million |
| Marius Hills radar tube, 30 km | 0.16–0.30 million | 1.6–2.9 million |
| Each GRAIL candidate, 1 km by 10 km to 4 km by 50 km, ten levels | 0.70–0.99 million | 14–20 million |

The small tube is a third of the rille's or collapses' width and height, the large tube their full size (Sauro et al.
2020). A GRAIL candidate's ten levels are the lowest in a vault as high as the collapses' depth-to-width ratio of
0.15–0.22 makes it, 150–880 m. The totals count the Marius Hills rille, its chains, radar tube and gravity low as one
tube, the Gruithuisen rille and chain as one and the Rima Sharp–Mairan gravity and radar candidates as one; rille K
lies on the edge of the Marius Hills gravity box and counts as its own tube, 0.40–0.70 and 5.6–9.8 million that leave
the totals with gravity if it is the same (*screen*).

### Before the air

The airless surface takes 160 mSv a year of cosmic rays at solar maximum and 416 at minimum, and up to 2,190 mSv in a
single solar particle event (Naito et al. 2020; Chang'E-4 measured 1,369 µSv a day, Zhang et al. 2020). The dose under
rock falls from 270 mSv a year at 1 m to 120 at 3 m and 20 at 5 m, and none is observable past 6 m (De Angelis et al.
2002, read from their figure for a 5 m layer of loose surface rock over basalt); the floor of a 43 m pit takes under 30
mSv a year, and the tube beyond 50–75 m from the skylight under 1 mSv (Naito et al. 2020). Rille roofs of 15–525 m
and the 20–25 m lava sheet over the Marius Hills pit (Haruyama et al. 2009) stand far past the 6 m.

| Roof of basalt | Mass per m² | Earth atmospheres |
|---|---:|---:|
| 5 m | 15–16 t | 1.5–1.6 |
| 10 m | 30–33 t | 2.9–3.2 |
| 20 m | 60–65 t | 5.8–6.3 |
| 100 m | 301–327 t | 29–32 |

At 3,010–3,270 kg/m³ of mare basalt (Kiefer et al. 2012), a roof holds 1 atm by its weight from 19–21 m and 0.5 atm
from 9–10 m (*screen*). A tube stays near −20 °C deep inside (Hörz 1985); beyond a skylight the shaded floor holds
about 290 K, near 17 °C, through the lunar day (Horvath et al. 2022), and the holes' albedo gave 4–19 °C (Haruyama et
al. 2010). Under the air the land's mean settles near the GCM's 22–23 °C at the sea's level (climate products), some
40 K above the tubes' rock, which follows by conduction: half the change reaches 20 m in 14–28 years, 50 m in 87–174
and 100 m in 350–700, for a thermal diffusivity of 0.5–1×10⁻⁶ m²/s (a design range).

Once the Open Moon's air presses 16.8 kPa at the surface, its column matches the 10.3 t/m² of Earth's whole
atmosphere (*screen*); from then on the open ground is shielded about as Earth's sea level is, and the tubes' shelter
matters most in the construction order's first steps: power and processing, protection, and extraction while the
surface is an industrial site
([industrial architecture](../../../engineering/reference/industrial_architecture/RECOVERED_INDUSTRIAL_ARCHITECTURE.md),
section 9). The build's workers live there and in the array's habitats until the gases accumulate
([population](../../../provisioning/population.md)).

### A lasting city

- **Daylight.** Heliostats track the Sun and send its beam down shafts into light guides. The beam reaching the mirror
  averages 36,800 lux over the equator's day and 23,500 at 45° (*screen*, from the haze product's clear-sky
  transmission and the calendar's sunlight above the air). For 250–500 lux on the floor with a third to a half of the
  light delivered (a design guess), mirrors cover 1.4–4.5% of the lit floor at the equator and 2.1–7.1% at 45°; the
  author's tube, with 7.3–9.3 times its footprint in floor, needs mirrors over a tenth to two-thirds of its footprint
  (*derived*). Through the night the tube lights itself as every home does, 11–45 kWh a person (daily life).
- **Air.** At 5–10 L/s of outdoor air a person, from ASHRAE 62.1's office rate with its floor allowance to England's
  Approved Document F, a city of a million moves 5,000–10,000 m³/s and one of 7 million 35,000–70,000 m³/s, for 2.5–17
  and 18–117 MW of fans at 300–1,000 Pa (a design range) (*screen*). The stack effect across a 100 m shaft with a 5 K
  difference is 3.9 Pa in the Open Moon's air at sea level, against 21 Pa on Earth (Straube 2007 gives its form), so
  the air moves by fans.
- **Floods.** Storms rain 17–87 mm an hour beneath their flashing cores (infrastructure review), so a pit's funnel
  150–300 m across takes 0.4–1.7 m³/s at the peak (*screen*) and its rim is raised. With the central ring, the lowest
  level of the author's tube fills in 3.8 years if every pump stops under waterlogged ground, 4 months under 500 m of
  sea and 34 days under 2 km; a breach of 1 m² admits 8–48 m³/s (*screen*).
- **Quakes.** About five shallow moonquakes a year shake the whole Moon (Nakamura et al. 1979), the largest above
  magnitude 5, and a site sees one above 4.5 within 100 km about once in 400 years, a rate like the central United
  States' (Oberst & Nakamura 1992). Eight of the 28 recorded lie within 30 km of young thrust-fault scarps (Watters et
  al. 2019), so tube cities keep clear of them. The tubes have stood since their maria formed, 3.2–3.7 billion years
  ago (Haruyama et al. 2010).

### How many, and when

The mapped tubes, each counted once, hold 2.3–4.1 million people at a third of their rilles' size and 29–52 million at
full size: 0.58–1.0 and 8.1–14 million on dry land, 0.93–1.7 and 12–22 million in places that straddle the coast, and
the rest under the sea (*screen*). With the four GRAIL candidates at the size gravity can see, 5–129 million, 1.3–34
million on dry land and 1.5–39 million straddling the coast. Against the Moon's 4–8 billion that is at most 3.2%, the
dry land's share at most 0.86% and 1.8% with the straddling places. The tubes' time is the build's first steps, before
the rain; once it falls, the central water table stands over every land pit's floor within years to four centuries,
so a tube kept as a city is sealed and pumped from then on. After the air and water arrive, those on dry ground become
cities among many, with pumps and membranes, and those under the sea hold enclosed communities like the undersea domes
of the heritage sites (conservation D6).

## The highlands

### The places

Highland cities stand in cool, bright air on the far side's high plateaus: terraces on crater walls above high lakes,
tower districts on ridges with gardens on their roofs, and the summit metropolis round its port, in the massif where
the rains are heaviest. Poleward of 30° the high ground is the Moon's most comfortable land and its driest.

### The high ground

From the 28% atlas less the rain-fed lakes (*screen*):

| Height above the sea | Dry land | Far side | Rain | Flashes per km² a year | 0–15° | 15–30° | 30–45° | 45–60° | 60–90° |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4–5 km | 1.77 million km² | 86% | 3.0 mm/day | 0.011 | 507,000 | 412,000 | 493,000 | 221,000 | 140,000 |
| 5–6 km | 1.02 million km² | 92% | 4.0 | 0.014 | 360,000 | 286,000 | 240,000 | 82,000 | 53,000 |
| 6–7 km | 625,000 km² | 96% | 5.9 | 0.021 | 302,000 | 214,000 | 76,000 | 15,000 | 18,000 |
| 7–8 km | 304,000 km² | 99% | 8.3 | 0.030 | 195,000 | 93,000 | 11,000 | 800 | 3,900 |
| 8–9 km | 141,000 km² | 100% | 10.8 | 0.039 | 118,000 | 22,000 | 200 | 40 | 300 |
| 9–10 km | 31,000 km² | 100% | 12.0 | 0.043 | 29,000 | 2,400 | – | – | – |
| above 10 km | 2,800 km² | 100% | 12.3 | 0.044 | 2,700 | 50 | – | – | – |

The highlands' runoff is 4,180 km³ a year within 15° of the equator, 717 at 15–30°, 81 at 30–45° and almost none
poleward (*screen*). Around the summit, dry land within 32, 100, 200 and 300 km covers 2,464, 18,700, 65,200 and
122,800 km², 9.5–7.8 km up, under 13.1–10.6 mm of rain a day, with lakes on 25–56% of the ground (*screen*).

### The air and its comfort

The climate products give the comfort of the land air by latitude and height, CM1's flat rings corrected for height
against the GCM's 3-day means (comfortable: air at 2 m of 18–26 °C with a dewpoint of at most 15 °C, strict, or 17–27 °C
and 16 °C, loose; hours of the 709-hour lunar day):

| Latitude | Height | CM1 strict–loose | GCM strict–loose |
|---|---|---:|---:|
| 0–30° | 0–0.5 km | 0–1 | 99–165 |
| 0–30° | 4–9 km | 16–54 | 229–425 |
| 30–60° | 0–0.5 km | 33–64 | 510–613 |
| 30–60° | 4–9 km | 141–248 | 550–649 |
| 60–90° | 4–9 km | 249–351 | 611–695 |

The 3-D CM1 box over the highlands at 44.7° S finds 82 strict hours below 1 km, 244 at 4–5 km, 262 at 5–6 km and 205
at 6–7 km, its days at 28–35 °C and nights at 17–19 °C (climate products). Both models put the most comfortable land
on high ground poleward of 30°, and the register carries comfort as the range between them.

### How tall

The engineering lattice model (four corner legs and braced faces sized at S690 steel's allowable stress, with one
member-buckling and one gust factor) is sized here for towers of 2,500 m² floor plates, storeys every 4 m and floors of
500 kg/m², the frame inside a closed facade. Each tower is the lightest frame that sways at most 1/500 of its height in
the service gust, keeps a buckling factor of at least 3 and stands no wider than half its height. The Moon's towers
stand on ground 6 km up in the design run's air with the summit tower study's 32.4 m/s design gust and 22.5 m/s service
gust; Earth's in the standard atmosphere with 50 m/s and 34.7 m/s (*screen*):

| Storeys | Height | Moon: frame steel per m² of floor | Earth | Moon: period, sway | Wind's share of the leeward leg, Moon / Earth | People per km² at a fifth of the ground | Heat at 1.51–4.52 kW a person |
|---|---:|---:|---:|---:|---:|---:|---:|
| 30 | 120 m | 2.2 kg | 10.7 kg | 2.8 s, 0.03 m | 26% / 11% | 85,700 | 129–387 W/m² |
| 60 | 240 m | 4.9 kg | 22.8 kg | 6.1 s, 0.14 m | 33% / 15% | 171,000 | 259–775 W/m² |
| 100 | 400 m | 9.0 kg | 40.6 kg | 9.9 s, 0.38 m | 36% / 19% | 286,000 | 431–1,290 W/m² |
| 200 | 800 m | 21 kg | 92 kg | 19 s, 1.5 m | 41% / 23% | 571,000 | 863–2,580 W/m² |
| 400 | 1.6 km | 59 kg | 222 kg | 21 s, 3.1 m | 52% / 24% | 1.14 million | 1,730–5,170 W/m² |

The frame here counts the legs and bracing only; the floors' own structure is in their 500 kg/m², and the
inhabited-volume budget gives the district's building 400 kg/m², 28 t a person at 70 m². On Earth the model's frame
weighs a fifth to two-fifths of the lateral steel of designed tube and diagrid towers, 17.8 kg/m² at 48 storeys against
their 44–73 and 38.7 at 96 against 98–215 (*screen*; Kheyroddin et al. 2018, after Khan's premium for height); its
Moon-to-Earth ratio is the result it carries. Steel's self-weight length is 25 km on the Moon and 4 km on Earth (summit
tower study). In the dense air the design gust presses 674 Pa at 6 km against Earth's 1,531 Pa, 0.44 of it, while weight
is 0.17 of Earth's, so wind leads the structure; from 400 storeys the windward legs pull in tension and anchor into the
ground (*screen*). The street between towers 50 m wide on a fifth of the ground is 62 m wide: at 30 storeys it sees at
least 25% of the sky from its middle, at 100 storeys at least 8% (*screen*). These are the view factors of an infinite
street canyon, cos(atan(2H/W)); the towers stand apart and 55% of each street side opens between them, so the street
sees more.

### People on the high ground

| On the ground | In the highlands | At 10,000 per km² | At 40,000 (the metropolis) | At 85,700 (the tower district) |
|---|---|---:|---:|---:|
| 2 billion | half | 100,000 km², 2.6% | 25,000 km², 0.6% | 11,700 km², 0.3% |
| 4 billion | three quarters | 300,000 km², 7.7% | 75,000 km², 1.9% | 35,000 km², 0.9% |
| 6 billion | all | 600,000 km², 15% | 150,000 km², 3.8% | 70,000 km², 1.8% |

Shares are of the 3.90 million km² of high ground (*screen*); the product also gives them against the ways-of-living
home, the ground poleward of 30° and the land outside kept shares of 30% and half.

### Water

Domestic water is 18–55 km³ a year per billion people (50–150 litres a day; ways of living). Lifting it takes 0.645 kWh
a cubic metre for each kilometre at 70% pump efficiency, 0.074 GW per km³ a year per kilometre (*screen*):

| Latitude | Lift to the nearest water, median and 90th percentile | Distance | Per billion from the nearest water | From the sea |
|---|---:|---:|---:|---:|
| 0–15° | 0.46, 1.5 km | 7.6, 21 km | 0.61–1.8 GW | 7.6–23 GW |
| 15–30° | 0.54, 2.8 km | 14, 47 km | 0.73–2.2 GW | 7.1–21 GW |
| 30–45° | 0.75, 4.8 km | 33, 108 km | 1.0–3.0 GW | 6.5–19 GW |
| 45–60° | 4.5, 5.6 km | 76, 146 km | 6.1–18 GW | 6.2–19 GW |
| 60–90° | 4.7, 6.0 km | 64, 136 km | 6.3–19 GW | 6.3–19 GW |

Recovering half, 90% or 99% of it cuts the lift twice, ten and a hundred times (*screen*). Within 30° the high ground's
own runoff supplies its people many times over; at 30–45° it covers 1.5–4.4 billion people's domestic water, or ten
times that at 90% recovery; poleward of 45° the comfortable high ground draws on the lakes and seas below.

### Getting about

The highland cities sit 4–10 km above the coasts. From the sky-fleet product (*screen*):

| Trip up from the coast | Sky boat, four seats | Sky ferry, 250 m | Regional ship, 350 m |
|---|---|---|---|
| 50 km | 4.2 kWh with the climb, 21 minutes | 0.57 kWh a passenger, 33 minutes | 1.2 kWh a passenger, 24 minutes |
| 100 km | 5.6 kWh, beyond its 67 km reach with the climb | 1.1 kWh, 67 minutes | 2.3 kWh, 48 minutes |
| 300 km | beyond its 150 km range | beyond its 180 km route | 7.0 kWh, 2.4 hours |
| 1,000 km | – | – | 23 kWh, 7.9 hours |

The sky boat's climb of 6 km takes 2.3 kWh at its 80% drive efficiency, a sixth of the same boat's climb on Earth and as
much as 83 km of cruise. Its battery holds 4.6 kWh for 150 km of cruise with its hover allowance of 0.46 kWh, which
includes the fleet's 500 m climb to its flight lane, and a 30% reserve, so with the climb it reaches 67 km; longer trips
up take a ferry or a regional ship. Ferry lines every two to three minutes carry 24,000–36,000 people an hour each way.

### The scattered communities

Beyond the highlands the communities are harbour towns on the seas and great lakes, wetland and canopy towns in the
rain belt, market villages among its farms, nearside towns under the fixed Earth, limb towns where it rises and sets,
the far side's dark-sky towns and the high latitudes' twilit ones, the tube cities on dry ground and the undersea domes
of the heritage sites. The other settings of the [ways-of-living study](../ways_of_living/README.md) cover 21.4 million
km² of home: coasts and great lakes 3.95 million, wet margins 3.58, the rain belt's forest edges and farms 4.79, the
nearside's middle latitudes 1.76, the limb 0.57, the far side's twilit-night places 1.75 and the high latitudes 5.00.
The tubes on dry land hold 0.6–14 million as mapped and up to 34 million with the gravity candidates, the undersea domes
48,600 (*screen*; ways of living):

| On the ground | In the highlands | Scattered | People per km² of home | Towns at 10,000 per km² |
|---|---:|---:|---:|---:|
| 2 billion | 50–75% | 0.5–1.0 billion | 23–47 | 0.2–0.5% of the homes |
| 4 billion | 50–75% | 1.0–2.0 billion | 47–93 | 0.5–0.9% |
| 6 billion | 50–75% | 1.5–3.0 billion | 70–140 | 0.7–1.4% |

## Assumptions

Stated once; the product's `assumptions` block holds them all.

- **Water table:** the three cases above, in uniform rock over a flat base the basalt's thickness below the lowest
  water, solved in two dimensions on the atlas's 7.6 km cells within 400 km of each place, the ground capping it; the
  rise from the sea's level under the land when the rain begins.
- **Sealing:** a grouted ring 5 m thick at Earth's 10⁻⁹–10⁻⁷ m/s, a sixth of that at lunar gravity; pumps at 70%.
- **Capacity (design guesses):** storeys of 4–5 m, half to 70% of the floor for homes, work and services; small and
  large tubes a third and the whole of the rille's or collapses' size; the GRAIL candidates at 1–4 km wide, 10–50 km
  long and the lowest ten levels of a vault 0.15–0.22 times as high as it is wide.
- **Rock:** mare basalt of 3,010–3,270 kg/m³; thermal diffusivity 0.5–1×10⁻⁶ m²/s (a design range).
- **Lasting city (design guesses):** a third to a half of the beam delivered from heliostat to floor; fans at
  300–1,000 Pa and 60%.
- **Towers:** the lattice model with 50 m plates, 4 m storeys and 500 kg/m² floors on ground 6 km up; Earth's service
  gust in the Moon's proportion to its design gust.
- **People:** 10,000, 40,000 and 85,700 per km² of settled land; ground populations of 2, 4 and 6 billion with half,
  three quarters or all in the highlands.

## Inputs and bindings

The run reads the committed products listed in its `INPUTS` and records their hashes: the atlas, drainage and
groundwater products; the CM1 comfort and highland-box products and the design run's global winds; the lightning
occurrence map; the haze product and the sea-appearance lighting calendar; the summit tower, summit tower form, sky
fleet, inhabited volume, ways-of-living and daily-life products. It reads the ignored grid products (the atlas and
drainage grids and the design climatology) and the geography domain's fetched inputs (LOLA's 16 px/deg grid, GRAIL's
GL0420A and the IAU gazetteer archive), and records their hashes when present. Main's corrected per-place product,
`research/studies/joint_synthesis/results/places.json`, is read in place since the joint integration and checked
against its sha256 (`9707e4d1…`).

## Reproduction and checks

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.land_cities.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q research/studies/land_cities
```

The tests check the product's regeneration, the producer, input and bound hashes and the shared constants; pressure,
conductivity and the stack effect against gravity; the mound's limits and the Moon-to-Earth mound ratios; the
two-dimensional solve against the strip's analytic mound and the free boundary of a capped strip, its rise against a
nonlinear one-dimensional run and its settling on the steady solve; the land water tables' cases, caps and floor
verdicts; the vault's levels and the estimate's arithmetic; seepage and pumping; the roof that holds air and the
conduction lag; the beam against the haze product; that every candidate takes one class; the capacity sums by kind and
class with the Marius Hills straddling the coast; the atlas cell on cell edges and a box across 0° E; that the high
ground's bands partition it; that the summit's heights flash most; that the lunar towers are lighter and their
densities consistent; the undersea residents, the sky boat's reach with its climb and the canyon's bound; the people,
water and grid helpers. The tests that regenerate the product skip cleanly without the ignored inputs. The independent
review and its fixes are in [checks.json](checks.json).

What this leaves for later work: the GRAIL candidates' own sizes and the paper tables this study could read only
through secondary sources (sources.json lists them); a groundwater model of the maria with layered basalt, the crust
beneath and the seas filling it from below, at the drainage grid's resolution; a structural model of tube roofs under
water pressure and moonquakes; the highland cities' wind over real terrain; and the distribution of people, which is
the author's to choose in step six.
