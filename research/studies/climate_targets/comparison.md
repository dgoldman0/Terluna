# The three target climates compared

17 October 2026 · `domain/climate`. The author is choosing the Moon's target climate among three, and may carry it as a
range or a distribution across them. The hot world is where the corrected GCM settles, 306.5 K over the 20-year test's
years 10–19 and 306.8–307.1 K at balance by its fitted feedback; the resort is the middle target, set at 301 K within
300–302 K; and the lock is the design temperature, 295 K, the design case itself. The hot world and the lock are the
GCM's output, two runs on two physics, and the resort is interpolated between them by the ground's temperature. Through
every correction of the day, from the heat weighted by area and read at 2 m (d3aff3c) to the model's own land (4916121)
and each mean's humidity (2191b43), their order has held: the hot world is the hottest, wettest and least dimmed, the
lock the coolest, driest and most dimmed, and the resort lies between them in every climate figure; the deep lakes'
overturns under the energetic regime are the one exception (below).

Each figure comes from the studies in this tree since the fourth joint integration of 19 October:
[what each target climate takes](README.md) and its product,
[results/climate_targets.json](results/climate_targets.json), with winter on the tower and the resort's best value (the
*targets*); [where light limits growth](../light_limitation/README.md) (the *light study*);
[the waters' oxygen](../lake_oxygen/README.md) on the corrected lake model, with
[the lakes by kind](../lake_kinds/README.md), [the seas' oxygen](../sea_oxygen/README.md) and
[the waters' methane](../lake_methane/README.md) at 301 K (the *lakes*); [the heat levers](../heat_levers/README.md)
(the *levers*); [atmosphere/rain](../../../atmosphere/rain/README.md) with
[the sound lab](../../../visualization/rain-audio/README.md) (the *rain*); and the [register](../../decisions.md).
[The resort at 301 K across the domains](../resort_301k/README.md) works out where they meet.

## The key figures

| | The hot world, 306.5 K | The resort, 301 K | The lock, 295 K | Evidence (source) |
|---|---|---|---|---|
| Dimming in all; the target's own share of it | 3.7–5.3%; none | 8.2–10.9%; 4.5–5.6% | 12.6–16.6%; 8.9–11.3% | Calculation on the test's feedback (targets) |
| Rain on land a day: every cell's land (the model's land cells) | 3.34 mm (2.94) | 2.82 mm (2.39 at 300 K) | 2.25 mm (1.98) | GCM output; the resort interpolated (targets) |
| Land under 0.5 mm a day, the same two counts | 18.9% (20.1%) | 22.5–28.1% (24.8–31.5%) | 38% (42%) | The same; the resort's range by the interpolation's form (coast water) |
| Land wet-bulb at 2 m: median; 99th percentile | 26.9; 31.5 °C | 22.0; 26.8 °C | 16.6; 21.8 °C | The same |
| Land's time with the 2 m wet-bulb above 28 °C | 31%, and 2.1% above 31 °C | under 1% | none | The same |
| Seas' surface, median | 34.5 °C | 29.5 °C | 24.1 °C | The same |
| Oxygen at saturation, at the seas' median warmth | 7.06 mg/L | 7.68 mg/L | 8.43 mg/L | Calculation (lakes); the resort interpolated (targets) |
| Deep crater lakes, 352 m: overturns in ten years (120 lunar cycles), Osborn's rule; Shih and Bouffard's at 0–15° | 29–30; 12 | 54–68; 8 at 300 K | 62–74; 11 | Lake model on the GCM's weather with the wind's mixing (lakes); the resort run at 301 K |
| Their deepest tenth, oligotrophic water: mean oxygen, Osborn's rule; Shih and Bouffard's | 4.2–4.9; 1.6 mg/L | 5.7–6.3; 3.5 mg/L | 6.8–7.4; 4.8 mg/L | The same |
| Winter on the tower: its length in sub-freezing mean air; the top's mean air; the open winter disk, 13.8 km up | none; 2.9 °C; 10.6 °C | 3.5 km; −3.0 °C; 5.6 °C | 9.9 km; −9.3 °C; 0.2 °C | GCM output over the summit; the resort interpolated (targets) |
| Growth lost through light to the target's own dimming | none | 0.31–1.64% | 0.62–3.3% | Calculation (targets, light study) |
| Margin under a saturated column, the added heat answered | 1.3–1.8 W/m² | 11.9–15.0 W/m² | 22.3–28.6 W/m² | Calculation against the column model's ceiling (targets) |

The dimming in all includes the 3.7–5.3% every target shares to answer the gases, the glow and the beamed power, which
costs each a further 0.26–1.56% of the growth through light and 2.0–7.3% of the land's rain (calculations). The resort's
wet-bulb 99th percentile stays under 28 °C across 300–302 K, by day and in 6.5-hour means, so its land passes 28 °C less
than 1% of the time; the product carries its rain on the model's own land cells at 300 K only. The interpolation moves
every climate figure linearly with the ground's temperature, so a range or a distribution reads from its ends. The
resort's deserts are given as a range because they turn on how each cell's rain is carried between the two runs: by a
fixed amount a kelvin they cover 22.5% of the land (24.8% on the model's own cells), by a fixed factor 28.1% (31.5%),
and the interpolated shares, 27.9% (30.8%), lie at the factor's end ([coast water](../coast_water/README.md)); a
design case at 301 K would settle it. Each kelvin warmer raises the freezing level over the summit about 1.2 km across
298–303 K and the ground's wet-bulb 99th percentile 0.84 K, takes 0.73–0.96 points off the dimming, 87–113 TW of
sunlight, adds 0.1 mm of rain a day on land with 0.9–1.8 points less desert by the interpolation's form, and costs the
seas 0.12 mg/L of oxygen and the margin 1.7–2.3 W/m² (targets; the freezing level read from the tower's rows, the
other rates as b83d542 corrected them).

## Each target's case

Here the climate figures are the GCM's output for the hot world and the lock and interpolations for the resort; the
dimming, its sunlight, the margin and the lakes' oxygen are calculations.

**The hot world** needs no dimming of its own, so it keeps the most rain, the least desert, the most light for growth
and the aerophytes, and the array's light for power. Its fresh seas hold 7.06 mg/L of oxygen at saturation, 10% more
than Earth's tropical seawater at 28 °C (lakes), and its comfort lives aloft, in air of 20–25 °C from about 8 to 12 km
up at 0.92–1.00 atm. The ground is hot and humid at every latitude: the air at 2 m has a median of 33.6 °C and a 99th
percentile of 39.3 °C, and the wet-bulb temperature passes 28 °C for 49% of the land's time within 15° of the equator,
and 31 °C, near where young, healthy adults can no longer keep up light activity (Vecellio et al. 2022), for 2.1%. Its
hottest days lie far above the register's pleasantly warm ones (25 September), and with the tower's top above freezing
it misses the register's wish for winter on the tallest sky tower (17 October). Its deep lakes turn over half as often
as the cooler targets' and keep 4.2–4.9 mg/L in their oligotrophic deepest tenth by Osborn's rule, 1.6 mg/L and a third
of the time anoxic by the energetic regime, with richer deep water low or anoxic. Its balance rests on unsaturated air,
1.3–1.8 W/m² under a saturated column's ceiling, and it still needs the shared 3.7–5.3% of dimming, without which the
added heat would warm it 4–7 K more.

**The resort**, at 301 K, keeps a winter on the tower and a warm, tropical ground: 3.5 km of the tower in sub-freezing
air, its top at −3.0 °C and the top's warmest 6.5-hour mean at −1.7 °C, about a kelvin inside the 302.25 K at which 2 km
of the tower stays below freezing in the mean; the land's wet-bulb at 2 m with a median of 22.0 °C and a 99th percentile
of 26.8 °C, 27.1 °C by day; the land's air at a median of 27.3 °C; and seas near 29.5 °C. It is the tropical resort Moon
the register names, with a fifth of the land kept as megaparks and gardens. Against 300 K it takes 0.7–1.0 points off
the dimming and adds 3.5% to the land's rain. It still takes 8.2–10.9% of the sunlight in all, 531–660 TW for its own
share, from a fleet whose tiles and photon authority already bind any dimmer. The open winter parks at 13.5–15 km lose
their winter above about 296 K; at 301 K winter holds the highest disk, 21.4 km up, and the crown, both sealed in the
present design, so open winter decks would move up to them. Its 352 m lakes turn over in 54–68 of 120 lunar cycles and
keep 5.7–6.3 mg/L in their oligotrophic deepest tenth, 3.5 mg/L by the energetic regime; its giants, 70% of the lake
area, go anoxic at depth when mesotrophic or richer. Mixing two physics, its figures are the least certain of the three
until a design case at the target settles them.

**The lock** is the design case, one physics and the GCM's own output, around which the plan's designs and experiences
have been built since 9 October. It is the mildest at the ground: the land's wet-bulb at 2 m has a median of 16.6 °C and
a 99th percentile of 21.8 °C and never passes 28 °C, and the air's 99th percentile is 27.8 °C, so its hottest days are
pleasantly warm. It keeps the tower's designed winter, 9.9 km of sub-freezing air with the open winter parks near
freezing, 0.2 °C at 13.8 km. Its waters hold the most oxygen, its oligotrophic and mesotrophic deep water staying oxic
under either mixing form, and its margin under a saturated column is the widest. It is also the costliest and the
driest: 12.6–16.6% of the sunlight in all and 1,052–1,340 TW for its own share, the dimmer growing from 3% of the
titania stack's light to 15–19%; 2.25 mm of rain a day on land with 38% of it desert; and 0.62–3.3% of the growth lost
through light, most of it from the light-sensitive life.

At the ground the rain falls from cloud bases 7.4, 6.4 and 5.1 km up in the hot world, the resort at 301 K and the lock,
and only drops that leave the cloud at 2.4, 2.0 and 1.6 mm or more survive the fall (rain, a calculation on independent
drops, its column at the repository's heat capacity). In the sound lab the Moon's heavy rain is 16–19 dB quieter than
Earth's and darker in every target, its drops striking at about a third of Earth's speed, and the hot world's is the
loudest (a rendering of an assumed transfer, not yet heard).

## Plant growth

Through light, each point of dimming costs 0.07–0.30% of the living Moon's growth (the light study and the targets'
`costs_per_percent`, a calculation on the canopy model and production by place): 0.31–1.64% for the resort's own
dimming, 0.62–3.3% for the lock's, and 0.26–1.56% in every target for the shared 3.7–5.3%. Through rain the drying costs
more. Rain-fed land grows 57–87% of the living Moon's growth by the water it transpires (the light study), and cooling
dries it: the land's rain falls from 3.34 mm a day in the hot world to 2.82 in the resort and 2.25 in the lock, and its
desert grows from 19% to 22.5–28% and 38%. The dimming's own cut at a held temperature takes a further 0.55% of the
rain a point in the GCM and 1.38% in Earth's models, bringing the resort's land to 2.60–2.75 mm and the lock's to
1.90–2.14.
The parts that follow the light grow 10–40% of the whole, and heavy dimming falls hardest on them: the aerophytes'
growth follows their light at 0.30–0.67 of each share taken, the lakes' and the river-fed coasts' at 0.38–1.0, the wet
margins' at 0.60–0.71, and the deep canopy loses 0.97% of its twilight photosynthesis for each 1% of light, against
0.39% under a high Sun. The light study's ledger is the design climate's.

## The lakes

The deep crater lakes, 352 m deep at the median, turn over far more often than first found: by the wind's energy below
the night's mixed layer, the share Earth's lakes pass to their stratified water (0.3%, 0.2–0.4%), the whole column turns
over in 62–74 of the 120 lunar cycles of ten years in the lock, 54–68 in the resort at 301 K and 29–30 in the hot world,
its deep water a mean 0.1–0.3 years from the air in the cooler targets and 0.4–1.4 in the hot world (lakes). That rests
on Osborn's rule for turning dissipation into mixing; Shih and Bouffard's energetic regime, the other end of the form,
turns the 0–15° lake over 11 times in the lock, 8 in the resort at 300 K (5 in its nutrient-poor water at 301 K) and 12
in the hot world, so that under this form the resort's lakes turn over least. Their deepest tenth keeps its oxygen in
nutrient-poor water: 6.8–7.4, 5.7–6.3 and 4.2–4.9 mg/L in the three by Osborn's rule, 4.8, 3.5 and 1.6 by the energetic
regime, the hot world's anoxic a third of the time. Mesotrophic deep water stays oxic in the cooler targets by Osborn's
rule (5.0–6.1 mg/L in the lock, 3.8–4.8 in the resort) and runs low or anoxic in the resort and the hot world by the
energetic regime; eutrophic deep water keeps 1.9–4.2 mg/L in the cooler targets and runs anoxic half the time in the hot
world, and new reservoirs run low in every target. Productive shallows, mixed through the day by the same share, stay
above 2 mg/L on almost every night: eutrophic ones fall below it on at most 1% of nights, new reservoirs on up to 3% of
the resort's and 4–19% of the hot world's, and the eutrophic ones on up to 12%, 38% and 70% of the three targets' nights
at the bracket's fixed dissipation.

At 301 K the lakes by kind put a tenth of the lake area in lakes the nights ventilate to their floors (403,700 km², from
226,500 to 629,700 by the bracket), a fifth in deep lakes to 1,000 m, whose floors stay oxic in nutrient-poor and
mesotrophic water, and seven tenths in 1,367 giants deeper than 1,000 m, which turn over once in years, keep 2.5–4.5
mg/L at depth to 1,925 m in oligotrophic water and go anoxic at depth when mesotrophic or richer; the 3,500 m giant goes
anoxic within thirteen years even oligotrophic. The held fresh seas renew within years by Osborn's rule, their
area-weighted deepest tenth at 4.1 mg/L nutrient-poor, 2.1 mesotrophic and 1.5 eutrophic, and the great basins of
Procellarum and Ingenii renew far less by the energetic regime (the seas' oxygen). Chemoclines do not hold under the
wind, and degassing does not pay.

The methane follows the water's richness more than its warmth. At 301 K nutrient-poor lakes and seas send 0.57 Tg a year
and, with the ledger's other sources cut by biology's 90%, hold the air at 1.01–1.75 ppm, under the held 1.87 ppm with
thin room, at the methane share's centre under Osborn's rule; at the share's 0.8 they need 0.5–2.6 Tg a year of removal,
and under Shih and Bouffard's form up to 0.77; mesotrophic waters send 7.8 Tg and need 4.7–6.9 Tg a year of removal at
the held pair, and the giants carry most of the richer waters' methane (the waters' methane). With the lakes alone
scaled from their 352 m lakes, nutrient-poor lakes send 0.37, 0.34 and 0.29 Tg a year in the lock, the resort at 300 K
and the hot world ([lake anoxia](../lake_anoxia/README.md)). The deep lakes stay nutrient-poor while people's dissolved
phosphorus stays within their room: with 8 billion people at the middle settling, 57–71 g a person a year while the
giants' floors stay oxic and 5–10 g if their depths go anoxic
([the resort at 301 K across the domains](../resort_301k/README.md)), and their remedies serve every target: oxygen
piped to the giants' depths or driven mixing, deep communities built for little oxygen, and shallows kept clear and
nutrient-poor or aerated.

## The levers that ease the shield

The levers give each saving in points, a point being 1% of the hot world's absorbed sunlight, 225.4 W/m², the unit of
the targets' dimming, on the corrected targets with the resort at 301 K and the glow's efficacy at 0.645–0.71. Ranked:

| Rank | Lever | Points saved | What it takes | Evidence |
|---|---|---|---|---|
| 1 | The glow's set-aside answers: the annulus absorber or the ultraviolet mirror | 0.62–0.83 | A film change and a check above the limb, or 60% more fleet mass | Screen on the fleet's products |
| 2 | Food made in orbit with the beams' power | 0.13–0.98 | Food plants in the array's industry; the food's carbon, 0.47–0.71 Gt a year, needs a return path | Arithmetic on the synthesis's ledgers and the food study |
| 3 | Methane and N2O at the low pair | 0.79–0.81 | With nutrient-poor waters, contactors through 0.4–0.8% of the air a year at the methane share's centre, up to 1.4% across its brackets; N2O decomposition 0.96–1.47 TW | Column forcing; design-guess costs; the waters traced at 301 K |
| 4 | White roofs and silver foliage | 0.07–0.39 | Little beyond adopted designs; cooler cities | Kernel arithmetic by latitude |
| 5 | Night-side cirrus thinning | 0 to 1.0, sign open | 12–240 kt of feldspar a year from high platforms; adds rain | Earth's models |
| 6 | The water top-up's fall | 0 to 1.6–2.3, first millennia only | Momentum binds any braking above the air | Energy and momentum accounting |
| 7 | 1.0 atm | 0.91–0.97 clear sky in the hot world, 0.46–0.71 with cloud masking; −0.11 in the lock | Reopens 1.2 atm; 17% less lift | Clear-sky columns, Earth's cloud masking |
| 8 | Brighter deserts | 0.17–1.94 | Covers wild land; takes rain four times as fast as dimming | Kernel arithmetic by latitude, Earth's models |
| 9 | CO2 below 400 ppm | 0.33–0.90 | 33–44% of the growth a point at 280 ppm | Column and leaf model |

The first four together save 1.6–3.0 points, against dimming in all of 3.7–5.3% for the hot world, 8.2–10.9% for the
resort and 12.6–16.6% for the lock. Taken with the case each belongs to, the glow's answers at the glow's efficacy and
food made in orbit at the beams' end, they save 43–57%, 17–31% and 11–21% of them and leave the resort 6.1–8.3%. They
can answer much of the plan's added heat, and bringing the Moon down to the resort or the lock stays mostly dimming. The
gas pairs carry their own overlap, the other levers act on separate terms, and no lever has yet run in a climate model.

## What stays open

- **The two physics in the interpolation.** The design case flew PlaSim's earlier longwave with liquid optics for all
  cloud and the test the refitted longwave with ice optics, so the resort mixes two models; a cooler Moon under the
  refit would keep some of its high cloud, in an amount unknown. Along the test's own warming the tower's top runs
  0.3–0.6 K colder at 300 K than the interpolation gives (targets, the GCM's output).
- **The GCM's high-cloud check.** The refit's high cloud carries about 5.6 W/m² of the hot world's net cloud warming
  (levers, arithmetic on the GCM's output). Its range from its water and its crystals' size, offline in minutes, may
  move every target's dimming more than seeding could, and the longwave refit's check at the settled state goes with it
  ([plan](../../plan.md), item 5).
- **Convection over water at lunar gravity.** PlaSim carries Miller et al.'s free convection over water at Earth's
  gravity, which at lunar gravity falls to 0.55 of that. In the lakes' slab check, a calculation, the Earth constant
  leaves the GCM's water 0.5–1.2 K cooler, so the seas' figures, their oxygen with them, may sit up to about a kelvin
  cool.
- **Runs at the chosen target.** A design case at the target, a 20-year GCM run of about 1.5 hours on four ranks,
  settles its climate and the tower's air with one physics; the glow pair from it, two more runs, fixes the glow's
  efficacy, now 0.645–0.71; the lake studies, which keep the targets they ran on until then, geography's lakes, the
  atlas globe and the studies that read the design climate then regenerate once. A range carried by its ends takes a
  design case at each end.
- **The author's choices.** The target, or the range or distribution to carry (the resort set at 301 K on 18 October);
  and the levers' questions: whether the glow's set-aside answers return if a cooler target is chosen, whether the
  beams' food may be made in orbit, with a return path for its carbon, 0.47–0.71 Gt a year, and whether the high cloud's
  range runs before the choice. The choice reopens the register's planned climate near 295 K (9 October), its hottest
  days pleasantly warm (25 September) and the 5% dimmer's 295 K.
