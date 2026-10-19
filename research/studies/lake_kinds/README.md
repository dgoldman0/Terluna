# The lakes by kind at the resort's 301 K, and what each kind needs

18 October 2026. The author set the resort at 301 K on 18 October, the climate-targets study's recommendation within
300–302 K ([register](../../decisions.md#the-target-climate-17-october-2026)), and, after
[the anoxic depths](../lake_anoxia/README.md) found the deep lakes' water going anoxic, asked for the lakes to be sorted
by what they need, out of a concern that lakes scattered over the whole Moon would each need tending. The same day the
author chose to drive the lake model's mixing below the mixed layer by the wind's own energy. This study runs the
[lake model](../../../biosphere/ecology/lake.py), as corrected after the independent review of 18 October and with that
wind-driven deep mixing, at 301 K on lakes whose floors lie at a ladder of depths from 5 m to 3,500 m, in each trophic
state; sorts every one of geography's 27,865 rain-fed lakes by kind from what the ladder shows; and sizes what each kind
needs. [census.py](census.py) writes the lake set's geometry, outlets, soda and groundwater, which no lake model sets,
to [results/lake_census.json](results/lake_census.json) (schema `terluna.research.lake-kinds-census/1`);
[run.py](run.py) writes [results/lake_kinds.json](results/lake_kinds.json) (schema `terluna.research.lake-kinds/1`) from
the runs of [ladder.py](ladder.py); [sources.json](sources.json) records each published value and how far it was read.
The ladder turns the wind's dissipation into mixing by Osborn's rule; the waters' oxygen study's
[mixing_forms.py](../lake_oxygen/mixing_forms.py) runs its 352 m lake and its 1,925 m giant at 301 K by Shih and
Bouffard's energetic regime, the other end of that range, which this study's figures stand beside.

## What it finds

**The resort at 301 K.** Built as the waters' oxygen study builds its resort, with climate_targets' share of the way at
301 K, 0.520, the seas' median is 29.6 °C and the land's air at 2 m 27.2 °C, against 28.7 and 26.1 °C at 300 K. The
352 m lakes there, run as lake_oxygen runs them, turn over in 54–68 of 120 lunar cycles against 55–68 at 300 K, and
their deepest tenth keeps up to 0.3 mg/L less oxygen in every trophic state: 5.7–6.3 mg/L oligotrophic, 3.8–4.8
mesotrophic, 1.8–3.0 eutrophic and 0.2–1.2 in new reservoirs.

**With the wind's mixing, the nights reach the floors of lakes down to about 170–240 m on most nights.** Each trophic
state run with its own physics, at the wind's central share of its energy below the mixed layer and Osborn's form for
the mixing, over the three bands:

| Floor depth | Nights mixing to the floor | Lunar cycles that overturn | Floor oxic (above 2 mg/L) | Deepest tenth, mg/L |
|---|---|---|---|---|
| 5–30 m, every state | 100% | 100% | 95–100% | 3.9–7.5 |
| 50 m | 96–100% | 100% | 94–100% | 3.2–7.0 |
| 100 m, oligotrophic; mesotrophic | 70–93%; 76–94% | 87–98%; 91–100% | 100% | 6.7–6.8; 5.6–6.0 |
| 100 m, eutrophic; new reservoir | 81–97% | 91–100% | 100%; 55–84% | 3.6–4.2; 2.0–2.9 |
| 352 m, oligotrophic; mesotrophic | 22–31%; 29–41% | 34–45%; 45–57% | 100% | 5.7–6.3; 3.8–4.8 |
| 352 m, eutrophic; new reservoir | 30–50% | 45–66% | 3–73%; 1–3% | 1.0–2.3; 0.0–0.2 |

The wind's share, spread over the little water a shallower lake holds below its mixed layer, keeps that water close to
the mixed layer's density, so most nights mix lakes to 100 m through to the floor, while the ordinary nights that stop
short of a 352 m lake's floor mix it to a median 20–23 m oligotrophic, 12–13 m mesotrophic and 6–7 m in richer, darker
water. The ventilated depth, the floor the nights reach on half of them, read along the rungs, is 168, 239 and 192 m in
the three bands for oligotrophic water, 201–284 m mesotrophic and 215 m to beyond 352 m eutrophic; at the share's 0.2%
and 0.4% the oligotrophic depth runs 159–218 and 196–278 m. A floor reached on three nights in four lies within
87–145 m, one reached on a night in four within 323 m to beyond 352 m.

**The kinds.** Every one of the 27,865 bodies by its deepest point against its band's ventilated depth for oligotrophic
water, and the edges at 50 and 1,000 m:

| Kind | Bodies | Area | Share of the area | Mean depth | Floor within the ventilated depth |
|---|---|---|---|---|---|
| Self-ventilated: deepest point within the ventilated depth | 22,384 | 403,700 km² | 10% | 39 m | all |
| The middle lakes: beyond it to 50 m | none, the ventilated depth lying below 50 m in every band | | | | |
| Deep: beyond the ventilated depth to 1,000 m | 4,114 | 776,100 km² | 20% | 202 m | 59% |
| Giants: deeper than 1,000 m | 1,367 | 2,743,000 km² | 70% | 1,070 m | 16% |

Four bodies in five are self-ventilated, holding a tenth of the area; seven tenths of the area lies in the giants. With
the ventilated depth taken where three nights in four reach the floor, the self-ventilated lakes hold 19,828 bodies and
227,500 km², and where one night in four does 24,190 bodies and 610,400 km². Sorted under their own ventilated depths,
mesotrophic lakes have 458,900 km² self-ventilated and eutrophic ones 505,800 km², and the share's bracket moves the
oligotrophic sorting between 361,700 and 448,600 km². The deep lakes, a median 368 m at their deepest point, turn over
in a third to two thirds of the lunar cycles: their floors stay oxic in oligotrophic and mesotrophic water, and in
eutrophic water and new reservoirs the deeper floors go anoxic between overturns.

The ventilated depth lies between the ladder's 100 m and 352 m rungs in every band and is interpolated in the log of the
depth across that gap, where many lakes have their deepest points; the census at those two rungs brackets the split,
the ladder's oligotrophic lakes at 100 m mixed to their floors on 70–93% of nights and those at 352 m on 22–31%:

| Band (its ventilated depth) | Deepest point within 100 m | Within the ventilated depth, self-ventilated | Within 352 m |
|---|---|---|---|
| 0–15° (168 m) | 9,618 bodies, 155,900 km² | 10,831 bodies, 262,700 km² | 12,162 bodies, 432,900 km² |
| 15–30° (239 m) | 4,649, 52,100 km² | 5,708, 107,700 km² | 6,083, 139,000 km² |
| 30–45° (192 m) | 1,605, 15,000 km² | 1,910, 27,830 km² | 2,107, 49,930 km² |
| 45–90°, at the 30–45° depth | 3,840, 3,528 km² | 3,935, 5,483 km² | 3,979, 7,929 km² |
| All | 19,712, 226,500 km² | 22,384, 403,700 km² | 24,331, 629,700 km² |

A rung near 200 m would narrow that bracket.

**The giants turn over once in years, and their richer depths go anoxic.** Each giant's fifty years in its own regime,
at 0–15° unless named:

| Giant | Lunar cycles that overturn, of 604 | Longest stretch without one | Deepest tenth, mg/L |
|---|---|---|---|
| 1,000 m, oligotrophic; mesotrophic | 80; 90 | 2.8; 2.4 years | 4.5; 1.0 |
| 1,925 m, oligotrophic (0.2% and 0.4%) | 40 (20 and 50) | 6.6 years (8.3 and 3.6) | 2.6–2.8 (0.9–1.3 and 3.5–3.6) |
| 1,925 m, mesotrophic (0.2% and 0.4%); eutrophic | 50 (35 and 65); 105 | 3.6 years (6.6 and 3.6); 1.4 | anoxic |
| 1,925 m at 15–30° and 30–45°, oligotrophic | 55 and 44 | 7.9 and 4.6 years | 3.8 and 2.5–2.6 |
| 1,925 m at 15–30° and 30–45°, mesotrophic | 65 and 70 | 3.9 and 3.3 years | 0.01–0.04 and anoxic |
| 3,500 m, oligotrophic; mesotrophic | 10; 30 | 9.6; 7.2 years | 0.5–0.7, anoxic within 13 years; anoxic |

An overturn mixes a giant through, but its oxygenated top is a small share of its volume, and the floor's and deep
water's demand empties the depths between. Oligotrophic giants to 1,925 m keep 2.5–4.5 mg/L at depth; the 3,500 m
giant goes anoxic within thirteen years; mesotrophic and richer giants are anoxic at depth. The deepest oligotrophic
giant had not settled in its 300 years of coarse passes and 100 more at the analysis step, nor had the 1,000 m
mesotrophic one and the 1,925 m oligotrophic one at 0.4% at the analysis step.

**A chemocline at 20 m does not hold under the wind's mixing.** The 30 and 50 m lakes at 0–15° with their water below
20 m made denser by 0.1 kg/m³, and the 30 m lake's mesotrophic and eutrophic water by 1 kg/m³, keep a mean 1–4% of that
contrast through their first year and 0–0.2% through their fiftieth: the nights mix them through to the floor from the
first months, and only new reservoirs' water below the step falls below 2 mg/L, within 0.25 and 1.2 years. Under
Osborn's form the wind's mixing leaves a lake of tens of metres no isolated depth held by such a step.

**Closed lakes and soda.** The drainage's 2,287 closed parts cover 32,150 km²: 2,105 of them self-ventilated (17,750
km², a mean 31 m deep), 164 deep (9,956 km²) and 18 within giants (4,444 km²). Along the sea chemistry's left-alone
path at each part's own turnover, the self-ventilated ones reach 1 g/kg in a median 10,450 years (6,985–12,420 by the
weathering's bracket) and Earth seawater's salinity in 312,000; the deep ones 95,600 and 2.9 million years. The 182
closed parts beyond the ventilated depth, 14,400 km², could hold a chemocline once their salinity outweighs a wet
spell's freshening of the mixed layer: the deep ones after 154 years at the cold spells' smaller margin and a
twentieth's freshening, 1,270 years at a hundredth's, and 6,970–306,000 years at the larger margin; the meromictic runs
above say a chemocline that weak would not hold against the wind's mixing.

**Springs.** At the 352 m oligotrophic lakes' mean bottom temperature, 27.2 °C, the groundwater's water is 0.051–0.139
kg/m³ denser than the lakes', its salts against the lakes' 0.040 g/kg, a density 0.18–0.50 K of warmth would cancel
(the product's spring_water row at 27.2 °C; the census's rows at 22, 25 and 28 °C give 0.22–0.61, 0.20–0.54 and
0.18–0.49 K). Entering the deep layer across the floor, it could
build the cold spells' smaller margin, 0.003 kg/m³, within the longest stretch without an overturn in 123 deep bodies,
4,075 km², and the larger, 0.03 kg/m³, in 15 bodies, 537 km²; the median resolved deep lake takes 407 years to the
smaller margin. No giant qualifies, as their overturns come within ten years.

**What each kind needs.** The needs at 301 K, with the oxygen study's remedies run on this study's 352 m lakes:

- **Self-ventilated lakes** (22,384 bodies): watched from orbit and the shore. Their floors are oxic through the
  record in every trophic state to 50 m and in all but new reservoirs to 100 m; productive ones still run low in the
  night, eutrophic 5 m lakes to a median 2.3–4.0 mg/L, below 2 mg/L on up to 36% of nights at 15–30°. They stay
  oligotrophic under 0.10 g of bioavailable phosphorus per m² a year with Earth's settling, 0.052 at the square root of
  g/6 and 0.032 at g/6, against a dissolved load of 0.018–0.048.
- **Deep lakes** (4,114 bodies): one profiling mooring at each deepest point and one more for each further 2,000 km²,
  4,138 moorings. Oligotrophic and mesotrophic ones keep oxic floors through the overturns; eutrophic ones and new
  reservoirs need their nutrient limits, 0.21, 0.095 and 0.049 g P m⁻² a year by the settling's bracket, or
  oxygenation at 10–40 kW per km² eutrophic and 14–58 in new reservoirs. Degassing does not pay at 301 K: the deep
  water holds too little methane, 0.25 kW of it per km² against 272 kW of pumping in eutrophic water.
- **Giants** (1,367 bodies): 2,061 moorings. Their depths go anoxic in mesotrophic and richer water and in the deepest
  oligotrophic ones, so deep-water life there needs oligotrophic water and, below that, oxygenation by design: 3.5–14 kW
  per km² oligotrophic. That rate taken for the deep kind and the giants alike, 1.2–5.0 GW oxygenates a tenth of their
  3.5 million km² together and 0.12–0.50 GW a hundredth. Their oligotrophic limit is 0.45, 0.20 and 0.090 g P m⁻² a
  year; where their floors go anoxic the retention can at most be cancelled, the limit then falling to the throughflow's
  0.018 g.
- **Built lakes.** A reservoir or park lake to 50 m keeps an oxic floor through the record with no running energy when
  oligotrophic (night minima 6.3–6.9 mg/L), and in its eutrophic first decade as a new reservoir keeps 94–100% of its
  floor oxic, its nights falling to 2.1–3.8 mg/L; a deep new reservoir consumes 1.31 t of oxygen per km² a day through
  that decade, which oxygenation delivers as 1.3–3.9 t, up to three times it as the devices' stirring raises the
  demand, for 14–58 kW per km².
- **Overturns in nutrient-poor lakes.** No overturn of an oligotrophic lake from 16 m to 352 m takes its surface below
  2 mg/L, in any band.

**The mixing's form.** All of this runs on Osborn's rule, K = 0.2 ε/N², which the lake model applies at every buoyancy
Reynolds number. In the near-neutral deep water the wind's share keeps, Re_b reaches 10³–10⁵, where Shih et al. (2005)
find the efficiency falling, K = 2ν Re_b^½, as Bouffard & Boegman (2013) carry it to lakes: Osborn's divided by 0.1
Re_b^½, 3–30 times less mixing there. The waters' oxygen study's [mixing_forms.py](../lake_oxygen/mixing_forms.py) runs
the 0–15° lake at 301 K under that form: the 352 m lake turns over in 5, 8, 14 and 13 of 120 lunar cycles oligotrophic,
mesotrophic, eutrophic and as a new reservoir, against 44–55 under Osborn's, and its deepest tenth keeps 3.5 mg/L
oligotrophic against 6.0, 1.0 mesotrophic (anoxic 48% of the time) against 4.2, and nothing in richer water; the 1,925 m
giant turns over twice in fifty years and goes anoxic at depth even oligotrophic, where it keeps 2.7 mg/L under
Osborn's, and its mesotrophic, eutrophic and new-reservoir runs turn over 15, 35 and 35 times, anoxic at depth as under
Osborn's. The ladder's shallower rungs were not run under it. The kinds' figures here therefore lean towards more
mixing: under the energetic regime the ventilated depth would be shallower, the deep lakes would overturn about once a
year, and the deep kind's and the giants' floors would go anoxic in mesotrophic and richer water and the giants' in
oligotrophic water too.

## How it is set up

**The resort at 301 K.** [setting.py](setting.py) builds the resort as the waters' oxygen study builds it: the lock's
weather, the design case's 3-day means, moved by the share of the way to the hot world's test that climate_targets
gives the target on the ground's temperature, (301 − 295.02)/(306.52 − 295.02) = 0.520 (0.433 at 300 K), with the air
carried from each GCM cell down to its lake along the lapse rate of the sky ladder interpolated the same way,
0.843 K/km between 0.9 and 5 km (0.828 at 300 K). Built at 300 K, it gives climate_targets' own share and ladder and
lake_oxygen's lapse rate exactly, and lake_oxygen's run_case() on it gives that study's resort runs, its three 352 m
lakes and the 0–15° pan, identical in every figure; this study's column runner gives the same 352 m lake.
`Setting301` is the one 301 K construction the other lake studies import. At 301 K the seas' median is 29.6 °C and the
land's air at 2 m 27.2 °C, against 28.7 and 26.1 °C at 300 K. Like lake_oxygen's resort, it takes the design case's
wind factor and mixes two physics.

**The ladder.** [ladder.py](ladder.py) runs the lake model
([biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py), as corrected after the review of 18 October and its
re-check, with the wind's deep mixing the author chose that day) through [columns.py](columns.py) on each band lake's
own hypsography, Mandel'shtam's, Gagarin's and Rabbi Levi's, scaled to floors at 5, 8, 12, 16, 20, 30, 50, 100 and 352
m. Below the mixed layer the dissipation is the share of the wind's energy flux at 10 m that reaches the stratified
water, spread over its mass, with Osborn's rule for the diffusivity: 0.3% of the flux, the central value of Earth's
lakes' 0.3 ± 0.1% (Wüest, Piepke & Van Senden 2000), with 0.2% and 0.4% as the bracket on the 352 m lakes' and the 1,925
m giant's oligotrophic and mesotrophic water. Each trophic state runs with its own physics, its water's clarity for the
light and its own respiration's products for the density, where lake_oxygen ran every deep lake with mesotrophic water;
each run still carries all 24 oxygen columns. Lakes to 50 m step every hour, as lake_oxygen's shallows, the deeper ones
every two hours, as its deep lakes. Each spins up as lake_oxygen's deep lakes do, from 1 K above the air's mean and
saturated, passing over the ten-year record in coarse steps until its deep water's temperature, oxygen, age and solute
all settle (to 0.005 K, 0.02 g/m³, 2% and 0.02 g/m³ a pass), or for twenty passes; the lakes of 100 m and more then pass
over it at their own step until the same test passes there, or for ten passes, and every lake then runs the record. The
settle test takes the age and the solute because a column whose deep water has stopped exchanging keeps ageing and
gathering respiration's products, which a test on its warmth and oxygen alone would pass; the second stage settles the
coarse steps' slightly different deep water at the step the record is analysed at, which the shallower rungs,
overturning within months, do within the analysed decade. The giants, floors at 1,000, 1,925 and 3,500 m, start from the
352 m column of the same band, trophic state and wind share carried down and oxygenated, as lake_oxygen starts its 1,925
m lake, pass over the record in the same coarse steps until they settle or for 300 years (thirty passes), and at the
analysis step until they settle there or for 100 years more, so that they reach their own regime, and then run fifty
years, five passes of the record, saved and checkpointed after each pass: at 0–15° the oligotrophic and mesotrophic
states at each depth and the eutrophic one at 1,925 m, at the other bands the oligotrophic and mesotrophic states at
1,925 m; a new reservoir is no giant. [meromictic.py](meromictic.py) holds the 30 and 50 m lakes at 0–15° meromictic for
fifty years: from the same oxygenated start, their water below 20 m carries an added step of the model's one density
tracer, respiration's solutes, worth 0.1 kg/m³ (and 1 kg/m³ for the 30 m lake's mesotrophic and eutrophic water). Every
column is saved every 23.8 hours, lake_anoxia's spacing, with its temperature, oxygen columns, debt, age, solutes and
labile pool, and each column's surface oxygen and the dissipation at every step.

**Reading a column.** [analysis.py](analysis.py) takes each night, sunset to sunrise at the lake, and its deepest
mixing; a layer is ventilated in a night when the mixed layer reaches its centre, and each layer's floor is the area the
basin loses across it, so the floor's shares are of the lake's bed in plan. A complete overturn is a lunar cycle in
which the mixed layer reaches the bottom layer, and consecutive ones make an episode. The ventilated depth is the floor
depth to which the nights mix a lake to its floor on half of them, a quarter and three quarters beside it, read along
the ladder's rungs and interpolated in the log of the depth between them: the wind's mixing keeps a shallower lake's
water below the mixed layer near the mixed layer's density, so the nights reach the floors of lakes far deeper than
they mix a 352 m lake on its ordinary nights, the nights that stop short of its floor, which are given beside it. An
isolated hypolimnion's years to its bubble point are lake_anoxia's (`isolated_years`, with its rule for the
hypolimnion's thickness) at each lake's own deep warmth.

**The census.** [census.py](census.py) takes geography's rain-fed lakes as lake_oxygen's lakes.py labels them on LOLA's
16 pixel/degree grid, the drainage summary's 27,865 bodies, each cell's depth its lake's level less the ground, and
writes [results/lake_census.json](results/lake_census.json) (schema `terluna.research.lake-kinds-census/1`). A body's
kind is set by its deepest point against its band's ventilated depth and the edges at 50 and 1,000 m; bodies poleward of
45° take the 30–45° band's lake. The census sorts every body under each trophic state's ventilated depth at the wind's
central share and, for the oligotrophic and mesotrophic states, at its bracket; its central sorting takes oligotrophic
water at the central share, the state the land's dissolved phosphorus holds the lakes in (lake_anoxia). A body with at
least half its cells in depressions that do not spill is closed; its closed cells are a lake of their own, and their
depth and their open-water evaporation and rain (geography's, from the climate run it routed) set their soda path.

**Soda.** The closed lakes run the sea chemistry's left-alone path of the fresh seas (natural sinks, kept open), each
at its own turnover: the sea chemistry's scaling takes the lake's depth over its evaporation against the waters'
1,908 years; the balance variant takes the inflow that the drainage balances against the evaporation, the evaporation
less the rain on the lake, carrying the land water's 1.41 meq/kg of alkalinity against the rivers' 0.75 (springs). A
fresher surface holds a chemocline once its density deficit passes the margin the cold spells need: a closed lake of
salinity S whose mixed layer a wet spell freshens by a share f holds one when 0.8 S f kg/m³ passes 0.003–0.03 kg/m³.

**Springs.** Groundwater comes from stage 2's mature water table at the central crust, re-solved by geography's own
solver from its saved heads; it reproduces the product's 114.6 km³ a year into the lakes and its 886 gaining lakes, and
its gross inflow is the springs study's 145.8 km³ a year. Bodies on its 15 km grid take their lake's inflow per area;
smaller ones take their band's median inflow per metre of shore along their own shore. The springs' water (springs
study, mature and central, by band) is compared with the lakes' (sea chemistry's held lakes, 0.040 g/kg) at the
model's 0.8 kg/m³ per g/kg and Lake Kivu's coefficient for its CO2, at the 352 m oligotrophic lakes' mean bottom
temperature, 27.2 °C. A spring-fed candidate is a lake beyond the
ventilated depth whose inflow, all of it entering the deep layer across the floor at the deep water's warmth, builds
the margin within the longest stretch without an overturn in the oligotrophic run of its band at the next depth at or
below its deepest point, the ladder's lakes and then the giants.

**The needs.** The remedies are lake_anoxia's, run on this study's three 352 m lakes at 301 K: oxygenation, degassing
pipes and circulation by Lorenzen & Fast's air per km², with Vollenweider's permissible load at each kind's mean depth
and the lakes' 1.78 m a year of throughflow under lake_anoxia's settling brackets, and with the retention cancelled, the
most anoxic floors can do at steady state; the oxygenation's bounds take the demand fixed by supply, with the warm
kinetics' transient beside them. The moorings are a proposed design at about the density of Lake Erie's oxygen
network.

## Assumptions and limits

- **The lake set is a first estimate.** Geography's drainage routed the runoff of an earlier climate run,
  A28_dim5_moon, at T21 over LOLA's 1.9 km cells; it regenerates once the design case for the chosen target settles,
  and the 301 K climate's stronger evaporation would shrink the closed lakes and move some of the 27,865 bodies across
  the kinds. Lakes smaller than a few cells, and depths within a cell, are below its reach. The drainage counts 2,495
  bodies with any closed cell; 208 of them are open lakes holding nested closed pits, 21,400 km² of pits in
  471,000 km² of open lake, which the census takes as open, a point for geography's drainage to settle.
- **One column for a lake.** Each lake is one column of the model through the GCM's 3-day means at T21 made hourly, at
  the band lake's cell. The Lake numbers, under 1 on many days whichever thermocline they are read at, say the wind
  stirs the lakes by basin-scale motions a column cannot hold; the column stands for them by the share of the wind's
  energy that Earth's lakes pass to their stratified water, spread evenly over the water below the mixed layer, where
  most of it is spent at the floor (Wüest, Piepke & Van Senden 2000), with Osborn's rule and the model's mixing
  efficiency of 0.2 against the 0.15 measured with it; the share's bracket stands on the 352 m lakes and a giant, and
  the drag the share was measured against is open (lake_oxygen). The form that turns the dissipation into mixing is
  Osborn's rule at every buoyancy Reynolds number, which leans towards more mixing: by Shih and Bouffard's energetic
  regime the 352 m lake and the 1,925 m giant at 301 K turn over a ninth to a quarter as often, the giant a twentieth to
  a third as often (mixing_forms), and the shallower rungs and the census were not run by it. The wind's sheltering by
  the crater rim, a small lake's short fetch, inflows and outflows, seiches and the sediment's heat are outside it.
  Shapes are the three band lakes' scaled: a giant with a broad, shallow-sided floor holds less of its water near the
  floor than Mandel'shtam's bowl, and 72% of the lake area lies in lakes over 1,000 km², whose median deepest points are
  1.7 km (to 10,000 km²) and 3.6 km (larger), so the deep and giant kinds rest on the 352 m lakes' nights and the
  giants' scaled columns.
- **The kinds by their deepest point.** A lake with a deep pit and wide shallows is sorted by the pit; its ventilated
  floor share, which the census carries for every body, says how much of its bed the nights reach.
- **Ten years of weather.** The deep water settles to the coldest spell in the GCM's record; longer records hold colder
  spells, and the design case's years 20–21 hold more overturns than years 10–19. The giants sit near the edge of
  overturning, where a few hundredths of a kelvin decide whether the coldest night reaches the floor.
- **The meromictic lakes.** The chemocline is a step in the model's one density tracer, respiration's solutes, which
  leave the mixed layer with the gas, as a lake whose surface a fresh inflow keeps fresh; its salt and its own heat
  balance are not modelled.
- **Soda by scaling.** The closed lakes run the sea chemistry's left-alone path at their own turnover, with its own
  limits: no floor of their own, no groundwater, and the rain on the lake counted only in the balance variant.
- **Springs.** The groundwater's temperature as it enters is not modelled; its density excess holds a deep layer only
  if it arrives within a few tenths of a kelvin of the deep water's warmth, and seepage on Earth enters mostly near the
  shore (McBride & Pfannkuch 1975), where the census spreads it over the floor.
- **Biology and remedies.** Demand comes from Earth's lakes for engineered communities not yet designed, fixed by
  their supply, the warm kinetics a transient beside it; the labile pool's e-folding is a design bracket; the remedies'
  energy and rates are Earth's, at the Moon's gravity, and the circulation's air is Earth's guideline with no plume
  computed; the moorings and the reservoirs' rule are proposed
  designs; one throughflow, the lake set's 1.78 m a year, serves every kind's phosphorus rule.
- **The targets as run.** The 301 K resort is built on the climate targets the lakes ran on, domain/ecology's at e9c12eb
  (schema /2), which lake_oxygen keeps byte for byte as climate_targets_as_run.json until the lakes rerun on the
  resort's 301 K design case. The merged targets' lapse rates over the GCM's own land mask, 1.2–1.8% lower, would cool
  these lakes' air at 301 K by 0.024 K at most and 0.015 K typically, the 301 K share unchanged, and their seas' median
  at 301 K, 0.12 K lower, would raise the saturation by 0.016 mg/L.

## The checks

[test_run.py](test_run.py) and [test_census.py](test_census.py) check the products' binding to their code and inputs
(the census, the ladder's index with every saved file's sha256, the lake model's); that climate_targets' own rule gives
the 301 K share and that the construction at 300 K reproduces the waters' oxygen study's resort exactly; that every run
keeps its heat and oxygen to round-off and no lake freezes; that the ladder holds every case, at the wind's share and at
its bracket where it was run, with the cap holding no interface back and the dissipation's median in proportion to the
share, and the giants passed over the record to their regime; that deeper floors are mixed through on fewer nights; that
the ventilated depth is where the nights stop reaching the floor along the ladder, three nights in four reaching no
deeper floor than one in four, and the 352 m lakes' ordinary nights deeper in clearer water, while the share moves the
overturns and leaves the ordinary night's depth; that the census is the drainage summary's lake set, its closed bodies
the mostly closed of the summary's, every body counted once in every table, and the groundwater stage 2's product
reproduced; the closed lakes' soda times against the sea chemistry's scaling; the springs denser than the lakes; that
the circulation, the isolated hypolimnion and the permissible load are lake_anoxia's at the same depths, with the
retention cancelled as anoxic floors cancel it; that the moorings add up; and that each overturn event is counted once.

## Running it

```sh
python -m research.studies.lake_kinds.census      # a few seconds: LOLA, the drainage, groundwater stage 2 re-solved
python -m research.studies.lake_kinds.ladder      # about 2.5 hours on two cores; a stopped run resumes where it stopped
python -m research.studies.lake_kinds.ladder --cache-only    # the index and the 300 K proof from the cache, seconds
python -m research.studies.lake_kinds.run         # a few minutes: the saved columns, the census and the needs
python -m research.studies.lake_kinds.run --cache-only       # the same, stopping unless ladder.json is current
python -m pytest research/studies/lake_kinds
```

ladder.py writes its cases to `research/runs/lake_kinds/cases` (ignored, on the research drive): for each, a summary in
lake_oxygen's form and its passes' columns, saved about daily in lake_anoxia's arrays (temperature, the 24 oxygen
columns and their debt, the water's age, the solutes and the labile pool) with the step-by-step series.
`research/runs/lake_kinds/ladder.json` indexes them with the sha256 of every file and the lake model's, and run.py
checks each against it. A case's key holds the code that runs it, the case and what the runs read from their inputs
(lake_oxygen's `Setting.values_fingerprint()`), so a product written again with the same values leaves the cache
standing; with `--cache-only` ladder.py computes no run, and run.py stops unless the index holds every case under the
key the current code and setting give it. [setting.py](setting.py)'s `Setting301` is the one 301 K construction the
other lake studies import.
