# How the fresh seas and lakes water their coasts

18 October 2026, on `domain/climate`. The author set the resort at 301 K that day
([register](../../decisions.md#the-target-climate-17-october-2026)), and the warmer Moon still keeps deserts. The author
hopes the shores' groundwater and plants offset some of them, asked whether anything had tested how well the fresh seas
and lakes water their coasts and islands under the Moon's own physics, and approved two steps: the groundwater within
plants' reach under the deserts by the seas, islands and lakes, and the dew and fog of the 354-hour night on the shores
against inland; with a plan, unrun, for a third, a two-dimensional CM1 slice from sea to dry land. This study takes
them at 301 K. [run.py](run.py) is steady groundwater flow and surface energy balances on the GCM's output and
geography's products ([results/coast_water.json](results/coast_water.json)); no climate or cloud model runs. It answers
part of [resources' question 16](../../../resources/questions.md), the shores' water table, Moon-wide at LOLA's
16 px/deg with a check on four coasts at 100 m; [the monthly tide, desert coasts at 100 m and the load's
balance](#the-monthly-tide-and-where-the-seas-load-settles) followed on 4 November 2026, and Jura Island's worked case
stays open.
An independent review the same day reframed the deserts' extent and the share within reach.

## What it finds

- **At 301 K interpolation puts 22.5–28% of the land in desert.** The resort is built as climate_targets builds it,
  0.5197 of the way from the design case (295.023 K on the ground) to the hot world's test (306.523 K), and this study
  works cell by cell, so that each cell's mask, rain and uptake agree. Each cell's rain moved linearly in the ground's
  temperature leaves 22.5% of the land under 0.5 mm a day counting every cell's land, and 24.8% on the model's own land
  cells, 5.11 million km² at 16 px/deg; moved by a constant factor for each kelvin it leaves 28.1% (31.5%), 6.42
  million km², close to climate_targets' 27.9% (30.8%), which interpolates the desert shares themselves. Two runs of
  mixed physics cannot choose between the forms. Of the 5.3 points between the linear form and the interpolated
  shares, 30–45° carries 3.2: it is 33.9% desert in the design case and its rain doubles from 1.52 to 3.01 mm a day
  in the hot world, so the linear form leaves 3.6% of it desert at 301 K against 17.8% by the shares (10.3% by the
  factor form). 45–60° carries 2.4: its rain rises from 0.24 to 1.22 mm a day, and it is 38.0% desert against 52.1%
  (62.1%). 15–30° adds 0.2 and 60–75° takes back 0.4, and all the land poleward of 60° is desert in both forms. The
  linear form's curve flattens near the resort, falling 0.9 points a kelvin over 301–302 K where the factor form falls
  1.8 and the shares 1.6, so how much a warmer resort shrinks the deserts depends on the form as well. The study
  carries three masks, the linear, the factor and the share-matched (the linear form's land under 0.82 mm a day, which
  gives the interpolated shares), and the belts' shares move by under a tenth between them.
- **The deserts hold 175,000 km of shore on the linear mask, half of it round small lakes.** At 16 px/deg the Moon's
  seas, islands and lakes have 1.27 million km of shore. Within the linear mask lie 33,800 km of sea coast, 27,400 km
  of island coast, 88,000 km round the common-level lakes (the atlas's small waters below the sea level, in closed
  hollows), 20,800 km round open rain-fed lakes and 5,000 km round closed ones, with 29 km² of desert behind each km
  of shore. The factor mask holds 206,000 km, the extra mostly round common-level lakes at 30–60°.
- **Their belts are a few hundred metres wide.** Plants reach the water table from 15–25 m above it, and the shore's
  water feeds it from the edge. In the linear mask's deserts at the central crust (5.2×10⁻⁴ m²/s at the ground):

  | Desert shore | Shore | Strip width beside the water: median (10th–90th percentile) | Belt per km of shore | Belt | Supply per km of shore |
  |---|---|---|---|---|---|
  | Sea coasts | 33,800 km | 186–300 m (56–93 to 576–788 m) | 254–362 m | 8,600–12,200 km² | 0.22–0.32 million m³ a year |
  | Island coasts | 27,400 km | 197–314 m (58–96 to 577–782 m) | 259–364 m | 7,100–10,000 km² | 0.24–0.33 million m³ |
  | Common-level lake shores | 88,000 km | 135–221 m (38–63 to 568–786 m) | 223–319 m | 19,700–28,100 km² | 0.18–0.25 million m³ |
  | Open lake shores | 20,800 km | 416–648 m (53–88 to 1,181–1,525 m) | 262–348 m | 5,400–7,200 km² | 0.10–0.13 million m³ |
  | Closed lake shores | 5,000 km | 309–493 m (53–89 to 1,155–1,491 m) | 391–552 m | 1,900–2,700 km² | 0.15–0.22 million m³ |
  | All | 175,000 km | | | 42,700–60,200 km² | 32–46 km³ a year in all |

  The ranges run from a reach of 15 m to 25 m. The strip widths are the belt on each pixel's slope from the
  shoreline, before the pixel's depth caps the strip; the belt per km of shore is the belt's area over the shore's
  length, after the caps and with the pixels inland. The slow crust (2.5×10⁻⁴ m²/s) narrows the sea coasts' median
  strip to 176–277 m and sandy ground (5.2×10⁻³) widens it to 197–326 m: most desert shores rise steeply enough that
  the ground's height, more than its supply, ends the belt. The sandy 25 m case widens the 90th percentile to 1.7 km.
- **About 1% of the desert land lies within reach at 16 px/deg, and 0.4–1.5% after the 100 m check.** The belts cover
  0.84–1.18% of the linear mask's desert land at the central crust and 0.72–1.57% across the crust's range, counting
  every cell's land; 0.59–0.83% (0.51–1.15%) on the model's own land cells; 0.78–1.10% (0.67–1.46%) on the factor
  mask; 0.85–1.21% (0.73–1.60%) on the share-matched deserts. By band at the central crust: 2.0–2.8% at 30–45°,
  1.0–1.5% at 45–60°, 0.45–0.64% at 60–75° and 1.7–2.3% at 75–90°, where the common-level lakes stand thickest, their
  area 14.9% of the land's. On four coasts mapped at 100 m, with one uptake for both grids, the 16 px/deg rule gives
  1.03–2.66 times the belts resolved there (1.10–2.13 at the central crust), above 1 in every case. Each case's share
  divided by its own range of ratios puts within reach **0.4–1.5% of the desert land** (0.40–1.07% at the central
  crust), and 0.36–1.41% on the factor mask, on the assumption that the ratios carry over to the deserts: the four
  coasts are sea coasts at 2–37°, outside the deserts, while two-thirds of the desert land lies beyond 60°, half its
  shore is common-level lakes, and its pixels are a half to a tenth as wide east-west. The 100 m belts converge only on
  western Procellarum, whose 200 m belts are 0.98–1.02 times its 100 m ones; on the other three coasts 200 m gives
  1.07–1.31 times the 100 m belts, so there the belts go on narrowing with finer spacing.
- **A kilometre of desert shore waters about a hundredth of the desert behind it.** At 301 K the plants' demand over
  the deserts is 2.1 mm a day, FAO-56's reference surface under the GCM's own sunlight, air and wind, of which the
  groundwater must supply 1.9 mm (0.71 m a year) after the deserts' 0.17 mm of rain; it is 1.65 mm a day in the design
  case and 2.54 in the hot world, below geography's design metre a year under the deserts' low Sun. A flat
  desert shore at the median uptake feeds a belt 0.84–1.08 km wide at the central crust (0.58–0.75 km in the slow
  crust, 2.6–3.4 km in sandy ground), 0.59–0.76 million m³ a year per km; the desert sea coasts, rising as they do,
  supply 0.22–0.32 million. The 29 km² of desert behind each km would use 22.5 million m³ a year if all of it were
  planted.
- **Closed lakes in the deserts give up area to green their rims.** The 860 closed desert lakes all lie at 75–90°.
  Round the median one the belt draws 1.9–2.9 times the lake's own evaporation at the central crust (1.7–4.4 across
  the range). The ratio leans high on two counts: the lake's evaporation is its row's sea cells' at 301 K, 0.9–1.5 mm a
  day under their marine air, where a lake in the deserts' dry air evaporates more; and the median lake is 0.59 km²,
  two pixels, where the strips' shore is crudest. Together the belts draw 0.76–1.08 km³ a year against the lakes'
  4.1 km³ of evaporation, most of which the large lakes carry: a small closed lake keeps its green rim only by
  shrinking or by more inflow.
- **Small islands can be green throughout.** Of the Moon's 2,474 islands, holding 383,600 km² of land, those under
  10 km² have 64–78% of their land within reach at the central crust, and 394–575 of the 1,123 lie wholly within it,
  their strips covering 95% of their land or more; islands of 10–100 km² 22–31%, larger ones 2–14%. The small islands
  average 4.4 km², one to four pixels, where the strips are crudest and the 100 m check has no island, and geography's
  screen keeps a long island wooded throughout only to about a kilometre across at the central crust; a one-pixel
  island at low latitude is 1.9 km across. Taken as discs of their own area whose highest pixel and the fall
  at their middle under uptake everywhere, Ea²/4T, stay within the reach, 189–387 of the small islands are green
  throughout. In the deserts, 170–253 of 507 small islands are green throughout by the strips and 90–170 by the discs,
  and island land as a whole is 4.5–6.3% within reach.
- **The 15 km water table under the deserts stays deep.** In geography's mature central state the water table stands
  at the ground on 0.37% of the rain mask's desert land and within 25 m of it on 2.8%, a median 774 m down: 1.7% at the
  ground at 30–45°, 0.7% at 45–60°, 0.2% at 60–75° and 0.1% beyond, with 5.3 km³ a year coming out. Its shallow land
  is low ground where the water table stands near the seas' level, only 22% of it in cells beside a sea or lake;
  the solver holds it there without the plants' draw, which from tens of kilometres the crust could not supply. In
  the young Moon, at 1,000 and at 10,000 years, the water table stands at the ground on 0.08% of the desert land and
  within 25 m on 0.5%, a median 1,320 m down.
- **Dew and fog give the desert shores almost nothing.** At 301 K the GCM's nights over the deserts hold the air at
  2 m 60% humid over the shore land of its water cells, 59% on its coastal land cells and 57% inland, the dew point
  8.6–9.5 K below the air. The dense, tall air sends back most of the ground's infrared, so the night's net loss is
  11–24 W/m² and a canopy cools only 0.5–1.7 K below the air: it gathers 0.001 mm a lunar night against 4.8–5.5 mm of
  rain a lunar day, and 0.012 mm with the air taken at the ground's own temperature. On the deserts' mean night a
  canopy begins to gather dew above 89.5% humidity inland, 90.5% on the coastal land cells and 97.2% under the water
  cells' marine air, whose smaller deficit cools it less. The night's means saturate at the ground over 0.4% of the
  desert land; with the 6.5-hour means' threefold fog, its water is at most 0.16 mm a lunar night with Earth's least
  deposition onto short vegetation, 2.1 cm/s, and 0.60 mm with its most, 8 cm/s, 12% of the deserts' rain.
- **Wet ground makes its own dew.** CM1's rings on the design ground run their land nights at 75–82% humidity, 3.2–4.6
  K above the dew point, and a canopy there gathers none. Over the earlier rings' ground evaporating 2 mm a day the
  nights run 96–97% humid and a canopy gathers 2.6–2.7 mm a lunar night, 3.3 mm on the equator's: dew of the wet
  ground's own water, returned.
- **T21 resolves no coast, and CM1 puts the fall-off within the first 6 km.** Each GCM cell's air is all land or all
  water, and at night the water cells' air is about as humid as the land's (60% against 59% and 57%). The Jura box, at
  6 km, holds the land's dew point 4.8 K below the sea's from the first column inland and level to 25 km and beyond,
  with an afternoon sea breeze of about 1 m/s and a night land breeze of 0.18–0.35 m/s offshore. The humidity falls
  off within the first few kilometres, finer than both models.

**What it means.** The shores' water comes from the ground. At 301 K it can green a rim of a few hundred metres round
the desert seas and lakes, and the smallest islands throughout: 0.4–1.5% of the desert land, drawing 32–46 km³ a year
from the seas and lakes through the crust. The night adds nothing to it, so the deserts stay 22.5–28% of the land,
as the form of the interpolation leaves them, with green rims.

## A two-dimensional slice: the plan

**What it would run.** Two two-dimensional CM1 slices, one from a desert sea coast and one from a desert lake shore,
each running from open water across a shore belt to dry inland: 150 km of water, the belt kept wet as wide as this
study's belts there (a few hundred metres to 2 km), and 300–400 km of land at the GCM's dryness, periodic with the
coast mirrored. The air at 301 K, held above 8 km to a 45–70° row of the design case at 301 K, the run the target's
choice already calls for (about 1.5 hours on four ranks); the Coriolis force of the latitude, a 27th of Earth's at 45°;
the GCM row's sinking imposed; the Sun's lunar day. A lunar day of spin-up and a lunar day analysed, as the rings ran;
CM1's horizontal stretching for the shore's spacing, a setup change for the land and sea laid along x.

**What it would answer.** How far the long day's sea breeze carries the sea's moisture inland when rotation barely
turns it; whether the breeze's front raises storms that rain on the coast; how strong the night's land breeze is, and
whether dew forms on the belt and beyond it; and whether a wet belt makes its own humid nights.

**What it would cost**, scaled from the rings' and the Jura box's own records on 6 threads (the rings take 273–342 s a
model day for 1,288 columns, a step of 69–75 s with 8–9 acoustic sub-steps; the box 965 s for 9,000 columns):

| Slice | Columns | Step | Wall hours each, 59 model days |
|---|---|---|---|
| 2 km throughout | 300 | 73 s | 1.3–1.6 |
| 500 m at the shore, stretched | 212 | 20 s | 3.4–4.2 |
| 1 km throughout | 600 | 40 s | 4.8–5.9 |
| 250 m at the shore, stretched | 263 | 10 s | 8.3–10.4 |

The finer slices' steps follow the horizontal Courant number with 25 m/s aloft. Only 500 m or finer resolves the
belts; the 2 km slice, about a tenth of the Jura box, resolves the breeze and the storms.

**Recommendation: run the pair at 500 m once the design case at 301 K exists and the coasts are being designed.** It
gives the coasts' weather: the sea breeze and the storms at its front, the land breeze and the belts' nights, which set
the resort's coastal comfort and the shore cities' climate, in about 7–8.5 hours for the pair. Steps 1 and 2 already
answer the deserts' water: the belts' reach is set by the ground, which no atmosphere changes; dew needs nights 90%
humid or more where both models give 57–82%; and the Jura box holds the land's dew point 4.8 K below the sea's from
the first 6 km. The slice would add to the deserts' water only if its breeze carried the sea's air tens of kilometres
inland.

## How it is set up

**The resort.** [extract.py](extract.py) reads the design case (`A28_corrected_970`) and the 20-year test
(`A28_corrected_970_lw_ice`) over years 10–19, 1,200 3-day means each, and the design case's two years written every
6.5 hours (`A28_corrected_970_eddies`). Each cell's rain, the plants' demand and each night figure are found in each
run and taken across cell by cell, a share 0.5197 of the way, the share climate_targets' resort band carries at 301 K:
linearly, and for the rain also by a constant factor for each kelvin, rd^(1−f) rt^f. Each mean's humidity is that of
its start and its end, and the day's means are those the Sun is up for more than half of, both as climate_targets
reads them. The share-matched mask takes the linear form's rain under the threshold midway between the last cell that
brings the land to climate_targets' share and the next.

**The shores.** [shores.py](shores.py) takes LOLA at 16 px/deg above the degree-200 GRAIL geoid, as geography's
drainage does, against the atlas's sea level: every pixel below it is common-level water, and the drainage product's
lake cells, with their levels and whether each spills, are the rain-fed lakes. Common-level bodies of at least the
atlas's listing size (5×10⁻⁴ of the Moon) are the 11 seas, 8.36 million km², and the other 5,697 the common-level
lakes, 2.26 million km²; the 27,781 rain-fed lakes cover 3.92 million km², 2,493 of them closed. The largest land body
is the main land and the others the islands, each with its highest pixel and its land's mean height. Each land pixel
within 12 km of water keeps its distance from the nearest water pixel, by distance transforms in bands of a degree
with their own east-west spacing, its height above that water's level, the water's class and its island; each pixel
beside the water keeps its edges facing each class.

**The belts.** [belts.py](belts.py) solves steady Dupuit flow under a belt fed from the water's edge, T h'' = E, the
head at the water's level at the shore and no flow at the belt's edge, where the deep water table stands flat. The
water table at the edge falls cW² below the water's level, c = E/2T, and the belt ends where the ground's height and
that fall reach the reach: sW + cW² = R on a shore rising at a slope s. A pixel beside the water holds a strip along
each water-facing edge as wide as the belt on the slope from the shoreline to its centre, at most its depth, the edges
counting π/4 as a staircase does for shores of every direction; pixels inland hold the share of their depth the belt
on their own slope from the shoreline covers. The uptake E is each cell's demand less its rain in the linear form,
for every mask. The transmissivity
takes groundwater stage 2's slow and central crust with the water table at the ground, 2.5 and 5.2×10⁻⁴ m²/s, and ten
times the central for sandy ground, as geography's shore screen does; the reach is that screen's 15–25 m, 3–5 m of
roots and 12–23 m of capillary rise at lunar gravity. At the deserts' 1.9 mm a day Gardner's sandy loam lifts the
water from 13.0 m and the finer soil from 26.0 m, so the reach holds at 301 K. An island lies wholly within reach
when its strips cover 95% of its land.

**The check at 100 m.** The same rules run on the four coasts geography resampled from LOLA's 256 px/deg to 100 m
(western Procellarum, the Ingenii coast, southern Nubium and the Smythii headland), over their inner squares 114–120
km across, in blocks of 200 and 400 m as well, and on the 16 px/deg pixels in the same squares, all with the uptake
of the tile centre's cell.

**The demand.** FAO-56's Penman-Monteith for the short reference surface (albedo 0.23, aerodynamic resistance 208/u2)
with the hourly convention for periods of day and night: surface resistance 50 s/m by day and 200 by night, the soil
taking a tenth of the net radiation by day and giving back half its loss by night, the psychrometric constant and the
air's density at the local pressure. The surface takes the GCM's sunlight at the ground and its sky's longwave, the
net longwave less the ground's own emission (emissivity 1 on land and 0.98 on water, as PlaSim's radiation takes
them). The wind at 2 m is the lowest level's, about 0.9 km up, times the factor by which a 3-day mean's vector
understates its 6.5-hour speeds (1.14 by day and 1.21 at night over land), the CM1 rings' 10 m to 1 km ratio on the
design ground (0.55 by day, 0.28 at night over land) and FAO-56's 0.748 from 10 m to 2 m. The 6.5-hour means give the
same demand as the 3-day means, 2.572 mm a day over the design case's land.

**The night.** ExoPlaSim forms no dew of its own: its surface flux can only raise the lowest level's humidity
(fluxmod.f90, line 664). [dew.py](dew.py) balances a grass-like canopy against the sky's longwave, half its loss drawn
back from the soil, taking heat and vapour from the air at 2 m through the aerodynamic resistance and its leaves'
boundary layers; an exposed leaf through its own boundary layer alone; and bare ground at the GCM's own temperature.
The air at 2 m is the lowest level's carried down the dry adiabat, which reads warm in calm nights, and a variant takes
it 0.3 K above the ground, as the rings' land nights run. Fog is the night's means saturated at the ground's
temperature, its water 0.01–0.5 g/m³ reaching the plants at 0.05–2.1 cm/s, from a 10 µm drop's settling at lunar
gravity to Earth's least deposition onto short vegetation, and at Earth's most, 8 cm/s. The rings' nights take each
ring's own surface pressure from its case. On an Earth night of 10 clear hours under an 87 W/m² deficit the canopy
gathers 0.21 mm at 90% humidity and 0.28 mm at 96%, between the 0.13 ± 0.10 mm a night typical of grassland and the
lysimeters' 0.5–0.6 mm most (Ritter et al. 2019).

**The water table.** Geography's 2 px/deg grids of the mature central state and the young states, each cell
attributed to its T21 cell, the deserts' cells chosen by the masks.

**The classes at T21.** The shore land in the model's water cells, under their marine air; the model's coastal land
cells, with a water cell among their eight neighbours; and its inland land cells. The 16 px/deg shore classes take the
air of the cells they lie in.

## Assumptions and limits

- **The resort mixes two physics.** The design case flew PlaSim's earlier longwave and the test the refitted one, and
  a design case at 301 K settles the deserts' extent, their nights and their demand.
- **The deserts' extent depends on the form of the interpolation.** Each cell's rain moved linearly gives 22.5%
  (24.8%), by a constant factor a kelvin 28.1% (31.5%), and the interpolated shares 27.9% (30.8%); the belts' shares
  move by under a tenth between them.
- **The belts are a screen at 16 px/deg.** Pixels 1.9 km across smooth cliffs and shore relief. The 100 m check
  ranges 1.03–2.66 on four coasts, none of them in a desert, and the 100 m belts converge on one of them; the
  corrected share takes the ratios to carry over to the deserts. Beaches, dunes, sediment and the developed soils,
  which set the shores' own permeability and their slopes at the water, are left out; the transmissivity is the
  crust's at the ground, with sandy ground its high end. The monthly tide, 3.7 m typical on the Near-side Sea, is left
  out, and so are the load the plants' draw leaves, the curvature of small islands, which feeds a round island better
  than a strip, and the interplay of belts round one lake.
- **The common-level lakes are the atlas's.** Every closed depression below the sea level is water in the atlas and
  the GCM; in the polar deserts their own evaporation would have to come through the crust, which the 15 km solver
  holds as fixed heads.
- **The lakes and the water table are the earlier design climate's.** Geography's lakes and the 15 km water table
  were built on `A28_dim5_moon`'s rain, with 42% of the land desert, and regenerate once a design case at the target
  exists.
- **The belts are the mature Moon's.** In the first millennia a frozen layer seals the polar crust, 1.5 km thick under
  the polar highlands and thawing over 18,000–250,000 years poleward of 60° at the design climate, and the shores'
  groundwater moves only through the thawed top.
- **The night rests on 3-day means in cells 170 km across.** The design case's 6.5-hour means give the deserts ten
  times the 3-day means' dew and three times their fog, both still small against the rain; the GCM's humidity
  is the lowest level's, 0.9 km up, mixed down, and over the sea the air near the water is moister.
- **The plants are FAO-56's reference grass.** Engineered plants with other resistances, leaves or roots shift the
  demand and the reach; the demand moves the flat shore's width as its inverse square root.

## The checks

[test_run.py](test_run.py) checks the belt's width against its limits and geography's screen (0.56–0.79 km at
5×10⁻⁴ m²/s with 10–20 m to spare, 1.8–2.5 km at ten times, 8 m under a kilometre-wide island), the steady flow
against finite differences of T h'' = E, belts.py's rules on a straight rising shore on a Cartesian grid against the
analytic belt, the shores' geometry (the banded distance transforms against brute force across the seam, the labels
joined across the seam and the poles, the edges' lengths and depths), the 16 px/deg pipeline from distances and edges
to belts on a synthetic strip of land at the equator against its own rule's analytic figure, the canopy's balance
closing and its dew rising with humidity, the dew on an Earth night against measured dew, the demand against FAO-56's
own daily equation, the dew point against saturation, the product's binding to its code and inputs, the basis
against climate_targets, the share-matched masks against climate_targets' shares, the case-matched correction, the
water table's grids against geography's product, the figures above against the product, and, with
TERLUNA_COAST_WATER_REGENERATE=1, the product's regeneration. These are runnable calculations with their inputs bound
by sha256. The rules agree within a tenth with the analytic belt on a synthetic shore at 100 m spacing, the local
coasts' sequence from 400 to 100 m converges on one coast of four, and the Earth dew check is the one comparison with
measurement.

## Running it

```sh
climate/gcm/.venv/bin/python -m research.studies.coast_water.extract   # about a minute; needs netCDF4
python -m research.studies.coast_water.shores                          # about twenty seconds, 1.5 GB
python -m research.studies.coast_water.run                             # about seven minutes, 2.4 GB
python -m pytest research/studies/coast_water                          # seconds; seven more minutes to regenerate
```

The extract writes `research/runs/coast_water/gcm_extract.npz` (578 MB) with a sidecar holding the sha256 of every
GCM file it read, and the shores `research/runs/coast_water/shores_16ppd.npz`, both ignored and on the research drive
(`/media/projectspace/terluna-research/coast-water/runs`). The run reads them, geography's water table grids in
`research/runs/resources/groundwater`, the four local coasts in `research/runs/sea_appearance/terrain` and the CM1
runs' records and cases in `climate/crm/runs`, and binds each by its sha256. [sources.json](sources.json) records
each published value and how far it was read.

## The monthly tide and where the sea's load settles

4 November 2026, on `domain/hydrosphere`. The seas rise and fall with the Earth's monthly tide, 3.7 m in a typical month over the Near-side Sea and 5.5–6.0 m in its southern bays ([geography's tide](../../../geography/README.md#the-monthly-tide)), and every shore's belt draws its water from them. This section puts that tide under the belts above, on the same 301 K basis and the same 16 px/deg shores, and follows the seas' dissolved load into the ground the plants draw from. [tide.py](tide.py) is a cross-shore groundwater model, [tide_run.py](tide_run.py) carries it to every shore, [tide_terrain.py](tide_terrain.py) brings in LOLA's 256 px/deg rows across desert sea coasts at 45–60°, and [load.py](load.py) follows the load ([results/coast_tide.json](results/coast_tide.json)).

### What it finds

- **The monthly tide reaches 45–140 m into the ground and arrives days late.** In linear theory (Ferris 1951) the 27.55-day line's swing falls to 1/e within 44–62 m in the slow crust, 63–89 m in the central crust and 121–137 m in sandy ground, across each ground's specific yield, and lags the sea by 5–10 days at 100 m in the crust and 3–4 days in sand. Lunar gravity enters through the transmissivity: the ground conducts a sixth as well as on Earth, and the month is 53 times Earth's semidiurnal period, so the monthly tide reaches about three times as far into the same ground as an Earth tide does. On a steep face the cross-shore model puts the swing at 75 m at 0.365 of the sea's, against Ferris's 0.367, and its lag at 4.39 days against 4.39.
- **On a sloping shore the tide raises the mean water table inland, towards the high water.** The rising sea enters the face over its whole wetted slope, and the falling sea leaves only through a seepage face, so the face takes the sea in faster than it lets it out. Without the plants' draw the model gives Nielsen's (1990) rise, εA/2 with ε = kA/slope, within 3% up to ε = 0.3; beyond it the rise saturates, at 0.42 of the tide's amplitude at ε = 1, 0.78 at 4 and 0.94 at 20, where the inland water table stands near the high water.
- **The plants' draw makes the high-water strip the belts' feed.** The monthly high water spreads A/slope inland, 80 m on a 1:50 shore at the Near-side Sea's typical half-range of 1.6 m and 270 m on 1:170, and the ground under that strip stays nearly saturated through the month, since the month's drainage reaches only a damping length into it. A belt drawing on the water table is then fed from near the high-water line. In the central crust with 15 m of reach, a draw of 2 mm a day and that tide, a belt 486 m wide on a 1:50 shore widens by 7% and one 694 m wide on 1:170 by 29%; on 1:14 it widens by 1.7% and on 1:4 by 0.5%, and its supply grows with its width. The slow crust's belts gain 0.6–43% across those slopes, and sandy ground's, wider to begin with, 0.2–11%.
- **Where the plants draw, the water table swings by decimetres to a metre.** At that tide the water table's monthly range is 0.8–0.9 m 100 m inland and about 0.2 m at 200 m on shores steeper than 1:15 in the central crust (1.5 m and 0.7 m in sandy ground), and stays near the sea's own range within the strip the high water floods. On the steepest shores the belt's edge lies within the tide's reach: on 1:4 the strip within reach at every low water is 57 m wide against the 60 m belt without the tide, so plants at its edge lose the water table for days each month.
- **Moon-wide the tide widens the desert belts by about 2%.** On the linear mask's deserts the belts grow from 42,730 to 43,620 km² at the central crust with 15 m of reach, by 2.1%, and their share of the desert land from 0.84% to 0.85%; across the six cases by 0.2% (sandy ground, 25 m) to 3.4% (slow crust, 15 m), and by 1.8% and 2.4% on the factor and share-matched masks. Only the sea and island coasts gain, by 5.7% and 5.6% at the central crust (0.6–9.5% across the cases). The desert sea coasts' typical monthly range is 2.8 m at the median (1.0–3.7 m), and their median rise is 0.17 m, 3.3 m at the 90th percentile, where the gentle pixels lie: the median strip barely moves, from 186 to 188 m, while the 90th percentile grows from 576 to 654 m. By band the share within reach moves from 2.01% to 2.07% at 30–45°, 1.05% to 1.09% at 45–60°, 0.45% to 0.46% at 60–75° and stays at 1.73% at 75–90°, where the common-level lakes' shores hold the belts. The belts draw 33.1 km³ a year against 32.2. The 100 m check's correction leaves the deserts' share within reach at 0.4–1.5%.

  | Desert shore, linear mask | Belt without the tide, central crust, 15 m | With it | Change, central (six cases) | Supply per km of shore, central |
  |---|---|---|---|---|
  | Sea coasts | 8,594 km² | 9,083 km² | +5.7% (+0.8 to +9.2%) | 0.221 to 0.234 million m³ a year |
  | Island coasts | 7,085 km² | 7,481 km² | +5.6% (+0.6 to +9.5%) | 0.237 to 0.251 million m³ |
  | Common-level lake shores | 19,664 km² | 19,667 km² | +0.02% (0 to +0.03%) | 0.176 million m³, unchanged |
  | Open and closed lake shores | 7,388 km² | 7,388 km² | under 0.01% | unchanged |
  | All | 42,730 km² | 43,620 km² | +2.1% (+0.2 to +3.4%) | 32.2 to 33.1 km³ a year in all |

- **The lakes carry their own tides, decimetres only in the large ones.** Each lake follows its own equilibrium tide with its volume held, tilting with the Earth's pull. The 50 rain-fed lakes over 10,000 km² reach a median greatest typical monthly range on their shores of 0.51 m, up to 1.38 m round the largest, 200,000 km² near 1° N, 232° E (Hertzsprung), and 1.17 m round the 141,000 km² lake near 4° S, 202° E (Korolev); lakes of 1,000–10,000 km² reach 0.17 m (up to 0.93 m), 100–1,000 km² 5 cm (up to 20 cm), 10–100 km² 1 cm and the 15,013 under 10 km² under 2.4 cm. 796 rain-fed lakes reach a decimetre somewhere on their shores, holding 2.94 of the rain-fed lakes' 3.92 million km², the smallest of them 172 km²; of the common-level lakes, all under the seas' threshold of 19,000 km², 493 reach a decimetre and the largest 0.61 m. In the deserts the common-level lake shores' median range is 5 cm and the rain-fed lakes' under 1 cm, so their belts move by under a part in a thousand.
- **At 100 m on desert sea coasts at 45–61° the tide adds 0–3% to the belts, a fifth of what the 16 px/deg screen gives the same squares.** Eight 12 km squares at the densest desert sea-coast degrees of the four new LOLA bands, with typical monthly ranges of 0.23–3.73 m, hold belts of 1.5–7.9 km² at 100 m in the central crust with 15 m of reach. The tide adds 0–2.9% to them, 1.5% over the eight, where the 16 px/deg strips in the same squares gain 0–16%, 6.7% over the eight: LOLA's rock at the water rises more steeply than the 2 km pixels' mean slope, so the high-water strip is narrow, and the sea coasts' Moon-wide gain is nearer 1–2% at this resolution. The 16 px/deg strips there are 1.3–2.5 times the 100 m belts, as on the coast-water study's four coasts. On 32 transects every 10 m the belts end a median 107 m from the water (32–1,104 m), the tide moves that edge by under a metre at the median, and the water table's half-range 100 m inland is 0.03–0.74 m. The squares stand for the desert sea coasts at 45–61°: the Near-side Sea's northern shores (Mare Frigoris at 60.5° N, 32.5° E and 50.6° N, 36.5° W; 48.7° N, 43.3° E; 50.6° N, 74.6° W), its north-eastern limb at 48.7° N, 84.5° E, the far side's 60.7° N, 97.5° E, and the South Pole–Aitken Sea's shore at 57.5–57.6° S, 128.6° and 147.5° E. They leave out the desert coasts poleward of 61°, where most of the desert sea coast lies and no 256 px/deg rows were fetched, the deserts at 30–45°, the islands' coasts and every lake shore. Each profile's drift over its last ten months is recorded: in four squares a central-crust transect had not settled when it stopped, its belt's edge still moving 18 m in those ten months at 60.7° N, 97.5° E (its mean water table 0.43 m), 3.2 m at 57.5° S, 147.5° E (0.21 m), 2.8 m at 50.6° N, 74.6° W (0.14 m) and 1.8 m at 60.5° N, 32.5° E (0.13 m); every other transect moved under 0.6 m and 6 cm.
- **The held seas lay calcite, then opal, at the roots.** The plants draw the shore's water from the capillary fringe and keep almost none of its salts, so each cubic metre they transpire leaves its load in the root zone. The held fresh sea, 0.84 mmol/kg of calcium with twice that in bicarbonate and 0.1 mmol/kg of silica at 0.14 g/kg, stands at the chemical divide of evaporating waters (Hardie and Eugster 1970): its calcite saturates at once in the air's CO2, where the sea already stands at three times saturation, and at a concentration of 1.9–2.9 times in a root zone of 1–3% CO2; amorphous silica follows at 20 times; all its calcium leaves as calcite and all its silica as opal before the water is gone, and nothing is left to salt the soil. The left-alone sea as it passed 0.14 g/kg carries magnesium and sodium besides, with more bicarbonate than twice its calcium: calcite and opal come first, magnesium clays where the water passes pH 9 (at 17 times in the air's CO2, at 470 times in 1%), and a sodium carbonate brine to the end.

  | Desert shore, central crust, 15 m | Draw a year | Hardpan from the held sea, solid; cemented soil at 30% porosity | From the left-alone sea, cemented; a sodic metre of root zone | With the delivered chloride: a saline root zone; halite a century |
  |---|---|---|---|---|
  | Sea coasts | 0.87 m | 2.9 mm a century; 9.8 mm | 13 mm a century; 380–1,900 years | 1.8 years; 40 kg/m² |
  | Island coasts | 0.92 m | 3.1 mm; 10.4 mm | 14 mm; 360–1,800 years | 1.7 years; 42 kg/m² |
  | Common-level lake shores | 0.79 m | 2.7 mm; 8.9 mm | 12 mm; 420–2,100 years | 1.9 years; 36 kg/m² |

  The draw is each class's supply over its belt area. The sodic layer's years take an exchange capacity of 5–25 cmol/kg in the root zone's metre; the saline root zone is its soil water at 4 dS/m. The rain-fed lakes' shores draw lake water, 0.04 g/kg in the held state, about 0.28 of the sea's load at their 0.37–0.39 m a year; closed lakes in the deserts turn to soda over 10⁴–10⁵ years (the [sea chemistry](../sea_chemistry/README.md)), and their rims salt with them. At the held sea's 9.8 mm a century the root zone cements 10 cm of soil in about a thousand years, within both groundwater calcretes' 10–3,000 mm a thousand years and pedogenic calcretes' 0.5–350 (Goudie 1985, through Fenske et al. 2025).
- **The tide and the rain carry little of the load back.** Where the high water lifts the water table within 3 m of the ground, the roots stand in groundwater that moves with the tide, and the soluble salts they leave there can go back to the sea with it: on the planar shores at the typical tide that strip holds 15–40% of the plants' draw across the grounds and slopes, 24–34% in the central crust with 15 m of reach. The calcite and opal stay where they form, so the tide can return at most that share of a soluble load and almost none of the held sea's, which is all hardpan. The rain flushes none in the deserts' mean climate: their 62 mm a year falls on belts drawing 870 mm, the plants use it, and nothing passes below the roots. Holding the root zone below the saline threshold with the delivered chloride would need 0.17 m a year of water passing below the roots, nearly three times the deserts' rain.

**What it means.** The tide gives the desert shores a little more water: about 2% more belt in all at 16 px/deg and less at 100 m, most of it on gentle sea coasts where the monthly high water spreads over the land, at the cost of a few dry days a month at the steepest belts' edges. The deserts' green rims stay where the ground's height ends them. The tide does not clean the shores. The held seas' calcium and silica stay at the roots as a calcrete hardpan, about a centimetre of cemented soil a century under every desert belt, which roots and people must break or outgrow over centuries to millennia; and if the seas carried the delivered water's chloride, the desert belts' root zones would turn saline within two years, with neither the tide nor the rain to flush them. Green desert shores then need that chloride cleaned from the water ([question 17](../../../resources/questions.md)) or about 0.17 m of water a year flushed below each belt's roots.

### How it is set up

**The basis.** [tide_run.py](tide_run.py) recomputes this study's 301 K basis with run.py's own functions: the three desert masks, the plants' uptake by T21 cell, the 16 px/deg shores and the six belt cases, and its belts without the tide reproduce run.py's. The plants' draw through the lunar day uses the same basis's split of the deserts' demand, 87% of it by day; no new climate input was needed, so nothing departs from the 18 October basis.

**The tide on each shore.** Each sea and island coast takes half the typical monthly range of the nearest node of geography's tide grid (the equilibrium tide with each sea's volume held and its own water's attraction, 2026–2045). Each lake takes its own equilibrium tide as a closed basin of fixed volume: the tide lines' forcing, recovered from the grid in geography's monomials (rank 15, the combinations constant or vanishing on the sphere left out; residual 9×10⁻⁸), less its mean over the lake, turned into a typical range by the ratio of the sea nodes' typical range to the rms of their lines, 2.84 (2.78–2.97). The lakes are labelled as shores.py labels them; their own water's attraction, a few per cent in the seas of their size, is left out. In Mare Fecunditatis geography's equilibrium tide stands in for the strait study's, due shortly, which can replace it there.

**The cross-shore model.** [tide.py](tide.py) solves unconfined Boussinesq flow across a shore by backward Euler on a uniform grid: nodes under the sea held at its level, so the waterline moves up and down the slope; nodes whose water table would rise above the ground seeping at the ground; the plants drawing the cell's uptake wherever the water table lies within the reach, nothing beyond. The transmissivity is this study's; the specific yield is 5–10% in the crust (Hunt 1996's Oahu lavas; Woessner and Poeter's fractured and vesicular basalt) and 21–27% in sandy ground (Johnson 1967's sand averages), their geometric means the central values. Each run repeats the forcing for at least three of the belt's slowest settling times and then until the window means' geometric tail falls under 3 mm. For each of the six cases, planar shores over slopes of 0.006–0.25, tide amplitudes of 0.25–2.8 m and draws of 1.2–3 mm a day, 288 runs, give the rise that the rule s W + c W² = R needs added to the reach to reproduce the model's belt; each 16 px/deg pixel takes it at its own slope, tide and draw, interpolated in the logs.

**The desert coasts at 100 m.** [tide_terrain.py](tide_terrain.py) fetched LOLA's 256 px/deg rows across the four one-degree bands where the linear mask's desert sea coasts at 45–60° are densest (48–49° N and 60–61° N in the eastern half, 50–51° N in the western, 57–58° S in the eastern), 24 MB each by exact byte range, their sha256s in [tide_inputs.json](tide_inputs.json). At the two densest desert sea-coast degrees of each band a 12 km square at 100 m is built as geography's coast terrain builds it; this study's fine belts run on it with and without each pixel's rise (the central crust with 15 m of reach and sandy ground with 25 m), beside the 16 px/deg strips in the same square, and the cross-shore model runs on four transects from the shoreline along the distance field's gradient, every 10 m from 400 m out to 2.5 km in, with the cell's uptake and the site's tide, to three of the belt's settling times and the geometric tail or at most 120 months, the run with the tide starting from the run without it raised by the planar table's rise.

**The load.** [load.py](load.py) concentrates the seas' waters step by step at the air's CO2 and at a root zone's 1% and 3% with the carbon cycle's Plummer and Busenberg constants and Davies activities, calcite at saturation, amorphous silica at its solubility and magnesium clays above pH 9 as the sea chemistry takes them, to the ionic strength of 0.5 where the activities end. The deposits follow the draw: the root zone gathers the draw times the water's concentration, every mineral of it in the end.

### Assumptions and limits

- **The shore's slope is the 16 px/deg pixel's.** The tide's effect depends on the slope over the waterline's own excursion, A/slope, metres to hundreds of metres; the screen takes the slope from the shoreline to the pixel's centre, about a kilometre, and holds it at 0.006 on the flattest 5% of sea coasts. Beaches and sediment, with their own slopes and permeability, are left out, as the waves and run-up are; on the 100 m squares the slopes at the water are the rock's.
- **The specific yield is a drained pore volume measured on Earth.** It holds at lunar gravity for a water table deeper than the soil's capillary length, 3–6 m in geography's Gardner soils at lunar gravity. Nearer the ground the fringe, six times taller than on Earth, keeps water that a falling water table would drain on Earth, lowering the apparent yield and letting the tide reach further; that is left out. Across the crust's 5–10% the swing at 100 m is 0.33 and 0.20 of the sea's.
- **The tide is one monthly line of the typical range.** The two leading lines, 27.55 and 27.21 days, beat over six years; a run with both at a Mare Frigoris node (typical range 3.3 m, lines of 0.58 and 1.44 m) gives the mean belt 0.8% narrower than one line of the typical range and the belt within reach at every low water 3.4% narrower, as the beat's strongest months lower the low water further. The fortnightly and 9.6-day lines are smaller and damp faster.
- **The rise is read from planar shores.** The 288 runs reproduce the analytic belt without the tide within 2%. Six direct runs off the grid, across the grounds and reaches, slopes of 0.01–0.15, tides of 0.5–2.4 m and draws of 1.3–2.5 mm a day, give the belt with the tide within 0.3% of the table's (412 m against 413, 772 against 770, 1,688 against 1,689). Fifteen sandy runs on the gentlest slopes stopped at 300 months with their far ends still moving under 0.5 mm a month.
- **The water is fresh throughout.** The seas, the groundwater and the fringe have one density, and no salt wedge forms, as this study takes it.
- **The lakes' tides are equilibrium tides,** without their own water's attraction or seiches.
- **The load's chemistry ends where the activities end.** Davies activities hold to an ionic strength of 0.5; the soda and halite past it are mass balance. Calcite forms at saturation in the root zone, where the sea chemistry's threefold supersaturation stands for calcifiers in open water. The sodic and saline thresholds, an exchangeable sodium share of 15% and soil water near 4 dS/m, are Earth's for crop plants and reference points for engineered ones; the soil's exchange capacity is a stated bracket. The root zone gathers the load over a metre; deeper or shallower roots spread or concentrate it.

### The checks

[test_tide.py](test_tide.py) checks Ferris's damping against its formula and against Earth's semidiurnal tide, the cross-shore model's swing and lag on a steep face against Ferris's, its rise at small ε against Nielsen's and its saturation below the high water, its belt without the tide against belts.width, its water balance over a month with the tide and the draw (what enters from the sea and the seepage face less what leaves, against what the plants draw, within 3%), the belts with a reach for each pixel against run.belts_for, a lake's tide vanishing at a point, the evaporation paths of the held and left-alone seas and the deposits' mass balance, and the product's binding to its code and inputs and its belts without the tide against this study's product. In the product the steep-face run puts the swing 75 m inland at 0.365 of the sea's against Ferris's 0.367 and its lag at 4.39 days against 4.39; the plants' draw through the lunar day, 87% of it by day, swings the water table about 0.3 m at 100–200 m without the tide and narrows the belt within reach at every hour by under 1%; a finite aquifer 100 m deep narrows the belts with and without the tide by under 1%, its finite-amplitude rise A²/4D 9 mm; and the crust's specific yield at 5% or 10% moves the belt with the tide by under 2%. These are runnable calculations and numerical checks; no measurement of a tidal water table or a phreatic hardpan on these shores exists to test them, and on Earth Raubenheimer et al. (1999) measured mean water tables above the sea across a beach of 0.03 slope, as the model gives them.

### Running it

```sh
python -m research.studies.coast_water.tide_terrain --download         # once, 96 MB of LOLA rows by byte range
python -m research.studies.coast_water.tide_run                        # 7 minutes, 1.8 GB, after the grids (20 minutes once)
python -m pytest research/studies/coast_water/test_tide.py             # about a minute
```

The run reads this study's GCM extract and shores (research/runs/coast_water), geography's tide grid (`geography/products/tides_28pct_4ppd.npz`), its drainage product and LOLA's 16 px/deg heights for the lakes, the sea chemistry's product and the four LOLA slices, and binds each by its sha256. The 288 planar runs are cached in `research/runs/coast_water/tide/numerics.json` under the hash of tide.py and the grid. [tide_inputs.json](tide_inputs.json) records each slice's source, byte range and sha256, and [sources.json](sources.json) each published value and how far it was read.
