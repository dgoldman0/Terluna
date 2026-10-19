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
16 px/deg with a check on four coasts at 100 m; the tide, the load's balance and Jura Island's worked case stay open.
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
