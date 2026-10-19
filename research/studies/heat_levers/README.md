# What could take pressure off the dimming

17 October 2026, brought to the corrected target climates on 19 October · `domain/resources`. The author asked whether
anything could reasonably reduce the greenhouse effect, or the heat the plan adds, to take pressure off the solar
shield's dimming, and which options may cost less. This study screens every lever the project's products suggest
against the three target climates of [what each target climate takes](../climate_targets/README.md), whose dimming in
all is 3.7–5.3% for the hot world near 307 K, 8.2–10.9% for the resort at 301 K, where the author set it on 18 October,
and 12.6–16.6% for the lock at 295 K. Each lever is given in W/m² of the Moon's balance and in points of dimming, a
point being 1% of the hot world's absorbed sunlight, 225.4 W/m², as the targets count the heat beyond the sunlight.
[run.py](run.py) is arithmetic on committed products ([results/heat_levers.json](results/heat_levers.json), schema
`terluna.research.heat-levers/2`); no model runs. The targets are the corrected product (schema `/5`) with the GCM
extract it binds. Published values and how far each was read are in [sources.json](sources.json). The array's studies
are read only.

## The methane and N2O pair

The held pair, 1.87 ppm of methane and 0.1 ppm of N2O, forces 3.58 W/m² all-sky in the deep column at 295 K
([reduced gases](../../../atmosphere/radiative_convective/README.md)), 1.59 points. Lower levels save:

| Pair, ppm | Saved, W/m² | Points | Methane contactors, air a year; fans | N2O decomposition | Years to reach from none |
|---|---|---|---|---|---|
| 0.3 and 0.03, the low pair | 1.78–1.82 | 0.79–0.81 | 1.0–1.4%; 0.07–0.10 TW | 0.96–1.47 TW | 66–90; 1,100–1,680 |
| 0.3 and 0.1, methane alone | 1.08–1.09 | 0.48 | 1.0–1.4%; 0.07–0.10 TW | as held | 66–90; as held |
| 1.87 and 0.03, N2O alone | 0.61 | 0.27 | as held | 0.96–1.47 TW | as held; 1,100–1,680 |
| 0.1 and 0.01, the grid's lowest | 2.52–2.61 | 1.12–1.16 | | | |
| 1.87 and 0.1, held | | | 0.03–0.09%; under 0.01 TW | 0.29–0.44 TW | 621–1,200; 3,650–5,610 |

The costs are the domain synthesis's ([ledgers](../domain_synthesis/README.md#where-they-are-held)), with biology at
the sources cutting 90% of the methane and 95% of the N2O. With biology cutting half of each, the low pair's methane
contactors pass 5.8–7.9% of the air a year and its N2O decomposition takes 9.6–14.7 TW. The low end of each saving
carries its own clear-sky trend through the columns to the hot world's ground: the low pair's saving falls from 2.99
W/m² at 290 K to 2.93 at 300 K, 2.2% less at 306.5 K, where the held pair's whole forcing, 5.94 W/m² at 290 K and 5.66
at 300 K, falls 5.7%.

- **The air has no OH**, so methane leaves only into soils and N2O not at all: every tonne biology leaves is taken
  from open air, and a lower level passes more air for the same tonne. The low pair's methane contactors move
  1.0–1.4×10⁹ kg of air a second, the scale of the build's whole delivery stream.
- **The removal's own heat is small.** At the low pair it is 0.03–0.04 W/m² (0.26–0.40 with half cuts), against the
  1.78–1.82 W/m² it saves.
- **The lakes tie the methane to the target.** Ecology's lake model ([the waters' oxygen](#sources), which the targets
  bind) finds the deep crater lakes' deepest tenth anoxic in mesotrophic and richer water in every target, and in
  oligotrophic water at 0.6–1.5 mg/L in the lock, 0.1–0.7 in the resort at 301 K and under 0.1 in the hot world, with
  productive shallows' floors anoxic most of the month. Anoxic water and sediment make methane, and in the synthesis's
  ledger the wild's lakes and wetlands make most of it (12.8–218 Tg a year, 58.9 central). Methane emissions rise with
  warmth at about 0.96 eV across scales (Yvon-Durocher et al. 2014): 1.9 times from the lock's sea median, 24.1 °C, to
  the resort's 29.5 °C and 3.5 times to the hot world's 34.5 °C. As an indicator, the wild sources the biology leaves,
  5.9 Tg a year, scaled so, with the farms' and herds' methane and the soils' uptake held, would raise the held pair's
  removal from 0.85–3.0 Tg a year to 6.4–8.6 Tg in the resort and 16–18 Tg in the hot world (0.20–0.26% and 0.48–0.54%
  of the air a year through the contactors), or leave the soils alone holding methane near 4.1–4.8 and 7.2–8.0 ppm,
  against the 2.2–2.9 ppm they hold on the unscaled sources. Ecology's tracer through the lake model's columns ([the
  waters' methane](#sources)) sizes the waters themselves at 301 K: with Earth's methanotrophs at their central rates,
  nutrient-poor lakes and seas send 0.57 Tg a year and hold the air at 1.01–1.75 ppm, under the held 1.87, and
  mesotrophic waters send 7.8 Tg and put it at 3.5–4.2 ppm. Lakes kept oligotrophic and oxygenated, by the piped oxygen
  or driven mixing the lake study names, are part of biology's cut.
- **Side effects.** None on health or growth at these levels. The low pair reaches its methane level within a century
  and keeps the contactors running from then; N2O's decomposition heats the air it treats by 400 K with 96% recovered
  (the synthesis's design guess).
- *Evidence: a runnable line-by-line column with Earth's cloud ratio; costs on stated design guesses; the warm waters'
  methane an indicator.*

## CO2

At 290 K, the only ground temperature the column's CO2 sweep covers, the 1.2 atm column sends out 1.99 W/m² more
infrared under a clear sky at 280 ppm than at 400 and absorbs 0.04 W/m² less sunlight, which also takes heat off: 2.04
W/m² in all, with 3.98 W/m² of infrared per doubling ([summary](../../../atmosphere/radiative_convective/README.md)).
Taken all-sky by AR6's effective forcing of doubled CO2 over the Earth control's own change, 3.93 over 4.13 W/m²:

| Level | Partial pressure | Saved, W/m² | Points | C3 leaf rate at 25 / 30 / 35 °C | Water-use efficiency |
|---|---|---|---|---|---|
| 280 ppm | 34.0 Pa | 1.94–2.04 | 0.86–0.90 | 0.71 / 0.67 / 0.62 | 0.70 |
| 350 ppm | 42.6 Pa | 0.74–0.78 | 0.33–0.35 | 0.89 / 0.87 / 0.85 | 0.88 |

The 350 ppm row carries the sunlight of the 280 ppm row by the log of the CO2 ratio, and the targets' grounds lie 5–17
K above the sweep's 290 K. The leaf rates are the [CO2 study](../atmospheric_co2/README.md)'s light-saturated C3 model
against 400 ppm, and the water-use efficiency follows the CO2 at a fixed internal ratio. Water holds most of the land's
growth ([light limitation](../light_limitation/README.md)), so the loss in carbon per litre is the one that counts: 280
ppm takes 30–38% of the growth for its 0.86–0.90 points, 33–44% a point, where a point of dimming costs 0.07–0.30% of
the growth through light and 0.55–1.38% of the land's rain. CO2 is also the upper air's main coolant: at 280 ppm the
warmest titania case's exobase rises 5 K ([middle atmosphere](../../../atmosphere/middle_atmosphere/README.md)). The
managed return that holds 400 ppm, kilns of 3.9–14 GW and an import of 0.12–0.33 Gt C a year at year 10,000
([carbon cycle](../carbon_cycle/README.md)), is small beside any of these costs.
*Evidence: a runnable column at one ground temperature and a textbook leaf model.*

## Surface pressure

The columns compared at the same ground temperature, with the titania stack's sunlight, Earth's humidity profile and a
clear sky: at 1.0 atm the air scatters less and absorbs 5.3–5.4 W/m² more sunlight, and sends out more infrared, 5.3
W/m² more at 290 K, 7.0 at 300 K and 9.1 at 310 K. Net, its balance falls 0.8 W/m² at 295 K, 1.8 at 301 K and 3.1–3.2
at 306.5–307.1 K. That column holds 400 ppm, 17% less CO2 partial pressure than the 1.2 atm column; keeping the plants'
48.6 Pa takes 480 ppm, which gives back 1.04 W/m² under a clear sky.

The clear sky is the upper end. Under cloud, part of the column's change stays hidden: taking the masking Earth's clouds
give the surface's reflection (Lenton & Vaughan's 0.48 over 0.73 at Earth's 67% cover, and carried to the GCM's
20–38%) for the column's whole net change, sunlight and infrared, and the CO2 kept at its all-sky ratio, 0.95:

| With the plants' CO2 kept | The lock, 295 K | The resort, 301 K | The hot world, 306.5–307.1 K |
|---|---|---|---|
| Clear sky, W/m² taken off | −0.25 (warmer) | 0.80 | 2.05–2.19 |
| Clear sky, points | −0.11 | 0.36 | 0.91–0.97 |
| With Earth's cloud masking, points | −0.21 to −0.12 | 0.10–0.26 | 0.46–0.71 |

- **It works only where the Moon is hot**, because the hot, wet column's infrared opacity is what the lower pressure
  thins; in the lock it warms the Moon under a clear sky and more under cloud.
- **Flight.** Buoyant lift falls 17% at every height, and the flight band's air lies 8–9 km lower. The register keeps
  1.2 atm for aerial craft.
- **Fire and breath.** Oxygen keeps Earth's partial pressure, so breathing is unchanged, and its share rises from 17.5%
  to 21%, so fires behave as on Earth; the port fire study's design fires would change.
- **Resources.** The air would need 20% less N2, 5.1×10¹⁷ kg of the import that takes 10²⁰–10²¹ kg of small bodies.
- *Evidence: clear-sky columns at prescribed humidity, with Earth's cloud masking carried over; no GCM at 1.0 atm.*

## The high cloud

In the hot world the cloud reflects 10.8 W/m² of sunlight and holds back 16.8 W/m² of infrared, 6.0 W/m² of net
warming, against 0.3 W/m² in the design case
([climate/gcm](../../../climate/gcm/README.md#where-the-corrected-moon-settles)). The refit's high cloud therefore
carries about 5.6 W/m² of net warming; the check year gave its ice cloud 5.7.

On Earth cirrus warm the planet by about 5 W/m² on net, and AR6 puts the most that thinning by seeding could take at
1–2 W/m², with low confidence. The model estimates span a change of sign:

| Study | Approach | Net effect, W/m² |
|---|---|---|
| Mitchell & Finnegan 2009 | First estimate, seeding by airliners | more negative than −2.8 |
| Storelvmo et al. 2013 | CAM, seeding concentration varied | down to −2.0 at an optimum; warming beyond it |
| Muri et al. 2014 | Ice falling twice as fast | −1.55, cooling 0.94 K |
| Gasparini et al. 2020 | Favourable seeding in ECHAM and CESM | −0.8 and −1.8 |
| Liu & Shi 2021 | CAM5, seeding by local conditions | −1.36; −2.0 kept out of low latitudes, where seeding warmed |
| Penner et al. 2015 | CAM5 | −0.74 only with weak updrafts and no growth on existing ice |
| Gasparini & Lohmann 2016 | ECHAM-HAM | no significant cooling |
| Tully et al. 2022 | ECHAM-HAM with P3 | no significant cooling (−0.36 at one concentration); up to +4.9 from overseeding |

- **Carried over** as Earth's share of the cirrus effect, up to two fifths, thinning the Moon's high cloud could take
  0–2.3 W/m², 0–1.0 points, and overseeding could warm it.
- **How a sky-ship world could do it.** The cloud sits in the GCM's top layers, 63–107 km up at 9–32 kPa, above the
  flight band and about the height of the high platforms near 70 km; at 0.16 g buoyant platforms float there readily.
  Mineral dust starts ice at 110–130% of ice saturation against clean air's 154–167% (climate/gcm), and the
  highlands' plagioclase is the Moon's own dust. The long night lets seeding work on the night side, where the
  cloud traps infrared and reflects no sunlight, as Gasparini et al. seeded at night only in ECHAM. Seeding 1–20
  particles a litre of 0.5 µm feldspar through the layer, renewed every ten days, takes 12–240 kt a year (*estimate*),
  against the film plant's 2.3–4.6 Gt.
- **Rain.** Thinning works through the infrared and strengthens the water cycle, where dimming weakens it
  (Kristjánsson et al. 2015).
- **The GCM's high cloud is itself the open check.** PlaSim makes cloud from relative humidity at Earth's thresholds, in
  coarse top layers, with crystal sizes carried past their tables; it has no ice nucleation, so it cannot represent
  seeding. The plan's item 3 gives the high cloud a range from its water and crystal size, offline in minutes. That
  range may move the hot world itself by more than seeding could, and comes first.
- *Evidence: a concept on Earth's literature, beyond what the project's models compute.*

## The heat the plan adds

### Beamed power

The beams put 33–157 TW of heat on the Moon, 0.87–4.14 W/m², 0.39–1.84 points. By the domain synthesis's array case
the loads are the residents' 11–33 TW, food 10–78, the sky towns' hydrogen remade 2.4–25, sky control 0.4–3.8 and close
computing 6.4 TW. Food is the part that need not be made on the Moon: made at the ring radius, its 11–84 TW with the
air's share would put 0.0005–0.004 W/m² on the Moon, and so take 0.29–2.22 W/m² off it, 0.13–0.98 points.

- **Its carbon arrives as food.** The beams' food power makes the blend's power-made and fully lamp-lit mixes whole and
  power-heavy's lamp-lit vegetables and fruit and microbial protein, 65% of the diet's energy
  ([food systems](../food_systems/README.md#the-authors-blend)). Made in orbit, that food brings 0.47 Gt C a year to
  the Moon as 8 billion people eat it (*estimate*, 0.25 kg C a person a day), and up to 0.71 if all of the food study's
  losses, 1.3–1.5 times what is eaten, fall after it arrives: 1.5–6.2 times the managed import of 0.12–0.33 Gt C a year
  that holds 400 ppm at year 10,000 ([carbon cycle](../carbon_cycle/README.md)). The food would replace that import and
  exceed it, so its carbon cycle closes only if the rest goes back to orbit or into burial on the Moon, the kilns
  resized with it; its phosphorus and nitrogen would arrive as well, which changes the phosphorus cycle's return path.
- **Costs.** Food plants in the array's industry and their mass; the food's fall heat is negligible; supply through
  the beams' reserves becomes food security. It departs from the chosen design in where the beamed food is made; the
  blend stays.
- **The rest stays.** Loads that must be on the Moon add heat whatever powers them, except the Moon's own sunlight,
  wind and rivers, which add none but cost land and storage through the 354-hour night.

### The water top-up's fall

Through the first thousand years the top-up adds 1.5–9.5 W/m², 0.68–4.20 points
([climate targets](../climate_targets/README.md#the-heat-beyond-the-sunlight)). From the exobase down, 1.22–1.72 MJ/kg
as the exobase lies 1.77–2.56 lunar radii out (resources' first screen), it is 1.5–5.9 W/m²; above it, at the high end,
3.5–5.3 W/m², 1.6–2.3 points.

- **Braking above the air meets the inflow's momentum.** Bringing 4.8–13×10⁷ kg/s to rest at the exobase from the
  1.43–1.74 km/s of its fall from the far terminal takes a steady force of 6.8–23×10¹⁰ N, which photons could carry
  only on a reflecting sail 790–2,700 times the Moon's disc.
  A tether in orbit needs as much mass going outward, or an anchored structure that hands the momentum to the Moon;
  elevators cross the ring stack twice a month ([resources](../../../resources/README.md), finding 11), and today's
  lunar elevator carries 584 t a year, 2.6–7.1×10⁹ times short. The energy shed above the exobase would be 49–199 TW.
- **Recovered and used on the Moon, it displaces the beams.** The fall carries 58–359 TW, more than the Moon's whole
  beamed intake of 33–157 TW, so recovery that ran the Moon's loads would spare the beams' heat, 0.87–4.14 W/m², the
  two heats no longer adding; no concept recovers this flow.
- **What stays practical is timing.** The top-up is a heat of the first millennia, 1.8–4.0 W/m² on to year 10,000,
  answered by temporary dimming whose light the array collects; in the fast case a slower top-up spreads the 12–28
  W/m² of its busiest centuries.
- *Evidence: solid energy and momentum accounting; no recovery hardware exists at this rate.*

### The glow

The GCM gives the glow's infrared an efficacy of 0.645 from the glow pair on the corrected radiation of 17 October and
0.71 from the held branch of 13 October. Rolling both halves 80° took 5.3–5.8 W/m² of absorbed-sunlight equivalent off
the as-built fleet's 9.2–10.0 (2.3–2.6 points) and leaves 3.89–4.24 W/m², 1.72–1.88 points
([glow](../solar_shield_array/glow.md), read only); the lower efficacy takes 0.16 points off the glow's dimming. The two
answers the author set aside on 13 October, with the plants' water design making up the glow's rain, would take more:
the annulus absorber 1.40–1.53 W/m² (0.62–0.68 points), 1.1 g/m² more on the annulus, which needs the atmosphere
domain's check of 250–400 nm above the limb and gives up two-thirds of the annulus film's light pressure, restored by
about 5.5% more reflection; and the outer faces' ultraviolet mirror 1.70–1.87 W/m² (0.75–0.83 points) at about 60%
more fleet mass. A glow pair at the chosen target measures the efficacy again.

## The Moon's own reflection

The deep air halves what the surface can do. In the GCM radiation check's clear-sky columns at the equator the Moon
reflects 0.34 more of the arriving sunlight for each unit of surface albedo, against Earth's 0.73 (Lenton & Vaughan
2009), and 0.28–0.31 under the GCM's cloud. The kernel falls toward the poles, where the Sun's slanting path crosses
more air: against the sea's reflection at the same latitude, the land columns give 0.31–0.32 at 36° and 0.28–0.29 at
58°. The deserts take the kernel at each row's latitude, poleward of 58° held there at the high end and carried on
along its fall at the low; the equator's kernel would read the hot world's deserts 15–24% high and the lock's 12–18%.

| Surface | Area | Brighter by | W/m² | Points |
|---|---|---|---|---|
| Roofs and pavements of the settled land | 0.26–0.74 million km², 0.7–1.9% of the Moon | 0.10–0.15 | 0.04–0.28 | 0.02–0.12 |
| Silver foliage in the hybrid parks | 4.68 million km², 12.3% | 0.02–0.05 | 0.13–0.60 | 0.06–0.26 |
| The hot world's deserts, mostly poleward of 45° | 11.1–11.7% of the Moon, 5.6–6.0% weighted by sunlight | 0.10–0.30 | 0.39–1.52 | 0.17–0.67 |
| The lock's deserts | 23.6–24.3%, 16.2–16.9% weighted by sunlight | 0.10–0.30 | 1.19–4.36 | 0.53–1.94 |

The deserts' areas span the extract's two counts on the model's land mask, every cell's land with its own rain and the
model's own land cells.

- **Roofs and foliage come with designs already adopted.** The dry gardens' silver foliage is in the author's parks
  design, and white roofs cool cities where people live: on Earth they cut the urban heat island by a third and the
  daily maximum by 0.6 °C (Oleson et al. 2010), a comfort lever for the hot world's cities.
- **The deserts cost rain and wild land.** On Earth brightening deserts takes 3.9% of the land's rain for each kelvin
  of cooling, against a sunshade's 1.0%, and weakens the monsoons (Irvine et al. 2011); here it would cover the wild
  deserts the author keeps.
- *Evidence: kernels by latitude from line-by-line clear-sky shares and Earth's cloud masking scaled by the GCM's cover;
  areas from the domain synthesis and the GCM's land mask and extract; brightening by design guess.*

## Against the shield's own cost

A point of dimming carries 111 TW of light, which collectors turn into 22–33 TW, so the dimming costs no energy. Its
costs lie in the fleet and the biosphere ([glow](../solar_shield_array/glow.md),
[constraints review](../solar_shield_array/constraints_review.md), read only):

- a dimmer absorbing in the window films returns 10% of its dimming as glow with both halves rolled, 25% as built;
- one reflecting from the window tiles adds about 0.022 a point to their light-pressure coefficient, 0.276 to about
  0.32 for two points, and under one sail loading each point takes from the annulus tiles' mass allowance, which no
  tested tile already meets;
- each point takes 0.07–0.30% of the living Moon's growth through light and 0.55–1.38% of the land's rain.

A lever is worth pursuing where it saves points at less than that: without taking rain or light, or without loading
the window tiles.

## What is worth pursuing

| Rank | Lever | Points saved | What it costs | Side effects | Evidence |
|---|---|---|---|---|---|
| 1 | The glow's set-aside answers: annulus absorber or ultraviolet mirror | 0.62–0.83 | A film change and a check above the limb, or 60% more fleet mass | Less of the glow's rain to make up | Screen on the fleet's products |
| 2 | Food made in orbit with the beams' power | 0.13–0.98 | Food plants in the array's industry | Its carbon, 1.5–6.2 times the import, needs a return path; nutrients from orbit | Arithmetic on the synthesis's ledgers |
| 3 | Methane and N2O at the low pair | 0.79–0.81 | Contactors through 1.0–1.4% of the air a year; 1.0–1.6 TW | None on health or growth; warm lakes push the other way | Column forcing; design-guess costs |
| 4 | White roofs and silver foliage | 0.07–0.39 | Little beyond adopted designs | Cooler cities | Kernel arithmetic |
| 5 | Night-side cirrus thinning | 0 to 1.0, sign open | 12–240 kt of feldspar a year from high platforms | Adds to the rain, which dimming takes | Earth's models; first the GCM's high-cloud range |
| 6 | The top-up's fall | 0 to 1.6–2.3, first millennia only | Momentum binds any braking above the air | Temporary dimming meanwhile | Energy and momentum accounting |
| 7 | 1.0 atm | 0.91–0.97 in the hot world under a clear sky, 0.46–0.71 with Earth's cloud masking; −0.11 in the lock | Reopens 1.2 atm; 17% less lift | 20% less N2 to import; Earth's fire | Clear-sky columns |
| 8 | Brighter deserts | 0.17–1.94 | Covering wild land | Takes the land's rain four times as fast as dimming | Kernel arithmetic, Earth's models |
| 9 | CO2 below 400 ppm | 0.33–0.90 | 33–44% of the growth a point at 280 ppm | Warmer upper air | Column and leaf model |

The first four together save 1.6–3.0 points: 30–81% of the hot world's dimming in all, 15–37% of the resort's and
10–24% of the lock's. The levers can answer much of the plan's added heat; bringing the Moon down to the resort or the
lock stays mostly dimming. The gas pairs carry their own overlap, and the others act on separate terms.

Questions for the author: whether the glow's set-aside answers come back if the resort or the lock is chosen; whether
food for the beams may be made in orbit, with a return path for the carbon it brings beyond the import; and whether the
GCM's high-cloud range runs before the target is chosen, since it moves every target's dimming.

## Assumptions and limits

- **Efficacy 1** for every infrared and ground-heat lever and the GCM's 0.645–0.71 for the glow, as the targets take
  them.
- **Columns** are clear-sky with prescribed humidity at 290–310 K; the gases' all-sky values carry Earth's ratios,
  CO2's AR6's, and the 1.0 atm column's change Earth's cloud masking of the surface's reflection. The CO2 sweep is at
  290 K only.
- **The methane indicator** scales the wild sources the biology leaves by the seasonal temperature dependence of
  emissions between sea medians, with the farms' and herds' methane and the soils' uptake held; substrate and oxygen
  also set the lakes' methane, and ecology's tracer gives the waters' own.
- **Cirrus thinning** takes Earth's share of the cirrus effect; the seeding mass assumes 0.5 µm feldspar particles
  through 63–107 km, renewed every ten days.
- **The albedo kernel** comes from the radiation check's land columns below 4 km against the sea's reflection at their
  latitudes, with Earth's cloud masking (0.48 over 0.73 at 67% cover) scaled by the GCM's 20–38% cover. No land column
  lies poleward of 58°. The cities and parks take the equator's kernel, which reads high where they lie poleward of 45°,
  and 0.6–1.0 of the Moon's mean sunlight. Over the lock's deserts the GCM's cover is 0.11, about half its Moon-wide
  0.20, which would raise their figures about 5%. The deserts by band come from the climate targets' extract (on the
  research drive, ignored) and the GCM's land mask.
- **Food's carbon** takes 0.25 kg C a person a day, the blend's powered share of the diet's energy, and at the high
  end the food study's losses after arrival.
- **Not done**: no lever runs in a climate model; the high cloud's range (climate) is other work.

## Running it

```sh
python -m research.studies.heat_levers.run          # under a second, after climate_targets' extract
python -m pytest research/studies/heat_levers
```

## Sources

Published values are in [sources.json](sources.json), with how far each was read. Ecology's waters' oxygen,
`research/studies/lake_oxygen` on `domain/ecology` at bd0097c, comes through the climate targets, which bind it, and the
product carries the deep lakes' figures. Ecology's waters' methane, `research/studies/lake_methane` on `domain/ecology`
at 4e4f710, was read with `git show`; the README quotes its numbers, and the product binds none of them.
