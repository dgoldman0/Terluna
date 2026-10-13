# Where light limits growth

13 October 2026. The author asked where and when light is the factor that limits growth on the Open Moon, if anywhere:
it bears on the gas pair, whose warming the dimming answers, and it is worth knowing in general. Light limits growth
where nothing else does first, so this screen takes each part of the living Moon in turn, reads from the project's
products what holds its growth there, and weighs how much of its growth follows the light.

`python -m research.studies.light_limitation.run` writes [results/light_limitation.json](results/light_limitation.json)
(schema `terluna.research.light-limitation/1`, with its producer, the sha256 of its seven inputs and of the GCM
climatology, its evidence statement and reading rule) in under a second. It reads the canopy's
[light response](../../../biosphere/canopy/README.md#the-light-response), computed for it on the plant model:
`python -m biosphere.canopy.light_response`, about a minute and a half on one thread. [sources.json](sources.json)
records each published value and how far it was read.

## What it finds

- **Water holds the land.** The rain supports 55–57% of the canopy's light-limited growth within 15° of the equator,
  31–32% at 15–30°, 15% at 30–45°, 2% at 45–60° and under 1% beyond (the food systems study's band basis, from
  ecology's crop model and the design climate). Rain-fed land grows 14.0–14.6 Gt C a year; with the wet margins the
  land grows 32–37% of the 47–55 Gt C its light would grow with water everywhere. The GCM's own soil water agrees:
  the design run's evaporation falls below its potential on 38% of the land within 15° of the equator, 76% at 15–30°,
  90% at 30–45° and all of it beyond. Temperature holds growth nowhere; the coldest three-day mean anywhere is 14 °C.
  On Earth, by comparison, water limits growth on about 40% of the vegetated land, temperature on 33% and radiation on
  27% (Nemani et al. 2003).
- **Light limits where water is ample.** On the wet margins of the lakes and rivers, 0.5–2 million km² in ecology's
  scenario, growth runs at the light's potential, 1.2–5.4 Gt C a year. Watered land would grow the same way: within
  30° of the equator the light binds first, and poleward of about 40° the sunlight's energy for transpiration binds
  first, so watered growth there follows the light one for one. The chosen food blend waters no fields and lights its
  lamp farms with power: its 0.68–1.86 million km² of fields are rain-fed, and its other food grows without sunlight.
- **In the waters,** light limits the lakes, which keep 96% of the phosphorus the land sends them, and the river-fed
  coasts. Their growth follows the light at 0.38–0.59 near the equator at fixed biomass by Talling's column (1957),
  more poleward, and up to one for one where the biomass grows with the light. The open seas are nutrient-poor by
  design, at 0.1 mg/m³ of chlorophyll; the production model's light term (Behrenfeld & Falkowski 1997) puts their
  response at 0.08–0.14.
- **In the sky,** the aerophytes' growth follows their light in the sky ecology's models, between a sparse canopy's
  photosynthesis, 0.30, and a deep canopy's growth, 0.67.
- **Through the lunar day,** light matters most at dawn, dusk and in the long twilight. The deep canopy loses 0.97% of
  its twilight photosynthesis for each 1% of light taken, 0.86% with the Sun below 10°, 0.62% at 10–30°, 0.45% at
  30–60° and 0.39% above 60°, where its top leaves are near saturation. The high Sun gives 44% of the cycle's
  photosynthesis and 36% of what a dimming takes. With water and nutrients ample, each 1% of light costs 0.46–0.65% of
  gross photosynthesis, more poleward and in deeper canopies, and 0.56–1.1% of growth, because upkeep does not shrink
  with the light. The night's store grows by 0.08–0.28% for each 1%.
- **What the dimming costs through light.** Light limits 10–40% of the living Moon's 16–26 Gt C a year of growth,
  and each 1% of dimming costs 0.07–0.30% of it. The dimming the plan now carries, 4.6–6.4%, costs 0.3–1.9%: the gases
  at the held pair take 1.69%, people's heat 0.26–2.0% with power from the array and close computing, and the fleet's
  glow 2.64% with both halves of each orbit rolled outside the service arc (research/array-industry, 8be6d1e). The
  GCM finds the glow's infrared needs 0.71 of its watts in dimming to hold the design temperature.
- **The rain may cost more than the light.** In Earth's models, offsetting quadrupled CO2 with 3.5–5.0% less sunlight
  held the temperature within 0.3 K and cut global rain by 4.5 ± 1.3% (GeoMIP G1; Tilmes et al. 2013), 0.76–1.38% for
  each 1% of dimming: the rain follows the energy the surface receives. If the Moon's rain answered the dimming for the
  gases and the glow the same way, the rain-fed land, 57–87% of the growth, would lose 0.43–1.2% of the living Moon's
  growth for each 1% of that dimming, 0.7–5.2% at the plan's 1.7–4.3%. People's heat is released at the ground, where
  the dimming that answers it takes about as much energy from the surface as the heat adds, so it is left out. The
  GCM run with these forcings and the dimming gives the Moon's own response
  ([gases, glow and dimming](../gases_glow_dimming/README.md)): the gases' dimming takes 0.55% of the rain for each 1%,
  about half Earth's, and the glow, whose infrared heats the upper air, takes 2.4% more with the temperature held.
- **The gas pair holds.** Through light, the middle pair's 1.69% costs 0.12–0.50% of the growth; the low pair would
  save 0.06–0.25% of it and Earth's present air would add 0.04–0.16%. Through the rain, at Earth's response, the low
  pair would save 0.4–1.0% and Earth's air would cost 0.2–0.6% more. Set against what the pairs move elsewhere, the
  low pair's contactors through 1.0–1.4% of the air from the first century and Earth's air's N2O plants a third the
  middle pair's, these leave the middle pair the balanced choice. The GCM finds the gases' rain response about half
  Earth's, which keeps it so.

## The land, band by band

| Band | Dry land, million km² | Rain supports | GCM land below its soil-water stress | Rain-fed growth, Gt C a year | With water ample, Gt C a year | Watered land held by |
|---|---|---|---|---|---|---|
| 0–15° | 5.44 | 55–57% | 38% (4.8 mm of rain a day) | 7.8–7.9 | 13.2–15.0 | light |
| 15–30° | 5.58 | 31–32% | 76% (2.7 mm) | 4.3–4.4 | 13.0–14.8 | light |
| 30–45° | 5.05 | 15% | 90% (1.3 mm) | 1.8–2.1 | 10.5–12.3 | light or its energy |
| 45–60° | 3.69 | 2.2–2.3% | all (0.19 mm) | 0.17–0.20 | 5.5–6.7 | the sunlight's energy |
| 60–75° | 2.71 | 0.6% | all (0.05 mm) | 0.03–0.04 | 3.5–4.3 | the sunlight's energy |
| 75–90° | 0.92 | 0.5% | all (0.04 mm) | 0.01 | 1.2–1.5 | the sunlight's energy |

The rain's share is the food systems study's, raised to the water need the design sunlight allows. The potential
takes the canopy model's clear-sky growth at the equator, 2.77–3.16 kg C m⁻² a year between the Earth-like plant and
one idling at a quarter of its upkeep through the night, times each band's light under the design clouds (0.89) or
the sunlight's cap on transpiration, whichever is lower. The GCM's stress is ExoPlaSim's own: its land evaporates at
its potential only while the soil holds 0.2 m of water, 0.4 of its 0.5 m capacity. The wet margins are placed in the
two wettest bands in proportion to their land and grow at the potential.

## Production by place

| Place | Growth, Gt C a year | Share of the living Moon's | Held by | Light elasticity | Basis |
|---|---|---|---|---|---|
| Rain-fed land | 14.0–14.6 | 57–87% | water | 0 | food systems band basis |
| Wet margins | 1.2–5.4 | 7–21% | light | 0.60–0.71 | ecology's scenario at the canopy's potential |
| Lakes | 0.39–1.57 | 2–6% | light | 0.38–1.0 | 3.92 million km² at 100–400 g C m⁻² (*design guess*) |
| River-fed coasts | 0.16–1.59 | 1–6% | light | 0.38–1.0 | 5–15% of the seas at 300–1,000 g C m⁻² (*design guess*) |
| Open seas | 0.45–1.51 | 3–6% | nutrients | 0.08–0.14 | 50–150 g C m⁻² (*design guess* on Field et al. 1998) |
| Sky | 0.02–0.83 | 0.1–3% | light | 0.30–0.67 | sky ecology at its reference light, 0.1–2% cover |
| All | 16.2–25.5 | | | 0.07–0.30 weighted | |

The elasticity is the share of growth lost for each share of light taken. Rain-fed land's is zero to first order:
its growth is the water it transpires times the plants' water-use efficiency, which the light does not set. The
weighted figure's ends pair the light-limited places' low production and response with the rest's high production,
and the reverse.

## When light limits

| Sun height | 90° | 60° | 30° | 10° | 0° | −10° | −20° |
|---|---|---|---|---|---|---|---|
| Deep canopy's photosynthesis lost per 1% of light | 0.35% | 0.38% | 0.52% | 0.79% | 0.94% | 0.99% | 1.00% |

The row is the canopy model's leaf area index 5 at 22 °C. Sparse canopies follow the light less: over the cycle at
the equator a leaf area index of 1 loses 0.30% of its photosynthesis for each 1% and one of 5 loses 0.47%; at 60° they
lose 0.49% and 0.64%. Growth follows the light more than photosynthesis does, since the upkeep of leaves, stems and
roots stays, and most for the plant that respires through the night as by day.

## What it leaves out

- **The rain's response on the Moon.** This screen takes Earth's G1 as its reference. Its land evaporation fell about
  10%, partly through CO2 closing the stomata, which the Moon's case lacks at 400 ppm; its fast response comes from
  CO2's absorption in the air, as the gases' does here. The GCM gives the Moon's own
  ([gases, glow and dimming](../gases_glow_dimming/README.md)): about half Earth's for the gases, and for the glow,
  whose infrared heats the upper air with about 4% reaching the ground, 0.46% of the rain for each W/m² of it with the
  temperature held.
- **Partial water limits.** Rain-fed land's light response is taken as zero. After rain, before the soil dries, plants
  follow the light for a while, so its true response is small and positive.
- **The waters' and the sky's production** are design guesses on Earth's waters and the sky ecology's ranges; the
  canopy model has no water, nutrients or acclimation to the long day.
- **Nutrients on land.** The engineered biosphere's nitrogen and phosphorus are taken as ample where the water is; the
  wet Moon's nitrogen cycle is ecology's next work.
- **The far side** differs from the near by 2–3% in the plant model and is not treated apart.

## Files

| File | Holds |
|---|---|
| [run.py](run.py) | The screen: the land by band, the wet margins, the waters, the sky, the cycle, the ledger, the dimming's cost through light and the rain's response from Earth's models |
| [test_run.py](test_run.py) | Its binding, regeneration and arithmetic |
| [sources.json](sources.json) | The published values and how far each was read |
| [results/light_limitation.json](results/light_limitation.json) | The product |

```
python -m biosphere.canopy.light_response          # the canopy with the light dimmed, about a minute and a half
python -m research.studies.light_limitation.run     # the screen
python -m pytest research/studies/light_limitation biosphere/tests/test_light_response.py
```
