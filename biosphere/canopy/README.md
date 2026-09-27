# Canopy photosynthesis

How much carbon a plant canopy can fix under the Open Moon's sky, set against the
same canopy under Earth's: the light at the ground from the illumination domain,
followed into the canopy by wavelength and into each layer's sunlit and shaded
leaves.

`python -m biosphere.canopy.fetch_inputs --download` restores the PROSPECT-D leaf
coefficients (pinned in [inputs.json](inputs.json), ignored by Git), and
`python -m biosphere.canopy.model` (about a minute and a half on one thread, with
`OPENBLAS_NUM_THREADS=1`) writes [results/canopy.json](results/canopy.json)
(schema `terluna.biosphere.canopy-photosynthesis/1`).

## Method

| File | Computes |
|---|---|
| [leaf_optics.py](leaf_optics.py) | Leaf reflectance and transmittance, 400–2500 nm, from PROSPECT-D (Féret et al. 2017): a pile of absorbing plates (Jacquemoud and Baret 1990) whose absorption comes from the leaf's chlorophyll, carotenoids, anthocyanins, water and dry matter |
| [radiation.py](radiation.py) | Light inside the canopy by wavelength: the direct beam, diffuse light in eight streams each way, scattering by bi-Lambertian leaves with a spherical angle distribution (backscatter after Sellers 1985, as in CLM5), a Lambertian soil, and the sunlit share of each layer |
| [leaf.py](leaf.py) | C3 leaf photosynthesis (Farquhar, von Caemmerer and Berry 1980) with Bernacchi et al. (2001) kinetics and the capacities of the [CO₂ study](../../research/studies/atmospheric_co2/README.md), whose light-saturated leaf it reproduces exactly; electron transport follows the absorbed photons |
| [model.py](model.py) | Reads the [surface-light product](../../illumination/surface_light/README.md), adds the sky's return of what the canopy reflects, and integrates over Sun heights, the solar cycle and latitude; leaf carbon over the cycle and the store it needs through the night (the reserve arithmetic of [long_night.py](../long_night.py)) |
| [plant.py](plant.py) | The whole plant through the lunar cycle, with twilight and the chosen climate's temperatures, for ways of getting through the night (below) |
| [fruit.py](fruit.py) | Fruit designed for the lunar day on that plant: how fast a large fruit can grow, when to set it, and how much of the plant's growth it can take (below) |

The base stand has a leaf area index of 5 in layers of 0.025, randomly placed leaves
with a spherical angle distribution, soil reflecting 10%, and PROSPECT-D green leaves
(N = 1.5, 40 µg/cm² of chlorophyll, 8 of carotenoids, 0.01 cm of water, 0.009 g/cm² of
dry matter). Rubisco capacity is 60 µmol m⁻² s⁻¹ at the top and falls as exp(−0.157 L)
with depth (Lloyd et al. 2010). A sunlit leaf's share of the direct beam is spread over
its orientations (uniform in the cosine for spherical leaves). Photons from 400 to
700 nm drive photosynthesis. Both worlds hold O₂ at Earth's partial pressure and
400 ppm of CO₂: 48.6 Pa on the Moon at 1.2 atm and 40.5 Pa on Earth. Leaves are at
22 °C by day and respire at 19 °C by night, the equatorial land air of the chosen
climate at its warmest and coldest; Earth gets the same temperatures, so the two
skies and the two CO₂ pressures are what differ.

The PROSPECT-D code matches the prosail package's implementation to 10⁻¹⁵. The canopy
radiation reproduces the exact direct and diffuse transmission of black leaves,
conserves energy to 10⁻⁹, counts sunlit leaf area exactly, and changes by under 0.2%
of the incident light when its layers are made four times thinner.

## Results

Clear sky, the equator, leaf area index 5. Rates are µmol CO₂ per m² of ground per
second.

| Sun height | Moon: light in / absorbed / fixed | Earth: light in / absorbed / fixed | Moon ÷ Earth, fixed |
|---|---|---|---|
| 90° | 1,523 / 1,373 / 33.0 | 2,238 / 1,994 / 28.4 | 1.16 |
| 30° | 614 / 574 / 22.4 | 1,037 / 981 / 18.3 | 1.22 |
| 10° | 162 / 150 / 9.7 | 287 / 269 / 8.8 | 1.09 |

Over the whole cycle the Moon's canopy absorbs 35% fewer photosynthetic photons than
Earth's and fixes 18% more carbon: 12.0 against 10.2 µmol m⁻² s⁻¹, or 12.5 against
10.6 g C per m² per 24 hours. Per photon absorbed it fixes 1.8 times as much.

Changing Earth's light into the Moon's one factor at a time, averaged over the cycle
at the equator:

| Step | Canopy photosynthesis | Change |
|---|---|---|
| Earth's light and CO₂ | 10.19 | |
| The Moon's CO₂ pressure (48.6 Pa) | 10.83 | +6% |
| Earth's spectrum and diffuse share, scaled to the Moon's photons | 9.06 | −16% |
| The Moon's spectrum, with Earth's diffuse share | 9.12 | +0.6% |
| The Moon's diffuse share: its own light | 12.03 | +32% |

- **The diffuse sky outweighs the missing photons.** Under Earth's direct beam the
  sunlit leaves at the top are saturated and the shaded ones starve; the Moon's
  sky, 35% diffuse overhead and 80% at 10°, spreads the light over the leaves that
  can use it. On light alone, at the same CO₂, the Moon's canopy fixes 11% more than
  Earth's at the equator and 12% more at 30° and 60°, where the diffuse gain is
  larger still (+36% and +47%).
- **Colour barely matters.** The Moon's redder light changes photosynthesis by 0.6%.
- **Deeper canopies gain most.** The Moon's advantage grows with leaf area: 9% at a
  leaf area index of 1, 14% at 2, 18% at 5 to 8.
- **Changes that help on Earth help less on the Moon.** Counting 700–750 nm photons
  adds 6% on the Moon and 11% on Earth, because the Moon's far-red is short. Paler
  leaves in the top third of the canopy (20 µg/cm² of chlorophyll) add 3% on the
  Moon and 8% on Earth, since the diffuse sky already spreads the light. A clumped
  canopy (index 0.7) adds 3% on the Moon and loses 2% on Earth. Leaf capacity moves
  both by about ±20% between 40 and 90 µmol m⁻² s⁻¹ at the top.
- **The night is the cost.** Over a lunar cycle at the equator the leaves fix
  369 g C/m², respire 40 by day and 33 by night, and keep 297 net (10.0 per
  24 hours, against Earth's 8.1). To last the 365 hours from dusk to the next
  morning's positive balance they need a store of 33 g C per m² of ground, thirty
  times Earth's 1.1 for one night, and 44 at a leaf area index of 8. That is about
  7 g of carbon, 17 g of carbohydrate, per m² of leaf, a large share of a leaf's
  own mass, so the store belongs in stems and roots. Leaves alone set the optimum
  leaf area index near 6 on both worlds.
- **Latitude.** At 30° the Moon's canopy fixes 11.5 g C/m² per 24 hours (Earth at
  equinox 9.7), and at 60° 7.9 (6.7).

## Limits

- Clear sky only, with the surface-light product's light, which understates the
  light below about 20° of Sun height; the Moon's totals are low by a few per cent
  for that reason (plant.py corrects it with the sky atlas).
- The leaf model is steady-state, with internal CO₂ fixed at 0.7 of ambient: no
  stomata, water supply or leaf energy balance, and no cost for CO₂'s slower
  diffusion in the denser air (the CO₂ study puts it at 5–7% at light saturation,
  which takes back most of the 6% gained from the higher partial pressure).
- No acclimation to 354 hours of continuous light, no slowing of photosynthesis as
  sugar stores fill, no photoinhibition. Whether leaves use two weeks of light as
  fully as the model assumes is the first biological question it leaves open.
- Leaf respiration only. Stems, roots and growth are not counted; with the leaves
  they respire roughly half of a forest's gross photosynthesis on Earth.
- One leaf type, a spherical angle distribution, random or uniformly clumped
  placement, and a neutral soil. The pale-leaf and far-red variants are scenarios.

## Through the night: the whole plant

`python -m biosphere.canopy.plant` (under a minute) writes
[results/plant.json](results/plant.json) (schema
`terluna.biosphere.plant-carbon-cycle/1`). Besides the canopy's inputs it reads two
generated products: the chosen climate's climatology with land air temperature by
local time (`climate/gcm/.venv/bin/python -m climate.gcm.climatology A28_dim5:15-24`)
and the sky atlases ([illumination/sky](../../illumination/sky/README.md)).

It follows the base stand through one solar cycle at the equator, 30° and 60°, on the
near and far sides, and treats the parts of the day the canopy model leaves out:

- **Temperature.** Leaves are at the land air temperature of the chosen climate
  through the lunar day (GCM run A28_dim5, 3-day means composited by hour angle).
- **Low Sun and twilight.** Below 30° of Sun height the diffuse light is scaled to
  the spherical sky atlas, and after sunset the canopy gets the atlas's twilight,
  cut by the shield's 6%, with the spectrum of the diffuse light at 1°.
- **The whole plant.** Stems and roots add maintenance respiration (Q10 = 2), set so
  that the same stand under Earth's sky, with the same temperatures over 24 hours,
  turns half of its gross photosynthesis into new tissue (carbon-use efficiency 0.5;
  Waring et al. 1998). New tissue costs a quarter of its carbon again. The store is
  the largest cumulative deficit over the cycle ([long_night.py](../long_night.py)).

The strategies are:

- **Earth-like:** respire in the dark as by day.
- **Idle:** slow the upkeep of leaves, stems and roots to a half, a quarter or a tenth
  whenever photosynthesis no longer covers it.
- **Regrow:** drop the canopy when dusk photosynthesis no longer covers its upkeep,
  and grow it again from sunrise over 5 days (90 g of dry leaf per m² of leaf, 47%
  carbon, plus growth respiration).

### Results

The chosen climate's land nights stay at 19–20 °C at every latitude and on both
sides, with days of 21–23.5 °C. Cold does not slow respiration anywhere.

The night is shorter than the Sun's geometry says. The Moon's tall sky stays sunlit
far above the ground after sunset. At the equator, photosynthetic light stays above
50 µmol m⁻² s⁻¹ for 10 hours after sunset, above 10 for 32 hours and above 1 for
58 hours, until the Sun is 29° down. Only 238 of the 354 hours of night are darker
than that; on Earth light
falls below 1 µmol m⁻² s⁻¹ 16 minutes after sunset. At 30° latitude 217 hours are
dark. At 60° the Sun never sinks more than 30° below the horizon, and only 34 hours
are.

Equatorial near side, per m² of ground; a cycle is the lunar month, and Earth's stand
has the same temperatures over 24 hours:

| Strategy | Growth per cycle | Per 24 hours | Carbon-use efficiency | Store for the night | Deficit lasts |
|---|---|---|---|---|---|
| Earth's stand, one day | 5.4 g C | 5.4 g C | 0.50 | 1.8 g C | 12 h |
| Earth-like | 229 g C | 7.8 g C | 0.56 | 46 g C | 336 h |
| Idle at half | 249 g C | 8.4 g C | 0.61 | 21 g C | 314 h |
| Idle at a quarter | 259 g C | 8.8 g C | 0.64 | 10 g C | 295 h |
| Idle at a tenth | 265 g C | 9.0 g C | 0.65 | 3.7 g C | 273 h |
| Regrow leaves | 28 g C | 0.9 g C | 0.07 | 206 g C | 467 h |

- **Every evergreen strategy closes its budget.** Even the Earth-like plant grows
  faster than the same stand on Earth. Its canopy fixes more (the diffuse sky, the
  CO₂ pressure and the twilight) while it respires the same per hour. Its store is
  25 times an Earth night's.
- **Idling shrinks the store more than it adds growth.** A quarter-rate night cuts
  the store from 46 to 10 g C/m² and adds 13% to growth. The trigger matters:
  idling only in full darkness, below 1 µmol m⁻² s⁻¹, leaves a store of 19.
- **Twilight is worth 14% of the Earth-like store and a quarter of the idling
  one.** Cutting the light at sunset raises them to 53 and 13.
- **Regrowth is the weakest strategy.** A new canopy costs 264 g C/m² every cycle
  and has to be built at dawn from reserves. Growth falls to 28 at the equator, 3 at
  30° and none at 60°, and the store rises to 206–232. Regrowing over 2 days leaves
  37; over 10 days, nothing. Idling stems and roots meanwhile adds only 11.
- **Near and far sides differ by 3–4% for evergreen plants.** The far side's cooler
  equatorial days (21 against 23.5 °C) lower photosynthesis and respiration together.
  Regrowth's thin margin shrinks by a third there (19 against 28). The sides differ
  more in their night cue: earthlight lights only the near side.
- **Latitude.** At 60° the Earth-like plant grows 135 g C/m² per cycle with a store
  of 35. Idling at a quarter there needs 5.6, because twilight fills most of the
  night.
- **Calibration.** With an Earth carbon-use efficiency of 0.4 (more upkeep for stems
  and roots), the Earth-like plant grows 197 with a store of 63, idling at a quarter
  239 with 14, and regrowth nothing. An efficiency of 0.6 leaves stems and roots
  almost no upkeep in this stand.
- **Storage capacity.** Stores of tens of grams of carbon per m² are small next to
  what Earth's forests carry as non-structural carbohydrates: hundreds of grams per
  m², several per cent of their dry mass (Martínez-Vilalta et al. 2016). For small
  herbaceous plants the same store is a large share of their mass, which is where
  storage organs come in.

What the model cannot test decides the matter:
- leaves that survive about 240 dark hours without dismantling themselves;
- metabolism that slows and restarts on cue;
- plants that keep using 354 hours of continuous light.

### Limits of the whole-plant model

- The strategies are carbon scenarios. The physiology they assume is untested:
  idling without dark-induced senescence, how fast metabolism slows and recovers,
  and how fast a canopy regrows.
- Stem and root upkeep is set from an Earth carbon-use efficiency. Lunar plants'
  biomass, tissues and nitrogen are unknown.
- Twilight comes from the atlas's exponential column and unfiltered sunlight scaled
  by the shield, with one spectrum. The sky is clear.
- The GCM temperatures are 3-day means, which smooth the coldest hours before dawn.
- Leaves are at the air temperature. There are no water limits, nutrients,
  reproduction or herbivores. The continuous-light physiology is as in the canopy
  model.

## Fruit designed for the lunar day

`python -m biosphere.canopy.fruit` (about a minute and a half) writes
[results/fruit.json](results/fruit.json) (schema `terluna.biosphere.fruit-carbon/2`).
It asks how fast a large fruit can grow on the whole-plant model's stand, when in the
cycle to set it, and how much of the plant's growth it can take while the plant keeps
enough for itself. The fruit's developmental program is the design variable. The runs
cover the equator, 30° and 60° on the near side and the far side's equator.

- **The fruit** has today's watermelon's size and make-up: 7.5 kg at harvest
  (7.1–8.3 kg in Noh et al. 2013), 8.6% dry matter and 42% carbon in that dry matter
  (USDA composition). Its program is the degree-days above 10 °C it needs from fruit
  set to harvest. The design case, the *day fruit*, needs the degree-days of one
  sunlit half at its site: 199 at the equator, 187 at 30°, 167 at 60° and 165 on the
  far side's equator. Today's watermelon needs 560, from growers' 35–45 days with
  25–35 °C days and 15–20 °C nights, and is kept as the baseline. At the equator the
  program is swept from 120 to 560 degree-days.
- **Its growth** follows a beta sigmoid (Yin et al. 2003) that is fastest at half the
  program (0.4–0.6 tested). That matches the third of the harvest weight Noh et al.
  found 15 days after fruit set.
- **Its cost** is its carbon, a quarter again in growth respiration, and upkeep of
  0.01 g CH₂O per g of dry matter per day at 25 °C with Q10 = 2. Crop models give
  that value to roots, and to every organ but the leaves in tropical crops (WOFOST,
  after Spitters et al. 1989 and Penning de Vries et al. 1989).
- **One cohort per cycle.** The plant sets fruit once each cycle. A fruit that takes
  longer than a cycle overlaps the next one, and the plant carries both.
- **Two ways of feeding it.** A *fed* fruit grows at its full rate throughout, drawing
  on the plant's store in the dark. A *daylight* fruit grows only on the plant's
  surplus while photosynthesis exceeds upkeep. Its development runs on in the dark,
  and the growth it misses there is lost.

### Results

The day fruit, set at sunrise, is ripe at sunset 14.8 days later. Its fastest growth
comes at local noon, when the plant's surplus peaks, and none of its growth falls
while the plant is in deficit. It reaches full size on the daylight surplus alone and
is off the plant before the night, so the plant's store stays at its own 46 g C per m²
(10 when it idles at a quarter of its daytime upkeep).

- **How fast a large fruit can grow.** At the equator the plant's daytime surplus
  would build 0.50 kg of fruit per m² of ground a day. A 7.5 kg day fruit grows at up
  to 0.87 kg a day, 3.5 times today's watermelon's fastest (0.25 kg a day at the
  equator's temperatures). Giant pumpkins, bred for more phloem into the fruit, gain
  about 15 kg a day at their peak (Savage et al. 2015), so what one fruit can take in
  leaves room for far larger fruits. Canopy area sets the size: at a 60% fruit share,
  each 7.5 kg fruit draws on about 2 m² of canopy.
- **The program is what to design.** Programs up to the length of the day cost
  nothing at night. Longer ones spill into it. At 230 degree-days (17.5 days) the fruit
  reaches 94% of full size on daylight alone, or needs 13 g C more store when fed; at
  260 (20 days) 83% or 32 more; at 300 (24.5 days) 70% or 56 more; today's
  watermelon, 45 days, 60% or 70 more. Shorter programs (120–160 degree-days, 9–12
  days) add nothing, because the carbon then limits the harvest.
- **When to set it.** At sunrise, with the harvest at sunset. The plant is in surplus
  from about half a day before sunrise to 0.3 days after sunset (15.5 days), which
  leaves a little slack at both ends.

How much can go to fruit: the day fruit at the equator on the near side. The fruit
share is the fruit's part of all new tissue carbon. The plant's own growth is also
given as the cycles it would take to build a whole new canopy (211.5 g C). The store
for the night is the plant's own at every share: 46 g C per m², or 10 when it idles.

| Fruit share | Harvest per m² per lunar cycle | Plant's own growth | A new canopy every | Harvest when the plant idles at night |
|---|---|---|---|---|
| None | | 229 g C | 0.9 cycles | |
| 30% | 1.9 kg | 158 g C | 1.3 cycles | 2.1 kg |
| 50% | 3.1 kg | 112 g C | 1.9 cycles | 3.5 kg |
| 60% | 3.7 kg | 89 g C | 2.4 cycles | 4.2 kg |
| 70% | 4.3 kg | 66 g C | 3.2 cycles | 4.9 kg |
| 80% | 4.9 kg | 44 g C | 4.8 cycles | 5.5 kg |

- **Against today's watermelon.** Today's watermelon takes 45 days and lives through a
  night, growing about a third of its mass in the dark whatever the set time. At a 60%
  share it yields 3.5 kg per m² per cycle and either needs a store of 116 g C or
  reaches 59% of full size on daylight alone. The day fruit yields 5% more and leaves
  the plant its own store.
- **Earth.** The same stand and temperatures on Earth, with today's watermelon at a
  60% share, grow 2.45 kg per m² in 29.5 days with a store of 3.8 g C. The day fruit
  on the Moon yields 51% more per month, and 71% more when the plant idles at night.
- **Latitude and side.** At a 60% share the day fruit yields 3.3 kg at 30°, 2.2 kg at
  60° and 3.6 kg on the far side's equator (3.8, 2.6 and 4.1 when the plant idles).
- **Sensitivity** (day fruit, equatorial near side, 60% share). A growth peak at 0.4
  or 0.6 of the program changes the harvest by under 1%. Halving the fruit's upkeep
  adds 2%. More stem and root upkeep (an Earth carbon-use efficiency of 0.4) gives
  3.2 kg, with the plant's own store at 63.

### What the design asks of the plant

The model sets these requirements; whether a plant can meet them is untested.

- **A program one lunar day long:** about 200 degree-days from set to ripe at the
  equator, a third of today's watermelon's. Cucumbers already grow from 5 to 30 cm,
  picking size, in 10–15 days, unripe (Wiechers et al. 2011).
- **Filling at up to 0.9 kg a day** for a 7.5 kg fruit: 3.5 times today's
  watermelon, and a seventeenth of a giant pumpkin's peak.
- **Setting on cue at sunrise.** Parthenocarpic cucumbers set fruit without
  pollination, which frees the set from pollinators and from seed development.
- **Cells ready by sunrise.** Growth this fast is mostly cell expansion. One way to
  get there is to make the fruit's cells in the ovary late in the night, a small draw
  on the store.

### Limits of the fruit model

- Size, composition and growth curve are today's watermelon's, measured under Earth's
  day; the day fruit's program is a design requirement.
- Development follows temperature alone, and every fruit that is triggered sets. For a
  fruit that lives through a night, whether it keeps developing through two weeks of
  darkness, whether it aborts when the night comes, and how its sugar and flesh turn
  out are open. Cucumber fruits starved early make up much of it later by expanding
  their cells more (Marcelis 1993); the daylight case allows none of that.
- The stand is the whole-plant model's (leaf area index 5, with stems and roots). A
  fruiting crop's own canopy, flowering, seeds, water and nutrients are not modelled.
