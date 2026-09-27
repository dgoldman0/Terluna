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
