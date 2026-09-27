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
  for that reason.
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
