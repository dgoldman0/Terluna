# Sunlight at the ground

What light reaches the ground of the Open Moon under a clear sky, by wavelength:
the photosynthetic photons, their colour, the share that arrives as diffuse sky
light, the red:far-red ratio plants sense, ultraviolet, and what a bright or dark
ground does to all of these.

`python -m illumination.surface_light.model` (about 6 minutes on two processes,
with `OPENBLAS_NUM_THREADS=1`) writes
[results/surface_light.json](results/surface_light.json) (schema
`terluna.illumination.surface-light/1`). `--atlas-dir illumination/sky/data`
adds a comparison with the sky solver's atlases when they have been generated.

## Method

The light is the atmosphere domain's radiation, run for the ground instead of for
heating:

- below 500 nm, the middle atmosphere's ultraviolet-visible calculation
  ([photolysis.py](../../atmosphere/middle_atmosphere/photolysis.py)) in 1-nm
  bins: Rayleigh scattering, O2, O3, NO2 and the other photolysed absorbers, and
  the O2 collision-induced bands;
- from 500 to 1000 nm, that model's near-infrared step: line-by-line H2O and O2
  absorption from HITRAN, the MT_CKD water continuum and collision-induced
  bands ([radiative_convective](../../atmosphere/radiative_convective/)), with
  ozone's Chappuis band, on every second level;
- the delta-two-stream solver with a pseudo-spherical direct beam, as both of
  those models use.

Columns: the design Moon (1.2 atm, a 294.9 K surface, the chosen climate's global
mean), the same air at 1.0 atm, and the middle atmosphere's Earth control (1 atm,
288 K, 277 DU of solved ozone). Each is the moist adiabat with Manabe–Wetherald
humidity up to the tropopause and the middle atmosphere's solved temperatures and
trace gases above it; for the Moon they are interpolated between its stored cases
at 288 and 298 K. The design Moon sits behind the titania film passing 5% less at
every wavelength (the author's choice of 2026-09-26). It is also run in
unfiltered sunlight, to separate what the air does from what the shield does and
to match the sky atlas; its ultraviolet is not that of an unshielded Moon, whose
air would grow ozone.

The ground is Lambertian. Over a ground of albedo A the downward light is the
light over a black ground divided by (1 − A R), where R is the atmosphere's
reflectance for diffuse light arriving from below; the direct beam does not
change. This is exact in the two-stream adding method (it matches direct
solutions to 10⁻¹⁴), so the product stores black-ground spectra and R, and any
ground follows, including one whose albedo changes with wavelength. The band
summaries are computed on the line-by-line grid for albedos 0–0.8.

## Results

Clear sky over ground of albedo 0.1, with the Sun overhead. Photons are
µmol m⁻² s⁻¹; photosynthetic photons are 400–700 nm.

| | Design Moon | Moon air, unfiltered sunlight | Moon, 1.0 atm | Earth control |
|---|---|---|---|---|
| Photosynthetic photons | 1,560 | 1,670 | 1,630 | 2,250 |
| Illuminance | 89,600 lux | 95,400 lux | 93,700 lux | 124,600 lux |
| Diffuse share of those photons (of blue) | 37% (65%) | 37% (65%) | 32% (58%) | 7% (12%) |
| Blue (400–500) ÷ red (600–700) photons | 0.58 | 0.60 | 0.61 | 0.80 |
| Red ÷ far-red, 660 ÷ 730 nm | 1.88 | 1.87 | 1.84 | 1.15 |
| Red ÷ far-red, 600–700 ÷ 700–750 nm | 2.66 | 2.64 | 2.66 | 2.18 |
| UV-A photons (315–400 nm) | 35 | | 40 | 200 |
| UV index | 0.11 | | 0.12 | 14.8 |
| Photosynthetic photons sent back to space | 33% | 33% | 30% | 13% |

- **Two thirds of Earth's light, a third of it from the sky.** Of the
  photosynthetic photons arriving at the top with the Sun overhead, 69% reach
  the design Moon's ground, 33% go back to space and 5% are absorbed in the
  air (Earth: 93%, 13%, 3%). The air alone passes 74% of Earth's photons; the
  shield (98% across 400–700 nm, then 5% less) brings it to 69%. As the Sun
  sinks the gap widens: 60% of Earth's at 30° of Sun height and 57% at 10°. The
  diffuse share grows from 37% overhead to 56% at 30° and 80% at 10° (Earth: 7%,
  12%, 30%).
- **The light is redder and short of far-red.** Rayleigh scattering at 7.7 times
  Earth's depth (0.756 at 550 nm) sends blue back to space: in the same
  sunlight the Moon's air brings to the ground 49% of what Earth's does at
  400–410 nm, 58% at 440–450, 75% at 540–550 and 79–87% in the red. The water
  vapour column (205 kg/m², ten times the Earth control's) and the O2 column cut
  deep bands into the far-red and near infrared: 46% of Earth's light at
  720–730 nm, 55% in the O2 A band at 760–770 nm and 3% at 940–950 nm. The
  red:far-red ratio in the 10-nm windows plants' phytochrome is usually measured
  by is 1.8 overhead and 2.3 at 30°, against Earth's 1.15 and 1.22; over the
  broader bands the shift is smaller (2.6 against 2.2).
- **Almost no ultraviolet.** Behind the titania film there is no UV-B at the
  ground and a sixth of Earth's UV-A; the UV index is 0.11 (Earth 14.8, which
  reproduces the middle atmosphere's own 15.0 for its Earth control).
- **A bright ground matters far more than on Earth.** The thick sky returns much
  of the light the ground reflects: its reflectance for light from below is 0.68
  at 400 nm, 0.56 at 450, 0.36 at 550 and 0.19 at 650 nm (Earth 0.21, 0.14,
  0.07, 0.035). Going from a ground of albedo 0.05 to 0.8 (fresh snow, salt,
  pale sand) raises the Moon's photosynthetic photons by 37% and its blue by
  73%, and the blue:red ratio from 0.57 to 0.83, close to Earth's; on Earth the
  same change adds 6%. At albedo 0.5 the gain is 19% and 34%.
- **The month-long day.** Integrated over one lunar day with the Sun on the
  celestial equator, the design Moon's ground gets 1,160 mol/m² of
  photosynthetic photons at the equator, 965 at 30° latitude and 472 at 60°,
  which averages 39, 33 and 16 mol/m² per 24 hours across the month. An Earth
  equinox day gives 60, 51 and 27. So a lunar day at the equator delivers about
  as many photons as 19–20 clear equinox days on Earth, in 14.8 days of
  continuous light followed by 14.8 days without.
- **Pressure.** At 1.0 atm the ground gets 5% more photosynthetic photons with
  the Sun overhead and a smaller diffuse share (32%).

### Against the sky solver

With unfiltered sunlight and ground of albedo 0.1, this model's illuminance
agrees with the spherical multiple-scattering atlases
([sky/](../sky/README.md)) to within 1–2% for the Moon with the Sun above 45°
and within 5% at 30°. Below that the two-stream's plane-parallel diffuse light
falls short: the atlas's Moon is 10% brighter at 20°, 35% at 10°, twice as
bright at 5° and seven times at 1°; Earth's agrees within 3% down to 7° and is
6% brighter at 5° and 20% at 3°. The direct beams agree within about 5%
throughout, so the difference is the diffuse light of a sky whose depth is a
large fraction of the curvature. Values below about 20° of Sun height here
therefore understate the light and its diffuse share; in the atlas the Moon's
ground is brighter than Earth's below about 7°. Scaling each height's photons by
the atlas's ratio raises the lunar day's total by 3% at the equator, 4% at 30°
and 12% at 60° (`checks.atlas_comparison`).

## Limits

- Clear sky only: no clouds, aerosols, haze or terrain. The chosen climate has
  17–25% cloud cover (climate/gcm), and in year 24 of its run the ground
  receives 89% of the clear-sky sunlight on average; how clouds change the
  colour and diffuse share is not computed.
- Global-mean columns: the water vapour, which sets the far-red, varies with
  place and weather.
- Two-stream transfer: its diffuse light is plane-parallel, so at low Sun it
  departs from the spherical sky solver (see the atlas comparison above).
- Ground albedos are neutral and Lambertian; vegetation, water and snow have
  their own spectra, which the reading rule accepts.
- The shield is the film at normal incidence; light through gaps in the swarm (a
  few 10⁻⁴ of sunlight) is left out.
- Physical light only. What plants do with it (action spectra, light
  saturation, photoinhibition, canopy light, the 354-hour day and night) is not
  modelled here.
