# Light under the water through the lunar day and night

4 November 2026. This study carries the solved sky's light down into the Open Moon's seas and lakes through the whole 709-hour lunar day: the direct beam and diffuse light for every solar elevation from noon through the long twilight, the near side's earthlight, through a flat surface and down the column of each design water, for latitudes 0–75°. It gives the euphotic, twilight and aphotic zones with explicit definitions, the sea and lake floor lit for benthic photosynthesis, the ultraviolet at depth and the depths animals can see to, with the same waters under the Earth control's sky beside them. [run.py](run.py) and [optics.py](optics.py) write [results/underwater_light.json](results/underwater_light.json) (schema `terluna.research.underwater-light/1`, with the sha256 of every input) in about a minute and a half, with the light tables and the floors' depth tables on the research drive (`/media/projectspace/terluna-research/underwater-light/`). The product also holds the spectral downwelling irradiance by 10 nm band at noon down to 200 m and photosynthetic-photon profiles at each latitude:

```sh
nice -n 10 env OPENBLAS_NUM_THREADS=1 python3 -m research.studies.underwater_light.run
```

## What it finds

**Photosynthetic light reaches 71 m in the open sea at the equator and 42 m at 75°, counted as a lunar-day dose.** At local noon the open sea keeps 1% of its photosynthetic photons to 77 m under an overhead Sun and to about 70 m at every latitude from 45° to 75°; the dose definition shortens the zone toward the poles because the whole day's light falls with latitude while the noon profile's shape changes little. The Moon's clear sky puts 1,500 µmol of photosynthetic photons per m² per second on the water at an overhead noon and 446 as the lunar-day mean at the equator (38.6 mol per m² per Earth day), against the Earth control's 2,240 and 687 at the same solar elevations. At the equator the productive coast's zone ends at 12.6 m by dose and 13.5 m at noon, the river mouth's at 2.3 and 2.5 m and the crater lakes' at 42 and 46 m. Dusk's chlorophyll, twice dawn's, moves the open sea's noon 1% depth at the equator between 74 and 81 m.

| Water, latitude | Noon 1% | Dose to 0.4 mol m⁻² d⁻¹ | Noon light to 0.1 lux | Disphotic end | Earth control: noon 1%; dose |
|---|---:|---:|---:|---:|---:|
| Open sea, 0° | 77 m | 71 m | 280 m | 689 m | 85; 85 m |
| Open sea, 45° | 70 | 60 | 241 | 610 | 72; 70 |
| Open sea, 75° | 71 | 42 | 214 | 585 | 68; 43 |
| Crater lake, 0° | 46 | 42 | 166 | 367 | 49; 48 |
| Crater lake, 45° | 41 | 35 | 142 | 324 | 41; 40 |
| Productive coast, 0° | 13.5 | 12.6 | 51 | 98 | 13.8; 14.1 |
| Productive coast, 45° | 12.5 | 10.8 | 45 | 88 | 12.4; 12.0 |
| River mouth, 0° | 2.5 | 2.3 | 7.9 | 14.2 | 2.5; 2.5 |

The Earth control's open sea keeps its light 8–14 m deeper at low latitudes: its sky is brighter and bluer, and the Moon's tall air both dims the beam and turns more of the light diffuse, whose slant paths fade faster. Toward the poles the two converge: under a low Sun most of the Moon's light arrives as diffuse skylight, less slanted on average than Earth's low beam, and the Moon's twilight adds 12% to the day's dose at 75°.

**The long twilight keeps the water lit for days after sunset, but it adds little to the dose.** The Moon's sky still puts 2,200 lux on the water with the Sun 6° below the horizon, 830 lux at 12° and 24 lux at 30°, where Earth's sky gives 2.8, 0.006 and 2×10⁻¹¹ lux. Under the open sea's water, 0.1 lux reaches 186 m at sunset, 148 m a day later (Sun 12° down at the equator), 79 m after 59 hours and 46 m after 79 hours, while Earth's water falls below 0.1 lux at every depth by the time the Sun is 12° down. The twilight brings 0.7% of the lunar day's photosynthetic photons at the equator, 1.5% at 45° and 12% at 75° (Earth: 0.02–0.4%).

**The 354 hours of night stay dark at low latitudes and never darken poleward of 65°.** Just below the surface the photosynthetic light stays under 1 µmol m⁻² s⁻¹ for 255 hours of each lunar day at the equator, 208 hours at 45° and 123 at 60°; from 65° poleward, where the Sun sinks no more than 25° below the horizon, the twilight never lets it fall that far. In the open sea, at the depth where the lunar-day dose falls to 0.4 mol m⁻² d⁻¹, the light stays below microphytobenthos' instantaneous compensation (2.8 µmol m⁻² s⁻¹) for 421 hours of each lunar day at the equator and 334 at 75°. A dose sets how much an organism can photosynthesise in a lunar day; whether it lives through the dark is a separate, open question that no Earth stonewort, seagrass or microbial mat answers.

**Clear water over the shallow seas lights about as much floor as the aquatic plan's meadows need; coastal water lights a fifth of it.** The floor depths are geography's 16 px/deg heights against the atlas's sea level: 235,700 km² of the seas' 10.62 million km² lie shallower than 30 m (the aquatic plan's 233,000 km²), and the rain-fed lakes cover 3.92 million km², 387,600 km² of it shallower than 30 m. The table gives the floor whose lunar-day mean light exceeds each threshold, if every floor lay under the named water, in thousands of km².

| Floor under | Seagrass, 5.1 mol m⁻² d⁻¹ | Charophytes, 4.0–1.3 | Coral reefs, 1.41 | Microphytobenthos, 0.4 |
|---|---:|---:|---:|---:|
| Seas, open-sea water | 140 | 171–322 | 313 | 498 |
| Seas, productive-coast water | 31 | 37–61 | 60 | 91 |
| Seas, river-mouth water | 6.4 | 7.4–12.1 | 11.8 | 17.1 |
| Lakes, crater-lake water | 207 | 237–362 | 354 | 483 |
| Lakes, productive-coast water | 85 | 96–142 | 139 | 193 |

The aquatic plan's 140,000–380,000 km² of stonewort meadow needs light like the open sea's over the shelves, where the charophytes' limit lies at 27–47 m at the equator, 20–39 m at 40° and 12–30 m at 60°; under productive-coast water the lit floor is 37,000–61,000 km², because the limit rises to 5.5–8.8 m at the equator. The lakes' clear water lights another 237,000–362,000 km² of floor, 92% of it within 30° of the equator, where the rain-fed lakes lie. In the seas, 0–30° holds 61–70% of the lit floor at every threshold. The shelves' clarity is the design's lever: the plankton and fines of a productive coast cut the meadows' lit floor to a fifth.

| Seas' floor lit above charophytes' 3% (open-sea water), thousand km² | 0–15° | 15–30° | 30–45° | 45–60° | 60–75° | 75–90° |
|---|---:|---:|---:|---:|---:|---:|
| Lit floor | 118 | 91 | 55 | 39 | 16 | 3.4 |
| Seas shallower than 30 m | 74 | 62 | 42 | 34 | 18 | 5.7 |

**The shield leaves no UV-B under the water, and by Setlow's DNA-damage spectrum a dose 830 times weaker than Earth's, carried by UV-A to three times Earth's depth.** Just below the surface at an overhead noon the Moon's water receives 10.3 W/m² of UV-A (315–400 nm) and no UV-B, against Earth's 63.9 and 2.5 W/m². Weighted by Setlow's DNA-damage spectrum, the dose is 3.2×10⁻⁴ W/m² against Earth's 0.27; 99.98% of it comes from 340–400 nm, where the action spectrum's analytic form levels off at 3×10⁻⁵ of its 300 nm value, and the part below 340 nm, where the spectrum falls steeply from its 300 nm value, is 5×10⁻⁸ W/m². That UV-A dose fades with depth as UV-A does: the lunar-day mean falls to 10% and 1% at 25 and 51 m in the open sea (Earth: 7.4 and 15 m), 9.8 and 20 m in the crater lakes, 2.2 and 4.6 m on a productive coast and 0.3 and 0.6 m at a river mouth. Setting pure water's ultraviolet absorption to zero, the bound below Mason et al.'s measurements, moves the open sea's 1% depth from 51 to 55 m; a steeper ultraviolet slope for dissolved organic matter below 350 nm changes nothing on the Moon, whose ultraviolet lies almost wholly above 350 nm, and shortens Earth's 1% depth from 15 to 12 m.

**Fish can see by day to about 690 m in the open sea at the equator, and full Earth lights the near side's water for sight to tens of metres.** Noon daylight falls to 0.1 lux, a level at which juvenile salmon stop feeding, at 280 m in the open sea, 166 m in a crater lake, 51 m on a productive coast and 7.9 m at a river mouth; it reaches the deep-sea fishes' limit at 689, 367, 98 and 14 m. Under the sub-Earth point, full Earth puts 5.4 lux on the water (the calendar's 5.44), which keeps 0.1 lux down to 62 m in the open sea, 43 m in a crater lake, 14 m on a productive coast and 2.2 m at a river mouth, and stays above the deep-sea limit to 435 m. With the Earth 25° up, 1.8 lux reaches the water and 0.1 lux reaches 36 m in the open sea. On the far side the night's only light is the twilight's.

**The season and the turned Moon move these depths by small, readable amounts.** The axis's 1.54° to the ecliptic moves the open sea's dose depth at 75° between 40 and 44 m through the year (42 m at equinox), and less at lower latitudes. Every result is tabulated by latitude: a place that the water's load moves by about 8° reads its new latitude in the zone tables. The floor areas by band use today's latitudes.

## Definitions

- **Euphotic zone, classic:** down to the depth where photosynthetic photons (400–700 nm) fall to 1% of their value just below the surface at local noon, with the water's mean chlorophyll; dawn's and dusk's beside it.
- **Euphotic zone, by dose:** down to the depth where the lunar-day mean photosynthetic light falls to 0.4 mol photons m⁻² per Earth day (4.6 µmol m⁻² s⁻¹ averaged over the 709 hours), the lowest benthic light at which Gattuso et al. (2006) found microphytobenthic communities' metabolism reported.
- **Twilight (disphotic) zone:** below the euphotic zone to the depth where noon daylight falls to the deep-sea fishes' limit. Warrant and Locket (2004) place the end of daylight for vision at about 1000 m in the clearest ocean; the limit here is the rod catch this study's optics give at 1000 m in Earth's clearest open ocean (chlorophyll 0.03 mg/m³, dissolved organic absorption 0.01 m⁻¹ at 440 nm, sea salinity) under an overhead Sun, 6.0×10⁻¹⁰ µmol m⁻² s⁻¹ weighted by a rod pigment peaking at 480 nm. Defined through the same optics, the limit moves with any common error in the attenuation much less than an absolute photon threshold would. The depth where noon light falls to 0.1 lux (0.01–1 lux, the visual feeding thresholds that Ono et al. (2010) review) marks where fish stop feeding by sight within it.
- **Aphotic zone:** below the disphotic end.
- **Benthic thresholds:** lunar-day mean photosynthetic light at the floor, in mol per m² per Earth day: seagrass 5.1 (Gattuso et al.'s median of 88 records), coral reefs 1.41 and microphytobenthos 0.4 (Gattuso et al.), charophytes 1.32–3.96 (3–9% of the sub-surface light at New Zealand charophytes' depth limits, from a search summary of Schwarz and colleagues (2002), converted with the Earth control's clear-sky equinoctial mean at 40°, 510 µmol m⁻² s⁻¹). Earth's values are reference points; the organisms are engineered.

[sources.json](sources.json) records every value with how far its source was read.

## Method

The light on the water comes from the solved spherical sky that the sea-appearance and optical-comfort studies read (`research/runs/optical_comfort/spherical/*_standard.npz`): direct and diffuse irradiance per spectral channel (360–830 nm) and per unit source, times the shielded Sun's energy, at every elevation from −90° to 90°. The Sun's arc is the present axis's at equinox, sin e = cos φ cos H over the 708.7-hour synodic day in half-hour steps. Earthlight uses the same per-source response at the Earth's elevation, weighted by the light calendar's Earth spectrum at full phase (Glenar et al.'s shape at Robinson et al.'s brightness), shared over each band's channels as the calendar's export shares it. The ultraviolet comes from the surface-light product at 1 nm (Sun 1–90°, two-stream), the case `moon_1.2atm/design` and the Earth control `earth_control/unfiltered`, with ground albedo 0.1 as the UV-biology study reads it.

At the surface the direct beam passes by Fresnel's equations for unpolarised light on flat water (n = 1.333), and the diffuse light is taken as an isotropic sky split into eight zenith streams, each with its own transmission (6.6% of an isotropic sky is reflected). Under the surface each stream fades with Lee et al.'s (2005) attenuation evaluated at its own zenith above the water, a form whose zenith term carries the refraction into the water, so the diffuse light's slant streams fade faster than an overhead beam and the profile bends as they die out. The water's constituents are the illumination domain's ([model.py](../../../illumination/water_column/model.py)) for the biosphere's design waters ([waters.json](../../../biosphere/living_water/waters.json)), with two Earth terms replaced in [optics.py](optics.py): pure water scatters as at 0.14 g/kg, 0.770 of the domain table's seawater at 38.4 g/kg (Zhang et al. 2009: seawater scatters 30% more than pure water, the salt's share close to linear in salinity), and pure water absorbs below 420 nm as the IOCCG protocol's composite gives (0.0147 m⁻¹ at 300 nm), where the domain's table holds Smith and Baker's (1981) bounds for natural water below 380 nm (0.141 m⁻¹ at 300 nm). Dissolved organic matter keeps the domain's slope of 0.0176 per nm into the ultraviolet. Each water's chlorophyll runs from its dawn value at sunrise to its dusk value at sunset and back through the night, mixed in the logarithm of the light.

The crater lakes take a design water of this study: chlorophyll 0.3 mg/m³ (0.1–1), dissolved organic absorption 0.05 m⁻¹ at 440 nm (0.02–0.3) and 0.2 g/m³ of developed fines (0–1), a rain-fed lake on a small catchment as poor in nutrients as the open sea, with a little more plankton and more organic matter and fines from its shores. Soda lakes would differ in ways the study does not compute: the nutrients a closed basin gathers can feed their plankton far beyond the open sea's, shrinking the euphotic zone toward the river mouth's few metres; carbonate precipitating in the water would raise backscattering; and the salt they accumulate would bring pure water's scattering back toward the domain table's seawater value, 30% above fresh water's at 38.4 g/kg, and raise the refractive index slightly. Their dissolved organic matter, for which the study has no value, would set their ultraviolet.

The floors' depths are built as geography's drainage builds its heights (`geography.topography.height_above_geoid` at 16 px/deg, degree-200 geoid): sea cells lie below the atlas's sea level of −1654.2 m and lake cells are the drainage product's, each lake's depth its level less the cell's height (1,420 of 1,252,355 lake cells stand above their lake's level and count as zero depth). Per 16 px/deg row the study keeps cumulative counts of floor shallower than each half-metre, so a threshold's lit floor is read from the limit depth at that row's latitude.

## Checks

- **Against the sea-appearance study, Sun 45° up.** Its own method, all light through Lee et al.'s attenuation at 45° in the domain's seawater, reproduces its 1% depths exactly (68.08, 12.40 and 2.34 m). This study gives 69.8, 12.5 and 2.35 m: fresh water's weaker molecular scattering deepens the open sea's by 1.3 m, and the surface's Fresnel transmission, the separate diffuse streams and the IOCCG absorption at 400–420 nm by 0.4 m; the coastal waters, whose particles dominate, move by under 1%.
- **The twilight** of the spherical sky agrees with the light calendar's converged point-source solve (182 scattering orders) within 1.2% from −40° to 90°.
- **Earthlight** at full phase on the water, 5.436 lux with the Earth overhead and 1.767 lux at 25°, matches the calendar's 5.437 and 1.774.
- **The floors:** 235,700 km² of sea shallower than 30 m against the aquatic plan's 233,000; 10.62 million km² of sea and 3.92 million km² of lakes, as the lake studies count them.
- 15 tests ([test_run.py](test_run.py)) cover the corrected optics, the Fresnel limits, a single overhead beam against Lee et al.'s attenuation, the DNA spectrum against TEMIS's table, the day's mean sunlight against its closed form, the area count and the product's checks.

## Limits

- The sky is clear, with ground albedo 0.1; clouds and haze, which the design climate has, would cut both the day's dose and the twilight. The surface is flat: waves' focusing and their spread of the transmitted beam are left out.
- The diffuse light is taken as isotropic for the surface and the streams; the real twilight sky brightens toward the horizon, which would make its light fade faster with depth than computed. Lee et al.'s attenuation was fitted to Earth's sunlit skies and is applied here stream by stream; the study has not compared it with a radiative-transfer solution in the water.
- The waters' constituents carry Earth's ocean optics on the biosphere's design guesses and this study's crater-lake guess; their ranges move every depth more than the method does.
- Inelastic scattering (Raman and fluorescence) is left out. In the clearest water it adds light deep in the column, most where daylight is faintest, and would deepen the disphotic end; bottom reflection is also left out.
- The charophytes' thresholds rest on a search summary converted with Earth's clear-sky light, which overstates Earth's real light and hence the thresholds; the deep-sea limit rests on a conventional 1000 m. Both are marked open in [sources.json](sources.json).
- The DNA weighting of the Moon's ultraviolet lies almost wholly beyond 340 nm, where the analytic form's level tail sets the weight; below 340 nm the dose is 5×10⁻⁸ W/m².
- Declination 0 throughout except the 75° sensitivity; the floors use today's latitudes and today's sea surface. The organisms' survival through the night is open.
