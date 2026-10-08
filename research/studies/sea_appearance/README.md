# How the Open Moon's seas look

This study works out how the seas look by day, through the long evening, under
earthlight and in the far side's twilit night. It couples the
[wave studies](../../../climate/waves/README.md), the
[monthly tide](../../../geography/README.md#the-monthly-tide), the solved sky and
cloud light of [illumination](../../../illumination/), the
[optical comfort](../optical_comfort/README.md) and
[cloud twilight](../cloud_twilight/README.md) studies, and the living water of the
[biosphere](../../../biosphere/). The plan below was adopted on 2026-10-04.

## What frames it

- **The evening lasts days.** The Sun sinks 0.51° an hour: its disk takes an
  hour to set, and reaching 18° below the horizon takes 35 hours. The tall
  atmosphere keeps the sky bright far below the horizon. The
  [solved sky](../../../illumination/sky/results/solved_sky.json) gives the ground
  2,210 lux with the Sun 6° down, 267 lux at 18° and 24 lux at 30°, where Earth's
  sky would be in full night. The register defines dawn and dusk by this light
  from the [lighting calendar](#the-lighting-calendar).
- **The nearside night is earthlit.** The Earth stays nearly fixed in the
  nearside sky and is near full through the night. Its 1.9° disk gives about
  7.8 lux above the air at full Earth, 1/16,000 of the sunlight there. The Earth's
  visual geometric albedo is 0.242 (Robinson et al. 2025), below the 0.367 long
  quoted from Danjon's earthshine work
  ([illumination/earthlight](../../../illumination/earthlight/README.md)).
- **The far side's nights are twilit.** The twilight wraps round the Moon: at
  latitude φ the Sun sinks at most 90° − φ below the horizon. Only seas near the
  far side's equator grow dark around local midnight, lit then by the stars and
  an airglow not yet estimated for this atmosphere.
- **Around the limb the Earth rises and sets with the libration.** Over
  Smythii–Marginis, where the wave studies' headland stands, it spends part of
  each month above the horizon and part below.
- **The horizon is close:** 2.6 km for an eye 2 m above the water, 18.6 km from
  100 m.

## Decisions

The author's, on 2026-10-04:

- **The seas are living.** The optical work uses best guesses with stated bases
  and ranges, labelled as design guesses: plankton pigments by kind of sea
  (nutrient-poor open water, productive coasts and river mouths) and through the
  15-day day and night, dissolved organic matter from the rain-fed lakes and
  rivers, regolith fines near the coasts, and bioluminescence in the dark-night
  seas where the [lunar-cycle ecology](../lunar_cycle_ecology/README.md) supports it.
  Earth's ocean optics are the reference and engineered life a design variable.
- **Cloud ice keeps the cloud scheme's own sizes.** The cloud-light model
  ([microphysics.py](../../../illumination/cloud_light/microphysics.py)) takes each
  cell's ice size from the scheme's distributions, which CM1's radiation cap does
  not touch; 140 µm serves as a check. The atmospheric-electricity work found the
  box's ice reaching 157 µm while CM1's radiation reads ice only up to 140 µm, so
  the saved cloud fields grew under anvils of large crystals that their radiation
  took as too thick. The study reports the ice sizes in the clouds it uses.
- **Four coasts for the renderings,** each showing what the others do not:
  western Oceanus Procellarum by Russell, the most exposed coast, where a full
  Earth stands about 14° up in the east-southeast over the open water at night;
  the eastern Smythii headland, where the Earth rises and sets over the western
  sea; southern Mare Nubium, with the largest tide; and the South Pole–Aitken
  coast by Mare Ingenii under the far side's twilit night, the more representative
  far-side night, chosen first once the calendar showed how far the twilight
  reaches.

## Building blocks

1. **The surface's reflection** ([illumination/water_surface](../../../illumination/water_surface/)):
   Fresnel reflection of the sky and the glitter of the Sun's and Earth's disks
   over a distribution of slopes. The resolved slopes come from the SWAN spectra,
   as in the optical comfort study's water slopes; the shorter waves, which carry
   much of the slope, come from the unified spectrum of Elfouhaily et al. (1997),
   whose gravity–capillary range holds gravity and surface tension explicitly,
   evaluated at lunar gravity.
2. **The water's colour:** absorption and scattering of seawater from published
   tables, the living constituents above and suspended regolith fines, mare basalt
   and highland, from Apollo soil spectra. The same optics give the light under
   water for the ecology work.
3. **The light sources:** the solved sky, whose cached scattering solution the
   calendar reads down to the Sun's lowest point; the earthlight, with a published Earth reflectance
   spectrum, and the earthlit sky from the sky solver; clouds from the cloud-light
   model; the stars, with airglow a labelled placeholder until the atmosphere
   domain estimates it.
4. **Whitecaps and foam** from the waves' breaking, at lunar gravity, where
   bubbles rise more slowly.

## Work order

1. **A lighting calendar.** Hour by hour through the month at six coasts (the four
   shores of the [shore month](../../../climate/waves/nearside.md#a-month-at-four-shores),
   the Smythii headland and a South Pole–Aitken coast): the Sun's and the Earth's
   positions, the Earth's phase, the light reaching the ground from the sky and
   from the Earth, which source dominates and which mode of vision applies. With
   it, a map of where the Earth stands over every sea and where it rises and sets.
   Done: [the lighting calendar](#the-lighting-calendar).
2. **The reflection and water-colour models.** Done: [the sea's slopes](#the-seas-slopes),
   [the surface's reflection](#the-surfaces-reflection) and
   [the waters' colour](#the-waters-colour).
3. **Results by regime** (day, the long evening, earthlit night, the far side's
   twilit night):
   radiance and colour by direction of view, the glitter paths and the contrast of
   the waves. Done for clear skies: [the seas by regime](#the-seas-by-regime).
4. **Scientific renderings of the four coasts** in [visualization](../../../visualization/),
   extending the [B1 spectral path tracer](../../../visualization/reference-renderer/),
   judged against photographs of seascapes. The renderer is built and its parts are tested;
   the full frames wait for a quieter machine. The 24 frames are
   [described scene by scene](#the-scenes-described) meanwhile.

Models belong to their domains (optics in illumination, sea states in climate,
the tide and the seas in geography, the living water's guesses in the biosphere);
this folder holds the coupled runners, their results and the write-up.

## The lighting calendar

[lighting.py](lighting.py) follows six coasts hour by hour through the shore month
of the wave studies, 7 February to 8 March 2038, and writes
[lighting_calendar.json](results/lighting_calendar.json); the
[figure](../../../visualization/sea-appearance/results/lighting_calendar.png) shows it.

**The light lasts days after sunset.** The solved sky's light on the ground, read
from its cached scattering solution, falls slowly as the Sun sinks:

| Sun below the horizon | 6° | 18° | 30° | 40° | 50° | 60° | 90° |
|---|---:|---:|---:|---:|---:|---:|---:|
| Clear-sky light on the ground | 2,200 lux | 266 lux | 24 lux | 3.6 lux | 0.58 lux | 0.088 lux | 0.001 lux |

An 18% grey surface stays in day vision (above 5 cd/m², CIE 191:2010) until the
Sun is 23.6° down, 46 hours after sunset at the equator. It would reach night
vision (below 0.005 cd/m²) only with the Sun 60° down, five days after sunset.
Away from the far side's equator the Sun never sinks that far: at latitude φ its
lowest point is 90° − φ below the horizon. The last scattering order adds under
0.05% of the light even with the Sun straight below.

**The nearside night is earthlit.** Only the Smythii headland, with the Earth
low, reaches night vision, for 16 hours.

| Coast | The Earth's elevation | Light from the Earth on the ground, brightest | Hours the Earth outshines the Sun's sky | Darkest moment |
|---|---:|---:|---:|---:|
| Western Oceanus Procellarum, by Russell | 6–21° | 0.94 lux | 108 | 0.48 lux |
| Southern Mare Imbrium, Montes Carpatus | 59–74° | 4.5 lux | 183 | 2.7 lux |
| Southern Mare Nubium, by Pitatus | 53–66° | 3.8 lux | 154 | 2.9 lux |
| Mare Nectaris, by Fracastorius | 43–58° | 3.5 lux | 165 | 2.3 lux |
| Eastern Smythii headland | −10.5 to +2.7° | 0.26 lux | 120 | 0.085 lux |
| South Pole–Aitken sea, by Mare Ingenii | below the horizon | none | 0 | 0.46 lux |

Each night lasts about 350 of the month's 709 hours. On the nearside the twilight
outshines the Earth until the Sun is 40–50° down; then the Earth lights the
night, between new Earth at local noon and full Earth at local midnight. Above
the air the brightest earthlight of the month is 7.2 lux. Its spectrum is bluer
than sunlight: chromaticity (0.295, 0.303) against (0.322, 0.332). At the Smythii
headland the Earth stands above the horizon for 213 hours and lights the night
faintly from near the horizon. The Ingenii coast never sees the Earth, yet its
darkest moment keeps 0.46 lux of twilight, as bright as the darkest moment at
western Procellarum and of a night under a full Moon on Earth.

**Where the Earth stands.** Every water cell at 1-degree spacing, followed every
six hours through 2026–2045
([earth_over_seas.npz](results/earth_over_seas.npz),
[map](../../../visualization/sea-appearance/results/earth_over_seas.png)): the Earth
is always up over 99% of the nearside sea and never over 94% of the South
Pole–Aitken sea and all of Moscoviense. It rises and sets with the libration over
93% of Smythii–Marginis, 81% of Humboldtianum and 99% of Orientale, and over the
polar seas, where it skims the horizon.

**Dawn and dusk.** The [register](../../decisions.md#light-and-time) now defines
them by this light. Practical dusk lasts until the light falls to the level at
which Earth's civil twilight ends, 2.98 lux in the same solver's Earth control
([optical comfort](../optical_comfort/README.md)), with the Sun 41.4° down: 82
hours after sunset at the equator and 98 hours at 30° latitude, and poleward of
48.6° a clear night never falls that low. Day vision holds for its first 46 hours
at the equator. Where the Earth stands high, its light holds the night near the
level: the darkest moments at the three high-Earth coasts are 2.3–2.9 lux.

**Boundaries.** The light is clear-sky, from a horizontally uniform solved column
without refraction; the twilight's paths cross the terminator through night-side
air this column does not hold. The earthlight is Glenar et al.'s (2019) model
spectrum of the whole Earth, scaled to the visual phase curve that Robinson et al.
(2025) fit to the observations, and carried through the air by wavelength. The
model's spectrum is 11% bluer, blue band to red, than Robinson et al.'s band
albedos, and the Earth's turning continents and clouds move its light by about
10%.
Starlight and airglow stand in as 0.001 lux. Modes of vision are for an 18% grey
surface; clouds come in the regime results.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.lighting
python visualization/sea-appearance/lighting.py
```

The runner reads the cached scattering solution from
`research/runs/optical_comfort/spherical/`, linked to the research drive.

## The sea's slopes

How the water reflects the sky, the Sun and the Earth follows from the
distribution of its surface slopes. [slopes.py](slopes.py) assembles them hour by
hour through the shore month at the calendar's six coasts and writes
[sea_slopes.json](results/sea_slopes.json) with its hourly series
([sea_slopes.npz](results/sea_slopes.npz)); the
[figures](../../../visualization/sea-appearance/) show them. Two parts add up:

- **The waves longer than about 1 m** (the wave runs' highest frequency, 0.497 Hz):
  their slope covariance from the wave runs' directional spectra. At the four
  nearside shores these come from the nearside run's restart files at the end of
  each 48-hour segment, which hold the full spectrum of every grid point and match
  SWAN's printed spectra to their printed precision; at the Smythii headland from
  the shore history's hourly spectra at its east face.
- **The shorter waves,** down through the capillary roll-off: the short-wave regime
  of the unified spectrum of Elfouhaily et al. (1997)
  ([short_waves.py](../../../illumination/water_surface/short_waves.py)), whose
  dispersion and capillary scales carry gravity and surface tension, driven by the
  friction velocity of the GCM's surface stress at each coast, interpolated as the
  wave coupling interpolates it. At Earth's gravity the spectrum stays within 0.006
  of Cox and Munk's clean-sea slopes from 5 to 14 m/s.

**Lunar gravity moves the capillary waves.** The gravity–capillary waves are
4.2 cm long against Earth's 1.7 cm, and the slowest ripples travel at 0.147 m/s
against 0.230. Short waves form once the friction velocity passes 0.054 m/s
(Earth's 0.085), and the same wind drives them harder. Under a 3 m/s wind the
waves shorter than 1 m have a mean square slope of 0.009 on the Moon and 0.001 on
Earth; under 5 m/s, 0.026 against 0.010.

**The winds are light, and the seas are often glassy.** The median 10 m wind at
the coasts is 2.2–3.5 m/s. For 6–29% of the month's hours the stress is too weak
for short waves, and only the longer waves tilt the water.

| Coast | Median 10 m wind | Hours without short waves | Waves longer than 1 m | Short waves | All waves, median (90th percentile) | Earth's clean sea under the same wind |
|---|---:|---:|---:|---:|---:|---:|
| Western Oceanus Procellarum | 3.1 m/s | 14% | 0.015 | 0.010 | 0.026 (0.042) | 0.019 |
| Southern Mare Imbrium | 3.3 m/s | 14% | 0.015 | 0.011 | 0.026 (0.031) | 0.020 |
| Southern Mare Nubium | 3.5 m/s | 6% | 0.015 | 0.012 | 0.026 (0.038) | 0.021 |
| Mare Nectaris | 2.6 m/s | 24% | 0.016 | 0.007 | 0.029 (0.039) | 0.016 |
| Eastern Smythii headland | 2.2 m/s | 29% | 0.012 | 0.004 | 0.016 (0.033) | 0.015 |
| South Pole–Aitken sea, by Mare Ingenii | 2.5 m/s | 22% | no wave run | 0.007 | — | 0.016 |

Values are mean square slopes: medians over the month, the longer waves' over the
times with spectra (15 at each nearside shore, hourly at the headland). At the
median the seas are 1.1–1.3 times as steep as Earth's clean sea under the same
wind. They are smoother than Earth's in the calm hours, down to a third as
steep, and up to twice as steep in the strongest winds. A mean square slope of
0.026 is an rms tilt of 9°.

**Boundaries.** The short-wave laws are Earth fits. At lunar gravity they hold if
the same balance of wind and capillarity sets the short waves, and the lunar seas'
u*/c_m, up to 1.8, stays inside the range the fit spans. The GCM's stress is a
170 km, three-hourly mean, without gusts, sea breezes or the surface films of a
living sea. Gusts would break the calm hours into patches of ripples. Films damp
capillary waves in light winds: Cox and Munk measured slicks at about half the
clean slope. The nearside shores have spectra every 48 hours; the South
Pole–Aitken sea has no wave run. The shape of the slope distribution comes with
the reflection model.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.slopes
python visualization/sea-appearance/slopes.py
```

## The surface's reflection

[reflection.py](../../../illumination/water_surface/reflection.py) reflects light
off a sea of Fresnel facets whose gradients follow a Gaussian with the slope
covariance above. Smith's (1967) shadowing hides facets from the observer and
blocks reflected rays on other waves. A small bright source, the Sun or the
Earth, appears as glitter by Cox and Munk's formula, averaged over its disk. The
sky reflects through the same kernel integrated over the facets the observer
sees. A ray sent into the water or onto another wave reflects once more off a
level surface. Seawater's index follows Quan and Fry (1995).

Checks: the visible facets' weights sum to one, a flat sea mirrors the sky, the
glitter of an overhead source returns Fresnel's 2%, and the glitter formula agrees
with direct sampling of the facets within 3%. Polarization, the slopes'
skewness and peakedness, foam, and multiple reflection beyond the second bounce
are left out. The renderings will trace multiple reflection on explicit waves.

## The waters' colour

[waters.py](waters.py) takes the biosphere's
[design guesses](../../../biosphere/living_water/waters.json) of what the seas
carry through the seawater optics of
[illumination/water_column](../../../illumination/water_column/README.md) to a
reflectance in each of the solved sky's channels. Multiplied by the clear-sky
daylight on the water with the Sun 45° up, it gives the light leaving the water
([water_colours.json](results/water_colours.json),
[figure](../../../visualization/sea-appearance/)). This is the water's own colour,
seen looking down past the reflected sky.

| Water | Chlorophyll | Light leaving the water, under the Moon's sky (Earth's) | Colour | Photosynthetic light falls to 1% at |
|---|---:|---:|---|---:|
| Open sea | 0.1 mg/m³ | 63 cd/m² (102) | deep blue | 68 m |
| Productive coast, mare fines | 2 mg/m³ | 156 cd/m² (241) | olive green | 9.5 m |
| Productive coast, high-titanium mare fines | 2 mg/m³ | 125 cd/m² (193) | darker olive | 8.9 m |
| Productive coast, highland fines | 2 mg/m³ | 317 cd/m² (489) | pale green | 12 m |
| River mouth | 5 mg/m³ | 186 cd/m² (283) | tan | 1.6 m |

**The lunar daylight warms every water's colour.** With the Sun 45° up the
Moon's clear sky puts 55,490 lux on the water against Earth's 85,200, and 951
against 1,529 µmol of photosynthetic photons per m² per second. The tall air takes
more blue out of the sunbeam, so each water's own light comes out greener and
yellower than under Earth's sky: the open sea at chromaticity (0.186, 0.237)
against (0.177, 0.202).

**The plankton's month changes the open sea little.** Dusk's chlorophyll, twice
dawn's by the design guess, turns the open sea slightly greener (y 0.246 against
0.225) and lifts the 1% depth from 65 to 72 m between dusk and dawn.

**The fines set the coasts.** At the same 3 g/m³, highland fines make a coast
twice as bright as mare fines and paler; high-titanium mare fines darken it.
Suspended at lunar gravity, fines settle six times slower than on Earth.

**Boundaries.** Earth's seawater optics and its nadir reflectance model without
bidirectional effects; the fines' optics through Hapke's isotropic and
equivalent-slab approximations for grains 4 µm across, from space-weathered Apollo
soils, which seawater would weather; the contents are design guesses. The
reflected sky, the glitter and the views across the water come in the results by
regime.

```sh
python -m illumination.water_column.fetch_inputs --download
OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.waters
python visualization/sea-appearance/waters.py
```

## The seas by regime

[regimes.py](regimes.py) puts the pieces together at the four rendering coasts
for six moments of the shore month: the Sun at its highest; 4° up in the
afternoon; 4°, 15° and 35° below the horizon; and the darkest hour. For an eye
2 m above the water it computes the colour and radiance of sky and sea in every
direction:
- the solved sky lit by the Sun and, weighted by its spectrum, by the Earth;
- the sky reflected by that hour's wave slopes;
- the glitter of both disks;
- the productive coast's own light.

[regimes.json](results/regimes.json) holds readouts toward the brighter source
and away from it. The panoramas themselves sit on the research drive;
[the figures](../../../visualization/sea-appearance/) show them. The sea fills
every direction below the horizon; the coastline is left to the renderings.

**By day the sea is darker than its sky and its waves show strongly.** With the
Sun 54–87° up the sky is a pale blue, 8,000–11,000 cd/m² at 30°, whitening to
about 9,500 at the horizon. The sea just below the horizon is a blue-grey
3,300–5,000 cd/m², 0.35–0.53 times the sky above it. The spread of reflected
light over the facets is 0.34–0.59 of its mean, so the waves stand out. The
Sun's glitter lies more than 30° below the horizon, near the observer's feet.

**A low Sun lays a deep orange glitter path that lasts for hours.** With the Sun
4° up, its glitter runs from the horizon to about 15° below it, deep orange
(chromaticity 0.59, 0.40). The Sun sinks half a degree an hour, so the path
stays for most of a day. On the glassy hours at the Smythii headland and southern
Mare Nubium (mean square slope 0.002–0.003) it is narrow and reaches
590,000–830,000 cd/m². On the rippled seas of western Procellarum and the Ingenii
coast (0.025–0.031) it is broader and reaches 100,000–114,000. Away from the Sun
the sea is 0.65–0.86 times the sky above it.

**For days after sunset the sea reflects a warm arch and the Moon's shadow.**
With the Sun 4° down the sky toward it glows orange at the horizon (about
1,000 cd/m², x 0.44), and the sea beneath reflects it at 0.56–0.92 of its
brightness. The Moon's own shadow lies on the opposite horizon as a grey band
(about 360 cd/m²). At 15° down the arch has narrowed to 250–260 cd/m² and the
eastern horizon to about 60. At 35° down it is an orange-red band at the horizon
(10.4–10.7 cd/m², x 0.53–0.54) under a yellow glow 10–15° up, and the ground has
9–12 lux. Where the Earth is up, as at western Procellarum, its glitter path lies
on the sea opposite the arch.

**The darkest hour differs coast by coast.**

| Coast | The Earth | Light on the ground | What the sea shows |
|---|---|---:|---|
| Southern Mare Nubium | 53° up, 90% lit | 2.9 lux | A blue earthlit sky (about 0.5 cd/m², chromaticity 0.27, 0.30) over a blue-grey sea, with a faint pink glow toward the hidden Sun |
| Western Oceanus Procellarum | 20° up, 58% lit | 0.48 lux | The Earth's yellow-orange glitter on a glassy sea, 2.5 cd/m², ten times its sky |
| Eastern Smythii headland | rising, 1.4° up, 38% lit | 0.085 lux | A deep-orange glitter of the rising Earth across a mirror-calm sea (mean square slope 0.0007), 5 cd/m² under a sky of 0.01–0.04 |
| South Pole–Aitken sea, by Mare Ingenii | 45° below the horizon | 0.46 lux | A golden twilight arch over the southern horizon all night, 0.44 cd/m² at its foot, reflected in the sea |

**Through this air the Sun and the Earth are yellow to orange.** The
column of air is seven times Earth's. The Sun's direct light is a pale yellow at
60° up (chromaticity 0.42, 0.42, with 41% of it reaching the ground). It is
orange at 15° (0.52, 0.44, 10%) and red-orange at 5° (0.58, 0.41, 2.7%). The
Earth is bluish-white above the air and a pale yellow disk at 60° up, orange near
the horizon. The sky, which takes the scattered blue, stays blue by day and
under a high Earth.

**Boundaries.** Clear sky only; clouds from the cloud-light model come next. The
coastline, foam, polarization and reflections beyond the second bounce are left
out. The panoramas are statistical, 1° by 0.5°, without wave texture. Below a
few cd/m² the eye sees colour ever more weakly, so the night colours are what a
long photographic exposure would record. The resolved waves come from the
nearest spectrum, up to 21 hours from the moment at the nearside shores; the
Ingenii coast borrows the nearside shores' median, since its sea has no wave run.
Every coast carries the productive coast's water.

```sh
NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 <environment>/bin/python -m research.studies.sea_appearance.regimes
python visualization/sea-appearance/regimes.py
```

The runner needs numba, as the optical comfort study's environment provides, and
writes its panoramas to `research/runs/sea_appearance/`.

## The scenes, described

[scenes.py](scenes.py) sets a camera and a view for each of the 24 frames and gathers
what the products say each one holds. It writes [scenes.json](results/scenes.json):
- where the Sun, the Earth and the stars stand in the frame;
- the colours and brightness of sky and sea from the regime panoramas;
- the coast's skyline from LOLA heights at the moment's tide;
- the hour's waves;
- the water's own colour.

[scenes.md](scenes.md) describes each frame in words and numbers, for picturing the
scenes before they are rendered. It is written from [scenes.template.md](scenes.template.md),
which holds the words and takes every number from the product.

Each frame's colours come from the regime results' model, evaluated for that frame's
directions with the scene's own wave slopes and the Sun and the Earth as its camera sees
them. Sizes in pixels follow the rectilinear lens.

**Two cameras moved off the study's points.** At the 118 m resolution of the coast's
terrain ([geography/coast_terrain.py](../../../geography/coast_terrain.py)):
- **S Nubium:** the wave grid's shore point falls on an island 4 by 6 km, so the
  camera floats 1.5 km off its west shore.
- **Smythii:** the headland's point lies on its east face, with the land between it and the
  western sea, so the camera floats 1 km off the west face. It takes the shore history's
  west-face waves.

**The coasts change what the open-sea panoramas showed.**
- At the Ingenii coast the darkest hour's twilight arch stands above a range 700 to
  1,300 m high, 10 to 21 km away, which hides the arch's lowest, reddest part.
- At S Nubium the darkest hour falls at low tide, 3.1 m below the mean level, with the
  island in front of the hidden Sun.
- Every coast reflects in the water below it, where the panoramas assumed open sky.

**Few stars show.** The night skies of these frames are 0.03 to 0.8 cd/m², and the thick
air dims stars near the horizon by two to four magnitudes. Only the darkest hour at Smythii,
under a low crescent Earth, records stars: Aldebaran and seven fainter ones. No frame shows
a star to the naked eye.

**Boundaries.**
- The coast's colour is a placeholder: bare Apollo soil, without beaches or surf. Land vegetation is
  omitted to isolate coastal optics; these shorelines do not represent the intended mature
  landscape. The water is modelled as a living, productive coast. Coastal vegetation and habitats
  remain to be supplied by the domain research.
- The sea's colours are the panoramas', statistical and without wave texture.
- Whitecaps are not modelled. The windiest frames, at 5 to 6 m/s, would show the first of them.

```sh
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 <environment>/bin/python -m research.studies.sea_appearance.scenes
python -m research.studies.sea_appearance.scenes --text     # rewrite scenes.md from the product only
```

The runner needs numba for the solved sky, as the results by regime do.


## Nubium at local midnight

[nubium_midnight.py](nubium_midnight.py) adds a separate
[Nubium midnight product](results/nubium-midnight.json). It selects the lowest
Sun on the hourly calendar, rather than SN-6's minimum total illuminance. At
2038-03-07 02:47 TT, the Sun is 60.6° below the horizon and Earth is 53.25° up,
98.6% illuminated and 1.804° across. The study station receives 3.59 lux.
The offshore viewpoint keeps the original Nubium location, with a portrait view
toward Earth: 65° horizontally, 87.4° vertically, pitched 16° upward.

The nearest SWAN restart is eleven hours later; its significant wave height is
0.626 m and peak wavelength 27.3 m. The scene retains the productive-water design
guess, clear molecular atmosphere and unsolved land ecology. This adds no new
climate, atmospheric or wave simulation and does not modify the original 24 frames.

The [visualization workflow](../../../visualization/sea-appearance/README.md#nubium-midnight-and-human-vision)
traces the sea in physical radiance, projects a historical NASA Earth texture,
and applies a documented human-vision display approximation. The wide view and
an Earth-centred view have different fields of view and adaptation assumptions.
Neither is an empirical validation of what an individual person would see.

```sh
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 <environment>/bin/python -m research.studies.sea_appearance.nubium_midnight
```


## Can an observer resolve Earth?

The [visibility screening](visibility.py) tests the Nubium midnight question with
observed spatial structure and the external [HDR-VDP-3 model](https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk2023hdrvdp3.pdf),
including its age-dependent CIE99 ocular scatter, local adaptation and neural
contrast response. The [result](results/earth-visibility.json) is a conditional
prediction of **detectable interior structure**, not a reconstruction of one
person's entire visual experience. The previous display's clipped white disk
is not supported as a visibility conclusion.

The clear 1.2-atm molecular atmosphere transmits **37.9%** of direct photopic
Earthlight at this elevation: normal illuminance falls from **6.762 to 2.564 lux**.
Earth remains **1.804 degrees** across, with a disk-mean luminance of
**3,294 cd/m²** before ocular scatter. Thus substantial atmospheric extinction
and visible planetary detail are compatible. The Earth-centred sky is
0.529 cd/m². Aerosols, local lunar clouds and refraction are not included.

The spatial input is a public NASA EPIC observation from 17 November 2015.
Its 443, 551 and 680 nm bands are tested separately as **relative achromatic
patterns**, each normalized to the study's photopic illuminance. The early test
product has inconsistent L1A/L1B metadata and unspecified absolute units; those
units are not used as an absolute calibration. Historical clouds and viewing
geometry do not predict the 2038 scene's weather. The
[input manifest](visibility_inputs.json) pins every downloaded file and credits
its provider. Raw external inputs and software remain in ignored run storage.

Each textured stimulus is compared with a counterfactual that removes only
interior texture. It preserves the limb, angular size and integrated normal
illuminance; the latter agrees within 2.7e-8 relative. The comparison therefore
cannot pass merely by distinguishing a disk from an empty sky.

| Test | Predicted detection of interior-texture removal |
|---|---|
| Earth, 551 nm pattern, observer ages 24, 50 and 70 | Saturated model prediction in all three cases |
| Earth, 443 and 680 nm patterns, age 24 | Saturated in both cases |
| Earth, interior contrast reduced to one quarter, age 70 | Saturated |
| Earth, 60 versus 120 pixels per degree; 6 versus 12 degree surround | Same qualitative conclusion; statistic not numerically converged |
| Terrestrial Moon control, ages 24 and 70 | Greater than 0.9999999999 |
| Moon control, quarter contrast, age 70 | 0.156 |
| Identical-image control | Exactly zero |

`P_det` is the model's probability of detecting the difference in its specified
task; it is **not** a fraction of observers who recognize Earth. The large
`C_max` values are nonlinear detection statistics, not literal multipliers of
physical contrast. Saturation does not establish precise confidence or resolve
individual coastlines. This is static, attended, achromatic detection; gaze
history, colour recognition and individual eye conditions remain open.

The Moon control uses quantitative LROC WAC 566-nm Hapke-normalized reflectance,
not a cosmetically enhanced lunar texture. Its mean is anchored at 5,000 cd/m²
within [Schmidt's measured bright near-full-Moon range](https://spaceweather.com/swpod2009/13jan09/Perigee_moon_2009_01_11_corr.pdf).
Its diameter is 0.518 degrees and prescribed surround 0.01 cd/m². The map's
reference geometry is g=i=60 degrees, e=0; opposition-dependent spatial contrast
is not solved. Unmapped polar caps use an explicit boundary extension outside
the tested interior. **This is a measurement-anchored consistency control, not
a calibrated coastal photograph or a new human observation experiment.**

Reproduction requires the existing solved sky and earthlight inputs, Python with
NumPy/SciPy/h5py/numba, GNU Octave 10.3.0, the image 2.16.1 and statistics 1.6.0
packages, and the external HDR-VDP-3.0.7 archive. Download the manifest URLs to an
ignored input directory and extract HDR-VDP there. The runner checks source
hashes and every extracted model file against the archive before running it.
The small [Octave compatibility function](octave_compat/dirac.m) supplies only
the numeric singular placeholders used by the published CIE Fourier expression;
the official function replaces its DC bin with one. No model parameters change.
DC, monotonicity and finite optical values are separately checked by the display
runner. CPU runs use one BLAS thread, no GPU and at most two numba threads.

```sh
# Activate the environment containing the Python and Octave dependencies first.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 python -m research.studies.sea_appearance.visibility --root "$PWD" --inputs research/runs/sea_appearance/visibility-inputs --out research/runs/sea_appearance/visibility --octave /path/to/octave-cli --hdrvdp research/runs/sea_appearance/visibility-inputs/hdrvdp-3.0.7
```

For isolated Octave installs, set `OCTAVE_HOME` to that installation and provide
`--image-package-list` and `--statistics-package-list` if the packages use
separate local lists. `stimuli.json`, `external-model.json`, `predictions.json`
and the MAT stimuli/maps retain the run's inputs and evidence. The compact
committed result includes source hashes, software versions and prediction
interpretation. Focused checks exercise light conservation, unchanged limbs,
angular integration, PDS sample decoding and rejection of modified model files.

The [revised display workflow](../../../visualization/sea-appearance/README.md#visibility-informed-nubium-reference)
uses these constraints while keeping display colour, tone mapping and image
generation distinct from the visibility finding.

## Nobili lake at night

The supplemental [Nobili scene](results/nobili-night.json) is a westward view from
0 N, 75.6278 E, with the eye 2 m above the lake in the selected 28% water atlas.
Native LOLA rows at 256 pixels/degree are joined across the equator, converted
to the degree-200 GRAIL geoid and sampled at 100 m. The west rim reaches 10.7
degrees at 9.3 km; the first land along that bearing is 2.35 km away. Sub-grid
terrain detail, land reflectance and mature vegetation remain placeholders.

NASA/JPL geometry for 21 June 2033, 09:07 TT puts Earth 14.707 degrees high,
1.980 degrees across and 65.62% illuminated. The Sun is 86.34 degrees below
the horizon. The molecular-atmosphere calculation gives 0.374 lux total on
a horizontal surface. This is an astronomical example, not a weather forecast
for that date. [nobili.py](nobili.py) and [its admission manifest](nobili_inputs.json)
record the exact sources and the match to the cloud experiment's solar longitude.

CM1 equatorial-ring snapshot 367, model day 45.75, column 381 supplies a
conditional cloud and wind state. The ring has flat terrain and no resolved
cross-ring dimension. A 20-km extrusion is a stated optical scenario. The
[cloud calculation](nobili_clouds.py) carries the measured-phase-curve Earth
irradiance through molecular air and all five condensate species. It uses a
uniform finite Earth disk for scattered light. [Results](results/nobili-cloud-light.json):
Earth-derived horizontal light is 0.320 lux, with 0.0097 lux standard error in
the diffuse part; the direct disk's sampled cloud transmission is 0.459–0.482.
Median per-direction Monte Carlo standard error is 7.1%; no paths truncated.
Local cloud shadows on the rim and terrain-coupled weather are not solved.

The recorded 10-m wind is 1.074 m/s and friction velocity 0.0738 m/s. Lake
water uses a conditional equilibrium spectrum: about 0.182 m significant
height and a 6.32-m peak wavelength at inverse wave age 0.84. A wave-age
1.5 sensitivity gives 0.068 m and 1.98 m. These are assumptions about wave
development, not a lake wave-history calculation. The new full-range spectrum
uses the correction to Elfouhaily equation 41 documented by
[Mobley's Light & Water supplement](https://misclab.umeoce.maine.edu/education/Light&Water/D:/PAPER/supp6.pdf):
the short-wave term includes the low-wavenumber suppression and peak factor.
The older short-wave-only calculations retain their existing baseline.

Run the study with the admitted NASA row slices and CM1 run available:

```sh
python -m research.studies.sea_appearance.nobili --crm-case <CM1-ring-equator-directory>
NUMBA_NUM_THREADS=1 python -m research.studies.sea_appearance.nobili_clouds --photons 1024
```

Raw large arrays stay in ignored `research/runs/sea_appearance/nobili/`. The
[visualization](../../../visualization/sea-appearance/README.md#nobili-lake) consumes
these results; its final photographic illustration is not a validated human-vision
prediction.

## Two additional cloud views

The selected scenes add an Earth-visible night and a twilight cloud view to the
Nobili lake image. [Selection and rejected candidates](results/cloud-scene-selection.json)
record why each is admitted. These are conditional examples from saved weather
experiments, not forecasts for the attached astronomical dates.

| Quantity | Fecunditatis: Earth and clouds | Eastern Smythii: twilight clouds |
|---|---|---|
| Observer | 0 N, 51.4427 E | 0 N, 92.2797 E |
| Example time, TT | 26 June 2033, 05:36 | 12 June 2033, 15:22 |
| Native water depth | 303 m | 1878 m |
| Earth centre / diameter | 42.645° / 1.937°; 99.30% lit | −7.667° / 1.886°; completely below the horizon |
| Sun elevation | −51.376°; deep night, not local midnight | −3.336°; bright twilight |
| Selected cloud-top altitude | 53.79 km | 59.79 km |
| Overhead all-condensate optical depth | 0.289 | 0.542 |
| Horizontal illumination | 3.240 lux | 2555 lux |
| 10-m wind / native upwind fetch | 1.10 m/s / 250.4 km | 4.43 m/s / 45.0 km |
| Conditional significant wave height / peak wavelength | 0.187 m / 6.56 m | 1.151 m / 33.04 m |

The cloud heights are model cell-centre **altitudes**, not atmospheric scale
heights or solid cloud-tower heights. The selected columns contain about 6 km
and 2 km of qualifying cloud layers in total. The broader saved sections supply
the optical calculation. Fine ice filaments in the photographs are generated
texture; the experiment does not resolve their three-dimensional shape.

[cloud_candidates.py](cloud_candidates.py) screens the equatorial experiment's
second cycle. Native LOLA/GRAIL heights reject a far-side camera that the coarse
climate mask placed over water. All five hydrometeors enter the opacity screen:
cloud top and a clear-looking Earth ray alone had admitted overcast candidates.
The final Earth screen retained 293 observer/cloud pairs at overhead optical
depth below 0.3. A strict twilight screen retained none; widening the solar
window and admitting optical depth below 1 retained 392 pairs. These counts
are neither independent weather events nor occurrence probabilities. The final
native-terrain search covers the available equatorial eastern strip, not the
whole world's possible weather.

The [Fecunditatis](results/fecunditatis-earth-clouds.json) and
[Smythii](results/smythii-twilight-clouds.json) scene products bind the source
terrain, dated geometry and solar-phase-matched CM1 snapshots. Both assume a
100-km cross-ring cloud extent. The finite-source transport uses 512 photons
per angular direction and 16 groups, with median raw pixel sampling errors of
10.8% and 19.0%. No photon paths truncated. A separate 65,536-photon ground run
refines each flux. Independent angular integration agrees within 0.77 and
0.90 combined standard errors; the reproducible
[Fecunditatis](results/fecunditatis-cloud-transport-audit.json) and
[Smythii](results/smythii-cloud-transport-audit.json) audits do not establish
empirical accuracy or angular convergence.

The Earth scene transmits 0.485–0.541 of the direct disk through clouds, in
addition to molecular extinction. Its unresolved solar cloud modulation is a
remaining material approximation: 0.482 lux of formal clear-sky twilight is
retained, about 15% of the total. The Smythii calculation transports sunlight;
its retained formal Earth contribution is only 0.000074 lux. The faint source
is never silently set to zero after a low-photon run finds no paths.

[coastal_waves.py](coastal_waves.py) measures the native upwind shore and applies
[Elfouhaily et al. (1997), equation 37](https://doi.org/10.1029/97JC00467), with
lunar gravity, through the domain's [fetch closure](../../../illumination/water_surface/fetch.py).
It supplies wave age to the existing full-range spectrum. This remains an
Earth empirical law, conditional on steady wind; duration, inherited swell,
coastal circulation and breaking are not forecasts. The spectral sources
already pinned for earlier studies remain unchanged.

To rebuild, restore the admitted external inputs from the two `*_inputs.json`
manifests and the existing sky/CM1 products. Fail on hash mismatch. The
following is the Earth case; substitute the Smythii manifest, scene and wave
filenames and use `--photons 0 --sun-photons 512` for its angular light run.

```sh
STUDY=research/studies/sea_appearance
RUN=research/runs/sea_appearance/fecunditatis-earth-clouds
python -m research.studies.sea_appearance.nobili --inputs "$STUDY/fecunditatis_earth_clouds_inputs.json" --crm-case <CM1-ring-equator-directory> --build-dir "$RUN" --out "$STUDY/results/fecunditatis-earth-clouds.json"
python -m research.studies.sea_appearance.coastal_waves --scene "$STUDY/results/fecunditatis-earth-clouds.json" --terrain-manifest "$STUDY/fecunditatis_earth_clouds_inputs.json" --out "$STUDY/results/fecunditatis-cloud-waves.json"
NUMBA_NUM_THREADS=1 python -m research.studies.sea_appearance.coastal_clouds --scene "$STUDY/results/fecunditatis-earth-clouds.json" --out "$RUN/light" --photons 512 --sun-photons 0 --fine
NUMBA_NUM_THREADS=1 python -m research.studies.sea_appearance.coastal_flux --scene "$STUDY/results/fecunditatis-earth-clouds.json" --cloud "$RUN/light/cloud-light.npz" --photons 65536 --out "$STUDY/results/fecunditatis-cloud-ground-light.json"
python -m research.studies.sea_appearance.coastal_transport_audit --cloud "$RUN/light/cloud-light.npz" --ground "$STUDY/results/fecunditatis-cloud-ground-light.json" --out "$STUDY/results/fecunditatis-cloud-transport-audit.json"
```

The [rendering workflow](../../../visualization/sea-appearance/README.md#additional-cloud-views)
retains calculated references, first generated passes, full prompts and image
reviews. The photographic illustrations have explicit residual differences;
they are not exact naked-eye predictions. Nobili remains the admitted hill-framed
night scene, with its thinner and lower clouds.
