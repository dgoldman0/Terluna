# Storms and the lunar day in a cloud-resolving model

The GCMs (`../gcm/`) diagnose clouds and convection from formulas fitted to
Earth, on cells about 170 km across. This folder runs CM1 (George Bryan, NCAR;
release cm1r22.0, MIT-style licence), a model that resolves convection and forms
its own clouds and rain, at lunar gravity. It shows the weather people would meet
through the month-long day: storms, rain, cloud, winds from the ground to the
flight band, and the air near the ground by day and by night.

## How CM1 is set up for the Moon

[cm1_run.py](cm1_run.py) downloads CM1 at a pinned hash, patches a fresh copy of
its source for every build and compiles it outside the repository
(`TERLUNA_CM1_HOME`, by default on the external drive):

- **Gravity** is a build-time constant used by the dynamics, by the RRTMG
  radiation (which turns pressure into the mass of each layer with it) and by
  the CAPE diagnostic. The `moon` builds use GM/R² = 1.6242 m/s²; the `earth_g`
  builds keep CM1's 9.81 m/s² for a gravity control.
- **Rain, snow and ice fall at lunar speeds.** In the Morrison scheme a particle
  regime V = A D^B has a Reynolds number growing as the Best number to the power
  (B + 1)/3, and the Best number is proportional to g, so at a fixed size
  V ∝ g^((B+1)/3). Cloud droplets (Stokes law, B = 2) fall 6 times slower,
  raindrops (B = 0.8) 2.9 times, graupel (B = 0.37) 2.3 times; the caps on mean
  fall speeds scale alike.
- **The Sun keeps the lunar day**, 29.53 Earth days, at the chosen 5% dimmer
  shield's 1,173.5 W/m² above the air (the GCM run `A28_dim5`). On a domain that wraps the Moon the hour angle
  grows eastward along x, so the terminator crosses the domain once a lunar day.
- **The air is the design air**: 400 ppm CO₂ and 17.5% O₂ at 1.2 atm (the GCM
  run's gases), no ozone (the titania shield blocks the ultraviolet that makes
  it), methane, nitrous oxide or halocarbons.
- **The soil's five layers are 5–80 cm thick** instead of 1–16 cm, reaching
  below the month-long day's thermal wave (about 0.6 m deep).
- **Land and water can be laid along x from a file**, and large-scale nudging
  of temperature, vapour and wind can be confined above a height, so the model's
  own boundary layer and convection stay free.

## The equatorial ring

Case `ring` is the equator as a two-dimensional ring 10,916 km round, in 1,816
columns 6.0 km wide and 111 levels to 150 km (100 m apart near the ground,
2 km apart above 27 km). The Sun crosses it once a lunar day, so every snapshot
holds all local times at once. The 28% scenario's seas and rain-fed lakes along
the equator (43% of the ring) are laid flat at sea level; the land is a
labelled placeholder (albedo 0.20, moisture availability 0.5, roughness 10 cm).
The seas hold the GCM design case's equatorial sea surface, 298.9 K. The air
starts from the design case's equatorial profile (the GCM run `A28_dim5`, years
15–24, sea-level seas within 10° of the equator), and above 8–16 km its mean
temperature, vapour and east–west wind are held to that profile over 3 days. The
ring has no north–south dimension, so this nudging stands in for what the
Hadley circulation carries away.

Physics: compressible, with the Morrison double-moment microphysics (graupel),
RRTMG radiation every 30 minutes with cloud interaction, the MM5 surface layer
and soil, fixed sea surface temperature, a simple first-order boundary-layer
scheme and damping above 115 km. The Rayleigh damping time and the boundary-layer
scheme's asymptotic length scale are Earth's values stretched by the ratio of
gravities, as the dynamics stretch.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp         # 2-D runs use OpenMP
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run ring --hours 5     # restartable segments of one model day
climate/gcm/.venv/bin/python -m climate.crm.cm1_run status ring
climate/gcm/.venv/bin/python -m climate.crm.ring_analysis ring --from-day 29.5
```

On eight threads a model day takes 7–10 wall minutes (a 40 s base step, which
CM1's adaptive stepping stretches to at most 80 s); the two lunar days below
took 8.6 wall hours in two runs with a rest between. The analysis writes its
summary to [../results/crm/ring_ring.json](../results/crm/ring_ring.json), the
fields a page or figure needs (the circulation by local time and height, the
heaviest storm's section and the rain map) to
[../results/crm/ring_ring_fields.npz](../results/crm/ring_ring_fields.npz), and
the full arrays to `products/ring_ring.npz`, which is not kept in Git.
[visualization/equator-weather](../../visualization/equator-weather/) draws them
as a page.

## What the ring shows

The results are for the second lunar day, model days 29.5–59.0 (237 snapshots),
after a first lunar day of spin-up. The ring took about 30 days to come into
balance: for 8 days evaporation only filled the air, then a burst of storms
released what it had stored.

**The day at the equator**, land and sea (hour angle 0° is noon; 10° is 0.82
Earth days, so sunrise to noon lasts 7.4 Earth days):

| Local time | Air at 2 m | Humidity | Wet-bulb | Wind at 10 m | Cloud | Rain |
|---|---|---|---|---|---|---|
| After midnight (−145°) | 21.3 / 24.4 °C | 99 / 87% | 21.2 / 22.7 °C | 0.1 / 0.7 m/s | 36 / 1% | — |
| Dawn (−95°) | 20.5 / 24.2 °C | 99 / 86% | 20.4 / 22.5 °C | 0.2 / 1.2 m/s | 18 / 1% | — |
| Mid-morning (−45°) | 29.0 / 25.7 °C | 67 / 77% | 24.4 / 22.9 °C | 2.6 / 2.8 m/s | 0 / 0% | — |
| Noon (+5°) | 32.7 / 27.0 °C | 57 / 76% | 25.9 / 23.8 °C | 2.8 / 1.7 m/s | 6 / 2% | — |
| Mid-afternoon (+45°) | 28.8 / 26.2 °C | 73 / 83% | 24.9 / 24.0 °C | 2.3 / 2.4 m/s | 31 / 22% | 0.29 / 0.14 mm/h |
| Late afternoon (+75°) | 24.9 / 24.8 °C | 87 / 89% | 23.3 / 23.4 °C | 2.8 / 3.6 m/s | 25 / 32% | 0.50 / 0.33 mm/h |
| Sunset (+95°) | 23.1 / 24.4 °C | 94 / 91% | 22.5 / 23.2 °C | 1.7 / 2.7 m/s | 21 / 19% | 0.04 / 0.11 mm/h |
| Midnight (175°) | 22.0 / 24.5 °C | 99 / 88% | 21.9 / 23.1 °C | 0.2 / 1.2 m/s | 9 / 1% | — |

Over land the night is calm and saturated. The air cools slowly, from 23 °C at
sunset to 20.5 °C at dawn a fortnight later, and in the second half of the night
fog or cloud topped about 200 m up covers a third of the land. After sunrise the
land warms for a week: by noon the air at 2 m is 33 °C, the ground 35 °C, the
humidity 57% and the wet-bulb 26 °C, the mixed layer is 12 km deep and cumulus
towers reach 16 km without raining. Storms fill the afternoon, from about 15° to
85° past noon (1 to 7 Earth days after noon), and cool the land air back to
25 °C. The seas keep the air above them between 24 and 27 °C and rain only in
the afternoon. Over the whole lunar day the air over land ranges from 19.7 to
33.8 °C (1st to 99th percentile) with the wet-bulb at most 28.3 °C; over the
seas, 22.8 to 28.2 °C with the wet-bulb at most 25.4 °C.

**Storms.** About three storms are raining somewhere on the ring at any time,
seven in ten of them over land. A storm's rain is 30 km wide at the median (a
tenth are wider than 96 km, the widest 246 km); its heaviest rain is 3.6 mm/h at
the median, 18 mm/h for the strongest tenth and 55 mm/h at most; its cloud
reaches 28 km at the median, 68 km for the tallest tenth and 82 km at most.
Updrafts reach 12–14 m/s, gusts at 10 m reach 11 m/s, and the outflow cools the
air near the ground by up to 7.7 K. Followed through the three-hourly snapshots,
storms last 6 hours at the median, 12 hours for the longest-lived tenth and up
to 36 hours, and drift slowly, between 2.7 m/s westward and 3.1 m/s eastward.
The largest become storm systems: several cells feeding an anvil about
1,000 km wide between 20 and 75 km. A place on land gets rain above 1 mm/h in
about two spells per lunar day, all in the afternoon, most lasting a few hours
and the longest a day: 61 mm per lunar day in all (about 750 mm per Earth year),
against 34 mm over the seas.

**The circulation that follows the Sun.** Near the ground the air parts about
45° before dawn and flows toward the afternoon on both sides: through the
morning at up to 4.6 m/s (at 1 km), and back from the evening and night at up
to 4.4 m/s. The two streams meet about 20° past noon, where the storms begin;
the air rises through the afternoon and sinks at about 0.6 cm/s (half a
kilometre per Earth day) everywhere else, which keeps the mornings and nights
clear. Above 10 km the wind blows eastward, toward later local times, at all
hours, strengthening with height; its mean is the GCM's, held by the nudging.

**Winds by height** (east–west speed; the ring has no north–south wind):

| Height | Median | 99th percentile | Highest |
|---|---|---|---|
| 10 m over land / sea | 1.0 / 1.7 m/s | 5.1 / 5.8 m/s | 9.3 / 10.9 m/s |
| 1 km | 2.4 m/s | 8.4 m/s | 14.1 m/s |
| 5 km | 2.3 m/s | 6.4 m/s | 13.0 m/s |
| 10 km | 1.7 m/s | 6.1 m/s | 12.7 m/s |
| 20 km | 2.6 m/s | 6.4 m/s | 13.4 m/s |
| 30 km | 3.8 m/s | 7.2 m/s | 14.2 m/s |
| 40 km | 5.0 m/s | 8.9 m/s | 20.6 m/s |
| 55 km | 6.7 m/s | 10.0 m/s | 16.4 m/s |
| 70 km | 7.9 m/s | 12.4 m/s | 17.9 m/s |

In the 35–45 km flight band the air is quiet: its vertical motion is under
0.6 m/s 99% of the time. Storm cloud reaches the band in 1–3% of the afternoon's
columns (0.4% of all land columns), with updrafts up to 14 m/s inside.

**Water.** A lunar air column holds about 215 mm of water, seven times Earth's
for the same humidity, since at 1.2 atm under a sixth of the gravity the column
has seven times Earth's mass per square metre. Evaporation, 2.3 mm/day on
average (3.1 over land, 1.4 over the seas), is Earth-like, so water stays in the
air about 90 days, against Earth's 9. In balance, rain is 1.7 mm/day, and the
nudging aloft removes 0.5 mm/day: the moisture the model's convection lifts
above the GCM's upper-air humidity. Rain comes in episodes of a few days: over
two-day spans the ring's mean rain ranges from 0.1 to 3.3 mm/day. The quietest
spans, in both lunar days, came as the afternoon crossed the ring's widest seas
(about 40% land).

**Against the GCM.** The ring rains less than the GCM's equator (2.1 against
4.7 mm/day over land, 1.1 against 1.7 over the seas) and has less cloud (16%
against 32% over land, 7% against 26% over the seas), and it gathers its rain
into the afternoon, where the GCM spreads it from mid-morning to evening. Part
of the difference is the ring's own: moisture cannot converge on its storms from
the north and south. The GCM's air temperatures are means over its lowest
layer, centred 0.9 km up: 1–2 K below the ring's 2 m air at night (19.4–20
against 20.5–22 °C) and about 9 K below it at noon (24 against 32.7 °C).

**What the ring leaves out.** It is two-dimensional, flat and at sea level,
with a placeholder land surface of fixed moisture availability, and the seas are
held at 25.7 °C. Above 8–16 km the mean temperature, vapour and east–west wind
are the GCM's. Two-dimensional storms organise into lines more readily than
real storms, 6-km columns (1 km at Earth scale) resolve deep storms but not the
smallest clouds, and lifetimes come from three-hourly snapshots.

## Rings near the poles

Cases `ring_80n` and `ring_80s` are circles of latitude, set up as the
equatorial ring with four differences:

- **The ring is shorter**, 2πR cos φ: 1,896 km in 312 columns at 80°, and
  3,734 km in 624 columns at 70°.
- **The Sun stays low.** It still crosses the ring once a lunar day, but climbs
  at most 10° at 80° and 20° at 70°. It keeps the equinox: the Moon's 1.5° tilt
  raises and lowers the daily sunlight at 80° by about a quarter over the year,
  alternately in the two hemispheres, and the rings leave that out.
- **The wind turns.** The Coriolis force of the latitude (5.2 × 10⁻⁶ s⁻¹ at 80°,
  an inertial period of 14 Earth days) acts on departures from the reference
  wind, whose balancing pressure gradient CM1 supplies (`lspgrad = 1`). The
  north–south wind is held to zero above 8–16 km along with the rest.
- **The reference is the GCM row nearest the latitude** (80.3° or 69.2°), and
  the placeholder land evaporates as readily as the GCM's land on that row.
  PlaSim's bucket lets land evaporate min(1, W / (0.4 × 0.5 m)) of what open
  water would, the same factor as CM1's moisture availability. On the
  equatorial ring's rows that factor is 0.61, against the ring's 0.5.

| Circle | Ring | Water along it | Land wetness | Sea | Land ground (sea level) |
|---|---|---|---|---|---|
| 80° N | 1,896 km | 8% | 0.02 | 20.2 °C | 17.8 °C |
| 80° S | 1,896 km | 43% | 0.03 | 20.3 °C | 17.1 °C |
| 70° N | 3,734 km | 12% | 0.05 | 21.3 °C | 18.7 °C |
| 70° S | 3,734 km | 37% | 0.22 | 21.6 °C | 17.1 °C |
| 45° N | 7,719 km | 37% | 0.41 | 23.1 °C | 20.1 °C |
| 45° S | 7,719 km | 33% | 0.46 | 24.4 °C | 19.2 °C |

Both hemispheres are run because they share a sky but not a ground. The GCM's
air over the two poles agrees to within a degree at every level, and the Sun
treats them alike over the year. The South Pole–Aitken basin's seas, though,
cover close to half of each southern circle and about a tenth of each northern
one. The northern land also stands lower (1.4 km against 2.1 km on average) and
is drier at 70°.

**The GCM's vertical wind.** A closed ring cannot rise or sink on average, but
the GCM's air sinks over both polar caps and rises near the equator: its
overturning between the equator and the poles. Cases `ring_80n_lsw`,
`ring_80s_lsw`, `ring_70n_lsw`, `ring_70s_lsw`, `ring_45n_lsw` and `ring_45s_lsw`
are the rings above, and two at 45° for the mid-latitudes, with the
GCM's mean vertical wind at their latitude imposed through CM1's own large-scale
vertical advection (`dolsw`), which carries temperature, vapour, condensate and
wind with it. The wind comes from the mass budget of the band around the ring's
row: the zonal and time mean of surface pressure times northward wind in each
layer, less its column mean (no net mass crosses a latitude circle), converges
above each layer interface and sinks through it. The layer interfaces are
PlaSim's own, with each full level midway between them; the output's `levp`
holds midpoints between full levels instead. The air sinks at 1–6.5 mm/s from
2 km to about 60 km at all four latitudes, fading to nothing by 80 km, and
rises at 1–4 mm/s between 11° N and 11° S. Sinking at 4 mm/s through a lapse
rate 1 K/km short of the dry adiabat warms the air by about 0.3 °C a day.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp               # after a patch change
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_80n_lsw
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run ring_80n_lsw --hours 3 --threads 4
climate/gcm/.venv/bin/python -m climate.crm.ring_analysis ring_80n_lsw --from-day 29.5
```

## Tilted rings, a candidate

A ring along a great circle tilted to the equator, short of the poles, keeps
the Sun sweeping it once a lunar day while it swings between its tilt's
latitude north and south, crossing the equator twice. Two such rings, tilted
about 35° and 65°, would cover every latitude to 65° in both hemispheres in
continuous runs, and let storms and air drift between the tropics and high
latitudes. They need four changes, each testable:

- **The Sun** from each column's latitude and longitude along the great circle,
  which spherical geometry gives exactly.
- **A Coriolis parameter that varies along the ring** and changes sign at the
  equator crossings. CM1 keeps a two-dimensional array for its beta-plane
  option. The patch would fill it along x and use it in the Coriolis and
  balancing-pressure terms, and a test would check the inertial period at each
  point.
- **Land, sea, sea temperatures and land wetness per stretch** of ring, from
  the GCM along the path.
- **One upper-air profile for the whole path.** That's defensible because the
  GCM's air is nearly uniform with latitude: temperatures agree within 1.5 °C
  at every level from the equator to 80°, and vapour aloft is about 20% lower
  toward the poles.

They stay two-dimensional, so storms still organise into lines. Each is a full
great circle at the equatorial ring's cost, about 8.6 hours of machine time for
two lunar days. Whether to run them waits on what the polar rings show.

## The gravity pair

If lunar convection were Earth's stretched six times in size and duration, a
run at Earth gravity with every length and time shrunk by g_moon/g_earth would
reproduce the lunar run level for level. The microphysics does not stretch:
drops, snow and graupel fall at speeds that follow gravity to other powers, and
rain forms at rates that do not depend on it. Cases `pair_moon` and `pair_earth`
set up that comparison over a sea at 298.9 K, without radiation, because
radiative transfer does not stretch either: at Earth gravity the same pressures
hold a sixth of the absorbing mass. [pair_analysis.py](pair_analysis.py)
compares the two level by level.

Neither design tried so far convects. Holding the domain-mean air at every
height to the GCM's profile kept the boundary layer at the GCM's mean humidity,
too dry to convect: no cloud in 20 days. Starting from the ring's own air over
its seas, with the air free below 8 km, the still sea warmed the air to its own
temperature within days and nothing lifted it: CAPE rose to 300 J/kg against
15 J/kg of inhibition, the wind at 10 m stayed near zero, and no cloud formed in
a week. Convection needs something to cool the air or lift it; the next design
prescribes the ring's own radiative cooling over its seas, stretched in time for
the Earth case, in both runs.
