# Storms and the lunar day in a cloud-resolving model

The GCMs (`../gcm/`) diagnose clouds and convection from formulas fitted to
Earth, on cells about 170 km across. This folder runs CM1 (George Bryan, NCAR;
release cm1r22.0, MIT-style licence), a model that resolves convection and forms
its own clouds and rain, at lunar gravity. It shows the weather people would meet
through the month-long day: storms, rain, cloud, winds from the ground to the
flight band, and the air near the ground by day and by night.

The [atmospheric-electricity study](../../research/studies/atmospheric_electricity/README.md)
(2026-10-02) reads these runs' storms: their graupel, cloud ice and snow, the air
where they meet supercooled water, and what slower settling at lunar gravity does
to them.

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

**A faster build** (the author's decision of 2026-09-29, in the
[decisions register](../../research/decisions.md); tested on 2026-09-30 with
[build_check.py](build_check.py),
[build_check_ring_equator.json](../results/crm/build_check_ring_equator.json)).
The builds used CM1's own gfortran settings, `-O2` with no processor target, so
the laptop's wider vector instructions and fused multiply-add went unused. Two
rebuilds of `moon_omp` were tested, each under a label of its own:
`-O3 -march=native` (`moon_omp_o3`), and a stricter variant without fused
multiply-adds or glibc's vector exp and log (`moon_omp_o3_strict`:
`-ffp-contract=off`, and `-nostdinc`, which drops the vector-maths declarations
gfortran pre-includes). Each ran one model day of `ring_equator` from a copy of
its last restart files, side by side with `moon_omp` and with two copies of
`moon_omp` nudged in the last bit, all five on one thread each. A compiler's
changes touch the last bit of every operation at once, so the nudge flips the
last bit of every value of the scalar restart file that is neither zero nor a
whole number (and the second nudge of every other one); a first try nudging a
single value hit fields the model does not carry forward, and one value could not
reach the whole ring within a day. Neither rebuild reproduced `moon_omp` bit for
bit, the strict variant included, and both drifted from it as the nudged runs
did, and less. After the day the air at 2 m differed from `moon_omp` by
0.28–0.32 °C (root mean square) in the rebuilds and 0.31–0.33 °C in the nudged
runs, with domain means within 0.04 °C. Over the day each field compared (air,
vapour, pressure and skin temperature at the surface, rain, wind at 10 m, and
potential temperature, vapour and wind aloft) differed by 0.69–0.97 of the
larger nudge's difference with `-O3 -march=native` and 0.82–0.99 with the strict
variant. By CM1's own timings `-O3 -march=native` ran 11–15% faster than
`moon_omp` in two checks, and the strict variant 8–10%; three runs of `moon_omp`
itself agreed within 0.7%. So every OpenMP build now compiles at
`-O3 -march=native` (`OMP_DEFAULT` in [cm1_run.py](cm1_run.py)), `moon_omp_o2`
keeps CM1's own `-O2` for comparisons, and the MPI builds keep `-O2`, untested.

**Threads and reproducibility** (found 2026-10-04). On one thread a run
repeats itself bit for bit; on several it does not always. Restarted from
`box_0e`'s day-10.5 restart on 4 threads, `moon_omp` repeated the one-thread
run once and once parted from it about 85 steps in (a domain total first,
then the state), and the electrified build, run the same way twice, parted
from itself about 50 steps in. Storm weather then grows the difference into
different storms within about an hour. A result that changes only in some
runs is a race between threads, not a change of summing order. It lies in
CM1's own threaded loops on these cases' paths, since every Terluna patch in
a threaded region writes only its own points and WRF-ELEC's threaded loops
write only their own slabs (checked by compiling them with no default
sharing). Each run is a valid realization; runs that differ in one setting
compare as realizations, not as the same storms. The one-thread checks above
and the supercell's restart check stand.
A case's executable changes only when its build is rebuilt (`build moon_omp`
replaces the one every such case links to), so the runner records the executable
and thread count of every segment.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_o3 moon_omp_o3_strict
climate/gcm/.venv/bin/python -m climate.crm.build_check ring_equator moon_omp_o3 moon_omp_o3_strict --threads 1    # about 25 minutes
```

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

## Tilted rings

Four rings run along great circles tilted 45° to the equator, in two mirrored
pairs:

- **A** crosses the equator heading north at 0° E and heading south at 180° E,
  and reaches 45° N over 90° E and 45° S over 270° E.
- **A′** is A mirrored across the equator. It crosses at the same two points,
  but reaches 45° S over 90° E and 45° N over 270° E, so at every longitude the
  pair stands the same distance from the equator on either side.
- **B and B′** are the same pair turned 90° in longitude, crossing the equator
  at 90° E and 270° E.

Together they reach every longitude at every latitude up to 45°, and the Sun
sweeps each of them once a lunar day, with local time advancing steadily
around the ring.

**Crossings.** Two great circles cross at two opposite points, so the four
rings cross at 12. Four lie on the equator, at 0°, 90°, 180° and 270° E, where
A meets A′ and B meets B′. The other eight lie at 35.3° N and 35.3° S, at 45°,
135°, 225° and 315° E, where each ring of one pair meets each ring of the
other. At a crossing two runs simulate the same place under the same Sun, each
in a vertical slice facing a different way. Where they agree, a result does
not depend on the slice's direction; where they disagree, the spread measures
the error of the two-dimensional approximation there. Agreement cannot reveal
errors every ring shares, such as the missing inflow from the sides. The
comparisons are statistical: composites by local time over a patch around each
crossing through the analysed lunar day, since independent runs do not line up
storm by storm. Every ring takes its local time from its longitude, so the two
rings at a crossing see the same Sun at the same moment.

**The equatorial ring at the crossings.** It passes through the four
equatorial crossings too. It joins those comparisons only after a rerun with
the GCM's mean vertical wind at the equator (rising at 1–4 mm/s); as run so
far it differs from the tilted rings in that forcing as well as in direction.

A tilted ring needs five things beyond a ring on a circle of latitude:

- **The Sun** from each column's latitude and longitude along the great
  circle, which spherical geometry gives exactly.
- **A Coriolis parameter that varies along the ring** and changes sign at the
  equator crossings. CM1's beta-plane option carries a two-dimensional array
  of it into the Coriolis terms; a patch fills it along x and makes it act on
  departures from the reference wind, as `lspgrad = 1` does on an f-plane.
- **Land, sea, sea temperatures and land wetness column by column** from the
  GCM along the path.
- **One upper-air profile for the whole path.** That is defensible because the
  GCM's air is nearly uniform with latitude: temperatures agree within 1.5 °C
  at every level from the equator to 80°, and vapour aloft is about 20% lower
  toward the poles.
- **The GCM's mean vertical wind at each column's latitude**, which along a
  tilted ring changes from rising in the tropics to sinking in the
  mid-latitudes, so CM1's large-scale vertical advection has to vary along x.

They stay two-dimensional, so storms still organise into lines. Each is a full
great circle the size of the equatorial ring; a pair run together on 4 threads
each takes about 10 hours (the A pair took 10.1–10.3).

**How the tilted rings are built** (cases `ring_a`, `ring_a_prime`, `ring_b`,
`ring_b_prime`; `var12` is the tilt). The Sun block takes each column's
latitude and the longitude it has gained since the ring's northward equator
crossing from spherical geometry. CM1 switches its beta plane off in
two-dimensional runs; a patch keeps it on for a tilted ring and fills its
Coriolis array with 2Ω sin(latitude) column by column, acting on departures
from the base-state wind (the rows' zonal-mean wind times cos(tilt)/cos(latitude),
averaged along the path). The GCM's mean vertical wind at each column's
latitude comes from a table, `terluna_wls2d.txt`, read into a per-column
version of CM1's large-scale vertical advection. Each column has its own
surface segment: water where the 28% scenario's seas and lakes cover half a
box 1° across, sea temperatures of its latitude, and on land the ground of the
nearest GCM land cell and a land-use row for its wetness (rows 20–23 and 25–30
become the placeholder land at wetnesses 0.05–0.90).

The patches were checked in three ways. The rebuilt CM1 reproduces the first
model day of `ring_80n_lsw` bit for bit. The same day, with its vertical wind
given as a table whose columns are all alike, is also bit-for-bit identical.
A six-hour test of ring A gives the Sun's elevation within 8 × 10⁻⁷ of the
formula and the Coriolis parameter within 3 × 10⁻¹² s⁻¹ of 2Ω sin(latitude).

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_a
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run ring_a --hours 16 --threads 4
climate/gcm/.venv/bin/python -m climate.crm.ring_analysis ring_a --from-day 29.5
climate/gcm/.venv/bin/python -m climate.crm.ring_crossings ring_a ring_a_prime --from-day 29.5
```

**What the A pair shows** (the second lunar day, days 29.5–59:
[ring_ring_a.json](../results/crm/ring_ring_a.json),
[ring_ring_a_prime.json](../results/crm/ring_ring_a_prime.json)). A crosses 72%
land and A′ 63%. Both are warm and humid from 45° S to 45° N alike. Over land
the air at 2 m has a median of 22–23 °C and a hottest 1% of 32–35 °C, and over
the seas a median of 24 °C. The dewpoint averages 21–22 °C in every 15° band
of latitude and never falls below 17 °C anywhere on either ring, so neither has
a comfortable hour. Over land it rains 1.3 mm/day on A and 2.2 on A′, over the
seas 0.8 and 1.1; cloud covers 14–16% of the land, and fog 4.5–4.7%, 94% of it
at night. The stretches between 30° and 45° rain 13–94 mm per lunar day over
land, against 9–12 mm on the rings at 45°, where the GCM's air sinks faster
(1.2–2.8 mm/s against 0.2–1.7 mm/s).

**At the crossings** ([ring_crossings.py](ring_crossings.py),
[crossings_ring_a_ring_a_prime.json](../results/crm/crossings_ring_a_ring_a_prime.json)).
Each ring's patch is its columns within 100 km of the crossing. At 0° E the
ground is wet land: A has 34 land columns there, A′ 25 and 9 of lake. At 180° E
it is a lake district: A has 21 land and 13 water columns, A′ 10 and 24. The
two rings agree on the air near the ground. Their lunar-day means differ by at
most 0.3 °C in temperature and 0.1 °C in dewpoint, and through the lunar day by
0.4–0.8 °C (root mean square over the 36° bins of local time). The largest
steady difference is 0.7 °C by day at 0° E (bootstrap range 0.1–1.3 °C). Winds
agree within 0.3 m/s, both at 10 m (means 1.1–2.4 m/s) and at 40 km in the
flight band (4.1–4.5 m/s). Rain, cloud and fog differ severalfold. At 0° E, A
brings 7 mm of rain per lunar day and A′ 70 mm; at 180° E, 28 and 8 mm over
land and 7 and 35 mm over the water. Cloud cover differs by up to 0.10 (23%
against 33% at 0° E), and fog by up to 0.06. Part of that is chance, since a
patch's rain in one lunar day comes from a few storms, but the bootstrap ranges
exclude zero at 0° E over land and at 180° E over the water. So a ring's
temperature, humidity and wind at a place hold whichever way its slice faces,
while its rain and cloud at a place depend on that direction by a factor of
several.

**Against the GCM at the crossings.** The GCM here is the corrected design case
(`A28_dim5_moon`, with its Sun, tilt and sunlight split corrected; the rings
themselves were forced with the GCM before its correction). Its air is its
lowest layer carried to the ground along a dry adiabat, but its vapour is that
of the layer itself, centred about 0.86 km up. There the rings' dewpoint is
about 2 °C below its value at 2 m (5 °C at midday, about 1 °C at night). At
that height the rings are more humid than the GCM over land in the lunar-day
mean: by 1.1–1.4 °C at 0° E (19.0–19.2 against 17.9 °C) and by 5.3–5.6 °C at
180° E (19.9–20.1 against 14.6 °C). The gap lies in the night. At 0° E the GCM
is 0.4 °C more humid by day and 2.8–3.2 °C drier at night; at 180° E it is
drier around the clock, most at night. Both slice directions share that
night-time excess, so it comes from what every ring leaves out or from the GCM.
The GCM also rains more at 0° E: 180 mm per lunar day on its land cells, most
of it in the early afternoon, against the rings' 7–70 mm; at 180° E it brings
78 mm, against their 8–28. Before its correction the GCM was 2.6–4.7 °C drier
than the rings at both crossings and rained 129–143 mm at both.

**Steep rings from the corrected GCM.** Three rings set up on 2026-09-29 take
their upper air, vertical wind, sea and ground temperatures and land wetness
from the corrected design case (`A28_dim5_moon`, years 20–29). Two are great
circles tilted 70° to the equator, heading north across it at 45° E
(`ring_70_45e`) and 135° E (`ring_70_135e`). Together they reach 70° on either
side, with 3,200–4,600 km of land in each band of latitude, and they cross each
other at 62.8° N, 180° E and 62.8° S, 0° E, at 83°, where the GCM puts the
start of comfortable air and CM1 has no answer yet. The third is the equator
again (`ring_equator`), set up as a great circle of no tilt so that, like the
tilted rings, each column takes its own surface, land wetness and ground
temperature, and the GCM's mean vertical wind rises through it (up to 4 mm/s).
It crosses each steep ring on the equator, at 45°, 135°, 225° and 315° E.

**The same ring with its highlands.** Every ring so far lays its land flat at
sea level, yet 44% of the Moon's land stands above 2 km and 17% above 4 km,
and about half the land along `ring_70_45e` lies above 2 km, up to 8.9 km.
`ring_70_45e_terrain` is that ring with its ground and raised lakes at their
real heights (the atlas's heights over 0.5° boxes, the land smoothed by two
passes of a 1-2-1 filter to slopes of at most 0.19, up to 8.3 km); the flat
ring is its twin, so their difference is what the highlands do. Ground and lake
temperatures are carried from sea level to their heights at the GCM's lapse
rate. CM1 allows none of the rings' upper-air hold, imposed vertical wind or
top damping toward the domain mean over terrain, so the terrain ring runs in
its own build (`moon_omp_terrain`): the upper air is held by comparing each
column with the reference at its own height, the imposed vertical wind is
advected on the sea-level column's metrics (understating vertical gradients by
under 6% over the highest ground), and the damping layer above 115 km damps
toward the initial state. A six-hour test gave surface pressure over ground at
5.6 km of 1,087 hPa against 1,211 hPa at sea level, as the Moon's 52-km scale
height implies, air there about 5 °C cooler, and slope winds building gently.

The four started together on 2026-09-29. After about two hours `ring_equator`
and `ring_70_135e` were paused, a sixth of the way through, until the twin
shows whether terrain matters enough that every ring needs it: if the twin's
low-lying stretches differ, flat against terrain, the other rings are rerun with
their ground heights; if only the high ground itself differs, the two resume
from their restart files and are read for low-lying land. In two dimensions air
cannot flow around a mountain, only over it, so the twin may overstate how much
highlands disturb the air beyond them.

That limit decided it. For eleven model days the terrain ring behaved as
expected: low-lying places within 0.4 °C of the flat ring in air (seas within
0.1 °C by day and 0.6 °C at night), dewpoints on low land about 0.8 °C lower,
high ground about 1 °C cooler per kilometre with dewpoints 2–6 °C lower, and
the ring's first comfortable hours, at night over land at 0.5–2 km and over the
raised lakes. Around day 11, as night reached much of its high ground, gales
broke out along the whole ring: a median wind at 10 m of 7–10 m/s against 1–2
in the flat ring, jets of 30–35 m/s through the lowest 2 km and 55 m/s aloft,
as strong over the seas and level land as on slopes, by day and by night, and
on slopes blowing downhill only half the time, so not cold air draining off the
highlands. Air forced over 8-km mountains in two dimensions is the likely cause,
perhaps helped by the steep slopes in CM1's terrain-following levels; highlands
the air can flow around very likely do not do this. The gales also cut its time
step from 80 s to about 23 s. It was stopped at day 13 (restart files kept); a
sound test of highlands needs three dimensions. The flat twin runs on alone.

**The plan for highlands** (the author's, 2026-09-29, recorded in the decisions
register). The rings stay flat and are read for low-lying land, so
`ring_equator` and `ring_70_135e` resume from their restart files and finish
beside `ring_70_45e`. The highlands' climate comes from the flat rings corrected
for height, about 1 °C cooler per kilometre with lower dewpoints as the terrain
ring's first eleven days showed, checked against the GCM's highland cells; that
covers how warm and humid high ground is and how often it is comfortable. The
highlands' own weather, winds up and down slopes, rain on the slopes facing the
wind and cold air pooling in valleys, needs three dimensions: a box the size of
`box_0e`, 385 km square, over a stretch of high ground on `ring_70_45e`'s path,
with its ground heights in the terrain build and the ring's day-night forcing.
That fits one overnight run on the laptop, and whether to run it is decided
after the flat ring's results. There is no budget for rented compute, so larger
highland domains stay out of reach.

**What the three flat rings show** (their second lunar day, days 29.6–59.0:
[ring_ring_equator.json](../results/crm/ring_ring_equator.json),
[ring_ring_70_45e.json](../results/crm/ring_ring_70_45e.json),
[ring_ring_70_135e.json](../results/crm/ring_ring_70_135e.json), the crossings
[crossings_ring_70_45e_ring_70_135e.json](../results/crm/crossings_ring_70_45e_ring_70_135e.json),
[crossings_ring_70_45e_ring_equator.json](../results/crm/crossings_ring_70_45e_ring_equator.json)
and [crossings_ring_70_135e_ring_equator.json](../results/crm/crossings_ring_70_135e_ring_equator.json),
and [ring_comfort.py](ring_comfort.py) with
[comfort_flat_rings.json](../results/crm/comfort_flat_rings.json)). All three
finished on 2026-09-30, in 11.5–12.3 wall hours each on 2–4 threads. Their
second lunar day has nearly settled: over its last third the land runs
0.1–0.5 °C warmer, with dewpoints 0.3–0.5 °C higher, than at the same local
times a lunar day before, and the seas 0.0–0.4 °C warmer.

Low-lying land is never comfortable. At sea level the land of all three rings
has 0–4 comfortable hours per lunar day in every band of latitude, and at most
94 on the loose band's terms. For 58–67% of the lunar day the air is in the
comfortable range but too humid, and for 30–42% it is too hot. Ranges span the
rings that cross each band:

| Latitude (land columns) | Sunlit air: median, hottest 5% | Night air, median | Dewpoint at 2 m, day and night | Hours too hot, too humid |
|---|---|---|---|---|
| 0–15° (1,330) | 27.7–29.7 °C, 33.2–33.9 °C | 22.0–22.4 °C | 21.4–23.4 °C | 240–290, 419–468 |
| 15–30° (262) | 29.0–29.8, 34.6–34.8 | 21.3–22.1 | 20.6–22.1 | 253–281, 427–455 |
| 30–45° (531) | 30.5–31.0, 35.6–35.8 | 20.8–21.7 | 19.0–20.9 | 284–299, 410–425 |
| 45–60° (565) | 29.8–30.2, 33.5–34.1 | 20.7–21.3 | 18.9–19.6 | 273–284, 422–434 |
| 60–70° (762) | 26.9–27.3, 30.5–31.3 | 20.0–20.9 | 18.1–19.0 | 211–217, 464–476 |

The rings agree where they cross. At 62.8° S, 0° E the steep pair's land
differs by 0.3 °C in air and 0.2 °C in dewpoint over the lunar day (bootstrap
ranges −0.1 to 0.7 and −0.5 to 0.9 °C), and at 62.8° N, 180° E by 0.7 and
1.7 °C (0.3 to 1.1 and 1.2 to 2.2 °C); there `ring_70_135e` has the only
comfortable hours at any crossing, 10 strict and 120 loose. Neither rains at
either. At the four equatorial crossings the steep rings run within 0.7 °C of
the equatorial ring in air and 0.5–1.0 °C drier in dewpoint, and they rain
2–30 times less, 5–84 mm per lunar day against 10–290 mm. As with the A pair, a
ring's temperature and humidity at a place hold whichever way its slice faces,
and its rain does not. Along its path the equatorial ring rains 3.9 mm a day
over land against 5.6 at the GCM's nearest land cells, and 3.4 against 3.1 over
the seas; the steep rings rain 0.1–0.2 mm a day over land against the GCM's
1.3–1.9, and 0.03–0.1 over the seas against 1.2–1.7.

Against the corrected GCM the gap lies over land. Where a ring and the GCM both
have sea at sea level, the ring's air is 1.8–2.4 °C warmer and its dewpoint at
the height of the GCM's lowest layer 0.8–2.5 °C higher. Where the GCM cell
nearest a land column is land below 1 km (580–730 m on average), the ring's air
is 2.4–4.7 °C warmer and that dewpoint 1.6–5.1 °C higher, and at 2 m the ring's
dewpoint stands a further 1.4–3.5 °C above it. The GCM counts 116–458
comfortable hours per lunar day at those cells, and the rings none.

Part of that gap poleward of 30° is the rings' own ground. Each land column
takes the land-use class nearest the GCM's wetness beneath it
(`WETNESS_CLASSES` in [cm1_run.py](cm1_run.py)), and the driest class has a
moisture availability of 0.05. The corrected GCM's ground poleward of 45° is
drier than that almost everywhere (median wetness 0.005–0.013, with 94–100% of
the columns below 0.05, and 36–39% at 30–45°). There the rings' land evaporates
0.3–0.5 mm a day against the GCM's 0.1–0.2; equatorward of 45° the two
evaporate within 15% of each other. The floor does not explain the rest: the
equatorial ring's land is as wet as the GCM's, evaporates less (2.9 against
3.7 mm a day), and is still 2.4 °C warmer and 3.0 °C more humid than the GCM's
low cells. Since 2026-09-30 the classes reach down to 0.003 (0.003, 0.01 and
0.02, in the land-use table's residential and commercial rows, which no case
uses and CM1 treats as any other land), so cases set up from then on take the
GCM's dry ground as it is; the three flat rings ran with the floor.

Corrected for height, comfort appears on high ground poleward of 30°. Over
model days 1–11 the terrain twin's land ran 1.0 °C cooler per kilometre of
height by day and by night, with dewpoints 1.0 °C lower per kilometre by day and
1.2 °C at night (fits through zero over its height classes above 0.5 km;
6.8–7.5 °C cooler and 6.1–7.9 °C drier at 7 km), and its land below 0.5 km ran
0.8–0.9 °C drier than the flat ring's. The GCM's own land cells change more
slowly with height: 0.7, 0.5 and 0.1 °C cooler per kilometre at 0–30°, 30–60°
and 60–90°, with lowest-layer dewpoints 0.4, 0.2 and 0.6 °C lower. So the check
does not bear out the twin's rate, and the GCM's comfort comes from air that is
cooler and drier at every height. The rings' land has a median height of
1.8–2.4 km, 42–55% of it above 2 km. With the twin's differences added at each
column's height, the comfortable hours per lunar day (strict), with the GCM's
at the nearest land cell in brackets, are:

| Latitude | Below 0.5 km | 0.5–2 km | 2–4 km | Above 4 km |
|---|---|---|---|---|
| 0–30° | 0 (99) | 0 (160) | 3 (233) | 16 (229) |
| 30–60° | 33 (510) | 38 (501) | 100 (516) | 141 (550) |
| 60–70° | 9 (626) | 50 (613) | 99 (638) | 249 (611) |

Over the rings' 3,450 land columns that is 42 hours per lunar day (103 on the
loose band's terms) against the GCM's 387 (490): 18% of the land has at least
100 hours against the GCM's 92%, and 5% at least 200 against 70%. The tropics
stay too humid at every height, with corrected dewpoints of 16–22 °C.

So CM1 and the GCM disagree about what comfort turns on, how warm and humid the
air over land is near the ground. The decisions register accepts the corrected
climate's drier land for the comfort the GCM counts; the rings find comfort
only on high ground poleward of 30°, at about a tenth of the GCM's hours. Open
for the author: a rerun of the steep rings with those drier classes, which
shows how much of the gap poleward of 30° is the rings' own ground; the 3-D box
over high ground, which the plan leaves to these results; and an independent
model for the rest of the gap (ROCKE-3D is open). The rain test, which the
author chose to run first on the faster build, found that the slow fall of rain
at lunar gravity explains none of the gap in warmth and humidity (below, "The
rain test"). With no affordable way left to tell which model is right, the
author paused the climate work on 2026-09-30 (decisions register), carrying
comfort as the range between the two.

```sh
climate/gcm/.venv/bin/python -m climate.crm.ring_analysis ring_70_45e --from-day 29.5    # likewise ring_equator and ring_70_135e
climate/gcm/.venv/bin/python -m climate.crm.ring_crossings ring_70_45e ring_70_135e --from-day 29.5    # and each steep ring with ring_equator
climate/gcm/.venv/bin/python -m climate.crm.ring_comfort ring_equator ring_70_45e ring_70_135e
```

## A three-dimensional box at a crossing

Case `box_0e` tests what the rings' two dimensions do to the air near the
ground, the cloud and the rain. It is ring A's patch where it crosses A′
(0° N, 0° E) given a second horizontal dimension: 64 × 64 columns at the rings'
6.0 km spacing, 385 km square and wrapping round on all four sides, on their
111 levels. It keeps what ring A has there:

- **the upper air and the wind along the ring**, the same reference profile and
  nudging, with x running along ring A (45° east of north);
- **land as wet as ring A's columns within 100 km** (0.86 on average, set as the
  0.90 class), all of it land;
- **the GCM's mean vertical wind there**, rising at up to 4.7 mm/s;
- **the Sun of the site**, the same across the box. At this size the box spans
  12.7° of longitude, a day of local time, which a box that wraps round cannot
  carry. On the equator there is no Coriolis force.

A box cannot make the planet-wide day–night circulation that the rings carry
along their length. Keeping only its longest waves (2,180 km and more), that
circulation carries air along the rings and lifts and sinks it. At the crossing
it rises through the afternoon (5 mm/s in the lowest 3 km, up to 23 mm/s at
3–8 km) and sinks through the night and morning (3–4 mm/s low down, 6–13 mm/s
above). In the lowest 3 km it warms and dries the night and early morning by up
to 0.5 K/day and 0.35 g/kg/day, moistens the late morning by up to 0.5 g/kg/day,
and cools the hours around noon and the afternoon by up to 0.8 K/day: as much as
the GCM's vertical wind does there. The box takes it as a heating and
moistening prescribed by local time and height (`var11 = 1`, read from
`terluna_lsadv.txt`): the mean of A and A′ over their second lunar day,
refitted with the day's first four harmonics, up to 16 km. A first run with
only the part along the rings, stopped at day 47, rained in the mornings and
stayed clear and hot through the afternoons. The box then differs from ring A at
the crossing in its second dimension, its uniform land and Sun, its size and
its surface wind, which lacks the circulation's own winds (1.1 against
2.1 m/s by day in the first run); the rings' circulation itself stays untested.

The patch was checked in two ways. The rebuilt CM1 reproduces the first day of
`ring_80n_lsw` bit for bit. In two three-hour runs of the box, a heating and a
moistening prescribed at 25–35 km, falling between two local-time bins, arrive
at 98% of the amounts expected (the nudging takes back the rest) and at 99.8%.
Two lunar days take about 7 hours on 8 threads, with 0.5 GB of memory and
about 23 GB of output and restart files.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e     # after ring_a and ring_a_prime
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e --hours 12 --threads 8
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e --from-day 29.5
```

**What the box shows** (the second lunar day, days 29.5–59, against the rings'
land within 100 km of the crossing:
[box_box_0e.json](../results/crm/box_box_0e.json)). Two lunar days took 6.6
hours on 8 threads. The box's day has the rings' shape: clear mornings, storms
from early afternoon with the most rain late in the afternoon, and a saturated
night. Its rain and daytime cloud fall between the two rings': 38 mm per lunar
day against 7 on A and 70 on A′, and 14% daytime cloud against 12% and 23%. It
is warmer and more humid than both around the clock, at the same relative
humidity (88% over the lunar day): by 1.6–1.9 °C in the air at 2 m (25.4
against 23.5–23.8 °C), 1.8 °C in the dewpoint at 2 m and 1.6–1.8 °C at the
height of the GCM's lowest layer, all far outside the day-to-day scatter. Its
nights are foggier: fog covers the box 56% of the night against 16–25% on the
rings, under low cloud up to about 1 km, which keeps it warmer. Two limits of
the box could run it warm. Its surface wind lacks the circulation's own winds
(0.65 against 1.1–1.3 m/s over the lunar day, and nearly calm at night). Its
prescribed heating cannot respond to its own warmth, where a warmer stretch of
ring drives a circulation that carries the heat away.

Against the corrected GCM the box stands further off than the rings. Its air
measured the GCM's way averages 25.5 °C against the GCM's 22.0 °C, and its
dewpoint at the GCM's layer 20.8 against 17.9 °C. It rains 38 mm per lunar day
against the GCM's 180, with 14% daytime cloud against 39% (before the GCM's
correction: 21.6 °C, 16.4 °C, 129 mm and 32%). No hour is comfortable in the
box or on either ring. So the rings' two dimensions do not explain why the
cloud-resolving model and the GCM disagree: a third dimension leaves rain and
cloud within the spread of the two slices and makes the air warmer. The
disagreement lies in what every CM1 case here shares (the placeholder land
surface, the large-scale forcing, the physics at lunar gravity) or in the GCM's
own convection and clouds.

## What else divides CM1 and the GCM at the crossing

With the rings' two dimensions ruled out, three more pieces were tested at the
same crossing: the land surface every CM1 case shares, the radiation of both
models, and where CM1's rain goes.

**The land surface** (cases `box_0e_small`, `box_0e_small_rough`,
`box_0e_small_soil` and `box_0e_small_gcm_land`: boxes 192 km square, a quarter
of `box_0e`, under the same forcing, 2.8–3.0 hours each on 4 threads). The GCM's
land (PlaSim as run) has a roughness length of 2 m everywhere against the
placeholder's 10 cm, a soil that holds heat about a quarter better, and the same
albedo of 0.20; at 0° E its bucket stays full. The small box with the
placeholder land runs 0.4 °C warmer and more humid than `box_0e`, with the same
rain and daytime cloud. Against it, over the second lunar day:

| Land given the GCM's | Air at 2 m | Dewpoint at the GCM's layer | Rain per lunar day |
|---|---|---|---|
| Roughness (2 m) | +0.05 °C | +0.20 °C | +1 mm |
| Soil | 0.00 °C | +0.07 °C | 0 mm |
| Roughness, soil and full wetness | −0.10 °C | +0.28 °C | +7 mm |

The ground divides its energy the same way in all four: 139–140 W/m² of net
sunlight and −6 to −10 W/m² of net infrared, 26–27 W/m² carried into the air as
heat and 101–105 W/m² evaporating water. The GCM's ground took the same
130 W/m² in all before its correction, carrying 8 W/m² into the air and
evaporating 122; corrected, it takes 121, carries 5 and evaporates 116. So the
land accounts for a few tenths of a degree of the gap in the air (4.2 °C against
the GCM before its correction, 3.9 °C after) and none of the gap in the dewpoint
(5.0 °C before, 3.5 °C after); rougher land halves the wind at 10 m and changes
little else.

**Radiation** ([radiation_check.py](radiation_check.py),
[radiation_check_box_0e_small.json](../results/crm/radiation_check_box_0e_small.json)).
The repository's line-by-line model, against which the GCM was calibrated,
runs on each model's own clear column at the crossing. CM1's radiation (RRTMG)
is within 2.5 points of it in the shares of sunlight reflected, absorbed in the
air and reaching the ground, near noon and with the Sun lower, and within
8 W/m² in each over the lunar day. The GCM's is not. On its own column, which
holds 195 kg/m² of water vapour, it reflects 26% of the clear-sky sunlight
arriving over the lunar day, as line-by-line does, while its air absorbs 23%
against line-by-line's 31% and its ground takes 51% against 43%. On the
328 W/m² that reaches its column that is 26 W/m² too little in the air and
25 W/m² too much in the ground, and the GCM shows the same 8 points in every
column of a wider check ([GCM README](../gcm/README.md), "How it splits the
sunlight"). The 328 W/m² is itself short. CM1, like the Moon, has 373 W/m²
arriving over a lunar day at 0° E; the GCM's Sun sweeps unevenly and has
Earth's tilt, which take about 30 and 15 W/m² ([GCM README](../gcm/README.md),
"How the Moon is set up"). The clear air's absorbed sunlight, 111 W/m² in CM1
against 76 in the GCM, divides as follows: the GCM's missing sunlight about
10 W/m², its split 29, CM1's wetter column (271 kg/m²) 9, and CM1's radiation
−13, of which 8 is RRTMG absorbing less than line-by-line and 5 is the Sun's
own spectrum it takes in place of the shield's. These are the GCM before
2026-09-29. With its Sun and tilt corrected and its sunlight split
recalibrated, 373 W/m² reaches its column at 0° E and its clear air takes 31.4%
of it, as line-by-line's does ([GCM README](../gcm/README.md), "Recalibrated").

**Where CM1's rain goes.** The microphysics keeps a water budget. Over the
second lunar day in the small box, 304 mm of water condensed; 52 mm evaporated
back from cloud, 197 mm evaporated from rain, snow and graupel on the way down,
and 34 mm reached the ground: 11% of what condensed, with 85% of the falling
water evaporating. The large box loses the same share; the rings lose 57–73%.
At lunar gravity raindrops fall about three times slower and graupel 2.3 times,
through a daytime mixed layer 8–10 km deep at 65–75% humidity, so most rain
evaporates back into the lower air. That keeps the lower air humid, and the
ground evaporates less into it. The GCM's rain scheme, fitted to Earth, brings
129 mm per lunar day to the ground here, about what its ground evaporates.

So neither the rings' two dimensions nor the land surface explains why the
cloud-resolving model and the GCM disagree. Three differences stood out: the
GCM's Sun, which brought 12% too little sunlight to 0° E; its radiation, which
put about 8% of the sunlight that did arrive into the ground where line-by-line
puts it into the air; and the fate of rain at lunar gravity, which CM1 follows
through the fall and the GCM's scheme takes from Earth. With the first two
corrected (`A28_dim5_moon`), the GCM at 0° E is 0.4 °C warmer and 1.5 °C more
humid at its lowest layer and rains 40% more. The box stays 3.5 °C warmer and
2.9 °C more humid, and the gap in rain widens: 180 mm per lunar day in the GCM
against 38 in the box. That left the fate of rain at lunar gravity as the
largest known difference. The three flat rings set up since from the corrected
run (above) keep the gap over land.

**The rain test** (the author's go-ahead of 2026-09-30:
`box_0e_small_earth_fall`,
[box_box_0e_small_earth_fall.json](../results/crm/box_box_0e_small_earth_fall.json)).
The small box ran again from its inputs as written, on a build where rain,
snow, graupel, ice and cloud droplets fall at Earth's speeds and everything
else, droplet activation included, keeps lunar gravity; the pair differs in the
fall alone, and both keep the forcing of the GCM before its correction. Over the
second lunar day, from CM1's running totals of its microphysics:

| | Lunar fall | Earth fall |
|---|---|---|
| Water condensed | 304 mm | 231 mm |
| Evaporated back from cloud | 52 mm | 53 mm |
| Evaporated from rain, snow and graupel on the way down | 197 mm | 108 mm |
| Reaching the ground | 34 mm | 60 mm |
| Share of the falling water evaporated | 85% | 64% |

With the fall at Earth's speeds 77% more rain reaches the ground, still a third
of the GCM's 180 mm. The air near the ground turns warmer and more humid: 1.0 °C
warmer at 2 m (26.8 against 25.8 °C) with a dewpoint 1.0 °C higher at 2 m and
at the height of the GCM's lowest layer (22.3 against 21.3 °C there, where the
GCM has 17.9 °C), at the same relative humidity (88%). Rain that falls quickly
leaves less cloud behind, 9% of the box by day against 14% and 75% at night
against 85%, so more sunlight reaches the ground, while the heat carried into
the air falls by 1.6 W/m² and the evaporation rises by 4.0 W/m². So the slow fall at lunar gravity
keeps the lower air cooler and less humid, and explains part of the gap in rain
and none of the gap in warmth and humidity. The rest lies in what every CM1
case shares (the large-scale forcing held from the GCM, the boundary layer, the
placeholder land) or in the GCM's own convection and clouds; an independent
model would show which.

```sh
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_small_rough --from-day 29.5 --reference box_0e_small
climate/gcm/.venv/bin/python -m climate.crm.radiation_check box_0e_small --from-day 29.5    # about a quarter of an hour on 8 cores
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_earth_fall
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_small_earth_fall    # box_0e_small's inputs as written
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_small_earth_fall --hours 5 --threads 4    # 2.3 hours
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_small_earth_fall --from-day 29.5 --reference box_0e_small
```

**The GCM's daytime cloud** (the author's go-ahead of 2026-10-01: GCM run
`A28_dim5_moon_thin_convective`,
[comfort_flat_rings_A28_dim5_moon_thin_convective.json](../results/crm/comfort_flat_rings_A28_dim5_moon_thin_convective.json)).
At 0° E by day the GCM has 42% cloud against the box's 14%, and the rain test
showed how directly cloud sets the warmth near the ground. Nearly all of the
GCM's daytime cloud there is convective: thunderstorm columns (5% at every
middle level) and their anvils (13–14% at 20–60 km), with 12 mm of convective
rain a day. Its humidity thresholds make little cloud by day, since the relative
humidity over tropical land averages 39–74% by level against thresholds of
85–98%. So the test thinned the convective cloud, whose cover is
0.245 + 0.125 ln(convective rain in mm per `convective_day_s`), by counting the
rain per 2,400 s in place of a day (cover near its floor of 0.05 even at 12 mm a
day). It branched from the corrected design case's last year with every input and
other setting checked identical, and ran 15 years.

Over years 5–14 the daytime cloud at 0° E fell to 26% and the sunlight reaching
the ground there by day rose from 246 to 272 W/m². The whole Moon warmed:
cloud cover 23.5% to 20.0%, planetary albedo 0.280 to 0.268, the surface 294.9
to 296.6 K, on its way to about 297.1 K. Near the ground the air warmed 1.5 °C
by day and by night at 0° E and along every ring, with dewpoints about as much
higher, which brings the GCM within 0.2–0.7 °C of CM1's air over the seas and
0.8–2.8 °C over low land. It warmed alike from the ground up: 1.5 °C at 0.9 km,
1.45–1.5 °C at 4 and 10 km and 1.6 °C at 18 km, over the globe and over the
tropical land, and 2.0–2.7 °C higher up. The air near the ground stays as warm
as before against the air above 8 km, which CM1 holds to the GCM's, and CM1
forced from this run would warm with it. So the GCM's daytime cloud sets how
warm the whole Moon is (its comfortable hours on the rings' land fall from 387
to 290 per lunar day with these 1.7 K). The warmth of the lowest air against the
air above, which is where the two models differ, stays as it was. That leaves the
lowest few kilometres: each model's boundary layer and convection, which the GCM
resolves with two layers (at 0.9 and 4 km), and CM1's nights, fogged 85% of the
time at 0° E against the GCM's 22% cloud, which is untested.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28 --set sunlight_scale=0.95 --set convective_day_s=2400 \
    --start-from climate/gcm/runs/A28_dim5_moon/model/MOST_REST.00029 --folder A28_dim5_moon_thin_convective --years 15
climate/gcm/.venv/bin/python -m climate.gcm.compare A28_dim5_moon:20-29 A28_dim5_moon_thin_convective:5-14
climate/gcm/.venv/bin/python -m climate.gcm.compare --levels A28_dim5_moon:20-29 A28_dim5_moon_thin_convective:5-14
climate/gcm/.venv/bin/python -m climate.crm.ring_comfort ring_equator ring_70_45e ring_70_135e --gcm A28_dim5_moon_thin_convective:5-14
```

**The GCM's low cloud at night** (the author's go-ahead of 2026-10-01: GCM run
`A28_dim5_moon_low_cloud`,
[comfort_flat_rings_A28_dim5_moon_low_cloud.json](../results/crm/comfort_flat_rings_A28_dim5_moon_low_cloud.json)).
CM1's nights at 0° E are fogged 85% of the time, against the GCM's 22% cloud,
and fog would keep the lowest air warm at night. The GCM makes cloud in its
lowest layer (the ground to about 1.8 km) only above 98.3% relative humidity.
The test lowered that threshold to 0.6 (the runner's `surface_cloud_rcrit`),
everything else as the design case, for 15 years from its last year. Cloud in
the lowest layer over tropical land rose from 4% to 32%: 36% in dark 3-day
means and 26% in sunlit ones, since the GCM's air there is barely more humid at
night than by day. The Moon cooled 5.0 K (294.9 to 289.9 K, settling near
289.0 K): cloud cover 23.5% to 43.2% and planetary albedo 0.280 to 0.309, the
clouds' effect on sunlight going from −8.5 to −16.0 W/m² and on the infrared
from +6.9 to +7.4 W/m², and the water vapour from 176 to 124 kg/m². Over the
tropical land the lowest layer cooled 4.9 °C in dark spells against 4.4 °C at
10 km, and 4.4 °C in sunlit spells against 4.45 °C. So more low cloud left the
GCM's lowest air at night slightly cooler against the air above, which runs
the wrong way to explain CM1's warm, muggy lowlands.

Both cloud tests leave the warmth of the lowest air against the air above within
half a degree of the design case's. In the GCM, cloud sets how warm the whole
Moon is; the gap with CM1 lies in how each model treats the lowest few
kilometres themselves: the surface's fluxes, the boundary layer's mixing and the
convection, which the GCM resolves with two layers, and the forcing CM1 takes
from the GCM (its air held from 8 km up, the vertical wind it is given). The
tests also show how much the GCM's climate rests on its cloud, fitted to Earth.
Settings that move its cloud cover between 20% and 43% move the Moon between
about 297 and 289 K. Its comfortable hours on the rings' land, 387 per lunar
day as designed, fall to 290 with the warmer Moon and to 121 with the cooler
one, whose air is too cool.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28 --set sunlight_scale=0.95 --set surface_cloud_rcrit=0.6 \
    --start-from climate/gcm/runs/A28_dim5_moon/model/MOST_REST.00029 --folder A28_dim5_moon_low_cloud --years 15
climate/gcm/.venv/bin/python -m climate.gcm.compare --levels A28_dim5_moon:20-29 A28_dim5_moon_low_cloud:5-14
climate/gcm/.venv/bin/python -m climate.crm.ring_comfort ring_equator ring_70_45e ring_70_135e --gcm A28_dim5_moon_low_cloud:5-14
```

## A box over the highlands

The author's plan for highlands (2026-09-29) left a three-dimensional box over
high ground to the flat rings' results, and the author chose to run it on
2026-10-01, with CM1 and the GCM still apart on the air near the ground. Case
`box_highland` is `box_0e`'s size, 64 by 64 columns 6.0 km apart (385 km
square) on 111 levels, centred where `ring_70_45e` crosses 44.7° S, 246.1° E on
the far side, its x axis along the ring (151° from north). The site was chosen
along the ring's path for high, varied land at mid-latitudes, where the flat
rings, corrected for height, put their comfort: the box's ground runs from the
floor of a small basin at sea level to 6.9 km, half of it between 1.3 and
5.2 km, with a mean of 3.1 km.

- **Ground.** Each column takes the atlas's height over its 0.25° neighbourhood
  ([atlas_heights](cm1_run.py)), smoothed by two passes of a 1-2-1 filter in x
  and y and blended toward the mean of the edge columns (2.9 km) over 50 km at
  each side, so the ground joins where the box wraps round
  ([box_heights](cm1_run.py)); the steepest slope is 0.18, as gentle as the
  terrain ring's. 3% of the columns hold water in the atlas, a basin's floor,
  laid as land at its level.
- **Surface.** All land at the site's wetness (0.06, class 0.05). Every column's
  ground and deep-soil temperature is carried from the ring's sea-level ground,
  293.7 K, at the GCM's lapse rate there, 0.46 °C per kilometre, as the terrain
  ring's were; the terrain build reads them from `terluna_surface2d.txt`, a
  patch to CM1's surface set-up for boxes over terrain.
- **Forcing.** The corrected GCM's upper air along `ring_70_45e`, the GCM's
  mean vertical wind at the site, and the heating and moistening of the ring's
  own day-night circulation there over its second lunar day. The terrain build
  places all three by the heights of the box's first column, a corner at 2.9 km
  on the blended edge, and its levels follow the ground. So the air is held from
  5.3 km above every column's ground (8.0 km above sea level at the corner),
  each level toward the reference at its own height, where the flat rings hold
  theirs from 8 km above the ground. The day-night forcing and the vertical wind
  reach every column 2.9 km higher above its ground than on the ring, and the
  ring's forcing in its lowest 2.9 km goes unused. The Coriolis force is the
  site's (f = −3.7 × 10⁻⁶ s⁻¹).

A first model day ran in 6.6 minutes on 8 threads at CM1's longest time step
(80 s), with the wind at 10 m under 1.3 m/s through the first night and the air
cooling 0.8 °C per kilometre of height by its end. The two lunar days run in
sessions of up to seven hours with an hour's rest between, about 6½ hours in
all on the faster build.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_terrain
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_highland
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_highland --hours 7 --threads 8    # restartable
```

**What the highland box shows** (its second lunar day, days 29.6–59.0:
[highland_analysis.py](highland_analysis.py),
[highland_box_highland.json](../results/crm/highland_box_highland.json)). The
run took 6.4 hours on 8 threads, its time step falling from 80 s to about 42 s
through the sunlit hours as convection grew. Over the last third of the second
lunar day its air ran 0.2 °C warmer and its dewpoint 0.1 °C lower than at the
same local times a lunar day before. The analysis takes the 2,116 columns beyond
the blended edges, from sea level to 6.9 km (mean 3.4 km): 768 valleys and
basins, 582 slopes and 766 ridges and summits.

These highlands are dry and sunny. Cloud covers 1% of the box (4% of the ground
above 6 km), no fog forms at night, and 0.2 mm of rain falls per lunar day
(1 mm above 6 km). Days are hot and nights mild: the sunlit air has a median of
30.0 °C (hottest 5%, 36.5 °C), the dark air 18.4 °C (coldest 5%, 15.3 °C), and
dewpoints are 15.5 °C by day and 13.9 °C at night. Of each lunar day 192 hours
are comfortable on the strict terms (308 on the loose), 262 too hot and 184 too
cold, and comfort rises with height:

| Ground | Columns | Sunlit air: median, hottest 5% | Dark air: coldest 5%, median | Dewpoint, day, night | Comfortable hours (flat ring corrected) |
|---|---|---|---|---|---|
| 0–1 km | 155 | 35.0, 38.4 °C | 15.3, 18.5 °C | 17.9, 14.5 °C | 134 (3) |
| 1–2 km | 450 | 33.6, 37.0 | 15.1, 18.4 | 17.2, 14.4 | 126 (71) |
| 2–3 km | 333 | 32.3, 35.7 | 16.1, 19.3 | 16.6, 14.6 | 160 (136) |
| 3–4 km | 273 | 30.7, 34.1 | 15.9, 19.1 | 15.8, 14.2 | 228 (170) |
| 4–5 km | 481 | 29.5, 32.8 | 15.3, 18.4 | 14.9, 13.5 | 243 (194) |
| 5–6 km | 297 | 27.9, 31.2 | 15.2, 17.8 | 14.2, 12.9 | 240 (175) |
| 6–7 km | 127 | 26.6, 29.9 | 14.5, 16.8 | 13.5, 12.2 | 199 (179) |

Valleys get 160 comfortable hours, slopes 194 and ridges 223. By day the air
cools 1.1 °C per kilometre of height and its dewpoint falls 0.6 °C; at night
only 0.35 and 0.4 °C, because cold air drains downhill. Through the night the
wind at 10 m blows down the slopes, 0.2–0.7 m/s and strongest just after
sunset, and valleys run 0.3 °C colder than ridges beyond what their heights
give. Through the sunlit hours it blows uphill at 1.2–1.7 m/s, 81–98% of the
time, the two switching about six Earth days before and after noon. What little
rain falls is three times heavier on slopes the wind blows up than on those it
blows down (0.011 against 0.004 mm a day).

Against the plan's height correction the box agrees on high ground, 228–243
comfortable hours at 3–6 km against the corrected flat ring's 170–194, and finds
far more comfort low down, 126–134 hours below 2 km against 3–71, since its
nights are cooler and drier than the flat ring's muggy ones. The twin's 1 °C per
kilometre holds by day; at night the air changes less with height. Against the
GCM's four land cells near the site (at 1.9–3.9 km: air 20.3 °C, lowest-layer
dewpoint 11.6 °C, 602 comfortable hours, 16 mm of rain per lunar day) the box
keeps CM1's warm, humid offset and rains almost nothing. Two of its settings
bear on these numbers. Its air is held from 5.3 km above the ground where the
flat rings' is held from 8 km, and every column takes the day-night forcing at
the corner's heights, which may keep its low ground drier and more comfortable
than a like-for-like ring would; the rerun below, with each column held and
forced at its own height, finds that it did, through the day-night forcing. Its
deep ground follows the GCM's 0.46 °C per kilometre, gentler than the air's
change by day, so high ground may run slightly warm.

```sh
climate/gcm/.venv/bin/python -m climate.crm.highland_analysis box_highland --from-day 29.5
```

**The box again, each column at its own height** (case
`box_highland_own_height`, the author's go-ahead of 2026-10-02). The terrain
build now holds and forces every column at its own height above sea level. The
hold to the GCM's air weighs each column's own height (from 8 km above sea
level, fully from 16 km): at each level it averages the held columns'
departures from the reference at their heights and applies that to each column
with its weight. The ring's day-night heating and moistening and the GCM's
vertical wind (a new mode, `var13 = 3`, staggered once to u and v points) are
interpolated to each column's own heights. The inputs are the first run's,
byte for byte. At start-up the model prints where each forcing lands, and the
first model day showed:

| Column | Ground | Held from | Vertical wind at its lowest w level | Day-night forcing from |
|---|---|---|---|---|
| Lowest | 0 m | 8.6 km above sea level | 0.015 mm/s | the ring's lowest level, 50 m |
| Corner | 2.9 km | 8.0 km | 0.38 mm/s | the ring's level at 2.9 km |
| Highest | 6.9 km | 8.1 km | 0.46 mm/s | the ring's level at 7.0 km |

The precomputed weights and threaded loops keep the cost near the first run's
(188 model seconds per wall second on the first day against 217). The stopped
terrain ring would now be held this way too if it were resumed, from 8 km above
sea level over every column, where it was held from 8 km above each column's
ground.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_terrain
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_highland_own_height
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_highland_own_height --hours 7 --threads 8    # restartable
```

**What the box shows at its own heights** (its second lunar day, days
29.6–59.0: [highland_analysis.py](highland_analysis.py),
[highland_box_highland_own_height.json](../results/crm/highland_box_highland_own_height.json);
both runs' summaries now carry their purpose and executable, and the first
run's states its departure). The run took 7.7 hours on 8 threads in two
sessions, at 185 model seconds per wall second (the first run took 6.4 hours at
220). It settled as the first did: over the last third of its second lunar day
its air ran 0.3 °C warmer and its dewpoint 0.2 °C lower than at the same local
times a lunar day before.

Holding and forcing each column at its own height changed the low ground and
little else. The air at 2 m stays within 0.4 °C of the first run's in every
height band, by day and by night (over the box a sunlit median of 30.1 °C and a
dark one of 18.5 °C), and so do the cloud (1% of the box) and the air's change
with height by day (1.15 °C cooler per kilometre). The low ground is more
humid, because the ring's day-night forcing now reaches it at its own height.
The ring's circulation cools its lowest kilometre at night by 0.04–0.18 °C a
day and moistens it by 0.07–0.28 g/kg a day, and in sunlight moistens it by
0.15–0.33 g/kg a day; at 2.9 km, where the first run took the low ground's
forcing from, it warms the air at night by 0.18 °C a day and moistens it by
0.05 g/kg. Over the low ground the air 0.5–1 km above sea level now runs
0.5–0.65 °C cooler and about 0.45 g/kg moister at night than in the first run,
and 0.5 g/kg moister by day. The GCM's vertical wind, which the first run also
took from 2.9 km, works the other way and more weakly: its ascent cooled that
air by about 0.09 °C a day at night. After sunrise the low ground's dewpoint at 2 m
climbs to 16.4 °C, where the first run's stayed near 13.4 °C, so the morning
passes through comfortable warmth in muggy air, and comfort falls low down:

| Ground | Comfortable hours, first run → own heights (flat ring corrected) | Too humid | Too cold | Rain, mm per lunar day |
|---|---|---|---|---|
| 0–1 km | 134 → 82 (3) | 92 → 133 | 182 → 197 | 0.01 → 0.02 |
| 1–2 km | 126 → 101 (71) | 103 → 126 | 189 → 193 | 0.01 → 0.03 |
| 2–3 km | 160 → 145 (136) | 131 → 145 | 134 → 133 | 0.02 → 0.03 |
| 3–4 km | 228 → 216 (170) | 74 → 83 | 142 → 138 | 0.10 → 0.06 |
| 4–5 km | 243 → 244 (194) | 37 → 39 | 181 → 176 | 0.24 → 0.12 |
| 5–6 km | 240 → 262 (175) | 18 → 17 | 226 → 204 | 0.60 → 0.19 |
| 6–7 km | 199 → 205 (179) | 9 → 11 | 305 → 296 | 0.99 → 0.26 |

Over the box 183 hours of each lunar day are comfortable on the strict terms
(192 in the first run) and 309 on the loose (308). Valleys and basins lose the
most, 134 hours against 160, while slopes and ridges keep theirs (194 and 223),
so comfort now climbs more steeply with height, from 82 hours below 1 km to 262
at 5–6 km. By day the dewpoint falls 0.8 °C per kilometre of height (0.6 in the
first run). At night the air cools 0.25 °C per kilometre (0.35), and valleys run
0.17 °C colder than ridges beyond what their heights give (0.28). Low cloud
gathers around sunrise over up to 17% of the low ground, with a trace of
drizzle. On the first lunar day, from the starting state, it grew into a deck
over 79% of the low ground, its base 150–650 m and its top up to 1.7 km above
the ground, which drizzled at up to 0.4 mm a day and cleared within three Earth
days of sunrise; the first run had only thin fog then, gone within an Earth day
of sunrise. From late morning the upslope wind is weaker, 1.4 m/s against
1.7 before noon and 0.8–0.9 against 1.2–1.3 in the afternoon, and it turns
downhill before sunset; just after sunset the downslope wind is stronger,
0.9 m/s against 0.7. Rain falls to 0.09 mm per lunar day (0.2), still most on
the highest ground.

Above 2 km the two runs' air agrees within 0.15 °C. Over the second lunar day,
moving the hold to 8 km above sea level changed the moisture of the air above
5 km by at most 0.35 g/kg: drier at 6–8 km over the low ground at night, where
the first run held it and this one does not, and moister at 7–9 km over the
high ground, weakly held now. So the first run's error reached the ground through the
day-night forcing over the low ground.

Against the plan's height correction the box now agrees at 1–3 km (101 and 145
hours against the corrected flat ring's 71 and 136) and still finds more
comfort in the basin below 1 km (82 against 3) and on high ground (216–262 hours
at 3–6 km against 170–194). The GCM's four land cells near the site (602
comfortable hours, 16 mm of rain per lunar day) stay far from both: the box
keeps CM1's warm, humid offset and rains almost nothing. Its deep ground still
follows the GCM's 0.46 °C per kilometre.

```sh
climate/gcm/.venv/bin/python -m climate.crm.highland_analysis box_highland_own_height --from-day 29.5
```

## What the storms hold for charging

The atmospheric-electricity study reads the storms of six runs over their second
lunar day with [mixed_phase_analysis.py](mixed_phase_analysis.py): the two
equatorial boxes, `box_0e_small`'s twin with Earth's fall speeds, and the three
flat rings. It takes Morrison's size distributions and fall-speed laws as the
runs used them and finds the charging zone, the air from 0 to −40 °C where
graupel, cloud ice or snow, and supercooled cloud water each reach 0.01 g/kg
(0.001 and 0.1 g/kg as checks):

| Run | Charging zone, m of air per m² of ground (snapshots with it) | Graupel aloft, kg/m², and its residence | Graupel in the zone, g/m³, mass-weighted diameter | Cloud ice in the zone, per litre | Updrafts in the zone, median / top tenth |
|---|---|---|---|---|---|
| `box_0e` | 22 (34 %) | 0.30, 4.9 h | 1.8–2.5, 5.4–6.9 mm | 3–7 | 3.1 / 10.3 m/s |
| `box_0e_small` | 22 (31 %) | 0.23, 4.4 h | 1.7–2.3, 5.0–6.9 mm | 2–10 | 3.4 / 9.9 m/s |
| `box_0e_small_earth_fall` | 11 (21 %) | 0.04, 2.2 h | 0.7–1.3, 5.3–6.2 mm | 3–7 | 4.5 / 12.5 m/s |
| `ring_equator` | 50 (98 %) | 0.10, 3.2 h | 0.5–0.9, 6.3–6.9 mm | 0.6–13 | 0.9 / 3.5 m/s |
| `ring_70_45e` | 3.6 (22 %) | 0.007, 2.8 h | 0.6–0.9, 6.8–7.0 mm | 0.2–2.8 | 1.0 / 3.6 m/s |
| `ring_70_135e` | 8.4 (35 %) | 0.017, 3.2 h | 0.5–1.1, 6.7–7.0 mm | 0.4–3.6 | 0.8 / 3.2 m/s |

In every run the 0 °C level stands near 25 km and −40 °C near 60 km, five to six
times their heights in Earth's tropics, and cloud tops reach 55–80 km. The
charging zone is shallow and warm: four-fifths of it lies between 0 and −10 °C
(30–32 km, 650–680 hPa), a fifth between −10 and −20 °C, almost none colder,
because the supercooled cloud water thins from 0.10–0.17 g/m³ near 0 °C to
0.03–0.05 g/m³ at −20 to −30 °C. Graupel is plentiful and large; a quarter of it
between 0 and −10 °C sits at the scheme's largest size (a mass-weighted diameter
of 8 mm), where the scheme limits growth. Cloud ice is scarce, and snow (0.5–1
g/m³) carries most of the ice that graupel meets. Graupel strikes cloud ice at
1.3–1.5 m/s and snow at 0.8–1.0 m/s. The zone forms from late morning to the
end of the afternoon, and in the boxes its patches are about 10 km wide, near the
grid's limit, and a median 6 km deep.

With Earth's fall speeds the same box holds a sixth as much graupel aloft, half
the charging zone and half the residence, and its graupel strikes ice at
2.9–3.3 m/s. The analysis also integrates graupel's collisions with ice and snow
over their size distributions under four laboratory charging laws, at the runs'
fall speeds and at Earth's; the study weighs the result against the air's
conductivity
([research/studies/atmospheric_electricity](../../research/studies/atmospheric_electricity/README.md)).

```sh
climate/gcm/.venv/bin/python -m climate.crm.mixed_phase_analysis box_0e --from-day 29.5
```

## Electrified storms

[cm1_elec.py](cm1_elec.py) builds CM1 with electrified storms, stage 2 of the
atmospheric-electricity study. WRF-ELEC's version of the NSSL two-moment
microphysics (MicroTed/wrf4-elec at commit e43041b, in the public domain under
the WRF notice; Mansell et al. 2005, 2010; Fierro et al. 2013) replaces CM1's
own copy of the same scheme. It carries charge on cloud droplets, rain, cloud
ice, snow, graupel and hail through every microphysical process. Charge
separates when graupel or hail rebounds from cloud ice or snow (non-inductive
charging; by default Saunders and Peck's law with Brooks et al.'s critical rime
accretion rate above −15 °C and none below −32.5 °C, WRF-ELEC's `isaund = 12`)
and when cloud droplets rebound from graupel polarised by the field (inductive
charging). The small ions' net charge attaches to the particles. Four Terluna
files in [fortran/](fortran/) and two of WRF-ELEC's complete it:

- [terluna_elec.F](fortran/terluna_elec.F) passes CM1's calls of the scheme to
  WRF-ELEC's driver with the arrays reordered from CM1's (i, j, k) to WRF's
  (i, k, j), and keeps the charges in CM1's passive tracers (C/kg), which CM1
  advects and mixes like any tracer, without its positivity limiter. It runs
  the steps as WRF-ELEC's own driver does: the microphysics once, without
  sedimentation, then the sedimentation in sub-steps of about 0.75 s, each
  followed by a solve for the electric field of the net charge and by
  lightning. Leakage, where a case asks, follows once a step. It keeps the
  vertical field for the next step's inductive charging, also across restarts
  (`terluna_ez_<time>.bin` beside CM1's restart files, as WRF keeps it in its
  own).
- [terluna_lightning.F](fortran/terluna_lightning.F) solves Poisson's equation
  for the potential by Fourier transform across the periodic domain and a
  tridiagonal solve on CM1's stretched levels, with the potential zero at the
  ground and at the model top. It carries WRF-ELEC's cylindrical discharge
  scheme (`light1d`; Ziegler and MacGorman 1994): where the field exceeds the
  breakdown field (284 kV/m at 1.225 kg/m³, scaled by air density and held
  within 50–180 kV/m, WRF-ELEC's `nssl_ibrkd = 4`), every column within 12 km
  that holds more than 0.1 nC/m³ joins the discharge, which removes shares of
  the charge above that level, set so that it carries comparable positive and
  negative charge (about 30 % of it), and hands them to the small ions with the
  opposite sign.
- WRF-ELEC's branched lightning (`lightmsz` in `module_discharge_msz.F`;
  MacGorman, Straka and Ziegler 2001), the builds' default, makes one flash at a
  time. It starts at a point where the field reaches 0.9 of the breakdown field
  in air colder than 275 K, drawn at random among such points, runs a channel up
  and down from it while the field stays above 0.8 of breakdown, and spreads
  each end through the connected cells whose potential differs from the
  channel's by more than its 200-V/m internal field allows. It neutralizes the
  charge in those cells, handing the opposite charge to the small ions. A
  downward channel strikes the ground when it reaches air warmer than −7 °C
  (about 5 km above Earth's ground), if the flash started where the potential
  lies beyond 50 MV of either sign with a vertical field above 10 kV/m in the
  matching direction, and the channel's end lies in charge of the matching
  sign. WRF-ELEC applies this rule by default (its ground height of −1); a
  height above zero replaces the temperature with that height, and a height of
  exactly zero leaves only a channel reaching the two lowest levels, as the Earth
  benchmarks below first ran (found 2026-10-05). A ground flash neutralizes
  charge of one sign in the cloud, the ground supplying the opposite charge.
  Each flash makes nitrogen oxides
  along its channel at Wang et al.'s (1998) yield per metre at the channel's
  pressure, times WRF-ELEC's factor of 0.1, reduced in proportion where the
  flash changes the charge by less than 1 nC/kg. After each flash the field is
  solved again, until no starting point is left or a chosen one gives no flash;
  after 120 flashes in a sub-step, 5 % of all charge goes to the small ions, as
  in WRF-ELEC. Its trilinear interpolation (`mlint2`) comes from WRF-ELEC's
  `module_boxmgsetup.F`.
- [terluna_branched.F](fortran/terluna_branched.F) runs `lightmsz` on one
  process; WRF-ELEC runs it under MPI. It draws the starting point as `lightmsz`
  does, with the same Park and Miller generator, seeded from the step so that a
  restarted run repeats its flashes, recentres the arrays on that point across
  the periodic domain (the channels stop at the arrays' edges), hands
  `lightmsz` that point alone and returns the charge, channels and nitrogen
  oxides to the domain. The patch to `lightmsz` replaces its few MPI calls
  outside its MPI blocks, takes the breakdown field's bounds from the build,
  recomputes that field at each call and sends its report to its own log.
- [terluna_screen.F](fortran/terluna_screen.F) ports WRF-ELEC's screening
  layers (`screen`; Ziegler et al. 1991), off unless a case asks for them.
  Clear air conducts better than cloud, so where cloud meets clear air the
  current the field drives into the edge does not match the current leaving it,
  and charge gathers there. Once a step, at each cloud point with cloud on one
  side and clear air beyond it on the other (along x, y or z), it takes the
  field normal to the edge inside and outside and gives the small ions the
  charge the two conduction currents leave in a step, capped at 0.25 nC/m³ of
  net charge there. The conductivity is WRF-ELEC's table of Earth's clear air
  with a tenth of it in cloud, or the clear air and cloud of
  [atmosphere/electricity](../../atmosphere/electricity/README.md) by height.
  Its vertical differences take the distances between CM1's scalar levels,
  where WRF-ELEC takes the depth of the layer below; the two agree on an even
  grid.
- The same file holds point discharge from the ground (`ground_discharge`),
  which WRF-ELEC does not have, off unless a case asks for it (`corona_v_m`,
  `var2`). Under a storm, plants and other grounded points go into corona once
  the field at the ground passes an onset, about 3 kV/m over dense vegetation
  and 5 kV/m on a barren ridge, and the ions they give off hold the field at the
  ground near 5–12 kV/m while it stands several times stronger a few hundred
  metres up (Standler and Winn 1979; Soula and Chauzy 1991). Once a step, over
  land, the field at the lowest level drives a current J = a(E² − E₀²) up from
  the ground, with a set by Standler and Winn's 1 nA/m² at 8 kV/m over a
  3-kV/m onset and E₀ scaled by the density of the air at the ground. The charge
  takes the sign of the ground's induced charge, joins the small ions' charge in
  the lowest layer and moves with the air from there; the small ions do not
  drift in the field in this scheme. A step gives off at most the charge that
  brings the field at the ground back to the onset, so where a storm charges
  faster than the points discharge, the field there settles near the onset, as
  measured. WRF-ELEC hands the small ions' charge in a cell at once to any
  cloud or precipitation particles there, so where rain or cloud reaches the
  lowest level the ground's ions go into it, and rain carries them back down.
  In a 10-kV/m field raindrops would take up such ions within about a minute,
  while the ions climb about 100 m, so the two are alike there; the scheme's
  ions only rise on the air's motion, so it may underplay screening in calm air.
  Every electrified run logs the field at the ground every ten steps
  (`terluna_ground.txt`): its largest value over land, and where no particles
  reach the lowest level, with the charge the ground gave off and the columns
  discharging, wet and dry. Without point discharge the field at the ground
  under the boxes' storms reaches 100–220 kV/m. In a 30-minute check from the
  coarse box's day-10.5 restart the ground gave off about 500 C every eight
  minutes, mostly negative, from 35–41 of the 4,096 columns at a time, at a
  mean field there of about 7–8 kV/m by the charge given off. Of those columns
  33–38 had cloud or rain at the lowest level, and there the strongest field at
  the ground stayed at 35–40 kV/m; where no particles reached the lowest level
  it was 3.6–3.8 kV/m, at the onset. Under rain this scheme screens the ground
  less than point discharge would, since its ions neither climb in the field
  nor wait to be taken up.

Lunar gravity enters the scheme as it enters Morrison's: the module's gravity is
the build's, the drag-law fall speeds (graupel and hail under drag laws, cloud
droplets under Stokes's law) follow it directly, and the fitted fall-speed laws
V = a D^b scale by (g/9.81)^((b+1)/3): rain with b = 0.8, snow with 0.42, cloud
ice with 0.55, and graupel and hail with each particle's exponent in WRF-ELEC's
default table (Milbrandt and Morrison 2013), whose coefficient also sets their
collection and ventilation. The scheme sets a particle's fall speeds and
coefficients only where that particle is present, and the scaling acts only
there; at Earth's 9.81 m/s² it leaves them as they are. WRF-ELEC's driver and
its sedimentation driver run in one thread, since WRF divides its work outside
the microphysics; the build gives each slab loop the OpenMP directive CM1's
copy has, with every per-slab array private, the charge totals summed across
threads, and what the driver sets before the loop shared.

The lunar boxes change two of WRF-ELEC's Earth settings (`cm1_elec.LUNAR`).
WRF-ELEC's 0.75-s sub-step lets graupel settle through about 1 % of the 500-m
layers of its storms in each; in the boxes' charging zone, 25–35 km up, the
layers are 2 km deep and graupel falls at 0.44 of Earth's speed, so a 6.8-s
sub-step keeps that share. On the Moon air warmer than −7 °C lies about 27–34
km up, so a downward channel strikes the ground when it comes within 5 km of it,
the height that rule stands at on Earth. The 50-MV and 10-kV/m conditions for a
ground strike, the 0.9 and 0.8 fractions of breakdown, the 200-V/m internal
field and the nitrogen oxide yield stay as WRF-ELEC calibrates them for Earth.
That height rule keeps WRF-ELEC's demand that the channel's end lie in charge of
the matching sign, and the lunar storms hold no charge 5 km up, so under it no
lunar flash could strike the ground; see "Ground strikes" below.

WRF-ELEC also caps the breakdown field at 180 kV/m. On Earth the cap applies
below about 4.5 km, beneath where flashes start. The lunar air is dense enough
for it to apply below about 36 km: there it sets the breakdown field below the
density scaling by 4 % at 34 km, 10 % at 30 km, 23 % at 20 km and 44 % at the
ground (180 against 322 kV/m). In `box_0e_elec`'s first lunar day 84 of its 162
flashes started at 32–34 km, where the cap lowered the field a flash starts at
by 4–7 %, and 152 ran channels below 36 km, where it lowers the field a channel
stops at by 10–19 % at their lower ends (24–30 km). The lunar boxes therefore
lift the cap (`lightning` 4; the author's decision, 2026-10-04), keeping the
density scaling at every height. WRF-ELEC's 50-kV/m floor stays: it applies
only above about 100 km in the lunar air, far above the clouds, near the model
top where the potential is held at zero.

Three faults in the first version were found and fixed. It discharged once a
6-second step, after the step's charging and all its sedimentation; WRF-ELEC
sub-steps the sedimentation and discharges after each sub-step, and without
that the field stood a median 3.2 times above breakdown at a step's first
discharge and each discharge removed a median 150 C. Two faults showed up as
results that changed with the number of threads. The scaling ran wherever
graupel might be, also where the scheme had set no graupel exponent, so it
multiplied by whatever value was left in memory there; and the directive made
private the fall-speed moments the driver sets before the loop (`infdo`),
leaving each thread's copy unset. A test keeps every private variable of both
directives from being set before their loops, and a build that fills unset
reals with signalling NaNs and stops at the first invalid operation located the
scaling fault. Restarted from its restart at 30 minutes, on 4 threads and on 2,
the benchmark storm below repeats its first run's output fields at 35 and 40
minutes and its 41 flashes in between bit for bit.

A case sets the electricity in CM1's namelist: `var3` the screening layers (0
none, 1 with WRF-ELEC's Earth conductivity, 2 with the conductivity file),
`var4` the height a downward
channel must reach to strike the ground (m; 0 for WRF-ELEC's −7 °C rule),
`var5` the sub-step (s; 0 for WRF-ELEC's 0.75 s), `var6` the charging (2
non-inductive, 3 with inductive), `var7` the charging law (12 Saunders and
Peck, cut off below −32.5 °C; 11 the same law as WRF-ELEC's test case sets it,
with a smooth critical rime accretion rate and charging below −32.5 °C; 1
Takahashi), `var8` the lightning (0 none, 1 cylinders, 3 branched; 2
and 4 the same with the breakdown field's scaling unbounded), `var9` leakage
through the conductivity of
[atmosphere/electricity](../../atmosphere/electricity/README.md) and `var10`
the cylinders' radius. Each run logs every flash or cylinder
call (`terluna_flashes.txt`: its kind, start, the field and breakdown field
there and the domain's largest field over breakdown before it, the area and the
levels where it neutralizes charge, the positive and negative charge it
neutralizes, the electrostatic energy before and after, its channel points and
nitrogen oxides),
`lightmsz`'s own report on the first 2,000 flashes of each run
(`terluna_msz.log`), and, every ten steps, the largest field, the domain's
charge and energy and the largest charging rates (`terluna_field.txt`);
WRF-ELEC's driver prints the domain's charging by collision pair every step.
[elec_analysis.py](elec_analysis.py) reads the logs and the charge tracers of
the snapshots.

The build leaves out CM1's own water budget from the NSSL scheme (the
condensation, evaporation and rain totals CM1's copy adds to its budget
output), the three-moment option, the activated CCN and IN arrays, terrain under
the field solver and `lightmsz`'s horizontal steps of the starting channel, which
WRF-ELEC leaves off as well. Charge leaves the air by lightning to the ground, by
falling to the ground on rain and hail, and, where a case turns it on, by
leakage, which relaxes the net charge at σ/ε₀ with the conductivity of clear air
or of cloud by height; the screening layers that conductivity gradients build at
cloud edges come with `var3`.

The solver, the lightning and the screening are checked against numpy and
against charge laid out by hand
([test_cm1_elec_field.py](../tests/test_cm1_elec_field.py),
[test_cm1_elec_branched.py](../tests/test_cm1_elec_branched.py),
[test_cm1_elec_screen.py](../tests/test_cm1_elec_screen.py)). The
potential matches a direct solve of the same difference equations on stretched
levels to one part in 10⁵, and the exact potential of a charged layer between
grounded plates to one part in 10³. A cylinder call removes the charge
WRF-ELEC's rule sets, only within the cylinders, and counts separate regions
across the domain's edges. A branched flash between two charge regions
neutralizes equal positive and negative charge in both and makes nitrogen
oxides; a flash at the domain's edge reaches both sides; flashes continue until
no starting point is left; a channel that runs down through a strong lower
positive region strikes the ground with negative charge, the same storm with
every charge reversed strikes with the same amount of positive charge, and
moving the ground-strike height from 3 to 5 km turns that flash from one in
cloud to one to ground. In an upward field a block of cloud gathers negative
screening charge on its top and positive on its base, each as the formula
gives it, none inside it, in clear air or on sides the field runs along, never
past the cap, and edges across the domain's wrap gather it too.

### The benchmark storm

`supercell_elec` is CM1's own supercell (Weisman and Rotunno 2000: the Weisman
and Klemp sounding with 14 g/kg of vapour at the ground, the quarter-circle
hodograph, a warm bubble, 1-km spacing, a domain 120 km square, two hours) at
Earth's gravity, with 600 CCN per cm³ and WRF-ELEC's defaults: its charging,
its 0.75-s sub-steps and its branched lightning
([elec_supercell_elec.json](../results/crm/elec_supercell_elec.json)).
`supercell_elec_cylinders` runs the same storm with the cylinders
([elec_supercell_elec_cylinders.json](../results/crm/elec_supercell_elec_cylinders.json)),
and `supercell_nssl` runs it with CM1's own copy of the NSSL scheme without
electricity. The electrified storm and CM1's own agree. At every half hour the
largest updraft is within 6 % (28, 46, 53 and 56 m/s against 28, 46, 56 and
56), reflectivity within 3.4 dBZ, and the domain's condensate within 9 % and at
each level holding much of it within 19 %; from 30 minutes on, graupel is
within 12 % and cloud-ice numbers within a factor of 1.5 (reflectivity once
differs by 7 dBZ, at 35 minutes), and the rain at the end is within a third.
Hail differs most, up to a factor of two once it
forms, as the two copies of the scheme are different versions. The reordering of
the arrays between CM1 and WRF-ELEC and the sub-stepped sedimentation therefore
leave the storm as CM1's own scheme makes it. The electricity, with its
sub-steps and flashes, triples the run's time: 56 minutes on 4 threads, beside
the lunar box on the other 4, against 18 for `supercell_nssl`.

No published run puts this scheme in CM1 or in WRF's idealized supercell, so the
benchmark sets the storm beside published runs of the same charging and
lightning schemes on similar storms (notes and pages in the study's
[sources](../../research/studies/atmospheric_electricity/sources.json)):

| | `supercell_elec` | Published runs |
|---|---|---|
| First flash | 28.5 min | 32.5 min (MacGorman et al. 2001, the same branched scheme), 27–28 min (Kuhlman et al. 2006), about 30 min (Ziegler et al. 2014), 30–40 min (Sun et al. 2023) |
| Flash rate | rising to 153 a minute at 75–80 min, about 110 a minute through the second hour | 34 a minute at the peak (Mansell 2014), 75 (Brothers et al. 2018, 1 km), 100–200 with peaks of 400–500 (Ziegler et al. 2014), 300–600 in a mature storm (Calhoun et al. 2014) |
| Kinds | all 7,134 in cloud, with the ground rule off (above); with it, 7,321 in cloud and 85 negative to ground (1.1 %), the first at 36.8 min | 99.67 % in cloud (Salinas et al. 2021, of Brothers et al.'s storm); a first ground flash at 34 min (Ziegler et al. 2014) |
| Charge per flash | median 9.4 C; 2.8–67 C from the 10th to the 90th percentile | 2.7–39.6 C (MacGorman et al. 2001) |
| Flash area | median 27 km², mean 79 km² | mean 160 km² at 1 km (Brothers et al. 2018, with a different branched scheme) |
| Net charge density | up to +4.7 and −4.8 nC/m³ | about 1 nC/m³ (Ziegler and MacGorman 1994), up to 4 (Sun et al. 2021) |
| Charge on one particle type | up to 45 (graupel), 48 (cloud ice) and 53 (hail) nC/m³, largely cancelling | up to 1.2 nC/m³, in a weaker storm (Sun et al. 2024) |
| Largest non-inductive charging | 3.2 nC/m³/s | 0.3–0.55 nC/m³/s (Kuhlman et al. 2006; Sun et al. 2023) |
| Inductive against non-inductive charging | a sixth of the domain's charge separated | about a tenth (two-moment runs) |

The storm electrifies and flashes when the published storms do, with their net
charge densities, flash rates and charge per flash. Its charging is several
times stronger and the charge on each particle type tens of times larger. Both
follow from its cloud ice: up to 2×10⁸ crystals per kilogram where droplets
freeze homogeneously in updrafts of 50 m/s, which CM1's own copy of the scheme
makes as well, and each crystal that bounces off graupel carries charge away.

At 30–35 minutes the storm's net charge lies as a normal tripole's does:
negative at 7.75 km (−24 to −25 °C) between positive above (8.25–8.75 km, −29
to −33 °C) and below (6.25 km), 14–38 C at the negative level. Snow, cloud ice
and cloud droplets carry the negative charge at 30 minutes, and cloud ice,
cloud droplets and graupel at 35. From 40 minutes, as the updrafts pass 50 m/s,
the main negative charge rises to 10.75–11.75 km (−47 to −55 °C), and from 75
minutes it lies mostly at 8.75–9.75 km (−33 to −40 °C), up to 340 C at one
level. Snow and cloud ice carry it, snow alone from 90 minutes, while graupel
and hail at the same levels carry positive charge. Saunders and
Peck's law charges graupel positively where it rimes fast in strong, wet
updrafts, and the published runs with this law report the same inverted charge
in such storms (Kuhlman et al. 2006; Mansell et al. 2005), as high CCN does (Sun
et al. 2023, 2024).

A flash starts a median 8.75 km up (5.75–9.75 km), where the field reaches a
median 107 kV/m against a breakdown field of 109 kV/m; it neutralizes charge
between a median 4.75 and 10.75 km over a median 27 km², and dissipates a
median 1.2 GJ (0.34–6.9 GJ) of electrostatic energy. A few span much of the
storm: 169 flashes neutralized more than 200 C each, and the largest 1,490 C
over 4,270 km². Before a sub-step's first
flash the domain's largest field stands a median 1.20 times the breakdown field
(1.03–1.59), before later flashes 1.12, and once 4.4. WRF-ELEC's sub-steps hold
it near breakdown: without them the first version's field stood a median 3.2
times above breakdown at a step's first discharge. Within one sub-step at most
12 flashes follow one another, against WRF-ELEC's cap of 120. The flashes make a
median 4.7 mol of nitrogen oxides each and 109,000 mol over the two hours, at
WRF-ELEC's yield along the model's 1-km channels; how that compares with the
yields observed per flash is open.

With the cylinders (`supercell_elec_cylinders`) the sub-steps also hold the
field near breakdown: a median 1.15 times it at a sub-step's first discharge
(1.02–1.35, at most 1.5), against 3.2 times in the first version, which
discharged once a 6-second step. A call removes a median 32 C of each sign
(14–124 C), against 150 C before, from 12-km cylinders reaching from near the
ground to 15 km, and dissipates a median 3.0 GJ. The 4,026 calls in 92 minutes,
up to 86 a minute, take 5,520 separate regions. Each call acts at once in every
cylinder around every point above breakdown, through the whole column, so its
charge and energy stand for several flashes together, and its counts follow the
trends of a storm's lightning (Fierro et al. 2013).

At the fine box's 2-km spacing (`supercell_elec_2km`, the same storm on a grid
twice as coarse, [elec_supercell_elec_2km.json](../results/crm/elec_supercell_elec_2km.json))
the scheme keeps its onset and rates: the first flash at 30 minutes, up to 125
flashes a minute and 90 a minute through the second hour (153 and 111 at 1
km), a median 15 C per flash over a median 44 km² (9.4 C and 27 km² at 1 km).
The field builds further past breakdown in places. Before a sub-step's first
flash it stands a median 1.22 times breakdown, as at 1 km, but 2.2 times at the
90th percentile (1.6 at 1 km), and it reached 680 kV/m at 9 km once (255 kV/m
at 1 km). Each grid point starts at most one flash a step, as in WRF-ELEC, and a
coarser grid holds fewer points in a charged volume; whether that holds the
flashes back there is not yet checked.

In WRF-ELEC's own test settings (`supercell_elec_wrf`, [elec_supercell_elec_wrf.json](../results/crm/elec_supercell_elec_wrf.json):
its 2-km grid over 84 km, 800 CCN per cm³, Saunders and Peck's law as its test
case sets it and the screening layers, run for two hours where WRF-ELEC's test
runs one) the storm flashes about twice as often, up to 298 a minute, a median 18
C per flash; its net charge reaches −19.6 nC/m³ and its field 620 kV/m. The
screening layers add 2–20 C a step at cloud edges. No output of WRF-ELEC's own
run of this case is at hand to set beside it. Neither 2-km storm strikes the
ground in two hours, but both ran with the ground rule off (above).

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build earth_g_omp_elec moon_omp_elec earth_g_omp
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup supercell_elec    # and supercell_elec_cylinders, supercell_nssl
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run supercell_elec --threads 4    # about an hour
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis supercell_elec --against supercell_nssl
```

### The electrified equatorial box

`box_0e_elec` runs `box_0e`'s inputs as they were written (its site, surface,
grid, starting air and forcing, from the GCM before its correction) with the
NSSL microphysics with hail in place of Morrison's, its CCN at `box_0e`'s 100
droplets per cm³, and WRF-ELEC's charging and branched lightning with the lunar
settings above, over two lunar days with three-hourly output and twelve-hourly
restarts. It runs from 2026-10-04. A first start, before the sub-steps and the
branched lightning, was stopped at day 4, before its storms had charged. Each
restart keeps its vertical field (`terluna_ez_<time>.bin`), so windows of the
storms can run again from it with frequent output.

**Ice sizes in the radiation, a note for the radiation and cloud-optics work.**
CM1's radiation (RRTMG) takes cloud ice sizes up to 140 µm, where its ice optics
table ends, and stops beyond. Lunar ice falls slowly and grows larger: the box
stopped at day 6.45, as its storms began, on ice of 157 µm. The electrified
build caps the ice size the radiation reads at 140 µm, as CM1's Morrison scheme
does in every earlier lunar run, and every ten steps `terluna_ice_optics.txt`
logs the ice mass above the cap, its mean size and the largest size where ice
reaches 0.001 g/kg. The capped ice reaches the radiation as optically thicker
than its size gives, by its size over 140 µm: anvils of large crystals reflect
too much sunlight by day and hold too much heat by night. Snow above 130 µm RRTMG
handles itself, cutting its mass by (130 µm / size)². The microphysics and the
electricity use the scheme's own sizes and do not see the cap. The box resumed
from its day-6 restart with this build.

From day 18, in its first lunar night, the box runs with the breakdown field's
cap lifted (`lightning` 4, executable `d2643310f7275844`), so its second lunar
day's storms have the density scaling at every height. A flash is the only
thing the cap changes, and none strikes between the first lunar day's last
(day 14.3) and the second's storms. Its first lunar day ran under the cap; those
storms run again without it from the day-10.5 restart (below). The box has no
point discharge from the ground on either day.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec --threads 8 --hours 30    # about half a day
```

### The fine box

`box_0e_elec_fine` is `box_0e_elec` at a third of its spacing, 2.0 km, over a
box a third as wide, 128 km, with the same 64 by 64 columns. CM1 cannot move a
run onto a finer grid, so it starts from the coarse run's output at a chosen day:
the air averaged over the coarse box's columns, by level (potential
temperature, vapour and wind), with the surface pressure, the air at 2 m and the
skin temperature of the land (the soil's layers between it and the deep
ground), under the Sun and the day-night forcing of that hour (`var18`). Its
nudging, vertical wind, land, CCN and electricity are the coarse box's, its first
time step a third of the coarse one's, and it writes output every 15 minutes and
restarts every 3 hours over two model days. It builds its own clouds from that
air. The day it starts from is chosen once the coarse run has shown its storms,
a few hours before a stormy window; the case refuses to set up until it is set.
A check from day 5 ran its first five model minutes cleanly. It ran on
2026-10-04 from the coarse run's day 10.75, 41° past local noon and four hours
before the coarse run's first flash (day 10.92), through days 11–12, when
`box_0e`'s storms peaked.

Its two model days ran in 2.8 hours on 4 threads beside the coarse run, and no
flash struck. Building storms from the averaged air took most of the first
day: graupel formed at day 11.2 and the storms charged from day 11.5. From day
11.9 they stood as strong as the coarse box's (updrafts to 21 m/s against 21,
graupel and hail per unit area within 3 % of the coarse box's on average), with
cloud tops to 68 km against 86. They held up to 1,505 C of
positive and 1,603 C of negative charge at day 12.55, more per unit area than
the coarse box held between its flashes, but the strongest field aloft reached
131 kV/m at 34 km (day 11.89), 73 % of the breakdown field there, where a flash
starts at 90 %. Over the same days the coarse box flashed 56 times, in bursts
from seven or eight storms over nine times the area, four or five of them after
day 11.9. At that rate per unit area the fine box would expect half a flashing
storm while its storms were mature and see none 55 % of the time, so
these two days leave open whether the finer grid changes the lightning. The box
is also small for these storms: lunar storms stand four to five times taller
than Earth's, so 128 km is like a 30-km box on Earth, room for one or two storms
at a time.

At the author's direction it ran on from its day-2 restart, with the cap lifted
as the coarse box now has it, toward day 14.75 and the coarse box's busiest
lightning: after day 12.75 the coarse box had 11 flashing storms and 106 of its
first lunar day's 162 flashes. At the coarse box's rate per unit area that gave
about a three-in-four chance of a storm that flashes. Its first two days hold
under either setting, since they made no flash. Until day 13.0 its storms
matched the coarse box's (updrafts to 19 m/s, graupel and hail per unit area
close to the coarse box's) and its field reached 120 kV/m at 34 km. Then its
deep convection collapsed: from day 13.25 its cloud tops fell from about 66 km
to 34–48 km, its graupel and hail a hundredfold, and its strongest field to
17–24 kV/m, while the coarse box's storms kept their tops at 78–88 km and
flashed 94 times between days 12.75 and 13.84. A box this small holds one or two
storms, and their outflow and anvils may stabilize all of its air with no
surroundings to renew it; that is not yet checked. The author stopped it at day
3.0 of its run (coarse day 13.75), 2026-10-04, after three days without a
flash. The 2-km lightning stays unmeasured; a box wide enough to hold several
storms is the way to measure it.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_fine    # after setting fine_from['day']
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_fine --threads 4 --hours 24
```

### The first lunar day again

A case with `restart_from` runs on from another case's restart: it takes that
case's inputs, restart files and kept vertical field, and its namelist with its
own electricity, and its progress starts at the restart (`setup_restart`).
`box_0e_elec_uncapped` runs `box_0e_elec`'s first lunar day's storms again,
days 10.5–14.5 from the day-10.5 restart, before the first flash, with the
breakdown field's cap lifted, beside the first run under it.
`box_0e_elec_uncapped_corona` adds point discharge from the ground, with an
onset of 3 kV/m at 1.225 kg/m³ (dense vegetation; 3.4 kV/m at the box's 1.39
kg/m³). The two share an executable, so the second differs from the first only
in the ground's discharge. They show what point discharge changes: the field at
the ground and how much charge the ground gives off under storms, and whether
the space charge alters the lightning above. Later runs take it up once that is
known. On several threads CM1 does not repeat itself bit for bit (above,
"Threads and reproducibility"): the uncapped rerun parted from the first run
about 50 steps after the restart and flashed first at day 10.842, against
10.919, so each pair compares two realizations of the same day. The cap's
effect shows as a shift in the flashes' heights, depths and charge, and point
discharge's as the field at the ground. They ran on 2026-10-04 after the fine
box, 3.5 and 2.6 hours on 4 threads.

Over days 10.5–14.5, all flashes in cloud:

| | First run (capped) | Cap lifted | Cap lifted, point discharge |
|---|---|---|---|
| Flashes | 162 | 123 | 86 |
| First flash | day 10.919 | 10.842 | 10.915 |
| Starting below 34.8 km | 84 (52 %) | 40 (33 %) | 28 (33 %) |
| Start height, median (90th percentile) | 33.8 km (37.8) | 35.8 (37.8) | 35.8 (49.8) |
| Charge per flash, median (90th) | 87 C (174) | 95 (201) | 103 (227) |
| Channel bottom, median (10th) | 27.8 km (23.9) | 27.8 (20.4) | 29.8 (22.1) |
| Charge neutralized in all | 16,300 C | 13,900 | 10,100 |
| Nitrogen oxides | 389 mol | 333 | 237 |

Lifting the cap moved the flashes' starts up, as the higher threshold below
36 km requires: a third of them started below 34.8 km, against half under the
cap. Each flash then neutralized about a tenth more charge, and there were a
quarter fewer. The channels did not stop higher, as a higher stopping field
would have them do. With one realization each, the weather alone could make
differences that size in counts, charge and depth; the shift in starting
height follows from the threshold itself.

Point discharge left the strongest field at the ground much as it was, a
median over ten-step intervals of 34 kV/m against 40 and at most 181 against
228, since the ground's ions go into the rain under storms: 58 % of the
discharging columns were wet. Where no particles reached the lowest level it
cut the strongest field from a median of 11.0 kV/m to 7.2 and from 17.5 to 12.5
at the 90th percentile; in 22 % of intervals it still passed 10 kV/m somewhere,
against 55 % without it. Without point discharge the field at the ground passed
10 kV/m somewhere in the box in 95 % of intervals and 50 kV/m in 43 %. The
ground gave off charge from a median of 431 of the 4,096 columns at a time
(1,022 at most), a median of 3.4 A over the box and up to 34 A, 1.5 million C
of positive and 0.58 million C of negative charge over the four days. The
box's own charge stayed at a few thousand coulombs of each sign (at most
12,600 C positive, against 7,900 without it), so nearly all of that charge
went back to the ground on the precipitation it attached to; the totals
measure the model's ion cycle more than charge the storms keep. The run made 86 flashes,
fewer and starting higher than the run without it (a tenth started at or
above 49.8 km). The weather alone could make a difference that size between
realizations, so its effect on the lightning stays unresolved.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_uncapped    # and box_0e_elec_uncapped_corona
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_uncapped --threads 4 --hours 8
```

### Ground strikes

No lunar flash struck the ground in any run until 2026-10-05, and the zero was
an artifact. About a third of the lunar flashes start where a ground strike
can: lightmsz's reports give a potential beyond 50 MV at the start (a mean of
250–460 MV, at most 1.2 GV, against Earth storms' tens to about a hundred MV)
with the field pointing toward the ground, nearly all at negative potential.
Two things then stopped them. WRF-ELEC's Earth-calibrated stopping field, 0.15
of breakdown, halts every downward channel 15–25 km up (none got lower than
15 km). And the lunar rule, a channel within 5 km of the ground, kept WRF-ELEC's
demand that the channel's end lie in charge of the matching sign, which the
storms do not have that low. The lunar storms hold their charge as Earth's do:
in the flashing storms of days 40.75–41.9 the main negative charge lies at
36–56 km (−8 to −32 °C, strongest at 40–42 km), weak positive charge spreads
through the warm rain below the 0 °C level at 27.8 km, a weak positive layer
sits at 62–70 km, and air of −7 °C lies at 33.8 km, just beneath the main
negative charge as on Earth. So WRF-ELEC's own temperature rule fits them.

`box_0e_elec_ground_rule` runs the coarse run's busiest storms again, days
40.5–42 from its day-40.5 restart with the cap lifted, under that rule, and
`supercell_elec_ground_rule` runs the Earth benchmark under it. The Earth storm
makes 85 negative ground strikes among 7,406 flashes (1.1 %), the first at 36.8
min against the 34 min of a published run of this kind of storm (Ziegler et al.
2014), rising to 1.5–1.9 a minute in the second hour, a median 30 C each (6–122
C, 10th to 90th percentile) and 2 GJ. The lunar window makes 21 negative ground
strikes among 144 flashes (15 %; the first run made 157 flashes over these days,
all in cloud), the first at day 40.79. Each brings down a median 195 C of
negative charge (72–430 C, at most 889 C) and releases a median 69 GJ (20–190
GJ, at most 437 GJ), starting 34–38 km up. Earth's ground strikes bring down
5–30 C; charges of hundreds of coulombs belong on Earth to the rare positive
strikes that set off sprites.

The rule makes this an upper bound: it counts a channel that reaches −7 °C as
reaching the ground, which on Earth leaves 3–5 km to cross and on the Moon about
34 km. WRF-ELEC's stopping field, which halted every channel 15–25 km up, makes
the lower bound of none. Between them lies the physics of a leader crossing 34
km of clear air from a cloud at hundreds of megavolts: its channel's internal
field of a few kV/m costs about 50–125 MV over that distance, leaving most of
the potential at its tip, and Earth's leaders cross 10–15 km of clear air (bolts
from the blue) or climb 70 km (gigantic jets) on far less. A rule for that
crossing is the next step. On several threads the window is a separate
realization of these days, but the change from none to 21 is the rule's.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_ground_rule    # and supercell_elec_ground_rule
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_ground_rule --threads 4 --hours 6
```

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
