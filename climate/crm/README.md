# Storms and the lunar day in a cloud-resolving model

The GCMs (`../gcm/`) diagnose clouds and convection from formulas fitted to
Earth, on cells about 170 km across. This folder runs CM1 (George Bryan, NCAR;
release cm1r22.0, MIT-style licence), a model that resolves convection and forms
its own clouds and rain, at lunar gravity. It shows the weather people would meet
through the month-long day: storms, rain, cloud, winds from the ground to the
flight band, and the air near the ground by day and by night.

[wave_coverage.py](wave_coverage.py) exports geographic 10 m wind vectors and
surface density from the three corrected rings, preserving their three-hourly
histories and independent realisations. The [coastal-wave study](../waves/coastal.md)
records their spatial coverage and the remaining basin forcing requirements.

The [atmospheric-electricity study](../../research/studies/atmospheric_electricity/README.md)
(2026-10-02) reads these runs' storms: their graupel, cloud ice and snow, the air
where they meet supercooled water, and what slower settling at lunar gravity does
to them. Its "Where the study stands" gathers what the runs show about the
storms and what stays open.

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
  fall speeds scale alike. The single-drop solver of
  [atmosphere/rain](../../atmosphere/rain/README.md#the-check-on-cm1s-fall-speeds)
  checks the rain's factor (17 October 2026): with gravity alone changed in the
  Moon's air, it slows a 0.2 mm drop 4.4 times and a 6 mm drop 1.9 times, and
  meets the 2.9 near 1.4 mm. So the drizzle here falls too fast and the large
  drops too slowly, and per metre of fall the drizzle loses about a third less
  water than it would and the millimetre drops about a tenth more. CM1 was not
  rerun.
- **Rain's size limits are Earth's in every build but one.** Morrison's breakup
  holds rain's mean diameter near 0.60 mm and caps it at 2.8 mm, and the lunar
  runs' rain sits at that equilibrium; the `moon_omp_rain_limits` build scales the
  limits from the namelist ([Rain limits at lunar gravity](#rain-limits-at-lunar-gravity)).
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
A search on 2026-10-06 ran `box_0e_elec` from its day-10.5 restart for
8,000 s on 4 threads, three times with each of several parts switched off.
The runs parted at random times 3,000–7,400 s in, always first in the
domain's total vertical momentum, and the plain Morrison box parts the same
way, so at least one race lies outside the electrification. With mass
conservation off, all three runs differed from one another; with radiation
off, all three agreed, and with the large-scale forcing off, the two that
finished agreed. With the race striking about half the runs, three runs agree
by luck about one time in eight and two runs one time in four, so neither part
is pinned down. The boundary-layer scheme cannot be switched off under the
box's turbulence setting (`cm1setup` 2 requires it), so the test without the
surface and boundary-layer schemes did not run. The search stopped there; the
scripts are in `atmospheric-electricity/smoke/race_bisect.py` on the data
drive.
A case's executable changes only when its build is rebuilt (`build moon_omp`
replaces the one every such case links to), so the runner records the executable
and thread count of every segment.

**CM1's simple boundary layer gave the air its density times the surface's fluxes** (found 2026-10-21 in the Near-side
Sea box's water budget; fixed with the author's go-ahead of that day). Every ring and box here runs CM1's simple
boundary-layer scheme (`ipbl` 2, the Louis-type scheme of Bryan and Rotunno 2009), 51 of the 59 runs on the research
drive. Its driver (`turb.F`, block `dobr09`) takes its tendencies from the routines that return the change of the
density-weighted field, rho times the change per second. CM1's ordinary turbulence path divides that density out
(`turbtend.F`), but this driver hands the tendencies to `solve1.F`, which adds them as changes of potential
temperature, the mixing ratios and the wind. So the air took rho times what the surface passed it:
- in the sea box, 1.35 times the sea's evaporation;
- in the earlier land box, 1.32 times the evaporation and 1.46–1.51 times the sensible heat.

The air at the ground weighs 1.38 kg/m³ in the Moon's 1.2 atm, and the surface drag was inflated by the same factor.
The water budget closes only with the inflated term, to 0.01 mm/day. The same code stands in the cm1r22.0 release and
in NCAR's current main branch (read on 21 October).

The build `moon_omp_rain_law_pbl` divides the density out where the driver's block ends (`PBL_DENSITY_PATCH` in
[cm1_run.py](cm1_run.py)). In its trials the boundary layer's water came to 0.995–1.000 times the evaporation. CM1's
averaging of the scheme's wind tendencies from scalar points, which this scheme fills at u and v points, is left as
CM1 has it. Every earlier ring and box put 30–50% more of its surface's heat and water into its air than its surface
fluxes say, and dragged its low wind harder. Their near-surface air carries that bias in its temperature, humidity,
fog and wind at 10 m, and so do the products read from them, among them the wave study's winds
([wave_coverage.py](wave_coverage.py)), the comfort tables ([ring_comfort.py](ring_comfort.py)) and the storms the
electricity study reads. A fixed rerun of a reference case has yet to size the effect.
[box_budget.py](box_budget.py) now reports the boundary layer's water and heat over the surface's for every box.

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

### The 45° rings on the design ground, with the curved air's light

The two 45° rings above took their forcing from the run before the corrected design
case, whose land at 45–60° evaporated about ten times as much as the design climate's,
so their humid nights describe ground kept moist. Cases `ring_45n_lsw_twilight` and
`ring_45s_lsw_twilight` run them again with what the setup now takes and with the
curved air's light. Against `ring_45n_lsw` and `ring_45s_lsw` they differ in:

- **the forcing**, the corrected design case `A28_dim5_moon`, years 20–29: seas of
  22.7 °C at 45° N and 23.2 °C at 45° S, land ground at sea level 21.7 °C and 21.6 °C;
- **the land's wetness**, the GCM's own on those rows: 0.012 at 45° N and 0.036 at
  45° S, against 0.41 and 0.46;
- **the water along the circles**, from the current geography product: 34% and 31%,
  against 37% and 33%;
- **the build**, `moon_omp_twilight` with `var20 = 1`. After RRTMG, every column takes
  the curved air's sunlight by its own solar elevation, below the horizon too, and the
  infrared the same geometry lets out, from `terluna_twilight.txt`, which the setup
  writes from [../results/gcm/twilight.json](../results/gcm/twilight.json) in its full
  form (the shells' widening as well as the curved light and escape; [the GCM's
  account](../gcm/README.md#the-curved-airs-light)). Each shell's share heats or cools
  the levels it overlaps, the infrared scaled with each level's temperature to the
  fourth power; the ground takes its sunlight in the column's clear share and loses its
  extra infrared from what reaches it. Shells above the model's 150 km are left out with
  the air they stand for.

The air, the sunlight, the grid and every other setting are the earlier rings'. In a
two-hour test of each, the columns 2–10° past sunset gain 0.016–0.036 K a day of
sunlight's heating, averaged over their levels, where RRTMG gives them none.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_twilight
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_45n_lsw_twilight     # and ring_45s_lsw_twilight
/media/projectspace/terluna-research/session-tools/twilight_rings.sh                 # both, 6 threads; rerun to resume
```

**Paused on 14 October**, the northern ring at day 9 of 59 and the southern one before it began. Two of their inputs
changed while they ran. The twilight ledger's first channels put the near-infrared part of its terms too high, and the
corrected ledger has replaced them ([the GCM's account](../gcm/README.md#the-curved-airs-light)). The GCM's design case
is running again with its corrected radiation, which changes the air above 8–16 km that the rings are held to. The
rings wait for both.

### The 45° rings on the corrected radiation's design case

Cases `ring_45n_corrected` and `ring_45s_corrected` run the 45° rings on the design case of the GCM's corrected
radiation, `A28_corrected_970` years 10–19 ([the GCM's account](../gcm/README.md#the-curved-airs-light)), and
`ring_45n_corrected_twilight` and `ring_45s_corrected_twilight` the same with the curved air's terms from the corrected
ledger. A case's `gcm` setting names the run whose reference profile, vertical wind, land wetness and sunlight its setup
takes. The builds keep the air they were compiled with, which the setup checks against the run's; it is the same, so
only the sunlight moves, to 1,198 W/m² from 1,173. Against the twilight rings above they differ in:

- **the forcing**: seas of 22.9 °C at 45° N and 23.4 °C at 45° S, land ground at sea level 22.0 °C and 21.9 °C, about
  0.2–0.3 K warmer, and the land's wetness 0.014 and 0.046 against 0.012 and 0.036;
- **the twilight table**, the corrected ledger's (its channels keep the heights of the absorbed sunlight);
- **the controls**, the same rings on the build without the twilight hook.

The water along the circles, the air, the grid and every other setting are the earlier rings'. Two model hours of each
ran clean before the queue. The four run one after another on 6 threads, in blocks of two hours with half an hour's rest
after each:

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_45n_corrected          # and the other three
/media/projectspace/terluna-research/session-tools/corrected_rings.sh                 # the four; rerun to resume
```

**The north pair** (`python -m climate.crm.ring_analysis ring_45n_corrected_twilight --from-day 29.5`, and the
control, [ring_ring_45n_corrected_twilight.json](../results/crm/ring_ring_45n_corrected_twilight.json) and
[ring_ring_45n_corrected.json](../results/crm/ring_ring_45n_corrected.json)), over the second lunar day:

- **The curved air's light changes the ground by hundredths of a kelvin.** Over land the air at 2 m is 0.015 K cooler
  at dusk twilight, 0.021 K cooler by day and 0.013 K warmer in the deep night, with the same sign and size over each
  half of the lunar day. Humidity moves by tenths of a percent and the air aloft by 0.03–0.08 K, warmer near 100 km at
  dusk and cooler there at night. Both rings are held above 8–16 km to the same corrected GCM state, which already
  carries the twilight's large-scale effect, so the pair measures its local effect through the day.
- **On the design ground the ring is a desert.** Against the old ring on the over-wet ground: no rain against
  0.31 mm/day over land, no fog against 4.1% of samples, cloud over 1.4% of the land against 6.9%, land evaporating
  0.18 mm/day against 1.97. The air over land reaches 37.8 °C at its warmest local time against 30.3 °C, and 39.2 °C in
  the hottest 1% of samples; it falls to 17.2 °C at its coolest. Nights run 75% humid with dewpoints near 15.9 °C
  against 96% and 19.5 °C, days 38% against 62%, the highest wet bulb is 25.5 °C, and the mixed layer reaches 13.2 km
  by day against 7.1 km. The fog-desert nights the old rings showed belonged to their wet ground.

**The south** resumed on 16 October with its twilight ring alone, from day 7. Its control is dropped: the north pair
settled the twilight's local effect, and the twilight ring gives the southern climate on the design ground, whose land
is about three times as wet as the north's (0.046 against 0.014). Over its second lunar day
([ring_ring_45s_corrected_twilight.json](../results/crm/ring_ring_45s_corrected_twilight.json)) its land evaporates
0.69 mm/day against the north's 0.18, its nights reach 94% humidity against 88%, with fog in 0.4% of land samples
against none, and its afternoons peak at 36.0 °C against 37.8; it rains 0.01 mm/day, the north's desert all but
unchanged. Both rings drift warm from the GCM's state, and the GCM's sideways transport accounts for about a third
of it ([below](#what-else-divides-cm1-and-the-gcm-at-the-crossing)), so their absolute warmth runs high and their
humidity with it.

## Jura Island

Jura Island is the arc of Montes Jura, with the craters Sharp and Louville, round the flooded Sinus Iridum: 58,357 km²
at 38.4–51.6° N and 25.4–46.9° W in the 28% scenario's atlas, its summit 2,514 m above the sea at 45.6° N, 39.4° W.
Small islets lie within 10–100 km of its coasts and La Condamine Island (519 km²) 32 km off; larger land begins about
250 km out. Lunokhod 1 and Luna 17 lie 865 m down 61 km off its coast, a proposed undersea heritage community
([conservation](../../research/studies/conservation/README.md#9-candidate-community-settings)).

Case `box_jura` is a 3-D box about the island's centre (45.0° N, 36.1° W), 100 by 90 columns at the 45° N ring's
5,993 m, 599 by 539 km, so 70–100 km of sea lie beyond its outer coasts before its periodic edges. The island stands at
its atlas heights (up to 2,135 m at the box's spacing, a mean of 630 m), the atlas's sea and islets round it at sea
level, held at the design case's 296.1 K; the land is the placeholder land at the GCM's wetness there (0.014). The
reference is the corrected radiation's design case on the 45° N row, held from 8 km above each column's ground, with
the row's mean vertical wind at each column's height, the Coriolis force of 45° and the day-night circulation of
`ring_45n_corrected_twilight`; the sunlight is 1,198 W/m². It carries the curved air's twilight from the corrected
ledger, so the box is complete, although the 45° N pair put its local effect at hundredths of a kelvin. The hook
places each shell's share by its column's own face heights, so a column on high ground takes the shells above its
ground.

Two pieces of the setup are new. A box can extend a ring along a circle of latitude, taking that ring's GCM row as its
reference (`parallel_box_site`, and `ring_crossings.ring_position` for the day-night forcing). The build
`moon_omp_terrain_sea` makes the two-dimensional surface's water columns sea for CM1's surface layer
(`SEA_COLUMNS_PATCH`): the terrain build keeps the domain's land flag, under which a water column would run the soil
with water's properties. `moon_omp_terrain_sea_twilight` adds the twilight hook to it. In two model hours of a test
run every water column was sea and every land column land, the sea held its temperature and every field stayed
finite. The box runs in two-hour blocks on 6 threads with half an hour's rest after each:

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_terrain_sea_twilight
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_jura
/media/projectspace/terluna-research/session-tools/box_jura.sh                        # rerun to resume
```

**What the box leaves out.**
- **The land.** It is the placeholder land, bare and at one dryness everywhere, and its slab soil keeps a fixed
  moisture availability, so rain never wets it. The shore belts the fresh sea feeds, the gardens and the plants are
  missing; the [island study](../../research/studies/jura_island/README.md) sizes the belts at 2–4.5% of the island.
- **The sea.** It holds one temperature, the bay and the open sea alike.
- **The air above.** From 8 km up it is held to the design case's mean over the whole 45° N row, and the sinking and
  the day-night circulation are imposed from the row and the flat-land ring, so the island cannot feed back on them.
  The GCM behind them keeps PlaSim's Earth-tuned clouds, and the base design case carries no reduced gases, glow or
  their dimming.
- **Resolution and size.** The 6 km spacing smooths the arc to 2,135 m against LOLA's 2,717 m at 16 px/deg and
  resolves sea-breeze fronts and storms only coarsely. The periodic edges put copies of the island 140–200 km beyond
  its outer coasts. Two lunar days give one analysed day.
- **The physics.** RRTMG takes the Sun's own spectrum scaled to 1,198 W/m², and the twilight table comes from one clear
  design column. The boundary layer is Earth's empirical scheme run at lunar gravity. Cloud drops have a fixed number,
  with no aerosols, lightning or tides, and earthshine is left out.

**What the box shows** (run on 15–16 October, 59 days in 16.4 wall hours on 6 threads;
`python -m climate.crm.jura_analysis box_jura --from-day 29.5`,
[jura_box_jura.json](../results/crm/jura_box_jura.json)), over the second lunar day, on the island's 1,572 columns
beyond the blended edges, a mean 645 m up:

- **No rain and no cloud over the island**, at any height or hour. The arc, 2,135 m high at the box's spacing, lifts the
  air as a 350 m hill would on Earth ([the island study](../../research/studies/jura_island/README.md#the-heights)). A
  trace of cloud gathers over the open sea before dawn, 0.2–0.5% of it.
- **Hot days everywhere, warm slopes at night.** By day the air at 2 m is 31.3 °C at the median and 34.9 °C at the 95th
  percentile, as warm 1.8 km up as at the shore (31.2 against 30.3 °C), since a dry mixed layer cools by only 1.6 K a
  kilometre. At night the low ground runs coolest, 18.1 °C at the 5th percentile below 100 m, and the slopes at
  0.5–1.5 km warmest, near 21 °C, as cold air drains to the low ground. The dewpoint stays near 14 °C over the land and
  19 °C over the bay and the sea.
- **Calm air and a weak sea breeze.** The median wind at 10 m over the land is 0.4 m/s. On the coasts the wind blows
  onshore at about 1.1 m/s around noon on the bay's side and 0.9 m/s on the open sea's, and offshore at 0.2–0.4 m/s
  through the night.
- **Comfort.** The island is comfortable for 249 hours of each lunar day on the strict terms and 343 on the loose,
  too hot for about 300. Low ground below 250 m does best, 318–320 hours, and the slopes at 0.5–1.5 km worst, 172–187,
  their nights too humid. The land facing the bay is mostly the arc's inner slopes, a mean 1,090 m up, with 170 hours;
  its low ground below 500 m, 39 columns, has 218 hours against 307 on the open sea's side, with muggier nights.
- **Beside the GCM and the flat ring.** The GCM's two land cells within 150 km, 0.5 and 1.1 km up, average 21.5 °C with
  a lowest-layer dewpoint of 15.6 °C, 211 comfortable hours and 11.4 mm of rain per lunar day: milder days, and rain
  where the box has none. The flat ring's land near the site, at sea level, averages 26.0 °C with a dewpoint of 19.1 °C
  and no comfortable hours; the island's air is drier.

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

**At 45° N on the corrected design case** (the author asked on 16 October for this and the land test below to be
prepared; [gcm_profiles.py](gcm_profiles.py),
[profiles_ring_45n_corrected.json](../results/crm/profiles_ring_45n_corrected.json)). The land surface was left out of
the gap on wet ground, where the land evaporated 101–105 W/m². On the design ground the 45° N ring is a desert, so the
comparison was taken there. Each model's air is averaged over the GCM's lowest four layers, 0–1.8, 1.8–6.5, 6.5–13.7
and 13.7–22.9 km above the ground, by local time in the GCM's 36-degree bins. The GCM's land on that row stands a mean
2.2 km up and the ring's at sea level, so over land the two compare through differences in height above the ground;
over the sea they compare directly.

- **The ground heats the air alike.** Near noon the land gives the air 195 W/m² of heat in CM1 against 208–225 in the
  GCM, and both evaporate little, 13–21 W/m².
- **The lowest air takes different shapes.** By day CM1 mixes its lowest 7 km to the dry adiabat, 1.6–1.7 K a kilometre
  between its two lowest layers, under a mixed layer 14–19.5 km deep in the early afternoon; through the night its air
  at 2 m falls up to 5 K below the layer 1.8–6.5 km up. The GCM stays stable by day, 1.4 K a kilometre, and its air
  near the ground stays 2.7 K above that layer all night.
- **Over the sea CM1's whole lower air is warmer.** At the same heights CM1 runs 4.3–7.8 K warmer than the GCM in the
  lowest layer, 5–6 K in the second and 2.2–2.8 K in the third, and within about a degree in the fourth, where it is
  held to the GCM. Its air at 2 m is warmer than the sea, which takes 8–12 W/m² of heat from it by day and evaporates
  13–52 W/m² against the GCM's 64–97.

So CM1's warm offset belongs to its whole free lower air, over the sea as over the land.

**The heat budget of the ring's lower air** (the author's go-ahead of 16 October;
[ring_heat_budget.py](ring_heat_budget.py),
[heat_budget_ring_45n_corrected.json](../results/crm/heat_budget_ring_45n_corrected.json)) finds why. From the
ring's snapshots, as ring means over the second lunar day: its radiation, as its potential temperature takes it
(CM1 stores its radiative tendencies at constant volume; [cm1_output.py](cm1_output.py)), the imposed sinking on each
column's own air, the nudging toward the GCM's reference above 8 km and the storage, with CM1's own transport,
turbulence, convection and condensation as the rest. From the GCM's 3-day means at 47° N and the rows either side:
the heat its circulation carries sideways, which a ring cannot have.

| Layer | Ring against GCM | Radiation | Imposed sinking | Nudging | Rest | Storage | GCM's sideways transport |
|---|---|---|---|---|---|---|---|
| 0–1.8 km | 25.6 against 19.8 °C | −0.23 | +0.07 | 0 | +0.21 | +0.055 | −0.29 (eddies −0.22) |
| 1.8–6.5 km | 24.1 against 18.3 °C | −0.23 | +0.07 | 0 | +0.18 | +0.017 | −0.08 |
| 6.5–13.7 km | 16.5 against 13.7 °C | −0.15 | +0.03 | −0.12 | +0.26 | +0.016 | −0.06 |
| 13.7–23 km | 6.5 against 5.9 °C | −0.11 | +0.06 | −0.01 | +0.08 | +0.009 | −0.06 |

K a day of temperature. Below 8 km the ring's air radiates 27.2 W/m² away and takes 7.5 W/m² from the imposed
sinking and 22.8 from its own turbulence, convection and condensation; the 3.1 W/m² left over still warms it. The
ground gives the air 33.2 W/m² of sensible heat, and the 10.4 W/m² between the two is the heat the convection carries
above 8 km, less what condensation adds below it net of the rain's evaporation. The GCM's circulation carries about
15 W/m² of heat sideways out of the same air, 8.5 W/m² of it from the lowest 1.8 km and mostly by its eddies, and
dries the lowest layer by 0.2 g/kg a day.

- **The ring drifts from the GCM's state.** It starts from the GCM's profile, the reference it is held to aloft, and
  without the sideways loss its lower air warms until radiation balances: the lowest 1.8 km from 20.0 to 26.0 °C over
  the first lunar day and on to 26.7 °C by day 59, still 0.055 K a day warmer near the ground, and 1.8–6.5 km from
  18.2 to 25.0 °C. At 6.5–13.7 km, partly held, it warms by 3.6 K, and at 13.7–23 km, held, by under a degree.
- **Part of the gap lies in how CM1 is driven.** The rings and boxes take the GCM's sinking and its upper air but
  none of the heat and vapour its circulation moves between latitudes. Their absolute warmth near the ground runs
  high, at 45° N by about 6 K over the second lunar day and still growing, about a third of it for that reason (the
  confirming test below); the shapes within a lunar day, the breezes and the contrasts between places rest on the
  terrain, the sea and the Sun, as before.
- **What it leaves out.** The GCM's sideways terms come from 3-day means, which miss its faster eddies, so they are a
  floor. Its own radiative heating is not in its output, so the budget does not close on the GCM's side.

**The confirming test** (case `ring_45n_corrected_sideways`, the author's go-ahead of 16 October) is
`ring_45n_corrected` with that sideways heat and vapour imposed by height, the same at every hour, through the
large-scale advection the boxes use for the day-night circulation (`var11 = 1`, one local-time bin in
`terluna_lsadv.txt` on the GCM's levels): −0.28 K and −0.20 g/kg a day in the lowest 2 km. Its inputs match the
control's byte for byte except that file and the switch. In a two-model-hour test its lowest 2 km ended 0.021 K
cooler and 0.015 g/kg drier than the control's, as the forcing gives. It ran on 16 October in two-hour blocks with
half-hour rests, 59 days in 5.7 wall hours.

```sh
climate/gcm/.venv/bin/python -m climate.crm.ring_heat_budget ring_45n_corrected --from-day 29.5       # seconds
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_45n_corrected_sideways
/media/projectspace/terluna-research/session-tools/sideways_ring.sh                  # 6 threads; rerun to resume
```

**What the confirming test found**
([heat_budget_ring_45n_corrected_sideways.json](../results/crm/heat_budget_ring_45n_corrected_sideways.json), ring
means over the second lunar day): the sideways transport takes about a third of the warm offset away.

| Layer | Control | With the sideways transport | GCM |
|---|---|---|---|
| 0–1.8 km | 25.6 °C | 23.5 °C | 19.8 °C |
| 1.8–6.5 km | 24.1 °C | 22.3 °C | 18.3 °C |
| 6.5–13.7 km | 16.5 °C | 15.4 °C | 13.7 °C |
| 13.7–23 km | 6.5 °C | 6.1 °C | 5.9 °C |

- **Part of the drift.** Its lowest 1.8 km warm from 20.0 to 24.2 °C over the first lunar day and to 24.6 °C by day
  59, against the control's 26.0 and 26.7, and still warm 0.050 K a day near the ground against 0.055. About 2 K
  cooler, that layer radiates 2.9 W/m² less and takes 5.5 W/m² more from its own turbulence, convection and
  condensation, 3.1 of it the ground's added sensible heat, which together balance the 8.5 W/m² the transport takes
  from it.
- **What is left.** About 4 K below 6.5 km. Each kelvin the lowest layer cools gives back about 3.9 W/m², so the
  3.7 K left there would need about 14 W/m² more taken from it, 1.7 times what the GCM's 3-day means carry sideways.
- **The faster eddies are not it.** On 17 October the design case ran two more years writing a mean every 6.5 hours
  (`A28_corrected_970_eddies`; [gcm_eddies.py](gcm_eddies.py),
  [gcm_eddies_A28_corrected_970_eddies.json](../results/crm/gcm_eddies_A28_corrected_970_eddies.json)). Its eddies
  take −0.208 K a day from the lowest layer and −0.100 from the next, and the same outputs averaged into 3-day means
  give −0.218 and −0.092: over the lowest two layers 13.4 W/m² against 13.1, so the 3-day means miss 2%. At 45° N the
  GCM's eddies are slow, as on a world that turns once a month, and its two years agree with the ten the ring used.
- **So the rest lies in the models' physics:** how each sheds the heat the ground gives it, by radiation, convection
  and turbulence. By day CM1 mixes its lower air dry to 14–19 km where the GCM's stays stable. The land test waits
  behind this, since over land the ground heats the air alike in the two models.
- **CM1's longwave matches the correlated-k model.** On 17 October 990 of the ring's own columns, cloud-free below
  25 km over the second lunar day, went through the atmosphere domain's correlated-k model, which follows the
  line-by-line model ([ring_radiation.py](ring_radiation.py),
  [radiation_ring_45n_corrected_sideways.json](../results/crm/radiation_ring_45n_corrected_sideways.json)). CM1
  stores its radiative tendencies at constant volume, cp/cvm or about 1.38 times RRTMG's own rate in this air
  ([cm1_output.py](cm1_output.py)), and hands RRTMG the ground's temperature as the air's at the ground. Read back to
  RRTMG's rate, on columns built as CM1 builds them, its longwave comes within 0.013 K a day of the correlated-k
  model in every layer: 0.008 K a day more cooling in the lowest 1.8 km (0.2 W/m²), 0.010 at 1.8–6.5 km (0.7 W/m²),
  0.013 at 6.5–13.7 km (1.3 W/m²) and 0.005 at 13.7–23 km (0.6 W/m²). Its sunlight came within 2.5 points of
  line-by-line at the crossing, with 8 W/m² less absorbed in the air
  ([radiation, above](#what-else-divides-cm1-and-the-gcm-at-the-crossing)). CM1's radiation is sound on both counts,
  so the gap lies in the GCM or in how the two models carry the ground's heat.
- **The GCM's side.** The same morning a year of the design case wrote PlaSim's heating by process at every level
  (`A28_corrected_970_heating`, with the runner's new `--heating-diagnostics`; [gcm_heating.py](gcm_heating.py),
  [gcm_heating_A28_corrected_970_heating.json](../results/crm/gcm_heating_A28_corrected_970_heating.json)). At 47° N
  its lowest 1.8 km take 42.8 W/m² from the ground, pass 17.6 W/m² up through turbulence and dry convective
  adjustment, and radiate 9.9 W/m² away, 12.8 by longwave less 2.9 of sunlight; the 17.3 W/m² left leaves with its
  circulation, about 8.3 sideways and the rest in its vertical motions. The processes' own sum matches the total its
  diabatic step hands the dynamics within 0.02 K a day. On 151 cloud-free 3-day mean columns there, PlaSim's
  longwave cools the lowest layer by 0.44 K a day where the correlated-k model gives 0.07, 10.6 W/m² more, and the
  next layer 11.4 W/m² more; CM1's own longwave comes within 0.2 and 0.7 W/m² of the correlated-k model in the same
  two layers.
- **So most of the gap is the GCM running cool.** PlaSim sheds about 10 W/m² more by longwave from its lowest layer
  than CM1 does, worth about 2.7 K of the 3.7 at the ring's 3.9 W/m² a kelvin. PlaSim over-cools its lower air
  against the line-by-line model, and CM1's RRTMG matches that model. PlaSim's water continuum was fitted to the
  line-by-line model's outgoing longwave at the top, which leaves its cooling by height free; since 17 October the
  GCM's longwave correction fits it by height ([the GCM's account](../gcm/README.md#the-longwave-by-height)). The
  rest of the gap, about a degree, lies in how the two models carry the ground's heat upward: by day CM1's resolved
  dry convection where the GCM's stable layer leaves part of it to its resolved circulation.
- **RRTMG's own gas optics match the correlated-k model, beyond its table's edge too, and so does CM1's call**
  (19 and 20 October;
  [rrtmg_edge.py](rrtmg_edge.py),
  [rrtmg_edge_ring_45n_corrected_sideways.json](../results/crm/rrtmg_edge_ring_45n_corrected_sideways.json)). RRTMG's
  tables start at 1,053.63 hPa, and the ring's ground stands at 1,214 hPa, so its lowest 7.5 km are extrapolated in
  ln p. CM1's own RRTMG longwave, compiled from the lunar build's objects into an offline driver
  ([rrtmg_offline/lw_columns.f90](rrtmg_offline/lw_columns.f90)), took the same 990 cloud-free columns on the
  correlated-k model's layers, with the build's gravity, cp, CO2 and O2 and the air at the ground at the ground's
  temperature. It comes within 0.02 K a day of the correlated-k model at every height (−0.016, −0.006, −0.013 and
  −0.007 K a day over the GCM's lowest four layers), the same within 0.01 K a day on twins scaled so that every layer
  lies inside the tables, and holding the absorption at the edge's values below it moves the lowest layers by
  0.003 K a day. The edge needs no patch. CM1's own call, its stored lwten read back to RRTMG's rate, comes within
  0.008 K a day of the same RRTMG on the same columns in every layer (+0.008, −0.003, 0.000 and +0.002). Read as
  stored, its lwten had seemed to cool 0.14, 0.11, 0.09 and 0.08 K a day more: the constant-volume convention accounts
  for 0.07–0.11 K a day of that in every layer, and in the lowest 1.8 km the air at the ground, which the earlier
  columns carried down from the layers above, for another 0.07. Its gases, gravity and the ground's emissivity (0.98
  over water, land-use 16) match the driver's within what would move it, and its own layers above the top on Earth's
  standard temperature profile, the ice aloft and its 30-minute radiation step leave less than 0.01 K a day between
  them in these layers.
- **The ring's ice high up came from its start.** Every column holds ice above 107 km, 1.1 g/m², which is why the
  radiation check took columns cloud-free below 25 km. Above the GCM's top level (130 km) the starting sounding and
  the nudging reference kept the top's vapour while the temperature kept falling, so 140–150 km began 6–15 times
  saturated over ice. The ice formed within three days and stayed at the same amount through both lunar days,
  falling and subliming that slowly, with no point reaching the 125% of ice saturation the Morrison scheme needs to
  make more. At 63–107 km the ring's ice comes and goes, in 10–33% of columns, in air held to the GCM's mean vapour
  at 76–79% of ice saturation; where it holds ice it holds 1.5–3.5 times the water the GCM gives its clouds there.
  Since 17 October the runner keeps the top level's share of ice saturation above it, its vapour falling with the
  temperature (`held_share_above`); cases set up before keep their inputs. The GCM's side of the high cloud is in
  [its account](../gcm/README.md#the-high-cloud-and-its-ice).
- **For the cases that follow.** The sideways transport is part of the GCM's forcing that a slice of it lacks, so
  forced rings and boxes are better driven with it. Their warmth near the ground ran 3.7 K above the GCM's at 45° N,
  most of it the GCM's longwave over-cooling, which its correction since 17 October answers; a ring on the design case
  will show what is left.

**The land test on dry ground** (case `ring_45n_corrected_gcm_land`) is `ring_45n_corrected` with the GCM's roughness
(2 m) and soil thermal inertia (THERIN 4.97), keeping the GCM's own wetness as the control does. A circle of latitude's
setup now passes a case's land settings to its one land row. Its inputs match the control's byte for byte except the
two placeholder rows of the land-use table, and a two-model-hour test ran clean with the land's roughness at 2 m
against the control's 0.1 m. It waits, since the warm offset over the sea lies outside the land.

```sh
climate/gcm/.venv/bin/python -m climate.crm.gcm_profiles ring_45n_corrected --from-day 29.5      # seconds
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring_45n_corrected_gcm_land
/media/projectspace/terluna-research/session-tools/gcm_land_ring.sh               # 6 threads; rerun to resume
```

## Rain limits at lunar gravity

Morrison's scheme limits its rain with Earth's values. Self-collection is offset by breakup above a mean diameter of
0.30 mm, CM1's `brkup_db`, by the factor 2 − exp(2300 m⁻¹ (D_m − 0.30 mm)) (Verlinde and Cotton 1993), so the mean
diameter D_m = 1/λ of rain's exponential spectrum settles at 0.60 mm; that mean diameter is capped at 2.8 mm, and the
rain's mean fall speeds at 9.1 m/s in the scheme's reference air, which the Terluna builds carry to the Moon by the
single factor 2.94 as 3.09 m/s. No runner sets `brkup_db`; the namelist every case starts from writes CM1's 300 µm.
The world-crafting team's check of 19 October raised the same point (their B28). The author approved patching and
testing the limits on 19 October.

**Where the rain sits** ([rain_sizes.py](rain_sizes.py), `rain_sizes_<case>.json`). Over the second lunar day of four
runs, the breakup equilibrium sets the size of nearly all the rain:

| | Small box | Its Earth-fall twin | 45° S ring | Equatorial ring |
|---|---|---|---|---|
| Rain's mass at mean diameters of 0.45–0.80 mm | 95% | 85% | 77% | 91% |
| Beyond the breakup cutoff, 0.30 mm | 98% | 95% | 92% | 97% |
| Held at the 2.8 mm cap | 0.05% | 0.0005% | 0.003% | 1.0% |
| Rain in the cells held at the cap, median; in all rain cells | 0.04; 0.66 mg/m³ | 0.06; 0.71 mg/m³ | 0.06; 1.4 mg/m³ | 0.04; 2.5 mg/m³ |
| Rain aloft, mean | 0.40 kg/m² | 0.14 kg/m² | 0.013 kg/m² | 0.54 kg/m² |

The 45° N ring on the corrected design case holds 0.07 g/m² of rain aloft, about 6,000 times less than the small box,
too little to test the limits on.

**The lunar limits** ([atmosphere/rain](../../atmosphere/rain/README.md#breakup-at-lunar-gravity),
`rain_breakup.json`). Colliding drops coalesce or break by their Weber number, and with the single-drop solver's
terminal speeds a lunar pair of drops collides as an Earth pair of the same size ratio when it is s times as large. At
equal air density s is the ratio of capillary lengths, 2.46; the Moon's air, 1.2 times as dense, raises it to 2.78
for Morrison's spectrum at the equilibrium. The cap follows the largest drop that holds together, which grows with the
capillary length alone, and the fall-speed cap that drop's terminal speed, (g/9.81)^¼ of Earth's:

| | CM1 (Earth) | Lunar, at the design air |
|---|---|---|
| Breakup cutoff | 0.30 mm | 0.84 mm |
| Breakup rate | 2,300 /m | 826 /m |
| Mean diameter where breakup balances | 0.60 mm | 1.68 mm |
| Cap on the mean diameter | 2.8 mm | 6.88 mm |
| Cap on the mean fall speeds, in the reference air | 9.1 m/s (the Terluna builds: 3.09 m/s) | 5.80 m/s |

At the lunar equilibrium Morrison's mass-weighted fall speed in the Moon's air at 2 m is 4.47 m/s by its power law
(4.21 m/s by the solver for the same spectrum); the builds' cap holds it at 2.70 m/s there, and the scaled cap, 5.07
m/s, leaves it free.

**The cap stays, rescaled.** It holds thin rain sorted by its fall, a tenth to a sixtieth of the median rain cell's,
1.2–1.6 km up. In rain that thin, breakup at the cap restores the size over about four hours, 1/(5.78 × 312 × W) s at
W = 0.04 mg/m³, while falling rain sorts its sizes across a 100 m layer in about a minute, and below 1e-8 kg/kg the
scheme leaves self-collection and breakup out. Without the cap that rain's mean diameter would have no bound. Where the
rain is thick breakup holds it: at the cap with 1 g/m³ it restores the size in half a second.

**The patch** (`RAIN_LIMITS_PATCHES` in [cm1_run.py](cm1_run.py); build `moon_omp_rain_limits`). It adds three entries
to the namelist's `&param3` beside CM1's `brkup_db`: `brkup_gscale`, 0 to keep CM1's limits or 1 to scale them with
gravity; `brkup_rhor`, the surface air's density over Earth's sea-level air at the same temperature, which in CM1,
whose air takes Earth's gas constant, is the surface pressure over 1013.25 hPa; and `brkup_rhoexp`, the exponent of
that ratio in the scale of the collisions. With `brkup_gscale = 1` and the build's gravity g, the breakup's two
lengths scale by √(9.81/g) · `brkup_rhor`^`brkup_rhoexp`, the cap by √(9.81/g), and the fall-speed cap becomes
9.1 (g/9.81)^¼ m/s. At its default it keeps CM1's own limits. Its header in `morrison.F` states what it changes, its
default and its source, so other CM1 series can take it; it applies after Terluna's core patch, which gives the
build's gravity and the fall-speed cap. Terluna's runner writes the lunar values for a case with `rain_limits='lunar'`
(`brkup_gscale = 1`, `brkup_rhoexp = 0.687` from the breakup product, `brkup_rhor` from the case's sounding: 1.1965 at
the small box's 1,212 hPa) and CM1's own for `rain_limits='earth'`. No other build takes the patch, so every existing
executable and run stays as it was.

**The build and the smoke test** (19 October). `moon_omp_rain_limits` built cleanly at `-O3 -march=native` in four
minutes on two jobs (executable `1dc13187bb93677a`); every source it leaves alone matches `moon_omp_o3`'s byte for byte.
From the small box's day-42.0 restart, with 2.4 kg/m² of rain aloft, each ran two model hours on one thread. At its
default the new build repeated `moon_omp_o3` bit for bit in every output, statistic and restart file. With the lunar
limits it ran normally and logged the limits above. After the first hour the rain's mass-median mean diameter stood at
1.57 mm against 0.60 mm at the default; after two the rain aloft stood at 2.3 kg/m² against 2.8, and 3.1 mm of rain had
reached the ground against 0.93 mm. That is the transient of limits switched on in mid-storm, in one realization; the
test below measures them settled. The script and its record are `crm-runs/smoke_rain_limits/` on the research drive.

**The proposed test** (cases `box_0e_small_rain_control` and `box_0e_small_lunar_rain`, defined and not set up or run).
The small box's inputs as written run two lunar days on the rain-limits build, once with CM1's limits written out and
once with the lunar ones, so the pair shares one executable and differs in its namelist's rain limits alone; both write
the budgets of potential temperature and vapour and the radiative heating. Of the four runs measured above, the box
holds the most rain aloft after the equatorial ring and resolves it in three dimensions, and its rain test of 30
September frames the comparison: 85% of its falling water evaporated. Its forcing is the GCM's before the correction, so
it shows how far the limits move CM1's answers in a known rainy box; the design climate at 301 K has no CM1 case yet.
Over the second lunar day, against the spread between the control and `box_0e_small` itself (the same physics built at
`-O2`, another realization):

- the rain reaching the ground per lunar day, and how often and how hard it rains (the share of ground under 1 mm/h
  or more);
- the rain's sizes: its mass by mean diameter, beyond the cutoff and the equilibrium and at the cap (`rain_sizes.py`),
  and the spectrum by height;
- the microphysics' water budget: water condensed, evaporated back from cloud, evaporated from rain on the way down
  and reaching the ground, and the share of the falling water that evaporates;
- cloud by day and by night, fog and cloud by height;
- the air at 2 m and at the GCM's lowest layer, its dewpoint and comfort;
- the column's heating: the microphysics', the radiation's and the boundary layer's heating and moistening by height
  from the budgets, and the energy at the ground.

Each run takes about 2.3–2.8 wall hours on 4 threads (the rain test's twin took 2.32 hours), so the pair takes about
5–6 hours one after the other with a half-hour rest between, or about 3 hours side by side on 8 threads; each writes
about 12 GB with the budgets, where the small box wrote 5.7 GB. The analysis takes minutes: `box_analysis` with the
control as reference and `rain_sizes` for each arm, and a short reading of the budgets that is still to be written. A
cheaper pair would start both from the small box's day-29.5 restart, in about half the time, with a first lunar day
built under Earth's limits.

```sh
python -m atmosphere.rain.breakup                                           # about a minute: rain_breakup.json
python -m climate.crm.rain_sizes box_0e_small --from-day 29.5               # about a minute
python -m climate.crm.cm1_run build moon_omp_rain_limits --jobs 2           # four minutes
# the test, once approved:
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_small_rain_control
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_small_lunar_rain
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_small_rain_control --hours 4 --threads 4
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_small_lunar_rain --hours 4 --threads 4
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_small_lunar_rain --from-day 29.5 \
    --reference box_0e_small_rain_control
```

**The test as run** (cases `box_0e_small_rain_control_day29` and `box_0e_small_lunar_rain_day29`, the cheaper pair, on
19–20 October). Both restarted from the small box's day-29.5 restart on the rain-limits build and ran its second lunar
day, about an hour each on 4 threads, writing the budgets of potential temperature and vapour as domain means every
three hours. The spread between two realizations of the same physics is the control against the small box itself
over the same days: the same state at day 29.5, compiled differently. Box means over days 29.6–59.0
([box_analysis.py](box_analysis.py),
[box_box_0e_small_lunar_rain_day29.json](../results/crm/box_box_0e_small_lunar_rain_day29.json),
[box_box_0e_small_rain_control_day29.json](../results/crm/box_box_0e_small_rain_control_day29.json);
[rain_sizes.py](rain_sizes.py); [box_budget.py](box_budget.py)), with the 90% ranges of the moving-block bootstrap:

| | Lunar limits | CM1's limits | Lunar against CM1's | Control against the small box |
|---|---|---|---|---|
| Air at 2 m | 27.2 °C | 25.7 °C | +1.47 K (+1.24 to +1.71) | −0.13 K |
| Dewpoint at 2 m | 24.1 °C | 23.2 °C | +0.92 K (+0.68 to +1.14) | −0.17 K |
| Air as the GCM measures it | 27.9 °C | 25.8 °C | +2.06 K (+1.75 to +2.36) | −0.10 K |
| Rain reaching the ground, per lunar day | 63 mm | 33 mm | +90% | −2% |
| Water condensed | 174 mm | 307 mm | −43% | +1% |
| Falling water evaporated on the way down | 48% | 86% | | 85.7% against 85.3% |
| Rain held aloft | 0.17 kg/m² | 0.41 kg/m² | | |
| Rain's mass-weighted mean diameter | 1.31 mm | 0.59 mm | | |
| Cloud cover | 36% | 47% | −11 points (−18 to −5) | −2 points |
| Fog, cloud below 200 m | 30% | 31% | −1 point (−10 to +7) | −4 points |
| Cloud fraction of the lowest level, about 100 m | 28% | 17% | +10 points | −3 points |
| Wind at 10 m | 0.53 m/s | 0.62 m/s | −0.09 m/s | 0.00 m/s |

The lunar limits move every measure by several times the spread. Their rain grows larger, 82% of its mass beyond the
lunar cutoff of 0.84 mm, and evaporates less on the way down, so twice as much reaches the ground from 43% less
condensation: less of the air's water cycles back through falling rain.

**Why the air near the ground warms**
([box_budget_box_0e_small_lunar_rain_day29.json](../results/crm/box_budget_box_0e_small_lunar_rain_day29.json); the
spread's in
[box_budget_box_0e_small_rain_control_day29.json](../results/crm/box_budget_box_0e_small_rain_control_day29.json)).
CM1's own budget closes within 6 W/m² in each of the GCM's lowest four layers once the box's prescribed day-night
heating is reckoned from its table, and its microphysics matches CM1's running water totals within 4%. The lunar arm
against the control, W/m², positive warming:

| Layer | Microphysics | Of it: rain's evaporation; condensation | Turbulence | Circulation | Radiation | Storage |
|---|---|---|---|---|---|---|
| 0–1.8 km | +8.2 | +8.2; −6.1 | −12.3 | +8.3 | −1.4 | +2.1 |
| 1.8–6.5 km | +13.9 | +25.1; −14.7 | −0.2 | −9.0 | −1.9 | +6.1 |
| 6.5–13.7 km | −1.8 | +40.3; −45.4 | +1.8 | +5.9 | −0.4 | +5.2 |
| 13.7–23 km | +15.3 | +48.4; −30.1 | +0.2 | −1.4 | +0.2 | +1.2 |
| Whole column | +29.8 | +135.0; −124.4 | −10.3 | −2.4 | −3.4 | +19.1 |

The spread's differences stay within 3.3 W/m² in every term and layer. Turbulence is the boundary layer and
dissipation, and circulation the box's own transport with its advection's implicit diffusion. Until 21 October the
tool put the implicit diffusion with the turbulence. CM1's slot for it also carries the imposed sinking's heating, and
its vertical advection has that heating taken out (`adv.F`), so the two columns differed from these by up to
11 W/m² (box_budget schema 2).

- **Over the whole column** the rain's evaporation cools 135 W/m² less and condensation heats 124 W/m² less; the net,
  about 30 W/m², is the heat of the extra rain that reaches the ground.
- **In the lowest 1.8 km** the rain's evaporation cools the air 8.2 W/m² less. The warmer air takes 6.6 W/m² less
  sensible heat from the ground, which evaporates 8.2 W/m² more, and its turbulence brings it 12.3 W/m² less in all;
  its own circulation brings 8.3 more.
- **So the air at 2 m runs 1.5 K warmer and 0.9 g/kg moister**, most by night (+1.6 K, the dewpoint +1.5 K), and the
  low cloud reaches the ground more often while the fog below 200 m holds at 30%.
- **Against the GCM** at 0° E, which forced this box before its corrections: the box's air as the GCM measures it,
  25.8 °C with CM1's limits against the GCM's 22.0 °C, becomes 27.9 °C with the lunar ones, and its rain at the ground,
  33 mm a lunar day against the GCM's 180 on its land cells, becomes 63. The GCM's own falling rain evaporates by one
  Earth-tuned constant on each layer's dryness (PlaSim's `gamma`, `rainmod.f90`), with no drop sizes or fall speeds,
  so how much of its rain evaporates on the way down is the question this test puts to the GCM.

```sh
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_small_lunar_rain_day29 --from-day 29.5 \
    --reference box_0e_small_rain_control_day29                               # about a minute each
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_small_rain_control_day29 --from-day 29.5 \
    --reference box_0e_small
climate/gcm/.venv/bin/python -m climate.crm.rain_sizes box_0e_small_lunar_rain_day29 --from-day 29.5
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_0e_small_lunar_rain_day29 \
    --reference box_0e_small_rain_control_day29                               # seconds
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_0e_small_rain_control_day29 --reference box_0e_small
```

**Other Earth-tuned constants in Morrison's rain**, changed only where the physics gives the value:

- The self-collection rate, 5.78 m³ kg⁻¹ s⁻¹ (Seifert and Beheng 2001), carries Earth's collision kernel, the volume
  drops sweep at their difference in fall speed with collision efficiency near one. At the scheme's own lunar speeds it
  would be 2.94 times smaller, by the solver's 2.5–3.4 times for drops of 0.5–3 mm. It sets how fast the rain reaches
  the breakup equilibrium, whose size is the breakup law's. The rain law below takes its lunar rate from the collision
  kernels, 1.98 times smaller at the lunar equilibrium, since Earth's large drops approach one speed.
- Autoconversion and accretion (Khairoutdinov and Kogan 2000: 1,350 q_c^2.47 N_c^−1.79 and 67 (q_c q_r)^1.15) are fits
  to Earth's bin model. At lunar gravity droplets settle 6 times slower and collide less efficiently, so rain forms
  more slowly by a factor that needs a lunar collision calculation; they stay.
- The floor on the mean diameter, 20 µm (`LAMMAXR`), bounds drizzle at the cloud's edge and has no scaling with
  gravity of its own; it stays.
- For the ice, drops shed by wet graupel are taken as 1 mm (`PRACG/5.2E-7`), snow's and graupel's mean diameters are
  capped at 2 mm (`LAMMINS`, `LAMMING`), and their fall-speed caps follow the single factor; listed.
- The NSSL scheme of the electrified runs carries Earth limits of its own: self-collection reduced above a mean volume
  diameter of 0.6 mm by exp(−2500 m⁻¹ (D − 0.6 mm)) and stopped above 2.0 mm, and the mass-weighted diameter capped at
  6 mm (`xvrmx0` with `imaxdiaopt = 3`), its active breakup off for two moments. Carried as Morrison's are, the
  collisional lengths would scale by about 2.8 and the cap by 2.46; it is not patched, since no NSSL run is planned.

**Since.** The lunar limits change the answers by several times the spread between realizations (above). On 21 October
the author chose for CM1's next runs, the pair of boxes on the design climate, the lunar limits with the size-dependent
fall law below. The similarity carries Earth's collision outcomes to the Moon by the Weber number; no lunar collision
has been measured.

## The rain law at lunar gravity

Morrison's scheme reads rain's fall law only through three moments of its exponential spectrum, and no single power law
follows the single-drop solver through drizzle and large drops alike: refitted to the solver, the closest power law
stands as far from it as the lunar builds' present one
([atmosphere/rain](../../atmosphere/rain/README.md#the-fall-law-for-cm1)). The patch therefore takes the moments from
a table of the solver's, and gives rain's self-collection its lunar rate beside them. The author approved building it
on 21 October, before the pair of boxes over the Near-side Sea and at 0° E.

**The patch** (`RAIN_LAW_PATCHES` in [cm1_run.py](cm1_run.py); build `moon_omp_rain_law`, which carries the rain
limits and the twilight hook with it). Two entries join `&param3`: `rain_fall_law`, 0 to keep the build's power law or
1 to read `terluna_rain_law.txt` from the run directory, and `rain_selfcoll`, rain's self-collection rate, CM1's
5.78 m³/kg/s when 0. With the table, three things come from it, read linearly in the log of the spectrum's slope:
rain's mass-weighted and number-weighted mean fall speeds, in its sedimentation and in its collection of snow and
graupel, times the scheme's density factor; and the ventilation integral of its evaporation, times that factor's
square root. Rain's collection of cloud ice keeps the power law, and the rain limits' cap on the mean fall speeds
stays. For a case with `rain_law='lunar'`, which requires `rain_limits='lunar'`, the runner writes the table and the
lunar rate, 2.92 m³/kg/s, from atmosphere/rain's product. The patch's header in `morrison.F` states what it changes,
its default and its source.

**The build and the smoke test** (21 October; `crm-runs/smoke_rain_law/` on the research drive). `moon_omp_rain_law`
built at `-O3 -march=native` in three minutes on two jobs (executable `721fba67fa38d856`). From the small box's
day-42.0 restart, with 2.4 kg/m² of rain aloft, each arm ran two model hours on one thread. With the law written out
at its default, the build parted from the rain-limits build in the last digits: after two hours the rain at the ground
differed in its sixth digit. Built both strictly (`-ffp-contract=off -nostdinc`, without fused multiply-adds or vector
maths; `strict.py`), the two matched bit for bit in every output, statistic and restart file. So the patch at its
default leaves the scheme as it was, and at `-O3` the compiler parts them: the patch's branches change which of
Morrison's loops it vectorises, and so where it calls glibc's vector exp and log. With the lunar law the build logged
its 301 table rows and the rate, and its reading of the table at a slope of 1,000/m matched the product's mirror to
seven digits. Over the two hours it held slightly more rain aloft (2.37 against 2.33 kg/m²) at a smaller mass-median
mean diameter (1.51 against 1.60 mm), and 8% less rain reached the ground (2.85 against 3.09 mm). This is the
transient of a law switched on in mid-storm, in one realization; the boxes below measure the law settled.

```sh
python -m atmosphere.rain.fall_law                                          # about a second: rain_fall_law.json
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_rain_law --jobs 2       # three minutes
python3 /media/projectspace/terluna-research/crm-runs/smoke_rain_law/smoke.py setup|run|compare
python3 /media/projectspace/terluna-research/crm-runs/smoke_rain_law/strict.py build|setup|run|compare
```

## Boxes on the design climate: the Near-side Sea and 0° E

The small box's 48% of falling rain evaporated, against PlaSim's evaporation constant gamma, rests on one box over
land, forced by the GCM before its corrections and without the heat and vapour its circulation carries sideways. The
GCM evaporates a larger share of its rain over the Near-side Sea than at that box's place, and the gap widens with
gamma ([the GCM's account](../gcm/README.md#the-evaporation-of-falling-rain-21-october)). So, at the author's word of
21 October, a pair of boxes on the design climate puts gamma to the test where the land and the sea differ: one over
the Near-side Sea at 5° S, 45° W in the southern Procellarum region, the author's site, open water 120 km from the
nearest land; its twin over the land at 0° E. Both take the lunar rain limits and the rain law above.

**Forcing from the GCM alone** ([gcm_site.py](gcm_site.py)). No ring passes either site on the design climate, so each
box takes its forcing from the GCM's column there: the cells within 5.6° of latitude and 8.5° of longitude of the site
that have its surface, six sea cells at sea level and four land cells 0.7–1.3 km up. It takes the GCM's mean air
there as the reference and the GCM's mean vertical wind, which CM1 applies to its own air through the column. It also
takes the heat and vapour the GCM's circulation brings by local time, prescribed through the table the rings'
day-night forcing used (`var11`), through the column. The box's air is held to the GCM's only above the top of the
GCM's convection, and the first two attempts differed here (below). These come from the GCM's own budget, the storage
less the physics, and not from differences of its winds and humidity. The Moon's column at 301 K holds about 220 kg/m²
of vapour, so its circulation moves water through a column about a hundred times as fast as the column keeps it, and
the net at a site is a residual of about a millimetre a day. At the Near-side Sea, centred and spectral differences of
the 3-day means brought the column +0.1 to +2.5 mm/day, where the GCM's rain less its evaporation takes 0.7 away.

- **The GCM's budget.** A branch of the cold-cloud run, `A28R_301k_cold_moist`, six years from its year-19 restart
  with its settings, writes PlaSim's physics heating and moistening of each layer as the model applies them
  (`--moisture-diagnostics` with `--heating-diagnostics`). The heating comes at each grid point (code 21) and is
  projected here onto the T21 harmonics, as the spectral step takes it; the moistening comes as that step applies it
  (code 40, whose slot is filled in the diabatic spectral step), and code 41 holds the horizontal diffusion after it.
- **The snapshots and the exact means.** The diagnostics are single time steps at the end of each 3-day mean, at ten
  times of day. So the vapour over the column is held to the GCM's exact means: the change in its precipitable water
  less the evaporation, plus the rain. The difference is spread over the levels in proportion to the vapour each
  holds. Over the six years it was 0.24 mm/day at the sea, where the snapshots gave −1.02 and the means −0.78, and
  0.48 at the land, −0.33 against +0.15.
- **The vertical wind.** PlaSim's postprocessor writes omega (`wap`) about a hundred times too weak at the equator
  against the meridional mass flux of the same winds. The wind here is the GCM's own continuity: the divergence and
  the surface pressure's gradient from the written winds and pressure projected onto the harmonics, summed down each
  column over its layers as PlaSim sums them. Over the equatorial band it gives ascent of 1–3 mm/s, as the band's
  meridional mass flux does. At the sea it rises 3.5 and 4.0 mm/s at 0.9 and 4.2 km and sinks 3.4 mm/s at 10 km; at
  the land it sinks 4.3, 1.0 and 2.3 mm/s at the same heights.
- **Local time.** The means cover 36° of hour angle each, centred on the site's hour angle from the subsolar
  longitude of the mean solar height; the snapshots fall 18° later. Each level's series is fitted with the day's mean
  and first four harmonics, the means' harmonics divided by the damping a 36° mean gives them.

CM1 applies the mean vertical wind to its own air, so the table holds the circulation's heating and moistening with
the mean wind's advection of the box's starting air added back, on CM1's own levels. At the start the two cancel to
the accuracy of their differences (centred in the table, third-order upstream in CM1), and the box's own advection
follows its air as it changes. The GCM holds its tendencies as layer values between half levels, so each layer's is
laid on CM1's levels by pressure, linear within the layer about its value, its slope the monotonized-central limit of
its neighbours' (flat in the top and bottom layers). Every layer keeps its total; the 2% of the column's mass above
CM1's 150 km top is left out. Over the day the GCM's circulation cools and
dries the lowest kilometre at both sites, by 0.12 K and 0.43 g/kg a day over the sea and 0.24 K and 0.52 g/kg a day
over the land: it carries the evaporation and the ground's heat away. The box cannot make the planet's day-night
circulation, whose wind turns with the Sun: over the sea the GCM's wind 0.9 km up blows 3.7 m/s toward the east
around noon and 3.8 m/s toward the west in the evening, while its mean over the day is 0.5 m/s. So the box's mean
wind is nudged toward the GCM's by local time, in 72 windows each a 36th of the lunar day (0.82 Earth days), the
lowest level's speed 1.4–4.3 m/s through the day.

**The first attempt** (kept as `crm-runs/box_near_side_sea_split`). It took the table on the GCM's levels, the mean
wind's advection added back along the GCM's coarse gradient between its levels, and held the box's wind to the GCM's
mean. Over the sea, CM1's own advection by the mean wind, acting on the box's well-mixed boundary layer and its own
lapse, took back less than the table added: between 1.8 and 6.5 km the large-scale terms warmed the box by 8.5 W/m²
where the GCM's circulation cools that layer by about 11. With its wind near the mean, the sea evaporated 32–61 W/m²
against the GCM's 102. Over its first 13.6 days the box's air there warmed 0.3 K a day, its air at 2 m rose above the
sea's temperature, the sea's evaporation fell to 25 W/m², its CAPE to nothing, its condensation level to 8–11 km, and
it did not rain. It was stopped at day 15 and the pair began again with the table on CM1's levels and the wind by
local time. With them the sea evaporated 78 W/m² by day 7; on the Moon, where the dry adiabat cools 1.6 K a
kilometre, the morning air at 73% relative humidity still condenses 4–5 km up.

**The second attempt: the nudging stood in for the convection** (kept as `crm-runs/box_near_side_sea_nudged`, its 59
days complete; [box_budget_box_near_side_sea_nudged.json](../results/crm/box_budget_box_near_side_sea_nudged.json)).
It took the table and the mean wind to 16 km, as the rings' day-night forcing does, with the layers' values
interpolated in height and the lowest held to the ground, and held its air to the GCM's from 8 km, fully from 16 km.
Over its second lunar day the sea box rained 0.06 mm/day against the GCM's 2.72 at the site. Its water and heat
budgets, set against the GCM's own budget at the site, account for the gap. CM1's vapour budget closes to
0.01 mm/day. In mm/day:

| Over the second lunar day | Box | GCM at the site |
|---|---|---|
| Water into the air from the surface | 2.99 (the surface flux 2.23, the boundary layer's density above) | 3.51 |
| Removed by the forcing: the table; the mean ascent | 2.13; 0.49 | 0.78, the circulation's net |
| Removed by the nudging, all of it at 8–16 km | 0.35 | none |
| Change in the column's water | −0.04 | 0 |
| Rain at the ground | 0.06 | 2.72 |

- **Above 16 km the nudging did the convection's work.** The air there holds 76% of the column's mass. Radiation cooled
  it by 85 W/m², the nudging put back 83.5 W/m², and the box's clouds added nothing. Over the whole column the nudging
  supplied 89 of the 107 W/m² radiation removed, so the box needed about 2 W/m² of latent heat, which is the 0.06
  mm/day it rained. In the GCM the latent heat of its rain does that work, 79 W/m². Its convection scheme forms 3.5 of
  the 4.0 mm/day formed at the site, 3.2 of it between 18 and 72 km. Its circulation takes 35.5 W/m² from the air above
  16 km, which the box's forcing left out.
- **The box's convection stayed shallow.** Its clouds stayed below 10 km, its condensate peaking at 7 km under 1%
  cover. At lunar gravity that is the trade-cumulus regime, and it formed 0.41 mm/day of condensate where the GCM forms
  4.0.
- **Below 16 km the water left through the forcing.** The mean ascent moistened the starting air by 0.80 mm/day but
  dried the box's own air by 0.49. The box's shallow clouds and boundary layer had carried water from below 5 km up to
  6–10 km, where the GCM's mean flow sinks. Interpolating the GCM's layer values in height left the table's column
  0.55 mm/day drier than the GCM's.
- **The earlier land box shows the same shape more mildly.** It was forced by the GCM before its corrections. Its
  nudging supplied 50 of the 90 W/m² its radiation removed and took 5.3 mm/day of water back from the 4.2 the mean
  ascent brought ([box_budget_box_0e_small_lunar_rain_day29.json](../results/crm/box_budget_box_0e_small_lunar_rain_day29.json)).
  It rained 63 mm a lunar day against 180 in the GCM that forced it.
- **Two published results frame it.** Randall and Cripe (1999, *JGR* 104, D20) set out the ways to force a column
  model: prescribe the whole advective tendency; prescribe the sideways tendencies and let the observed vertical
  motion act on the model's own profiles; or relax the model toward observed soundings. These boxes combine all three,
  and the second, the mean ascent, turned against them. In the EUROCS study (Derbyshire and others 2004, *QJRMS* 130),
  two cloud-resolving models gave deep, raining convection in moister free tropospheres and only shallow convection in
  the driest. The study names mixing between cloud and surrounding air as the mechanism behind many differences among
  convection schemes. PlaSim's Kuo scheme lifts an undiluted parcel from the surface and rains a share of the column's
  moisture supply wherever that parcel is buoyant (`rainmod.f90`).

So the second attempt's 53% of falling water evaporated over that day (2.3 mm on the way down, 2.0 mm reaching the
sea, from 12.1 mm condensed; [box_box_near_side_sea_nudged.json](../results/crm/box_box_near_side_sea_nudged.json))
came from clouds below 10 km. Its rain fell a fraction of the way the GCM's falls, so the share tests nothing about
gamma. The small box's 48% stands with the same doubt.

**The pair as it runs now** (the author's go-ahead of 21 October). Three changes:
1. **The build** `moon_omp_rain_law_pbl` divides the density out of the boundary layer's tendencies.
2. **The forcing** reaches CM1's top, each GCM layer's total kept (above).
3. **The box's temperature and vapour** are held to the GCM's only above the top of the GCM's convection.

The convection top is the GCM's lowest level above which its convection forms less than 1% of its convective rain:
72 km at both sites, the sea's convection forming 99.7% of that rain at or below it. The air is held from the GCM's
next level up, 93 km, fully from its top level at 129 km, where CM1's own damping begins at 115.5 km. The GCM resolves
no tropopause to hold from: its temperature falls through its top level, to 166 K at 129 km at both sites, and its
lapse never falls below 0.8 K/km. The convection's top stands in for the tropopause. The wind is still nudged by
local time throughout. The mean ascent now acts above 16 km too: 3.6 mm/s at 40 km over the sea and 6.7 mm/s over the
land.

**The boxes** (cases `box_near_side_sea` and `box_0e_design`, on the build `moon_omp_rain_law_pbl`). They have the
small box's size and spacing, 32 × 32 columns 6.0 km apart. The sea is held at the GCM's 302.0 K. The land is the
placeholder land at the GCM's wetness, 0.83 (the 0.90 class, as box_0e's), at the ground's own height, with heights
above it. Both start at local sunrise and run two lunar days, the second analysed. They carry the curved air's
twilight from the 301 K ledger and write the heat and vapour budgets.

**The trials** (21 October; `crm-runs/trial_site_boxes/` on the research drive). Each case, copied, ran three model
hours on two threads, in about half a minute. Every trial ran normally and wrote its four outputs, with every field
the analyses read present and finite and the 22 fields of their heat, vapour and radiation budgets.
- **The first trials' logs** show the rain limits, the rain law's 301 rows and rate, the table, the vertical wind
  falling to zero above 16 km and the twilight's 167 elevations. With the split's fixes, the land box's showed the
  table's 111 levels and all 72 nudging windows read.
- **The second trials**, with the three changes, ran in 29 and 26 s. The boundary layer's water came to 0.995–1.000
  times the evaporation, against 1.32–1.35 before. Its heat came to 0.90–1.19 times the sensible heat, which is only
  3–4 W/m² at sunrise. The namelists hold the air from 93.1 and 92.8 km, and the vertical wind reads through to the top.

The pair runs one box after the other on four threads, with half an hour's rest between.

**The sea box with the three changes** (21 October; its second lunar day, days 29.5–59.0;
[box_budget_box_near_side_sea.json](../results/crm/box_budget_box_near_side_sea.json),
[box_box_near_side_sea.json](../results/crm/box_box_near_side_sea.json)). It ran its 59 days in 2.0 hours on four
threads. Its water budget closes to 0.06 mm/day, and its boundary layer passes on 0.998 times the sea's evaporation.
- **Rain.** The box rained 0.85 mm/day over the second lunar day, 31% of the GCM's 2.72, against 0.06 with the
  nudging. Its first lunar day barely rained, 0.02 mm/day: its upper air had first to cool enough for deep convection to
  start.
- **Convection.** Cloud water and ice now reach about 70 km and peak at 52 km, the heights where the GCM forms its rain.
  The column's radiation cools it by 50 W/m², against 107 in the nudged run, under the deep cloud.
- **Evaporation.** The sea evaporates 2.69 mm/day, against the GCM's 3.51.
- **Drift.** Without the nudging the upper air cools. At 40 km it went from 262 to 254 K over the two lunar days and at
  60 km from 238 to 224.2 K, still falling 0.1–0.2 K a day at the end, so the second lunar day is not yet the box's own
  equilibrium. (Corrected on 24 October: this entry first gave 226.5 K at 60 km; the domain diagnostics and the 3-D
  output at day 59.0 both give 224.2 K.)
- **The share.** 61% of the water falling from its clouds evaporates on the way down: 39.2 mm, against 25.2 mm reaching
  the sea, from 80.7 mm condensed. On PlaSim's pooled shares for this sea (49.7% at gamma 0.02, 75.8% at 0.04) that is
  gamma about 0.027, provisional on the drift and on the timing below; `rain_evaporation.py`'s fit follows with the land
  box.

**The gap follows the time of day.** The box takes the GCM's own circulation as prescribed forcing, and at this site
that circulation lifts the air by day and lets it sink by night. The box's rain follows it. Over the second lunar day,
by quarter:

| Second lunar day | Prescribed heating above 16 km | Prescribed water | Inhibition | Box rain | GCM rain |
|---|---|---|---|---|---|
| Sunrise to noon | −37 W/m² | −4.6 mm/day | 55 J/kg | 0.00 mm/day | 0.92 mm/day |
| Noon to sunset | −119 W/m² | +3.3 mm/day | 54 J/kg | 1.19 mm/day | 3.50 mm/day |
| Sunset to midnight | +97 W/m² | −2.1 mm/day | 3 J/kg | 2.04 mm/day | 3.54 mm/day |
| Midnight to sunrise | +93 W/m² | −4.5 mm/day | 3 J/kg | 0.17 mm/day | 2.93 mm/day |

- **Afternoon.** The rising air cools the upper air and moistens the column, and the box builds deep storms.
- **Evening.** They peak after sunset, spending what the afternoon built.
- **Night.** After midnight the sinking air warms the upper air by about 95 W/m² and dries the column, and the storms
  die.
- **Morning.** Through the morning the column dries further, the inhibition stands at 55 J/kg, the condensation level
  sits at 6.2 km, and nothing rains.

PlaSim's convection is a Kuo scheme (`rainmod.f90`). Wherever a parcel lifted undiluted from the surface is buoyant, it
turns the column's whole moisture supply into rain, the surface's evaporation included. So after midnight the GCM's
sea evaporates 3.4 mm/day and its scheme rains 2.9 of it while the circulation dries the column. That is the closure's
known behaviour:
- Kuo assumed that all the moisture supply makes cloud wherever the air is conditionally unstable (as reviewed by
  Villalba-Pradas and Tapiador 2022, *GMD* 15).
- Its moisture-convergence closure has been criticized for reversing the causality between convection and convergence
  (Emanuel, Neelin and Bretherton 1994, *QJRMS* 120; Arakawa 2004, *J. Climate* 17).
- Most single-column convection schemes start rain quickly where cloud-resolving models develop it through the day
  (Guichard and others 2004, *QJRMS* 130).
- ECMWF replaced its Kuo scheme with Tiedtke's mass-flux scheme in 1989.

So at least part of the gap is PlaSim's rain in the hours when resolved convection cannot run. The box's 61% is the
share of its own afternoon and evening storms. A third of the GCM's rain falls in the night and morning, when the box
makes almost none, so gamma fitted to the box measures how fast those storms' rain evaporates. It cannot settle whether
PlaSim should rain at those hours. The land box runs with the same changes as the Kuo-forced twin; over land the ground
heats and cools with the Sun, where the sea box's water is held at one temperature. PlaSim's replacement scheme is
planned in [the GCM's account](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october).

**The land box with the three changes** (21 October; its second lunar day;
[box_budget_box_0e_design.json](../results/crm/box_budget_box_0e_design.json),
[box_box_0e_design.json](../results/crm/box_box_0e_design.json)). It ran its 59 days in 2.3 hours.
- **Closure.** Its water budget closes to 0.04 mm/day, and its boundary layer passes on 0.997 times the ground's
  evaporation.
- **Rain.** It rained 3.28 mm/day over the second lunar day, 92% of the GCM's 3.58 on the same four land cells, and
  1.84 mm/day over its first.
- **Its water.** The ground evaporates 2.86 mm/day, against the GCM's 3.42, and the prescribed day-night forcing brings
  it 1.37 mm/day. Its column was still gaining 0.60 mm/day.
- **The share.** 61% of the water falling from its clouds evaporates on the way down.

Over land the two models rain at the same hours (mm/day):

| Second lunar day | Land box | GCM |
|---|---|---|
| Sunrise to noon | 2.87 | 1.91 |
| Noon to sunset | 8.59 | 11.06 |
| Sunset to midnight | 1.53 | 0.96 |
| Midnight to sunrise | 0.07 | 0.21 |

Both rain mostly through the afternoon and barely after midnight, because over land the ground's evaporation follows
the Sun and leaves Kuo's scheme little to rain at night. The sea box's gap is the sea's own: there the evaporation
goes on through the night, Kuo's scheme rains it, and the sinking air keeps CM1's convection down.

**What the pair says of gamma** ([rain_evaporation.json](../results/gcm/rain_evaporation.json)). Both boxes evaporate
61% of their falling water. PlaSim at one gamma evaporates far more over the sea than on land: 75.8% against 58.2% at
0.04. So the sea box's 61% comes at gamma 0.027, while the land box's lies beyond PlaSim's land share at 0.04, about
0.046 if the shares are extended in log gamma. No single constant fits both places, which points at the form of
PlaSim's evaporation of falling rain rather than its constant alone. Both shares are provisional: each box drifts
still, the sea box's upper air cooling and the land box's column gaining water.

**The near-ground gap, with the boundary layer fixed.** The rings of September ran 2.4–4.7 °C warmer than the GCM
near the ground over land and 1.8–2.4 °C over the seas, which made the GCM look cool. The corrected boxes' air at
2 m runs 1.15 K warmer than the GCM's over land and 0.23 K over the sea, at the same cells over the second lunar
day. So most of that gap was CM1's boundary layer passing on its density times the surface's heat and water. A
corrected rerun of a reference ring would show how much of it the rings carried.

**The sea box's third lunar day** (24 October; case `box_near_side_sea_day3`, days 59–88.5;
[box_days_box_near_side_sea.json](../results/crm/box_days_box_near_side_sea.json),
[box_budget_box_near_side_sea_day3.json](../results/crm/box_budget_box_near_side_sea_day3.json),
[box_box_near_side_sea_day3.json](../results/crm/box_box_near_side_sea_day3.json)). The sea box ran on from its
day-59 restart on the same build, inputs and settings, in 1.3 hours on three threads, its outputs numbered on from
474. The runner's restart setup now copies the rain law's table and carries the wind's windows by local time on by
whole lunar days: CM1 reads at most 100, so the run takes the 37 that span its days, each with the profile of the
first two lunar days' window at the same hour (`restart_windows`; the case records runner `0fd2316a921208a1`). A
3-hour trial from the restart read the 37 windows, the rain law's 301 rows and the restart's running totals, and
wrote every field the analyses read (`trial_sea_day3.py` in `session-tools/agents/climate` on the research drive).

| The sea box | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| Rain, second lunar day | 0.00 | 1.19 | 2.04 | 0.17 | 0.85 mm/day |
| Rain, third lunar day | 0.00 | 2.23 | 2.97 | 0.14 | 1.33 mm/day |
| CM1's inhibition, third lunar day | 37 | 35 | 3 | 2 J/kg | |

- **Rain.** The afternoon's rain rises by 87% and the evening's by 46%, and the box rains 49% of the forcing GCM's
  2.72 mm/day at the site, against 31% over the second day. The nights and mornings stay quiet, 0.14 and
  0.00 mm/day, under CAPE of 1,270–1,490 J/kg.
- **Closure.** The water budget closes to 0.04 mm/day, the boundary layer passing on 0.998 times the sea's
  evaporation of 2.98 mm/day, and the column loses 0.37 mm/day. The heat budget, with the storage taken between the
  day's sunrises, closes to 4 W/m² of the 43 W/m² the radiation removes; with box_budget's storage, a linear trend
  over the day that the day's own cycle tilts, it leaves 22 W/m² (`heat_storage_check.py` beside the trial).
- **Drift.** The upper air has nearly settled. Between the third day's sunrises the domain-mean air at 40 km went
  from 253.5 to 253.0 K and at 60 km from 224.2 to 224.5 K, and the 16–60 km band's mean by mass fell 0.5 K, against
  2.2 K over the second day. The lunar day's mean fell 1.5 K at 40 km and 1.4 K at 60 km from the second day, against
  4.0 and 6.1 K from the first to the second; the box started from the GCM's 261.9 and 237.7 K. Below 16 km the air
  cooled 0.5 K.
- **The share.** 60.4% of the water falling from the box's clouds evaporates on the way down (60.7 mm, against
  39.7 mm reaching the sea, from 133.5 mm condensed): gamma 0.027 on PlaSim's pooled shares, as over the second day.

So the second lunar day's pattern holds as the upper air settles: the sea box rains in the afternoon and evening and
barely after midnight, where the GCM with the new scheme rains 1.45 mm/day.

**The new scheme on the boxes' own air** (24 October; [sbm_offline.json](../results/gcm/sbm_offline.json), the
account in [the GCM's README](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october)). PlaSim's
simplified Betts–Miller scheme was run at once on each box's domain-mean column every 3 hours of its second lunar day,
laid on PlaSim's ten levels, at reference humidities of 0.6–0.8, parcel mixing of 0–16 per e-folding of pressure and
inhibition limits of 10 and 30 J/kg. Wherever it convects it rains tens to hundreds of mm/day, and on the sea box's
air it convects most after midnight: at rh 0.7 without mixing on every column, while the box rains 0.17 mm/day
there with 1,420 J/kg of CAPE and 3 J/kg of inhibition. Mixing raises its rain until the diluted parcel loses its
buoyancy, which on PlaSim's levels comes at every hour together at a mixing of 16, and no setting brings the sea's
nights near the box's while keeping its afternoons and evenings. On the land box's air the scheme convects on none of
the afternoons and on every night, although in the GCM it rains over that land at the box's hours. The box's
convection spends the vapour its columns gather, so the scheme's rate at once measures the vapour left above its
reference, most of it after the box's convection has stopped. Only the GCM can show when the scheme would rain in a
model running it.

**Why the sea box's nights stay quiet** (24 October; [box_onset.py](box_onset.py),
[box_onset_box_near_side_sea_day3.json](../results/crm/box_onset_box_near_side_sea_day3.json),
[box_onset_box_0e_design.json](../results/crm/box_onset_box_0e_design.json)). The diagnosis covers the sea box's third
lunar day and, for contrast, the land box's second. At every 3-hourly output, parcels rise through the domain-mean
sounding:
- **The parcels.** They start from the lowest level and from the lowest 0.5 and 1 km mixed, keep their moist static
  energy and drop their condensate. They rise undiluted or mixing in surrounding air at 1–8 per e-folding of pressure,
  as the scheme's `eps_bm` counts it. The lunar column's scale height near the ground is about 48 km, so 8 is
  0.17/km. If lunar convection is Earth's stretched by the ratio of the gravities, as the gravity pair tests, 8 is
  Earth's 1/km (Neale and others 2008 used 1e-3/m).
- **The forcing's dCAPE** (Xie and Zhang 2000). This is the CAPE's change over an hour under the day-night table and
  the mean ascent alone, applied above 1 km.
- **The 3-D output.** The strongest updrafts, the cloud tops, and CM1's own CAPE, inhibition and boundary-layer depth
  in each column.

By quarter, for the parcel from the lowest 0.5 km (undiluted unless marked):

| | Rain, mm/day | Boundary layer | CAPE | CAPE at eps 8 | Inhibition | LCL; LFC | LFC inside the boundary layer | Free-troposphere dCAPE | Strongest updraft, mean of outputs |
|---|---|---|---|---|---|---|---|---|---|
| Sea, sunrise–noon | 0.00 | 3.7 km | 1,646 J/kg | 10 | 39.5 J/kg | 6.4; 12.7 km | 6% of columns | +1.2 J/kg/h | 1.8 m/s |
| Sea, noon–sunset | 2.23 | 2.3 km | 2,008 | 37 | 38.1 | 7.5; 9.6 | 1% | +3.7 | 10.1 |
| Sea, sunset–midnight | 2.97 | 4.9 km | 1,876 | 16 | 1.4 | 4.7; 5.1 | 58% | −1.6 | 6.8 |
| Sea, midnight–sunrise | 0.14 | 6.1 km | 1,415 | 10 | 0.4 | 5.9; 6.2 | 59% | −1.6 | 1.5 |
| Land, sunrise–noon | 2.91 | 3.8 km | 1,900 | 46 | 99.7 | 5.9; 11.8 | 20% | +4.6 | 5.8 |
| Land, noon–sunset | 8.59 | 4.0 km | 4,593 | 35 | 9.0 | 4.3; 5.6 | 37% | +14.9 | 16.5 |
| Land, sunset–midnight | 1.53 | 0.2 km | 2,195 | 78 | 73.1 | 2.7; 9.0 | 0% | −4.3 | 2.2 |
| Land, midnight–sunrise | 0.08 | 0.1 km | 473 | 1 | 148.4 | 3.3; 15.2 | 0% | −4.1 | 2.3 |

- **Over land the quiet nights belong to the boundary layer.** After sunset it collapses to about 100 m, and the parcel
  meets 136–157 J/kg of inhibition with its LFC at 14–16 km. Convection starts each morning once the deepening layer,
  6.9 km, holds the LFC in 56% of the columns and the inhibition falls to 10–12 J/kg. Then the box resolves
  low-level updrafts of 2.6–3.0 m/s (the 99.9th percentile below 16 km).
- **After midnight the sea box meets those conditions and does not convect.** Its boundary layer is 6.1 km deep and
  holds the LFC in 59% of the columns; the parcel's inhibition is 0.4 J/kg and its CAPE 1,415 J/kg. Its resolved
  low-level updrafts reach 1.0–1.25 m/s, and its cloud stays shallow, with tops below 16 km in 14–25% of the columns.
- **The sea box convects on the CAPE the afternoon's circulation builds.** Its storms start at 1,822 J/kg against
  65 J/kg of inhibition, as the free-troposphere dCAPE runs at +3.9 to +4.2 J/kg/h. They peak at 2,609 J/kg and spend
  it, to 1,173 J/kg by midnight. Through the night the CAPE stays at 1,350–1,460 J/kg, and in the morning it rises
  to 1,630–1,680 J/kg under growing inhibition.
- **The forcing alone does not separate evening from night.** Its free-troposphere dCAPE is −0.7 to −2.4 J/kg/h
  through both the rainy evening and the quiet night. CAPE at eps 8 is 18–63 J/kg as the storms start and peak, 5–16
  at night and in the morning, and 0–1 through the evening's decay.

So over the sea a trigger on whether boundary-layer air reaches its LFC would favour the nights, and the dCAPE treats
the evening and the night alike. The CAPE the circulation builds marks the box's convection, carried into the evening
by the storms the afternoon starts. The box's own grid is a candidate cause of the quiet nights. Its night boundary
layer is as deep as its 6-km spacing is wide, so the boundary-layer scheme carries the layer's eddies, which the grid
cannot resolve. Petch and others (2002, *QJRMS* 128) found deep convection's onset delayed on grids too coarse for such
eddies. A finer box would show whether the sea's nights convect.
[The GCM's account](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october) tests a CAPE trigger on the
GCM's own columns.

**The sea box's night on finer grids** (25 October; cases `box_near_side_sea_night_6km`, `_3km` and `_2km`;
[box_resolution.py](box_resolution.py),
[box_resolution_near_side_sea_night.json](../results/crm/box_resolution_near_side_sea_night.json)). Three boxes put
that candidate to the test, at the author's go-ahead of 25 October. Each starts from the settled box's own night, and
they differ from one another in their grid alone.
- **The start.** Each takes box_near_side_sea_day3's air at day 84.5 (its output 677), averaged over its columns, and
  CM1 adds its random perturbations of up to 0.25 K at every point. That hour lies 40° of hour angle past midnight on
  the third lunar day, 14 minutes into a window of the wind by local time. There the source box rained 0.16 mm/day
  under 1,305 J/kg of CM1's CAPE and 0.6 J/kg of inhibition, its boundary layer 6.0 km deep and holding the LFC in
  61% of its columns.
- **The forcing, held at that hour for two Earth days** (`setup_night`). The clock that times the Sun and the day-night
  table (`var19`) is stopped. The table keeps its values at the hour angle −139.9°, warming the air above 16 km by
  87 W/m² and drying the column by 4.7 mm/day. The Sun stays 49.6° below the horizon, beneath the twilight table's
  lowest elevation (−40°), while the twilight's infrared acts as at every hour. The wind is nudged toward that window's
  for the whole run. The hold on the air above the GCM's convection, the mean vertical wind, the sea at 302.0 K, the
  rain law and limits, the build `moon_omp_rain_law_pbl` and its physics are the source's. The cases record runner
  `a0e3d9887788dbea`.
- **The grids.** The control is the box's own grid, 6 km on 32 × 32 columns (192 km). The finer boxes keep domains
  wide enough for deep cells: 3 km on 64 × 64 columns (192 km) and 2 km on 64 × 64 (128 km). All keep the 111 levels
  to 150 km and CM1's boundary-layer scheme. Their first time steps are divided as the spacing, and CM1's adaptive
  step settles near 75, 39 and 26 s. They ran in 4, 34 and 49 minutes on four threads.
- **What they write.** The 3-D output comes every 6 hours without the budgets' fields (the domain diagnostics keep
  their means), the domain diagnostics every half hour, CM1's statistics every 15 minutes and restarts every
  12 hours. The three runs wrote 1.33 GB and their smoke tests 0.34 GB; the box's third lunar day wrote 5.6 GB.
- **The smoke tests** (`trial_night.py` in `session-tools/agents/climate-resolution` on the research drive). Each
  case, copied, ran half a model hour on four threads in 3–28 s and wrote its 3-D output, diagnostics, statistics and
  a restart, with every field the analysis reads present and finite. Their logs show the one wind window, the
  table's 36 bins, the twilight's 167 elevations, the rain law's 301 rows and the stopped clock.

**The rule, set before the runs** (`CONTROL` and `READING`). Over its second day the control must reproduce the
source's quiet night: rain under 0.3 mm/day, the strongest updraft below 60 km averaging under 3 m/s, the 99.9th
percentile of w below 16 km under 2 m/s in every snapshot, and cloud tops above 40 km in under 5% of the columns on
average over its snapshots. Otherwise the method fails and the finer grids go unread. A grid convects when it rains
at least 0.5 mm/day, a third of the GCM's 1.45 after midnight, with its strongest updraft averaging 5 m/s or deep
tops in 5% of its columns.

**The control reproduces the quiet night, and the finer grids stay as quiet.** Over each box's second day, beside the
source box over the same hours of its own run (48 diagnostics and 4 snapshots in each box, 8 of each in the source):

| Second day | Source box, days 85.5–86.5 | 6 km, the control | 3 km | 2 km |
|---|---|---|---|---|
| Rain | 0.03 mm/day | 0.07 | 0.05 | 0.04 |
| Strongest updraft below 60 km: mean; highest | 1.28; 1.51 m/s | 1.45; 1.81 | 1.88; 2.18 | 1.96; 2.36 |
| w's 99.9th percentile below 16 km; at 16–60 km | 1.02; 0.10 m/s | 1.11; 0.13 | 1.35; 0.17 | 1.45; 0.21 |
| Columns with cloud tops below 16 km | 36% | 48% | 47% | 43% |
| Highest top of the cloud below 40 km: mean; highest | 10.5; 25.8 km | 8.6; 8.6 | 8.9; 9.3 | 8.7; 9.3 |
| Columns with ice standing apart at 50–52 km | 0.3% | 0.8% | 1.6% | 3.0% |
| CM1's CAPE; inhibition | 1,292; 1.6 J/kg | 1,343; 0.8 | 1,320; 0.7 | 1,313; 0.6 |
| The 0.5 km parcel of the mean sounding: CAPE; inhibition; CAPE at eps 8 | 1,440; 0.5; 13 J/kg | 1,500; 0.5; 13 | 1,476; 0.6; 13 | 1,468; 0.7; 13 |
| Boundary layer | 6.1 km | 6.1 | 6.0 | 6.0 |

- **The finer grids resolve stronger eddies, and their cloud stays shallow.** The 99.9th percentile of w below 16 km
  rises from 1.11 m/s at 6 km to 1.35 at 3 km and 1.45 at 2 km, and the strongest updraft from 1.45 to 1.96 m/s, while
  the boundary layer stays 6.0–6.1 km deep. The cloud rising from below tops out at 8.6–9.3 km on every grid, no column
  holds cloud at 16–40 km, and the rain falls slightly as the grid refines. The two halves of the second day differ by
  under 0.1 m/s in the strongest updraft on every grid, so none was building toward an onset when the runs ended.
- **The cloud above 40 km stands apart.** A thin layer of ice, up to 0.03 g/kg at 50–52 km, comes and goes over
  0–11% of the columns with no cloud beneath it at 16–40 km. The source box carries it at night too, so it forms in
  place. The rule's share of cloud tops above 40 km counts it, and averaged over the second day it stays under 5% on
  every grid.
- **CAPE stays high and inhibition nil on every grid.** The mean sounding's parcel reaches its LFC at 6.2–6.5 km,
  and in 47–57% of the columns CM1's own LFC lies inside the boundary layer. Mixing at eps 8 leaves the parcel 13 J/kg
  on every grid, as in the source's quiet hours; the source's storms start and peak with 18–63 J/kg.
- **Drift.** From the homogeneous start each box's boundary layer reaches 6.1 km within its first hour, and its
  shallow cloud covers 29–34% of the columns by hour 12. Over the two days the air at 60 km cools 1.37 K on every
  grid, against 1.33 K in the source over the same hours. The lowest 1.77 km cools 0.03–0.06 K and dries
  0.23–0.25 g/kg, against 0.10 K and 0.32 g/kg in the source, whose forcing grows drier toward morning (5.3 mm/day by
  the end of its two days).
- **Water.** The boundary layer passes on 0.998–1.003 times the sea's evaporation of 3.05–3.08 mm/day. The vapour
  budget's rest over the second day, 0.13–0.17 mm/day, matches the vapour CM1's microphysics puts back above 55 km,
  0.160 mm/day, where the held table dries air that holds none; the source's at those hours are 0.157 and 0.135. Less
  that refill, the budgets close to 0.03 mm/day.

So under the rule set before the runs none of the three convect, and the quiet sea nights hold up across the grids.
The rings' and boxes' 6 km grid passes this test at the sea box's night: at 3 and 2 km the boxes resolve stronger
eddies under the same CAPE and the same nil inhibition, and their cloud stays below 10 km. The GCM's next candidate is
the CAPE threshold near 100 J/kg that
[its account](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october) points to (`--sbm-cape`, off by
default), for a branch at the author's word. The test is one realization of each grid over two Earth days, from one
hour's mean air with that hour's forcing held. CM1's boundary-layer scheme stays on at every spacing, as in every lunar
case, so the test measures the grid within that set-up. The 2 km box is 128 km wide against the others' 192.

**What ends the night's cloud: its entrainment under the night's sounding** (25 October, at the author's go-ahead;
[box_entrainment.py](box_entrainment.py),
[box_entrainment_near_side_sea.json](../results/crm/box_entrainment_near_side_sea.json), 22 seconds, no new runs). With
the grid ruled out, the box's own clouds show what stops them. Every 3-D output of the sea box's third lunar day is
sampled level by level, as Siebesma and Cuijpers (1995, *J. Atmos. Sci.* 52) sampled large-eddy trade cumulus: 59
snapshots 3 hours apart in each quarter of the lunar day, each of 1,024 columns, and the held-night boxes' 7 snapshots
from hour 12 (1,024 columns at 6 km, 4,096 at 3 and 2 km).
- **The samples.** Cloud is air holding at least 0.01 g/kg of cloud water and ice. The cloudy updraft is cloud that
  rises, and the core is cloudy updraft lighter than its level's mean, counting the weight of every water species.
  Each sample's environment is the rest of its level.
- **The entrainment.** Total water qt and moist static energy h, which condensation leaves unchanged, give a bulk
  plume's fractional entrainment ε = −(dφ_c/dz)/(φ_c − φ_e) from each quarter's pooled profiles; the detrainment is ε
  less the mass flux's fractional growth with height. Rain that forms takes qt away and freezing adds to h, so the
  night's warm cloud, which barely rains, keeps both, while in the afternoon's deep cells neither is kept exactly. The
  bulk rate is the mass-flux-weighted mean of ε below 16 km, with the 5–95% range of 300 resamplings of the snapshots.
- **The plume.** A plume with the mean air of the lowest 0.5 km (box_onset's parcel) rises undiluted to its
  condensation level, then mixes in a share 1 − exp(−ε Δz) of the surrounding air's h and water between levels of
  each quarter's mean sounding, dropping its condensate as it forms (with its condensate kept as weight, in brackets
  below). Its cap is the first height above its level of free convection (LFC) where its buoyancy vanishes.
- **The tops.** Each column's cloud top is the top of its lowest unbroken cloud below 40 km, which leaves out the ice
  that forms in place near 50 km.

| | Night | Morning | Afternoon | Evening |
|---|---|---|---|---|
| Rain | 0.14 mm/day | 0.00 | 2.23 | 2.97 |
| Cloudy columns; their tops' median and 90th percentile | 22%; 7.3, 7.9 km | 0.9%; 6.7, 7.3 | 10%; 22.1, 27.8 | 10%; 7.3, 25.8 |
| Cloudy updraft's entrainment by qt (5–95%); by h | 1.34 (1.20–1.53); 1.10 × 10⁻³ m⁻¹ | 0.50; 0.32 | 0.33 (0.30–0.35); 0.22 | 0.31 (0.29–0.34); 0.24 |
| Core's entrainment by qt (5–95%); by h | 0.57 (0.55–0.61); 0.18 × 10⁻³ m⁻¹ | too few cores | 0.29 (0.25–0.31); 0.16 | 0.27 (0.25–0.29); 0.18 |
| Undiluted parcel's CAPE below 10 km; at 10–16 km; above | 5; 32; 1,378 J/kg | 0.6; 12; 1,632 | 0.1; 54; 1,944 | 19; 84; 1,733 |
| Plume at the quarter's core rate by h: inhibition; LFC; cap | 0.5 J/kg; 6.2; 8.3 (6.9) km | no cores | 39; 10.0; 26.8 (24.0) | 1.5; 5.3; 7.2 (6.6) |
| Plume at 1 × 10⁻⁴ m⁻¹: inhibition; LFC; cap | 0.5 J/kg; 6.2; 8.8 (7.2) km | 54; 17.4; 31.4 (28.6) | 38; 10.0; 32.9 (28.3) | 1.5; 5.3; 27.8 (23.9) |

**At night the boundary layer's air meets a layer of almost no buoyancy, and mixing ends its cloud there.** Nearly all
of the night's 1,415 J/kg of undiluted CAPE lies above 10 km. Between its LFC at 6.2 km and 10 km the parcel gains 5
J/kg, at most 0.0013 m/s² of buoyancy (0.2 K), and carrying its condensate it loses even that at 8.0 km. The cloudy
updrafts mix in their surroundings at 1.1–1.3 × 10⁻³ m⁻¹ and their buoyant cores at 0.18–0.57 × 10⁻³ m⁻¹. At these rates
the plume stops between 6.9 and 8.3 km or never becomes buoyant, where the columns' clouds end: half by 7.3 km, 90% by
7.9 km. The cloud itself shows the same end. Its rising air is 0.05–0.10 K heavier than its level, and its mass flux
falls fivefold from 6.2 to 7.9 km, 0.019 to 0.004 kg/m²/s, as it detrains 1.0–2.4 × 10⁻³ m⁻¹ into a layer that cloud
covers 6–10% of. Above 10 km the parcel's buoyancy grows to 0.012 m/s² by 16 km in air that warms 0.51 K/km in potential
temperature, less stable than the parcel's own moist adiabat, and that deep layer holds the CAPE the night's diluted
cloud never reaches.

**In the afternoon the same mixing carries cloud deep once the ascent has lifted it through the inhibition below
10 km.** There the parcel is 0.9–1.4 K colder than its surroundings at 4–8 km, and the plume crosses 38–41 J/kg of
inhibition. The cloudy updrafts cross that layer as forced cloud, 0.06–0.24 K heavier than their levels, rising at
0.6–0.7 m/s under the prescribed ascent. Above 10 km the parcel's buoyancy grows to 0.017 m/s² by 16 km (54 J/kg
between 10 and 16 km, against the night's 32). The updrafts' mass flux grows more than tenfold from 6.7 to 12 km, and a
plume at the afternoon's core rate by h (0.16 × 10⁻³ m⁻¹) stops at 26.8 km, or 24.0 km with its condensate, beside the
columns' tops: half by 22.1 km, 90% by 27.8 km. At the core rate by qt (0.29 × 10⁻³ m⁻¹) it stops at 12.6 km, or never
becomes buoyant with its condensate, so the two rates bracket the afternoon's tops.

**At a common rate the two soundings diverge.** At 1 × 10⁻⁴ m⁻¹, below every measured rate, the plume
stops at 8.8 km at night and at 32.9 km in the afternoon; at 2 × 10⁻⁴ m⁻¹ at 8.3 and 24.3 km. The night's core rate
by h carries a plume to 25.3 km through the afternoon's air, and the afternoon's to 8.4 km through the night's. So the
sounding sets the difference, and the night's mixing deepens it: its cloudy updrafts mix four to five times as fast
per metre as the afternoon's, its cores up to twice as fast.

**The night's air above the boundary layer, and what the sinking does to it.** Between 6 and 10 km the potential
temperature rises 0.74–0.95 K/km at night against 0.13–0.17 K/km in the afternoon. At 10–14 km the relative humidity
is 0.73–0.77 against 0.82–0.89, and at 10–12 km the saturation deficit h* − h is 5.9–6.2 kJ/kg against 2.8–3.4. The
stable layer stands on a cooler boundary layer: at 1 km the night's potential temperature is 283.6 K against the
afternoon's 285.5, while at 10 km the two match, 287.3 and 287.4 K. The prescribed sinking keeps that shape. After
midnight it cools the lowest 6 km by 18.7 W/m² and dries them by 2.9 mm/day against the sea's evaporation, and it warms
6–10 km by 11.6 W/m² and 10–16 km by 25.4 W/m², where radiation cools them by 13.4 and 21.6 W/m². So the sinking's
warming holds the air above the cloud at the afternoon's temperature while the boundary layer beneath cools. In the
afternoon the forcing cools 6–16 km by 15.2 W/m² and moistens it by 2.2 mm/day, and radiation adds 1.5 W/m² at 6–10 km
and takes 1.8 at 10–16 km.

**Against BOMEX, the trade-wind case.** BOMEX's cumulus grow under subsidence reaching 0.65 cm/s at 1.5 km, radiative
cooling of 2 K/day and drying near the sea, below an inversion at 1.48–2.0 km. Its LES clouds mix in 1–3 × 10⁻³ m⁻¹,
0.8–1.7 × 10⁻³ m⁻¹ by the bulk-plume method this diagnosis uses. The night box has the regime's parts: sinking that
warms the air above the cloud and dries the boundary layer against the sea's evaporation, a cloud layer spreading
beneath a stable layer, little rain, and entrainment and detrainment near 10⁻³ m⁻¹. Its lid is a layer where the
parcel's buoyancy almost vanishes, beneath a deep layer holding nearly all of its CAPE, where BOMEX's is an inversion.
If lunar convection is Earth's stretched by the ratio of the gravities, as the gravity pair tests, a BOMEX-like cloud
would mix 0.13–0.28 × 10⁻³ m⁻¹ by the bulk method on the Moon; this inference from scaling alone puts the night's
cloudy updrafts at four to ten times that and its cores at 0.6–4 times.

**The grid.** On the held-night boxes the rate stays put as the spacing shrinks. Their cloudy updrafts mix at 2.18,
2.03 and 2.12 × 10⁻³ m⁻¹ by qt (1.52, 1.36 and 1.58 by h) at 6, 3 and 2 km, and at 2.21, 2.15 and 2.22 by qt over the
second days alone, while half their columns' cloud tops lie by 7.9 km on every grid. Their buoyant cloud occupies only
the cloud's topmost level, too thin for a core profile on any grid. They mix faster than the source's night (1.34 ×
10⁻³ m⁻¹), their cloud a single thin deck at 6.2–7.9 km under the held forcing. All keep the same 111 levels,
0.5–0.6 km apart in the cloud layer, so the night's cloud spans four or five of them, and CM1's boundary-layer scheme;
a finer vertical grid stays untested.

**The evidence.** The rates are measurements of the model, the bulk mixing of the cloud that the 6, 3 and 2 km boxes
resolve, from one realization each; the resampling ranges treat 3-hourly snapshots as independent and so understate
the spread. Romps (2010) found that the bulk-plume method gives at most half the entrainment measured directly in LES,
and the night's core rates by qt and h differ threefold, since a core chosen for its buoyancy keeps its h. The caps
are a calculation on mean soundings. The mechanism, a layer of almost no buoyancy between the boundary layer and the
CAPE, is an inference from the caps, the cloud's own buoyancy and mass flux, and the soundings, and BOMEX is an
analogy from the literature. The evening marks a limit of the bulk plume: a fresh plume in the evening's mean air stops
at 5.8–7.9 km at the measured rates, or never becomes buoyant, yet the evening rains 2.97 mm/day from storms the
afternoon set going. Their cores mix 0.07–0.23 × 10⁻³ m⁻¹ by qt between 6.7 and 17 km while rising at 0.7–4.7 m/s,
more slowly than the evening's bulk rate, which the cloud near its base sets.
[The GCM's account](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october) traces why its scheme rains
after midnight all the same.

**The sources and how they were read.** A helper agent retrieved the literature on 25 October, within a bounded search;
the publishers' pages for the AMS and Wiley papers refused it.
- Romps (2010, *J. Atmos. Sci.* 67, 1908–1927), read in full from the author's copy (davidromps.com): BOMEX's
  bulk-plume and direct entrainment, and the factor of two.
- Siebesma and others (2003, *J. Atmos. Sci.* 60, 1201–1219), not read. BOMEX's forcing comes from MicroHH's
  implementation of the case (`cases/bomex/bomex_input.py` on GitHub), and its LES entrainment of 1–3 × 10⁻³ m⁻¹ from
  Siebesma's 2005 GCSS presentation (atmos.washington.edu, `~breth/GCSS/Arrowhead-200506`), both read.
- Siebesma and Cuijpers (1995), not read. Its method is taken as Romps (2010) and Villalba-Pradas and Tapiador (2022,
  *Geosci. Model Dev.* 15, 3447, read) describe it, cores entraining 1.5–2 × 10⁻³ m⁻¹ and the updraft about
  3 × 10⁻³ m⁻¹ in BOMEX, five to ten times Tiedtke's (1989) 3 × 10⁻⁴ m⁻¹ for shallow convection.
- De Rooy and others (2013, *QJRMS* 139) and Derbyshire and others (2004, *QJRMS* 130), not read.

**Why the settled box runs colder aloft** (25 October, at the author's go-ahead; [box_upper_air.py](box_upper_air.py),
[box_upper_air_near_side_sea.json](../results/crm/box_upper_air_near_side_sea.json), 8 seconds, no new runs). Over its
third lunar day the box is 7.7 K colder than Kuo's GCM, the run that forced it, at 42 km and 10.0 K at 56 km. Its
convection holds the upper air on a diluted moist adiabat, which sets most of that gap; the GCM's moister air near the
ground sets most of the rest, and the radiation codes add little.
- **Against each model's own moist adiabat.** Kuo's air at 19–56 km sits 0.3–1.3 K below the undiluted
  pseudo-adiabat of its own lowest level, and the rh 0.7 branch's 1.0–1.6 K above its own. The box's sits 3.0–6.3 K
  below the adiabat of its lowest 0.5 km. A plume held at zero buoyancy while it entrains air of the box's own
  humidity matches the box's air between 0.5 and 1 × 10⁻⁴ m⁻¹ (a calculation on its mean sounding, after the
  zero-buoyancy plume of Singh and O'Gorman 2013, not read): at 1 × 10⁻⁴ its column stays within 1.7 K of the box's
  at 19–56 km, and at 0.5 × 10⁻⁴ the box runs 1.7–2.6 K colder than it. That is 2.4–4.8 per e-folding of pressure,
  below the afternoon cores' bulk rates under 16 km.
- **The air near the ground.** The GCM's lowest level holds 13.5 g/kg of vapour against 11.8 in the box's lowest
  0.5 km, 4.6 kJ/kg more moist static energy, which sets its adiabat 3.1 K warmer at 42 km and 4.0 K at 56 km.
- **Radiation on one column.** On the GCM's air, CM1's clear-sky longwave (its RRTMG through rrtmg_edge's driver) cools
  35–48 km as PlaSim's does, 27.6 against 27.2 W/m² (29.7 and 29.8 with each code's curved-air escape), and 48–63 km 3.8
  W/m² more (2.7 with the escapes). The box absorbs 3.4 W/m² less sunlight at 48–63 km and 3.8 more at 35–48 km.
  PlaSim's water-vapour copy shifts by −3.9 and +4.9 W/m² between the two columns, so that difference follows the box's
  drier air at 48–63 km (0.22 against 0.48 g/kg). On the box's own colder, drier air RRTMG cools 35–63 km 17–18 W/m²
  more than PlaSim, whose fitted correction holds on the GCM's columns; RRTMG agreed within 2% with the correlated-k
  model aloft on the 45° N ring's columns.

| Lunar-day means, W/m² (48–63 km; 35–48 km) | Box | Kuo GCM | rh 0.7 GCM |
|---|---|---|---|
| Temperature | 233.9; 252.9 K | 243.8; 260.8 | 245.1; 262.7 |
| Longwave, sunlight | −24.8, 11.5; −30.5, 18.1 | −24.1, 14.9; −22.9, 14.3 | −27.0, 16.7; −25.0, 16.2 |
| Convection (box: microphysics, resolved transport, turbulence) | 19.1; 10.5 | 20.0; 17.5 | 9.3; 11.2 |
| Forcing or circulation | −6.9; −2.4 | −10.8; −8.9 | 1.0; −2.4 |
| Storage | −1.2; −1.8 | 0; 0 | 0; 0 |

**How the budgets settle.** At 48–63 km the box's convection heats as much as Kuo's, while its radiation cools
4.1 W/m² more and its forcing 3.9 less. On the box's start air its forcing cools that layer by 11.3 W/m², as the GCM's
circulation does (10.8), so the drift weakened it: the mean ascent there cools the box's less stable air less. The
rh 0.7 branch heats the layer by 9.3 W/m² and still sits 1.3 K warmer than Kuo. The convective heating's size follows
the radiation and the circulation, and the temperature follows the parcel to which the convection holds the air. So
the explanation through mixing fits, at 5.3–6.3 K of the gap at 42–56 km, with radiation adding a few W/m² at
48–63 km. The GCM's moister boundary layer gives 3–4 K, which a mixing setting in the scheme would leave in place.
The evidence is one realization of the box, the GCM's single-step heating averaged over the lunar day, and radiation
on lunar-day mean columns: RRTMG on the box's mean column cools 56 km 18% less than the mean of the box's own
clear-sky heating. The sunlight copy leaves out the gas bands and the curved air's twilight and serves for the
difference between columns only.

```sh
/media/projectspace/terluna-research/session-tools/run_column_budget.sh       # the GCM branch, and gamma's years
climate/gcm/.venv/bin/python -m climate.crm.gcm_site A28R_301k_cold_moist 0 5 -5 -45 sea      # about ten seconds
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_rain_law_pbl                 # three minutes
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_design
python3 /media/projectspace/terluna-research/crm-runs/trial_site_boxes/trial.py run|check
/media/projectspace/terluna-research/session-tools/run_site_boxes.sh             # the pair; rerun to resume
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_nudged --from-day 29.5   # the second attempt
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_day3    # the third lunar day, from day 59
/media/projectspace/terluna-research/session-tools/agents/climate/run_sea_day3.sh    # 1.3 hours; rerun to resume
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_day3 --from-day 59
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_near_side_sea_day3 --from-day 59
climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea box_near_side_sea_day3     # its lunar days
climate/gcm/.venv/bin/python -m climate.gcm.sbm_offline           # the new scheme on the boxes' air, 15 seconds
climate/gcm/.venv/bin/python -m climate.crm.box_onset box_near_side_sea_day3 59 88.5    # why the nights stay quiet
climate/gcm/.venv/bin/python -m climate.crm.box_onset box_0e_design 29.5 59               # the land, for contrast
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_night_6km     # and _3km, _2km
climate/gcm/.venv/bin/python /media/projectspace/terluna-research/session-tools/agents/climate-resolution/trial_night.py run|check box_near_side_sea_night_6km
/media/projectspace/terluna-research/session-tools/agents/climate-resolution/run_night.sh 6km   # 4 min; 3km 34, 2km 49
climate/gcm/.venv/bin/python -m climate.crm.box_resolution          # the three grids beside the source's two nights
climate/gcm/.venv/bin/python -m climate.crm.box_entrainment         # the clouds' entrainment and the scheme's nights
climate/gcm/.venv/bin/python -m climate.crm.box_upper_air           # why the box runs colder aloft than the GCM
```

### The sea box with Earth's boundary-layer length (27 October)

CM1's simple boundary layer (`ipbl` 2) mixes with K = l² √max(S² − N²/Pr, small), its length
l = (l_inf⁻² + (κz)⁻²)^−½, and `case_settings` stretches the asymptotic length l_inf by the gravities from Earth's 75 m
to 453 m, as it stretches the Rayleigh layer's time. No experiment here has tested that scaling, and at the top of the
sea's night boundary layer the scheme carries most of the heat and vapour the box exchanges there
([the GCM's account](../gcm/README.md#the-single-column-with-cm1s-own-convection-27-october)). On 27 October the author
approved a run of the sea box's third lunar day with Earth's length, to see whether its deep, evenly moist night layer
rests on the stretched one ([box_length.py](box_length.py),
[box_length_near_side_sea.json](../results/crm/box_length_near_side_sea.json)).

**The case** (`box_near_side_sea_day3_linf75`). It runs box_near_side_sea_day3 again from box_near_side_sea's day-59
restart, sunrise to sunrise over days 59–88.5, on the build `moon_omp_rain_law_pbl` with the day-3 run's 37 wind
windows, inputs, forcing and physics; every input file is the day-3 run's byte for byte except the namelist, which
differs in l_inf (75.0 m) and in what it writes. l_inf is a namelist value, so the build is the day-3 run's. The domain
diagnostics keep the day-3 run's 3-hour cadence, so the two compare output for output; the 3-D output comes every 6
hours without the budgets' fields, as the held nights' does; and the run keeps only its newest two restarts besides the
one it starts from. It writes 18 3-D fields where its source wrote 45, so CM1 puts all its output under a second output
configuration, the files suffixed `.02` and numbered from 1 (cm1.F, `nconfig`); [box_lid.py](box_lid.py) and box_length
read either.

**The runner.** `cm1_run` gains two keys a case may set: `restart_namelist`, entries set over the source's namelist when
a case runs on from a restart, and `keep_restarts`, under which the runner removes, after each segment, the restart sets
the case wrote beyond its newest few, never the one it started from. No earlier case sets either (a test checks it), so
every case already set up or run sets up and runs as before. The new cases record runner `f3ca175b41d3a536`. The held
night with the same length (`box_near_side_sea_night_6km_linf75`, `setup_night` as for box_near_side_sea_night_6km)
differs from its 453 m twin in l_inf alone.

**The smoke test** (`trial_linf75.py` and `trial_linf75_checks.py` in `session-tools/agents/climate-pbl-tests` on the
research drive). A copy of the set-up case ran six model hours from the restart on six threads. CM1 ended normally, its
log showing l_inf 75 m, the second output configuration and 18 3-D fields; with the 37 windows in place it wrote
diagnostics at 59.125 and 59.25 days, the 3-D output at 59.25 and a restart, every field the analyses read present and
finite, and its running domain rain at the restart equal to the day-3 run's. The copy was removed after the checks.

**The held night first** (`box_near_side_sea_night_6km_linf75`). Two Earth days of the held night at 6 km with Earth's
length ran in 3.6 minutes on six threads. Over the second day, beside its 453 m twin:

| Held night at 6 km, second day | l_inf 453 m | 75 m |
|---|---|---|
| Mixed depth; CM1's hpbl, km | 6.48; 6.10 | 6.44; 6.10 |
| Vapour drop 0.31–6.2 km, g/kg | 0.24 | 0.25 |
| θv rise 6.2–8.2 km, K | 1.35 | 1.45 |
| Rain, mm/day | 0.07 | 0.16 |
| Heat through 6.45–6.99 km, W/m²: the scheme; all parts | −9.6; −11.2 | −0.7; −2.2 |
| Vapour through 6.45–6.99 km, W/m²: the scheme; all parts | 11.6; 15.3 | 0.9; 5.3 |
| The scheme's khv there, m²/s | 10.2 | 0.8 |

At Earth's length the scheme's diffusivity at the lid falls thirteenfold and its exchange there by 92%, while the
resolved parts with the implicit diffusion barely change (−1.6 and −1.5 W/m² of heat, 3.8 and 4.3 W/m² of vapour). Over
those two days the layer keeps its depth and its even moisture, and its lid gains 0.1 K. The held night starts from the
453 m box's own air, so it shows how the lid's exchange answers the length at once and leaves the layer's own evolution
to the lunar day.

**The rule, set before the run** (`CRITERIA` in box_length.py). After midnight the run with Earth's length is robust if
its mixed depth lies within 1 km of the stretched run's 6.48 km and its 0.31–6.2 km vapour drop within 1 g/kg of
0.35 g/kg; it depends on the length if the depth falls by more than 1.5 km or the drop grows by more than 1 g/kg; it is
ambiguous otherwise.

**The lunar day.** Both runs over days 59–88.5, read output for output: 59 domain diagnostics in each quarter, and
29–30 3-D outputs at the 6-hourly times both wrote. Each cell gives the stretched length (453 m) first and Earth's
(75 m) second. The depth, drop and rise are read on the GCM's full levels at the site, 0.314, 6.201 and 8.238 km, as
the comparison with the GCM read them; CM1's own boundary-layer depth (hpbl) is the mean of its columns'
bulk-Richardson heights.

| Quarter | Mixed depth, km | Vapour drop, g/kg | θv rise 6.2–8.2 km, K | CM1's hpbl, km | Rain, mm/day |
|---|---|---|---|---|---|
| Sunrise–noon | 3.25; 3.19 | 0.64; 0.51 | 1.28; 1.44 | 3.69; 3.87 | 0.00; 0.00 |
| Noon–sunset | 2.64; 2.78 | 1.08; 1.08 | 0.22; 0.24 | 2.31; 2.60 | 2.27; 2.35 |
| Sunset–midnight | 5.58; 5.41 | 2.14; 2.33 | 1.25; 1.31 | 4.91; 4.85 | 2.96; 3.00 |
| Midnight–sunrise | 6.48; 6.15 | 0.35; 0.53 | 1.44; 1.70 | 6.07; 5.78 | 0.16; 0.28 |

The cloud below 40 km (liquid and ice of at least 0.01 g/kg) in the 3-D outputs, and what each run carries through
6.45–6.99 km (W/m², upward; the scheme's part, the resolved part with the advection's implicit diffusion, and all):

| Quarter | Columns with cloud | Their top: median, 90th percentile, km | Columns with cloud above 16 km | Heat: scheme, resolved, all | Vapour: scheme, resolved, all | khv, m²/s |
|---|---|---|---|---|---|---|
| Sunrise–noon | 0.7%; 1.8% | 9.8, 13.8; 9.4, 9.7 | 0.1%; 0.1% | −1.2, 0.1, −1.2; 0.0, 0.5, 0.5 | 0.8, 0.4, 1.2; −0.2, −1.1, −1.4 | 13.3; 0.7 |
| Noon–sunset | 10.2%; 9.5% | 17.2, 31.4; 18.7, 29.5 | 6.8%; 6.5% | −0.5, 0.0, −0.6; −0.1, 1.6, 1.5 | 3.3, 33.1, 36.4; 0.4, 31.7, 32.1 | 31.7; 1.5 |
| Sunset–midnight | 10.1%; 11.4% | 10.8, 20.8; 8.8, 17.8 | 3.0%; 1.3% | −7.2, −2.3, −9.5; −0.4, −2.3, −2.7 | 19.0, 69.1, 88.0; 1.5, 73.1, 74.6 | 10.0; 0.5 |
| Midnight–sunrise | 21.9%; 45.8% | 7.4, 10.8; 8.9, 12.5 | 0.7%; 8.2% | −10.1, −0.6, −10.7; −0.7, −0.1, −0.8 | 14.1, 8.5, 22.6; 1.2, 9.6, 10.8 | 11.3; 0.8 |

**Under the rule the night layer is robust.** After midnight, with Earth's length, the layer stands 6.15 km deep
against 6.48 and its vapour drops 0.53 g/kg against 0.35, inside both margins.

**What the length does.** At Earth's length the scheme's diffusivity at the lid falls 13–21-fold in every quarter, and
its exchange there with it: after midnight it carries 0.75 W/m² of heat down and 1.2 W/m² of vapour up, against 10.1
and 14.1. The resolved motion does not take that exchange over. Through the lid after midnight the run with Earth's
length moves 0.8 W/m² of heat against 10.7 and 10.8 W/m² of vapour against 22.6. The layer below changes little: its
lid is 0.26 K stronger, its air 0.2–0.4 g/kg moister up to 6.2 km and 0.5 g/kg drier at 8.2 km. The night's
cloud changes most: it covers 46% of the columns after midnight against 22%, its median top stands at 8.9 km against
7.4, 8.2% of the columns reach above 16 km against 0.7%, and it rains 0.28 mm/day against 0.16, a fifth of the GCM's
1.45 with the new convection. By day the two runs differ little: in the afternoon the layer's depth, its drop and the
rain lie within 0.14 km, 0.01 g/kg and 0.08 mm/day of each other, and the evening's storms rain alike, 3.00 against
2.96 mm/day, their tops lower at Earth's length (median 8.8 km against 10.8).

So the sea's deep, evenly moist night layer holds at Earth's length, and the scheme's mixing at the layer's top bears
mostly on the night's cloud. In CM1 the exchange the scheme carries at the lid can almost vanish while the layer's depth
moves 0.33 km and its lid 0.26 K. In the single column of
[the GCM's account](../gcm/README.md#the-single-column-with-cm1s-own-convection-27-october), where CM1's convective
tendencies are prescribed from the 453 m run and cannot adjust, the absence of the same exchange doubled the lid; the
column's doubled lid therefore overstates what that exchange does in a model whose convection answers its absence (an
interpretation).

**Evidence and limits.** The two runs are one realization each, started from one restart, so their differences mix the
length's effect with two realizations drifting apart; the stretched run's second and third lunar days differ in their
afternoon rain by 87%, so rain differences of a few tenths of a mm/day lie within that spread. The fluxes are CM1's own
diagnostics, each output one time step's mean over the box's columns, the scheme's part its `thfs` and `qvfs`. Both
runs use CM1's scheme on a 6 km grid, so the test measures the box's sensitivity to that scheme's length within CM1; no
observation stands for either length. The cloud tops come from the 3-D outputs both runs wrote every 6 hours, 118 of
them, against 236 diagnostics.

**Cost.** The lunar day ran in 51 minutes on six threads, in 30 segments, and its folder holds 1.1 GB with its two
kept restarts and the starting one; the held night ran in 3.6 minutes and holds 0.15 GB; the smoke test's copy was
removed after its checks. The day-3 run took 1.3 hours on three threads and holds 5.6 GB with every restart and the
budgets' 3-D fields.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_night_6km_linf75   # and _day3_linf75
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_near_side_sea_night_6km_linf75 --hours 0.5 --threads 6
PYTHONPATH=. climate/gcm/.venv/bin/python /media/projectspace/terluna-research/session-tools/agents/climate-pbl-tests/trial_linf75.py
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_near_side_sea_day3_linf75 --hours 3 --threads 6   # 51 min
climate/gcm/.venv/bin/python -m climate.crm.box_length       # both lengths by quarter, and the held nights; seconds
climate/gcm/.venv/bin/python -m climate.crm.box_lid          # the held nights' lid split by part
```

### The sea box forced by the fifteen-level GCM (27 October)

Every sea box above took its forcing from the Kuo run on ten levels (`A28R_301k_cold_moist`), while the night rain it
was set against came, from 22 October on, from GCM runs with the Betts–Miller scheme and, from 26 October, on fifteen
levels, whose circulation differs from Kuo's. On 27 October the author chose to test whether CM1's quiet sea nights
survive the fifteen-level GCM's own circulation. The case `box_near_side_sea_l15` is box_near_side_sea with one change,
the run whose column forces it, read by a rule set before it ran.

**The forcing run.** The fifteen-level baseline (`A28R_301k_cold_sbm_rh07_deep200_l15cap`) had two years. It was run
again from its own start on the shared ExoPlaSim source as it now stands, for four years
(`A28R_301k_cold_sbm_rh07_deep200_l15cap_cwsource`; `runs.sh` and `settling.py` in `session-tools/option3` on the
research drive). Its first two years repeat the baseline's in all 61 output fields and both restarts, so the source's
later patches, switched off, leave the model as it was. Over the four years it settles: the ground cools from 298.10 to
297.53, 297.15 and 296.87 K, and the imbalance at the top shrinks from −4.12 to −3.27, −2.63 and −1.90 W/m². A line
through the imbalance against the ground reaches zero at 295.7 K, and an exponential approach fitted to the ground
levels off at 296.3 K with an e-folding time of 2.7 years. The sea site's rain after midnight, 2.31, 1.85, 2.33 and
2.62 mm/day in years 0–3, does not follow that drift. Years 2 and 3 fall inside 1.4–2.8 mm/day, the first two years'
range widened by their difference, which was the rule set before the run, so the box takes all four years.

**What the forcing brings** (`forcing_audit.py` beside them), against the Kuo run's six years at the same six cells:
- **The surface and the reference column.** The sea stands at 301.02 K against 302.03, the column 0.7–1.3 K cooler in
  potential temperature below 10 km and 1.4–1.5 g/kg drier below 7 km.
- **The mean vertical wind.** It rises 2.6 mm/s at 1 km against 3.5 and 2.1 against 4.0 at 4 km, and it sinks from
  about 5.5 km, 1.3 mm/s at 7 km, where Kuo's still rises.
- **The held air.** The GCM's convection forms its rain up to 55 km against 72, so the box's air is held from 71.7 km,
  fully from 92.4 km, against 93.1 and 129.3.
- **The circulation by local time.** After midnight it warms the column by 88 W/m² against 57 and dries it by
  4.3 mm/day against 3.3; through the afternoon it moistens it by 2.3 mm/day against 4.5. Over the lunar day it takes
  2.0 mm/day of vapour from the site's column against 0.8, as the GCM rains 1.78 mm/day there against 2.72 and
  evaporates 3.75 against 3.50.

Every other input follows from the fifteen-level run by the setup the Kuo-forced box had. The rain limits differ only
in the breakup's reference density, 1.1991 against 1.199 kg/m³.

**The runner.** A case run on from a case that was itself run on from a restart now takes its wind windows from the
case the chain began from (`window_root`). A run-on case numbers its windows from its own restart while
`restart_windows` counts phases from the run's start, so the fourth lunar day's windows would otherwise have come one
window, 0.82 days, out of phase. Every case set up before ran on from a case that began its chain, so its windows are
unchanged; a test checks both. `box_days` reads a fourth lunar day, 88.5–118 days. The first two cases record runner
`73493ca056a4fa5e`, the fourth day `c8d7968f708640a0`.

**The trial** (`trial_option3/trial.py` in `crm-runs`, after the site boxes' trial). Three model hours ran in 28 s on
two threads. Every field the analyses read was present and finite, the namelist held the air from 71.7 and 92.4 km,
and the boundary layer passed on 0.999 times the sea's evaporation.

**The runs.** Two lunar days from sunrise ran in 1.8 hours on four threads, the third from the day-59 restart in 1.0
hour and the fourth from the day-88.5 restart in 1.1 hours (`box_runs.sh` and `box_day4.sh` beside the GCM's). Rain by
quarter, in mm/day ([box_days_box_near_side_sea_l15.json](../results/crm/box_days_box_near_side_sea_l15.json)):

| | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| Fifteen-level forcing, second lunar day | 0.00 | 0.95 | 1.03 | 0.50 | 0.62 |
| Fifteen-level forcing, third lunar day | 0.00 | 1.16 | 1.20 | 1.62 | 1.00 |
| Fifteen-level forcing, fourth lunar day | 0.00 | 1.01 | 1.06 | 1.88 | 0.99 |
| Kuo forcing, third lunar day | 0.00 | 2.23 | 2.97 | 0.14 | 1.33 |
| The fifteen-level GCM at the site, four years | 0.19 | 1.64 | 2.96 | 2.28 | 1.77 |

**The rule, set before the run**, on the box's rain after midnight over its third lunar day, against the GCM's
2.28 mm/day: at most 0.5 mm/day would have left CM1's quiet night standing under the fifteen-level GCM's own
circulation, and the GCM's night rain to its moist physics; at least 1.0 means that the quiet night came with the Kuo
forcing. The third day rains 1.62 mm/day after midnight, and the fourth, run because the box was still changing,
1.88.

**The night after midnight** ([box_onset_box_near_side_sea_l15_day3.json](../results/crm/box_onset_box_near_side_sea_l15_day3.json)
and [_day4](../results/crm/box_onset_box_near_side_sea_l15_day4.json), against box_near_side_sea_day3's), over the third
and fourth lunar days:
- **A shallower, unstable night.** The boundary layer is 3.6 and 4.5 km deep against 6.1 under the Kuo forcing (the
  GCM's 4.3). The condensation and free-convection levels stand at 3.3–3.7 and 3.7–3.9 km against 5.8 and 6.6, and
  CM1's CAPE at 1,770–1,880 J/kg against 1,265, with 1.0–1.3 J/kg of inhibition against 1.7.
- **Strong updrafts in a few columns.** About 2% of the columns have tops above 16 km under either forcing, but here
  their updrafts are strong: 6.5–6.7 m/s at 18–19 km at the most, and 2.4–2.7 m/s at the 99.9th percentile between 16
  and 60 km against 0.3.
- **The cloud layer.** At 5.2–11.8 km on the night budget's bands (`read_box4.py` beside the runs) it condenses 1.86
  and 1.71 mm/day and re-evaporates 83% of that on both days, against 1.14 mm/day and 81%. The prescribed circulation
  warms the layer by 31 W/m² against 18, as the GCM's own circulation does there (28 W/m²). Its heat budget leaves
  9–11 W/m² unclosed, against under 1 for the Kuo-forced box.

**The evening.** Before midnight the box rains 1.20 and 1.06 mm/day against the GCM's 2.96. Its storms are held down
by 18–26 J/kg of inhibition, against 3 under the Kuo forcing, with a free-convection level at 9.1–9.4 km against 5.3.
Its rain therefore peaks after midnight, where the GCM's peaks in the evening.

**Budgets and drift.** Over the third and fourth lunar days the water budget closes to 0.07 and 0.04 mm/day, and the
boundary layer passes on 0.998 times the sea's evaporation of 2.76 and 3.03 mm/day
([box_budget_box_near_side_sea_l15_day3.json](../results/crm/box_budget_box_near_side_sea_l15_day3.json) and
[_day4](../results/crm/box_budget_box_near_side_sea_l15_day4.json)). The forcing takes 2.46–2.48 mm/day from the column
against the GCM's own circulation's 1.99, a gap of about 0.5 against the Kuo-forced box's 1.26; after midnight it
takes 4.9–5.0 mm/day against 4.3. 61–62% of the box's falling water evaporates, gamma 0.027 on PlaSim's pooled shares,
as under the Kuo forcing. The box is still settling. Between the third day's sunrises its domain-mean air cooled by
0.4 K at 3 km, 0.8 K at 6 km, 1.4 K at 10 km, 0.8 K at 40 km and 1.0 K at 60 km; over the fourth by 0.3, 0.5, 0.5,
0.7 and 0.4 K, its vapour falling 0.5–0.6 g/kg below 6 km.

**Why the evening parts** (`evening.py` beside the runs, taken up once the fourth day had shown the night holding). The
same parcel, from the lowest 0.5 km mixed and undiluted, was raised through both models' columns: the box's domain means
at each 3-hourly output and the GCM's own 3-day means at the six cells, laid on the box's levels by height. PlaSim's
Betts–Miller scheme, with the run's settings (reference humidity 0.7, no mixing, no inhibition limit, the 200 hPa depth
test and the shallow cap), was applied at once to both on the run's fifteen levels. By quarter, the box's third and
fourth lunar days and the GCM's four years, the GCM's rain here from its 3-day means binned by the site's hour angle:

| | CM1's rain, mm/day | Parcel on CM1's air: CAPE; inhibition, J/kg | The scheme on CM1's air: deep share; rain at once, mm/day | The GCM's rain (convective), mm/day | Parcel on the GCM's columns: CAPE; inhibition, J/kg |
|---|---|---|---|---|---|
| Noon–sunset | 1.16; 1.01 | 1,240 and 750; 56 and 33 | 36% and 25%; 54 and 24 | 1.71 (1.43) | 118; 51 |
| Sunset–midnight | 1.20; 1.06 | 940 and 980; 27 and 20 | 37% and 36%; 109 and 73 | 3.02 (2.81) | 27; 23 |
| Midnight–sunrise | 1.65; 1.91 | 2,230 and 2,050; 1.2 and 0.7 | none on either day | 2.18 (1.49) | 28; 20 |

- **The GCM's columns hold little CAPE,** 27–118 J/kg on its 3-day means against 750–2,230 J/kg in CM1's air, as a
  Betts–Miller scheme keeps its column near its reference.
- **On CM1's own air the scheme convects in the afternoon and evening.** It goes deep at a quarter to a third of those
  outputs, raining tens of mm/day at once, where CM1 rains about 1 mm/day under 20–56 J/kg of inhibition. Split by 15°
  of hour angle ([below](#the-late-evening-on-finer-grids-and-where-the-models-part-28-october)), those outputs fall
  where CM1's own storms rain, at 60–120°, and none falls in CM1's lull from 135° to midnight.
- **After midnight it never convects deep.** There CM1 convects through about 1 J/kg and rains 1.65–1.91 mm/day, while
  the scheme stays shallow at every output: on the third day its cloud was too thin for the depth test at 26 outputs
  and the column drier than its reference at 33.

On the same air, then, the scheme keeps resolved convection's timing through the evening and leaves it after
midnight. Read by quarter alone, these outputs were first taken to show the scheme convecting through the inhibition
that held CM1's storms down; the split by hour angle corrects that reading, and puts the evening's parting in the two
models' air. The GCM's night rain comes from its own columns, held near its reference, so the two models' nights come
close in amount by different routes.

**What this changes.** The sea's night gap that the accounts since 24 October chased came largely with the Kuo run's
circulation. Under the fifteen-level GCM's own circulation CM1 rains through the night at 71–82% of the GCM's rate, so
at this site the fifteen-level GCM's night rain is close to what resolved convection makes of the same large-scale
conditions. Several readings that rested on the Kuo-forced box are qualified:
- CM1's 6 km evenly moist night layer, against which
  [the boundary-layer comparison](../gcm/README.md#plasims-boundary-layer-mixing-at-lunar-gravity-27-october) and the
  nonlocal scheme's test set PlaSim, belongs to the Kuo-forced box; under the fifteen-level forcing the layer is
  3.6–4.5 km deep.
- The night water budget's pointer to PlaSim's night-cloud physics
  ([the GCM's account](../gcm/README.md#the-seas-night-water-by-process-and-height-27-october)) and the cloud water's
  premise compared PlaSim's night with that box's.
- On 24 October the circulation's night alone, warmer aloft and drier under the new GCM, was read as showing that the
  Kuo forcing did not explain the quiet nights. The day's forcing turns out to matter more: under the fifteen-level
  circulation the afternoon is fed less, the evening's storms stay weaker, and the night keeps its CAPE (an
  interpretation of the quarters above).

The Kuo GCM's own comparison stands: under its circulation CM1 rains 0.14 mm/day after midnight where Kuo rains 2.9.
The models now part in the late evening, from 135° to midnight, where the GCM's low air is warmer and moister than the
box's ([28 October](#the-late-evening-on-finer-grids-and-where-the-models-part-28-october)). Open, at the author's
word: settling the fifteen-level configuration with the evening's difference stated, its longwave refit and a 15-year
run; a physical account of the evening before any GCM test, taken up on 28 October; and the land box forced by the
same GCM.

**Evidence and limits.** One realization on the 6 km grid, forced by four years of a GCM still settling, and the box
itself still settling over its fourth lunar day. The prescribed circulation includes the GCM's response to its own
convection, and a box forced by it cannot judge whether that circulation is the Moon's; it shows what resolved
convection does under it. The cloud layer's heat budget leaves 9–11 W/m² unclosed after midnight. The scheme applied at
once to another model's air meets that air away from its own equilibrium, so its rain at once measures how far the air
stands from the reference ([the offline test of 24 October](../gcm/README.md#plasims-convection-at-the-moons-long-day-21-october)),
and the GCM's 3-day means smooth the instants the scheme fires on. The GCM's quarter rain comes from night_water's
fitted 3-day means at its six cells, CM1's from its domain diagnostics every 3 hours.

**Cost.** The GCM's four years ran in 31 minutes on four cores and hold 0.9 GB. The box's four lunar days ran in
3.9 hours on four threads and hold 23 GB with every restart; the trial's copy holds 0.1 GB. The read-outs take
seconds to minutes.

```sh
S=/media/projectspace/terluna-research/session-tools/option3
$S/runs.sh l15_four_years                                                      # the GCM, 31 minutes
climate/gcm/.venv/bin/python $S/settling.py                                    # the years' rule
climate/gcm/.venv/bin/python $S/forcing_audit.py                               # the forcing beside the Kuo run's
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_l15     # _day3 and _day4 in turn
$S/box_runs.sh                                                                 # two lunar days, then the third
$S/box_day4.sh                                                                 # the fourth
climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea_l15 box_near_side_sea_l15_day3 box_near_side_sea_l15_day4
climate/gcm/.venv/bin/python -m climate.crm.box_onset box_near_side_sea_l15_day4 88.5 118     # and _day3 59 88.5
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_l15_day4 --from-day 88.5
climate/gcm/.venv/bin/python $S/read_box4.py --no-products                     # the rule and the night's budget
climate/gcm/.venv/bin/python $S/evening.py                                     # why the evening parts
```

### The late evening on finer grids, and where the models part (28 October)

Under the fifteen-level forcing the two models part in the evening (above). Split by 15° of hour angle
(`evening_bins.py` beside the runs), the evening holds two spans. From 60° to 120°, through sunset, both rain, the GCM
three times as much at first and about as much by the end. From 135° to midnight the box rains 0.1–0.2 mm/day where
the GCM rains 2.3–3.6:

| Hour angle | CM1's rain, third; fourth lunar day, mm/day | CM1's inhibition, J/kg | PlaSim's scheme on CM1's air: outputs convecting deep | The GCM's rain (convective), mm/day |
|---|---|---|---|---|
| 60–75° | 0.88; 0.83 | 52; 33 | 100%; 60% | 2.46 (2.09) |
| 75–90° | 2.12; 1.93 | 61; 43 | 100%; 100% | 4.90 (4.20) |
| 90–105° | 3.19; 2.68 | 54; 39 | 100%; 100% | 4.67 (4.27) |
| 105–120° | 2.87; 2.45 | 44; 31 | 100%; 100% | 2.73 (2.63) |
| 120–135° | 0.64; 0.72 | 30; 21 | 20%; 10% | 1.89 (1.85) |
| 135–150° | 0.20; 0.21 | 20; 13 | none | 2.82 (2.74) |
| 150–165° | 0.11; 0.09 | 9.1; 4.8 | none | 2.29 (2.00) |
| 165–180° | 0.10; 0.12 | 1.2; 0.4 | none | 3.64 (3.29) |

Applied at once to the box's air with the run's settings, PlaSim's scheme convects deep where the box's own storms
rain, while its columns hold 30–60 J/kg of inhibition on average, and stays shallow in the lull. The quarters above
averaged the two spans together, and their reading, that the scheme convects through the inhibition that holds the
box's storms down, does not survive the split: on the box's air the scheme keeps resolved convection's timing through
the evening and leaves it only after midnight, where it stays shallow while the box rains. From 135° to midnight,
then, the models part because their air differs. Two questions follow: whether the box's lull comes from its grid, and
how the GCM's late-evening air differs from the box's.

**The lull on finer grids** (cases `box_near_side_sea_l15_evening_6km`, `_3km` and `_2km`;
[box_resolution.py](box_resolution.py) with `--test evening`;
[box_resolution_near_side_sea_l15_evening.json](../results/crm/box_resolution_near_side_sea_l15_evening.json)). In
the lull the box's boundary layer is 5.5 km deep, nearly as deep as its 6 km columns are wide, and the free-convection
level of its mean sounding stands near 9 km. Coarse grids delay convection's onset under such a layer (Petch and
others 2002), so three boxes put the lull to the test the night passed on 25 October (above).
- **The start.** Each box takes box_near_side_sea_l15_day4's air at day 107.625 (its output 862), 142° of hour angle
  in the fourth lunar day's lull, averaged over its columns, with CM1's random perturbations of up to 0.25 K; the mean
  air carries no cloud. Over the next Earth day the source rained 0.11 mm/day under 1,000 J/kg of CM1's CAPE and
  10 J/kg of inhibition.
- **The forcing, held at that hour for two Earth days** (`setup_night`, as for the night). The clock that times the Sun
  and the day-night table stops at 142.03°, with the Sun 51.75° below the horizon, and the held table dries the column
  by 1.18 mm/day. The wind is nudged toward the window of 140–150° for the whole run. The hold on the air above the
  GCM's convection, the mean vertical wind, the sea at 301.02 K, the rain law and limits, the build and its physics are
  the source's. The cases record runner `58208d6dfe0b435b`.
- **The grids** are the night's: 6 km on 32 × 32 columns (192 km), the control; 3 km on 64 × 64 (192 km); and 2 km on
  64 × 64 (128 km), all on 111 levels to 150 km with CM1's boundary-layer scheme, their adaptive steps near 75, 39 and
  26 s. They ran in turn in 4, 30 and 43 minutes on four threads and wrote 1.33 GB with the night's output settings.
  They had no smoke tests of their own; their set-up is the night's, whose smoke tests passed.

**The rule, set before the runs** (`EVENING_CONTROL` and `EVENING_READING`). Over its second day the control must
reproduce the source's lull: rain under 0.4 mm/day, twice the source's at 135–150°; the strongest updraft below 60 km
averaging under 3 m/s; and the 99.9th percentile of w below 16 km under 2 m/s in every snapshot. The source carries the
sunset storms' ice above 40 km over much of its domain, which a box started from mean air lacks, so no share of cloud
tops enters. Otherwise the method fails and the finer grids go unread. A grid convects when it rains at least
0.93 mm/day, a third of the GCM's 2.8 at 135–150°, with its strongest updraft averaging 5 m/s or w's 99.9th percentile
at 16–60 km averaging 1.5 m/s, below the 2.4–2.7 m/s of the source's storms after midnight. Should the finer grids
convect, the lull comes from the grid; should they stay quiet, it holds across the grids, and CM1 does not support the
GCM's late-evening rain.

**The control reproduces the lull, and the finer grids stay as quiet.** Over each box's second day, beside the source
box over the same hours of its run (48 diagnostics and 4 snapshots in each box, 8 of each in the source):

| Second day | Source box, days 108.625–109.625 | 6 km, the control | 3 km | 2 km |
|---|---|---|---|---|
| Rain | 0.09 mm/day | 0.14 | 0.16 | 0.17 |
| Strongest updraft below 60 km: mean; highest | 1.00; 1.15 m/s | 0.92; 1.05 | 1.39; 1.64 | 1.54; 1.86 |
| w's 99.9th percentile below 16 km; at 16–60 km | 0.82; 0.12 m/s | 0.74; 0.15 | 1.04; 0.19 | 1.15; 0.22 |
| CM1's CAPE; inhibition | 1,309; 3.3 J/kg | 1,232; 5.5 | 1,204; 6.2 | 1,188; 6.7 |
| The 0.5 km parcel of the mean sounding: CAPE; inhibition; CAPE at eps 8 | 1,464; 3.2; 67 J/kg | 1,372; 6.2; 59 | 1,340; 7.1; 59 | 1,322; 7.6; 58 |
| That parcel's free-convection level | 7.7 km | 8.7 | 8.9 | 9.0 |
| Columns whose own free-convection level lies in the boundary layer | 30% | 17% | 12% | 9% |
| Columns with cloud at 16–40 km; mean top of the cloud below 40 km | 63%; 24.7 km | 63%; 25.4 | 61%; 25.3 | 61%; 25.4 |
| Boundary layer | 5.5 km | 5.5 | 5.5 | 5.5 |

- **Stronger eddies, no deep cells.** The 99.9th percentile of w below 16 km rises from 0.74 m/s at 6 km to 1.04 at
  3 km and 1.15 at 2 km, and the strongest updraft from 0.92 to 1.54 m/s, in a boundary layer 5.5 km deep on every
  grid. Aloft the 99.9th percentile stays at 0.15–0.22 m/s, under a tenth of what the source's storms reach after
  midnight.
  The rain rises from 0.14 to 0.17 mm/day as the grid refines. Between the second day's halves the strongest updraft
  moves by under 0.05 m/s on every grid.
- **The air recharges at one pace on every grid.** Held at that hour, the sea evaporates 2.93–2.98 mm/day into each
  box and the boundary layer passes all of it on; the held table takes 1.18 mm/day and the rain 0.14–0.17, so the
  vapour below 6.5 km rises by about 0.45 g/kg a day. CM1's CAPE climbs from 890–940 J/kg over the first half-day to
  1,250–1,300 over the last, and its inhibition falls from 11–13 to 4–5 J/kg, alike on the three grids. The source,
  whose forcing moves on toward midnight, reaches 1.9 J/kg over the same hours and rains after midnight.
- **Cloud aloft.** Cloud with tops at 22–25 km forms from the mean air, which carried none, over 89–90% of the columns
  in the first day and covers 61–63% over the second, as in the source; the updrafts there are the weak ones above. No
  column holds cloud above 40 km.
- **Drift and water.** Over the two days the air at 1 km warms 0.18–0.22 K and at 16 km 0.25–0.30 K, and at 40 km it
  cools 0.37–0.38 K on every grid, against −0.02, +0.17 and −0.42 K in the source over the same hours. The vapour
  budget's rest over the second day is −0.06 to −0.08 mm/day, against 0.01 in the source.

So under the rule set before the runs none of the three convect, and the box's quiet late evening holds up across the
grids. At 3 and 2 km the boxes resolve stronger eddies in the same boundary layer, recharge at the same pace and form
no deep cell in two days. On the box's own air CM1 does not support the GCM's late-evening rain.

**Where the two models' air parts** (`late_evening_columns.py`, `low_air.py` and `sea_surface.py` beside the runs).
The GCM's 3-day means at the six cells were laid on the box's heights, as above, and set beside the box's domain means
over its third and fourth lunar days.
- **The late evening** (135° to midnight: 58 of the box's outputs, 360 of the GCM's means). The GCM's air is warmer
  through the whole column, by 1.0–1.7 K below 8 km, 2.3 K at 10 km and 5.0–7.6 K at 16–30 km. Its lowest 4 km hold
  0.25–0.97 g/kg more vapour, the most near the ground (11.92 against 10.95 g/kg at 0.3 km), and above 10 km it stands
  at 67–70% of saturation, where the box's air stands at 83–98%.
- **PlaSim's scheme on each mean column**, with the run's settings on its fifteen levels. On the GCM's, the parcel from
  the lowest level passes 20 J/kg of inhibition and is then buoyant without a break from 5.4 to 55 km, so the scheme
  convects deep. On the box's, the parcel is buoyant only at the level 2.5 km up, inside the boundary layer, and
  nowhere above its condensation level, so the scheme stays shallow. The box's own 0.5 km parcel holds 1,100–1,500 J/kg
  above a free-convection level at 7.7–9.3 km, 2–4 km above the boundary layer's top, and mixing at eps 8 leaves it
  50–70 J/kg. The scheme goes deep on only 7% of the GCM's individual 3-day means, whose averaging smooths the instants
  it fires on, while the GCM rains 2.69 mm/day convectively at those hours.
- **Through the lunar day**, by quarter, with the GCM's wind at 10 m taken from the 3-day means of its components,
  which understates its speed, and the box's as the mean of its columns' speeds:

| Quarter | Air at 0.3 km: GCM; box | The sea: GCM; box | Evaporation: GCM; box | Wind at 10 m: GCM; box | The box's mean wind at 0.3 km; its target's |
|---|---|---|---|---|---|
| Sunrise–noon | 298.83 K, 12.36 g/kg; 298.29, 11.17 | 300.77; 301.02 K | 3.61; 2.74 mm/day | 2.75; 1.23 m/s | 1.42; 3.23 m/s |
| Noon–sunset | 300.24, 12.11; 299.62, 9.59 | 301.32; 301.02 | 4.54; 3.03 | 2.81; 1.57 | 1.64; 3.10 |
| Sunset–midnight | 299.66, 11.93; 298.06, 10.34 | 301.16; 301.02 | 3.50; 3.22 | 2.16; 1.41 | 1.48; 2.21 |
| Midnight–sunrise | 298.72, 12.29; 297.35, 12.33 | 300.83; 301.02 | 3.37; 2.60 | 2.40; 1.27 | 1.37; 2.49 |

  The GCM's low air is 0.5–1.6 K warmer in every quarter and 1.2–2.5 g/kg moister from sunrise to midnight, its vapour
  level with the box's after midnight. Its sea, within 0.3 K of the box's held 301.02 K, evaporates 9–50% more, most in
  the afternoon, under a wind at 10 m 1.5–2.2 times the box's.
- **The box's wind falls short of its forcing.** Its domain-mean wind at 0.3 km runs at 44–67% of its target's speed
  through most of the lunar day and reaches 92% only in the late evening, where the target is weak and steady
  (1.6–1.7 m/s at 0.3–1 km). The nudging relaxes the box's mean wind over 3 days, a time scale carried over from the
  rings, where it held the air to a fixed profile. Here the target changes with each window, every 0.82 days, as the
  GCM's wind turns with the Sun, faster than such a relaxation follows (an interpretation; the box's momentum budget is
  not read here). The Kuo-forced box was nudged the same way, unmeasured here.

**What it shows** (an interpretation of the runs above). The box's lull is a recharge. Its sunset storms end with
20–30 J/kg of inhibition in its columns and a free-convection level near 9 km over a boundary layer 5.5 km deep; the
sea then refills the low air, the inhibition decays to nil by midnight, and after midnight the box convects again from
a free-convection level at 3.7–3.9 km (above). On the held hour the recharge runs at one pace on 6, 3 and 2 km, so the
grid does not set the lull's length. The GCM enters the evening with warmer, moister low air, and on that air its
scheme keeps convecting through the late evening. Its late-evening rain thus rests on low air that CM1, under the same
circulation, does not build. Two sources of that difference stand out, neither measured yet. On PlaSim's side, its
sea's transfer uses Earth's free-convection constant, which on 26 October accounted for about a quarter of its extra
evaporation, and its wind is stronger. On the box's side, a wind at about half the GCM's through the morning and
afternoon may leave its sea evaporating less, and its air drier, than the GCM's circulation implies. Which side
dominates decides whether the late-evening rain is PlaSim's alone.

**What it changes.**
- The reading of 27 October, that PlaSim's scheme on the box's air convects through the inhibition that holds the box's
  storms down, is withdrawn (corrected above).
- Under the rule set before the runs, CM1 on 6, 3 and 2 km grids does not support the GCM's rain from 135° to midnight
  on the box's own air: that rain leans toward PlaSim's making, pending the account of the low air.
- The box's short wind qualifies every comparison of the sea's evaporation with the GCM's, among them the reading of 26
  October that a lunar free-convection constant removes a quarter of the GCM's extra evaporation, an excess measured
  against the Kuo-forced box, nudged the same way.

Open, at the author's word: a box whose wind follows the GCM's, nudged fast enough for windows of 0.82 days, would show
how much of the low air's difference is the box's; the free-convection constant over the sea is the author's call
(26 October); the items open above stand.

**Evidence and limits.** Each grid is one realization over two Earth days, started from one hour's mean air with that
hour's forcing held; the 2 km box is 128 km wide against the others' 192, and CM1's boundary-layer scheme stays on at
every spacing. The held window is one where the box's wind nearly meets its target, but the air it starts from carries
the day's shortfall. The account sets the GCM's 3-day means, which smooth what its scheme fires on, against the box's
3-hourly domain means, and the GCM's wind taken from its components' means understates its speed, so the ratio of the
winds is a lower bound. The scheme applied at once to another model's air measures how far that air stands from its
reference, not a rain rate. The quarters' vapour is a mixing ratio, the late evening's profile specific humidity. The
box was still settling over its fourth lunar day (above).

**Cost.** The three boxes ran in 76 minutes on four threads and wrote 1.33 GB; each read-out takes seconds to two
minutes.

```sh
S=/media/projectspace/terluna-research/session-tools/option3
climate/gcm/.venv/bin/python $S/evening_bins.py                                # the evening by 15 degrees
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_l15_evening_6km   # and _3km, _2km
$S/run_evening.sh                                                              # 6, 3 and 2 km in turn, 76 minutes
climate/gcm/.venv/bin/python -m climate.crm.box_resolution --test evening      # the rule's reading
climate/gcm/.venv/bin/python $S/late_evening_columns.py                        # the late evening's air in both
climate/gcm/.venv/bin/python $S/low_air.py                                     # the air at 0.3 km by quarter
climate/gcm/.venv/bin/python $S/sea_surface.py                                 # the sea, its fluxes, the box's wind
```

### The sea box with its wind following the GCM's, and the gust over water at lunar gravity (28 October)

The late evening's account (above) found two things in the box's own set-up that bear on its evaporation, and so on
the low air where the two models part. Its mean wind ran at 44–67% of the speed its forcing aims for through most of
the lunar day. And over water its surface layer adds a gust with no gravity in it. Two boxes test them side by side,
each the box's third and fourth lunar days again from its day-59 restart, read by rules set before they ran
(`rules.txt` beside the runs, parts A and D).

**The wind's own time scale** (`WIND_TAU_PATCH`, the build `moon_omp_rain_law_pbl_wind`). CM1 relaxes the domain-mean
wind and the held air toward their targets on one time scale, `lsnudge_tau`, 3 days in every case, a value the rings
took for a profile held fixed. The boxes' wind targets change with each window of local time, every 0.82 days, faster
than a 3-day relaxation follows. With `var10` above 0 the build relaxes the wind alone on `var10` seconds, the held air
keeping `lsnudge_tau`; with `var10` at CM1's default of 0 it runs as `moon_omp_rain_law_pbl` does. Three-hour trials
from the day-59 restart (`trial.py` beside the runs) checked it: with `var10` 0 the new build repeats the old one bit
for bit in all 16 files the trial writes, and with `var10` at 7200 s the box's mean wind at 0.3 km rose from 1.44 to
2.49 m/s within the three hours, against a target of 3.24, and its evaporation from 2.66 to 3.15 mm/day. The boxes take
2 hours, a tenth of a window, the time scale on which cloud-resolving runs forced by observations relax their mean wind
(Khairoutdinov, Randall and DeMott 2005).

**The gust over water** (`GUST_PATCH`, the build `moon_omp_rain_law_pbl_wind_gust`). CM1's surface layer (`sfclay`,
`sfcmodel` 2) adds a gust to the wind in its bulk formulas: over land Beljaars's (1995) convective velocity of the
boundary layer's eddies, vconvc (g/T_s · h · B)^1/3 with vconvc 1, at the build's gravity; over water the old MM5
formula √Δθ_v m/s, with no gravity in it. Over the lunar sea, whose air stands 3–4 K below the sea in virtual
temperature, that gust was 1.6–2.1 m/s and raised the wind in the bulk formulas 1.6–1.8 times
([the GCM's account](../gcm/README.md#the-free-convection-constant-over-four-years-and-the-two-models-light-wind-exchange-28-october)).
The build takes Beljaars's branch over water too; the COARE algorithm's gust over Earth's ocean has the same form
(Fairall and others 2003). In a three-hour trial the gust fell from 2.0 to 0.61 m/s, the surface layer's wind from 3.18
to 2.62 m/s and the evaporation from 3.15 to 2.75 mm/day, every field finite.

**The cases** (`box_near_side_sea_l15_wind_day3`, `_wind_day4`, `_wind_gust_day3` and `_wind_gust_day4`). The third
days' inputs equal box_near_side_sea_l15_day3's byte for byte, the namelist differing in `var10` alone, and the gust
cases differ from the wind cases in their build alone; each fourth day runs on from its own third day's restart. The
first case records runner `1c3858b609e083ea`, set up before the gust's cases were added, and the others
`c2184716516f7265`.

**The rules, set before the runs.**
- **Method.** Over the fourth lunar day the box's domain-mean wind at 0.3 and 1 km averages at least 0.85 of its
  target's speed in every quarter; for the gust box, its surface layer's gust also stays at most 1.0 m/s in every
  quarter. Otherwise the box is not read as a test of the change.
- **Evaporation (the wind box).** G = (E_wind − E_old)/(E_GCM − E_old), with E the lunar day's mean surface vapour flux,
  E_old 3.03 mm/day (the fourth day as first run) and E_GCM 3.75. At least 0.5: the box's short wind accounts for most
  of the gap; at most 0.25: the gap lies mostly in the two models' surface exchange.
- **The late evening (both).** Rain over the fourth day's 135–180° of at least 0.93 mm/day, a third of the GCM's 2.8:
  CM1 supports part of the GCM's late-evening rain; at most 0.4: the lull holds.
- **The night (both).** Rain after midnight over the fourth day of at least 1.0 mm/day: the box still rains at night as
  the GCM does; at most 0.5: the quiet night returns (the thresholds of 27 October).

**The runs.** The two boxes ran side by side on three threads each, two lunar days in 2.6 hours (`run_wind.sh`). Over
the fourth lunar day, by quarter from sunrise, beside the box as first run and the GCM at the site over four years
with Earth's and with the lunar free-convection constant (`wind_read.py` beside the runs;
[box_days_box_near_side_sea_l15_wind_gust_day3.json](../results/crm/box_days_box_near_side_sea_l15_wind_gust_day3.json)
and the wind box's, with their onset and budget products):

| Fourth lunar day | Box as first run | Wind following the GCM's | Wind and lunar gust | GCM, Earth constant | GCM, lunar constant |
|---|---|---|---|---|---|
| Mean wind at 0.3 km against its target | 0.42, 0.53, 0.66, 0.55 | 0.89, 0.95, 1.03, 0.85 | 0.89, 0.96, 1.03, 0.86 | | |
| Gust over the sea, m/s | 1.9, 1.7, 2.0, 2.1 | 1.7, 1.5, 1.9, 2.0 | 0.63, 0.48, 0.67, 0.67 | | |
| Evaporation; daily mean, mm/day | 2.79, 3.30, 3.32, 2.72; 3.03 | 3.41, 4.01, 3.36, 2.89; 3.42 | 3.28, 3.92, 3.04, 2.61; 3.21 | 3.61, 4.54, 3.50, 3.37; 3.75 | 3.34, 4.67, 3.20, 3.12; 3.58 |
| Air at 0.3 km, K | 298.18, 299.74, 298.02, 297.25 | 298.80, 300.07, 298.39, 297.58 | 298.54, 299.90, 298.12, 297.24 | 298.83, 300.24, 299.66, 298.72 | 298.90, 300.35, 299.67, 298.69 |
| Vapour at 0.3 km, g/kg | 11.08, 8.90, 10.11, 12.06 | 11.18, 9.18, 10.85, 12.36 | 10.82, 8.81, 10.42, 11.93 | 12.36, 12.11, 11.93, 12.29 | 11.99, 11.62, 11.77, 12.08 |
| CM1's inhibition, J/kg | 11, 30, 19, 1.0 | 5.7, 24, 8.9, 1.1 | 6.1, 25, 12, 1.5 | | |
| Rain; daily mean, mm/day | 0.00, 1.01, 1.06, 1.91; 1.00 | 0.23, 1.06, 1.26, 2.32; 1.22 | 0.15, 0.97, 1.16, 2.30; 1.14 | 0.15, 1.71, 3.02, 2.18; 1.78 | 0.20, 1.18, 2.40, 1.72; 1.38 |
| Rain at 60–75°, 75–90°, 90–105°, 105–120°, 120–135° | 0.83, 1.93, 2.68, 2.45, 0.72 | 0.86, 1.67, 2.96, 2.57, 0.48 | 0.52, 1.48, 2.64, 3.44, 0.46 | 2.47, 4.90, 4.67, 2.73, 1.89 | 1.38, 4.13, 2.28, 1.84, 2.03 |
| Rain at 135–150°, 150–165°, 165–180°; the late evening | 0.21, 0.09, 0.13; 0.14 | 0.17, 0.09, 1.28; 0.49 | 0.17, 0.09, 0.06; 0.11 | 2.82, 2.29, 3.64; 2.93 | 2.94, 2.86, 2.47; 2.76 |

The vapour is a mixing ratio; the GCM's columns are its 3-day means at the six cells laid on the box's heights.

**By the rules.**
- **The box with the wind and the lunar gust passes both method checks** (its wind 0.86–1.03 of its target, its gust
  0.48–0.67 m/s), so it is read. Its late evening rains 0.11 mm/day: the lull holds. Its night rains 2.30 mm/day after
  midnight: it still rains at night as the GCM does (2.28 with Earth's constant, 1.72 with the lunar one).
- **The box with the wind alone misses its method check narrowly:** its wind after midnight averages 0.847 of its
  target's at 0.3 km (0.88 at 1 km), under the 0.85 its rule asked of every quarter, so it is not read as a test of the
  wind. Its measures stand as observations: an evaporation of 3.42 mm/day, half the gap to the GCM (G 0.54); a late
  evening of 0.49 mm/day, its storms returning at 165–180°; a night of 2.32. Over its third lunar day, which passed the
  same check, G was 0.33 and the late evening 0.44 mm/day.

**What the two show** (an interpretation of the runs above).
- **The wind made part of the evaporation gap, and the gust gave part of it back.** With the GCM's wind the box
  evaporates 3.42 mm/day against 3.03; with the lunar gust as well, 3.21. The Earth gust had stood in for about half of
  what the short wind withheld.
- **The late evening stays quiet with lunar surface physics.** The GCM's wind with the Earth gust brings the box's
  storms back before midnight; with the lunar gust the lull holds through 180°, as in the box as first run. The GCM's
  rain over those hours, 2.8–2.9 mm/day with either constant, stays unsupported.
- **The night holds in every box.** After midnight the three boxes rain 1.9–2.3 mm/day from the same regime: a
  boundary layer 4.2–4.5 km deep, 1,670–1,830 J/kg of CAPE released through 1–1.5 J/kg of inhibition, a
  free-convection level at 3.9–4.2 km and the 99.9th percentile of w at 16–60 km at 2.7–3.8 m/s. With the GCM's wind
  it rains 2.3 mm/day, as the GCM with Earth's constant does.
- **The low air's difference now lies above the surface.** The corrected box evaporates within 2–16% of the GCM with the
  lunar constant in every quarter, yet its air at 0.3 km stays 1.2–2.8 g/kg drier by day and 1.5 K cooler in the evening
  and the night. The difference that feeds the GCM's late-evening rain lies less in the sea's supply than in how each
  model carries heat and vapour up from its lowest layers, where CM1 mixes them through a boundary layer far deeper than
  PlaSim's, as
  [the boundary-layer comparison](../gcm/README.md#plasims-boundary-layer-mixing-at-lunar-gravity-27-october) found.
  The test of that reading the same day did not bear it out: the GCM carries vapour out of its lowest kilometre as CM1
  does, and the gap lies in how much vapour each column holds
  ([the next section](#the-sea-box-on-the-corrected-build-its-wind-pulled-from-the-gcms-lowest-level-up-28-october)).

**What it changes.** Every CM1 run over water so far, the rings and boxes among them, took MM5's gust over the sea, and
every box forced by local time took the 3-day relaxation of its wind; both bear on their evaporation, in opposite
directions. The box with both corrections is the nearest CM1 stands to the GCM's own forcing with lunar surface physics:
on it the night's agreement with the GCM holds and the late evening's gap remains. The two corrections became the
defaults for the boxes and rings that follow, the wind's pull starting at the GCM's lowest level (register, 28
October; [the next section](#the-sea-box-on-the-corrected-build-its-wind-pulled-from-the-gcms-lowest-level-up-28-october));
PlaSim took Beljaars's gust for a four-year test
([the GCM's account](../gcm/README.md#beljaarss-gust-in-plasims-surface-transfer-over-four-years-28-october)); and
the vertical carrying of the low air's heat and vapour, tested next, put the gap in the columns' vapour instead.

**Evidence and limits.** One realization of each box over two lunar days from the day-59 restart, whose first day
carries the box's earlier state, forced by four years of a GCM still settling. The boxes hold the sea at 301.02 K, so
unlike the GCM's mixed layer it cannot warm as its evaporation changes. The gust takes CM1's land setting, vconvc 1,
where Beljaars suggested 1.2 and COARE takes 1.25, a gust 0.1–0.15 m/s larger. The quarters' measures are means of
the domain diagnostics every 3 hours; the gust is inferred from the domain means of the surface layer's wind and the
resolved wind at 50 m. The GCM's columns are 3-day means.

**Cost.** The builds took 4 and 2 minutes; the four trials 18–35 s each. The two boxes ran side by side in 2.6 hours on
six threads and hold 22.4 GB with every restart; the trials' copies hold 0.2 GB.

```sh
W=/media/projectspace/terluna-research/session-tools/wind_fc
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_rain_law_pbl_wind moon_omp_rain_law_pbl_wind_gust --jobs 4
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_near_side_sea_l15_wind_day3   # and _wind_gust_day3
climate/gcm/.venv/bin/python $W/trial.py run && climate/gcm/.venv/bin/python $W/trial.py check   # bit check and trials
$W/run_wind.sh wind & $W/run_wind.sh wind_gust &                               # side by side, 2.6 hours
climate/gcm/.venv/bin/python $W/wind_read.py                                   # the rules' readings, both boxes
climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea_l15_wind_gust_day3 box_near_side_sea_l15_wind_gust_day4
climate/gcm/.venv/bin/python -m climate.crm.box_onset box_near_side_sea_l15_wind_gust_day4 88.5 118   # and _wind_day4
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_l15_wind_gust_day4 --from-day 88.5
```

### The sea box on the corrected build: its wind pulled from the GCM's lowest level up (28 October)

The wind's own time scale and Beljaars's gust over water became the defaults for the boxes and rings that follow
(register, 28 October), and the wind took one more change first. Below the GCM's lowest level, 313 m over the sea site,
the boxes' wind targets are that level's wind held to the ground (`gcm_site.wind_windows`), and a ring's targets are
held the same way. On 3 days the pull hardly reached the box's surface layer. On 2 hours it set the lowest levels
against their own drag: in the gust box's fourth day the 50 m wind stood at 0.92–0.98 of the 0.3 km wind (domain-mean
vector winds, by quarter), against 0.90–0.92 in the box as first run, although a stronger wind should leave the air
near the ground further behind, the air less unstable and the sea rougher.

**The build** (`WIND_ALOFT_PATCH`, the build `moon_omp_rain_law_pbl_wind_gust_aloft`, `CORRECTED_BUILD`). With `var8`
above 0 it pulls the wind only at levels from `var8` m up and leaves the levels below to the box's own surface layer and
turbulence; with `var8` at CM1's default of 0 it runs as the gust build. A case on it takes the corrected wind unless its
namelist sets otherwise (`corrected_wind`): `var10` 2 hours and `var8` the GCM's lowest level from the case's reference
profile, in a fresh set-up and in a run on from a restart alike. Over the sea site that level stands at 312.7 m, so
CM1's levels at 50, 150 and 250 m are the box's own.

**The trials**, each read by a rule set before it ran (`rules.txt` part A beside the runs; `trial_aloft.py`). From the
day-59 restart on three threads:
- With `var8` 0 the build repeats the gust build bit for bit in all 16 files a three-hour trial writes.
- With `var8` 312.7 m the 50 m wind falls to 0.92–0.93 of the 0.3 km wind, against 1.00–1.02 with the pull to the
  ground.
- The three-hour rule also asked that the wind at 0.3 and 1 km come within 10% of the gust build's. At 0.3 km it fell
  12% short, so that rule failed as written: the free levels are mixed with the pulled ones, the surface's drag reaches
  both, and three hours from a 1.4 m/s start are still far from steady. A 24-hour pair, read by a rule set before it
  ran, narrowed the gap to 8% at 0.3 km and 2% at 1 km and held there from day 59.6, the 50 m wind staying below the
  full pull's: it passed. Over that first day both builds level off below their targets, the full pull at 0.82–0.84 of
  its target at 0.3 km and the pull from 313 m at 0.74–0.76, and the evaporation runs 9% lower with the levels below
  left free.

**The box days** (`box_near_side_sea_l15_aloft_day3` and `_aloft_day4`, from the day-59 restart; their inputs the gust
box's byte for byte but for `var8`). Day 3 ran on two threads beside the GCM, in 1.6 hours, and day 4 on three, in 1.2
hours. They were read on the fourth lunar day by the rules of `rules.txt` part B, set before they ran. Those rules
check the wind at the pulled levels, 0.5 and 1 km, with 0.3 km now in the layer the box shapes itself. The table compares
the fourth day with the gust box's and with the GCM that takes the same gust over four years
([the GCM's account](../gcm/README.md#beljaarss-gust-in-plasims-surface-transfer-over-four-years-28-october)); the
products are
[box_days_box_near_side_sea_l15_aloft_day3.json](../results/crm/box_days_box_near_side_sea_l15_aloft_day3.json) and
the fourth day's onset and budget.

| Fourth lunar day, by quarter from sunrise | Wind pulled to the ground, lunar gust | Wind pulled from 313 m, lunar gust | GCM with Beljaars's gust |
|---|---|---|---|
| Mean wind at 0.5 km; at 1 km, against its target | 0.89, 0.96, 1.03, 0.87; 0.89, 0.92, 1.02, 0.89 | 0.87, 0.93, 1.03, 0.84; 0.88, 0.91, 1.02, 0.88 | |
| Mean wind at 0.3 km against its target | 0.89, 0.96, 1.03, 0.86 | 0.85, 0.92, 1.03, 0.81 | |
| The 50 m wind against the 0.3 km wind | 0.97, 0.94, 0.94, 0.98 | 0.92, 0.90, 0.93, 0.92 | |
| Gust over the sea, m/s | 0.63, 0.48, 0.67, 0.67 | 0.62, 0.49, 0.68, 0.66 | |
| Evaporation; daily mean, mm/day | 3.28, 3.92, 3.04, 2.61; 3.21 | 3.13, 3.79, 3.09, 2.43; 3.11 | 3.23, 4.87, 2.88, 2.71; 3.41 |
| Air at 0.3 km, K | 298.54, 299.90, 298.12, 297.24 | 298.40, 299.84, 298.15, 297.25 | 299.05, 300.59, 299.78, 298.85 |
| Vapour at 0.3 km, g/kg | 10.82, 8.81, 10.42, 11.93 | 10.68, 8.71, 10.25, 11.88 | 11.87, 11.59, 11.62, 12.02 |
| Rain; daily mean, mm/day | 0.15, 0.97, 1.16, 2.30; 1.14 | 0.03, 1.17, 1.03, 1.89; 1.03 | 0.14, 1.46, 2.81, 1.15; 1.40 |
| Rain at 60–75°, 75–90°, 90–105°, 105–120°, 120–135° | 0.52, 1.48, 2.64, 3.44, 0.46 | 0.60, 1.37, 3.07, 2.12, 0.47 | 2.27, 4.20, 4.02, 2.73, 2.20 |
| Rain at 135–150°, 150–165°, 165–180°; the late evening | 0.17, 0.09, 0.06; 0.11 | 0.30, 0.10, 0.06; 0.15 | 2.04, 3.70, 2.20; 2.62 |

**By the rules.** The fourth day misses its wind check narrowly. After midnight its wind at 0.5 km averages 0.840 of its
target, under the 0.85 asked of every quarter (1 km holds at 0.876). Its gust (0.49–0.68 m/s) and its near-ground
profile (the 50 m wind at 0.90–0.93 of the 0.3 km wind) pass. Like the wind box of the morning, it is not read as a
test, and its measures stand as observations. The late evening rains 0.15 mm/day and the night 1.89 mm/day after
midnight. The third day, which passed the same checks, rains 0.13 and 1.91. With the pull left off the lowest levels
the evaporation falls 3% from the gust box's and the night's rain 18%, while the late evening stays quiet. Against the
GCM on the same surface physics, the box rains 1.03 mm/day over the lunar day against 1.40: half the GCM's from noon
to midnight (1.10 against 2.13 mm/day), and 0.15 against 2.62 in the late evening. After midnight it rains 1.6 times
the GCM's 1.15 mm/day, at evaporation within 11%.

**How the two models carry vapour up from the sea** (`rules.txt` part D, set before either run's budgets were read;
`carry.py`, `carry_direct.py`). This morning's account read the low air's gap as a sign that PlaSim carries heat and
vapour up from its lowest layers more slowly than CM1. The test set that against the vapour each model's turbulence and
convection carry up through 0.6 and 1.2 km, as a share of the evaporation, and the fall of vapour from 0.3 to 1.2 km. In
the box the flux came first from CM1's vapour budget by process. That budget returns the evaporation over the whole
column within 1.6% in every quarter, as the rule asked, but not over the lowest 1.2 km, where its terms sum to −1.3 mm/day
against a change of −0.03. A direct measure, the resolved eddy flux ρ⟨w′q′⟩ from the 3-D output every 3 hours plus the
boundary-layer scheme's flux, stands beside it.

| Fourth day / four years, by quarter | CM1, direct | CM1, from its budget | GCM with the gust |
|---|---|---|---|
| Share of the evaporation carried up through 0.6 km | 0.91, 0.76, 0.88, 0.79 | 1.27, 1.19, 1.34, 1.18 | 0.90, 0.94, 0.95, 0.88 |
| The same through 1.2 km | 0.58, 0.43, 0.66, 0.43 | 0.97, 0.75, 1.16, 0.76 | 0.59, 0.80, 0.74, 0.57 |
| Vapour's fall from 0.3 to 1.2 km, g/kg | 0.49, 0.56, 0.24, 0.16 | | 0.36, 0.49, 0.18, 0.22 |
| Vapour at 0.3, 1.2, 2.5, 4 and 6 km by day (sunrise to noon), g/kg | 10.7, 10.2, 9.9, 9.8, 9.2 | | 11.9, 11.5, 10.6, 10.2, 7.5 |
| The same, noon to sunset | 8.7, 8.2, 8.1, 8.2, 8.4 | | 11.6, 11.1, 9.8, 9.4, 8.7 |
| Potential temperature at 0.3 and 4 km, sunrise to noon; sunset to midnight, K | 282.9, 283.3; 282.6, 282.7 | | 283.4, 284.7; 284.1, 284.3 |

By the rule the slow-carrying reading holds nowhere: on either measure of CM1's flux the GCM carries a like share
through the lowest kilometre with a vapour fall no steeper than CM1's (0.74–0.87 of it by day). From noon to sunset
the falls agree within 13% while the air at 0.3 km differs by 2.9 g/kg, so there the difference lies in the column's
whole vapour; in the morning and evening the falls sit 26–28% apart, the GCM's the smaller, with gaps of 1.2–1.4 g/kg;
after midnight the low air agrees within 0.14 g/kg. The profiles say the same: by day the GCM is 1.2–2.9 g/kg moister
from 0.3 to 1.2 km alike, less so higher, and drier at 6 km in the morning, while CM1's vapour stands nearly uniform up
to 4 km. The low air's gap is not in how either model carries vapour out of its lowest layers; it lies in how much
vapour each column holds.

**A lead on the column's vapour.** The box's day-night forcing takes the GCM's circulation by local time, less the
part the GCM's mean vertical wind carries, which CM1 is to apply itself on the box's own profile. Over the lunar day the
GCM's circulation exports 1.99 mm/day of vapour from the site's column (held to its exact budget), and the mean wind's
part the table leaves out comes to about +0.48 mm/day, so the table removes about 2.5. In the box the table removes
2.45 mm/day (box_budget) and CM1's imposed ascent gives −0.10, not +0.48. The box therefore exports about half a
millimetre a day more vapour than the GCM it stands for, about a quarter of the GCM's export, which would dry its whole
column. Why CM1's large-scale vertical advection gives so little back is checked in the next section.

**What it changes.** The corrected build stands as the default, its near-ground wind the box's own. On the same surface
physics the two models rain alike in amount over the lunar day and differently in time: the GCM in the evening, CM1
after midnight. The low air's gap lies in the column's vapour, with the box's forcing a named suspect, rather than in
the boundary layer's carrying.

**Evidence and limits.** One realization of each box day from one restart, forced by the baseline GCM's four years,
whose light-wind term differs from the gust GCM's; the GCM's numbers are its four years at the six cells. CM1's
direct flux takes resolved snapshots every 3 hours and the boundary-layer scheme's tendency, and neglects density
fluctuations; why CM1's budget does not close over its lowest kilometre is not traced. The GCM's flux comes from its
single-step budgets at the end of each 3-day mean, and below 1.2 km its convection condenses little. The heat budget
of the box does not close over the column (its resolved terms leave several W/m² against a surface flux of a few), so
the heat flux is not read; the potential temperatures are.

```sh
W=/media/projectspace/terluna-research/session-tools/gust_aloft
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_rain_law_pbl_wind_gust_aloft --jobs 4
python3 $W/trial_aloft.py run && python3 $W/trial_aloft.py check                # bit check and the three-hour trial
python3 $W/trial_aloft.py run gust24 aloft24 && python3 $W/trial_aloft.py check24   # the 24-hour pair
$W/run_aloft.sh                                                                  # day 3 on two threads, day 4 on three
climate/gcm/.venv/bin/python $W/aloft_read.py                                    # part B
climate/gcm/.venv/bin/python $W/carry.py && climate/gcm/.venv/bin/python $W/carry_direct.py   # part D
climate/gcm/.venv/bin/python $W/box_low_budget.py && climate/gcm/.venv/bin/python $W/gcm_low_budget.py
climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea_l15_aloft_day3 box_near_side_sea_l15_aloft_day4
climate/gcm/.venv/bin/python -m climate.crm.box_onset box_near_side_sea_l15_aloft_day4 88.5 118
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_l15_aloft_day4 --from-day 88.5
```

### The box's forcing export checked (28 October, evening)

The lead above asked why CM1's large-scale vertical advection gives back so little. The site forcing builds its table
from the GCM's circulation by local time and adds back the GCM's mean vertical wind's advection of the box's starting
air, the GCM's mean air (`gcm_site.to_cm1`); CM1 then applies the same mean wind to its own air. The two cancel while
the box's air is its starting air. Afterwards their sum is a departure term, the mean wind acting on the box's departure
from its start. The check was read by rules set before it was computed (`rules.txt` beside `forcing_check.py`):

- **The table** rebuilt from the GCM run's four years matches the box's `terluna_lsadv.txt` to two parts in a million.
- **The excess.** The box's forcing removes 0.45 mm/day (third day) and 0.55 (fourth) more vapour from the column than
  the GCM's circulation, which removes 1.99. The departure term reproduces both to 0.01 mm/day, and CM1's own imposed
  ascent matches the term recomputed from the box's vapour to 0.01.
- **Where.** The GCM's mean wind over the site rises through the lowest 5 km, fastest near 2.2 km (0.39 cm/s), and sinks
  slowly from 6 to 14 km: it draws air into the column below about 2.5 km and sends it out between 2.5 and 6 km. In
  CM1's advective form the air it draws in carries the box's own vapour. On the fourth day the box's air from 0.25 to
  1 km stands 1.5–1.7 g/kg below its start; in flux form the inflow levels take 0.83 mm/day and the outflow levels,
  drier too, give back 0.25. The layer below 1.77 km holds 0.62 of it.
- **When.** Through the first lunar day the term swung with the box's adjustment (+0.16, +0.42, +0.14 and −0.97 mm/day
  by quarter), and with the term still small the box's air at 0.2–0.7 km fell 0.8 g/kg below its start within the
  first week. From the second day it took 0.44, 0.47 and 0.57 mm/day, 43 kg/m² over three lunar days, while the
  column's vapour fell from 259 to 217 kg/m². At the fourth day's evaporation and rain, the column would have gained
  0.21 mm/day without it, against the 0.36 it lost.
- **Heat.** The same term warms the box by 12.0 W/m² on the fourth day, 12.3 of it above 23 km, where the box's air has
  fallen 7–10 K below its start under the GCM's weak mean ascent; there it restores.

**The reading.** The forcing's mean-wind part is what Randall and Cripe (1999, J. Geophys. Res. 104, 24,527) call
horizontal advective forcing, the horizontal advection prescribed and the vertical reckoned on the model's own sounding;
its hour-by-hour part is their revealed forcing, the whole advection prescribed. Over the days to weeks of field
experiments they found the two comparable, preferred revealed forcing for its simplicity, and noted that a forced
column whose sounding errors carry it away cannot find its way back. Over four lunar days the advective form here
turned the low-level deficit the box opened in its first week into a steady drying of the whole column, because the
air the GCM's mean wind draws in is the box's own and drier than the GCM's. Coupled columns under the
weak-temperature-gradient approximation carry the air entering a column at the upwind value, the neighbouring
column's or the reference's, and the air leaving at the column's own (Daleu, Woolnough and Plant 2012, J. Atmos. Sci.
69, 3683, their equation 9). In that form the box's inflow would carry the GCM's air and its drier outflow would export
less, so the departure would restore the column rather than drain it.

**What it changes.** The afternoon's reading stands where it placed the gap: in the column's vapour, not in the
boundary layer's carrying. That vapour, though, has drifted with the forcing's form since the first lunar day, and on
the fourth day the box's low air stands 1.6 g/kg below the GCM's mean air, about the low air's whole gap over the
lunar day; the comparisons made on this forcing carry that drift. How much of the gap CM1's physics would keep with its
inflow at the GCM's humidity needs a box run in that form. Boxes forced from here on should take a form that does not
build on the box's departure: the inflow at the GCM's air, recommended, or the mean wind's part prescribed with the
rest. The choice is the author's.

**Evidence and limits.** Domain means every 3 hours from one realization of each box day. The departure term takes
centred differences on the box's levels, as the table's added-back term does; CM1's own term (third-order upstream)
matches it to 0.01 mm/day. The flux-form split interpolates the density to the w levels. Nothing has been run in
another form.

```sh
F=/media/projectspace/terluna-research/session-tools/forcing_check
PYTHONPATH=. climate/gcm/.venv/bin/python $F/forcing_check.py       # rules F0-F5; forcing_check.json and .txt
PYTHONPATH=. climate/gcm/.venv/bin/python $F/low_air_by_day.py      # the low air against its start by lunar day
```

### The sea box with the mean wind's inflow at the GCM's air (28 October, evening)

Boxes forced from here on draw the GCM's air in where its mean wind converges (register, 28 October), and one lunar
day on the fifteen-level GCM's forcing sizes how far the box's vapour drifted with its own air drawn in.

**The build** (`INFLOW_PATCH`, the build `moon_omp_rain_law_pbl_wind_gust_aloft_inflow`). CM1's large-scale vertical
advection (`wsub`) is the flux form less the scalar times the mean wind's mass convergence, div, so the air the
convergence draws in carries the box's own mean. With `var7` = 1 the build carries the reference's air in where div
is positive and the box's own out where it is negative, the upwind form of coupled columns (Daleu, Woolnough and Plant
2012, their equation 9). Each level's domain-mean potential temperature and vapour are drawn toward the nudging
reference, the GCM's mean air and the box's starting sounding, at the rate max(div, 0)/ρ₀ from CM1's own base state.
The term vanishes while the box's air is the reference, so the table's added-back advection still cancels at the
start. The rates go to `terluna_inflow_rate.txt`, from which `box_budget` reckons the term at each output, as it does
the table; with `var7` at 0 the build runs as the aloft build. It is the default build from here on
(`CORRECTED_BUILD`): a box forced by the GCM's column at its site takes `var7` 1 from `corrected_inflow` unless its
namelist sets otherwise, and a ring on it keeps `var7` 0 and runs as the aloft build. The aloft days keep the build they
ran on.

**The trials**, read by rules set before the build finished (`rules.txt` part G beside `trial_inflow.py`), three model
hours from the day-59 restart on three threads:
- With `var7` 0 the build repeats the aloft build's trial bit for bit in all 16 files the trial writes.
- With `var7` 1 the rates CM1 writes match the rates reckoned from the GCM's vertical wind and the box's density to
  0.983–1.020 at the 48 levels compared, and stand at zero where the flow diverges. The largest, 2.6 × 10⁻⁶ per second
  at 450 m, renews the air there in about 4½ days; the inflow reaches 2.5 km.

**The day** (`box_near_side_sea_l15_inflow_day3`, from the day-59 restart, its inputs the corrected third day's byte
for byte but for `var7`; 0.9 hours on four threads), read by the rules of `rules.txt` part H, set before the build
finished; the products are
[box_days_box_near_side_sea_l15_inflow_day3.json](../results/crm/box_days_box_near_side_sea_l15_inflow_day3.json) and
[box_budget_box_near_side_sea_l15_inflow_day3.json](../results/crm/box_budget_box_near_side_sea_l15_inflow_day3.json).

| Third lunar day, by quarter from sunrise | The box's own air drawn in (corrected build) | The GCM's air drawn in |
|---|---|---|
| Mean wind at 0.5 km against its target | 0.888, 0.953, 1.034, 0.857 | 0.876, 0.954, 1.033, 0.849 |
| Evaporation; daily mean, mm/day | 2.82, 2.99, 2.91, 2.26; 2.74 | 2.72, 2.76, 2.62, 2.10; 2.55 |
| Air at 0.3 km, K | 298.56, 299.70, 298.27, 297.48 | 298.78, 299.70, 298.71, 298.28 |
| Vapour at 0.3 km, g/kg | 11.52, 10.49, 10.67, 12.42 | 11.61, 11.18, 11.28, 12.52 |
| Rain; daily mean, mm/day | 0.00, 1.28, 1.02, 1.91; 1.05 | 0.00, 1.00, 0.82, 0.61; 0.61 |
| Rain in the late evening (135–180°), mm/day | 0.13 | 0.09 |
| The column's vapour from the day's start to its end, kg/m² | 247.4 to 227.8 | 247.4 to 238.5 |
| Vapour against the start at 0.2–0.7, 2–5 and 5–10 km, g/kg | −0.75, +0.05, +0.90 | −0.38, +0.34, +1.22 |
| Heat stored by the column; brought by the inflow, W/m² | −20.4; none | −10.6; +22.5 |

**By the rules.** The vapour budget closes with the inflow reckoned (its rest 0.05 mm/day against the corrected day's
0.075), and the inflow does what it was built to do: it gives back 0.147 mm/day where CM1's advection at the inflow
levels takes 0.148, and what remains of the forcing's excess over the GCM's circulation, 0.31 mm/day, is the outflow
levels' part (0.33). The column's drying falls from 0.67 to 0.30 mm/day, 0.37 less where the rule asked 0.3. The day
misses its wind check by 0.001: after midnight its wind at 0.5 km averages 0.849 of its target, where the corrected
day's passed at 0.857. By the rule its measures, the drying's with them, stand as observations.

**What the day shows**, as observations:
- **The drift halves.** The column ends 10.7 kg/m² moister than on the corrected day, its air at 0.2–0.7 km 0.38 g/kg
  below the GCM's mean air rather than 0.75, and 0.6–0.7 g/kg moister at 0.3 km from noon to midnight. What drying
  remains comes from the outflow: the box's air at 2–10 km stands moister than the GCM's (+0.34 and +1.22 g/kg), and the
  air leaving there carries more vapour out, a term that pulls the column back toward the GCM's.
- **The rain falls with it.** The day rains 0.61 mm/day against 1.05, after midnight 0.61 against 1.91, and the late
  evening stays quiet (0.09). The water budget says why: the rain is the evaporation and the forcing less the column's
  change, and on the corrected day the column's drawdown, 0.67 mm/day, fed it. With the GCM's air drawn in the drawdown
  is 0.30, the forcing takes 0.15 mm/day less, and the evaporation falls 0.19 under the moister low air.
- **The evaporation now leads.** Once the column stops drying, the box rains its evaporation, 2.55 mm/day, less the
  forcing's export, which lies between the GCM's 1.99 and today's 2.29: 0.3–0.6 mm/day, where the GCM rains 1.78 from
  3.75 mm/day of evaporation (the gust GCM evaporates 3.41). The two models' rain at the site, and its timing, rest
  first on that gap.
- **Heat.** The inflow brings 22.5 W/m², 14.5 of it at 16–60 km, where the box's air stands 5.5–8.4 K below the GCM's
  mean air. The box still cools, by 10.6 W/m² against 20.4, and its condensation gives 12 W/m² less as it rains less.

**What it changes.** Boxes forced by the GCM's column take the inflow from here on. The old form's drawdown fed part of
the rain the earlier comparisons read on this forcing, the night's 1.9 mm/day after midnight among it; with the GCM's
air drawn in, the night rains 0.6 mm/day. Where the two models part at the site now turns first on the box's
evaporation, a third below the GCM's, then on its remaining cooling.

**Evidence and limits.** One realization of one lunar day from one restart, forced by the fifteen-level GCM's four
years without the gust, read as observations. The inflow draws in the GCM's mean air, not its air hour by hour, at
rates from CM1's base state. The rain once the column stops drying is arithmetic on this day's terms, not a run to
equilibrium.

```sh
I=/media/projectspace/terluna-research/session-tools/inflow
climate/gcm/.venv/bin/python -m climate.crm.cm1_run build moon_omp_rain_law_pbl_wind_gust_aloft_inflow --jobs 4
PYTHONPATH=. climate/gcm/.venv/bin/python $I/trial_inflow.py run          # part G: the two three-hour trials
PYTHONPATH=. climate/gcm/.venv/bin/python $I/trial_inflow.py check
$I/run_inflow.sh                                                          # the day on four threads; rerun to resume
PYTHONPATH=. climate/gcm/.venv/bin/python $I/inflow_read.py               # part H
climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea_l15_inflow_day3
climate/gcm/.venv/bin/python -m climate.crm.box_budget box_near_side_sea_l15_inflow_day3 \
    --reference box_near_side_sea_l15_aloft_day3 --from-day 59
```

### Where the two models' evaporation parts over the sea (28 October, night)

With the GCM's air drawn in, the box's rain at the site turns first on its evaporation, a third below the forcing GCM's
(2.55 against 3.75 mm/day). The trace asks whether that gap lies in the air each model holds or in how each turns air
into evaporation, by rules set before it was computed (`rules.txt` beside `evap_trace.py`). PlaSim's own transfer
(`surface_gust.surface_transfer`, its `surflx` step for step, with its Charnock roughness iterated to its own stress and
its implicit evaporation step) is applied to the GCM's own 3-day means at the site's sea cells and to the box's
domain-mean air at the GCM's lowest level, 312.7 m, averaged over the same tenths of the lunar day; CM1's own transfer
at its 50 m level is read from the box's diagnostics.

| Over the Near-side Sea | PlaSim on the GCM's air | PlaSim on the box's air, Earth's term | PlaSim on the box's air, the gust | CM1's own, at 50 m |
|---|---|---|---|---|
| Evaporation, mm/day (inflow day; fourth day) | 3.43 (it wrote 3.75) | 3.28; 3.94 | 1.72; 2.04 | 2.55; 3.11 |
| Transfer ρCU, 10⁻³ kg/m²/s | 5.44 | 4.91 | 2.63 | 3.85 |
| Its stability factor | 1.96 | 2.35 | 1.11 | about 1.23 |
| Wind in the formula, m/s | 3.17 | 2.50 | 2.55 | 2.45 |
| Humidity difference, g/kg | 7.48 | 8.03 | 8.03 | 7.43 |
| Roughness, m | 1.0 × 10⁻⁴ | 5.7 × 10⁻⁵ | 5.8 × 10⁻⁵ | 1.1 × 10⁻⁴ |

The columns' last three are the inflow day's; CM1's stability factor is its transfer coefficient, 1.13 × 10⁻³, over
its neutral value at 50 m from its own roughness lengths.

**By the rules.** PlaSim's transfer reproduces the GCM's own evaporation from its 3-day means to 8.5% (with the gust,
the gust GCM's to 14%), within the 15% the method check allowed; the means smooth the wind, which accounts for the
shortfall. The gap lies in the surface transfer: on the inflow day the air's part is 0.15 of the 0.88 mm/day gap
between PlaSim on the GCM's air and CM1 (17%) and the transfer's 0.73 (83%); on the fourth day the box's drier air
would draw more under PlaSim's transfer, and the transfer's part exceeds the whole gap.

**What the trace shows.**
- **Earth's free-convection term carried the forcing GCM's extra evaporation.** On the same air PlaSim's transfer with
  it evaporates 29% more than CM1's surface layer. Its interpolation over unstable water roughly doubles PlaSim's
  stability factor, where CM1's similarity raises its coefficient about 1.23 times neutral.
- **With Beljaars's gust, as the design run will take it, PlaSim's transfer falls a third below CM1's** on the same air.
  Its stability factor stands at 1.11 at the 313 m level and its roughness at half CM1's, from a smaller u*. The box's
  wind and vapour change less between 50 and 313 m (about 1.08–1.09 times) than the log profile PlaSim's coefficient
  assumes from its roughness (about 1.13–1.14), which suggests the level sits near the top of the surface layer or
  above it, an inference the trace does not test.
- **CM1's transfer sits between PlaSim's two forms**, and the GCM makes up part of any shortfall with its own air:
  with the gust its lowest level dries by 1.9 g/kg, so its evaporation falls 9% where its transfer falls 32%.

**What it changes.** The gap the box's rain turns on is a property of PlaSim's surface transfer at its lowest level,
not of either model's air. The design run takes the gust, which brings the sea's transfer from 29% above CM1's to about
a third below it on the same air. A next step, proposed: set PlaSim's transfer from its 313 m level beside
Monin–Obukhov similarity integrated to the same height on the box's and the GCM's states, a column calculation of
minutes; if they part by more than a tenth, the design configuration would take similarity over water beside the gust.

**The similarity, set beside it (29 October).** CM1's own surface layer over water, carried to the GCM's lowest level,
gives CM1's written evaporation at its own 50 m within 0.3% and the same evaporation at 313 m with CM1's wind. On the
same air and the gust's wind at 313 m, PlaSim's transfer evaporates 0.72 of the similarity's on the box's inflow day and
0.82 on the gust GCM's own air, most of the difference in light wind over warm water, so the design configuration takes
similarity over water beside the gust ([Similarity over water in
PlaSim](../gcm/README.md#similarity-over-water-in-plasim-29-october)).

**Evidence and limits.** The GCM's side rests on 3-day means at six cells over four years, the box's on domain means
every 3 hours of one lunar day each; PlaSim's formula on a domain-mean wind leaves out the spread of the box's own
winds, which CM1's column-by-column fluxes keep and the gust is meant to stand for. PlaSim's roughness is iterated on
the mean state.

```sh
E=/media/projectspace/terluna-research/session-tools/evaporation
PYTHONPATH=. climate/gcm/.venv/bin/python $E/evap_trace.py       # rules EV1-EV3; evap_trace.json
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
there; at Earth's 9.81 m/s² it leaves them as they are. Until 2026-10-06 rain
and cloud ice ventilated at Earth's speeds: rain's ventilation, where it
evaporates (NUCOND) and in the gather-scatter routine, takes rain's fall-speed
coefficient from the module's constants rather than from the scaled speeds, and
cloud ice's takes a Reynolds number fitted on Earth. A review found it, and both
now scale as the fall speeds do; before, lunar rain evaporated 1.3–1.5 times too
fast for its drop sizes (a slope diameter of 0.1–0.5 mm) and the larger ice
crystals deposited up to 15 % too fast, in every NSSL run before
`box_0e_elec_corrected`'s restart from day 0. Snow, graupel and hail ventilate
from their scaled speeds, and Morrison's scheme scales its coefficients
themselves, so neither changed. WRF-ELEC's driver and
its sedimentation driver run in one thread, since WRF divides its work outside
the microphysics; the build gives each slab loop the OpenMP directive CM1's
copy has, with every per-slab array private, the charge totals summed across
threads, and what the driver sets before the loop shared.

The lunar boxes change two of WRF-ELEC's Earth settings (`cm1_elec.LUNAR`).
WRF-ELEC's 0.75-s sub-step lets graupel settle through about 1 % of the 500-m
layers of its storms in each; in the boxes' charging zone, 25–35 km up, the
layers are 2 km deep and graupel falls at 0.44 of Earth's speed, so a 6.8-s
sub-step keeps that share; the driver rounds it to whole sub-steps of the model
step, which in the boxes is one sub-step of 8 s. On the Moon air warmer than −7 °C lies about 27–34
km up, so a downward channel strikes the ground when it comes within 5 km of it,
the height that rule stands at on Earth. The 50-MV and 10-kV/m conditions for a
ground strike, the 0.9 and 0.8 fractions of breakdown, the 200-V/m internal
field and the nitrogen oxide yield stay as WRF-ELEC calibrates them for Earth.
That height rule keeps WRF-ELEC's demand that the channel's end lie in charge of
the matching sign, and the lunar storms hold no charge 5 km up, so under it no
lunar flash could strike the ground; see "Ground strikes" below. From
2026-10-05 the lunar boxes take WRF-ELEC's own rule (the author's decision), an
upper bound on their ground strikes until a rule for a leader's crossing of the
air below the storm replaces it; every run before then used the height rule.

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

A case sets the electricity in CM1's namelist: `var2` the onset of point
discharge from the ground (V/m at sea-level density; 0 for none), `var3` the
screening layers (0 none, 1 with WRF-ELEC's Earth conductivity, 2 with the
conductivity file), `var4` the height a downward channel must reach to strike
the ground (m; below zero for WRF-ELEC's rule, air warmer than −7 °C, its
default −1; at zero only the two lowest levels count), `var5` the sub-step (s;
0 for WRF-ELEC's 0.75 s), `var6` the charging (2 non-inductive, 3 with
inductive), `var7` the charging law (12 Saunders and Peck, cut off below
−32.5 °C; 11 the same law as WRF-ELEC's test case sets it, with a smooth
critical rime accretion rate and charging below −32.5 °C; 1 Takahashi),
`var8` the lightning (0 none, 1 cylinders, 3 branched; 2 and 4 the same with
the breakdown field's 180-kV/m cap lifted), `var9` leakage
through the conductivity of
[atmosphere/electricity](../../atmosphere/electricity/README.md) and `var10`
the cylinders' radius. Each run logs every flash or cylinder
call (`terluna_flashes.txt`: its kind, start, the field and breakdown field
there and the domain's largest field over breakdown before it, the area and the
levels where it neutralizes charge, the positive and negative charge it
neutralizes, the electrostatic energy before and after, its channel points and
nitrogen oxides), the area and levels being those where the flash hands
charge to the small ions (before 2026-10-05, only `lightmsz`'s channel points,
which it marks where a flash takes more than its threshold `scth` and which a
weak flash high in thin air could miss everywhere: 13 of `box_0e_elec`'s
flashes, of 5–72 C starting at 50–68 km, logged no reach),
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
realizations, so its effect on the lightning stays unresolved. Both runs left
the ions' charge in place, as nothing conducted it away, and that leftover
charge made most of the field away from the storms; with leakage it falls to a
median 0.2 kV/m, below the onset ("Leakage", below).

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
strikes that set off sprites. With leakage, which clears the charge the runs
otherwise leave in the air, the window makes 39 ground strikes among 161
flashes, a median 110 C each and at most 624 C ("Leakage", below).

The rule makes this an upper bound: it counts a channel that reaches −7 °C as
reaching the ground, which on Earth leaves 3–5 km to cross and on the Moon about
34 km. WRF-ELEC's stopping field, which halted every channel 15–25 km up, makes
the lower bound of none. Between them lies the physics of a leader crossing 34
km of clear air from a cloud at hundreds of megavolts: its channel's internal
field of a few kV/m costs about 50–125 MV over that distance, leaving most of
the potential at its tip, and Earth's leaders cross 10–15 km of clear air (bolts
from the blue) or climb 70 km (gigantic jets) on far less. "The leader's
crossing" below adds a rule for that crossing: at the literature's internal
field it removes up to a fifth of the strikes this rule counts. On several threads the window is a separate
realization of these days, but the change from none to 21 is the rule's.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_ground_rule    # and supercell_elec_ground_rule
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_ground_rule --threads 4 --hours 6
```

### The charging law

Every run above charges by Saunders and Peck's law. Stage 1's
one-dimensional estimate, from the Morrison storms' graupel, found that law
nearly stopped for slow lunar graupel, which rimes below its threshold, and
Takahashi's law charging two hundred to five hundred times faster. On
2026-10-05 both windows ran again under Takahashi's law (WRF-ELEC's `isaund`
1, with its crystal-size and impact-speed factor, on its table, whose 0 to
−30 °C and 0.01–30 g/m³ of cloud water hold the lunar charging zone), each
changing only the law: `box_0e_elec_takahashi_first` beside
`box_0e_elec_uncapped` (days 10.5–12, from before the first flash, with that
run's 5-km ground rule) and `box_0e_elec_takahashi` beside
`box_0e_elec_ground_rule` (days 40.5–42, WRF-ELEC's ground rule). The other
runs of the same days under Saunders and Peck show how far a realization moves
each number (`elec_analysis laws`, `climate/results/crm/elec_charging_laws.json`):

| Days | Run | Law | Flashes (to ground) | Charge separated, kC (non-inductive + inductive) | Hail's share of the non-inductive | In-cloud flash: median charge, energy, start | Main negative charge |
|---|---|---|---|---|---|---|---|
| 10.5–12 | `box_0e_elec` (cap) | Saunders and Peck | 13 (0) | 18.0 + 2.8 | 84 % | 80 C, 55 GJ, 35.8 km | 42–44 km, −15 to −18 °C |
| | `box_0e_elec_uncapped` | Saunders and Peck | 13 (0) | 22.0 + 3.7 | 86 % | 69 C, 55 GJ, 35.8 km | 34–44 km, −7 to −17 °C |
| | `box_0e_elec_uncapped_corona` | Saunders and Peck | 6 (0) | 13.5 + 2.4 | 82 % | 69 C, 54 GJ, 35.8 km | 38 km, −10 °C |
| | `box_0e_elec_takahashi_first` | Takahashi | 14 (0) | 25.1 + 4.8 | 37 % | 95 C, 106 GJ, 37.8 km | 48–56 km, −23 to −33 °C |
| 40.5–42 | `box_0e_elec` (5-km rule) | Saunders and Peck | 157 (0) | 122.1 + 22.4 | 88 % | 88 C, 70 GJ, 37.8 km | 38–44 km, −11 to −18 °C |
| | `box_0e_elec_ground_rule` | Saunders and Peck | 144 (21) | 116.8 + 21.2 | 88 % | 106 C, 79 GJ, 35.8 km | 38–44 km, −9 to −17 °C |
| | `box_0e_elec_takahashi` | Takahashi | 139 (10) | 111.8 + 28.2 | 43 % | 101 C, 119 GJ, 43.8 km | 50–58 km, −24 to −35 °C (8 of 12 outputs) |

The main negative charge is the level holding most negative charge at the
three-hourly outputs that hold more than 300 C of it, leaving out the decaying
storms' rain below 0 °C.

**The rate.** Takahashi's law makes as many flashes and separates as much
charge as Saunders and Peck's, within the spread of the realizations, and its
storms flash no sooner: 10.0 hours after day 10.5 against 8.2–10.1, and 4.4
hours after day 40.5 against 4.3. Stage 1's two hundred to five hundred times
does not appear because NSSL's hail does 82–88 % of the charging under
Saunders and Peck: it falls faster than graupel and rimes above the 0.1 g m⁻²
s⁻¹ below which the law gives no charge. The Morrison runs carried no hail
(their one rimed ice category ran as graupel). Under Takahashi's law graupel
separates four to six times more charge than under Saunders and Peck's and
hail two-fifths to four-fifths as much, and the totals meet.

**The sign.** At the lunar storms' cloud water, 0.03–0.17 g/m³, Takahashi's
table gives rimed ice positive charge at every temperature; it turns negative
only colder than about −25 °C at 0.2 g/m³, −15 °C at 0.5 g/m³ and −10.5 °C
at 1 g/m³. Every collision in both Takahashi windows charged graupel and hail
positively (the negative parts come to at most 1 C against 25–111 kC), and
the ice crystals and snow negatively. The rising ice carries the main
negative charge to 48–58 km (−23 to −35 °C) in 15 of the 23 three-hourly
outputs, above positive charge at 32–36 km (−3 to −9 °C) in about half of
them and lower elsewhere, and 19 of the outputs have no upper positive layer:
an inverted dipole. Under Saunders and Peck the storms hold the normal
arrangement: in 18 of the second window's 23 outputs the main negative charge
sits at 36–44 km (−8 to −18 °C), with positive charge above it at 48–72 km
and below. The laboratory measurements at the graupel's impact speeds point the same
way as Takahashi's table: at 1.8 m/s Ávila et al. (2013) found the rimer charging
positively at every temperature (−7 to −13 °C) and cloud water they used
(0.05–0.5 g/m³), and at 1.2 m/s Pradeep Kumar et al. (2024) found it positive
from −7 to −18 °C. Negative charging appears only at faster impacts: at 1.8
m/s below −17 °C (Pradeep Kumar et al.) and at 2.4 m/s below −10 °C (Ávila et
al.).

**The flashes.** Under Takahashi's law the flashes start higher (a median
37.8–43.8 km against 35.8–37.8), reach higher (channels to a median 48–50 km
against 40–42) and release one and a half to two times the energy for a similar or
somewhat larger charge; the largest in-cloud flash released 1,964 GJ against 848 GJ
under Saunders and Peck, and the second window's flashes made 484 mol of
nitrogen oxides against 350–389. The window's 10 ground strikes, by WRF-ELEC's
rule and so an upper bound, brought down a median 204 C (37–839 C, 10th to
90th percentile), the largest 1,489 C and 2.4 TJ, starting at a median 42 km
and up to 62 km. Each law here is one realization of the window; the
difference in structure lies far outside the realizations' spread, the
difference in ground strikes, 10 against 21, within reach of it.

So the charging law sets the lunar storms' polarity more than their rate.

**Impact speeds.** In the charging zone (−5 to −30 °C) of `box_0e_elec`'s
stormy windows (days 10.5–12 and 40.5–42, three-hourly output), NSSL's graupel
falls at a mass-weighted 1.4–3.9 m/s (10th to 90th percentile; median 1.9–2.4),
3.7–4.4 mm across at 410–450 kg/m³, and its hail at 4.8–7.7 m/s (median
5.9–6.6), 12–13 mm across at 710–740 kg/m³: under lunar gravity hail stays
aloft long enough to grow large. Cloud ice and snow fall at 0.25–0.49 m/s, so
hail, which does 82–88 % of the charging under Saunders and Peck, strikes them
at 4.4–7.6 m/s (median 5.4–6.5), within the 1.2–11 m/s Earth's laboratories
have measured, and graupel at 0.95–3.6 m/s, where the slow-impact studies'
1.2–1.8 m/s lie.
(`elec_analysis.py impact` applies NSSL's own laws to the output: Milbrandt and
Morrison's coefficients at each particle's density for graupel and hail, the
adjusted Ferrier law for cloud ice and Ferrier's for snow, with CM1's shape
parameters, 0 for graupel and 0.5 for hail, and the build's gravity factors.)

Most of the hail in the charging zone, 80–82 % of its mass, meets too little
cloud water to charge at all: its rime accretion rate stays under the 0.1
g m⁻² s⁻¹ below which WRF-ELEC separates no charge. The hail that charges sits
at −5 to −16 °C in 0.02–0.18 g/m³ of cloud water, and Saunders and Peck's law
charges it both ways: positively where it rimes faster than the critical rate,
which above −15 °C is Brooks et al.'s (0.1–1.5 g m⁻² s⁻¹), and negatively
where it rimes slower. The model's charging totals show the split. In the
second lunar day's window hail took +49 and −45 kC in collisions with snow
and +6 and −8 kC with cloud ice; in the first lunar day's, +8 and −5 kC and
+1.2 and −0.9 kC. Takahashi's table charges all of it positively at that
cloud water, and that turns the storms' polarity over. The measurement that
would decide between the laws is the sign of the charge rimed ice takes at
4–8 m/s in 0.02–0.2 g/m³ of cloud water between −5 and −16 °C, conditions
Earth's laboratories can reach.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_takahashi      # and box_0e_elec_takahashi_first
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_takahashi --threads 4 --hours 6
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis laws
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis impact    # elec_impact_speeds.json
```

### Leakage

None of the runs above conducts charge away, and the charge that evaporating
and subliming particles hand to WRF-ELEC's small ions stays for weeks. In
`box_0e_elec` the ions held 1.9–3.1 kC of positive and 2.0–3.7 kC of negative
charge through the first lunar night and 2.6–5.2 and 0.7–3.0 kC through the
second, half of it below about 30 km and nine-tenths below about 40 km
(`climate/results/crm/elec_box_0e_elec.json`), and in the second lunar day the box's net charge
rose to +10.7 kC at day 40 (at day 40.5, 7.0 kC on snow and 4.1 kC on the
ions), its negative partner having fallen to the ground on precipitation; it
stood near +2 kC for the rest of the run. WRF-ELEC gives the ions' charge to
each new cloud at once. The lunar clear air relaxes net charge in 15 minutes
at 35 km and 35 minutes at the ground (stage 1's conductivity with 100 aerosol
particles per cm³), cloud in about a day.

On 2026-10-05 both windows ran again with leakage (`var9` 1: each cell's net
charge relaxing at σ/ε₀, with the cloud's conductivity where cloud water and
ice reach 0.01 g/kg and the clear air's elsewhere), each changing only that:
`box_0e_elec_leakage_first` beside `box_0e_elec_uncapped` and
`box_0e_elec_leakage` beside `box_0e_elec_ground_rule`
(`climate/results/crm/elec_leakage.json`). Within 33 minutes of day 40.5 the
box's net charge fell from 11.4 kC to 1.2 kC.

| Days | Run | Leakage | Flashes (to ground) | Charge separated, kC | In-cloud flash: median charge, energy | Net charge at the outputs | Field at the ground, median and largest (kV/m) | Away from particles, median and largest (kV/m) |
|---|---|---|---|---|---|---|---|---|
| 10.5–12 | `box_0e_elec_uncapped` | none | 13 (0) | 22.0 + 3.7 | 69 C, 55 GJ | −1.9 to +0.3 kC | 35, 228 | 7.5, 16 |
| | `box_0e_elec_uncapped_corona` | none; point discharge | 6 (0) | 13.5 + 2.4 | 69 C, 54 GJ | −1.0 to −0.2 kC | 21, 123 | 5.0, 15 |
| | `box_0e_elec` (cap) | none | 13 (0) | 18.0 + 2.8 | 80 C, 55 GJ | −1.4 to −0.3 kC | | |
| | `box_0e_elec_leakage_first` | on | 11 (0) | 15.9 + 2.2 | 54 C, 47 GJ | −0.25 to +0.03 kC | 1.4, 39 | 0.23, 3.4 |
| 40.5–42 | `box_0e_elec_ground_rule` | none | 144 (21) | 116.8 + 21.2 | 106 C, 79 GJ | +8.0 to +14.8 kC | 81, 243 | 21, 36 |
| | `box_0e_elec` (5-km rule) | none | 157 (0) | 122.1 + 22.4 | 88 C, 70 GJ | +7.8 to +13.3 kC | | |
| | `box_0e_elec_leakage` | on | 161 (39) | 124.3 + 20.6 | 98 C, 81 GJ | −1.0 to +4.8 kC | 4.7, 48 | 1.6, 19 |

The field at the ground is the largest over land in each ten-step interval
(the coarse run predates that log); away from particles, over the columns
whose lowest level no particles reach.

The leftover charge was an artifact of the missing conduction, and in the
second lunar day it held the field at the ground near 20 kV/m all over the
box, over a hundred times Earth's fair-weather field. With leakage the ions
hold at most 175 C of either sign between storms, the field at the ground
away from storms falls to a median 0.2–1.6 kV/m, below point discharge's
lunar onset of 3.4 kV/m, and the strongest under the storms falls from
228–243 kV/m to 39–48 kV/m, since leakage also acts on the charge rain
carries below cloud base, where the clear air's conductivity applies. The
storms are left much as they were: as many flashes, as much charge separated,
flashes of a similar size starting at the same heights, and while active the
same arrangement of charge (the main negative charge at 38–44 km with
positive charge above and below it). Between storms the outputs hold less
charge, some of it in decaying anvils at 54–60 km.

The window with leakage makes 39 negative ground strikes among 161 flashes,
against 21 among 144 without it: 70 of its flash starts meet WRF-ELEC's
conditions for a negative ground flash against 54, and more of those reach
air of −7 °C. They are smaller: a median 110 C (61–278 C, 10th to 90th
percentile), at most 624 C, and a median 54 GJ. The start potentials do not
show why there are more; with one realization each, part of the difference
may be the weather's.

So the lunar storms' lightning stands without leakage, but the field at the
ground does not, nor the point-discharge comparison of "The first lunar day
again": their fields away from storms were the leftover charge's. The author
made leakage the lunar runs' default on 2026-10-05; the runs above keep the
settings they ran with. The Takahashi windows ran without leakage as well; their arrangement of charge
follows from the law's sign, while their ground strikes carry the same
leftover charge as `box_0e_elec_ground_rule`'s.

Since 2026-10-05 the leakage takes the conductivity column with the aerosol
grown by the air's humidity (atmosphere/electricity README, "Humidity"). Its
ions attach to the particles as in the continuum regime, which overstates the
attachment as the air thins, so the clear air's conductivity aloft is the low
end of what the particles allow: corrected for the transition regime it is
5–12 % higher at 34 km, 12–25 % at 50 km and 23–45 % at 70 km, and the ground
and cloud are hardly changed. `box_0e_elec_corrected` takes the column as it
stands; the correction waits for stage 3's rebuild of the column with the
regional aerosol (the author's decision, 2026-10-06).

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_leakage      # and box_0e_elec_leakage_first
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_leakage --threads 4 --hours 6
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis leakage
```

### The leader's crossing

WRF-ELEC's ground-strike rule counts a downward channel that reaches air of
−7 °C in charge of the matching sign as reaching the ground: 3–5 km short of it
on Earth, 34 km in the lunar air. Its stopping field, 0.15 of breakdown, would
halt every lunar channel 15–25 km up. The literature on leaders (notes on the
data drive, `leader_crossing.md`) stops a leader by neither: a leader carries on
while the potential difference between its tip and the air ahead lasts
(Lalande et al. 2002, their potential budget, Eq. 9; Mazur and Ruhnke's model,
Mazur 2002), and thermalized leaders cross fields of a few hundred V/m
(Gallimberti et al. 2002). What a long crossing costs is the drop along the
channel, its internal field times its length: about 1 kV/m for a thermalized
lightning leader at sea-level density (Mansell 2000; Boggs et al. 2018, who
scale it with density as an assumption they say has not been well studied),
1–10 kV/m in Lalande et al.'s range, near 100 kV/m in a channel that does not
heat. The streamer zone at the tip needs only 0.2–0.4 MV.

`var1` sets that internal field (V/m at 1.225 kg/m³, scaled by the air's
density; `leader_v_m` in a case, 0 for none). Where it is set, a downward
channel that meets WRF-ELEC's rule strikes only if the leader can carry on to
the ground (`terluna_crosses` in `terluna_lightning.F`): the channel holds the
potential where the flash started, as Mansell's model takes it, less the drop
along it from the starting height, and at every level below the tip and at the
ground's zero the tip must stay beyond the air's potential, on the channel's
side, by what its streamer zone needs, 0.225 MV for a positive leader
(Lalande et al.) and 0.4 MV for a negative one, raised 1.3 % per g/m³ of water
vapour as the positive streamer field is (Allen and Boutlendj 1991). Humidity
thus enters, as the streamer field's density scaling cancels in a potential,
but it moves a fraction of a megavolt against drops and potentials of tens to
hundreds. The path lightmsz draws is vertical, so the crossing is the column
below the channel's start.

`box_0e_elec_leader_1kv` and `box_0e_elec_leader_10kv` run
`box_0e_elec_leakage`'s window (days 40.5–42, leakage on) with the crossing at
1 and 10 kV/m, each changing only that (2026-10-05;
`climate/results/crm/elec_leader_crossing.json`):

| Internal field | Flashes (to ground) | Ground strikes' charge: median (10th–90th), largest | Energy: median, largest |
|---|---|---|---|
| none (WRF-ELEC's rule alone) | 161 (39, 24 %) | 110 C (61–278), 624 C | 54 GJ, 302 GJ |
| 1 kV/m | 150 (30, 20 %) | 200 C (69–547), 927 C | 80 GJ, 306 GJ |
| 10 kV/m | 229 (28, 12 %) | 81 C (59–289), 350 C | 47 GJ, 332 GJ |

Each run is a separate realization of the window. Applied to the same
strikes, the 39 of `box_0e_elec_leakage` (from lightmsz's report of each
flash's starting potential, a median 307 MV, and the 32 km of sea-level-density
path below its start, so 32 MV per kV/m), the crossing passes 35 at 1 kV/m, 32
at 3 kV/m, 30 at 5 kV/m, 17 at 10 kV/m and 1 at 20 kV/m; for the 21 of
`box_0e_elec_ground_rule`, 21, 19, 18, 12 and 2. That test leaves out the air's
potential on the way down, which the model's test includes. At the literature's
1–3 kV/m, then, the 34 km of clear air removes up to a fifth of the ground
strikes WRF-ELEC's rule counts, and about a fifth of the lunar flashes
strike the ground; the crossing takes half of them only near 10 kV/m. The
crossing does not decide whether a lunar storm strikes the ground; WRF-ELEC's
Earth-calibrated conditions at the start still do (a potential beyond 50 MV, a
field above 10 kV/m, charge of the matching sign at −7 °C), as do processes the
scheme leaves out: whether a downward end forms and leaves the cloud, whether
the cloud end holds the channel's potential up through the 0.2–0.4 s the
crossing takes at Earth's leader speeds, and current cutoff in a long, branched
channel.

The author made the crossing at 1 kV/m the lunar runs' default on 2026-10-05;
the runs above keep the settings they ran with.

`supercell_elec_leader_1kv` runs the Earth benchmark with the crossing at 1
kV/m beside `supercell_elec_ground_rule`: its flash log is the same to the byte,
7,406 flashes and 85 ground strikes, since the 3–5 km below Earth's −7 °C level
cost a few megavolts against starting potentials beyond 50 MV. (The benchmark
on 4 threads thus repeated itself bit for bit, where the lunar box does not;
"Threads and reproducibility".)

Humidity elsewhere (notes `humidity_corona_breakdown.md`). Point discharge's
onset takes no humidity term: corona inception's threshold moves under 0.15 %
per g/m³ (Ortéga et al. 2007), other measurements disagree in sign, and
vegetation sets the onset more (3 kV/m over dense bushes, 5 kV/m on a barren
ridge). Its linear scaling with density is likely too strong, the square root
of the density ratio being better supported, a difference of at most 0.35
kV/m at the lunar ground, inside the vegetation's range. The breakdown field
takes none either: the runaway threshold follows the air's electron density
(Dwyer 2007), which water vapour at the charging zone's 0.2–0.6 mol % changes
by under 0.1 %, and the humidity-sensitive processes of initiation in cloud are
in WRF-ELEC's Earth calibration already, the lunar charging zone holding
slightly less vapour than Earth's at the same temperature (0.63 against 0.68
mol % at −6 °C).

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_leader_1kv   # and _10kv, supercell_elec_leader_1kv
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_leader_1kv --threads 4 --hours 6
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis leader
```

### The main run

`box_0e_elec_corrected` is stage 2's main run: `box_0e_elec`'s two lunar days
again from day 0 under the lunar rules from the start (the breakdown field's
cap lifted, WRF-ELEC's ground rule with the leader's crossing at 1 kV/m,
leakage through the humidity-corrected conductivity column) and with NSSL's
rain and cloud ice ventilating at their lunar fall speeds. It ran on
2026-10-07 and 08, 27.2 hours on 6 threads with one executable throughout
(`climate/results/crm/elec_box_0e_elec_corrected.json`). Its first attempt,
with the old ventilation, ran to day 7 and is kept as
`box_0e_elec_corrected_before_ventilation_fix`; the two runs' strongest
updrafts agree within about 1 m/s day by day.

The storms charged from day 7, the field reached 100 kV/m on day 9, and the
first flash came at day 10.66. The box flashed 530 times on 12 model days,
199 times in the first lunar day and 331 in the second, the last at day 43.8,
none in the lunar nights. The busiest hour held 45 flashes and the median
hour with any 4. Over the box's 148,000 km² that is 0.022 flashes per km² per
year, about a hundredth of Earth's average (2.7).

| | In cloud | Negative to ground |
|---|---|---|
| Flashes | 424 | 106 (20 %) |
| Start height, median (10th–90th percentile) | 35.8 km (33.8–39.8) | 35.8 km (33.8–37.8) |
| Charge neutralized or lowered, median (10th–90th), largest | 108 C (48–218), 703 C | 149 C (63–434), 1,025 C |
| Energy released, median (10th–90th), largest | 86 GJ (42–185), 758 GJ | 66 GJ (23–162), 325 GJ |
| Channels in cloud, median | 29.8–41.8 km over 144 km² | 37.8–43.8 km over 72 km² |

No positive ground strike occurred. The storms held at most 7.1 kC of
positive and 5.4 kC of negative charge and 4.2 TJ of electrostatic energy,
and the strongest field was 226 kV/m. The original run (`box_0e_elec`, under
its own rules) made 504 flashes, all in cloud under its 5-km rule, with a
median 90 C and 76 GJ. Runs of the same days that differ only by chance moved
flash counts and median charges by 10–20 % (the windows above), so differences
of that size between the two runs say nothing of the changed rules.

Hail did 87 % of the non-inductive charging. In collisions with snow it took
+221 kC and −135 kC, with cloud ice +26 and −23 kC; graupel took +34 and −18
kC with snow and +5.5 and −4.7 kC with ice; inductive charging added +60 and
−7 kC on graupel and +17 and −6 kC on snow. Both lunar days held the normal
arrangement. In the outputs where the main negative region held at least 20 C
(21 in the first lunar day, 27 in the second) it sat at a median 40–42 km
(−13 to −14 °C), carried 90–96 % by snow, with positive charge beneath it at
about 32 km (−4 °C) in most outputs and above it at 58–64 km (−35 to −44 °C)
in 5 of the first lunar day's outputs and 19 of the second's.

With leakage the ions held at most 24 C of positive and 34 C of negative
charge through the lunar nights, against 1.9–5.2 kC in the original run, and
the box ended its two lunar days with no net charge. Within six hours of a
flash the field at the ground over land reached a median of 4.9 kV/m at its
strongest point in each ten-step interval (25.5 kV/m at the 90th percentile,
61.5 at most), and 0.95 kV/m where no particles reached the lowest level (9.1
at most). The cloud ice the radiation reads was cut at 140 µm in a third of
its mass, at a mean radius of 150 µm above the cap, so the radiation sees that
ice 7 % optically thicker than its size gives.

Beside Morrison's storms at the same site and forcing (`box_0e`, the second
lunar day; `climate/results/crm/box_box_0e_elec_corrected.json` against
`box_box_0e.json`), the NSSL run's air at 2 m was 1.0 °C warmer and its dew
point 0.8 °C higher at the same humidity, and it rained 55 % more by day (0.16
against 0.11 mm/h). Its nights were much clearer near the ground: cloud at the
ground covered 29 % of the box's columns through the night against 56 %, and
cloud filled 7 % of the lowest kilometre against 23 %. Above the lowest
kilometre the two schemes' cloud agrees within a few hundredths of cover. The
Open Moon's night fog, which the aerosol study takes from the Morrison rings,
therefore depends on the microphysics scheme by about a factor of two.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_corrected
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_corrected --threads 6 --hours 36
climate/gcm/.venv/bin/python -m climate.crm.elec_analysis box_0e_elec_corrected
climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e_elec_corrected --from-day 29.5
```

### Storm lives

[storm_lives.py](storm_lives.py) follows each storm through two windows of the
main run written every ten minutes, `box_0e_elec_corrected_storms_first` (days
11.5–13.0) and `box_0e_elec_corrected_storms_second` (days 40.5–42.0), each run
on 2026-10-08 from the main run's restart under its settings and executable
(0.7 and 0.8 hours on 6 threads;
`climate/results/crm/storms_box_0e_elec_corrected_storms_first.json` and
`_second.json`). A storm's core is the columns holding at least 1 kg/m² of
graupel and hail, joined across the box's periodic edges and followed from
output to output by its overlap. CM1's threads make each window a new
realization of its days. The windows' flashes follow the main run's for the
first hours and then part from them: 68 against the main run's 129 on the
first window's days, 125 against 174 on the second's. The windows therefore
give the storms' lives, and the main run their numbers.

| | Days 11.5–13.0 | Days 40.5–42.0 |
|---|---|---|
| Storm cores followed (whole lives within the window) | 68 (41) | 71 (31) |
| Storms that flashed; their flashes (to ground) | 6; 68 (18) | 8; 125 (31) |
| First flash after the core appeared, median (range) | 1.05 h (0.73–1.30) | 0.83 h (0.67–1.54) |
| First flash after the core held 10 C, median (range) | 0.38 h (0.23–0.56) | 0.50 h (0.19–1.13) |
| Time from first to last flash, median (range) | 0.57 h (0.13–1.44) | 0.37 h (one flash to 7.2) |
| Flashes a minute in its busiest ten minutes, median (range) | 0.25 (0.1–0.6) | 0.35 (0.1–0.8) |
| Life of a storm that flashed, median (range) | 7.2 h (5.7–8.8) | 7.7 h (5.0–12.2) |
| Life of a storm that did not, median (90th percentile) | 1.7 h (5.0) | 1.7 h (5.8) |
| Its core after the last flash, median (range) | 5.1 h (4.5–6.5) | 5.9 h (3.9–8.6) |

Every flash started inside a core. The first window's track 14 shows the
sequence. Its core appeared at hour 3.3 of the window with an updraft of 7 m/s.
Charge appeared at hour 4.2, and by the first flash, at 4.55, the updraft had
strengthened to 19 m/s, the graupel and hail in the charging zone (−5 to
−30 °C) had grown from 0.07 to 12.5 million tonnes and the top had risen from
40 to 64 km. The storm flashed 26 times in 1.2 hours, 7 of them to ground and at
most five in ten minutes, while the updraft held at 17–21 m/s. The flashes
stopped as the updraft weakened, from 17 m/s at hour 5.8 to 8 m/s twenty
minutes later. The heaviest rain beneath it, 46 mm/h, came 40 minutes after
the last flash. Its charge, at most 1,410 C positive and 1,210 C negative,
drained over the next three hours, and its graupel fell out over five and a
half, the core ending at hour 11.3 with its anvil still at 66 km.

Most storms that flashed had their updraft, their graupel and hail in the
charging zone and their flashes peak together, within half an hour of one
another (medians of 10 minutes or less). Storms that met other cells flashed in
bursts: the second window's longest-flashing storm (its track 48) took in seven
other cores, two of which had flashed, and flashed at hours 25.2–26.2,
28.8–29.7 and 32.5. The heaviest rain beneath a storm follows its most frequent
flashes by a median 0.75–0.9 hours. The lives in the table count the storms
each window holds whole, 4 and 3 of those that flashed.

The storms that flashed were the strongest and deepest: updrafts of 12–36 m/s
at their peak (medians 19 and 21 m/s) against a median 6 m/s for the rest, 3–49
million tonnes of graupel and hail in the charging zone (medians 10 and 19)
against 0.2, and tops at 66–98 km against medians of 54–56 km. The ranges
overlap: some storms without a flash reached 18–23 m/s, 11–18 million tonnes
and 78–88 km. Output by output, a storm's flashes follow its updraft
(Spearman's ρ 0.50–0.54 with the volume rising faster than 5 m/s and with the
strongest updraft) more closely than its graupel and hail (0.29–0.45), and
hardly at all its cloud ice and snow, its top or its area (−0.13 to 0.09). With
cores at 0.5 or 2 kg/m² the medians of the first flash's delay stay at
0.83–1.05 hours and of the charge's at 0.38–0.50, while the flashing storms'
median lives run 5.8–11 hours, longest at the lower threshold, which joins more
of a storm's cells and its trailing precipitation.

Earth's benchmark supercell (`supercell_elec`, output every five minutes, open
edges; `climate/results/crm/storms_supercell_elec.json`), tracked the same way,
flashed within the five minutes in which its core appeared and up to 151
times a minute. At its peak it held 16 million tonnes of graupel and hail in
its charging zone, within the lunar flashing storms' range. The lunar storms
take about an hour from core to first flash and flash 190–1,500 times less
often at their busiest. Both comparisons cross grids, 6 km against 1 km: on
the benchmark a coarser grid gives fewer, larger flashes, and the lunar box at
2 km charged its storms to 73 % of breakdown without a flash ("The fine box"),
so how much of the delay and of the rate the grid sets is open.

```sh
climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup box_0e_elec_corrected_storms_first   # and _second
climate/gcm/.venv/bin/python -m climate.crm.cm1_run run box_0e_elec_corrected_storms_first --threads 6 --hours 8
climate/gcm/.venv/bin/python -m climate.crm.storm_lives box_0e_elec_corrected_storms_first box_0e_elec_corrected_storms_second supercell_elec
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
