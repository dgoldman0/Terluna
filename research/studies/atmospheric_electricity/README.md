# Atmospheric electricity on the Open Moon

The author set this study's goal on 2026-10-02: to understand the Open Moon's
atmospheric electricity about as well as the project understands its other
subjects of similar importance. That covers the charge its storms build and the
fields they make; its lightning (how often, of what kinds, how long and how
energetic, and the nitrogen it fixes); the conductivity of its air from the
ground to the conducting upper atmosphere; its fair-weather field and global
circuit; and its transient luminous events (sprites, halos, jets and elves).

The study draws on cloud microphysics in [climate/crm](../../../climate/crm/README.md),
the air column and its ionization in [atmosphere/](../../../atmosphere/README.md),
and the particles the protection lets through in
[protection/](../../../protection/README.md). Models and diagnostics live in those
domains; this folder holds the comparison between them and what it means.

## What sets the answer

Three features of the Open Moon pull in different directions.

**Particles fall slowly.** CM1's Morrison scheme, patched for lunar gravity,
scales each particle's fall speed at a given size by (g/g⊕)^((b+1)/3), where b is
the exponent of its speed law: graupel falls at 0.44 times Earth's speed, snow
0.43, cloud ice 0.30, rain 0.34 and cloud droplets 0.17. Graupel therefore meets
ice crystals at about 1 m/s. In laboratory experiments the charge a graupel
particle gains in each bouncing collision with an ice crystal rises steeply with
impact speed, and collisions occur in proportion to the speed, so the same
particles separate charge between about five and thirty times more slowly than
on Earth, depending on how steep that rise is. How charge transfer behaves at
these slow impacts is the largest uncertainty.

**Storms are deep and slow.** The air's scale height is about six times Earth's,
so convection and its mixed-phase layers reach about six times as high. A storm
six times larger in every direction makes the same electric field from a sixth
of the charge density, and a slower storm gives charge longer to build. Whether
lunar storms are as much wider as they are deeper is measured in the boxes.

**The air conducts little.** The column holds 7.25 times Earth's mass of air (1.2
atm at a sixth of the gravity), so most of the troposphere lies deeper in the
air than Earth's sea level. The nucleonic and electromagnetic cascades of cosmic
rays die away above it, and muons carry the ionization there. The Moon has no
global magnetic field, and the cosmic rays whose muons reach the lower air pass
the regional magnets' fields of hundreds of nanotesla almost unbent. Weakly
ionized air lets storm charge build longer before it leaks away; over land, the
ground's radioactivity adds ionization near the surface.

These effects are of similar size, so the answer turns on the particle
populations, the storms' shape and lifetime, and the conductivity, which stage 1
measures.

## Stage 1: the storms already run, and the conducting air

The CM1 runs saved graupel, cloud ice, snow and rain mass and number, cloud
water, vertical wind, potential temperature and pressure every three model hours
over two lunar days: the equatorial boxes `box_0e` (385 km square) and
`box_0e_small` (192 km), `box_0e_small`'s twin with Earth's fall speeds
(`box_0e_small_earth_fall`), and the flat rings `ring_equator`, `ring_70_45e`
and `ring_70_135e`. Each run's second lunar day gives:

- graupel, cloud ice and snow mass and number by height and temperature, with
  sizes from the scheme's exponential distributions (graupel at 400 kg/m³);
- the charging zone, the air between 0 and −40 °C where graupel, cloud ice or
  snow and supercooled cloud water meet: its depth per unit of ground by
  temperature, its share of the time, its course through the lunar day, and the
  width of its storms in the boxes, with the thresholds varied;
- updrafts through the charging zone, the mass-weighted fall speeds of graupel
  and ice and their difference, and the graupel's bulk residence time (its mass
  aloft over the rate it falls through the melting level);
- the non-inductive charging rate from published laboratory laws applied to the
  model's particles and impact speeds, with the speed dependence varied across
  the laboratory range, and its integral up each column, the storm's generator
  current density.

The fall-speed pair shows what slower settling alone does to each of these.

Alongside, an ionization and conductivity column for the Open Moon's air in
[atmosphere/](../../../atmosphere/README.md): cosmic-ray ionization by depth at
zero geomagnetic cutoff, at solar minimum and maximum, including the muons below
Earth's sea-level depth; ionization from the ground over land; ion recombination
and attachment to aerosol and cloud droplets; and ion mobility at the air's
density. It gives the conductivity σ and the charge relaxation time ε₀/σ by
height, in clear air and in cloud. Titan, whose low-gravity nitrogen atmosphere
holds about ten times Earth's column of air, is the measured analogue: the
Huygens probe measured its conductivity on the way down.

A storm can make lightning where its generator current density exceeds the
current the air conducts at the breakdown field (σ times that field, which
scales with air density). Comparing the two, and the time to reach breakdown
with the storms' lifetimes, shows which charging cases produce lightning, on the
Moon and in the Earth-fall twin.

### What stage 1 found

The Open Moon has lightning. Its air holds storm charge so well that some
storm columns reach the breakdown field under every charging law tried; how
many, and how fast, turns on how much charge slow, lightly rimed graupel
separates, which the laboratory has barely measured
([stage1.py](stage1.py), [results/stage1.json](results/stage1.json); 2026-10-03).

**The storms** ([climate/crm](../../../climate/crm/README.md), "What the storms
hold for charging"). The 0 °C level stands near 25 km and −40 °C near 60 km,
five to six times their heights in Earth's tropics. The charging zone is
shallow and warm: four-fifths of it lies between 0 and −10 °C, a fifth between
−10 and −20 °C, almost none colder, because the supercooled cloud water thins
from 0.10–0.17 g/m³ to 0.03–0.05 g/m³. Graupel is plentiful (1.7–2.5 g/m³ in the
boxes), 5–7 mm across by mass, and stays aloft for 4–5 hours; cloud ice is
scarce, 2–10 crystals per litre, and snow carries most of the ice graupel
meets. Graupel strikes ice at 1.3–1.5 m/s and rimes at 0.1–0.2 g m⁻² s⁻¹, the
bottom of the range the laboratory measured.

**Charging** ([charging.py](../../../atmosphere/electricity/charging.py)), in the
equatorial box (385 km) and with the same particles falling at Earth's speeds,
with the fall-speed pair's total charge separation at lunar against Earth's fall
speeds:

| Law (laboratory impact speeds) | Charging in the zone, pC m⁻³ s⁻¹ | At Earth's fall speeds | Storm-column current, top tenth / strongest, nA/m² | Fall-speed pair, lunar ÷ Earth's |
|---|---|---|---|---|
| Saunders and Peck (3–14 m/s) | 0.0003–0.0015 | 0.04–0.10 | 0.036 / 0.32 | 0.10 |
| Takahashi (9 m/s) | 0.27–0.71 | 0.89–2.0 | 6.9 / 93 | 1.35 |
| Pradeep Kumar et al. 2024 (1.2 m/s) | 0.02–0.09 | 0.05–0.21 | 1.1 / 15 | 1.29 |
| Ávila et al. 2013 (1–3 m/s) | 0.0006–0.0014 | 0.0014–0.003 | 0.024 / 0.33 | 1.24 |

Saunders and Peck's law nearly stops because most slow graupel rimes below its
threshold of 0.1 g m⁻² s⁻¹. Takahashi's charge grows only linearly with impact
speed, and the two experiments near lunar impact speeds bracket what slow bounces
carry, 0.2 fC (Ávila et al.) to 6–8 fC (Pradeep Kumar et al.). Slower settling on
its own cuts total charging tenfold under Saunders and Peck and raises it by a
quarter to a third under the other laws, because six times more graupel stays
aloft.

**The conducting air** ([conductivity.py](../../../atmosphere/electricity/conductivity.py),
[muons.py](../../../atmosphere/electricity/muons.py)). Below Earth's sea-level
depth only muons and their decay electrons ionize the air. Traced with MCEq
through the Open Moon's column, they make 0.05–0.06 ion pairs per cm³ per second
at the charging zone (2,700–4,100 g/cm², 30–50 km), against 12 in an Earth storm's
charging zone at 6 km, and 0.04 at the ground, against 2 over Earth's seas. Clear
air there conducts 1.2–1.9 × 10⁻¹⁴ S/m without aerosol and 0.10–0.13 × 10⁻¹⁴ S/m with
1000 particles per cm³ swollen by the air's humidity (0.2–0.3 × 10⁻¹⁴ S/m before
the swelling and the attachment's pressure scaling were added on 2026-10-05;
atmosphere/electricity README, "Humidity"). Inside cloud, where droplets take up the scarce ions, it
conducts 1 × 10⁻¹⁶ S/m, so separated charge leaks away over about a day; in an
Earth storm, over half an hour. The solar cycle leaves this unchanged below about
90 km, since the muons' parent cosmic rays lie far above the energies it
modulates. The same ion chain reproduces Earth's measured fair-weather
conductivity (Gringel 1978) within 2–8 % from 5 to 30 km, and MCEq Earth's
sea-level muons within 6–13 % above 10 GeV/c; it falls 30 % short at 1–10 GeV/c,
energies that never reach the lunar troposphere.

**Breakdown.** The runaway threshold at the zone's air density is 178–188 kV/m.
Lunar cloud conducts 0.016–0.019 nA/m² at that field, against 0.8 nA/m² in Earth's
storm clouds, whose charging currents of roughly 20–1000 nA/m² reach breakdown in
seconds to minutes. In the lunar storms:

| Law | Storm columns reaching breakdown, box (equatorial ring) | Time to breakdown in the box, top-tenth column / strongest |
|---|---|---|
| Saunders and Peck | 21 % (2 %) | 16 h / 1.4 h |
| Takahashi | 67 % (28 %) | 4 min / 17 s |
| Pradeep Kumar et al. | 62 % (15 %) | 24 min / 2 min |
| Ávila et al. | 13 % (2 %) | 31 h / 1.3 h |

With Earth's fall speeds (the twin), Saunders and Peck's law brings 55 % of
storm columns to breakdown, the top tenth in 45 minutes. If Takahashi's law or
the 1.2 m/s measurements hold, most lunar storms make lightning within minutes,
as Earth's do. (Stage 2 overturns this: in the NSSL storms hail does
most of Saunders and Peck's charging, and Takahashi's law changes the storms'
polarity and leaves their flash rate as it was; item 6 below.) If Saunders and Peck's law or Ávila's small charges hold, only the
strongest storm columns do, after one to sixteen hours, which may outlast a
single convective cell; the three-hourly output cannot tell. The rings' weaker
updrafts, in two dimensions, put fewer columns over the threshold.

The comparison is a one-dimensional bound: all separated charge forms one layer,
and it leaks through cloud. Stage 2 replaces it with charge carried on the
particles, a field solver and discharges, and follows each cell's life. The
decisive laboratory quantity is the charge per bounce at impact speeds of
1–1.5 m/s and rime accretion rates of 0.1–0.3 g m⁻² s⁻¹, within reach of existing
wind tunnels.

```sh
climate/gcm/.venv/bin/python -m climate.crm.mixed_phase_analysis box_0e          # each storm run
/media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons moon   # and earth
climate/gcm/.venv/bin/python -m atmosphere.electricity.conductivity
climate/gcm/.venv/bin/python -m research.studies.atmospheric_electricity.stage1
```

## Stage 2: an electrified storm in CM1

CM1 r22 includes the NSSL two-moment microphysics with Mansell's electrification
hooks, the collision terms that charging uses, behind a switch CM1 leaves off.
It has no charge-transfer law, charge transport, field solver or discharge
scheme, and it fixes gravity at 9.8 m/s² in its own code.

The electrified builds `earth_g_omp_elec` and `moon_omp_elec` (2026-10-03;
[climate/crm](../../../climate/crm/README.md), "Electrified storms") put
WRF-ELEC's NSSL module, which carries charge through every microphysical process
([S8]), in place of CM1's copy of the scheme. Terluna's code passes CM1's arrays
to WRF-ELEC's driver in WRF's order and keeps the charges in CM1's tracers. As
WRF-ELEC's driver does, it runs the sedimentation in sub-steps, each followed by
a solve for the field and by lightning: WRF-ELEC's branched flashes
(MacGorman, Straka and Ziegler 2001 [S27]), in cloud and to ground, with the
charge each neutralizes, its channels and the nitrogen oxides it makes, or its
cylinders. Leakage through stage 1's conductivity and WRF-ELEC's screening
layers at cloud edges are switches. The scheme's fall speeds take lunar gravity as
Morrison's do.

1. The builds, their checks (the field solver and the lightning against numpy
   and against charge laid out by hand, the results against the number of
   threads and across restarts) and an Earth benchmark storm, CM1's own
   supercell, compared with published runs of the same charging and lightning
   schemes ([S7], [S16], [S26]–[S30]; below). Done on 2026-10-03, and again on
   2026-10-04 with the sub-steps and the branched flashes, at 2 km and in WRF-ELEC's
   own test settings.
2. `box_0e_elec`: `box_0e`'s inputs as written, with the NSSL microphysics,
   charging and branched lightning from the start, over two lunar days with
   three-hourly output and twelve-hourly restarts. Its sub-steps are 6.8 s,
   scaled to its 2-km layers and slower graupel, and a downward channel strikes
   the ground when it comes within 5 km of it. It compares the NSSL storm with
   Morrison's at the same site and forcing, and gives the first lunar charge
   structure and lightning. Running from 2026-10-04 at the author's go-ahead; a
   first start without the sub-steps was stopped at day 4, before its storms had
   charged.
3. `box_0e_elec_fine`, the same site at 2 km over a box 128 km square, started
   from the coarse run's averaged air a few hours before a stormy window and run
   for two model days with output every 15 minutes, to resolve the storm cells
   and draw the flashes on a finer grid. Set up on 2026-10-04 at the author's
   request and run the same day at the author's go-ahead, from the coarse run's
   day 10.75, four hours before its first flash. Its storms took most of the
   first model day to grow from the averaged air, then matched the coarse box's
   in updrafts and graupel and charged to about 1,500 C of each sign, but no
   flash struck: the field peaked at 73 % of breakdown. At the coarse box's rate
   of flashing storms per unit area, a box this size sees none in that time
   55 % of the time, so the effect of the finer grid on the lightning is
   still open (climate/crm README, "The fine box"). At the author's direction
   it ran on toward the coarse box's busiest lightning, but its deep convection
   collapsed after coarse day 13.25 while the coarse box's carried on, and the
   author stopped it at coarse day 13.75 after three days without a flash. The
   2-km lightning stays unmeasured; it needs a box wide enough for several
   storms.
4. The first lunar day's storms again from the coarse run's day-10.5 restart
   (climate/crm README, "The first lunar day again"). WRF-ELEC caps the
   breakdown field at 180 kV/m, a cap that on Earth applies only below about
   4.5 km; in the lunar air it applies below about 36 km, where half the first
   lunar day's flashes started. The author lifted it for the lunar runs on
   2026-10-04, and the coarse run runs without it from day 18;
   `box_0e_elec_uncapped` gives its first lunar day without it.
   `box_0e_elec_uncapped_corona` adds point discharge from the ground (Standler
   and Winn 1979), which WRF-ELEC lacks and the author left to the agent's
   judgement: without it the field at the ground under the storms reaches
   100–220 kV/m, where on Earth the ground's plants and points hold it near
   5–12 kV/m. The pair shows what that discharge changes before other runs take
   it up. A 30-minute check held the ground's field at the onset where no
   particles reach the lowest level, but not under rain, where WRF-ELEC hands
   the ground's ions to the rain at once and the scheme's ions do not climb in
   the field. Both ran on 2026-10-04. Without the cap a third of the flashes
   started below 34.8 km, against half under it, with 123 flashes against 162
   and a tenth more charge each. Point discharge cut the dry ground's
   strongest field from a median 11 to 7 kV/m but left the strongest under
   storms near 34 kV/m (40 without it). On several threads CM1 does not repeat
   itself bit for bit, so each pair compares realizations, and the counts lie
   within what the weather alone could change (climate/crm README).
5. Ground strikes (climate/crm README, "Ground strikes"). No lunar flash struck
   the ground until 2026-10-05, an artifact: the lunar ground rule kept
   WRF-ELEC's demand for matching charge 5 km up, where the storms have none,
   and the Earth benchmarks ran with WRF-ELEC's rule switched off. The storms
   hold their charge as Earth's do, with −7 °C just beneath the main negative
   charge, so WRF-ELEC's own rule fits them. Under it the busiest lunar storms
   (days 40.5–42) make 21 negative ground strikes among 144 flashes, a median
   195 C and 69 GJ each (up to 889 C and 437 GJ), and the Earth benchmark makes
   1.1 % ground strikes from 36.8 min, as published runs do. The rule takes a
   channel at −7 °C, 34 km up, to reach the ground, so this is an upper bound;
   WRF-ELEC's stopping field, which halts every channel 15–25 km up, gives the
   lower bound of none. A rule for a leader's crossing of that gap comes next.
6. The charging law (climate/crm README, "The charging law"). Both windows,
   days 10.5–12 and 40.5–42, ran again on 2026-10-05 under Takahashi's law,
   each changing only the law, with the runs' three-hourly output. The storms
   make as many flashes as under Saunders and Peck's law (14 against 6–13,
   and 139 against 144–157), separate as much charge and flash no sooner.
   Stage 1's two hundred to five hundred times came from the Morrison storms'
   graupel; in the NSSL storms hail, which falls faster and rimes above
   Saunders and Peck's threshold, does 82–88 % of that law's charging. The law
   sets the storms' polarity instead. At the lunar storms' 0.03–0.17 g/m³ of
   cloud water, Takahashi's table charges rimed ice positively at every
   temperature, and the storms take an inverted dipole: the main negative
   charge at 48–58 km (−23 to −35 °C) above positive charge at 32–36 km, where
   Saunders and Peck's law gives the normal arrangement with the main negative
   charge at 38–44 km. The laboratory measurements at 1.2–1.8 m/s also charge
   rimed ice positively ([S12], [S13]), but they bear only on the graupel: the
   storms' hail, which does most of the charging, falls at 4.8–7.7 m/s and
   strikes ice and snow at 4.4–7.6 m/s. Where it charges, at −5 to −16 °C in
   0.02–0.18 g/m³ of cloud water, Saunders and Peck's law charges it both ways
   and Takahashi's table positively, and the polarity rests on that
   disagreement. Under Takahashi's law the flashes start
   higher and release one and a half to two times the energy (a median 106–119
   GJ against 54–79), and the second window's ten ground strikes, an upper
   bound by WRF-ELEC's rule, brought down a median 204 C, the largest 1,489 C.
7. Leakage (climate/crm README, "Leakage"). With nothing to conduct it, the
   charge evaporating cloud leaves on WRF-ELEC's small ions stayed for weeks:
   1.9–3.7 kC of each sign through the first lunar night and 0.7–5.2 kC
   through the second, and 11.4 kC net at day
   40.5, which held the field at the ground near 20 kV/m all over the box.
   Both windows ran again on 2026-10-05 with leakage through stage 1's
   conductivity, each changing only that. The leftover charge went within the
   hour, and the field at the ground away from storms fell to a median 0.2–1.6
   kV/m and under them from 228–243 to 39–48 kV/m at most. The lightning stayed
   as it was (11 flashes against 6–13, 161 against 144–157), while the busy
   window's ground strikes rose to 39 against 21, smaller ones (a median 110 C,
   at most 624 C; item 8 adds the leader's crossing). The fields at the ground of item 4's point-discharge
   comparison were the leftover charge's. The author made leakage the lunar
   runs' default the same day.
8. Humidity and the leader's crossing (2026-10-05; climate/crm README, "The
   leader's crossing"; atmosphere/electricity README, "Humidity"). The lunar
   air is humid, 84 % relative humidity in clear air at the ground and 52–64 %
   from 10 to 50 km. Its aerosol now swells with that humidity (κ-Köhler,
   [S41], with κ 0.3 and 0.1–1.0 as a bracket, [S45]) and takes up ions by a coefficient that scales with pressure ([S42]),
   which lowers clear air's conductivity by 12 % at the ground and 22 % at 34
   km. That scaling is the continuum regime's and makes the values aloft the
   low end, 5–12 % low at 34 km and 23–45 % at 70 km; they stand until stage 3
   rebuilds the column (author, 2026-10-06). The ions' mobility and
   recombination take no humidity term ([S43],
   [S44]), nor do point discharge's onset ([S51]) or the breakdown field
   ([S52]). A ground strike can now need the leader to cross the 34 km below
   the −7 °C level: it carries on while its tip keeps the potential its
   streamer zone needs beyond the air's ([S46], [S47]), the channel losing its
   internal field times its length, about 1 kV/m at sea-level density for a
   thermalized leader ([S48], [S49]; 1–10 kV/m), and humidity raising the
   streamer zone's need by 1.3 % per g/m³ ([S50]), a fraction of a megavolt.
   At 1–3 kV/m the crossing passes 82–100 % of the ground strikes WRF-ELEC's
   rule counts, and about a fifth of the lunar flashes strike the ground (30
   of 150 at 1 kV/m); it takes half of them near 10 kV/m (28 of 229) and
   nearly all at 20 kV/m. On Earth's benchmark storm at 1 kV/m it changes
   nothing: the same 85 ground strikes among 7,406 flashes. The author made
   the crossing at 1 kV/m the lunar runs' default the same day. So the ground strikes' number turns on WRF-ELEC's
   Earth-calibrated conditions at the start and on what the scheme leaves out
   (a downward end leaving the cloud, the cloud holding the channel's
   potential through the 0.2–0.4 s crossing, current cutoff), and on the
   internal field above about 5 kV/m.
9. The aerosol near the ground ([open_moon_aerosol](../open_moon_aerosol/README.md),
   2026-10-05), built from the ecology register's landscapes and organisms and the
   climate runs: the seas, forests, fog-desert plains and polar lands give the
   ground a conductivity of about 2.2–2.7 × 10⁻¹⁵ S/m (0.6–4.7 across its cases),
   charge relaxing in about an hour, 2–4 % of that in fog, and 80–160 cloud nuclei
   per cm³ at 0.3 % supersaturation by day, about the storm runs' 100.
10. Windows of `box_0e_elec`'s second lunar day run again from its restarts with
   output every few minutes, to follow each storm's life.

It gives the storms' charge structure, the field by height, flash rates and
types (within cloud and to ground), flash extent, the charge each flash moves
and the nitrogen oxides it makes.

### Thunder

[thunder.py](thunder.py) ([results/thunder.json](results/thunder.json),
2026-10-05) traces thunder from the 465 flashes the box made with the cap lifted
(its second lunar day and the first lunar day's rerun) to the ground, through the
box's own air at the hours that flashed and in the design air's composition, in
eight directions ([atmosphere/electricity/thunder.py](../../../atmosphere/electricity/thunder.py)).
A flash's sound is 0.18 % of the electrostatic energy it releases (Holmes et
al. 1971; tried at a tenth and ten times that), spread along its channel, in a
spectrum peaked at Few's (1969) frequency.

A typical flash (78 GJ, its channel 30–42 km up) peaks at 27 Hz, a
large one (the top tenth by energy, 232 GJ, 30–54 km) at 19 Hz. Directly
below, its thunder arrives 87 s after the flash and lasts about 35 s (69 s for
the large flash): 89 dB at 31.5 Hz and 85 dB at 63 Hz in its loudest second,
62–64 dBA, a deep rumble with almost nothing above 250 Hz. Across all 465 flashes it
reaches 56–68 dBA there (10th to 90th percentile). Eastward it falls to 57 dBA
at 25 km, 50 at 50 km, 42 at 100 km and 31 at 200 km; in the other directions
the silent zone begins before 200 km. The air cools only 1.0 K per
km over the lowest 30 km, so sound bends upward gently, and the silent zone
beyond the ray that grazes the ground starts 175–225 km out; the box's weak
eastward wind (up to 7 m/s aloft) carries it about 30 km farther east than
west. Within that edge the thunder stays above the threshold of hearing; the
edge, not the flash's strength, sets the range. It stays above an ordinary
daytime background of 45 dBA to 54–113 km (median 80 km), and at a tenth or ten
times the acoustic share to about 30 or about 180 km.

The same calculation for Earth's benchmark flash (1.2 GJ over 4.75–10.75 km)
in the US standard atmosphere puts its thunder at 76.5 dBA directly below,
arriving after 15 s and lasting 17 s, audible to 41 km and above 45 dBA to
31 km, where thunder on Earth is seldom heard beyond 15–25 km. The calculation
leaves out wind gusts and turbulence near the ground, soft ground and the
storm's own noise, so its ranges run long by a similar factor on the Moon,
which would put lunar thunder above a daytime background to roughly 50 km and
audible in quiet to 100–150 km. Under the storm itself lunar thunder is about
14 dB quieter on the A-weighted scale than Earth's, being four to five times
farther overhead, but it is far deeper and longer, and it carries several times
as far. Thirteen of the flashes, starting at 50–68 km on the second lunar day,
logged no channel though they neutralized 5–72 C; they are placed at their
starting points, and why the log missed their channels is not yet traced.

## Stage 3: the global circuit and transient luminous events

The fair-weather field and global circuit follow from stage 2's storm currents,
weighted by how often each kind of storm occurs over the Moon (from the GCM and
the rings), the column resistances from stage 1, and a conducting upper
boundary from the ionization column carried upward. Transient luminous events
follow from the charge that large flashes move: a quasi-electrostatic
calculation for sprites and halos [S4], an electromagnetic-pulse calculation for
elves, and a discharge-propagation treatment for jets [S5], with thresholds set
by the gas density at the Open Moon's heights.

## How the work is done

Matched boxes differ only in the treatment under test: site, surface, forcing,
grid and starting state stay the same, and comparisons use each run's second
lunar day. Grid spacing and output frequency are checked on short windows from
saved restarts; every run kept its twelve-hourly restarts. Runs fit the laptop:
8 cores and 31 GB shared with other work, in sessions of hours with rests, and
no rented compute. CM1 and the GCM disagree about the air near the ground over
land ([climate/crm](../../../climate/crm/README.md)); each stage carries that as a
range of storm environments.

Products follow [shared/README.md](../../../shared/README.md). Sources and what was
read of them are in [sources.json](sources.json).
