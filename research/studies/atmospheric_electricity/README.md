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
air there conducts 1.2–1.9 × 10⁻¹⁴ S/m without aerosol and 0.2–0.3 × 10⁻¹⁴ S/m with
1000 particles per cm³. Inside cloud, where droplets take up the scarce ions, it
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
as Earth's do. If Saunders and Peck's law or Ávila's small charges hold, only the
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
to WRF-ELEC's driver in WRF's order, keeps the charges in CM1's tracers, solves
for the field every step, and runs WRF-ELEC's cylindrical lightning; leakage
through stage 1's conductivity is a switch. The scheme's fall speeds take lunar
gravity as Morrison's do.

1. The builds, their checks (the field solver and discharges against numpy, the
   results against the number of threads) and an Earth benchmark storm, CM1's
   own supercell, compared with published runs of the same charging scheme
   ([S7], [S16], [S26]–[S30]; below). Done.
2. `box_0e_elec`: `box_0e`'s inputs as written, with the NSSL microphysics and
   charging from the start, over two lunar days with three-hourly output and
   twelve-hourly restarts. It compares the NSSL storm with Morrison's at the
   same site and forcing, and gives the first lunar charge structure and
   lightning. A multi-hour run, waiting for the author's go-ahead.
3. Windows of `box_0e_elec`'s second lunar day run again from its restarts with
   output every few minutes, to follow each storm's life, under Takahashi's law,
   with leakage, and with the unbounded breakdown field.

It gives the storms' charge structure, the field by height, flash rates and
types (within cloud and to ground), flash extent, the charge each flash moves
and the nitrogen oxides it makes.

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
