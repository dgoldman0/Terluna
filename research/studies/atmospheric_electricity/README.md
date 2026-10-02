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

## Stage 2: an electrified storm in CM1

CM1 r22 includes the NSSL two-moment microphysics with Mansell's electrification
hooks, the collision terms that charging uses, behind a switch CM1 leaves off.
It has no charge-transfer law, charge transport, field solver or discharge
scheme, and it fixes gravity at 9.8 m/s² in its own code.

1. Patch the NSSL scheme for lunar gravity (fall speeds, drag and collection)
   and run it beside the Morrison box at the same site and forcing, to see how
   much the storm itself changes with the scheme.
2. Add charge carried on each particle type, a laboratory charging law, a
   Poisson solver for the field and a discharge scheme, following WRF-ELEC
   (Mansell's electrification in WRF; [S2], [S3]) where its licence allows.
   Reproduce an Earth benchmark storm against WRF-ELEC first.
3. Run the equatorial box electrified from saved restarts, with output frequent
   enough to follow each storm's life.

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
