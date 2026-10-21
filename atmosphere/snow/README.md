# Snow

The current [repair and literature review](repair.md) preserves these earlier tables while testing an independent conservative solver and physical constraints on geometry and breakup.

The [21 October analysis](analysis.md) extends this work to snow-conditioned cloud motion, corrected aggregate geometry, full size distributions and lunar snowstorm optics. The tables below retain the 20 October baseline; the source-definition corrections and current limitations are recorded in the new analysis. The inherited fractional-particle bookkeeping fails an arbitrary-area invariance requirement; read the [numerical audit](numerical_audit.md) before interpreting either set of size estimates.

Snow in the Open Moon's sky winter: how snowflakes fall, grow and melt at a sixth of Earth's gravity, in the air about
the freezing level, beside the raindrop solver ([../rain](../rain/README.md)). The Moon's winter weather is built up in
the sky, and on 20 October the author asked how large snowflakes can grow there and how they fall, as the second half
of the work on the high cloud's crystals ([register](../../research/decisions.md), 20 October).

| File | What it holds |
|---|---|
| [flakes.py](flakes.py) | The snowflake's mass, area and fall law, its growth from the vapour, its melting, the winter columns and the product |
| [results/snow.json](results/snow.json) | Fall speeds by size, growth across the winter cloud, crossing times and melting depths, Moon and Earth |
| [aggregation.py](aggregation.py) | A column of falling snow: super-particles that grow, stick together and break, and the snow they bring to the freezing level |
| [results/aggregation.json](results/aggregation.json) | The earlier arriving-snow sizes by supply and crystal numbers, retained as a baseline |
| [population.py](population.py), [results/population.json](results/population.json) | Geometry controls and number, mass, area and flux distributions, including cloud startup |
| [growth_depth.py](growth_depth.py), [results/growth_depth.json](results/growth_depth.json) | Shorter and warmer continuously saturated growth columns |
| [numerical_audit.py](numerical_audit.py), [numerical_audit.md](numerical_audit.md) | Fractional-weight bookkeeping, arbitrary-area controls and collision-batch diagnostics |
| [sectional.py](sectional.py), [repair.py](repair.py), [results/repair.json](results/repair.json) | Conservative independent solver and completed domain/refinement diagnostics, with explicit adequacy flags |
| [kernel_diagnostics.py](kernel_diagnostics.py), [results/kernel_diagnostics.json](results/kernel_diagnostics.json) | Fixed-geometry encounter scaling and the installed chipping rule's mass limits |
| [geometry_literature.md](geometry_literature.md), [strength_literature.md](strength_literature.md) | Targeted primary evidence for evolving geometry, interception, ice contacts and aggregate failure |
| [turbulence.py](turbulence.py), [turbulence.md](turbulence.md) | Archived cloud motions conditioned on snow and growth environments; separate flake-scale loading scenarios |
| [habits.py](habits.py), [habits.md](habits.md) | Crystal-form context, source-definition corrections and representative aggregate geometries |
| [assessment.py](assessment.py) | Audit-linked evidence metadata on the unchanged numerical products, with separate provenance |
| [optics.py](optics.py), [optics.md](optics.md) | Snow extinction, height-resolved lunar illumination, thin-path radiance and multiple-scattering controls |

```sh
python -m atmosphere.snow.flakes        # seconds: results/snow.json
python -m atmosphere.snow.aggregation   # about 20 minutes: results/aggregation.json
```

## The physics

The baseline treats a snowflake as an aggregate of ice crystals, with m = 0.0185 D^1.9 kg and projected area 0.35 of the circle across D. The [source review](habits.md) found that the mass coefficient uses a mean image dimension rather than Dmax, and that the cited area fits use an ellipse rather than that circle. The 2026 source is Seifert and co-authors, using Leinonen's aggregation method. These baseline closures are preserved for comparison; the new population analysis uses the corrected dimension convention and explicit area sensitivities. It falls by Heymsfield and
Westbrook's (2010) law, in which gravity enters through the modified Best number alone,
X* = 8 m g ρ / (π Ar^½ η²), with Re = (δ0²/4)((1 + 4 X*^½ / (δ0² C0^½))^½ − 1)², δ0 = 8.0, C0 = 0.35 and the speed
η Re / (ρ D).

It grows from the vapour by diffusion to a capacitance of 0.25 D (Westbrook, Hogan and Illingworth 2008), ventilated as
Hall and Pruppacher (1976) give it for ice, against Murphy and Koop's (2005) saturation over ice and supercooled water.
Below the freezing level it melts by the heat it takes from the air less what its wet surface's evaporation carries
away (Pruppacher and Klett 2010, chapter 16), at the same capacitance and ventilation. The air's viscosity,
conductivity and vapour diffusivity are the raindrop solver's (`drops.transport`) for each world's air.

A winter column is a straight temperature profile through the freezing level, its pressure hydrostatic:

| Column | Freezing level | Pressure there | Air there | Lapse |
|---|---|---|---|---|
| Moon: the resort's sky ladder at 301 K | 29.2 km above sea level | 659 hPa | 0.83 kg/m³ | 0.92 K a km (fitted over 20–40 km) |
| Earth: a winter column | about 1.5 km | 850 hPa | 1.08 kg/m³ | 6.5 K a km |

The Moon's ladder is the climate targets' (research/studies/climate_targets), read as the rain products read it
([targets.py](../rain/targets.py)). The cloud runs from a top temperature down to the freezing level at liquid
saturation, a mixed-phase cloud where snow grows fastest.

## What it finds

**Snow drifts down.** Each world's flakes in the air at its freezing level, m/s:

| Flake | 0.2 mm | 0.5 mm | 1 mm | 2 mm | 5 mm | 1 cm | 2 cm | 5 cm |
|---|---|---|---|---|---|---|---|---|
| Moon | 0.12 | 0.22 | 0.31 | 0.40 | 0.51 | 0.57 | 0.61 | 0.64 |
| Earth | 0.54 | 0.82 | 1.02 | 1.19 | 1.35 | 1.43 | 1.47 | 1.50 |

Lunar flakes fall at a quarter to two-fifths of Earth's speed; the largest drift down at about 0.6 m/s, 2.3 km an hour.

**The winter cloud is deep, and snow spends hours in it.** The air cools only 0.92 K a km near the Moon's freezing
level, so the same span of cloud temperatures stands seven times deeper than on Earth and holds about five times the
air. From each cloud top down to the freezing level:

| Cloud top | Moon: depth; air | Moon: a 5 mm flake's crossing | Earth: depth; air | Earth: crossing |
|---|---|---|---|---|
| −10 °C | 10.9 km; 8,275 kg/m² | 5.8 hours | 1.5 km; 1,543 kg/m² | 0.31 hours |
| −15 °C | 16.3 km; 11,851 kg/m² | 8.5 hours | 2.3 km; 2,226 kg/m² | 0.45 hours |
| −20 °C | 21.7 km; 15,089 kg/m² | 11.1 hours | 3.1 km; 2,856 kg/m² | 0.59 hours |

**So snow grows far more.** A 0.2 mm crystal crossing the cloud gains from the vapour alone 2.2, 12.5 and 43 mg on the
Moon from the three tops, against 0.010, 0.025 and 0.049 mg on Earth: 200–900 times the mass, as a growing crystal
grows faster the larger it is. With five times the air in the cloud, a flake also has about five times the cloud's
small crystals to sweep up at one mixing ratio. A single flake sweeping a cloud that never runs out grows without
bound; how large flakes grow in a whole population, with a breakup law, follows below.

**Snow lasts far below the freezing level.** Below it the air warms only 0.92 K for each km a flake falls, so flakes
melt slowly. Where all of a flake's ice has melted, metres below the freezing level:

| Flake | Moon, air at 90% | Moon, at 70% | Earth, at 90% | Earth, at 70% |
|---|---|---|---|---|
| 2 mm | 1,045 | 2,849 | 228 | 450 |
| 5 mm | 1,137 | 2,942 | 266 | 489 |
| 1 cm | 1,212 | 3,019 | 297 | 522 |
| 2 cm | 1,294 | 3,103 | 332 | 557 |

In drier air a flake's own evaporation keeps its surface below freezing for longer. The sky winter's snowfall is a
band kilometres deep: snow forming in a cloud above the freezing level reaches decks 1–3 km below it.

## How large flakes grow

On 20 October the author asked for a breakup law, to say how large the flakes grow. [aggregation.py](aggregation.py)
follows a column of falling snow. Crystals of 0.1–0.4 mm enter at a −15 °C cloud top, 10 or 100 a litre, and fall at
their own speeds to the freezing level, as super-particles each standing for many real ones (Shima and co-authors
2009). They grow from the vapour as far as the cloud condenses it for a set snowfall, spread through the cloud's
depth. They meet the others at their height, sweeping the area of their projected areas (Mitchell 1988) at their
difference in speed, and a meeting sticks with Lin, Farley and Orville's efficiency. One that does not stick breaks
fragments off the smaller particle by Phillips, Yano and Khain's (2017) law for snow, which Sotiropoulou and
co-authors set out (Atmospheric Chemistry and Physics Discussions, acp-2020-786, appendix B): up to 100 fragments of
0.1 mm a collision, more the harder it is. The snow arriving at the freezing level, its sizes by mass and the largest
flake landing at least once an hour on each square metre:

| Snowfall | Crystals | Moon: median; 99% below; largest hourly | Earth: median; 99% below; largest hourly |
|---|---|---|---|
| 0.3 mm/h | 10 a litre | 10; 34; 50 mm | 0.7; 1.2; 1.5 mm |
| 0.3 mm/h | 100 a litre | 18; 63; 86 mm | 0.5; 1.4; 2.0 mm |
| 1 mm/h | 10 a litre | 51; 152; 197 mm | 0.9; 1.8; 2.3 mm |
| 1 mm/h | 100 a litre | 61; 179; 225 mm | 0.9; 2.5; 5.4 mm |
| 3 mm/h | 10 a litre | 141; 387; 518 mm | 0.9; 1.8; 2.2 mm |
| 3 mm/h | 100 a litre | 147; 468; 600 mm | 1.7; 9.2; 20 mm |

The snowfall is in millimetres of water an hour. On Earth 10 crystals a litre take up only 0.54 mm/h of the cloud's
vapour in the two heavier cases, the rest left unused; the arriving snow also counts the crystals' own mass entering
the top, about 0.4 mm/h at 100 a litre on Earth.

- **The breakup law does not cap the flakes.** With it and without it the sizes agree within 10–25%, and not always
  smaller with it. Phillips's law chips small fragments off the smaller particle when a meeting does not stick, which
  feeds the cloud more small crystals but splits no large flake. On both worlds the flakes' size is set by the
  snowfall's water and the time the snow takes to fall.
- **Lunar flakes grow 15–90 times the size of Earth's at the same settings:** a centimetre or two in light
  snow, 5–15 cm in moderate snow and tens of centimetres in heavy snow. The smallest crystals take 1–1.5 days to fall
  through the Moon's 16 km of cloud against about an hour on Earth, in a model whose sticking probability depends on temperature and whose encounters use differential settling alone.
- **How far to trust the sizes.** The model's Earth flakes run small against Earth's real snow, where flakes of
  several millimetres to a couple of centimetres are common, so it aggregates too little, most likely for the spread of
  fall speeds among like flakes, which it leaves out; that alone does not determine the sign of the lunar size bias. On the other side it leaves
  out large, tenuous flakes breaking in turbulence and finite cloud lifetimes,
  and Brown and Francis's mass law carried far past the sizes it was measured at. By that law a 20 cm lunar flake
  weighs about 0.9 g, a fifth of a kilogram of ice to each cubic metre it spans: a loose lace of crystals.

## Assumptions and limits

- The laws for mass, area, fall, capacitance and ventilation were fitted on Earth; they are carried to the Moon through
  the quantities they depend on (gravity, the air's density, viscosity and diffusivity), not refitted.
- Single-crystal seeds use the aggregate mass–dimension law. A habit-specific seed mass and capacitance could change both their dimensions and their subsequent mass gain.
- The columns are straight profiles and the clouds idealised: liquid saturation from top to freezing level, no
  updraughts. The population calculation includes collision chipping; large-aggregate fracture, riming, the melting flake's collapse and sublimation are left out.
- The Moon's ladder is the target study's runs interpolated at 301 K; the design years' own profile may set the
  freezing level a few hundred metres apart.

## Next

The same fall and growth physics gives the high cloud's crystals at lunar gravity for the GCM
([plan](../../research/plan.md), item 5). The [new analysis](analysis.md) measures the archived cloud environments and calculates flake-scale load scenarios. The independent [sectional reference](repair.md) now provides consistent particle accounting, but its first three full-column controls fail the finite size-domain check. Numerical convergence, a coupled moist trajectory, a spread of fall speeds among like flakes and a measured mechanical failure law remain needed to predict a survival-limited size distribution.
