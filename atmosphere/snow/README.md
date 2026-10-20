# Snow

Snow in the Open Moon's sky winter: how snowflakes fall, grow and melt at a sixth of Earth's gravity, in the air about
the freezing level, beside the raindrop solver ([../rain](../rain/README.md)). The Moon's winter weather is built up in
the sky, and on 20 October the author asked how large snowflakes can grow there and how they fall, as the second half
of the work on the high cloud's crystals ([register](../../research/decisions.md), 20 October).

| File | What it holds |
|---|---|
| [flakes.py](flakes.py) | The snowflake's mass, area and fall law, its growth from the vapour, its melting, the winter columns and the product |
| [results/snow.json](results/snow.json) | Fall speeds by size, growth across the winter cloud, crossing times and melting depths, Moon and Earth |

```sh
python -m atmosphere.snow.flakes        # seconds: results/snow.json
```

## The physics

A snowflake is an aggregate of ice crystals. Its mass follows Brown and Francis (1995), m = 0.0185 D^1.9 kg for a
maximum dimension D in metres, and its projected area is 0.35 of the circle across D, the area ratio that large
aggregates of dendrites and plates approach (0.33–0.41; Leinonen and co-authors 2026). It falls by Heymsfield and
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
small crystals to sweep up at one mixing ratio. Followed through, sweeping alone grows lunar flakes without bound, so
their size is set by what breaks them apart, which this estimate leaves out. Earth's largest flakes reach a few
centimetres; lunar flakes meet at about two-fifths of Earth's speeds, with about a sixth of the energy for the same
flakes, so larger ones are plausible. How large waits on a breakup law.

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

## Assumptions and limits

- The laws for mass, area, fall, capacitance and ventilation were fitted on Earth; they are carried to the Moon through
  the quantities they depend on (gravity, the air's density, viscosity and diffusivity), not refitted.
- Growth from the vapour takes the aggregate's mass law, which gives a single crystal too large a size for its mass,
  so the masses are the robust figures.
- The columns are straight profiles and the clouds idealised: liquid saturation from top to freezing level, no
  updraughts. Breakup, riming, the melting flake's collapse and sublimation below the freezing level are left out.
- The Moon's ladder is the target study's runs interpolated at 301 K; the design years' own profile may set the
  freezing level a few hundred metres apart.

## Next

A breakup law for how large flakes grow, from the energy of their collisions at lunar fall speeds. The same fall and
growth physics then gives the high cloud's crystals at lunar gravity for the GCM ([plan](../../research/plan.md),
item 5).
