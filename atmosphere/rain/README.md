# Rain

Single raindrops in the Open Moon's air: how fast each size falls, with its drag and its shape, how it evaporates on
the way down, and which sizes reach the ground from a cloud base, with tools for weighing drop-size spectra. The
product gives these in the Moon's air at 2 m in the three target climates beside Earth's air, and checks the fall
speeds CM1 gives its rain.

The model comes from the author's rain work of 11–13 October 2026, done in ChatGPT and recovered on 17 October (the
assessment of that day on the research drive, `integration-history-review-2026-10-17/assessment.md`, finding F02). Its
files match the original archive (`open-moon-rain-generator.zip`, SHA-256 `f44c3edf…`) byte for byte, and its solver
matches the hash its checks recorded on 11 October. It comes in cleaned up: the code reads the repository's constants
and air, its prose follows the project's standards, and its weather is rebound from the old CM1 box to the target
climates. The rain's sound, which the same work built on this physics, is rendered from the drops this domain gives
arriving at the ground by the audio lab [visualization/rain-audio](../../visualization/rain-audio/README.md).

| File | Holds |
|---|---|
| [drops.py](drops.py) | The single-drop solver: the air column, equilibrium shape, drag, terminal speed, stability screen, ventilated evaporation and the fall from cloud base to ground |
| [distributions.py](distributions.py) | Cloud-base gamma spectra, the ledger of what reaches the ground, the collision-encounter diagnostic and arrival fluxes |
| [product.py](product.py) | The product [results/rain_drops.json](results/rain_drops.json), schema `terluna.atmosphere.rain-drops/1` |
| [arrivals.py](arrivals.py) | The product [results/rain_arrivals.json](results/rain_arrivals.json), schema `terluna.atmosphere.rain-arrivals/1`: the drops arriving at the ground |
| [rates.py](rates.py) | The product [results/rain_rates.json](results/rain_rates.json), schema `terluna.atmosphere.rain-rates/1`: the GCM's rain rates on each target's land |
| [sources.json](sources.json) | The equations' and data's sources, with how far each was read |
| [vendor/rainprops/](vendor/rainprops/) | Loftus's rainprops at commit `e2df9c2` (MIT licence): the code the solver adapts, with the Gunn–Kinzer fall speeds, unmodified, each file pinned in its `PROVENANCE.json` |

```sh
python -m atmosphere.rain.product                 # about three and a half minutes on one thread
python -m atmosphere.rain.arrivals                # about two minutes
climate/gcm/.venv/bin/python -m atmosphere.rain.rates   # seconds; reads the GCM's output with netCDF4
python -m pytest atmosphere/tests/test_rain_*.py
```

## The physics

A drop is an oblate spheroid in its equilibrium shape (Green 1975), which balances surface tension against the
pressure of the air it falls through and depends on gravity through the capillary length. It falls under gravity and
buoyancy against the drag of a sphere (Clift and Gauvin 1970) raised for the spheroid's surface (Loth 2008), and
trades vapour and heat with the air by diffusion, ventilated by its fall (Pruppacher and Klett 2010). The solver
follows a drop's height, its speed over the ground, its mass and its temperature, with the air-relative speed setting
drag and ventilation. The air's viscosity, conductivity and vapour diffusivity come from kinetic theory for each gas,
mixed by Wilke's rule (Reid, Prausnitz and Sherwood 1977), so the Moon's air, with less oxygen in a denser column,
takes its own values. These are the equations of Loftus and Wordsworth (2021) as their code implements them.

The port changes three things from the recovered solver, each named in its code. The saturation pressure is the
repository's (Wagner and Pruss 2002), where the recovered copy dropped a digit from one coefficient. The surface
tension takes the IAPWS coefficient of 235.8 mN/m and critical temperature, where rainprops carries 235 mN/m, which
read 0.3% low at 20 °C. The gas constant and gravities come from [shared/constants.json](../../shared/constants.json).
Together they move fall speeds by at most 0.08% (at 6 mm on Earth) and fall times by at most 0.03%.

## Checks

| Check | Result | Evidence state |
|---|---|---|
| Earth's fall speeds: 34 Gunn–Kinzer measurements, 0.1–5.8 mm, in air at 20 °C and 50% humidity | 0.200 m/s root mean square and 7.43% at most (at 0.1 mm); the solver runs 2.6–6.4% slow from 0.6 to 2.4 mm, most at 1 mm | Empirical validation on Earth |
| The recovered solver's saved speeds, drops and drips, on its own inputs | Within 0.08% in speed and 0.03% in fall time | Reproduction of the original |
| Stokes limit: drops of 2–10 µm | Within 0.5% of (ρw − ρa) g D² / 18μ, falling as g and as D² | Physical limit |
| Inertial regime: rigid 6 mm spheres | Slowed by lunar gravity 2.47–2.49 times, against √(9.81/1.62) = 2.46 | Physical limit |
| The stability screen | Grows 2.46 times at lunar gravity, as g^(−1/2) (Loftus and Wordsworth 2021) | Physical scaling |
| Released at rest over 0.1 m and 1 m | Below the vacuum free-fall speed and the terminal speed | Physical bound |
| Saturated isothermal air over 0.1 m; the water lost; tighter integration | No evaporation or cooling; the water a drop loses equals its evaporation rate integrated along its fall within 10⁻⁶; status unchanged and evaporated share within 10⁻⁴ | Numerical checks |

The lunar values are this physics carried to lunar gravity and air, and the evaporation rests on its laboratory
ventilation law: both are runnable calculations.

## Results

The air at 2 m is each target's median over land in the
[climate targets](../../research/studies/climate_targets/README.md), whose 2 m state is the GCM's lowest level,
0.88–0.91 km up, carried down to the ground: their 2 m temperature is that level's on the dry adiabat, 1.42–1.47 K
warmer in the medians, and their extractor computes the 2 m wet-bulb temperature from that level's specific humidity,
capped at saturation. The GCM's output holds that humidity as one step at each 3-day window's boundary, which
`domain/climate` is correcting. Each target here pairs the two medians, and the humidity inverts the extractor's
wet-bulb relation at 1.2 atm, where the extractor used each cell's own surface pressure; that moves the vapour pressure
by 0.6%. Started from the lowest level itself, the parcel's cloud base comes 21–30 m lower (`column_sensitivity`).
The Moon's air is the design air at 1.2 atm, with Earth's oxygen partial pressure (21,227 Pa) and 400 ppm of CO2
(`earthlike_air` in [thermodynamics.py](../radiative_convective/thermodynamics.py)). Earth's air is at 1 atm, with the
same temperature and relative humidity.

| | Hot world | Resort | Lock |
|---|---|---|---|
| Air at 2 m over land: median temperature; relative humidity | 33.5 °C; 58% | 26.1 °C; 62% | 20.4 °C; 67% |
| Air density: Moon; Earth | 1.362; 1.138 kg/m³ | 1.400; 1.170 kg/m³ | 1.429; 1.195 kg/m³ |

### Fall speed

In still air at 2 m, in m/s:

| Diameter | Moon, hot world | Moon, resort | Moon, lock | Earth, the resort's air |
|---|---|---|---|---|
| 0.2 mm | 0.155 | 0.156 | 0.157 | 0.710 |
| 0.5 mm | 0.556 | 0.556 | 0.556 | 2.03 |
| 1 mm | 1.158 | 1.152 | 1.149 | 3.79 |
| 2 mm | 2.14 | 2.12 | 2.11 | 6.31 |
| 3 mm | 2.92 | 2.90 | 2.88 | 7.83 |
| 6 mm | 4.43 | 4.39 | 4.36 | 9.24 |
| 10 mm | 5.21 | 5.16 | 5.11 | above the screen |
| Stability screen | 16.2 mm | 16.3 mm | 16.4 mm | 6.6 mm |

Drops on the Moon also keep rounder: a 6 mm drop is 0.88 as tall as it is wide, against 0.65 on Earth. The hot
world's warmer air is 5% less dense than the lock's, so its large drops fall 2% faster. At 1.16 atm, the GCM's global
mean surface pressure and close to the 1.155 atm its air implies over land, the resort's drops fall 0.2% (0.1 mm) to
1.6% (6 mm) faster than at 1.2 atm (`pressure_sensitivity`).

### Survival below the cloud

Each drop leaves the parcel cloud base at its terminal speed and falls through still air. The column below the cloud is
hydrostatic, cools at the dry adiabat g/cp and keeps the 2 m vapour mole fraction, so the cloud base is where it
saturates. Under lunar gravity that air cools 1.6 K per km as it rises, six times more slowly than on Earth, so air of
the same temperature and humidity saturates six times higher. The column stands for the daytime convective mixed layer
under a raining cell, which CM1's box finds 8–10 km deep at 65–75% humidity
([climate/crm](../../climate/crm/README.md)); the GCM's mean air over land, through day and night, cools only
0.76–0.92 K per km between 0.9 and 5 km and stands 2.8–4.0 K warmer than the column at 5 km (`column_sensitivity`). The
rain neither cools nor moistens the column: rain falling through air its own evaporation has moistened would lose less,
and a cloud base below the parcel height would shorten the fall. Gravity holds its surface value, 0.85% above its value
at 7.4 km. The dry adiabat takes the gas-transport model's heat capacity, which carries argon, 1% of the air, as
nitrogen and so stands 0.8% above the repository's (`thermodynamics.py`); with the repository's, the cloud bases sit
49–78 m (1.0%) lower and the smallest surviving drop 0.3–0.6% smaller (`column_sensitivity`).

| | Hot world | Resort | Lock |
|---|---|---|---|
| Cloud base: Moon; Earth | 7.42 km; 1.23 km | 6.09 km; 1.01 km | 4.96 km; 0.82 km |
| Smallest drop reaching the ground: Moon; Earth | 2.42 mm; 0.75 mm | 1.93 mm; 0.61 mm | 1.53 mm; 0.50 mm |
| Water reaching the ground, mean drop 0.75 mm: Moon; Earth | 0%; 20–26% | 0–0.3%; 36–40% | 0.3–1.7%; 53–55% |
| The same, mean drop 1.5 mm | 1.7–4.6%; 65–66% | 6.9–12%; 76–77% | 18–24%; 84–86% |
| The same, mean drop 2.5 mm | 19–22%; 83–86% | 35–37%; 89–91% | 52%; 93–94% |
| Encounters on the way down, central spectrum: Moon; Earth | 36; 1.3 | 26; 1.1 | 22; 0.9 |

A 1 mm drop leaves the cloud base at 1.2 m/s and evaporates after 56–62 minutes, 2.8–3.1 km below it. The spectra are
the rain work's bounded set: gamma concentrations at the cloud base with mass-weighted mean diameters of 0.75, 1.5 and
2.5 mm, each range spanning shape parameters 0 and 3, on 64 classes from 0.1 to 6 mm; on twice as many classes the
resort's shares move by at most 0.014 percentage points (`grid_check`). With independent drops the share reaching the
ground is the same at any rain rate. The encounters are the expected number of drops each surviving
drop's path meets (geometric kernel, efficiency 1), weighted by the water reaching the ground, for the central
spectrum (1.5 mm, shape 3) at 5 mm/h; they grow in proportion to the rain rate. On the Moon all the water reaching the
ground travels in drops that meet more than one other, so collisions and coalescence would reshape the spectrum, and
these shares are a filter on independent drops. CM1's box loses 85% of its falling water to evaporation
([climate/crm](../../climate/crm/README.md)), within the range above.

### The drops arriving at the ground

`rain_arrivals.json` keeps the ground's side of the same falls, class by class, for every bounded spectrum: each
surviving class's size and speed at impact and its drops per square metre per second for each mm/h of rain on the
ground. For the central spectrum (1.5 mm, shape 3), at 1 mm/h on the ground:

| | Hot world | Resort | Lock |
|---|---|---|---|
| Drops arriving per m² per second: Moon; Earth | 151; 294 | 213; 308 | 203; 304 |
| Their median size and speed: Moon | 1.23 mm, 1.40 m/s | 1.03 mm, 1.18 m/s | 1.09 mm, 1.25 m/s |
| The same: Earth | 0.93 mm, 3.61 m/s | 0.91 mm, 3.52 m/s | 0.93 mm, 3.54 m/s |
| Rain leaving the cloud base per mm/h on the ground: Moon; Earth | 58; 1.5 mm/h | 14; 1.3 mm/h | 5.5; 1.2 mm/h |

The Moon's drops arrive shrunk to sizes near Earth's, at about a third of the speed, and fewer of them carry the same
water. Its rain rate at the ground asks a heavy rain of the cloud base in this column, where coalescence and the rain's
own moistening of the air would let more water through. A 3 mm drop released at rest 1 m up, as from a leaf, strikes at
1.59–1.60 m/s on the Moon and 4.01 m/s on Earth; drops of 8–12 mm leaving the cloud base reach the Moon's ground with
85–97% of their water at 4.8–5.4 m/s, and exceed Earth's stability screen. The airs are the rain product's, so their
humidity inherits the climate targets' pairing of the GCM's 2 m temperature with its lowest level's humidity, which the
GCM's output holds as a snapshot at each 3-day window's boundary (`domain/climate` is correcting it).

### Rain rates on the land

`rain_rates.json` reads the GCM's rain, a true mean over each output window, on the land of the two runs behind the
climate targets, each cell weighed by its land area, and sets each target between them as the targets do. A mean is
wet at 0.1 mm/h or more. The design case's two years written every 6.5 hours (`A28_corrected_970_eddies`) give the
weather within the 3-day means, and each target's 6.5-hour figures take the design case's ratios of the two:

| Land rain | Hot world | Resort | Lock |
|---|---|---|---|
| Mean | 3.34 mm/day | 2.72 mm/day | 2.25 mm/day |
| 3-day means: wet share of the time; wet median; 99th percentile | 24%; 0.39; 2.48 mm/h | 21%; 0.35; 2.36 mm/h | 18%; 0.31; 2.28 mm/h |
| 6.5-hour means: the same | 12.5%; 0.71; 5.50 mm/h | 10.9%; 0.62; 5.25 mm/h | 9.6%; 0.56; 5.05 mm/h |

The 3-day means smooth the rain into spells nearly twice as common and about half as hard: in the design case the land
is wet in 18.4% of its 3-day means and 9.7% of its 6.5-hour means, and the 6.5-hour wet median is 1.8 times the 3-day
one, its 99th percentile 2.2 times and its largest 2.3 times (14.3 against 6.3 mm/h over the same two years). The
6.5-hour years averaged eleven at a time match the design case's own 3-day output over its years 10–19 within 0.8% at
the wet share and every quantile to the 99th. Every figure is a mean over a T21 cell, about 170 km across: rain at a
point under a shower falls harder and for less time.

### The check on CM1's fall speeds

CM1's Morrison and NSSL schemes slow every raindrop at lunar gravity by one factor, (g / 9.81)^((B + 1)/3) with
B = 0.8, which is 2.94 ([climate/crm](../../climate/crm/README.md)). The solver, with gravity alone changed in the
Moon's air at 2 m, slows each size by its own factor. At the resort's air:

| Diameter | Slowing by gravity: solver | Rigid sphere | CM1 | CM1's lunar law over the solver | Its Earth law over the solver on Earth | Water lost per metre, CM1 over the solver |
|---|---|---|---|---|---|---|
| 0.1 mm | 5.26 | 5.26 | 2.94 | 3.39 | 2.05 | 0.32 |
| 0.2 mm | 4.37 | 4.37 | 2.94 | 1.75 | 1.25 | 0.64 |
| 0.5 mm | 3.44 | 3.46 | 2.94 | 1.03 | 0.91 | 0.98 |
| 1 mm | 3.08 | 3.13 | 2.94 | 0.86 | 0.85 | 1.10 |
| 2 mm | 2.76 | 2.92 | 2.94 | 0.81 | 0.89 | 1.13 |
| 3 mm | 2.49 | 2.78 | 2.94 | 0.82 | 0.99 | 1.11 |
| 6 mm | 1.93 | 2.48 | 2.94 | 0.95 | 1.46 | 1.03 |

The hot world and the lock differ from these by under 3% (`cm1_check`). At the review's air, 20 °C and 80% humidity,
the solver gives 4.35, 3.43, 3.07, 2.76, 2.49 and 1.93 at 0.2, 0.5, 1, 2, 3 and 6 mm. The assessment's rerun of the
recovered solver gave the same to 3 mm and 1.92 at 6 mm (1.9243), which the port's surface tension moves to 1.9252.

The two agree at 1.37–1.39 mm. Small drops fall in the viscous regime, where speed follows gravity itself (Stokes
drops fall 6.04 times slower); large ones fall in the inertial regime, where a rigid sphere's speed follows the square
root of gravity, and under lunar gravity a large drop also flattens less, which cuts its drag further. One exponent
cannot follow both. Against the Earth law it starts from, CM1's lunar drizzle falls about 1.5 times too fast and its
6 mm drops about 1.5 times too slowly.

CM1's whole law also carries its Earth fit. Morrison's V = 842 D^0.8 with its density factor gives Earth's 0.2 mm drops
23–27% more speed than the solver, its 1–2 mm drops 11–16% less and its 6 mm drops 45–47% more. In the Moon's air at
2 m, CM1's drizzle falls 1.72–1.79 times the solver's speed at 0.2 mm and 3.3–3.5 times at 0.1 mm, the two meet at
0.52–0.56 mm, CM1's 1–4 mm drops fall 13–19% slower than the solver's, and its 6 mm drops 5% slower. CM1 also caps the
rain's mean fall speed at 2.66–2.74 m/s near the ground, the speed of a 2.7 mm drop in the solver, where the largest
lunar drops reach 4.4–5.2 m/s.

For CM1's rain, the time spent below the cloud sets the evaporation. With the same ventilation law, a 0.2 mm drop at
CM1's speed loses 35–38% less water per metre of fall than at the solver's, and a 0.1 mm drop 68–69% less, so CM1's
drizzle survives longer falls than it would; its 1–3 mm drops lose 9–13% more per metre than they would, so its main
rain evaporates somewhat too fast. How the two offset in CM1's 85% loss depends on how its rain divides between
drizzle and millimetre drops; a CM1 run with a size-dependent law would measure it. CM1 was not rerun.

## Assumptions and limits

- **The drops.** Pure liquid water in the equilibrium shape, held during acceleration, with near-ambient liquid
  properties; the warm-water domain runs from 0 to 40 °C at the ground with the cloud base above freezing. Argon, 1% of
  the air, takes nitrogen's transport properties. The ventilation law was measured below a Reynolds number of 2,600,
  which lunar drops stay under to about 7 mm and Earth's to about 4.5 mm.
- **Independent drops.** Each drop falls alone; the encounter diagnostic shows that on the Moon collisions,
  coalescence and breakup would reshape the spectrum. The stability screen is a size scale for breakup.
- **The canopy and the impacts.** The recovered work's impact descriptors and canopy water ledger
  (`physical_event_interface.py`) wait for `biosphere/canopy`, whose interception, storage and drip release they need;
  until then the rain's sound takes an assumed share of the water as drips (rain_arrivals.json gives their falls).
