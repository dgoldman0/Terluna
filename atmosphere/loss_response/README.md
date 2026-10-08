# Loss response to protection

How fast the finished 1.2 atm atmosphere loses mass as a function of its UV
transmission, the share of sunlight below 175 nm that the optical shield lets
through, and of how much of the solar wind reaches the air. It sizes the two
protection functions against a total loss budget for the protection design study
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).
[tides.py](tides.py) finds how much Earth's tide raises escape.
[traced.py](traced.py) reads the middle atmosphere's limb tracing and gives the
state the upper air settles to under the light the shield lets through.
[absorption.py](absorption.py) finds how far out the optical shield must reach
and what the sunlit exosphere beyond it produces. [exosphere.py](exosphere.py)
finds what that exosphere loses, as ions and as the fragments of broken
molecules, with and without a lunar magnetosphere, and the UV transmission each
budget then allows. [fate.py](fate.py) follows the escaping air into the space
around Earth and the Moon. [cycle.py](cycle.py) reads FISM2's daily record band by
band from 1947 and gives each calendar year's measured sunlight, for the
solar-cycle mean that R1 reads, and the year around each solar maximum since
cycle 19, the strongest of which is the loss response's solar maximum.
[oxygen.py](oxygen.py) follows the oxygen atoms and hydrogen the glow and the
traced light make above the base through the whole column, from the ground to the
exobase, and finds what of them escapes. [gravity_waves.py](gravity_waves.py)
estimates how strongly breaking gravity waves mix the upper air, on Earth and on
the Moon, which decides how much of that oxygen reaches the exobase. All nine are
screening models; their methods and limits are in their docstrings.

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # about 4 minutes; resumes if stopped
python -m atmosphere.loss_response.cycle        # FISM2's daily record by band, a few seconds
OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.limb_heat   # then the limb tracing, about 30 CPU minutes (cached)
python -m atmosphere.loss_response.model        # then the loss response, about a minute
python -m atmosphere.loss_response.absorption   # then the protected radius, a few minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.exosphere   # then the exosphere's losses, about 20 minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.fate    # where the escaping air goes, about 15 s
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.oxygen  # atomic oxygen and hydrogen, about 40 minutes in three processes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.gravity_waves   # the waves' mixing, a minute
python -m pytest atmosphere/tests/test_tides.py atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py atmosphere/tests/test_exosphere.py atmosphere/tests/test_fate.py atmosphere/tests/test_cycle.py atmosphere/tests/test_oxygen.py atmosphere/tests/test_gravity_waves.py
```

It needs the middle atmosphere's stored results, its limb tracing
([results/limb_heat.json](../middle_atmosphere/results/limb_heat.json)), the WHI
2008 spectrum (`python -m atmosphere.middle_atmosphere.fetch_inputs --download`)
and the FISM2 spectra (`python -m atmosphere.middle_atmosphere.fetch_limb_inputs
--download`). The product is
[results/loss_response.json](results/loss_response.json) (schema
`terluna.atmosphere.loss-response/3`), with the full sweep in
[results/euv_sweep.csv](results/euv_sweep.csv). Its field
`allowed_leak_fraction` is the allowed UV transmission from ultraviolet-driven
escape and the solar-wind scalings alone. The exosphere's losses, and the allowed
transmissions with them included, are in
[results/exosphere_loss.json](results/exosphere_loss.json) (schema
`terluna.atmosphere.exosphere-loss/2`). The tidal multipliers are
in [results/tidal_escape.json](results/tidal_escape.json) (schema
`terluna.atmosphere.tidal-escape/1`), and the protected radius and the
exosphere's ion production in
[results/absorption_radius.json](results/absorption_radius.json) (schema
`terluna.atmosphere.absorption-radius/2`). The cloud the escaping air forms is in
[results/escape_fate.json](results/escape_fate.json) (schema
`terluna.atmosphere.escape-fate/1`), and each year's measured sunlight in
[results/solar_cycle.json](results/solar_cycle.json) (schema
`terluna.atmosphere.solar-cycle/2`). The oxygen atoms and hydrogen are in
[results/oxygen_escape.json](results/oxygen_escape.json) (schema
`terluna.atmosphere.oxygen-escape/1`), and the waves' mixing in
[results/gravity_waves.json](results/gravity_waves.json) (schema
`terluna.atmosphere.gravity-wave-mixing/1`); the oxygen step needs the radiative inputs
(`python -m atmosphere.radiative_convective.fetch_inputs --download`) for the
middle atmosphere's photolysis.

## How it works

The ultraviolet branch takes the heat that reaches the upper air from the
middle atmosphere's limb tracing
([limb_heat](../middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07),
read through [traced.py](traced.py)). It counts the light below 175 nm that
passes the shield's window stack and the ring fleet's 4 µm annulus film, a UV
transmission f through gaps over the whole aperture of the protected radius,
and unfiltered sunlight beyond the aperture, each along its slant paths through
air swollen by the heat it takes, and the sky's interplanetary Lyman-alpha glow
(about 1,000 rayleigh, 4.1×10⁻⁶ W/m² onto a surface facing open sky), which
reaches the upper air from every direction. The state is the first heat,
counting up from zero, that the swollen air returns; where none lies in the
tables, the air runs away. The author adopted this count on 7 October
([decisions](../../research/decisions.md)) in place of a quarter of the light
reaching the disk. Solar maximum is FISM2's year around the strongest maximum
since 1978, cycle 21's of 1979–80, every band at its own measured rise and the
glow at Lyman-α's ([cycle.py](cycle.py)). Two stress cases sit beside it: cycle
19's year of 1957–58, the strongest on record, whose spectrum FISM2 builds from
the 10.7 cm radio flux alone, and the escape model's 2.5 on the ultraviolet with
FISM2's mean rise of the X-rays over the last three maxima
([escape.py](../middle_atmosphere/escape.py)). The films are the design's behind
both shields, and the protected radius is the ring fleet's 4 lunar radii unless
stated. The thermal column turns the heat into an exobase and a
molecular Jeans loss, over base temperatures from the middle atmosphere behind
the titania stack and behind the 200-nm edge, for each treatment of the upper
air. It takes each state's heat where the tracing puts it in height: each part
of the light by the limb tables' record of where it heats, mixed by band at the
state's activity, and the glow where O2 absorbs it just above the base
([traced.py](traced.py), `model.state_config`), as the author approved on
7 October in place of the column's fixed 'middle' shape. Beyond the column's domain the energy-limited expression of the September
feasibility report bounds the loss. The solar-wind branch combines pickup,
capped by mass loading at the solar wind's own mass flux through a disc three
lunar radii across (0.29 kg/s), with sputtering by returning pickup ions and by
precipitating protons, over stated parameter ranges.

Earth's tide lowers the barrier a molecule must climb to leave the Moon: in the
frame turning with the Moon's orbit, a molecule can slip out of the Moon's Hill
sphere through the necks around the L1 and L2 points, 58,000 km from the Moon,
while slower than the Moon's own escape speed. The tidal step launches test
molecules from the exobase in the Earth–Moon three-body problem and counts how
many of those with energy enough to pass L1 actually leave; most that do take
about five days. Each species' Jeans loss is multiplied by the result for its
exobase radius and temperature.

A protected radius covers the thermosphere's heating where unfiltered sunlight
beyond its aperture, crossing the collisional air below the exobase at a grazing
angle, gives at most a tenth of the heat of the state the air settles to with
that radius; the smallest such radius is the heating radius. The absorption step
keeps its band-mean limb deposit as a cross-check. Above the exobase the air is
collisionless. Chamberlain's
exosphere, built from the exobase density, temperature and radius, gives the
molecules on ballistic flights, which rise and fall back, and the escaping
ones, which the Jeans loss already counts. Sunlight ionizes a molecule in about
19 days at quiet Sun and 8 days at solar maximum, and the step counts the ions
made from ballistic molecules in sunlight outside the shadow.

The exosphere step takes the same exosphere for N2 and O2 separately, with the
O2 share at the exobase either well mixed, as the column has it (17.5%), or
separated by diffusion above a homopause at 0.01 Pa (3–4%). Outside the shadow
sunlight breaks molecules apart faster than it ionizes them: N2 in 10–26 days
and O2 in 4–5 days, on the WHI spectrum and within 10% of Heays et al. (2017).
The atoms fly apart with 0.1–1.7 eV each, above the 0.1–0.2 eV an atom needs
to leave from two to five lunar radii, and escape unless their path meets the
exobase. The solar wind's protons also take electrons from exospheric
molecules and turn them into ions, inside the optical shadow as well, since the
wind closes in behind an absorbing screen over about eight screen radii. The
ring fleet's screen is about 20,000 km out, under three screen radii from the
Moon, so the Moon may sit in its unrefilled wake; a sensitivity takes the wind
as absent on paths within the wake's core, 2.6 lunar radii behind a 4-radius
screen. With no
magnetosphere the wind carries the ions off, and those that curve back into the
air sputter out more molecules. A lunar dipole holds the wind off at its
stand-off, 10 lunar radii for the September design's magnets
(1.5×10²¹ A·m²). Inside it, ions on open field lines leave down the tail, and
ions on closed field lines drain along the field into the air and stay, or are
carried to the boundary by the convection the wind drives, or recombine into
fast atoms that mostly escape, in proportion to the three rates.

## Results

Rederived on 7 October 2026 with the traced count, FISM2's X-ray cycle and the
ring fleet's protected radius of 4 lunar radii, the column taking each state's
heat where the tracing puts it in height. The design point adds the mean over
the solar cycle measured band by band.

**Earth's tide roughly doubles escape.** Across the loss response's sweep the
tide multiplies molecular escape by 1.8–4.9; at the allowed transmissions it is
1.8–2.7. A uniform lowering of the barrier, the feasibility baseline's factor,
gives more there, because it counts every molecule above the lowered barrier
as lost, while many of them are aimed away from the necks and fall back.

**Where the heat lands matters as much as how much.** The light through gaps
heats the thermosphere at three depths: Lyman-α within a few e-folds of
pressure above the base, the far ultraviolet in the lower thermosphere and the
extreme ultraviolet in its upper half
([the limb tracing](../middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07)).
Heat laid low is conducted down to the base with little warming of the air
above. Placed where the tracing puts it, the same heat leaves the exobase
26–51 K cooler than the column's 'middle' shape gives, and the
ultraviolet-driven loss 16–3,600 times smaller, the factor largest for the
coolest air.

**The sky's glow heats the air just above the base.** O2 takes the glow's
Lyman-α within a few e-folds of pressure above the 0.3 Pa base, where its heat
is conducted away: with the UV transmission at 10⁻⁶ it raises the exobase by
6–9 K at quiet Sun, 10–14 K at solar maximum and in the escape model's stress
case, and 12–18 K in cycle 19's year. Behind the titania stack, with the films'
traced X-rays, the floor loss is 5×10⁻¹⁰–0.002 kg/s at quiet Sun and
7×10⁻⁹–0.02 kg/s at solar maximum (3×10⁻⁸–0.04 in cycle 19's year and
5×10⁻⁹–0.01 in the escape model's). Behind the 200-nm edge, whose upper air
starts warmer, it is 0.01–0.3 kg/s at quiet Sun and 0.1 kg/s at solar maximum
(0.3 and 0.09 kg/s), and with all near-infrared heating at any maximum the air
runs away with no gaps at all.

**At solar maximum the ultraviolet rises far less than the escape model's 2.5,
the glow and the shortest light more.** FISM2's daily record from 1947
([cycle.py](cycle.py)) puts the light from 10 to 175 nm, weighted by energy, at
1.25–1.43 times the WHI 2008 quiet week in the year around each maximum of
cycles 20 to 25, at 1.64 around cycle 19's of 1958, the strongest on record, and
at 1.14 over solar cycles 23 and 24; the far ultraviolet from 122 to 175 nm, most
of the energy, rises 1.07–1.59 at maximum. Lyman-α, which sets the sky's glow,
rises 1.31–1.88 at maximum against the escape model's 1.5, the X-rays below 10 nm
4.6–9.6 against its 6.0, and the light from 10 to 30 nm up to 3.8 against its
2.5. The tables here take cycle 21's year, the strongest since 1978, when FISM2
adds the MgII and Lyman-α proxies to the 10.7 cm radio flux, as solar maximum,
whole; cycle 19's year, whose spectrum rests on the radio flux alone, and the
escape model's factors are stress cases. The
protection design point gives each state's mean over the measured years of
cycles 23 to 25, which R1, a long-term average, reads.

**Ultraviolet-driven escape alone lets the gaps pass about a tenth of a
percent.** Largest allowed UV transmission from ultraviolet-driven escape
and the solar-wind scalings, across the three upper-air treatments, with the
solar-wind loss at the low end of its range. Each budget's cycle time is the
atmosphere's mass (3.1×10¹⁸ kg, the feasibility baseline's hydrostatic column)
over the loss rate, the time the loss and its resupply take to replace the
whole atmosphere:

| Shield and Sun | 1 kg/s (100 billion years) | 10 kg/s (10 billion years) | 100 kg/s (1 billion years) |
|---|---|---|---|
| Titania, quiet Sun | 0.058–0.18% | 0.081–0.22% | 0.081–0.26% |
| Titania, solar maximum (cycle 21) | 0.025–0.11% | 0.032–0.12% | 0.032–0.14% |
| Titania, cycle 19's year (stress) | 0.018–0.086% | 0.018–0.10% | 0.018–0.11% |
| Titania, the escape model's 2.5 (stress) | 0.016–0.062% | 0.021–0.071% | 0.021–0.079% |
| 200-nm edge, quiet Sun | 0.0092–0.035% | 0.014–0.049% | 0.014–0.049% |
| 200-nm edge, solar maximum (cycle 21) | 0.011–0.012% (none for the warmest upper air) | 0.012% (none for the warmest) | 0.012% (none for the warmest) |
| 200-nm edge, cycle 19's year (stress) | 0.0012–0.0013% (none for the warmest) | 0.0012–0.0013% (none for the warmest) | 0.0012–0.0013% (none for the warmest) |
| 200-nm edge, the escape model's 2.5 (stress) | 0.0078–0.0081% (none for the warmest) | 0.0087–0.0088% (none for the warmest) | 0.0087–0.0088% (none for the warmest) |

From 10 to 100 kg/s the allowed transmission rises little, and for the warmest
upper air not at all: the runaway caps it before the loss does. The measured
maximum allows the titania stack 1.5–1.8 times what the escape model's 2.5 does
and the 200-nm edge 1.4 times, its far ultraviolet being weaker. Cycle 19's year
allows the cooler titania cases a quarter to two-fifths more than the 2.5 and
the 200-nm edge a sixth as much, its glow and shortest light being stronger.

**With a 4-radius shadow the air runs away near the shadow's edge.** As the heat
swells the air, its exobase approaches the aperture's edge and the unfiltered
light beyond it heats the thermosphere more, until the heat outgrows itself.
Behind the titania stack that begins at a transmission of 0.14–0.26% with
collisional upper air (solar maximum to quiet Sun; 0.11% in cycle 19's year and
0.079% in the escape model's), 0.10–0.19% in LTE (0.076% and 0.060%) and
0.032–0.081% with all near-infrared heating (0.018% and 0.021%), and behind the
200-nm edge at 0.012–0.049% (0.0013% and 0.0087%). The loss just short of it is
71–85 kg/s, 36–47 kg/s and 2–7 kg/s behind the titania stack and 1–4 kg/s
behind the edge, down to 57, 26, 1 and 0.4 kg/s in cycle 19's year. These onsets come from
tables traced through air the column heats with its middle shape, which swells
more than the traced shape's: at each onset the traced column's own exobase
still lies at 3.2–3.9 lunar radii, so the onsets lean early.

**The solar-wind branch as screened.** Its range is 0.03 kg/s (low),
0.19 kg/s (central) and 2.1 kg/s (high); at the high end sputtering by
returning pickup ions dominates. With pickup capped at 0.29 kg/s it barely
changes the allowed transmission for budgets of 10 kg/s and more, and at 1 kg/s
its high end alone exceeds the budget. The ion production below exceeds that cap
many times over, and the exosphere step finds that those ions are carried off;
its losses replace the capped pickup and the ion sputtering, and the table above
holds for the ultraviolet-driven part alone.

**For the thermosphere's heating the shield must reach just beyond the
exobase.** With its edge inside the exobase, the unfiltered light beyond the
aperture heats the thermosphere along slant paths and the air swells into it
until it runs away. At the allowed transmissions the heating needs 3.5–4 lunar
radii for 1 and 10 kg/s and 4–4.5 for 100 kg/s behind the titania stack, and 4
behind the 200-nm edge: the first tabulated radii beyond the exobase of the
swollen air the tables trace. The ring fleet's 4 lunar radii cover the heating up
to 10 kg/s and in all but one case at 100 kg/s. The extreme ultraviolet itself
is absorbed lower, with unit optical depth at grazing incidence at 2.0–2.5 lunar
radii, and the thin air above that still takes up too much.

**The sunlit exosphere beyond the shadow makes ions faster than the budget.**
At the allowed transmissions the exobase sits at 2.6–3.9 lunar radii, where
gravity is weak. Outside a shadow sized for heating, its ballistic molecules
become ions at 2–5 kg/s for the 1 kg/s budget, 8–20 kg/s for 10 kg/s and
9–55 kg/s for 100 kg/s, counted to a third of the Hill radius; counting to the
whole Hill radius adds up to three fifths. A shadow that held production to a
tenth of the budget would need 4–10 lunar radii counted to a third of the Hill
radius, 4–13 to half of it, or 4–24 to the whole. The exosphere step follows
these ions.

**At those transmissions the sunlit exosphere overruns every budget.** Where
ultraviolet-driven escape already fills the budget, the sunlit exosphere
outside the ring fleet's 4-radius shadow adds a central 5–174 kg/s with no
magnetosphere (cycle times of 20 to 0.6 billion years on its own) and
1.7–108 kg/s with the September magnets (58 to 0.9 billion years). The solar
wind's charge exchange alone makes up to 5–21 kg/s of ions there. Most of it
comes from the first few scale heights above the exobase.

**Without a magnetosphere the ions leave.** In the flowing solar wind the
motional electric field pulls an ion about 40,000 times harder than lunar
gravity at three lunar radii, and gravity holds an ion only where the flow
stops, at an obstacle's surface. The exosphere's ions are carried off, as at
comets, whose ion production likewise far exceeds the wind's mass flux through
the region that makes it: the loaded wind slows and its interaction region
grows to tens of lunar radii. A tenth to a half curve back into the air and
sputter out 1–10 molecules each. The solar wind's charge exchange with the dense
exosphere near the exobase makes ions inside the optical shadow too; in the
central estimate the design point's exosphere then loses 0.6–1.9 kg/s in five of
its six cases, which alone uses a 1 kg/s budget.

**The ring fleet's wake may shelter the inner exosphere.** The ring fleet's
screen, about 20,000 km out, lies under three screen radii from the Moon, and
the wind closes in behind an absorbing screen over about eight. If the wind is
absent on paths within the wake's core at the Moon (2.6 lunar radii behind a
4-radius screen), the charge exchange of the cool upper air nearly vanishes: at
the design point's standard level the three coolest cases lose 0.064–0.074 kg/s
without a magnetosphere, against 0.65–0.90 kg/s with the wind everywhere. Where
the exobase lies outside the core the wake changes less. The refill length is an
assumption a plasma model has to test.

**A lunar magnetosphere holds the wind off and keeps part of the ions.** With
the September magnets, of the ions made outside the 4-radius shadow at the
allowed transmissions, 28–42% drain along closed field lines into the air and
stay, 32–33% leave along the open field lines over the poles, and
recombination into fast atoms and the wind-driven convection take the rest; the
warmer the upper air, the farther its exosphere reaches and the less the
magnetosphere keeps. The magnetosphere also keeps the solar wind, and its
charge exchange, off the exosphere inside the stand-off. For most upper air a
much weaker dipole does that: at the design point's standard level the three
coolest cases lose under 1 kg/s from their exosphere without one, LTE at solar
maximum and all near-infrared heating at quiet Sun need 2–3×10¹⁹ A·m² (a
stand-off near 2.7 lunar radii, a fiftieth of the September moment), and no
moment up to 10²³ A·m² holds the warmest at solar maximum to 1 kg/s; in cycle
19's year it finds no state within the 4-radius shadow at that level, which no
magnet changes. Dipoles that hold the wind off only just above the exobase fall in the
range where hybrid simulations find a weak field can raise ion escape (Egan et
al. 2019), which this step leaves out.

**Magnets do not hold the fragments.** Sunlight breaks molecules apart in the
sunlit exosphere faster than it ionizes them, and the atoms are neutral. The
escaping fragments are 69–81% as much as the ions made, and with a
magnetosphere they are most of what the sunlit exosphere loses. Only the shadow
removes them.

**Atomic oxygen made above the base matters only for warm upper air with weak
mixing.** Behind the titania stack nothing splits O2 below the base, but above it
the glow makes 7–14 kg/s of oxygen atoms and the light through gaps 3–11 kg/s
more at the standard level, from quiet Sun to the solar maxima
([oxygen.py](oxygen.py)). With the middle atmosphere's eddy mixing, scaled for
the Moon, almost all of it goes back down as odd oxygen and builds 1.6–2.9 DU of
ozone in the middle atmosphere. The homopause then lies at 3–6×10⁻⁵ Pa and the
atoms are 0.1–1.3% of the gas at the exobase. They lose 0.002–0.006 kg/s with
collisional upper air and in LTE at the solar cycle's mean spectrum, and
0.12 kg/s averaged over the cycle's years with all near-infrared heating
(0.06 kg/s with the September magnets), 0.34 kg/s in its worst year. Most of it
leaves as ions the solar wind makes from the sunlit exosphere's atoms. Atoms made
within a collision length or so of the exobase leave hot, since photolysis and
recombination give them more than the escape energy there; they add
0.001–0.002 kg/s over the cycle (at most 0.011 kg/s under the escape model's 2.5, letting
them through two collisions), half the loss of the coolest case. Earth's own
mixing, unscaled, puts the homopause at 1–3×10⁻³ Pa: the atoms are then 2–11%
of the exobase gas and lose 0.035–0.12 kg/s in the cooler cases and 2.1 kg/s
with all near-infrared heating (1.0 kg/s with the magnets). Behind the 200-nm
edge, whose middle atmosphere already holds atomic oxygen under about 300 DU of
ozone, the atoms add 2–4 kg/s. Hydrogen from water and H2 escapes at
0.001–0.003 kg/s.

**Breaking gravity waves mix the upper air about as strongly as the Moon-scaled
profile.** How strongly the air above the base mixes decides the atoms
([gravity_waves.py](gravity_waves.py)). Gravity waves launched from half the
tropopause pressure with Earth's 4 mPa, the calculation calibrated so that it
gives Earth's measured eddy diffusion of 1.1×10⁶ cm²/s at 96 km (Swenson et al.
2021) and its homopause at 110 km, put the Moon's homopause for atomic oxygen at
1–2×10⁻⁵ Pa, a little above the Moon-scaled profile's. The lunar thermosphere's
buoyancy frequency falls to about 10⁻³ per second, so short, fast waves turn back
there, while waves a few hundred kilometres long, which a troposphere six times
deeper than Earth's favours, rise almost undamped and break near the exobase.
With the waves' mixing above the base the warmest case's atoms lose 0.10 kg/s
over the cycle (0.05 kg/s with the magnets; a trillion years) and the cooler
cases 0.002–0.006 kg/s; with Earth's wavelengths and a third of the calibration,
the low end of Earth's measured mixing, the warmest case's atoms lose 0.45 kg/s
(0.22). Titan, the nearest analogue, a 1.5-bar nitrogen atmosphere at
1.35 m/s², mixes less: its measured homopause eddy diffusion, 2×10⁷–10⁸ cm²/s
(Yelle et al. 2008; Bell et al. 2014), carried to the Moon is about a tenth of
the Moon-scaled mixing, and with it the atoms lose 0.77 kg/s (0.37 with the
magnets; 130 billion years). Titan's troposphere is driven by about a hundredth
of the sunlight, so its waves may be the weaker. Both put the Moon's mixing well
above Earth's own, so for the warmest air the atoms lose 0.1–0.8 kg/s.

**With the exosphere counted, the shield has to reach beyond the exobase.**
Largest allowed UV transmission for the titania stack with the sunlit
exosphere's central loss added, by protected radius, with no magnetosphere,
with no magnetosphere but the ring fleet's wake, and with the September
magnets. The ranges span the three treatments of the upper air and quiet Sun to
solar maximum; each radius carries its own traced state, and
a radius counts only where it covers the thermosphere's heating, light beyond its
aperture giving at most a tenth of the heat. Cases that cannot meet the budget at
any transmission are counted out of six. The ring fleet's tiles cover 4 lunar
radii; the other radii are sensitivities, and each lunar radius beyond 4 would
need about 23% more rings:

| Budget (cycle time) | Magnetosphere | 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.0056–0.099% (2 cannot) | 0.011–0.10% (2 cannot) | 0.012–0.10% (2 cannot) | 0.012–0.10% (2 cannot) |
| | none, the ring fleet's wake | 0.016–0.11% (2 cannot) | 0.017–0.14% (1 cannot) | 0.015–0.15% | 0.024–0.16% |
| | September | 0.012–0.15% (1 cannot) | 0.0091–0.17% | 0.018–0.18% | 0.025–0.18% |
| 10 kg/s (10 billion years) | none | 0.024–0.15% (1 cannot) | 0.021–0.19% | 0.031–0.21% | 0.037–0.21% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.022–0.19% | 0.034–0.22% | 0.042–0.23% |
| | September | 0.024–0.15% (1 cannot) | 0.028–0.21% | 0.037–0.22% | 0.043–0.23% |
| 100 kg/s (1 billion years) | none | 0.024–0.15% (1 cannot) | 0.032–0.25% | 0.066–0.28% | 0.081–0.30% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.032–0.25% | 0.067–0.29% | 0.084–0.30% |
| | September | 0.024–0.15% (1 cannot) | 0.032–0.25% | 0.073–0.29% | 0.085–0.30% |

Without a magnetosphere the two cases with all near-infrared heating cannot meet
1 kg/s at any radius up to 10 lunar radii; the ring fleet's wake, as the
sensitivity takes it, leaves only the one at solar maximum at 4 lunar radii, and
none from 6. With the September magnets every case meets 1 kg/s at 4 lunar
radii, the warmest at solar maximum only below 0.0091% (0.0026% in cycle 19's
year and 0.0054% under the escape model's 2.5). Without a magnetosphere the
stress cases alone allow 0.0053–0.040% (cycle 19's year) and 0.0078–0.032% (the
2.5) for 1 kg/s at 4 lunar radii, one of the three cases in each unable to meet
it. Behind the 200-nm
edge no case meets 1 kg/s without a magnetosphere at any radius; with the
September magnets two of six do at 4 lunar radii and five from 6, and five of
six meet 10 kg/s from 4 lunar radii with or without the magnets.

## Limits

The thermal column is molecular only, with no infrared cooling. The tidal step
uses the circular restricted three-body problem: the Moon's eccentric, inclined
orbit and the Sun's tide (under 1% of Earth's at the Hill radius) are left out,
and a molecule past two Hill radii counts as lost. Exobases beyond six lunar
radii, reached only at UV transmissions of a percent or more, take the
six-radius multiplier, which understates the tide there. The tide is applied to the
column's escape without feeding back on its structure, and the energy-limited
bound leaves it out. The exosphere in the absorption step is two-body and
collisionless, with no satellite particles and no depletion along each flight;
its ion production is an upper bound on pickup loss. The escaped cloud uses
assumed intact lifetimes and leaves out photon pressure, the Sun's tide within
the Earth–Moon region, and plasma. The heat comes from the limb tracing's
spherically symmetric, global-mean tables, interpolated between their heats, so
the onset of a runaway is found to the tables' spacing. The tables are traced
through air the column heats with its 'middle' shape, which swells more than
air heated where the tracing puts it, so the states, the heating radius and
the runaway onsets lean conservative; a retrace through air heated with the
traced shape would refine them. The absorption step's band-mean limb deposit is a
cross-check. The ring fleet's wake takes the wind to close linearly over eight
screen radii, an assumed length. The totals here leave out atomic oxygen and
hydrogen, which the oxygen step gives apart; its eddy mixing above the base sets
how much oxygen reaches the exobase, and the breaking gravity waves and Titan's
measured mixing bracket it. The waves leave out mean winds, tides and the day-night
circulation, and take Earth's wave sources for the lunar troposphere's, which no
model here has yet supplied. The base temperatures come from a global-mean,
radiative-only middle atmosphere. The glow's brightness is an assumption, and
it is taken to follow solar Lyman-α over the cycle (1.55 at solar maximum, 1.88
in cycle 19's year and 1.5 in the escape model's stress case); the interplanetary hydrogen that scatters it is ionized faster
at solar maximum, so the glow may rise less, which this leaves out. The solar-wind branch is a set of scalings with
assumed ranges; for comparison, Venus and Mars lose ions to the solar wind at
roughly a tenth of a kilogram a second, and scaling Titan's plasma-driven losses
to the solar wind at 1 AU gives a few kilograms a second (literature magnitudes,
not yet source-checked). The exosphere step's ion fates with a magnetosphere are a
comparison of timescales for draining, convection and recombination, with the
rates for N2+ applied to both ions; the open-field share, the convection
efficiency and the draining speed are ranges, and cusp entry and its sputtering,
the magnetosphere's own plasma pressure, Earth's magnetotail passages and the
extra escape a weak field can drive are left out. Its charge exchange runs on
straight paths as an upper bound, with a cross-section not re-read for this
step, and takes 0.1–1 of it; electron-impact ionization is left out. It counts
the exosphere to a third of the Hill radius, beyond which the tidal step's
escape takes over, and leaves out the photochemistry of the collisional air just
below the exobase in the sunlit ring. The allowed transmissions with the
exosphere use its central loss. The loss budget itself is the author's to
choose.
