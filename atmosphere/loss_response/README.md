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
the Moon, which decides how much of that oxygen reaches the exobase.
[infrared.py](infrared.py) gives the thermal column what it radiates with in the
infrared: CO2 above the base with its band's escape from HITRAN's lines, and the
base's atomic oxygen and NO. All ten are screening models; their methods and
limits are in their docstrings.

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # about 4 minutes; resumes if stopped
python -m atmosphere.loss_response.cycle        # FISM2's daily record by band, a few seconds
OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.limb_heat   # then the limb tracing, about 75 CPU minutes (cached)
python -m atmosphere.loss_response.model        # then the loss response, about two minutes
python -m atmosphere.loss_response.absorption   # then the protected radius, a few minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.exosphere   # then the exosphere's losses, about 20 minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.fate    # where the escaping air goes, about 15 s
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.oxygen  # atomic oxygen and hydrogen, about 40 minutes in three processes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.gravity_waves   # the waves' mixing, a minute
python -m pytest atmosphere/tests/test_tides.py atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py atmosphere/tests/test_exosphere.py atmosphere/tests/test_fate.py atmosphere/tests/test_cycle.py atmosphere/tests/test_oxygen.py atmosphere/tests/test_gravity_waves.py atmosphere/tests/test_infrared.py
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
`terluna.atmosphere.gravity-wave-mixing/1`). The infrared cooling and the oxygen
step need the radiative inputs
(`python -m atmosphere.radiative_convective.fetch_inputs --download`): HITRAN's
CO2 lines, and the middle atmosphere's photolysis.

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
7 October in place of the column's fixed 'middle' shape. Since 8 October the
column also radiates in the infrared ([infrared.py](infrared.py)): CO2's 15 µm
band at the middle atmosphere's 400 ppm, with that model's rates out of local
thermodynamic equilibrium and the band's escape from HITRAN's lines, and the
oxygen atoms' fine-structure lines and NO's 5.3 µm band from the base's
atomic oxygen and NO carried up, each in excess of what air at the base
temperature emits, so that air left at the base temperature neither cools nor
warms, as before. Beyond the column's domain the energy-limited expression of the September
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
heat where the tracing puts it in height, and again on 8 October with the column
radiating in the infrared. The design point adds the mean over the solar cycle
measured band by band.

**CO2 radiates most of the heat the upper air takes.** Warmed as the column
had it without infrared cooling, the air above the base would radiate 400–3,500
times the heat it takes through CO2's 15 µm band, though the band is still thick
there (about 3% of its photons leave straight up from the base). So the glow's
heat, laid within a few e-folds of pressure above the base, leaves as infrared,
heat laid higher is conducted down to where CO2 radiates it, and almost none
reaches the base. The air stays within a few kelvin of the base temperature for
three or four e-folds of pressure, and the column needs about three times the
heat to swell to the shadow's edge: the warmest titania case's exobase reaches
four lunar radii at 2.4×10⁻⁵ W/m², against 9×10⁻⁶ without the cooling, in the
limb tables' 'middle'-heated air. Behind the titania stack no light splits O2
below the base, so the base carries no atomic oxygen or NO; behind the 200-nm
edge its 2×10⁻⁴ of oxygen atoms radiate about 2% of the cooling. The atoms the
glow and the gaps make above the base are left out, so the cooling here is a
lower estimate.

**Earth's tide roughly doubles escape.** Across the loss response's sweep the
tide multiplies molecular escape by 1.6–4.9; at the allowed transmissions it is
1.8–2.7. A uniform lowering of the barrier, the feasibility baseline's factor,
gives more there, because it counts every molecule above the lowered barrier
as lost, while many of them are aimed away from the necks and fall back.

**Where the heat lands matters as much as how much.** The light through gaps
heats the thermosphere at three depths: Lyman-α within a few e-folds of
pressure above the base, the far ultraviolet in the lower thermosphere and the
extreme ultraviolet in its upper half
([the limb tracing](../middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07)).
Heat laid low is radiated with little warming of the air above. Placed where
the tracing puts it, the same heat leaves the exobase 17–42 K cooler than the
column's 'middle' shape gives, and the ultraviolet-driven loss 6–3,600 times
smaller, the factor largest for the coolest air.

**The sky's glow heats the air just above the base, and CO2 radiates it.** O2
takes the glow's Lyman-α within a few e-folds of pressure above the 0.3 Pa base:
with the UV transmission at 10⁻⁶ it raises the exobase by 0.1–0.3 K at quiet
Sun and at most 0.8 K in cycle 19's year. Behind the titania stack, with the
films' traced X-rays, the floor loss is 6×10⁻¹¹–6×10⁻⁴ kg/s at quiet Sun and
10⁻¹⁰–9×10⁻⁴ kg/s at solar maximum (2×10⁻¹⁰–0.001 in cycle 19's year and
10⁻¹⁰–8×10⁻⁴ in the escape model's). Behind the 200-nm edge, whose upper air
starts warmer, it is 0.004–0.1 kg/s at quiet Sun and 0.006–0.1 kg/s at solar
maximum (0.007–0.2 and 0.006–0.1 kg/s); with all near-infrared heating the air
there settles at every maximum.

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

**Ultraviolet-driven escape alone lets the gaps pass a few tenths of a
percent.** Largest allowed UV transmission from ultraviolet-driven escape
and the solar-wind scalings, across the three upper-air treatments, with the
solar-wind loss at the low end of its range. Each budget's cycle time is the
atmosphere's mass (3.1×10¹⁸ kg, the feasibility baseline's hydrostatic column)
over the loss rate, the time the loss and its resupply take to replace the
whole atmosphere:

| Shield and Sun | 1 kg/s (100 billion years) | 10 kg/s (10 billion years) | 100 kg/s (1 billion years) |
|---|---|---|---|
| Titania, quiet Sun | 0.20–0.52% | 0.29–0.64% | 0.35–0.79% |
| Titania, solar maximum (cycle 21) | 0.12–0.31% | 0.17–0.38% | 0.21–0.46% |
| Titania, cycle 19's year (stress) | 0.099–0.27% | 0.14–0.32% | 0.17–0.38% |
| Titania, the escape model's 2.5 (stress) | 0.074–0.20% | 0.11–0.23% | 0.13–0.27% |
| 200-nm edge, quiet Sun | 0.058–0.14% | 0.12–0.22% | 0.12–0.24% |
| 200-nm edge, solar maximum (cycle 21) | 0.031–0.078% | 0.061–0.13% | 0.061–0.13% |
| 200-nm edge, cycle 19's year (stress) | 0.024–0.065% | 0.043–0.11% | 0.043–0.11% |
| 200-nm edge, the escape model's 2.5 (stress) | 0.020–0.053% | 0.038–0.083% | 0.038–0.083% |

These are 2.8–5.5 times what the column allowed without infrared cooling behind
the titania stack, and 2.6–85 times behind the 200-nm edge (the most in cycle
19's year), whose warmest case now settles at every maximum. From 10 to 100 kg/s the allowed
transmission rises by a fifth behind the titania stack and not at all behind
the 200-nm edge: the runaway caps it before the loss does. The measured maximum
allows 1.5–1.7 times what the escape model's 2.5 does, its far ultraviolet being
weaker, and cycle 19's year allows 1.1–1.4 times as much as the 2.5.

**With a 4-radius shadow the air runs away near the shadow's edge.** As the heat
swells the air, its exobase approaches the aperture's edge and the unfiltered
light beyond it heats the thermosphere more, until the heat outgrows itself.
Behind the titania stack that begins at a transmission of 0.48–0.82% with
collisional upper air (solar maximum to quiet Sun; 0.40% in cycle 19's year and
0.27% in the escape model's), 0.36–0.61% in LTE (0.29% and 0.20%) and
0.21–0.35% with all near-infrared heating (0.17% and 0.13%), and behind the
200-nm edge at 0.13–0.24% (0.11% and 0.082%), with all near-infrared heating
0.061–0.12% (0.043% and 0.038%). These are three to eleven times the onsets
without infrared cooling, and 85 times behind the 200-nm edge in cycle 19's
year. The loss just short of an onset is already 140–190 kg/s, 96–120 kg/s and
37–38 kg/s behind the titania stack and 5–16 kg/s behind the edge, down to 125,
78, 27 and 3.3 kg/s in the stress cases, so behind the titania stack the
budgets up to 10 kg/s bind before the runaway does. These onsets come
from tables traced through air the column heats with its middle shape, which
swells more than the traced shape's: at each onset the traced column's own
exobase still lies at 3.4–3.9 lunar radii, so the onsets lean early.

**The runaway is a latch, and the infrared cooling puts its threshold beyond
every maximum on record.** Every case has two regimes, with the cooling as
without it. Below a threshold heat the air settles to its compact state; above
it the exobase lies far enough past the shadow's edge that the unfiltered
sunlight beyond the aperture heats the outer air more than the heat that swells
it, and the air swells past the tables' top, beyond six lunar radii. At quiet
Sun with no gaps the threshold lies at 3.9–8.4×10⁻⁵ W/m² behind the titania
stack, the exobase then at 5.0–5.8 lunar radii, and at 1.9–3.0×10⁻⁵ W/m² behind
the 200-nm edge (4.5–4.7 lunar radii), against compact states of
2.1–2.8×10⁻⁶ W/m². The second regime comes from the geometry: beyond about 4.5
lunar radii the outer air's take of that light doubles with every half to
three-quarters of a lunar radius the exobase gains, while the heat the column
needs to put the exobase there rises by about a quarter. The cooling raises
that heat and so moves the thresholds outward and up three and a half times,
from 0.5–1.9×10⁻⁵ W/m² near 4.2–4.6 lunar radii in the uncooled tables, and
leaves the second regime in place. A maximum that carried the air past its threshold
would still leave it swollen after the Sun quiets. At the standard level none
does: the compact states of every maximum on record, cycle 19's year and the
escape model's 2.5 included, take a tenth to three-fifths of the threshold heat
at the same activity, the least margin behind the 200-nm edge with all
near-infrared heating. In transmission, the gaps would have to pass 0.35–0.82%
at quiet Sun behind the titania stack to tip the air over (0.12–0.24% behind
the edge). The thresholds come from the same tables as the onsets and lean
early in the same way. What the column leaves out could move them either way:
the oxygen atoms the glow and the gaps make above the base would strengthen the
cooling, since an oxygen atom excites CO2 about a thousand times faster than N2
does; in the swollen regime the unfiltered light beyond the aperture would split
much of the outer air's O2, and its atoms would radiate at 63 µm but also lighten
the air and escape more readily.

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
until it runs away. At the allowed transmissions the heating needs 3.5 lunar
radii for 1 kg/s, 3.5–4 for 10 kg/s and 4 for 100 kg/s behind the titania stack,
and 3.5–4 behind the 200-nm edge: the first tabulated radii beyond the exobase
of the swollen air the tables trace. The ring fleet's 4 lunar radii cover the
heating at every budget. The extreme ultraviolet itself is absorbed lower, with
unit optical depth at grazing incidence at 2.0–2.6 lunar radii, and the thin air
above that still takes up too much.

**The sunlit exosphere beyond the shadow makes ions faster than the budget.**
At the allowed transmissions the exobase sits at 2.5–3.8 lunar radii, where
gravity is weak. Outside a shadow sized for heating, its ballistic molecules
become ions at 2–8 kg/s for the 1 kg/s budget, 7–23 kg/s for 10 kg/s and
13–71 kg/s for 100 kg/s, counted to a third of the Hill radius; counting to the
whole Hill radius adds up to three fifths. A shadow that held production to a
tenth of the budget would need 4–10 lunar radii counted to a third of the Hill
radius, 4–14 to half of it, or 5–24 to the whole; each lunar radius beyond the
ring fleet's 4 would take about 23% more rings. The exosphere step follows these
ions.

**At those transmissions the sunlit exosphere overruns every budget.** Where
ultraviolet-driven escape already fills the budget, the sunlit exosphere
outside the ring fleet's 4-radius shadow adds a central 4.6–190 kg/s with no
magnetosphere (cycle times of 22 to 0.5 billion years on its own) and
1.5–114 kg/s with the September magnets (66 to 0.9 billion years). The solar
wind's charge exchange alone makes up to 5–20 kg/s of ions there. Most of it
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
central estimate the design point's exosphere then loses 0.6–1.8 kg/s at the
standard level, most or all of a 1 kg/s budget.

**The ring fleet's wake may shelter the inner exosphere.** The ring fleet's
screen, about 20,000 km out, lies under three screen radii from the Moon, and
the wind closes in behind an absorbing screen over about eight. If the wind is
absent on paths within the wake's core at the Moon (2.6 lunar radii behind a
4-radius screen), the charge exchange of the cool upper air nearly vanishes: at
the design point's standard level the three coolest cases lose 0.063–0.065 kg/s
without a magnetosphere, against 0.60–0.79 kg/s with the wind everywhere. Where
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
much weaker dipole does that: at the design point's standard level the
collisional and LTE cases lose under 1 kg/s from their exosphere without one at
every maximum on record, and the warmest needs about 3×10¹⁹ A·m² (a stand-off
near 2.7 lunar radii, a fiftieth of the September moment) at quiet Sun and at
every maximum, cycle 19's year and the escape model's 2.5 included. Dipoles that hold the wind off only just above the exobase fall in the
range where hybrid simulations find a weak field can raise ion escape (Egan et
al. 2019), which this step leaves out.

**Magnets do not hold the fragments.** Sunlight breaks molecules apart in the
sunlit exosphere faster than it ionizes them, and the atoms are neutral. The
escaping fragments are 69–81% as much as the ions made, and with a
magnetosphere they are most of what the sunlit exosphere loses. Only the shadow
removes them.

**Atomic oxygen made above the base matters only for warm upper air with weak
mixing.** Behind the titania stack nothing splits O2 below the base, but above it
the glow makes 7–15 kg/s of oxygen atoms and the light through gaps 3–9 kg/s
more at the standard level, from quiet Sun to the solar maxima
([oxygen.py](oxygen.py)). With the middle atmosphere's eddy mixing, scaled for
the Moon, almost all of it goes back down as odd oxygen and builds 1.6–2.8 DU of
ozone in the middle atmosphere. The homopause then lies at 2–4×10⁻⁵ Pa and the
atoms are 0.06–0.6% of the gas at the exobase. They lose 0.002–0.003 kg/s with
collisional upper air and in LTE at the solar cycle's mean spectrum, and
0.021 kg/s averaged over the cycle's years with all near-infrared heating
(0.007 kg/s with the September magnets), 0.035 kg/s in its worst year; the
infrared cooling, which keeps that case's air cooler and more compact, cut its
atoms' loss about sixfold. Most of it
leaves as ions the solar wind makes from the sunlit exosphere's atoms. Atoms made
within a collision length or so of the exobase leave hot, since photolysis and
recombination give them more than the escape energy there; they add
0.001–0.003 kg/s over the cycle (at most 0.009 kg/s under the escape model's 2.5, letting
them through two collisions), about half the loss of the coolest case. Earth's own
mixing, unscaled, puts the homopause at 0.7–1.2×10⁻³ Pa: the atoms are then
2–12% of the exobase gas and lose 0.020–0.15 kg/s in the cooler cases and
0.46 kg/s with all near-infrared heating over the cycle (0.13 kg/s with the
magnets). Behind the 200-nm edge, whose middle atmosphere already holds atomic
oxygen under about 300 DU of ozone, the atoms add 0.6–6 kg/s. Hydrogen from water
and H2 escapes at about 0.0025 kg/s.

**Breaking gravity waves mix the upper air about as strongly as the Moon-scaled
profile.** How strongly the air above the base mixes decides the atoms
([gravity_waves.py](gravity_waves.py)). Gravity waves launched from half the
tropopause pressure with Earth's 4 mPa, the calculation calibrated so that it
gives Earth's measured eddy diffusion of 1.1×10⁶ cm²/s at 96 km (Swenson et al.
2021) and its homopause at 110 km, put the Moon's homopause for atomic oxygen at
about 10⁻⁵ Pa, a little above the Moon-scaled profile's. The lunar thermosphere's
buoyancy frequency falls to about 10⁻³ per second, so short, fast waves turn back
there, while waves a few hundred kilometres long, which a troposphere six times
deeper than Earth's favours, rise almost undamped and break near the exobase.
With the waves' mixing above the base the warmest case's atoms lose 0.017 kg/s
over the cycle (0.006 kg/s with the magnets, which brings that case to
0.050 kg/s, two trillion years) and the cooler
cases 0.002–0.003 kg/s; with Earth's wavelengths and a third of the calibration,
the low end of Earth's measured mixing, the warmest case's atoms lose 0.079 kg/s
(0.023). Titan, the nearest analogue, a 1.5-bar nitrogen atmosphere at
1.35 m/s², mixes less: its measured homopause eddy diffusion, 2×10⁷–10⁸ cm²/s
(Yelle et al. 2008; Bell et al. 2014), carried to the Moon is about a tenth of
the Moon-scaled mixing, and with it the atoms lose 0.15 kg/s (0.044 with the
magnets, bringing that case to 0.088 kg/s, a trillion years). Titan's troposphere is driven by about a hundredth
of the sunlight, so its waves may be the weaker. Both put the Moon's mixing well
above Earth's own, so for the warmest air the atoms lose 0.02–0.15 kg/s.

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
| 1 kg/s (100 billion years) | none | 0.056–0.28% (2 cannot) | 0.072–0.30% (2 cannot) | 0.076–0.30% (2 cannot) | 0.077–0.31% (2 cannot) |
| | none, the ring fleet's wake | 0.086–0.32% (2 cannot) | 0.038–0.42% | 0.087–0.49% | 0.11–0.52% |
| | September | 0.024–0.43% | 0.070–0.48% | 0.10–0.51% | 0.11–0.53% |
| 10 kg/s (10 billion years) | none | 0.061–0.47% | 0.11–0.56% | 0.13–0.60% | 0.14–0.62% |
| | none, the ring fleet's wake | 0.061–0.47% | 0.11–0.57% | 0.14–0.65% | 0.16–0.69% |
| | September | 0.089–0.47% | 0.13–0.60% | 0.15–0.67% | 0.17–0.69% |
| 100 kg/s (1 billion years) | none | 0.089–0.47% | 0.19–0.73% | 0.23–0.85% | 0.26–0.87% |
| | none, the ring fleet's wake | 0.089–0.47% | 0.19–0.74% | 0.23–0.87% | 0.26–0.90% |
| | September | 0.089–0.47% | 0.21–0.76% | 0.24–0.88% | 0.26–0.90% |

Without a magnetosphere the two cases with all near-infrared heating cannot meet
1 kg/s at any radius up to 10 lunar radii; the ring fleet's wake, as the
sensitivity takes it, lets every case meet it from 4 lunar radii. With the
September magnets every case meets 1 kg/s at 4 lunar radii, the warmest at solar
maximum below 0.070% (0.057% in cycle 19's year and 0.043% under the escape
model's 2.5). Without a magnetosphere the stress cases alone allow 0.059–0.13%
(cycle 19's year) and 0.050–0.10% (the 2.5) for 1 kg/s at 4 lunar radii, one of
the three cases in each unable to meet it. Behind the 200-nm edge no case meets
1 kg/s without a magnetosphere at any radius; with the September magnets four of
six do at 4 lunar radii and all six from 6, and all six meet 10 kg/s from 4
lunar radii with or without the magnets.

## Limits

The thermal column is molecular. Its infrared cooling takes CO2 at the middle
atmosphere's 400 ppm as one two-level band whose photons escape with the mean of
their chances up and down (a level's exchange with the others is not solved
band by band), air at the base temperature as in radiative balance, and the
base's atomic oxygen and NO carried up; the oxygen atoms made above the base,
which would strengthen it, and CO2's near-infrared heating above the base are
left out. The tidal step
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
