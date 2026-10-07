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
around Earth and the Moon. All six are screening models; their methods and
limits are in their docstrings.

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # about 4 minutes; resumes if stopped
OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.limb_heat   # then the limb tracing, about 2 CPU hours (cached)
python -m atmosphere.loss_response.model        # then the loss response, about a minute
python -m atmosphere.loss_response.absorption   # then the protected radius, a few minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.exosphere   # then the exosphere's losses, about 20 minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.fate    # where the escaping air goes, about 15 s
python -m pytest atmosphere/tests/test_tides.py atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py atmosphere/tests/test_exosphere.py atmosphere/tests/test_fate.py
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
`terluna.atmosphere.escape-fate/1`).

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
reaching the disk. Solar maximum takes 2.5 on the ultraviolet and FISM2's
measured rise of the X-rays below 10 nm by band
([escape.py](../middle_atmosphere/escape.py)). The films are the design's behind
both shields, and the protected radius is the ring fleet's 4 lunar radii unless
stated. The thermal column turns the heat into an exobase and a
molecular Jeans loss, over base temperatures from the middle atmosphere behind
the titania stack and behind the 200-nm edge, for each treatment of the upper
air. Beyond the column's domain the energy-limited expression of the September
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
ring fleet's protected radius of 4 lunar radii.

**Earth's tide roughly doubles escape.** Across the loss response's sweep the
tide multiplies molecular escape by 1.8–4.9; at the allowed transmissions it is
1.8–2.6. A uniform lowering of the barrier, the feasibility baseline's factor,
gives 2.6–4.5 there, because it counts every molecule above the lowered barrier
as lost, while many of them are aimed away from the necks and fall back.

**The glow sets a floor the shield cannot remove.** With the UV transmission at
10⁻⁶ the glow alone raises the exobase by 31–42 K at quiet Sun and 48–61 K at
solar maximum. Behind the titania stack, with the films' traced X-rays, the
floor loss is 3×10⁻⁷–0.1 kg/s at quiet Sun and 2×10⁻⁵–1.5 kg/s at solar
maximum. Behind the 200-nm edge, whose upper air starts warmer, it is
0.5–6.4 kg/s at quiet Sun and 7–8 kg/s at solar maximum, and with all
near-infrared heating at solar maximum the air runs away with no gaps at all.

**Traced, ultraviolet-driven escape alone lets the gaps pass a few hundredths
of a percent.** Largest allowed UV transmission from ultraviolet-driven escape
and the solar-wind scalings, across the three upper-air treatments, with the
solar-wind loss at the low end of its range. Each budget's cycle time is the
atmosphere's mass (3.1×10¹⁸ kg, the feasibility baseline's hydrostatic column)
over the loss rate, the time the loss and its resupply take to replace the
whole atmosphere:

| Budget (cycle time) | Titania, quiet Sun | Titania, solar maximum | 200-nm edge, quiet Sun | 200-nm edge, solar maximum |
|---|---|---|---|---|
| 1 kg/s (100 billion years) | 0.025–0.15% | 0.024–0.046% (none for the warmest upper air) | 0.0058–0.0064% (none for the warmest) | none |
| 10 kg/s (10 billion years) | 0.055–0.18% | 0.010–0.059% | 0.0053–0.033% | 0.0016–0.0019% (none for the warmest) |
| 100 kg/s (1 billion years) | 0.081–0.25% | 0.021–0.076% | 0.014–0.049% | 0.0087–0.0088% (none for the warmest) |

That is a third or less of what a count over the disk allowed (0.088–0.45% for
1 kg/s behind the titania stack at quiet Sun). Solar maximum with the warmest
upper air is the binding case.

**With a 4-radius shadow the air runs away near the shadow's edge.** As the heat
swells the air, its exobase approaches the aperture's edge and the unfiltered
light beyond it heats the thermosphere more. The largest steady loss is
147 kg/s with collisional upper air, 103 kg/s in LTE and 36 kg/s with all
near-infrared heating behind the titania stack, and 28–29 kg/s behind the
200-nm edge (19 kg/s with all near-infrared heating at quiet Sun); in each case
the exobase then sits at 3.9–4.1 lunar radii. Past that the air runs away, so a
budget looser than those losses asks a wider shadow.

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
until it runs away. At the allowed transmissions the heating needs
2.75–3.5 lunar radii for 1 kg/s, 3.5 for 10 kg/s and 4 for 100 kg/s behind the
titania stack, and 3.5–4 behind the 200-nm edge: the first tabulated radii
beyond the exobase (the tables run from 2 to 10 lunar radii in steps of a
quarter to a half radius). The ring fleet's 4 lunar radii cover the heating at
every budget up to 100 kg/s, at the edge for 100 kg/s. The extreme ultraviolet
itself is absorbed lower, with unit optical depth at grazing incidence at
2.0–2.7 lunar radii, and the thin air above that still takes up too much.

**The sunlit exosphere beyond the shadow makes ions faster than the budget.**
The exosphere is large because the exobase sits at 2.7–4.1 lunar radii at
these transmissions, where gravity is weak. Outside a shadow sized for heating,
its ballistic molecules become ions at 4–18 kg/s for the 1 kg/s budget,
12–53 kg/s for 10 kg/s and 27–121 kg/s for 100 kg/s, counted to a third of the
Hill radius. Counting to the whole Hill radius adds 4–7% at 1 kg/s and a
quarter to a half at 100 kg/s. A shadow that held production to a tenth of the
budget would need a radius of 6–13 lunar radii counted to half the Hill
radius, or 7–25 counted to the whole. The exosphere step follows these ions.

**At those transmissions the sunlit exosphere overruns every budget.** At the
allowed transmissions of the table above, where ultraviolet-driven escape
already fills the budget, the sunlit exosphere outside the ring fleet's
4-radius shadow adds a central loss of 5–273 kg/s with no magnetosphere (cycle
times of 20 to 0.4 billion years on its own) and 1.7–168 kg/s with the
September magnets (58 to 0.6 billion years). The solar wind's charge exchange
alone makes up to 6–23 kg/s of ions there. Most of it comes from the first few
scale heights above the exobase.

**Without a magnetosphere the ions leave.** In the flowing solar wind the
motional electric field pulls an ion about 40,000 times harder than lunar
gravity at three lunar radii, and gravity holds an ion only where the flow
stops, at an obstacle's surface. The exosphere's ions are carried off, as at
comets, whose ion production likewise far exceeds the wind's mass flux through
the region that makes it: the loaded wind slows and its interaction region
grows to tens of lunar radii. A tenth to a half curve back into the air and
sputter out 1–10 molecules each. The solar wind's charge exchange with the dense
exosphere near the exobase makes ions inside the optical shadow too; in the
central estimate it costs 0.7–6 kg/s at the design point's exobases, which
alone uses a 1 kg/s budget.

**The ring fleet's wake may shelter the inner exosphere.** The ring fleet's
screen, about 20,000 km out, lies under three screen radii from the Moon, and
the wind closes in behind an absorbing screen over about eight. If the wind is
absent on paths within the wake's core at the Moon (2.6 lunar radii behind a
4-radius screen), the charge exchange of the coolest upper air nearly
vanishes: at the design point's standard level the three coolest cases lose
0.07–0.31 kg/s without a magnetosphere, against 0.77–1.3 kg/s with the wind
everywhere. Where the exobase lies outside the core, with warmer upper air or at
the transmissions that fill a budget, the wake changes little. The refill
length is an assumption a plasma model has to test.

**A lunar magnetosphere holds the wind off and keeps part of the ions.** With
the September magnets, of the ions made outside the 4-radius shadow at the
allowed transmissions, 25–41% drain along closed field lines into the air and
stay, 32–33% leave along the open field lines over the poles, and
recombination into fast atoms and the wind-driven convection take the rest; the
warmer the upper air, the farther its exosphere reaches and the less the
magnetosphere keeps. The magnetosphere also keeps the solar wind, and its
charge exchange, off the exosphere inside the stand-off. For the coolest upper
air a much weaker dipole does that: at the design point's standard level the
coolest case loses under 1 kg/s from its exosphere without one, and the next
two need about 2×10¹⁹ A·m² (a stand-off near 2.3 lunar radii, just above their
exobase, a seventy-fifth of the September moment). No moment up to
10²³ A·m² holds the warmer three to 1 kg/s. Dipoles that hold the wind off only
just above the exobase fall in the range where hybrid simulations find a weak
field can raise ion escape (Egan et al. 2019), which this step leaves out.

**Magnets do not hold the fragments.** Sunlight breaks molecules apart in the
sunlit exosphere faster than it ionizes them, and the atoms are neutral. The
escaping fragments are 67–81% as much as the ions made, and with a
magnetosphere they are most of what the sunlit exosphere loses. Only the shadow
removes them.

**With the exosphere counted, the shield has to reach beyond the exobase.**
Largest allowed UV transmission for the titania stack with the sunlit
exosphere's central loss added, by protected radius, with no magnetosphere,
with no magnetosphere but the ring fleet's wake, and with the September
magnets. The ranges span the three treatments of the upper air and quiet Sun to
solar maximum; each radius carries its own traced state, and a radius counts
only where it covers the thermosphere's heating, light beyond its aperture
giving at most a tenth of the heat. Cases that cannot meet the budget at any
transmission are counted out of six:

| Budget (cycle time) | Magnetosphere | 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.007–0.052% (4 cannot) | 0.011–0.057% (4 cannot) | 0.0005–0.059% (3 cannot) | 0.0006–0.060% (3 cannot) |
| | none, the ring fleet's wake | 0.012–0.065% (3 cannot) | 0.006–0.11% (2 cannot) | 0.013–0.13% (1 cannot) | 0.023–0.13% (1 cannot) |
| | September | 0.003–0.11% (2 cannot) | 0.005–0.13% (1 cannot) | 0.018–0.14% (1 cannot) | 0.023–0.15% (1 cannot) |
| 10 kg/s (10 billion years) | none | 0.0015–0.14% (1 cannot) | 0.022–0.14% (1 cannot) | 0.0012–0.15% | 0.006–0.16% |
| | none, the ring fleet's wake | 0.002–0.14% (1 cannot) | 0.023–0.15% (1 cannot) | 0.003–0.17% | 0.009–0.18% |
| | September | 0.017–0.15% (1 cannot) | 0.027–0.17% (1 cannot) | 0.005–0.18% | 0.009–0.18% |
| 100 kg/s (1 billion years) | none | 0.024–0.15% (1 cannot) | 0.012–0.22% | 0.022–0.25% | 0.032–0.26% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.012–0.22% | 0.022–0.25% | 0.034–0.27% |
| | September | 0.024–0.15% (1 cannot) | 0.016–0.23% | 0.025–0.26% | 0.034–0.27% |

The case that cannot meet 1 kg/s with the September magnets at any radius up to
10 lunar radii is the warmest upper air at solar maximum. Without a
magnetosphere the solar wind's charge exchange keeps three or four of the six
cases above 1 kg/s at any radius, and the ring fleet's wake, as the sensitivity
takes it, leaves one to three. Behind the 200-nm edge no case meets 1 kg/s
without a magnetosphere, and with the September magnets two of six do, from
8 lunar radii; 10 kg/s holds for most cases only from 10 lunar radii, and
100 kg/s from 6.

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
the onset of a runaway is found to the tables' spacing; the column spreads that
heat with its default shape, while the traced light lands between its low and
middle shapes, and the losses differ many times between them (the middle
atmosphere's limb tracing). The absorption step's band-mean limb deposit is a
cross-check. The ring fleet's wake takes the wind to close linearly over eight
screen radii, an assumed length. Atomic-oxygen escape, photochemical escape and hydrogen
from water are left out. The base temperatures come from a global-mean,
radiative-only middle atmosphere. The glow's brightness and its solar-cycle
factor (1.5) are assumptions. The solar-wind branch is a set of scalings with
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
