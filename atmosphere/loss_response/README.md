# Loss response to protection

How fast the finished 1.2 atm atmosphere loses mass as a function of its UV
transmission, the share of sunlight below 175 nm that the optical shield lets
through, and of how much of the solar wind reaches the air. It sizes the two
protection functions against a total loss budget for the protection design study
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).
[tides.py](tides.py) finds how much Earth's tide raises escape.
[absorption.py](absorption.py) finds how far out the optical shield must reach
and what the sunlit exosphere beyond it produces. [exosphere.py](exosphere.py)
finds what that exosphere loses, as ions and as the fragments of broken
molecules, with and without a lunar magnetosphere, and the UV transmission each
budget then allows. [fate.py](fate.py) follows the escaping air into the space
around Earth and the Moon. All five are screening models; their methods and
limits are in their docstrings.

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # about 4 minutes; resumes if stopped
python -m atmosphere.loss_response.model        # then the loss response, about 10 s
python -m atmosphere.loss_response.absorption   # then the protected radius, about 10 s
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.exosphere   # then the exosphere's losses, about 3 minutes
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.fate    # where the escaping air goes, about 15 s
python -m pytest atmosphere/tests/test_tides.py atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py atmosphere/tests/test_exosphere.py atmosphere/tests/test_fate.py
```

It needs the middle atmosphere's stored results and the WHI 2008 spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`). The product
is [results/loss_response.json](results/loss_response.json) (schema
`terluna.atmosphere.loss-response/2`), with the full sweep in
[results/euv_sweep.csv](results/euv_sweep.csv). Its field
`allowed_leak_fraction` is the allowed UV transmission from ultraviolet-driven
escape and the solar-wind scalings alone. The exosphere's losses, and the allowed
transmissions with them included, are in
[results/exosphere_loss.json](results/exosphere_loss.json) (schema
`terluna.atmosphere.exosphere-loss/1`). The tidal multipliers are
in [results/tidal_escape.json](results/tidal_escape.json) (schema
`terluna.atmosphere.tidal-escape/1`), and the protected radius and the
exosphere's ion production in
[results/absorption_radius.json](results/absorption_radius.json) (schema
`terluna.atmosphere.absorption-radius/1`). The cloud the escaping air forms is in
[results/escape_fate.json](results/escape_fate.json) (schema
`terluna.atmosphere.escape-fate/1`).

## How it works

The ultraviolet branch deposits heat above the 0.3 Pa base from three sources:
the transmitted share f of sunlight below 175 nm, the titania film's own hard
X-rays, and the sky's interplanetary Lyman-alpha glow (about 1,000 rayleigh,
4.1×10⁻⁶ W/m² onto a surface facing open sky), which reaches the upper air from
every direction. The thermal column turns the heat into an exobase and a
molecular Jeans loss, over base temperatures from the middle atmosphere behind
the titania stack and behind the 200-nm edge, for each treatment of the upper
air. Beyond the column's domain the energy-limited expression of the September
feasibility report bounds the loss. The solar-wind branch combines pickup,
capped by mass loading at the solar wind's own mass flux through a disc three
lunar radii across (0.29 kg/s), with sputtering by returning pickup ions and by
precipitating protons, over stated parameter ranges.

The transmitted share counts as a quarter of the light that reaches the disk,
all of it above the base. The middle atmosphere's slant-path tracing
([limb_heat](../middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07),
7 October) follows light through gaps over the whole aperture and finds it heats
the thermosphere 2.5–5 times as much, as it reaches the air above the limb.
It also finds the X-rays' solar-maximum factor from FISM2: 4–6 below 10 nm and
13–24 for the light the titania films pass, where this branch takes 2.5 and 100.
Whether this branch takes either is the author's choice.

Earth's tide lowers the barrier a molecule must climb to leave the Moon: in the
frame turning with the Moon's orbit, a molecule can slip out of the Moon's Hill
sphere through the necks around the L1 and L2 points, 58,000 km from the Moon,
while slower than the Moon's own escape speed. The tidal step launches test
molecules from the exobase in the Earth–Moon three-body problem and counts how
many of those with energy enough to pass L1 actually leave; most that do take
about five days. Each species' Jeans loss is multiplied by the result for its
exobase radius and temperature.

The protected radius comes from the thermal column's profile at each budget's
allowed transmission. Sunlight passing outside a shadow cylinder of radius R
crosses the collisional air below the exobase at a grazing angle and leaves part
of its ultraviolet there; the radius is where that deposit equals a tenth of the
allowed transmission. Above the exobase the air is collisionless. Chamberlain's
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
molecules and turn them into ions, inside the optical shadow as well. With no
magnetosphere the wind carries the ions off, and those that curve back into the
air sputter out more molecules. A lunar dipole holds the wind off at its
stand-off, 10 lunar radii for the September design's magnets
(1.5×10²¹ A·m²). Inside it, ions on open field lines leave down the tail, and
ions on closed field lines drain along the field into the air and stay, or are
carried to the boundary by the convection the wind drives, or recombine into
fast atoms that mostly escape, in proportion to the three rates.

## Results

**Earth's tide roughly doubles escape.** Across the loss response's sweep the
tide multiplies molecular escape by 1.4–4.5; at the allowed transmissions it is
1.8–2.6. A uniform lowering of the barrier, the feasibility baseline's factor,
gives 2.6–4.5 there, because it counts every molecule above the lowered barrier
as lost, while many of them are aimed away from the necks and fall back.

**The glow sets a floor the shield cannot remove.** With the UV transmission at
10⁻⁶ the glow alone raises the exobase by 31–42 K at quiet Sun and 47–61 K at
solar maximum. Behind the titania stack the floor loss is 3×10⁻⁷–0.09 kg/s at
quiet Sun and up to 0.7 kg/s at solar maximum. Behind the 200-nm edge, whose
upper air starts warmer, it is 0.5–6 kg/s at quiet Sun and 2.7–21 kg/s at solar
maximum.

**Ultraviolet-driven escape alone lets the shield pass a few tenths of a
percent.** Largest allowed UV transmission from ultraviolet-driven escape and
the solar-wind scalings, across the three upper-air treatments, with the
solar-wind loss at the low end of its range. Each budget's cycle time is the atmosphere's mass
(3.1×10¹⁸ kg, the feasibility baseline's hydrostatic column) over the loss rate,
the time the loss and its resupply take to replace the whole atmosphere:

| Budget (cycle time) | Titania, quiet Sun | Titania, solar maximum | 200-nm edge, quiet Sun | 200-nm edge, solar maximum |
|---|---|---|---|---|
| 1 kg/s (100 billion years) | 0.088–0.45% | 0.004–0.15% | 0.027–0.031% (none for the warmest upper air) | none |
| 10 kg/s (10 billion years) | 0.21–0.64% | 0.050–0.21% | 0.028–0.15% | 0.026–0.028% (none for the warmest upper air) |
| 100 kg/s (1 billion years) | 0.51–0.97% | 0.17–0.36% | 0.29–0.44% | 0.081–0.14% |

Solar maximum with the warmest upper air is the binding case. The ozone-forming
edge cannot meet 1 kg/s at solar maximum even with a perfect shield, nor 10 kg/s
in its warmest case, and its figures leave out atomic-oxygen escape, which adds
about half again.

**The solar-wind branch as screened.** Its range is 0.03 kg/s (low),
0.19 kg/s (central) and 2.1 kg/s (high); at the high end sputtering by
returning pickup ions dominates. With pickup capped at 0.29 kg/s it barely
changes the allowed transmission for budgets of 10 kg/s and more, and at 1 kg/s
its high end alone exceeds the budget. The ion production below exceeds that cap
36–470 times, and the exosphere step finds that those ions are carried off; its
losses replace the capped pickup and the ion sputtering, and the table above
holds for the ultraviolet-driven part alone.

**For the thermosphere's heating the shield must reach to just under the
exobase.** Ultraviolet grazing the
thermosphere outside the shadow deposits a tenth of the allowed transmission
when the protected radius is within about 0.1 lunar radii of the exobase:
2.6–3.1 lunar radii for 1 kg/s, 3.1–3.7 for 10 kg/s and 3.9–4.6 for 100 kg/s. The
aperture at 78,000 km then has 0.78–1.04, 1.04–1.5 and 1.6–2.3 times the
September design's area. The September three lunar radii cover five of the
eight 1 kg/s cases and none at 10 kg/s or more. The extreme ultraviolet itself
is absorbed lower, with unit optical depth at grazing incidence at 2.1–2.9
lunar radii, and the thin air above that still takes up too much.

**The sunlit exosphere beyond the shadow makes ions faster than the budget.**
The exosphere is large because the exobase sits at three to five lunar radii,
where gravity is weak. Outside a shadow sized for heating, its ballistic
molecules become ions at 10–30 kg/s for the 1 kg/s budget, 22–58 kg/s for 10 kg/s
and 47–134 kg/s for 100 kg/s, counted to a third of the Hill radius. Counting to
the whole Hill radius adds 2–3% at 1 kg/s and about half at 100 kg/s. A shadow
that held production to a tenth of the budget would need a radius of 9–14 lunar
radii counted to half the Hill radius, or 17–26 counted to the whole. The
exosphere step follows these ions.

**At those transmissions the sunlit exosphere overruns every budget.** At the
allowed transmissions of the table above, where ultraviolet-driven escape
already fills the budget, the sunlit exosphere outside a heating-sized shadow
makes ions at 11–136 kg/s, sheds fragments that escape at 7–91 kg/s, and meets
the solar wind, whose charge exchange makes up to 6–27 kg/s more ions. With no
magnetosphere the central loss it adds is 25–303 kg/s (cycle times of 4 to 0.3
billion years on its own); with the September magnets it is 11–193 kg/s (9 to
0.5 billion years). Most of it comes from the first few scale heights above the
exobase, just outside a shadow that reaches only to the exobase.

**Without a magnetosphere the ions leave.** In the flowing solar wind the
motional electric field pulls an ion about 40,000 times harder than lunar
gravity at three lunar radii, and gravity holds an ion only where the flow
stops, at an obstacle's surface. The exosphere's ions are carried off, as at
comets, whose ion production likewise far exceeds the wind's mass flux through
the region that makes it: the loaded wind slows and its interaction region
grows to tens of lunar radii. A tenth to a half curve back into the air and
sputter out 1–10 molecules each. The solar wind's charge exchange with the dense
exosphere near the exobase makes ions that no optical shadow removes; in the
central estimate it costs 0.7–3 kg/s at the design point's exobases, which alone
uses a 1 kg/s budget.

**A lunar magnetosphere holds the wind off and keeps part of the ions.** With
the September magnets, of the ions made outside a heating-sized shadow, 25–75%
drain along closed field lines into the air and stay, 20–37% leave along the
open field lines over the poles, and recombination into fast atoms and the
wind-driven convection take the rest; the warmer the upper air, the farther its
exosphere reaches and the less the magnetosphere keeps. The magnetosphere also
keeps the solar wind, and its charge exchange, off the exosphere inside the
stand-off. For cold upper air a much weaker dipole does that: with a protected
radius of 4 lunar radii, about 3×10¹⁹ A·m² (a stand-off near 2.7 lunar radii, a
fiftieth of the September moment) brings the exosphere's loss at the design
point's cold exobases below 0.3 kg/s. Dipoles that hold the wind off only just
above the exobase fall in the range where hybrid simulations find a weak field
can raise ion escape (Egan et al. 2019), which this step leaves out.

**Magnets do not hold the fragments.** Sunlight breaks molecules apart in the
sunlit exosphere faster than it ionizes them, and the atoms are neutral. The
escaping fragments are 55–80% as much as the ions made, and with a magnetosphere
they are most of what the sunlit exosphere loses. Only the shadow removes them.

**With the exosphere counted, the shield has to reach beyond the exobase.**
Largest allowed UV transmission for the titania stack with the sunlit
exosphere's central loss added, by protected radius, with no magnetosphere and
with the September magnets. The ranges span the three treatments of the upper
air and quiet Sun to solar maximum, and a protected radius counts only where it
also covers the thermosphere's heating (a tenth of the whole heating, as in the
design study). Cases that cannot meet the budget at any transmission are
counted out of six:

| Budget (cycle time) | Magnetosphere | 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.025–0.14% (4 cannot) | 0.0008–0.17% (3 cannot) | 0.003–0.17% (3 cannot) | 0.003–0.17% (3 cannot) |
| | September | 0.016–0.33% (2 cannot) | 0.017–0.39% (1 cannot) | 0.057–0.43% (1 cannot) | 0.002–0.45% |
| 10 kg/s (10 billion years) | none | 0.006–0.41% (1 cannot) | 0.074–0.49% (1 cannot) | 0.015–0.54% | 0.035–0.57% |
| | September | 0.055–0.48% (1 cannot) | 0.007–0.55% | 0.031–0.60% | 0.047–0.64% |
| 100 kg/s (1 billion years) | none | 0.004–0.58% | 0.057–0.82% | 0.10–0.89% | 0.15–0.94% |
| | September | 0.004–0.58% | 0.077–0.88% | 0.12–0.93% | 0.16–0.96% |

The case that cannot meet 1 kg/s with the September magnets short of 10 lunar
radii is the warmest upper air at solar maximum, whose exobase sits near three
lunar radii even with a perfect shield. Without a magnetosphere the solar wind's
charge exchange keeps three or four of the six cases above 1 kg/s at any radius.
Behind the 200-nm edge no case meets 1 kg/s without a magnetosphere, and with
the September magnets only two of six do, from 6 lunar radii; 10 kg/s needs the
magnets and 6–8 lunar radii for most cases, and 100 kg/s holds from 4–5 lunar
radii.

**The escaping air makes a short-lived cloud around Earth.** Of the molecules
with the energy to pass L1, 65% leave the Moon's Hill sphere onto orbits around
Earth near the Moon's. Left alone, two-thirds of those would still orbit Earth
four months later, but sunlight destroys them first: N2 is ionized in 8–19 days
by this repository's own rate, and O2 is dissociated within days. For intact
lifetimes of 3 to 20 days, the steady cloud holds 150–1,400 tonnes for each kg/s
of escape. Of the escaped molecules, 91–99.5% end in space by sunlight, 0.5–4.5%
fall back on the Moon, up to 4% leave the Earth–Moon region, and almost none
reach Earth's atmosphere. With N2's lifetime the cloud reaches 140–155° ahead of
and behind the Moon along its orbit, spans 200,000 to about 1 million km from
Earth, and lies typically 40,000 km (up to 170,000 km) above or below the orbit's
plane: a thick, nearly closed ring. Near the Moon's orbit it holds 0.1–0.35
molecules per cm³ for each kg/s of escape, below the solar wind's 5 protons per
cm³. Collecting a tenth of the escape from it would take a collector 80,000–
170,000 km across, so recovering escaped air has to happen before it leaves.
About 13% of the cloud's destruction happens inside Earth's magnetotail, the
most that Earth's magnetosphere could hold for a while as ions; the rest ionizes
in the solar wind, which carries it off.

**Without a shield the air is lost fast.** The energy-limited loss with the full
ultraviolet reaching the air is 335–80,000 kg/s at quiet Sun (η 0.01–0.3,
absorption at 1.5–3 lunar radii). Deliveries slower than that would mostly be
lost, so the optical shield has to be working before the atmosphere accumulates.

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
the Earth–Moon region, and plasma. The band cross-sections in the limb are
averages. Atomic-oxygen escape, photochemical escape and hydrogen
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
