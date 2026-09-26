# Loss response to protection

How fast the finished 1.2 atm atmosphere loses mass as a function of its UV
transmission, the share of sunlight below 175 nm that the optical shield lets
through, and of how much of the solar wind reaches the air. It sizes the two
protection functions against a total loss budget for the protection design study
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).
[tides.py](tides.py) finds how much Earth's tide raises escape.
[absorption.py](absorption.py) finds how far out the optical shield must reach
and what the sunlit exosphere beyond it produces. All three are screening
models; their methods and limits are in their docstrings.

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # about 4 minutes; resumes if stopped
python -m atmosphere.loss_response.model        # then the loss response, about 10 s
python -m atmosphere.loss_response.absorption   # then the protected radius, about 10 s
python -m pytest atmosphere/tests/test_tides.py atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py
```

It needs the middle atmosphere's stored results and the WHI 2008 spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`). The product
is [results/loss_response.json](results/loss_response.json) (schema
`terluna.atmosphere.loss-response/2`), with the full sweep in
[results/euv_sweep.csv](results/euv_sweep.csv). Its field
`allowed_leak_fraction` is the allowed UV transmission. The tidal multipliers are
in [results/tidal_escape.json](results/tidal_escape.json) (schema
`terluna.atmosphere.tidal-escape/1`), and the protected radius and the
exosphere's ion production in
[results/absorption_radius.json](results/absorption_radius.json) (schema
`terluna.atmosphere.absorption-radius/1`).

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

**The shield may pass a few tenths of a percent.** Largest allowed UV
transmission, across the three upper-air treatments, with the solar-wind loss at
the low end of its range. Each budget's cycle time is the atmosphere's mass
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
36–470 times, so these figures hold only if most of those ions stay with the
Moon.

**The shield must reach to just under the exobase.** Ultraviolet grazing the
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
radii counted to half the Hill radius, or 17–26 counted to the whole. How much
of the production escapes depends on how the solar wind meets the ionized
exosphere, which none of these steps models.

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
its ion production is an upper bound on pickup loss. The band cross-sections in
the limb are averages. Atomic-oxygen escape, photochemical escape and hydrogen
from water are left out. The base temperatures come from a global-mean,
radiative-only middle atmosphere. The glow's brightness and its solar-cycle
factor (1.5) are assumptions. The solar-wind branch is a set of scalings with
assumed ranges; for comparison, Venus and Mars lose ions to the solar wind at
roughly a tenth of a kilogram a second, and scaling Titan's plasma-driven losses
to the solar wind at 1 AU gives a few kilograms a second (literature magnitudes,
not yet source-checked). The loss budget itself is the author's to choose.
