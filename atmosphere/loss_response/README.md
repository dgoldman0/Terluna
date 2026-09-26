# Loss response to protection

How fast the finished 1.2 atm atmosphere loses mass as a function of how much
light below 175 nm the optical shield lets through and how much of the solar
wind reaches the air. It sizes the two protection functions against a total loss
budget for the protection design study
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).
A second step, [absorption.py](absorption.py), finds how far out the optical
shield must reach and what the sunlit exosphere beyond it produces. Both are
screening models; their methods and limits are in their docstrings.

```sh
python -m atmosphere.loss_response.model        # writes results/ in about 10 s
python -m atmosphere.loss_response.absorption   # then the protected radius, about 10 s
python -m pytest atmosphere/tests/test_loss_response.py atmosphere/tests/test_absorption.py
```

It needs the middle atmosphere's stored results and the WHI 2008 spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`). The product
is [results/loss_response.json](results/loss_response.json) (schema
`terluna.atmosphere.loss-response/1`), with the full sweep in
[results/euv_sweep.csv](results/euv_sweep.csv). The protected radius and the
exosphere's ion production are in
[results/absorption_radius.json](results/absorption_radius.json) (schema
`terluna.atmosphere.absorption-radius/1`).

## How it works

The ultraviolet branch deposits heat above the 0.3 Pa base from three sources:
the leaked fraction f of sunlight below 175 nm, the titania film's own hard
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

The protected radius comes from the thermal column's profile at each budget's
allowed leak. Sunlight passing outside a shadow cylinder of radius R crosses
the collisional air below the exobase at a grazing angle and leaves part of its
ultraviolet there; the radius is where that deposit equals a tenth of the
allowed leak. Above the exobase the air is collisionless. Chamberlain's
exosphere, built from the exobase density, temperature and radius, gives the
molecules on ballistic flights, which rise and fall back, and the escaping
ones, which the Jeans loss already counts. Sunlight ionizes a molecule in about
19 days at quiet Sun and 8 days at solar maximum, and the step counts the ions
made from ballistic molecules in sunlight outside the shadow.

## Results

**The glow sets a floor the shield cannot remove.** With the leak at 10⁻⁶ the
glow alone raises the exobase by 31–42 K at quiet Sun and 47–61 K at solar
maximum. Behind the titania stack the floor loss stays small, 6×10⁻⁸–0.03 kg/s at
quiet Sun and up to 0.3 kg/s at solar maximum. Behind the 200-nm edge, whose upper
air starts warmer, it is 0.2–2.5 kg/s at quiet Sun and 1.1–9.4 kg/s at solar
maximum.

**The shield may leak a few tenths of a percent.** Largest allowed fraction of
sunlight below 175 nm, across the three upper-air treatments, with the
solar-wind loss at the low end of its range:

| Budget | Titania, quiet Sun | Titania, solar maximum | 200-nm edge, quiet Sun | 200-nm edge, solar maximum |
|---|---|---|---|---|
| 1 kg/s | 0.13–0.49% | 0.018–0.17% | 0.06% (none for the warmest upper air) | none |
| 10 kg/s | 0.27–0.72% | 0.075–0.24% | 0.08–0.21% | 0.002–0.05% |
| 100 kg/s | 0.74–1.2% | 0.26–0.45% | 0.50–0.65% | 0.17–0.23% |

Solar maximum with the warmest upper air is the binding case. The ozone-forming
edge cannot meet a 1 kg/s budget at solar maximum even with a perfect shield,
and its figures leave out atomic-oxygen escape, which adds about half again.

**The solar-wind branch as screened.** Its range is 0.03 kg/s (low),
0.19 kg/s (central) and 2.1 kg/s (high); at the high end sputtering by
returning pickup ions dominates. With pickup capped at 0.29 kg/s it barely
changes the allowed leak for budgets of 10 kg/s and more, and at 1 kg/s its high
end alone exceeds the budget. The ion production below exceeds that cap 40–600
times, so these figures hold only if most of those ions stay with the Moon.

**The shield must reach to just under the exobase.** Ultraviolet grazing the
thermosphere outside the shadow deposits a tenth of the allowed leak when the
protected radius is within about 0.1 lunar radii of the exobase: 2.7–3.2 lunar
radii for 1 kg/s, 3.3–4.0 for 10 kg/s and 4.2–5.2 for 100 kg/s. The aperture at
78,000 km then has 0.84–1.1, 1.2–1.7 and 1.9–2.8 times the September design's
area. The September three lunar radii cover four of the eight 1 kg/s cases and
none at 10 kg/s or more. The extreme ultraviolet itself is absorbed lower, with
unit optical depth at grazing incidence at 2.1–3.0 lunar radii, and the thin
air above that still takes up too much.

**The sunlit exosphere beyond the shadow makes ions faster than the budget.**
The exosphere is large because the exobase sits at three to five lunar radii,
where gravity is weak. Outside a shadow sized for heating, its ballistic
molecules become ions at 13–36 kg/s for the 1 kg/s budget (13–42 times the Jeans
loss), 27–80 kg/s for 10 kg/s and 59–172 kg/s for 100 kg/s, counted to a third of
the Hill radius. Counting to the whole Hill radius adds about 5% at 1 kg/s and 68–77% at
100 kg/s. A shadow that held production to a tenth of the budget would need a
radius of 11–15 lunar radii counted to half the Hill radius, or 22–29 counted to
the whole. How much of the production escapes depends on how the solar wind
meets the ionized exosphere, which neither step models.

**Without a shield the air is lost fast.** The energy-limited loss with the full
ultraviolet reaching the air is 335–80,000 kg/s at quiet Sun (η 0.01–0.3,
absorption at 1.5–3 lunar radii). Deliveries slower than that would mostly be
lost, so the optical shield has to be working before the atmosphere accumulates.

## Limits

The thermal column is molecular only, with no infrared cooling and no tides; it
flags exobases beyond 3.5 lunar radii for tidal review. The feasibility
baseline's tidal barrier factor, a restricted-three-body sensitivity,
multiplies the molecular escape at the allowed leaks by 2.5–4.2, so the allowed
leaks are upper limits until tides are included. The exosphere is two-body and
collisionless, with no satellite particles and no depletion along each flight;
its ion production is an upper bound on pickup loss. The band cross-sections in
the limb are averages. Atomic-oxygen escape,
photochemical escape and hydrogen from water are left out. The base
temperatures come from a global-mean, radiative-only middle atmosphere. The
glow's brightness and its solar-cycle factor (1.5) are assumptions. The
solar-wind branch is a set of scalings with assumed ranges; for comparison,
Venus and Mars lose ions to the solar wind at roughly a tenth of a kilogram a
second, and scaling Titan's plasma-driven losses to the solar wind at 1 AU gives
a few kilograms a second (literature magnitudes, not yet source-checked). The
loss budget itself is the author's to choose.
