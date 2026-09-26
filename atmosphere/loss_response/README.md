# Loss response to protection

How fast the finished 1.2 atm atmosphere loses mass as a function of how much
light below 175 nm the optical shield lets through and how much of the solar
wind reaches the air. It sizes the two protection functions against a total loss
budget for the protection design study
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).
It is a screening model; the method and its limits are in the docstring of
[model.py](model.py).

```sh
python -m atmosphere.loss_response.model     # writes results/ in about 10 s
python -m pytest atmosphere/tests/test_loss_response.py
```

It needs the middle atmosphere's stored results and the WHI 2008 spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`). The product
is [results/loss_response.json](results/loss_response.json) (schema
`terluna.atmosphere.loss-response/1`), with the full sweep in
[results/euv_sweep.csv](results/euv_sweep.csv).

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

**The solar wind matters only for tight budgets.** Its screening range is
0.03 kg/s (low), 0.19 kg/s (central) and 2.1 kg/s (high); at the high end
sputtering by returning pickup ions dominates. For budgets of 10 kg/s and more
it barely changes the allowed leak. At 1 kg/s the high end alone exceeds the
budget, so a solar-wind system would be needed there.

**Without a shield the air is lost fast.** The energy-limited loss with the full
ultraviolet reaching the air is 335–80,000 kg/s at quiet Sun (η 0.01–0.3,
absorption at 1.5–3 lunar radii). Deliveries slower than that would mostly be
lost, so the optical shield has to be working before the atmosphere accumulates.

## Limits

The thermal column is molecular only, with no infrared cooling and no tides; it
flags exobases beyond 3.5 lunar radii for tidal review. Atomic-oxygen escape,
photochemical escape and hydrogen from water are left out. The base
temperatures come from a global-mean, radiative-only middle atmosphere. The
glow's brightness and its solar-cycle factor (1.5) are assumptions. The
solar-wind branch is a set of scalings with assumed ranges; for comparison,
Venus and Mars lose ions to the solar wind at roughly a tenth of a kilogram a
second, and scaling Titan's plasma-driven losses to the solar wind at 1 AU gives
a few kilograms a second (literature magnitudes, not yet source-checked). The
loss budget itself is the author's to choose.
