# Protection architecture

The design study for the Open Moon's protection: the optical shield that stops
extreme ultraviolet and sets the climate, charged-particle protection against
solar-wind stripping, and later protected destinations and Earth–Moon corridors.
It holds the requirements, the interfaces between lanes and the architecture
trade. Each lane keeps its own models, and results cross between them as data
products. The scope is the Moon and its atmosphere.

[requirements.md](requirements.md) sets out what the protection must achieve,
with sources and status, and what each loss budget asks of it.

## Where the parts live

| Lane | Holds for this design |
|---|---|
| This folder | Requirements, interfaces, the architecture trade and its record |
| [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | Loss of the finished atmosphere against ultraviolet leak and solar-wind exposure: the physics that sizes both functions |
| [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Upper-air temperature, ozone and surface ultraviolet behind each shield |
| [climate/gcm](../../../climate/gcm/README.md) | The climate's response to dimming: about 1.9 K per 1% of sunlight, for an even cut |
| [protection](../../../protection/README.md) | Hardware: the film's optics and spectral products, the September reference design and its report, magnets, and the module catalogue |
| [engineering](../../../engineering/README.md) | Industry: the reconstructed 19 September architecture, and the network model's vocabulary with a supply ledger for the protection hardware |
| [illumination](../../../illumination/) | The shield's effects on Earth and on the sky, to build |

## Interfaces

| Product | Producer | Read by | Status |
|---|---|---|---|
| Shield transmission T(λ), 200 nm–5 µm | `protection/spectra/shield_transmission.json` | atmosphere, climate | Exists |
| Film transmission, 0.1–210 nm | `protection/spectra/stack_short_wave.json` | atmosphere | Exists |
| Loss response | `atmosphere/loss_response/results/loss_response.json` | this study | Exists, as a screening model |
| Protected radius and exospheric ion production | `atmosphere/loss_response/results/absorption_radius.json` | this study | Exists, as a screening model |
| Climate response to dimming | climate/gcm runs `A28`, `A28_dim3`, `A28_dim5`, `A28_dim8` | this study | Exists for an even cut |
| Module catalogue: each unit's mass, power, materials, lifetime and failure behaviour | `protection/modules/catalogue.json` | engineering's network model | Exists for the September reference design |
| Supply ledger: replacement, propellant and fresh material over time | `engineering/network/results/protection_supply_ledger.json` | this study | Exists for the September reference design |
| Effects on Earth: reflections, crossings of the full Moon and the Sun | illumination | this study | To build |

## Working rules

Numbers come from this repository's models. The reconstructions in
`protection/reference/historical/` and `engineering/reference/industrial_architecture/`
were written on 2026-09-26 from conversation records; they are sources of
ideas, and a value from them enters a calculation only after an independent
check. Protection hardware is a formation of replaceable units, so ageing and
damage appear as replacement rates in the supply ledger. Decisions live in
[research/decisions.md](../../decisions.md) and are reopened only for good
reasons, stated.

## Next

The loss budget is designed across 1–100 kg/s until the designs show which
rate works best (the author's decision of 2026-09-26); requirements.md shows
what each budget asks. The protected radius is computed: the shield reaches to
just under the exobase, 2.7–5.2 lunar radii across the range. Two losses the
loss response leaves out are large enough to change the allowed leaks, Earth's
tidal lowering of the escape barrier and the escape of ions made in the sunlit
exosphere, and the second decides whether charged-particle protection is needed
at every budget. Next come those two losses in the loss response, with a plasma
estimate of the ions' escape; the trade of how the shield is held, against
requirement S6, with the module catalogue and supply ledger scoring each
option; and the assessment of effects on Earth.
