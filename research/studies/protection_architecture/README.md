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
| [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | Loss of the finished atmosphere against UV transmission and solar-wind exposure: the physics that sizes both functions |
| [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Upper-air temperature, ozone and surface ultraviolet behind each shield |
| [climate/gcm](../../../climate/gcm/README.md) | The climate's response to dimming: about 1.9 K per 1% of sunlight, for an even cut |
| [protection](../../../protection/README.md) | Hardware: the film's optics and spectral products, the UV transmission a swarm of cells can hold, the September reference design and its report, magnets, and the module catalogue |
| [engineering](../../../engineering/README.md) | Industry: the reconstructed 19 September architecture, and the network model's vocabulary with a supply ledger for the protection hardware |
| [illumination](../../../illumination/) | The shield's effects on Earth and on the sky, to build |

## Interfaces

| Product | Producer | Read by | Status |
|---|---|---|---|
| Shield transmission T(λ), 200 nm–5 µm | `protection/spectra/shield_transmission.json` | atmosphere, climate | Exists |
| Film transmission, 0.1–210 nm | `protection/spectra/stack_short_wave.json` | atmosphere | Exists |
| Loss response | `atmosphere/loss_response/results/loss_response.json` | this study | Exists, as a screening model |
| Protected radius and exospheric ion production | `atmosphere/loss_response/results/absorption_radius.json` | this study | Exists, as a screening model |
| UV transmission a swarm of cells can hold, by design level | `protection/transmission/results/uv_transmission.json` | this study | Exists, from assumed design choices |
| Design point: loss, cycle time, protected radius and ions at the swarm's transmission | `research/studies/protection_architecture/results/design_point.json` | this study | Exists, as a screening model |
| Climate response to dimming | climate/gcm runs `A28`, `A28_dim3`, `A28_dim5`, `A28_dim8` | this study | Exists for an even cut |
| Module catalogue: each unit's mass, power, materials, lifetime and failure behaviour | `protection/modules/catalogue.json` | engineering's network model | Exists for the September reference design |
| Supply ledger: replacement, propellant and fresh material over time | `engineering/network/results/protection_supply_ledger.json` | this study | Exists for the September reference design |
| Effects on Earth: reflections, crossings of the full Moon and the Sun | illumination | this study | To build |

## Design point

`python -m research.studies.protection_architecture.run` (about 40 s, with
`OPENBLAS_NUM_THREADS=1`) takes the UV transmission a swarm can hold at three
levels of design choices and evaluates the loss response and its absorption
step there. The swarm holds 3–4×10⁻⁵ with tight choices, about 2×10⁻⁴ with
standard ones and 2.5×10⁻³ with relaxed ones, most of it from failed cells
waiting to be covered. At the tight and standard levels the sky's Lyman-alpha
glow supplies 76–98% of the upper air's heat. Behind the titania stack at the
standard level the loss is 0.03–4.6 kg/s with the solar-wind range and Earth's
tide, a cycle time of 22 billion years or longer, and the shield must reach
1.9–3.1 lunar radii. Behind the 200-nm edge the loss is 0.8–40 kg/s and the
radius 3.0–4.1 lunar radii. The ions the exosphere makes outside the shield, 2–38 kg/s behind the
titania stack, are the largest open term. The requirements hold the details.

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

The loss budget is designed across 1–100 kg/s, atmospheric cycle times of
about 100 billion to 1 billion years, until the designs show which rate works
best (the author's decision of 2026-09-26); requirements.md shows what each
budget asks. Earth's tide is in the loss response and roughly doubles molecular
escape. The design point puts a swarm with standard choices at the bottom of the
range or below it for the modelled channels, so the loss now turns on the escape
of ions made in the sunlit exosphere, which also decides whether
charged-particle protection is needed. Next comes an estimate of how many of
those ions escape, from published results for bodies in flowing plasma and from
how far out the September design's magnets would hold off the solar wind; then
the trade of how the shield is held, against requirement S6, with the module
catalogue and supply ledger scoring each option; and the assessment of effects
on Earth.
