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
| [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | Loss of the finished atmosphere against UV transmission, protected radius and solar-wind exposure, with and without a lunar magnetosphere: the physics that sizes both functions |
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
| The sunlit exosphere's losses (ions, fragments of broken molecules, the solar wind's charge exchange) with and without a lunar magnetosphere, and the allowed transmissions with them | `atmosphere/loss_response/results/exosphere_loss.json` | this study | Exists, as a screening model |
| UV transmission a swarm of cells can hold, by design level | `protection/transmission/results/uv_transmission.json` | this study | Exists, from assumed design choices |
| Design point: loss with the exosphere, cycle time, protected radius and magnetic moment at the swarm's transmission | `research/studies/protection_architecture/results/design_point.json` | this study | Exists, as a screening model |
| Climate response to dimming | climate/gcm runs `A28`, `A28_dim3`, `A28_dim5`, `A28_dim8` | this study | Exists for an even cut |
| Module catalogue: each unit's mass, power, materials, lifetime and failure behaviour | `protection/modules/catalogue.json` | engineering's network model | Exists for the September reference design |
| Supply ledger: replacement, propellant and fresh material over time | `engineering/network/results/protection_supply_ledger.json` | this study | Exists for the September reference design |
| Effects on Earth: reflections, crossings of the full Moon and the Sun | illumination | this study | To build |

## Design point

`python -m research.studies.protection_architecture.run` (about 15 minutes,
with `OPENBLAS_NUM_THREADS=1` and three single-threaded processes for the solar
cycle) takes the UV transmission a swarm can hold at three levels of design
choices and evaluates the loss response and its absorption and exosphere steps
there, at the ring fleet's protected radius of 4 lunar radii. It uses the traced
count and the heat placed where the tracing puts it, both as the author approved
on 7 October, at quiet Sun, at solar maximum (FISM2's year around cycle 21's
maximum, with cycle 19's year and the escape model's 2.5 as stress cases) and
over solar cycles 23 and 24 year by year from FISM2's daily record. The swarm holds
3–4×10⁻⁵ with tight choices, about 2×10⁻⁴ with standard ones and 2.5×10⁻³ with
relaxed ones, most of it from failed cells waiting to be covered. Traced along
its slant paths, light through gaps heats the upper air 2.5–5 times a count over
the disk; the sky's Lyman-alpha glow supplies 75–93% of the heat at the tight
level and 49–70% at the standard level, but it lands just above the base and
barely warms the exobase. Averaged over the solar cycle, the standard level at
4 lunar radii loses 0.68–0.97 kg/s with no magnetosphere with collisional upper
air and in LTE (144 to 101 billion years), almost all of it the solar wind's;
0.064–0.093 kg/s if the ring fleet's wake holds the wind off; and
10⁻⁵–0.003 kg/s with the September magnets. With all near-infrared heating it
loses 3.1–3.5 kg/s with no magnetosphere (32 to 28 billion years), 2.3–2.7 kg/s
in the wake and 0.75–0.95 kg/s with the magnets (131 to 104 billion years). At
the relaxed level the air runs away in most cases. The requirements hold the
details.

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

The [2026-10-05 primary magnetic comparison](../solar_shield_array/magnetic_architecture.md)
now carries this study's stored atmospheric requirements into finite regional
coil and upstream-source budgets. It read the design point from before the
traced count of 7 October, and its rerun waits on the closure trajectories it
needs, which this machine lacks. Larger regional paths merit development;
the quiet-wind lower-moment cases need a plasma/storm check before installation
is sized. The four-radius optical footprint and magnetic requirement are kept
separate. No loss budget or primary architecture has been selected.

The loss budget is designed across 1–100 kg/s, atmospheric cycle times of
about 100 billion to 1 billion years, until the designs show which rate works
best (the author's decision of 2026-09-26); requirements.md shows what each
budget asks. Earth's tide and the sunlit exosphere's losses are now in the loss
response. The exosphere step shows that the protected radius has to reach
beyond the exobase, because the sunlit exosphere just outside it loses ions and
the fragments of broken molecules, and that a 1 kg/s budget needs
charged-particle protection against the solar wind's charge exchange; for the
cooler upper air a dipole far weaker than the September magnets would do. Next
comes the trade of how the shield is held, against requirement S6, with the
module catalogue and supply ledger scoring each option, now across protected
radii of 3–10 lunar radii and with the magnets' mass and power; then a check of
the ion fates against published plasma simulations of small magnetospheres and
loaded exospheres, since the design now leans on them; the assessment of
effects on Earth; and the swarm's own assumptions (formation precision, failure
and cover times, hole size).
