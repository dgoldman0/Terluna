# Protection requirements

What the Open Moon's protection systems must achieve, where each requirement
comes from, and its status. Decisions are the author's
([research/decisions.md](../../decisions.md)). Derived requirements follow from
them through the calculation named. Open items wait on the author or on a
calculation. Every number here comes from this repository's models; the
reconstructed records in `protection/reference/historical/` and
`engineering/reference/industrial_architecture/` are sources of ideas only.

## The functions

The **optical shield** stops extreme and far ultraviolet before it heats and
ionizes the upper air, which drives escape, and dims the sunlight that reaches
the surface to set the climate. **Charged-particle protection** limits the
solar wind's stripping of the upper air by ion pickup and sputtering.
**Protected destinations and Earth–Moon corridors** come later, as their own
requirement for radiation dose along actual routes, secondary to protecting the
lunar atmosphere and surface. The Earth–Moon–L1–L2 megastructure concept is
outside this design.

## Environment to hold

| ID | Requirement | Source | Status |
|---|---|---|---|
| E1 | Surface pressure 1.2 atm (121,590 Pa), with 17.5% oxygen in the dry air | decisions.md (Sep 4); the oxygen share is the reports' starting choice | Fixed |
| E2 | A climate near 294 K: the shield passes 5% less sunlight than the titania stack at every wavelength, 1,173 W/m² at the top of the air | decisions.md; [climate/gcm](../../../climate/gcm/README.md) | Fixed; the adjustable range is open |
| E3 | An exobase near 250 K preferred, with 260 K a candidate limit | decisions.md (Sep 4) | Fixed as a preference; the loss budget R1 is the binding measure |
| E4 | A surface radiation dose of at most 0.027 mSv/day | decisions.md (Sep 4) | Fixed; the air column (about 7,900 g/cm²) is expected to provide it, pending a particle-transport calculation |

## Retention

| ID | Requirement | Source | Status |
|---|---|---|---|
| R1 | A total loss budget: the long-term average loss the protection must hold. 1 kg/s loses 1% of the atmosphere per 10⁹ years, 10 kg/s about 10%, 100 kg/s about one atmosphere | September feasibility report, section 5; decisions.md (the author, 2026-09-26) | Designed across 1–100 kg/s until the designs show which rate works best |
| R2 | The budget is shared among ultraviolet-driven escape, solar-wind stripping, and the channels not yet in the loss response: Earth's tidal lowering of the escape barrier, the escape of ions made in the sunlit exosphere, atomic oxygen, photochemical escape and hydrogen from water | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | Allocation open |

What each budget asks of the protection, from the loss response. The ranges span
the three treatments of the upper air and quiet Sun to solar maximum; solar
maximum with the warmest upper air is the binding case.

| Budget | Titania stack: allowed leak | 200-nm edge: allowed leak | Exobase at that leak | Protected radius for heating | Ions made outside that radius |
|---|---|---|---|---|---|
| 1 kg/s | 0.018–0.49% | none at solar maximum | 2.8–3.2 lunar radii | 2.7–3.2 lunar radii | 13–36 kg/s |
| 10 kg/s | 0.075–0.72% | 0.002–0.21% | 3.3–4.0 lunar radii | 3.3–4.0 lunar radii | 27–80 kg/s |
| 100 kg/s | 0.26–1.2% | 0.17–0.65% | 4.3–5.2 lunar radii | 4.2–5.2 lunar radii | 59–172 kg/s |

A tighter budget needs a tighter shield but a smaller one. The upper air expands
with the heat the allowed leak deposits, and the shield has to cover it to just
under the exobase; the September design's three lunar radii suffice for four of
the eight 1 kg/s cases. The ozone-forming edge makes every budget harder and
cannot meet 1 kg/s at solar maximum, because its upper air starts warmer and the
sky's Lyman-alpha glow alone then drives 1–9 kg/s. Without the optical shield
the loss is 300–80,000 kg/s.

The allowed leaks are upper limits, for two reasons. Earth's tide lowers the
escape barrier: the feasibility baseline's barrier factor multiplies the
molecular escape at these leaks by 2.5–4.2. And the exosphere beyond the shadow
makes ions at 0.5–42 times the Jeans loss, which the solar wind can carry off.
The loss response caps that pickup at 0.29 kg/s; how many of the ions escape
needs a plasma model. If most escape, an optical shield of practical size
cannot hold 1 or 10 kg/s, and holds 100 kg/s only with about the leak now
allowed for 10 kg/s. Holding the ions to a tenth of the budget with the shadow
alone would take a radius of 11–29 lunar radii, so the choice then lies with
charged-particle protection.

## Optical shield

| ID | Requirement | Source | Status |
|---|---|---|---|
| O1 | The leak of sunlight below 175 nm, averaged over time and over the protected region, stays within the allowed fraction for R1. Gaps, pinholes, edges, off-normal light and outages all count | loss response | Derived; the value follows R1 and O3 |
| O2 | The protected radius reaches to just under the exobase of the design case, where ultraviolet deposited in the thermosphere outside it is a tenth of the allowed leak: 2.7–3.2 lunar radii for 1 kg/s, 3.3–4.0 for 10 kg/s and 4.2–5.2 for 100 kg/s. At 78,000 km the aperture then has 0.84–1.1, 1.2–1.7 and 1.9–2.8 times the September design's area | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), absorption step; [protection/report.md](../../../protection/report.md), section 2 | Derived |
| O3 | The ultraviolet cut-off: the titania stack (no ozone, surface UV index 0.1) or an ozone-forming edge near 200 nm | [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Open: the author's choice |
| O4 | The window passes the climate's sunlight (E2), evenly unless a pattern is chosen and run through the climate model. At 78,000 km an even cut needs the dimmer to reach 2,100 km from the axis, or 2,400 km to include the twilight band | E2; this study | Derived |
| O5 | Light the shield rejects never reaches Earth's night side: specular reflections are steered at least about 5° off the plane in which Earth moves as seen from the shield, and Sun-facing surfaces keep diffuse reflectance near 0.1% | this study | Derived; to confirm in illumination |
| O6 | The optical shield works before the atmosphere accumulates | feasibility report, section 12; loss response | Derived |
| O7 | Effects on Earth are assessed and bounded: reflections, and the shield crossing the full Moon and the Sun as seen from Earth | protection/report.md, section 6 | Required assessment |

## Charged-particle protection

| ID | Requirement | Source | Status |
|---|---|---|---|
| C1 | Solar-wind loss stays within its share of R1 | loss response, with its absorption step | Open. With pickup capped at 0.29 kg/s it is needed only near 1 kg/s; if most of the 13–172 kg/s of ions made outside the shadow escape, it is needed at every budget. A plasma model decides between them |
| C2 | If needed, it grows in stages with the atmosphere: seed stations with local cavities, then more nodes whose cavities merge into a protected volume | the 2025 Lunashield staging, known from a 2026-09-26 reconstruction | Concept |
| C3 | No superconducting planetary ring | decisions.md (Sep 4) | Fixed |

## The system

| ID | Requirement | Source | Status |
|---|---|---|---|
| S1 | Formations of replaceable units, never single megastructures; ageing and damage are replacement rates | decisions.md (the author, 2026-09-26) | Fixed |
| S2 | Operation for at least 10⁹ years through replacement, within a supply ledger that counts propellant, replacement units and fresh material | decisions.md (Sep 4) | Fixed; the ledger exists for the September reference design ([engineering/network](../../../engineering/network/README.md)) |
| S3 | Construction within 500 years of 2026, with protection working first (O6) | decisions.md (Sep 4) | Fixed |
| S4 | All protection logistics follow the six traffic-safety rules | decisions.md (Sep 19, confirmed 2026-09-26) | Fixed |
| S5 | The recorded protection architecture: a solar-filter complex with an industrial hub, magnetic protection, destinations and corridors later, and Earth–Sun L1/L2 hubs for freight, power and industry | decisions.md (Sep 9) | Fixed |
| S6 | Holding the optical shield costs less material over its life than the atmosphere it saves. The September reference fails this: 2.8×10⁵ kg/s of propellant against 300–80,000 kg/s of loss without a shield | this study | Derived |

## Open decisions for the author

The loss budget (R1), designed across 1–100 kg/s for now; the ultraviolet
cut-off (O3); the dimmer's form (an even cut or one taken from the near
infrared, fixed or adjustable) and its range;
the use of the ring outside the window (transparent film or power collector);
and whether and when to build charged-particle protection (C1).
