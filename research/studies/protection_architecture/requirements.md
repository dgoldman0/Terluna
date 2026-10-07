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
| R1 | A total loss budget: the long-term average rate at which the atmosphere escapes, which resupply must match. Each budget has an atmospheric cycle time, the atmosphere's mass (3.1×10¹⁸ kg) over the loss rate: the time the loss and its resupply take to replace the whole atmosphere. 1 kg/s is a cycle of about 100 billion years and a resupply of about 32,000 tonnes a year; 10 kg/s, 10 billion years; 100 kg/s, 1 billion years | September feasibility report, section 5; decisions.md (the author, 2026-09-26); the mass from the feasibility baseline's hydrostatic column (`research/baselines/feasibility/reference.json`) | Designed across 1–100 kg/s until the designs show which rate works best |
| R2 | The budget is shared among ultraviolet-driven escape (with Earth's tide), the sunlit exosphere's losses (the ions it makes, the fragments of broken molecules and the solar wind's charge exchange), sputtering, and the channels not yet in the loss response: atomic oxygen, photochemistry below the exobase and hydrogen from water | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), with its exosphere step | Allocation open |

What each budget asks of the protection, from the loss response with its
exosphere step. The allowed UV transmission is the largest share of sunlight
below 175 nm that may reach the protected region once the sunlit exosphere's
losses (the ions it makes, the fragments of broken molecules and the solar
wind's charge exchange, central estimate) are added to ultraviolet-driven
escape. It depends on how far the protected region reaches and on whether a
magnetosphere holds the solar wind off. For the titania stack, with the ranges
spanning the three treatments of the upper air and quiet Sun to solar maximum,
and counting the cases of those six that cannot meet the budget at any
transmission:

| Budget (cycle time) | Magnetosphere | Protected radius 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.025–0.14% (4 cannot) | 0.0008–0.17% (3 cannot) | 0.003–0.17% (3 cannot) | 0.003–0.17% (3 cannot) |
| | September magnets | 0.016–0.33% (2 cannot) | 0.017–0.39% (1 cannot) | 0.057–0.43% (1 cannot) | 0.002–0.45% |
| 10 kg/s (10 billion years) | none | 0.006–0.41% (1 cannot) | 0.074–0.49% (1 cannot) | 0.015–0.54% | 0.035–0.57% |
| | September magnets | 0.055–0.48% (1 cannot) | 0.007–0.55% | 0.031–0.60% | 0.047–0.64% |
| 100 kg/s (1 billion years) | none | 0.004–0.58% | 0.057–0.82% | 0.10–0.89% | 0.15–0.94% |
| | September magnets | 0.004–0.58% | 0.077–0.88% | 0.12–0.93% | 0.16–0.96% |

Solar maximum with the warmest upper air is the binding case: its exobase sits
near three lunar radii even with a perfect shield, and it meets 1 kg/s only with
the September magnets and a protected radius of 10 lunar radii. Without a
magnetosphere the solar wind's charge exchange keeps three or four of the six
cases above 1 kg/s at any radius. The ozone-forming edge makes every budget
harder: without a magnetosphere it cannot meet 1 kg/s, and with the September
magnets two of six cases do, from 6 lunar radii; 10 kg/s needs the magnets and
6–8 lunar radii for most cases, and 100 kg/s holds from 4–5 lunar radii.
Without the optical shield the loss is 300–80,000 kg/s.

The protected radius has to reach beyond the exobase. At the transmissions
ultraviolet-driven escape alone would allow (0.004–0.97% for the titania stack),
where that escape already fills the budget, a shadow sized only for the
thermosphere's heating, just under the exobase, leaves a sunlit exosphere that
adds 25–303 kg/s with no magnetosphere and 11–193 kg/s with the September
magnets; most of that comes from its first few scale heights above the exobase.
Earth's tide is included throughout and multiplies molecular escape by 1.8–2.6.

**The design point.** A formation of replaceable cells can hold a UV
transmission of 3–4×10⁻⁵ with tight design choices, about 2×10⁻⁴ with standard
ones and 2.5×10⁻³ with relaxed ones
([protection/transmission](../../../protection/transmission/README.md)); covering
failed cells quickly is the main lever (O8). The sky's Lyman-alpha glow heats
the upper air as much as a transmission of about 10⁻³, so at the tight and
standard levels it supplies 76–98% of the heat and the transmission barely
moves the loss. Behind the titania stack at the standard level the
ultraviolet-driven loss, with Earth's tide, is 7×10⁻⁷–2.5 kg/s, and the exobase
sits at 1.9–3.2 lunar radii (179–267 K). With the protected radius covering only
the heating (1.9–3.1 lunar radii), the sunlit exosphere beyond it raises the
central total to 5.7–89 kg/s with no magnetosphere and 2.1–47 kg/s with the
September magnets, cycle times of 1.1–17 and 2.1–48 billion years. With the
protected radius at 4 lunar radii and the September magnets the central total
falls to 0.00004–0.2 kg/s in the four cooler cases (cycle times of 500 billion
years or longer), about 1 kg/s with all near-infrared heating at quiet Sun (100
billion years) and 16–18 kg/s with it at solar maximum (about 6 billion years).
Without a magnetosphere the same radius leaves 0.7–1.7 kg/s in the cooler cases
(60–140 billion years). Behind the 200-nm edge the central total at 4 lunar
radii with the September magnets is 3.7–49 kg/s (27 to 2 billion years), where
that radius covers the heating
([results/design_point.json](results/design_point.json), from
[run.py](run.py)). These figures count a gap's light over the disk; traced
along its slant paths it heats the air 2.5–5 times as much (O1). The author
adopted the traced count and FISM2's X-ray factors on 7 October, and these
figures are to be rederived with them.

## Optical shield

| ID | Requirement | Source | Status |
|---|---|---|---|
| O1 | The UV transmission, the share of sunlight below 175 nm that reaches the protected region, averaged over time and over that region, stays within the allowed fraction for R1. Gaps, pinholes, edges, off-normal light and outages all count | loss response with its exosphere step; [protection/transmission](../../../protection/transmission/README.md) | Derived. A swarm can hold 3×10⁻⁵–2×10⁻⁴ (O8). With a protected radius of 4 lunar radii and the September magnets the standard level meets 1 kg/s in the four cooler titania cases and sits at about 1 kg/s with all near-infrared heating at quiet Sun. With that heating at solar maximum the standard level meets 10 kg/s from about 6 lunar radii, and 1 kg/s would take 10 lunar radii and 2×10⁻⁵, below even the tight level. Below about 10⁻³ the sky's glow heats the upper air more than the transmission does. These levels count a gap's light over the disk. Traced along its slant paths above the limb ([limb heat](../../../atmosphere/middle_atmosphere/README.md#slant-paths-above-the-limb-2026-10-07), 7 October) it heats the air 2.5–5 times as much, the glow's crossover falls to about 5×10⁻⁴, and at the standard level and solar maximum the molecular loss alone is 30–36 kg/s with all near-infrared heating while the cooler titania cases stay under 1 kg/s ([the X-ray check](../solar_shield_array/integrated_comparison.md#the-x-ray-check)). On 7 October the author decided that O1 counts gaps by the traced heat and that the X-rays take FISM2's measured solar-cycle factors ([decisions.md](../../decisions.md)); the levels above are to be rederived with both |
| O2 | The protected radius covers the thermosphere's heating, which takes it to just under the exobase, and reaches far enough beyond the exobase that the sunlit exosphere's losses stay within their share of R1. At the design point, with the September magnets, the exosphere's central loss falls to 1 kg/s at 2.0–3.1 lunar radii in the four cooler cases, 3.9 with all near-infrared heating at quiet Sun and 8.5–8.7 at solar maximum; to 0.1 kg/s at 2.2–4.6, 7.0–7.1 and 11 lunar radii. At 78,000 km apertures of 4, 6 and 10 lunar radii have 1.7, 3.7 and 10 times the September design's area | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), absorption and exosphere steps; [protection/report.md](../../../protection/report.md), section 2 | Derived |
| O3 | The ultraviolet cut-off: the titania stack (no ozone, surface UV index 0.1) or an ozone-forming edge near 200 nm | [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Open: the author's choice |
| O4 | The window passes the climate's sunlight (E2), evenly unless a pattern is chosen and run through the climate model. At 78,000 km an even cut needs the dimmer to reach 2,100 km from the axis, or 2,400 km to include the twilight band | E2; this study | Derived |
| O5 | Light the shield rejects never reaches Earth's night side: specular reflections are steered at least about 5° off the plane in which Earth moves as seen from the shield, and Sun-facing surfaces keep diffuse reflectance near 0.1% | this study | Derived; to confirm in illumination |
| O6 | The optical shield works before the atmosphere accumulates | feasibility report, section 12; loss response | Derived |
| O7 | Effects on Earth are assessed and bounded: reflections, and the shield crossing the full Moon and the Sun as seen from Earth | protection/report.md, section 6 | Required assessment |
| O8 | Neighbouring cells overlap by several times their relative position error (50 m in the design levels, 1% more film); cells are patched or replaced at least every 10–20 years against micrometeoroid holes; a failed cell is covered within a day to a week by spares or overlays kept in the formation | [protection/transmission](../../../protection/transmission/README.md) | Derived; cell failure rates, the time to cover and the film's hole size are assumptions to test |

## Charged-particle protection

| ID | Requirement | Source | Status |
|---|---|---|---|
| C1 | Solar-wind loss stays within its share of R1 | loss response, with its absorption and exosphere steps | Derived from screening, to confirm with a plasma model. Without a magnetosphere the solar wind carries off the ions the sunlit exosphere makes, and its charge exchange with the dense exosphere near the exobase adds 0.7–3 kg/s (central) at the design point that no optical shadow removes, so 1 kg/s needs charged-particle protection. A lunar dipole whose stand-off clears the dense exosphere removes the charge exchange and keeps a quarter to three quarters of the ions it encloses: about 3×10¹⁹ A·m² does it for the design point's cooler upper air, and the September magnets (1.5×10²¹ A·m², stand-off 10 lunar radii) with margin. Stand-offs just above the exobase fall where simulations find a weak field can raise escape. Magnets do not hold the neutral fragments, which set the protected radius (O2) |
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
| S7 | Outflows stay out of the escape region. Everything the protection systems release or redirect (propellant exhaust with its broad edges and charge-exchange ions, unionized propellant, rejected light, shed material and failed units) is traced to where it ends up, including later passages near Earth, and the energy and mass it delivers inside the protected region stay within a stated share of R1. A held screen's exhaust leaves into the Moon's hemisphere. By the September report's energy diagnostic (a tenth of the deposited energy doing escape work against 0.94 MJ/kg at three lunar radii), each kg/s of budget allows about 9.4 MW deposited: 4×10⁻⁷ of the zoned held screen's 23 TW of jet power for 1 kg/s and 4×10⁻⁶ for 10 kg/s. Returning ions sputter 1–10 molecules each | the author, 2026-10-06; [protection/report.md](../../../protection/report.md), section 8; [atmosphere/loss_response](../../../atmosphere/loss_response/README.md) | Adopted by the author on 2026-10-07 ([decisions.md](../../decisions.md)); the share of R1 is open. As designed, the zoned held screen exceeds it through its direct plumes and its unionized gas. The September magnets (C1) would hold off its ions but pass the fast atoms that charge exchange makes and the unionized gas, so a 1 kg/s budget would need gridded ion thrusters and capture of 81–98% of that gas ([exhaust isolation](../solar_shield_array/integrated_comparison.md#where-the-held-screens-exhaust-goes)) |

## Open decisions for the author

The loss budget (R1), designed across 1–100 kg/s (cycle times of 100 billion to 1 billion years) for now; the ultraviolet
cut-off (O3); the dimmer's form (an even cut or one taken from the near
infrared, fixed or adjustable) and its range;
the use of the ring outside the window (transparent film or power collector);
the protected radius (O2); whether and when to build charged-particle
protection, which a 1 kg/s budget needs (C1); and the share of the budget left
to outflows under S7.
