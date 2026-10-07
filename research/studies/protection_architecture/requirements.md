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
| R2 | The budget is shared among ultraviolet-driven escape (with Earth's tide), the sunlit exosphere's losses (the ions it makes, the fragments of broken molecules and the solar wind's charge exchange), sputtering, and the oxygen atoms and hydrogen the upper air's photochemistry makes, which leave by Jeans escape, as ions from the sunlit exosphere and hot from near the exobase | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), with its exosphere and oxygen steps | Allocation open. Behind the titania stack at the standard level the oxygen atoms add 0.002–0.006 kg/s with collisional upper air and in LTE at the solar cycle's mean spectrum, with the Moon-scaled mixing or the breaking gravity waves' (gravity_waves.py), and 0.10–0.12 kg/s with all near-infrared heating over the cycle's years (0.05–0.06 with the September magnets); hot atoms are 0.001–0.002 kg/s of this. Weaker mixing raises the warmest case's atoms to 0.45 kg/s (the waves at the low end of Earth's measured mixing) and 0.77 kg/s (Titan's measured mixing carried to the Moon, 0.37 with the magnets); Earth's own unscaled mixing, which neither supports, would give 2.1 kg/s (1.0). Hydrogen adds 0.001–0.003 kg/s; behind the 200-nm edge the atoms add 2–4 kg/s |

What each budget asks of the protection, from the loss response with its
exosphere step. The allowed UV transmission is the largest share of sunlight
below 175 nm that may reach the protected region once the sunlit exosphere's
losses (the ions it makes, the fragments of broken molecules and the solar
wind's charge exchange, central estimate) are added to ultraviolet-driven
escape. It depends on how far the protected region reaches and on whether a
magnetosphere holds the solar wind off. The ultraviolet-driven escape is the
traced count the author adopted on 7 October: the transmission passes gaps over
the whole aperture of the protected radius, and its light, the films' X-rays
(with FISM2's measured solar cycle) and the unfiltered light beyond the
aperture follow their slant paths through air swollen by the heat it takes
([the X-ray check](../solar_shield_array/integrated_comparison.md#the-x-ray-check)).
The thermal column takes that heat where the tracing puts it in height, with the
sky's glow just above the base, as the author approved the same day. For the
titania stack, with the ranges spanning the three treatments of the upper air and
quiet Sun to solar maximum (FISM2's year around cycle 19's maximum; the escape
model's 2.5, the stress case, is given apart), and counting the cases of those six
that cannot meet the budget at any transmission:

| Budget (cycle time) | Magnetosphere | Protected radius 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.031–0.099% (3 cannot) | 0.0053–0.10% (2 cannot) | 0.006–0.10% (2 cannot) | 0.006–0.10% (2 cannot) |
| | none, the ring fleet's wake | 0.009–0.11% (2 cannot) | 0.017–0.14% (1 cannot) | 0.008–0.15% | 0.015–0.16% |
| | September magnets | 0.012–0.15% (1 cannot) | 0.0026–0.17% | 0.011–0.18% | 0.016–0.18% |
| 10 kg/s (10 billion years) | none | 0.024–0.15% (1 cannot) | 0.012–0.19% | 0.022–0.21% | 0.028–0.21% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.013–0.19% | 0.024–0.22% | 0.032–0.23% |
| | September magnets | 0.024–0.15% (1 cannot) | 0.018–0.21% | 0.027–0.22% | 0.033–0.23% |
| 100 kg/s (1 billion years) | none | 0.024–0.15% (1 cannot) | 0.018–0.25% | 0.049–0.28% | 0.065–0.30% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.018–0.25% | 0.050–0.29% | 0.068–0.30% |
| | September magnets | 0.024–0.15% (1 cannot) | 0.018–0.25% | 0.055–0.29% | 0.069–0.30% |

Without a magnetosphere the solar wind's charge exchange keeps the two cases with
all near-infrared heating above 1 kg/s at any radius up to 10 lunar radii. The
ring fleet's screen, about 20,000 km out, may leave the Moon in its unrefilled
wake; if the wind is absent within that wake's core (2.6 lunar radii behind a
4-radius screen, a sensitivity resting on an assumed refill length), only the
warmest at solar maximum is left at 4 lunar radii, and none from 6. With the
September magnets every case meets 1 kg/s at 4 lunar radii, the warmest at solar
maximum only below a transmission of 2.6×10⁻⁵ (5.4×10⁻⁵ in the stress case).
The ozone-forming edge makes every budget harder: without a magnetosphere it
cannot meet 1 kg/s at any radius, and with the September magnets two of six
cases do at 4 lunar radii, three at 6 and five at 10; five of six meet 10 kg/s
from 4 lunar radii. Without the optical
shield the loss is 300–80,000 kg/s.

The protected radius has to reach beyond the exobase. With its edge inside the
exobase, unfiltered sunlight beyond the aperture heats the thermosphere along
slant paths, and the air swells into that light until it runs away; the heating
alone needs the shadow to reach just beyond the exobase, to 3.5–4.5 lunar radii
at the transmissions ultraviolet-driven escape alone allows (0.018–0.26% for
the titania stack with a 4-radius shadow, 0.016% in the stress case). At those
transmissions, where that escape already fills the budget, the sunlit exosphere
outside the ring fleet's 4-radius shadow adds 5–151 kg/s with no magnetosphere
and 1.7–89 kg/s with the September magnets (up to 174 and 108 in the stress
case); most of that comes from its first few scale heights above the exobase.
With a 4-radius shadow the air runs away from a transmission of 0.11–0.26% with
collisional upper air (solar maximum to quiet Sun; 0.079% in the stress case),
0.076–0.19% in LTE (0.060%) and 0.018–0.081% with all near-infrared heating
(0.021%), and from 0.0012–0.049% behind the 200-nm edge (0.0087%), as its
exobase nears the shadow's edge. The tables find those
onsets in air that swells more than air heated where the tracing puts it, so
they lean early. Earth's tide is included throughout and multiplies molecular
escape by 1.8–2.7.

**The design point.** A formation of replaceable cells can hold a UV
transmission of 3–4×10⁻⁵ with tight design choices, about 2×10⁻⁴ with standard
ones and 2.5×10⁻³ with relaxed ones
([protection/transmission](../../../protection/transmission/README.md)); covering
failed cells quickly is the main lever (O8). Traced along its slant paths, light
through gaps heats the upper air 2.5–5 times what a count over the disk gives
it, and the sky's Lyman-alpha glow heats it as much as a transmission of about
5×10⁻⁴: at the tight level the glow supplies 75–93% of the heat, at the
standard level 49–70%. The glow's heat lands within a few e-folds of pressure
above the base and is conducted away, and the far ultraviolet through gaps lands
in the lower thermosphere, so placed where the tracing puts it the heat leaves
the exobase 26–51 K cooler than the column's 'middle' shape gives. These are the
figures at the ring fleet's protected radius of 4 lunar radii, rederived on
7 October with the traced count, FISM2's X-ray cycle and the traced heating shape
([results/design_point.json](results/design_point.json), from [run.py](run.py)).

R1 is a long-term average, so the design point also gives each state's mean
over solar cycles 23 and 24 (1997–2019), each calendar year's sunlight measured
band by band in FISM2's daily record, as the author chose on 7 October
([atmosphere/loss_response/cycle.py](../../../atmosphere/loss_response/cycle.py)).
The ultraviolet from 10 to 175 nm rises only 1.25–1.41 times at the last three
maxima, weighted by energy, where the escape model takes 2.5 for its solar
maximum. Each year is taken as a steady state, so the mean leans high against an
upper air that follows the cycle with a lag.

Averaged over the cycle, at the standard level behind the titania stack and with
4 lunar radii, the central total is 0.68–0.97 kg/s with no magnetosphere with
collisional upper air and in LTE (144 to 101 billion years), almost all of it the
solar wind's charge exchange and sputtering; 0.064–0.093 kg/s in the ring fleet's
wake; and 10⁻⁵–0.003 kg/s with the September magnets. With all near-infrared
heating it is 3.1–3.5 kg/s with no magnetosphere (32 to 28 billion years),
2.3–2.7 kg/s in the wake and 0.75–0.95 kg/s with the magnets (131 to 104 billion
years); its largest year, 2000, reaches 6.6–7.9, 5.8–7.2 and 2.5–3.2 kg/s. A
6-radius shadow brings that case's mean to 2.1–2.3 kg/s with no magnetosphere and
0.38–0.47 kg/s in the wake. At the tight level the means are 0.63–0.86 kg/s with
no magnetosphere in the cooler cases and 1.7–1.8 kg/s with all near-infrared
heating, which the wake brings to 0.79–0.85 and the magnets to 0.13–0.15 kg/s.
At solar maximum, FISM2's year around cycle 19's maximum of 1958, the warmest
case's air runs away at the standard level: its heating needs a 5-radius shadow,
and at 6 lunar radii it loses 7.7 kg/s (3.6 with the magnets). The cooler cases
lose 0.78–1.3 kg/s there with no magnetosphere and 0.0002–0.04 kg/s with the
magnets, about as in the stress case, where the warmest loses 26 kg/s (15 with
the magnets) and runs away with holes four times the particle size. At the tight
level the warmest at solar maximum loses 3.8 kg/s (1.0 with the magnets; 2.7 and
0.59 in the stress case). The measured years of cycles 23 to 25 reach at most
7.9 kg/s. The oxygen
atoms the upper air's photochemistry makes add little to the cooler cases
(R2), and to the warmest 0.10–0.12 kg/s over the cycle with the Moon-scaled
mixing or the breaking gravity waves' (0.05–0.06 with the magnets, which brings
that case to 0.80–1.0 kg/s); with Titan's measured mixing carried to the Moon
they add 0.77 kg/s (0.37 with the magnets, bringing it to 1.1–1.3 kg/s). At the
relaxed level the air runs away in most cases and years. Behind the 200-nm edge
the standard level runs away in the years around maximum (in every year with all
near-infrared heating), and the tight level
averages 3.1–3.4 kg/s with no magnetosphere in the cooler cases (0.69–0.84 kg/s
with the magnets), its warmest case running away in eleven of the 23 years.

## Optical shield

| ID | Requirement | Source | Status |
|---|---|---|---|
| O1 | The UV transmission, the share of sunlight below 175 nm that reaches the protected region, averaged over time and over that region, stays within the allowed fraction for R1. Gaps, pinholes, edges, off-normal light and outages all count | loss response with its exosphere step; [protection/transmission](../../../protection/transmission/README.md) | Derived on the traced count, which the author adopted on 7 October with FISM2's measured X-ray cycle ([decisions.md](../../decisions.md); [the X-ray check](../solar_shield_array/integrated_comparison.md#the-x-ray-check)), the heat placed where the tracing puts it: light through gaps over the whole aperture follows its slant paths and heats the air 2.5–5 times a count over the disk. A swarm can hold 3×10⁻⁵–2×10⁻⁴ (O8). With the ring fleet's 4-radius shadow and averaged over the solar cycle, the standard level keeps every titania case within 10 kg/s with no magnetosphere (at most 3.5 kg/s, 7.9 in the worst year) and within 1 kg/s with the September magnets (at most 0.95 kg/s). At solar maximum, FISM2's year around cycle 19's maximum, the warmest case runs away at the standard level and the cooler ones lose 0.78–1.3 kg/s with no magnetosphere; in the stress case, the escape model's 2.5, the warmest loses 26 kg/s and runs away with larger holes. The relaxed level runs away in most cases. Below about 5×10⁻⁴ the sky's glow heats the upper air more than the transmission does, so low that it raises the exobase only 6–18 K |
| O2 | The protected radius covers the thermosphere's heating, which takes it just beyond the exobase (with its edge inside the exobase the unfiltered light beyond the aperture runs the air away), and reaches far enough beyond the exobase that the sunlit exosphere's losses stay within their share of R1. At the design point's standard level, with the September magnets, the exosphere's central loss falls to 1 kg/s at 1.9–2.2 lunar radii in the three coolest cases and 2.7–3.2 in LTE at solar maximum and with all near-infrared heating at quiet Sun, and to 0.1 kg/s at 2.1–2.5 and 3.5–4.5 lunar radii. With that heating at solar maximum the heating alone needs 5 lunar radii, and at 6 the case loses 7.7 kg/s (3.6 with the magnets); in the stress case its loss falls to 1 kg/s at 8.2 and to 0.1 kg/s at 10.9 lunar radii. At 20,000 km apertures of 6 and 10 lunar radii have 2.2 and 6.2 times the ring fleet's 4-radius area, though the ring fleet's tile mass grows with the aperture's width, about a fifth for each added lunar radius | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), absorption and exosphere steps; [protection/report.md](../../../protection/report.md), section 2 | Derived |
| O3 | The ultraviolet cut-off: the titania stack (no ozone, surface UV index 0.1) or an ozone-forming edge near 200 nm | [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Open: the author's choice |
| O4 | The window passes the climate's sunlight (E2), evenly unless a pattern is chosen and run through the climate model. At 78,000 km an even cut needs the dimmer to reach 2,100 km from the axis, or 2,400 km to include the twilight band | E2; this study | Derived |
| O5 | Light the shield rejects never reaches Earth's night side: specular reflections are steered at least about 5° off the plane in which Earth moves as seen from the shield, and Sun-facing surfaces keep diffuse reflectance near 0.1% | this study | Derived; to confirm in illumination |
| O6 | The optical shield works before the atmosphere accumulates | feasibility report, section 12; loss response | Derived |
| O7 | Effects on Earth are assessed and bounded: reflections, and the shield crossing the full Moon and the Sun as seen from Earth | protection/report.md, section 6 | Required assessment |
| O8 | Neighbouring cells overlap by several times their relative position error (50 m in the design levels, 1% more film); cells are patched or replaced at least every 10–20 years against micrometeoroid holes; a failed cell is covered within a day to a week by spares or overlays kept in the formation | [protection/transmission](../../../protection/transmission/README.md) | Derived; cell failure rates, the time to cover and the film's hole size are assumptions to test |

## Charged-particle protection

| ID | Requirement | Source | Status |
|---|---|---|---|
| C1 | Solar-wind loss stays within its share of R1 | loss response, with its absorption and exosphere steps | Derived from screening, to confirm with a plasma model. Without a magnetosphere the solar wind carries off the ions the sunlit exosphere makes, and its charge exchange with the dense exosphere near the exobase leaves the exosphere losing 0.65–2.0 kg/s (central) in the five cases that keep a state at the design point's standard level (the sixth runs away), so 1 kg/s needs charged-particle protection or the ring fleet's wake. The ring fleet's screen may hold the wind off the inner exosphere in its unrefilled wake: as a sensitivity, the three coolest cases then lose 0.064–0.074 kg/s with no magnetosphere, a plasma question still to test. A lunar dipole whose stand-off clears the dense exosphere removes the charge exchange and keeps 28–42% of the ions it encloses at these exobases: at the standard level the three coolest cases lose under 1 kg/s from their exosphere without one, LTE at solar maximum and all near-infrared heating at quiet Sun need about 3×10¹⁹ A·m² (a stand-off of 2.7 lunar radii), and the warmest at solar maximum finds no state within the 4-radius shadow, which no magnet changes (in the stress case no moment up to 10²³ A·m² holds it to 1 kg/s). Averaged over the solar cycle, the September magnets (1.5×10²¹ A·m², stand-off 10 lunar radii) keep every titania case within 1 kg/s at the standard level. Stand-offs just above the exobase fall where simulations find a weak field can raise escape. Magnets do not hold the neutral fragments, which set the protected radius (O2) |
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
protection, which a 1 kg/s budget needs unless the ring fleet's wake holds the
solar wind off (C1); and the share of the budget left to outflows under S7.
