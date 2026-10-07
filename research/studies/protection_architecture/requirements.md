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
magnetosphere holds the solar wind off. The ultraviolet-driven escape is the
traced count the author adopted on 7 October: the transmission passes gaps over
the whole aperture of the protected radius, and its light, the films' X-rays
(with FISM2's measured solar cycle) and the unfiltered light beyond the
aperture follow their slant paths through air swollen by the heat it takes
([the X-ray check](../solar_shield_array/integrated_comparison.md#the-x-ray-check)).
For the titania stack, with the ranges spanning the three treatments of the
upper air and quiet Sun to solar maximum, and counting the cases of those six
that cannot meet the budget at any transmission:

| Budget (cycle time) | Magnetosphere | Protected radius 3 lunar radii | 4 lunar radii | 6 lunar radii | 10 lunar radii |
|---|---|---|---|---|---|
| 1 kg/s (100 billion years) | none | 0.007–0.052% (4 cannot) | 0.011–0.057% (4 cannot) | 0.0005–0.059% (3 cannot) | 0.0006–0.060% (3 cannot) |
| | none, the ring fleet's wake | 0.012–0.065% (3 cannot) | 0.006–0.11% (2 cannot) | 0.013–0.13% (1 cannot) | 0.023–0.13% (1 cannot) |
| | September magnets | 0.003–0.11% (2 cannot) | 0.005–0.13% (1 cannot) | 0.018–0.14% (1 cannot) | 0.023–0.15% (1 cannot) |
| 10 kg/s (10 billion years) | none | 0.0015–0.14% (1 cannot) | 0.022–0.14% (1 cannot) | 0.0012–0.15% | 0.006–0.16% |
| | none, the ring fleet's wake | 0.002–0.14% (1 cannot) | 0.023–0.15% (1 cannot) | 0.003–0.17% | 0.009–0.18% |
| | September magnets | 0.017–0.15% (1 cannot) | 0.027–0.17% (1 cannot) | 0.005–0.18% | 0.009–0.18% |
| 100 kg/s (1 billion years) | none | 0.024–0.15% (1 cannot) | 0.012–0.22% | 0.022–0.25% | 0.032–0.26% |
| | none, the ring fleet's wake | 0.024–0.15% (1 cannot) | 0.012–0.22% | 0.022–0.25% | 0.034–0.27% |
| | September magnets | 0.024–0.15% (1 cannot) | 0.016–0.23% | 0.025–0.26% | 0.034–0.27% |

Solar maximum with the warmest upper air is the binding case: it cannot meet
1 kg/s at any protected radius up to 10 lunar radii, even with the September
magnets, and it meets 10 kg/s from 5 lunar radii with them and 6 without.
Without a magnetosphere the solar wind's charge exchange keeps three or four of
the six cases above 1 kg/s at any radius. The ring fleet's screen, about
20,000 km out, may leave the Moon in its unrefilled wake; if the wind is absent
within that wake's core (2.6 lunar radii behind a 4-radius screen, a
sensitivity resting on an assumed refill length), one to three cases are left.
The ozone-forming edge makes every budget harder: without a magnetosphere it
cannot meet 1 kg/s at any radius, and with the September magnets two of six
cases do, from 8 lunar radii; 10 kg/s holds for most cases only from 10 lunar
radii, and 100 kg/s from 6. Without the optical shield the loss is
300–80,000 kg/s.

The protected radius has to reach beyond the exobase. With its edge inside the
exobase, unfiltered sunlight beyond the aperture heats the thermosphere along
slant paths, and the air swells into that light until it runs away; the heating
alone needs the shadow to reach just beyond the exobase, to 2.75–4 lunar radii
at the transmissions ultraviolet-driven escape alone allows (0.010–0.25% for
the titania stack with a 4-radius shadow). At those transmissions, where that
escape already fills the budget, the sunlit exosphere outside the ring fleet's
4-radius shadow adds 5–273 kg/s with no magnetosphere and 1.7–168 kg/s with the
September magnets; most of that comes from its first few scale heights above
the exobase. With a 4-radius shadow the air runs away above 147 kg/s with
collisional upper air, 103 kg/s in LTE and 36 kg/s with all near-infrared
heating (19–29 kg/s behind the 200-nm edge), its exobase then at the shadow's
edge, so a looser budget asks a wider shadow. Earth's tide is included
throughout and multiplies molecular escape by 1.8–2.6.

**The design point.** A formation of replaceable cells can hold a UV
transmission of 3–4×10⁻⁵ with tight design choices, about 2×10⁻⁴ with standard
ones and 2.5×10⁻³ with relaxed ones
([protection/transmission](../../../protection/transmission/README.md)); covering
failed cells quickly is the main lever (O8). Traced along its slant paths, light
through gaps heats the upper air 2.5–5 times what a count over the disk gives
it, and the sky's Lyman-alpha glow heats it as much as a transmission of about
5×10⁻⁴: at the tight level the glow supplies 75–93% of the heat, at the
standard level 49–70%. These are the figures at the ring fleet's protected
radius of 4 lunar radii, rederived on 7 October with the traced count and
FISM2's X-ray cycle ([results/design_point.json](results/design_point.json),
from [run.py](run.py)).

At the relaxed level the air runs away in ten of the twelve titania cases, its
exobase reaching the shadow's edge. At the standard level behind the titania
stack the ultraviolet-driven loss, with Earth's tide, is 4×10⁻⁶–32 kg/s, with
the exobase at 1.95–3.8 lunar radii (187–301 K), and the heating needs a shadow
reaching 2.0–4.0 lunar radii; with all near-infrared heating at solar maximum
and holes four times the particle size the air runs away at 4 and needs 5. With the September
magnets the central total at 4 lunar radii is 0.0001–0.08 kg/s with
collisional upper air and in LTE at quiet Sun (cycle times over a trillion
years), 2.9–4.2 kg/s in LTE at solar maximum (34 to 23 billion years), 2.8–3.3
kg/s with all near-infrared heating at quiet Sun (35 to 30 billion years) and
129 kg/s with it at solar maximum (0.8 billion years). Without a magnetosphere
the same radius leaves 0.77–1.3 kg/s in the coolest three cases (130 to 76
billion years), 6.8–9.0 kg/s in the next two (14 to 11 billion years) and
193 kg/s in the warmest (0.5 billion years). If the ring fleet's wake holds the wind off its
2.6-radius core, the coolest three lose 0.07–0.31 kg/s without a magnetosphere
(1.4 trillion to 320 billion years); where the exobase lies outside that core
the wake changes little. At the tight level the central totals with the
September magnets are 3×10⁻⁵–0.2 kg/s in the four cooler cases (490 billion
years or longer), 0.84–0.92 kg/s
with all near-infrared heating at quiet Sun (117 to 107 billion years) and
19–23 kg/s with it at solar maximum (5.2 to 4.3 billion years). Behind the
200-nm edge the standard level runs away in eight of twelve cases at 4 lunar
radii, and the rest lose 12–15 kg/s with the September magnets (8 to 6.5 billion
years).

R1 is a long-term average, and the design point now also gives each state's mean
over solar cycles 23 and 24 (1997–2019), each calendar year's sunlight measured
band by band in FISM2's daily record, as the author chose on 7 October
([atmosphere/loss_response/cycle.py](../../../atmosphere/loss_response/cycle.py)).
The ultraviolet from 10 to 175 nm rises only 1.25–1.41 times at the last three
maxima, weighted by energy, where the escape model takes 2.5. With the column's
middle heating shape, the standard level at 4 lunar radii averages 0.84–1.6 kg/s
with no magnetosphere in the collisional and LTE cases (117 to 61 billion years)
and 0.0016–0.15 kg/s with the September magnets. With all near-infrared heating it
averages 24–27 kg/s with no magnetosphere and 13–16 kg/s with the magnets (4.1 to
3.6 and 7.5 to 6.3 billion years); its largest year, 2000, reaches 84 and 53 kg/s.

## Optical shield

| ID | Requirement | Source | Status |
|---|---|---|---|
| O1 | The UV transmission, the share of sunlight below 175 nm that reaches the protected region, averaged over time and over that region, stays within the allowed fraction for R1. Gaps, pinholes, edges, off-normal light and outages all count | loss response with its exosphere step; [protection/transmission](../../../protection/transmission/README.md) | Derived on the traced count, which the author adopted on 7 October with FISM2's measured X-ray cycle ([decisions.md](../../decisions.md); [the X-ray check](../solar_shield_array/integrated_comparison.md#the-x-ray-check)): light through gaps over the whole aperture follows its slant paths and heats the air 2.5–5 times a count over the disk. A swarm can hold 3×10⁻⁵–2×10⁻⁴ (O8). With the ring fleet's 4-radius shadow and the September magnets the standard level meets 1 kg/s with collisional upper air and in LTE at quiet Sun, gives 2.9–4.2 kg/s in LTE at solar maximum and 2.8–3.3 kg/s with all near-infrared heating at quiet Sun, and 129 kg/s or a runaway with that heating at solar maximum, where even 6 lunar radii leave 64–77 kg/s. The tight level meets 1 kg/s in every case but that one (19–23 kg/s). Below about 5×10⁻⁴ the sky's glow heats the upper air more than the transmission does |
| O2 | The protected radius covers the thermosphere's heating, which takes it just beyond the exobase (with its edge inside the exobase the unfiltered light beyond the aperture runs the air away), and reaches far enough beyond the exobase that the sunlit exosphere's losses stay within their share of R1. At the design point's standard level, with the September magnets, the exosphere's central loss falls to 1 kg/s at 2.1–2.7 lunar radii in the three coolest cases, 4.9–5.3 in LTE at solar maximum and with all near-infrared heating at quiet Sun, and 11 with that heating at solar maximum; to 0.1 kg/s at 2.3–3.5 and 9.0–9.4 lunar radii, and beyond 11.5 in the warmest. At 20,000 km apertures of 6 and 10 lunar radii have 2.2 and 6.2 times the ring fleet's 4-radius area | [atmosphere/loss_response](../../../atmosphere/loss_response/README.md), absorption and exosphere steps; [protection/report.md](../../../protection/report.md), section 2 | Derived |
| O3 | The ultraviolet cut-off: the titania stack (no ozone, surface UV index 0.1) or an ozone-forming edge near 200 nm | [atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/README.md) | Open: the author's choice |
| O4 | The window passes the climate's sunlight (E2), evenly unless a pattern is chosen and run through the climate model. At 78,000 km an even cut needs the dimmer to reach 2,100 km from the axis, or 2,400 km to include the twilight band | E2; this study | Derived |
| O5 | Light the shield rejects never reaches Earth's night side: specular reflections are steered at least about 5° off the plane in which Earth moves as seen from the shield, and Sun-facing surfaces keep diffuse reflectance near 0.1% | this study | Derived; to confirm in illumination |
| O6 | The optical shield works before the atmosphere accumulates | feasibility report, section 12; loss response | Derived |
| O7 | Effects on Earth are assessed and bounded: reflections, and the shield crossing the full Moon and the Sun as seen from Earth | protection/report.md, section 6 | Required assessment |
| O8 | Neighbouring cells overlap by several times their relative position error (50 m in the design levels, 1% more film); cells are patched or replaced at least every 10–20 years against micrometeoroid holes; a failed cell is covered within a day to a week by spares or overlays kept in the formation | [protection/transmission](../../../protection/transmission/README.md) | Derived; cell failure rates, the time to cover and the film's hole size are assumptions to test |

## Charged-particle protection

| ID | Requirement | Source | Status |
|---|---|---|---|
| C1 | Solar-wind loss stays within its share of R1 | loss response, with its absorption and exosphere steps | Derived from screening, to confirm with a plasma model. Without a magnetosphere the solar wind carries off the ions the sunlit exosphere makes, and its charge exchange with the dense exosphere near the exobase costs 0.7–6 kg/s (central) at the design point's standard level, so 1 kg/s needs charged-particle protection. The ring fleet's screen may hold the wind off the inner exosphere in its unrefilled wake: as a sensitivity, the three coolest cases then lose 0.07–0.31 kg/s with no magnetosphere, a plasma question still to test. A lunar dipole whose stand-off clears the dense exosphere removes the charge exchange and keeps a quarter to two fifths of the ions it encloses at these exobases: at the standard level the coolest case loses under 1 kg/s without one, the next two need about 2×10¹⁹ A·m² (a stand-off of 2.3 lunar radii, just above their exobase), and no moment up to 10²³ A·m² holds the warmer three to 1 kg/s, where the September magnets (1.5×10²¹ A·m², stand-off 10 lunar radii) leave 2.8–129 kg/s in all. Stand-offs just above the exobase fall where simulations find a weak field can raise escape. Magnets do not hold the neutral fragments, which set the protected radius (O2) |
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
