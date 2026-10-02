# Decisions

The author's decisions about the Open Moon (Terluna), with where each is recorded
and whether it stands. A later decision replaces an earlier one, and the
replaced entry names what replaced it. Modelling choices, reference values and
results appear here only where the author adopted them.

Sources:

- **July**: the recovered canon of the July 2026 knowledge bundle
  ([archive/july_2026_knowledge_bundle/](../archive/july_2026_knowledge_bundle/README.md);
  its file `01_RECOVERED_CANON_AND_CORRECTIONS.md` and `data/decision_log.csv`).
- **Sep 4**: the user requirements of the September feasibility report
  ([baselines/feasibility/report.md](baselines/feasibility/report.md), section 1).
- **Sep 9**: the September protection report
  ([protection/report.md](../protection/report.md), section 3).
- **Sep 19**: the 19 September 2026 industry discussion, known through a
  reconstruction written on 2026-09-26
  ([engineering/reference/industrial_architecture/](../engineering/reference/industrial_architecture/README.md),
  `RECOVERED_INDUSTRIAL_ARCHITECTURE.md`, section 8). The author confirmed its
  six traffic-safety rules on 2026-09-26.
- **Author**: decisions the author stated in working sessions, with the date.
- The conservation study keeps its own decisions, D1–D11, in
  [its README](studies/conservation/README.md#4-decisions); they are not repeated here.

## The world

| Decision | Source | Status |
|---|---|---|
| The world is Earth's own Moon, transformed in place, and the project is called Terluna. | July | Stands. The paper ensemble's working title is *Constructing and Sustaining an Open Moon*. |
| The Moon keeps its orbit, gravity, synchronous rotation and 29.53-day solar cycle. Life and culture adapt to the long day and night, storing resources across them. | July | Stands |
| Construction within 500 years of 2026, and operation for at least 10⁹ years. | Sep 4 | Stands |
| An open atmosphere with no pressure dome over the Moon. | Sep 4 | Stands. Local enclosed places, such as the undersea Tranquility community (conservation D6), are separate. |

## Light and time

| Decision | Source | Status |
|---|---|---|
| Practical dawn and dusk each last about five hours, inside the 29.53-day cycle. This replaced an earlier 30-hour figure. At the equator five hours is about 2.54° of solar elevation. | July | Stands; settled on 2026-09-26 against the ensemble plan's nominal 12.5-hour sunset, a design case. It still needs a photometric definition for the actual atmosphere. |

## Atmosphere

| Decision | Source | Status |
|---|---|---|
| Surface pressure 1.2 atm (121,590 Pa). | Sep 4 | Stands; the design pressure of the climate runs. The July bundle's 80 kPa was its own reference value. |
| An exobase near 250 K is preferred, with 260 K a candidate limit. | Sep 4 | Stands. The solved upper air behind the titania film gives 142–193 K before the sky's Lyman-alpha glow is counted ([middle atmosphere](../atmosphere/middle_atmosphere/README.md)). |
| Water, nitrogen and oxygen come from the Solar-System-wide resource operation. | Conservation D4 | Stands; replaced the September baseline's oxygen from lunar rock. |
| Atmospheric layers take descriptive names, never L1–L5, which stay reserved for Lagrange points. | July | Stands |
| "Gravity wave" means an atmospheric buoyancy wave; write "atmospheric gravity wave" on first use. | July | Stands |
| A band at 35–45 km was placed as the top of the main troposphere, provisionally, pending a climate model. | July | Settled on 2026-09-26: the band stays as a named flight band, and its physics comes from the models. The 1-D column puts the tropopause near 86 km (288 K surface); in the 5% climate the band holds air at −10 to −20 °C and 0.58–0.48 atm ([climate/gcm](../climate/gcm/README.md)). |

## Water and climate

| Decision | Source | Status |
|---|---|---|
| Standing water covers 28% of the surface. | [shared/scenarios/water.json](../shared/scenarios/water.json); conservation D10 | Stands; selected on 2026-09-25 within the author's 25–35% range. |
| The shield passes 5% less sunlight than the titania stack at every wavelength, and the Moon settles near 295 K (294.9 K with the GCM's Sun, tilt and sunlight split corrected; 294 K before). | [climate/gcm](../climate/gcm/README.md) | Stands; chosen on 2026-09-26 and kept on 2026-09-29 after the GCM's correction. |
| The corrected climate's drier land and smaller lakes are accepted for the comfort it brings. Rain-fed lakes cover 10.3% of the Moon (11.6% before the GCM's correction), 42% of the land gets under 0.5 mm of rain a day, and the GCM's comfortable hours over land rise, most of all poleward of 45°. | Author, 2026-09-29; [climate/gcm](../climate/gcm/README.md), [geography](../geography/README.md) | Stands. CM1's check disagrees: its three flat rings (2026-09-30) find no low-lying land comfortable, and corrected for height 42 comfortable hours per lunar day on their land against the GCM's 387 at the same places, with comfort only on high ground poleward of 30°. Which model is right about the air near the ground over land is open ([climate/crm](../climate/crm/README.md)); comfort is carried as the range between them while the climate work pauses (Research practice). |
| The CM1 rings run flat and are read for low-lying land; the highlands' climate comes from their results corrected for height (about 1 °C cooler per kilometre with lower dewpoints, from the terrain ring's first eleven days, checked against the GCM's highland cells). Rings with terrain are not used: in two dimensions air forced over 8-km mountains turned into gales along the whole ring. The highlands' own weather (winds up and down slopes, rain on the slopes facing the wind, cold air pooling in valleys) waits for a possible 3-D box, 385 km square, over high ground on `ring_70_45e`'s path, forced by that ring and decided after its results. | Author, 2026-09-29; [climate/crm](../climate/crm/README.md) | Stands. The flat rings and the height correction were done on 2026-09-30. The 3-D box ran on 2026-10-01 at the author's go-ahead, at 44.7° S, 246.1° E, and again on 2026-10-02 with each column held and forced at its own height, the first run having taken its forcing from one corner column at 2.9 km: dry and sunny, with 183 comfortable hours per lunar day rising with height from 82 below 1 km to 262 at 5-6 km, muggy mornings and some cloud at sunrise on the low ground, slope winds and mild cold pooling at night ([climate/crm](../climate/crm/README.md), "A box over the highlands"). |
| Investigate atmospheric electricity in stages: compare saved CM1 graupel and mixed-phase diagnostics, use matched 3-D boxes for follow-up and rings for context, evaluate WRF-ELEC for charging/fields/lightning while retaining CM1 as the baseline, then add conductivity, global-circuit and TLE models. ROCKE-3D remains a separate background-climate check. | Author, 2026-10-02; [atmospheric-electricity plan](studies/atmospheric_electricity/README.md) | Stands as the agreed research direction. The first task is the output inventory and controlled comparison design. Graupel changes remain unmeasured; electrical software selection/adaptation and numerical runs remain future work. The wider GCM–CM1 discrepancy and existing compute limits remain open and binding, respectively. |

## Protection

| Decision | Source | Status |
|---|---|---|
| The protection architecture: a solar-filter complex with an industrial hub; magnetic protection; protection that grows around occupied destinations and Earth–Moon traffic corridors, secondary to protecting the lunar atmosphere and surface; Earth–Sun L1/L2 hubs for freight, power and industry. | Sep 9, carrying the original 2025 architecture | Stands |
| No superconducting planetary ring. | Sep 4 | Stands |
| The protection is designed across a range of total loss budgets, 1–100 kg/s, until the designs show which rate works best. Each budget is stated with its atmospheric cycle time, the time the loss and the resupply that matches it take to replace the whole atmosphere: about 100 billion years at 1 kg/s, 10 billion at 10 kg/s and 1 billion at 100 kg/s. | Author, 2026-09-26 | Stands |
| A surface radiation dose of at most 0.027 mSv/day. | Sep 4 | Stands |
| Protection hardware is a formation of replaceable units, never a single megastructure. Ageing and damage are handled by replacing units, supplied over the long term. | Author, 2026-09-26 | Stands |
| The protection design covers the Moon and its atmosphere. Destinations and corridors come later, and the older Earth–Moon–L1–L2 megastructure concept stays outside it. | Author, 2026-09-26 | Stands |

## Transport and safety

The author confirmed these rules on 2026-09-26 and regards them as set, to be reopened only for good reasons.

| Decision | Source | Status |
|---|---|---|
| Passive failures miss planets: nominal transfer paths and their failure dispersions avoid inhabited bodies unless capture has been verified. | Sep 19 | Stands |
| Earth never lies behind the lunar catcher, so a missed lunar capture cannot become an Earth-impact trajectory. | Sep 19 | Stands |
| Major braking happens remotely: high-energy packets approach the inhabited Moon only after capture and deceleration. | Sep 19 | Stands |
| Packet energy is capped. Mature dense traffic uses packets of roughly 10⁶–10⁸ kg in preference to 10¹²–10¹⁵ kg units, so each failure's consequences stay bounded; the exact cap awaits risk optimization. | Sep 19 | Stands |
| Planetary-scale surveillance, tracking, interception and tug capability, and assigned arrival corridors are part of the infrastructure. | Sep 19 | Stands |
| No civilization-threatening kinetic payload enters a planetary intercept corridor before verified capture. | Sep 19 | Stands |

## Research practice

| Decision | Source | Status |
|---|---|---|
| Simulation work fits the author's laptop (8 cores and 31 GB shared with the author's own work, runs of hours with rests between); there is no budget for rented compute. | Author, 2026-09-29 | Stands |
| Before any further CM1 run, CM1 is rebuilt for speed and tested: gfortran at `-O3` for this laptop's processor (`-march=native`), and a milder variant without fused multiply-adds or glibc's vector maths that may reproduce the current build bit for bit. A rebuild is used for every later run if it runs faster and matches the current build, byte for byte or within the current build's own spread after a one-bit nudge. The flat rings running on 2026-09-29 (`ring_70_45e`, `ring_equator`, `ring_70_135e`) finish on the current build. | Author, 2026-09-29; [climate/crm](../climate/crm/README.md) | Stands. Tested on 2026-09-30: neither rebuild matched the current build byte for byte, and both stayed within the spread of two nudges of it (the last bit of every value of its scalar restart file, since a single value's nudge could not reach the whole ring in a day). `-O3 -march=native`, 11–15% faster, is now every OpenMP build's setting; the MPI builds keep `-O2`, untested. |
| The climate work pauses. CM1 and the GCM disagree on how warm and humid the air over land is near the ground, and no affordable option settles which is right: the cheaper runs refine CM1's answer without an independent check, and a second GCM (ROCKE-3D) would take weeks of laptop sessions. Comfort is carried as the range between the two models (on the CM1 rings' land, 42 comfortable hours per lunar day in CM1 corrected for height against 387 in the GCM), and both put the most comfortable land poleward of 30°. The corrected GCM design case and the products built from it stand. The work reopens with an independent GCM or new evidence on daytime cloud over a slowly rotating world. | Author, 2026-09-30; [climate/crm](../climate/crm/README.md) | Stands. On 2026-10-01, at the author's go-ahead, two GCM runs tested its cloud: convective cloud thinned toward CM1's warmed the whole Moon by 1.7 K, and low cloud added toward CM1's foggy nights cooled it by 5.0 K, each about alike from the ground to 18 km, so the GCM's cloud does not explain the gap. What remains is how each model treats the lowest few kilometres. The runs also put the GCM's own climate between about 289 and 297 K across those cloud settings, and its comfortable hours on the rings' land between 121 and 387. |

## Life and people

| Decision | Source | Status |
|---|---|---|
| Aerial and gliding life is a major feature of the world. | July | Stands |
| High-altitude life is sparse: microbes or engineered radiation-hard films. | July | Stands |
| Long-endurance high-altitude platforms serve farside astronomy and nearside earthshine photometry. | July | Stands |
