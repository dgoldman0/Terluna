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
| Dawn and dusk follow the light on the ground under the solved clear sky. Practical dusk, while ordinary outdoor activity needs no artificial light, runs from sunset until the light falls to the level at which Earth's civil twilight ends: 2.98 lux in the same solver, with the Sun 6° down on Earth. The Moon reaches it with the Sun 41.4° below the horizon, 82 hours after sunset at the equator and 98 hours at 30° latitude; poleward of 48.6° a clear night never falls that low. Day vision for an 18% grey surface (CIE 191:2010) holds for the first 46 hours at the equator, to a depression of 23.6°. Dawn mirrors dusk. Where the Earth stands high, its light keeps the night above the level: 5–6 lux at the darkest at three nearside coasts with the Earth 43–74° up. | July; revised 2026-10-04 at the author's request, from the sea-appearance study's [lighting calendar](studies/sea_appearance/README.md#the-lighting-calendar) | Stands. The July canon gave about five hours (2.54° of solar motion at the equator), replacing a 30-hour figure, and was held on 2026-09-26 against the ensemble plan's nominal 12.5-hour sunset pending a photometric definition for the actual atmosphere. On 2026-10-04 the author asked that this row follow the current understanding. The threshold carries civil twilight's meaning to the Moon through the solver's Earth control ([optical comfort](studies/optical_comfort/README.md)). Clear sky, from a horizontally uniform column without refraction; clouds come in the sea-appearance study's results by regime. |
| How the seas look is studied by day, through the long evening, under earthlight and in the far side's twilit night, in this order: a lighting calendar of six coasts with a map of where the Earth stands over the seas, the surface's reflection and the water's colour, results by regime, then renderings of four coasts (western Oceanus Procellarum by Russell, the eastern Smythii headland, southern Mare Nubium and the South Pole–Aitken coast by Mare Ingenii). The far-side rendering shows the twilit night first, as the more representative case. Cloud ice keeps the cloud scheme's own sizes, with 140 µm as a check. | Author, 2026-10-04; [research/studies/sea_appearance](studies/sea_appearance/README.md) | Stands. The calendar found the twilight wrapping round the Moon: at latitude φ the Sun sinks at most 90° − φ below the horizon, so only seas near the far side's equator grow dark around local midnight, and the author chose the twilit night for the far-side coast. The calendar also gave dawn and dusk their photometric definition (the row above). |

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
| The wave work turns to the main seas and their water levels: next come the nearside sea's waves through the 60-day global atmospheric record, a tide product in geography, and one synthesis of the wave studies. The nearside sea runs at a 300 s timestep, checked against 150 s over its stormiest days. | Author, 2026-10-03; [climate/waves/review.md](../climate/waves/review.md) | Done on 2026-10-03. Through the second lunar cycle the nearside sea's mean Hs is 1.39 m (Smythii–Marginis 0.88 m), Hs of at least 1 m covers 68% of its area and time, and its largest Hs is 4.9 m; the seas travel toward the equator, and the most exposed coasts take up to 731 W/m on average ([climate/waves/nearside.md](../climate/waves/nearside.md)). The tide product is in [geography](../geography/README.md#the-monthly-tide) and the synthesis in [climate/waves/README.md](../climate/waves/README.md). |
| The review of the wave studies stands as their assessment. SWAN's and SWASH's lunar handling holds, and the runs forced by the GCM carry lunar air–sea momentum transfer: ExoPlaSim's open-sea roughness follows Charnock's law at lunar gravity. The studies leave out the seas beyond Smythii–Marginis (which holds 2.7% of the water), the monthly tide (a typical monthly range of 3.7 m over the nearside sea, 2.1 m in the South Pole–Aitken sea and 0.55 m in Smythii–Marginis, from the geography tide product), coasts of regolith whose sediment moves far more readily at lunar gravity, and how the waves look (lunar seas move at 0.41 times the speed of Earth seas of the same wavelength). Their largest uncertainties are the 1° basin grid, the forcing and the single lunar cycle. | Author, 2026-10-03 (asked for the analysis to be locked in); [climate/waves/review.md](../climate/waves/review.md), [results/review_checks.json](../climate/waves/results/review_checks.json) | Stands. The geography tide product (2026-10-03) replaced the review's first tide estimate; it also finds Mare Fecunditatis exchanging water with the nearside sea through a narrow strait every 9.5 days, whose geometry and friction set its fortnightly tide: twice the equilibrium with bed friction for slow flows, 0.8–1.05 times it with the friction of the metre-per-second flows the strait would carry. The review's other recommendations stay open: a Janssen drag-and-growth run, a 0.5° basin parent, further cycles, a wider frequency band, the seas' composition, sediment and beach evolution, and figures drawn over relief. |

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
| Data products record the shared constants their code reads by name and value instead of pinning `shared/constants.json` or `constants.py` by hash, and their tests compare those values with the current ones. Adding a constant then leaves every product current, while a change to a value a product used shows up. Workflow tools such as DVC and Snakemake track a step's parameters the same way. | Author, 2026-10-04; [shared/provenance.py](../shared/provenance.py), [shared/README.md](../shared/README.md#data-products) | Stands. The optical, sky and cloud-twilight producers follow it; regenerated on 2026-10-04, every data value reproduced and the sky's embedded numerical checks were refreshed with the current check code. The wave, tide and coastal products still record the constants file's hash as a record of their run, without a test of it. |

## Life and people

| Decision | Source | Status |
|---|---|---|
| Aerial and gliding life is a major feature of the world. | July | Stands |
| High-altitude life is sparse: microbes or engineered radiation-hard films. | July | Stands |
| Long-endurance high-altitude platforms serve farside astronomy and nearside earthshine photometry. | July | Stands |
| The seas are living: optical studies use stated best guesses with ranges for plankton pigments, dissolved organic matter and bioluminescence, drawn from the biosphere and lunar-cycle ecology work and labelled as design guesses, beside the regolith fines. | Author, 2026-10-04; [research/studies/sea_appearance](studies/sea_appearance/README.md) | Stands |
