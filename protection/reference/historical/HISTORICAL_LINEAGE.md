# Historical solar-shield / protection lineage

This file deliberately preserves superseded branches because they contain ideas that may still be useful. Values labeled **conversation-recovered** came from prior project conversations and are not original artifact bytes.

## 2025-08 — shared magnetosphere / plasma-channel ideas — HISTORICAL

Early work explored a superconducting-loop / artificial plasma-magnetic channel in the Earth–Moon region to reduce solar-wind interaction, sputtering and charged-particle exposure. A separate speculative plasma-funnel idea considered collecting ionized terrestrial/exospheric gases and guiding them toward lunar collection facilities. Those gas-supply rates were far below later 500-year atmosphere-build requirements and are not part of the current protection baseline.

Reusable idea: magnetic/plasma infrastructure can be staged around destinations and routes. Superseded claim: it cannot substitute for a physical EUV/X-ray photon filter.

## 2025-09-26 — Lunashield-L1 mini-magnetosphere — HISTORICAL / charged-particle subsystem

Conversation-recovered baseline for `Lunashield-L1`:

- plasma-inflated mini-magnetosphere at a lunar–solar L1-style location;
- dual REBCO seed coils, ~5 m radius / ~6 m stack, ~0.05–0.1 T, ~0.4–0.8 MA-turns, 30–50 K, order 10–50 MJ stored energy;
- 6–12 helicon/ECR Ar/He plasma sources, roughly 5–50 kW each, target 10^19–10^21 ions/s during boost;
- whistler/lower-hybrid RF current drive: ~0.2–1 MW typical, ~0.5–2 MW peak;
- 3–5 MWe fission/Brayton power;
- ~8–12 MWth waste heat and ~1,500–2,500 m² radiator scale;
- ~300–600 kWe Hall electric propulsion;
- ~23–40 t dry, ~25–46 t wet depending configuration;
- target boundary field ~70–150 nT against ~1–10 nPa solar-wind pressure, with a nominal few-lunar-radius plasma standoff;
- argon hold/boost consumption ~0.1–2 kg/day in the century-buildout notes.

Century program recovered from the same discussion:

| phase | interval | recovered scale |
|---|---|---|
| I Foundations | 2025–2040 | simulations and demos; ~$4–7B CAPEX |
| II Pathfinder | 2040–2052 | ~100–300 kWe class, ~100–300 km cavity; ~$2–4B CAPEX |
| III Block-0 | 2052–2065 | 3 MWe fission, 0.5–1.0 MW RF, 100–300 kW plasma, 100–300 kW EP, dual 5 m HTS coils; 23–35 t dry / 26–40 t wet; ~$12–20B CAPEX, ~$150–300M/yr OPEX |
| IV Block-1 / redundancy | 2065–2085 | 5 MWe, RF 1–2 MWe peak, EP 300–600 kWe, hot spare; ~$6–10B per upgraded unit, ~$10–16B fleet addition, ~$250–500M/yr two-unit OPEX; argon recycling/ISRU pilot |
| V maturity / sustainment | 2085–2125 | reactor swaps ~10–12 y; ~$2–4B refit per 10-y cycle; ~$200–400M/yr OPEX; ≥30% argon ISRU share, cislunar EP propellant |

Recovered century roll-up: **$45–75B P50 / $70–120B P90**, with operations roughly **$15–29B** and EP propellant ≤~100 kg/year in those notes.

Why it was superseded as the main atmospheric shield: it acts on charged particles, while the current retention problem is dominated by direct high-energy photons/upper-atmosphere energetics. It remains potentially relevant as a solar-wind subsystem or experiment.

## 2025-09-29 — EL1-SH / Earth–Sun L1-L2 super-hub metascreen — HISTORICAL

Conversation-recovered optical concept:

- EUV-focused metascreen around Earth–Sun L1;
- old target ~40% EUV attenuation while retaining ≥98% visible/IR;
- materials discussed: Mo/Si multilayers, Zr/Si filters and SiN_x windows;
- transparent-PV / waveguide harvesting capped at ~1–2% of total solar irradiance so the screen remained mostly transparent;
- free-flying 5–50 m tiles, 50–300 m semi-rigid rafts and kilometer-scale veil assemblies;
- 1–10 km² raft units inside ~1,000 km² logical veil blocks;
- UHMWPE Hoytether primary lattice, Vectran/CNT companion elements, ALD Al2O3/SiO2 protective jacket and CNT conductor wraps;
- optical/RF mesh and laser backhaul;
- tile photon-sail tabs plus electrospray/FEEP control;
- halo stationkeeping discussion around ~1–3 m/s/year;
- L1↔L2 and L1↔Moon/Venus logistics used invariant-manifold / low-energy transfer concepts, with scheduled insertion/capture rather than treating L1/L2 as rigid points.

The old ~40% EUV rejection target is far weaker than the 2026 atmosphere-retention design and is retained only as a materials/control lineage.

The associated BOM is reconstructed separately in `recovered_historical_specs/EL1_SH_BOM_recovered.csv`.

## 2025-09-27/29 — peak L1 hub infrastructure — HISTORICAL

A recovered LaTeX/report discussion titled variants of **“Earth--Sun L1 Hub at Peak: Power, Industry, and Cislunar Commerce”** used a mature-hub scale near **300 TW**. Conversation-recovered values included roughly 6.9×10^5 km² photovoltaic collector area, ~0.6–1.5 Gt array mass, ~3.1 MN photon force, ~52 GW electric-propulsion trim, ~9 TW electronics heat and ~1.7×10^4 km² radiator area. Industrial output was discussed at ~0.5–2 Gt/year, with very large autonomous workcell and docking counts.

The original LaTeX artifact was not recovered in the Library in this session. These values are included as historical infrastructure-scale context, not current shield requirements.

## 2026-04 — aeronomy / filter placement branch — HISTORICAL, partly superseded

Known code artifact: `lunar_atmosphere_sim.py` (original bytes unrecovered). Conversation history describes it as a 1-D spherical hydrostatic N2/O2/Ar model with tunable EUV/FUV/X-ray shielding, exobase and Jeans-escape calculations, an old `shield_factor=0.03` reference, and an optional `sweep_shield_factor()` study. Its earlier 250/260 K lifetime framing was later corrected by the September coupled work.

April engineering estimates treated ~99.7–99.9% high-energy attenuation (optical depth ~7) as a useful target and explored 1–100 g/m² shield areal mass. Lunar-disk area was ~9.5×10^12 m², implying ~9.5×10^9–9.5×10^11 kg at those areal masses before expanded atmospheric coverage. A maintenance scaling used roughly `3 kg/s × area_multiplier × (sigma / 1 g/m²) × (replacement_rate / 1%/yr)`.

Placement geometry derived the same basic finite-Sun relation later used in September: shield radius grows approximately as `R + 0.00465 d`. Examples in that discussion included ~3,500 km-class diameter close to the Moon, ~3,570 km at 10,000 km, ~3,940 km at 50,000 km, ~4,400 km at 100,000 km and ~6,700 km near a ~340–350 thousand km lunar-sunward scale, for lunar-disk protection only. The September model superseded these with a 3-lunar-radius atmospheric target and explicit force/mass closure.

## September 4, 2026 — integrated feasibility baseline — SEP-2026 ARCHIVE

`Lunar_Terraforming_Model.zip` and `Lunar_Terraforming_Feasibility.pdf` established the spherical atmospheric inventory and early optical/trajectory comparisons. The model already warned that a screen near Earth–Sun L1 does not automatically follow the Moon, and that screen mass alone can hide a severe stationkeeping/propellant problem.

The complete package is preserved in `raw_sources/` and extracted under `extracted/Lunar_Terraforming_Model/`.

## September 9, 2026 — concrete protection design — SEP-2026 ARCHIVE / current foundation

This is the pivot to the current design: 78,000 km actively held Moon-following titania/silica aperture, explicit power/propellant mass feedback, four regional anchored magnetic installations, component-level maintenance and manufacturing budgets, and explicit separation of photon and charged-particle jobs.

The complete original package is in `extracted/Lunar_Protection_Model/`.

## September 20, 2026 — environment-research import — HISTORICAL/CURRENT INTERFACE

`Terluna_environment_research.zip` and `.patch` preserve atmosphere spectral-interface, thermal-column, climate-cycle and biosphere long-night code plus research screening tables. They are useful because they show the shield being treated as an environmental boundary condition rather than an isolated engineering object.

## September 24–25, 2026 — spectral and atmospheric coupling — CURRENT

The `restructure/lanes` branch added:

- spectral products for the stored titania/silica stack and idealized 200/220/230/240/310 nm edge filters;
- climate response to shield-filtered sunlight;
- a middle-atmosphere O-H-N photochemical model used to study ozone and upper-air temperature under those filters;
- exact short-wave evaluation of the stored stack from hard X-rays through EUV/FUV.

These changes move the central question from “what nominal shield factor should we assume?” to “what spectrum does the physical aperture actually deliver, including defects, and what atmosphere results?”
