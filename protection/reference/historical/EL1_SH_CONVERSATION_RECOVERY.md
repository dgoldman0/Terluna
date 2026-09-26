# EL1-SH / Earth–Sun L1–L2 Super Hub — conversation recovery

**Status:** conversation-recovered technical record. Original 2025 draft / LaTeX / BOM file bytes remain unavailable. The numbers below are preserved from the September 27–29, 2025 project conversations and the later source-recovery pass.

## Identified original titles

- **Earth–Sun L1 Super Hub (EL1-SH) — Concept Draft v0.1**
- **Century Plan v0.1**
- variants of **Earth--Sun L1 Hub at Peak: Power, Industry, and Cislunar Commerce**

[Two private conversation identifiers listed here as search keys were removed on import to this repository, 2026-09-26.]

## Optical / metascreen branch

- EUV-focused metascreen around Earth–Sun L1.
- Historical target: ~40% EUV attenuation while retaining ≥98% visible/IR transmission.
- Recovered effective-area figure: ~3.8×10^13 m².
- Materials discussed: Mo/Si multilayers, Zr/Si filters, SiN_x windows.
- Transparent-PV / waveguide harvesting kept small; recovered BOM discussion used a default harvest fraction `c = 0.001`, while other notes capped harvesting around 1–2% of total irradiance to preserve optical transparency.
- Free-flying 5–50 m tiles; 50–300 m semi-rigid rafts; kilometer-scale veil assemblies.
- 1–10 km² raft units inside ~1,000 km² logical veil blocks.
- UHMWPE Hoytether primary lattice; Vectran/CNT companion elements; ALD Al2O3/SiO2 protective jacket; CNT conductor wraps.
- Optical/RF mesh and laser backhaul.
- Photon-sail tabs plus electrospray/FEEP control.
- Halo-stationkeeping discussion: order 1–3 m/s/year.
- L1↔L2 and L1↔Moon/Venus logistics used invariant-manifold / low-energy-transfer concepts with scheduled insertion and capture.

## Recovered magnetic / plasma module inventory

The existing reconstructed BOM in this bundle preserves the following counts and descriptions from the 2025 conversation:

| Item | Recovered count | Recovered description |
|---|---:|---|
| `MS-VEIL-1k` | ~38,000 | ~1,000 km² logical EUV-veil block; ~1,000 t each in the old note |
| `MS-CTRL-50` | ~800 | control node for roughly 50 veil blocks |
| `SW-SEED-10GJ` | ~40,000–80,000 | segmented REBCO HTS seed-coil module; order 10^10 J; local tens-to-hundreds of µT |
| `SW-PLAS-50/200MW` | ~5,000–12,000 | 50–200 MW RF/ECR plasma-drive module with H/Ar injection |
| `SW-MINI-NODE` | 6–12 | larger mini-magnetosphere / plasma node |
| `GC-RING-STN` | 200–400 | historical GCR-corridor ring-station family |
| `GC-LUN-LOBE` | 200–600 | historical lunar-lobe magnetic station |
| `PW-MW-5GW` | 500–2,000 | 5 GW-class microwave power/beaming module |
| `PW-LAS-1GW` | 1,000–5,000 | 1 GW-class laser power/beaming module |
| `TH-RAD-10M` | 1,000–3,000 | large thermal-radiator module family |
| `RB-MICRO` | 10^8–10^9 | micro-robot / inspection-repair swarm element |
| `RB-RAFT-TUG` | 50,000–150,000 | raft tug / deployment and servicing vehicle |

Recovered geometry/operations notes also used an Earth–Moon corridor of radius ~2,000 km and length ~400,000 km, a GCR-corridor field target around 50–100 µT, and a lunar “umbrella” concept aiming for >95% coverage with a 99% storm mode. Those are historical design targets, not demonstrated shielding performance.

## Peak hub scale recovered from the 2025 report discussion

A mature-hub discussion used a scale near **300 TW** and included roughly:

- photovoltaic collector area ~6.9×10^5 km²;
- array mass ~0.6–1.5 Gt;
- photon force ~3.1 MN;
- electric-propulsion trim ~52 GW;
- electronics heat ~9 TW;
- radiator area ~1.7×10^4 km²;
- industrial output ~0.5–2 Gt/year.

A separate full L1–L2 super-hub cost model in that conversation gave a mid-case of about **$15–28T (2025 USD)**, discussed in the context of a broader **$20–40T** program scale.

## Relation to the current architecture

The reusable pieces are the modular construction logic, servicing / logistics architecture, optical-swarm control concepts, and local particle-protection staging. The old ~40% EUV target and magnetic/plasma assumptions are historical and were superseded by the 2026 coupled atmosphere/protection work.
