# Recovery and provenance status

This record separates accessible calculations from historical descriptions. Earlier whole-project archive reviews remain conversation artifacts; their model-run claims are not inherited as new test results. The current [check record](checks.json) records this persistence pass.

| Material | Recovery condition | Disposition |
|---|---|---|
| April retention simulator and report | Already in `atmosphere/sim/`; code recovered and inspected in the previous audit | Preserve untouched. Parameterized energy/escape assumptions need reconciliation with September. It is no longer an unrecovered-code item. |
| September feasibility package | Original ZIP recovered again on 2026-09-26 inside the solar-shield dump; its hash matches | Verbatim model and all nine CSV tables imported, plus selected original reference fields in JSON. The full [report](baselines/feasibility/report.md) and its four figures were imported verbatim on 2026-09-26; the PDF builder stays in the source archive. Every member is indexed. |
| September protection package | Original ZIP recovered again on 2026-09-26 inside the solar-shield dump; its hash matches | Verbatim model/verifier and compact tables imported, plus selected original reference fields in JSON. Five external inputs pinned and restorable. The full design report was imported on 2026-09-26 as [protection/report.md](../protection/report.md), verbatim except for two private conversation links removed from its first source. Large grids and the PDF are indexed but omitted. |
| Solar-shield recovery dump, 2026-09-26 | Prepared in the author's other project space; all 136 members indexed | Holds both September archives, extracted copies, recovery notes and conversation-recovered reconstructions. The reconstructions of Lunashield-L1 (2025), the Earth–Sun L1 super hub EL1-SH (2025) and the protection lineage are imported verbatim to [protection/reference/historical/](../protection/reference/historical/). They are labelled reconstructions from conversations; their original files remain unrecovered. |
| September 2025 four component drafts and earlier root | Recovered through prior chat descriptions; complete manuscript files remain unlocated in this pass | Research questions and source leads, with historical scenario assumptions. No complete-source reading is claimed. |
| July 2026 Terluna knowledge bundle, `Terluna_Moon_Terraforming_Knowledge_Bundle_v1.0.zip` (435,651 bytes; recorded SHA-256 `cce7a34d101a741ffce4cb7ff8958a6a80c7d4772ba57eac4879dc989442a799`) | Located as a Project ZIP; raw-byte materialization was blocked again for the 2026-09-26 dump | Unread. Its creation conversation describes 52 files: decisions and corrections, atmosphere and climate, hydrology, soil, biosphere, settlement, infrastructure, failure analysis, 85 questions and answers, calculations and sources. Its 80 kPa scenario differs from the current 1.2 atm target. The highest-priority recovery target. |
| Earlier April 2026 atmosphere script, `lunar_atmosphere_sim.py` | Described in conversations; bytes not found in the Library or GitHub | A 1-D spherical N2/O2/Ar atmosphere with a tunable EUV/FUV/X-ray shield factor (reference 0.03), exobase and Jeans escape, and a shield-factor sweep. The recovered April simulator in `atmosphere/sim/` has no shield-factor sweep, so this appears to be a separate, earlier file. |
| Original EL1-SH and super-hub documents (2025) | Titles and values known from conversations; LaTeX, BOM and report files unrecovered | *Earth–Sun L1 Super Hub (EL1-SH), Concept Draft v0.1*, *Century Plan v0.1*, and variants of *Earth–Sun L1 Hub at Peak: Power, Industry, and Cislunar Commerce*. Reconstructed values are in [protection/reference/historical/](../protection/reference/historical/). |
| `plasma_scaling_corridor_vs_region.csv` (2025) | Described in conversations; file unrecovered | A corridor (~2,000 km radius, ~400,000 km long) against a regional plasma volume (~200,000 km radius); reconstructed notes are in [protection/reference/historical/](../protection/reference/historical/). |
| September megaforest/aphotic discussion and images | Earlier archive review recovered concepts and inspected illustrations | Proposed dimensions and ecosystem arrangements, not measurements or executable validation. |
| Legacy training data | Existing `archive/training_data/` files | Synthetic/fictional and correction records; not empirical or simulation evidence. |
| Accepted root and companion seeds | Already committed under `ensemble/papers/` | Working paper prose. Their selection does not admit their provisional sources or validate scenarios. |

## Known conflicts retained for research

Pressure/composition cases differ; the 250/260 K universal lifetime cliff is superseded; filter rejection percentages can use incompatible spectral definitions; optical blocking and magnetic particle effects need separate treatment. Angular twilight arithmetic is not illuminance. Large air volume is not supported habitable capacity. Low-gravity mechanical opportunity is not complete life-cycle evidence.

Older transport budgets and oxygen-process estimates mix different system boundaries unless explicitly reconciled. Thermal waste must be located where it is deposited. Small annual replacement can accumulate into a large resource commitment. Sagan K≈1.17, preferred cargo packets and cultural practices remain conditional scenarios.

## Original source archives

- `Lunar_Terraforming_Model.zip`, Library file `file_0000000085d082309e5433f1aa02c37e`.
- `Lunar_Protection_Model.zip`, Library file `file_00000000155881f59331ebf972a6add8`.
- Unrecovered bytes: `Terluna_Moon_Terraforming_Knowledge_Bundle_v1.0.zip`, previously listed Library file `file_00000000387071f58518961fba45ebbf`.
- `Terluna_solar_shield_full_dump_2026-09-26.zip`, 5,513,169 bytes, SHA-256 `e8d5fa05800f6917f27e5813546fe639f5be8562466b84cdc9d659daaba10bea`, held by the author; it contains both archives above.

Library identifiers are provenance, not portable download URLs. [provenance.json](provenance.json) contains independently calculated archive hashes, member sizes and dispositions. No original third-party research papers are bundled or newly admitted by this organization pass.
