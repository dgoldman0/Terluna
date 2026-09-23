# Shoreline checkpoint history

Records from the shoreline prototype's development, grouped by checkpoint. They
document what was built and checked at each stage; none of them describes the
current state. The commit that added each checkpoint is the place to reproduce
it exactly.

| Checkpoint | Date | Commit | What changed |
|---|---|---|---|
| [01-prototype](01-prototype/) | 2026-09-21 | `c010049`, `4877d8e` | First-person shoreline prototype built on the month-of-light sky data |
| [02-m1](02-m1/) | 2026-09-21 | `4c7affa` | Continuous land and seabed, optical water depth, reflections, calibrated noon exposure (Three.js r180) |
| [03-landscape](03-landscape/) | 2026-09-22 | `caa4caa` | One landscape field for geometry, collision, materials, placement and water; nested terrain; surface-water ledger |
| [04-r186](04-r186/) | 2026-09-22 | `032765c` | Three.js r180 → r186 migration, connected ponding, WebGPU/TSL renderer lab |
| [05-cloud-columns](05-cloud-columns/) | 2026-09-22 | `3f353e6` | Revision 04: clouds from atmospheric column soundings |

Revision 04's column methods now live with the model in
[atmosphere/column](../../../atmosphere/column/METHODS.md), and its rendering notes in
[../cloud-rendering.md](../cloud-rendering.md). The accuracy program that followed
(A1–A3) is in [illumination/references](../../../illumination/references/), and the
B1 reference renderer in [visualization/reference-renderer](../../../visualization/reference-renderer/).
