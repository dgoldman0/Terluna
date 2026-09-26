# Industrial architecture of 19 September 2026, reconstructed

A reconstruction of the industrial design worked out with the author on
19 September 2026. The recovery pass in the author's other project space wrote
these files on 2026-09-26 from the conversation record; they are imported
verbatim from `Terluna_industry_energy_logistics_recovery_2026-09-26.zip` (see
[research/provenance.json](../../../research/provenance.json)). No file from
19 September itself survives. The documents mark which statements they recall
from the conversation and which quantities they recomputed, and every number in
`results/` was computed on 2026-09-26 by `calculations.py`. It is input for the
integrated industrial model the engineering domain lists as next work.

| File | Holds |
|---|---|
| [RECOVERED_INDUSTRIAL_ARCHITECTURE.md](RECOVERED_INDUSTRIAL_ARCHITECTURE.md) | The reconstructed architecture: power ledgers, bootstrap and power sequence, heat rejection, source industry, manufacturing closure, the packet "matter stream", cislunar capture and buffering, traffic safety, construction order, and a specification for the network model |
| [INDUSTRIAL_ARCHITECTURE.json](INDUSTRIAL_ARCHITECTURE.json) | The same reconstruction in structured form |
| [HISTORY_RECOVERY_NOTES.md](HISTORY_RECOVERY_NOTES.md) | Dated notes recalling the 19 September discussion |
| [CALCULATION_PROVENANCE.md](CALCULATION_PROVENANCE.md), [calculations.py](calculations.py), [results/](results/) | Formulas and assumptions, a dependency-free script, and its outputs: fusion and fission fuel at about 500 PW, collector and radiator areas, packet rates and kinetic energies, arrival heat |
| [EVIDENCE_LEDGER.csv](EVIDENCE_LEDGER.csv) | The provenance of each claim: archived executable, recalled design history, or recomputed |

`python engineering/reference/industrial_architecture/calculations.py` rewrites
`results/` in place; on import it reproduced every file byte for byte.

The author confirmed the six traffic-safety rules in section 8 as decisions on
2026-09-26, and [research/decisions.md](../../../research/decisions.md) records
them on that basis. K≈1.17 (about 500 PW) is a civilization-scale scenario, not
the Open Moon's own load. Source choices follow the conservation study's later
principle: common icy bodies, with Saturn's rings, Titan's atmosphere and the
ocean worlds left untouched.
