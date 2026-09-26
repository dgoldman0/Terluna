# Recovered industrial architecture (19 September 2026)

The industrial design worked out with the author on 19 September 2026, recovered
from the author's other project space on 2026-09-26 and imported verbatim from
`Terluna_industry_energy_logistics_recovery_2026-09-26.zip` (see
[research/provenance.json](../../../research/provenance.json)). It is design
history: the recovered architecture and the arithmetic behind it, labelled in
the files as recovered notes or as quantities recomputed there. It is input for
the integrated industrial model the engineering domain lists as next work.

| File | Holds |
|---|---|
| [RECOVERED_INDUSTRIAL_ARCHITECTURE.md](RECOVERED_INDUSTRIAL_ARCHITECTURE.md) | The architecture: power ledgers, bootstrap and power sequence, heat rejection, source industry, manufacturing closure, the packet "matter stream", cislunar capture and buffering, traffic safety, construction order, and a specification for the network model |
| [INDUSTRIAL_ARCHITECTURE.json](INDUSTRIAL_ARCHITECTURE.json) | The same architecture in structured form |
| [HISTORY_RECOVERY_NOTES.md](HISTORY_RECOVERY_NOTES.md) | The dated design notes from the 19 September discussion |
| [CALCULATION_PROVENANCE.md](CALCULATION_PROVENANCE.md), [calculations.py](calculations.py), [results/](results/) | Formulas and assumptions, a dependency-free script, and its outputs: fusion and fission fuel at about 500 PW, collector and radiator areas, packet rates and kinetic energies, arrival heat |
| [EVIDENCE_LEDGER.csv](EVIDENCE_LEDGER.csv) | The provenance of each claim: archived executable, recovered design history, or recomputed |

`python engineering/reference/industrial_architecture/calculations.py` rewrites
`results/` in place; on import it reproduced every file byte for byte.

The six traffic-safety rules in section 8 are the author's decisions and are
recorded in [research/decisions.md](../../../research/decisions.md). K≈1.17
(about 500 PW) is a civilization-scale scenario, not the Open Moon's own load.
Source choices follow the conservation study's later principle: common icy
bodies, with Saturn's rings, Titan's atmosphere and the ocean worlds left
untouched.
