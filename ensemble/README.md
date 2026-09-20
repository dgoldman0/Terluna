# Constructing and Sustaining an Open Moon

Five-paper ensemble by Daniel S. Goldman. **Current stage: all five working seeds selected for continued development.**

The root uses seed-04; the four companions use seed-02. The author accepted the companion seeds on 20 September 2026 with further revision expected as research proceeds. Their text and working evidence notes are preserved unchanged. Full scholarly clearance remains pending.

The root develops a scientifically grounded dream of a living, inhabited Moon. Four standalone companions examine its physical environment, biosphere, human life, and engineering. Biological viability and ecology belong together; the human paper develops ways of living, culture, settlement, and use of aerial space. Each companion develops its own whole-world argument.

## Start here

Read the [editorial charter](planning/Open_Moon_Editorial_Charter_control-01.md), then the [detailed ensemble outline](planning/Open_Moon_Ensemble_Outline_plan-01.md). The [paper manifest](planning/paper_manifest.json) records the permanent short IDs and complete UUIDs. [AGENTS.md](AGENTS.md) supplies workspace instructions for later assisted work.

| Paper ID | Working title | Current seed |
|---|---|---|
| `cef466e6-d9f8` | [Constructing and Sustaining an Open Moon](papers/cef466e6-d9f8/README.md) | [seed-04](papers/cef466e6-d9f8/seed-04.md) |
| `7c83e7a3-5955` | [The Physical Environment of an Open Moon](papers/7c83e7a3-5955/README.md) | [seed-02](papers/7c83e7a3-5955/seed-02.md) |
| `74a946e2-8f98` | [Establishing and Sustaining a Lunar Biosphere](papers/74a946e2-8f98/README.md) | [seed-02](papers/74a946e2-8f98/seed-02.md) |
| `babdd5f6-4920` | [Living on an Open Moon](papers/babdd5f6-4920/README.md) | [seed-02](papers/babdd5f6-4920/seed-02.md) |
| `155b039f-a0e8` | [Engineering and Industry for Lunar Terraforming](papers/155b039f-a0e8/README.md) | [seed-02](papers/155b039f-a0e8/seed-02.md) |

## Layout and provenance

- `planning/` preserves all 14 files from the original plan-01 package byte for byte, including the locked charter, outline, source-control originals, hashes, review status, and twilight diagnostic.
- `papers/<paper-id>/` provides one workspace per paper with its selected seed and current-stage README. The root's import record is retained; the engineering workspace also contains its seed's calculation diagnostic.
- `seed_records/companion-seeds-02/` preserves the companion source and review records, arithmetic, and exact-file import map. See its [README](seed_records/companion-seeds-02/README.md).
- `tools/validate_setup.py` checks the planning snapshot, IDs, workspace links, and diagnostic arithmetic using only the Python standard library.
- `repository_setup.json` records the original planning import and repository base revision.

The archived editorial originals retain examples from their earlier projects for provenance. The Open Moon charter defines their application here. Historical planning and seed-package status files describe their original rounds; current workspace status is recorded in the READMEs and import records.

The existing [atmosphere simulator](../atmosphere/sim/sim1.py), [simulator report](../atmosphere/sim/report.md), and [training data](../training_data/) remain separate, unchanged project history. Their presence does not establish source admission or scientific validation for this ensemble. Inspect and reconcile their assumptions before use.

## Check the workspace

From the repository root:

```sh
python3 ensemble/tools/validate_setup.py
python3 ensemble/planning/twilight_diagnostic.py
```

A successful check establishes file integrity, identity consistency, workspace linkage, and the recorded arithmetic only. Scientific source admission, manuscript review, and final-PDF inspection remain separate gates.

## Continuing the work

Preserve the five roles and permanent IDs. Working titles and outlines can evolve through recorded decisions. Keep plan-01 and control-01 intact; introduce versioned successors for authorized revisions. Follow the complete charter before research, drafting, or editing. The five seeds are accepted working drafts, subject to further revision. Full manuscripts and paper PDFs remain pending.
