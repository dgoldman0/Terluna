# Instructions for work under /ensemble

Read `planning/Open_Moon_Editorial_Charter_control-01.md` in full, then `planning/Open_Moon_Ensemble_Outline_plan-01.md`, `planning/paper_manifest.json`, and the relevant paper workspace before making substantive changes. The charter is authoritative; this file is an entry point.

## Scope

Maintain exactly five papers unless the author explicitly changes the structure. Preserve the UUID-derived paper IDs across titles and revisions. The root is the scientifically grounded dream paper. The biosphere companion owns organism-level viability and ecology together, including relevant human biological constraints. The human companion owns life, culture, settlement, and spatial/temporal design. Engineering remains one companion within the whole-world ensemble.

The author has authorized the root-paper seed and selected revision 04 for repository inclusion. Its exact conversation artifact is preserved at `papers/cef466e6-d9f8/seed-04.md`; read the workspace README and `seed-import-01.json` for current status. The author has also accepted revision 02 of all four companion seeds as working starting points, with further revision expected. Each companion workspace contains `seed-02.md`; current status and exact-file provenance are in its README and `seed_records/companion-seeds-02/import.json`. Preserve each companion's independent whole-world argument rather than narrowing it to the root seed's illustrative setting. Further drafting or expansion follows subsequent author instructions. Do not fill empty workspaces with invented abstracts, results, citations, or review approvals.

## Provenance and change control

Treat `planning/` as the immutable imported plan-01 snapshot. Record later authorized changes in versioned successor documents with their scope and reasons; retain the original hashes and historical review statuses. Project-specific examples in `planning/source_controls/` are provenance, not additional lunar claims or instructions overriding the charter.

Legacy files elsewhere in this repository, older chat-derived material, and model outputs require inspection and reconciliation before scientific use. Keep scenario assumptions, model boundaries, evidence status, and dependencies explicit. A newer claim or successful numerical run does not by itself establish physical validity.

## Writing and evidence

Apply all G1-G7 and E1-E4 requirements and the ten separate prose passes. Review argument order, inheritance, audience, cadence, precise action, living/human presence, affirmative development, redundancy, and the opening-to-conclusion arc. Avoid defensive framing and reader-management. Preserve technical meaning while revising style.

Archive complete original sources and record full readings and claim locators before final citation. Distinguish sourced observations, calculations, conditional models, designs, cultural possibilities, and fictional archive material. Do not transfer old source-clearance or review labels into new work.

Record actual review activity and exact artifact hashes. Leave unperformed checks DEFERRED or BLOCKED. Never weaken a gate or prefill PASS to make a release appear complete. Final PDF review applies to the exact rendered build.

## Validation and repository boundaries

Run `python3 ensemble/tools/validate_setup.py` after workspace changes. This checks the imported snapshot and structure, not scientific or editorial adequacy. A future build/release checker must preserve that distinction.

Keep task changes within `/ensemble` unless broader edits are requested. Preserve the existing simulator, reports, and training data. Do not publish credentials, private account data, or source documents without appropriate redistribution rights.
