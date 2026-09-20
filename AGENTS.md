# Terluna research and paper workspace

The repository contains shared subject research alongside the five-paper manuscript ensemble. Start with `research/README.md`, `research/status.json` and the relevant topic README. `ensemble/AGENTS.md` and its complete editorial charter govern manuscript work.

## Current direction

Constructing and Sustaining an Open Moon is the next full manuscript. Support it with several regional, vertical, biological and human possibilities; the root's opening landscape is one example. Preserve exactly five paper roles and their established UUID identities unless the author changes them. The biosphere paper owns biology and ecology together; the human paper owns potential human life and society.

## File and evidence practice

Use ordinary stable filenames and Git history for new working material. Keep the already locked `ensemble/planning/` snapshot and selected seed files intact. The author has authorized the top-level research organization, so shared research belongs in the relevant folders outside `/ensemble`.

Historical imported model implementations and reference tables are byte-pinned in `research/provenance.json`. Preserve them as reproducible baselines; new scientific development should be explicitly distinguished and tested. Do not copy the same monolithic model into several topic folders.

A runnable calculation, numerical convergence, an input hash, a proposed design and empirical validation are different evidence states. State assumptions and boundaries near outputs. Never fill an empty topic with invented results or mark future work complete. Never silently substitute data after a failed download or hash mismatch.

External spectral bytes are restored through `protection/fetch_inputs.py`; they are ignored by Git. Review redistribution rights before adding outside datasets or full papers. Private credentials and unrelated personal information never belong here. Legacy synthetic training records are not empirical evidence.

## Checks

Run `python research/check.py` after changes affecting this organization or its imported files. It reports missing protection inputs explicitly; `--require-inputs` makes that state a nonzero exit. Results normally go to ignored `research/runs/`. Run the existing `python3 ensemble/tools/validate_setup.py` when paper workspace structure or links change.

Report actual checks and unperformed work. Full scholarly source reading, substantive editorial review, and exact-build PDF review remain separate from numerical tests. Changes outside the current task require author authorization.
