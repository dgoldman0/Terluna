# Commits

Use informative commit messages. Multiple phrases for a commit message.

# System Usage Limits

Be mindfulof system resource availability esp. when other tests are running in other sessions as well.

# Terluna research and paper workspace

The repository holds shared research in subject domains, cross-domain studies, scientific visualization, an immersive experience and the five-paper manuscript ensemble. Start with `research/README.md`, `research/status.json` and the relevant domain README. `ensemble/AGENTS.md` and its complete editorial charter govern manuscript work.

## Where work goes

Decide the lane before writing anything:

| Work | Home |
|---|---|
| A model, simulation or reference calculation | Its domain: `atmosphere/`, `climate/`, `biosphere/`, `geography/`, `geology/`, `illumination/`, `protection/`, `engineering/`, `habitation/`, `humanities/`, `provisioning/`, `resources/` |
| A study that couples several domains | `research/studies/<study>/`, with its runner, results and write-up together |
| Rendering that shows what a model computed (plots, labs, reference renderers, viewers) | `visualization/` (see `visualization/AGENTS.md`) |
| The explorable world: world definitions and the exports the Unreal game reads (the game itself is github.com/dgoldman0/terluna-game; the web engine here is frozen) | `immersion/` (see `immersion/AGENTS.md`) |
| Constants, scenarios and the data-product convention | `shared/` |
| Deprecated material | `archive/` |

Science never lives in the experience. If immersion work needs a physical quantity no domain provides, add it to the domain, with tests, and have the experience read its product. Results cross lanes as data products (a versioned schema, the producer and model hash, an evidence statement and the rule for reading it), not as imports. `shared/README.md` lists which lanes may import which; `shared/check_layers.py` enforces it. Constants come from `shared/constants.json`, not new literals.

## Current direction

Constructing and Sustaining an Open Moon is the next full manuscript. Support it with several regional, vertical, biological and human possibilities; the root's opening landscape is one example. Preserve exactly five paper roles and their established UUID identities unless the author changes them. The biosphere paper owns biology and ecology together; the human paper owns potential human life and society.

Immersive environments come from the domain research (biology, ecology, geology, hydrology) as it matures, not from the manuscripts' narrative scenes. Until then the experience builds general capability and uses placeholders that are labelled as placeholders.

## File and evidence practice

Use ordinary stable filenames and Git history for new working material. Keep the already locked `ensemble/planning/` snapshot and selected seed files intact.

Historical imported model implementations and reference tables are byte-pinned in `research/provenance.json`. Preserve them as reproducible baselines; new scientific development should be explicitly distinguished and tested. Do not copy the same monolithic model into several topic folders.

A runnable calculation, numerical convergence, an input hash, a proposed design and empirical validation are different evidence states. State assumptions and boundaries near outputs. Never fill an empty topic with invented results or mark future work complete. Never silently substitute data after a failed download or hash mismatch.

External spectral bytes are restored through `protection/fetch_inputs.py`; they are ignored by Git. Review redistribution rights before adding outside datasets or full papers, and record their source, hash and credit. Private credentials and unrelated personal information never belong here. Legacy synthetic training records are not empirical evidence.

## Checks

For ordinary local changes, run the checks relevant to the changed work before committing. `make check` selects them from staged, unstaged and untracked files relative to `HEAD`; `make check-plan` explains the selection without running tests. For work already committed, use `make check CHECK_ARGS='--base <revision-before-the-work>'`. A clean working tree alone does not verify earlier commits.

The selector follows Python imports, file references and recorded data-product dependencies into downstream consumers. Prose and planning edits receive structural checks and any tests that read them. Model and input changes receive affected tests; an unrecognized dependency falls back to the owner's suite, and shared foundations or unknown repository-wide changes select the full suite. Review the printed scope: dynamic paths and external data can require additional checks. `--include <paths>` adds affected inputs to the Git changes. The [shared checking guide](shared/README.md#checks) describes the commands and limits.

Use `make check-full` for substantial branch integrations, changes whose effects cannot be bounded, or an explicitly requested complete verification. It retains the layer check, all configured Python and JavaScript tests, the calendar build and tests, historical provenance reproduction and ensemble validation. A local commit is not by itself a reason to repeat that suite. `research/check.py` reports missing protection inputs explicitly; `--require-inputs` makes that state a nonzero exit. Study runs write to ignored `research/runs/`.

A successful check may cover later commits while its relevant source, tests, inputs (including ignored data), dependencies and execution environment remain unchanged. Record the command, scope and result; recording those results or splitting checked work into commits does not itself require another full run. The selector does not cache passing test results, so it cannot silently reuse a pass after an external input changes.

Report actual checks and unperformed work. Full scholarly source reading, substantive editorial review, and exact-build PDF review remain separate from numerical tests. Changes outside the current task require author authorization.
