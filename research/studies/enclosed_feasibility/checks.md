# Engineering feasibility: verification record

Assessment and checks: 13 November 2026, existing `domain/provisioning` worktree, based on revision `a8af8761`.

## Calculation and provenance

The producer ran with one numerical thread:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/studies/enclosed_feasibility/feasibility.py
```

It verified all six pinned inputs before computing and wrote `results/feasibility.json`. The product records its producer hash and input identities. The calculation extracts requirements from the existing staged physical product; it does not rerun or alter that producer. All 93 case identities and their existing admission results are retained: 81/81 primary, 6/6 poor-source and 4/6 faster-rollout cases. All nine portfolios have year-200 scenario ranges. Their ranges are not probability intervals.

The physical source product, enclosure assumptions, reclamation inputs and shared constants remain unchanged. Resource-distribution and preservation source identities are separately recorded in `feasibility_sources.json` at their existing resources revision. They inform the review; this study does not perform new source-release evaluations.

## Automated checks

The current repository affected-check selector was run against the provisioning worktree. The installed wrapper, `/tmp/terluna-feasibility/scoped_check.py`, loads `shared/check_changed.py` from the main checkout and sets its root to the existing provisioning checkout. It implements the current scoped-check procedure without changing the older branch's build files.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python /tmp/terluna-feasibility/scoped_check.py
```

Result: **13 tests passed**, comprising the study's ten tests and three selected check-selector regression tests. The layer check examined **1,668 files with zero violations**. JSON/Python structural validation passed. The selector reported the renderer as lacking automated coverage; its documented rendering and visual inspection were performed separately.

The numerical controls cover ideal circular-aperture encircled power at the first diffraction zero; inversion and aperture/range scaling; thermal balance and fourth-power temperature scaling; traffic probability limits; repair's all-new and material break-even limits; rejection of changed input hashes; retained case counts and checkpoints; delivered-power and transmitter heat balances; launch mass and storage energy; pressure-force accounting; and the distinction between opening demand, construction averages and operating peaks. The coefficient audit and controls preserve the original meaning of 200 W/m² as delivered electricity and reproduce the existing 72-hour reserve after discharge efficiency and the protected margin. These tests establish the behaviour of the bounded calculations, not empirical validation of their engineering assumptions.

## Source and artifact review

The source register records the actual reading scope. Sources include NASA's reported regolith experiment, CaRD prototype, life-support baseline tables, ISS water recovery, dust mitigation, LOXSAT preparation, a partial-gravity paper abstract, the OSIRIS-REx sample announcement, Caltech's mission report, a DOE development description and a vendor design-review announcement. Programme targets and proposed missions remain distinguished from completed demonstrations. This was a targeted primary-source review, not a complete scholarly review of every subsystem or all papers behind the linked reports.

The report's main numerical tables were compared with the generated product. Substantive review checked the U boundary; startup versus operating requirements; capacity margins already present in the producer; the difference between cargo-loss illustrations and planetary-impact risk; recovery versus complete ecological closure; and the dependence of repair savings on restored life. The report preserves architectural decisions and identifies recommendations as proposals. No financial returns, currency definitions or source permissions were changed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python visualization/enclosed_feasibility/plot.py
```

The renderer produced PNG and PDF figures plus source/renderer hashes under ignored `research/runs/enclosed_feasibility/`. The PNG was inspected at full rendered size for legibility, units, legend placement, plotted values and evidence captions. The source-product hash in the rendering record matches the current result. The PDF was exported but did not receive a separate exact-PDF review. Generated figures remain outside Git; their renderer is committed.

The repository-wide full suite was not requested or run. New trajectory integration, structural finite-element analysis, radiation-dose calculations, whole-ecosystem validation, flight hardware trials, endurance qualification and a dated initial deployment simulation remain unperformed. Their absence is part of the study's feasibility conclusion, not a failed retrieval or a silently substituted result.
