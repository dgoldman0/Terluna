# Verification of material potential and reclamation

This record distinguishes numerical implementation, source inspection and engineering qualification. The study is a conditional resource-distribution and physical-reclamation assessment; no permitted reserve, operating mine or fully closed industrial ecology has been demonstrated.

## Inputs and source review

The full JPL snapshot, MITHNEOS table, NEOMOD orbital grid, prior Lunar Prospector product, cost coefficients and committed provisioning reclamation export are checked against their recorded SHA-256 identities. Raw external data remain outside Git. The reclamation producer is provisioning commit `210a20bd2cf16ec951f02a01f6ee3758b497df03`; its copied compressed export has SHA-256 `f6b43f00ea585af509e0f9e4ec3545ccf1956ed8c04348dd3404cf965bb54f2c`.

Source reading was bounded to the sections recorded in `population_sources.json` and `lunar_potential_sources.json`. It includes the population-count and albedo/composition methods, the USGS meteorite analogue fits and printed program, and Hayne's density formulation and constraints. The rendered USGS program page was inspected; its apparent percentage-conversion and density-expression defects are recorded as printed-code issues, without claiming to reproduce the original execution. The USGS 2023 lunar-report PDF remained unavailable with HTTP 403; only the obtained official publication abstract informs the study. No numerical claim is derived from the unobtained report. SWRI grid retrieval did not validate the TLS chain in this environment; the paper-linked public source and exact bytes are recorded, and the restoration tool requires the pinned hash.

## Numerical controls

Population checks cover mutually exclusive known/unknown contributions, retained science exclusions, missing-characterization treatment, paired nested pools, monotone population branches, weight-percent conversion, physical composition bounds and provenance. The producer retains 2,048 draws and records the compound-moment approximation's numerical error separately from geological uncertainty. Its worst estimated single-body water-mean numerical standard error is 5.57%; no claim of asymptotically exact tail probabilities is made.

Lunar checks compare the analytic density integral with quadrature, verify spherical-area fractions, preserve six-oxide covariance and assay nonclosure, and check source hashes, depth factors and partial cost/grade ordering. The regional sampling retains 16,384 draws. An initial integer-weight array failure in the synthetic sampling test was corrected by explicitly using floating-point weights before the producer was run.

Comparison checks cover recovery exactly once, raw-feed versus product units, finite-sample frequency intervals, paired joint thresholds, exact dated demand and seed inclusion, stock-equivalent sustaining years, source identities, all 36 cases at three checkpoints, and availability monotonicity. The 81-case baseline ranges remain distinct from the six orderbooks rerun under six reclamation variants. A whole-pool grade sensitivity is checked independently of the inherited selected 10% source grade.

The provisioning producer reports 107 passing tests and one optional skip, with zero layer violations. After prose/provenance cleanup, its 11 focused tests passed again, including fresh reproduction of a complete combined case. Six control cases reproduce the frozen physical baseline exactly. These checks establish numerical accounting and the stated equipment/workforce screens, separately from engineering validation of recovery or restored service life.

## Review and figure

A second reader checked the stock/grade distinction, marginal versus joint probability thresholds, lunar density/area units, depth scope, partial cost scope and baseline/combined demand figures. The completed synthesis was read in context against the products. The resource summaries point to the broader population assessment while preserving the earlier finite-catalogue calculations.

The three-panel figure is rendered by `visualization/source-access/potential.py` into ignored `research/runs/source_access/potential/`. It displays model water quantiles, annual fresh feed and conditional stock-sufficiency frequencies. Its manifest records input and renderer hashes and the plotted water quantiles. The PNG was visually inspected for units, legibility, legends and the distinction between model intervals, the additional-water-allowance ceiling and the theoretical full-inventory access point.

## Selected repository checks

The existing resources checkout predates the change-selected Makefile target. The current main-branch `shared/check_changed.py` is run through `/tmp/terluna_resource_scope_check.py`, setting only its repository root to this resources worktree. This applies the current selector without changing unrelated tooling on the domain branch. The selector chose the layer check and all 14 resources test files, including the new population, lunar-potential and comparison tests. It identified the renderer as requiring its documented run and visual inspection.

Command: `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python /tmp/terluna_resource_scope_check.py`. Result: **145 passed, one skipped** in 15.86 seconds; **1,643 files checked with zero layer violations**. The skip is the opt-in bound-water regeneration test; its environment switch was not enabled. The final renderer was run separately and its PNG inspected. Study-local links and `git diff --check` also passed. No numerical source, input or test changed after this selected run.

Full-repository verification is not requested for this bounded domain change. Full scholarly reading of every cited paper, empirical industrial qualification, recurring cargo trajectories and an integrated financial rerun remain separate work.
