# Verification of the enclosed surface power comparison

Assessment date: 13 November 2026. Worktree: existing `domain/provisioning`; starting revision `fd431684`. The work adds a separate power comparison and qualifies the existing engineering synthesis. The inherited physical and financial producers and their numerical products are unchanged.

## Numerical and accounting controls

The final power comparison was regenerated with one numerical thread:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python provisioning/enclosed_power/model.py
```

The product retains SHA-256 identities for its producer, engineering inputs, physical industry product, constants and comparison assumptions. All input hashes matched. The producer rejects a mismatch; no substitute input was used. Product SHA-256: `6bd746a53e1e0046866b048dce929c897962a16cc5e86e9aae48254d7d7746c6`. Producer SHA-256: `9fa048b1e051885904df6c8e6f07d160c2a0c2d487d75d1a305b9b3f22960c4d`.

The fourteen power tests cover the ordinary-night analytical solution; conservation under load shifting; hybrid energy accounting; one conserved hydrogen/oxygen/water inventory and its tanks; thermal retention and local procurement of bulk storage; beam capacity scaling and its aperture floor; persistent inventories and replacement dates; polar numerical resolution; product/input identity; all 81 primary year-200 load mappings; the fission radiator mass check; and rejected invalid options/import fractions.

Across the 18 principal architecture/bundle cases and 24 site/scheduling cases, the largest absolute cyclic residual is **3.49 × 10⁻¹² kWh per kW of mean load**. Increasing the synthetic polar profile from 4,096 to 8,192 samples changes required photovoltaic capacity by **0.0263%** and leaves its calculated storage span unchanged. These checks establish the represented energy/accounting implementation, not the accuracy of a real site's illumination or hardware performance.

## Repository checks

The current scoped selector was loaded from the main worktree and directed at the provisioning worktree through `/tmp/terluna-feasibility/scoped_check.py`. This uses the project's current changed-file procedure without modifying the existing branch's build files. Commands:

```bash
python /tmp/terluna-feasibility/scoped_check.py --plan
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python /tmp/terluna-feasibility/scoped_check.py
```

The selected scope comprises the provisioning domain tests, the new power tests, five affected shared-check tests, lane enforcement and historical provenance reproduction. The source-register JSON has no executable consumer, so the selector conservatively includes the provisioning owner's suite. It identifies the new renderer for manual artifact inspection.

The first run reports **137 tests passed and one skipped**, with no lane violations across **1,671 source files** and historical provenance **PASS**. Its enclosing process exits with status 143 after the test summary; that invocation alone is not recorded as a successful complete command. A repeat with output captured in `/tmp/terluna-enclosed-power/selected-checks.log` resolves the runner status below.

Final selected-check result: **PASS, exit status 0**. The repeated run reports **137 passed, one skipped in 60.55 seconds**, followed by the selector's successful-completion message. Lane and historical-provenance checks pass again. The initial driver interruption is resolved by this complete repeat; its cause was not established.

The skipped test extracts the existing climate simulation's 2.3 GB of output when `TERLUNA_RERUN_EXTRACT` is explicitly enabled. This power study neither changes nor consumes that extraction. No repository-wide full suite was requested or run.

## Sources and substantive review

The [source register](enclosed_power_sources.json) identifies eight primary agency/research sources, the technical sections consulted and what each can establish. NTRS metadata was checked against source titles. The 1991 regolith-storage paper is represented by its abstract, not a claimed full reading. An ESA page retrieval failed and was excluded from evidence; no numerical input depends on it.

The report was read in context against the input manifest and generated tables. Review checked the vacuum/Open Moon boundary, the separate meanings of U and its economic zone, mean versus specified electrical capacity, reserve independence, storage inventories, local thermal-medium mass, imported equipment fractions, future cost assumptions, development bounds and checkpoint stocks versus campaign spending. The engineering synthesis was reviewed to ensure beam-specific qualification is no longer presented as necessary for every enclosed-phase power architecture.

The review retains a substantive objection in the 2006 lunar-fission paper: high-temperature regolith heat distribution and recovery were judged impractical in that assessment. The scalar thermal-retention assumption does not answer that objection. Reactor masses from paper concepts likewise do not demonstrate complete qualified lunar installations.

## Render inspection

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python visualization/enclosed_power/plot.py
```

The renderer generated PNG and PDF artifacts and a hash record in ignored `research/runs/enclosed_power/`. The final PNG was inspected at its full 2,400 × 1,728 resolution. Legends were moved clear of the bars and numerical labels. Beam curves use calculated costs from the producer, and their crossover marks agree with the tabulated thresholds. The PDF was exported but not separately inspected page by page.

## Work required before selecting and financing a utility

No actual terrain/ephemeris dispatch, probabilistic outage calculation, integrated hot-store conduction model, radiation-transport/shielding design, reactor fuel-cycle account, detailed installation estimate or physical hardware trial is supplied by these numerical checks. The solar thermal design especially needs a temperature-dependent recovery and containment calculation. Each supply needs a dated commissioning sequence and a complete cost account.

The study does not rerun the economy or attribute its comparative capacity-cost differences to investor returns. A selected generating/storage design must first replace the surface-supply cohorts and propagate materials, freight, maintenance, reserve requirements and any changed production schedule through the physical and financial producers. Remote source and depot power remain separate obligations. The existing common development allowance is retained once.
