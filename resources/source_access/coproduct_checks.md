# Coproduct assessment: evidence and checks

## Source and input review

The resource-code audit found a bulk account for water, native iron and major lunar oxides, without explicit precious-metal and rare-earth output. The extension retains that calculation and adds a separate element account. The original population, lunar-potential, reclamation and source-access products remain unchanged.

Six public PDFs were retrieved, hashed and read at the scopes recorded in `coproduct_sources.json`. The CI abundance table was checked against its mass-concentration column; the KREEP table uses its “This Work” column. The retained deviations describe the published compilation or analytical procedure, not mine-grade probability. The iron and chondrite examples from Day et al. use ng/g in the source and are converted to ppm by dividing by 1,000. Only selected La and Sc rows from the lunar-soil table are used. Full-table visual inspection, full scholarly review of all cited original measurements and independent laboratory validation were not performed.

The NASA abundance PDF was retrieved successfully after the initial web failure. The Cannon et al. full text remained inaccessible (403 from publisher and mirror); indexed abstract/introduction material supplied context only, and no numerical inputs depend on it. The unresolved reading limit is retained in the source register. External full papers are not committed.

An independent comparison with the committed provisioning product at `adeaa8ea` verified the compressed product hash and **all 1,458 physical values** in the 81-case, three-checkpoint export. The source producer and its earlier source revision are recorded separately from the checkout containing the product. Each exported row uses six annual physical fields. Fresh feed is distinct from retained tails, total generic metal and virgin generic metal.

The two named metal-rich asteroid records were extracted from the retained JPL snapshot without changing its bytes. Both fail the existing wider geometry limits in semimajor axis and eccentricity. The protection evaluator leaves their release unresolved, and the existing Bennu, Ryugu, Itokawa and Eros holds remain active. A site-release calculation supplies no mining authorization.

## Numerical verification

Run from the resources worktree:

```sh
python -m resources.source_access.coproducts
```

The calculation writes `results/coproducts.json.gz`, including 243 physical-case checkpoints. Eleven focused tests cover:

- ppm-to-tonnage conversion and conservation of an element across untreated feed, stream capture, separation residue, finishing residue and refined output;
- zero output for absent treatment or failed recovery, without deleting contained inventory;
- rejection of invalid masses, grades and recovery fractions;
- assay transcription controls, ng/g conversion, disjoint platinum-group and REY totals, the NdPr subset, and missing measurements remaining unknown;
- common geological-transfer factors, stage recovery applied once and monotonic treatment-share scaling;
- incremental-circuit break-even arithmetic;
- exact saved-product reproduction, producer and input identities;
- all annual outputs using fresh feed once and retaining the inherited feasibility flags;
- square-kilometre conversion against the original lunar column-mass product;
- named-target orbital screening and preservation holds.

An initial focused run identified an incorrect hand-entered test expectation for reference year-200 platinum-group output. The independent calculation gives 46.005230885 tonnes/year; the expectation was corrected from 46.006430873. The production calculation and saved result did not change.

## Selected repository checks

This existing resources checkout predates the current change-selected Makefile target. The main checkout's `shared/check_changed.py` was run through `/tmp/terluna-coproducts/scoped_check.py`, with its root set to the resources checkout. The wrapper changes no repository tooling.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python /tmp/terluna-coproducts/scoped_check.py
```

The selector chose the layer check and `resources/source_access/tests/test_coproducts.py`. Result: **11 passed**, **1,643 files checked with zero layer violations**. Pytest reported one non-failing cache-write warning because the sibling checkout is read-only under the ordinary sandbox. The result log is `/tmp/terluna-coproducts/selected-check.log`. Final prose and check-record edits do not change the checked model, inputs, tests or dependencies.

The completed report was read in context against its numerical result and source qualifications. Current resource summaries and the research status index now point to the element-specific assessment. Study-local file links, the serialized research index and `git diff --check` were checked separately. These activities establish accounting, source identity and document consistency; they do not establish an extractable reserve, a functioning refinery, a statistical asteroid-grade distribution or incremental commercial returns. Full-repository verification, coupled financial reruns and rendered-publication review were not performed for this bounded resource extension.

## First-eight-year Earth-delivery extension

The extension uses the existing **12,500-resident opening at year 8**. Its exporter copies 81 combined-recovery seed accounts from the same pinned provisioning product. An independent comparison verified all **486 physical values**: total, base and late-procurement asteroid/lunar feed. Each total equals its base and late components; year-8 operating flows are excluded. No provisioning model or accepted growth schedule changed.

The assumed annual feed and refinery ramps, Earth allocation, fractional delivery lags and return loss are recorded in `startup_delivery_inputs.json`. Seven new tests check once-only seed allocation; fractional-year arrival and complete mass balance; absent refining, zero export allocation and arrival after the cutoff; preservation of the opening population and feed; refinery-capacity scaling with independent lunar selection; platinum-group totals; and exact saved-result reproduction with input hashes. The full calculation was run, generating 243 case/coverage comparisons.

Command: `PYTHONDONTWRITEBYTECODE=1 PYTEST_ADDOPTS='-p no:cacheprovider' OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python /tmp/terluna-coproducts/scoped_check.py`, with the existing wrapper setting the checker root to the resources checkout. The selector chose `resources/source_access/tests/test_startup_deliveries.py` and the layer check. Result: **7 passed**, **1,646 files checked with zero layer violations**. The log is `/tmp/terluna-startup-coproducts/resources-check.log`. No numerical code, input, test or dependency changed after this run.

The new report section was read against the result, including cumulative versus annual units, distinct refinery scales, hypothetical lunar feed selection, retained material and delivery after year 8. Study-local links, JSON parsing and whitespace checks passed. The USGS nickel price is contextual evidence from the search-indexed primary table; direct full-PDF retrieval returned 403. That reading limit remains recorded in the startup inputs, and no prospective commodity-price or revenue model uses the observation.

This is a timing and mass-accounting estimate. Source deployment, the peak throughput implied by the assumed feed ramp, coproduct refinery capital and staffing, added feed from their construction, flight dispatch, product specifications and buyer demand have not been solved jointly. The existing local-power alternatives also remain outside the inherited feed calculation. These qualifications apply to every delivery estimate; passing the numerical tests establishes neither a funded export programme nor a financial return.
