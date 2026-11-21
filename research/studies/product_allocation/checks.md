# Product allocation verification

This record covers a research-informed product and processing-location screen. It does not establish a qualified refinery, an optimized global production mix, a sales commitment, extraction permission or revised total returns.

## Sources and admitted records

The resource assays, dated physical feed and recovery assumptions are copied from resources revision `84886850`; the local/beamed power comparison is taken from provisioning revision `adeaa8ea`. The existing purchase account and manufacturing comparisons are admitted from humanities revision `a63cd71b`. Full file identities appear in `allocation_inputs.json` and `inputs/basis.json`.

`export_inputs.py` checks each upstream file against its stated Git revision before admitting it. Re-exporting to `/tmp/terluna-product-routing/reproduced-inputs` produced exact byte matches for all four saved resource/power input records. The exporter reads data and source documents; it does not import another lane's model implementation.

Direct HTTPS retrieval obtained the intended USGS Mineral Commodity Summaries 2026 report and rare-earth chapter after the web reader returned HTTP 403. Their complete SHA-256 identities are in `allocation_sources.json`; the earlier retrieval failure was resolved without changing source year or substituting a report. Selected price, production, purity, lot-size and unit notes were read. The platinum-group and rare-earth price tables were additionally inspected as rendered PDF pages. The local text extractor's font-weight warnings are recorded; the selected numerical values were legible in the source tables.

NASA's optical-fibre article and manufacturing slides, ESA's regolith electrolysis and metal-printing accounts, and NASA's research-to-Earth-production example were read for the claims used. The 2025 NASA presentation reports are not treated as independent industrial replication. The separate NASA materials survey was downloaded as background and supplies no numerical input or claim of complete scholarly review. Full source PDFs remain outside the repository; the repository contains selected factual transcriptions and provenance.

## Calculation and tests

The runner completed 21 product screens, 54 concentrate-refining comparisons, nine fabrication mass comparisons, 36 power/energy comparisons, eighteen joint-circuit cost ceilings and eight annual feed-treatment opportunities. The opportunities retain the original dated feed and do not create commercial receipts. Missing assays and unpriced elements remain explicit; oxide prices are not multiplied by elemental assay mass.

The repository's current dependency selector was loaded from main revision `4059f8bb` and directed to the existing humanities checkout with `/tmp/terluna-earth-purchases/check.py`. Descriptive study filenames resolve the selector's otherwise ambiguous suffix matches to unrelated files named `inputs.json`, `sources.json` and `model.py`. The final plan covered sixteen changed paths and selected the allocation tests and layer check.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_ADDOPTS='-p no:cacheprovider' OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python /tmp/terluna-earth-purchases/check.py
```

Result: **10 tests passed** in the selected allocation suite; **1,695 files checked with zero layer violations**. The tests reproduce the saved results; independently reconcile the Earth and off-Earth refining cash expressions; verify a recovery-driven preference reversal; check added fabrication mass, absent assays and oxide units; distinguish the inherited freight floor from cheaper future sensitivities; verify joint cost accounting and annual feed scaling; and retain the existing customer purchase totals.

Thirteen reported price benchmarks, six fabrication penalties, four refining thresholds and four joint-cost figures were compared with the saved inputs and results at their displayed precision, using conventional half-up rounding at the $452.50 threshold. All matched. Annual potential outputs, production-scale comparisons, power costs and the year-200 customer target were also checked against the generated records.

The report was read continuously as a product-selection argument, including the location table, empirical evidence, numerical comparisons, allocation sequence and preservation conditions. Its recommendations remain distinct from qualified manufacturing performance and from adopted investment decisions. This was a Markdown and source review; no PDF of the new report was rendered. No full repository suite, new trajectory simulation, metallurgy experiment or customer survey was performed. The later edits record these checks and refine prose, without changing numerical implementation, inputs or results.

Result SHA-256: `ba7cee408c1d98143a92cc008b0240ae1547fbe5c1b5ed9dc6bea4d7aafa28fe`.

## Reproduction

From the humanities checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m research.studies.product_allocation.allocation --output /tmp/product_allocation.json
PYTHONDONTWRITEBYTECODE=1 python -m pytest research/studies/product_allocation/tests/test_allocation.py -q
```

Rebuilding the imported records also requires checkouts containing their recorded revisions:

```sh
python research/studies/product_allocation/export_inputs.py \
  --resources-root /path/to/resources-checkout \
  --provisioning-root /path/to/provisioning-checkout \
  --humanities-root /path/to/humanities-checkout \
  --output /tmp/product_allocation_inputs
```

A changed or missing upstream revision fails admission rather than silently substituting a source. New cost, market or recovery evidence should be admitted explicitly, with the affected comparison and financial integration rerun.
