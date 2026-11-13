# Verification of the complex returns appraisal

12 November 2026. These checks establish source identity, reproducible arithmetic and internal consistency. They do not establish product qualification, market demand, procurement sufficiency or a global welfare valuation.

## Numerical and structural checks

- `python3 -m pytest -q research/studies/complex_returns/tests/test_returns.py`: **16 passed**. Controls include saved-product reproduction and snapshot hashes; an independent annuity formula; agreement with the provisioning external balance; exclusion of GDP from cash; the two-year growth pipeline, capacity ceiling and cap; marginal expansion accounting; loss of GDP and cash when demand binds; separate financial and socioeconomic switching values; reserves as cash restrictions rather than economic damage; and complementary capital-provider claims. The incremental commercial comparison is positive while the corresponding whole-campaign comparison remains negative.
- `python3 /tmp/check-complex-returns.py`: selected layer check and the 16 study tests passed. This wrapper imports the current main `shared/check_changed.py` and sets its root to the isolated humanities checkout. The humanities branch predates the current changed-file Makefile; no branch tooling was replaced.
- After moving the renderer into the visualization lane, `python3 shared/check_layers.py` checked **1,643 files with zero violations**. The final selector plan still selects the study tests and layer check. It identifies the renderer as uncovered by automated tests; its actual rendering and visual review are recorded below.
- Final JSON parsing, local-link checks and `git diff --check` passed. The record and final prose do not change tested numerical inputs.

## Source and substantive review

The two copied provisioning products are byte-identical to the versions at `146a558275f486c9fd9a5cf79908b7c6fdf21888`. Their original repository paths, complete SHA-256 identities and schema versions are recorded in `returns_inputs.json`. They are frozen input snapshots, not copied physical model implementations. A bounded audit also verified the ten identities consumed by the upstream economy model; that is source inspection, not a rerun of its numerical suite.

One bounded research review examined primary appraisal methods and another audited the current provisioning account. Subsequent accounting reviews checked the actual implementation. Their findings led to removal of initial commercial resource reserves from the economic cost proxy, separation of the domestic growth invoice from overlapping labour and inputs, and correction of GDP diagnostics for lost sales and construction output. The study was then read substantively in context, including baseline selection, investment scopes, stakeholder standing, counterfactuals, tables, growth, long horizons and conclusions. This is internal review, not independent scholarly peer review.

`returns_sources.json` records targeted source-reading scope and the unsuccessful World Bank/EIB retrievals. Those unavailable documents are not used as evidence. No full literature review or terrestrial nonmarket price transfer is claimed.

## Figure review

`MPLCONFIGDIR=/tmp/terluna-roi-mpl XDG_CACHE_HOME=/tmp/terluna-roi-cache python3 visualization/complex-returns/plot.py` produced the PNG, SVG and hash manifest in ignored `research/runs/complex_returns/figures/`. The PNG was inspected at full figure scale for readable titles, axes and legends. Curves and bars were checked against the reviewed JSON: population plateaus, 100-year funding NPVs, proportional versus capped-demand external balances and socioeconomic switching values. Both image formats carry product and renderer hashes. Generated images are not committed.

## Remaining verification limits

The full repository suite was not run because the new calculation reads frozen data and has no physical-model consumers. Upstream engineering tests were not rerun. No PDF manuscript was typeset or reviewed. The actual procurement schedule, asset-cohort replacement timing, empirical demand, negotiated capital claims, complete resource opportunity costs and measured nonmarket outcomes remain open. The centuries-long scenarios are not validated forecasts or an extension of UK discounting guidance.
