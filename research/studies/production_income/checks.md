# Verification of the year-200 production and income study

The accounting consumer and its study were checked on the existing `domain/humanities` branch against the revised financial/wealth baseline at `d06172c0`. The provisioning source remains commit `a8af87614e9c30b6fa4879440765f7e74fe094be`. These checks establish source identity, numerical implementation and reproducibility; they do not establish empirical probabilities or qualify a two-century development programme.

## Source and account checks

The extractor verified SHA-256 and schema identity for all six admitted financial, physical, commercial and ownership products. It matched each commercial orderbook to its physical case and matched the 18 detailed asset accounts to both their physical-case and recovery-intervention hashes. All 81 cases contain one year-200 row at the same population, 6,468,289.238 under the source's continuous planning convention.

The annual cash plan reproduces the saved revised gross receipts and operating cash. For the 63 cases without the saved detailed finance-flow export, the unchanged financing engine reconstructs the equity control and U-credit path; year-200 debt is checked against the committed checkpoint. This preserves supplier-held U interest as outside earned income and separately retains local interest. The maximum external interest represented at year 200 is 27.866 billion USD2041; the local-interest total is zero in these cases.

The asset approximation uses the 18 available cohort exports. Its median depreciation/replacement ratio is 1.275156; leaving each case out gives relative errors from −11.870% to +13.115%. This is an internal comparison against engineering estimates, not empirical validation of economic depreciation.

The admitted SNA document's local bytes matched SHA-256 `47fa283670e5f85bffc8df7e44f41db45e0bb2f00cead758db6435e709999d1f`. Targeted source reading and its boundaries are recorded in [income_sources.json](income_sources.json). No full-manual reading or outside dataset redistribution is claimed.

## Numerical controls and tests

The saved run uses 8,192 draws per case, or 663,552 paired GDP/GNI evaluations. Across these draws, the largest production/expenditure difference is **0.015625 USD2041**, the largest production/income difference **0.00390625 USD2041**, and the largest GNI allocation difference **0.0078125 USD2041**. These are floating-point differences on trillion-dollar accounts.

The focused suite passed **12 tests**. It covers imported capital cancelling between investment and imports, in-kind compensation not creating profit, depreciation remaining inside gross income, owner distributions and restored-core value not creating earned income, signed foreign losses, outside income and ownership boundaries, scenario weights, source tampering, case coverage, independent account identities, saved-product identities, and separation of portfolio choices. An initial exact-equality assertion for a trillion-dollar capital increment encountered floating-point rounding; its replacement checks the increment to 0.01 USD2041 absolute and 1e-12 relative tolerance. The underlying accounts required no numerical correction.

The current repository check selector was run with its root set to this existing worktree:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python /tmp/terluna_scope_check.py
```

That wrapper uses main's `shared/check_changed.py` and changes only its `ROOT`, because the older worktree Makefile predates the current scoped checking procedure. The selected checks were the layer check, historical provenance reproduction, the new test file and the affected shared lane-rule test. **All 13 selected tests passed**; the layer check inspected **1,689 files with zero violations**. Historical provenance reported integrity PASS for **111 retained files**, with illumination, feasibility and protection numerical checks PASS. Full repository verification was not requested or run.

The manifests use study-specific names so the dependency selector resolves the new inputs to this consumer. An initial generic `inputs.json` filename caused ambiguous basename references to unrelated models; the plan was inspected before execution, and the final selection above follows the actual dependency scope. No selector implementation was changed.

## Sampling, figures and reading

An independent run used seed 2041201 and **16,384 draws per case**, twice the saved sample size. Across 648 comparisons of the 10th, 50th and 90th percentiles—three priors, nine portfolios, two selection conditions and four GDP/GNI measures—the maximum relative difference was **0.6756%** and the median difference **0.1088%**. This sampling check supports the displayed numerical precision; uncertainty in priors and institutions remains substantially larger.

Both generated PNG figures were inspected at full rendered size. The nine-portfolio figure shows separate GDP and ownership-dependent GNI intervals; the diversified comparison shows cumulative probability under each prior. Labels, units, selection conditions and probability qualifications are visible. All **54 plotted percentile triplets** were compared with the saved product and matched exactly. The SVG files and source-hash manifest are generated alongside the PNGs. The final metadata-only result update was rendered again with the same numerical values.

The study was read in full for its scope, table interpretation, account definitions and treatment of uncertainty. The reading checked the distinction between the intended financing strategy and source-path selection, the institutional assumption behind in-kind benefits, depreciation inside gross measures, signed foreign losses, actual ownership versus counterfactual endpoint shares, and the absence of a completion probability. This was substantive author-facing review by the implementing agent; no independent scholarly review or forecast calibration is claimed.

## Artifact identity and remaining work

The extracted annual input has SHA-256 `0b587c16ae83ed5a31d7fa5ef9a2b4f95200a6b85818b1bf7a031d905e707c3a`. The saved distribution product has SHA-256 `506fe8fe881333f6a3cd3484107f3047bee635909a52dcfad552e835d2081f84`. The result records producer, input, evidence and interpretation metadata. Figures are derived artifacts in ignored `research/runs/production_income/figures/` and are not committed.

Unperformed work includes empirical elicitation or calibration of the scenario probabilities, a full asset-cohort export for all 81 cases, identification of resource-provision recipients, national accounts for portfolio-equity ownership, complete fiscal and inventory accounts, returns on external financial assets, and a stochastic investment/financing history. Geological and engineering qualification remain with their producing domains. These omissions bound the interpretation of the completed conditional distribution study.
