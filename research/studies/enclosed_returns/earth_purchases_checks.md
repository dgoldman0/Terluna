# Cumulative purchase verification

The purchase account is a decomposition of existing commercial scenarios. It introduces no new orders, refinery assumptions, prices, customer budgets or financial returns. The source economics revision is `2b7df218` on the existing `domain/humanities` branch. The input manifest records complete SHA-256 identities for orders, pricing, acceptance rules, the market producer and revised financial results.

The runner completed **81 alternatives** with checkpoints at years 100 and 200. Cumulative and annual invoices reconcile to **486 financial checkpoints**, covering all three financing modes without adding their duplicate sales. The saved result reproduces exactly. Independent forward acceptance and demand-limit calculations reproduce the premium quantities recovered from invoices.

The current repository check selector was loaded from main revision `4059f8bb` and directed to the humanities worktree using `/tmp/terluna-earth-purchases/check.py`. Its printed plan covered eight changed paths and selected the purchase tests and layer check. The command was:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_ADDOPTS='-p no:cacheprovider' OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python /tmp/terluna-earth-purchases/check.py
```

Result: **7 tests passed**; **1,693 files checked with zero layer violations**. The tests cover full saved-product reproduction and financial reconciliation, rejected commodity offers, the propellant example, non-overlapping mass categories, inclusive checkpoint dates, premium acceptance and customer limits, gradual versus mature-rate cumulative accounting, and a missing checkpoint error.

The report was read continuously against the input definitions and generated totals. Its quantity tables were compared with the saved result; all displayed physical and monetary entries match their stated rounding. Goods, propellant and services remain separate, and the text identifies the missing elemental allocation. The existing source and return products remain unchanged. No rendered PDF, new assay review, market validation or engineering simulation was performed for this arithmetic extension. Recording these checks changes no tested implementation or numerical inputs.

Result SHA-256: `1721b41981f53c3b587c6be64b3c0c518644af4375a59c4822753cb57141c1f4`.

Reproduction within the humanities checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m research.studies.enclosed_returns.earth_purchases --output /tmp/earth_purchases.json.gz
PYTHONDONTWRITEBYTECODE=1 python -m pytest research/studies/enclosed_returns/tests/test_earth_purchases.py -q
```
