# Figures for the founding complex returns study

The renderer reads the reviewed [returns product](../../research/studies/complex_returns/results/returns.json) and shows population paths, funding NPV, external balances and socioeconomic switching values. It introduces no new economic model.

```sh
MPLCONFIGDIR=/tmp/terluna-roi-mpl XDG_CACHE_HOME=/tmp/terluna-roi-cache python3 visualization/complex-returns/plot.py
```

PNG and SVG figures and a source-hash manifest are written to ignored `research/runs/complex_returns/figures/`. Each image includes the source product and renderer hashes. The study states the scope, units and conditional assumptions; figures preserve those limits. Generated images are not committed.

The companion `staged_plot.py` reads the reviewed staged-finance product and
writes `staged.png`, `staged.svg` and `staged_manifest.json` to the same ignored
figure directory. It shows service invoices, industrial capacity, capital-recovery
NPV and sensitivity cases. Product and renderer SHA-256 identities are embedded;
all scenario assumptions belong to the producing study.

`strategy_plot.py` reads the integrated export-strategy product and produces
`strategies.png`, `strategies.svg` and `strategies_manifest.json` in the same
ignored directory. It compares whole-enterprise NPV, incremental operating
contributions, appraisal horizons and imposed enclosed-cap operation.

`consolidated.py` reads the [consolidated return product](../../research/studies/complex_returns/results/consolidated.json) and verifies its producer and every declared dependency hash before plotting. It compares the latest U package's sponsor cash NPV across four horizons, middle and strong demand, and real discount rates of 3% and 7%. The plotted account deducts initial seed and monetary-reserve funding, includes the conditional annual release of interest savings to outside saving assets, and excludes terminal book claims. Reproduce it from the repository root with:

```sh
MPLCONFIGDIR=/tmp/terluna-returns-mpl XDG_CACHE_HOME=/tmp/terluna-returns-cache python3 visualization/complex-returns/consolidated.py
```

The PNG, SVG and `.manifest.json` files use the ignored output stem `research/runs/complex_returns/consolidated-returns`. Both images embed the product, producer, renderer and dependency identities; the manifest also records every displayed source value and the image hashes. The renderer imports no economic model and performs only selection and display-unit conversion. The figure states the conditional evidence and keeps sponsor returns distinct from global welfare.
