# Figures for the founding complex returns study

The renderer reads the reviewed [returns product](../../research/studies/complex_returns/results/returns.json) and shows population paths, funding NPV, external balances and socioeconomic switching values. It introduces no new economic model.

```sh
MPLCONFIGDIR=/tmp/terluna-roi-mpl XDG_CACHE_HOME=/tmp/terluna-roi-cache python3 visualization/complex-returns/plot.py
```

PNG and SVG figures and a source-hash manifest are written to ignored `research/runs/complex_returns/figures/`. Each image includes the source product and renderer hashes. The study states the scope, units and conditional assumptions; figures preserve those limits. Generated images are not committed.
