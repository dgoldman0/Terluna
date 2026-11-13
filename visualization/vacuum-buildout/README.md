# Enclosed-build and returns figures

These renderers read reviewed products from `research/studies/vacuum_buildout/`; all modeling and source assumptions remain in that study. `render.py` draws capacity, advance construction, source feed and equipment utilization. `returns.py` draws the strong-demand year-100/year-200 sponsor returns, external-equity requirement, operating deficit and discount-rate sensitivity.

```sh
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/vacuum-buildout/render.py
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/vacuum-buildout/returns.py
```

PNG, SVG and source/renderer hash manifests are written to ignored `research/runs/vacuum_buildout/`. Inspect the resulting figures and compare their values with the frozen products. The plots establish neither sufficient mineral reserves nor committed funding.
