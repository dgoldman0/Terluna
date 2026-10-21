# Snow analysis figures

`plot.py` reads the versioned atmosphere snow population, turbulence and optics
products and draws a four-panel scientific overview. It imports no domain model
and performs only selection, unit conversion and normalization of stored moments.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/plot.py
```

The figure compares arriving-mass distributions for three imposed snowfall
rates; number, mass and extinction shares in the heavy case; resolved vertical
motions in two archived full-cycle weather boxes; and direct and diffuse light
at 35 km. The population panels use the corrected-dimension, representative
dendrite-area scenario at seed 0 and the 100/L input flux normalization. The
wind bars describe conditioned cell-times, not flake trajectories or resolved
small-scale turbulence. Direct normal and diffuse horizontal illuminance remain
separate because their receiving-plane conventions differ.

The renderer requires all three current products. It checks their schemas,
checks that the optics product hashes the population being plotted, and compares
stored total moments with the sums behind each cumulative curve. An open-ended
overflow bin remains in the normalization and manifest; the figure annotates
its share and draws no invented finite upper size for that tail. The inherited
particle algorithm has unresolved numerical sensitivity; these distributions
are exploratory, including their large-size tail. Extinction
shares equal projected-area shares under the optics product's size-independent
extinction efficiency; the renderer checks the corresponding optical total.

`generated/overview.png`, `generated/overview.svg` and
`generated/manifest.json` are ignored by Git. Image metadata records the source
and product hashes. The manifest additionally records every plotted value,
selection rule, upstream evidence statement and output hash. Rebuild after a
product change; the chart itself does not select a preferred population law or
supply missing results.
