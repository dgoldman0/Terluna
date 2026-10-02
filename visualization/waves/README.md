# Wind-wave comparison

`plot.py` reads `climate/waves/results/waves.json`, schema
`terluna.climate.wind-waves/1`, and displays its Earth and lunar gravity
experiments at 100 km fetch. The three columns use the three wind magnitudes
in the product. The top row shows significant wave height and the bottom row
shows mean period. Wind duration runs from 1 to 48 hours on a logarithmic axis.
Lines connect the sampled outputs without introducing intermediate model runs.

Solid lines show the AGROW sensitivity patch, which scales an arbitrary
wind-wave onset frequency with gravity; dashed lines show stock SWAN. Both
use uniform steady wind, flat 1000 m water, initially calm conditions and no
incoming swell or currents. The comparison is a numerical experiment.
Wind percentile labels select forcing magnitudes; the resulting waves do
not acquire those percentiles. Consult the source product for numerical
checks and the remaining physical evidence limits.

Early growth is sensitive to the timestep. The numerical resolution checks
cover only the 48-hour lunar endpoint at the middle wind magnitude; they do
not establish convergence of the other endpoints or the entire time histories.

Run from the repository root after producing the climate data:

```sh
python visualization/waves/plot.py
```

The generated `results/waves.png` and `results/waves.json` stay outside Git.
The sidecar records the product and renderer hashes, image hash, source
evidence, numerical checks and exact plotted rows. The renderer checks the
schema, units, complete plotted time series and the stated scenario conditions.
It imports no climate model and performs no wave calculation.
