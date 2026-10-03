# Wind-wave comparison

`plot.py` reads `climate/waves/results/waves.json`, schema
`terluna.climate.wind-waves/1`, and displays its Earth and lunar gravity
experiments at 100 km fetch. The three columns use the three wind magnitudes
in the product. The top row shows significant wave height and the bottom row
shows mean period. Wind duration runs from 1 to 48 hours on a logarithmic axis.
Lines connect the sampled outputs without introducing intermediate model runs.

Solid lines show the AGROW sensitivity patch, which scales an arbitrary
wind-wave onset frequency with gravity; dashed lines show stock SWAN. Both
use uniform steady wind, flat 1000 m water and initially calm conditions,
with incoming swell and currents set to zero. Wind percentile labels select
forcing magnitudes. The source product records the numerical checks and
remaining physical uncertainties.

Early growth is sensitive to the timestep. The numerical resolution checks
cover the 48-hour lunar endpoint at the middle wind magnitude. Earlier outputs
retain the reported timestep sensitivity; other wind cases need further checks.

Run from the repository root after producing the climate data:

```sh
python visualization/waves/plot.py
```

The generated `results/waves.png` and `results/waves.json` stay outside Git.
The sidecar records the product and renderer hashes, image hash, source
evidence, numerical checks and exact plotted rows. The renderer checks the
schema, units, complete plotted time series and the stated scenario conditions.
Wave calculations come from the climate data product.

## Basin and slope pilot

`python visualization/waves/pilot.py` reads
`climate/waves/results/pilot.json` and writes `results/pilot.png` with a
provenance sidecar. The two maps show the actual coarse SWAN nodes at hour 12;
the lower panels show the offshore node's hourly values and the three assumed
slope profiles. Each map cell displays its computed value. Source and image
hashes are recorded, and the displayed arrays are checked against the product.
The climate domain's `results/pilot_checks.json` records the resolution
differences and their implications for geographic predictions.

## Coastal resolution and wind coverage

`python visualization/waves/coastal.py` reads the coastal, CM1 ring-coverage
and GCM coverage products. It writes `results/coastal.png` and
`results/wind_coverage.png`, each with a hash sidecar. `--coverage-only` draws
the recovered atmospheric records independently of the coastal simulation.

The coastal maps display native wave-grid values and differences at common
wet nodes. The coverage figure shows the three CM1 paths, one equatorial
surface-wind history, and the GCM's final saved regional mean at its original
model level. Captions distinguish their time sampling and height. Array checks
compare plotted values with the products; source hashes identify the inputs.

## Waves through changing weather

`python -m visualization.waves.weather` reads
`climate/waves/results/weather.json` and writes `results/weather.png` with a
hash sidecar. It displays the recovered month of basin-mean atmospheric
stress, the selected week's wave-height map and histories, and mean period
at the offshore reference point. The input-speed curve is labelled as the
stress-equivalent SWAN input. Array checks compare the plotted values exactly
with the climate product. The climate study records timestep, prior-history
and remaining geographic sensitivities.

## The second lunar cycle

`python -m visualization.waves.cycle` reads
`climate/waves/results/cycle.json` and writes `results/cycle.png` with its
hash sidecar. Maps show the second cycle's largest local significant height
and the share of time above 1 m. Histories show basin and reference-point
heights, the area above four height thresholds, mean period at the reference
point, and the atmospheric surface stress supplied to the basin.

The climate product defines the spin-up and reporting intervals. Threshold
durations integrate its hourly heights; the renderer uses those supplied
statistics. Array checks cover every plotted map and curve. The caption
states the single-cycle scope and the remaining weak-wave and coastal
resolution questions.

## Energy approaching the shore

`python -m visualization.waves.shore` reads
`climate/waves/results/shore.json` and its hashed coastal case and terrain
products. It writes `results/shore.png` with a provenance sidecar. The three
maps show native flooded-rock depth, computed significant wave height and
net wave-energy transport for the selected eastward phase of the second
solar cycle. Arrows retain the computed transport direction and relative
magnitude, with the latitude–longitude geometry accounted for.

The renderer checks all three displayed arrays against their products.
Station labels connect the maps to the domain report's separate incoming
and outgoing powers. The caption states the stationary weather phase and
the geographic and physical assumptions; the domain report retains the
refinement and boundary checks.

With all three weather phases present, the command also writes
`results/shore_phases.png` and its sidecar. These maps share a power scale
and the same 237 m wave grid. Equal-length arrows show the local direction
of net wave-energy transport; the colours compare its magnitude.

## Changing coastal waves and shoreline motion

`python -m visualization.waves.coastal_history` reads the domain's
`terluna.climate.coastal-history/1` product and its hashed array archive.
`results/coastal_history.png` shows the complete reported cycle of station
heights and incoming power. `results/coastal_history_maps.png` shows the
six-hourly mean and maximum heights across the coast. The captions retain
the distinct sampling intervals, first-cycle spin-up and terrain assumptions.

`python -m visualization.waves.runup` reads the selected cases from
`climate/waves/results/runup.json`. `results/runup.png` displays the actual bed
and computed water surface near each case's largest recorded run-up, followed
by its complete reported waterline history. The upper panels use the nearest
one-second surface sample; the lower panels retain the 0.1 s waterline record.
Rock and assumed beach geometries carry separate labels.

Both renderers verify the displayed arrays and write source, renderer and
image hashes beside the ignored figures. The domain products own the
calculations and the finite-record interpretation of the run-up percentiles.
