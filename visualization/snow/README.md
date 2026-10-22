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

## Condensation-growth numerical controls

`growth.py` compares the frozen first-order and limited growth-transport
products with their independent continuous-mass characteristic references:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/growth.py
```

Three panels show the cumulative share of remaining ice mass at 96, 144 and
216 mass pivots, with identical logarithmic mass axes. Each grid uses its own
independently grown projected seed population as the reference. The common
shaded band spans the independently grown original continuous seed support;
projecting the initial seeds onto each grid changes that grid's reference
support slightly. The fourth panel shows number-weighted mass variance divided
by squared mean mass. The limited scheme reduces numerical spreading, while
its remaining discrepancy is visible on every grid.

These are closed condensation-only box controls at fixed air conditions. The
132.589 h duration stresses the growth transport numerically; it is not a flake
residence time. The controls include no collisions, settling or new seed
injection and do not predict a lunar snowstorm. Every displayed run uses the
same 30 s time-step cap and the same represented mass domain.

The renderer checks product schemas, the limited product's binding to the
first-order product, upstream source and input hashes, settings, grids and
endpoint times. It reconstructs number, mass, mean mass and variance from the
stored arrays and checks them against stored moments. CDFs are discrete sums at
stored masses, normalized by all remaining ice mass, including the tails outside
the view. The common display window excludes less than 0.01% of every curve;
the exact fractions and the negligible upper-boundary losses remain in the
manifest. The renderer imports no domain code and computes no growth physics.

`generated/growth_comparison.svg`, `generated/growth_comparison.png` and
`generated/growth_manifest.json` are ignored by Git. Image metadata records
input and renderer hashes. The manifest preserves upstream identities, evidence
statements, reading rules, plotted arrays, reconstructed moments and output
hashes. Regenerate the figure when its inputs change.

## Collision-distribution controls

`collision.py` compares an imposed collector capturing small particles from a
replenished bath under the frozen snow kernel:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/collision.py
```

The two panels show number CDFs after 1.36 and 4.54 minutes, with identical axes.
The independent reference samples exact waiting times and retains continuous
masses; the three fixed-pivot grids use matrix exponentials, so a timestep error
cannot explain their broad distributions. The shaded reference band is a 95%
DKW number-CDF band at each displayed time, not a joint two-time confidence band.
This control starts with a 15.86 cm collector; it does not establish the formation
or survival of such a flake.

The renderer verifies the product schema, producer and input identities,
comparison times, collector-number budget and means reconstructed from native
arrays. The common view retains at least 99.99% of each number distribution;
all omitted tails remain in the normalization and are recorded in the manifest.
The code imports no domain physics. Generated PNG, SVG and manifest files use
the `collision_comparison` and `collision_manifest` stems under `generated/`
and remain outside Git. Their provenance records input, renderer and output
hashes, plotted arrays, display exclusions and upstream evidence limits.

## Published Earth profile comparison

`earth.py` compares the pinned author-model reruns with the two published raw
Earth profiles. Supply the directory holding the hash-verified external archive
members; the renderer keeps those bytes outside the repository:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/earth.py \
  --archive-dir /tmp/lunar_snow_earth_population
```

Two rows show deposition plus settling, then the case with aggregation. Four
columns show number and mass concentrations and their velocity moments, with
matching axes between rows. Every line retains the actual 20 m bins and archive
missing-value gaps. The source layer is shaded. The plotted rerun means use 31
snapshots from model hours five through ten, inclusive; the archive's exact
producer and averaging times remain unidentified.

The renderer checks saved producer identities, embedded run/build records, the
shared executable identity and archive-member hashes. It independently
reconstructs the saved profile-difference metrics before plotting. Its manifest
records the plotted values, input and output hashes, library versions and
upstream evidence limits. The comparison reproduces author calculations under
declared settings; it supplies no independent population or empirical
validation. `Fn` and `Fm` are cell-average velocity moments, not measured
boundary crossings. Generated `earth_comparison.png`, `earth_comparison.svg`
and `earth_manifest.json` remain ignored under `generated/`.

## Within-cell collision candidates

`moments.py` compares the candidate collision distributions with an analytic
boundary test and the frozen exact-event bath reference:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/moments.py \
  --entropy-bath atmosphere/snow/results/collision_entropy_bath.json \
  --entropy-boundary atmosphere/snow/results/collision_entropy_boundary.json
```

Omitting both optional paths reproduces the frozen Beta-only diagnostic.

The upper panels show the constant-rate Poisson control when the expected mass
reaches the 191-cell reference face, then its distribution errors before and
after that crossing. Observation times remain fixed relative to this reference
face. The Beta reconstruction retains accurate first two mass moments while
its number CDF differs by 18.4 percentage points near the face. The discrepancy
is a difference in cumulative collector number, not an error of 18.4% in mass or
diameter. Its later agreement does not clear this failed distribution test.
The entropy reconstruction reduces the reference-face CDF error to about
0.185 percentage points in its completed boundary control.

The lower panels show the separate physical snow-kernel bath at 4.54 minutes.
A close-up spans six reference standard deviations on each side of the empirical
mean; the horizontal coordinate is departure from that mean as a percentage of
initial collector mass. The reference spread is about 0.0629% of initial mass.
The broad view retains the much larger spreading of the old fixed-pivot grids.
Every curve keeps the original collector normalization, including tails outside
the displayed windows. Grid, quadrature-order and timestep labels identify the
saved controls.

Each entropy product contains one completed control out of five requested.
The 191-cell bath refinement and the boundary control with half the maximum timestep each reached
the 600 s wall-time limit; the remaining three controls in each matrix were not
reached. Both full-matrix acceptance flags remain false. The figure displays
these incomplete-matrix counts and plots only completed saved distributions;
it supplies no uncomputed refinement curves.

The renderer imports no domain model. It verifies product schemas and pinned
producer/input identities, matches physical conditions, reconstructs native
moments and stable variance, checks number and mass ledgers, and independently
reproduces saved CDF differences. The reference's unbiased sample variance
(`ddof=1`) remains distinct from the candidate's population variance. Candidate
CDF values and the analytic Poisson reference come directly from the domain
products; connecting continuous-CDF query values is display interpolation. The
manifest retains plotted arrays, excluded tails, source and output identities,
upstream evidence statements, acceptance decisions and incomplete controls.

The shaded 95% DKW band belongs only to the empirical reference at the displayed
physical-bath endpoint. It describes sampling uncertainty; the analytic Poisson
control has no Monte Carlo band. These imposed-collector tests omit deposition,
settling and breakup. Endpoint compatibility supplies neither empirical ice
validation nor a lunar size distribution, and does not complete the full-column
collision repair. Generated `moment_comparison.png`, `moment_comparison.svg` and
`moment_manifest.json` remain ignored under `generated/`.

## Finite-time capture controls

`capture.py` reads the new Poisson-mixture and size-dependent birth-process
products, together with the frozen exact-event reference:

```sh
OPENBLAS_NUM_THREADS=1 python visualization/snow/capture.py
```

Four panels show analytic CDF errors around each grid's own first cell face,
the final physical-bath CDFs, their discrepancies from the independent event
histories at both nonzero observation times, and each physical control's recorded
runtime. Six variations in each matrix test timestep, quadrature and mass-grid
choices; shifted faces belong to the analytic matrix and ladder capacity to the
physical matrix. Analytic curves use a declared 0.005 CDF threshold. The separate
physical reference has a 95% DKW band at each time, uniform in mass but without
a joint guarantee across times or controls. It describes sampling uncertainty,
not a numerical convergence bound. The physical curves retain the original
collector normalization, with display exclusions recorded in the manifest.

The renderer verifies schemas and producer/input hashes, checks physical
settings and times, reconstructs native moments and stable variance, and checks
number and donor-water ledgers. It independently reproduces physical CDF
differences from the saved event histories and analytic differences from the
stored exact Poisson CDF. Independent gamma-function checks of the latter belong
to the domain tests. The renderer requires both complete, accepted matrices;
failed or missing controls cannot become unlabelled curves. It imports no domain
model and calculates no new snow physics.

The tests start with an imposed 15.86 cm collector and replenish every captured
donor. They establish neither collector formation nor a closed storm water
budget. Aggregate survival and predictive optics remain outside the calculation.
Recorded runtimes describe the initial two-worker execution and its single-worker
continuation; they are not universal performance bounds. The initial incomplete
matrix remains a separately bound product, and the final product identifies
every retained or repeated control. Generated `capture_comparison.png`,
`capture_comparison.svg` and `capture_manifest.json` remain ignored under
`generated/`. Their metadata and manifest bind renderer, inputs and outputs,
plotted arrays, display limits, sampling conventions and upstream evidence.


## Full prescribed snow column

`event_column.py` reads the completed [column control matrix](../../atmosphere/snow/results/event_column.json). Its four panels show lunar volume dependence, mass quantiles, the separate timestep/height controls and the generic Earth transfer control. Bands and bars show one standard error across whole independent realizations. The figure labels the failed numerical refinement and the absence of physical validation.

```sh
python -m visualization.snow.event_column
python -m pytest -q visualization/snow/tests/test_event_column.py
```

The generated `event_column_comparison.png`, `event_column_comparison.svg` and `event_column_manifest.json` live under ignored `generated/`. The reader verifies producer and input identities, recomputes ledger residuals, arrival-block totals, fluxes, CDF means and standard errors, checks quantiles against their native CDF brackets, and independently reconstructs each refinement decision. It rejects incomplete matrices and altered scientific quantities. Failed numerical refinements are displayed rather than admitted as forecasts. The [full-column report](../../atmosphere/snow/event_column_results.md) and [verification record](../../atmosphere/snow/results/event_column_checks.json) document the actual runs and artifact inspection.

## Volume-dependence diagnosis

`column_diagnosis.py` reads the primary 40-run diagnostic matrix and the separate eight-run sticking-only follow-up:

```sh
python -m visualization.snow.column_diagnosis
python -m pytest -q visualization/snow/tests/test_column_diagnosis.py
```

Its panels show lower-cloud empty-time fractions, the allocation of the prescribed 3 mm/h vapour supply, upper mass quantiles under process removals, and changes in the mean mass CDF when the statistical volume increases fourfold. The first two panels retain all three complete-model volumes; the process comparisons use the first two. Occupancy and water accounting describe growth-stage observations during gathering. Bands and quantile error bars show one standard error across four whole realizations. The CDF panel's five-percentage-point line tests that metric alone; overall acceptance additionally requires flux and all three reported quantiles to change by at most 5%.

The reader requires both complete matrices and verifies the follow-up's binding to the primary product. It checks producer and input hashes, reconstructs the water allocation and its deposition ledger, empty-time profiles, quantile means and standard errors, mean mass CDFs and comparison decisions. Its corruption tests reject altered scientific fields. Removing a process changes subsequent particle histories, so these comparisons do not partition independent physical contributions. They establish no calibrated giant-flake geometry or snowstorm optics.

Generated `column_diagnosis.png`, `column_diagnosis.svg` and `column_diagnosis_manifest.json` remain ignored under `generated/`. The manifest binds both input products, the renderer, verified upstream files and both image outputs. The [diagnosis report](../../atmosphere/snow/column_diagnosis_results.md) records the numerical findings and actual inspection; [methods](../../atmosphere/snow/column_diagnosis_methods.md) distinguish the observations from continuous residence times and from the separate offline averaging calculation.
