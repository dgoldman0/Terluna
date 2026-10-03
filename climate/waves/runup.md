# Individual waves on the Smythii shore

SWASH resolves the water surface, breaking and shoreline motion on explicit
cross-shore profiles. The forcing comes from the directional coastal spectra
in the [shore exposure](shore.md) and [evolving coast](coastal_history.md)
studies. The [result product](results/runup.json) keeps every completed case,
its forcing, numerical choices, finite sampling window and checks.

## Results for the evolving coast's strongest incoming states

The western and southern profiles use the hour-855 spectra, about 6.09 Earth
days into the second cycle. The eastern profile uses hour 1,061, about
14.68 days into that cycle. These saved hours maximize incoming power at
their respective stations. Each selected spectrum is held steady during its
wave-resolving calculation, with a 30-minute reported waterline history.

| Profile | Incoming normal Hs | Incoming mean period | 98th time percentile of run-up | Largest run-up in the window |
|---|---:|---:|---:|---:|
| Western rock, 0.25 m grid | 1.68 m | 12.37 s | 2.07 m | 2.40 m |
| Southern rock, 0.5 m grid | 1.04 m | 13.29 s | 1.31 m | 1.68 m |
| Eastern rock, six layers, 0.5 m grid | 0.85 m | 8.62 s | 0.67 m | 0.76 m |
| Western 1:20 beach, 1 m grid | 1.68 m | 12.37 s | 1.05 m | 1.10 m |
| Western 1:50 beach, 0.5 m grid | 1.68 m | 12.37 s | 0.98 m | 1.08 m |
| Western 1:100 beach, 1 m grid | 1.68 m | 12.37 s | 0.87 m | 0.90 m |

Run-up elevations are vertical metres above still water. The incoming heights
in this table describe the directionally projected spectrum. The coastal
history reports the full directional significant height at the station.
The eastern result also carries the incident-boundary sensitivity described
below. The gentler beaches retain their measured grid and initialization
uncertainty.

The western rock produces a concentrated breaking reach. Its seaward edge
lies about 10.6 m from the still-water shoreline in 2.82 m water, with a local
resolved Hs of 1.87 m. The active reach spans 17 m, extending onto rock above
still water as waves run up. The same incoming sea on the 1:20 beach begins
its sampled breaking reach about 60.6 m offshore, in 3.03 m water. That reach
spans 71 m. The assumed underwater slope therefore changes both the spread
of breaking and the elevation reached at shore.

On the refined 1:50 beach, the breaking reach starts about 149 m offshore in
2.98 m water, with a local Hs of 1.50 m, and spans 152 m. Its 1-to-0.5 m
grid comparison changes the upper time percentile from 1.029 to 0.979 m,
a 5.19% difference relative to the finer result. The maximum changes from
1.103 to 1.080 m. The upper-percentile comparison exceeds the 5% bound;
the beach's reported height and breaking geometry retain that sensitivity.

Halving the western rock spacing changes its upper time percentile by 1.1%
and its maximum by 3.6%. A second phase realization at the original spacing
gives a 2.66 m maximum, compared with 2.31 m for the original phases, while
the upper time percentile changes from 2.09 to 2.15 m. Adding a dimensionless
bottom-friction coefficient of 0.005 gives 2.05 m for that percentile;
reducing the breaking-onset coefficient from 0.6 to 0.5 gives 2.05 m.
These cases show greater stability in the upper waterline percentile than in
the largest individual excursion. More realizations are needed for a
probability distribution of extremes.

The gentlest beach develops slow shoreline excursions within the finite
window: its upper time percentile is 0.73 m in the first half and 0.89 m in
the second, with 15 detected peaks across the full 30 minutes. Longer sampling
and a finer 1:100 grid remain necessary to assess that response.

The [profile figure](../../visualization/waves/results/runup.png) shows actual
computed water surfaces near the largest excursions and each complete
reported waterline history. The cycle-long arrival and duration statistics
remain in the [coastal history](coastal_history.md); a continuous shoreline
history through those changing spectra is a further calculation.

## Geography and the question each profile answers

[shore_profiles.py](../../geography/shore_profiles.py) samples the pinned
LOLA/GRAIL flooded-rock terrain along the direction of decreasing depth at
each coastal station. The source raster has 118.45 m spacing. Bilinear samples
every 2 m provide a continuous numerical bed; the geographic information
retains the source spacing. The wave solver then uses a finer computational
grid to resolve motion across that assumed bed.

| Station | Starting depth | Distance to still-water shore | Mean submerged slope |
|---|---:|---:|---:|
| Western face | 22.38 m | 85.1 m | about 1:3.8 |
| Eastern face | 49.04 m | 318.7 m | about 1:6.5 |
| Southern shore | 13.55 m | 65.6 m | about 1:4.8 |

These are steep rock profiles at the atlas's fixed water level,
−1,654.2 m above the geoid. Smooth 1:20, 1:50 and 1:100 slopes provide explicit
beach scenarios. Their eventual occurrence depends on sediment supply,
transport and shoreline evolution. All profiles include a 100 m reach of
constant depth between the wavemaker and the start of the slope, and land
extending at least 15 m above still water.

## Coupling the models

The one-dimensional boundary spectrum is

`E₁ᴅ(f) = ∫ E(f, θ) max(n · uθ, 0) dθ`,

where `n` points toward shallower water. This retains the shoreward projection
of the incoming energy transport at the source depth. It reduces the incoming
directional sea to normal incidence. Alongshore variation, oblique-wave
transformation and currents require a two-dimensional wave-resolving model.

The code conservatively integrates the piecewise-linear SWAN spectrum onto
a uniform 1/3,600 Hz frequency grid. Each spectral bin receives its integrated
variance. The transfer records its incoming significant height, removed
variance and group-velocity-weighted power error. The main cases reject a
vertical-layer cutoff that would remove over 1% of incoming variance. Two
explicit eastern sensitivity cases use an identical narrower band at four
and six layers, removing 4.52% of incoming variance under a recorded 5%
guard. These cases isolate the layer and wavemaker choices at a common input.
The spectral
record and upstream case manifest are hash checked; gravity and water density
must agree across the coupling.

The transfer retains the frequency band resolved by SWAN. Each selected
history case records the incoming variance near both edges of that source
band. The upstream short-wave coverage limit therefore accompanies its
shoreline response.

The boundary uses first-order wave components with a recorded random seed.
The phases repeat after 3,600 s. This allows the initialization check to
compare the same forcing phases one full period later. The model develops
nonlinear waves and low-frequency motion within the domain; bound harmonics
and incoming infragravity motion at the offshore boundary remain additional
forcing choices.

## Solver and lunar settings

The unmodified **SWASH 12.01** source comes from the
[official distribution](https://swash.sourceforge.io/download/download.htm).
The source archive has SHA-256
`c8a2e149e4c1583f3f5a10fa6e9fa110e6db7d1eb41978420a4a387877ac397f`.
It is licensed under GNU GPL version 3 or later and credits Delft University
of Technology. [swash_build.py](swash_build.py) retains source, compiler and
executable hashes on the research drive.

The source audit follows configured gravity through dispersion, the spectral
boundary and the breaking detector. Gravity-scaling controls exercise the
pressure solution and shoreline motion. The inspected routines are
`SwashReadInput`, `SwashBCshortwave`, `SwashBCspecfile`, `SwashBreakPoint`,
`SwashRunupHeight`, `SwashCheckPrep`, `SwashBounCond` and `SwashBCStokeswave`;
their source hashes are in the build record. The
[SWASH manual](https://swash.sourceforge.io/online_doc/swashuse/swashuse.html)
describes the boundary, pressure and output choices used here.

| Setting | Choice |
|---|---|
| Gravity | `GM/R²` from shared constants, 1.62421887656 m/s² |
| Water density | 1,025 kg/m³, the existing seawater scenario |
| Vertical resolution | Four layers initially; six-layer control |
| Pressure | Non-hydrostatic Keller-box scheme; ILU preconditioner |
| Advection | Momentum-conserving upwind treatment |
| Incoming boundary | Weakly reflective spectral wavemaker; 120 s ramp |
| Bottom friction | Zero in the reference cases; explicit sensitivity cases retain their coefficient |
| Breaking | Onset coefficient 0.6; persistence coefficient 0.3 |
| Wetting and drying | Minimum wet depth 0.005 m; run-up diagnostic depth 0.01 m |
| Outputs | Waterline every 0.1 s; native-node surface and breaking flag every 1 s |

The density is recorded in the shared [wave scenario](../../shared/scenarios/waves.json).

For up to four uniform Keller-box layers, SWASH constructs incident velocities
from its discrete first-order Stokes solution. Above four layers it uses the
continuum hyperbolic-cosine velocity profile. The manual describes the
resulting boundary-height adjustment at large depth-to-wavelength ratios.
A four-to-six-layer comparison consequently changes both vertical resolution
and the incident velocity formulation. The full eastern spectrum needs six
layers to meet the 1% retained-variance bound; its absolute shoreline response
therefore receives particular scrutiny in the common-band cases.

With identical incoming spectrum bytes, the eastern four-layer case gives
a 0.718 m upper time percentile and a 0.841 m maximum. Six layers give
0.672 and 0.760 m. The 6.75% percentile difference and 10.66% maximum
difference exceed the respective 5% and 10% bounds. At six layers alone,
restoring the wider input band changes the upper percentile by only 0.0018 m
and the maximum by 0.0004 m in these realizations. The combined vertical and
wavemaker sensitivity is therefore material to this eastern result.

The default ILUD pressure preconditioner failed in the initial four-layer
spectral diagnostic. ILU completed both the flat-depth diagnostic and the
steep profile. The retained failed runs document that solver choice.

## What the reported quantities mean

Run-up is the vertical elevation of the moving waterline above still water.
The **98th time percentile** is exceeded during 2% of the sampled time.
The separate **98th peak percentile** describes detected local maxima,
separated by at least 0.35 mean wave periods and with 1 cm prominence.
The maximum is the largest value in that finite realization. Each case
records its reporting duration and the number of detected peaks.

A native model node belongs to the reported breaking reach when its breaking
flag is active in at least 1% of the one-second output samples. The reported
depth range uses the still-water bed; negative depths lie above still water
and can become wet during run-up. The reach spans the active nodes and may
contain gaps. The flag indicates the model's bulk
breaking treatment. Overturning crests, spray and entrained air require
additional fluid modelling.

## Independent controls

At lunar, fourfold lunar and Earth gravity, the linear-wave tests compare
amplitude and phase against the dispersion relation. The largest amplitude
error is 1.76%, and the largest phase error is 0.289%. Scaling gravity by four
and the clock by one half preserves the linear waveform to 0.000908% of the
imposed amplitude.

Still water stays level at both 0.5 and 0.25 m grid spacing. The waterline
diagnostic has a constant rest offset of 4.07 and 0.892 cm respectively from
its wet-depth and interpolation convention. These offsets accompany the
reported run-up values. The nonlinear gravity control, including breaking
and shoreline motion, changes the 98th time percentile by 1.20% and the
largest profile-height comparison by 0.189%. Individual waterline samples
differ by up to 0.277 m in that nonlinear comparison.

A closed-wall control separates incident and reflected waves at three gauges.
Its reflected-to-incident amplitude ratio is 1.0102, with a 0.776% fit residual.
Together these are executable and numerical checks. Physical calibration of
lunar wave growth, breaking, friction and shoreline turbulence remains open.

A separate control drives a flat 22.38 m deep flume with the transferred
stationary western spectrum. Three gauges at 25, 50 and 75 m cover a complete
3,600 s phase period after 900 s of initialization. The specified Hs is
1.719 m and mean period 12.736 s. At 0.5 m spacing and four layers, gauge
heights are 1.659, 1.639 and 1.629 m; their mean periods are 13.200, 13.356
and 13.537 s. The maximum differences, 5.22% in height and 6.28% in period,
exceed the 5% input-control bounds. Six layers at the same spacing give
5.62% and 6.74%. The controls have zero breaking flags. Their differences
include boundary reconstruction, numerical propagation, nonlinear evolution
and residual reflection; the shoreline estimates retain this input uncertainty.
Halving the four-layer spacing to 0.25 m gives heights of 1.666, 1.645 and
1.634 m and mean periods of 13.095, 13.269 and 13.466 s. Its maximum height
difference falls to 4.93%, while its 5.73% period difference still exceeds
the bound. The small improvement leaves the irregular-input control open.

For the earlier stationary hour-852 spectrum, changing the rock grid from
0.5 to 0.25 m moves the 98th time percentile from 2.179 to 2.136 m and the
maximum from 2.471 to 2.541 m. Increasing from four to six vertical layers
at 0.5 m gives a 98th time percentile of 2.175 m. These cases establish the
numerical reference for the later changing-coast inputs.

Extending the constant-depth wavemaker reach from 100 to 200 m changes the
rock case's upper time percentile from 2.179 to 2.104 m and its maximum from
2.471 to 2.558 m. The change includes the different travelling phases in
the finite sampling window. This records the sensitivity to the artificial
offshore boundary's position alongside the spatial and vertical controls.

The beach initialization check uses the identical boundary phases one or two
3,600 s periods later. On the 1:50 beach, extending initialization from
600 to 4,200 s changes the mean waterline from 0.409 to 0.482 m. Extending
it again to 7,800 s gives 0.482 m, while the upper time percentile changes
from 0.985 to 0.994 m. The main beach cases therefore use 4,200 s of
initialization followed by a 1,800 s reporting window. On the 1:100 slope,
the shorter initialization produces a clear transient; its results remain
in the numerical record alongside the longer calculation.

For the settled 1:50 reference beach, halving the horizontal spacing from
1 to 0.5 m changes the upper time percentile from 0.985 to 1.048 m and the
maximum from 1.056 to 1.151 m. The sampled breaking reach extends from 117
to 161 m. This comparison exceeds the 5% upper-percentile bound; it motivates
the separate 0.5 m beach calculation for the evolving-history peak. The
1:20 and 1:100 beach cases retain their 1 m numerical spacing.

## Reproduction and storage

Run commands from the repository root. The established wave-run link points
to `/media/projectspace/terluna-research/wave-runs`; source archives, builds,
spectra, surface histories and model logs stay there. The repository holds
the runners and compact evidence products.

```bash
python -m geography.shore_profiles
python -m climate.waves.swash_build
python -m climate.waves.runup controls
python -m climate.waves.runup case \
  --source research/runs/waves/shore/h852_s1_d36/sites.spc \
  --name calibration_west_rock_dx0p5_k4_ilu
python -m climate.waves.runup_controls
python -m climate.waves.runup case --source <coastal-sites.spc> --hour <saved-hour> \
  --station west_face --kind rock --dx 0.5 --step 0.025 --name <case-name>
python -m climate.waves.runup case --source <coastal-sites.spc> --hour <saved-hour> \
  --station west_face --kind slope --slope 50 --dx 1 --step 0.05 \
  --seconds 6000 --spinup 4200 --name <beach-case-name>
python -m climate.waves.runup_checks --select <case-names>
python -m visualization.waves.runup
```

Case identities include executable and input hashes. Repeated commands reuse
completed matching runs, and changed inputs require a new case directory.
Each SWASH run uses one CPU and a 2 GiB address-space cap. The figure displays
the computed bed, water surface and complete reported waterline history,
with exact array checks and a provenance sidecar.

The original hour-855 1:50 worker ended during initialization. Its partial
directory and interruption record remain on the research drive; the complete
1 m replacement is `history_west_slope50_dx1`. The
[verification record](shore_history_verification.json) links these retained
diagnostics, all completed controls and the repository checks.
