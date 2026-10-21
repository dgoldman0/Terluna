# Snow in the aerial winter: growth, motion and light

This analysis answers the 21 October question about turbulence, the kinds and sizes of flakes, and the light in a lunar snowstorm. It combines a revised population calculation with diagnostics of existing cloud simulations and the existing solved lunar sky. The atmospheric and optical diagnostics are conditional calculations. The inherited particle algorithm also has an identified numerical inconsistency, described below, so its size distributions remain exploratory. These results do not establish the frequency of such storms or validate the survival of decimetre aggregates.

The main distinction is between **snow in moist growth regions and snow falling through subsaturated air**. Moderate horizontal wind can carry a large aggregate along with the surrounding air. Its growth, vertical displacement and mechanical loading depend on moisture, updraughts, collisions and velocity differences across its structure. The existing simulations show much more vertical motion in the places where mixed-phase growth is possible than in cold air as a whole.

The subsequent [numerical repair and targeted physics review](repair.md) adds a conservative independent solver and distinguishes supported fractal geometry from unresolved contact-network strength. It preserves the distributions below as the audited earlier calculation; the first three repaired controls still fail their finite size-domain check and are reported separately. A further [growth-transport repair](growth_repair.md) isolates and substantially reduces artificial spreading during vapour deposition. Its numerical benchmarks, full-column controls and Earth comparison limits are separate from the historical distributions below. The latest [collision and moist-exposure analysis](collision_and_exposure.md) independently measures substantial collision-grid spreading and adds time-dependent paths through the archived storm; these findings further qualify the historical large-flake tail.

## The environments the flakes encounter

[turbulence.py](turbulence.py) reads 690 archived CM1 snapshots: all saved second-cycle outputs from two three-dimensional boxes, and a selected active storm window. It selects air between −20 and 0 °C at 20–60 km and then conditions the statistics on snow, ice supersaturation and cloud liquid. The fractions below give equal weight to sampled cell-times; they are not probabilities along a flake's trajectory.

| Environment | Morrison box: fraction with absolute vertical wind above 0.6 m/s | Corrected electrical box |
|---|---:|---:|
| All selected cold air | 0.62% | 0.34% |
| Snow mixing ratio at least 0.01 g/kg | 11.8% | 13.5% |
| That snow plus ice supersaturation | 37.8% | 44.0% |
| That snow plus cloud liquid at least 0.01 g/kg | 88.2% | 89.3% |

In the last row, median vertical wind is upward at 3.24 and 3.01 m/s; these liquid-bearing cells account for only 3.07% and 2.34% of snow-containing cell-times. Thus many growing flakes can be carried upward, even though a freely settling aggregate falls through still air at about 0.6 m/s. Horizontal wind in all snow-containing cells has medians of 3.53 and 3.90 m/s, with 99th percentiles of 6.13 and 6.65 m/s. These are results from the archived boxes' older forcing. They establish neither a Moon-wide wind climatology nor the weather of a newly simulated 301 K resort climate.

Moisture is equally consequential. In those two boxes, 75.8% and 64.5% of snow-containing cells have humidity below 95% of ice saturation. Snow can therefore occupy a deep falling curtain without growing throughout its depth. Among snowy columns, median maximum contiguous snow-bearing depth is 10 and 6 km, while the median ice-supersaturated depth is zero. Its 90th percentile is 4 km and its 99th percentile 16 km. The approximately 2 km cells and three-hourly full-cycle outputs cannot establish thin-layer structure or the lifetime of a moist path; the selected storm window has ten-minute output.

The ideal population column below instead supplies liquid-saturated air continuously over 16.3 km. This is a useful growth scenario, whose occurrence and persistence require additional climate evidence. Subsaturation could remove mass and delicate branches; ascent, recycling and turbulence-enhanced encounters could lengthen growth and increase aggregation. The archived statistics do not determine their net effect on the resulting size distribution.

## What turbulence can establish about survival

The resolved cloud motions are much larger in scale than a snowflake. The runs use diagnostic vertical mixing (`ipbl=2`) and have no prognostic LES turbulence-energy field to export. Enabling `output_tke` alone would not provide dissipation at flake scales. [The turbulence note](turbulence.md) records the code audit and the complete conditioning definitions.

A separate sensitivity imposes dissipation ε from 10⁻⁵ to 10⁻² m²/s³ and estimates the inertial-range velocity difference as `(ε D)^(1/3)`, with unit prefactor. The terrestrial mixed-phase reference range of 10⁻⁴–10⁻³ gives 0.025–0.053 m/s across 15 cm and 0.039–0.084 m/s across 60 cm. These differences are substantially smaller than the horizontal wind carrying the snow, while still imposing unequal aerodynamic loads on a delicate structure. No frequency is assigned to these scenarios in lunar clouds.

The load calculation reports an adverse aligned-slip pressure increment, `ρ V δu + ρ δu²/2`, assuming an unchanged drag coefficient. This includes the cross term around the settling slip, but does not resolve the loads on individual branches. Converting that loading into a survival threshold requires the aggregate's contact network, orientation, flexibility and bond strengths. The reviewed laboratory observations include aerodynamic shedding by artificial aggregates during release, so a quiet ambient wind cannot by itself establish intact survival. The current analysis calculates loading scales and identifies the necessary mechanical inputs; it supplies no invented breakup diameter.

## Crystal forms and population geometry

The large objects are aggregates of many crystals. Mapping terrestrial laboratory growth examples onto the temperature sequence in the 301 K reference column identifies conditions for plates and dendritic growth near −15 °C at 45.5 km, hollow columns around −7 °C at 36.8 km, needles and columns around −5 °C at 34.7 km, and plates near −2 °C at 31.4 km. An aggregate can carry crystals grown in several of those environments. Supersaturation, nucleation, pressure and riming influence the forms, so temperature alone cannot assign percentages of each habit. [The morphology review](habits.md) gives source identities and the relevant limits.

That review found two definition errors in the original geometry. Brown–Francis's coefficient 0.0185 applies to a mean image dimension; the published maximum-dimension conversion is `m = 0.0121 Dmax^1.9` in kilograms and metres. The quoted 2026 projected-area fits use a circumscribing ellipse; the settling law here needs area relative to a circle across the maximum dimension. Converted representative area ratios are 0.245 for needle aggregates, 0.313 for plates and 0.260 for dendrites.

The new population runs preserve the old geometry as a control, isolate the mass-definition correction, and use those three representative areas as sensitivities. They retain the inherited growth, sticking and small-fragment production equations. Holding an area ratio constant across every size and applying the aggregate mass law to both small seeds and very large flakes remain extrapolations. A case labelled “dendrites” describes that geometry choice, rather than a prediction that every crystal is dendritic.

## Population calculation and numerical controls

[population.py](population.py) exports the number, mass, projected area and settling flux in 96 logarithmic size bins through 10 m plus an explicitly open-ended tail bin, both inside 15 cloud layers and at the freezing level. At the boundary it converts each arriving particle's flux into concentration using that particle's own fall speed. The binning conserves the independent area moment, which the optical calculation needs. The specified 10 or 100 crystals per litre normalizes the injected number flux using a representative 0.2 mm crystal; it is not a separately imposed concentration for every injected size.

The cases vary water supply, seed concentration, geometry, gravity and random seed. A coupled numerical refinement doubles the collision layers, raises particle sampling and halves the maximum step. A mass ledger records injection, vapour deposition, arrivals, remaining ice and the tiny fragment mass discarded by the inherited particle-count rule. Three-hour startup windows record when snow first arrives and how its distribution develops; they do not simulate a cloud decaying or a flake passing through dry layers.

The [numerical audit](numerical_audit.md) found an additional limitation in the inherited collision bookkeeping. Injection uses fractional particle weights, but an exhausted collision pair is split by rounding half its weight down to an integer, and newly produced fragments are retained only when their weight reaches one. Rescaling the arbitrary represented horizontal area therefore changes the normalized particle population. The original super-particle method assumes integer multiplicities; applying those integer operations to fractional weights changes its premise. This is a demonstrated numerical defect, separate from the unknown strength of giant flakes.

A second concern arises when one sampled pair represents several collisions in a step: one sticking decision governs the whole batch. For a batch of `k` encounters with sticking efficiency `E`, its sticking-count variance is `k²E(1−E)`, compared with `kE(1−E)` for independent decisions. Random-seed and resolution controls quantify sensitivity within that implementation. They cannot establish convergence to the intended collision process. The audit preserves the baseline, tests arbitrary-area scaling and records collision batch sizes; it does not silently substitute a new algorithm.

Six shallow-column controls quantify the arbitrary-area effect: changing represented area from 1 to 256 m² at fixed physical inputs changes freezing-level mass and projected area by about −0.55% for the selected seed with fragmentation. Those millimetre-scale controls have at most one encounter per sampled pair and do not exercise the deep-column batching problem. Their paired differences include subsequent stochastic divergence; they are not an ensemble estimate of bias.

An instrumented replay reproduces the heavy, 100/L calculation exactly, including its mass ledger. The first crossing of 15 cm occurs at 6.07 model hours: a roughly 14.4 cm collector receives thirteen 1.63 cm donors under one sticking decision. The run reaches 10 m at 12.38 hours and 1 km at 38.13 hours, before roundoff-scale positive weight residuals first appear at 65.11 hours. Those residuals then amplify the late catastrophe: one merger involving weights of about 2.22 × 10⁻¹⁶ per square metre crosses a billion metres while carrying 0.403 kg/m² of represented ice.

The astronomical tail is physically impossible. Its timing shows that tiny rounding residuals do not explain the initial runaway; collision sampling and indefinitely extrapolated particle laws already permit nonphysical growth. The exports retain those particles, their mass and projected area so the failure remains visible. Clipping them would leave the earlier centimetre-to-decimetre batching untouched. Physical interpretation requires a consistent particle algorithm, demonstrated numerical convergence, and justified aggregate geometry and breakup treatment; this audit supplies no physical upper diameter.

The 20 October tables remain in [README.md](README.md) as the earlier calculation. Their “largest hourly” size is a cumulative arrival-rate threshold of one particle per square metre per hour, rather than a hard maximum or a usual flake size. The new product keeps that statistic beside mass percentiles and the full concentration distributions.

### Deep-column distribution results

The six central runs below use the corrected maximum-dimension mass law and representative dendrite area in the prescribed −15 to 0 °C saturated column. These are outputs of the inherited algorithm under the limitations just identified, not an established lunar snowfall distribution. Sizes are weighted by arriving mass. The imposed supply is the vapour-deposition budget; incoming seed ice adds water, while incomplete deposition and storage in the column can reduce the arriving flux.

| Imposed supply | Input flux normalization | Arriving water | Mass median | 90% of arriving mass below | 99% below |
|---|---:|---:|---:|---:|---:|
| 0.3 mm/h | 10/L | 0.30 mm/h | 1.75 cm | 3.6 cm | 5.8 cm |
| 0.3 mm/h | 100/L | 0.35 mm/h | 2.33 cm | 5.0 cm | 8.0 cm |
| 1 mm/h | 10/L | 0.95 mm/h | 8.37 cm | 15.8 cm | 22.5 cm |
| 1 mm/h | 100/L | 0.99 mm/h | 9.27 cm | 17.4 cm | 25.2 cm |
| 3 mm/h | 10/L | 2.70 mm/h | 19.70 cm | 36.4 cm | 157.1 cm |
| 3 mm/h | 100/L | 2.77 mm/h | 20.49 cm | 52.0 cm | Numerical runaway |

The final row puts 3.04% of arriving mass into model particles larger than 10 m, with its 99th mass percentile in an absurd numerical tail. That failure prevents reading the upper percentiles as physical flake sizes. In the same raw distribution, particles above 10 m account for 43.4% of the freezing-level optical extinction. Keeping an unchanged below-1-m component is useful for locating the optical contamination, but it neither repairs the smaller particles nor redistributes the excluded mass as a fracture model would.

The moderate-supply runs still place substantial calculated mass in decimetre aggregates. The remaining question is how much of that result persists under a consistent collision algorithm and a cloud history supported by the climate. Neither a single mass median nor an Earth comparison resolves those tests.

The geometry controls further expose the instability. The legacy-geometry moderate and heavy summaries exactly reproduce the stored 20 October results. With the corrected mass convention but the legacy area ratio, the heavy run puts 93.0% of arriving mass above 10 m; the representative plate-area case puts 43.2% there, compared with 3.04% for dendrite area and 0.099% for needle area. Those differences cannot be read as physical habit-dependent snow forecasts. They show how strongly the extrapolated geometry and numerical population interact.

Across three random seeds at 100/L, the moderate mass median ranges from 9.07 to 9.27 cm and its 99th mass percentile from 24.4 to 25.2 cm. The heavy mass median ranges from 20.04 to 20.49 cm, while its 99th percentile ranges from metres to the astronomical tail. Agreement among seeds in one part of a distribution does not establish accuracy of the underlying algorithm.

The coupled refinement uses 60 rather than 30 collision layers, twelve rather than eight injected packets per step, and a 15 rather than 30 second maximum step. For the heavy, 100/L, seed-0 case it changes the mass median from 20.49 to 19.46 cm (−5.0%) and the 90th percentile from 52.0 to 32.5 cm (−37.5%). Its 99th percentile is 1.05 m, and its raw maximum remains about 1.01 km. The arriving mass fraction above 10 m falls from 3.04% to 0.557%. This is substantial upper-tail sensitivity, not numerical convergence. The control changes three discretization settings together and uses one seed; it neither isolates their individual effects nor establishes an error bound.

The startup histories also distinguish initial-crystal crossing time from the transit of growing snow. At 100/L with either 1 or 3 mm/h supplied, the corrected deep column first delivers snow in the 9–12 hour window. Its initial 0.2 mm crystal would take 44.24 hours at the cloud-top speed, but growth increases its settling speed. The quoted post-fill distributions are sampled from about 66.4 to 132.7 model hours; that numerical sampling interval is not a measured cloud lifetime. The histories do not follow recycling in updraughts or passage through dry gaps.

### Shorter and warmer growth columns

Four additional runs shorten the continuously saturated path while keeping the same total supplied water, 100/L flux normalization and representative dendrite aggregate area. They inherit the numerical limitations above; their role is to test how strongly the assumed growth environment matters.

| Cloud top | Saturated depth | Imposed water supply | Arriving water | Mass median span | 99% of arriving mass below |
|---|---:|---:|---:|---:|---:|
| −5 °C | 5.43 km | 1 mm/h | 0.80 mm/h | 2.65 cm | 10.2 cm |
| −5 °C | 5.43 km | 3 mm/h | 1.66 mm/h | 8.59 cm | 30.5 cm |
| −10 °C | 10.85 km | 1 mm/h | 0.94 mm/h | 6.48 cm | 20.3 cm |
| −10 °C | 10.85 km | 3 mm/h | 2.40 mm/h | 16.37 cm | 49.7 cm |

For comparison, the −15 °C, 16.28 km case at 1 mm/h has a 9.27 cm mass median. Changing cloud top changes the temperature history as well as the available growth depth, and fixing total supply increases condensation per unit depth in the shorter cases. Holding the dendrite-area closure fixed isolates part of that comparison; it does not establish dendritic crystal production at a −5 °C cloud top. The [growth-depth product](results/growth_depth.json) records the complete distributions and settings. These imposed paths are still continuously saturated, whereas the archived weather contains dry gaps.

## Light through the snow

[optics.py](optics.py) evaluates snow-only extinction from the number-weighted projected area, using the large-particle approximation `β = 2 Σ n A`. Extinction derived from the particle populations inherits their numerical limitations; agreement with optical transport controls does not validate the supplied size distribution. Optical distances and beam transmission refer to a homogeneous snow-only medium with the time-mean distribution; exponentiating mean extinction does not give the mean transmission of fluctuating snowfall. For these porous geometries, area per unit mass scales only as `D^0.1`: increasing aggregate diameter fourfold increases area per mass by about 15%. Large, sparse-looking flakes can consequently still obscure distant views, particularly because slow settling raises the airborne ice concentration needed to carry a given mass flux.

The heavy, 100/L raw population gives a 5% direct-beam distance of 44.2 m at the freezing level; the numerical refinement changes it to 76.8 m. Considering only the unchanged below-1-m component gives 81.6 and 78.4 m, respectively. This comparison locates much of the opacity difference in the unstable tail. The component calculation is a diagnostic subtraction, not a predicted storm after fragmentation.

The illumination comes from the existing spherical molecular-sky solution at snow altitude. At 35 km, a Sun on the local horizontal still provides 8,119 lux normal to its beam and 5,855 lux of diffuse light onto an upper-facing horizontal surface. At −10° the corresponding values are 164 and 1,533 lux; at −20° the direct beam is blocked while diffuse light remains about 216 lux. Direct normal and horizontal illuminance are different quantities and are not added in this comparison. The optical column has a 294.9 K surface, while the snow population uses the 301 K interpolation; the mismatch is retained in the evidence statement.

A thin-path calculation shows how the reddened direct beam can make a forward view warmer than a side view illuminated mainly by diffuse skylight. An independent multiple-scattering slab calculation shows how diffuse illumination persists when the unscattered beam is strongly attenuated. These calculations do not yet resolve a finite storm jointly with its molecular atmosphere and liquid cloud. They predict neither the facets of individual aggregates nor halos and glints. The separately labelled full-Earth example is a brightness proxy with solar colour, rather than an Earthlight spectrum.

Supercooled droplets may dominate the appearance inside the growth cloud. An imposed liquid-water content of 0.1 g/m³ with 10 μm effective radius gives a droplet-only 5% beam distance of about 199 m; its extinction adds to the snow. A community beneath that cloud can see falling aggregates under a substantially clearer view than one embedded in its liquid droplets. [The optical note](optics.md) gives the assumptions, directional results and light table.

## Reproduction and evidence

Run the producers from the `domain/atmosphere` worktree:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.population
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.growth_depth
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.numerical_audit
python -m atmosphere.snow.assessment
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.turbulence
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.optics
python -m pytest atmosphere/tests/test_snow*.py -q
python visualization/snow/plot.py
```

The [assessment stage](assessment.py) attaches the audit warning directly to the population and growth-depth JSON products. It preserves their numerical results and original producer identities, records each original product hash and the numerical-payload hash, and identifies the separate assessor and audit. Re-running either raw producer requires repeating this stage before downstream use.

The population caches case outputs in ignored `research/runs/snow_population`, binding the case settings, model files, input products and named constants. The turbulence product records the actual archived field selections and their hashes. The optical producer rejects missing archives and mismatched source hashes. Numerical conservation, sampling sensitivity and transport controls establish implementation evidence; empirical validation and a mechanically justified large-aggregate survival limit remain outstanding.

The [four-panel overview](../../visualization/snow/generated/overview.png) compares the raw distributions, optical-area weighting, conditioned winds and incident light. Its [renderer and provenance rules](../../visualization/snow/README.md) preserve the off-scale tail in the normalization and identify the separate input calculations. The PNG was visually inspected for labels, scaling, clipping and qualification placement; the SVG was checked as XML with embedded provenance, without a separate SVG rendering review.

The final focused suite passed **65 tests**. Layer checks covered tracked and new source files with no violations, the research integrity and numerical checks passed, and the legacy moderate and heavy summaries reproduced the stored baseline exactly. [The check record](results/checks.json) gives commands, hashes and scope. These checks establish implementation and provenance evidence; they do not resolve the demonstrated failure of the heavy size distribution. The full `make check` suite was not run, and no commit was made.
