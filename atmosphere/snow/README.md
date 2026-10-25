# Snow

The snow work asks four things of the Open Moon's sky winter: how large snowflakes grow, how they fall, whether the
largest survive the cloud's turbulence, and how a snowstorm looks under lunar light. The author opened these questions
on 20 and 21 October, as the second half of the work on the high cloud's crystals
([register](../../research/decisions.md), 20 and 21 October). How flakes fall, how deep the snow zone stands and how
far below the freezing level snow lasts are calculated in the sections below. The size of lunar snowflakes is open.
The calculations made between 20 and 24 October put typical flakes anywhere from below a millimetre to tens of
centimetres across, and the spread follows from three inputs: the water a storm supplies to growing ice, the number of
crystals that share it, and the population solver that follows their collisions.

On 24 October the author made three decisions about the course of the work. The crystal-growth model is frozen as it
stood that day, and its errors against Earth's laboratory crystals enter every lunar result as a stated uncertainty. A rule sized to the laboratory data replaces the 1,320-point accuracy rule. The work then moves to the
lunar answer, in the order of the plan below.

## The plan, in order

1. **Check the frozen model by the data-sized rule** (below), naming every measured crystal that falls outside it.
2. **Carry the Earth errors through the lunar column.** Run the frozen model at the edges of its Earth error band in
   mass growth, axis growth and fall speed, through the archived lunar storm's humidity and the representative
   moisture cases, and set the change in lunar flake size beside the spread from moisture supply and crystal number.
   This step tests whether the crystal model's Earth errors matter for the Moon.
3. **Give lunar flake sizes as a range** over moisture supply and crystal number.
4. **Join the flakes, checked first on Earth.** With measured crystals as input and established aggregation
   efficiencies, the aggregation calculation is compared with the sizes of observed Earth snowflakes, then run through
   the lunar column.
5. **Estimate survival in turbulence**: whether the largest aggregates hold together in the motions of the cloud that
   grows them ([register](../../research/decisions.md), 21 October).

Steps 4 and 5 start from the crystals of step 3 and carry the uncertainty measured in step 2.

## The frozen crystal model

The model has two parts, both frozen at their 24 October state.

The **improved McSnow component** calls McSnow's own deposition, habit, geometry and fall-speed routines (Welss and
co-authors 2024; pinned revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`) with one change: the extra oblate
mean-ventilation term is removed ([growth trace](benchmark_growth_trace.md)). It runs at any temperature, pressure and
gravity, and every lunar column and storm calculation since 23 October uses it ([lunar column](mcsnow_lunar_growth.md),
[lunar storm](mcsnow_weather.md)).

The **two-temperature calibration** fits growth controls to the laboratory histories at −5.3 °C, where hollow columns
grow, and separately at −14.4 °C, where dendrites grow ([column density](mcsnow_column_density.md),
[growth limits](mcsnow_growth_limits.md); commit `ded734e`). It holds no rule for other temperatures. It shows how
closely the component can be brought to the laboratory, and it cannot run through a lunar column, which spans every
temperature from the cloud top to the freezing level.

Before the freeze, the component's −14.4 °C dendrite carried 26% too much mass at 30 minutes against Takahashi and
co-authors' (1991) fitted laboratory history, with a span 35% too wide and a thickness 17.5% too thin; its −5.3 °C
column grew too short, 31% below the fitted length at the fitted mass ([shape diagnosis](mcsnow_shape_diagnosis.md)).
The calibration brings the warm column within 9.73% of the fitted mass and dimensions from 3 to 30 minutes, and the
cold dendrite within 4.75% at 30 minutes with its 3-minute mass about 17% low; the warm column's fall speed stays 21–23%
high at 30 minutes. These agreements are fits to the same histories ([joint audit](mcsnow_joint_acceptance.md)). Both
parts share three recorded limits. Crystals started from elongated seeds run away. The shape and density rules ignore humidity
([recent-science check](mcsnow_science_check.md)). And the laboratory crystals grew for at most 30 minutes, against the
18–54 hours a lunar crystal spends in the cloud.

## The data-sized rule

The rule compares the model with each crystal Takahashi and co-authors measured: 169 crystals grown in a vertical wind
tunnel at water saturation, whose masses and axes Fuchs and co-authors archived in 2025
([measurements](results/takahashi_measurements.json)). Each crystal is compared at its own age and temperature. It
passes on mass and on each measured axis when the model lies within the scatter between repeat crystals grown at the
same temperature and age. The archive holds no fall speeds crystal by crystal, so the model's fall speed is compared
with the paper's fall-speed relation (its Table 4, each line the average of about four sets of measurements) within the
paper's stated measurement accuracy of ±3 cm/s. The rule replaces the [1,320-point rule](mcsnow_joint_acceptance.md),
which remains as the record of that audit and the candidates it tested.

## Files and commands

The snow solver sits beside the raindrop solver ([../rain](../rain/README.md)). The table lists the principal files;
the history below links the rest.

| File | What it holds |
|---|---|
| [scale_size_results.md](scale_size_results.md), [scale_size_diagnosis.py](scale_size_diagnosis.py), [scale_size_sources.json](scale_size_sources.json) | Gravity/depth attribution, water-per-crystal and aggregation controls, active-source audit and limits of physical interpretation |
| [scale_lunar_results.md](scale_lunar_results.md), [scale_lunar_methods.md](scale_lunar_methods.md), [scale_lunar_comparison.py](scale_lunar_comparison.py) | Complete SCALE/custom lunar comparison, numerical refinements, sticking and atmospheric-feedback contrasts |
| [flakes.py](flakes.py) | The snowflake's mass, area and fall law, its growth from the vapour, its melting, the winter columns and the product |
| [results/snow.json](results/snow.json) | Fall speeds by size, growth across the winter cloud, crossing times and melting depths, Moon and Earth |
| [aggregation.py](aggregation.py) | A column of falling snow: super-particles that grow, stick together and break, and the snow they bring to the freezing level |
| [results/aggregation.json](results/aggregation.json) | The earlier arriving-snow sizes by supply and crystal numbers, retained as a baseline |
| [population.py](population.py), [results/population.json](results/population.json) | Geometry controls and number, mass, area and flux distributions, including cloud startup |
| [growth_depth.py](growth_depth.py), [results/growth_depth.json](results/growth_depth.json) | Shorter and warmer continuously saturated growth columns |
| [numerical_audit.py](numerical_audit.py), [numerical_audit.md](numerical_audit.md) | Fractional-weight bookkeeping, arbitrary-area controls and collision-batch diagnostics |
| [sectional.py](sectional.py), [repair.py](repair.py), [results/repair.json](results/repair.json) | Conservative independent solver and completed domain/refinement diagnostics, with explicit adequacy flags |
| [column_diagnosis_methods.md](column_diagnosis_methods.md), [column_diagnosis_results.md](column_diagnosis_results.md) | Passive water and occupancy accounting; 48 complete process-removal and refinement controls |
| [event_column.py](event_column.py), [event_column.cpp](event_column.cpp), [event_column_benchmarks.py](event_column_benchmarks.py) | Continuous-mass collision events, donor depletion, complete prescribed columns and independent numerical controls |
| [event_column_results.md](event_column_results.md), [results/event_column.json](results/event_column.json) | All 28 full-column results, whole-realization uncertainties and failed lunar volume refinement |
| [kernel_diagnostics.py](kernel_diagnostics.py), [results/kernel_diagnostics.json](results/kernel_diagnostics.json) | Fixed-geometry encounter scaling and the installed chipping rule's mass limits |
| [geometry_literature.md](geometry_literature.md), [strength_literature.md](strength_literature.md) | Targeted primary evidence for evolving geometry, interception, ice contacts and aggregate failure |
| [turbulence.py](turbulence.py), [turbulence.md](turbulence.md) | Archived cloud motions conditioned on snow and growth environments; separate flake-scale loading scenarios |
| [habits.py](habits.py), [habits.md](habits.md) | Crystal-form context, source-definition corrections and representative aggregate geometries |
| [assessment.py](assessment.py) | Audit-linked evidence metadata on the unchanged numerical products, with separate provenance |
| [optics.py](optics.py), [optics.md](optics.md) | Snow extinction, height-resolved lunar illumination, thin-path radiance and multiple-scattering controls |

```sh
python -m atmosphere.snow.flakes        # seconds: results/snow.json
python -m atmosphere.snow.aggregation   # about 20 minutes: results/aggregation.json
```

## The physics

The baseline treats a snowflake as an aggregate of ice crystals, with m = 0.0185 D^1.9 kg and projected area 0.35 of the circle across D. The [source review](habits.md) found that the mass coefficient uses a mean image dimension rather than Dmax, and that the cited area fits use an ellipse rather than that circle. The 2026 source is Seifert and co-authors, using Leinonen's aggregation method. These baseline closures are preserved for comparison; the new population analysis uses the corrected dimension convention and explicit area sensitivities. It falls by Heymsfield and
Westbrook's (2010) law, in which gravity enters through the modified Best number alone,
X* = 8 m g ρ / (π Ar^½ η²), with Re = (δ0²/4)((1 + 4 X*^½ / (δ0² C0^½))^½ − 1)², δ0 = 8.0, C0 = 0.35 and the speed
η Re / (ρ D).

It grows from the vapour by diffusion to a capacitance of 0.25 D (Westbrook, Hogan and Illingworth 2008), ventilated as
Hall and Pruppacher (1976) give it for ice, against Murphy and Koop's (2005) saturation over ice and supercooled water.
Below the freezing level it melts by the heat it takes from the air less what its wet surface's evaporation carries
away (Pruppacher and Klett 2010, chapter 16), at the same capacitance and ventilation. The air's viscosity,
conductivity and vapour diffusivity are the raindrop solver's (`drops.transport`) for each world's air.

A winter column is a straight temperature profile through the freezing level, its pressure hydrostatic:

| Column | Freezing level | Pressure there | Air there | Lapse |
|---|---|---|---|---|
| Moon: the resort's sky ladder at 301 K | 29.2 km above sea level | 659 hPa | 0.83 kg/m³ | 0.92 K a km (fitted over 20–40 km) |
| Earth: a winter column | about 1.5 km | 850 hPa | 1.08 kg/m³ | 6.5 K a km |

The Moon's ladder is the climate targets' (research/studies/climate_targets), read as the rain products read it
([targets.py](../rain/targets.py)). The cloud runs from a top temperature down to the freezing level at liquid
saturation, a mixed-phase cloud where snow grows fastest.

## What it finds

**Snow drifts down.** Each world's flakes in the air at its freezing level, m/s:

| Flake | 0.2 mm | 0.5 mm | 1 mm | 2 mm | 5 mm | 1 cm | 2 cm | 5 cm |
|---|---|---|---|---|---|---|---|---|
| Moon | 0.12 | 0.22 | 0.31 | 0.40 | 0.51 | 0.57 | 0.61 | 0.64 |
| Earth | 0.54 | 0.82 | 1.02 | 1.19 | 1.35 | 1.43 | 1.47 | 1.50 |

Lunar flakes fall at a quarter to two-fifths of Earth's speed; the largest drift down at about 0.6 m/s, 2.3 km an hour.

**The winter cloud is deep, and snow spends hours in it.** The air cools only 0.92 K a km near the Moon's freezing
level, so the same span of cloud temperatures stands seven times deeper than on Earth and holds about five times the
air. From each cloud top down to the freezing level:

| Cloud top | Moon: depth; air | Moon: a 5 mm flake's crossing | Earth: depth; air | Earth: crossing |
|---|---|---|---|---|
| −10 °C | 10.9 km; 8,275 kg/m² | 5.8 hours | 1.5 km; 1,543 kg/m² | 0.31 hours |
| −15 °C | 16.3 km; 11,851 kg/m² | 8.5 hours | 2.3 km; 2,226 kg/m² | 0.45 hours |
| −20 °C | 21.7 km; 15,089 kg/m² | 11.1 hours | 3.1 km; 2,856 kg/m² | 0.59 hours |

**So snow grows far more.** A 0.2 mm crystal crossing the cloud gains from the vapour alone 2.2, 12.5 and 43 mg on the
Moon from the three tops, against 0.010, 0.025 and 0.049 mg on Earth: 200–900 times the mass, as a growing crystal
grows faster the larger it is. With five times the air in the cloud, a flake also has about five times the cloud's
small crystals to sweep up at one mixing ratio. A single flake sweeping a cloud that never runs out grows without
bound; how large flakes grow in a whole population, with a breakup law, follows below.

**Snow lasts far below the freezing level.** Below it the air warms only 0.92 K for each km a flake falls, so flakes
melt slowly. Where all of a flake's ice has melted, metres below the freezing level:

| Flake | Moon, air at 90% | Moon, at 70% | Earth, at 90% | Earth, at 70% |
|---|---|---|---|---|
| 2 mm | 1,045 | 2,849 | 228 | 450 |
| 5 mm | 1,137 | 2,942 | 266 | 489 |
| 1 cm | 1,212 | 3,019 | 297 | 522 |
| 2 cm | 1,294 | 3,103 | 332 | 557 |

In drier air a flake's own evaporation keeps its surface below freezing for longer. The sky winter's snowfall is a
band kilometres deep: snow forming in a cloud above the freezing level reaches decks 1–3 km below it.

## How large flakes grow (20 October, superseded)

**Superseded.** The sizes in this section, a centimetre or two in light snow to tens of centimetres in heavy snow, come
from one solver, the original super-particle column, and later work superseded them. The numerical audit of 21 October
found that this solver's fractional-weight bookkeeping changes its results with the arbitrary area it represents
([numerical audit](numerical_audit.md)). An independent event solver under the original lunar forcing never converged:
its median grew from 12.3 to 18.5 to 23.5 cm with each fourfold increase in sampled volume
([event column](event_column_results.md)). SCALE-SDM under the same forcing gave medians near 1.6 mm
([SCALE lunar column](scale_lunar_results.md)). The improved crystal component grew 3.14 mm crystals under the same
3 mm/h allowance, and reached a 0.72 mm median peak through an archived lunar storm's humidity
([lunar column](mcsnow_lunar_growth.md), [lunar storm](mcsnow_weather.md)). Step 3 of the plan replaces these sizes with
a range over moisture supply and crystal number.

On 20 October the author asked for a breakup law, to say how large the flakes grow. [aggregation.py](aggregation.py)
follows a column of falling snow. Crystals of 0.1–0.4 mm enter at a −15 °C cloud top, 10 or 100 a litre, and fall at
their own speeds to the freezing level, as super-particles each standing for many real ones (Shima and co-authors
2009). They grow from the vapour as far as the cloud condenses it for a set snowfall, spread through the cloud's
depth. They meet the others at their height, sweeping the area of their projected areas (Mitchell 1988) at their
difference in speed, and a meeting sticks with Lin, Farley and Orville's efficiency. One that does not stick breaks
fragments off the smaller particle by Phillips, Yano and Khain's (2017) law for snow, which Sotiropoulou and
co-authors set out (Atmospheric Chemistry and Physics Discussions, acp-2020-786, appendix B): up to 100 fragments of
0.1 mm a collision, more the harder it is. The snow arriving at the freezing level, its sizes by mass and the largest
flake landing at least once an hour on each square metre:

| Snowfall | Crystals | Moon: median; 99% below; largest hourly | Earth: median; 99% below; largest hourly |
|---|---|---|---|
| 0.3 mm/h | 10 a litre | 10; 34; 50 mm | 0.7; 1.2; 1.5 mm |
| 0.3 mm/h | 100 a litre | 18; 63; 86 mm | 0.5; 1.4; 2.0 mm |
| 1 mm/h | 10 a litre | 51; 152; 197 mm | 0.9; 1.8; 2.3 mm |
| 1 mm/h | 100 a litre | 61; 179; 225 mm | 0.9; 2.5; 5.4 mm |
| 3 mm/h | 10 a litre | 141; 387; 518 mm | 0.9; 1.8; 2.2 mm |
| 3 mm/h | 100 a litre | 147; 468; 600 mm | 1.7; 9.2; 20 mm |

The snowfall is in millimetres of water an hour. On Earth 10 crystals a litre take up only 0.54 mm/h of the cloud's
vapour in the two heavier cases, the rest left unused; the arriving snow also counts the crystals' own mass entering
the top, about 0.4 mm/h at 100 a litre on Earth.

- **The breakup law does not cap the flakes.** With it and without it the sizes agree within 10–25%, and not always
  smaller with it. Phillips's law chips small fragments off the smaller particle when a meeting does not stick, which
  feeds the cloud more small crystals but splits no large flake. On both worlds the flakes' size is set by the
  snowfall's water and the time the snow takes to fall.
- **Lunar flakes grow 15–90 times the size of Earth's at the same settings:** a centimetre or two in light
  snow, 5–15 cm in moderate snow and tens of centimetres in heavy snow. The smallest crystals take 1–1.5 days to fall
  through the Moon's 16 km of cloud against about an hour on Earth, in a model whose sticking probability depends on temperature and whose encounters use differential settling alone.
- **How far to trust the sizes.** The model's Earth flakes run small against Earth's real snow, where flakes of
  several millimetres to a couple of centimetres are common, so it aggregates too little, most likely for the spread of
  fall speeds among like flakes, which it leaves out; that alone does not determine the sign of the lunar size bias. On the other side it leaves
  out large, tenuous flakes breaking in turbulence and finite cloud lifetimes,
  and Brown and Francis's mass law carried far past the sizes it was measured at. By that law a 20 cm lunar flake
  weighs about 0.9 g, a fifth of a kilogram of ice to each cubic metre it spans: a loose lace of crystals.

## Assumptions and limits

- The laws for mass, area, fall, capacitance and ventilation were fitted on Earth; they are carried to the Moon through
  the quantities they depend on (gravity, the air's density, viscosity and diffusivity), not refitted.
- Single-crystal seeds use the aggregate mass–dimension law. A habit-specific seed mass and capacitance could change both their dimensions and their subsequent mass gain.
- The columns are straight profiles and the clouds idealised: liquid saturation from top to freezing level, no
  updraughts. The population calculation includes collision chipping; large-aggregate fracture, riming, the melting flake's collapse and sublimation are left out.
- The Moon's ladder is the target study's runs interpolated at 301 K; the design years' own profile may set the
  freezing level a few hundred metres apart.


## History, 21–24 October 2026

Each step of the investigation was recorded here as it happened, newest first. These paragraphs are the record of what
each step computed and the limits it found; where their conclusions differ from the standing at the top of this page,
the standing holds.

The [recent-science check](mcsnow_science_check.md) confirms that the frozen tuned component lacks the observed humidity-dependent transition between compact and complex growth. In 60 short response probes, a 60 µm warm crystal retains an added-volume density near 248 kg/m³ when ice supersaturation falls from 5.24% to 0.26%, although its mass uptake falls about twentyfold. All 24 sampled baseline states reproduce exactly, and the maximum timestep refinement change is 0.0567%. This is a qualitative failure of moisture-dependent morphology, with no new pointwise laboratory error assigned. The author's condition for continuing to freeze and aggregation is unmet; the fit and earlier evidence remain preserved.

The [joint accuracy audit](mcsnow_joint_acceptance.md) now requires every tracked physical error to be at most 10% and no greater than its value in the frozen column-density baseline. That baseline fails 426 of 1,320 comparisons, including a 40.39% warm three-minute fall-speed error. Coupled drag, compact-growth and state-dependent-density experiments improve selected results but fail individual non-regression checks; no replacement model is selected. Constant cold deposition density and a single monotone warm density transition cannot reproduce all the derived local-density targets within 10%. The earlier 9.73% warm mass-and-dimension result remains a limited diagnostic fit; current lunar results remain conditional.

The [warm-column density experiment](mcsnow_column_density.md) reduces the largest thirty-minute mass-or-dimension error from 23.4% to 9.33%, while exactly preserving all 18 checked cold states. Across quarter-minute samples from three to thirty minutes and both sides of the six-minute reference-fit join, the warm maximum is 9.73% for the spherical seed and fitted restart. The fit uses evolving deposition density, a 51.1 µm compact-stage limit and explicit early-age constraints. These diagnostic controls remain temperature-specific: warm thirty-minute speed error rises to 23%, and other initial shapes retain severe failures, including a 4:1 seed runaway. The seventeen fits and their rejected alternatives are retained; no physical replacement or lunar transfer is adopted.

The [independent early/later growth test](mcsnow_growth_limits.md) reduces the cold dendrite's training RMS log error from 0.183 to 0.027. At the reserved thirty-minute age, mass and both dimensions are within 4.8% for both the spherical seed and fitted restart. Warm-column RMS improves from 0.082 to 0.074, while length remains 23.4% and 17.2% high. Constrained refits show that shape alone is insufficient: density and branching matter strongly in the cold case, while restoring original transport and mean ventilation leaves warm thirty-minute mass 128% and 93% high. Widened bounds expose nonunique early transition parameters; unknown seed shape and supplied-state rate errors remain material. A further refit without the assumed warm compact cutoff gives RMS log error 0.098 versus 0.074 with the cutoff, and predicts an unobserved early aspect ratio as low as 0.394 in the sampled history. Without the cutoff, a fixed-mass 4:1 seed gives 89.1 times fitted late length, versus 4.21 times with it. Fifteen searches/refinements, 19,901 parameter-batch evaluations and 48 numerical checks are recorded. These temperature-specific fits diagnose model limitations; no physical replacement or lunar aggregation law is adopted.

The [automatic growth tuner](mcsnow_tuning.md) completes fourteen searches over six diagnostic controls. Separate three-minute restarts can fit mass and both axes closely, but their settings fail when applied from tiny spherical seeds. The shared fit reduces its squared-log training score by 25.5% while retaining large restart errors; fixed-control refits expose compensating parameters. These searches identified early shape evolution and growth allocation as unresolved, with fitted-history and geometry uncertainties retained. No calibrated replacement or full lunar aggregate distribution is adopted.

The [Earth shape and collection diagnosis](mcsnow_shape_diagnosis.md) adds 48 trajectories, 148 fixed-state probes, 36 directional controls and 1,098 native pair evaluations. At the supplied thirty-minute dendrite state, uptake is 3.7% high while thickness growth is 71.8% low. Matching mass leaves shape errors, and a closer spherical-seed result under unit directional ventilation fails the fitted restart. Prescribed-property collection passes independent implementation controls, while predicted crystal properties alter selected pair kernels by factors 0.267–5.218. These are conditional comparisons; the growth law and full lunar aggregate distribution remain unvalidated.

The [archived-lunar-humidity test](mcsnow_weather.md) follows the improved component through 35 hours of ten-minute CM1 weather. From 4 µm spherical seeds, median peak dimension is 0.724 mm versus 5.986 mm in the matched saturation control; 253 of 354 selected paths later reach a tiny sublimation residual. The −20 to 0 °C layers hold a median 25.03 mm of vapour inventory as liquid-water equivalent, but snowy cells have median ice-relative humidity of 89.47%. Rare highly elongated columns reach the declared 1 m censor, and weather cadence affects the results. Twenty executions retain both seed sizes, humidity controls and numerical sensitivities. These are dilute probes of the older A28_dim5 storm, with no population depletion or atmospheric heat feedback.

The [controlled lunar transfer](mcsnow_lunar_growth.md) now follows the improved plate-growth component through the existing 16.28 km snow zone. From 4 µm spherical seeds, the prescribed 3 mm/h allowance gives 3.14 mm crystals after 54.18 hours; tenfold and hundredfold lower incoming number fluxes give 9.51 and 29.08 mm. Continuous saturation without that allowance gives 14.75 cm in 18.07 hours but requires about 2,847 times the imposed water allowance at the original flux. A 0.2 mm seed gives a much more compact final crystal, exposing strong shape-history dependence. All 68 numerical controls complete or reach their declared censor; these conditional trajectories extend the improved component without establishing a lunar snow distribution.

The [growth-rate and shape trace](benchmark_growth_trace.md) now separates the refined archive discrepancy into excess uptake and a smaller residual shape difference. Sixty executions show that restarting from the archived three-minute state retains most of the error. Removing only the extra oblate mean-ventilation term reduces the ten-minute IGF2+ mass error from 49.10% to 9.21%, while preserving almost the same shape at equal mass. Its thirty-minute width is 5.64 mm against a fitted 4.17 mm, with 26.0% excess mass. The improvement does not establish a joint growth repair; native drag remains unchanged.

The [archived-crystal drag trace](benchmark_drag_trace.md) now verifies the selected state through native geometry, three native drag call paths and the authors' own Python helper. Their speeds agree; the saved archive speed remains unexplained. The plotted mass-and-axis comparison still fails when speed is excluded. A solid-envelope substitution approaches two saved speeds while changing mass and failing the HP control, so no drag or growth repair is adopted.

The [Earth benchmark reconstruction](benchmark_reproduction.md) now records 709 native-source controls and duration tests. The recovered author outputs remain unreproduced: saturation, timestep refinement and the documented historical ventilation switch do not resolve the mismatch, and selected archived plate states also give different native fall speeds. Exact water saturation retains late plate runaway; disabling both mean-ventilation operations gives a closer width but deficient mass. The [handoff plan](benchmark_reproduction_plan.md) continues from this checkpoint; no new lunar growth law is adopted.

The [moving-crystal repair](earth_growth_assessment.md#sealed-gas-needs-a-finite-water-treatment) removes a reproduced singularity from sealed voxel air without inventing a water source. The repaired fine Earth history completes with 92.6% of fitted mass but only 33.6% of fitted length. [Resolved initial-flux and attachment controls](earth_growth_assessment.md#local-resolution-and-attachment-limits-narrow-the-remaining-failure) locate the remaining discrepancy in the tested distribution of growth and shape; neither initial-rim adjustment nor uniformly faster basal attachment repairs it. A separate [humidity-inference test](earth_growth_assessment.md#chamber-humidity-can-couple-the-two-parameter-inferences) identifies a possible source of incompatible surface parameters, without adopting a fitted correction.

The [radius and parameter audit](earth_growth_assessment.md#radius-conventions-and-incompatible-step-energy-inferences) rules out a radius-convention change as a complete warm-growth repair in 84 conditional calculations. It also finds that the microscopy-inferred step energy and the broad-facet nucleation barrier cannot be used as one equilibrium parameter without further reconciliation; neither source is silently rescaled.

The [surface-transfer comparison](earth_growth_assessment.md#surface-transport-changes-where-the-incoming-water-becomes-ice) conserves incoming water while testing whether adjoining prism faces can feed a thin tip. Assumed micrometre migration changes length growth substantially with little change in uptake; nanometre migration barely helps. The unmeasured prism transport, later sheath geometry and incomplete growth-history match prevent adoption. Direct step measurements and new nitrogen-pressure evidence sharpen the remaining tests.

The [new facet-law tests](earth_growth_assessment.md#new-facet-laws-do-not-yet-close-the-earth-growth-test) compare published 2023 kinetics with all recorded chamber shapes and rims. They remove the inherited cold symmetry but reject an apparent two-axis match because its prescribed uptake fails local supply conditions. Source location, evolving profiles and surface transport remain unresolved.

The [spatial diffusion control](earth_growth_assessment.md#spatial-diffusion-resolves-finite-rims-and-local-growth) now resolves prescribed cavities and local attachment, with analytic sphere, published cylinder and independent conservation checks. It exposes a substantial thin-ring approximation error; evolving morphology and Earth validation remain open.

The [Earth growth assessment](earth_growth_assessment.md) isolates the late dendritic runaway and narrows the remaining transfer problem. A 240-case ventilation matrix reduces the revised spherical-start crystal from 195.39 to 4.68 mm against a fitted 4.17 mm, but fitted restarts and other habits retain errors. Another 288 humidity and seed-history executions show strong sensitivity to unreported early shape. Prescribing measured dimensions and fall speeds brings conditional mass uptake much closer to the Earth fits; this tests transport while leaving autonomous shape, branching and density evolution unresolved. The reference lunar snow zone overlaps terrestrial reduced-pressure conditions. No absolute lunar size distribution or storm optical forecast is adopted.

The [explicit face-transport controls](earth_growth_assessment.md#explicit-face-transport-rejects-two-tempting-shortcuts) distinguish a hollow-geometry failure from numerical error. Common-corner solid growth fails its own uniform-face condition for the warm column; prescribed thin rims change the rates strongly but do not jointly match mass and both axes. A separate [cold control](earth_growth_assessment.md#the-cold-control-exposes-an-inherited-symmetry-assumption) finds that inherited equal-face parameters force equal axis speeds, contradicted by two measured windows within their reported error bounds. The next repair requires facet-specific kinetics and spatially meaningful hollow growth; no lunar size law is adopted.

The [hollow-crystal assessment](earth_growth_assessment.md#published-hollow-crystal-geometry-and-new-time-series-constraints) reproduces published 2026 geometry routines and adds two measured substrate histories with rim widths. A humidity step widens the median rim by 3.58 times, with explicit measurement-error bounds. Together with the warm column's near-correct uptake but incorrect axis-growth ratio, this supplies a joint geometry and growth test. The complete 2026 growth equations remain unretrieved, and no new growth law is adopted.

The [native McSnow growth comparison](mcsnow_growth.md) records 85 Earth process executions across two source versions, two growth packages and two initial histories. It finds selected improvements but no general repair: all ten dendritic fitted-restart executions fail after uncontrolled growth, while the refined revised spherical-start crystal grows from a fairly close ten-minute result to 195.39 mm at thirty minutes against a fitted 4.17 mm. Mass and both axes must be tested through the full exposure; numerical completion and a plausible short-time size do not establish physical accuracy. Source analysis identifies coupled directional and mean ventilation as a specific control target. No lunar growth closure, size distribution or storm optical forecast is adopted.

The [five-temperature growth comparison and ventilation control](growth_benchmarks.md) find habit-dependent errors in native SCALE: columns are too short and plates too wide relative to published fitted histories, with high endpoint masses. Removing only directional ventilation reduces the −14.4 °C dendrite from 27.66 to 15.81 mm, still above the fitted 4.17 mm; the coupled growth treatment remains unresolved. A [light model and environmental review](growth_options.md) identifies McSnow growth alternatives and a 2026 aggregate-geometry implementation, recovers 169 laboratory measurements with source and data-quality limits, and locates lower-pressure and natural-cloud controls. These results support a process-level comparison before another physical lunar-size or optical forecast.

The [known-limitations review and new controls](scale_limitations.md) identify two separate problems with absolute flake size. Experimental offset contacts raise the lunar mass median from about 1.59 to 2.74 mm and D99 from 3.11 to 7.82 mm; the lunar timestep screen passes, while the Earth distribution screen still fails. More seriously, a direct −14.4 °C Earth single-crystal comparison gives 27.66 mm at 30 minutes against a measured fitted 4.17 mm when initialized from the fitted thin crystal, but 0.96 mm from a sphere with the same initial mass and span. Refinement changes are small compared with that discrepancy. The native axis-growth rule and its ventilation correction retain a strong dependence on starting shape, making individual-crystal growth the next repair target. Neither the original nor experimental absolute lunar sizes are adopted as a physical forecast.

The [size-constraint diagnosis](scale_size_results.md) separates gravity, cloud depth, seed flux and water supply within SCALE. The matched model produces larger lunar snow: approximately 1.59 mm mass-median dimension versus 0.744–0.746 mm in the most refined shallow Earth-gravity controls. Their bulk sizes change by less than 0.8% in the final timestep check, but the exact size-distribution criterion still fails. Lower gravity and greater depth each increase aggregation. The small absolute sizes arise within a shared source and particle treatment: the imposed water ceiling provides about 0.081 mg per injected crystal, while surviving aggregates remain comparatively compact. This rejects a lunar-specific suppression in the tested configuration. That diagnosis did not isolate contact geometry or validate individual-crystal growth; the newer controls above test those two issues.

The [full lunar-column comparison](scale_lunar_results.md) now completes ten SCALE executions and four additional custom controls over the original 16.28 km, 132.59-hour, continuously supplied scenario. SCALE's native prescribed-air cases give mass-median dimensions near 1.6 mm and 99th percentiles near 3.1 mm, with particle, timestep and vertical refinements passing the declared screens. Matching the original stronger sticking law gives approximately 3.56 and 7.47 mm after a shorter-step check removes collision omissions. Custom controls matching ballistic settling and omitting chipping still give 14–18 cm medians. Similar water fluxes therefore accompany very different particle populations. The [methods](scale_lunar_methods.md) distinguish matched forcing from changed particle laws; the result establishes model dependence and a usable numerical alternative, while the original volume-limit diagnosis, physical geometry validation and storm optics remain separate questions.

The [SCALE aggregation controls](analysis.md#scale-sdm-aggregation-controls) now add 19 declared runs, including ice boxes that form aggregates, an analytic liquid-coagulation comparison, time/particle/seed refinements and passive collision counters. Bulk ice statistics are stable over the tested refinements while upper-tail uncertainty remains. The counters find no omitted collisions in the observed ice case; all 789 compared arrays remain identical. The [runner](scale_aggregation.py), [product](results/scale_aggregation.json) and [source record](scale_aggregation_sources.json) distinguish numerical evidence from a reproduced Earth snowstorm and a coupled lunar prediction.

The [SCALE-SDM feasibility trial](analysis.md#scale-sdm-finite-water-column-trial) now completes twelve finite-water Earth/lunar controls with native vapor/heat feedback, evolving ice shape and settling. Total water closes, while the small precipitating fraction remains sensitive to particle count and seed; their saved particle snapshots contain no ice aggregates. The [build](scale_build.py), [column adapter](scale_column.py) and [numerical product](results/scale_column.json) establish a reusable framework, with published Earth population comparison, aggregate fracture and a weather-forced lunar distribution still open. The [existing-solver assessment](analysis.md#existing-solvers-for-the-next-column) and [source register](solver_reuse_sources.json) retain the alternative candidates and their evidence limits.

The current [repair and literature review](repair.md) preserves these earlier tables while testing an independent conservative solver and physical constraints on geometry and breakup.

The [volume-dependence diagnosis](column_diagnosis_results.md) completes 48 full lunar columns. It measures reduced vapour uptake from a sparse population whose demand fluctuates across the local supply cap: deposition rises from 2.199 to 2.498 mm/h across the three statistical volumes under the same 3 mm/h allowance. Removing all collisions passes the tested volume screen. Removing growth, chipping, or both leaves failed upper-tail checks; the sticking-only 99th percentile changes by 22.5%, with appreciable four-realization uncertainty. This identifies a water-allocation contribution and a separate tail concern, while the coupled distribution and physical lunar sizes remain unresolved.

The [full-column event study](event_column_results.md) now completes 512 independent closed-box runs and 28 coupled Earth/lunar controls. The event solver passes the reference comparison and closes the water and number ledgers, while the lunar distribution fails volume convergence: mean mass medians are 12.3, 18.5, 23.5 cm and mean 99th percentiles are 19.8, 36.3, 74.1 cm across successive fourfold volume increases. These are numerical diagnostics under the inherited growth, geometry and chipping laws. They leave the physical 15–60 cm tail and predictive storm optics unresolved.

The earlier [finite-time capture controls](collision_jump_results.md) pass six analytic and six prescribed-bath comparisons. The full column instead extends the independent event reference with donor depletion, continuous growth and the matched upwind settling generator; it is a separate numerical study. [The method note](event_column_methods.md) records its volume, splitting and fractional-fragment approximations.

The [21 October analysis](analysis.md) extends this work to snow-conditioned cloud motion, corrected aggregate geometry, full size distributions and lunar snowstorm optics. The tables above retain the 20 October baseline; the source-definition corrections and current limitations are recorded in the new analysis. The inherited fractional-particle bookkeeping fails an arbitrary-area invariance requirement; read the [numerical audit](numerical_audit.md) before interpreting either set of size estimates.

### Next steps recorded on 20–22 October

The same fall and growth physics gives the high cloud's crystals at lunar gravity for the GCM
([plan](../../research/plan.md), item 5). The [new analysis](analysis.md) measures the archived cloud environments and calculates flake-scale load scenarios. The independent [sectional reference](repair.md) now provides consistent particle accounting, but its first three full-column controls fail the finite size-domain check. Numerical convergence, a coupled moist trajectory, a spread of fall speeds among like flakes and a measured mechanical failure law remain needed to predict a survival-limited size distribution.

The [implicit facet component](facet_coupling.py) and [432-state audit](results/facet_coupling_audit.json) now solve attachment, surface humidity and latent heating together under explicit ventilation controls. Numerical balances close; autonomous habit and branching remain under assessment in the [Earth growth assessment](earth_growth_assessment.md).

The [subsequent Earth trajectory controls](earth_growth_assessment.md#earth-trajectories-identify-a-remaining-shape-failure) complete 204 integrations and isolate a column-shape failure even when fitted intermediate states are supplied: conditional mass growth is close while relative long-axis growth remains too slow. A newer empirical humidity screen identifies where compact growth gives way to substantial hollowing or branching; no replacement density law is selected.
