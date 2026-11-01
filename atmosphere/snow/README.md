# Snow

The snow work asks four things of the Open Moon's sky winter: how large snowflakes grow, how they fall, whether the largest survive the cloud's turbulence, and how a snowstorm looks under lunar light. The questions opened on 20 and 21 October as the second half of the work on the high cloud's crystals ([register](../../research/decisions.md), 20 and 21 October), and the answers to the first three now enter the climate model, where PlaSim's melted rain takes the size lunar snow melts into and its snow sublimates by a closure built from the lunar flakes (register, 29 and 30 October). How flakes fall, how deep the snow zone stands and how far below the freezing level snow lasts are calculated in the sections below.

Without any breakup, crystals that grow and aggregate together carry the 90th percentile of flake span at the lunar freezing level to 4.1–4.9 cm in light snow of 0.3 mm/h and to 12.6–12.7 cm at 3 mm/h with a hundredth of the crystals ([growth and aggregation](mcsnow_growth_aggregation.md)), in pairs of superparticle numbers that converge, on the lunar McSnow build corrected on 25 October and rerun in full on 29 October ([diffusivity](mcsnow_diffusivity.md)). A comparison of loads with Earth's largest observed flakes caps the flakes that break at about 5 cm for ordinary flakes and 8–15 cm for the largest ([survival](turbulence_survival.md)), and a calculation of the loads on a falling flake's junctions puts the lunar equivalent of Earth's 6 cm flakes at 11.7–13.5 cm ([breakup network](breakup_network.md)). The sizes reached without breakup and the cap meet in moderate snow. Light snow stays below the cap when it forms at −15 °C and reaches it when it forms near −30 °C. In heavy snow at the reference crystal number, and in snow of 1 mm/h from tops near −30 °C, McSnow without breakup grows aggregates metres across, so since 30 October its lunar build caps them at 14 cm, the load analogy's lunar equivalent of Earth's largest observed flakes, in place of the breakup law no Earth measurement yet supports. The capped columns are the reference from which PlaSim's snow takes its sublimation below a cloud base ([sub-cloud test](mcsnow_subcloud.md), [closure](results/subcloud_closure.json)). No snow reaches any ground: below a freezing level at 22–26 km the flakes sublimate, melt and fall the last 15–25 km as rain ([melting](snow_melting.md)).

On 24 October the author made three decisions about the course of the work. The crystal-growth model is frozen as it
stood that day, and its errors against Earth's laboratory crystals enter every lunar result as a stated uncertainty. A
rule sized to the laboratory data replaces the 1,320-point accuracy rule. The work then moves to the lunar answer, in
the order of the plan below.

The first three steps of the plan are done. Judged crystal by crystal against the scatter between repeat crystals, the
two-temperature calibration passes on mass and axes, with 5 of its 93 comparisons outside where chance gives about 4,
and fails on the warm column's fall speed ([data-sized check](data_sized_check.md)). The improved component, which every
lunar run uses, fails the check: 83 of its 332 mass and axis comparisons fall outside, and it falls 1.1–2.0 times faster
than the laboratory crystals at most temperatures. The band carried through the lunar column is that component's
measured error. At the band's edges the component moves single-crystal size by a factor of 1.4–1.6 wherever water limits
growth, against a factor of 9.3 from a hundredfold change in crystal number and 8.2 between archived storm humidity and
continuous saturation ([Earth error band](earth_error_band.md)). Only under unlimited saturation, which would need about
2,800 times the water the column allows, do the errors compound to a factor of 5.4. Across light to heavy snow and three
decades of crystal number, the crystal arriving at the freezing level follows the water each crystal receives, from
0.4 mm to 14 cm, with the crystal model's errors moving sizes above 2 mm by 15–31% and the sub-millimetre crystals by
up to 41% ([size range](lunar_size_range.md)). Above a centimetre those single crystals extrapolate far past any
measured crystal, and aggregation and breakup, the next two steps, decide what falls.

Steps 4 and 5 gave the first flakes. McSnow's aggregation, fed with crystals of measured size and run with Connolly and
co-authors' (2012) sticking law, reproduces the median size of Earth's snowflakes within 16% of two observed relations
from 0.3 to 9 mm/h, while its largest flakes outgrow the observed ones tenfold in the heaviest snow
([aggregation](mcsnow_aggregation.md)). Fed with full-size lunar crystals at the cloud top, the same aggregation carried
the 90th percentile of flake span past 10–14 cm within 1.5 km of the top at 3 mm/h, within 4.5 km at 1 mm/h and in the
lower half of the column at 0.3 mm/h. Under the same turbulence a lunar flake carries no more load at its junctions than
an Earth flake until it is 2.2–2.6 times larger, so Earth's ordinary 2 cm flakes and its 4–6 cm storm aggregates
correspond to lunar flakes of about 5 cm and 8–15 cm ([survival](turbulence_survival.md)). That analogy assumes that
breakup limited Earth's largest flakes. A first contact-network calculation on 30 October
([breakup network](breakup_network.md)) finds the bending moment at a flake's most loaded junction growing faster than
its weight, so that a junction failing in bending, as a narrow ice neck does, puts the lunar equivalent of Earth's
4–6 cm flakes at 7.8–11.7 cm with the upper crystals in the lower ones' shadow and 9.0–13.5 cm without, against the
analogy's 9.4–14.1 cm for Mitchell's aggregates; a breaking flake keeps about two thirds of its mass in one piece.

On 25 October the author approved three further steps, and all three are done. No law for whole-flake breakup has Earth
support ([breakup literature](breakup_literature.md)). The measurements count the fragments chipped from centimetre
flakes struck by dense graupel at 10⁻⁷–10⁻⁵ J, or record the speed, about 2 m/s, at which flakes blown against a snow
surface begin to fracture. Collisions in free fall are gentler than any of these: lunar flakes close on crystals and on
one another at 0.5 m/s or less, with energies near 10⁻⁸ J and about 3 × 10⁻⁸ J in strong turbulence. That finding stands
beside the analogy's premise and questions whether anything breaks lunar flakes aloft, so the sizes reached without
breakup form the second bound. Grown from 4 µm seeds while they aggregate, the flakes reach the analogy's span far lower
than full-size crystals did. On the corrected build, rerun in full on 29 October, the 90th percentile at 3 mm/h passes
10 cm about 6.5 km below the −15 °C top and 14 cm about 8.5–9.5 km down, against 1.5 km for both with full-size
crystals; with a hundredth of the crystals it passes 10 cm 11.5–12.5 km down and arrives at 12.6–12.7 cm; and light snow
arrives at 4.1–4.9 cm, both sizes converged under four times the superparticles. On the earlier build these were about 8
and 10 km, 14.5 km, 10.5–10.9 cm and 4.2–4.5 cm, and the crystal model's Earth error band moves those depths by up to 4
km. As crystals gather into flakes they take up less of a steady supply, 46–52% of it at 3 mm/h (40% on the earlier
build), and a cloud would hold the rest as supercooled water. Along the archived storm's paths most crystals sublimate
before the freezing level, and only the paths held at saturation with dilute crystals give a converged size, about 13–15
cm on the median path. The design-climate storm boxes hold half their snow within 18–20 km of the freezing level and
nine-tenths within 30–35 km ([snow zone](snow_zone_depth.md)). Light snow from the reference column's −15 °C top
therefore melts before it reaches the size limit. Light snow forming near −30 °C, higher in that snow zone, falls twice
as far and arrives with a median span of 5.3–6.5 cm, a 90th percentile of 12–23 cm and a 99th of 28–44 cm, from runs at
two superparticle numbers that agree on the median within 21% and part nearly twofold at the 90th percentile, a tail a
finer run would settle. It reaches the size limit, so whether light snow does depends on where in the snow zone it
forms.

Below the freezing level the flakes sublimate, then melt, and then fall as rain for 15–25 km, which the author asked to
follow on 25 October with sublimation carried all the way down ([melting](snow_melting.md)). In the design-climate storm
boxes the freezing level stands 22–26 km up, and the air under falling snow is 93–97% humid there and warms about 1 K
for each kilometre below it, so melting begins 0.3–0.55 km down. Flakes of 2 cm finish melting about 1.6 km below the
freezing level and flakes of 10–20 cm 3.7–5.2 km below it, three to four times as deep as the flakes of 20 October.
Their 9–15 mm drops hold together under lunar gravity and bring 81–93% of their water to the surface, while the drops of
flakes up to 2 cm evaporate on the way. Light snow of 0.3 mm/h therefore sends 10–15% of its water to the surface, 6–7%
with a tenth of the crystals, and snow whose flakes reach 11–13 cm 56–65%, against the 39–40% CM1's own storms deliver
of their falling water; no snow reaches the ground, even at the highest summit. Sublimation takes 2–8% of the snow over
land and 9–23% over the sea, and at 2.83 MJ/kg against 0.33 for melting it draws about as much heat from the air over
the sea as the melting does in the heavier snow, and more than twice as much in light snow. On Earth the same code
begins and ends the melting of ordinary flakes where Matsuo and Sasyo (1981) computed and where observed melting layers
lie. Those shares take the corrected runs' flakes, as every table of the melting record now does
([diffusivity](mcsnow_diffusivity.md)); the earlier build's flakes gave 11–13% and 55–60%.

## The plan, in order

1. **Check the frozen model by the data-sized rule** (below), naming every measured crystal that falls outside it.
   Done 24 October ([results](data_sized_check.md)).
2. **Carry the Earth errors through the lunar column.** Run the frozen model at the edges of its Earth error band in
   mass growth, axis growth and fall speed, through the archived lunar storm's humidity and the representative
   moisture cases, and set the change in lunar flake size beside the spread from moisture supply and crystal number.
   This step tests whether the crystal model's Earth errors matter for the Moon. Done 24 October
   ([results](earth_error_band.md)).
3. **Give lunar flake sizes as a range** over moisture supply and crystal number. Done for single crystals on
   24 October ([results](lunar_size_range.md)); flakes follow from step 4.
4. **Join the flakes, checked first on Earth.** With measured crystals as input and established aggregation
   efficiencies, the aggregation calculation is compared with the sizes of observed Earth snowflakes, then run through
   the lunar column. Done 24 October ([results](mcsnow_aggregation.md)).
5. **Estimate survival in turbulence**: whether the largest aggregates hold together in the motions of the cloud that
   grows them ([register](../../research/decisions.md), 21 October). Done 24 October as a load analogy
   ([results](turbulence_survival.md)).

Steps 4 and 5 start from the crystals of step 3 and carry the uncertainty measured in step 2. The 25 October stage
followed:

6. **A breakup law in place of the load analogy.** No law with Earth support exists for whole flakes breaking in
   collisions or under drag and turbulence; the analogy stays as the cap, with its limits stated. Done 25 October
   ([results](breakup_literature.md)). A first contact-network calculation of junction loads followed on 30 October
   ([results](breakup_network.md)): in idealised aggregates bending decides, the lunar breaking span comes out 1.95–2.25
   times Earth's where the analogy gives 2.35, and the larger piece of a breaking flake keeps a median 0.64–0.67 of its
   mass. It rests on geometric drag rules and an Earth calibration, with no measured junction strength.
7. **Growth and aggregation together**, through the reference column's water supply and along the archived storm
   paths, at the crystal-number cases and the error band's edges. Done 25 October
   ([results](mcsnow_growth_aggregation.md)).
8. **Light snow and melting**: how deep the storm boxes' snow zone stands above the freezing level, and whether light
   snow reaches the size limit before it melts. Done 25 October ([results](snow_zone_depth.md)).
9. **Below the freezing level**: how far the flakes fall before they melt and their water evaporates, with sublimation
   carried from the cloud to the ground, the snow's water and latent heat by layer, and the rain each ground height
   receives. Done 25 October ([results](snow_melting.md)).
10. **The lunar builds' vapour diffusivity**: McSnow's fit given back its 1 atm reference, which the lunar column's
   surface-pressure patch had moved, and a minimal set of runs repeated on the corrected build. Done 25 October
   ([results](mcsnow_diffusivity.md)). The remaining eighteen runs, held on 26 October so the compute could go to the
   climate model's levels and convection scheme, were repeated on 29 October, and the growth-and-aggregation and
   melting products now read the corrected runs.
11. **Snow below a cloud base, for the climate model**: McSnow's column with the air below a cloud base held over ice, the reference for replacing PlaSim's snow sublimation with the lunar flakes' own. A closure that takes the share of the snow's flux sublimated per metre by depth below the snow's top from reference columns matched by the snow's rate at the cloud base, and raises it as the flakes shrink, follows McSnow within 0–9% over both design boxes in snow of 1 and 0.3 mm/h. Carried to the boxes' own range, with references from −30 °C beside those from −15 °C and McSnow's aggregates capped at 14 cm, it holds on nine of the eleven resolved box runs, the land's 1 mm/h snow from −30 °C included, and overstates the loss only in 0.1 mm/h snow from −15 °C, by 1.48 and 1.65 of McSnow's e-foldings. Done 29–30 October ([results](mcsnow_subcloud.md)), and adopted on 30 October with that miss stated as its limit; its tables are the product PlaSim's build reads ([closure](results/subcloud_closure.json)).

What the snow work leaves open:

- **Where flakes break.** No measurement of whole flakes breaking, flake on flake or under sustained load, supports a breakup law, and the 14 cm cap stands in its place; with the references capped at 10 cm the closure's sublimation rises by 14–15% in the land's 1 mm/h snow from −30 °C and by 0.6% or less elsewhere. The junction-load calculation, which puts the lunar equivalent of Earth's 6 cm flakes at 11.7–13.5 cm, between the cap and its 10 cm sensitivity, rests on idealised aggregates, geometric drag rules and an Earth calibration, with no measured junction strength.
- **The tails.** The 90th and 99th percentiles of light snow from −30 °C move nearly twofold and by 58% when the superparticles are halved. Heavy snow at the reference crystal number runs away to metres without the cap and gathers just under it with the cap, 89–120 mm across at the boxes' cloud bases.
- **The closure's limits.** It overstates the loss in 0.1 mm/h snow from −15 °C by about five points of the base flux; four capped box runs, which lose 19–25% of their flux below the cloud base, stay unresolved; a top above 32.56 km takes the −30 °C family untested; and k and q keep the 301 K design ladder's temperatures ([sub-cloud test](mcsnow_subcloud.md)).
- **A storm whose vapour responds to its crystals**, in place of the reference columns' prescribed supply.
- **Riming.** McSnow's superparticles only fall, so its riming needs an updraft added and checked on Earth before it can give the graupel's sublimation, which PlaSim takes meanwhile from CM1's graupel (register, 29 October).
- **Wet flakes sticking together below the freezing level**, which the melting calculation leaves out ([melting](snow_melting.md)).
- **The storms' optics** under lunar light, beyond the 21 October extinction study ([optics](optics.md)).
- **The cap's cost.** With the cap set, McSnow computes every superparticle's geometry afresh at each step, which costs up to 3% of a run's time; trimming that cost comes before the next long McSnow batch.
- **The cold cloud.** The same fall and growth physics gives the high cloud's crystals at lunar gravity ([below](#the-high-clouds-crystals)), and the GCM carries them ([cold_cloud.py](../../climate/gcm/cold_cloud.py)); their aggregation, heterogeneous ice nuclei and the cooling of gravity waves are not yet modelled.

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

The calibration is the Earth benchmark result. The lunar runs, and the uncertainty carried with them, come from the
improved component, the only part of the model that runs through a lunar column; its Earth errors are the larger of the
two, so the band it sets is the wider.

The component diffuses vapour by Pruppacher and Klett's fit, 0.211 cm²/s at 0 °C and 1 atm, rising as the 1.94 power of
temperature and falling inversely with pressure (McSnow's `diffusivity_vapour`). It stands 3.4% below Massman's (1998)
measured fit at −5 °C and 4.6% below it at −30 °C, where Massman's fit is extended below its 0–100 °C data, and that gap
belongs to the component's stated uncertainty. In the lunar builds the patch that sets McSnow's surface pressure to the
freezing level's, 659 hPa, also sets that fit's 1 atm reference, so every lunar McSnow run before 25 October diffused
vapour at 0.651 of the fit at every height. Where diffusion limits growth, lunar crystals then grew at 0.68–0.79 of the
fit's rate, slowest in the coldest air; the Earth builds keep the 1 atm reference, so the calibration and the data-sized
check are unaffected. With the reference restored in the lunar build ([diffusivity](mcsnow_diffusivity.md)), light snow
at the reference flux, from −15 °C or near −30 °C, moves by under 5% at the median and the 90th percentile, while at 3
mm/h the crystals take up 24–37% more of the supply and the flakes arrive 11–26% larger at the median and the 90th
percentile, reaching 10 cm 2–3 km higher in the column.

Before the freeze, the component's −14.4 °C dendrite carried 26% too much mass at 30 minutes against Takahashi and
co-authors' (1991) fitted laboratory history, with a span 35% too wide and a thickness 17.5% too thin; its −5.3 °C
column grew too short, 31% below the fitted length at the fitted mass ([shape diagnosis](mcsnow_shape_diagnosis.md)).
The calibration brings the warm column within 9.73% of the fitted mass and dimensions from 3 to 30 minutes, and the cold
dendrite within 4.75% at 30 minutes with its 3-minute mass about 17% low; the warm column's fall speed stays 21–23% high
at 30 minutes. These agreements are fits to the same histories ([joint audit](mcsnow_joint_acceptance.md)). Both parts
share three recorded limits. Crystals started from elongated seeds run away. The shape and density rules ignore humidity
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

Repeat crystals scatter by 0.193 in ln mass and by 0.159 and 0.183 in the ln equatorial and polar axes, so the rule's
two-sigma bands span about 1.4–1.5 times either way. Under it the calibration has 5 of 93 mass and axis comparisons
outside, close to the 4 that chance gives, and its −5.3 °C fall speed fails from 12 minutes, up to 6.1 cm/s fast. The
component has 83 of 332 comparisons outside and misses the speed accuracy at 55 of 78 ages, mostly by falling too fast
([data-sized check](data_sized_check.md)).

## Files and commands

The snow solver sits beside the raindrop solver ([../rain](../rain/README.md)). The table lists the principal files;
the index of the 21–24 October investigation below links the rest. Material the work kept in the laptop's /tmp until 24 October is copied to the
project-space drive and reached through the ignored link `research/runs/snow_tmp_copy`;
[results/tmp_relocation.json](results/tmp_relocation.json) maps each /tmp path a committed record names to its copy and
hash, and lists the 52 transient paths the copy does not hold. New work keeps nothing it needs in /tmp.

Thirty-three products of 21–24 October were built with the rain solver as it stood before 25 October, whose
Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998) measured value, which the solver has taken since:
the collision benchmarks and their checks (`collision_entropy_*`, `collision_jump_*`, `collision_moment_*`,
`collision_reference`), the column diagnoses (`column_diagnosis*`), the event column (`event_column*`,
`event_closed_boxes`), the SCALE comparisons (`scale_crystal_growth`, `scale_growth_benchmarks`,
`scale_growth_controls`, `scale_lunar_comparison`, `scale_separation`, `scale_size_diagnosis`),
`crystal_transport_replay`, `geometry_continuation`, `growth_depth`, `growth_repair`, `growth_transport`,
`kernel_diagnostics`, `limited_growth`, `numerical_audit`, `optics`, `population`, `repair` and `turbulence`. The
current standing rests on none of them. They keep their values; their records name `atmosphere/rain/drops.py` at that
solver's hash, kept byte for byte as
[atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and their identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
The flake product of 20 October and the aggregation baseline were rebuilt with the corrected solver.

| File | What it holds |
|---|---|
| [data_sized_check.md](data_sized_check.md), [data_sized_check.py](data_sized_check.py), [results/data_sized_check.json](results/data_sized_check.json) | The frozen model against each measured crystal: repeat-crystal scatter, every verdict, the crystals outside |
| [earth_error_band.md](earth_error_band.md), [earth_error_band.py](earth_error_band.py), [results/earth_error_band.json](results/earth_error_band.json) | The component's Earth errors as three multipliers, their sizing, parity with the frozen component, and lunar sizes at the band's edges |
| [lunar_size_range.md](lunar_size_range.md), [lunar_size_range.py](lunar_size_range.py), [results/lunar_size_range.json](results/lunar_size_range.json) | Single-crystal sizes at the freezing level over water supply and crystal number, with the band's edges |
| [mcsnow_aggregation.md](mcsnow_aggregation.md), [mcsnow_aggregation.py](mcsnow_aggregation.py), [results/mcsnow_aggregation.json](results/mcsnow_aggregation.json) | McSnow aggregation from measured crystals checked against observed Earth snowflake sizes, and the lunar column's depths of outgrowth |
| [turbulence_survival.md](turbulence_survival.md), [turbulence_survival.py](turbulence_survival.py), [results/turbulence_survival.json](results/turbulence_survival.json) | Lunar flake spans with the same junction load as Earth's largest observed flakes, over turbulence and mass laws |
| [mcsnow_growth_aggregation.md](mcsnow_growth_aggregation.md), [mcsnow_growth_aggregation.py](mcsnow_growth_aggregation.py), [results/mcsnow_growth_aggregation.json](results/mcsnow_growth_aggregation.json) | Crystals that grow and aggregate together in the reference column and along the storm paths: depths of outgrowth, arriving sizes, convergence, the band's edges and the water the snow takes |
| [snow_zone_depth.md](snow_zone_depth.md), [snow_zone_depth.py](snow_zone_depth.py), [results/snow_zone_depth.json](results/snow_zone_depth.json) | How far above the freezing level the design-climate storm boxes hold their snow, and light snow against the size limit |
| [snow_melting.md](snow_melting.md), [melting.py](melting.py), [melting_environment.py](melting_environment.py), [results/snow_melting.json](results/snow_melting.json) | How far below the freezing level flakes melt and their drops evaporate in the storm boxes' air; sublimation, melting and evaporation by layer with their latent heat; the Earth check; the storm paths' losses above the freezing level |
| [mcsnow_diffusivity.md](mcsnow_diffusivity.md), [mcsnow_diffusivity.py](mcsnow_diffusivity.py), [results/mcsnow_diffusivity.json](results/mcsnow_diffusivity.json) | McSnow's vapour diffusivity in the lunar builds: the corrected build, its source differences, and a minimal set of runs repeated on it against the old runs |
| [mcsnow_subcloud.md](mcsnow_subcloud.md), [mcsnow_subcloud.py](mcsnow_subcloud.py) | McSnow's column with the air below a cloud base held over ice and its aggregates capped at 14 cm: the two declared changes, the design boxes' sub-cloud air, the cases, and the five tests of the climate model's snow sublimation against it |
| [melting_sources.json](melting_sources.json) | The 25 October melting sources: the quantities taken from each and how each was read |
| [breakup_literature.md](breakup_literature.md), [breakup_sources.json](breakup_sources.json) | The 25 October search for an Earth-supported breakup law: the measurements found, their quantities and access, and the collision energies aloft |
| [breakup_network.md](breakup_network.md), [breakup_network.py](breakup_network.py), [results/breakup_network.json](results/breakup_network.json) | Junction loads in model aggregates under drag and weight: how the most loaded junction's force and moment grow with span, the lunar breaking span they imply, and the largest remnant |
| [subcloud_closure.py](subcloud_closure.py), [results/subcloud_closure.json](results/subcloud_closure.json) | The adopted sub-cloud closure for PlaSim's snow: the capped references' tables by depth, the rule between the families, a reference implementation and its test vectors |
| [tmp_relocation.py](tmp_relocation.py), [results/tmp_relocation.json](results/tmp_relocation.json) | Each /tmp path the committed records name, mapped to its hash-checked copy on the project-space drive |
| [scale_size_results.md](scale_size_results.md), [scale_size_diagnosis.py](scale_size_diagnosis.py), [scale_size_sources.json](scale_size_sources.json) | Gravity/depth attribution, water-per-crystal and aggregation controls, active-source audit and limits of physical interpretation |
| [scale_lunar_results.md](scale_lunar_results.md), [scale_lunar_methods.md](scale_lunar_methods.md), [scale_lunar_comparison.py](scale_lunar_comparison.py) | Complete SCALE/custom lunar comparison, numerical refinements, sticking and atmospheric-feedback contrasts |
| [flakes.py](flakes.py) | The snowflake's mass, area and fall law, its growth from the vapour, its melting, the winter columns and the product |
| [results/snow.json](results/snow.json) | Fall speeds by size, growth across the winter cloud, crossing times and melting depths, Moon and Earth |
| [aggregation.py](aggregation.py) | A column of falling snow: super-particles that grow, stick together and break, and the snow they bring to the freezing level |
| [results/aggregation.json](results/aggregation.json) | The earlier arriving-snow sizes by supply and crystal numbers, retained as a baseline |
| [cirrus.py](cirrus.py) | A slab of cold air rising through a layer of the upper air: its haze freezing into crystals that grow and fall out, the Moon against Earth, for the GCM's cold cloud ([cold_cloud.py](../../climate/gcm/cold_cloud.py)) |
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

cirrus.py has no product of its own; [cold_cloud.py](../../climate/gcm/cold_cloud.py) runs it on the GCM's layers.

## The physics

A snowflake is an aggregate of ice crystals. Its mass follows Brown and Francis (1995), m = 0.0185 D^1.9 kg, and its projected area is 0.35 of the circle across D. It falls by Heymsfield and Westbrook's (2010) law, in which gravity enters through the modified Best number alone,
X* = 8 m g ρ / (π Ar^½ η²), with Re = (δ0²/4)((1 + 4 X*^½ / (δ0² C0^½))^½ − 1)², δ0 = 8.0, C0 = 0.35 and the speed
η Re / (ρ D).

The [source review](habits.md) of 21 October found Brown and Francis's D to be the mean of two orthogonal image dimensions, which on the maximum dimension gives m = 0.0121 Dmax^1.9 and spans about 25% larger at a given mass, and found the area ratios cited for the flake, Seifert and co-authors' (2026), taken against a circumscribing ellipse: against the circle across Dmax their aggregates of 1,024 crystals give 0.25–0.31. The flake model keeps the original closures, and the population analysis of 21 October carried the corrected ones as sensitivities.

The flake grows from the vapour by diffusion to a capacitance of 0.25 D (Westbrook, Hogan and Illingworth 2008), ventilated as
Hall and Pruppacher (1976) give it for ice, against Murphy and Koop's (2005) saturation over ice and supercooled water.
Below the freezing level it melts by the heat it takes from the air less what its wet surface's evaporation carries
away (Pruppacher and Klett 2010, chapter 16), at the same capacitance and ventilation. The air's viscosity,
conductivity and vapour diffusivity are the raindrop solver's (`drops.transport`) for each world's air, the
diffusivity Massman's (1998) measured fit since 25 October.

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

**So snow grows far more.** A 0.2 mm crystal crossing the cloud gains from the vapour alone 2.5, 14.8 and 53 mg on the
Moon from the three tops, against 0.011, 0.029 and 0.058 mg on Earth: 200–900 times the mass, as a growing crystal
grows faster the larger it is. With five times the air in the cloud, a flake also has about five times the cloud's
small crystals to sweep up at one mixing ratio. A single flake sweeping a cloud that never runs out grows without
bound; how large flakes grow in a whole population, with a breakup law, follows below.

**Snow lasts far below the freezing level.** Below it the air warms only 0.92 K for each km a flake falls, so flakes
melt slowly. Where all of a flake's ice has melted, metres below the freezing level:

| Flake | Moon, air at 90% | Moon, at 70% | Earth, at 90% | Earth, at 70% |
|---|---|---|---|---|
| 2 mm | 1,119 | 3,119 | 237 | 487 |
| 5 mm | 1,209 | 3,211 | 274 | 525 |
| 1 cm | 1,284 | 3,287 | 306 | 557 |
| 2 cm | 1,364 | 3,370 | 340 | 592 |

In drier air a flake's own evaporation keeps its surface below freezing for longer. The sky winter's snowfall is a
band kilometres deep: snow forming in a cloud above the freezing level reaches decks 1.1–3.4 km below it. The [melting
record](snow_melting.md) of 25 October replaces these depths: it carries the flake's sublimation before and while it
melts, Mitra's melting flake, the large flakes and the storm boxes' air.

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
| 0.3 mm/h | 10 a litre | 10; 34; 52 mm | 0.7; 1.3; 1.4 mm |
| 0.3 mm/h | 100 a litre | 18; 64; 96 mm | 0.6; 1.5; 2.2 mm |
| 1 mm/h | 10 a litre | 51; 151; 195 mm | 0.9; 1.9; 2.4 mm |
| 1 mm/h | 100 a litre | 61; 184; 233 mm | 0.8; 3.1; 4.5 mm |
| 3 mm/h | 10 a litre | 142; 436; 572 mm | 0.9; 1.9; 2.4 mm |
| 3 mm/h | 100 a litre | 145; 526; 662 mm | 1.8; 7.2; 18 mm |

The snowfall is in millimetres of water an hour. On Earth 10 crystals a litre take up only 0.61–0.62 mm/h of the cloud's
vapour in the two heavier cases, the rest left unused; the arriving snow also counts the crystals' own mass entering the
top, about 0.4 mm/h at 100 a litre on Earth.

- **The breakup law does not cap the flakes.** With it and without it the medians agree within 11% and the upper sizes
  within 40%, and not always smaller with it. Phillips's law chips small fragments off the smaller particle when a
  meeting does not stick, which feeds the cloud more small crystals but splits no large flake. On both worlds the
  flakes' size is set by the snowfall's water and the time the snow takes to fall.
- **Lunar flakes grow 15–90 times the size of Earth's at the same settings:** a centimetre or two in light
  snow, 5–15 cm in moderate snow and tens of centimetres in heavy snow. The smallest crystals take 1–1.5 days to fall
  through the Moon's 16 km of cloud against about an hour on Earth, in a model whose sticking probability depends on temperature and whose encounters use differential settling alone.
- **How far to trust the sizes.** The model's Earth flakes run small against Earth's real snow, where flakes of
  several millimetres to a couple of centimetres are common, so it aggregates too little, most likely for the spread of
  fall speeds among like flakes, which it leaves out; that alone does not determine the sign of the lunar size bias. On the other side it leaves
  out large, tenuous flakes breaking in turbulence and finite cloud lifetimes,
  and Brown and Francis's mass law carried far past the sizes it was measured at. By that law a 20 cm lunar flake
  weighs about 0.9 g, a fifth of a kilogram of ice to each cubic metre it spans: a loose lace of crystals.

## The high cloud's crystals

[cirrus.py](cirrus.py) follows a slab of cold air rising through a layer of the upper air, for the GCM's cold cloud
([plan](../../research/plan.md), item 5). The slab starts at ice saturation at the layer's bottom and rises at a steady
pressure velocity, cooling on its dry adiabat. Its haze of sulphate droplets freezes at Koop and co-authors' (2000)
homogeneous rate; the crystals grow from the vapour with the gas-kinetic limit of small crystals, at a deposition
coefficient of 0.5, fall at Heymsfield and Westbrook's speed and leave through the slab's base. A crystal is a solid
ice sphere up to 97 µm, where Brown and Francis's mass law meets the sphere's, and takes that law beyond, its area
ratio falling to the aggregates' 0.35 at 1 mm.

Two worlds rising at the same pressure velocity cool their air alike and condense the same vapour per kilogram. On the
Moon a crystal falls more slowly, in proportion to gravity while it is small, and a pascal of the slab is six times
deeper, so its crystals, though they grow larger, stay ten to fifteen times longer. [cold_cloud.py](../../climate/gcm/cold_cloud.py) runs the slabs
through the GCM's top three layers in the design years, at the ascent at which Earth's slab holds the GCM's own Earth
water path (0.053–0.058 Pa/s, about 2 cm/s on Earth):

| Layer of the GCM | Slab's temperatures | Earth: ice path; radius; ice fallen out | Moon: ice path; radius; fallen out |
|---|---|---|---|
| Above 107 km | 180 to 148 K | 0.32 g/m²; 1.4 µm; 1% | 2.0 g/m²; 1.4 µm; 0% |
| 81–107 km | 208 to 169 K | 2.1 g/m²; 18 µm; 97% | 77 g/m²; 31 µm; 30% |
| 63–81 km | 234 to 205 K | 6.1 g/m²; 51 µm; 98% | 258 g/m²; 82 µm; 92% |

- **Earth's cirrus loses nearly all its ice; the Moon's keeps much of it.** Below 107 km Earth's crystals fall out of
  the slab within hours of forming, so it holds only a few per cent of the ice it has condensed. The Moon's
  stay, so per kilogram of air its cloud holds six to seven times the ice and, with six times the air in each pascal,
  37–42 times the path.
- **Fewer, larger crystals.** Crystals that stay keep the supersaturation down, so the haze freezes less often: 37
  crystals a litre at 81–107 km against Earth's 255, and radii 1.6–1.7 times Earth's.
- **The top layer freezes its haze at once on both worlds.** Its air cools to 148 K, nearly every droplet freezes,
  and crystals of 1.4 µm hardly fall on either world: the same ice per kilogram, and six times the path on the Moon.

The factors move little with the deposition coefficient (0.1 to 1) or the haze's number (a fifth to five times):
water 30–42 times, radius 1.6–2.3 times. At the same ascent Earth's slab gives radii 1.45–1.6 times Kristjánsson and
Mitchell's, the GCM's Earth crystals; the GCM takes their radii times the factors, so that bias largely cancels. Left
out: aggregation, which on the Moon brings crystals together more slowly but for longer, so that a crystal's chance of
joining another over its life should change little (an estimate, not run); heterogeneous ice nuclei; gravity waves in
the cooling, which raise the number of crystals on Earth; and ice falling in from above.

## Assumptions and limits

- The laws for mass, area, fall, capacitance and ventilation were fitted on Earth; they are carried to the Moon through
  the quantities they depend on (gravity, the air's density, viscosity and diffusivity), not refitted.
- Single-crystal seeds use the aggregate mass–dimension law. A habit-specific seed mass and capacitance could change both their dimensions and their subsequent mass gain.
- The columns are straight profiles and the clouds idealised: liquid saturation from top to freezing level, no
  updraughts. The 20 October population column includes collision chipping; large-aggregate fracture, riming, the melting flake's collapse and sublimation are left out.
- The Moon's ladder is the target study's runs interpolated at 301 K; the design years' own profile may set the
  freezing level a few hundred metres apart.

## The investigation of 21–24 October 2026

Before 24 October the calculations put typical flakes anywhere from below a millimetre to tens of centimetres across, a spread carried by the water a storm supplies to growing ice, the number of crystals that share it, and the population solver that follows their collisions. The investigation tested population solvers, a second microphysics model and the crystal model against Earth's laboratory crystals until the crystal model was frozen on 24 October. Each step keeps its own write-up; where a write-up's conclusions differ from the standing at the top of this page, the standing holds.

| Strand | Write-ups | What they established |
|---|---|---|
| The 20 October population column and its solvers | [analysis](analysis.md), [numerical audit](numerical_audit.md), [sectional reference](repair.md), [growth repair](growth_repair.md), [collision and exposure](collision_and_exposure.md), [capture controls](collision_jump_results.md), [event column](event_column_results.md) and its [methods](event_column_methods.md), [volume diagnosis](column_diagnosis_results.md) | The original super-particle column's fractional-weight bookkeeping changes its results with the area it represents. An independent event solver under the same forcing closes its water and number ledgers, while its lunar sizes fail to converge with volume: their mass medians grow from 12.3 to 18.5 to 23.5 cm with each fourfold increase in sampled volume. A conservative sectional reference fails its finite size-domain check in heavy lunar snow, and 48 full columns trace part of the volume dependence to vapour allocation in a sparse population. The 20 October sizes are superseded |
| SCALE-SDM and published Earth columns | [SCALE lunar column](scale_lunar_results.md) and its [methods](scale_lunar_methods.md), [size diagnosis](scale_size_results.md), [known limitations](scale_limitations.md), [aggregation controls](analysis.md#scale-sdm-aggregation-controls), [Earth deposition and settling](welss_column.md), [Earth aggregation rerun](welss_aggregation.md) | Under the original lunar forcing SCALE-SDM gives mass medians near 1.6 mm, against 0.74 mm in its matched Earth configuration, since lower gravity and the deeper cloud each raise aggregation; custom controls under the same water flux give 14–18 cm, so the size followed the particle laws. Its single-crystal growth at −14.4 °C reaches 27.66 mm in 30 minutes from the fitted thin crystal and 0.96 mm from a sphere of the same mass and span, against the laboratory's fitted 4.17 mm |
| The crystal model against Earth's laboratory crystals | [growth benchmarks](growth_benchmarks.md), [growth options](growth_options.md), [McSnow growth](mcsnow_growth.md), [Earth growth assessment](earth_growth_assessment.md) with its facet, diffusion and surface-transfer components, [benchmark reconstruction](benchmark_reproduction.md) and its [plan](benchmark_reproduction_plan.md), [drag trace](benchmark_drag_trace.md), [growth trace](benchmark_growth_trace.md), [shape diagnosis](mcsnow_shape_diagnosis.md), [tuning](mcsnow_tuning.md), [state growth](mcsnow_state_growth.md), [growth limits](mcsnow_growth_limits.md), [column density](mcsnow_column_density.md), [joint audit](mcsnow_joint_acceptance.md), [recent-science check](mcsnow_science_check.md) | Removing McSnow's extra oblate mean-ventilation term cuts the archived benchmark's ten-minute mass error from 49.10% to 9.21% and gives the improved component every lunar run uses; the two-temperature calibration (`ded734e`) brings the warm column within 9.73% and the cold dendrite within 4.75% at 30 minutes. Shape, branching and density stay tied to the unreported starting shape, crystals from elongated seeds run away, and the shape and density rules lack the humidity-dependent transition between compact and complex growth. The 1,320-point rule these were judged by gave way to the data-sized rule |
| The crystal component in the lunar column | [lunar column](mcsnow_lunar_growth.md), [lunar storm](mcsnow_weather.md) | Under the 3 mm/h allowance the component grows 3.14 mm crystals from 4 µm seeds in 54 hours, and 9.51 and 29.08 mm at a tenth and a hundredth of the incoming crystals; through the archived storm's humidity the median peak is 0.724 mm, and 253 of 354 paths later sublimate to a tiny residual. The size range of 24 October supersedes these as single-crystal results |
| Geometry, strength, motions and light | [habits](habits.md), [geometry literature](geometry_literature.md), [strength literature](strength_literature.md), [turbulence](turbulence.md), [optics](optics.md) | The source review places the mass law's coefficient on a mean image dimension and the cited area fits on an ellipse; the archived cloud motions are conditioned on snow and growth environments, with flake-scale loading scenarios beside them; the optics compute snow extinction under height-resolved lunar light on a 294.9 K column, apart from the 301 K snow column |
