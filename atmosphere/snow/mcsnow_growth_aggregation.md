# Lunar snowflakes when crystals grow and aggregate together

Grown from 4 µm ice spheres while they settle and aggregate, lunar snowflakes reach the load analogy's survivable span
of 10–14 cm far lower in the snow zone than the full-size crystals of 24 October did. In the reference column at 3 mm/h
with the reference crystal flux, the 90th percentile of span passes 10 cm about 8 km below the −15 °C cloud top and 14
cm about 10 km down; full-size crystals passed both within 1.5 km. With a hundredth of the crystals at 3 mm/h it passes
10 cm only 14.5 km down and reaches the freezing level at 10.5–10.9 cm, converged; a single run at 1 mm/h gives the same
depth and 12 cm. Light snow of 0.3 mm/h never reaches 10 cm: with the reference flux it arrives with a 90th percentile
of 4.2–4.5 cm, converged, and a single run with a tenth of the crystals gives 4.8 cm. Formed near −30 °C instead and
falling 32.6 km, the same light snow arrives with a median of 5.2–5.5 cm and a 90th percentile of 9.7–12.7 cm, from two
superparticle numbers that agree on the median within 4%. Nothing breaks the flakes in these runs, so where they
converge their sizes are the uncapped bound that stands beside the analogy's cap of about 8–15 cm. The two meet where
crystals are few; light snow stays below the cap without breaking when it forms at −15 °C and reaches it when it forms
near −30 °C. In heavy snow with plentiful
crystals, and under continuous saturation at the reference flux, the lower column grows metre-scale aggregates whose
size falls by half or more when the superparticles are quadrupled; those sizes are unconverged and carry no conclusion.

The flakes also limit the water they take. Each cell's supply is shared among the ice it holds, and when crystals gather
into a few large flakes the population presents far less surface to the vapour. At 3 mm/h the snow arriving at the
freezing level carries 1.2 mm/h with the reference flux and 0.6 mm/h with a hundredth of it, 40% and 20% of the supply;
at 1 and 0.3 mm/h it carries 73–90%. The single crystals of plan step 3 took 95–100% of their share. In a cloud, vapour
left behind by the ice would condense as supercooled water, which these runs, without riming, leave out.

## The calculation

The 24 October aggregation runs fed the lunar column with crystals already grown to the size they reach at the
freezing level, so aggregation worked on full-size crystals from the cloud top down
([aggregation](mcsnow_aggregation.md)). Here McSnow grows the crystals while they settle and collide. Crystals enter
the top tenth of the column as 4 µm ice spheres, grow from the vapour by the frozen improved component (McSnow's own
deposition, habit and fall-speed routines with the extra oblate mean-ventilation term omitted), and aggregate with
the same kernel, Böhm collision efficiency, Connolly and co-authors' (2012) sticking and Mitchell's (1996) aggregate
geometry as the 24 October runs. Nothing breaks the flakes: no breakup law with Earth support exists
([breakup literature](breakup_literature.md)).

Two settings carry the storm's water supply.

**The reference lunar column** is the 301 K ladder from −15 °C down to the freezing level, 16.28 km at 0.92 K/km
under lunar gravity, held at water saturation. In the water-limited cases a supply of 3, 1 or 0.3 mm/h is spread
evenly over the column's depth, and in each 20 m cell the ice present shares that cell's supply in proportion to the
growth each particle requests. For a column of identical crystals this gives each crystal its share of the supply for
every metre it falls, the rule of the single-crystal column program of plan step 3 ([size range](lunar_size_range.md));
with aggregation, flakes and the crystals they sweep up compete for the same supply. In the saturated cases nothing
limits the water. The crystal-number cases are the incoming fluxes of the earlier runs: the reference flux of 10,232
crystals per square metre per second, a tenth and a hundredth of it.

**The storm paths** are the 354 paths through the archived A28_dim5 storm that the error band followed with single
probe crystals ([Earth error band](earth_error_band.md), [lunar storm](mcsnow_weather.md)). Each path carries its own
column of 500 m cells from 46 km down to 24 km, below every freezing level in the archive, whose air follows the
storm at the probe's horizontal position every ten minutes, and after the probe ends moves with the archived wind at
its last height. A puff of 4 µm spheres fills the probe's release cell once; the crystals then grow, aggregate and
settle within their path's column until they reach air at 0 °C, where they are counted as arriving, or the archive
ends. The columns see the archived vapour, or the same air held at water saturation. Their crystal-number cases are
concentrations in the release cell: the mean concentration of the reference column at 3 mm/h, 123 crystals a litre,
a tenth and a hundredth of it. Crystals that sublimate to the mass of a 0.2 µm sphere are removed, and single crystals
that grow to a metre, the probes' own censor, are removed and recorded, as the probes were.

## Checks of the build

Changes to the pinned McSnow source join the three patches of the aggregation runs, each off by default: the improved
component's ventilation change with the error band's capacitance and axis multipliers, deposition substeps limited to
a 1% mass change, the cell water budget, forced cells that settle within sub-columns, removal at a sublimation
residual and a size censor, and binary snapshots of the superparticles. With every new option at its default each
build reproduces the 24 October runs byte for byte: an Earth aggregation run in all 19 of its output files, and a lunar
run in its final superparticle state with its profile output trimmed. Each later build also reproduces this module's
first column in its final state and in every snapshot field except one: the fall speed of the two seeds created in a
snapshot's step, which is unset there and carries no weight in the analysis.

With collisions switched off, the water-limited column at 3 mm/h and the reference flux delivers crystals of 3.18 mm
(10th–90th percentile 2.95–3.40 mm) and 0.084 mg to the freezing level, against 3.14 mm and 0.081 mg from the
single-crystal column program. The improved component inside McSnow's population model and the cell budget therefore
reproduce the single-crystal results of plan step 3 within 1.3% in span and 2.8% in mass.

## The reference column

Sizes are maximum dimensions weighted by mass flux, pooled over the last 12 hours of each run, when the column is
steady: its statistics at the final time and 12 hours before agree within the sampling scatter. Depths are below the
column top; the column is 16.28 km deep, and 32.56 km for light snow from −30 °C.

| Case | Snow at the freezing level | 90th percentile passes 10 cm | passes 14 cm | Arriving: median | 90th | 99th | largest | Superparticle check |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 3 mm/h, reference flux | 1.2 mm/h | 7.5–8.5 km | 9.5–10.5 km | unconverged | | | | ×4: agrees to 9.5 km |
| 3 mm/h, a tenth of the crystals | 1.0 mm/h | 8.5 km | 9.5 km | unconverged | | | | not run |
| 3 mm/h, a hundredth of the crystals | 0.6 mm/h | 14.5 km | not passed | 62–63 mm | 105–109 mm | 141–192 mm | 24–29 cm | ×4: agrees throughout |
| 1 mm/h, reference flux | 0.73 mm/h | 14.5 km | not passed | 65 mm | 119 mm | 190 mm | 31 cm | not run |
| 0.3 mm/h, reference flux | 0.23–0.27 mm/h | not passed | not passed | 22–24 mm | 42–45 mm | 58–66 mm | 80–85 mm | ×4: agrees throughout |
| 0.3 mm/h, a tenth of the crystals | 0.25 mm/h | not passed | not passed | 23 mm | 48 mm | 70 mm | 10 cm | not run |
| 0.3 mm/h from −30 °C, reference flux | 0.23–0.26 mm/h | 23.5–25.5 km | 23.5–25.5 km | 52–55 mm | 97–127 mm | 16–25 cm | 30–34 cm | ×½: median within 4%, 90th percentile agrees to 10.5 km |
| Saturated, reference flux | 1.6 mm/h | 2.5–3.5 km | 2.5–3.5 km | unconverged | | | | not run; 30–90 superparticles a kilometre below 13 km |
| Saturated, a hundredth of the crystals | 0.6 mm/h | 13.5 km | not passed | 61 mm | 103 mm | 139 mm | 17 cm | not run |

The ranges give the base run and the run with four times the superparticles, and for light snow from −30 °C the base run
and the run with half of them. Sizes in the rows without a check are single runs.

The crystal model's Earth error band, at its two combined edges, moves the 3 mm/h reference-flux depths: the edge
with 16% less capacitance, 18% less axis growth and 44% slower fall passes 10 cm at 4.5 km and 14 cm at 6.5 km, and the
opposite edge at 8.5 km for both, against 8.5 and 9.5 km at the band's centre in the same pooled reading. The slower
edge keeps its crystals longer in the cloud, so they meet more often and gather sooner, and it delivers 0.7 mm/h of
snow against 1.2 mm/h.

Growth and aggregation together move every depth of outgrowth down the column, because crystals that start as seeds
spend the first kilometres growing before they meet often enough to aggregate. In the 24 October runs the 90th
percentile passed both 10 and 14 cm within 1.5 km of the top at 3 mm/h with all three crystal numbers, within 4.5 km at
1 mm/h, and at 7.5 and 12.5 km at 0.3 mm/h. Here the 90th percentile at 3 mm/h and the reference flux grows from 4 mm
half a kilometre below the top to 25 mm at 2.5 km and 95 mm at 4.5 km, and holds near 10 cm through the middle of the
column before passing 14 cm. At 1 mm/h it reaches 25 mm at 4.5 km, 62 mm at 8.5 km and 115 mm at 14.5 km; at 0.3 mm/h, 7
mm at 4.5 km, 13 mm at 8.5 km and 34 mm at 14.5 km. With a hundredth of the crystals at 3 mm/h each crystal receives a
hundred times the water, and the column holds fewer, larger flakes that meet less often and reach 10 cm only 14.5 km
down.

## Convergence

Three cases were rerun with four times the superparticles at the same sampled volume, and light snow from −30 °C with
half of them; the others are single runs, marked so in the table above. Agreement of the 90th percentile by depth ends
at the first two successive kilometre bins that differ by more than 25%.

| Case | Superparticles in the column | 90th percentile agrees down to | Arriving 90th percentile, base and refined |
|---|---:|---:|---:|
| 3 mm/h, reference flux, ×4 | 5,256 and 20,861 | 9.5 km | 3.0 m and 1.2 m |
| 3 mm/h, a hundredth of the crystals, ×4 | 5,170 and 20,819 | the freezing level | 109 mm and 105 mm |
| 0.3 mm/h, reference flux, ×4 | 12,963 and 52,291 | the freezing level | 42 mm and 45 mm |
| 0.3 mm/h from −30 °C, ×½ | 30,299 and 15,188 | 10.5 km | 127 mm and 97 mm |

At 3 mm/h and the reference flux the onset depths hold, 7.5–8.5 km for 10 cm and 9.5–10.5 km for 14 cm, while below
9.5 km the metre-scale tail shrinks with more superparticles: the arriving 90th percentile falls from 3.0 to 1.2 m and
the median from 30 to 14 cm. The same tail appeared in the 24 October runs and shrank under refinement there. With a
hundredth of the crystals the whole column agrees within 18% at every kilometre, and in light snow every kilometre but
one agrees within 29%, with the arriving median at 22 and 24 mm, the 99th percentile at 58 and 66 mm and the largest
flakes at 80 and 85 mm. Light snow from −30 °C, run at half its superparticles, agrees down to 10.5 km and keeps its
arriving median within 4%, at 52 and 55 mm. Below 10.5 km its kilometre bins scatter in both runs, and halving the
superparticles moves the arriving 90th percentile from 127 to 97 mm and the 99th from 25 to 16 cm, so a finer run would
settle those two.

## Along the storm paths

| Humidity, release concentration | Paths delivering snow to 0 °C | Their arriving 90th percentile across paths: 10th, 50th, 90th percentile | Share over 10 cm | Parcel ice at its peak, median path |
|---|---:|---:|---:|---:|
| Archived, 1.2 a litre | 28 of 354 | 0.9 mm, 32 mm, 0.10 m | 11% | 0.005 kg/m² |
| Archived, 12 a litre | 63 | 15 mm, 0.12 m, 22 m | 56% | 0.05 kg/m² |
| Archived, 123 a litre | 77 | 0.11 m, 19 m, 38 m | 92% | 0.55 kg/m² |
| Saturated, 1.2 a litre | 299 | 0.6 mm, 0.11 m, 6 m | 53% | 0.86 kg/m² |
| Saturated, 12 a litre | 296 | 3 mm, 10 m, 26 m | 70% | 1.6 kg/m² |
| Saturated, 123 a litre | 276 | 0.07 m, 25 m, 34 m | 89% | 3.0 kg/m² |

In the archived humidity most crystals sublimate before they reach the freezing level, as the single probes did: with
1.2 crystals a litre only 28 of the 354 paths bring snow down, and the median of their flakes' 90th percentile is 32 mm.
The parcels' ice then stays well inside the storm's own: at its peak the median path holds 0.005 kg/m², and the
90th-percentile path 0.26 kg/m², against a median 0.36 kg/m² of ice and snow that the storm itself carries in the same
column, averaged over the archive. With ten and a hundred times more crystals, more paths deliver snow and their flakes
run away, and at 123 crystals a litre the parcels' ice exceeds the storm's. Held at water saturation, the paths take up
more than the storm's own ice at every concentration, 2.4 times it at 1.2 crystals a litre.

Four times the superparticles settle only one of these figures. Held at saturation with 1.2 crystals a litre, the median
path's arriving 90th percentile moves from 105 to 119 mm and the share of paths over 10 cm from 53% to 54%, while the
tail across paths, flakes of tens of metres, persists. In the archived humidity the few paths that deliver snow do not
converge: with 1.2 crystals a litre the median moves from 32 to 45 mm and the share over 10 cm from 11% to 26%, and with
12 crystals a litre the median moves from 0.12 m to 2.8 m. The storm paths therefore give one result: in this storm's
air held at saturation, a dilute puff of crystals grows into flakes whose 90th percentile reaches the freezing level at
about 11–12 cm on the median path, while the 90th-percentile path runs away to metres. The runaway flakes form where a
puff of crystals released together in strongly supersaturated air grows and gathers with nothing to break it or to
deplete its vapour.

## Limits

Nothing breaks the flakes, so each size is the size reached under aggregation alone; the McSnow aggregation that
reproduces the median of Earth's snow also overgrows its largest flakes at 9 mm/h tenfold, 44 cm against the 4–6 cm
observed ([aggregation](mcsnow_aggregation.md)). Riming is off, so vapour that the ice does not take stays unused. The
reference column is a straight profile without updraughts, and its supply is spread evenly with height without
responding to the snow. The storm paths' air does not respond to their crystals, and each path's crystals start as one
puff in one cell. The crystal model is the frozen improved component, whose Earth errors move single-crystal sizes by a
factor of 1.4–1.6 where water limits growth ([Earth error band](earth_error_band.md)).

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation build
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation parity
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation growth_parity
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation column --case q3_f1 [--variant sp4] [--edge combined_small]
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation paths --case saturated_c0.01 [--variant sp4]
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_growth_aggregation assemble
```

The [runner](mcsnow_growth_aggregation.py) builds McSnow from the hash-checked source copy, names each build by the
hash of its patched sources, and writes the compact [product](results/mcsnow_growth_aggregation.json). Builds,
namelists, forcing tables, snapshots and logs are in `research/runs/snow_agent_2/growagg/`, on the project-space drive.
