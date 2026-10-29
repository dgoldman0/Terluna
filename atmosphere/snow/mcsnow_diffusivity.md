# McSnow's vapour diffusivity in the lunar builds

The lunar McSnow builds of 24 and 25 October diffused vapour at 0.651 of McSnow's own fit, and restoring the fit's
reference makes their crystals take up more of the vapour and their flakes larger. McSnow takes Pruppacher and Klett's
fit, Dv = 0.211 cm²/s × (T / 273.15 K)^1.94 × ps / p, in which ps, McSnow's constant surface pressure, stands for the
fit's 1 atm reference. To anchor the lunar column at the freezing level, the patch of
[mcsnow_aggregation.py](mcsnow_aggregation.py) set ps to that level's pressure, 65,925 Pa, and so set the reference as
well: every lunar run built from it diffused vapour at 65,925 / 101,325 = 0.651 of the fit, at every height. The fit
itself stands 3.4% below Massman's (1998) measured value at −5 °C and 4.6% below it at −30 °C, the frozen component's
stated uncertainty.

Rerun with the fit's reference restored, light snow barely moves, whether it forms at −15 °C or near −30 °C, while every
case in which the crystals take up less than the supply, or grow without a supply limit, grows faster and larger:

| Case, base run | Snow at the freezing level | 90th percentile passes 10 cm; 14 cm | Arriving: median | 90th | 99th | largest |
|---|---:|---:|---:|---:|---:|---:|
| 0.3 mm/h, reference flux | 0.23 → 0.23 mm/h | not passed → not passed | 22 → 21 mm | 42 → 41 mm | 58 → 59 mm | 80 → 77 mm |
| 3 mm/h, a hundredth of the crystals | 0.58 → 0.79 mm/h | 14.5 → 11.5 km; not passed | 62 → 69 mm | 109 → 126 mm | 192 → 197 mm | 29 → 43 cm |
| 3 mm/h, reference flux | 1.25 → 1.56 mm/h | 8.5 → 6.5 km; 9.5 → 8.5 km | unconverged | | | |
| Held at saturation, a hundredth of the crystals | 0.59 → 0.79 mm/h | 13.5 → 11.5 km; not passed | 61 → 71 mm | 103 → 130 mm | 139 → 250 mm | 17 → 50 cm |
| 0.3 mm/h from −30 °C, reference flux | 0.26 → 0.26 mm/h | 25.5 → 23.5 km; 25.5 → 32.5 km | 55 → 53 mm | 127 → 124 mm | 246 → 276 mm | 34 → 29 cm |

Sizes are by mass flux, pooled over the last 12 hours; depths are below the column top. In light snow at the reference
flux the crystals could take up more vapour than the supply gives them, so the supply sets their growth and the
correction moves the arriving median and 90th percentile by under 5%; from −30 °C the 99th percentile and the largest
flake move by 12–13%, inside the 27–41% that halving the superparticles moved them, and the kilometre bins that set the
14 cm depth scatter as before. At 3 mm/h the crystals took up only part of the supply, at the rate diffusion allowed,
and the corrected crystals take up 24–37% more of it: 52% of the supply at the reference flux and 26% with a hundredth
of the crystals, against about 40% and 20%. The flakes then reach 10 cm 2–3 km higher in the column and arrive 11–26%
larger at the median and the 90th percentile. Along the archived storm's paths at a tenth of the reference
concentration, 62 paths deliver snow against 63, sublimation still removes 16% of the ice grown, and 15% more ice grows,
so 51% of it arrives at 0 °C against 37%, with 33% leaving at the single-crystal censor against 46%.

Melted through the response tables of the [melting record](snow_melting.md)'s median columns, with each range spanning
land and sea, the corrected flakes of light snow from −15 °C finish within 0.05 km of before and send 10% of their water
to the surface against 11%, and those from −30 °C finish within 0.13 km and send 49–51% against 50–51%. With a hundredth
of the crystals at 3 mm/h, 99% of the ice has gone 3.83–3.95 km down against 3.64–3.75 km, and 61–63% of the water
reaches the surface against 56–57%; held at saturation, 4.0–4.2 km against 3.4–3.5 km, and 63–65% against 55–56%. In the
heavier snow sublimation's share falls by 0.8–1.3 points over the sea and 0.3–0.5 over land, and in light snow it moves
by 0.6 points or less.

These are the single base runs compared case by case. On 29 October the rest were repeated on the corrected build:
the other reference-column cases (3 mm/h at a tenth of the crystals, 1 mm/h, 0.3 mm/h at a tenth, saturation at the
reference flux), the superparticle refinements (×4 at 3 mm/h with the reference flux and with a hundredth of it and at
0.3 mm/h, ×½ for the deep column), the error band's two edges at 3 mm/h and the eight other storm-path runs, eighteen
runs in all, every one returning cleanly, with the growth parity and the parity and reproduction checks after them. The
products of [growth and aggregation](mcsnow_growth_aggregation.md) and the [melting record](snow_melting.md) now read
the corrected runs, and their records quote them. The lunar columns of the 24 October aggregation runs, which the
growth runs superseded, are not repeated; they stand as the old build's record, and the corrected build no longer
reproduces their lunar parity run, where the Earth build still reproduces its own in all 19 files.

## The correction

[mcsnow_diffusivity.py](mcsnow_diffusivity.py) gives the fit its own constant, p_dv = 101,325 Pa, in the lunar build
alone, in McSnow's `diffusivity_vapour` and in the same line of its condensation routine. The corrected lunar sources
differ from the old lunar build in those two lines, their two use lists and the new constant, and in nothing else; the
Earth build is the old one, executable `087804b7…`, so the calibration and the data-sized check are untouched. The
corrected runs repeat the old runs' recorded configurations, seeds and forcing on the corrected build
(`research/runs/snow_agent_3/dv_fix`), and the comparison reads both with the same statistics: sizes by mass flux
pooled over each column's last 12 hours, the depths below the top where the 90th percentile passes 10 and 14 cm, and
for the storm path its arrivals and the column budget of the [melting record](snow_melting.md).

## Which runs carry the factor

Every lunar run built from `mcsnow_aggregation.patched_sources('moon')` carries it: the lunar columns of the aggregation
runs of 24 October ([aggregation](mcsnow_aggregation.md)) and every lunar column, storm path, parity and reproduction
run of growth with aggregation on 25 October ([growth and aggregation](mcsnow_growth_aggregation.md)), and through them
the flakes of the [melting record](snow_melting.md); since 29 October those two records read the corrected runs, and
only the superseded aggregation columns of 24 October carry the factor. The single-crystal lunar runs of 23–24 October
patch only gravity and keep ps at 101,325 Pa, so their diffusivity is the fit's: the lunar column and storm of the
improved component ([lunar column](mcsnow_lunar_growth.md), [lunar storm](mcsnow_weather.md)), the [Earth error
band](earth_error_band.md) and the [size range](lunar_size_range.md). The Earth builds keep the 1 atm reference.

## Reproduction

`python -B -m atmosphere.snow.mcsnow_diffusivity build` compiles the corrected lunar build and records its source
differences, `run --case <name>` repeats one recorded run on it, `compare` reads the old and corrected runs and
writes [results/mcsnow_diffusivity.json](results/mcsnow_diffusivity.json), and `assemble` writes the growth-and-
aggregation and melting products from the corrected runs, repeating the growth parity and the parity checks first
(the Earth build stays the old cache's). The eighteen runs of 29 October went four at a time through
`session-tools/snow_reruns/run_reruns.py` on the research drive, about 3.9 hours of one core in all.
`atmosphere/tests/test_snow_mcsnow_diffusivity.py` checks the source differences, the old factor and this record's table
against that product.
