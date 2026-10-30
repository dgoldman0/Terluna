# McSnow's lunar column with the air below a cloud base held over ice

PlaSim's snow sublimates below a cloud base by gamma, an Earth-calibrated rate measured against water saturation, and
replacing it with the lunar flakes' own sublimation needs a reference in which those flakes fall out of a cloud into the
air the design boxes hold below it. McSnow's idealised profile cannot give one as it stands. It holds its liquid-water
region between the heights h1 and h2 at water saturation and dries the air below h1 at a fixed lapse of relative
humidity, 0.05% per metre, so the snow work's reference columns set h1 to the freezing level, and their snow never meets
subsaturated air above it. [mcsnow_subcloud.py](mcsnow_subcloud.py) adds one declared change to the corrected lunar
build ([diffusivity](mcsnow_diffusivity.md)), off by default: with `atmo_nml`'s `ssat_below` at −1 or above, the air
below h1 holds that supersaturation over ice. The cloud above h1 is then the reference columns' own, and the snow that
leaves it sublimates in air as dry over ice as the boxes'. Against these columns, a closure that takes the share of the
snow's flux sublimated per metre from three reference columns by depth below the snow's top, matches them by the snow's
flux at the cloud base and raises the rate as the flakes shrink follows McSnow within 0–9% over both boxes, in snow of
1 and 0.3 mm/h; without the last two terms it loses a quarter to a third of McSnow's rate in the light snow. Carried
to the boxes' own range, with references from −30 °C as well, it holds in light snow from −30 °C and at the dry end of
the boxes' sub-cloud air, overstates the loss where nearly all of a light snow from −15 °C sublimates, and meets
McSnow's own limit in snow of 1 mm/h from −30 °C, whose flakes grow metres across without a breakup law.

## The change

The change touches thirteen lines of two sources. In `mo_atmo_nml.f90` it declares `ssat_below`, adds it to the
namelist, sets its default to −999 and prints it when set. In `mo_atmo.f90`, where `init_atmo_ideal` fills each height,
the air at or below h1 takes a vapour pressure of (`ssat_below` + 1) times saturation over ice, with its
supersaturation and relative humidity to match, and the surface vapour pressure, from which McSnow builds the column's
dry pressure, stays at water saturation, as with the cloud reaching the ground. The column's pressure therefore stays
the reference columns'. Off, the changed build gives the corrected build's `col_sat_f0.01_base` back byte for byte in
its snapshots and every kept table (29 October); the Earth build is the old cache's and stays untouched.
`atmosphere/tests/test_snow_mcsnow_subcloud.py` checks the thirteen lines, the default, the Earth sources and each
case's configuration against the run it repeats.

## The boxes' air

Each box's cloud base is the lowest height above its freezing level where its mean humidity over ice reaches 0.98, and
its sub-cloud air the mean humidity over ice below that, weighted by the air's mass, over days 29.5–59.0 of the design
boxes (`session-tools/snow_sublimation/subcloud_air.py` on the research drive, by the rule written with it). The sea box
(`box_near_side_sea`) has its cloud base 3,957 m above its freezing level at 23.11 km and 0.9469 below it, the land box
(`box_0e_design`) 8,254 m above 25.61 km and 0.9604. The cases hold the air below h1 at −0.0531 and −0.0396 over ice.

## The cases

Each case repeats a reference column's recorded configuration and seed, with h1 and `ssat_below` set from a box. McSnow
gives each cell the water supply over the column's depth, per metre, so a column's own supply gives each cloud cell its
share.

| Case | Column repeated | Box | Supply | Snapshots |
|---|---|---|---|---|
| `sub_sea`, `sub_land` | `col_q1_f1_base`, 1 mm/h, 96 h | sea, land | raised by the column's depth over the cloud's | 6 h |
| `sub_sea_b`, `sub_land_b` | `col_q1_f1_base` | sea, land | its own | 6 h |
| `sub_sea_c`, `sub_land_c` | `col_q1_f1_base` | sea, land | its own | 30 min |
| `sub_sea_l`, `sub_land_l` | `col_q0.3_f1_base`, 0.3 mm/h, 132 h | sea, land | its own | 30 min |
| `ref_q0.3_f1`, `ref_q1_f1`, `ref_q3_f1` | the three references at the reference crystal flux, 132, 96 and 72 h | none | its own | 30 min |
| `ref_q0.1_f1` | 0.1 mm/h from −15 °C, 180 h, built as the snow work builds its columns | none | its own | 30 min |
| `rdeep_q0.1_f1`, `rdeep_q0.3_f1`, `rdeep_q1_f1`, `rdeep_q3_f1` | 0.1, 0.3, 1 and 3 mm/h from −30 °C, 228, 168, 120 and 96 h (`rdeep_q0.3_f1` repeats `col_deep_q0.3_f1_base`) | none | its own | 30 min |
| `sdeep_sea_q1`, `sdeep_land_q1`, `sdeep_sea_q0.1`, `sdeep_land_q0.1` | `rdeep_q1_f1`, `rdeep_q0.1_f1` | sea, land | its own | 30 min |
| `sub_sea_q0.1`, `sub_land_q0.1` | `ref_q0.1_f1` | sea, land | its own | 30 min |
| `sub_sea_dry`, `sdeep_land_dry` | `col_q1_f1_base`, `rdeep_q1_f1` | sea, land at the dry end (−0.0967, −0.1091) | its own | 30 min |
| `smid_sea_q1` | 1 mm/h from −22.5 °C, halfway in height between the two tops, 108 h | sea | its own | 30 min |

The first pair raised the supply on a misreading of McSnow's, which is already per metre, and gave each cloud cell 1.32
(sea) and 2.03 (land) times its share; it stands as the record of that reading. The snapshot writer works on copies of
the superparticles, so writing every 30 minutes leaves a run itself unchanged, which the second test checks byte for
byte at every sixth hour.

## The first test

The first test (rules set on 29 October at 20:27, before its runs) took PlaSim's closure as k by depth below the snow's
top from the moderate column, the share of the snow's flux a unit of ice subsaturation sublimates per metre by the snow
work's flake law, and asked whether the share of the cloud base's flux lost by the lowest 0.5 km layer fell within
0.80–1.25 of McSnow's. On the first pair the cloud carried 1.77 and 1.42 times the moderate column's flux at its base,
and the closure lost 0.886 (sea) and 1.776 (land) of McSnow's share. On the corrected pair the cloud carried 1.142 and
0.909 times it, and the closure lost 0.725 and 0.869 of McSnow's share, 0.269 against 0.371 over the sea and 0.426
against 0.490 over land. Nothing was adopted.

The readings added after that test's rules, which change none of its verdicts, find two things behind it. McSnow's
share made the rules' way, from three snapshots six hours apart and two layers, scatters by 0.11–0.37 of itself between
consecutive 18-hour windows of one run, more than the band's width on two of the four runs, and the snapshot fluxes of
the moderate column's lowest layer range over 0.37–0.89 mm/h. Within that scatter, the closure lost less than McSnow on
the corrected pair by every reading tried, 0.63–0.88. The closure also leaves out two pieces of physics. Flakes that
sublimate shrink, and a smaller flake loses a larger share of its mass per metre: by the flake law, a population keeping
a share r of its mass has k raised by r^−q, q 0.24–0.27, while the moderate column's flakes grow all the way down. And
lighter snow keeps smaller flakes: at the reference crystal flux, k stands 3.3 times higher 7 km below the top and 2.7
times higher 13 km below it at 0.3 mm/h than at 1 mm/h, a dependence a test at 1 mm/h alone cannot see.

## The second test

The second test (rules set on 29 October at 21:08, before its runs) gave the closure both missing terms and sampled
McSnow finely enough to judge them. Three references at the reference crystal flux, 0.3, 1 and 3 mm/h, give k and its
shrinking exponent q by depth, and a straight line through each one's flux below the case's cloud base. The snow's flux
at the cloud base sets weights between the references, linear in the logarithm of the flux, and with them k, q and the
line. Below the cloud base k is raised by [(F/F(h1)) / (F̂(z)/F̂(h1))]^−q, the flakes' mass against the matched
reference's at the same height, since the reference's flakes keep growing where the case's sublimate. The statistic is
D = ln(F_base/F_low), the e-foldings the flux loses from the 0.5 km layer above the cloud base to the lowest 0.5 km,
which measures the rate in light snow where a share saturates. McSnow's D is pooled over each run's last half, 97 or 133
snapshots 30 minutes apart, and its error taken from eight consecutive blocks. Each band combines the earlier tests' 25%
tolerance with twice that error.

Every six-hourly snapshot of the five reruns came out byte for byte as the original run's, so the finer snapshots leave
the runs unchanged. McSnow's relative errors came to 0.021–0.077, all four cases resolved, and the closure held on each:

| Box run | McSnow's D | Plain | Shrinking only | Rate only | Both | Band |
|---|---:|---:|---:|---:|---:|---:|
| Sea, 1 mm/h (`sub_sea_c`) | 0.373 | 0.876 | 0.926 | 0.859 | **0.908** | 0.762–1.312 |
| Land, 1 mm/h (`sub_land_c`) | 0.620 | 0.888 | 1.022 | 0.917 | **1.060** | 0.776–1.289 |
| Sea, 0.3 mm/h (`sub_sea_l`) | 0.912 | 0.359 | 0.379 | 0.857 | **1.001** | 0.792–1.263 |
| Land, 0.3 mm/h (`sub_land_l`) | 2.301 | 0.239 | 0.275 | 0.677 | **1.024** | 0.797–1.255 |

The columns give the closure's D over McSnow's; plain is the first test's closure. At the moderate rate the plain
closure holds once McSnow is sampled enough, so the first test's misses came from its sampling; in light snow it loses
a quarter to a third of McSnow's rate, the snow's rate carries most of the correction and the shrinking the rest. Layer
by layer below the cloud base the closure follows McSnow's flux within 0.04 of the base flux, and the flake law applied
to each run's own flakes gives 0.93–1.04 of McSnow's D, so the law follows McSnow's own sublimation. The clouds carry
their references' flux at the base, 0.96–1.02 of it. The flakes' flux-weighted size at the cloud base runs from 55 mm
(sea, 1 mm/h) to 8 mm (land, 0.3 mm/h), and the sublimation below takes 245–708 W/m² for each mm/h arriving at the base.
By the test's adoption rule PlaSim's snow then sublimates by this closure in place of gamma, built with tests and a
smoke day before the design run.

The references start at −15 °C, 16.28 km above the freezing level, and the closure is tested for snow from such a top.
The design boxes' snow tops stand a median 29–41 km above the freezing level, at −34 °C over the sea and −50 °C over
land ([snow zone](snow_zone_depth.md)), beyond the references and beyond McSnow's habit table, which ends at −30 °C.
Snow that forms higher grows longer and sublimates more slowly: near the freezing level k stands at 2.0 × 10⁻³ per metre
31 km below the top of the 0.3 mm/h column from −30 °C, against 3.6 × 10⁻³ 16 km below the top of the column from
−15 °C. How PlaSim takes the deeper tops is the build's open choice.

## The third test

The third test (rules set on 29 October at 22:47, before its runs) carried the closure to the boxes' own range: a second
family of references from −30 °C, the coldest top McSnow's habit table allows, at 0.1, 0.3, 1 and 3 mm/h; light snow of
0.1 mm/h from −15 °C; the dry end of the air the boxes' falling ice meets, the flux-weighted tenth percentile of its
humidity over ice (0.903 sea, 0.891 land); and a top halfway in height between the families, which tests the rule
between them, ln k, q and the flux line interpolated at the same height by the top's height. The statistic is D to the
lowest layer that keeps 5% of the base flux, and a case whose D changes between the two halves of its pool by more than
six of its errors is unresolved. The new runs' directories sit on the research drive; the deep rerun repeats
`col_deep_q0.3_f1_base` byte for byte at every sixth hour.

Five of the nine box runs resolved. Four did not, and the reason lies in McSnow: from −30 °C, and from halfway, snow of
1 mm/h falls far enough to grow metre-scale aggregates without a breakup law, flux-weighted 0.7–3.2 m across at the
cloud base and up to 4.6–7.6 m, with the largest 1% of a layer's superparticles carrying 25–79% of its flux. Such flakes
lose almost nothing below the cloud base, and a handful of superparticles set the flux's noise, so McSnow's D sits at
zero within errors of 0.08–0.17. The 3 mm/h column from −15 °C grows the same, 1.6 m across on arrival. The snow work's
load analogy puts the survivable span at 5 cm for ordinary flakes and 8–15 cm for the largest; the 0.1 mm/h snow from
−30 °C stays well inside it (18–22 mm at the cloud base), and the 0.3 mm/h snow from −30 °C reaches its upper edge.

On the five resolved runs the closure holds within 10% in light snow from −30 °C (1.055 sea, 1.100 land, where McSnow
loses 61% and 82% of the flux) and at the dry end (1.077), and overstates the loss by 1.48 and 1.65 in 0.1 mm/h snow
from −15 °C, where McSnow loses 94%. The flake law applied to each run's own flakes follows McSnow within 4% in all
five, so the miss lies in the shrinking term: it takes every flake as losing the same share of its mass, while where
most of the flux sublimates the smallest flakes vanish whole and leave the flux to larger ones. By its rule nothing was
adopted.

## The fourth test

The fourth test (rules set on 30 October at 02:45, offline on the nine resolved runs of the second and third tests)
replaced the shrinking term with the matched references' own flakes, carried down from the cloud base by the flake law,
each losing mass at its own rate and gone when it has none left. It mends the light snow from −15 °C, 1.07 and 1.13 of
McSnow's D against 1.48 and 1.65, and misses where the snow falls through the land's 8 km below the cloud base at 1 and
0.3 mm/h, 1.54 and 1.37 where the shrinking closure gave 1.06 and 1.02: the carried flakes stop aggregating at the cloud
base, while McSnow's keep merging into larger flakes that lose less per metre. Each closure leaves out one of the two
processes, the shrinking closure the vanishing of the smallest flakes and the carried flakes the aggregation. As shares
of the cloud base's flux sublimated, the shrinking closure stays within 0.03 of McSnow's except in the lightest snow
from −15 °C, where it sublimates 98–99% against 94%, and the carried flakes within 0.06 except at 1 mm/h over land,
0.62 against 0.46. The carried flakes were chosen after these runs showed the miss and fit no parameter; the choice
between the two is open.

## Reproduction

`OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud build` compiles the changed lunar build
(`research/runs/snow_agent_3/subcloud`, executable `b5c6f8c1…`), `parity` repeats `col_sat_f0.01_base` on it with
`ssat_below` off, and `run --case <name>` runs one case. The rules, scripts and outcomes of the tests are kept on the
research drive in `session-tools/snow_sublimation`: `rules_subcloud.txt` with `subcloud_test.py` and the added readings
(`subcloud_narrowness.py`, `subcloud_pooled9.py`, `subcloud_shrink.py`), and `rules_subcloud2.txt` with
`subcloud_test2.py` and `run_subcloud2.sh`, and `rules_subcloud3.txt` with `subcloud_test3.py`, `run_subcloud3.sh`,
`box_snow_rates.py`, `reference_convergence.py` and `subcloud3_noise.py`, and `rules_subcloud4.txt` with
`subcloud_test4.py`. `run --case <name> --smoke` runs a case's first 6 model hours, and `--store DIR` places a new run's
directory on another disk, linked from the cache.
