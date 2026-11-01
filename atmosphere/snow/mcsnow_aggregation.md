# Snowflakes from measured crystals: an Earth check, then the lunar column

McSnow's aggregation, fed with crystals of the size Takahashi and co-authors measured and run with the sticking
efficiency of Connolly and co-authors (2012), reproduces the median size of Earth's snowflakes. At 0.3, 1, 3 and 9 mm/h
of snowfall its median-volume melted diameter lies within 16% of both Gunn and Marshall's (1958) and Sekhon and
Srivastava's (1970) observed relations. In the heaviest snow its largest flakes outgrow the observed ones: at 9 mm/h the
90th-99th percentiles reach tens of centimetres where the CASP II storm produced 4–6 cm aggregates, because McSnow has
no law for whole aggregates breaking apart. On the Moon the same aggregation carries flakes past every span that could
plausibly survive within the first one to three kilometres of the 16 km snow zone. Lunar snowflake size is therefore set
by breakup, and plan step 5 estimates where breakup begins ([survival](turbulence_survival.md)).

## The calculation

McSnow (Welss and co-authors 2024; Brdar and Seifert 2018), at the pinned revision every snow calculation uses, follows
superparticles that settle and collide in a one-dimensional column. Two superparticles in a cell meet at a rate set by
their encounter area, the square of the sum of their spans, times the difference of their fall speeds, with Böhm's
collision efficiency; a meeting sticks with the chosen sticking efficiency. The column is held at ice saturation, so
crystals neither grow nor sublimate and aggregation alone acts on them. Crystals enter in the top tenth of the column at
a steady rate set by the snowfall, with masses drawn about their mean with a 30% spread, close to the spread between
repeat laboratory crystals ([data-sized check](data_sized_check.md)).

| | Earth check | Lunar column |
|---|---|---|
| Column | −15 °C to the freezing level at 6.2 K/km, 2.42 km, Earth gravity, 101 kPa at freezing | The 301 K ladder from −15 °C, 16.28 km at 0.92 K/km, lunar gravity, 65.9 kPa at freezing |
| Crystals | 28 µg, 2.75 mm across: Takahashi's −14.4 °C history at 20 minutes, with the measured crystals' mass–span law | The frozen component's crystals arriving at the lunar freezing level (plan step 3), with their mass–span and area laws |
| Snowfall | 0.3, 1, 3 and 9 mm/h | 3 mm/h at the reference crystal flux, and the cases below |
| Aggregates | McSnow's mixed-aggregate geometry (Mitchell 1996) | The same |
| Sticking | McSnow's four temperature laws | Connolly and co-authors (2012), with the others as sensitivities |

The crystals' area law is McSnow's own sector-branched crystal law (Mitchell 1996), and their mass law is fitted to the
21 measured crystals within 0.5 °C of −14.4 °C (m = 0.0124 D^2.21 in SI units). Three changes to the pinned source make
these inputs possible, each recorded as a patch: the monomer mass and area laws read from the namelist, a fall-speed
multiplier that carries the crystal model's fall-speed band into settling and collisions together, and the lunar
gravity, freezing-level pressure and lapse rate. With default values the patched Earth build reproduces the unpatched
rebuild of the authors' aggregation case byte for byte in all 25 output files over an hour.

## Earth: medians that match, tails that run long

Median-volume melted diameter of the snow arriving in the lowest 300 m, against the observed relations at the same
snowfall:

| Snowfall | Gunn and Marshall | Sekhon and Srivastava | Connolly (option 2) | Pruppacher and Klett (1) | Half Connolly (3) | Karrer (4) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.3 mm/h | 0.81 mm | 0.93 mm | 0.92 mm | 1.03 mm | 0.81 mm | 1.15 mm |
| 1 mm/h | 1.44 mm | 1.60 mm | 1.39 mm | 1.74 mm | 1.07 mm | 2.25 mm |
| 3 mm/h | 2.44 mm | 2.63 mm | 2.74 mm | 2.86 mm | 1.72 mm | 3.52 mm |
| 9 mm/h | 4.13 mm | 4.31 mm | 4.78 mm | 11.2 mm | 3.06 mm | 89 mm |

Connolly's sticking law stays within 0.87–1.16 of both relations at every snowfall; Gunn and Marshall fitted their
relation above 1 mm, so its 0.3 mm/h value is an extrapolation. Two further seeds at 1 mm/h give 1.32 and 1.37 mm. The
flakes' maximum spans with that law have mass-weighted medians of 4.0, 7.2, 19 and 42 mm across the four snowfalls, with
99th percentiles of 9, 17, 43 and 436 mm. Up to 3 mm/h these are the familiar sizes of Earth snow. At 9 mm/h the tail
reaches 44 cm at its 99th percentile and 60 cm at its largest flake, ten times the 4–6 cm aggregates observed in the
CASP II storm, and the two stickier laws run away altogether. The observed relations average many storms and crystal
types, while each run uses one crystal and one column, so the agreement is a consistency check of the aggregation, and
the heavy-snow tail marks where breakup, absent from McSnow, takes over on Earth.

## The Moon: aggregation outruns survival

Without a breakup law, the lunar runs grow aggregates that no snowflake could survive, metres across by mid-column, so
their final sizes say little. The telling quantity is how far below the cloud top the flakes first pass a span that
could plausibly survive. Plan step 5 places that span near 10–15 cm, the lunar equivalent of Earth's 4–6 cm storm
aggregates. The table gives the depth below the top at which the mass-weighted 90th percentile of span first exceeds
10 cm and 14 cm, read from 1 km layers of the final snapshot after 72 hours; the column is 16.28 km deep.

| Lunar case | 90th percentile passes 10 cm | passes 14 cm |
|---|---:|---:|
| 3 mm/h, reference crystal flux | 1.5 km | 1.5 km |
| 3 mm/h, a tenth of the crystals | 1.5 km | 1.5 km |
| 3 mm/h, a hundredth of the crystals | 0.5 km | 1.5 km |
| 1 mm/h, reference flux | 3.5 km | 4.5 km |
| 0.3 mm/h, reference flux | 7.5 km | 12.5 km |
| 0.3 mm/h, a tenth of the crystals | 10.5 km | not passed; median 7–8 cm, 90th percentile 11–12 cm near the freezing level |
| 3 mm/h, reference flux, fall speed ×0.555 / ×1.159 | 0.5 / 1.5 km | 0.5 / 2.5 km |
| 3 mm/h, reference flux, crystal span ×0.85 / ×1.18 | 1.5 / 1.5 km | 2.5 / 1.5 km |
| 3 mm/h, a hundredth, fall speed ×0.555 / ×1.159 | 0.5 / 1.5 km | 1.5 / 1.5 km |
| 3 mm/h, a hundredth, crystal span ×0.85 / ×1.18 | 0.5 / 1.5 km | 0.5 / 1.5 km |
| 3 mm/h, reference flux, sticking laws 1, 3 and 4 | 1.5, 2.5, 1.5 km | 1.5, 2.5, 1.5 km |
| 3 mm/h, reference flux, crystals entering through the whole column | 4.5 km | 4.5 km |

In every case the flakes reach the edge of the survivable range before the freezing level, and snowfall decides how
soon: within a kilometre or two of the top at 3 mm/h, within four or five at 1 mm/h, and only in the lower half of the
column at 0.3 mm/h. In the lightest case, with a tenth of the crystals, aggregation alone delivers flakes with a median
of 7–8 cm and a 90th percentile of 11–12 cm to the freezing level. The crystal model's fall-speed band, its largest
measured error, moves the 3 mm/h depths by about a kilometre: flakes that fall more slowly stay longer and crowd closer,
and pass the range sooner. The band's crystal-size edges and the other sticking laws move them by a kilometre, and
letting crystals enter through the whole column by three. Snowfall moves the 14 cm depth from 1.5 km at 3 mm/h to 12.5
km at 0.3 mm/h.

These depths hold under refinement where they are measured. At 3 mm/h, one, five and twenty-five times as many
superparticles give the 10 cm depth as 1.5 km each time and the 14 cm depth as 1.5, 1.5 and 2.5 km. Lower down the
column is unconverged: the 90th percentile 2–3 km below the top falls from 3.7 m to 2.2 m to 0.27 m as superparticles
are added, and metre-scale aggregates still form further down at every resolution.

Below those depths the runs keep aggregating with nothing to stop them; real flakes would break. If breakup holds
flakes near the span at which they fail, the flakes reaching the freezing level are those that survive, which plan
step 5 places at several centimetres to about 15 cm. That inference joins this numerical result to the load analogy,
and a breakup law would replace it.

## Limits

Crystals enter at their full size near the top, which lengthens their aggregation path; on the Moon they reach that size
only at the freezing level (plan step 3), so the runs bound aggregation from above, and the uniform-injection case moves
their entry through the column. Vapour growth, latent heating and depletion are switched off. McSnow has no
whole-aggregate breakup, so neither world's largest flakes are limited. The lunar column's lower part is unconverged in
superparticle number: more superparticles push the metre-scale flakes deeper, while the spans reached in the top three
kilometres, which carry the conclusion, hold at every resolution tested.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation build
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation parity
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation earth
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation moon --part 0
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation moon --part 1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_aggregation assemble
```

The [runner](mcsnow_aggregation.py) builds McSnow from the hash-checked source copy on the project-space drive and
writes the compact [product](results/mcsnow_aggregation.json); builds, namelists, raw output and logs are in
`research/runs/snow_agent/mcsnow_agg/`.
