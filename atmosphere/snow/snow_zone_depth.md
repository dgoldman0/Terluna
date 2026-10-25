# How deep the snow zone stands, and whether light snow reaches the size limit

The design-climate storm boxes hold their snow far above the freezing level. Weighted by the snow each column holds,
half of it lies within 18–20 km of the freezing level and nine-tenths within 30–35 km, and the snow's top stands a
median 29–41 km above it, at −34 °C over the Near-side Sea and −50 °C over land. Light snow of 0.3 mm/h does not reach
the load analogy's survivable span of 10–14 cm before the freezing level. Through the reference column, which spans the
lowest 16.3 km of that snow zone from −15 °C, its flakes arrive with a 90th percentile of 4.2–4.5 cm, a 99th of 5.8–6.6
cm and the largest at 8–8.5 cm, converged under four times the superparticles. The light snow that forms higher, in the
boxes' snow between −15 °C and −30 °C, has not yet been followed to the freezing level. Below the freezing level the
flakes melt within 1–3 km while they fall, and these runs follow no aggregation there.

## The storm boxes

The two CM1 boxes on the design climate, the cold-cloud run at 298.5 K, are 192 km squares forced by the GCM's column at
their sites: over land at 0° N, 0° E, and over the Near-side Sea at 5° S, 45° W. Each saves its fields every three hours
for 59 days on 111 levels to 150 km, 2 km apart through the snow zone. For every saved column the calculation finds the
freezing level, where the air first falls through 0 °C going up, and then the snow above it: the highest level holding
at least 10⁻⁵ kg of snow per kilogram of air (10⁻⁶ as a sensitivity), the temperature there, and the heights above the
freezing level below which half and nine-tenths of the column's snow lie. Columns are counted alone and weighted by the
snow mass they hold above the freezing level. CM1's snow is the snow category of its Morrison double-moment
microphysics, the aggregates and large crystals of the boxes' own clouds; the calculation reads it as the boxes save it
and simulates nothing.

## The snow zone

| | Land, 0° N 0° E | Near-side Sea, 5° S 45° W |
|---|---:|---:|
| Columns holding snow above the freezing level | 45% of 484,352 | 35% of 484,352 |
| Freezing level, median (weighted by snow) | 25.7 km (26.1 km) | 23.4 km (23.2 km) |
| Height above the freezing level below which half the snow lies, median (weighted) | 16.2 km (20.2 km) | 17.4 km (18.3 km) |
| The same for nine-tenths of the snow | 33.1 km (34.9 km) | 32.4 km (29.6 km) |
| Top of the snow above the freezing level, median (weighted) | 41.2 km (51.2 km) | 29.1 km (36.8 km) |
| Temperature at the snow's top, median (weighted) | −50 °C (−64 °C) | −34 °C (−45 °C) |
| Snow above the freezing level in a snowy column, median | 0.71 kg/m² | 0.28 kg/m² |

The heights hold when the snow threshold drops tenfold: the snow-weighted half and nine-tenths heights move by less than
0.2 km, while the top of the thinnest snow rises 4–16 km. The heaviest tenth of the snowy columns keep their snow at the
same heights, half within 18–21 km and nine-tenths within 28–37 km. The snow zone of these storms therefore spans about
20 km above the freezing level for half its snow and 30–35 km for most of it. The reference column of the snow work,
16.3 km from 0 to −15 °C on the ladder, covers the lower half of that snow; the boxes' own −15 °C stands at 15–16 km.

## Light snow

In the reference column, 16.28 km from −15 °C down to the freezing level at water saturation, light snow of 0.3 mm/h
grown from 4 µm seeds while it aggregates never carries its 90th percentile of span past 10 cm. With the reference
crystal flux its flakes arrive with a median of 22–24 mm, a 90th percentile of 42–45 mm, a 99th of 58–66 mm and the
largest at 80–85 mm, the pairs from the base run and the run with four times the superparticles; a single run with a
tenth of the crystals gives 23, 48, 70 and 104 mm ([growth and aggregation](mcsnow_growth_aggregation.md)). The 90th
percentile grows steadily down the column, 7 mm at 4.5 km below the top, 13 mm at 8.5 km and 34 mm at 14.5 km, so the
flakes that reach the freezing level are still growing. The load analogy places breakup near 10–14 cm, two to three
times the 90th percentile and above the largest flakes, so light snow from the reference column melts before it reaches
the size limit. Below the freezing level the flakes melt within 1–3 km while they fall, and they would aggregate further
there as wet, sticky flakes, which these runs do not follow.

The boxes' snow reaches higher than the reference column's top: half of it lies above 18–20 km over the freezing level
and the −15 °C level stands 15–16 km up in their air. Light snow that forms near −30 °C would fall through about twice
the depth before it melts, through air where Connolly and co-authors' sticking efficiency is 0.09–0.14, against 0.64
at −15 °C. Whether those flakes reach the limit has not yet been computed.

## Limits

The boxes' snow is CM1's Morrison scheme's own, so its heights describe where these storms make and carry snow, while
the flake sizes come from McSnow's separate calculation in idealised columns ([growth and
aggregation](mcsnow_growth_aggregation.md)). The columns follow the ladder's lapse rate of 0.92 K/km, which puts −15 °C
16.3 km and −30 °C 32.6 km above the freezing level; in the boxes' mean air −15 °C stands 15–16 km and −30 °C 26.5–27.6
km above it, at 1.1 K/km. The crystal model's habit table ends at −30 °C, and no column starts colder. Melting is
followed only as the single-flake melting depths of the flake model (README, "Snow lasts far below the freezing level");
flakes that aggregate while they melt, wet and sticky, are outside these runs.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.snow_zone_depth
```

The [runner](snow_zone_depth.py) reads the boxes' archived fields from `climate/crm/runs` on the project-space drive,
four fields from each saved file, and writes the compact [product](results/snow_zone_depth.json) with one hash for
each box over every field it read.
