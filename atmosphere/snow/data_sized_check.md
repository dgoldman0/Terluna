# The frozen crystal model against each measured crystal

The two-temperature calibration passes the data-sized rule on mass and axes. Of its 93 comparisons with individual
crystals at −5.3 and −14.4 °C, 5 fall outside the scatter between repeat crystals, where chance alone would place about
4. Its fall speed passes at every −14.4 °C age and fails at −5.3 °C from 12 minutes on, where the hollow column falls
up to 6.1 cm/s faster than the published relation allows. The improved McSnow component, which runs every lunar
calculation, fails the rule: 83 of its 332 mass and axis comparisons fall outside, against about 15 by chance, and its
fall speed misses the ±3 cm/s accuracy at 55 of 78 ages, mostly by falling too fast.

## The rule and the scatter it uses

Takahashi and co-authors (1991) grew 169 crystals in a vertical wind tunnel at water saturation and weighed and
measured each once, at one age; Fuchs and co-authors archived the values in 2025
([measurements](results/takahashi_measurements.json)). The experiments ran at ten nominal temperatures, and 134 of the
crystals lie within 0.5 °C of one of them. Crystals of one nominal temperature and one age are repeats. At 34 such
pairs two to five repeats exist, and their spread measures how much crystals grown under the same conditions differ.
Pooling the spread over all pairs, with each pair's deviations taken from its own mean, gives one standard deviation
per quantity in natural logarithms:

| Quantity | Pooled scatter (ln) | Degrees of freedom | Two-sigma band, model/crystal |
|---|---:|---:|---|
| Mass | 0.193 | 54 | 0.68–1.47 |
| Equatorial axis | 0.159 | 54 | 0.73–1.37 |
| Polar axis | 0.183 | 23 | 0.69–1.44 |

A crystal passes on a quantity when the model at its age and nominal temperature lies within two of these standard
deviations; about 4.6% of repeat crystals lie outside that band by chance. Every crystal uses the pooled value, ages
with a single crystal included: a pair of repeats measures the spread poorly, and pooling assumes that the relative
spread is similar at every age and temperature. The scatter differs by temperature, from 0.04 to 0.27 in ln mass, on 2–11 degrees of freedom each. Using each
temperature's own scatter where it rests on at least five degrees of freedom changes 1 of the calibration's 93
verdicts and 28 of the component's 332. The authors inferred thin planar thickness from mass and basal area, so the
plates, sectors and dendrites carry no polar axis crystal by crystal.

The archive holds no fall speeds crystal by crystal. The paper's Table 4 gives fall speed against age at each
temperature, each line the average of about four sets of measurements, and states the measurement accuracy as
±3 cm/s. The model's speed is compared with that relation at every age where crystals exist, from 3 minutes to the end
of each relation's interval, and passes within ±3 cm/s. This test covers the average speed; the spread of speeds
between crystals is untested.

The model values come from the frozen executables without change. The calibration's runs reproduce its recorded
baseline exactly at both temperatures, and halving their timestep changes no value by more than 0.011%. The
component's values come from its retained traces at 1/128 s, which differ from the 1/32 s traces by at most 0.092%.

## The two-temperature calibration

| Temperature | Quantity | Crystals | Outside | Median model/crystal |
|---|---|---:|---:|---:|
| −5.3 °C | Mass | 17 | 1 | 0.982 |
| −5.3 °C | Equatorial axis | 17 | 0 | 0.983 |
| −5.3 °C | Polar axis | 17 | 0 | 1.017 |
| −14.4 °C | Mass | 21 | 2 | 1.017 |
| −14.4 °C | Equatorial axis | 21 | 2 | 0.998 |

Three crystals fall outside. Crystal 33, grown at −5.8 °C and measured at 5 minutes, weighs 1.7 times the model's
crystal. Crystals 94 and 95, grown at −13.9 and −14.0 °C and measured at 10 and 30 minutes, are much smaller than the
model and than their repeats: crystal 94 spans 0.69 mm where the model spans 1.38 mm, and crystal 95 spans 2.41 mm
where the model spans 3.99 mm and crystal 114, its repeat at −14.6 °C, spans 3.77 mm. All three sit near the edge of
their temperature window.

Fall speed passes at all eleven −14.4 °C ages, with the largest difference −1.6 cm/s at 30 minutes. At −5.3 °C it
passes to 10 minutes and fails at 12, 15, 20, 25 and 30 minutes, where the model falls 3.3–6.1 cm/s fast, 21–26% above
the relation. The calibration fitted these mass and axis histories and left fall speed out of its objective, so its
agreement on mass and axes is a fit to these data, while its warm fall speed is a measured failure.

## The improved component

| Temperature | Mass: outside of n; median | Equatorial: outside of n; median | Polar: outside of n; median | Speed: ages outside; largest difference |
|---|---|---|---|---|
| −3.7 °C | 2 of 14; 1.255 | 2 of 14; 1.239 | 4 of 13; 0.760 | 7 of 9; +10.2 cm/s |
| −5.3 °C | 9 of 17; 1.485 | 0 of 17; 0.955 | 3 of 17; 0.993 | 10 of 10; +22.4 cm/s |
| −8.6 °C | 4 of 10; 1.338 | 0 of 10; 0.913 | 9 of 10; 1.903 | 7 of 7; +18.7 cm/s |
| −10.6 °C | 0 of 12; 1.148 | 1 of 12; 1.064 | 1 of 8; 0.841 | 6 of 6; +5.6 cm/s |
| −12.2 °C | 6 of 15; 1.454 | 6 of 15; 1.281 | — | 6 of 9; +5.0 cm/s |
| −14.4 °C | 4 of 21; 0.924 | 11 of 21; 0.897 | — | 3 of 11; −5.8 cm/s |
| −16.5 °C | 5 of 20; 1.177 | 10 of 20; 1.251 | — | 4 of 9; −8.5 cm/s |
| −18.2 °C | 2 of 9; 1.310 | 3 of 9; 1.240 | — | 2 of 7; +3.6 cm/s |
| −20.1 °C | 1 of 8; 0.934 | 0 of 8; 0.915 | 0 of 8; 1.083 | 6 of 6; +14.3 cm/s |
| −22.0 °C | 0 of 8; 0.848 | 0 of 8; 0.962 | 0 of 8; 1.106 | 4 of 4; +21.3 cm/s |

The component grows too much mass at most temperatures, by a median 1.3–1.5 times at −5.3, −8.6, −12.2 and
−18.2 °C, and its thick plates at −8.6 °C come out nearly twice as thick as the measured ones. Its −14.4 °C dendrite
shows why a median can hide a failure: the span runs 0.55–0.68 of the measured crystals' over the first seven
minutes and 1.3–2.3 times theirs at 25 and 30 minutes, so the median sits near one while half the crystals fall
outside. The component's fall speed is the largest error. It falls 1.1–2.0 times faster than the relation at most
temperatures and ages, and its dendrites and sectors at −14.4 and −16.5 °C fall too slowly after 10–20 minutes. These errors, measured
crystal by crystal, size the band carried through the lunar column ([Earth error band](earth_error_band.md)).

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.data_sized_check run
```

The [runner](data_sized_check.py) reads the committed measurements, verifies each retained component trace against the
hash in the [shape diagnosis](results/mcsnow_shape_diagnosis.json), and runs the calibration's frozen executable through
the [science check's](mcsnow_science_check.py) native interface. The compact [product](results/data_sized_check.json)
holds every crystal's verdict, the scatter estimates and the numerical checks; raw executable inputs and outputs are in
`research/runs/snow_agent/data_sized_check/`, on the project-space drive.
