# The crystal model's Earth errors in the lunar column

The improved McSnow component's measured Earth errors move lunar crystal size by a factor of 1.4–1.6 wherever water
limits growth. Water supply and crystal number move it much further: a hundredfold change in crystal number changes the
water-limited crystal by a factor of 9.3, archived storm humidity and continuous saturation differ by a factor of 8.2,
and the cases together span a factor of 203, from a 0.73 mm median peak in the archived storm to 147 mm under unlimited
saturation. Only where growth runs unchecked do the errors compound, to a factor of 2.3 in the saturated storm and 5.4
in the unlimited column, a case that needs about 2,800 times the water the column allows.

## How the band is set

The improved McSnow component, which every lunar run uses, fails the [data-sized check](data_sized_check.md): 83 of
its 332 mass and axis comparisons with Takahashi and co-authors' (1991) crystals fall outside the scatter between
repeat crystals, and it falls 1.1–2.0 times faster than those crystals at most temperatures. The two-temperature
calibration passes the same check on mass and axes, but it holds no rule between its two temperatures and cannot run
through a lunar column. The band therefore carries the component's own measured error, at ten temperatures from −3.7
to −22 °C. That range covers the lunar column, from −15 °C to the
freezing level, and the storm paths' release levels between −20 and 0 °C; the multipliers act at every temperature,
including the colder air some paths drift into. Each kind of error becomes one multiplier on the component's own
physics, in a separate build that leaves the frozen component untouched:

| Error | Multiplier | What it changes |
|---|---|---|
| Mass growth | The crystal's capacitance | The vapour flux to the crystal, at every size and temperature |
| Axis allocation | The ventilated inherent growth factor | How new ice divides between the polar and equatorial axes |
| Fall speed | The terminal speed | Settling, the per-crystal water allowance that follows it, and ventilation |

At each laboratory temperature, each multiplier is sized so that the perturbed Earth crystal (a 4 µm sphere at
101 kPa and water saturation, as in the check) has zero median log error against the measured crystals there: mass
against the crystals' masses, axes against their aspect ratios where both axes were measured, and fall speed against
the published relation at the crystals' ages. The smallest and largest multiplier over the ten temperatures form the
band's edges. The component's errors differ in sign and size between temperatures, while each multiplier acts at every
temperature at once, so an edge applies the component's worst measured error of one kind throughout the column.

| Temperature | Mass (capacitance) | Axes (growth factor) | Fall speed |
|---|---:|---:|---:|
| −3.7 °C | 0.871 | 1.165 | 0.726 |
| −5.3 °C | 0.838 | 1.004 | 0.604 |
| −8.6 °C | 0.836 | 0.818 | 0.652 |
| −10.6 °C | 0.925 | 1.074 | 0.733 |
| −12.2 °C | 0.853 | — | 0.784 |
| −14.4 °C | 1.043 | — | 0.929 |
| −16.5 °C | 0.919 | — | 1.159 |
| −18.2 °C | 0.862 | — | 0.827 |
| −20.1 °C | 1.043 | 0.952 | 0.607 |
| −22.0 °C | 1.109 | 0.967 | 0.555 |
| **Band edges** | **0.836–1.109** | **0.818–1.165** | **0.555–1.159** |

The planar crystals at −12.2 to −18.2 °C carry no measured polar axis, so they set no axis multiplier. Growth feeds on
size, so small multipliers remove large mass errors: a 16% smaller capacitance removes the 48% excess mass at −5.3 °C.
The fall-speed multipliers depart furthest from one, because the component falls too fast at nine of the ten
temperatures, by up to 1.8 times for the −22 °C thick plates.

With every multiplier at one, the perturbed builds reproduce the frozen component exactly: the two retained column runs
byte for byte, and the retained archived-vapour storm run in its peaks, terminations, states, ages, deposits and
positions. The Earth library at factor one agrees with the component's retained 1/128 s traces to within 0.33% at every
crystal age, the difference between its adaptive steps and the traces' fixed ones.

## Lunar crystal sizes at the band's edges

Every case starts from 4 µm spheres. The column cases follow one crystal down the 16.28 km reference column at water
saturation, sharing a 3 mm/h water supply among the incoming crystal flux of the earlier lunar runs or a tenth or a
hundredth of it, or growing without a supply limit; sizes are maximum dimensions at the freezing level. The storm cases
follow 354 crystals through the archived lunar storm's humidity, temperature and winds, or through the same paths held
at water saturation; sizes are each path's peak dimension.

| Case | Unperturbed | Band, smallest–largest | Factor |
|---|---:|---:|---:|
| Column, 3 mm/h shared by the reference crystal flux | 3.14 mm | 2.66–3.70 mm | 1.39 |
| Column, a tenth of the crystals | 9.51 mm | 7.89–11.35 mm | 1.44 |
| Column, a hundredth of the crystals | 29.1 mm | 23.7–35.1 mm | 1.48 |
| Column, unlimited saturation | 147 mm | 73–398 mm | 5.42 |
| Archived storm humidity, median peak | 0.727 mm | 0.553–0.863 mm | 1.56 |
| Storm paths at continuous saturation, median peak | 5.99 mm | 4.56–10.70 mm | 2.35 |

Where the column's water allowance limits growth, each crystal deposits its share of the supply whatever its growth
rate, so the mass and fall-speed errors barely move its size. At the reference crystal flux they change the arriving
crystal by less than 0.1% and at most 2.4%, although the slowest fall-speed edge stretches the descent from 54 to
96 hours. The axis error remains, because it reshapes the same mass: putting more of the ice into the polar axis makes
the crystal 13% smaller across, and putting less there makes it 17% larger. In the archived storm, growth is limited by
the time each crystal spends in moist air; the mass edges move the median peak by −13% and +9%, and the combined edges
by −24% and +19%.

Where nothing limits the water, the errors compound through size and residence. Under unlimited saturation a crystal
that grows faster becomes larger, is ventilated more strongly and grows faster still, and a crystal that falls more
slowly stays longer in the cloud: the combined edges span 73 to 398 mm around the unperturbed 147 mm. In the saturated
storm paths, the edge that puts more ice into the polar axis sends 40 of the 354 paths to the 1 m size censor as
extreme needles, the elongation runaway already recorded for the component, and their 90th-percentile peak sits at
that censor.

The answer to the doubt that opened this step therefore depends on water. For crystals that share a finite supply,
which every storm imposes, the crystal model's Earth errors amount to a correction of tens of percent, while the water
each crystal receives sets its size across a factor of about 200 in the cases computed here. The errors would matter at
the factor-of-five level only in growth that no modelled storm supplies. Halving the column timestep changes no arriving
size by more than 0.032%.

## Limits

The multipliers measure errors over at most 30 minutes of growth at Earth gravity and pressure; lunar crystals grow for
18–96 hours at a sixth of the gravity, and the component's errors there are unmeasured. The calculations follow single
crystals in prescribed air, with no depletion of the vapour, latent heating of the air or collisions. The archived
storm is the older A28_dim5 case, saved every ten minutes, rather than the 301 K design climate. The sizes are maximum
dimensions of single crystals, which aggregation (plan step 4) turns into flakes.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.earth_error_band run
```

The [runner](earth_error_band.py) builds the band versions of the column program and the storm library from the pinned
McSnow source with the improved component's patches, sizes the edges, checks parity, runs the column and storm cases and
writes the compact [product](results/earth_error_band.json). Builds, sizing, raw runs and logs are in
`research/runs/snow_agent/earth_error_band/`, on the project-space drive.
