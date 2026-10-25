# Lunar snow crystal sizes over water supply and crystal number

The crystal arriving at the lunar freezing level follows the water each crystal receives. Across snow supplies of
0.1–10 mm/h of water and incoming crystal fluxes from the reference flux of the earlier lunar runs down to a thousandth
of it, the water per crystal runs from 0.0027 to 271 mg and the arriving crystal from 0.40 to 142 mm across. Each size
depends on the ratio of supply to crystal number alone: crystals take up 95–100% of their share up to 27 mg, and 77–91%
at the largest shares, where the descent ends before the crystal has gathered all of it. The crystal model's measured
Earth errors move sizes above 2 mm by −15%
to −25% and +18% to +31%, and the sub-millimetre crystals by up to −41%. These are single crystals; the flakes they
form are plan step 4.

## The calculation

The column is the 301 K reference ladder of the earlier lunar runs ([lunar column](mcsnow_lunar_growth.md)): 16.28 km
from −15 °C down to the freezing level, at water saturation, with lunar gravity and no vertical air motion. A steady
supply Q is shared among a steady incoming flux F of 4 µm spheres, so each crystal may deposit at most Q/F over the full
descent, in proportion to the distance it has fallen. The grid crosses five supplies (0.1, 0.3, 1, 3 and 10 mm/h) with
seven fluxes (1, 0.3, 0.1, 0.03, 0.01, 0.003 and 0.001 times the reference 10,232 crystals per square metre per second).
Every point runs with the frozen improved McSnow component and at the two combined edges of its Earth error band, the
edge combinations that most reduce and most increase the arriving size in this column
([Earth error band](earth_error_band.md)). All 105 runs complete at the freezing level and keep their deposition within
the allowance.

| Water per crystal | Example: supply; flux | Arriving crystal | Band | Descent |
|---:|---|---:|---:|---:|
| 0.0027 mg | 0.1 mm/h; reference | 0.40 mm | 0.23–0.47 mm | 133 h |
| 0.027 mg | 1 mm/h; reference | 1.82 mm | 1.55–2.14 mm | 76 h |
| 0.081 mg | 3 mm/h; reference | 3.14 mm | 2.66–3.70 mm | 54 h |
| 0.27 mg | 10 mm/h; reference | 5.60 mm | 4.70–6.65 mm | 40 h |
| 0.81 mg | 3 mm/h; a tenth | 9.51 mm | 7.89–11.35 mm | 32 h |
| 2.7 mg | 1 mm/h; a hundredth | 17.0 mm | 14.0–20.4 mm | 26 h |
| 8.1 mg | 3 mm/h; a hundredth | 29.1 mm | 23.7–35.1 mm | 22 h |
| 27 mg | 10 mm/h; a hundredth | 52.2 mm | 42.3–63.3 mm | 20 h |
| 81 mg | 3 mm/h; a thousandth | 87.4 mm | 70.2–107.5 mm | 19 h |
| 271 mg | 10 mm/h; a thousandth | 142 mm | 106–186 mm | 18 h |

Points with equal water per crystal give the same crystal: 3 mm/h shared by the reference flux and 0.3 mm/h shared by
a tenth of it both give 3.14 mm. Between 0.027 and 27 mg the arriving span grows as the water per crystal to the power
0.486, close to the square root that suits thin crystals whose mass grows with the square of their span. A crystal that
receives more water also grows larger sooner, falls faster and spends less time in the cloud, so the descent shortens
from 133 hours for the smallest crystals to 18 hours for the largest. At 271 mg the crystal deposits 210 mg of its share
before reaching the freezing level; under unlimited saturation the same column gives 147 mm and 232 mg, the most this
crystal can take up in the time it has.

## Reading the range

The grid does not say which part of it lunar storms occupy. The archived lunar storm holds most snowy cells below ice
saturation, at a median 89% relative humidity, and its crystals peak at a median 0.73 mm before most sublimate
([lunar storm](mcsnow_weather.md)); along the same paths held at water saturation the median peak reaches 6.0 mm. A
storm sets its supply and its crystal number together, through its updrafts, cooling and ice nucleation, and a coupled
calculation is needed to place it on this grid.

Above about a centimetre the component carries single crystals far past any measured: Takahashi and co-authors'
laboratory crystals reached about 4 mm in 30 minutes. There the arriving size measures the water each crystal would
gather, while aggregation (plan step 4) and breakup in turbulence (step 5) decide what actually falls.

## Limits

The supply is imposed: it neither responds to the crystals nor varies in time and height, and each crystal receives it
in proportion to the distance it falls. The column is a straight profile at water saturation without updraughts. The
crystals start as 4 µm spheres of one size and grow alone, without collisions or depletion of one another's vapour.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B -m atmosphere.snow.lunar_size_range
```

The [runner](lunar_size_range.py) reuses the band's column program and writes the compact
[product](results/lunar_size_range.json); raw runs are in `research/runs/snow_agent/earth_error_band/column_runs/`.
