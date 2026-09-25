# Climate

The original [scalar thermal-response table](reference/thermal_response.csv) and [September baseline](../research/baselines/feasibility/model.py) are preserved.

[cycle.py](cycle.py) adds a finite-volume latitude–longitude energy-balance model. It uses an equal-area spherical grid, moving sunlight, prescribed land/water heat capacities and a conservative diffusive heat-transport operator. An exact linear interval propagator and a periodic-boundary solve remove arbitrary initial-temperature/spin-up choices.

## Interpretation

This is an uncalibrated thermal screen. Synthetic dry, concentrated-water and distributed-water maps are compared. Albedo, atmospheric participation, transport and a linear outgoing-radiation law are inputs. The global mean follows analytically from those inputs; a temperate mean is not independent evidence that a 1.2-atm lunar atmosphere produces it.

Moisture, clouds, winds, latent heat, freezing feedbacks, actual terrain and photochemistry remain outside the calculation. Outputs reaching freezing or far outside the reference temperature range are flagged. The reported 273–313 K fraction is a chosen grid-time diagnostic, not a habitability or liquid-water certification.

The [54-case output](../research/studies/environment_screens/results/climate_budget.csv) closes its discrete global energy balance and periodic condition. Timestep tests pass, while a 6x24 to 12x48 spatial comparison changes some coarse-cell averaged temperatures by about 2 K; spatial resolution uncertainty remains. [Findings](../research/findings.md) explain the exact checks and boundaries.

## The month-long day in one column (2026-09-25)

[lunar_day.py](lunar_day.py) follows a single column of air over land or sea
through the 29.5-day lunar day and night, at the equator, 60° and 80°, and an
Earth control through its 24-hour day. The column is the middle-atmosphere
equilibrium at a 288 K surface ([atmosphere/middle_atmosphere](../atmosphere/middle_atmosphere/)),
regridded to 10 m layers at the ground. Thermal radiation comes from the
correlated-k scheme and sunlight from a line-by-line table for the case's shield
and composition. The ground trades heat and water with the lowest layer by bulk
formulas with a stability correction, in a prescribed 5 m/s wind:
- land conducts heat into 10 m of soil and resists evaporation like
  well-watered vegetation (100 s/m);
- the sea is a 50 m mixed layer.

Unstable air mixes upward, and condensing vapour releases its heat and falls
out. Above 3 km, humidity is capped at the reference profile.

The rest of the Moon enters as a relaxation toward the global-mean
equilibrium, standing in for the circulation:
- **Above 3 km, temperature relaxes over about half a day.** The Moon turns
  once a month, so its Rossby radius is far larger than the Moon. Gravity waves
  of about 50 m/s, as fast as on Earth, then erase temperature differences in
  the free atmosphere in about half a day. This is the weak-temperature-gradient
  regime of Earth's tropics and of slowly rotating planets, and its
  single-column form carries away the heat that convection and rain deliver
  aloft.
- **Vapour relaxes over 20 days.**

The base case uses 12 hours, and the variants 3 hours and 5 days. The soil and
sea are brought to their periodic state between days. The sea is held to 288 K
on average, and the model reports the heat that currents would have to carry to
hold it there. The albedo is 0.13 throughout.

```sh
python -m climate.lunar_day                                  # every case, into results/lunar_day/
python -m climate.lunar_day --only moon_titania_equator_land
```

A case takes 20–40 minutes on one core. A stopped run resumes: finished cases
are skipped, and unfinished ones restart from their last simulated day. Results
are in [results/lunar_day/](results/lunar_day/): `lunar_day.csv` (one row per
case) with `lunar_day.json` (schema, producer hash, evidence and reading rule),
and `series/<case>.csv` for the recorded day. All values come from the last of 8
lunar days (60 Earth days for the control).

| Case | Ground, K (day mean) | Air at 5 m, K | Warmest, hours from noon | At noon: sunlight absorbed / warm air / evaporation (W/m²) |
|---|---|---|---|---|
| Equator, land (base) | 280.1–296.9 (287.7) | 280.7–294.4 | +9 | 645 / 202 / 264 |
| Equator, land, 200-nm edge | 280.3–297.0 (287.8) | 280.9–294.4 | +7 | 653 / 206 / 265 |
| Equator, sea | 287.3–288.7 (288.0) | 287.2–288.9 | +132 | 645 / −3 / 14 |
| 60°, land | 279.6–292.8 (285.6) | 280.3–292.3 | +29 | 265 / 54 / 125 |
| 60°, sea | 287.7–288.3 (288.0) | 287.5–288.2 | +127 | 265 / 0 / 16 |
| 80°, land | 276.2–283.4 (279.4) | 277.0–283.4 | +77 | 69 / 9 / 17 |
| Equator, land, 2 m/s wind | 279.9–298.0 (288.0) | 280.6–294.5 | +15 | 645 / 208 / 281 |
| Equator, land, 8 m/s wind | 280.4–296.2 (287.5) | 280.9–294.2 | +16 | 645 / 190 / 254 |
| Equator, land, long stable tail | 280.1–296.9 (287.7) | 280.7–294.4 | +9 | 645 / 202 / 264 |
| Equator, land, weak night mixing | 278.6–296.9 (287.1) | 278.9–294.4 | +9 | 645 / 202 / 264 |
| Equator, drying soil (500 s/m) | 280.1–298.4 (288.2) | 280.7–294.6 | +8 | 645 / 320 / 76 |
| Equator, arid (10,000 s/m) | 279.6–298.8 (288.2) | 280.2–294.5 | +9 | 645 / 362 / 4 |
| Equator, land, exchange 3 h | 279.8–295.0 (286.8) | 280.4–292.3 | −20 | 645 / 228 / 225 |
| Equator, land, exchange 5 days | 284.4–301.5 (292.6) | 285.0–300.1 | +55 | 645 / 150 / 345 |
| Earth, equator, land | 293.1–305.8 (298.5) | 294.7–303.8 | +1 | 945 / 215 / 387 |
| Earth, equator, arid | 293.4–310.0 (301.4) | 295.8–306.6 | +1 | 945 / 372 / 10 |
| Earth, equator, sea | 288.0 (288.0) | 287.6–288.2 | +4 | 945 / −1 / 20 |

What the column shows:

- **At the equator, moist ground swings 17 K over the month.** The ground
  runs 280–297 K and the air 281–294 K. The warmest time is about 9 hours after
  noon; the coldest is sunrise, after the 15-day night. At noon the ground
  absorbs 645 W/m², a third less than on Earth under the same Sun, because the
  heavier air scatters and absorbs more. Evaporation takes 264 W/m² and warm air
  202 W/m². Convection carries the heat through 12 km of air, and 1.2 mm of dew
  forms over the night.
- **The free atmosphere caps the day.** Low gravity makes the Moon's
  dry-adiabatic lapse rate only 1.6 K/km. So air mixed 12 km deep ends near the
  temperature of the free atmosphere above, which the circulation holds near the
  global mean. The exchange time is the model's largest uncertainty:
  - between 3 hours and 5 days, the peak runs 295.0–301.5 K and the day mean
    286.8–292.6 K;
  - the global-mean profile is itself an assumption, since the sunlit convecting
    half may set a warmer free atmosphere.

  A GCM must settle both.
- **Dry ground hardly matters on the Moon.** Cutting evaporation almost to
  nothing raises the peak ground temperature by 2 K and leaves the air
  unchanged. On Earth the
  same change raises the peak by 4 K and the air by 3 K. The arid cases keep the
  base case's soil, its 0.13 albedo and its humid air aloft. A real desert of
  bright sand under dry air would swing more.
- **The sea hardly swings; currents decide its level.** The 50 m mixed layer
  varies by 1.4 K over the month. Holding it at 288 K under these clear skies
  needs currents to remove 147 W/m² at the equator and 30 W/m² at 60°. Clouds and
  the circulation decide where it really settles.
- **Higher latitudes are cooler but not frozen.** Land at 60° runs 280–293 K,
  and at 80° 276–283 K. The free atmosphere, held at the global mean, keeps even
  80° above freezing in this model. The poles are where the relaxation is least
  certain, so a GCM must test this.
- **Wind and mixing move the extremes by 1–2 K.**
  - Calm air (2 m/s) adds 1.1 K at the peak, and 8 m/s wind removes 0.7 K.
  - Weak night mixing lowers the minimum by 1.5 K and keeps the air at the ground
    saturated for 156 of the 354 night hours, which means dew or fog.
  - The form of the stable tail makes no difference.
  - The 200-nm edge instead of the titania stack changes the ground by 0.1 K.
- **The Earth control reproduces humid equatorial land.** The ground runs
  293–306 K and the air 295–304 K, a 9 K daily range. At noon the ground receives
  389 W/m² of downward infrared, and evaporation averages 151 W/m². Observed
  equatorial rainforest has air near 300 K, a daily range of about 8–10 K and
  downward infrared near 400 W/m². The model's day varies by about half a kelvin
  from one day to the next; every lunar case repeats within 0.1 K.
- **So the month-long day swings the lowlands only somewhat more than one humid
  tropical day on Earth**, spread over a month. That holds under these
  assumptions: clear skies, a fixed wind, and a free atmosphere held near the
  global mean.

The single column has no clouds, no winds of its own, no terrain, and neither
soil moisture nor vegetation that responds. Its exchange with the rest of the
Moon is a stated relaxation. It shows which processes set the swing at a place
and how far each assumption moves it. It is not a forecast of lunar weather;
that needs the GCM.

References: Louis (1979), Boundary-Layer Meteorology 17, 187; Manabe and
Wetherald (1967), J. Atmos. Sci. 24, 241; Sobel and Bretherton (2000), J.
Climate 13, 4378; Mills and Abbot (2013), Astrophys. J. Lett. 774, L17.

## Three-dimensional preparation

[gcm/](gcm/) holds the boundary files (land fraction, elevation and ocean depth
for 25% and 35% water on common GCM grids), the planet parameters, the
experiment plan and a restartable runner for ExoPlaSim, whose radiation is
calibrated against the line-by-line model. In ExoPlaSim the design case settles at
295.5 K, with the air near the ground within about a kelvin of that from equator
to pole. PlaSim's clouds, which cool it by only about 4 W/m², are the main caveat.

## Next step

The atmosphere domain's [radiative–convective column](../atmosphere/radiative_convective/)
now computes clear-sky outgoing longwave radiation and absorbed sunlight against
surface temperature, line by line, for the Moon and an Earth control. Those are
the quantities this screen prescribes through A, B and the albedo. Using them
here is the next step, but it needs an explicit cloud assumption, because the
column is clear-sky. Water/ice treatment and measured terrain follow. The present temperature/irradiance interface can supply biological requirement calculations, with all atmospheric and trait assumptions retained.

Run `python -m pytest climate` and `python -m research.studies.environment_screens.run` from the repository root. This folder supplies no GCM weather prediction.
