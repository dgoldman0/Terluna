# Radiative–convective column

A line-by-line, clear-sky energy budget for a prescribed atmospheric column: how
much thermal radiation leaves the top of the atmosphere (OLR) and how much sunlight
is absorbed (ASR), as functions of surface temperature. The Open Moon and an Earth
control run through the same code, so every lunar number has an Earth
counterpart. Equations, data, tests and limits are in [METHODS.md](METHODS.md).

| File | Holds |
|---|---|
| [thermodynamics.py](thermodynamics.py) | Saturation vapour pressure, the non-dilute moist adiabat, and the column with GM/r² gravity |
| [spectroscopy.py](spectroscopy.py) | HITRAN lines (multi-mesh Voigt), the MT_CKD 4.3 water continuum, collision-induced absorption, Rayleigh scattering |
| [optics.py](optics.py) | Layer optical depths, in parallel, with a cached node table for air-broadened CO2 and O2 |
| [longwave.py](longwave.py), [shortwave.py](shortwave.py) | Thermal fluxes; delta-scaled two-stream solar fluxes with a pseudo-spherical direct beam |
| [climate.py](climate.py) | Scenarios and the OLR/ASR budget at a surface temperature |
| [run.py](run.py) | Named sweeps, written to `results/` as a data product with schema, hashes and evidence statement |
| [analysis.py](analysis.py) | Balance temperatures for stated cloud effects, runaway limits and sensitivities from the sweep table |
| [validation.json](validation.json) | Cross-code (PyRADS), Monte Carlo and convergence records |
| [inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py) | The external spectroscopic inputs (about 176 MB), with URLs, sizes, hashes and credits |

## Run

From the repository root:

```sh
python -m atmosphere.radiative_convective.fetch_inputs --download   # once; hash-checked, ignored by Git
python -m pytest atmosphere/tests/test_radiative_convective.py
python -m atmosphere.radiative_convective.run --list
python -m atmosphere.radiative_convective.run --only earthlike_humidity
```

A column takes about a minute for each of the thermal and solar parts on eight
cores. The CO2 node table (`cache/`, ignored) is built on first use. The runner
appends rows as they finish and resumes after an interruption.

## First results (2026-09-24)

These come from [results/inverse_climate.csv](results/) (67 columns; the schema,
producer hashes and evidence statement are in `inverse_climate.json`; derived
numbers are in `summary.json`). Settings: clear sky, 400 ppm CO2, stratosphere at
200 K, surface albedo 0.13, Manabe–Wetherald humidity unless marked saturated.
Balance temperatures add an assumed, constant cloud effect C to the absorbed
sunlight.

| | Earth control, 1 atm | Moon, 1.0 atm | Moon, 1.2 atm |
|---|---|---|---|
| OLR at 290 K (W/m²) | 276.7 | 233.1 | 227.8 |
| Clear-sky planetary albedo at 290 K | 0.160 | 0.267 | 0.284 |
| Water vapour column at 290 K (kg/m²) | 23.5 | 144 | 136 |
| Tropopause height at 290 K (km) | 12.4 | 78.3 | 75.1 |
| Balance with no clouds (K) | 293.7 | 300.5 | 301.1 |
| Balance with C = −20 W/m², about Earth's observed value (K) | 285.3 | 287.4 | 286.7 |
| Balance change per 10 W/m² of C, between C = −10 and −30 (K) | 4.3 | 7.0 | 7.6 |
| Saturated-column OLR limit (W/m²) | 299.5 | 239.4 | 239.0 |
| Clear-sky ASR of that saturated column (W/m²) | 296.9 | 263.8 | 260.0 |
| OLR change per CO2 doubling at 290 K (W/m²) | −4.1 | −4.0 | −4.0 |
| Surface net longwave loss at 290 K (W/m²) | 87 | 25 | 25 |

What these calculations show:

- **The two lunar effects nearly cancel.** Lunar gravity puts about six times
  more water vapour and CO2 over each unit of pressure, so OLR is about 45 W/m²
  lower at the same surface temperature. The same deep column scatters sunlight
  back (Rayleigh optical depth 0.75 at 550 nm against Earth's 0.10), which
  raises the clear-sky albedo from 0.16 to 0.28. The two almost offset each
  other: with Earth-like clouds, humidity and CO2, the balance temperature is
  within about 2 K of the Earth control's.
- **The Moon is more sensitive.** The net budget changes less per kelvin, so the
  same cloud, albedo or CO2 change moves the lunar surface temperature 1.6–1.8
  times as far as Earth's.
- **The margin against a runaway greenhouse is smaller.** The saturated OLR
  limit falls from 299 to 239 W/m². A fully saturated, cloud-free Moon absorbs
  about 21–24 W/m² more than that limit, so it has no balance, while the Earth
  control stays just below its own. Clouds and dry regions (not modelled here)
  must supply that margin. The Manabe–Wetherald columns show no runaway up to
  320 K.
- **CO2 is a modest control.** It acts about equally per doubling on both
  bodies, about 4 W/m². Between 150 and 1000 ppm it shifts the lunar balance by
  roughly ±4 K around the 400-ppm case.
- **The stratosphere matters.** At 290 K, a 150 K instead of 200 K stratosphere
  lowers lunar OLR by 7 W/m², which warms the surface. The stratospheric
  temperature, and whether any ozone heats it, is a real uncertainty for the
  Moon.
- **Clear-sky nights lose heat slowly.** The water-rich lower air is nearly
  opaque in the infrared, so at 290 K the ground loses about 25 W/m² net
  (Earth: 87). The ground stops cooling radiatively about 5 K below the air
  temperature, against about 17 K on Earth. How far nights actually cool needs
  a model of the ground and the air just above it.
- **1.0 and 1.2 atm hardly differ in climate.** Balance temperatures agree to
  within 1–2 K across the cloud effects tried, and the runaway limit is the same.
- **The spectral shield is not yet in the sunlight.** These runs use unfiltered
  sunlight. The protection model's titania–silica stack passes about 90% of
  sunlight and none below about 330 nm ([uv_absorption.csv](../../protection/results/uv_absorption.csv),
  `results.json`). An approximate version of that filter (its stored ultraviolet
  points and 96.3% above 400 nm) lowers the Moon's clear-sky ASR at 290 K from
  243.5 to 230.4 W/m². That is about 13 W/m², roughly 9–10 K cooler, and it
  leaves no ultraviolet to form ozone. The full spectral transmission belongs in
  the solar calculation before any balance temperature is quoted.

None of this is a climate prediction: clouds, the month-long day and night, the
poles, circulation and ozone are all outside the calculation.

## What it is and is not

It is a radiation calculation for columns whose temperature and humidity are
prescribed. It is not a solved climate, a cloud model or a habitability result.
Clouds, ozone, day and night, the poles and circulation are all outside it.
Comparing ASR with OLR tells whether a clear column at a given surface
temperature would warm or cool, and how strongly lunar gravity changes that
compared with Earth.
