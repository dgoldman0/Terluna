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
| [ck.py](ck.py), [validate_ck.py](validate_ck.py) | Correlated-k thermal radiation (16 bands × 16 g-points) built from the same spectroscopy, with ozone; its line-by-line check writes `results/ck_validation.json` |
| [trace_gases.py](trace_gases.py) | Forcing of inert fluorinated trace gases (SF6, CF4, NF3) on the Earth control and the Moon, with band saturation; writes `results/trace_gases.json` |
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

These come from [results/inverse_climate.csv](results/) (67 columns, plus 54 of them
repeated behind each of three shields; the schema,
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
- **The spectral shield costs 7–10 K.** Sunlight filtered through each shield
  in the protection domain's product
  ([shield_transmission.json](../../protection/spectra/shield_transmission.json))
  loses 9% in the titania–silica stack (309 of 340 W/m² reaches the Moon) and
  3.5–5% in the idealised edge filters. At 290 K the lunar clear-sky ASR falls
  from 243.5 to 230.4 W/m² (1.2 atm, titania). The balance temperatures below
  are for the 1.2 atm Moon. The 200-nm edge's ASR leaves out the ozone that
  shield lets form.

  | Cloud effect C (W/m²) | Unfiltered | Titania stack | 200-nm edge | 310-nm edge |
  |---|---|---|---|---|
  | 0 | 301.1 K | 291.9 K | 294.5 K | 294.1 K |
  | −20 | 286.7 K | 277.0 K | 279.5 K | 279.0 K |
  | −40 | 272.1 K | 263.5 K | 265.9 K | 265.4 K |

  The saturated column's excess over the runaway limit shrinks from 21 W/m² to
  9 W/m² (titania) or 11 W/m² (edges). It is still positive, so clouds and dry
  regions must still supply the margin. The 1.0-atm Moon differs from these
  temperatures by −0.6 to +1.7 K, warmer the stronger the cloud cooling.

None of this is a climate prediction: clouds, the month-long day and night, the
poles and circulation are all outside the calculation. A single column through
the month-long day and night is in [climate/lunar_day.py](../../climate/lunar_day.py). Ozone and a
self-consistent stratosphere are added in [../middle_atmosphere](../middle_atmosphere/),
which moves these balance temperatures: behind the titania stack the
stratosphere is colder than the 200 K assumed here, and with C = −20 W/m² the
balance rises to about 284 K.

## Correlated-k scheme (2026-09-24)

[ck.py](ck.py) reproduces the line-by-line thermal fluxes for the calculations
that need many radiation calls. The numbers are in [validation.json](validation.json)
and `results/ck_validation.json`:

- OLR is within 0.4 W/m² for Earth and Moon columns at 288–290 K, and 1.4 W/m²
  for a saturated 330 K lunar column.
- Surface downward longwave is within 0.8–2.9 W/m².
- Ozone 9.6-µm heating is within 0.08 K/day at every level.
- CO2 15-µm heating is within about 5% up to 30 Pa (Earth column) and 5 Pa
  (lunar column). Above that it is off by 20–40%, and several times in the
  model's top layer: 16 g-points do not resolve the line peaks that only the
  uppermost layers see.
- Above a few pascals, collisions no longer keep the CO2 15-µm band in local
  thermodynamic equilibrium. `CKLongwave.fluxes_nonlte` treats the band
  (500–820 cm⁻¹) as a two-level system above 50 Pa, quenched by N2, O2 and
  atomic oxygen (López-Puertas & Taylor 2001), and solves its source function
  together with the radiation field on the correlated-k optical depths. With
  frequent collisions it returns the equilibrium fluxes (within 10⁻⁴ K/day);
  at 0.3 Pa on the Moon it cuts the cooling about four-fold.

## Trace greenhouse gases (2026-09-24)

[trace_gases.py](trace_gases.py) asks how much warming inert, chlorine- and
bromine-free fluorinated gases could add. Each gas is modelled as a well-mixed
absorber, grey across its strongest band, and its forcing is computed line by
line on columns at 288 K. The IPCC AR6 radiative efficiency fixes each band's
mean cross-section through the Earth control. All-sky values assume the
all-sky adjusted forcing is 0.6–0.9 of the clear-sky instantaneous forcing.
Results are in `results/trace_gases.json`.

| Gas (band) | Earth, W/m² per ppb | Moon 1.2 atm, W/m² per ppb (thin limit) | Most one gas can give on the Moon (all-sky) | ppb for 2 W/m² on the Moon |
|---|---|---|---|---|
| SF6 (925–955 cm⁻¹) | 0.57 | 3.5 | 3.7–5.6 W/m² | 0.7–1.0 |
| NF3 (890–920 cm⁻¹) | 0.20 | 1.3 | 4.0–6.1 W/m² | 1.9–2.4 |
| CF4 (1270–1290 cm⁻¹) | 0.10 | 0.46 | 0.7–1.1 W/m² | not reachable |

- A part per billion does five to seven times more on the Moon than on Earth,
  because each ppb is a six-to-seven-fold thicker absorber over each square
  metre.
- For the same reason each band saturates within a few ppb. One gas cannot
  supply more than about 4–6 W/m², however much is added. SF6, NF3 and CF4
  together, in their separate bands, give roughly 8–13 W/m², about 6–10 K at
  the Moon's sensitivity. The grey-band curve is a lower bound once a band is
  saturated, because real bands keep absorbing weakly in their wings and hot
  bands.
- CF4's band sits where the Moon's deeper water-vapour column already absorbs,
  so it is nearly useless there.
- These gases are inert and harmless to breathe at ppb levels, and carry no
  chlorine or bromine to destroy ozone. On Earth they last centuries to tens of
  millennia (SF6 3,200 years, NF3 570, CF4 50,000). Behind a shield that blocks
  the far ultraviolet that breaks them up high in the atmosphere they would
  last longer still: a lever that is effectively permanent once pulled.

## What it is and is not

It is a radiation calculation for columns whose temperature and humidity are
prescribed. It is not a solved climate, a cloud model or a habitability result.
Clouds, ozone, day and night, the poles and circulation are all outside it.
Comparing ASR with OLR tells whether a clear column at a given surface
temperature would warm or cool, and how strongly lunar gravity changes that
compared with Earth.
