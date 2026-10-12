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
| [trace_gases.py](trace_gases.py) | Forcing of inert fluorinated trace gases (SF6, CF4, NF3) on the Earth control and the Moon from measured cross-sections; writes `results/trace_gases.json` |
| [validation.json](validation.json) | Cross-code (PyRADS), Monte Carlo and convergence records |
| [inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py) | The external spectroscopic inputs (about 176 MB), with URLs, sizes, hashes and credits |

## Run

From the repository root:

```sh
python -m atmosphere.radiative_convective.fetch_inputs --download   # once; hash-checked, ignored by Git
python -m atmosphere.radiative_convective.fetch_added_inputs --download   # CH4 and N2O lines, the same way
python -m pytest atmosphere/tests/test_radiative_convective.py
python -m atmosphere.radiative_convective.run --list
python -m atmosphere.radiative_convective.run --only earthlike_humidity
python -m atmosphere.radiative_convective.reduced_gases      # methane and N2O, 15-23 minutes on six workers
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

## Trace greenhouse gases (2026-09-25)

[trace_gases.py](trace_gases.py) asks how much warming inert, chlorine- and
bromine-free fluorinated gases could add. Their forcing is computed line by line
on columns at 288 K from measured infrared cross-sections:
- SF6 and NF3 from the PNNL database (Sharpe et al. 2004), at 278–323 K, with
  every band from 560–600 cm⁻¹ up;
- CF4's strong band (1255–1290 cm⁻¹) from NCAR (Massie et al. 1991), at
  203–293 K.

These are the sets HITRAN serves, taken from AER's public copy for LBLRTM and
pinned in [inputs.json](inputs.json). Each layer uses the cross-sections at its
temperature, and the coldest set below that. Every level uses the spectra as
measured, without adjusting for pressure (760 Torr for SF6 and NF3).

All-sky values scale the clear-sky curves by the factor that brings the Earth
control onto the IPCC AR6 radiative efficiency: 0.56 for SF6, 0.53 for NF3 and
0.33 for CF4. The factor absorbs clouds, stratospheric adjustment and the
overlapping gases these columns lack (N2O, CH4 and O3). CF4's factor is the
smallest because its band lies under N2O and CH4. Results are in
`results/trace_gases.json`, beside the earlier grey-band estimates.

| Gas | Earth, W/m² per ppb | Moon 1.2 atm, W/m² per ppb (thin limit) | Moon at 1 / 10 / 100 / 1000 ppb | ppb for 2 / 5 / 10 / 20 W/m² on the Moon |
|---|---|---|---|---|
| SF6 | 0.57 | 3.3 | 1.4 / 3.5 / 8.5 / 17.9 W/m² | 2.3 / 24 / 150 / over 1000 |
| NF3 | 0.20 | 1.3 | 1.1 / 5.1 / 12.7 / 24.1 W/m² | 2.0 / 9.4 / 49 / 450 |
| CF4 | 0.10 | 0.45 | 0.18 / 0.46 / 0.68 / 0.70 W/m² | not reachable |

The Moon's values are all-sky estimates.

- **A part per billion does five to six times more on the Moon than on Earth.**
  Each ppb is a six-to-seven-fold thicker absorber over each square metre. In
  the thin limit the ratio is 5.8 for SF6 and 6.4 for NF3. For CF4 it is 4.6,
  because the Moon's deeper water-vapour column partly covers its band.
- **The strongest band saturates within a few ppb, but the rest keep adding.**
  The measured spectra's wings, hot bands and weaker bands go on absorbing. At
  1 ppm, SF6 gives 17.9 W/m² and NF3 24.1 W/m², and both are still rising. The
  earlier grey-band model kept only each gas's strongest band, so it capped
  every gas at 4–6 W/m². That cap was an artefact of the model.
- **Tens of watts per square metre are within reach.** The balance runs give
  about 1.6 W/m² per K, so 10 W/m² is roughly 6 K before feedbacks. That takes
  about 150 ppb of SF6 or 50 ppb of NF3. The two gases' bands overlap near
  900–950 cm⁻¹, so their forcings would not simply add; they were not computed
  together.
- **CF4 remains a weak lever.** Its strong band gives at most about 0.7 W/m² on
  the Moon. The data cover only that band, so its weaker bands would add a
  little. With N2O present, as in the photochemistry cases, it would give less.
- **These gases are inert and long-lived.** They are harmless to breathe at ppb
  levels and carry no chlorine or bromine to destroy ozone. On Earth they last
  centuries to tens of millennia (SF6 3,200 years, NF3 570, CF4 50,000). Behind
  a shield that blocks the far ultraviolet that breaks them up high in the
  atmosphere they would last longer still: a lever that is effectively
  permanent once pulled.

## Methane and N2O in the deep column (2026-10-12)

[reduced_gases.py](reduced_gases.py) answers ecology's question 3, which the author made the next work on 12
October: how far methane and N2O warm the Moon. It adds HITRAN's CH4 and N2O lines (0–3,500 cm⁻¹, pinned in
[inputs_added.json](inputs_added.json) beside inputs.json, which other products bind by hash) to the clear-sky
columns at 290–300 K and records the drop in outgoing longwave, the instantaneous forcing at the top, beside the net
change at the tropopause. A gas's optical depth is computed once per column for 1 ppm and scaled. The Earth control
sets the method against AR6's effective forcing for the rise from 1750 to 2019: this column gives 0.94 W/m² of
clear-sky forcing for methane against AR6's 0.54, and 0.32 for N2O against 0.21. Those ratios, 0.58 and 0.65, carry
Earth's clouds, adjustments and methane's absorption of sunlight into the Moon's all-sky estimates. Results are in
`results/reduced_gases.json`, a run of about 38 minutes on six workers, the levels below Earth's 1750 air and their
pairs included.

The Moon's 1.2 atm column at 295 K, each gas from none, the way the planned climate has it:

| Gas | Clear sky, W/m² | All-sky estimate, W/m² | Warming at 1.0–1.6 W/m²/K |
|---|---|---|---|
| Methane, 0.01 ppm | 0.24 | 0.14 | 0.09–0.14 K |
| Methane, 0.03 ppm | 0.52 | 0.30 | 0.19–0.30 K |
| Methane, 0.1 ppm | 1.06 | 0.61 | 0.38–0.61 K |
| Methane, 0.3 ppm (held, by the author's rule) | 1.83 | 1.05 | 0.65–1.05 K |
| Methane, 0.73 ppm (Earth, 1750) | 2.7 | 1.6 | 1.0–1.6 K |
| Methane, 1.87 ppm (Earth, 2019) | 4.0 | 2.3 | 1.4–2.3 K |
| Methane, 3.4 ppm (the screen's lowest steady state) | 5.1 | 2.9 | 1.8–2.9 K |
| Methane, 21.6–29 ppm (its central steady states) | 9.5–10.3 | 5.5–5.9 | 3.4–5.9 K |
| Methane, 230 ppm (its highest) | 16.8 | 9.7 | 6.0–9.6 K |
| N2O, 0.001 ppm | 0.09 | 0.06 | 0.04–0.06 K |
| N2O, 0.003 ppm | 0.23 | 0.15 | 0.09–0.15 K |
| N2O, 0.01 ppm | 0.59 | 0.39 | 0.24–0.39 K |
| N2O, 0.03 ppm (held, by the author's rule) | 1.28 | 0.83 | 0.52–0.83 K |
| N2O, 0.1 ppm | 2.71 | 1.76 | 1.1–1.8 K |
| N2O, 0.27 ppm (Earth, 1750) | 4.6 | 3.0 | 1.9–3.0 K |
| N2O, 0.33 ppm (Earth, 2019) | 5.1 | 3.4 | 2.1–3.3 K |
| N2O, 0.55 ppm (about a thousand years of the wild's rise) | 6.5 | 4.2 | 2.6–4.2 K |
| N2O, 1 ppm | 8.4 | 5.5 | 3.4–5.5 K |
| N2O, 4.2 ppm (about ten thousand years) | 14.2 | 9.3 | 5.8–9.2 K |
| Both at 0.1 and 0.01 ppm | 1.58 | 0.97 | 0.60–0.96 K |
| Both held, 0.3 and 0.03 ppm | 2.86 | 1.75 | 1.09–1.75 K |
| Both at Earth's 1750 levels | 6.35 | 3.89 | 2.4–3.9 K |
| Both at Earth's 2019 levels | 7.7 | 4.7 | 2.9–4.7 K |
| Both at the central steady state (21.6–29 and 0.55 ppm) | 12.9–13.5 | 7.9–8.3 | 4.9–8.3 K |
| Both at the highest (230 and 4.2 ppm) | 23.4 | 14.4 | 8.9–14.3 K |

The warming takes provisioning's response of 1.0–1.6 W/m²/K on the all-sky estimates. The column's own clear-sky
response near the design case, 1.38 W/m²/K with the water-vapour feedback and the dimmed titania sunlight (the
inverse-climate rows at 290 and 300 K), gives 9.3–9.8 K for the central steady state from the clear-sky forcing.

- **A part per million does more on the Moon.** Each is about a seven times thicker absorber over each square metre.
  At Earth's 2019 levels methane forces 4.0 W/m² on the Moon's clear sky against 2.5 on the Earth control, and N2O
  5.1 against 2.7. The bands saturate, so the forcing grows nearly with the logarithm of the amount: over a hundred
  times Earth's present methane (230 ppm) forces about four times as much.
- **The planned climate carries neither gas.** The design GCM run's radiation holds CO2 as its only trace gas
  ([climate/gcm](../../climate/gcm/README.md)), so whatever methane and N2O the Moon keeps adds to it. Earth's present
  levels would add 4.7 W/m² all-sky, and the screen's central steady state 7.9–8.3 W/m², within the range of the
  ring fleet's infrared glow (5.7–9.9 W/m²).
- **N2O is the larger lever.** It forces more per ppm than methane, and with no sink it keeps building: 9.3 W/m² at
  the 4.2 ppm ten thousand years bring. Methane levels off where the soils take what is made.
- **The two overlap.** Their bands share the region near 1,300 cm⁻¹, so together at the central steady state they
  force 13.5 W/m² clear-sky against 16.8 for the two apart.
- **Below Earth's 1750 air the forcing falls fast.** The bands come out of saturation, and the forcing falls about as
  the square root of the amount: a tenth of Earth's present N2O, 0.03 ppm, forces 0.83 W/m² all-sky against 3.35, and
  a sixth of its methane, 0.3 ppm, 1.05 against 2.32. The author's rule of 12 October holds the gases at whichever
  levels keep things closest to the already planned world, and the planned climate carries neither, so the domain
  synthesis holds them at that pair: 1.75 W/m², answered by 0.60% added dimming
  ([domain synthesis](../../research/studies/domain_synthesis/README.md#methane-and-n2o)).
- **What it leaves out.** Clouds enter only through AR6's Earth ratio; there is no stratospheric adjustment, no
  absorption of sunlight by the gases in the Moon's own columns and no chemistry, such as the stratospheric water that
  methane makes where OH exists.

## What it is and is not

It is a radiation calculation for columns whose temperature and humidity are
prescribed. It is not a solved climate, a cloud model or a habitability result.
Clouds, ozone, day and night, the poles and circulation are all outside it.
Comparing ASR with OLR tells whether a clear column at a given surface
temperature would warm or cool, and how strongly lunar gravity changes that
compared with Earth.
