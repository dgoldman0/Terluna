# Middle atmosphere

The temperature, ozone and surface ultraviolet of a column when radiation,
photochemistry and sunlight are mutually consistent, for the Earth control and
for the Moon behind each candidate spectral shield. It supplies the base of the
thermal column (`../thermal_column.py`), which carries the calculation up to the
exobase. Equations, data and limits are in [METHODS.md](METHODS.md).

| File | Holds |
|---|---|
| [photolysis.py](photolysis.py) | Ultraviolet–visible sunlight in the column (202–850 nm): actinic fluxes, photolysis rates, heating, surface UV index |
| [chemistry.py](chemistry.py) | O–H–N photochemistry (14 species, JPL 19-5 rates) with eddy diffusion, solved to steady state |
| [equilibrium.py](equilibrium.py) | Radiative–convective equilibrium at a fixed surface temperature, iterated with the chemistry, optionally with CO2 cooling out of local thermodynamic equilibrium; linear day–night response |
| [run.py](run.py) | The Earth control and the lunar cases; writes `results/` as a data product |
| [balance.py](balance.py) | Balance surface temperatures solved from each case's runs at several surface temperatures; writes `results/balance.json` |
| [escape.py](escape.py) | Exobase temperature and molecular loss with each case's upper temperature as the thermal column's base; heating by ultraviolet leaking through a filter; bounds on atomic-oxygen loss |
| [lbl_check.py](lbl_check.py) | Re-solves a finished case's upper air with line-by-line CO2 cooling, to measure the correlated-k error there |
| [benchmarks.py](benchmarks.py) | Line-by-line thermal and solar fluxes on the equilibrium profiles, for checking a 3-D model's radiation; writes `results/radiation_benchmarks.json` |
| [limb_heat.py](limb_heat.py) | Short-wave sunlight that passes the solar shield's window, annulus films and gaps, traced through the air along slant paths above the limb: where it heats, the state the heated air settles to, its loss, and the gap transmission each budget allows; writes `results/limb_heat.json` |
| [inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py) | Ultraviolet cross-sections and quantum yields (JPL recommendations via the MPI-Mainz atlas), hash-checked, not redistributed |
| [limb_inputs.json](limb_inputs.json), [fetch_limb_inputs.py](fetch_limb_inputs.py) | X-ray and extreme-ultraviolet absorption of N2, O2 and O (CXRO atomic scattering factors, Leiden photoabsorption database) and FISM2 daily X-ray spectra for the solar cycle, hash-checked, not redistributed |

## Run

From the repository root, after the radiative–convective inputs are restored:

```sh
python -m atmosphere.middle_atmosphere.fetch_inputs --download
python -c "from atmosphere.radiative_convective import ck; ck.build()"   # correlated-k tables, once (hours)
python -m atmosphere.middle_atmosphere.run --list
python -m atmosphere.middle_atmosphere.run
python -m atmosphere.middle_atmosphere.balance
python -m atmosphere.middle_atmosphere.lbl_check moon_1.2atm_edge_200nm moon_1.2atm_titania_stack
python -m atmosphere.middle_atmosphere.escape
python -m atmosphere.middle_atmosphere.benchmarks
python -m atmosphere.middle_atmosphere.fetch_limb_inputs --download
OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.limb_heat   # after atmosphere.loss_response.tides
```

Results are written in [results/](results/): `middle_atmosphere.csv` (one row per
case) with `middle_atmosphere.json` (schema, producer hashes, evidence),
`profiles/<case>.csv`, `balance.json`, `lbl_check.json`, `escape_coupling.csv`,
`radiation_benchmarks.json` and `limb_heat.json`. A case takes 7–15 minutes on four cores, and
10–20 minutes on three with non-equilibrium CO2 cooling. The limb tracing takes
about 50 CPU minutes and keeps each case's heat tables in the ignored `cache/limb_heat/`.

## Results (2026-09-24)

Every case holds a 288 K surface under clear skies. The settings are 400 ppm
CO2, Manabe–Wetherald humidity, 330 ppb N2O and 0.53 ppm H2 at the ground, and
lunar eddy mixing 36 times Earth's, unless a row says otherwise. The balance
column assumes a cloud effect of −20 W/m², about Earth's. It is solved: the
same columns were run at other surface temperatures (268–308 K for the Moon's
four base cases, 278 and 298 K for the rest) and the surplus interpolated
between them (`balance.json`).

| Case | Ozone (DU) | Tropopause | Stratopause | T at 1 Pa (LTE) | Subsolar UV index | Surplus at 288 K (W/m²) | Balance, C = −20 |
|---|---|---|---|---|---|---|---|
| Earth control | 277 | 215 K, 25.7 kPa | 256 K, 82 Pa | 198 K | 15.0 | 24.0 | 289.6 K |
| Moon 1.2 atm, titania stack | 0 | 179 K, 17.3 kPa | none | 169 K | 0.1 | 14.8 | 284.6 K |
| Moon 1.0 atm, titania stack | 0 | 183 K, 14.8 kPa | none | 170 K | 0.1 | 13.8 | 284.1 K |
| Moon 1.2 atm, 200-nm edge | 298 | 194 K, 22.8 kPa | 264 K, 29 Pa | 225 K | 2.4 | 22.8 | 289.9 K |
| Moon 1.0 atm, 200-nm edge | 285 | 198 K, 19.5 kPa | 265 K, 31 Pa | 224 K | 3.2 | 22.7 | 289.5 K |
| Moon 1.2 atm, 220-nm edge | 136 | 194 K, 22.8 kPa | 247 K, 29 Pa | 211 K | 4.4 | 22.6 | 289.8 K |
| Moon 1.2 atm, 230-nm edge | 96 | 194 K, 22.8 kPa | 234 K, 39 Pa | 198 K | 5.7 | 22.4 | 289.6 K |
| Moon 1.2 atm, 240-nm edge | 69 | 194 K, 22.8 kPa | 209 K, 90 Pa | 182 K | 7.2 | 22.3 | 289.5 K |
| Moon 1.2 atm, 310-nm edge | 138 | 194 K, 22.8 kPa | none | 175 K | 1.1 | 21.5 | 289.0 K |
| Moon 1.2 atm, unfiltered sunlight | 317 | 194 K, 22.8 kPa | 265 K, 29 Pa | 226 K | 2.3 | 32.1 | 296.3 K |
| 1.2 atm, 200-nm edge, no N2O | 558 | 202 K, 26.3 kPa | 264 K, 39 Pa | 220 K | 1.4 | 25.2 | 291.5 K |
| 1.0 atm, 200-nm edge, no N2O | 494 | 206 K, 22.4 kPa | 264 K, 36 Pa | 219 K | 2.0 | 25.3 | 290.9 K |
| 1.2 atm, 200-nm edge, Earth-rate mixing | 154 | 194 K, 22.8 kPa | 263 K, 29 Pa | 221 K | 4.9 | 21.9 | 289.3 K |
| 1.0 atm, 200-nm edge, Earth-rate mixing | 141 | 191 K, 17.0 kPa | 266 K, 24 Pa | 226 K | 6.9 | 20.3 | 288.2 K |
| 1.2 atm, 200-nm edge, 5 ppm stratospheric water | 393 | 210 K, 30.2 kPa | 256 K, 39 Pa | 213 K | 1.7 | 30.1 | 294.0 K |

Balance temperatures, solved, for a range of cloud effects C (W/m²):

| Case | C = 0 | C = −10 | C = −20 | C = −30 | C = −40 | Net budget change per K |
|---|---|---|---|---|---|---|
| Earth control | 297.4 K | 293.5 K | 289.6 K | 285.5 K | 281.3 K | 2.46 W/m² |
| Moon 1.2 atm, titania stack | 297.9 K | 291.2 K | 284.6 K | 278.1 K | 271.9 K | 1.61 W/m² |
| Moon 1.0 atm, titania stack | 296.6 K | 290.4 K | 284.1 K | 277.9 K | 272.0 K | 1.72 W/m² |
| Moon 1.2 atm, 200-nm edge | 302.3 K | 296.7 K | 289.9 K | 283.4 K | 277.2 K | 1.64 W/m² |
| Moon 1.0 atm, 200-nm edge | 300.7 K | 295.3 K | 289.5 K | 283.5 K | 277.4 K | 1.72 W/m² |

The variants were run only at 278 and 298 K, which brackets their balance for
three or four of these cloud effects. Their budgets change by 1.4–1.7 W/m²
per K, as their base cases' do; every bracketed balance is in `balance.json`.

What the calculations show:

- **Behind the titania stack the Moon has no ozone and a cold, dry upper
  atmosphere.** No light below about 340 nm reaches the air, so neither O2 nor
  O3 is split. Above an 86-km, 179 K tropopause the air stays at 150–177 K,
  and the stratosphere holds only 24 ppb of water. The surface UV index is 0.1:
  essentially no UV-B, including the UV-B skin needs to make vitamin D.
- **A shield that passes 200–242 nm grows an Earth-like ozone column, much
  higher up.** The 200-nm edge gives 298 DU at 1.2 atm and 285 DU at 1.0 atm.
  The six-fold deeper O2 column absorbs the ozone-making light at lower
  pressures. The ozone number density peaks near 5 kPa (130 km) and the mixing
  ratio near 100 Pa (320 km), and a 264 K stratopause sits near 30 Pa (410 km).
- **Surface ultraviolet stays low even with that ozone.** The subsolar UV index
  is 2.4, against 15 for the Earth control. At 1.2 atm the air's Rayleigh
  optical depth is about 7 at 310 nm, so most UV-B is scattered back out before
  it reaches the ground.
- **Where the filter's edge sits matters twice.** Ozone falls from 298 DU (edge
  at 200 nm) to 69 DU (240 nm). Edges between 220 and 240 nm also block the
  light that destroys N2O harmlessly, so more N2O meets O(¹D) and becomes NOx
  (up to 118 ppb, against 37 ppb for the 200-nm edge), which destroys more
  ozone. The 310-nm edge still builds 139 DU from its assumed 0.1% leakage
  below the edge. It also blocks ozone's own Hartley-band photolysis, so that
  ozone lives long; the number rests on the leakage.
- **Biology and mixing set the ozone as much as the shield does.** Without N2O
  there is no NOx and the column reaches 558 DU (494 DU at 1.0 atm). With
  Earth-rate eddy mixing (slower per scale height on the Moon) it is 154 DU
  (141 DU at 1.0 atm). With 5 ppm of
  stratospheric water it is 393 DU: more HOx ties NOx up as HNO3. Behind any of
  these shields nothing destroys NOy aloft (NO photolysis needs light below
  191 nm), so it builds to 30–120 ppb.
- **The self-consistent upper atmosphere warms the climate estimate.** Behind
  the titania stack, the cold stratosphere lowers outgoing longwave by about
  9 W/m² compared with the fixed 200 K stratosphere of the radiative–convective
  sweep. With Earth-like clouds the balance is then near 284 K rather than
  277 K. Every ozone-forming shield balances at 289.0–289.9 K, and the same
  method gives the Earth control 289.6 K (Earth: 288 K). The stratosphere's
  chemistry moves the balance more than the filter's edge does: 5 ppm of
  stratospheric water gives 294.0 K and no N2O 291.5 K. Unfiltered sunlight
  gives 296.3 K. All of this still hangs on the cloud effect, which only a 3-D
  model can supply.
- **Ozone falls as the surface warms.** Behind the 200-nm edge at 1.2 atm the
  column drops from 423 DU with a 268 K surface to 234 DU at 308 K (1.0 atm:
  351 to 216 DU). A warmer surface raises the tropopause and carries more water
  into the stratosphere, so the subsolar UV index rises from 1.5 to 3.3.
- **The month-long day would swing the ozone-heated upper air by up to about
  60 K from radiation alone.** The linear response at the equator, with no
  circulation, is about 30 K at 100 Pa and at 1 Pa and peaks near 40 Pa
  behind the 200-nm edge, where radiative relaxation takes 4–12 days. Behind the
  titania stack the swing is 1–3 K below 10 Pa and about 10 K near 1 Pa. The
  same method gives the Earth control a 4.7 K daily swing near its stratopause,
  comparable to the few-kelvin diurnal variation observed there. Swings this
  large drive strong tides and day–night winds.
- **The whole air column barely cools over the long night.** Its heat
  capacity, 8×10⁷ J/m²/K (a 19 m layer of water), limits cooling by the mean
  outgoing longwave to about 3.5 K over a 15-day night (4.2 K at 1.0 atm). The
  ground and the air just above it are followed through the month-long day in
  [climate/lunar_day.py](../../climate/lunar_day.py) (see the
  [climate README](../../climate/README.md)).

### Upper air out of local thermodynamic equilibrium

Above a few pascals, collisions are too rare to keep CO2's 15-µm band in local
thermodynamic equilibrium (LTE), so the cases above overstate its cooling
there. The `_nonlte` cases treat the band as a two-level system above 50 Pa
(METHODS.md). One question stays open: how much of the near-infrared sunlight
that CO2 absorbs up there becomes heat. Either all of it does (upper bound), or
only the part that collisions take before the molecules re-emit it
(collisional, lower bound).

| Case | 10 Pa | 3 Pa | 1 Pa | 0.3 Pa | 0.1 Pa (model top) |
|---|---|---|---|---|---|
| Earth control, LTE | 234 K | 210 K | 198 K | 192 K | 199 K |
| Earth control, non-LTE, collisional | 233 K | 207 K | 198 K | 194 K | 199 K |
| Earth control, non-LTE, all near-infrared heats | 235 K | 214 K | 211 K | 228 K | 256 K |
| US Standard Atmosphere 1976 | 232 K | 212 K | 198 K | 188 K | 188 K |
| Moon 1.2 atm, titania stack, LTE | 163 K | 170 K | 169 K | 163 K | 168 K |
| Moon 1.2 atm, titania stack, collisional | 156 K | 155 K | 146 K | 142 K | 158 K |
| Moon 1.2 atm, titania stack, all heats | 163 K | 172 K | 175 K | 192 K | 243 K |
| Moon 1.2 atm, 200-nm edge, LTE | 253 K | 235 K | 225 K | 200 K | 186 K |
| Moon 1.2 atm, 200-nm edge, collisional | 253 K | 234 K | 222 K | 200 K | 187 K |
| Moon 1.2 atm, 200-nm edge, all heats | 254 K | 237 K | 229 K | 215 K | 239 K |

- **The Earth control favours the collisional bound.** With it the column is
  within 6 K of the US Standard Atmosphere from 10 to 0.3 Pa (about 65–88 km).
  Letting all the absorbed near-infrared heat the air makes 0.3 Pa 40 K too
  warm. The real mesosphere is also shaped by waves and circulation that this
  column lacks, so the agreement supports the bound; it does not validate it.
- **Behind the titania stack the upper air is colder still.** With the
  collisional bound, 0.3 Pa is 142 K, 21 K below the LTE estimate. The upper
  bound is 192 K.
- **Behind the 200-nm edge little changes.** Ozone heating dominates there: the
  collisional bound matches LTE within 3 K, and the upper bound is 15 K warmer at
  0.3 Pa.
- The 1.0 atm columns agree with these within 2 K. The top-of-atmosphere budget
  moves by less than 0.1 W/m², so the surface results above stand.
- The model top at 0.1 Pa absorbs the line-core sunlight meant for everything
  above it, and there correlated-k and line-by-line differ by 30–40 K
  (`lbl_check.json`). Treat 0.1 Pa temperatures as bounds, not values.

### Exobase

The thermal column (`../thermal_column.py`) takes each case's temperature at
0.3 Pa as its base; correlated-k and line-by-line agree there within about 1 K.
Heat deposited above the base is swept. For a filter that passes 0.1% of
sunlight below its edge, the heat is estimated from the WHI 2008 solar spectrum
below 175 nm: 1.7×10⁻⁶ W/m² for the quiet Sun and 4.2×10⁻⁶ W/m² near solar
maximum (METHODS.md). For the titania stack the film's own transmission is
computed from published optical constants
([protection/spectra/stack_short_wave.json](../../protection/spectra/stack_short_wave.json)),
so its 0.1% rows stand for light that bypasses the film through gaps in the
aperture.

| 1.2 atm, upper air | Base | Exobase at 0 / 10⁻⁶ / 3×10⁻⁶ / 10⁻⁵ W/m² | 0.1% leak: quiet Sun / solar maximum | Film alone: quiet Sun / solar maximum | Loss at solar maximum: N2 and O2 + atomic O |
|---|---|---|---|---|---|
| Titania stack, collisional | 142 K | 142 / 157 / 188 / 299 K | 167 / 207 K | 142 / 142 K | 3×10⁻⁵ kg/s + none |
| Titania stack, LTE | 163 K | 163 / 179 / 209 / 316 K | 189 / 229 K | 163 / 164 K | 3×10⁻³ kg/s + none |
| Titania stack, all heats | 192 K | 192 / 207 / 239 / 313 K | 218 / 259 K | 192 / 192 K | 0.5 kg/s + none |
| Titania stack, all heats, 280 ppm CO2 | 197 K | 197 / 212 / 244 / 310 K | 223 / 264 K | 197 / 197 K | 1.1 kg/s + none |
| Titania stack, all heats, 150 ppm CO2 | 206 K | 206 / 222 / 253 / 304 K | 233 / 271 K | 206 / 206 K | 3.9 kg/s + none |
| 200-nm edge, collisional | 200 K | 200 / 216 / 247 / 308 K | 226 / 266 K | | 1.7 + 0.8 kg/s |
| 200-nm edge, LTE | 200 K | 200 / 216 / 248 / 308 K | 227 / 267 K | | 1.8 + 0.8 kg/s |
| 200-nm edge, all heats | 215 K | 215 / 231 / 261 / 296 K | 242 / 274 K | | 10 + 2.8 kg/s |
| Earlier assumption, 180 K at 0.1 Pa | 180 K | 180 / 194 / 224 / 317 K | | | |

- **Blocking 99.9% below 175 nm keeps the exobase near the 250 K target, but
  not through solar maximum.** Behind the 200-nm edge (collisional bound) the
  exobase is 226 K for the quiet Sun and 266 K near solar maximum.
  Interpolating, blocking 99.95% would hold it near 233 K at solar maximum
  (248 K for the upper bound). Behind the titania stack the same 0.1% leak
  gives 167–207 K.
- **The titania stack's film lets almost no heat through.**
  Its 10 µm of silica and 1 µm of titania pass only hard X-rays shorter than
  about 1.5 nm. From 5 nm to 200 nm they transmit less than 10⁻²⁰. The
  heat that reaches the upper air is at most 1×10⁻⁹ W/m² for the quiet Sun
  and 2×10⁻⁸ W/m² at solar maximum, where FISM2's measured spectra raise
  these X-rays about twentyfold (`escape.xray_cycle`). That is 1,700 and 220
  times less than a 0.1% leak, so the film alone holds the exobase within half
  a kelvin of its base, at 142–192 K. Beyond
  that, light passing gaps, pinholes and edges of the aperture warms it (a
  design number, not a material one), and so does the sky's own Lyman-alpha.
- **The sky's Lyman-alpha is not in these tables; the loss response adds it.**
  Interplanetary hydrogen glows at about 1,000 rayleigh at 1 AU
  (protection/report.md, section 9), about 4×10⁻⁶ W/m² onto a surface facing
  open sky, and reaches the upper air from every direction, so no Sun-facing
  shield blocks it. [atmosphere/loss_response](../loss_response/README.md)
  deposits it above the base (O2 absorption over the column above 0.3 Pa, taken
  with gravity at the base radius, for light from the whole sky): 2.1–2.8×10⁻⁶
  W/m² at quiet Sun. That raises the exobase by 31–42 K at quiet Sun and 47–61 K
  at solar maximum (a solar-cycle factor of 1.5 assumed), for every shield.
- **Some loss channels are handled downstream or not yet included.** Earth's
  tide lowering the escape barrier is in the
  [loss response](../loss_response/README.md): test molecules in the Earth–Moon
  three-body problem give 2.6 times the nitrogen Jeans flux at the September
  feasibility report's reference exobase, where the report's uniform barrier
  factor gives about 4.5
  ([research/baselines/feasibility/report.md](../../research/baselines/feasibility/report.md),
  sections 4–5). Solar-wind sputtering and ion pickup are screened there too;
  the report gives 0.4–40 kg/s over a three-lunar-radius cross-section, as a
  budget translation. Hydrogen from water is not yet included.
- **Less CO2 warms the upper air, within the film's margin.** CO2 is the upper
  air's main coolant. In the warmest treatment (all absorbed near-infrared
  heats), cutting it from 400 to 280 ppm raises the base by 5 K and to 150 ppm
  by 15 K, and the exobase follows one for one: behind the film alone it stays
  at 197–206 K. Only the 0.1% leak at solar maximum passes 250 K, as it already
  does at 400 ppm. 150 ppm (18 Pa) lies below the plant-growth floor
  ([research/studies/atmospheric_co2](../../research/studies/atmospheric_co2/)), so any CO2 level
  biology accepts is safe here.
- **The losses stay small either way.** At 267 K, 2.6 kg/s removes about 0.3% of
  an atmosphere of 2.8×10¹⁸ kg in 100 million years.
- **Atomic oxygen adds about half again behind the 200-nm edge.** O2 photolysis
  there leaves 200 ppm of O at the base. With the chemistry's eddy mixing
  continued upward, the homopause is at 4–9×10⁻⁵ Pa, O makes up about 1% of
  the gas at the exobase, and its escape is about half the molecular loss. The
  bounds at solar maximum are wide: 0.02 kg/s if eddy mixing reached the
  exobase, 33 kg/s if the oxygen separated from the base up. The escaping oxygen's 63-µm cooling, which
  would lower the exobase, is not included. Behind the titania stack there is
  no atomic oxygen.
- **At these leaks the base carries through almost one for one.** Across the
  three treatments the base spans 142–215 K and the exobase at the
  solar-maximum leak 207–274 K, so the upper air matters about as much as the
  leak. Only at 10⁻⁵ W/m² does the deposited heat dominate (296–316 K).
- All of these are global means with no infrared cooling in the thermal column.
  The day–night swing of the upper air (above) is tens of kelvin, and the
  circulation it drives is not included.

### Slant paths above the limb (2026-10-07)

The heat counts above put a quarter of the light that reaches the disk above
the base. [limb_heat.py](limb_heat.py) traces the sunlight that passes the
solar shield along its actual paths through the air, for the ring fleet's
aperture: a screen at 20,000 km covering four lunar radii, with the window
stack over the disk, a light film over the annulus, gaps that pass a grey share
of the band over the whole aperture, and open sky beyond it. It takes the loss
response's six cases, Lyman-α glow and Earth's tide, finds the state the heated
air settles to and the gap transmission each budget allows, and writes
[results/limb_heat.json](results/limb_heat.json). It keeps each ray's deposit,
so the product also gives the heat for protected radii of 2 to 10 lunar radii,
which the loss response reads ([traced.py](../loss_response/traced.py)): on
7 October the author made this traced count O1's and took FISM2's measured rise
of the X-rays at solar maximum, which these tables use. Methods and limits are
in [METHODS.md](METHODS.md).

- **Gaps heat the thermosphere 2.5–4.2 times the escape model's count.** A
  grey transmission through gaps over the aperture lights the air above the
  limb as well as the disk, and the tall air takes the extreme ultraviolet along slant paths out
  to near the exobase. Around the quiescent air the factor is 2.5 behind the
  titania stack with collisional upper air, 2.7 in LTE, 3.2 with all
  near-infrared heating and 3.8–4.2 behind the 200-nm edge. As the heat
  swells the air it grows, to 3.1–3.9 at the heats where the cooler titania
  cases lose 1–100 kg/s. On rays tangent between the base and the exobase,
  [absorption.py](../loss_response/absorption.py), with band-mean cross
  sections, finds the same heating within 5%.
- **The films pass X-rays near 1 nm, which the thermosphere stops in a ring
  just above the limb.** The 4 µm film passes 0.75–1.35 nm, the 2 µm film
  0.85–2.55 nm and the window stack 0.55–0.95 nm. Absorbing them takes a
  column near 10²³ m⁻², so only rays tangent between the limb and about 1.5
  lunar radii give them up: 7–17% of what passes the annulus, none of it below
  the base. Over the annulus the 4 µm film heats the thermosphere 1.1–2.5 times
  what the escape model's film count gives it over the disk.
- **The films' X-rays rise about twentyfold from quiet Sun to solar maximum.**
  The escape model used to multiply every X-ray below 10 nm by 100 at solar
  maximum. FISM2's daily spectra, the year around each of the last three maxima
  over the WHI 2008 quiet week, give by band 53 below 0.5 nm falling to 3.6 at
  5–10 nm. Weighted by the light each passes, the mean of the three maxima is 18
  for the 4 µm film (13–21 across them), 16 for the 2 µm film, 20 for the window
  stack and 6.0 for the whole band below 10 nm (4.6–6.9), which sets a gap's
  X-rays. Single days at maximum run up to about twice the yearly mean.
- **With those factors the 4 µm film adds a tenth to a fifth of the glow's
  heat.** At solar maximum it heats the thermosphere by 3.5–5.2×10⁻⁷ W/m²
  behind the titania stack (2.5–6.0×10⁻⁷ across the three maxima), 11–14% of
  the 3.1–3.8×10⁻⁶ W/m² the glow gives there, and 7–8×10⁻⁷ W/m² behind the
  200-nm edge. With a factor of 100 it would add 2.0–2.9×10⁻⁶ W/m², about as
  much as the glow. The 2 µm film adds 1.3–1.9×10⁻⁶ W/m², and the window stack
  over the annulus 2–3×10⁻⁸ W/m².

The states with the 4 µm film and the ring fleet's 4-radius shadow, as
molecular loss with Earth's tide (kg/s), against the escape model's disk count
at the same transmission. Solar maximum is the mean of the three maxima, with
the range across them.

| 1.2 atm, titania stack | Quiet Sun, standard level (2×10⁻⁴): traced / disk count | Solar maximum, no gaps | Tight level (3×10⁻⁵) | Standard level | Disk count, standard level |
|---|---|---|---|---|---|
| Collisional upper air | 4×10⁻⁶ / 8×10⁻⁷ | 2×10⁻⁵ (1–2×10⁻⁵) | 4×10⁻⁵ (3–5×10⁻⁵) | 0.003 (0.0025–0.004) | 6×10⁻⁵ |
| LTE | 0.0013 / 0.0003 | 0.004 (0.003–0.005) | 0.009 (0.007–0.011) | 0.42 (0.32–0.49) | 0.009 |
| All near-infrared heats | 0.61 / 0.16 | 1.5 (1.2–1.7) | 2.8 (2.2–3.1) | 33 (29–36) | 2.1 |

- **The cooler titania cases stay small; the warmest swells.** At the
  standard level and solar maximum the exobase sits at 2.2 lunar radii
  (collisional, 230–234 K) and 2.7 (LTE, 264–269 K), and with all
  near-infrared heating at 3.8–3.9 (300–302 K). Behind the 200-nm edge the
  standard level gives 3.5–3.7 kg/s at quiet Sun with collisional or LTE upper
  air. At solar maximum there, and with all near-infrared heating at quiet Sun
  as well, the heat the swollen air takes stays above the heat that swells it
  up to the tables' top, where the exobase passes six lunar radii. Beyond four
  lunar radii the thermosphere reaches unfiltered light outside the aperture.
- **The gaps a budget allows fall to a quarter or a third.** For the same
  allowance of molecular loss, traced gaps may pass 0.22–0.32 of what the disk
  count allows in the cooler titania cases, at quiet Sun and at solar maximum
  alike, and 0.12–0.27 with all near-infrared heating. For 1 kg/s at solar
  maximum: 4.7×10⁻⁴ (collisional) and 2.4×10⁻⁴ (LTE), against the disk count's
  1.5×10⁻³ and 9×10⁻⁴; with all near-infrared heating, the film and glow alone
  exceed 1 kg/s, and 10 kg/s allows 1.0×10⁻⁴. These
  allowances leave out the solar wind and the exosphere step's losses, so they
  compare the two counts at the same allowance.
- **The tables keep the light by band and by height.** Fifteen bands of the
  light below 175 nm are traced apart, so any solar spectrum weights them, and
  each part of the light's heat is followed in log pressure above the base. At
  the last three maxima FISM2's daily record puts the light from 10 to 175 nm at
  1.25-1.41 times the quiet week, weighted by energy, where the escape model
  takes 2.5: the extreme ultraviolet below 121 nm rises 1.6-1.9, Lyman-α 1.3-1.6
  and the far ultraviolet from 122 to 175 nm, which carries most of the energy,
  1.07-1.18 ([loss_response/cycle.py](../loss_response/cycle.py)). Through gaps
  the far ultraviolet and Lyman-α heat the lower thermosphere and the extreme
  ultraviolet its upper half.
- **Where the heat lies matters as much as how much.** The thermal column
  spreads the heat with a fixed shape in log pressure. The traced light lands
  at a mean depth 0.38–0.43 of the way from the base to the exobase, at about
  the depth of the disk's own light, between the column's 'low' shape (0.25)
  and its 'middle' shape (0.5), which the loss response uses and the table
  above follows. With 'low' the losses fall 14–550 times. Both counts pass
  through the same column, so the shape moves their losses together.
