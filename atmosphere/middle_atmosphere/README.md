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
| [equilibrium.py](equilibrium.py) | Radiative–convective equilibrium at a fixed surface temperature, iterated with the chemistry; linear day–night response |
| [run.py](run.py) | The Earth control and the lunar cases; writes `results/` as a data product |
| [escape.py](escape.py) | Exobase temperature and molecular loss with each case's upper temperature as the thermal column's base |
| [lbl_check.py](lbl_check.py) | Re-solves a finished case's upper air with line-by-line CO2 cooling, to measure the correlated-k error there |
| [inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py) | Ultraviolet cross-sections and quantum yields (JPL recommendations via the MPI-Mainz atlas), hash-checked, not redistributed |

## Run

From the repository root, after the radiative–convective inputs are restored:

```sh
python -m atmosphere.middle_atmosphere.fetch_inputs --download
python -c "from atmosphere.radiative_convective import ck; ck.build()"   # correlated-k tables, once (hours)
python -m atmosphere.middle_atmosphere.run --list
python -m atmosphere.middle_atmosphere.run
python -m atmosphere.middle_atmosphere.lbl_check moon_1.2atm_edge_200nm moon_1.2atm_titania_stack
python -m atmosphere.middle_atmosphere.escape
```

Results are written in [results/](results/): `middle_atmosphere.csv` (one row per
case) with `middle_atmosphere.json` (schema, producer hashes, evidence),
`profiles/<case>.csv`, `lbl_check.json` and `escape_coupling.csv`. A case takes
7–15 minutes on four cores.

## Results (2026-09-24)

Every case holds a 288 K surface under clear skies. The settings are 400 ppm
CO2, Manabe–Wetherald humidity, 330 ppb N2O and 0.53 ppm H2 at the ground, and
lunar eddy mixing 36 times Earth's, unless a row says otherwise. The balance
temperature linearises the top-of-atmosphere surplus with the
radiative–convective model's slope near 288 K (1.29–1.40 W/m²/K for the Moon,
2.34 for Earth) and an assumed cloud effect of −20 W/m², about Earth's.

| Case | Ozone (DU) | Tropopause | Stratopause | T at 1 Pa | Subsolar UV index | Surplus at 288 K (W/m²) | Balance, C = −20 |
|---|---|---|---|---|---|---|---|
| Earth control | 277 | 215 K, 25.7 kPa | 256 K, 82 Pa | 198 K | 15.0 | 24.0 | 289.7 K |
| Moon 1.2 atm, titania stack | 0 | 179 K, 17.3 kPa | none; 150–177 K above | 169 K | 0.1 | 14.8 | 284.0 K |
| Moon 1.0 atm, titania stack | 0 | 183 K, 14.8 kPa | none | 170 K | 0.1 | 13.8 | 283.6 K |
| Moon 1.2 atm, 200-nm edge | 298 | 194 K, 22.8 kPa | 264 K, 29 Pa | 225 K | 2.4 | 22.8 | 290.2 K |
| Moon 1.0 atm, 200-nm edge | 285 | 198 K, 19.5 kPa | 265 K, 31 Pa | 224 K | 3.2 | 22.7 | 289.9 K |
| Moon 1.2 atm, 220-nm edge | 136 | 194 K | 247 K | 211 K | 4.4 | 22.6 | 290.0 K |
| Moon 1.2 atm, 230-nm edge | 96 | 194 K | 234 K | 198 K | 5.7 | 22.4 | 289.9 K |
| Moon 1.2 atm, 240-nm edge | 69 | 194 K | 209 K | 182 K | 7.2 | 22.3 | 289.8 K |
| Moon 1.2 atm, 310-nm edge | 139 | 194 K | none; 191 K near the tropopause | 175 K | 1.1 | 21.5 | 289.2 K |
| Moon 1.2 atm, unfiltered sunlight | 317 | 194 K | 265 K | 226 K | 2.3 | 32.1 | 297.3 K |
| 200-nm edge, no N2O | 558 | 202 K | 264 K | 220 K | 1.4 | 25.2 | 292.0 K |
| 200-nm edge, Earth-rate mixing | 154 | 194 K | 263 K | 221 K | 4.9 | 21.9 | 289.5 K |
| 200-nm edge, 5 ppm stratospheric water | 393 | 210 K, 30.2 kPa | 257 K | 213 K | 1.7 | 30.1 | 295.8 K |

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
  there is no NOx and the column reaches 558 DU. With Earth-rate eddy mixing
  (slower per scale height on the Moon) it is 154 DU. With 5 ppm of
  stratospheric water it is 393 DU: more HOx ties NOx up as HNO3. Behind any of
  these shields nothing destroys NOy aloft (NO photolysis needs light below
  191 nm), so it builds to 30–120 ppb.
- **The self-consistent upper atmosphere warms the climate estimate.** Behind
  the titania stack, the cold stratosphere lowers outgoing longwave by about
  9 W/m² compared with the fixed 200 K stratosphere of the radiative–convective
  sweep. With Earth-like clouds the balance is then near 284 K rather than
  277 K. Every ozone-forming shield balances near 290 K, and the same method
  gives the Earth control 289.7 K (Earth: 288 K). All of this still hangs on the
  cloud effect, which only a 3-D model can supply.
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
  outgoing longwave to about 3.5 K over a 15-day night (4.2 K at 1.0 atm). What
  the ground and the air just above it do needs a surface model or a GCM.

### Exobase

The thermal column (`../thermal_column.py`) takes each case's temperature at
0.3 Pa as its base; correlated-k and line-by-line agree there within about 1 K.
Heat deposited above the base by extreme-ultraviolet leakage is swept, as in
the earlier screens.

| Base (1.2 atm) | Base temperature | Exobase at 0 / 10⁻⁶ / 3×10⁻⁶ / 10⁻⁵ W/m² | Molecular loss at 3×10⁻⁶ W/m² |
|---|---|---|---|
| Titania stack | 163 K | 163 / 179 / 210 / 316 K | 2×10⁻⁴ kg/s |
| 200-nm edge | 200 K | 200 / 216 / 248 / 308 K | 0.3 kg/s |
| Earlier assumption, 180 K at 0.1 Pa | 180 K | 180 / 194 / 224 / 317 K | 5×10⁻³ kg/s |

- A shield that passes 200–242 nm makes the base of the thermosphere about 37 K
  warmer than the titania stack does. At low leakage that raises the exobase by
  20–40 K and molecular loss several hundred-fold. The molecular loss stays
  negligible for an atmosphere of 2.8×10¹⁸ kg.
- The same shield leaves atomic oxygen at 250 ppm at the model top and rising.
  The thermal column holds only N2 and O2, so the escape of those lighter atoms
  is the open question for this shield. Behind the titania stack there is none.
- The top layer at 0.1 Pa absorbs the line-core sunlight meant for everything
  above it. There, correlated-k and line-by-line differ by 30–40 K
  (`lbl_check.json`: 168 against 127 K behind titania, 186 against 155 K behind
  the 200-nm edge), and local thermodynamic equilibrium fails. Treat 0.1 Pa
  temperatures as bounds, not values.
