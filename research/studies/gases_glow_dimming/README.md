# The gases, the glow and the dimming in the GCM

13 October 2026. The author asked for the GCM to run again with the gases, the fleet's glow and the dimming that
answers them before the gas pair and the dimming's cost are decided. The light limitation study had found that in
Earth's models a dimming that holds the temperature also takes rain, and that most of the Moon's land grows with its
rain. Three runs branch from the design case's year-29 restart ([climate/gcm](../../../climate/gcm/README.md),
`A28_dim5_moon`) and fly 20 years each on the same model. Their years 10–19 are compared with the design case's
years 20–29.

| Run | The air | Infrared from space at the top of the air | Dimming added |
|---|---|---|---|
| `A28_gases` | the middle pair, carried as 755 ppm of CO2 | none | 1.69%, the gases' |
| `A28_gases_glow` | the same | the glow and the returned infrared: 10.17 W/m² by day, 0.52 by night | 4.22%: the gases' and 2.53% matching the infrared watt for watt |
| `A28_gases_glow_held` | the same | the same | 3.52%: the gases' and 1.83% |

`python -m research.studies.gases_glow_dimming.run` writes
[results/gases_glow_dimming.json](results/gases_glow_dimming.json) (schema `terluna.research.gases-glow-dimming/1`,
with its producer, the sha256 of its three inputs and of the climatologies, run records and extract it reads, its
evidence statement and reading rule) in about a second. It reads each run's climatology
(`python -m climate.gcm.climatology A28_gases:10-19` and so on, with the GCM environment) and an extract of the layers
and the first 3-day means (`climate/gcm/.venv/bin/python -m research.studies.gases_glow_dimming.extract`, seconds).
`python -m research.studies.gases_glow_dimming.glow_absorption` writes
[results/glow_absorption.json](results/glow_absorption.json), the spectral check below, in a few seconds once the
atmosphere domain's correlated-k tables are cached.

## What it finds

- **The gases' dimming holds the climate and takes little rain.** With the middle pair and its 1.69% the Moon stays
  at the design temperature, 0.14 K warmer, and its clouds hold. The rain falls 0.9% over land and sea alike, 1.7% on
  land within 15° of the equator and less poleward. That is 0.55% of the rain for each 1% of dimming, about half the
  0.76–1.38% of Earth's GeoMIP G1 models (Tilmes et al. 2013). The rain-fed land loses 0.12 Gt C a year, 0.5–0.7%
  of the living Moon's growth, beside the light study's 0.1–0.5% through light.
- **The middle pair stands.** Scaled by their warming, the pairs cost the growth 0.2–0.4% (low), 0.5–0.7% (middle)
  and 0.6–1.0% (Earth's present air) through the rain, and 0.06–0.24%, 0.12–0.50% and 0.16–0.66% through light. The
  low pair saves 0.3–0.6% of the growth in all for contactors through 1.0–1.4% of the air a year from the first
  century; Earth's air adds 0.2–0.4% with N2O decomposition about three times the middle pair's
  ([domain synthesis](../domain_synthesis/README.md#where-they-are-held)). The middle pair stays the balanced choice.
- **The glow heats the upper air.** About 4% of its infrared reaches the ground in the GCM, and 2.0–2.5% in the
  atmosphere domain's correlated-k model under clear sky. The films are cold, 150–190 K, so their glow falls in the
  far infrared and the CO2 band, which the moist air takes high up. The infrared they return goes back to the air that
  sent it. Three-quarters of the two is absorbed above half the surface pressure, about 34 km up.
- **Matched watt for watt, the glow cools the Moon.** With the 2.53% that matches its 5.34 W/m² in the energy balance,
  the ground ends 1.26 K cooler and is still cooling 0.16 K a decade. The rain falls 7.6%, and 11–14% on land at
  30–60°.
- **The glow needs 0.7 of its watts in dimming.** The two glow runs differ by 0.70% of dimming and 1.20 K, 1.72 K for
  each 1%. The dimming that holds the design temperature with the glow is 1.80% for its infrared, 0.71 of watt for
  watt, and 0.68 to hold the gases' run. With the 0.21% for the returned sunlight the glow takes about 2.0% of the
  design sunlight, where the glow study put 2.7% at an efficacy of 1. The gases and the glow together take 3.49% where
  watt for watt they took 4.22%.
- **Held, the glow still takes rain.** In the held run the ground stays within 0.06 K of the design case while the
  upper air warms 0.3 K about 55 km up, 0.7 K at 70 km, 1.2 K at 90 km and 2.1 K in the top layer. The warmer upper air
  steadies the column: the convective rain, nine-tenths of the Moon's, falls 3.9% while the large-scale rain holds. At
  the design temperature the rain falls 3.4%, 4.0% over land: 2.4% and 3.1% of it are the glow's, 0.46% of the rain
  for each W/m² of its infrared. The high cloud's cover falls by a fifth, from 14.3% to 11.3%, and the low cloud holds.
- **What it costs the living Moon.** At the design temperature the rain-fed land loses 0.55–0.58 Gt C a year,
  2.3–3.4% of the living Moon's growth, and the plan's 3.49% of dimming costs 0.24–1.03% more through light. The glow's
  part through the rain, 1.8–2.6%, is more than three times the gases'. Each W/m² taken off the glow's infrared returns
  about 0.46% of the rain; the glow study's further answers (the annulus absorber, the ultraviolet mirror, the wider
  gap) are the levers.
- **The plants make it up.** The author set the glow's further answers aside on 13 October, with the plants' water
  design to make up the rain the glow takes. Rain-fed growth is the water the plants pass times the carbon they fix
  per litre, so plants that fix 3.2% more carbon per litre, or take a 3.2% larger share of what the land evaporates,
  restore the glow's part, and 4.1% the whole plan's. The dry bands need more: 2.0% within 15° of the equator, 4.0% at
  15–30°, 5.6% at 30–45° and 11% at 45–60°. CAM plants on Earth fix 2–6 times as much carbon per litre as C3 plants
  (Hartzell et al. 2018), so the design has room to spare.

## The runs against the design case

| | Design case, years 20–29 | Gases | Glow, watt for watt | Glow, held |
|---|---|---|---|---|
| Dimming added | | 1.69% | 4.22% | 3.52% |
| Surface temperature | 294.90 K | +0.14 K | −1.26 K | −0.05 K |
| Its trend over years 10–19 | | +0.06 K a decade | −0.16 K a decade | −0.05 K a decade |
| Rain | 2.24 mm/day | −0.9% | −7.6% | −3.5% |
| Rain over land | 2.04 mm/day | −0.9% | −7.6% | −4.2% |
| Evaporation over land | 1.65 mm/day | −0.9% | −8.4% | −4.2% |
| Cloud cover; high deck; low deck | 23.5%; 14.3%; 5.3% | 23.5%; 14.4%; 5.3% | 20.1%; 11.6%; 5.4% | 19.6%; 11.3%; 5.4% |
| Water vapour | 176 kg/m² | +1.5% | −7.1% | +0.5% |
| Convective rain, nine-tenths of the rain | | −1.0% | −8.2% | −3.9% |
| Rain-fed growth | 14.0–14.6 Gt C a year | −0.12 | −1.14 to −1.20 | −0.57 to −0.60 |

| Land rain by band | 0–15° | 15–30° | 30–45° | 45–60° | 60–75° | 75–90° |
|---|---|---|---|---|---|---|
| Design case, mm/day | 4.78 | 2.70 | 1.30 | 0.19 | 0.05 | 0.04 |
| Gases | −1.7% | +0.1% | +0.3% | +1.1% | +3.3% | +1.3% |
| Glow, watt for watt | −7.0% | −6.5% | −11.3% | −14.0% | +4.5% | +4.6% |
| Glow, held | −3.9% | −3.5% | −5.9% | −9.2% | +4.1% | −3.6% |

| What the plants make up, by band | 0–15° | 15–30° | 30–45° | 45–60° | All rain-fed land |
|---|---|---|---|---|---|
| Evaporation over land, the glow's part | −2.0% | −3.9% | −5.3% | −9.8% | |
| More carbon per litre that restores it | +2.0% | +4.0% | +5.6% | +10.9% | +3.2% |
| The same for the whole plan, gases and glow | +3.4% | +4.1% | +6.4% | +10.5% | +4.1% |

The poleward bands' small rains change by a few per cent either way and carry little growth.

The cloud decks are the climatology's: the low deck the most cover among the layers below half the surface pressure,
the high deck above it. The poleward bands rain little, so their changes move little water. The held point, at 1.80%
for the glow, lies 0.03% beyond the held run and is found by carrying the two glow runs' difference that far.

| The air by layer (change from the design case) | 0.9 km | 10 km | 28 km | 40 km | 54 km | 71 km | 92 km | 129 km |
|---|---|---|---|---|---|---|---|---|
| Design case | 292.9 K | 285.6 K | 268.8 K | 255.9 K | 239.2 K | 219.2 K | 196.6 K | 170.0 K |
| Gases | +0.11 | +0.10 | +0.18 | +0.19 | +0.18 | +0.15 | +0.21 | −0.20 |
| Glow, watt for watt | −1.12 | −1.04 | −1.30 | −1.51 | −1.48 | −1.06 | −0.12 | +0.93 |
| Glow, held | −0.05 | 0.00 | +0.06 | +0.08 | +0.26 | +0.71 | +1.23 | +2.13 |

The heights are the design case's layers (sigma 0.983 to 0.038), from its layer temperatures at lunar gravity.

## Where the energy goes in the first month

Each branch against the one it differs from by that forcing alone, over its first ten 3-day means, before the climate
responds (W/m², downward positive):

| | Gases, against the design year on the patched build | Glow, watt for watt, against gases | Glow, held, against gases |
|---|---|---|---|
| Top of the air: sunlight; infrared; net | −3.45; +3.33; −0.12 | −5.22; +5.07; −0.15 | −3.76; +5.10; +1.34 |
| Ground: sunlight; infrared; net radiation | −2.04; +0.43; −1.61 | −2.98; +0.20; −2.77 | −2.20; +0.22; −1.99 |
| The air gains | 1.49 | 2.62 | 3.33 |

Every forcing here heats the air and takes radiation from the ground. Air heated directly needs less of the heat that
condensing rain releases, so the rain falls at once, before the temperature moves: Earth's models find the same fast
response to CO2 (Tilmes et al. 2013). The gases move about 1.5 W/m² that way and the glow about 2.6–3.3 more, which is
why the glow takes rain even with the temperature held.

## The spectral check

[glow_absorption.py](glow_absorption.py) takes the glow and the returned infrared through the middle atmosphere's
clear titania-stack columns at 288 and 298 K by the atmosphere domain's correlated-k model (16 bands, 16 g-points),
interpolated to 294.9 K. The glow is grey at the films' temperatures and arrives evenly from the sky; the returned
infrared comes back at the angle and g-point at which it left, as a specular film returns it.

| Infrared from space | Reaches the ground | Absorbed below 34 km | 34–62 km | Above 62 km |
|---|---|---|---|---|
| Glow, films at 150 K | 0.4% | 10% | 38% | 52% |
| Glow, films at 190 K | 1.2% | 17% | 38% | 44% |
| Glow, films at 300 K (the deepest stacks' inner plates) | 3.2% | 31% | 33% | 33% |
| The Moon's returned infrared | 5.2% | 44% | 33% | 17% |
| Both, as the run carries them (films at 150–190 K) | 2.0–2.5% | 21–26% | 36–37% | 36–40% |
| The GCM, all sky, from the first means | 3.8–4.1% | | | |

The boundaries are PlaSim's half levels at half and a quarter of the surface pressure. The two models agree that the
ground gets almost none of it.

## How it is set up

- **The gases as CO2.** PlaSim's longwave radiation has water vapour, CO2 and ozone only. A calibration branch with
  1,022 ppm in place of 400 gives 5.27 W/m² at the top of the air over the first 3-day mean and 5.30 over ten, 5.63
  W/m² for each e-fold, so the middle pair's 3.58 W/m² is 755 ppm (`--co2-ppm`).
- **The glow as infrared from space.** A patch to PlaSim's longwave scheme adds infrared to the downward flux at the
  top of its air, 10.17 W/m² where the Sun is up and 0.52 where it is down (`--top-infrared`), and PlaSim carries it
  down with its broadband transmissivities. Those are the recommended answer's glow, 3.60 W/m² (6.72 by day and 0.48
  by night), and returned infrared, 1.75 W/m² (3.46 and 0.04), from the
  [glow study](../solar_shield_array/glow.md), read from research/array-industry at a59f5a9. The returned sunlight,
  0.61 W/m², and the dimming that answers it cancel and are left out. At zero the patch changes nothing: the design
  year flown again on the patched build is identical to the unpatched one.
- **The dimming as less sunlight.** The design dimmer passes 0.95 of the titania stack's sunlight; the runs pass
  0.933907, 0.909867 and 0.916522.
- **The held run's dimming** came from the watt-for-watt run's cooling and the 1.9 K for each 1% of sunlight that the
  earlier dimmer runs gave ([climate/gcm](../../../climate/gcm/README.md), "A dimmer shield"). It landed 0.05 K below
  the design temperature.

## What it leaves out

- **One model.** T21 with 10 layers and a 50 m mixed-layer sea; the rain's response and the glow's efficacy are this
  model's. Earth's G1 models spread by a factor of two in their rain response.
- **PlaSim's broadband longwave** takes the glow down as if it had the air's own spectrum. The correlated-k check
  agrees on what reaches the ground; PlaSim's own absorption by height is not diagnosed, and its top layer, which warms
  most, holds all the air above about 107 km, the top 7.6% of its mass.
- **The clouds** are PlaSim's diagnostic, Earth-tuned scheme. The high cloud's loss under a warmer upper air
  follows from its humidity threshold; a cloud-resolving model could test it.
- **Where the glow lands.** The runs spread it evenly over each hemisphere; the glow study maps it by place and time.
- **People's heat** is released at the ground and is not in these runs.
- **The growth** follows the land's evaporation band by band, as the light study's water-limited land does; the
  timing of the rain within the lunar day and the soils' water beyond that are not followed.

## Next work

- The glow's efficacy went into the glow study at the integration of 17 October: its recommended answer takes 2.0% of
  the design sunlight, within 1% of the domain synthesis's convention, and its whole table of answers follows
  ([glow](../solar_shield_array/glow.md)). The rerun on the GCM's final radiation waits on the target climate.
- Build the water design into ecology's plants: the carbon per litre and the transpired share by band above.

## Files

| File | Holds |
|---|---|
| [run.py](run.py) | The comparison: the runs' surface means, layers and first means, the glow's dimming that holds the temperature, the rain and the growth there, what the plants make up, and the pairs |
| [extract.py](extract.py) | The layers and first 3-day means read from the GCM's output |
| [glow_absorption.py](glow_absorption.py) | The spectral check of where the glow and the returned infrared are absorbed |
| [test_run.py](test_run.py) | Their binding, regeneration and arithmetic |
| [sources.json](sources.json) | The published values and code read |
| [results/gases_glow_dimming.json](results/gases_glow_dimming.json), [results/glow_absorption.json](results/glow_absorption.json) | The products |

```
climate/gcm/.venv/bin/python -m research.studies.gases_glow_dimming.extract
python -m research.studies.gases_glow_dimming.glow_absorption
python -m research.studies.gases_glow_dimming.run
python -m pytest research/studies/gases_glow_dimming
```
