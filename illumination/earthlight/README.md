# Earthlight

The Earth's light at the Moon by wavelength. [model.py](model.py) takes the
spectral shape of the whole Earth from Glenar et al. (2019), a fit to NASA Virtual
Planetary Laboratory models of the Earth seen from the Moon, and scales it at each
phase angle to the Earth's visual brightness as Robinson et al. (2025) fit it to
ground and spacecraft photometry from 5° to 144° of phase.

That visual brightness is lower than the value long quoted. Robinson et al. give
the Earth a visual geometric albedo of 0.242 and a phase integral of 1.22; the
0.367 of Allen's tables, which [shared/constants.json](../../shared/constants.json)
still carries, rests on Danjon's earthshine work, extrapolated to full phase
before the Moon's opposition effect was known. A full Earth at its mean distance
puts 7.8 lux on the Moon above the air, against 13.5 lux for a grey Lambert
sphere of albedo 0.367, and its light is bluer than sunlight: chromaticity
(0.295, 0.303) against (0.322, 0.332).

The spectral shape holds within the model's 10–12% near 0.3 µm and its phase
curve to 60°; beyond the 144° of Robinson et al.'s data the visual scale stays at
its 144° value and the model's curve carries the Earth to dark at new Earth. The
model's spectrum is 11% bluer, blue band to red, than Robinson et al.'s band
albedos, and the turning continents and clouds move the Earth's light by about
10% through a day.

```sh
python -m illumination.earthlight.fetch_inputs --download
```

[inputs.json](inputs.json) records the table's source and its CC BY 4.0 licence.
