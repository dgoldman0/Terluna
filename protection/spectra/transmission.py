"""Spectral transmission of candidate spectral shields, written as a data product.

    python -m protection.spectra.transmission          # writes shield_transmission.json

Three shields, each a transmission T(lambda) applied to all sunlight reaching the Moon:

- titania_stack: the September protection design (protection/model.py): anti-
  reflection layers, 1 um titania and 10 um silica at the stored optimised
  thicknesses, computed with that model's own normal-incidence thin-film code
  from 210 nm (its silica dispersion limit) to 2730 nm. Titania is opaque below
  ~330 nm. Beyond 2730 nm the 2500-2730 nm mean is held, since the model
  neglects silica absorption there; below 210 nm it is opaque.
- edge_200nm: an idealised filter, not a designed coating: 99.9% blocking below
  200 nm (the extreme and far ultraviolet that heats the upper atmosphere), a
  smooth rise to full pass by 210 nm, the titania stack's transmission above
  420 nm and its 420-nm value through the near ultraviolet, so that only the
  ultraviolet differs from the current design. Sunlight that splits O2
  (200-242 nm) reaches the atmosphere, so an ozone layer can form.
- edge_310nm: the same idealisation with the edge at 310 nm: no light that
  splits O2, part of the UV-B and all UV-A, with no ozone layer to rely on.

The two edges are scenarios for the atmosphere; whether a coating can realise
them (hafnia and alumina have band edges near 220 and 150 nm) is open design work.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PROTECTION = HERE.parent
SCHEMA = 'terluna.protection.shield-transmission/1'


def titania_stack(wavelength_nm):
    from protection import model
    ar = json.loads((PROTECTION / 'results' / 'results.json').read_text())['optics']['AR_layers_um']
    w = np.asarray(wavelength_nm, dtype=float)
    t = np.zeros_like(w)
    ok = (w >= 210.0) & (w <= 2730.0)
    _, t[ok] = model.optical_stack(w[ok] / 1000.0, np.array(ar))
    tail = (w > 2730.0)
    ref = (w >= 2500.0) & (w <= 2730.0)
    t[tail] = t[ref].mean()
    return t


JOIN_NM = 420.0   # above this the idealised filters copy the titania stack exactly


def edge(wavelength_nm, stack, edge_nm, blocked=1e-3, width_nm=10.0):
    """Idealised edge filter: the stack's transmission above JOIN_NM, its value at
    JOIN_NM through the near ultraviolet, and a smooth edge to `blocked` below edge_nm."""
    w = np.asarray(wavelength_nm, dtype=float)
    plateau = float(np.interp(JOIN_NM, w, stack))
    passing = np.where(w >= JOIN_NM, stack, plateau)
    x = np.clip((w - edge_nm) / width_nm, 0.0, 1.0)
    s = x * x * (3 - 2 * x)
    return blocked + (passing - blocked) * s


def compute():
    wavelength = np.round(np.arange(200.0, 5000.0 + 0.25, 0.5), 3)
    stack = titania_stack(wavelength)
    return wavelength, {'titania_stack': stack,
                        'edge_200nm': edge(wavelength, stack, 200.0),
                        'edge_310nm': edge(wavelength, stack, 310.0)}


def main():
    wavelength, shields = compute()
    code = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in (HERE / 'transmission.py', PROTECTION / 'model.py')}
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='protection', files=code, stored_design='protection/results/results.json optics.AR_layers_um'),
        evidence=('titania_stack is the stored September design evaluated with its own thin-film code (normal incidence, '
                  'effective-medium mixed layers, no silica absorption, no coating ageing); edge_200nm and edge_310nm '
                  'are idealised filters, not designs. Transmission applies to direct sunlight; gaps in the aperture '
                  'and off-normal incidence are not included.'),
        reading_rule='Multiply the top-of-atmosphere solar spectral irradiance by T(lambda); interpolate linearly in wavelength; '
                     'below 200 nm every shield is treated as blocking at least 99.9%.',
        units=dict(wavelength='nm (vacuum)', transmission='fraction'),
        join_nm=JOIN_NM,
        wavelength_nm=wavelength.tolist(),
        transmission={k: [round(float(x), 6) for x in v] for k, v in shields.items()})
    out = HERE / 'shield_transmission.json'
    out.write_text(json.dumps(product, separators=(',', ':')) + '\n')
    for k, v in shields.items():
        sel = lambda lo, hi: (wavelength >= lo) & (wavelength <= hi)
        print(f"{k:14s} T(250)={v[sel(249.9,250.1)][0]:.2e} T(300)={v[sel(299.9,300.1)][0]:.2e} T(350)={v[sel(349.9,350.1)][0]:.3f} "
              f"T(400)={v[sel(399.9,400.1)][0]:.3f} T(550)={v[sel(549.9,550.1)][0]:.3f} T(1000)={v[sel(999.9,1000.1)][0]:.3f}")


if __name__ == '__main__':
    main()
