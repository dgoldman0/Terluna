"""Fixed-column UV-window screen using stored spectral transfer, no new equilibrium.

python -m research.studies.uv_biology.run
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from biosphere.ecology.screen import erythema, govardovskii
from shared.constants import AVOGADRO, PLANCK, SPEED_OF_LIGHT
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'illumination/surface_light/results/surface_light.json'
SOURCES = ROOT / 'research/studies/aerial_ecology/sources.json'
SCHEMA = 'terluna.research.uv-biology/1'
OUT = HERE / 'results/uv_biology.json'
WINDOWS = {'reference': None, 'uvb_295_315_1pct': (295,315,0.01),
           'uvb_295_315_5pct': (295,315,0.05), 'uvb_295_315_10pct': (295,315,0.1),
           'long_uvb_310_315_10pct': (310,315,0.1), 'uv_310_400_95pct': (310,400,0.95)}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def spectrum(product, case, sun):
    c = product['cases'][case]
    r = np.asarray(product['columns'][c['column']]['reflectance_from_below'])
    s = c['spectra'][sun]
    return (np.asarray(s['direct_w_m2_nm']) + np.asarray(s['diffuse_black_w_m2_nm'])) / (1 - 0.1*r)


def window_spectrum(w, baseline, unfiltered, window):
    result = baseline.copy()
    if window:
        lo, hi, transmission = window
        band = (w >= lo) & (w < hi)
        # Absolute fraction of unfiltered incoming sunlight inside the band.
        result[band] = transmission * unfiltered[band]
    return result


def metrics(w, flux):
    width = np.diff(w)
    if not np.allclose(width, 1):
        raise ValueError('this stored product requires 1 nm bins')
    photons = flux * w*1e-9 / (PLANCK*SPEED_OF_LIGHT*AVOGADRO) * 1e6
    bands = {f'{lo}_{hi}': float(photons[(w >= lo)&(w < hi)].sum())
             for lo, hi in ((280,315),(295,305),(310,315),(315,350),(350,400))}
    return dict(uv_index=float(40*np.sum(flux*erythema(w))),
                band_photons_umol_m2_s=bands,
                bee_uv_relative_catch=float(np.sum(photons*govardovskii(w,344))),
                uv_energy_w_m2=float(flux[(w >= 280)&(w < 400)].sum()))


def run():
    d = json.loads(SOURCE.read_text())
    if d['schema'] != 'terluna.illumination.surface-light/1':
        raise ValueError(d['schema'])
    w = np.asarray(d['centres_nm'])
    rows = {}
    for sun in ('90.0', '30.0', '10.0'):
        ref = spectrum(d, 'moon_1.2atm/design', sun)
        bare = spectrum(d, 'moon_1.2atm/unfiltered', sun)
        earth = metrics(w, spectrum(d, 'earth_control/unfiltered', sun))
        rows[sun] = dict(earth=earth, cases={})
        for name, window in WINDOWS.items():
            flux = window_spectrum(w, ref, bare, window)
            m = metrics(w, flux)
            m['added_uv_energy_w_m2'] = m['uv_energy_w_m2'] - metrics(w, ref)['uv_energy_w_m2']
            m['bee_uv_catch_vs_earth'] = m.pop('bee_uv_relative_catch') / earth['bee_uv_relative_catch']
            m['unchanged_below_295nm'] = bool(np.array_equal(flux[w < 295], ref[w < 295]))
            rows[sun]['cases'][name] = m
    paths = [Path(__file__), ROOT/'biosphere/ecology/screen.py']
    out = dict(schema=SCHEMA,
        producer=dict(files={str(p.relative_to(ROOT)):digest(p) for p in paths},
                      inputs={str(p.relative_to(ROOT)):digest(p) for p in (SOURCE, SOURCES)}, constants=constants_used(paths)),
        evidence='Linear spectral reweighting of the solved, almost ozone-free lunar reference column; '
                 'multiple scattering is already in its transfer. No new ozone, living-atmosphere chemistry, '
                 'film design, ring thermal balance or exobase calculation is performed.',
        reading_rule='Window fractions are absolute incoming transmission, not percentages of the old film. '
                     'Outside each band all filtering is retained, including wavelengths below the stored '
                     '202 nm limit. Values are horizontal irradiance at ground albedo 0.1. Photon bands are '
                     'unweighted proxies, not vitamin-D, UVR8 or DNA-damage action spectra.',
        windows=WINDOWS, by_sun_degrees=rows,
        adoption=dict(current_spectrum_retained=True, added_uv_adopted=False,
                      required=['current ring geometry and spectral stack', 'live-atmosphere photochemistry',
                                'vertical actinic biological doses', 'non-LTE exobase solar-cycle feedback',
                                'film infrared and compensating climate dimming']))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    return out


if __name__ == '__main__':
    run()
