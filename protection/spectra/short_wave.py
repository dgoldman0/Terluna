"""Short-wave transmission of the stored titania-silica shield, from hard X-rays to the far ultraviolet.

    python -m protection.spectra.short_wave          # writes stack_short_wave.json

The shield product (transmission.py) starts at 200 nm and treats shorter light as blocked. This
evaluates the stored September design from 0.1 to 210 nm: the same anti-reflection layers, 1 um of
titania and 10 um of silica at the stored thicknesses, through protection/model.py's own
normal-incidence thin-film recursion and effective-medium mixing, with published optical constants:

- below 24.8 nm (above 50 eV), both oxides from CXRO atomic scattering factors (Henke, Gullikson and
  Davis 1993), the independent-atom estimate the design's own X-ray table uses and trusts from 50 eV;
- silica from 24.8 nm up: fused silica fitted over 24.8 nm-125 um (Franta et al. 2016);
- titania from 120 nm up: the design's own film data (Siefke et al. 2016). Between 24.8 and 120 nm no
  measured set was found, so the CXRO estimate is used there, outside its trusted range and without
  its real part (the tables give none below about 30 eV). The 10 um of silica is opaque across that
  range by itself, so the estimate does not change the result.

Densities and molar masses are the design table's (protection/model.py): silica 2.2 g/cm^3, titania
3.9 g/cm^3. The transmission is for direct sunlight at normal incidence through intact film. Gaps in
the aperture, pinholes, edges and scattered light are not included; they, not the film, decide how
much extreme ultraviolet the shield lets through.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np

from shared.constants import AVOGADRO, CLASSICAL_ELECTRON_RADIUS, ELEMENTARY_CHARGE, PLANCK, SPEED_OF_LIGHT

HERE = Path(__file__).resolve().parent
PROTECTION = HERE.parent
SOURCES = PROTECTION / 'sources'
SCHEMA = 'terluna.protection.stack-short-wave/1'
HC_EV_NM = PLANCK * SPEED_OF_LIGHT / ELEMENTARY_CHARGE * 1e9
CXRO_LIMIT_NM = 24.8                     # 50 eV; shorter wavelengths take both oxides from CXRO
SILICA = dict(formula={'si': 1, 'o': 2}, molar_mass_kg=0.0600843, density_kg_m3=2200.0)
TITANIA = dict(formula={'ti': 1, 'o': 2}, molar_mass_kg=0.079866, density_kg_m3=3900.0)
INPUTS = ('o.nff', 'si.nff', 'ti.nff', 'SiO2_Franta.yml', 'TiO2_Siefke.yml')


def _table(name):
    """Rows of three numbers (wavelength in um, n, k) from a refractiveindex.info file."""
    rows = []
    for line in (SOURCES / name).read_text().splitlines():
        parts = line.split()
        if len(parts) == 3:
            try:
                rows.append([float(v) for v in parts])
            except ValueError:
                pass
    return np.array(rows)


def titania_table_start_nm():
    """First wavelength of the design's measured titania data."""
    from protection import model
    return float(model.TI[0, 0] * 1000.0)


def cxro_index(wavelength_nm, formula, molar_mass_kg, density_kg_m3):
    """Complex refractive index n + ik of a compound from CXRO atomic scattering factors.

    n = 1 - (r_e lambda^2 / 2 pi) N sum(f1 + i f2) over its atoms. f2 is interpolated in log-log,
    as in the design's X-ray table; f1 is kept only where every atom has it."""
    w = np.asarray(wavelength_nm, dtype=float)
    energy = HC_EV_NM / w
    f1, f2, have_f1 = np.zeros_like(w), np.zeros_like(w), np.ones(w.shape, dtype=bool)
    for element, count in formula.items():
        e, a1, a2 = np.loadtxt(SOURCES / f'{element}.nff', skiprows=1).T
        f2 += count * np.exp(np.interp(np.log(energy), np.log(e), np.log(a2)))
        known = a1 > -9000.0
        have_f1 &= energy >= e[known][0]
        f1 += count * np.interp(np.log(energy), np.log(e[known]), a1[known])
    scale = CLASSICAL_ELECTRON_RADIUS * (w * 1e-9) ** 2 / (2 * np.pi) * density_kg_m3 / molar_mass_kg * AVOGADRO
    return 1.0 - scale * np.where(have_f1, f1, 0.0) + 1j * scale * f2


def indices(wavelength_nm):
    """Silica and titania refractive indices, from the sources the module notes give for each range."""
    from protection import model
    w = np.asarray(wavelength_nm, dtype=float)
    fs = _table('SiO2_Franta.yml')
    at = np.clip(w / 1000.0, fs[0, 0], fs[-1, 0])
    franta = np.interp(at, fs[:, 0], fs[:, 1]) + 1j * np.interp(at, fs[:, 0], fs[:, 2])
    silica = np.where(w < CXRO_LIMIT_NM, cxro_index(w, **SILICA), franta)
    start = titania_table_start_nm()
    titania = np.where(w < start, cxro_index(w, **TITANIA),
                       model.ti_nk(np.clip(w, start, model.TI[-1, 0] * 1000.0) / 1000.0))
    return silica, titania


def design_transmission(wavelength_nm, silica, titania, ar):
    """The stored design's layer sequence, built as protection/model.py's optical_stack builds it."""
    from protection import model
    w_um = np.asarray(wavelength_nm, dtype=float) / 1000.0
    porous = model.mix_n(np.ones_like(silica), silica, 0.5)
    mixed = model.mix_n(silica, titania, 0.5)
    layers = [porous, silica, mixed, titania, mixed, silica, porous]
    thickness_um = [ar[0], ar[1], ar[2], 1.0, ar[3], 10.0, ar[4]]
    return model.stack_rt(w_um, layers, thickness_um)[1]


def stored_ar_layers():
    return np.array(json.loads((PROTECTION / 'results' / 'results.json').read_text())['optics']['AR_layers_um'])


def compute():
    wavelength = np.round(np.arange(0.1, 210.0 + 1e-9, 0.05), 3)
    silica, titania = indices(wavelength)
    t = design_transmission(wavelength, silica, titania, stored_ar_layers())
    # The silica alone, by Beer-Lambert, where titania rests on the CXRO estimate.
    gap = (wavelength >= CXRO_LIMIT_NM) & (wavelength < titania_table_start_nm())
    silica_alone = np.exp(-4 * np.pi * silica.imag[gap] * 10.0e3 / wavelength[gap])
    return wavelength, t, float(silica_alone.max())


def main():
    wavelength, t, silica_gap_max = compute()
    at = lambda nm: float(t[int(np.argmin(abs(wavelength - nm)))])
    code = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in (HERE / 'short_wave.py', PROTECTION / 'model.py')}
    inputs = {n: hashlib.sha256((SOURCES / n).read_bytes()).hexdigest()[:16] for n in INPUTS}
    summary = {
        'transmission_at_nm': {f'{nm:g}': at(nm) for nm in (0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0,
                                                           121.6, 150.0, 200.0)},
        'max_transmission_2.5_to_210_nm': float(t[wavelength >= 2.5].max()),
        'silica_alone_max_transmission_24.8_to_120_nm': silica_gap_max,
    }
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='protection', files=code, inputs=inputs,
                      stored_design='protection/results/results.json optics.AR_layers_um; 1 um titania, 10 um silica'),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('Multiply the top-of-atmosphere solar spectral irradiance by T(lambda), interpolating linearly '
                      'between the 0.05 nm samples. For the titania stack this replaces the "at least 99.9% blocked '
                      'below 200 nm" of shield_transmission.json, for the film only: light that passes gaps, '
                      'pinholes or edges of the aperture is not included.'),
        units=dict(wavelength='nm (vacuum)', transmission='fraction of direct sunlight at normal incidence'),
        ranges=[dict(nm=f'0.1-{CXRO_LIMIT_NM}', silica='CXRO', titania='CXRO'),
                dict(nm=f'{CXRO_LIMIT_NM}-{titania_table_start_nm():.1f}', silica='Franta et al. 2016',
                     titania='CXRO, independent-atom absorption only (outside its trusted range)'),
                dict(nm=f'{titania_table_start_nm():.1f}-210', silica='Franta et al. 2016', titania='Siefke et al. 2016')],
        summary=summary,
        wavelength_nm=wavelength.tolist(),
        transmission=[float(f'{x:.4g}') for x in t])
    (HERE / 'stack_short_wave.json').write_text(json.dumps(product, separators=(',', ':')) + '\n')
    for k, v in summary['transmission_at_nm'].items():
        print(f'T({k} nm) = {v:.3e}')
    print('max T beyond 2.5 nm', f"{summary['max_transmission_2.5_to_210_nm']:.2e}",
          '| silica alone, 24.8-120 nm, at most', f'{silica_gap_max:.2e}')


if __name__ == '__main__':
    main()
