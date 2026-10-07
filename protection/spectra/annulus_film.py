"""Ultraviolet transmission of light annulus films, from soft X-rays to the near infrared.

    python -m protection.spectra.annulus_film          # writes annulus_film.json

Beyond the climate window the shield's aperture, out to four lunar radii, has to stop the ultraviolet and pass
visible light (research/decisions.md): its rays cross the upper air and exosphere above the limb and go on toward
Earth, which must not see a deep shadow. The window carries the 26 g/m^2 climate stack, 1 um of titania on 10 um of
silica with anti-reflection layers. This evaluates thinner films of the same two oxides for the annulus, titania
sunward on a silica support, bare and with the stack's stored anti-reflection and matching layers, through
protection/model.py's normal-incidence thin-film recursion: below 210 nm with short_wave's optical constants (CXRO
below 24.8 nm, Franta silica, Siefke titania from 120 nm and the CXRO estimate for titania between), above it with
model.py's (Siefke titania, Malitson silica without absorption).

Each film is weighted by the quiet-Sun WHI 2008 spectrum below 202 nm and TSIS-1 above it. Requirement O1 counts
the sunlight below 175 nm; its transmission is given with and without the soft X-rays below 2.5 nm, which no film of
this kind stops and the window stack also passes. The heat that light deposits in the upper air follows the
atmosphere domain's film_heat (all of it absorbed above the base, an upper bound for X-rays, at quiet Sun and at
solar maximum), beside its heat for light that bypasses the film at O1's swarm levels. Reflected and absorbed shares of the whole spectrum give the light
pressure on the film at normal incidence, (2R + A) times the solar constant over c: reflected light pushes twice, absorbed
light once.

The transmission is for direct sunlight at normal incidence through intact film; oblique light crosses more
absorber and is stopped better. Gaps, pinholes, edges and scattered light are not included; with the film's own
transmission far below them, they decide what the annulus lets through. The support's strength, handling and
ageing are not modelled.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np

from protection.spectra import short_wave as sw

HERE = Path(__file__).resolve().parent
PROTECTION = HERE.parent
SCHEMA = 'terluna.protection.annulus-film/1'
TITANIA_UM = (0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0)
SILICA_UM = (0.5, 1.0, 2.0, 4.0, 10.0)
LEVELS = (1e-5, 1e-6, 1e-8)
SPLIT_NM = 210.0
BANDS_NM = {'xray_0.1_2.5': (0.1, 2.5), 'euv_2.5_121': (2.5, 121.0), 'lyman_alpha_121_122.5': (121.0, 122.5),
            'fuv_122.5_175': (122.5, 175.0), '175_200': (175.0, 200.0), '200_242': (200.0, 242.0),
            '242_320': (242.0, 320.0), '320_400': (320.0, 400.0), 'visible_400_700': (400.0, 700.0),
            'near_infrared_700_2500': (700.0, 2500.0)}
DENSITY_G_UM = dict(silica=2.2, titania=3.9, porous=1.1, mixed=3.05)   # g/m^2 per micrometre
INPUTS = sw.INPUTS+('TSIS1_HSRS_stride100.csv',)


def grid():
    """Wavelengths (nm): short_wave's 0.05 nm steps to 210 nm, then 0.5 nm steps to 2,500 nm."""
    return np.r_[np.round(np.arange(0.1, SPLIT_NM, 0.05), 3), np.arange(SPLIT_NM, 2500.0+1e-9, 0.5)]


def indices(wavelength_nm):
    """Silica and titania indices: short_wave's below 210 nm, protection/model.py's above."""
    from protection import model
    w = np.asarray(wavelength_nm, float)
    low = w < SPLIT_NM
    silica, titania = np.zeros(w.shape, complex), np.zeros(w.shape, complex)
    if low.any():
        silica[low], titania[low] = sw.indices(w[low])
    if (~low).any():
        silica[~low] = model.si_n(w[~low]/1000.0)
        titania[~low] = model.ti_nk(w[~low]/1000.0)
    return silica, titania


def film(wavelength_nm, silica, titania, titania_um, silica_um, coated):
    """Reflectance and transmittance from the sunward side: titania on silica, bare or with the stack's layers."""
    from protection import model
    w_um = np.asarray(wavelength_nm, float)/1000.0
    if not coated:
        return model.stack_rt(w_um, [titania, silica], [titania_um, silica_um])
    ar = sw.stored_ar_layers()
    porous = model.mix_n(np.ones_like(silica), silica, 0.5)
    mixed = model.mix_n(silica, titania, 0.5)
    return model.stack_rt(w_um, [porous, silica, mixed, titania, mixed, silica, porous],
                          [ar[0], ar[1], ar[2], titania_um, ar[3], silica_um, ar[4]])


def areal_mass(titania_um, silica_um, coated):
    """Film areal mass (g/m^2) from the design table's densities; mixed layers half silica, half titania."""
    mass = DENSITY_G_UM['titania']*titania_um+DENSITY_G_UM['silica']*silica_um
    if coated:
        ar = sw.stored_ar_layers()
        mass += DENSITY_G_UM['porous']*(ar[0]+ar[4])+DENSITY_G_UM['silica']*ar[1]+DENSITY_G_UM['mixed']*(ar[2]+ar[3])
    return float(mass)


def spectrum(wavelength_nm):
    """Solar spectral irradiance (W m^-2 nm^-1): quiet-Sun WHI 2008 below 202 nm, TSIS-1 HSRS above."""
    from atmosphere.middle_atmosphere.escape import whi_quiet_sun
    w = np.asarray(wavelength_nm, float)
    whi_w, whi_f = whi_quiet_sun()
    hsrs = np.loadtxt(sw.SOURCES/'TSIS1_HSRS_stride100.csv', delimiter=',', skiprows=1)
    return np.where(w < 202.0, np.interp(w, whi_w, whi_f), np.interp(w, hsrs[:, 0], hsrs[:, 1]))


def power(w, f, lo, hi):
    sel = (w >= lo) & (w < hi)
    return float(np.trapezoid(f[sel], w[sel]))


def weighted(w, f, x, lo, hi):
    sel = (w >= lo) & (w < hi)
    return float(np.trapezoid(x[sel]*f[sel], w[sel])/np.trapezoid(f[sel], w[sel]))


def heat(w, f, t):
    """Heat of the upper air (W/m^2 of lunar surface, global mean) from the film's transmission below 175 nm, by the
    escape model's film_heat: all of it absorbed above the base with its heating efficiency, quiet Sun and solar
    maximum (FISM2's measured rise by band below 10 nm and 2.5 for the rest). The middle atmosphere's limb tracing
    gives the heat along the annulus's slant paths."""
    from atmosphere.middle_atmosphere import escape
    sel = w < 175.0
    passed = t[sel]*f[sel]
    scale = np.where(w[sel] < 10.0, escape.xray_scale(w[sel], escape.XRAY_SOLAR_MAXIMUM), escape.SOLAR_MAXIMUM)
    factor = 0.25*escape.HEATING_EFFICIENCY
    return dict(quiet=float(factor*np.trapezoid(passed, w[sel])), maximum=float(factor*np.trapezoid(passed*scale, w[sel])))


def evaluate(w, f, r, t, mass):
    a = 1.0-r-t
    uv = dict(below_175=weighted(w, f, t, 0.1, 175.0), from_2_5_to_175=weighted(w, f, t, 2.5, 175.0),
              max_2_5_to_175=float(t[(w >= 2.5) & (w < 175.0)].max()),
              max_10_to_175=float(t[(w >= 10.0) & (w < 175.0)].max()))
    bands = {name: weighted(w, f, t, lo, hi) for name, (lo, hi) in BANDS_NM.items()}
    passed = {name: power(w, f*t, lo, hi) for name, (lo, hi) in BANDS_NM.items() if hi <= 200.0}
    bol = dict(R=weighted(w, f, r, 0.1, 2500.0), A=weighted(w, f, a, 0.1, 2500.0), T=weighted(w, f, t, 0.1, 2500.0))
    pressure = 2*bol['R']+bol['A']
    return dict(areal_mass_g_m2=mass, uv_transmission=uv, band_transmission=bands, transmitted_W_m2=passed,
                upper_air_heat_W_m2=heat(w, f, t), bolometric=bol,
                visible_reflectance=weighted(w, f, r, 400.0, 700.0), pressure_coefficient=pressure,
                pressure_coefficient_per_g_m2=pressure/mass)


def gap_heat():
    """The escape model's heat from light that bypasses the film at the swarm levels of O1 (O8), for comparison."""
    from atmosphere.middle_atmosphere import escape
    return {f'{level:g}': dict(quiet=escape.leakage_heat(level),
                               maximum=escape.leakage_heat(level, activity=escape.SOLAR_MAXIMUM,
                                                           xray_activity=escape.XRAY_SOLAR_MAXIMUM))
            for level in (3e-5, 2e-4)}


def compute():
    w = grid()
    f = spectrum(w)
    silica, titania = indices(w)
    rows = []
    for coated in (False, True):
        for t_ti in TITANIA_UM:
            for t_si in SILICA_UM:
                r, t = film(w, silica, titania, t_ti, t_si, coated)
                rows.append(dict(titania_um=t_ti, silica_um=t_si, coated=coated,
                                 **evaluate(w, f, r, t, areal_mass(t_ti, t_si, coated))))
    window = [x for x in rows if x['coated'] and x['titania_um'] == 1.0 and x['silica_um'] == 10.0][0]
    lightest = {}
    for level in LEVELS:
        for coated in (False, True):
            ok = [x for x in rows if x['coated'] == coated and x['uv_transmission']['from_2_5_to_175'] <= level
                  and x['uv_transmission']['max_10_to_175'] <= 10*level]
            if ok:
                best = min(ok, key=lambda x: x['areal_mass_g_m2'])
                lightest[f"{level:g}_{'coated' if coated else 'bare'}"] = {
                    k: best[k] for k in ('titania_um', 'silica_um', 'areal_mass_g_m2', 'pressure_coefficient')}
    return w, f, rows, window, lightest


def main():
    w, f, rows, window, lightest = compute()
    from atmosphere.middle_atmosphere import escape, fetch_inputs, fetch_limb_inputs
    code = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
            for p in (HERE/'annulus_film.py', HERE/'short_wave.py', PROTECTION/'model.py', Path(escape.__file__))}
    inputs = {n: hashlib.sha256((sw.SOURCES/n).read_bytes()).hexdigest()[:16] for n in INPUTS}
    for path in (fetch_inputs.path('whi2008_ref_solar_irradiance_ver2.dat'),
                 *(fetch_limb_inputs.path(n) for n in (escape.FISM2_QUIET, *escape.FISM2_MAXIMA.values()))):
        inputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='protection', files=code, inputs=inputs,
                      stored_design='protection/results/results.json optics.AR_layers_um'),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('A film meets a level when its WHI-weighted transmission between 2.5 and 175 nm is at or below it '
                      'and no wavelength between 10 and 175 nm passes more than ten times that. Transmissions are '
                      'shares of the band\'s sunlight; the pressure coefficient times 1361/c W/m^2 is the light pressure '
                      'at normal incidence. The window row is the climate stack, as the dynamics models carry it the '
                      'stack and dimmer redirect 13.8% of the light (pressure coefficient 0.276 if all reflected).'),
        units=dict(wavelength='nm (vacuum)', thickness='micrometre', areal_mass='g/m^2'),
        designs=dict(titania_um=list(TITANIA_UM), silica_um=list(SILICA_UM), coatings=['bare', 'stored layers']),
        band_power_W_m2={name: power(w, f, lo, hi) for name, (lo, hi) in BANDS_NM.items()},
        gap_heat_W_m2=gap_heat(),
        window_stack=window, lightest_meeting_level=lightest, films=rows)
    (HERE/'annulus_film.json').write_text(json.dumps(product, indent=1)+'\n')
    for x in rows:
        print(f"{'AR' if x['coated'] else '  '} Ti {x['titania_um']:4.2f} Si {x['silica_um']:4.1f} um "
              f"{x['areal_mass_g_m2']:6.2f} g/m2  T<175 {x['uv_transmission']['below_175']:.1e} "
              f"T2.5-175 {x['uv_transmission']['from_2_5_to_175']:.1e} max {x['uv_transmission']['max_10_to_175']:.1e} "
              f"vis {x['band_transmission']['visible_400_700']:.3f} R {x['bolometric']['R']:.3f} "
              f"A {x['bolometric']['A']:.3f} 2R+A {x['pressure_coefficient']:.3f}")
    print('lightest', json.dumps(lightest, indent=1))


if __name__ == '__main__':
    main()
