"""Eddy mixing of the upper air by breaking gravity waves, on Earth and on the Open Moon: a screening estimate.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.gravity_waves   # writes results/gravity_waves.json

The middle atmosphere scales Earth's eddy diffusion by (g_Earth / g_Moon)^2 = 36 (run.py), which keeps the mixing time
per scale height. That is the scaling of Lindzen's (1981) saturated gravity-wave mixing, K = k c^4 / (2 H N^3), with
Earth's waves unchanged: the buoyancy frequency N falls and the scale height H grows six-fold at a sixth of the
gravity. How strongly the air above the 0.3 Pa base mixes decides how much atomic oxygen reaches the exobase
(atmosphere/loss_response/oxygen.py), so this step follows the waves themselves:

- each wave of a spectrum (phase speeds c, both directions alike with no mean wind, and horizontal wavelengths sharing
  the flux evenly) rises from half the tropopause pressure with its share of the launch momentum flux. Its vertical
  wavenumber follows the full dispersion relation, m^2 = k^2 (N^2 / w^2 - 1) - 1/(4 H^2) with w = k c, and the wave
  turns back where m^2 falls below zero: in the lunar thermosphere, where N falls to about 1e-3 per second, waves
  short and fast enough that k c passes N cannot rise. The Coriolis parameter is left out, as it may be on the slowly
  turning Moon;
- molecular viscosity and heat conduction take the flux at (nu + kappa) (k^2 + m^2) over the vertical group velocity;
  where the flux would pass the saturation flux rho k c^2 / (2 m), at which the wave overturns, the wave breaks and
  holds there, and the turbulence holding it diffuses at D = [(1/H_rho + d ln m/dz) c_g / (k^2 + m^2) - (nu + kappa)]
  / 2, Lindzen's k c^4 / (2 H N^3) for a hydrostatic wave in an isothermal column, summed over the waves;
- the constituent eddy diffusion is that sum times one factor, the waves' intermittency over their effective Prandtl
  number, set so that the same calculation through Earth's standard atmosphere (US Standard Atmosphere 1976) gives
  Earth's measured eddy diffusion at its homopause;
- on the Moon the same spectrum runs through the loss response's column at a state, ground to exobase, with gravity
  falling with radius, and its wavelengths either stretched by the ratio of gravities, as the climate domain's CM1
  runs stretch Earth's lengths for a troposphere six times deeper, or left at Earth's.

Left out: mean winds and their filtering of the spectrum, the waves' own driving of winds, tides and the day-night
circulation, which also carry the air up and down, the wave sources the lunar troposphere actually has, and the
breakdown of the slowly varying wave picture where a vertical wavelength nears the scale height.
"""
from __future__ import annotations
import functools
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from shared.constants import BOLTZMANN, MOON_GM, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'gravity_waves.json'
SCHEMA = 'terluna.atmosphere.gravity-wave-mixing/1'

GAS_CONSTANT = 8.314462618        # J/(mol K)
# The launched spectrum, Earth's: 4 mPa from half the tropopause pressure, about 15 km on Earth (Alexander and Dunkerton
# 1999, JAS 56, 4167, who cite observed means of 3-6 mPa outside the tropics), Gaussian in phase speed with a
# 30 m/s e-folding (half-width at half-maximum 25 m/s; their broad case's 60 m/s puts Earth's homopause near 119 km
# under the same calibration, against the observed 100-110), and shared evenly among horizontal wavelengths of 50 to
# 800 km (Holton and Zhu 1984 take about 200 km, Alexander and Dunkerton 10 to 1000 km).
PHASE_SPEEDS_M_S = tuple(float(c) for c in range(5, 81, 5))
PHASE_WIDTH_M_S = 30.0
WAVELENGTHS_KM = (50.0, 100.0, 200.0, 400.0, 800.0)
LAUNCH_FLUX_PA = 4e-3
LAUNCH_OVER_TROPOPAUSE = 0.5

# US Standard Atmosphere 1976 (NOAA, NASA and USAF): sea-level values, the geopotential layers to 84.852 km' (86 km
# geometric) and the defined thermosphere above.
US76_EARTH_RADIUS_KM = 6356.766
US76_MOLAR = 0.0289644
US76_LAYERS = ((0.0, -6.5), (11.0, 0.0), (20.0, 1.0), (32.0, 2.8), (47.0, 0.0), (51.0, -2.8), (71.0, -2.0),
               (84.8520, None))
US76_SURFACE = (101325.0, 288.15)
US76_TROPOPAUSE_KM = 11.019      # the standard's tropopause (11 km' geopotential)
# Earth's measured constituent eddy diffusion, which calibrates the raw diffusion: 1.1e6 cm^2/s, flat from 80 to
# 96 km in the equatorial annual mean from SABER and SCIAMACHY atomic oxygen (Swenson et al. 2021, JGR Atmospheres,
# doi:10.1029/2021JD035343), taken at the top of that range. Other measures run from (1-3)e5 in the upper mesosphere
# (Strobel 1989) to 1.8e6 (Lubken 1997) and more; the calibrated scheme puts Earth's homopause at 110 km.
EARTH_HOMOPAUSE_KM = 96.0
EARTH_HOMOPAUSE_K_CM2_S = 1.1e6
# Lunar waves are stretched by the ratio of gravities, as the climate domain's CM1 runs stretch Earth's lengths.
STRETCH = STANDARD_GRAVITY / MOON_SURFACE_GRAVITY
# Mean molar mass above 86 km (g/mol at km), which reproduces the standard's tabulated pressures within 2% to 150 km.
US76_MOLAR_ABOVE = ((86.0, 28.9522), (100.0, 28.40), (110.0, 27.27), (120.0, 26.20), (150.0, 24.10))


def us76_temperature(z_km):
    """US Standard Atmosphere 1976 temperature (K) at geometric heights (km), up to 1000 km."""
    z = np.atleast_1d(np.asarray(z_km, float))
    r = US76_EARTH_RADIUS_KM
    h = r * z / (r + z)                                        # geopotential height, km'
    t = np.empty_like(z)
    base_t = US76_SURFACE[1]
    for (h0, lapse), (h1, _) in zip(US76_LAYERS[:-1], US76_LAYERS[1:]):
        sel = (h >= h0) & (h < h1)
        t[sel] = base_t + lapse * (h[sel] - h0)
        base_t += lapse * (h1 - h0)
    t[(z >= 86.0) & (z < 91.0)] = 186.8673
    sel = (z >= 91.0) & (z < 110.0)
    t[sel] = 263.1905 - 76.3232 * np.sqrt(1.0 - ((z[sel] - 91.0) / -19.9429) ** 2)
    sel = (z >= 110.0) & (z < 120.0)
    t[sel] = 240.0 + 12.0 * (z[sel] - 110.0)
    sel = z >= 120.0
    xi = (z[sel] - 120.0) * (r + 120.0) / (r + z[sel])
    t[sel] = 1000.0 - (1000.0 - 360.0) * np.exp(-0.01875 * xi)
    # Between 84.852 km' and 86 km the layer above 71 km' continues.
    return t


def us76_column(top_km=150.0, step_km=0.25):
    """Earth's standard column (surface first): heights (m), pressures (Pa), temperatures (K), gravity (m/s^2) and
    mean molar mass (kg/mol), integrated hydrostatically from sea level."""
    z_km = np.arange(0.0, top_km + 1e-9, step_km)
    t = us76_temperature(z_km)
    g = STANDARD_GRAVITY * (US76_EARTH_RADIUS_KM / (US76_EARTH_RADIUS_KM + z_km)) ** 2
    above_z, above_m = zip(*US76_MOLAR_ABOVE)
    molar = np.where(z_km <= above_z[0], US76_MOLAR, np.interp(z_km, above_z, above_m) * 1e-3)
    dz = np.diff(z_km) * 1e3
    integrand = molar * g / (GAS_CONSTANT * t)
    log_p = np.log(US76_SURFACE[0]) - np.r_[0.0, np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dz)]
    return dict(z_m=z_km * 1e3, p_pa=np.exp(log_p), t_k=t, g_m_s2=g, molar_kg_mol=molar)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def viscosity(t):
    """Dynamic viscosity of air (Pa s): Sutherland's law as the US Standard Atmosphere 1976 gives it."""
    return 1.458e-6 * t ** 1.5 / (t + 110.4)


def conductivity(t):
    """Thermal conductivity of air (W/(m K)), as the US Standard Atmosphere 1976 gives it."""
    return 2.64638e-3 * t ** 1.5 / (t + 245.4 * 10.0 ** (-12.0 / t))


def properties(col, cp=1005.0):
    """Buoyancy frequency (s^-1), density (kg/m^3), inverse density scale height and d ln N/dz (1/m), and kinematic
    viscosity and thermal diffusivity (m^2/s) at a column's levels (z_m, p_pa, t_k, g_m_s2, molar_kg_mol)."""
    z, p, t, g = (np.asarray(col[k], float) for k in ('z_m', 'p_pa', 't_k', 'g_m_s2'))
    molar = np.asarray(col['molar_kg_mol'], float) * np.ones_like(t)
    rho = p * molar / (GAS_CONSTANT * t)
    n2 = g / t * (np.gradient(t, z) + g / cp)
    n = np.sqrt(np.maximum(n2, 1e-12))
    return dict(n=n, rho=rho, inv_h=-np.gradient(np.log(rho), z), dlnn=np.gradient(np.log(n), z),
                nu=viscosity(t) / rho, kappa=conductivity(t) / (rho * cp))


def spectrum(stretch=1.0, flux_pa=None, speeds_m_s=None, width_m_s=None, wavelengths_km=None):
    """The waves launched: (phase speed m/s, horizontal wavenumber 1/m, launch momentum flux Pa) for each phase speed
    and horizontal wavelength, the flux Gaussian in phase speed (both directions alike) and shared evenly among the
    wavelengths, which `stretch` multiplies."""
    flux_pa = LAUNCH_FLUX_PA if flux_pa is None else flux_pa
    c = np.asarray(PHASE_SPEEDS_M_S if speeds_m_s is None else speeds_m_s, float)
    weight = np.exp(-(c / (PHASE_WIDTH_M_S if width_m_s is None else width_m_s)) ** 2)
    lam = np.asarray(WAVELENGTHS_KM if wavelengths_km is None else wavelengths_km, float) * 1e3 * stretch
    share = flux_pa * weight / weight.sum() / lam.size
    return [(float(ci), 2 * math.pi / float(li), float(fi)) for ci, fi in zip(c, share) for li in lam]


def vertical_wavenumber(k, c, n, inv_h):
    """m^2 (1/m^2) of a wave of horizontal wavenumber k and intrinsic phase speed c: k^2 (N^2 / w^2 - 1) - 1/(4 H^2),
    with w = k c; where it is negative the wave cannot propagate and turns back."""
    w = k * c
    return k ** 2 * (n ** 2 / w ** 2 - 1.0) - 0.25 * inv_h ** 2


def diffusion(col, waves, launch_pa, cp=1005.0):
    """The turbulent diffusion (m^2/s) breaking waves need at each level, summed over the waves, and for each wave the
    pressure where it first breaks and where it turns back (Pa; None if it does not).

    A wave's momentum flux F holds against molecular damping, which takes it at (nu + kappa) (k^2 + m^2) over its
    vertical group velocity c_g = w m / (k^2 + m^2 + 1/(4 H^2)), until it reaches the saturation flux
    rho k c^2 / (2 m), where its horizontal wind equals c and it overturns. Held there, it needs turbulence (equal
    on momentum and heat) of D = [(1/H_rho + d ln m/dz) c_g / (k^2 + m^2) - (nu + kappa)] / 2, which is Lindzen's
    k c^4 / (2 H N^3) for a hydrostatic wave in an isothermal column. A wave stops where m^2 turns negative."""
    z, p = np.asarray(col['z_m'], float), np.asarray(col['p_pa'], float)
    q = properties(col, cp)
    start = int(np.argmax(p <= launch_pa))
    total = np.zeros_like(p)
    out = []
    for c, k, f0 in waves:
        m2 = vertical_wavenumber(k, c, q['n'], q['inv_h'])
        m = np.sqrt(np.maximum(m2, 1e-30))
        group = k * c * m / (k ** 2 + m2 + 0.25 * q['inv_h'] ** 2)
        damping = (q['nu'] + q['kappa']) * (k ** 2 + m ** 2) / group           # flux lost per metre
        saturation = 0.5 * q['rho'] * k * c ** 2 / m
        need = 0.5 * ((q['inv_h'] + np.gradient(np.log(m), z)) * group / (k ** 2 + m ** 2) - (q['nu'] + q['kappa']))
        flux, broke, turned = f0, None, None
        if m2[start] <= 0:
            out.append(dict(phase_speed_m_s=c, wavelength_km=2e-3 * math.pi / k, launch_flux_pa=f0, breaks_at_pa=None,
                            turns_at_pa=float(p[start])))
            continue
        flux = min(flux, saturation[start])
        for j in range(start + 1, p.size):
            if m2[j] <= 0:
                turned = float(p[j])
                break
            flux *= math.exp(-0.5 * (damping[j - 1] + damping[j]) * (z[j] - z[j - 1]))
            if flux >= saturation[j]:
                flux = saturation[j]
                if need[j] > 0:
                    total[j] += need[j]
                    broke = float(p[j]) if broke is None else broke
        out.append(dict(phase_speed_m_s=c, wavelength_km=2e-3 * math.pi / k, launch_flux_pa=f0, breaks_at_pa=broke,
                        turns_at_pa=turned))
    return total, out


def earth():
    """Earth's standard column and the waves' raw turbulent diffusion through it (m^2/s), launched from half the
    standard's tropopause pressure."""
    col = us76_column()
    tropopause = float(np.exp(np.interp(US76_TROPOPAUSE_KM * 1e3, col['z_m'], np.log(col['p_pa']))))
    raw, waves = diffusion(col, spectrum(), LAUNCH_OVER_TROPOPAUSE * tropopause)
    return col, raw, waves


@functools.lru_cache(maxsize=None)
def calibration():
    """The factor on the raw diffusion (the waves' intermittency over their effective Prandtl number) that gives
    Earth's measured eddy diffusion at its homopause, and the raw value there (m^2/s)."""
    col, raw, _ = earth()
    at = float(np.interp(EARTH_HOMOPAUSE_KM * 1e3, col['z_m'], raw))
    return EARTH_HOMOPAUSE_K_CM2_S * 1e-4 / at, at


def lunar_column(col, planet_gm, planet_radius_m, molar_kg_mol):
    """The wave module's view of a column of the loss response (levels: z_m, p_pa, t_k), with gravity falling with
    radius."""
    z = np.asarray(col['z_m'], float)
    return dict(z_m=z, p_pa=np.asarray(col['p_pa'], float), t_k=np.asarray(col['t_k'], float),
                g_m_s2=planet_gm / (planet_radius_m + z) ** 2, molar_kg_mol=molar_kg_mol)


def eddy_diffusion(col, tropopause_pa, stretch=STRETCH, factor=None):
    """Constituent eddy diffusion (cm^2/s) at a column's levels from the waves launched at half its tropopause
    pressure, Earth's spectrum with its wavelengths stretched by `stretch`, calibrated on Earth."""
    factor = calibration()[0] if factor is None else factor
    raw, _ = diffusion(col, spectrum(stretch), LAUNCH_OVER_TROPOPAUSE * tropopause_pa)
    return factor * raw * 1e4


def homopause(p_pa, t_k, eddy_cm2_s):
    """The homopause for atomic oxygen (Pa): the level above the highest at which eddy diffusion still matches its
    molecular diffusion in N2 (Banks and Kockarts 1973); None if the eddies never do, or do at the top (levels
    surface first)."""
    p, t, eddy = (np.asarray(x, float) for x in (p_pa, t_k, eddy_cm2_s))
    molecular = 9.69e16 * t ** 0.774 / (p / (BOLTZMANN * t) * 1e-6)
    mixed = np.nonzero(eddy >= molecular)[0]
    if not mixed.size or mixed[-1] == p.size - 1:
        return None
    return float(p[mixed[-1] + 1])


def main(argv=None) -> int:
    from atmosphere.loss_response import model as lr, oxygen as ox     # read late: oxygen imports this module
    from atmosphere.middle_atmosphere import run as middle_run
    factor, raw_at = calibration()
    col_e, raw_e, waves_e = earth()
    keep = slice(0, None, 8)
    out_earth = dict(z_km=(col_e['z_m'][keep] / 1e3).tolist(), p_pa=col_e['p_pa'][keep].tolist(),
                     t_k=col_e['t_k'][keep].tolist(), eddy_cm2_s=(factor * raw_e[keep] * 1e4).tolist(),
                     homopause_pa=homopause(col_e['p_pa'], col_e['t_k'], factor * raw_e * 1e4), waves=waves_e)
    lunar = {}
    transmission = ox.standard_transmission()
    activity = ox.activity_named('cycles_23_24_mean')
    for treatment in lr.TREATMENTS:
        st = ox.state('titania_stack', treatment, activity, transmission)
        case = ox.middle_case('titania_stack', treatment)
        layers, summary = ox.middle_layers(case)
        col = ox.column(case, layers, summary, st['profile'], ox.Setup())
        moon = lunar_column(col, MOON_GM, col['planet'].radius_m, case.air().molar_mass)
        p, t = col['p_pa'], col['t_k']
        middle = case.mixing().profile(p, col['tropopause_pa'])
        row = dict(p_pa=p.tolist(), t_k=t.tolist(), z_km=(col['z_m'] / 1e3).tolist(),
                   moon_scaled_cm2_s=middle.tolist(), earth_cm2_s=(middle / middle_run.MOON_KZZ).tolist())
        for name, stretch in (('waves_stretched', STRETCH), ('waves_earth_wavelengths', 1.0)):
            raw, waves = diffusion(moon, spectrum(stretch), LAUNCH_OVER_TROPOPAUSE * col['tropopause_pa'])
            k = factor * raw * 1e4
            row[f'{name}_cm2_s'] = k.tolist()
            row[f'{name}_turned_back'] = sum(w['turns_at_pa'] is not None for w in waves) / len(waves)
            above = p < lr.BASE_PA
            row[f'{name}_homopause_pa'] = homopause(p[above], t[above], k[above])
        for name in ('moon_scaled', 'earth'):
            row[f'{name}_homopause_pa'] = homopause(p, t, np.asarray(row[f'{name}_cm2_s']))
        lunar[treatment] = row
        print(treatment, {k: (f'{v:.3g}' if isinstance(v, float) else v) for k, v in row.items()
                          if k.endswith('homopause_pa') or k.endswith('turned_back')}, flush=True)
    files = [Path(__file__), Path(ox.__file__)]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={str(f.relative_to(ROOT)): digest(f) for f in files},
                      constants=constants_used(files),
                      products={name: digest(ROOT / name) for name in ox.PRODUCTS}),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('Eddy diffusion in cm^2/s at the levels given. earth is the US Standard Atmosphere 1976 with the '
                      'waves calibrated to Earth\'s homopause; lunar[treatment] is the titania stack\'s column at the '
                      'standard level and the solar cycle\'s mean spectrum, with the waves\' mixing (wavelengths '
                      'stretched by the ratio of gravities, or Earth\'s), the middle atmosphere\'s Moon-scaled mixing '
                      'and Earth\'s for comparison, and the homopause each gives for atomic oxygen (the waves\' above '
                      'the 0.3 Pa base, where the oxygen step takes them).'),
        spectrum=dict(phase_speeds_m_s=PHASE_SPEEDS_M_S, phase_width_m_s=PHASE_WIDTH_M_S, wavelengths_km=WAVELENGTHS_KM,
                      launch_flux_pa=LAUNCH_FLUX_PA, launch_over_tropopause=LAUNCH_OVER_TROPOPAUSE, stretch=STRETCH),
        calibration=dict(factor=factor, raw_at_homopause_m2_s=raw_at, earth_homopause_km=EARTH_HOMOPAUSE_KM,
                         earth_homopause_k_cm2_s=EARTH_HOMOPAUSE_K_CM2_S),
        earth=out_earth, lunar=lunar)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
