"""The cloud water PlaSim's radiation sees, layer by layer, recomputed from a year of ExoPlaSim output.

    climate/gcm/.venv/bin/python -m climate.gcm.cloud_water climate/gcm/runs/A/model/MOST.00039.nc

PlaSim carries no cloud water: its condensation and convection schemes rain out excess vapour at once.
For the clouds' radiation it diagnoses their water separately (rainmod.f90, mkclouds) with CCM3's
formula (Kiehl et al. 1996, NCAR/TN-420+STR, pp. 49-50): a density of 0.21 g/m3 at the ground falling
off as exp(-z/h), h = 700 ln(1 + PW) m, where z is the layer's height and PW the precipitable water in
kg/m2, a scale fitted to Earth. This module mirrors the formula and the runner's cloud_water choices
(exoplasim_run.CLOUD_WATER) on a year of output. The output holds 3-day means of temperature, humidity
and cloud cover, so the numbers show how the formula behaves, not the model's instantaneous values.
"""
import sys
import numpy as np

from climate.gcm.exoplasim_run import CLOUD_WATER, cloud_water_namelist

SURFACE_DENSITY_KG_M3 = 0.00021     # CCM3's cloud water density at the ground, as in PlaSim
SCALE_M = 700.0                     # CCM3's scale height per unit ln(1 + PW), as in PlaSim
MAX_PATH_G_M2 = 1000.0              # PlaSim's radiation caps a layer's water path here
GASCON = 287.0                      # PlaSim's dry-air gas constant in the formula


def half_levels(sigma):
    """PlaSim's layer boundaries (0 at the top to 1 at the ground) from its full levels, which lie midway
    between them. (The output's 'levp' holds midpoints between full levels instead.)"""
    half = [0.0]
    for s in sigma:
        half.append(2.0 * float(s) - half[-1])
    if abs(half[-1] - 1.0) > 2e-3:
        raise ValueError('the full levels do not lie midway between boundaries running from 0 to 1')
    half[-1] = 1.0
    return np.array(half)


def cloud_water_path(ta, q, ps, sigma, gravity, choice='plasim'):
    """In-cloud liquid water path of each layer (g/m2) as PlaSim's radiation computes it, for
    temperature (K) and specific humidity shaped (time, level, lat, lon) from the top down and surface
    pressure (Pa) shaped (time, lat, lon)."""
    g_ref, scale = cloud_water_namelist(choice, gravity)
    g = g_ref if g_ref > 0 else gravity
    sigma = np.asarray(sigma, dtype=float)
    half = half_levels(sigma)
    dsig = np.diff(half)[None, :, None, None]
    pw = (q * dsig).sum(axis=1) * ps / gravity
    h = SCALE_M * np.log1p(pw * gravity / g)
    z = np.empty_like(ta)
    top = np.zeros_like(ps)
    for j in range(len(sigma) - 1, 0, -1):
        dz = ta[:, j] * GASCON / g * np.log(half[j + 1] / half[j])
        z[:, j] = top + 0.5 * dz
        top = top + dz
    z[:, 0] = top + 0.5 * ta[:, 0] * GASCON / g * np.log(half[1] / sigma[0])
    mixing_ratio = scale * SURFACE_DENSITY_KG_M3 * np.exp(-z / h[:, None]) * GASCON * ta / (sigma[None, :, None, None] * ps[:, None])
    return np.minimum(MAX_PATH_G_M2, 1000.0 * mixing_ratio * ps[:, None] / gravity * dsig), z


def optical_depth(path_g_m2):
    """PlaSim's cloud optical depth for a layer's water path (Stephens 1978 fit, radmod.f90)."""
    return 2.0 * np.log10(path_g_m2 + 1.5) ** 3.9


def main(argv=None) -> int:
    import json
    import netCDF4
    from climate.gcm.exoplasim_run import PRODUCTS
    path = (argv or sys.argv[1:])[0]
    gravity = json.loads((PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']['gravity_m_s2']
    with netCDF4.Dataset(path) as d:
        ta, q, cover = (np.asarray(d[k][:], dtype=float) for k in ('ta', 'hus', 'cl'))
        ps = np.asarray(d['ps'][:], dtype=float) * 100.0
        sigma = np.asarray(d['lev'][:], dtype=float)
        lat = np.asarray(d['lat'][:], dtype=float)
    _, weights = np.polynomial.legendre.leggauss(lat.size)            # Gaussian weights of the rows
    area = (weights[np.argsort(np.argsort(np.sin(np.radians(lat))))][:, None] * np.ones(ta.shape[-1]))[None]
    def mean(x, w=None):
        w = area * (np.ones_like(x) if w is None else w)
        return float((x * w).sum() / w.sum())
    results = {c: cloud_water_path(ta, q, ps, sigma, gravity, c) for c in CLOUD_WATER}
    print(f'{path}: in-cloud water path (g/m2) and optical depth of each layer, weighted by cloud cover')
    print('  sigma  cover  height km' + ''.join(f'  {c:>18}' for c in CLOUD_WATER))
    for j, s in enumerate(sigma):
        c = cover[:, j]
        cells = ''.join(f'  {mean(results[k][0][:, j], c):9.1f} {mean(optical_depth(results[k][0][:, j]), c):7.1f} ' for k in CLOUD_WATER)
        print(f'  {s:.3f}  {mean(c):.3f}  {mean(results["plasim"][1][:, j]) / 1000:8.1f} {cells}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
