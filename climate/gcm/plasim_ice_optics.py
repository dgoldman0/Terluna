#!/usr/bin/env python3
"""Ice optics for PlaSim's cold cloud: what the runner's --ice-optics patch (ICE_OPTICS_EDITS) gives a cloud layer
colder than -38 °C, mirrored here beside PlaSim's own liquid-water optics.

    climate/gcm/.venv/bin/python -m climate.gcm.plasim_ice_optics table
    climate/gcm/.venv/bin/python -m climate.gcm.plasim_ice_optics write

PlaSim gives every cloud the optics of liquid water: Stephens' (1978) fit of reflection to water path in sunlight, and
a mass absorption of 0.1 m²/g in the infrared. Below -38 °C (235.15 K) cloud water freezes by itself, so cloud there is
ice. The patch gives such cloud CAM3's ice optics, as RRTMG carries them: Ebert and Curry's (1992) extinction,
single-scattering albedo and asymmetry in sunlight (RRTMG's iceflag 1), and their grey infrared absorption,
0.005 + 1/r m²/g (iceflag 0). The effective radius r comes from the layer's temperature through Kristjánsson and
Mitchell's table for hexagonal columns, 180-274 K, as CAM3 and WRF's RRTMG carry it (reicalc). Below 180 K, where the
table ends, r follows the table's own trend and shrinks by 5.5% per kelvin. Ebert and Curry fitted crystals of 10 µm
and larger, and the table gives that size at 190 K, so colder cloud rests on their formulas carried to smaller crystals.

PlaSim's two sunlight bands take a lone ice layer's reflection and transmission from the atmosphere domain's
delta-scaled PIFM two-stream (shortwave.layer_properties): its visible band (below 0.75 µm) from Ebert and Curry's
visible band, its near-infrared band as the sum over their other four bands, each weighted by its share of the TSIS-1
sunlight beyond 0.75 µm. As in PlaSim, light that has passed cloud higher up counts as diffuse. The cloud's water path
stays PlaSim's own (the runner's cloud_water setting); only how that water meets light changes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-ice-optics/1'
ICE_BELOW_K = 235.15                    # -38 °C: cloud water freezes by itself below this
DIFFUSIVITY = 1.66                      # PlaSim's and RRTMG's factor from vertical to flux optical depth
LIQUID_KAPPA_M2_G = 0.1                 # PlaSim's acllwr
TSWR = (0.077, 0.065, 0.0055)           # PlaSim's tswr1-3, its fit of liquid cloud reflection
MU_FLOOR = 1e-4                         # the Sun's cosine is held above this at the terminator

# Kristjánsson and Mitchell's ice effective radius (µm) for hexagonal columns at 180, 181, ... 274 K, as CAM3 and
# WRF's RRTMG (module_ra_rrtmg_lw.F, retab) carry it.
RETAB = np.array([
    5.92779, 6.26422, 6.61973, 6.99539, 7.39234, 7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,
    10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,
    20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, 27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,
    31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, 34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,
    38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, 42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,
    50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, 65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,
    93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, 124.954, 130.630, 136.457, 142.446, 148.608, 154.956,
    161.503, 168.262, 175.248, 182.473, 189.952, 197.699, 205.728, 214.055, 222.694, 231.661, 240.971, 250.639])

# Ebert and Curry's ice optics as RRTMG's shortwave carries them (module_ra_rrtmg_sw.F, abari-fbari), for their bands
# above 14300, 7700-14300, 5300-7700, 4000-5300 and below 4000 cm-1: extinction per g/m2 of ice A + B/r, 1 - omega
# C + D r, asymmetry E + F r, with r in µm.
EC_A, EC_B = 3.448e-03, 2.431
EC_C = np.array([1.000e-05, 1.100e-04, 1.240e-02, 3.779e-02, 4.666e-01])
EC_D = np.array([0.000e+00, 1.405e-05, 6.867e-04, 1.284e-03, 2.050e-05])
EC_E = np.array([7.661e-01, 7.730e-01, 7.865e-01, 8.172e-01, 9.595e-01])
EC_F = np.array([5.851e-04, 5.665e-04, 7.204e-04, 7.463e-04, 1.076e-04])
NIR_EDGES_CM = (1e4 / 0.75, 7700.0, 5300.0, 4000.0)      # PlaSim's band 2 from 0.75 µm, then Ebert and Curry's
# Their four near-infrared bands' shares of the TSIS-1 sunlight beyond 0.75 µm (nir_weights), as the patch has them.
NIR_WEIGHTS = (0.6164, 0.2283, 0.0811, 0.0742)


def ice_radius(t):
    """Ice effective radius (µm) at temperature t (K): the table between whole kelvins, its own trend below 180 K."""
    t = np.asarray(t, float)
    i = np.clip(np.floor(t - 179.0).astype(int), 1, 94)
    c = t - (179.0 + i)
    inside = RETAB[i - 1] * (1.0 - c) + RETAB[i] * c
    below = RETAB[0] * np.exp(np.log(RETAB[1] / RETAB[0]) * (t - 180.0))
    return np.where(t < 180.0, below, inside)


def nir_weights():
    """Ebert and Curry's near-infrared bands' shares of the TSIS-1 sunlight in PlaSim's band 2 (beyond 0.75 µm); the
    last takes the sunlight beyond 5 µm as well."""
    from atmosphere.radiative_convective import shortwave as sw
    nu = np.linspace(2000.0, NIR_EDGES_CM[0], 200001)
    spec, beyond = sw.solar_spectrum(nu)
    edges = NIR_EDGES_CM + (2000.0,)
    parts = []
    for hi, lo in zip(edges[:-1], edges[1:]):
        m = (nu >= lo) & (nu <= hi)
        parts.append(float(np.trapezoid(spec[m], nu[m])))
    parts[-1] += float(beyond)
    return tuple(p / sum(parts) for p in parts)


def band_optics(r):
    """Per g/m2 extinction, single-scattering albedo and asymmetry in Ebert and Curry's five bands at radius r."""
    tau = np.full(5, EC_A + EC_B / r)                      # the same in every band
    omega = 1.0 - EC_C - EC_D * r
    g = np.minimum(EC_E + EC_F * r, 0.9999)
    return tau, omega, g


def pifm(tau, omega, g, mu):
    """A lone layer's direct-beam reflection and total transmission, and its diffuse reflection and transmission."""
    from atmosphere.radiative_convective.shortwave import layer_properties
    a = [np.atleast_2d(np.asarray(x, float)) for x in (tau, omega, g, max(mu, MU_FLOOR))]
    r_dif, t_dif, r_dir, t_dir, eb = layer_properties(*a)
    return float(r_dir[0, 0]), float(t_dir[0, 0] + eb[0, 0]), float(r_dif[0, 0]), float(t_dif[0, 0])


def ice_sw(water_g_m2, t, mu):
    """What the patch gives a fully cloudy ice layer under a clear sky: PlaSim's band 1 direct and diffuse reflection,
    and band 2's direct reflection and transmission, then diffuse."""
    r = float(ice_radius(t))
    tau, omega, g = band_optics(r)
    r1, _, r1s, _ = pifm(water_g_m2 * tau[0], omega[0], g[0], mu)
    band2 = np.zeros(4)
    for w, k in zip(NIR_WEIGHTS, range(1, 5)):
        band2 += w * np.array(pifm(water_g_m2 * tau[k], omega[k], g[k], mu))
    return dict(r1=r1, r1s=r1s, r2=band2[0], t2=band2[1], r2s=band2[2], t2s=band2[3])


def ice_kappa(t):
    """Ebert and Curry's grey infrared absorption of ice (m²/g) at the radius the temperature gives."""
    return 0.005 + 1.0 / ice_radius(t)


def emissivity(water_g_m2, kappa):
    return 1.0 - np.exp(-DIFFUSIVITY * kappa * water_g_m2)


def liquid_sw(water_g_m2, mu):
    """PlaSim's own cloud optics in sunlight (swr, Stephens 1978) for the same cases as ice_sw; its diffuse values are
    its direct formulas at a Sun cosine of 0.5."""
    tau = 2.0 * np.log10(water_g_m2 + 1.5) ** 3.9
    ln = np.log(1000.0 / tau)
    out = {}
    for key, m in (('', max(mu, MU_FLOOR)), ('s', 0.5)):
        b1 = TSWR[0] * np.sqrt(m)
        b2 = TSWR[1] * np.sqrt(m) / np.log(3.0 + 0.1 * tau)
        om = min(0.9999, 1.0 - TSWR[2] * m * m * ln)
        un = 1.0 - om
        uz = un + 2.0 * b2 * om
        u = np.sqrt(uz / un)
        ex = np.exp(min(25.0, tau * np.sqrt(uz * un) / m))
        rr = (u + 1.0) ** 2 * ex - (u - 1.0) ** 2 / ex
        out['r1' + key] = 1.0 - 1.0 / (1.0 + b1 * tau / m)
        out['t2' + key] = 4.0 * u / rr
        out['r2' + key] = (u * u - 1.0) / rr * (ex - 1.0 / ex)
    return out


def disk_mean(f):
    """Mean of a reflection over the sunlit disk (weights 2 mu0)."""
    from atmosphere.radiative_convective.shortwave import disk_nodes
    mu, w = disk_nodes(6)
    return float(sum(2.0 * m * wt * f(m) for m, wt in zip(mu, w)))


def comparison(waters=(0.5, 1.5, 4.7, 15.0), temps=(165.0, 180.0, 189.0, 200.0, 217.0, 230.0)):
    """Disk-mean reflection in each band and infrared emissivity of a fully cloudy layer, ice at each temperature
    against PlaSim's liquid optics, for water paths in g/m2."""
    rows = []
    for wp in waters:
        liq = dict(water_g_m2=wp, visible_reflection=disk_mean(lambda m: liquid_sw(wp, m)['r1']),
                   near_ir_reflection=disk_mean(lambda m: liquid_sw(wp, m)['r2']),
                   near_ir_absorption=disk_mean(lambda m: 1.0 - liquid_sw(wp, m)['r2'] - liquid_sw(wp, m)['t2']),
                   emissivity=float(emissivity(wp, LIQUID_KAPPA_M2_G)))
        ice = []
        for t in temps:
            ice.append(dict(t_k=t, radius_um=float(ice_radius(t)),
                            visible_reflection=disk_mean(lambda m: ice_sw(wp, t, m)['r1']),
                            near_ir_reflection=disk_mean(lambda m: ice_sw(wp, t, m)['r2']),
                            near_ir_absorption=disk_mean(lambda m: 1.0 - ice_sw(wp, t, m)['r2'] - ice_sw(wp, t, m)['t2']),
                            emissivity=float(emissivity(wp, ice_kappa(t)))))
        rows.append(dict(liquid=liq, ice=ice))
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('action', choices=('table', 'write'))
    args = parser.parse_args(argv)
    rows = comparison()
    for row in rows:
        liq = row['liquid']
        print(f"{liq['water_g_m2']:5.1f} g/m2  liquid: visible {liq['visible_reflection']:.4f}, near-IR "
              f"{liq['near_ir_reflection']:.4f} reflected {liq['near_ir_absorption']:.4f} absorbed, emissivity "
              f"{liq['emissivity']:.3f}")
        for ice in row['ice']:
            print(f"        ice at {ice['t_k']:5.1f} K, r {ice['radius_um']:5.1f} um: visible {ice['visible_reflection']:.4f}, "
                  f"near-IR {ice['near_ir_reflection']:.4f} / {ice['near_ir_absorption']:.4f}, emissivity "
                  f"{ice['emissivity']:.3f}")
    if args.action == 'write':
        weights = nir_weights()
        out = dict(schema=SCHEMA, ice_below_k=ICE_BELOW_K, near_ir_weights=[round(w, 4) for w in weights],
                   radius_below_180k_per_k=float(np.log(RETAB[1] / RETAB[0])), comparison=rows,
                   evidence=('Ebert and Curry (1992) ice optics and Kristjánsson and Mitchell\'s radius by temperature as '
                             'CAM3 and RRTMG carry them, in the atmosphere domain\'s delta-scaled PIFM two-stream, '
                             'against PlaSim\'s liquid-water optics on the same water paths; the radius follows the '
                             'table\'s trend below 180 K and the optics their formulas below 10 µm.'),
                   reading_rule=('A fully cloudy layer alone over a black ground: reflection and absorption are means '
                                 'over the sunlit disk (weights 2 mu0) in PlaSim\'s visible band (below 0.75 µm) and '
                                 'near-infrared band, and emissivity is 1 - exp(-1.66 kappa W) for water path W.'),
                   producer=dict(domain='climate', files={f'gcm/{Path(__file__).name}':
                                                          hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}))
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / 'ice_optics.json').write_text(json.dumps(out, indent=1) + '\n')
        print('near-infrared weights', out['near_ir_weights'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
