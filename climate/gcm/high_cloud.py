#!/usr/bin/env python3
"""How far the GCM's high cloud could move the design climate through its water and its crystals' size: a first-order
range from a run's own years.

    climate/gcm/.venv/bin/python -m climate.gcm.high_cloud A28R_301k_s0897 --years 5-14

PlaSim diagnoses its cloud water with Earth's CCM3 formula (cloud_water.py) and, with --ice-optics, gives cloud colder
than -38 °C the ice optics of plasim_ice_optics.py, the crystals' radius set by the temperature through Kristjánsson and
Mitchell's table. Neither has a lunar measure: where CM1's 45° N ring made ice of its own it held 1.5-3.5 times
PlaSim's water, and in a sixth of Earth's gravity crystals can grow larger before they fall. For the run's top three
layers this takes, over the given years, the cover-weighted in-cloud water path and temperature, a fully cloudy
layer's disk-mean reflection of the shielded sunlight and its infrared emissivity as the patch gives them, and their
change when the water or the radius is scaled. A change dR in a layer's reflection changes the sunlight the planet
reflects by cover dR S (1 - a)^2, with S the sunlight arriving at the top of the air and a the planetary albedo below,
as for a thin layer over a reflecting column; a change de in its emissivity cuts the outgoing infrared by
cover de (F - sigma T^4), with F the clear-sky outgoing infrared and T the layer's temperature. The layers count apart,
without overlap, and the changes are first order: the cloud's amount, the vapour and the circulation do not respond,
and the sunlight is the mean over day and night. The ground's change is the net over the 301 K branches' own feedback
(compare --settle). Writes results/gcm/high_cloud_<run>.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.gcm import plasim_ice_optics as io

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-high-cloud-range/1'
SIGMA_SB = 5.670374419e-8
LAYERS = 3                                   # the top three, above about 63 km
WATER_SCALES = (0.5, 1.5, 2.5, 3.5)          # CM1's own ice held 1.5-3.5 times PlaSim's water at 63-107 km
RADIUS_SCALES = (0.5, 2.0)                   # smaller crystals in the coldest air, larger at lunar gravity
FEEDBACK_W_M2_K = (1.20, 1.56)               # the 301 K branches with and without the clear share (compare --settle)


def visible_share(shield=None):
    """The share of the sunlight below 0.75 µm, PlaSim's first band, for a shield's transmission (None: none)."""
    from atmosphere.radiative_convective import shortwave as sw
    nu = np.linspace(2000.0, 50000.0, 480001)
    spec, beyond = sw.solar_spectrum(nu, shield=shield)
    total = float(np.trapezoid(spec, nu)) + float(beyond)
    m = nu >= io.NIR_EDGES_CM[0]
    return float(np.trapezoid(spec[m], nu[m])) / total


def layer_optics(water_g_m2, radius_um, visible, nir_weights=io.SHIELD_NIR_WEIGHTS):
    """A fully cloudy ice layer's disk-mean reflection of direct sunlight, both bands weighted by their shares, and
    its infrared emissivity, for a water path and a crystal radius."""
    tau, omega, g = io.band_optics(radius_um)

    def reflection(mu):
        r1 = io.pifm(water_g_m2 * tau[0], omega[0], g[0], mu)[0]
        r2 = sum(w * io.pifm(water_g_m2 * tau[k], omega[k], g[k], mu)[0] for w, k in zip(nir_weights, range(1, 5)))
        return visible * r1 + (1.0 - visible) * r2
    return io.disk_mean(reflection), float(io.emissivity(water_g_m2, 0.005 + 1.0 / radius_um))


def run_cloud(run, years):
    """The top layers' cover, cover-weighted in-cloud water path and temperature over the years, area-weighted, with the
    sunlight at the top, the planetary albedo and the clear-sky outgoing infrared."""
    import netCDF4
    from climate.gcm.cloud_water import cloud_water_path
    from climate.gcm.exoplasim_run import PRODUCTS
    gravity = json.loads((PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']['gravity_m_s2']
    cfg = json.loads((HERE / 'runs' / run / 'progress.json').read_text())['configuration']
    choice = cfg.get('cloud_water', 'earth_path')
    sums = dict(cover=np.zeros(LAYERS), area=0.0, water=np.zeros(LAYERS), t=np.zeros(LAYERS), sun_in=0.0, sun_out=0.0,
                olr_clear=0.0)
    for year in years:
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            ta, q, cl = (np.asarray(d[k][:], float) for k in ('ta', 'hus', 'cl'))
            ps = np.asarray(d['ps'][:], float) * 100.0
            sigma = np.asarray(d['lev'][:], float)
            lat = np.asarray(d['lat'][:], float)
            rst, rsut, rltcs = (np.asarray(d[k][:], float) for k in ('rst', 'rsut', 'rltcs'))
        area = np.cos(np.radians(lat))[None, :, None] * np.ones(ta.shape[-1])
        water, _ = cloud_water_path(ta, q, ps, sigma, gravity, choice)
        for j in range(LAYERS):
            w = cl[:, j] * area
            sums['cover'][j] += w.sum()
            sums['water'][j] += (water[:, j] * w).sum()
            sums['t'][j] += (ta[:, j] * w).sum()
        a = np.broadcast_to(area, rst.shape)
        sums['area'] += a.sum()
        sums['sun_in'] += ((rst - rsut) * a).sum()          # net less outgoing (negative upward): what arrives
        sums['sun_out'] += (-rsut * a).sum()
        sums['olr_clear'] += (-rltcs * a).sum()
    return dict(cover=(sums['cover'] / sums['area']).tolist(), water_g_m2=(sums['water'] / sums['cover']).tolist(),
                t_k=(sums['t'] / sums['cover']).tolist(), sunlight_w_m2=sums['sun_in'] / sums['area'],
                planetary_albedo=sums['sun_out'] / sums['sun_in'], olr_clear_w_m2=sums['olr_clear'] / sums['area'],
                cloud_water=choice, sigma=sigma[:LAYERS].tolist())


def effect(cloud, visible, water_scale=1.0, radius_scale=1.0):
    """The change at the top of the air (W/m2, positive warms) when the top layers' water and crystal radius are
    scaled, all alike or each by its own (a sequence, the top first), layer by layer and in all."""
    rows = []
    layers = len(cloud['cover'])
    waters = list(water_scale) if np.ndim(water_scale) else [water_scale] * layers
    radii = list(radius_scale) if np.ndim(radius_scale) else [radius_scale] * layers
    for c, w, t, ws, rs in zip(cloud['cover'], cloud['water_g_m2'], cloud['t_k'], waters, radii):
        r = float(io.ice_radius(t))
        r0, e0 = layer_optics(w, r, visible)
        r1, e1 = layer_optics(w * ws, r * rs, visible)
        sun = c * (r1 - r0) * cloud['sunlight_w_m2'] * (1.0 - cloud['planetary_albedo']) ** 2
        infrared = c * (e1 - e0) * (cloud['olr_clear_w_m2'] - SIGMA_SB * t ** 4)
        rows.append(dict(radius_um=r, reflection=r0, emissivity=e0, scaled_reflection=r1, scaled_emissivity=e1,
                         sunlight_w_m2=-sun, infrared_w_m2=infrared, net_w_m2=infrared - sun))
    net = sum(x['net_w_m2'] for x in rows)
    return dict(water_scale=water_scale, radius_scale=radius_scale, layers=rows, net_w_m2=net,
                ground_k=[net / f for f in FEEDBACK_W_M2_K])


def analyse(run, years):
    from illumination.surface_light import model as surface
    shield, _ = surface.sunlight('design')
    visible = visible_share(shield)
    cloud = run_cloud(run, years)
    cases = [effect(cloud, visible, water_scale=s) for s in WATER_SCALES]
    cases += [effect(cloud, visible, radius_scale=s) for s in RADIUS_SCALES]
    cases += [effect(cloud, visible, water_scale=min(WATER_SCALES), radius_scale=max(RADIUS_SCALES)),
              effect(cloud, visible, water_scale=max(WATER_SCALES), radius_scale=min(RADIUS_SCALES))]
    return dict(schema=SCHEMA, run=run, years=list(years), visible_share=visible, cloud=cloud,
                base=effect(cloud, visible)['layers'], cases=cases, feedback_w_m2_k=list(FEEDBACK_W_M2_K),
                evidence=('First-order radiative changes of the GCM\'s own top three layers of cloud, its water and '
                          'crystal radius scaled through the ice optics the GCM runs with, from the run\'s own means; '
                          'no cloud, vapour or circulation response, no overlap, day and night averaged. A range for '
                          'judging the high cloud, not a run.'),
                reading_rule=('cloud holds the top layers\' area-weighted cover and their cover-weighted in-cloud '
                              'water path and temperature over the years, the top first. Each case scales every top '
                              'layer\'s water or radius; sunlight_w_m2 and infrared_w_m2 are its changes at the top of '
                              'the air, positive warming, and ground_k the net over each feedback, K.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('--years', default='5-14', help='first-last')
    args = parser.parse_args(argv)
    first, last = (int(x) for x in args.years.split('-'))
    out = analyse(args.run, range(first, last + 1))
    out['producer'] = dict(domain='climate', files={f'gcm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'plasim_ice_optics.py',
                                                              HERE / 'cloud_water.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'high_cloud_{args.run}.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    c = out['cloud']
    print(f"visible share {out['visible_share']:.3f}; sunlight {c['sunlight_w_m2']:.1f} W/m2, albedo "
          f"{c['planetary_albedo']:.3f}, clear-sky OLR {c['olr_clear_w_m2']:.1f} W/m2")
    for j, b in enumerate(out['base']):
        print(f"layer {j + 1}: cover {c['cover'][j]:.3f}, water {c['water_g_m2'][j]:.2f} g/m2, {c['t_k'][j]:.1f} K, "
              f"r {b['radius_um']:.1f} um: reflection {b['reflection']:.4f}, emissivity {b['emissivity']:.3f}")
    for case in out['cases']:
        print(f"water x{case['water_scale']}, radius x{case['radius_scale']}: net {case['net_w_m2']:+.2f} W/m2 "
              f"(sunlight {sum(x['sunlight_w_m2'] for x in case['layers']):+.2f}, infrared "
              f"{sum(x['infrared_w_m2'] for x in case['layers']):+.2f}); ground "
              + ' to '.join(f'{k:+.2f}' for k in sorted(case['ground_k'])) + ' K')
    print(f'wrote {path.relative_to(HERE.parent.parent)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
