#!/usr/bin/env python3
"""The GCM's cold cloud at lunar gravity: factors on its water and its crystals' radius for each of its top three
layers, from slabs of rising air on the Moon and on Earth (atmosphere/snow/cirrus.py).

    climate/gcm/.venv/bin/python -m climate.gcm.cold_cloud A28R_301k_s0897 --years 5-14    # about 20 minutes, 4 workers

PlaSim's cloud water at the design case is Earth's in-cloud water path at each pressure (the runner's cloud_water
'earth_path', CCM3's formula counted as on Earth), and the ice optics take Kristjánsson and Mitchell's radii, measured
in Earth's cirrus. Both hold Earth's gravity. Each top layer's bounds come from the run's sigma levels and mean surface
pressure over the years, its bottom temperature from the run's mean temperatures interpolated in ln sigma. A slab of
air rises through the layer from ice saturation at its bottom, on the Moon and on Earth at the same pressure velocity
omega and slab depth: the same cooling and condensation per kilogram, the dynamics equal in pressure. The top layer's
slab stops at its full level, as the air above it warms upward and holds little ascent. The factors are the Moon's
in-cloud ice path (the slab's mean ice mixing ratio times the layer's pressure thickness over g) and effective radius
over Earth's.

The central case for each layer takes the omega at which Earth's slab of 4 kPa holds the GCM's own water path there,
interpolated in log omega over the grid, so that the factors start from the dynamics PlaSim's Earth values stand
for; where no omega on the grid gives that water (the top layer, whose formula water exceeds what its air holds at
ice saturation), the layer below's. The range takes every slab of 2-8 kPa at omegas within a factor of three of the
central one, and the deposition coefficient (0.1 and 1 against 0.5) and the haze's number (a fifth and five times) at
it. The first-order change these factors make at the top of the air and the ground comes from high_cloud.effect, with
the run's own cloud. Writes results/gcm/cold_cloud_<run>.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-cold-cloud-factors/1'
LAYERS = 3
OMEGAS = (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3)          # Pa/s: Earth's 0.3-10 cm/s near 200 hPa
SLABS = (2000.0, 4000.0, 8000.0)                            # Pa: Earth's 0.6-2.5 km of cirrus near 200 hPa
CENTRAL_SLAB = 4000.0
ALPHAS = (0.1, 1.0)                                         # beside the central 0.5
HAZE_SCALES = (0.2, 5.0)
CODE = ('climate/gcm/cold_cloud.py', 'climate/gcm/high_cloud.py', 'climate/gcm/plasim_ice_optics.py',
        'atmosphere/snow/cirrus.py', 'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py')


def run_column(run, years):
    """The run's sigma levels, its area-weighted mean surface pressure and temperature profile over the years."""
    import netCDF4
    ps_sum = area_sum = 0.0
    t_sum = None
    for year in years:
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            sigma = np.asarray(d['lev'][:], float)
            lat = np.asarray(d['lat'][:], float)
            ps = np.asarray(d['ps'][:], float) * 100.0
            ta = np.asarray(d['ta'][:], float)
        w = np.cos(np.radians(lat))[None, :, None] * np.ones(ps.shape)
        ps_sum += float((ps * w).sum())
        area_sum += float(w.sum())
        t = np.array([(ta[:, j] * w).sum() for j in range(sigma.size)])
        t_sum = t if t_sum is None else t_sum + t
    return dict(sigma=sigma.tolist(), surface_pa=ps_sum / area_sum, t_k=(t_sum / area_sum).tolist())


def layer_bounds(column):
    """Each top layer's bottom and top sigma, its pressure thickness, the slab's bottom temperature and end."""
    from climate.gcm.cloud_water import half_levels
    sigma = np.asarray(column['sigma'])
    half = half_levels(sigma)
    t = np.asarray(column['t_k'])
    ps = column['surface_pa']
    rows = []
    for j in range(LAYERS):
        bottom = half[j + 1]
        k = j if bottom < sigma[j] else j + 1                   # the full levels either side of the bottom
        lo, hi = (k - 1, k) if k >= 1 else (0, 1)
        t_bottom = t[lo] + (t[hi] - t[lo]) * math.log(bottom / sigma[lo]) / math.log(sigma[hi] / sigma[lo])
        top = sigma[j] if j == 0 else half[j]                  # the top layer's slab stops at its full level
        rows.append(dict(layer=j + 1, sigma_bottom=float(bottom), sigma_top=float(half[j]), p_bottom_pa=bottom * ps,
                         p_end_pa=top * ps, layer_pa=(half[j + 1] - half[j]) * ps, t_bottom_k=float(t_bottom)))
    return rows


def _worlds():
    from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
    from atmosphere.radiative_convective.thermodynamics import DryAir, earthlike_air
    from atmosphere.snow import cirrus
    return dict(moon=cirrus.Air(MOON_SURFACE_GRAVITY, earthlike_air(121590.0, 400.0)),
                earth=cirrus.Air(STANDARD_GRAVITY, DryAir(101325.0)))


def _slab(job):
    from atmosphere.snow import cirrus
    world, layer, omega, slab, alpha, haze = job
    w = _worlds()[world]
    r = cirrus.slab(w, layer['t_bottom_k'], layer['p_bottom_pa'], layer['p_end_pa'], omega, slab, alpha=alpha,
                    haze_per_kg=cirrus.HAZE_PER_KG * haze)
    r.update(world=world, layer=layer['layer'], haze_scale=haze,
             path_g_m2=1e3 * r['ice_kg_kg'] * layer['layer_pa'] / w.gravity)
    return r


def _interpolate(xs, ys, x):
    lx, ly = np.log(np.asarray(xs)), np.log(np.asarray(ys))
    return float(np.exp(np.interp(math.log(x), lx, ly)))


def central_omega(earth_paths, target):
    """The omega (log-interpolated) at which Earth's slab holds the target path, or None when no omega reaches it."""
    paths = np.asarray(earth_paths)
    for i in range(len(OMEGAS) - 1):
        a, b = paths[i], paths[i + 1]
        if min(a, b) <= target <= max(a, b) and a > 0 and b > 0 and a != b:
            f = math.log(target / a) / math.log(b / a)
            return float(math.exp(math.log(OMEGAS[i]) + f * math.log(OMEGAS[i + 1] / OMEGAS[i])))
    return None


def analyse(run, years, workers=4):
    from climate.gcm import high_cloud as hc
    from illumination.surface_light import model as surface
    cloud_product = RESULTS / f'high_cloud_{run}.json'
    high = json.loads(cloud_product.read_text())
    column = run_column(run, years)
    layers = layer_bounds(column)
    jobs = [(w, layer, om, sl, 0.5, 1.0) for layer in layers for om in OMEGAS for sl in SLABS for w in ('earth', 'moon')]
    with ProcessPoolExecutor(workers) as ex:
        grid = list(ex.map(_slab, jobs))
    key = {(r['world'], r['layer'], r['omega_pa_s'], r['slab_pa']): r for r in grid}
    centrals = []
    for layer, water in zip(layers, high['cloud']['water_g_m2']):
        earth = [key[('earth', layer['layer'], om, CENTRAL_SLAB)]['path_g_m2'] for om in OMEGAS]
        centrals.append(central_omega(earth, water))
    for j in range(LAYERS):                                     # a layer no omega fits takes the one below's
        if centrals[j] is None:
            below = next((c for c in centrals[j + 1:] if c is not None), None)
            centrals[j] = below
    sens = [(w, layer, om, CENTRAL_SLAB, a, 1.0) for layer, om in zip(layers, centrals) for a in ALPHAS
            for w in ('earth', 'moon')]
    sens += [(w, layer, om, CENTRAL_SLAB, 0.5, h) for layer, om in zip(layers, centrals) for h in HAZE_SCALES
             for w in ('earth', 'moon')]
    sens += [(w, layer, om, CENTRAL_SLAB, 0.5, 1.0) for layer, om in zip(layers, centrals) for w in ('earth', 'moon')]
    with ProcessPoolExecutor(workers) as ex:
        extra = list(ex.map(_slab, sens))
    rows = []
    from climate.gcm import plasim_ice_optics as io
    for j, (layer, om) in enumerate(zip(layers, centrals)):
        n = layer['layer']
        pairs = lambda rs: [(e, m) for e, m in zip(rs[0::2], rs[1::2])]
        at = [r for r in extra if r['layer'] == n]
        centre_e, centre_m = at[-2], at[-1]
        factors = dict(water=centre_m['path_g_m2'] / centre_e['path_g_m2'],
                       radius=centre_m['effective_radius_um'] / centre_e['effective_radius_um'])
        spread = [(key[('earth', n, o, s)], key[('moon', n, o, s)]) for o in OMEGAS for s in SLABS
                  if om / 3.0 <= o <= om * 3.0]
        spread += pairs(at[:-2])
        water = [m['path_g_m2'] / e['path_g_m2'] for e, m in spread]
        radius = [m['effective_radius_um'] / e['effective_radius_um'] for e, m in spread]
        rows.append(dict(layer, gcm_water_g_m2=high['cloud']['water_g_m2'][j], gcm_t_k=high['cloud']['t_k'][j],
                         gcm_radius_um=float(io.ice_radius(high['cloud']['t_k'][j])), central_omega_pa_s=om,
                         central=dict(earth=centre_e, moon=centre_m), factors=factors,
                         water_range=[min(water), max(water)], radius_range=[min(radius), max(radius)]))
    shield, _ = surface.sunlight('design')
    visible = hc.visible_share(shield)
    first_order = hc.effect(high['cloud'], visible, water_scale=[r['factors']['water'] for r in rows],
                            radius_scale=[r['factors']['radius'] for r in rows])
    return dict(schema=SCHEMA, run=run, years=list(years), column=column, layers=rows, grid=grid, sensitivity=extra,
                first_order=first_order,
                inputs={f'climate/results/gcm/{cloud_product.name}':
                        dict(path=f'climate/results/gcm/{cloud_product.name}',
                             sha256=hashlib.sha256(cloud_product.read_bytes()).hexdigest()[:16])},
                evidence=('Offline slab calculations of the cold cloud on both worlds at equal dynamics in pressure, '
                          'with homogeneous freezing, kinetically limited growth and fall-out; no aggregation, '
                          'heterogeneous nuclei, waves or ice falling in from above. Factors on the GCM\'s Earth '
                          'values, and a first-order radiative estimate from the run\'s own cloud; not a run.'),
                reading_rule=('layers holds each top layer (1 the top) with its bounds, the central omega, the slabs '
                              'there on both worlds (central), the factors (Moon over Earth: in-cloud ice path and '
                              'effective radius) and their ranges over slabs of 2-8 kPa at omegas within a factor of '
                              'three, the deposition coefficient and the haze. first_order is high_cloud.effect with '
                              'the central factors layer by layer: W/m2 at the top of the air, positive warming, and '
                              'the ground through the 301 K branches\' feedback.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('--years', default='5-14', help='first-last')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args(argv)
    first, last = (int(x) for x in args.years.split('-'))
    out = analyse(args.run, range(first, last + 1), args.workers)
    out['producer'] = dict(domain='climate', files={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
                                                    for p in CODE})
    path = RESULTS / f'cold_cloud_{args.run}.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    for r in out['layers']:
        e, m = r['central']['earth'], r['central']['moon']
        print(f"layer {r['layer']} ({r['p_end_pa']/1e3:.1f}-{r['p_bottom_pa']/1e3:.1f} kPa, from {r['t_bottom_k']:.1f} K): "
              f"omega {r['central_omega_pa_s']:.3f} Pa/s; Earth {e['path_g_m2']:.2f} g/m2, {e['effective_radius_um']:.1f} um "
              f"(GCM {r['gcm_water_g_m2']:.2f}, {r['gcm_radius_um']:.1f}); Moon {m['path_g_m2']:.1f} g/m2, "
              f"{m['effective_radius_um']:.1f} um; water x{r['factors']['water']:.1f} ({r['water_range'][0]:.1f}-"
              f"{r['water_range'][1]:.1f}), radius x{r['factors']['radius']:.2f} ({r['radius_range'][0]:.2f}-"
              f"{r['radius_range'][1]:.2f})")
    f = out['first_order']
    print(f"first order: {f['net_w_m2']:+.1f} W/m2 at the top of the air; ground "
          + ' to '.join(f'{k:+.1f}' for k in sorted(f['ground_k'])) + ' K')
    print(f'wrote {path.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
