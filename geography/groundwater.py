"""Water the porous crust would take up: a first estimate of the groundwater term in the water inventory.

    python -m geography.groundwater            # results/groundwater.json

GRAIL showed the lunar highland crust is porous: a bulk density of 2550 kg/m3, implying an average
porosity of 12% to depths of at least a few kilometres (Wieczorek et al. 2013), falling to about 4% at
20 km (as summarised by Wiggins et al. 2022). Water delivered to the Moon fills that pore space where it
can reach it, beneath the seas and below the water table on land, in addition to the standing water of
the seas and lakes. This module fits an exponential porosity profile, phi(z) = phi0 exp(-z/d), to the
two constraints, integrates the pore volume down to a range of saturation depths, and sets it beside
the seas (the atlas) and the rain-fed lakes (the drainage estimate) as global equivalent layers and
masses.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.geography.groundwater/1'
WATER_DENSITY = 1000.0                    # kg/m3, as the other geography products
MEAN_POROSITY_UPPER = (0.12, 4000.0)      # average porosity and the depth it averages over (m)
POROSITY_AT_DEPTH = (0.04, 20000.0)       # porosity and its depth (m)
SATURATION_DEPTHS_M = (1000.0, 2000.0, 4000.0, 10000.0, 20000.0)
EVIDENCE = ('A first estimate from two published constraints on the highland crust\'s porosity (GRAIL): an '
            'exponential profile through them, taken to hold everywhere and to be filled completely below the '
            'saturation depth chosen. Not modelled: the maria\'s basalt fill, which is less porous; pores that are '
            'closed or unconnected; water bound into new minerals; the water table\'s shape on land; crustal '
            'loading; and how fast the water gets there, which depends on a permeability unknown within orders of '
            'magnitude.')
READING_RULE = ('Read the groundwater rows as the extra water the crust would hold once saturated to that depth, '
                'beneath seas and land alike, and the totals as the inventory to deliver for the standing water '
                'to stay at the scenario share. The saturation depth is a question of time and permeability, not '
                'settled here.')


def fit_profile(mean_upper=MEAN_POROSITY_UPPER, at_depth=POROSITY_AT_DEPTH):
    """(phi0, d) of phi = phi0 exp(-z/d) whose mean over the upper layer and value at depth match."""
    (mean, depth_mean), (phi_deep, depth_deep) = mean_upper, at_depth

    def mismatch(d):                       # phi0 from the deep value, then its upper mean against the target
        phi0 = phi_deep * math.exp(depth_deep / d)
        return phi0 * d * (1.0 - math.exp(-depth_mean / d)) / depth_mean - mean

    lo, hi = 1000.0, 200000.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mismatch(mid) > 0:              # too porous near the top: a longer decay lowers phi0
            lo = mid
        else:
            hi = mid
    d = 0.5 * (lo + hi)
    return phi_deep * math.exp(depth_deep / d), d


def pore_column_m(depth_m, phi0, d):
    """Pore volume per unit area (m of water) from the surface to depth_m."""
    return phi0 * d * (1.0 - math.exp(-depth_m / d))


def results():
    atlas = json.loads((HERE / 'results' / 'atlas.json').read_text())
    drainage = json.loads((HERE / 'results' / 'drainage.json').read_text())
    area = 4 * math.pi * MOON_RADIUS ** 2
    phi0, d = fit_profile()
    seas_kg = atlas['water_mass_kg']
    lakes_kg = drainage['lakes']['volume_km3'] * 1e9 * WATER_DENSITY
    layer = lambda kg: kg / WATER_DENSITY / area
    rows = []
    for depth in SATURATION_DEPTHS_M:
        column = pore_column_m(depth, phi0, d)
        crust_kg = column * area * WATER_DENSITY
        total = seas_kg + lakes_kg + crust_kg
        rows.append(dict(saturation_depth_m=depth, crust_layer_m=round(column, 1), crust_kg=float(f'{crust_kg:.3e}'),
                         total_kg=float(f'{total:.3e}'), total_over_seas=round(total / seas_kg, 2)))
    return dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        porosity_profile=dict(phi0=round(phi0, 4), e_folding_m=round(d), constraints=dict(
            mean_upper=dict(porosity=MEAN_POROSITY_UPPER[0], over_m=MEAN_POROSITY_UPPER[1]),
            at_depth=dict(porosity=POROSITY_AT_DEPTH[0], depth_m=POROSITY_AT_DEPTH[1])),
            sources=['Wieczorek et al. 2013, Science 339, 671 (12% to at least a few km; bulk density 2550 kg/m3)',
                     'Wiggins et al. 2022, Nature Communications, doi:10.1038/s41467-022-32445-3 (12% over the upper '
                     '~4 km of the highland crust, 4% at 20 km)']),
        share=atlas.get('share'),
        standing_water=dict(seas_kg=float(f'{seas_kg:.3e}'), seas_layer_m=round(layer(seas_kg), 1),
                            rain_fed_lakes_kg=float(f'{lakes_kg:.3e}'), lakes_layer_m=round(layer(lakes_kg), 1),
                            lakes_from='results/drainage.json (runoff from climate run A)'),
        crust=rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    data = results()
    data['producer'] = dict(domain='geography', files={'groundwater.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'groundwater.json').write_text(json.dumps(data, indent=1) + '\n')
    p = data['porosity_profile']
    print(f"porosity {p['phi0']:.3f} at the surface, e-folding {p['e_folding_m']/1000:.1f} km; seas "
          f"{data['standing_water']['seas_layer_m']} m, lakes {data['standing_water']['lakes_layer_m']} m (global layers)")
    for r in data['crust']:
        print(f"  saturated to {r['saturation_depth_m']/1000:4.0f} km: crust {r['crust_layer_m']:6.0f} m; total {r['total_kg']:.2e} kg "
              f"= {r['total_over_seas']} x the seas")
    return 0


if __name__ == '__main__':
    sys.exit(main())
