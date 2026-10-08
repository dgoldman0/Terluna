"""Reference three-radius aperture, finer time/Sun grids and spatial gap maps."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.fleet import (fleet_acceleration, coverage_snapshot,
    square_vertices, projected_polygons, sun_points)
from protection.dynamics.cycling import sail_basis
from protection.dynamics.optical import arriving_ray_vectors, unit
from .fleet_run import read_raw, environment, identities, digest
from .fleet_analyze import summarize_coverage
from .cycling_search import ROOT, HERE
from .cycling_analysis import compact_series_json


def spatial_map(state, normal, sample, side, radius, suns=64, width=97):
    """Pointwise finite-Sun interception map, with actual spatial coordinates.

    This supplemental map uses the same perspective square polygons as the
    area quadrature. Only uneclipsed dates are admitted. It is a visualization
    sample; exact polygon areas supply the coverage statistic.
    """
    import shapely
    q, v = state[:, :3], state[:, 3:]
    sun, _ = arriving_ray_vectors(sample['positions']['sun'], sample['positions']['earth'],
        sample['sun_v'], sample['earth_v'], sample['sun_a'], sample['earth_a'])
    u = unit(sun); _, b, c = sail_basis(sun)
    tau = np.maximum(q@u, 0)/K.SPEED_OF_LIGHT
    q = q-(v+sample['moon_v'])*tau[:, None]
    r = np.linalg.norm(q-u*(q@u)[:, None], axis=-1)
    keep = (q@u > 0) & (r < radius+side+np.maximum(q@u, 0)*K.SUN_RADIUS/np.linalg.norm(sun))
    vertices = square_vertices(q[keep], normal[keep], np.broadcast_to(sun, q[keep].shape), side-110.)
    axis = np.linspace(-radius, radius, width)
    x, y = np.meshgrid(axis, axis); mask = x*x+y*y <= radius**2
    points = shapely.points(np.column_stack([x[mask], y[mask]]))
    hits = np.zeros(len(points), dtype=int)
    for source in sun_points(sun, suns, .173):
        polys = projected_polygons(vertices, source, u, b, c)
        merged = shapely.union_all(polys)
        hits += shapely.covers(merged, points)
    field = np.full(x.shape, -1.)
    field[mask] = hits/suns
    return dict(axis_km=np.round(axis/1000, 4).tolist(),
        intercepted_solar_fraction=field.tolist(), outside_disk_value=-1,
        sampled_points=len(points), solar_directions=suns,
        points_with_zero_filtering_fraction=float(np.mean(hits == 0)),
        points_with_full_sampled_sun_filtering_fraction=float(np.mean(hits == suns)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--products', nargs='+', required=True)
    parser.add_argument('--names', nargs='*')
    parser.add_argument('--cadence-hours', type=float, default=3.)
    parser.add_argument('--suns', type=int, default=32)
    parser.add_argument('--output', default='fleet_refinement.json')
    a = parser.parse_args()
    cases, inputs = {}, {}
    for filename in a.products:
        path = HERE/'results'/filename
        p = json.loads(path.read_text())
        cases.update(p['cases']); inputs[filename] = digest(path)
    if a.names:
        cases = {n: cases[n] for n in a.names}
    results = {}
    for name, case in cases.items():
        raw, metadata = read_raw(name)
        t = raw['t']; states = raw['state']; side = metadata['config']['side_km']*1000
        env = environment(metadata['config']['days'])
        indices = np.arange(0, len(t), round(a.cadence_hours*3600/(t[1]-t[0])))
        if indices[-1] != len(t)-1:
            indices = np.r_[indices, len(t)-1]
        variants = {}
        for target in (3, 4):
            for label in ('nominal', 'collision_and_shadow_screened'):
                ids = case['selections'][label]['member_ids']
                records = []
                for i in indices:
                    state = states[i, ids]
                    d = fleet_acceleration(env, t[i], state, 4*K.MOON_RADIUS, side, diagnostics=True)
                    records.append(coverage_snapshot(state, d['normal'], env.at(t[i]), side,
                        target*K.MOON_RADIUS, sun_count=a.suns, rotation=.173))
                summary = summarize_coverage(records, t[indices])
                variant = dict(summary=summary, records=records,
                    mean_ray_coverage_change_from_coarse_4R=summary['ray_coverage']['mean']-
                        case['selections'][label]['summary']['ray_coverage']['mean'] if target == 4 else None)
                if label == 'nominal':
                    worst = int(np.argmin([r['ray_coverage'] for r in records])); i = indices[worst]
                    d = fleet_acceleration(env, t[i], states[i], 4*K.MOON_RADIUS, side, diagnostics=True)
                    convergence = {}
                    for n in (16, 64, 256):
                        convergence[str(n)] = coverage_snapshot(states[i], d['normal'], env.at(t[i]),
                            side, target*K.MOON_RADIUS, sun_count=n, rotation=.381)
                    variant['worst_date_days'] = float(t[i]/K.JULIAN_DAY)
                    variant['solar_quadrature_at_worst_date'] = convergence
                    if records[worst]['natural_eclipse_fraction'] > 1e-8:
                        raise ValueError('Spatial map requires an uneclipsed date')
                    variant['spatial_map_at_worst_date'] = spatial_map(states[i], d['normal'],
                        env.at(t[i]), side, target*K.MOON_RADIUS)
                variants[f'{target}R_{label}'] = variant
                print(json.dumps(dict(name=name, target=target, subset=label, summary=summary)), flush=True)
        results[name] = dict(variants=variants, times_days=(t[indices]/K.JULIAN_DAY).tolist())
        raw.close()
    sources = {**identities(), 'research/studies/solar_shield_array/fleet_refine.py': digest(__file__),
        'research/studies/solar_shield_array/fleet_analyze.py': digest(HERE/'fleet_analyze.py'),
        'research/studies/solar_shield_array/cycling_analysis.py': digest(HERE/'cycling_analysis.py')}
    product = dict(schema='terluna.research.solar-shield-fleet-refinement/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources), input_products=inputs),
        evidence='Actual simultaneous finite-square ray unions on the saved ephemeris histories; numerical aperture comparison and grid sensitivity',
        reference_requirement='protection/report.md section 2: established protected radius of three lunar radii; finite-Sun margin is added at the aperture',
        target_radii_m={'3R': 3*K.MOON_RADIUS, '4R': 4*K.MOON_RADIUS},
        trajectory_control_target_radii=4,
        reading_rule='The three-radius requirement is inherited. Both apertures use the same propagated fleet and four-radius service controller; no three-radius trajectory re-optimization is claimed. Large squares remain geometry probes with centre-force and structural limitations. Maps are point samples; polygon unions determine area.',
        parameters=vars(a), cases=results)
    (HERE/'results'/a.output).write_text(compact_series_json(product))


if __name__ == '__main__':
    main()
