"""Compare saved angular cloud radiance with an independent ground-flux run.

Angular quadrature uses the raw photon groups, before display smoothing.
Errors describe Monte Carlo sampling, not angular discretization or model bias.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summary(values, unit):
    return {f'mean_{unit}': float(values.mean()),
            f'standard_error_{unit}': float(values.std(ddof=1)/np.sqrt(values.size))}


def audit(cloud, ground):
    refinement = json.loads(ground.read_text())
    if refinement['cloud_map_sha256'] != digest(cloud):
        raise ValueError('Ground estimate belongs to another cloud map')
    with np.load(cloud) as z:
        elevation = np.radians(z['elevation_deg'])
        azimuth = np.radians(z['azimuth_deg'])
        cloudy, clear = z['cloudy_groups'][..., 1], z['clear_groups'][..., 1]
        meta = json.loads(str(z['metadata']))
    if refinement['input_scene_sha256'] != meta['input_scene_sha256']:
        raise ValueError('Ground estimate and angular map describe different scenes')
    da = (np.r_[azimuth[1:], azimuth[0]+2*np.pi]
          - np.r_[azimuth[-1]-2*np.pi, azimuth[:-1]])/2
    mu = np.r_[0, np.sin(elevation), 1]

    def flux(values):
        padded = np.concatenate([values[:1], values, values[-1:]], axis=0)
        return np.trapezoid((padded*da[None, :, None]).sum(axis=1)*mu[:, None], mu, axis=0)

    cf, qf = flux(cloudy), flux(clear)
    source_names = meta['source_map_scope']
    ground_mean = sum(refinement['components'][name]['diffuse_horizontal_lux'] for name in source_names)
    ground_se = np.sqrt(sum(refinement['components'][name]['diffuse_standard_error_lux']**2 for name in source_names))
    angular = summary(cf, 'lux')
    combined_se = np.hypot(angular['standard_error_lux'], ground_se)
    record = dict(schema='terluna.research.coastal-transport-audit/1',
        cloud_map_sha256=digest(cloud), ground_product_sha256=digest(ground),
        source_map_scope=source_names, cloudy_angular_flux=angular,
        clear_angular_flux=summary(qf, 'lux'), paired_angular_flux_change=summary(cf-qf, 'lux'),
        independent_ground_diffuse=dict(mean_lux=ground_mean, standard_error_lux=float(ground_se)),
        difference_in_combined_standard_errors=float((cf.mean()-ground_mean)/combined_se),
        patches={})
    for label, er, ar in [('cloud_bank', (15, 35), (260, 280)),
                          ('upper_sky', (55, 80), (240, 300))]:
        emask = (elevation >= np.radians(er[0])) & (elevation <= np.radians(er[1]))
        amask = (azimuth >= np.radians(ar[0])) & (azimuth <= np.radians(ar[1]))
        c = cloudy[emask][:, amask].mean(axis=(0, 1))
        q = clear[emask][:, amask].mean(axis=(0, 1))
        record['patches'][label] = dict(elevation_range_deg=er, azimuth_range_deg=ar,
            cloudy_mean_cd_m2=float(c.mean()), clear_mean_cd_m2=float(q.mean()),
            paired_change=summary(c-q, 'cd_m2'))
    record.update(producer_sha256=digest(__file__),
        reading_rule='Flux integrates raw angular samples with periodic azimuth weights and endpoint continuation; patch means weight samples equally. Ground flux uses an independent seed. Excludes direct beams and retained formal sources from this comparison.',
        evidence='Agreement of two numerical estimators within the same optical scenario; not empirical validation.',
        limits=['Sampling errors omit angular discretization, unresolved cross-ring morphology and particle-optics uncertainty.',
                'A matching scene/cloud hash is required; no substitute data are admitted.'])
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cloud', type=Path, required=True)
    p.add_argument('--ground', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    record = audit(a.cloud, a.ground)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k: record[k] for k in ['cloudy_angular_flux', 'independent_ground_diffuse', 'difference_in_combined_standard_errors']}, indent=2))


if __name__ == '__main__':
    main()
