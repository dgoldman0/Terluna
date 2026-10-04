"""Gaze and finite matte scenes illuminated by the solved spherical sky product.

OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.spherical
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .angular import AngularSky
from .directional import SUNS, angular_rows, scene_rows
from .run import CASES, HERE, ROOT, digest

INPUT = ROOT/'illumination/sky/results/solved_sky.json'
BASELINE = HERE/'results/native_sky.json'
SCHEMA = 'terluna.research.optical-comfort-spherical-sky/1'
WORLDS = dict(moon='moon_1.2atm',earth='earth_control')


class SolvedAtlas:
    def __init__(self,product):
        if product['schema'] != 'terluna.illumination.solved-spherical-sky/1':
            raise ValueError('Unsupported solved sky product')
        self.records = {}
        for world,name in WORLDS.items():
            record = product['worlds'][name]['atlas']
            path = ROOT/record['path']
            if digest(path) != record['sha256']:
                raise ValueError(f'Solved sky archive hash differs: {path}')
            with np.load(path,allow_pickle=False) as archive:
                header = json.loads(str(archive['meta']))
                if header['schema'] != 'terluna.illumination.solved-spherical-sky-atlas/1' or header['world'] != name:
                    raise ValueError('Unsupported solved sky archive')
                self.records[world] = {k:archive[k].copy() for k in ('suns','elevations','azimuths','xyz')}

    def get(self,world,sun_deg):
        a = self.records[world]
        index = np.flatnonzero(a['suns']==sun_deg)
        if len(index)!=1:
            raise ValueError('An exact stored solar elevation is required')
        sky = AngularSky(np.sin(np.radians(a['elevations'])),np.radians(a['azimuths']),a['xyz'][index[0],:,:,1])
        return sky,dict(atlas_world=WORLDS[world],sun_deg=sun_deg,reconstructed_diffuse_lux=sky.horizontal,
                        interpretation='Angular shape and horizontal flux from the same solved spherical atmosphere')


def light_table(product):
    """Adapt the product's flux names to the existing scene integrator interface."""
    cases = {}
    for world,name in WORLDS.items():
        summaries = {str(r['sun_deg']):{'0.1':dict(direct_illuminance_lux=r['direct_horizontal_lux'],
                                                illuminance_lux=r['total_horizontal_lux'])}
                     for r in product['worlds'][name]['samples'] if r['sun_deg']>0}
        cases[CASES[world]] = dict(summaries=summaries)
    return dict(cases=cases)


def build():
    product = json.loads(INPUT.read_text())
    baseline = json.loads(BASELINE.read_text())
    atlas = SolvedAtlas(product)
    settings = baseline['numerical_settings']
    data = light_table(product)
    views,audits,error = angular_rows(data,atlas)
    for row in views:
        row['flux_quality'] = 'solved_column_spherical_multiple_scattering'
    scenes = scene_rows(data,atlas,int(np.log2(settings['paths_per_scramble'])),settings['maximum_surface_bounces'])
    comparisons = []
    for old,new in zip(baseline['angular_views'],views,strict=True):
        keys = ('world','sun_deg','gaze_deg','azimuth_from_sun_deg')
        if any(old[k]!=new[k] for k in keys):
            raise ValueError('Gaze comparison order differs')
        comparisons.append(dict(**{k:new[k] for k in keys},ambient_relative_change=new['ambient_lux']/old['ambient_lux']-1,
                                 total_relative_change=new['total_lux']/old['total_lux']-1))
    scene_changes = []
    for old,new in zip(baseline['scenes'],scenes,strict=True):
        keys = ('world','sun_deg','scene')
        if any(old[k]!=new[k] for k in keys):
            raise ValueError('Scene comparison order differs')
        scene_changes.append(dict(**{k:new[k] for k in keys},
                                  probe_total_relative_change={k:v['total_lux']/old['probes'][k]['total_lux']-1
                                                               for k,v in new['probes'].items()},
                                  step_contrast_change=new['tasks']['step_signed_michelson']-old['tasks']['step_signed_michelson']))
    files = [HERE/f for f in ('spherical.py','directional.py','angular.py','scenes.py','run.py')]
    probes = [p for s in scenes for p in s['probes'].values()]
    beam_projection = [dict(world=world,sun_deg=row['sun_deg'],relative_normal_beam_difference=
                            row['direct_horizontal_lux']/np.sin(np.radians(row['sun_deg']))/row['direct_normal_lux']-1)
                       for world,name in WORLDS.items() for row in product['worlds'][name]['samples'] if row['sun_deg'] in SUNS]
    return dict(schema=SCHEMA,producer=dict(lane='research',runner='research/studies/optical_comfort/spherical.py',
                                          files={str(p.relative_to(ROOT)):digest(p) for p in files},
                                          inputs={str(p.relative_to(ROOT)):digest(p) for p in (INPUT,BASELINE)}),
                evidence='Conditional gaze and finite matte-scene calculations with one consistent spherical sky '
                         'for each solved atmospheric column. The lunar source includes the chosen shield and 5% dimming.',
                reading_rule='Absolute directional radiance and flux originate in solved_sky.json. The existing '
                             'scene quadrature normalizes each angular interpolant to that solution’s horizontal '
                             'diffuse illuminance, obtained from the same interpolant. Scenes place the direct beam '
                             'at the solar centre and infer normal illuminance from horizontal flux; the projection '
                             'difference from the finite-disk beam is recorded. Geometry, reflectance and Sobol '
                             'scrambles match native_sky.json.',
                units=dict(illuminance='lux',luminance='cd m^-2',angles='degrees',changes='dimensionless'),
                numerical_settings=settings,atlas_audit=audits,angular_views=views,scenes=scenes,
                direct_beam_projection_audit=beam_projection,
                changes_from_hybrid=dict(angular=comparisons,scenes=scene_changes),
                checks=dict(maximum_angular_coefficient_change_64_to_128=error,
                            maximum_direct_beam_projection_relative_difference=max(abs(r['relative_normal_beam_difference']) for r in beam_projection),
                            maximum_scene_nested_sample_relative_change=max(p['nested_sample_relative_change'] for p in probes),
                            maximum_scene_scramble_relative_difference=max(p['scramble_relative_difference'] for p in probes),
                            maximum_sampled_tail_upper_fraction=max(p['tail_upper_fraction'] for p in probes)))


def main():
    product = build()
    path = HERE/'results/spherical_sky.json'
    path.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps(product['checks'],indent=2))


if __name__=='__main__':
    main()
