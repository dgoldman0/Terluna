"""Repeat directional and finite-scene screens using the local native sky atlases.

OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.local_sky
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .angular import AngularSky, PackedAtlas, WORLD_ATLAS, Y_WEIGHTS
from .directional import angular_rows, scene_rows
from .run import HERE, INPUT, ROOT, digest, read_input

SCHEMA = 'terluna.research.optical-comfort-native-sky/1'
BASELINE = HERE / 'results/directional_scenes.json'
MANIFEST = HERE / 'local_inputs.json'


class NativeAtlas:
    def __init__(self, directory=None):
        directory = Path(directory or ROOT / 'illumination/sky/data')
        expected = json.loads(MANIFEST.read_text())['native_atlases']
        self.maps, self.records, self.inputs = {}, {}, {}
        for name in WORLD_ATLAS.values():
            path = directory / (name + '_atlas.npz')
            actual = digest(path)
            if actual != expected[path.name]['sha256']:
                raise ValueError(f'Native atlas hash differs from the admitted local input: {path}')
            self.inputs[name] = dict(file=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name,
                                     sha256=actual)
            with np.load(path, allow_pickle=False) as a:
                self.records[name] = {key: a[key].copy() for key in
                                      ('suns', 'elevs', 'azs', 'direct', 'direct_horizontal', 'diffuse')}
                self.records[name]['metadata'] = json.loads(str(a['metadata']))
                rgb = a['rgb'].astype(float)
            r = self.records[name]
            if (rgb.shape != (len(r['suns']), 81, 65, 3) or
                    not np.isfinite(rgb).all() or np.any(np.diff(r['suns']) <= 0) or
                    not np.allclose(r['elevs'], np.sign(np.linspace(-1, 1, 81)) *
                                    np.linspace(-1, 1, 81) ** 2 * 90) or
                    not np.allclose(r['azs'], np.linspace(0, 180, 65))):
                raise ValueError('Unsupported native atlas grid')
            self.maps[name] = rgb @ Y_WEIGHTS
            if np.any(self.maps[name] < 0):
                raise ValueError('Negative native photopic radiance')

    def get(self, world, sun_deg):
        name = WORLD_ATLAS[world]
        r = self.records[name]
        indices = np.flatnonzero(r['suns'] == sun_deg)
        if len(indices) != 1:
            raise ValueError('An exact stored solar elevation is required')
        i = indices[0]
        sky = AngularSky(np.sin(np.radians(r['elevs'][40:])), np.radians(r['azs']), self.maps[name][i, 40:])
        stored = float(r['diffuse'][i] @ Y_WEIGHTS)
        return sky, dict(atlas_world=name, sun_deg=sun_deg, stored_diffuse_lux=stored,
                         reconstructed_diffuse_lux=sky.horizontal,
                         reconstructed_over_stored=sky.horizontal / stored,
                         source_archive_sha256=self.inputs[name]['sha256'])


def build():
    native, packed, data = NativeAtlas(), PackedAtlas(), read_input()
    baseline = json.loads(BASELINE.read_text())
    # Fixed seeds and geometry make this a paired change of angular resolution.
    settings = baseline['numerical_settings']
    power = int(np.log2(settings['paths_per_scramble']))
    angles, audits, error = angular_rows(data, native)
    scenes = scene_rows(data, native, power, settings['maximum_surface_bounces'])
    recovery = []
    for world, name in WORLD_ATLAS.items():
        full, r, w = native.maps[name], native.records[name], packed.data['worlds'][name]
        if not np.array_equal(r['suns'], w['suns']):
            raise ValueError('Native and packed solar grids differ')
        old_hash = next(p['sha256'] for p in packed.data['provenance'] if p['file'] == name + '_atlas.npz')
        shared = full[:, ::2, ::2]
        difference = abs(shared - packed.maps[name])
        lit = shared > 1e-6
        recovery.append(dict(world=world, native_sha256=native.inputs[name]['sha256'],
                             packed_recorded_archive_sha256=old_hash,
                             archive_bytes_identical=native.inputs[name]['sha256'] == old_hash,
                             beam_and_flux_tables_exactly_equal=all(np.array_equal(r[k], w[k]) for k in
                                                                    ('direct', 'direct_horizontal', 'diffuse')),
                             maximum_shared_sample_relative_difference=float(np.max(difference[lit] / shared[lit])),
                             shared_sample_absolute_difference_over_sum=float(difference.sum() / shared.sum()),
                             atmosphere=r['metadata']['atmosphere'], quality=r['metadata']['quality']))
    view_changes = []
    for old, new in zip(baseline['angular_views'], angles, strict=True):
        keys = ('world', 'sun_deg', 'gaze_deg', 'azimuth_from_sun_deg')
        if any(old[k] != new[k] for k in keys):
            raise ValueError('Baseline angular rows differ')
        view_changes.append(dict(**{k: new[k] for k in keys},
                                 ambient_relative_change=new['ambient_lux'] / old['ambient_lux'] - 1,
                                 total_relative_change=new['total_lux'] / old['total_lux'] - 1))
    scene_changes = []
    for old, new in zip(baseline['scenes'], scenes, strict=True):
        keys = ('world', 'sun_deg', 'scene')
        if any(old[k] != new[k] for k in keys):
            raise ValueError('Baseline scene rows differ')
        scene_changes.append(dict(**{k: new[k] for k in keys},
                                  probe_total_relative_change={k: p['total_lux'] / old['probes'][k]['total_lux'] - 1
                                                               for k, p in new['probes'].items()},
                                  step_contrast_change=new['tasks']['step_signed_michelson'] -
                                                       old['tasks']['step_signed_michelson']))
    files = [Path(__file__), MANIFEST, *[HERE / f for f in ('directional.py', 'angular.py', 'scenes.py', 'run.py')]]
    probes = [p for s in scenes for p in s['probes'].values()]
    return dict(schema=SCHEMA, producer=dict(lane='research', runner=str(Path(__file__).relative_to(ROOT)),
                files={str(p.relative_to(ROOT)): digest(p) for p in files},
                inputs={str(p.relative_to(ROOT)): digest(p) for p in (BASELINE, INPUT)},
                native_inputs=native.inputs),
                evidence='Local full-resolution angular-shape transfer and paired finite matte scenes. '
                         'The original exponential unfiltered atmospheric profiles remain in use for sky direction; '
                         'current solved-column shielded photopic fluxes set its amplitude. '
                         'This isolates archive packing and angular resolution within the existing hybrid model.',
                reading_rule='Same scene geometry, fluxes, gaze grid, interpolation rule and fixed Sobol seeds as '
                             'directional_scenes.json. Native upper sky has 41 elevation by 65 azimuth samples. '
                             'Signed RGB is converted directly to photopic Y. Local and packed archive hashes '
                             'differ; array comparisons below establish their numerical correspondence.',
                units=dict(illuminance='lux', luminance='cd m^-2', angles='degrees', changes='dimensionless'),
                numerical_settings=settings, recovery_audit=recovery, atlas_audit=audits,
                angular_views=angles, scenes=scenes, angular_changes=view_changes, scene_changes=scene_changes,
                checks=dict(maximum_angular_coefficient_change_64_to_128=error,
                            maximum_ambient_relative_change=max(abs(r['ambient_relative_change']) for r in view_changes),
                            maximum_scene_total_relative_change=max(abs(v) for r in scene_changes
                                                                    for v in r['probe_total_relative_change'].values()),
                            maximum_step_contrast_absolute_change=max(abs(r['step_contrast_change']) for r in scene_changes),
                            maximum_scene_nested_sample_relative_change=max(p['nested_sample_relative_change'] for p in probes),
                            maximum_scene_scramble_relative_difference=max(p['scramble_relative_difference'] for p in probes),
                            maximum_sampled_tail_upper_fraction=max(p['tail_upper_fraction'] for p in probes)))


def main():
    product = build()
    path = HERE / 'results/native_sky.json'
    path.write_text(json.dumps(product, indent=2, allow_nan=False) + '\n')
    print(json.dumps(product['checks'], indent=2))


if __name__ == '__main__':
    main()
