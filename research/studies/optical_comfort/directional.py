"""Directional sky, physical shade and visual-task extension of the optical study.

python -m research.studies.optical_comfort.directional
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import scipy

from .angular import ATLAS, PackedAtlas, LightField, direction, sky_quadrature
from .run import CASES, HERE, INPUT, ROOT, digest, eye_illuminance, read_input
from .scenes import PROBES, SOBOL_BOUNCES, probe, scene_set, with_step

SCHEMA = 'terluna.research.optical-comfort-directional-scenes/1'
SUNS = (90., 60., 30., 20., 10., 5., 1.)
GAZES = (-90., -75., -60., -45., -30., -15., 0., 15., 30.)
AZIMUTHS = (0., 45., 90., 135., 180.)
SEEDS = (17, 43)


def field(data, atlas, world, h):
    s = data['cases'][CASES[world]]['summaries'][str(float(h))]['0.1']
    sky, audit = atlas.get(world, h)
    return LightField(sky, h, s['direct_illuminance_lux'],
                      s['illuminance_lux'] - s['direct_illuminance_lux']), audit


def angular_rows(data, atlas):
    rows, audits, errors = [], [], []
    angles = [(g, a) for g in GAZES for a in AZIMUTHS]
    normals = np.array([direction(g, a) for g, a in angles])
    for world in ('moon', 'earth'):
        for h in SUNS:
            light, audit = field(data, atlas, world, h)
            audits.append(dict(world=world, **audit))
            estimates = []
            for order in (64, 128):
                rays, weights = sky_quadrature(order)
                coefficients = (light.sky(rays) * weights) @ np.maximum(rays @ normals.T, 0)
                estimates.append(coefficients)
            errors.extend(abs(estimates[0] - estimates[1]).tolist())
            for (g, a), n, coefficient in zip(angles, normals, estimates[1]):
                sky = float(coefficient * light.diffuse_h)
                ground = .1 * (light.direct_h + light.diffuse_h) * (1 - n[2]) / 2
                direct = light.dni * max(float(n @ light.sun), 0)
                uniform = eye_illuminance(light.direct_h, light.diffuse_h, .1, h, g, a)
                rows.append(dict(world=world, sun_deg=h, landscape_albedo=.1,
                                 gaze_deg=g, azimuth_from_sun_deg=a,
                                 sky_lux=sky, ground_lux=float(ground), sun_direct_lux=direct,
                                 ambient_lux=sky + ground, total_lux=sky + ground + direct,
                                 uniform_ambient_lux=uniform['sky_lux'] + uniform['ground_lux'],
                                 uniform_total_lux=uniform['total_lux'],
                                 sky_eye_over_horizontal_diffuse=float(coefficient),
                                 flux_quality='two_stream_low_sun_caution' if h <= 20 else 'two_stream_screen'))
    return rows, audits, max(errors)


def scene_rows(data, atlas, power=17, max_bounces=18):
    rows = []
    for scene in scene_set():
        stepped = with_step(scene)
        for h in (90., 30.):
            for world in ('moon', 'earth'):
                light, _ = field(data, atlas, world, h)
                probes = {}
                for name, (position, normal) in PROBES.items():
                    geometry = stepped if name.startswith('step_') else scene
                    estimates = [probe(geometry, light, position, normal, power, seed, max_bounces)
                                 for seed in SEEDS]
                    values = {k: float(np.mean([e[k] for e in estimates])) for k in estimates[0]}
                    values.update(
                        nested_sample_relative_change=abs(values['coarse_total_lux'] / values['total_lux'] - 1),
                        scramble_relative_difference=abs(estimates[0]['total_lux'] - estimates[1]['total_lux'])
                                                     / values['total_lux'],
                        tail_upper_fraction=values['sampled_tail_upper_lux'] / values['total_lux'])
                    probes[name] = values
                tread = .2 * probes['step_tread']['total_lux'] / math.pi
                riser = .2 * probes['step_riser']['total_lux'] / math.pi
                work = probes['workplane']['total_lux']
                rows.append(dict(world=world, sun_deg=h, scene=scene.name, probes=probes,
                                 tasks=dict(step_tread_luminance_cd_m2=tread,
                                            step_riser_luminance_cd_m2=riser,
                                            step_signed_michelson=(tread - riser) / (tread + riser),
                                            work_background_luminance_cd_m2=.2 * work / math.pi,
                                            work_mark_luminance_cd_m2=.1 * work / math.pi,
                                            work_mark_weber=-.5)))
            print(f'Integrated {scene.name}, Sun {h:g} degrees', flush=True)
    return rows


def build(power=17, max_bounces=18):
    data, atlas = read_input(), PackedAtlas()
    angles, audits, quadrature_error = angular_rows(data, atlas)
    scenes = scene_rows(data, atlas, power, max_bounces)
    producer_files = [HERE / f for f in ('directional.py', 'angular.py', 'scenes.py', 'run.py', 'sources.json')]
    # These files document the packed payload's axes, source separation and Y preservation.
    interpretation_files = [ROOT / p for p in ('immersion/bake/sky_atlas.py',
                            'illumination/sky/render_data.py', 'illumination/sky/colour_matching.py',
                            'illumination/sky/METHODS.md')]
    all_probes = [p for r in scenes for p in r['probes'].values()]
    return dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/optical_comfort/directional.py',
                      files={str(p.relative_to(ROOT)): digest(p) for p in producer_files},
                      input_interpretation_files={str(p.relative_to(ROOT)): digest(p) for p in interpretation_files},
                      inputs={str(INPUT.relative_to(ROOT)): dict(schema=data['schema'], sha256=digest(INPUT)),
                              str(ATLAS.relative_to(ROOT)): dict(format='packed-atlas version 1', sha256=digest(ATLAS))},
                      runtime=dict(numpy=np.__version__, scipy=scipy.__version__)),
        evidence='Angular-shape transfer and prescribed finite-scene transport. The packed atlas computes '
                 'an unfiltered exponential proxy atmosphere over albedo 0.1. Its normalized sky shape is '
                 'combined with current solved-column direct and diffuse photopic fluxes, separately for each world. '
                 'This is a hybrid sensitivity, not a shielded spherical atmosphere solution. '
                 'Neutral Lambertian local scenes include occlusion and repeated diffuse reflection. '
                 'Task luminance and contrast are optical signals, without a psychophysical visibility threshold.',
        reading_rule='Positive z is up; Sun is at azimuth 0 along +x; azimuth 90 is +y. '
                     'Packed upper-sky RGB is converted with Y=(0.2126,0.7152,0.0722), then bilinearly '
                     'interpolated in sin(view elevation) and symmetric azimuth. The exact horizontal '
                     'integral of that interpolant normalizes each frame before applying current diffuse lux. '
                     'The Sun is a separate collimated beam. No interpolation in solar elevation. '
                     'All boundary fields use landscape albedo 0.1; finite pale patches and walls never '
                     'change this atmospheric boundary. Ground reflection is recalculated by scene transport, '
                     'rather than taken from the atlas lower hemisphere. Scene rows are means of two fixed '
                     'scrambled Sobol sequences. ambient excludes the beam reaching the eye directly but includes '
                     'reflected sunlight. The step is inserted only for its two surface probes; other probes '
                     'describe the same setting without that object.',
        units=dict(length='metres', illuminance='lux', luminance='cd m^-2', angles='degrees',
                   reflectance='dimensionless, photopic and spectrally neutral'),
        assumptions=dict(landscape_albedo=.1, sun_elevations_scene_deg=[90., 30.],
                         local_atmospheric_transport='No extinction or in-scattering between local surfaces; fixed sky boundary.',
                         sky_limit='Packed resolution is 41 x 33 over the full elevation span and symmetric half azimuth; '
                                   '21 elevation samples cover the upper sky. Quantization and interpolation checks '
                                   'do not validate the source model or bound unresolved angular structure.',
                         upstream_limit='Clear aerosol-free sky, old exponential optical profiles, unfiltered sunlight, '
                                        'open upstream whole-sphere energy residual up to 4.6%. Current two-stream '
                                        'diffuse flux is deficient at low Sun; samples at <=20 degrees are flagged.',
                         materials='Opaque, nonemissive Lambertian rectangles, same reflectance on both faces; '
                                   'no gloss, transmission, polarization, textures, rough water or canopy model.',
                         task='Neutral mark reflectance 0.1 on a 0.2 workplane; infinitesimal measurement patch at '
                              '(0,-1,0.8). A separate 2 x 1 x 0.15 m raised block has reflectance 0.2. '
                              'Its tread and viewer-facing riser are sampled at their centres. '
                              'Task contrast excludes ocular scatter, adaptation and angular-resolution limits.'),
        scene_geometry=[s.record() for s in scene_set()],
        step_geometry=with_step(scene_set()[0]).record()['rectangles'],
        probe_geometry={name: dict(position_m=p, normal=n) for name, (p, n) in PROBES.items()},
        numerical_settings=dict(paths_per_scramble=2 ** power, nested_coarse_paths=2 ** (power - 2),
                                scramble_seeds=SEEDS, maximum_surface_bounces=max_bounces,
                                sobol_dimensions=2 * (SOBOL_BOUNCES + 1), angular_quadrature_orders=[64, 128]),
        atlas_audit=audits, angular_views=angles, scenes=scenes,
        checks=dict(maximum_angular_coefficient_change_64_to_128=quadrature_error,
                    maximum_packed_horizontal_relative_difference=max(abs(a['reconstructed_over_stored'] - 1) for a in audits),
                    maximum_scene_nested_sample_relative_change=max(p['nested_sample_relative_change'] for p in all_probes),
                    maximum_scene_scramble_relative_difference=max(p['scramble_relative_difference'] for p in all_probes),
                    maximum_sampled_tail_upper_fraction=max(p['tail_upper_fraction'] for p in all_probes),
                    numerical_interpretation='Nested refinement and two-scramble differences are diagnostics, '
                                             'not confidence intervals. Tail bounds are conservative per surviving '
                                             'path and then sampled; they are not a certified global quadrature bound.'),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--power', type=int, default=17)
    parser.add_argument('--bounces', type=int, default=18)
    parser.add_argument('--output', type=Path, default=HERE / 'results/directional_scenes.json')
    args = parser.parse_args()
    product = build(args.power, args.bounces)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, allow_nan=False) + '\n')
    print(json.dumps(product['checks'], indent=2))
    print(args.output)


if __name__ == '__main__':
    main()
