"""Independent angular integrals, occlusion limits and scene energy accounting."""
import json
import math

import numpy as np
import pytest

from .angular import AngularSky, LightField, PackedAtlas, direction, sky_quadrature, uniform_sky
from .directional import field
from .run import HERE, ROOT, digest, eye_illuminance, read_input
from .scenes import Rectangle, Scene, cosine_directions, intersect, probe, scene_set, with_step


def test_exact_sky_normalization_for_linear_sine_and_azimuth_field():
    mu = np.array([0., .07, .3, .8, 1.])
    az = np.array([0., .3, 1.2, math.pi])
    sky = AngularSky(mu, az, (1 + 2 * mu[:, None]) * (1 + .2 * az[None, :]))
    expected = 2 * (math.pi + .1 * math.pi ** 2) * (1 / 2 + 2 / 3)
    assert sky.horizontal == pytest.approx(expected, rel=1e-14)
    rays, weights = sky_quadrature(128)
    assert np.sum(sky(rays) * rays[:, 2] * weights) == pytest.approx(1, rel=2e-6)


def test_angular_sky_symmetry_and_domain():
    sky, _ = PackedAtlas().get('moon', 30)
    assert sky(direction(10, 45)[None])[0] == pytest.approx(sky(direction(10, -45)[None])[0])
    with pytest.raises(ValueError):
        sky(direction(-10, 0)[None])
    with pytest.raises(ValueError):
        PackedAtlas().get('moon', 75)
    with pytest.raises(ValueError):
        AngularSky([0, 1], [0, math.pi], np.zeros((2, 2)))


def test_uniform_sky_eye_integral_matches_analytic_tilted_hemispheres():
    light = LightField(uniform_sky(), 30, 60, 40)
    for e, a in [(-60, 15), (-30, 180), (0, 90), (30, 45), (90, 0)]:
        calculated = light.open_plane(direction(e, a), .3)
        expected = eye_illuminance(60, 40, .3, 30, e, a)
        assert calculated['total_lux'] == pytest.approx(expected['total_lux'], rel=8e-5)


def test_packed_sky_reconstructs_stored_horizontal_flux():
    atlas = PackedAtlas()
    for world in ('moon', 'earth'):
        for h in (90, 60, 30, 20, 10, 5, 1):
            _, audit = atlas.get(world, h)
            assert abs(audit['reconstructed_over_stored'] - 1) < .004


def test_shape_transfer_keeps_current_horizontal_flux_and_beam():
    atlas, data = PackedAtlas(), read_input()
    for world in ('moon', 'earth'):
        light, _ = field(data, atlas, world, 90)
        upward = light.open_plane(direction(90, 0))
        assert upward['sky_lux'] == pytest.approx(light.diffuse_h, rel=3e-5)
        assert upward['sun_direct_lux'] == light.direct_h


def test_ground_and_two_sided_rectangle_intersections():
    scene = Scene('test', patch_reflectance=.8,
                  rectangles=(Rectangle('roof', 2, 3, (-3, -3), (3, 3), .3),))
    origins = np.array([[0, 0, 1], [0, 0, 1], [5, 0, 1], [0, 0, 4.]])
    rays = np.array([[0, 0, 1], [0, 0, -1], [0, 0, -1], [0, 0, -1.]])
    t, n, rho, ids = intersect(scene, origins, rays)
    assert t == pytest.approx([2, 1, 1, 1])
    assert rho == pytest.approx([.3, .8, .1, .3])
    assert ids.tolist() == [1, 0, 0, 1]
    assert n[:, 2] == pytest.approx([-1, 1, 1, 1])


def test_cosine_sampling_rotates_to_the_requested_hemisphere():
    # These samples give a cosine distribution whose mean cosine is 2/3.
    u = (np.arange(10000) + .5) / 10000
    samples = np.stack([u, np.mod(u * 199, 1)], axis=-1)
    n = direction(-35, 80)
    rays = cosine_directions(samples, n)
    assert np.linalg.norm(rays, axis=1) == pytest.approx(np.ones(len(u)), abs=1e-14)
    assert np.min(rays @ n) >= 0
    assert np.mean(rays @ n) == pytest.approx(2 / 3, rel=1e-6)


def test_open_scene_transport_recovers_uniform_sky_and_lambertian_ground():
    light = LightField(uniform_sky(), 30, 60, 40)
    for e, a in [(0, 90), (-30, 180), (90, 0)]:
        r = probe(Scene('open'), light, [0, 0, 1.6], direction(e, a), power=15)
        expected = eye_illuminance(60, 40, .1, 30, e, a)['total_lux']
        assert r['total_lux'] == pytest.approx(expected, rel=.0015)
        assert r['sampled_tail_upper_lux'] == 0


def test_black_canopy_matches_independent_rectangle_view_factor():
    scene = Scene('black', ground_reflectance=0, patch_reflectance=0,
                  rectangles=(Rectangle('roof', 2, 3, (-3, -3), (3, 3), 0),))
    light = LightField(uniform_sky(), 90, 60, 40)
    r = probe(scene, light, [0, 0, 0], [0, 0, 1], power=16)
    # Point-to-square projected solid angle / pi, half width = height = 3 m.
    factor = 4 / math.pi / math.sqrt(2) * math.atan(1 / math.sqrt(2))
    assert r['total_lux'] == pytest.approx(40 * (1 - factor), rel=.001)
    assert r['sun_direct_lux'] == 0
    assert r['sun_reflected_lux'] == 0


def test_pale_roof_needs_illuminated_ground_to_reflect_from_its_underside():
    # With a perfectly black infinite ground, an opaque horizontal roof's
    # lower face has no source. Its reflectance cannot change the ground lux.
    light = LightField(uniform_sky(), 90, 60, 40)
    results = []
    for rho in (0, .8):
        scene = Scene('roof', ground_reflectance=0, patch_reflectance=0,
                      rectangles=(Rectangle('roof', 2, 3, (-3, -3), (3, 3), rho),))
        results.append(probe(scene, light, [0, 0, 0], [0, 0, 1], power=13)['total_lux'])
    assert results[0] == results[1]


def test_multiple_reflections_are_positive_and_component_ledgers_close():
    scene = scene_set()[-1]
    light = LightField(uniform_sky(), 90, 60, 40)
    shallow = probe(scene, light, [0, 0, 1.6], [0, -1, 0], power=14, max_bounces=1)
    deep = probe(scene, light, [0, 0, 1.6], [0, -1, 0], power=14, max_bounces=12)
    assert deep['total_lux'] > shallow['total_lux']
    assert deep['ambient_lux'] == pytest.approx(deep['sky_unobstructed_lux'] + deep['sky_reflected_lux']
                                               + deep['sun_reflected_lux'])
    assert deep['ambient_lux'] == pytest.approx(deep['sky_unobstructed_lux'] + deep['ground_in_view_lux']
                                               + deep['structures_in_view_lux'])
    assert deep['sampled_tail_upper_lux'] < shallow['sampled_tail_upper_lux']


def test_scene_transport_is_linear_in_illumination():
    scene = scene_set()[2]
    first = probe(scene, LightField(uniform_sky(), 30, 60, 40), [0, 0, 1.6], [0, -1, 0], power=12)
    second = probe(scene, LightField(uniform_sky(), 30, 30, 20), [0, 0, 1.6], [0, -1, 0], power=12)
    for key in first:
        assert second[key] == pytest.approx(first[key] / 2, abs=1e-12)


def test_step_casts_a_shadow_and_does_not_self_occlude_outgoing_probes():
    scene = with_step(Scene('step'))
    ray = direction(30, 90)
    # A low point behind the block looks into it; a point above the tread clears it.
    t = intersect(scene, np.array([[0, -2.1, .01], [0, -1.5, .150001]]),
                  np.broadcast_to(ray, (2, 3)))[0]
    assert np.isfinite(t[0])
    assert not np.isfinite(t[1])


def test_scene_rejects_nonconservative_reflectance():
    with pytest.raises(ValueError):
        Scene('bad', patch_reflectance=1)


def test_stored_directional_product_ledgers_and_hashes():
    product = json.loads((HERE / 'results/directional_scenes.json').read_text())
    for section in ('files', 'input_interpretation_files'):
        for path, expected in product['producer'][section].items():
            assert digest(ROOT / path) == expected
    for path, item in product['producer']['inputs'].items():
        assert digest(ROOT / path) == item['sha256']
    assert len(product['angular_views']) == 630
    assert len(product['scenes']) == 24
    for row in product['scenes']:
        for p in row['probes'].values():
            assert p['total_lux'] == pytest.approx(p['ambient_lux'] + p['sun_direct_lux'])
            assert p['ambient_lux'] == pytest.approx(p['sky_unobstructed_lux'] + p['sky_reflected_lux']
                                                    + p['sun_reflected_lux'])
        tasks = row['tasks']
        a, b = tasks['step_tread_luminance_cd_m2'], tasks['step_riser_luminance_cd_m2']
        assert tasks['step_signed_michelson'] == pytest.approx((a - b) / (a + b))
        assert tasks['work_mark_luminance_cd_m2'] / tasks['work_background_luminance_cd_m2'] - 1 == -.5


@pytest.mark.parametrize('world', ['moon', 'earth'])
def test_stored_enclosed_scene_reproduces(world):
    product = json.loads((HERE / 'results/directional_scenes.json').read_text())
    stored = next(r for r in product['scenes'] if r['scene'] == 'screened_canopy_pale'
                  and r['world'] == world and r['sun_deg'] == 30)['probes']['step_riser']
    settings = product['numerical_settings']
    light, _ = field(read_input(), PackedAtlas(), world, 30)
    scene = with_step(scene_set()[-1])
    estimates = [probe(scene, light, [0, -1, .075], [0, 1, 0],
                       power=int(math.log2(settings['paths_per_scramble'])), seed=seed,
                       max_bounces=settings['maximum_surface_bounces'])
                 for seed in settings['scramble_seeds']]
    for key in estimates[0]:
        assert stored[key] == pytest.approx(np.mean([e[key] for e in estimates]), rel=1e-12, abs=1e-9)
