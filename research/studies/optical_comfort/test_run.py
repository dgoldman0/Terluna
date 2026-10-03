"""Independent photometric limits and evidence boundaries for the optical screen."""
import json
import math

import numpy as np
import pytest

from research.studies.optical_comfort import run as m


def test_eye_plane_recovers_horizontal_up_and_down_limits():
    up = m.eye_illuminance(60, 40, 0.3, 30, gaze_deg=90)
    down = m.eye_illuminance(60, 40, 0.3, 30, gaze_deg=-90)
    assert up['total_lux'] == pytest.approx(100)
    assert down['total_lux'] == pytest.approx(30)
    assert down['direct_lux'] == 0
    assert m.eye_illuminance(60, 40, 0.3, 30, azimuth_from_sun_deg=180)['total_lux'] == pytest.approx(35)


def test_uniform_radiance_is_independent_of_plane_direction():
    # Uniform sky and equally luminous ground give E = pi L in every plane.
    for elevation in (-90, -30, 0, 30, 90):
        assert m.eye_illuminance(0, math.pi, 1, 30, gaze_deg=elevation)['total_lux'] == pytest.approx(math.pi)


def test_tilted_eye_against_independent_full_sphere_quadrature():
    # Integrate radiance over the sphere, clipping rays behind the eye plane.
    # Equal weights in azimuth and Gauss weights in z represent solid angle.
    z, weights = np.polynomial.legendre.leggauss(240)
    az = (np.arange(960) + 0.5) * 2 * math.pi / 960
    radial = np.sqrt(1 - z ** 2)
    rays = np.stack(np.broadcast_arrays(radial[:, None] * np.cos(az),
                                        radial[:, None] * np.sin(az), z[:, None]), axis=-1)
    luminance = np.where(z > 0, 40 / math.pi, 30 / math.pi)
    for elevation in (-60, -30, 30):
        e = math.radians(elevation)
        normal = np.array([math.cos(e), 0, math.sin(e)])
        diffuse = np.sum(np.maximum(rays @ normal, 0) * luminance[:, None] * weights[:, None]) * 2 * math.pi / 960
        sun = np.array([math.cos(math.radians(30)), 0, math.sin(math.radians(30))])
        direct = 120 * max(0, float(sun @ normal))
        calculated = m.eye_illuminance(60, 40, 0.3, 30, elevation, 0)
        assert calculated['total_lux'] == pytest.approx(diffuse + direct, rel=2e-5)


def test_eye_rejects_negative_or_nonfinite_light():
    for value in (-1, float('nan'), float('inf')):
        with pytest.raises(ValueError):
            m.eye_illuminance(value, 40, 0.3, 30)


@pytest.mark.parametrize('shape', list(m.SKY_SHAPES.values()))
def test_sky_factors_against_solid_angle_quadrature(shape):
    # Independent elevation/azimuth integral over the forward upper quarter-sphere.
    x, w = np.polynomial.legendre.leggauss(100)
    e = (x + 1) * math.pi / 4
    az = x * math.pi / 2
    a, b = shape
    radiance = a + b * np.sin(e)
    vertical = np.sum(w * radiance * np.cos(e) ** 2) * math.pi / 4 * np.sum(w * np.cos(az)) * math.pi / 2
    horizontal = 2 * math.pi * np.sum(w * radiance * np.sin(e) * np.cos(e)) * math.pi / 4
    assert vertical / horizontal == pytest.approx(m.sky_vertical_factor(a, b), rel=1e-12)


def test_masks_endpoints_and_horizon_weight():
    assert m.sky_fraction_below(0) == 0
    assert m.sky_fraction_below(90) == pytest.approx(1)
    assert m.sky_fraction_below(45) == pytest.approx(0.5 + 1 / math.pi)


def test_glare_domain_and_unit_scaling():
    assert m.veiling_point(100, 10, 70) == pytest.approx(20)
    for theta in (0.1, 1, 30, 90):
        with pytest.raises(ValueError):
            m.veiling_point(100, theta, 40)
    for transmittance in (0.15, 0.5, 1.0):
        assert m.contrast_retained(1000 * transmittance, 100 * transmittance) == pytest.approx(10 / 11)


def test_annulus_against_fine_midpoint_quadrature():
    # This catches degrees/radians confusion and a missing eye-plane cosine.
    edges = np.linspace(math.radians(5), math.radians(25), 100001)
    theta = (edges[:-1] + edges[1:]) / 2
    reference = np.sum(20 / np.degrees(theta) ** 2 * 2 * math.pi * np.sin(theta) * np.cos(theta) * np.diff(edges))
    assert m.annular_veil_per_luminance(70) == pytest.approx(reference, rel=1e-9)


def test_fresnel_normal_grazing_and_brewster_limits():
    n = m.WATER_REFRACTIVE_INDEX
    assert m.fresnel_water(1) == pytest.approx(((n - 1) / (n + 1)) ** 2)
    assert m.fresnel_water(0) == pytest.approx(1)
    # At Brewster incidence the parallel reflectance is zero.
    brewster = math.atan(n)
    expected = 0.5 * ((1 - n * n) / (1 + n * n)) ** 2
    assert m.fresnel_water(math.cos(brewster)) == pytest.approx(expected)


def test_equatorial_residence_and_no_interval_above_daily_maximum():
    assert m.interval_hours(0, 30) == pytest.approx(m.SYNODIC_MONTH_DAYS * 2)
    assert m.interval_hours(0, 90, period_hours=24) == pytest.approx(6)
    assert m.interval_hours(40, 60, latitude_deg=60) == 0


def test_inputs_schema_and_component_budget(tmp_path):
    bad = tmp_path / 'wrong.json'
    bad.write_text(json.dumps({'schema': 'other'}))
    with pytest.raises(ValueError):
        m.read_input(bad)
    rows = m.surface_rows(m.read_input())
    for r in rows:
        assert r['horizontal_total_lux'] == pytest.approx(r['horizontal_direct_lux'] + r['horizontal_diffuse_lux'])
        assert r['vertical_ambient_lux'] == pytest.approx(r['vertical_ground_lux'] + r['vertical_sky_lux'])
        if r['sun_deg'] <= 20:
            assert r['transfer_quality'] == 'two_stream_low_sun_caution'
    moon = next(r for r in rows if r['world'] == 'moon' and r['sun_deg'] == 90 and r['landscape_albedo'] == 0.1)
    earth = next(r for r in rows if r['world'] == 'earth' and r['sun_deg'] == 90 and r['landscape_albedo'] == 0.1)
    assert moon['horizontal_total_lux'] < earth['horizontal_total_lux']
    assert moon['vertical_ambient_lux'] > earth['vertical_ambient_lux']
    assert moon['vertical_ambient_lux'] == pytest.approx(20327.13)


def test_low_sun_proxy_keeps_components_and_albedo_separate():
    data = m.read_input()
    for r in m.low_sun_proxy(data):
        assert r['landscape_albedo'] == 0.1
        assert r['proxy_total_lux'] == pytest.approx(r['proxy_direct_lux'] + r['proxy_diffuse_lux'])
        assert r['source_atlas_hash']
    earth = next(r for r in m.low_sun_proxy(data) if r['world'] == 'earth' and r['sun_deg'] == 1)
    assert earth['proxy_total_lux'] == pytest.approx(149.075 + 1519.38, rel=3e-6)


def test_local_reflectance_does_not_change_scene_lighting():
    patches = m.local_surfaces(m.surface_rows(m.read_input()))
    patches = [p for p in patches if p['world'] == 'moon' and p['sun_deg'] == 90 and p['landscape_albedo'] == 0.1]
    assert len(patches) == 4
    for patch in patches:
        assert patch['horizontal_cd_m2'] / patch['patch_reflectance'] == pytest.approx(89602.6 / math.pi)


def test_downward_gaze_reverses_the_noon_comparison():
    gazes = m.gaze_scenes(m.surface_rows(m.read_input()))
    for a in (0.1, 0.8):
        def value(world, gaze):
            return next(r['total_lux'] for r in gazes if r['world'] == world and r['sun_deg'] == 90
                        and r['landscape_albedo'] == a and r['gaze_deg'] == gaze and r['azimuth_from_sun_deg'] == 0)
        assert value('moon', 0) > value('earth', 0)
        assert value('moon', -60) < value('earth', -60)


def test_stored_scene_products_reproduce_and_provenance_matches():
    stored = json.loads((m.HERE / 'results/optical_comfort.json').read_text())
    rows = m.surface_rows(m.read_input())
    for key, calculated in [('surface_scenes', rows), ('gaze_scenes', m.gaze_scenes(rows))]:
        assert len(stored[key]) == len(calculated)
        for actual, expected in zip(stored[key], calculated):
            assert set(actual) == set(expected)
            for field, value in expected.items():
                if isinstance(value, (float, int)):
                    assert actual[field] == pytest.approx(value, rel=1e-12, abs=1e-9)
                else:
                    assert actual[field] == value
    for path, digest in stored['producer']['files'].items():
        assert m.digest(m.ROOT / path) == digest
    assert stored['producer']['inputs'][str(m.INPUT.relative_to(m.ROOT))]['sha256'] == m.digest(m.INPUT)
