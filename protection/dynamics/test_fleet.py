import numpy as np
import pytest

from shared import constants as K
from .fleet import (safe_normals, beam_margins, square_vertices, projected_polygons,
                    coverage_snapshot, swept_conflicts, conflict_free_subset)
from .test_cycling import static_sample


def test_earth_guard_moves_the_reflected_solar_cone_off_earth():
    sample = static_sample()
    sample['positions']['earth'] *= -1
    q = np.array([[15e6, 0., 0.]])
    state = np.column_stack([q, np.zeros_like(q)])
    n, rays, changed = safe_normals(state, sample, 4*K.MOON_RADIUS, 10000.)
    assert changed[0]
    assert np.min(beam_margins(q, rays['sun'], n, sample, 4*K.MOON_RADIUS, 10000.)) > 0
    assert np.dot(n[0], -rays['sun'][0]) > 0
    assert np.linalg.norm(n[0]) == pytest.approx(1)


def test_finite_square_perspective_and_roll():
    pytest.importorskip('shapely')
    source = np.array([K.AU, 0., 0.]); q = np.array([[15e6, 0., 0.]])
    n = np.array([[-1., 0., 0.]])
    v = square_vertices(q, n, np.array([source]), 10000.)
    assert np.linalg.norm(v[0, 1]-v[0, 0]) == pytest.approx(10000.)
    poly = projected_polygons(v, source, np.array([1., 0., 0.]),
                               np.array([0., 1., 0.]), np.array([0., 0., 1.]))[0]
    assert poly.area == pytest.approx((10000.*K.AU/(K.AU-15e6))**2)


def test_overlap_is_a_union_and_seams_reduce_actual_shadow():
    pytest.importorskip('shapely')
    sample = static_sample(); q = np.array([[15e6, 0., 0.]])
    state = np.column_stack([q, np.zeros_like(q)])
    n = np.array([[-1., 0., 0.]])
    a = coverage_snapshot(state, n, sample, 1e6, 4*K.MOON_RADIUS, sun_count=4, seam=0, placement=0)
    b = coverage_snapshot(np.repeat(state, 2, axis=0), np.repeat(n, 2, axis=0), sample,
                          1e6, 4*K.MOON_RADIUS, sun_count=4, seam=0, placement=0)
    c = coverage_snapshot(state, n, sample, 1e6, 4*K.MOON_RADIUS, sun_count=4, seam=1000, placement=0)
    assert a['ray_coverage'] == pytest.approx(b['ray_coverage'])
    assert b['excess_intersection_area_fraction'] == pytest.approx(a['ray_coverage'])
    assert c['ray_coverage']/a['ray_coverage'] == pytest.approx(.999**2)


def test_finite_sun_removes_the_umbra_of_an_isolated_ten_km_tile():
    pytest.importorskip('shapely')
    sample = static_sample(); q = np.array([[15e6, 0., 0.]])
    state = np.column_stack([q, np.zeros_like(q)])
    result = coverage_snapshot(state, np.array([[-1., 0., 0.]]), sample,
                               10000., 4*K.MOON_RADIUS, sun_count=32)
    assert result['ray_coverage'] > 0
    assert result['all_sampled_sun_coverage'] == 0


def test_swept_guard_catches_collision_between_safe_endpoints():
    state = np.zeros((2, 2, 6))
    state[0, :, 0] = [-100., 100.]; state[1, :, 0] = [100., -100.]
    conflicts, _ = swept_conflicts([0., 1.], state, 1., acceleration_bound=0., clearance=1.)
    assert (0, 1) in conflicts
    assert conflicts[0, 1] == pytest.approx(0.)
    assert len(conflict_free_subset(2, conflicts)) == 1


def test_separated_tiles_can_shadow_each_other_without_colliding():
    state = np.zeros((2, 2, 6)); state[:, :, 0] = [15e6, 16e6]
    conflicts, shadows = swept_conflicts([0., 1.], state, 10000., acceleration_bound=0.,
        sun_directions=np.array([[1., 0., 0.], [1., 0., 0.]]))
    assert not conflicts
    assert (0, 1) in shadows
