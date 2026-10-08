import numpy as np
import pytest
import shapely
from scipy.spatial.transform import Rotation

from .active_formation import lattice, solar_frame
from .natural_pattern import (rectangle_union, mutual_visibility,
    neighbour_exclusion, transported_roll, local_coverage)
from .test_cycling import static_sample


def test_rectangle_union_matches_independent_polygon_boolean():
    rng = np.random.default_rng(49021)
    lo = rng.uniform(-1, 1, (12, 8, 2))
    hi = lo+rng.uniform(0, 2, lo.shape)
    active = rng.random((12, 8)) > .2
    expected = []
    for low, high, keep in zip(lo, hi, active):
        shapes = shapely.box(low[keep, 0], low[keep, 1], high[keep, 0], high[keep, 1])
        expected.append(shapely.union_all(shapes).area)
    np.testing.assert_allclose(rectangle_union(lo, hi, active), expected, atol=1e-12)
    np.testing.assert_allclose(rectangle_union(lo, hi, np.zeros_like(active)), 0.)


def test_shadow_neighbour_reduction_matches_all_tiles_and_detects_omission():
    sample = static_sample(); frame = solar_frame(sample)
    offsets, neighbour, valid, _ = lattice(5, 8700., 3000.)
    state = np.c_[np.array([15e6, 0., 0.])+offsets@frame.T, np.zeros((25, 3))]
    reduced, _ = mutual_visibility(state, frame, sample, neighbour, valid)
    every = np.tile(np.arange(25), (25, 1))
    all_valid = ~np.eye(25, dtype=bool)
    full, _ = mutual_visibility(state, frame, sample, every, all_valid)
    np.testing.assert_allclose(reduced, full, atol=1e-12)
    assert reduced.min() < 1 and reduced.max() == pytest.approx(1.)
    assert neighbour_exclusion(state, frame, sample, neighbour, valid) > 0
    # Force a non-neighbour tile into the solar cylinder. The guard must reject.
    state[24, :3] = state[0, :3]+frame[:, 0]*1000
    assert neighbour_exclusion(state, frame, sample, neighbour, valid) < 0


def test_transported_roll_removes_prescribed_spin_without_changing_plane():
    times = np.arange(21.)
    orientations = Rotation.from_euler('zy', np.c_[times*.1, times*.01]).as_matrix()
    spline = transported_roll(times, orientations)
    frames = spline(times).as_matrix()
    np.testing.assert_allclose(frames[:, :, 2], orientations[:, :, 2], atol=1e-12)
    assert np.linalg.norm(spline(times, 1), axis=1).max() < .02
    assert np.isfinite(spline(times, 2)).all()


def test_local_finite_sun_geometry_finds_small_pattern_edge_failure():
    sample = static_sample()
    centre = np.array([15e6, 0., 0., 0., 0., 0.])
    # A central tile covers central rays, but cannot cover the solar limb's
    # ~70 km displacement at a 15 Mm intercept distance.
    result = local_coverage(centre[None, :], sample, centre, radius=1000.)
    assert result['all_sampled_sun_coverage'] == 0
    assert result['worst_missed_ray'] is not None
    assert result['minimum_source_coverage'] == pytest.approx(0., abs=1e-12)
