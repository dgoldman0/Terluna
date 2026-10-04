"""The nearside-sea runner: a grid across the 0-degree meridian, periodic forcing and restartable segments."""
import numpy as np
import pytest

from climate.waves import nearside as N
from climate.waves.weather import apply_weights


def test_circular_crop_keeps_a_margin_on_both_sides_and_crosses_the_seam():
    occupied = np.zeros(360, bool)
    occupied[100:150] = True
    assert N.circular_crop(occupied).tolist() == list(range(99, 151))
    occupied = np.zeros(360, bool)
    occupied[275:] = occupied[:73] = True
    crop = N.circular_crop(occupied)
    assert crop[0] == 274 and crop[-1] == 73 and len(crop) == 160


def test_periodic_weights_interpolate_across_the_seam():
    lon = np.arange(0.0, 360.0, 5.625)
    lat = np.array([-10.0, 0.0, 10.0])
    water = np.ones((3, len(lon)), bool)
    field = np.broadcast_to(np.cos(np.radians(lon)), (3, len(lon)))[None]
    target_lon = np.array([[-2.0, 0.5, 357.0]])
    indices, weights, fallback = N.periodic_weights(lon, lat, water, target_lon, np.zeros_like(target_lon))
    assert not fallback.any() and np.allclose(weights.sum(axis=-1), 1)
    expected = np.interp(np.mod(target_lon, 360), np.r_[lon, 360.0], np.r_[np.cos(np.radians(lon)), 1.0])
    assert np.allclose(apply_weights(field, indices, weights)[0], expected)


def test_segments_cover_the_record_with_a_final_short_segment():
    bounds = N.segment_bounds()
    assert bounds[0] == (0, 48) and bounds[-1] == (1392, 1419) and len(bounds) == 30
    assert all(b == c for (_, b), (c, _) in zip(bounds[:-1], bounds[1:]))


def test_segment_input_keeps_the_smythii_physics_and_resumes_from_a_hotfile():
    sea = dict(body=433, lon=np.arange(-3.0, 3.0), lat=np.arange(-2.0, 2.0))
    first = N.segment_input(sea, 0, 48, step_s=300, frequency_intervals=48, water_density=1025.0, hotstart=False)
    later = N.segment_input(sea, 48, 96, step_s=300, frequency_intervals=48, water_density=1025.0, hotstart=True)
    for line in ("GEN3 KOMEN DRAG WU AGROW", "BREAKING CONSTANT 1.0 0.73", "PROP BSBT",
                 "CIRCLE 36 0.00496872696556 0.496872696556 48", "SET GRAV 1.62421887656 RHO 1025"):
        assert line in first and line in later
    assert "INIT ZERO" in first and "INIT HOTSTART SINGLE 'initial.hot' UNFORMATTED" in later
    # A resumed segment starts its output one hour (three for spectra) after its hotfile.
    assert "OUTPUT 20000103.010000 1 HR" in later and "OUTPUT 20000103.030000 3 HR" in later
    assert later.rstrip().endswith("HOTFILE 'end.hot' UNFORMATTED\nSTOP")


def test_the_nearside_sea_loads_across_the_seam():
    if not N.ATLAS.exists():
        pytest.skip("geography atlas product not restored")
    try:
        sea = N.load_sea()
    except FileNotFoundError:
        pytest.skip("IAU gazetteer input not restored")
    assert sea["lon"][0] < 0 < sea["lon"][-1] and np.allclose(np.diff(sea["lon"]), 1.0)
    assert sea["wet"].sum() == 7053
    assert all(r["depth_m"] >= N.REFERENCE_DEPTH_M and r["offset_km"] < 150 for r in sea["references"])


def test_time_above_matches_the_cycle_episode_count():
    from climate.waves import nearside_analysis as A
    from climate.waves.cycle import exceedance
    rng = np.random.default_rng(5)
    t = np.arange(60.0)
    heights = np.abs(rng.normal(1.0, 0.7, size=(60, 12)))
    for threshold in (0.5, 1.0, 2.0):
        expected = [exceedance(t, heights[:, k], threshold)["hours"] for k in range(12)]
        assert np.allclose(A.time_above(t, heights, threshold), expected)


def test_coastal_normals_point_to_the_land():
    from climate.waves import nearside_analysis as A
    wet = np.zeros((6, 8), bool)
    wet[1:5, 1:7] = True
    wet[2:4, 6] = False          # a notch of land on the east side
    nodes, normals = A.coastal_nodes(dict(wet=wet))
    rows, cols = np.nonzero(wet)
    lookup = {(r, c): n for r, c, n in zip(rows[nodes], cols[nodes], normals)}
    assert np.allclose(lookup[(2, 5)], (1, 0))            # land to the east
    assert np.allclose(lookup[(1, 3)], (0, -1))           # the southern edge (row 0 is south)
    assert (3, 3) not in lookup                           # an interior node
    assert np.allclose(np.hypot(*normals.T), 1)


def test_energy_mean_weights_periods_by_wave_energy():
    from climate.waves.nearside_analysis import energy_mean
    t = np.array([0.0, 1.0, 2.0])
    hs = np.array([[1.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    period = np.array([[10.0, 7.0], [10.0, 7.0], [16.0, 7.0]])
    # Trapezoids of E = Hs^2 and E T: (1 + 1)/2 + (1 + 4)/2 = 3.5 and (10 + 10)/2 + (10 + 64)/2 = 47.
    assert np.allclose(energy_mean(t, hs, period), [47 / 3.5, 0.0])


def test_landmark_near_names_the_closest_feature_with_its_distance():
    from climate.waves.nearside_analysis import landmark_near
    features = [dict(name="A", lon=0.0, lat=0.0), dict(name="B", lon=10.0, lat=0.0)]
    name, km = landmark_near(features, 8.0, 0.0)
    assert name == "B" and km == pytest.approx(N.MOON_RADIUS * np.radians(2.0) / 1000)
