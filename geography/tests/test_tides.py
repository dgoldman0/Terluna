"""The tide product's parts: ephemeris, potential, harmonics, self-attraction, volume and dynamics."""
import math

import numpy as np

from geography import tide_dynamics as D
from geography import tides as T
from shared.constants import MOON_SURFACE_GRAVITY


def test_ephemeris_matches_jpl_horizons():
    check = T.check_ephemeris()
    assert check["samples"] > 1600
    # The physical librations, left out, stay under 0.05 degrees.
    for key in ("sub_earth_longitude_error_max_deg", "sub_earth_latitude_error_max_deg",
                "sub_solar_longitude_error_max_deg", "sub_solar_latitude_error_max_deg"):
        assert check[key] < 0.06, key
    assert check["earth_distance_error_max_km"] < 15


def test_monomials_reproduce_the_degree_2_and_3_potential():
    rng = np.random.default_rng(3)
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    points = rng.normal(size=(50, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    c2, c3 = 12.9, 0.058
    cosine = points @ direction
    direct = c2 * (1.5 * cosine**2 - 0.5) + c3 * (2.5 * cosine**3 - 1.5 * cosine)
    expanded = T.monomial_values(points) @ T.monomial_coefficients(direction, c2, c3)
    assert np.allclose(expanded, direct)


def test_harmonics_round_trip_on_a_coarse_grid():
    lat = np.arange(89.0, -90.0, -2.0)
    lon = np.arange(1.0, 360.0, 2.0)
    harmonics = T.Harmonics(lat, lon, 20)
    rng = np.random.default_rng(4)
    a = rng.normal(size=(21, 21)) * np.tri(21)
    b = rng.normal(size=(21, 21)) * np.tri(21)
    b[:, 0] = 0
    a2, b2 = harmonics.analyse(harmonics.synthesise(a, b))
    assert max(np.abs(a2 - a).max(), np.abs(b2 - b).max()) < 0.02 * np.abs(a).max()


def global_ocean(depth=1000.0, step=2.0):
    lat = np.arange(90 - step / 2, -90, -step)
    lon = np.arange(step / 2, 360, step)
    labels = np.ones((len(lat), len(lon)), dtype=int)
    return dict(labels=labels, depth=np.full(labels.shape, depth), lat=lat, lon=lon)


def test_self_attraction_of_a_degree_2_load_follows_its_factor():
    grid = global_ocean()
    response = T.Response(grid, [1], nmax=20, water_density=1000.0)
    la = np.radians(grid["lat"])[:, None] * np.ones((1, len(grid["lon"])))
    pattern = 1.5 * np.sin(la)**2 - 0.5
    attraction = response.self_attraction(pattern)
    factor = 3 * 1000.0 / (T.moon_mean_density() * 5)
    assert np.allclose(attraction, factor * pattern, atol=0.01 * factor)


def test_each_sea_keeps_its_volume():
    grid = global_ocean()
    grid["labels"] = np.zeros_like(grid["labels"])
    grid["labels"][10:30, 10:60] = 1
    grid["labels"][50:70, 100:150] = 2
    response = T.Response(grid, [1, 2], nmax=20, water_density=1000.0)
    la, lo = np.meshgrid(np.radians(grid["lat"]), np.radians(grid["lon"]), indexing="ij")
    forcing = np.cos(la)**2 * np.cos(2 * lo)
    height, iterations = response.solve(forcing)
    assert iterations > 1
    for mask in response.masks:
        assert abs(np.sum(height[mask] * response.area[mask])) < 1e-12 * np.sum(response.area[mask])
    assert np.all(height[~response.wet] == 0)


def test_lines_are_named_by_their_lunar_arguments():
    rates = T.argument_rates()
    named = T.identify(1 / 27.55455, rates, 1e-4)
    assert named["argument"] == "M'" and abs(named["argument_frequency_cpd"] - 1 / 27.55455) < 1e-5
    assert T.identify(1 / 31.812, rates, 1e-4)["argument"] == "2D -M'"


def rectangle(cells_x=24, cells_y=6, depth=400.0, step=1.0):
    lat = np.arange(90 - step / 2, -90, -step)
    lon = np.arange(step / 2, 360, step)
    labels = np.zeros((len(lat), len(lon)), dtype=int)
    rows = slice(int(90 / step) - cells_y // 2, int(90 / step) + cells_y // 2)
    labels[rows, 10:10 + cells_x] = 1
    return labels, np.full(labels.shape, depth), lat, lon


def test_the_gravest_seiche_of_a_rectangular_basin():
    labels, depth, lat, lon = rectangle()
    sea = D.SeaGrid(labels, depth, lat, lon, 1)
    periods, _ = D.seiche_periods(sea, 2)
    length = 24 * math.radians(1.0) * D.MOON_RADIUS * math.cos(math.radians(1.0))
    analytic = 2 * length / math.sqrt(MOON_SURFACE_GRAVITY * 400.0) / 3600
    assert abs(periods[0] / analytic - 1) < 0.03


def test_slow_forcing_settles_at_equilibrium():
    labels, depth, lat, lon = rectangle()
    sea = D.SeaGrid(labels, depth, lat, lon, 1)
    forcing = (np.cos(sea.lon) * np.cos(sea.lat)).astype(complex)
    result = D.compare(sea, forcing, 1 / 10000.0)
    assert result["departure_max_relative"] < 1e-3 and abs(result["rms_amplitude_ratio"] - 1) < 1e-3


def test_the_main_body_drops_cells_joined_only_at_a_corner():
    mask = np.zeros((10, 10), bool)
    mask[2:5, 2:5] = True
    mask[5, 5] = True          # touches the block only diagonally
    main = D.largest_component(mask)
    assert main.sum() == 9 and not main[5, 5]
