"""Tests of the trial close-up look's own pieces (fine_look.py): the grain is bounded and has no seam, the dune crests
run along the parallels, the drainage count gathers downhill, and a window is the four tiles round its point; and of
fine_tiles.py's water-only tile: land keeps the land's colour, deep water the open water's, the shallows the coast's."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine_look  # noqa: E402
import fine_tiles  # noqa: E402


def points(lat, lon):
    la, lo = np.meshgrid(lat, lon, indexing='ij')
    return fine_look.sphere_points(la, lo)


def test_grain_is_bounded_and_has_no_seam_at_the_first_meridian():
    lat = np.linspace(-10, 10, 64)
    n = fine_look.fbm(points(lat, np.array([-1e-4, 1e-4, 359.9999])), 6.0, 11)
    assert 0.0 <= float(n.min()) and float(n.max()) <= 1.0
    assert np.abs(n[:, 0] - n[:, 1]).max() < 1e-3 and np.abs(n[:, 0] - n[:, 2]).max() < 1e-3


def test_dune_crests_run_along_the_parallels():
    step = 0.03                                             # degrees: about 0.9 km
    c = fine_look.crests(points(10.0 - step * np.arange(200), 40.0 + step * np.arange(200)), 3.5, 20.0, 53)
    assert -1.0 <= float(c.min()) and float(c.max()) <= 1.0
    north = np.abs(np.diff(c, axis=0)).mean()
    east = np.abs(np.diff(c, axis=1)).mean()
    assert north > 2.0 * east                               # the pattern changes faster going north than going east


def test_drainage_gathers_down_a_valley_and_stops_at_water():
    y, x = np.mgrid[0:20, 0:21].astype(np.float32)
    surf = 10.0 * np.abs(x - 10.0) + 2.0 * (20.0 - y)       # a valley along column 10, falling to the last row
    wet = np.zeros(surf.shape, bool)
    wet[19] = True
    acc = fine_look.drainage(surf, wet)
    assert acc[18, 10] > acc[5, 10] > acc[5, 3]
    assert acc[18, 10] > 100                                # most of the 19 x 21 dry cells pass the valley's foot


def test_a_window_is_the_four_tiles_round_its_point():
    tiles = fine_look.window_tiles(14.0, 2.0)
    assert len(set(tiles)) == 4
    rows = sorted({t[0] for t in tiles})
    cols = [t[1] for t in tiles[:2]]
    assert rows[1] == rows[0] + 1
    top = fine_tiles.LAT_TOP - fine_tiles.TILE_DEG * rows[0]
    assert top >= 14.0 >= top - 2 * fine_tiles.TILE_DEG
    assert cols == [fine_tiles.TILE_COLS - 1, 0]            # the window crosses the first meridian


def test_a_water_only_tile_colours_land_deep_water_and_shallows():
    import appearance
    T = fine_tiles.TILE
    r0 = (90 - fine_tiles.LAT_TOP) * fine_tiles.PPD                  # the first tile's corner
    rows, cols = r0 + T, T
    height = np.zeros((rows, cols), np.int16)                        # half metres
    level = np.zeros((rows, cols), np.int16)
    water = np.zeros((rows, cols), np.uint8)
    water[r0:, :200] = 1                                             # deep water: 100 m
    height[r0:, :200] = -200
    water[r0:, 200:300] = 1                                          # shallows: 1 m
    height[r0:, 200:300] = -2
    height[r0:, 300:] = 100                                          # land
    land = np.full((2880, 5760, 3), 0.25, np.float32)
    waters = dict(open_sea=np.array([0.01, 0.02, 0.05]), productive_coast=np.array([0.02, 0.06, 0.04]),
                  river_mouth=np.array([0.2, 0.15, 0.1]))
    saved = dict(fine_tiles._G)
    try:
        fine_tiles._G.update(height=height, level=level, water=water, land=land, plume=np.zeros((2880, 5760), np.float32),
                             waters=waters, rotation=np.eye(3))
        t = fine_tiles.water_tile(0, 0)
    finally:
        fine_tiles._G.clear()
        fine_tiles._G.update(saved)
    shade = lambda c: np.round(appearance.srgb_encode(np.asarray(c, float)) * 255).astype(np.uint8)  # noqa: E731
    assert t['wet'][:, :300].all() and not t['wet'][:, 300:].any()
    assert (t['colour'][10, 100] == shade(waters['open_sea'])).all()
    assert (t['colour'][10, 250] == shade(waters['productive_coast'])).all()
    assert (t['colour'][10, 400] == shade([0.25, 0.25, 0.25])).all()
