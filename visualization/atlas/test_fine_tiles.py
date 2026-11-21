"""Tests of the close views' tiles (fine_tiles.py): the tile-window versions of the page's shading and sky agree with
appearance.py's, and the 16 px/deg colour is read where it lies."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appearance  # noqa: E402
import fine_tiles  # noqa: E402


def test_normals_match_the_page_build():
    rng = np.random.default_rng(3)
    h = np.cumsum(rng.normal(0, 200, (180, 360)), axis=1).astype(np.float32)
    lat = 90.0 - (np.arange(180) + 0.5)
    assert np.array_equal(fine_tiles.normals(h, lat, 1), appearance.relief_normals(h, 1))


def test_sky_matches_the_page_build_away_from_the_edges():
    rng = np.random.default_rng(4)
    h = (rng.normal(0, 3000, (180, 360))).astype(np.float32)
    water = np.zeros(h.shape, bool)
    water[80:90, 100:120] = True
    lat = 90.0 - (np.arange(180) + 0.5)
    mine = fine_tiles.sky(h.copy(), lat, 1, water)
    page = appearance.sky_view(h.copy(), 1, water)
    inner = (slice(20, 160), slice(20, 340))
    assert np.allclose(mine[inner], page[inner], atol=1e-5)


def test_sample16_reads_cell_centres():
    field = np.arange(2880 * 5760, dtype=np.float32).reshape(2880, 5760) % 997
    i, j = 1234, 4321
    lat = np.array([[90.0 - (i + 0.5) / 16.0]], np.float32)
    lon = np.array([[(j + 0.5) / 16.0 - 180.0]], np.float32)
    assert abs(float(fine_tiles.sample16(field, lat, lon)[0, 0]) - field[i, j]) < 1e-3
