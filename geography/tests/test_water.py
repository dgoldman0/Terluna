"""Geoid synthesis, level curves and the depression tree on synthetic and real topography."""
from __future__ import annotations
import math
import unittest
import numpy as np

from geography import topography as tp, basins as bs, fetch_inputs


def inputs_present():
    try:
        return fetch_inputs.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


def synthetic(ppd=2):
    lat = 90.0 - (np.arange(180 * ppd) + 0.5) / ppd
    lon = (np.arange(360 * ppd) + 0.5) / ppd
    grid = tp.Grid(lat, lon, ppd)
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing='ij')

    def distance(lat0, lon0):
        c = np.sin(la) * math.sin(math.radians(lat0)) + np.cos(la) * math.cos(math.radians(lat0)) * np.cos(lo - math.radians(lon0))
        return np.degrees(np.arccos(np.clip(c, -1, 1)))
    return grid, distance


class GeoidTests(unittest.TestCase):
    def test_legendre_orthonormal(self):
        x, w = np.polynomial.legendre.leggauss(300)
        u = np.sqrt(1 - x * x)
        for m in (0, 1, 13, 120):
            p = tp.legendre_column(m, 150, x, u)
            gram = (p * w) @ p.T * (0.5 if m == 0 else 0.25)
            self.assertLess(np.abs(gram - np.eye(gram.shape[0])).max(), 1e-10)

    def test_p20_closed_form(self):
        t = np.linspace(-1, 1, 11)
        p = tp.legendre_column(0, 2, t, np.sqrt(1 - t * t))
        self.assertTrue(np.allclose(p[2], math.sqrt(5) * (3 * t * t - 1) / 2))

    def test_cell_areas_cover_the_sphere(self):
        grid, _ = synthetic(3)
        self.assertAlmostEqual(grid.cell_area_m2.sum() / bs.AREA, 1.0, delta=1e-12)


class StorageTests(unittest.TestCase):
    def test_level_curve_of_a_flat_bowl(self):
        grid, distance = synthetic()
        area = grid.cell_area_m2
        bowl = distance(0.0, 180.0) < 30.0
        height = np.where(bowl, -1000.0, 0.0)
        frac, vol = bs.level_curve(height, area, [-500.0])
        cap = area[bowl].sum()
        self.assertAlmostEqual(frac[0], cap / bs.AREA, delta=1e-12)
        self.assertAlmostEqual(vol[0] / (500.0 * cap), 1.0, delta=1e-12)
        level = bs.level_for_volume(height, area, 250.0 * cap)
        self.assertAlmostEqual(level, -750.0, delta=1e-3)

    def test_two_bowls_spill_at_the_saddle(self):
        grid, distance = synthetic()
        area = grid.cell_area_m2
        deep = distance(0.0, 100.0) < 20.0
        shallow = distance(0.0, 150.0) < 20.0
        saddle = (distance(0.0, 125.0) < 6.0)
        height = np.full(deep.shape, 1000.0)
        height[deep] = -2000.0
        height[shallow] = -800.0
        height[saddle & ~deep & ~shallow] = 300.0
        records = bs.depression_tree(height, area, grid.lat_deg, grid.lon_deg)
        joins = bs.major_merges(records, 1e-3)
        level, parts = joins[0]
        self.assertAlmostEqual(level, 300.0)
        vols = sorted(d.volume_m3 for d in parts)
        self.assertAlmostEqual(vols[0] / ((300.0 + 800.0) * area[shallow].sum()), 1.0, delta=1e-9)
        self.assertAlmostEqual(vols[1] / ((300.0 + 2000.0) * area[deep].sum()), 1.0, delta=1e-9)
        sizes, _ = bs.seas(height, area, 0.0)
        self.assertEqual(len(sizes), 2)
        sizes, _ = bs.seas(height, area, 500.0)
        self.assertEqual(len(sizes), 1)


@unittest.skipUnless(inputs_present(), 'topography inputs not restored (python -m geography.fetch_inputs --download)')
class LunarTests(unittest.TestCase):
    def test_geoid_shape_and_topography(self):
        height, geoid, grid = tp.height_above_geoid(4)
        area = grid.cell_area_m2
        self.assertAlmostEqual((geoid * area).sum() / area.sum(), 0.0, delta=0.5)
        # The Moon's flattening (J2) dominates: high equator, low poles.
        self.assertGreater(geoid[np.abs(grid.lat_deg) < 1].mean(), 150.0)
        self.assertLess(geoid[np.abs(grid.lat_deg) > 89].mean(), -250.0)
        i = np.unravel_index(np.argmin(height), height.shape)
        self.assertLess(grid.lat_deg[i[0]], -60.0)                 # South Pole-Aitken floor
        self.assertTrue(170.0 < grid.lon_deg[i[1]] < 200.0)
        self.assertTrue(-9000.0 < height.min() < -7500.0 and 9000.0 < height.max() < 11000.0)


if __name__ == '__main__':
    unittest.main()
