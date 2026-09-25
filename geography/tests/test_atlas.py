"""The atlas pieces on synthetic topography, the gazetteer reader, and checks against the real inputs."""
from __future__ import annotations
import csv
import io
import math
import struct
import tempfile
import unittest
import zipfile
from pathlib import Path
import numpy as np

from geography import topography as tp, basins as bs, atlas, nomenclature as nm, fetch_inputs

HERE = Path(__file__).resolve().parents[1]


def inputs_present():
    try:
        return fetch_inputs.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


def synthetic(ppd=4):
    lat = 90.0 - (np.arange(180 * ppd) + 0.5) / ppd
    lon = (np.arange(360 * ppd) + 0.5) / ppd
    grid = tp.Grid(lat, lon, ppd)
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing='ij')

    def distance(lat0, lon0):
        c = np.sin(la) * math.sin(math.radians(lat0)) + np.cos(la) * math.cos(math.radians(lat0)) * np.cos(lo - math.radians(lon0))
        return np.degrees(np.arccos(np.clip(c, -1, 1)))
    return grid, distance


def dbf_bytes(fields, rows):
    """A minimal dBase III table: fields are (name, type, length, decimals)."""
    header_len = 32 + 32 * len(fields) + 1
    record_len = 1 + sum(f[2] for f in fields)
    out = struct.pack('<BBBBIHH20x', 3, 126, 9, 25, len(rows), header_len, record_len)
    for name, kind, length, dec in fields:
        out += struct.pack('<11sc4xBB14x', name.encode(), kind.encode(), length, dec)
    out += b'\r'
    for row in rows:
        out += b' '
        for (name, kind, length, dec), value in zip(fields, row):
            text = f'{value:>{length}.{dec}f}' if kind == 'N' else f'{value:<{length}}'
            out += text.encode('utf-8')[:length].ljust(length)
    return out + b'\x1a'


class LevelTests(unittest.TestCase):
    def test_level_for_area_of_a_flat_bowl(self):
        grid, distance = synthetic(2)
        area = grid.cell_area_m2
        bowl = distance(0.0, 180.0) < 30.0
        height = np.where(bowl, -1000.0, 0.0)
        cap = area[bowl].sum() / bs.AREA
        self.assertEqual(bs.level_for_area(height, area, 0.5 * cap), -1000.0)
        self.assertEqual(bs.level_for_area(height, area, 1.01 * cap), 0.0)


class BodyAndIslandTests(unittest.TestCase):
    def setUp(self):
        self.grid, distance = synthetic(4)
        self.area = self.grid.cell_area_m2
        self.bowl = distance(10.0, 90.0) < 25.0
        self.peak = distance(10.0, 90.0) < 3.0
        self.height = np.where(self.bowl, -2000.0, 500.0)
        self.height[self.peak] = 800.0

    def test_a_bowl_holds_one_sea_of_known_depth(self):
        water, lab, stats = atlas.bodies(self.height, self.area, 0.0, self.grid)
        self.assertEqual(len(water), 1)
        wet = self.bowl & ~self.peak
        self.assertAlmostEqual(water[0]['area_km2'], self.area[wet].sum() / 1e6, delta=1.0)
        self.assertEqual(water[0]['mean_depth_m'], 2000)
        self.assertEqual(water[0]['max_depth_m'], 2000)
        self.assertAlmostEqual(stats['water_share'], self.area[wet].sum() / bs.AREA, delta=1e-4)

    def test_a_central_peak_is_an_island_with_its_summit(self):
        land, lab, stats = atlas.islands(self.height, self.area, 0.0, self.grid)
        self.assertEqual(len(land), 1)
        self.assertEqual(land[0]['summit_m'], 800)
        self.assertAlmostEqual(land[0]['area_km2'], self.area[self.peak].sum() / 1e6, delta=1.0)
        self.assertAlmostEqual(stats['main_land_share'], self.area[~self.bowl].sum() / bs.AREA, delta=1e-4)

    def test_water_names_follow_their_centres(self):
        features = [dict(name='Mare Testium', code='ME', diameter_km=500.0, lat=10.0, lon=95.0),
                    dict(name='Lacus Siccus', code='LC', diameter_km=50.0, lat=-40.0, lon=0.0),
                    dict(name='Mons Insulae', code='MO', diameter_km=20.0, lat=10.0, lon=90.0)]
        water, _, _ = atlas.bodies(self.height, self.area, 0.0, self.grid, features)
        self.assertEqual(water[0]['water_names'], ['Mare Testium'])
        land, _, _ = atlas.islands(self.height, self.area, 0.0, self.grid, features)
        self.assertEqual(land[0]['features'], ['Mons Insulae'])

    def test_a_sea_across_the_seam_is_one_body(self):
        grid, distance = synthetic(4)
        height = np.where(distance(0.0, 0.0) < 20.0, -1000.0, 100.0)
        water, _, _ = atlas.bodies(height, grid.cell_area_m2, 0.0, grid)
        self.assertEqual(len(water), 1)


class NomenclatureTests(unittest.TestCase):
    def test_reads_the_gazetteer_table(self):
        fields = [('name', 'C', 30, 0), ('clean_name', 'C', 30, 0), ('diameter', 'N', 12, 3), ('center_lon', 'N', 12, 4),
                  ('center_lat', 'N', 12, 4), ('type', 'C', 20, 0), ('code', 'C', 4, 0), ('approval', 'C', 20, 0)]
        rows = [('Mare Testium', 'Mare Testium', 500.0, 350.0, 10.0, 'Mare, maria', 'ME', 'Adopted by IAU'),
                ('Crater X', 'Crater X', 12.5, 20.0, -5.0, 'Crater, craters', 'AA', 'Dropped')]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / nm.ARCHIVE
            with zipfile.ZipFile(path, 'w') as archive:
                archive.writestr('MOON.dbf', dbf_bytes(fields, rows))
                archive.writestr('MOON.cpg', 'UTF-8')
            adopted = nm.features(path)
            everything = nm.features(path, adopted_only=False)
        self.assertEqual([f['name'] for f in adopted], ['Mare Testium'])
        self.assertAlmostEqual(adopted[0]['lon'], -10.0)
        self.assertAlmostEqual(adopted[0]['diameter_km'], 500.0)
        self.assertEqual(len(everything), 2)


@unittest.skipUnless(inputs_present(), 'geography inputs not restored')
class RealTopographyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.height, _, cls.grid = tp.height_above_geoid(4, 200)
        cls.area = cls.grid.cell_area_m2

    def test_coarse_levels_match_the_fine_level_table(self):
        with open(HERE / 'results' / 'water_levels.csv') as handle:
            fine = {round(float(r['area_fraction']), 2): float(r['level_m']) for r in csv.DictReader(handle)}
        for share in (0.25, 0.35):
            self.assertAlmostEqual(bs.level_for_area(self.height, self.area, share), fine[share], delta=10.0)

    def test_the_near_side_maria_form_one_sea_at_the_scenario_share(self):
        features = nm.features()
        level = bs.level_for_area(self.height, self.area, atlas.scenario_share())
        water, _, _ = atlas.bodies(self.height, self.area, level, self.grid, features)
        names = set(water[0]['water_names'])
        for mare in ('Oceanus Procellarum', 'Mare Imbrium', 'Mare Serenitatis', 'Mare Crisium', 'Mare Nectaris'):
            self.assertIn(mare, names)

    def test_landing_heights_match_published_values(self):
        # Chang'e 4 and Chang'e 3 landing sites, published at -5,935 m and -2,640 m on the 1737.4 km sphere.
        _, sphere = atlas.point_heights(np.array([-45.444, 44.121]), np.array([177.599, -19.512]))
        self.assertAlmostEqual(float(sphere[0]), -5935.0, delta=40.0)
        self.assertAlmostEqual(float(sphere[1]), -2640.0, delta=40.0)


if __name__ == '__main__':
    unittest.main()
