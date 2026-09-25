"""ExoPlaSim runner: surface files, the sea mask, pruning and process bookkeeping (no model run)."""
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np

from climate.gcm import exoplasim_run as er


class RunnerTests(unittest.TestCase):
    def test_surface_file_layout(self):
        field = np.arange(32 * 64, dtype=float).reshape(32, 64) / 7.0
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'f.sra'
            er._write_sra(path, 172, field)
            lines = path.read_text().splitlines()
        self.assertEqual([int(v) for v in lines[0].split()], [172, 0, 20090101, 0, 64, 32, 0, 0])
        self.assertEqual(len(lines), 1 + 32 * 64 // 8)
        values = np.array([float(v) for line in lines[1:] for v in line.split()]).reshape(32, 64)
        np.testing.assert_allclose(values, field, atol=5e-4)       # three decimals, north to south from 0 E

    def test_sea_mask_covers_the_water_share(self):
        rng = np.random.default_rng(3)
        land = rng.uniform(0.0, 1.0, (32, 64))
        area = np.cos(np.radians(np.linspace(-85.0, 85.0, 32)))[:, None] * np.ones(64)
        mask = er.land_mask(land, area, 0.25)
        sea_share = float(((1 - mask) * area).sum() / area.sum())
        self.assertGreaterEqual(sea_share, 0.25)
        self.assertLess(sea_share - 0.25, float(area.max() / area.sum()))
        # Every sea cell has less land than every land cell.
        self.assertLessEqual(land[mask == 0].max(), land[mask == 1].min())

    def test_spectrum_tail_follows_the_sun(self):
        # A blackbody "measurement" that stops at 2.73 um must continue as the same blackbody, not as a
        # constant: a constant tail would move most of the energy into ExoPlaSim's near-infrared band.
        h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
        planck = lambda um: 1.0 / ((um * 1e-6) ** 5 * np.expm1(h * c / (um * 1e-6 * k * 5772.0)))
        w_nm = np.linspace(202.0, 2730.0, 3000)
        grid = np.concatenate([np.geomspace(0.2, 0.75, 1025)[:-1], np.geomspace(0.75, 100.0, 1024)])
        f = er.filtered_spectrum(grid, w_nm, planck(w_nm / 1000.0) / 1000.0, [0.0, 1e6], [1.0, 1.0])
        inside = grid >= 0.202
        np.testing.assert_allclose(f[inside], planck(grid[inside]), rtol=2e-3)
        energy = lambda sel: float(np.trapezoid(f[sel], grid[sel]))
        # A 5772 K blackbody puts 54% of its energy below 0.75 um (lambda T = 4330 um K).
        self.assertAlmostEqual(energy(grid <= 0.75) / energy(grid > 0), 0.54, delta=0.01)

    def test_pruning_keeps_recent_and_every_tenth_year(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for y in range(25):
                (work / f'MOST_REST.{y:05d}').touch()
                (work / f'MOST.{y:05d}.nc').touch()
            er._prune(work, 24)
            rest = sorted(int(p.name.split('.')[1]) for p in work.glob('MOST_REST.*'))
            out = sorted(int(p.name.split('.')[1]) for p in work.glob('MOST.?????.nc'))
        self.assertEqual(rest, [9, 19, 22, 23, 24])
        self.assertEqual(out, [9, 19, 20, 21, 22, 23, 24])

    def test_partial_year_is_cleared(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for name in ('plasim_output', 'plasim_status', 'MOST.00003', 'MOST_DIAG.00003', 'MOST.00002.nc',
                         'MOST_REST.00002'):
                (work / name).touch()
            er._clean_partial(work, 3)
            self.assertEqual(sorted(p.name for p in work.iterdir()), ['MOST.00002.nc', 'MOST_REST.00002'])

    def test_process_lookup_finds_this_process(self):
        self.assertIn(os.getpid(), er._members(os.getpgrp()))
        self.assertIn(os.getpid(), er._using(Path(os.getcwd())))
        self.assertTrue(er._alive(os.getpid()))


if __name__ == '__main__':
    unittest.main()
