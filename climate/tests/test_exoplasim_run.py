"""ExoPlaSim runner: surface files, the sea mask, settings, pruning and process bookkeeping (no model run)."""
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

    def test_cloud_water_settings(self):
        self.assertEqual(er.cloud_water_namelist('plasim', 1.62), (0.0, 1.0))           # PlaSim unchanged
        gref, scale = er.cloud_water_namelist('earth_path', 1.62)
        self.assertEqual(gref, er.STANDARD_GRAVITY)
        self.assertAlmostEqual(scale, 1.62 / er.STANDARD_GRAVITY)                        # Earth's water path
        self.assertEqual(er.cloud_water_namelist('full_column', 1.62), (er.STANDARD_GRAVITY, 1.0))
        with self.assertRaises(ValueError):
            er.cloud_water_namelist('thick', 1.62)

    def test_patches_apply_in_order_with_distinct_markers(self):
        markers = [marker for _, marker, _ in er.PATCHES]
        self.assertEqual(len(markers), len(set(markers)))
        # The convective patch edits text the cloud-water patch wrote, so it must come after it.
        order = [marker for name, marker, _ in er.PATCHES if name == 'rainmod.f90']
        self.assertEqual(order[0], 'Terluna: gravity for the heights below')
        # The Sun's clock extends the namelist and broadcast lines the calibration patch wrote, the shortwave
        # path extends the Sun's, and the infrared from space extends the shortwave path's.
        radmod = [marker for name, marker, _ in er.PATCHES if name == 'radmod.f90']
        self.assertEqual(radmod, ['Terluna: calibration multiplier', "Terluna: the Sun's clock", 'Terluna: calibrated absorptivity',
                                  'Terluna: infrared from space at the top of the air'])

    def test_default_settings_are_the_corrected_physics(self):
        self.assertEqual(er.MODEL['cloud_water'], 'earth_path')
        self.assertEqual(er.MODEL['convective_day_s'], 86400.0)
        self.assertEqual(er.MODEL['surface_cloud_rcrit'], 0.0)       # PlaSim's own threshold near the ground
        self.assertEqual(er.MODEL['sunlight_scale'], 1.0)                 # the shield product as it is
        self.assertEqual(er.MODEL['sun_clock'], 'synodic')               # an even Sun, one turn per 29.53 days
        self.assertEqual((er.MODEL['fixed_orbit'], er.MODEL['whole_lunar_days']), (1, 1))
        self.assertEqual(er.SUN_CLOCKS.index('plasim'), 0)                # PlaSim's own clock is the default off
        # Since 2026-10-14 the curved air and the CO2 and O2 bands are the standard, with the multipliers refitted
        # with them.
        import inspect
        for f in (er.run, er.configuration):
            params = inspect.signature(f).parameters
            self.assertEqual((params['twilight'].default, params['gas_bands'].default), ('full', True))
        refit = er.ROOT / 'climate' / 'results' / 'gcm' / 'shortwave_fit_A28_dim5_lunar.json'
        if refit.exists():
            import json
            scale = json.loads(refit.read_text())['forms']['scale']
            self.assertEqual((er.MODEL['sw_water_absorption'], er.MODEL['rayleigh_scale']),
                             (round(scale['k'], 3), round(scale['rayleigh_scale'], 3)))

    def test_clear_sky_output(self):
        off = er.output_variables({**er.MODEL, 'clear_sky': 0})
        on = er.output_variables({**er.MODEL, 'clear_sky': 1})
        self.assertEqual(off, er.OUTPUT)
        self.assertEqual(on[:len(off)], off)
        self.assertTrue({'rstcs', 'rltcs'} <= set(on))

    def test_patch_markers_mark_patched_text_only(self):
        # A marker must appear in what a patch adds and never in what it replaces, or it would apply twice.
        for _, marker, edits in er.PATCHES + er.TWILIGHT_PATCHES:
            self.assertTrue(any(marker in new for _, new in edits))
            self.assertFalse(any(marker in old for old, _ in edits))

    def test_the_twilight_patch_stays_out_of_the_shared_model(self):
        markers = {marker for _, marker, _ in er.PATCHES}
        self.assertFalse(markers & {marker for _, marker, _ in er.TWILIGHT_PATCHES})
        # It extends the namelist and broadcast lines the infrared from space wrote, and keeps PlaSim's lines short.
        old = ''.join(o for o, _ in er.TWILIGHT_EDITS)
        self.assertIn('irtopnight', old)
        for _, new in er.TWILIGHT_EDITS:
            self.assertTrue(all(len(line) <= 132 for line in new.splitlines()))
        self.assertEqual(er.TWILIGHT_MODES, ('curved', 'full', 'path'))
        added = ''.join(new for _, new in er.TWILIGHT_EDITS)
        self.assertIn('call twilight_term', added)                      # the sunlight, at the end of swr
        self.assertIn('call twilight_lw', added)                        # the infrared, at the end of lwr

    def test_the_gas_bands_follow_the_twilight_patch_and_the_fitter(self):
        from climate.gcm import plasim_shortwave as ps
        # They extend the module, namelist and broadcast lines the twilight patch wrote and the line the calibrated
        # absorptivity wrote, and keep PlaSim's lines short.
        twilight = ''.join(new for _, new in er.TWILIGHT_EDITS)
        anchors = [old for old, _ in er.GAS_EDITS]
        for old in anchors[:3]:
            self.assertIn(old, twilight)
        self.assertIn(anchors[6], ''.join(new for _, new in er.SHORTWAVE_EDITS))
        for _, new in er.GAS_EDITS:
            self.assertTrue(all(len(line) <= 132 for line in new.splitlines()))
        # ZGASQ's columns are the fitter's bands in order: CO2 and O2 beside the water, then O2 in the visible.
        self.assertEqual(er.GAS_ORDER, tuple(name for name, _, _ in ps.GAS_BANDS))
        self.assertEqual([band for _, _, band in ps.GAS_BANDS], [2, 2, 1])
        # Earth's air: O2 is 23.1% of its mass and CO2 at 400 ppm 0.061%.
        earth = er.gas_mass_fractions(dict(pN2=0.7808, pO2=0.2095, pAr=0.00934, pCO2=0.0004))
        self.assertAlmostEqual(earth['O2'], 0.2314, places=3)
        self.assertAlmostEqual(earth['CO2'] * 1e4, 6.08, places=1)

    def test_the_longwave_correction_follows_the_twilight_and_gas_patches(self):
        # It extends the module, namelist, broadcast and lwr lines the twilight and gas patches wrote, keeps PlaSim's
        # lines short, and runs after the curved air's escape.
        written = ''.join(new for _, new in er.TWILIGHT_EDITS + er.GAS_EDITS)
        for old, _ in er.LW_CORRECTION_EDITS:
            self.assertIn(old, written)
        for _, new in er.LW_CORRECTION_EDITS:
            self.assertTrue(all(len(line) <= 132 for line in new.splitlines()))
        added = ''.join(new for _, new in er.LW_CORRECTION_EDITS)
        self.assertLess(added.index('call twilight_lw'), added.index('call lw_correction'))

    def test_ice_optics_follow_the_longwave_patch(self):
        # They extend the module, namelist and broadcast lines the longwave patch wrote and the subroutine it added,
        # keep PlaSim's lines short, and act only where NICEOPT asks, after PlaSim has set the liquid optics.
        written = ''.join(new for _, new in er.LW_CORRECTION_EDITS)
        anchors = [old for old, _ in er.ICE_OPTICS_EDITS]
        for old in anchors[:3] + anchors[5:]:
            self.assertIn(old, written)
        for _, new in er.ICE_OPTICS_EDITS:
            self.assertTrue(all(len(line) <= 132 for line in new.splitlines()))
        sw, lw = er.ICE_OPTICS_EDITS[3][1], er.ICE_OPTICS_EDITS[4][1]
        self.assertLess(sw.index('endwhere'), sw.index('if (niceopt == 1) call ice_cloud_sw'))
        self.assertLess(lw.index('endif'), lw.index('if (niceopt == 1 .and. clgray <= 0.) call ice_cloud_lw'))
        markers = [marker for _, marker, _ in er.TWILIGHT_PATCHES]
        self.assertEqual(markers[3], 'Terluna: ice optics for cold cloud')

    def test_the_twilight_clear_share_follows_the_patches_it_extends(self):
        # It extends the module, namelist and broadcast lines and the subroutines the twilight, longwave and ice patches
        # wrote, keeps PlaSim's lines short, and comes last.
        written = ''.join(new for _, _, edits in er.TWILIGHT_PATCHES[:-1] for _, new in edits)
        for old, _ in er.TWILIGHT_CLEAR_EDITS:
            self.assertIn(old, written)
        for _, new in er.TWILIGHT_CLEAR_EDITS:
            self.assertTrue(all(len(line) <= 132 for line in new.splitlines()))
        self.assertIs(er.TWILIGHT_PATCHES[-1][2], er.TWILIGHT_CLEAR_EDITS)
        self.assertIn('ntwlclear', ''.join(new for _, new in er.TWILIGHT_CLEAR_EDITS))

    def test_the_twilight_clear_share_patch_leaves_the_model_as_it_was_when_off(self):
        # The patched twilight subroutines, compiled beside the originals on a stand-in module: with NTWLCLEAR 0 they give
        # the originals' fluxes bit for bit; with 1 each layer's sunlight comes in the cover clear above it and each
        # layer's infrared leaves through the longwave transmissivity of the clouds above it.
        import shutil
        import subprocess
        if shutil.which('gfortran') is None:
            self.skipTest('no gfortran')
        sw_new, lw_new = er.TWILIGHT_EDITS[-2][1], er.TWILIGHT_EDITS[-1][1]
        cut = lambda text, name: text[text.index(f'      subroutine {name}\n'):text.index(f'      end subroutine {name}\n')
                                      + len(f'      end subroutine {name}\n')]
        original = cut(sw_new, 'twilight_term') + cut(lw_new, 'twilight_lw')
        patched = original
        for old, new in er.TWILIGHT_CLEAR_EDITS:
            if old in patched:
                patched = patched.replace(old, new)
        self.assertNotEqual(patched, original)
        renamed = patched.replace('twilight_term', 'clear_term').replace('twilight_lw', 'clear_lw')
        module = '''
      module radmod
      implicit none
      integer, parameter :: NHOR = 3, NLEV = 3, NLEP = 4, NTWLMAX = 2
      real, parameter :: PI = 3.14159265358979323846
      real :: gsol0 = 1200.0, gdist2 = 1.0
      integer :: nclouds = 1, ntwl = 2, ntwlclear = 0
      real :: twlel(NTWLMAX) = (/-10.0, 10.0/), twlg(NTWLMAX) = (/0.01, 0.03/)
      real :: twlq(NLEV,NTWLMAX) = reshape((/0.002, 0.004, 0.006, 0.003, 0.005, 0.007/), (/NLEV, NTWLMAX/))
      real :: twlw(NLEV+1) = (/4.0, 2.0, 1.0, 0.5/), twlt(NLEV+1) = (/200.0, 250.0, 290.0, 300.0/)
      real :: gmuz(NHOR), dcc(NHOR,NLEV), dt(NHOR,NLEP), twltau(NHOR,NLEV)
      real :: dswfl(NHOR,NLEP), dfd(NHOR,NLEP), dftu(NHOR,NLEP), dlwfl(NHOR,NLEP)
      end module radmod
'''
        driver = '''
      program drive
      use radmod
      implicit none
      integer :: j, k
      gmuz(:) = (/0.05, -0.02, 0.1/)
      dcc(1,:) = (/0.0, 0.0, 0.0/)
      dcc(2,:) = (/0.5, 0.0, 0.2/)
      dcc(3,:) = (/0.3, 0.6, 0.1/)
      twltau(1,:) = (/1.0, 1.0, 1.0/)
      twltau(2,:) = (/0.8, 1.0, 0.9/)
      twltau(3,:) = (/0.85, 0.5, 0.95/)
      do j = 1, NHOR
       dt(j,:) = (/205.0, 252.0, 288.0, 301.0/) + 2.0 * j
      enddo
      do k = 0, 2
       dswfl(:,:) = 0.
       dfd(:,:) = 0.
       dftu(:,:) = 0.
       dlwfl(:,:) = 0.
       if (k == 0) then
        call twilight_term
        call twilight_lw
       else
        ntwlclear = k - 1
        call clear_term
        call clear_lw
       endif
       do j = 1, NHOR
        print '(16es26.17)', dswfl(j,:), dfd(j,:), dftu(j,:), dlwfl(j,:)
       enddo
      enddo
      end program drive
'''
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'clear.f90').write_text(module + original + renamed + driver)
            subprocess.run(['gfortran', '-fdefault-real-8', '-O0', '-o', 'clear', 'clear.f90'], cwd=tmp, check=True)
            out = np.loadtxt(subprocess.run(['./clear'], cwd=tmp, check=True, capture_output=True,
                                            text=True).stdout.split('\n')).reshape(3, 3, 4, 4)
        before, off, on = out
        self.assertTrue(np.array_equal(before, off))                          # NTWLCLEAR 0: the model as it was
        self.assertTrue(np.array_equal(on[0], before[0]))                     # a cloud-free column is unchanged
        dcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.2], [0.3, 0.6, 0.1]])
        tau = np.array([[1.0, 1.0, 1.0], [0.8, 1.0, 0.9], [0.85, 0.5, 0.95]])
        for j in (1, 2):
            above = np.r_[1.0, np.cumprod(1.0 - dcc[j])[:-1]]                  # the cover clear above each layer
            self.assertTrue(np.allclose(np.diff(on[j][0]), np.diff(before[j][0]) * above, rtol=1e-12, atol=1e-14))
            self.assertEqual(on[j][0][-1], before[j][0][-1])                  # the ground's kept the column's clear share
            escape = lambda f: np.r_[np.diff(f), -f[-1]]                       # each layer's, top first, then the ground's
            lw_above = np.r_[1.0, np.cumprod(tau[j])]
            self.assertTrue(np.allclose(escape(on[j][2]), escape(before[j][2]) * lw_above, rtol=1e-12, atol=1e-14))
            self.assertTrue(np.allclose(on[j][3], on[j][2]))                  # dlwfl takes the same upward change

    def test_the_output_asks_only_for_what_plasim_writes(self):
        # The surface stress joins the output; the 10 m wind, which PlaSim 3.4.2 never writes, stays out, since the
        # postprocessor would leave it out of the file without an error.
        self.assertTrue({'tauu', 'tauv'} <= set(er.OUTPUT))
        self.assertFalse({'uas', 'vas'} & set(er.OUTPUT))
        self.assertTrue({180, 181} <= set(er.PLASIM_CODES) and not {165, 166} & set(er.PLASIM_CODES))
        import re
        source = er.TWILIGHT_PACKAGE / 'exoplasim'
        if not (source / 'pyburn.py').exists():
            self.skipTest('the twilight copy of ExoPlaSim is not at hand')
        names = {m.group(2): int(m.group(1))
                 for m in re.finditer(r'"(\d+)":\["(\w+)"', (source / 'pyburn.py').read_text().replace(' ', ''))}
        self.assertEqual((names['tauu'], names['tauv'], names['uas'], names['vas']), (180, 181, 165, 166))
        for name in er.OUTPUT:
            self.assertTrue(name in er.DERIVED or names[name] in er.PLASIM_CODES, name)
        # The codes are those PlaSim's output routine writes on its raw output unit, 40.
        text, codes = (source / 'plasim' / 'src' / 'outmod.f90').read_text(), set()
        for m in re.finditer(r'call\s+write(?:gp|sp)\(', text):
            depth, i, args, cur = 1, m.end(), [], ''
            while depth:
                c = text[i]
                depth += (c == '(') - (c == ')')
                if depth == 1 and c == ',':
                    args.append(cur.strip())
                    cur = ''
                elif depth:
                    cur += c
                i += 1
            args.append(cur.strip())
            if args[0] == '40' and args[2].isdigit():
                codes.add(int(args[2]))
        self.assertEqual(codes, set(er.PLASIM_CODES))

    def test_new_runs_take_the_301_k_ledger_in_the_clear_share(self):
        import inspect
        for f in (er.run, er.configuration):
            params = inspect.signature(f).parameters
            self.assertEqual((params['twilight_ledger'].default, params['twilight_clear_share'].default), ('301k', True))
            self.assertEqual(params['longwave_fit'].default, '301k')
        self.assertEqual(er.TWILIGHT_PRODUCT, er.TWILIGHT_LEDGERS['301k'])
        self.assertEqual(er.LONGWAVE_PRODUCT, er.LONGWAVE_FITS['301k'])
        from climate.gcm import plasim_longwave as pl
        self.assertEqual(er.LONGWAVE_FITS['301k'], pl.REFIT)
        self.assertEqual(er.LONGWAVE_FITS['295k'], pl.CORRECTION)
        self.assertEqual(er.SLANT_PRODUCT.name, 'twilight.json')               # the calibration's slant factors

    def test_the_twilight_table_reads_as_the_patch_expects(self):
        ledgers = [path for path in er.TWILIGHT_LEDGERS.values() if path.exists()]
        if not ledgers:
            self.skipTest('no twilight ledger built')
        for ledger in ledgers:
            with self.subTest(ledger=ledger.name):
                self._twilight_table(ledger)

    def _twilight_table(self, ledger):
        import json
        product = json.loads(ledger.read_text())
        term, cone, lw = product['gcm_term'], product['gcm_term_cone'], product['longwave']
        tables = {}
        for mode in ('curved', 'full'):
            with tempfile.TemporaryDirectory() as d:
                rows = er.twilight_table(Path(d) / 'twilight.dat', mode, ledger).read_text().strip().split('\n')
            n = int(rows[0])
            self.assertEqual((n, len(rows)), (len(term['elevation_deg']), n + 3))
            table = np.array([[float(x) for x in r.split()] for r in rows[1:n + 1]])
            self.assertEqual(table.shape, (n, 2 + len(term['layers'])))      # elevation, ground, each layer
            self.assertTrue(np.all(np.diff(table[:, 0]) > 0) and n <= 400)  # rising, inside the patch's table
            np.testing.assert_allclose(table[:, 1], term['ground'], rtol=1e-6, atol=1e-12)
            escape, temperature = (np.array(r.split(), float) for r in rows[n + 1:n + 3])
            self.assertEqual(escape.size, len(term['layers']) + 1)          # each layer, then the ground
            self.assertTrue(np.all(temperature > 100.0))
            tables[mode] = table, escape
        # The full table adds the shells' widening to both the sunlight and the infrared (the product keeps six digits).
        np.testing.assert_allclose(tables['full'][0][:, 2:] - tables['curved'][0][:, 2:], np.array(cone['layers']).T,
                                   rtol=1e-4, atol=2e-9)
        np.testing.assert_allclose(tables['full'][1][:-1] - tables['curved'][1][:-1], lw['gcm_layers']['cone'],
                                   rtol=1e-5, atol=1e-8)

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
