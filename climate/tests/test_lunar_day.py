"""Single-column lunar-day model: soil conduction, surface exchange and stability checks."""
from pathlib import Path
import tempfile
import unittest
import numpy as np

from climate import lunar_day as ld


class LunarDayTests(unittest.TestCase):
    def test_tridiagonal(self):
        rng = np.random.default_rng(4)
        n = 12
        a, c = -rng.uniform(0, 1, n), -rng.uniform(0, 1, n)
        b = 3 + rng.uniform(0, 1, n)
        a[0] = c[-1] = 0.0
        m = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        d = rng.uniform(-1, 1, n)
        np.testing.assert_allclose(ld._tridiagonal(a, b, c, d), np.linalg.solve(m, d), rtol=1e-12)

    def test_stability_function(self):
        self.assertAlmostEqual(float(ld.stability(0.0, 'short')), 1.0)
        self.assertLess(float(ld.stability(0.5, 'short')), float(ld.stability(0.5, 'long')))
        self.assertGreater(float(ld.stability(-0.5, 'short')), 1.0)

    def test_soil_conduction_conserves_heat(self):
        model = ld.Model.__new__(ld.Model)
        model.case = ld.Case('stub', 'earth_control')
        edges = np.concatenate([[0.0], 0.01 * (1.35 ** np.arange(25) - 1) / 0.35])
        model.dz_soil = np.diff(edges * 10.0 / edges[-1])
        model.soil = np.linspace(290.0, 280.0, model.dz_soil.size)
        s = model.case.surface
        before = float((model.soil * model.dz_soil).sum()) * s.soil_heat_capacity_j_m3_k
        flux, dt = 50.0, 600.0
        model._conduct(flux, dt, 0.0)
        after = float((model.soil * model.dz_soil).sum()) * s.soil_heat_capacity_j_m3_k
        self.assertAlmostEqual((after - before) / (flux * dt), 1.0, delta=1e-9)

    def test_vapour_condensed_aloft_heats_the_air(self):
        class Column:
            p = np.linspace(1.0e5, 4.0e4, 7)
            bl = np.array([True, True, False, False, False, False])
            t_ref = np.linspace(290.0, 240.0, 6)
            mass = -np.diff(p) / 9.81
        c = Column()
        c.x_ref = 0.5 * ld.saturation_mole_fraction(c.t_ref, np.sqrt(c.p[1:] * c.p[:-1]))
        above = ~c.bl
        case = ld.Case('stub', 'earth_control', relaxation_days=1e12)
        half_life = case.vapour_relaxation_days * ld.DAY_S * np.log(2.0)
        for factor, dt, expected in ((2.0, 1.0, 1.0), (0.5, half_life, 0.75)):
            model = ld.Model.__new__(ld.Model)
            model.case, model.c = case, c
            model.cp_air, model.condensed = 1004.0, 0.0
            model.t, model.x = c.t_ref.copy(), factor * c.x_ref
            t0, x0 = model.t.copy(), model.x.copy()
            model.relax(dt)
            # Vapour beyond the reference humidity condenses at once; drier air regains vapour only at the
            # exchange rate (half the deficit in a half-life).
            np.testing.assert_allclose(model.x[above], expected * c.x_ref[above], rtol=1e-9)
            rained = 0.622 * np.maximum(x0 - c.x_ref, 0.0)[above]
            # Condensed vapour releases its latent heat where it condenses; arriving vapour brings none.
            np.testing.assert_allclose(model.cp_air * (model.t - t0)[above],
                                       ld.th.latent_heat(t0[above]) * rained, rtol=1e-9, atol=1e-12)
            self.assertAlmostEqual(model.condensed, float((rained * c.mass[above]).sum()), delta=1e-9)

    def test_checkpoint_restores_the_saved_day(self):
        class Column:
            t_ref = np.linspace(290.0, 200.0, 5)
            x_ref = np.full(5, 1e-3)
            table_layers = np.ones((3, 5))
            table_ground = np.ones(3)

        def model(case):
            m = ld.Model.__new__(ld.Model)
            m.case, m.c = case, Column()
            m.t, m.x, m.soil = Column.t_ref + 1.5, Column.x_ref * 2.0, np.linspace(295.0, 288.0, 25)
            m.ocean_heat, m.day_means = -12.5, [290.0, 290.4]
            return m
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'stub.npz'
            saved = model(ld.Case('stub', 'earth_control'))
            saved._save(path, 2)
            fresh = model(ld.Case('stub', 'earth_control'))
            fresh.t, fresh.x, fresh.soil, fresh.ocean_heat, fresh.day_means = 0 * fresh.t, 0 * fresh.x, 0 * fresh.soil, 0.0, []
            self.assertEqual(fresh._resume(path), 2)
            for name in ('t', 'x', 'soil'):
                np.testing.assert_array_equal(getattr(fresh, name), getattr(saved, name))
            self.assertEqual((fresh.ocean_heat, fresh.day_means), (saved.ocean_heat, saved.day_means))
            # A checkpoint from another case (or other code or column) is ignored.
            self.assertEqual(model(ld.Case('stub', 'earth_control', days=9))._resume(path), 0)
            self.assertEqual(model(ld.Case('stub', 'earth_control'))._resume(Path(tmp) / 'missing.npz'), 0)

    def test_spin_up_moves_ground_to_periodic_state(self):
        class Column:
            day_s = 29.5 * 86400.0
            ts_ref = 288.0
        land = ld.Model.__new__(ld.Model)
        land.case, land.c = ld.Case('stub', 'earth_control'), Column()
        land.soil = np.linspace(295.0, 285.0, 25)
        mean = np.linspace(292.0, 285.0, 25)
        wave = land.soil - mean
        land.spin_up(land.soil.copy(), mean)
        # Every layer's day-mean moves to the skin's; the daily wave about it is kept.
        np.testing.assert_allclose(land.soil - wave, np.full(25, 292.0))

        sea = ld.Model.__new__(ld.Model)
        sea.case, sea.c = ld.Case('stub', 'earth_control', surface=ld.SEA), Column()
        sea.ocean_heat = 0.0
        start = np.array([288.0])
        sea.soil = np.array([289.0])
        sea.spin_up(start, np.array([288.6]))
        cap = ld.SEA.mixed_layer_m * 1025.0 * 3990.0
        self.assertAlmostEqual(sea.ocean_heat, -cap * 1.0 / Column.day_s)
        self.assertAlmostEqual(float(sea.soil[0]), 289.0 + 288.0 - 288.6)
