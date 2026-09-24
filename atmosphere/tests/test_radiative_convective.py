"""Radiative-convective column: thermodynamics, radiative transfer and spectroscopy checks.

Tests that need the external spectroscopic inputs skip, with that reason, when
atmosphere/radiative_convective/fetch_inputs.py has not restored them.
"""
from __future__ import annotations
import math
import unittest
import numpy as np
from scipy.special import expn, voigt_profile

from shared.constants import STEFAN_BOLTZMANN, SOLAR_CONSTANT
from atmosphere.radiative_convective import thermodynamics as th, longwave as lw, shortwave as sw
from atmosphere.radiative_convective import fetch_inputs
from atmosphere.radiative_convective.spectroscopy import C2


def inputs_present():
    try:
        return fetch_inputs.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


NEEDS_INPUTS = unittest.skipUnless(inputs_present(), 'spectroscopic inputs not restored (fetch_inputs --download)')


class ThermodynamicsTests(unittest.TestCase):
    def test_saturation_reference_points(self):
        self.assertAlmostEqual(float(th.saturation_pressure(373.124)) / 101325.0, 1.0, delta=2e-5)
        self.assertAlmostEqual(float(th.saturation_pressure(300.0)), 3536.8, delta=0.5)      # IAPWS-95
        self.assertAlmostEqual(float(th.saturation_pressure(253.15)), 103.3, delta=0.3)      # over ice
        self.assertAlmostEqual(float(th.saturation_pressure(273.16)), 611.657, delta=0.01)
        below = float(th.saturation_pressure(273.16 - 1e-9))
        self.assertAlmostEqual(below, 611.657, delta=0.01)

    def test_latent_heats(self):
        self.assertAlmostEqual(float(th.latent_heat(273.16)) / 2.5008e6, 1.0, delta=2e-3)
        self.assertAlmostEqual(float(th.latent_heat(273.15 - 1e-6)) / 2.834e6, 1.0, delta=3e-3)

    def test_dilute_adiabat_matches_textbook_form(self):
        air = th.earthlike_air(101325.0, 400.0)
        t, p = 260.0, 60000.0
        es, lv = float(th.saturation_pressure(t)), float(th.latent_heat(t))
        r = air.gas_constant / th.R_VAPOUR * es / (p - es)
        dilute = air.gas_constant / air.cp * (1 + lv * r / (air.gas_constant * t)) / (
            1 + lv * lv * r / (air.cp * th.R_VAPOUR * t * t))
        self.assertAlmostEqual(th.adiabat_slope(t, p, air) / dilute, 1.0, delta=0.01)

    def test_pure_vapour_limit_is_clausius_clapeyron(self):
        air = th.earthlike_air(101325.0, 400.0)
        t = 350.0
        p = float(th.saturation_pressure(t)) * (1 + 1e-9)
        self.assertAlmostEqual(th.adiabat_slope(t, p, air), th.R_VAPOUR * t / float(th.latent_heat(t)), delta=1e-6)

    def test_column_mass_and_gravity(self):
        air = th.earthlike_air(121590.0, 400.0)
        p = np.geomspace(121590.0, 1.0, 200)
        t = np.full_like(p, 250.0)
        x = np.zeros_like(p)
        flat = th.column_from_levels(th.MOON, air, p, t, x, constant_gravity=True)
        self.assertAlmostEqual(flat['layer_mass_kg_m2'].sum(), (p[0] - p[-1]) / th.MOON.surface_gravity, delta=1e-6)
        curved = th.column_from_levels(th.MOON, air, p, t, x)
        self.assertGreater(curved['layer_mass_kg_m2'].sum(), flat['layer_mass_kg_m2'].sum())
        scale = th.GAS_CONSTANT / air.molar_mass * 250.0 / th.MOON.surface_gravity
        z_e = float(np.interp(1.0, np.log(p[0] / p), flat['z_m']))   # one scale height up
        self.assertAlmostEqual(z_e / scale, 1.0, delta=1e-3)

    def test_built_column_structure(self):
        spec = th.ColumnSpec(th.MOON, th.earthlike_air(121590.0, 400.0), 288.0)
        col = th.build_column(spec)
        self.assertTrue(np.all(np.diff(col['p_pa']) < 0))
        self.assertAlmostEqual(col['t_k'][0], 288.0)
        above = col['p_pa'] <= col['tropopause_pa'] * (1 + 1e-12)
        self.assertTrue(np.allclose(col['t_k'][above], 200.0))
        self.assertTrue(np.allclose(col['x_h2o'][above], col['tropopause_h2o']))
        self.assertTrue(np.all(np.diff(col['t_k'][~above]) <= 1e-9))


class LongwaveTests(unittest.TestCase):
    def test_planck_integral(self):
        nu = np.linspace(0.5, 12000.0, 400000)
        total = lw.integrate(nu, math.pi * lw.planck(nu, 288.0))
        self.assertAlmostEqual(total / (STEFAN_BOLTZMANN * 288.0**4), 1.0, delta=1e-4)

    def test_isothermal_equilibrium(self):
        nu = np.linspace(100.0, 2500.0, 2000)
        tau = np.abs(np.random.default_rng(0).normal(size=(12, nu.size))) * 3
        up, down = lw.fluxes(nu, tau, np.full(13, 250.0), 250.0)
        self.assertTrue(np.allclose(up, math.pi * lw.planck(nu, 250.0)[None, :], rtol=1e-12))

    def test_grey_slab_matches_exponential_integral(self):
        nu = np.array([500.0, 800.0])
        for tau0 in (0.05, 0.7, 3.0):
            up, down = lw.fluxes(nu, np.full((1, 2), tau0), np.array([220.0, 220.0]), 290.0, n_angles=12)
            t3 = 2 * expn(3, tau0)
            expected = math.pi * (lw.planck(nu, 290.0) * t3 + lw.planck(nu, 220.0) * (1 - t3))
            # Twelve Gauss angles reproduce 2 E3(tau) to about 1e-5.
            self.assertTrue(np.allclose(up[-1], expected, rtol=2e-5))
            scale = math.pi * lw.planck(nu, 220.0)
            self.assertTrue(np.allclose(down[0], scale * (1 - t3), rtol=0, atol=2e-5 * scale.max()))

    def test_linear_source_matches_subdivision(self):
        nu = np.array([700.0])
        t_levels = np.array([290.0, 230.0])
        up, down = lw.fluxes(nu, np.full((1, 1), 2.0), t_levels, 290.0)
        n = 4000
        tl = np.linspace(290.0, 230.0, n + 1)
        # Linear in optical depth means linear in Planck radiance, not temperature.
        b = np.linspace(lw.planck(nu, 290.0)[0], lw.planck(nu, 230.0)[0], n + 1)
        tl = C2 * nu[0] / np.log1p(lw._B_FACTOR * nu[0]**3 / b)
        upf, downf = lw.fluxes(nu, np.full((n, 1), 2.0 / n), tl, 290.0)
        self.assertAlmostEqual(up[-1, 0] / upf[-1, 0], 1.0, delta=1e-6)
        self.assertAlmostEqual(down[0, 0] / downf[0, 0], 1.0, delta=1e-6)


class ShortwaveTests(unittest.TestCase):
    rng = np.random.default_rng(3)

    def test_conservative_and_absorbing_limits(self):
        tau = 10 ** self.rng.uniform(-3, 2, 500)
        g = self.rng.uniform(0, 0.9, 500)
        mu = self.rng.uniform(0.05, 1, 500)
        rd, td, rr, tr, e = sw.layer_properties(tau, np.ones_like(tau), g, mu)
        self.assertLess(np.max(abs(rd + td - 1)), 1e-6)
        self.assertLess(np.max(abs(rr + tr + e - 1)), 1e-6)
        rd, td, rr, tr, e = sw.layer_properties(tau, np.zeros_like(tau), g, mu)
        self.assertEqual(rd.max(), 0.0); self.assertEqual(rr.max(), 0.0); self.assertEqual(tr.max(), 0.0)
        self.assertTrue(np.allclose(td, np.exp(-2 * tau)))

    def test_layer_equals_two_half_layers(self):
        tau = 10 ** self.rng.uniform(-3, 2, 500)
        om, g, mu = self.rng.uniform(0, 1, 500), self.rng.uniform(0, 0.9, 500), self.rng.uniform(0.05, 1, 500)
        full = sw.layer_properties(tau, om, g, mu)
        rda, tda, rra, tra, ea = sw.layer_properties(tau / 2, om, g, mu)
        d = 1 - rda * rda
        pair = (rda + tda * tda * rda / d, tda * tda / d, rra + tda * (ea * rra + tra * rda) / d,
                ea * tra + (tra + ea * rra * rda) * tda / d, ea * ea)
        for a, b in zip(full, pair):
            self.assertLess(np.max(abs(a - b)), 1e-10)

    def test_column_energy_conservation(self):
        nlay, nnu = 20, 50
        tau_abs = 10 ** self.rng.uniform(-4, 0, (nlay, nnu))
        tau_sca = 10 ** self.rng.uniform(-3, 0, (nlay, nnu))
        radius = 1737400.0 + np.linspace(300e3, 0.0, nlay + 1)
        for mu0 in (0.05, 0.4, 1.0):
            direct, diffuse, up = sw.column_fluxes(tau_abs, tau_sca, np.zeros_like(tau_sca), 0.2, mu0, radius)
            net = direct + diffuse - up
            absorbed_layers = net[:-1] - net[1:]
            surface = net[-1]
            self.assertTrue(np.all(absorbed_layers > -1e-12))
            self.assertTrue(np.allclose(absorbed_layers.sum(0) + surface + up[0], mu0, rtol=1e-12))
            self.assertTrue(np.allclose(surface, 0.8 * (direct[-1] + diffuse[-1]), rtol=1e-12))

    def test_transparent_planet(self):
        nu = np.linspace(2000.0, 30000.0, 200)
        spectrum = np.full(nu.size, 0.05)
        zero = np.zeros((3, nu.size))
        res = sw.planetary_mean(nu, spectrum, zero, zero, zero, 0.3)
        self.assertAlmostEqual(res['reflected'] / res['incident'], 0.3, delta=1e-9)
        self.assertAlmostEqual(res['net'][-1] / res['incident'], 0.7, delta=1e-9)

    def test_slant_factors(self):
        r = 1737400.0 + np.linspace(200e3, 0.0, 11)
        self.assertTrue(np.allclose(sw.slant_factors(r, 1.0)[-1], 1.0))
        flat = 1e12 + np.linspace(200e3, 0.0, 11)
        self.assertTrue(np.allclose(sw.slant_factors(flat, 0.3)[-1], 1 / 0.3, rtol=1e-5))
        # In spherical shells the path through upper layers is shorter than 1/mu0.
        self.assertTrue(np.all(sw.slant_factors(r, 0.3)[-1] < 1 / 0.3))


@NEEDS_INPUTS
class SpectroscopyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from atmosphere.radiative_convective import spectroscopy as sp
        cls.sp = sp
        cls.h2o = sp.load_lines('H2O', 1e-28)
        cls.co2 = sp.load_lines('CO2', 1e-28)
        cls.sums = sp.PartitionSums(np.concatenate([cls.h2o.gid, cls.co2.gid]))

    def brute(self, grid, lines, t, p, p_self, pedestal):
        sp = self.sp
        lines = lines.select(grid.start - 25, grid.stop + 25)
        s = (lines.strength * self.sums.ratio(lines.gid, t) * np.exp(-sp.C2 * lines.elower * (1 / t - 1 / 296))
             * (-np.expm1(-sp.C2 * lines.nu / t)) / (-np.expm1(-sp.C2 * lines.nu / 296)))
        pa, ps = p / 101325, p_self / 101325
        c = lines.nu + lines.delta_air * (pa - ps)
        gl = (296 / t) ** lines.n_air * (lines.gamma_air * (pa - ps) + lines.gamma_self * ps)
        gg = lines.nu * np.sqrt(1.380649e-23 * t * 6.02214076e23 * 1e3) / 299792458.0 / np.sqrt(lines.mass)
        out = np.zeros(grid.size)
        for k in range(len(lines)):
            d = grid.nu - c[k]
            m = np.abs(d) <= 25
            base = voigt_profile(25, gg[k], gl[k]) if pedestal else 0.0
            out[m] += s[k] * (voigt_profile(d[m], gg[k], gl[k]) - base)
        return out

    def test_multimesh_kernel_matches_direct_sum(self):
        grid = self.sp.Grid.span(1590.0, 1610.0, 0.01)
        for lines, t, p, ps, ped in ((self.h2o, 288.0, 1e5, 1500.0, True), (self.h2o, 220.0, 100.0, 0.0, True),
                                     (self.h2o, 330.0, 8e4, 4e4, True)):
            fast = self.sp.line_absorption(grid, lines, self.sums, t, p, ps, pedestal=ped)
            ref = self.brute(grid, lines, t, p, ps, ped)
            rel = np.abs(fast - ref) / np.maximum(ref, 1e-3 * ref.max())
            self.assertLess(rel.max(), 0.01)
            self.assertAlmostEqual(fast.sum() / ref.sum(), 1.0, delta=1e-3)
        grid = self.sp.Grid.span(660.0, 675.0, 0.01)
        for t, p in ((250.0, 3e4), (200.0, 10.0)):
            fast = self.sp.line_absorption(grid, self.co2, self.sums, t, p)
            ref = self.brute(grid, self.co2, t, p, 0.0, False)
            rel = np.abs(fast - ref) / np.maximum(ref, 1e-3 * ref.max())
            self.assertLess(rel.max(), 0.01)

    def test_pressure_shift_past_the_selection_margin(self):
        # A line selected just inside the 25 cm-1 margin whose pressure shift at 3 atm carries its centre
        # beyond the cutoff contributes nothing, and must not index outside the far-wing mesh.
        sp = self.sp
        grid = sp.Grid.span(1000.0, 1010.0, 0.01)
        fields = [(975.0005, 1005.0), (1e-20, 1e-20), (0.08, 0.08), (0.4, 0.4), (100.0, 100.0), (0.7, 0.7),
                  (-0.05, -0.05)]
        two = sp.LineList('H2O', *(np.array(v) for v in fields), np.array([1, 1]), np.array([18.010565] * 2))
        fast = sp.line_absorption(grid, two, self.sums, 296.0, 3 * 101325.0)
        ref = self.brute(grid, two, 296.0, 3 * 101325.0, 0.0, False)
        self.assertTrue(np.all(np.isfinite(fast)))
        self.assertAlmostEqual(float(fast.sum()) / float(ref.sum()), 1.0, delta=1e-3)

    def test_mt_ckd_matches_aer_example(self):
        # AER MT_CKD_H2O 4.3 run_example (p=1013 mb, T=300 K, h2o_frac=0.00990098),
        # self and foreign coefficients at 497, 550 and 603 cm^-1, (c) AER.
        nu = np.array([497.0, 550.0, 603.0])
        s, f = self.sp.WaterContinuum().cross_sections(nu, 300.0, 101300.0, 0.00990098)
        self.assertTrue(np.allclose(s[0], 3.05589121e-23, rtol=1e-5))
        self.assertTrue(np.allclose(f[0], 2.41129709e-23, rtol=1e-5))
        self.assertTrue(np.all(s > 0) and np.all(f > 0))

    def test_cia_first_dataset_takes_precedence(self):
        # Two O2-O2 datasets cover 20000-29837 cm^-1; only the first is used there,
        # and the second supplies the wavenumbers beyond it.
        table = self.sp.CollisionInduced('hitran_cia_O2-O2_2024.cia')
        first = next(b for b in table.bands if b['lo'] < 21000 < b['hi'] and '32' in b['refs'])
        second = next(b for b in table.bands if '54' in b['refs'])
        s = next(s for s in first['sets'] if s['t'] == 253.0)
        peak = float(s['nu'][np.argmax(np.where(s['nu'] > 20000, s['k'], 0))])
        k = table.coefficient(np.array([peak]), 253.0)[0]
        self.assertAlmostEqual(k / np.interp(peak, s['nu'], s['k']), 1.0, delta=1e-9)
        s2 = second['sets'][0]
        beyond = float(s2['nu'][np.argmax(np.where(s2['nu'] > first['hi'] + 50, s2['k'], 0))])
        self.assertGreater(table.coefficient(np.array([beyond]), s2['t'])[0], 0.0)
        self.assertTrue(any(lo <= 25000 <= hi for lo, hi in second['exclude']))

    def test_rayleigh_optical_depth_of_earth(self):
        column = 101325.0 * 6.02214076e23 / (0.0289644 * 9.80665) * 1e-4
        tau = self.sp.rayleigh_cross_section(np.array([1e4 / 0.55]), 3.6e-4)[0] * column
        self.assertAlmostEqual(tau, 0.0973, delta=0.0006)     # Bodhaine et al. (1999), Table 3

    def test_solar_spectrum_total(self):
        nu = np.linspace(2000.0, 49500.0, 200000)
        spec, longward = sw.solar_spectrum(nu)
        total = np.trapezoid(spec, nu) + longward
        self.assertAlmostEqual(total / SOLAR_CONSTANT, 1.0, delta=2e-3)
        self.assertLess(longward / SOLAR_CONSTANT, 0.01)


@NEEDS_INPUTS
class NodeTableTests(unittest.TestCase):
    def test_table_matches_direct(self):
        from atmosphere.radiative_convective import optics, spectroscopy as sp
        cfg = optics.OpticsConfig(nu_min=640.0, nu_max=700.0, processes=1)
        optics._load(cfg)
        state = optics._STATE
        for p, t in ((712.0, 231.0), (23500.0, 263.0)):
            i0, wi, j0, wj = optics._brackets(p, t)
            optics._map(optics._node, [(cfg.grid, cfg, 'CO2', i0 + a, j0 + b) for a in (0, 1) for b in (0, 1)], cfg)
            direct = sp.line_absorption(cfg.grid, state['lines']['CO2'], state['sums'], t, p)
            table = optics._tabulated(cfg, 'CO2', p, t)
            rel = np.abs(table - direct) / np.maximum(direct, 1e-3 * direct.max())
            self.assertLess(rel.max(), 0.015)
            self.assertAlmostEqual(table.sum() / direct.sum(), 1.0, delta=0.005)


if __name__ == '__main__':
    unittest.main()
