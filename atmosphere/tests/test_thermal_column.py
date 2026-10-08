"""Molecular thermal column: conduction, hydrostatics and Jeans-escape boundary tests."""
from __future__ import annotations
import math, unittest
import numpy as np
from atmosphere.thermal_column import *


class ColumnTests(unittest.TestCase):
    def test_conduction_flux_from_independent_temperature_gradient(self):
        result,profile,_=solve_column(3e-6)
        r=profile['radius_R']*R; t=profile['temperature_K']
        flux=-(r/R)**2*(ColumnConfig().conductivity_coefficient*t)*np.gradient(t,r)
        self.assertLess(np.max(abs(flux[2:-2]-profile['conductive_flux_per_surface_W_m2'][2:-2])),1e-10)
    @classmethod
    def setUpClass(cls):
        cls.cold=solve_column(0.)
        cls.warm=solve_column(3e-6)
    def test_upper_temperature_is_solved(self):
        self.assertLess(self.cold[0]['exobase_temperature_k'],200)
        self.assertGreater(self.warm[0]['exobase_temperature_k'],220)
        self.assertLess(self.warm[0]['exobase_temperature_k'],250)
    def test_knudsen_boundary(self):
        self.assertAlmostEqual(self.warm[1]['knudsen'][-1],1,places=6)
    def test_energy_boundary(self):
        self.assertLess(abs(self.warm[0]['net_energy_boundary_residual_W_m2']),1e-11)
        self.assertLess(self.warm[0]['max_scaled_BC_residual'],1e-6)
    def test_spherical_hydrostatic_integral(self):
        _,profile,_=self.warm
        from scipy.integrate import cumulative_trapezoid
        _,rs,_,_=lower_boundary(ColumnConfig())
        x=np.log(.1/profile['pressure_Pa'])
        delta=-R*rs/GM*cumulative_trapezoid(profile['temperature_K'],x,initial=0)
        inverse=1/profile['radius_R']
        self.assertLess(np.max(abs(inverse-inverse[0]-delta)),2e-7)
    def test_grid_convergence(self):
        a=self.warm[0];b=solve_column(3e-6,ColumnConfig(mesh_points=161,tolerance=2e-7))[0]
        self.assertAlmostEqual(a['exobase_temperature_k'],b['exobase_temperature_k'],delta=2e-4)
        self.assertLess(abs(a['molecular_loss_kg_s']/b['molecular_loss_kg_s']-1),1e-4)
    def test_boundary_dependence_exposed(self):
        b=solve_column(1e-6,ColumnConfig(lower_temperature_k=210))[0]
        a=solve_column(1e-6)[0]
        self.assertGreater(b['exobase_temperature_k'],a['exobase_temperature_k'])
    def test_high_heating_domain_flag(self):
        r=solve_column(1e-5)[0]
        self.assertIn('KINETIC_ESCAPE_SENSITIVITY_REQUIRED',r['domain_flags'])
    def test_jeans_energy_matches_maxwell_integral(self):
        from scipy.integrate import quad
        for lam in [1.,10.,25.]:
            number=quad(lambda y:y*math.exp(-y),lam,np.inf,epsabs=1e-20)[0]
            energy=quad(lambda y:y*(y-lam)*math.exp(-y),lam,np.inf,epsabs=1e-20)[0]
            self.assertAlmostEqual(energy/number,(lam+2)/(lam+1),places=9)
    def test_invalid_heat(self):
        with self.assertRaises(ValueError):solve_column(-1)
    def test_molecular_species_flux(self):
        r=self.warm[0]
        self.assertGreater(r['N2_loss_kg_s'],r['O2_loss_kg_s'])
        self.assertAlmostEqual(r['molecular_loss_kg_s'],r['N2_loss_kg_s']+r['O2_loss_kg_s'])


class TracedShapeTests(unittest.TestCase):
    """A heating shape given as the share of the heat below each log pressure above the base."""
    cfg = ColumnConfig(lower_pressure_pa=0.3, lower_temperature_k=160.)

    def traced(self, x, f):
        import dataclasses
        return dataclasses.replace(self.cfg, heating_shape='traced', heating_log_pressure=tuple(x),
                                   heating_fraction=tuple(f))

    def test_the_middle_shape_traced_reproduces_it(self):
        middle = solve_column(1e-5, self.cfg)[0]
        length = math.log(.3 / middle['exobase_pressure_pa'])
        x = np.linspace(0, length, 400)
        shaped = solve_column(1e-5, self.traced(x, heating_cdf(x / length, 'middle')))[0]
        self.assertAlmostEqual(shaped['exobase_temperature_k'], middle['exobase_temperature_k'], delta=.01)
        self.assertLess(abs(shaped['molecular_loss_kg_s'] / middle['molecular_loss_kg_s'] - 1), 1e-3)

    def test_heat_laid_at_the_base_barely_warms_the_exobase(self):
        cold = solve_column(0., self.cfg)[0]['exobase_temperature_k']
        based = solve_column(1e-5, self.traced([0., .5, 30.], [0., 1., 1.]))[0]['exobase_temperature_k']
        middle = solve_column(1e-5, self.cfg)[0]['exobase_temperature_k']
        self.assertLess(based - cold, .1 * (middle - cold))

    def test_a_traced_shape_must_be_a_cumulative_share_from_the_base(self):
        for x, f in (([0., 1.], [0., .5]), ([.5, 1.], [0., 1.]), ([0., 1., 2.], [0., .7, .6]), ([0.], [0.])):
            with self.assertRaises(ValueError):
                solve_column(1e-6, self.traced(x, f))


class InfraredCoolingTests(unittest.TestCase):
    """CO2, O and NO radiating in excess of what air at the base temperature emits."""
    base = ColumnConfig(lower_pressure_pa=0.3, lower_temperature_k=170.)
    columns = tuple(float(v) for v in np.geomspace(1e10, 1e24, 57))

    def cooled(self, co2=4e-4, o=0., no=0., escape=None, **kwargs):
        import dataclasses
        x = np.arange(0., 40.001, .5)
        escape = escape or tuple(float(v) for v in 1 / (1 + np.asarray(self.columns) / 1e15))
        return dataclasses.replace(self.base, coolant_log_pressure=tuple(x), co2_mole_fraction=(co2,) * x.size,
                                   o_mole_fraction=(o,) * x.size, no_mole_fraction=(no,) * x.size,
                                   co2_escape_column_cm2=self.columns, co2_escape=escape, **kwargs)

    def test_without_coolants_the_column_is_unchanged(self):
        plain = solve_column(1e-5, self.base)[0]
        empty = solve_column(1e-5, self.cooled(co2=0.))[0]
        self.assertAlmostEqual(empty['exobase_temperature_k'], plain['exobase_temperature_k'], delta=.01)
        self.assertLess(abs(empty['molecular_loss_kg_s'] / plain['molecular_loss_kg_s'] - 1), 1e-3)
        self.assertNotIn('infrared_cooling_W_m2', plain)
        self.assertIn('NO_IR_COOLING', plain['domain_flags'])

    def test_air_left_at_the_base_temperature_neither_cools_nor_warms(self):
        summary, profile, _ = solve_column(0., self.cooled())
        self.assertLess(np.ptp(profile['temperature_K']), 1e-6)
        self.assertLess(abs(summary['infrared_cooling_W_m2']), 1e-12)

    def test_the_heat_leaves_by_radiation_conduction_and_escape(self):
        q = 1e-5
        summary, profile, _ = solve_column(q, self.cooled())
        j = np.array([summary['N2_loss_kg_s'], summary['O2_loss_kg_s']]) / AREA
        advected = j @ (3.5 * RS_SPECIES * summary['lower_temperature_k'] - GM / (summary['lower_radius_R'] * R))
        into_base = -(summary['lower_conductive_flux_W_m2'] + advected)
        out = summary['infrared_cooling_W_m2'] + into_base + summary['upper_energy_to_infinity_W_m2']
        self.assertLess(abs(out / q - 1), 1e-6)
        self.assertGreater(summary['co2_cooling_W_m2'], .9 * q)          # 400 ppm of CO2 takes nearly all of it
        self.assertAlmostEqual(summary['infrared_cooling_W_m2'], summary['co2_cooling_W_m2']
                               + summary['o_cooling_W_m2'] + summary['no_cooling_W_m2'], delta=2e-3 * q)
        self.assertIn('IR_COOLING_PRESCRIBED_COMPOSITION', summary['domain_flags'])

    def test_cooling_keeps_the_column_more_compact(self):
        for q in (3e-6, 1e-5):
            plain = solve_column(q, self.base)[0]
            cooled = solve_column(q, self.cooled())[0]
            self.assertLess(cooled['exobase_radius_R'], plain['exobase_radius_R'])
            self.assertLess(cooled['molecular_loss_kg_s'], plain['molecular_loss_kg_s'])

    def test_thin_co2_follows_the_middle_atmospheres_two_level_rates(self):
        from atmosphere.radiative_convective import ck
        cfg = self.cooled(o=1e-2, escape=(1.0,) * len(self.columns))
        t, p = np.array([150., 250., 400.]), np.array([1e-3, 1e-2, 1e-1])
        got = infrared_cooling(p, t, np.log(.3 / p), np.full(3, 1.2 * R), cfg, lower_boundary(cfg)[3], parts=True)
        molecules = 1 - 4e-4 - 1e-2           # the N2 and O2 that deactivate CO2 beside the oxygen atoms
        layer = dict(layer_t_k=t, layer_p_pa=p, layer_x_h2o=np.zeros(3),
                     dry_fractions=dict(N2=.8246 * molecules, O2=.1754 * molecules))
        nonlte = ck.NonLTE()
        self.assertEqual(nonlte.einstein_a, CO2_BAND['einstein_a'])
        eps = nonlte.epsilon(layer, o_mixing=np.full(3, 1e-2))
        n = p / (KB * t)
        theta = CO2_BAND['theta_k']
        excess = 1 - np.expm1(theta / t) / np.expm1(theta / 170.)
        thin = 4e-4 * n * KB * theta * 1.5 * 2 * np.exp(-theta / t) * excess * eps / (1 + eps)
        np.testing.assert_allclose(got['co2'], thin, rtol=1e-3)

    def test_oxygen_lines_cool_in_equilibrium_when_thin(self):
        # So little oxygen that its lines are thin both ways; with more, the column below the level absorbs the
        # downward half.
        t, p = np.array([250., 400.]), np.array([1e-6, 1e-6])
        cfg = self.cooled(co2=0., o=1e-7)
        got = infrared_cooling(p, t, np.log(.3 / p), np.full(2, 3 * R), cfg, lower_boundary(cfg)[3], parts=True)['o']
        thick = infrared_cooling(p, t, np.log(.3 / p), np.full(2, 3 * R), self.cooled(co2=0., o=1e-3),
                                 lower_boundary(cfg)[3], parts=True)['o']
        np.testing.assert_allclose(thick / got * 1e-7 / 1e-3, .5, rtol=.05)
        z = 5 + 3 * np.exp(-HC_K * 158.265 / t) + np.exp(-HC_K * 226.977 / t)
        lines = (3 * np.exp(-HC_K * 158.265 / t) / z * 8.91e-5 * KB * HC_K * 158.265
                 * (1 - np.expm1(HC_K * 158.265 / t) / np.expm1(HC_K * 158.265 / 170.))
                 + np.exp(-HC_K * 226.977 / t) / z * 1.75e-5 * KB * HC_K * 68.712
                 * (1 - np.expm1(HC_K * 68.712 / t) / np.expm1(HC_K * 68.712 / 170.)))
        np.testing.assert_allclose(got, 1e-7 * p / (KB * t) * lines, rtol=5e-3)

    def test_doppler_escape_matches_direct_quadrature(self):
        from scipy.integrate import quad
        from scipy.special import expn
        for tau in (0., .1, 1., 10., 300.):
            direct = quad(lambda x: math.exp(-x * x) / math.sqrt(math.pi) * expn(2, tau * math.exp(-x * x)),
                          -8, 8, points=[0.], limit=200)[0]
            self.assertAlmostEqual(float(doppler_escape(tau)), direct, delta=2e-3 * max(direct, 1e-3))

    def test_cooling_inputs_are_checked(self):
        import dataclasses
        bad = (dict(co2_mole_fraction=(4e-4,)), dict(coolant_log_pressure=(1., 2.) + (3.,) * 79),
               dict(co2_escape=(0.5,) * 56 + (0.9,)), dict(cooling_multiplier=-1.))
        for change in bad:
            with self.assertRaises(ValueError):
                solve_column(1e-6, dataclasses.replace(self.cooled(), **change))


if __name__ == "__main__":
    unittest.main()
