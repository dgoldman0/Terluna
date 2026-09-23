"""Analytic, conservation, regression and invalid-input tests for wind screening."""
import json
import math
from dataclasses import replace
import unittest
import numpy as np
from scipy.linalg import solve, eigh
from atmosphere.lower_air import LowerAir, AirConfig, A1Profile, R_MOON, GM_MOON
from climate.canopy_flow import CanopyFlow, CanopyConfig
from biosphere.megaforest_mechanics import (StaticTree, TreeConfig, RootConfig,
    root_capacity, beam_matrices, distributed_beam_load)


class AirTests(unittest.TestCase):
    def test_surface_and_spherical_isothermal(self):
        a = LowerAir()
        z = np.array([0., 100., 500., 5000.])
        p = a.sample(z)
        expected = a.config.pressure_pa*np.exp(-GM_MOON*(1/R_MOON-1/(R_MOON+z))/a.gas_constant/288)
        np.testing.assert_allclose(p['pressure_pa'], expected, rtol=2e-14)
        np.testing.assert_allclose(p['density_kg_m3'], p['pressure_pa']/a.gas_constant/288)
        self.assertAlmostEqual(p['pressure_pa'][0], 121590)

    def test_hydrostatic_derivative(self):
        for lapse in (0., .0015, -.002):
            a = LowerAir(AirConfig(lapse_k_per_m=lapse))
            z = np.array([20., 250., 1000., 4000.])
            deriv = (a.sample(z+.02)['pressure_pa']-a.sample(z-.02)['pressure_pa'])/.04
            s = a.sample(z)
            np.testing.assert_allclose(deriv, -s['density_kg_m3']*s['gravity_m_s2'], rtol=2e-8)

    def test_humidity_and_composition(self):
        a = LowerAir().sample(0)['density_kg_m3']
        self.assertLess(LowerAir(AirConfig(vapor_kg_per_kg_dry=.01)).sample(0)['density_kg_m3'], a)
        self.assertGreater(LowerAir(AirConfig(oxygen_mole_fraction=.21)).sample(0)['density_kg_m3'], a)

    def test_air_rejects_invalid_values(self):
        for kwargs in ({'pressure_pa':0}, {'oxygen_mole_fraction':1.1}, {'vapor_kg_per_kg_dry':-.1},
                       {'temperature_k':float('nan')}, {'lapse_k_per_m':.1}):
            with self.assertRaises(ValueError): LowerAir(AirConfig(**kwargs))
        for z in (-1, 5001, float('nan')):
            with self.assertRaises(ValueError): LowerAir().sample(z)

    @staticmethod
    def state():
        return dict(schema='open-moon-atmospheric-profile/1', world='moon', regime='test_fixture',
            planet=dict(radius=1737400.,g0=1.62), upper=dict(top_m=1000.),
            rows=[dict(z=0.,p=121590.,T=294.,r=.008), dict(z=1000.,p=119300.,T=292.,r=.006)])

    def test_a1_contract(self):
        state = self.state()
        a = A1Profile(state)
        expected_p = math.sqrt(121590*119300)
        expected_tv = 293*(1+.007/(287.05/461.5))/(1+.007)
        sample = a.sample(500.)
        self.assertAlmostEqual(sample['pressure_pa'], expected_p, places=8)
        self.assertAlmostEqual(sample['density_kg_m3'], expected_p/(287.05*expected_tv), places=12)
        state['rows'][0]['p'] = 1  # Adapter keeps its own immutable numerical state.
        self.assertGreater(a.sample(0)['pressure_pa'], 100000)
        with self.assertRaises(ValueError): a.sample(1001)
        with self.assertRaises(FileNotFoundError): A1Profile.from_file('/no/such/a1/export.json')

    def test_a1_invalid_schema_and_domain(self):
        for change in ('schema', 'pressure', 'z', 'boundary'):
            s = self.state()
            if change == 'schema': s['schema'] = 'unknown'
            if change == 'pressure': s['rows'][1]['p'] = 130000
            if change == 'z': s['rows'][1]['z'] = 0
            if change == 'boundary': s['upper']['top_m'] = 2000
            with self.assertRaises(ValueError): A1Profile(s)


class CanopyTests(unittest.TestCase):
    def test_zero_drag_constant_density_is_linear(self):
        class ConstantAir:
            def sample(self, z): return {'density_kg_m3':np.ones_like(np.asarray(z, float))}
        flow = CanopyFlow(ConstantAir(), CanopyConfig(leaf_area_index=0))
        z = np.linspace(0,300,33)
        np.testing.assert_allclose(flow.ratio(z), z/300, atol=1e-12)

    def test_momentum_balance_and_boundary(self):
        flow = CanopyFlow(LowerAir())
        self.assertLess(flow.diagnostics['momentum_relative_residual'], 1e-6)
        np.testing.assert_allclose(flow.ratio([0,300]), [0,1], atol=1e-10)
        self.assertTrue(np.all(np.diff(flow.ratio(np.linspace(0,600,300))) >= 0))
        self.assertAlmostEqual(flow.ratio(300-1e-5), flow.ratio(300+1e-5), places=5)

    def test_drag_mixing_and_scaling(self):
        a = LowerAir()
        low = CanopyFlow(a, CanopyConfig(leaf_area_index=5))
        high = CanopyFlow(a, CanopyConfig(leaf_area_index=22))
        mixed = CanopyFlow(a, CanopyConfig(leaf_area_index=22, mixing_k_over_uh=.12))
        self.assertLess(high.ratio(150), low.ratio(150))
        self.assertGreater(mixed.ratio(150), high.ratio(150))
        np.testing.assert_allclose(high.sample([50,150,300],20), 2*high.sample([50,150,300],10))
        np.testing.assert_allclose(high.sample([50,150,300],0), 0)

    def test_canopy_invalid_inputs(self):
        for cfg in (CanopyConfig(leaf_area_index=-1),CanopyConfig(mixing_k_over_uh=0),
                    CanopyConfig(above_displacement_fraction=.95),CanopyConfig(height_m=6000)):
            with self.assertRaises(ValueError): CanopyFlow(LowerAir(),cfg)
        with self.assertRaises(ValueError): CanopyFlow(LowerAir()).sample(30,-1)


class BeamTests(unittest.TestCase):
    def test_cantilever_tip_load(self):
        length, rigidity, force = 20., 1e8, 1000.
        z = np.linspace(0,length,21)
        k,_ = beam_matrices(z,np.full(20,rigidity),np.zeros(20))
        load = np.zeros(len(k)); load[-2]=force
        u = solve(k[2:,2:],load[2:],assume_a='pos')
        self.assertAlmostEqual(u[-2]/(force*length**3/(3*rigidity)),1,places=8)
        self.assertAlmostEqual(u[-1]/(force*length**2/(2*rigidity)),1,places=8)

    def test_uniform_distributed_load(self):
        length, rigidity, q = 30., 2e8, 300.
        z = np.linspace(0,length,31)
        k,_ = beam_matrices(z,np.full(30,rigidity),np.zeros(30))
        load = distributed_beam_load(z,np.full(30,q))
        u = np.r_[0.,0.,solve(k[2:,2:],load[2:],assume_a='pos')]
        reactions = k@u-load
        self.assertAlmostEqual(u[-2]/(q*length**4/(8*rigidity)),1,places=7)
        self.assertAlmostEqual(abs(reactions[0])/(q*length),1,places=8)
        self.assertAlmostEqual(abs(reactions[1])/(q*length**2/2),1,places=8)

    def test_euler_buckling(self):
        length, rigidity, force = 20., 1e8, 10000.
        z=np.linspace(0,length,31)
        k,g=beam_matrices(z,np.full(30,rigidity),np.full(30,force))
        last=len(k)-3
        beta=eigh(g[2:,2:],k[2:,2:],subset_by_index=[last,last],eigvals_only=True)[0]
        expected=math.pi**2*rigidity/(4*length**2)
        self.assertLess(abs((force/beta)/expected-1), 1e-6)

    def test_load_shape_rejection(self):
        with self.assertRaises(ValueError): distributed_beam_load([0,1,2],[1])
        with self.assertRaises(ValueError): distributed_beam_load([0,1],[float('nan')])

    def test_reject_invalid_beam(self):
        with self.assertRaises(ValueError): beam_matrices([0,1],[0],[1])
        with self.assertRaises(ValueError): beam_matrices([1,0],[1],[1])


class TreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.rigid=StaticTree(LowerAir())

    def test_zero_wind_and_quadratic_load(self):
        zero=self.rigid.evaluate(0)
        self.assertEqual(zero['base_total_moment_nm'],0)
        self.assertEqual(zero['tip_displacement_m'],0)
        a,b=self.rigid.evaluate(10),self.rigid.evaluate(20)
        for key in ('base_wind_moment_nm','base_total_moment_nm','tip_displacement_m','total_drag_n'):
            self.assertAlmostEqual(b[key]/a[key],4,places=9)

    def test_reconfiguration_reduces_crown_drag(self):
        flex=StaticTree(LowerAir(),replace(TreeConfig(),vogel_exponent=-.8))
        a,b=flex.evaluate(10),flex.evaluate(20)
        self.assertAlmostEqual(b['crown_force_n']/a['crown_force_n'],2**1.2,places=10)
        self.assertLess(a['root_utilization'],self.rigid.evaluate(10)['root_utilization'])
        # At the floor the asymptotic quadratic drag law returns.
        c,d=flex.evaluate(300),flex.evaluate(600)
        self.assertAlmostEqual(d['crown_force_n']/c['crown_force_n'],4,places=10)

    def test_thresholds_and_domain_labels(self):
        env=self.rigid.envelope()
        self.assertEqual(env['first_mode'],'root')
        self.assertTrue(env['first_threshold_within_displacement_domain'])
        self.assertFalse(env['stem_threshold_within_domain'])
        for mode in ('root','stem','crown','domain'):
            u=env[mode+'_threshold_m_s']
            if u is not None:
                self.assertAlmostEqual(self.rigid.evaluate(u)[mode+'_utilization'],1,places=7)

    def test_right_censoring(self):
        env=self.rigid.envelope(1)
        self.assertIsNone(env['first_threshold_m_s'])
        self.assertEqual(env['first_mode'],'right_censored')
        json.dumps(env,allow_nan=False)

    def test_root_serial_bottleneck(self):
        root=RootConfig()
        result=root_capacity(root,TreeConfig(),1.62)
        self.assertEqual(result['capacity_nm'],min(result['soil_capacity_nm'],result['tissue_capacity_nm']))
        weak=root_capacity(replace(root,root_tensile_area_to_stem_area=.001),TreeConfig(),1.62)
        self.assertEqual(weak['bottleneck'],'root_tissue')
        for g in (1.62,9.80665):
            r=root_capacity(root,TreeConfig(),g)
            self.assertGreater(r['capacity_nm'],0)
        self.assertAlmostEqual(root_capacity(root,TreeConfig(),9.80665)['soil_weight_moment_nm']/
                               root_capacity(root,TreeConfig(),1.62)['soil_weight_moment_nm'],9.80665/1.62)

    def test_wet_mass_changes_prestress_not_primary_drag(self):
        dry=StaticTree(LowerAir(),replace(TreeConfig(),retained_water_mm=0))
        wet=self.rigid
        self.assertGreater(wet.crown_mass,dry.crown_mass)
        self.assertGreater(wet.buckling_utilization,dry.buckling_utilization)
        self.assertAlmostEqual(wet.evaluate(10)['base_wind_moment_nm'],dry.evaluate(10)['base_wind_moment_nm'])
        self.assertGreater(wet.evaluate(10)['base_gravity_moment_nm'],dry.evaluate(10)['base_gravity_moment_nm'])

    def test_buckling_rejected(self):
        weak=StaticTree(LowerAir(),replace(TreeConfig(),young_modulus_pa=1e7))
        self.assertEqual(weak.envelope()['status'],'SELF_WEIGHT_BUCKLING')

    def test_grid_convergence(self):
        fine=StaticTree(LowerAir(),replace(TreeConfig(),elements=240))
        coarse=self.rigid.envelope(); reference=fine.envelope()
        for mode in ('root','stem','domain'):
            self.assertLess(abs(coarse[mode+'_threshold_m_s']/reference[mode+'_threshold_m_s']-1), .015)
        self.assertLess(abs(self.rigid.buckling_utilization/fine.buckling_utilization-1), .01)

    def test_invalid_tree_root_wind(self):
        for t in (replace(TreeConfig(),base_diameter_m=-1), replace(TreeConfig(),vogel_exponent=-3),
                  replace(TreeConfig(),crown_base_fraction=.999999)):
            with self.assertRaises(ValueError): StaticTree(LowerAir(),t)
        with self.assertRaises(ValueError): StaticTree(LowerAir(),roots=replace(RootConfig(),plate_depth_m=0))
        with self.assertRaises(ValueError): self.rigid.evaluate(float('nan'))


class ScopeTests(unittest.TestCase):
    def test_immersion_output_rejected(self):
        from research.run_megaforest_wind import run, ROOT
        with self.assertRaises(ValueError): run(ROOT/'immersion'/'forbidden-wind-output')
        self.assertFalse((ROOT/'immersion'/'forbidden-wind-output').exists())


if __name__ == '__main__': unittest.main()
