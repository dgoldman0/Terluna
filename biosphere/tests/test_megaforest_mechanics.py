"""Single-tree beam, anchorage and wind-envelope tests."""
import json
import math
from dataclasses import replace
import unittest
import numpy as np
from scipy.linalg import solve, eigh
from atmosphere.lower_air import LowerAir
from biosphere.megaforest_mechanics import StaticTree, TreeConfig, RootConfig, root_capacity, beam_matrices, distributed_beam_load


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


if __name__ == "__main__":
    unittest.main()
