"""Periodic latitude-longitude energy-balance screen tests."""
from __future__ import annotations
import unittest
import numpy as np
from climate.cycle import ClimateConfig, laplacian, solve_climate


class ClimateTests(unittest.TestCase):
    def test_transport_conservative(self):
        x,L=laplacian(6,24)
        self.assertLess(np.max(abs(L.sum(axis=0))),1e-12)
        self.assertLess(np.max(abs(L-L.T)),1e-12)
        self.assertLess(np.linalg.eigvalsh(L).max(),1e-10)
    def test_global_mean_from_analytic_energy_balance(self):
        c=ClimateConfig();r,_=solve_climate(c)
        expected=273.15+(c.solar_w_m2*c.transmission*(1-c.albedo)/4-c.olr_at_273_w_m2)/c.olr_slope_w_m2_k
        self.assertAlmostEqual(r['global_mean_K'],expected,places=9)
        self.assertLess(r['max_discrete_energy_residual_W_m2'],1e-7)
        self.assertLess(r['periodic_residual_K'],1e-8)
    def test_timestep_convergence(self):
        a,_=solve_climate(ClimateConfig(steps=360));b,_=solve_climate(ClimateConfig(steps=720))
        self.assertAlmostEqual(a['minimum_K'],b['minimum_K'],delta=.01)
    def test_olr_is_an_assumption(self):
        a,_=solve_climate();b,_=solve_climate(ClimateConfig(olr_at_273_w_m2=220))
        self.assertAlmostEqual(a['global_mean_K']-b['global_mean_K'],10,places=8)
    def test_freezing_is_flagged(self):
        a,_=solve_climate(ClimateConfig(olr_at_273_w_m2=260))
        self.assertTrue(a['freezing_flag'])
    def test_water_damping(self):
        a,_=solve_climate(ClimateConfig(layout='dry'))
        b,_=solve_climate(ClimateConfig(layout='distributed',water_depth_m=50))
        self.assertLess(b['maximum_K']-b['minimum_K'],a['maximum_K']-a['minimum_K'])


if __name__ == "__main__":
    unittest.main()
