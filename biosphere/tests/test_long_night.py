"""Periodic carbon-reserve and aquatic-oxygen requirement tests."""
from __future__ import annotations
import unittest
import numpy as np
from biosphere.long_night import *


class CarbonTests(unittest.TestCase):
    def brute(self,a,dt):
        return max([0.]+[-sum(np.tile(a,2)[s:s+k])*dt for s in range(len(a)) for k in range(1,len(a)+1)])
    def test_circular_deficit_matches_bruteforce(self):
        rng=np.random.default_rng(7)
        for n in range(2,21):
            a=rng.normal(size=n);a+=max(0.,-.9*a.mean())
            r=periodic_storage_requirement(a,.2)
            self.assertAlmostEqual(r['worst_single_cycle_deficit'],self.brute(a,.2),places=11)
    def test_negative_balance_has_no_finite_solution(self):
        r=periodic_storage_requirement([1,-2],1)
        self.assertFalse(r['cycle_balance_feasible']);self.assertIsNone(r['minimum_capacity'])
    def test_square_wave_analytic(self):
        r=periodic_storage_requirement([2,2,-.25,-.25],1)
        self.assertAlmostEqual(r['minimum_capacity'],.5)
    def test_periodic_sufficiency_and_insufficiency(self):
        a,dt,_=carbon_trace();required=periodic_storage_requirement(a,dt)['minimum_capacity']
        success=simulate_store(a,dt,required+1e-9,5)
        fail=simulate_store(a,dt,.9*required,5)
        self.assertLess(success['unmet_demand'],1e-8)
        self.assertGreater(fail['cycles'][-1]['unmet'],.01)
        self.assertLess(abs(success['ledger_residual']),1e-8)
    def test_dim_edges_increase_required_storage(self):
        a,dt,_=carbon_trace(night_demand_fraction=.25)
        self.assertGreater(periodic_storage_requirement(a,dt)['minimum_capacity'],.25*PERIOD_DAYS/2)
    def test_carbon_grid_convergence(self):
        vals=[]
        for n in [720,1440,2880]:
            a,dt,_=carbon_trace(steps=n);vals.append(periodic_storage_requirement(a,dt)['minimum_capacity'])
        self.assertLess(abs(vals[2]-vals[1]),2e-4)
    def test_infeasible_initial_store(self):
        with self.assertRaises(ValueError):simulate_store([1,-1],1,2,initial=3)
    def test_uniform_zero_net(self):
        self.assertEqual(periodic_storage_requirement([0,0],1)['minimum_capacity'],0)


class OxygenTests(unittest.TestCase):
    def test_no_exchange_analytic(self):
        self.assertEqual(oxygen_dark_end(8,8,.5,0,10),3)
    def test_dark_solution_vs_small_steps(self):
        y=8.;dt=.0005
        for _ in range(20000):y+=dt*(-.5+.1*(8-y))
        self.assertAlmostEqual(y,oxygen_dark_end(8,8,.5,.1,10),delta=1e-4)
    def test_periodic_and_mass_balance(self):
        r=oxygen_periodic()
        self.assertLess(abs(r['periodic_residual']),1e-8)
        self.assertLess(abs(r['mass_balance_residual']),1e-8)
    def test_infeasible_demand_is_exposed(self):
        r=oxygen_periodic(respiration_g_m3_day=4,mean_day_production_g_m3_day=.1)
        self.assertFalse(r['prescribed_aerobic_demand_feasible'])
    def test_oxygen_capacity_is_enforced(self):
        r=oxygen_periodic(exchange_per_day=.01,respiration_g_m3_day=2.,mean_day_production_g_m3_day=6.)
        self.assertLessEqual(r['maximum_g_m3'],9.6)
        self.assertGreater(r['unmet_respiration_g_m3'],0)
        self.assertGreater(r['outgassing_above_cap_g_m3'],0)
    def test_convergence(self):
        a=oxygen_periodic(steps=720);b=oxygen_periodic(steps=1440)
        self.assertAlmostEqual(a['minimum_g_m3'],b['minimum_g_m3'],delta=.002)


if __name__ == "__main__":
    unittest.main()
