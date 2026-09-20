from __future__ import annotations
import math,unittest
from dataclasses import replace
import numpy as np
from atmosphere.thermal_column import *
from atmosphere.spectral_interface import Band,deposited_heat,historical_solar_bands,gap_energy_bounds,audit_optical_files
from biosphere.long_night import *
from climate.cycle import ClimateConfig,laplacian,solve_climate

class SpectralTests(unittest.TestCase):
    def test_missing_input_rejected(self):
        with self.assertRaises(ValueError):deposited_heat([], (1,120))
        b=Band(1,120,1e-3,None,.9,.1,2,'incomplete test')
        with self.assertRaises(ValueError):deposited_heat([b],(1,120))
    def test_gap_and_overlap_rejected(self):
        for start in [49,51]:
            bs=[Band(1,50,.001,.01,1,.1,2,'test'),Band(start,120,.001,.01,1,.1,2,'test')]
            with self.assertRaises(ValueError):deposited_heat(bs,(1,120))
    def test_energy_ledger(self):
        q,rows=deposited_heat([Band(1,120,.004,.001,.8,.1,3,'synthetic unit test')],(1,120))
        self.assertAlmostEqual(q,7.2e-7)
        self.assertAlmostEqual(rows[0]['heat_W_m2']+rows[0]['other_energy_W_m2'],rows[0]['absorbed_W_m2'])
    def test_historical_band_total_and_units(self):
        self.assertAlmostEqual(sum(b.irradiance_w_m2 for b in historical_solar_bands()),.00464)
    def test_historical_bands_do_not_invent_response(self):
        with self.assertRaises(ValueError):deposited_heat(historical_solar_bands(),(.1,118.))
    def test_gap_energy_without_intraband_interpolation(self):
        r=gap_energy_bounds(historical_solar_bands(),(1239.841984/50,120.181141))
        self.assertAlmostEqual(r['lower_fraction'],.375)
        self.assertAlmostEqual(r['upper_fraction'],3.79/4.64)
    def test_changed_optical_archive_rejected(self):
        with self.assertRaises(ValueError):audit_optical_files('wavelength,irradiance\n202,1\n300,1\n','0.12 1.2 0.1\n0.3 1.3 0.2\n')
    def test_invalid_transmission(self):
        with self.assertRaises(ValueError):
            deposited_heat([Band(1,120,.01,1.1,1,.1,3,'test')],(1,120))

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

if __name__=='__main__':unittest.main()
