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


if __name__ == "__main__":
    unittest.main()
