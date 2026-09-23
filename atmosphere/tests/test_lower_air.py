"""Lower-air profile and A1 profile-contract tests."""
import math
import unittest
import numpy as np
from atmosphere.lower_air import LowerAir, AirConfig, A1Profile, R_MOON, GM_MOON


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


if __name__ == "__main__":
    unittest.main()
