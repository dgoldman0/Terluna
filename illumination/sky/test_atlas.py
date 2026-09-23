"""Invariants of the generated clear-sky atlases in data/*_atlas.npz.

Numerical and data properties only; these tests do not validate weather or climate.
Generate the atlases with calculate_full_cycle.py (numba required).
"""
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
WORLDS=('earth','moon','moon_no_ozone')


def load_atlases(folder):
 """Load the generated atlases, or skip when they have not been built."""
 missing=[k for k in WORLDS if not (folder/f'{k}_atlas.npz').exists()]
 if missing:raise unittest.SkipTest('sky atlases not generated ('+', '.join(missing)+'); run illumination/sky/calculate_full_cycle.py')
 return {k:np.load(folder/f'{k}_atlas.npz') for k in WORLDS}


class AtlasTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.maps=load_atlases(ROOT/'data')
 def test_endpoints(self):
  for d in self.maps.values():self.assertEqual(d['suns'][0],-90);self.assertEqual(d['suns'][-1],90)
 def test_angle_sort(self):
  for d in self.maps.values():self.assertTrue(np.all(np.diff(d['suns'])>0));self.assertEqual(len(d['suns']),174)
 def test_matching_angle_grid(self):
  for d in self.maps.values():np.testing.assert_array_equal(d['suns'],self.maps['earth']['suns'])
 def test_finite_maps(self):
  for d in self.maps.values():
   for k in ('rgb','direct','direct_horizontal','diffuse','sh','cloud_direct','cloud_diffuse'):self.assertTrue(np.all(np.isfinite(d[k])),k)
 def test_photopic_nonnegative(self):
  for d in self.maps.values():
   for k in ('rgb','direct','direct_horizontal','diffuse'):self.assertGreaterEqual(float(np.min(d[k]@np.array([.2126,.7152,.0722]))),-1e-5)
 def test_night_direct_zero(self):
  for d in self.maps.values():self.assertLess(np.max(np.abs(d['direct'][d['suns']< -1])),1e-10)
 def test_noon_is_bright(self):
  for d in self.maps.values():self.assertGreater((d['direct_horizontal'][-1]+d['diffuse'][-1])@np.array([.2126,.7152,.0722]),80000)
 def test_solar_beam_noon_vertical(self):
  for d in self.maps.values():np.testing.assert_allclose(d['direct'][-1],d['direct_horizontal'][-1],rtol=1e-4)
 def test_noon_azimuth_symmetry(self):
  for d in self.maps.values():
   sky=d['rgb'][-1,42:];spread=np.ptp(sky,axis=1);self.assertLess(float(np.max(spread/np.maximum(np.abs(sky).max(1),1e-8))),.003)
 def test_extinction_no_ozone(self):
  d=self.maps['moon'];z=self.maps['moon_no_ozone'];self.assertGreaterEqual(float(z['direct'][-1]@np.array([.2126,.7152,.0722])),float(d['direct'][-1]@np.array([.2126,.7152,.0722])))
if __name__=='__main__':unittest.main(verbosity=2)
