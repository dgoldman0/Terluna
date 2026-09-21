"""Numerical and data invariants. These tests do not validate weather or climate."""
import unittest,json,math
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parent
class FullCycleTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.maps={k:np.load(ROOT/'data'/f'{k}_atlas.npz') for k in ('earth','moon','moon_no_ozone')}
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
 def test_quarter_cycle_angles(self):
  a=np.rad2deg(np.arcsin(np.cos(np.linspace(0,1,5)*2*np.pi)));np.testing.assert_allclose(a,[90,0,-90,0,90],atol=1e-10)
 def test_daylight_halves(self):
  phase=np.linspace(0,1,1000,endpoint=False)+.00001;angles=np.arcsin(np.cos(phase*2*np.pi));self.assertEqual(np.sum(angles>0),500)
 def test_clock_ratio(self):self.assertAlmostEqual((29.53059*24*3600)/86400,29.53059)
 def test_clear_weather_identity(self):
  for d in self.maps.values():
   F=d['diffuse'];H=d['direct_horizontal'];C=0.;tau=0.;T=1/(1+.75*.15*tau);B=math.exp(-tau/.5);np.testing.assert_array_equal(F*(1-C+C*T)+H*C*max(0,T-B),F)
 def test_cloud_slab_bounds(self):
  for tau in np.linspace(0,40,81):
   T=1/(1+.75*.15*tau);self.assertTrue(0<T<=1)
   for mu in [.035,.1,.5,1.]:self.assertLessEqual(math.exp(-tau/mu),T+1e-15)
 def test_fog_identity_and_bounds(self):
  d=np.array([0.,10,100,1000]);np.testing.assert_array_equal(np.exp(-d*0),np.ones(4));self.assertTrue(np.all(np.diff(np.exp(-3.912/220*d))<0))
 def test_landscape_assets(self):
  for name in ('albedo','normal','geometry','shadow'):
   im=Image.open(ROOT/'data'/'landscape'/f'{name}.png');self.assertEqual(im.size,(3072,1536));self.assertEqual(im.mode,'RGB')
 def test_kind_palette(self):
  im=np.array(Image.open(ROOT/'data'/'landscape'/'geometry.png'));self.assertTrue(set(np.unique(im[:,:,2])).issubset({0,1,2,3,4,5}))
 def test_depth_codec(self):
  d=np.geomspace(1,19000,1000);q=np.floor(np.log2(1+d)/16*65535);v=2**(q/65535*16)-1;self.assertLess(float(np.max(abs(v-d)/d)),.00035)
 def test_wave_dispersion_ratio(self):self.assertAlmostEqual(math.sqrt(1.62*.616)/math.sqrt(9.80665*.616),math.sqrt(1.62/9.80665))
if __name__=='__main__':unittest.main(verbosity=2)
