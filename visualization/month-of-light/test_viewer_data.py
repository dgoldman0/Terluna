"""Display-model and landscape-asset invariants of the month-of-light viewer.

The viewer's weather, fog, water and clock expressions are display approximations
(see METHODS.md); these tests check their stated identities and bounds only.
"""
import unittest,math
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parent
SKY=ROOT.parents[1]/'illumination'/'sky'
WORLDS=('earth','moon','moon_no_ozone')


def load_atlases(folder):
 """Load the generated atlases, or skip when they have not been built."""
 missing=[k for k in WORLDS if not (folder/f'{k}_atlas.npz').exists()]
 if missing:raise unittest.SkipTest('sky atlases not generated ('+', '.join(missing)+'); run illumination/sky/calculate_full_cycle.py')
 return {k:np.load(folder/f'{k}_atlas.npz') for k in WORLDS}


class DisplayModelTests(unittest.TestCase):
 def test_quarter_cycle_angles(self):
  a=np.rad2deg(np.arcsin(np.cos(np.linspace(0,1,5)*2*np.pi)));np.testing.assert_allclose(a,[90,0,-90,0,90],atol=1e-10)
 def test_daylight_halves(self):
  phase=np.linspace(0,1,1000,endpoint=False)+.00001;angles=np.arcsin(np.cos(phase*2*np.pi));self.assertEqual(np.sum(angles>0),500)
 def test_clock_ratio(self):self.assertAlmostEqual((29.53059*24*3600)/86400,29.53059)
 def test_cloud_slab_bounds(self):
  for tau in np.linspace(0,40,81):
   T=1/(1+.75*.15*tau);self.assertTrue(0<T<=1)
   for mu in [.035,.1,.5,1.]:self.assertLessEqual(math.exp(-tau/mu),T+1e-15)
 def test_fog_identity_and_bounds(self):
  d=np.array([0.,10,100,1000]);np.testing.assert_array_equal(np.exp(-d*0),np.ones(4));self.assertTrue(np.all(np.diff(np.exp(-3.912/220*d))<0))
 def test_depth_codec(self):
  d=np.geomspace(1,19000,1000);q=np.floor(np.log2(1+d)/16*65535);v=2**(q/65535*16)-1;self.assertLess(float(np.max(abs(v-d)/d)),.00035)
 def test_wave_dispersion_ratio(self):self.assertAlmostEqual(math.sqrt(1.62*.616)/math.sqrt(9.80665*.616),math.sqrt(1.62/9.80665))


class WeatherOnAtlasTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.maps=load_atlases(SKY/'data')
 def test_clear_weather_identity(self):
  for d in self.maps.values():
   F=d['diffuse'];H=d['direct_horizontal'];C=0.;tau=0.;T=1/(1+.75*.15*tau);B=math.exp(-tau/.5);np.testing.assert_array_equal(F*(1-C+C*T)+H*C*max(0,T-B),F)


class LandscapeAssetTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  if not (ROOT/'data'/'landscape'/'geometry.png').exists():raise unittest.SkipTest('landscape textures not generated; run make_landscape.py')
 def test_landscape_assets(self):
  for name in ('albedo','normal','geometry','shadow'):
   im=Image.open(ROOT/'data'/'landscape'/f'{name}.png');self.assertEqual(im.size,(3072,1536));self.assertEqual(im.mode,'RGB')
 def test_kind_palette(self):
  im=np.array(Image.open(ROOT/'data'/'landscape'/'geometry.png'));self.assertTrue(set(np.unique(im[:,:,2])).issubset({0,1,2,3,4,5}))
if __name__=='__main__':unittest.main(verbosity=2)
