"""Fast analytic/numerical unit checks; no large atmosphere solve required."""
import math, unittest
import numpy as np
from solver import EARTH,MOON,Atmosphere,densities,solar_columns,solar_table,ray_path,weights_path,angular_quadrature,bracket,interpolate
from render_data import cmf,to_rgb,spherical_harmonics

class GeometryAndTransport(unittest.TestCase):
    def test_vertical_gas_column(self):
        for a in [EARTH,MOON]:
            gas,_=solar_columns(a.radius_m,1.,a.radius_m,a.top,a.scale_height_m,a.ozone_scale,1200)
            self.assertAlmostEqual(gas/(a.scale_height_m*(1-math.exp(-10))),1.,delta=2e-5)
    def test_vertical_ozone_column(self):
        for a in [EARTH,MOON]:
            _,oz=solar_columns(a.radius_m,1.,a.radius_m,a.top,a.scale_height_m,a.ozone_scale,1200)
            self.assertAlmostEqual(oz/(15000*a.ozone_scale),1.,delta=3e-5)
    def test_earth_tangent_asymptotic(self):
        a=EARTH;gas,_=solar_columns(a.radius_m,0.,a.radius_m,a.top,a.scale_height_m,a.ozone_scale,2000)
        asym=math.sqrt(math.pi*a.radius_m*a.scale_height_m/2)*(1+3*a.scale_height_m/(8*a.radius_m))
        self.assertAlmostEqual(gas/asym,1.,delta=3e-5)
    def test_body_shadow(self):
        a=EARTH;x=solar_columns(a.radius_m,-.1,a.radius_m,a.top,a.scale_height_m,a.ozone_scale)
        self.assertEqual(x,(-1.,-1.))
    def test_vertical_boundary(self):
        a=MOON;sh=np.linspace(a.radius_m,a.top,20)
        p,g,c,s=ray_path(a.radius_m+1.7,-1.,a.radius_m,a.top,sh,a.scale_height_m,a.ozone_scale)
        self.assertTrue(g);self.assertAlmostEqual(p[:,6].sum(),1.7,places=5)
    def test_upper_ray(self):
        a=MOON;sh=np.linspace(a.radius_m,a.top,20)
        p,g,c,s=ray_path(a.radius_m+1.7,1.,a.radius_m,a.top,sh,a.scale_height_m,a.ozone_scale)
        self.assertFalse(g);self.assertAlmostEqual(p[:,6].sum(),a.top-a.radius_m-1.7,places=5)
    def test_segment_conservation(self):
        a=MOON;sh=a.radius_m+np.linspace(0,1,81)**2*(a.top-a.radius_m)
        for mu in [0.,.1,.6,1.]:
            p,g,c,s=ray_path(a.radius_m+1.7,mu,a.radius_m,a.top,sh,a.scale_height_m,a.ozone_scale)
            w,T=weights_path(p,np.array([1e-5,4e-5]),np.zeros(2))
            np.testing.assert_allclose(w.sum(axis=0)+T,1.,rtol=2e-14)
    def test_absorbing_segment_account(self):
        p=np.zeros((1,7));p[0,4]=1.;p[0,5]=1.;p[0,6]=100.
        w,T=weights_path(p,np.array([.01]),np.array([.03]))
        self.assertAlmostEqual(T[0],math.exp(-4),places=13)
        self.assertAlmostEqual(w[0,0],(1-math.exp(-4))*.25,places=13)
    def test_vacuum(self):
        p=np.ones((4,7));w,T=weights_path(p,np.zeros(3),np.zeros(3))
        np.testing.assert_array_equal(w,0);np.testing.assert_array_equal(T,1)
    def test_solar_half_disk_in_vacuum(self):
        a=EARTH;t=solar_table(np.array([a.radius_m]),np.array([0.]),a.radius_m,a.top,a.scale_height_m,a.ozone_scale,np.zeros(1),np.zeros(1),np.ones(1),np.deg2rad(a.sun_radius_deg),12,30)
        self.assertAlmostEqual(t[0,0,0,0],.5,places=12)
    def test_solar_disk_full_and_hidden(self):
        a=EARTH;t=solar_table(np.array([a.radius_m]),np.deg2rad(np.array([-2.,2.])),a.radius_m,a.top,a.scale_height_m,a.ozone_scale,np.zeros(1),np.zeros(1),np.ones(1),np.deg2rad(a.sun_radius_deg),12,30)
        self.assertEqual(t[0,0,0,0],0.);self.assertAlmostEqual(t[0,1,0,0],1.,places=12)
    def test_disk_contact_monotone(self):
        a=EARTH;t=solar_table(np.array([a.radius_m]),np.deg2rad(np.linspace(-.4,.4,81)),a.radius_m,a.top,a.scale_height_m,a.ozone_scale,np.zeros(1),np.zeros(1),np.ones(1),np.deg2rad(a.sun_radius_deg),12,30)
        self.assertTrue(np.all(np.diff(t[0,:,0,0])>=-1e-15))
    def test_column_ratio(self):
        self.assertAlmostEqual(MOON.scale_height_m*MOON.density_scale/EARTH.scale_height_m,1.2*9.80665/1.62)
    def test_rayleigh_phase_normalization(self):
        x,w=np.polynomial.legendre.leggauss(12)
        self.assertAlmostEqual(np.sum(3/(16*np.pi)*(1+x*x)*w)*2*np.pi,1.,places=13)
    def test_rayleigh_isotropic_field(self):
        J=4*np.pi;K=4*np.pi/3
        self.assertAlmostEqual(3/(16*np.pi)*(J+K),1.,places=14)
    def test_quadrature_full_sphere(self):
        mu,w,cp,sp,dp=angular_quadrature(20,16)
        self.assertAlmostEqual(w.sum()*len(cp)*dp,4*np.pi,places=12)
        self.assertAlmostEqual(np.sum(w*np.maximum(mu,0))*len(cp)*dp,np.pi,places=12)
    def test_bilinear_interpolation(self):
        rg=np.array([0.,1.]);ag=np.array([0.,1.]);t=np.array([[[[0.]],[[2.]]],[[[1.]],[[3.]]]])
        self.assertAlmostEqual(interpolate(t,rg,ag,.3,.4,0,0),1.1)
    def test_cie_peak(self):
        x=cmf(np.array([555.]))[0];self.assertGreater(x[1],.99);self.assertLess(x[1],1.01)
    def test_equal_energy_color(self):
        lam=np.arange(380.,801.,1.);xyz=np.trapezoid(cmf(lam),lam,axis=0)
        self.assertLess(np.max(xyz)/np.min(xyz),1.015)
    def test_lunar_clock(self):
        rate=360/(29.53059*24)
        self.assertAlmostEqual((32/60)/rate,1.0499765,places=5)
        self.assertAlmostEqual(6/rate,11.812236,places=5)
    def test_ozone_is_bounded(self):
        for h in np.linspace(0,500000,100):
            _,o=densities(h,MOON.scale_height_m,MOON.ozone_scale)
            self.assertGreaterEqual(o,0.);self.assertLessEqual(o,1.)

if __name__=='__main__':unittest.main(verbosity=2)
