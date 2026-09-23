"""Three-dimensional porous-drag patch flow tests."""
import unittest
from dataclasses import replace
import numpy as np
from climate.forest_patch_flow import FlowConfig, PatchGrid


class FlowTests(unittest.TestCase):
    def setUp(self):self.grid=PatchGrid(FlowConfig(nx=24,ny=16,nz=12))

    def test_empty_wind_tunnel_is_uniform(self):
        r=self.grid.solve(np.zeros(self.grid.shape))
        np.testing.assert_allclose(r.centered[...,0],1.,atol=1e-13)
        np.testing.assert_allclose(r.centered[...,1:],0.,atol=1e-13)
        self.assertEqual(r.diagnostics['steps'],1)

    def test_projection_removes_gradient(self):
        g=self.grid;x,y,z=np.meshgrid(*g.centers,indexing='ij')
        p=np.cos(np.pi*(x+6)/12)*np.cos(np.pi*(z)/3)
        v=g.faces()
        for a in range(3):
            s=[slice(None)]*3;s[a]=slice(1,-1)
            v[a][tuple(s)]+=np.diff(p,axis=a)/g.spacing[a]
        v,_=g.project(v)
        for actual,expected in zip(v,g.faces()):np.testing.assert_allclose(actual,expected,atol=1e-12)

    def test_projection_random_field(self):
        g=self.grid;rng=np.random.default_rng(42)
        v=[a+rng.normal(0,.1,a.shape) for a in g.faces()];g._boundary(v)
        v,_=g.project(v)
        self.assertLess(np.max(abs(g.divergence(v))),1e-12)

    def test_incompatible_flux_rejected(self):
        v=self.grid.faces();v[0][0]=2
        with self.assertRaises(ValueError):self.grid.project(v)

    def test_negative_or_nan_drag_rejected(self):
        for bad in (-1.,np.nan):
            with self.assertRaises(ValueError):self.grid.solve(np.full(self.grid.shape,bad))

    def test_drag_slab_pressure_balance(self):
        g=self.grid;b=2.
        r=g.solve(np.full(g.shape,b))
        np.testing.assert_allclose(r.centered[...,0],1.,atol=1e-12)
        np.testing.assert_allclose(np.diff(r.pressure,axis=0)/g.spacing[0],-b/2,atol=1e-11)

    def test_porosity_deflects_flow_and_closes_momentum(self):
        g=self.grid;x,y,z=np.meshgrid(*g.centers,indexing='ij');b=np.zeros(g.shape)
        b[(abs(x)<1)&(abs(y)<1)&(z>.3)&(z<1.2)]=3
        r=g.solve(b)
        self.assertLess(r.centered[b>0,0].mean(),.95)
        self.assertGreater(np.max(abs(r.centered[...,2])),.01)
        self.assertLess(r.diagnostics['momentum_relative_residual'],1e-4)
        self.assertLess(r.diagnostics['max_divergence'],1e-11)

    def test_nonconvergence_raises(self):
        g=PatchGrid(replace(self.grid.config,max_steps=4,tolerance=1e-15))
        b=np.zeros(g.shape);b[8:12,6:9,1:5]=4
        with self.assertRaises(RuntimeError):g.solve(b)


if __name__ == "__main__":
    unittest.main()
