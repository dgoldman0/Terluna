"""One-dimensional canopy momentum-balance tests."""
import unittest
import numpy as np
from atmosphere.lower_air import LowerAir
from climate.canopy_flow import CanopyFlow, CanopyConfig


class CanopyTests(unittest.TestCase):
    def test_zero_drag_constant_density_is_linear(self):
        class ConstantAir:
            def sample(self, z): return {'density_kg_m3':np.ones_like(np.asarray(z, float))}
        flow = CanopyFlow(ConstantAir(), CanopyConfig(leaf_area_index=0))
        z = np.linspace(0,300,33)
        np.testing.assert_allclose(flow.ratio(z), z/300, atol=1e-12)

    def test_momentum_balance_and_boundary(self):
        flow = CanopyFlow(LowerAir())
        self.assertLess(flow.diagnostics['momentum_relative_residual'], 1e-6)
        np.testing.assert_allclose(flow.ratio([0,300]), [0,1], atol=1e-10)
        self.assertTrue(np.all(np.diff(flow.ratio(np.linspace(0,600,300))) >= 0))
        self.assertAlmostEqual(flow.ratio(300-1e-5), flow.ratio(300+1e-5), places=5)

    def test_drag_mixing_and_scaling(self):
        a = LowerAir()
        low = CanopyFlow(a, CanopyConfig(leaf_area_index=5))
        high = CanopyFlow(a, CanopyConfig(leaf_area_index=22))
        mixed = CanopyFlow(a, CanopyConfig(leaf_area_index=22, mixing_k_over_uh=.12))
        self.assertLess(high.ratio(150), low.ratio(150))
        self.assertGreater(mixed.ratio(150), high.ratio(150))
        np.testing.assert_allclose(high.sample([50,150,300],20), 2*high.sample([50,150,300],10))
        np.testing.assert_allclose(high.sample([50,150,300],0), 0)

    def test_canopy_invalid_inputs(self):
        for cfg in (CanopyConfig(leaf_area_index=-1),CanopyConfig(mixing_k_over_uh=0),
                    CanopyConfig(above_displacement_fraction=.95),CanopyConfig(height_m=6000)):
            with self.assertRaises(ValueError): CanopyFlow(LowerAir(),cfg)
        with self.assertRaises(ValueError): CanopyFlow(LowerAir()).sample(30,-1)


if __name__ == "__main__":
    unittest.main()
