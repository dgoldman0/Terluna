"""Conditional 1-D mean-canopy momentum balance; no weather or gust prediction.

Solve d(rho K dU/dz)/dz = rho Cd a_front U|U|/2 with U(0)=0, U(Hc)=Uref.
K/(Uref Hc) is a prescribed positive mixing parameter, not a turbulence closure.
LAI is converted to frontal area density explicitly by a projection factor.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_bvp, simpson


@dataclass(frozen=True)
class CanopyConfig:
    height_m: float = 300.0
    leaf_area_index: float = 22.0
    crown_base_fraction: float = .4
    projection_factor: float = .3
    element_drag_coefficient: float = 1.0
    mixing_k_over_uh: float = .04
    above_displacement_fraction: float = .65
    above_roughness_fraction: float = .1

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite canopy parameters required')
        if self.height_m <= 0 or self.leaf_area_index < 0 or self.mixing_k_over_uh <= 0:
            raise ValueError('Invalid canopy size, leaf area or mixing')
        if not 0 <= self.crown_base_fraction < 1 or not 0 <= self.projection_factor <= 1 or self.element_drag_coefficient < 0:
            raise ValueError('Invalid crown or frontal drag parameters')
        if (not 0 <= self.above_displacement_fraction < 1
                or not 0 < self.above_roughness_fraction < 1-self.above_displacement_fraction):
            raise ValueError('Invalid above-canopy log-profile geometry')


class CanopyFlow:
    def __init__(self, air, config=CanopyConfig(), tolerance=1e-7):
        config.validate()
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Positive BVP tolerance required')
        self.config = config
        h, k = config.height_m, config.mixing_k_over_uh
        rho0 = float(air.sample(0)['density_kg_m3'])
        def density(x):
            return air.sample(x*h)['density_kg_m3']/rho0
        def area(x):
            s = np.clip((x-config.crown_base_fraction)/(1-config.crown_base_fraction), 0, 1)
            return 6*s*(1-s)*config.leaf_area_index*config.projection_factor/(1-config.crown_base_fraction)
        def rhs(x, y):
            r = density(x)
            return np.vstack((y[1]/r, r*.5*config.element_drag_coefficient*area(x)*y[0]*abs(y[0])/k))
        def bc(a, b):
            return np.array([a[0], b[0]-1])
        mesh = np.unique(np.r_[np.linspace(0, 1, 121), config.crown_base_fraction])
        self.solution = solve_bvp(rhs, bc, mesh, np.vstack((mesh, np.ones_like(mesh))),
                                  tol=tolerance, max_nodes=10000)
        if not self.solution.success:
            raise RuntimeError('Canopy momentum BVP failed: '+self.solution.message)
        x = np.linspace(0, 1, 4001)
        u, stress = self.solution.sol(x)
        if np.min(u) < -1e-8 or np.max(u) > 1+1e-8:
            raise RuntimeError('Unphysical canopy velocity solution')
        drag = .5*config.element_drag_coefficient*density(x)*area(x)*u*abs(u)
        integral = float(simpson(drag, x=x))
        net_stress = float(k*(stress[-1]-stress[0]))
        self.diagnostics = dict(status='CONDITIONAL_STEADY_MEAN_FLOW',
            max_collocation_residual=float(max(self.solution.rms_residuals)),
            dimensionless_drag_integral=integral, dimensionless_stress_difference=net_stress,
            momentum_relative_residual=abs(integral-net_stress)/max(abs(integral), 1e-12),
            canopy_config=asdict(config), turbulence_closure='prescribed_constant_eddy_viscosity',
            pressure_gradient=False, gusts=False, climate_prediction=False)

    def ratio(self, z):
        z = np.asarray(z, float)
        if not np.all(np.isfinite(z)) or np.any(z < 0):
            raise ValueError('Finite nonnegative wind heights required')
        c = self.config
        x = z/c.height_m
        inside = self.solution.sol(np.minimum(x, 1))[0]
        outside = np.log((np.maximum(x, 1)-c.above_displacement_fraction)/c.above_roughness_fraction)
        outside /= np.log((1-c.above_displacement_fraction)/c.above_roughness_fraction)
        return np.maximum(0, np.where(x <= 1, inside, outside))

    def sample(self, z, reference_speed_m_s):
        if not np.isfinite(reference_speed_m_s) or reference_speed_m_s < 0:
            raise ValueError('Finite nonnegative imposed reference speed required')
        return reference_speed_m_s*self.ratio(z)
