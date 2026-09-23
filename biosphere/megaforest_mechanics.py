"""Conditional static wind screening for tapered megatrees, in SI units.

Euler--Bernoulli beam elements include linearized geometric stiffness from
self-weight. Crown streamlining is a bounded empirical response law. Anchorage
and branch thresholds are explicit idealized requirements, not calibrated tree
predictions. No gust resonance, fatigue, growth, hydraulics or post-damage solve.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.optimize import brentq


@dataclass(frozen=True)
class TreeConfig:
    height_m: float = 300.0
    base_diameter_m: float = 300.0/35
    tip_diameter_fraction: float = .12
    crown_base_fraction: float = .45
    crown_radius_m: float = 54.0
    crown_frontal_fraction: float = .30
    crown_drag_coefficient: float = .8
    stem_drag_coefficient: float = 1.0
    vogel_exponent: float = 0.0
    reconfiguration_start_m_s: float = 5.0
    minimum_drag_fraction: float = .2
    young_modulus_pa: float = 8e9
    bending_strength_pa: float = 20e6
    compression_strength_pa: float = 20e6
    wood_density_kg_m3: float = 650.0
    crown_wood_to_stem_mass: float = .25
    leaf_area_index: float = 22.0
    leaf_mass_kg_m2_leaf: float = .30
    epiphyte_mass_kg_m2_footprint: float = 35.0
    retained_water_mm: float = 30.0
    branch_count: int = 32
    branch_diameter_fraction: float = .18
    branch_span_to_crown_radius: float = .6
    branch_strength_pa: float = 20e6
    elements: int = 60

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite tree inputs required')
        for k in ('height_m', 'base_diameter_m', 'crown_radius_m', 'young_modulus_pa',
                  'bending_strength_pa', 'compression_strength_pa', 'wood_density_kg_m3',
                  'branch_strength_pa', 'branch_span_to_crown_radius', 'reconfiguration_start_m_s'):
            if getattr(self, k) <= 0:
                raise ValueError(k+' must be positive')
        if not 0 < self.tip_diameter_fraction <= 1 or not 0 <= self.crown_base_fraction < 1:
            raise ValueError('Invalid taper or crown base')
        if not 0 < self.branch_diameter_fraction <= 1 or not 0 <= self.crown_frontal_fraction <= 1:
            raise ValueError('Invalid branch/crown geometry')
        if not -1.5 <= self.vogel_exponent <= 0 or not 0 < self.minimum_drag_fraction <= 1:
            raise ValueError('Invalid bounded drag-response law')
        for k in ('crown_drag_coefficient', 'stem_drag_coefficient', 'crown_wood_to_stem_mass',
                  'leaf_area_index', 'leaf_mass_kg_m2_leaf', 'epiphyte_mass_kg_m2_footprint', 'retained_water_mm'):
            if getattr(self, k) < 0:
                raise ValueError(k+' must be nonnegative')
        if (not isinstance(self.elements, int) or not 8 <= self.elements <= 500
                or not isinstance(self.branch_count, int) or self.branch_count < 1):
            raise ValueError('Invalid beam resolution or branch count')


@dataclass(frozen=True)
class RootConfig:
    plate_radius_m: float = 30.0
    plate_depth_m: float = 4.0
    soil_density_kg_m3: float = 1600.0
    cohesion_pa: float = 10000.0
    friction_angle_deg: float = 30.0
    pore_pressure_ratio: float = .2
    weight_lever_fraction: float = .3
    shear_mobilized_fraction: float = .25
    shear_lever_fraction: float = .5
    root_tensile_area_to_stem_area: float = .1
    root_tensile_strength_pa: float = 8e6
    root_lever_fraction: float = .5

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite root inputs required')
        if min(self.plate_radius_m, self.plate_depth_m, self.soil_density_kg_m3,
               self.root_tensile_strength_pa, self.root_tensile_area_to_stem_area) <= 0:
            raise ValueError('Positive root geometry/material required')
        if self.cohesion_pa < 0 or not 0 <= self.friction_angle_deg < 80:
            raise ValueError('Invalid soil strength')
        for k in ('pore_pressure_ratio', 'weight_lever_fraction', 'shear_mobilized_fraction',
                  'shear_lever_fraction', 'root_lever_fraction'):
            if not 0 <= getattr(self, k) <= 1:
                raise ValueError('Invalid root coefficient '+k)


def root_capacity(root, tree, gravity):
    """Two serial load-path bottlenecks: root tissue and mobilized root-soil plate.

    tau=c+sigma_eff*tan(phi); mobilized area and moment arms are prescribed.
    Root-tissue and soil capacities are NOT added. The plate weight and the
    selected shear contribution form one explicitly idealized soil mechanism.
    """
    root.validate()
    if not np.isfinite(gravity) or gravity <= 0:
        raise ValueError('Positive finite gravity required')
    r, d = root.plate_radius_m, root.plate_depth_m
    soil_mass = math.pi*r*r*d*root.soil_density_kg_m3
    normal = .5*root.soil_density_kg_m3*gravity*d*(1-root.pore_pressure_ratio)
    tau = root.cohesion_pa+normal*math.tan(math.radians(root.friction_angle_deg))
    weight_moment = root.weight_lever_fraction*soil_mass*gravity*r
    shear_moment = root.shear_mobilized_fraction*math.pi*r*r*tau*root.shear_lever_fraction*r
    tissue_moment = (root.root_tensile_area_to_stem_area*math.pi*tree.base_diameter_m**2/4
                     *root.root_tensile_strength_pa*root.root_lever_fraction*r)
    soil_moment = weight_moment+shear_moment
    if min(tissue_moment, soil_moment) <= 0:
        raise ValueError('Anchorage requires a positive resisting mechanism')
    return dict(capacity_nm=min(tissue_moment, soil_moment), soil_capacity_nm=soil_moment,
                tissue_capacity_nm=tissue_moment, soil_weight_moment_nm=weight_moment,
                soil_shear_moment_nm=shear_moment, effective_shear_pa=tau,
                mobilized_soil_mass_kg=soil_mass,
                bottleneck='root_tissue' if tissue_moment < soil_moment else 'root_soil_plate')


def beam_matrices(z, rigidity, axial_force):
    """Consistent 2-DOF/node beam stiffness and compressive geometric stiffness."""
    z = np.asarray(z, float)
    rigidity, axial_force = np.asarray(rigidity, float), np.asarray(axial_force, float)
    n = len(z)-1
    if (n < 1 or rigidity.shape != (n,) or axial_force.shape != (n,)
            or not np.all(np.isfinite(z)) or not np.all(np.isfinite(rigidity))
            or not np.all(np.isfinite(axial_force)) or np.any(np.diff(z) <= 0)
            or np.any(rigidity <= 0) or np.any(axial_force < 0)):
        raise ValueError('Invalid beam grid, rigidity or axial force')
    k, kg = np.zeros((2*n+2, 2*n+2)), np.zeros((2*n+2, 2*n+2))
    for i, length in enumerate(np.diff(z)):
        l = length
        ke = rigidity[i]/l**3*np.array([[12,6*l,-12,6*l], [6*l,4*l*l,-6*l,2*l*l],
                                       [-12,-6*l,12,-6*l], [6*l,2*l*l,-6*l,4*l*l]])
        ge = axial_force[i]/(30*l)*np.array([[36,3*l,-36,3*l], [3*l,4*l*l,-3*l,-l*l],
                                            [-36,-3*l,36,-3*l], [3*l,-l*l,-3*l,4*l*l]])
        k[2*i:2*i+4, 2*i:2*i+4] += ke
        kg[2*i:2*i+4, 2*i:2*i+4] += ge
    return k, kg


def distributed_beam_load(z, force_per_m):
    z, force_per_m = np.asarray(z, float), np.asarray(force_per_m, float)
    if (z.ndim != 1 or len(z) < 2 or force_per_m.shape != (len(z)-1,)
            or not np.all(np.isfinite(z)) or not np.all(np.isfinite(force_per_m))
            or np.any(np.diff(z) <= 0)):
        raise ValueError('Invalid distributed load grid or values')
    out = np.zeros(2*len(z))
    for i, (l, q) in enumerate(zip(np.diff(z), force_per_m)):
        out[2*i:2*i+4] += q*np.array([l/2, l*l/12, l/2, -l*l/12])
    return out


def _tail(values):
    return np.r_[np.cumsum(values[::-1])[::-1], 0.0]


class StaticTree:
    """Fixed-base, small-displacement intact-tree screen under imposed winds."""
    def __init__(self, air, tree=TreeConfig(), roots=RootConfig(), wind_ratio=None):
        tree.validate()
        roots.validate()
        self.tree, self.roots = tree, roots
        self.z = np.linspace(0, tree.height_m, tree.elements+1)
        self.dz = np.diff(self.z)
        self.mid = (self.z[1:]+self.z[:-1])/2
        t = tree
        diameter = lambda z: t.base_diameter_m*(1-(1-t.tip_diameter_fraction)*z/t.height_m)
        self.diameter, dm = diameter(self.z), diameter(self.mid)
        self.area = math.pi*self.diameter**2/4
        self.section_modulus = math.pi*self.diameter**3/32
        air_mid = air.sample(self.mid)
        self.rho = air_mid['density_kg_m3']
        self.gravity = air_mid['gravity_m_s2']
        self.stem_mass_segments = t.wood_density_kg_m3*math.pi*dm**2/4*self.dz
        self.stem_mass = float(self.stem_mass_segments.sum())
        footprint = math.pi*t.crown_radius_m**2
        self.crown_mass = (t.crown_wood_to_stem_mass*self.stem_mass+footprint*(
            t.leaf_area_index*t.leaf_mass_kg_m2_leaf+t.epiphyte_mass_kg_m2_footprint+t.retained_water_mm))
        # 1 mm water = 1 kg/m2. LAD is not used as drag area.
        s = (self.mid/t.height_m-t.crown_base_fraction)/(1-t.crown_base_fraction)
        crown_shape = np.where((s > 0) & (s < 1), np.sqrt(np.maximum(0, 4*s*(1-s))), 0)
        fractions = crown_shape*self.dz
        if fractions.sum() <= 0:
            raise ValueError('Beam grid does not resolve the specified crown')
        fractions /= fractions.sum()
        self.crown_mass_segments = self.crown_mass*fractions
        self.crown_width = 2*t.crown_radius_m*crown_shape*t.crown_frontal_fraction
        self.stem_width = dm
        self.weights = (self.stem_mass_segments+self.crown_mass_segments)*self.gravity
        self.axial = _tail(self.weights)
        axial_mid = (self.axial[1:]+self.axial[:-1])/2
        k, kg = beam_matrices(self.z, t.young_modulus_pa*math.pi*dm**4/64, axial_mid)
        kr, gr = k[2:, 2:], kg[2:, 2:]
        last = len(kr)-1
        self.buckling_utilization = float(eigh(gr, kr, subset_by_index=[last,last], eigvals_only=True)[0])
        self.factor = None if self.buckling_utilization >= 1 else cho_factor(kr-gr, lower=True)
        self.anchorage = root_capacity(roots, tree, float(air.sample(0)['gravity_m_s2']))
        self.wind_ratio = np.ones_like(self.mid) if wind_ratio is None else np.asarray(wind_ratio(self.mid), float)
        if self.wind_ratio.shape != self.mid.shape or not np.all(np.isfinite(self.wind_ratio)) or np.any(self.wind_ratio < 0):
            raise ValueError('Invalid prescribed wind shape')
        span = t.crown_radius_m*t.branch_span_to_crown_radius
        self.branch_gravity_moment = float(np.sum(self.crown_mass_segments*self.gravity))*span/(2*t.branch_count)
        self.branch_section_modulus = math.pi*(t.base_diameter_m*t.branch_diameter_fraction)**3/32
        self.span = span

    def evaluate(self, speed_m_s, include_profile=False):
        if not np.isfinite(speed_m_s) or speed_m_s < 0:
            raise ValueError('Finite nonnegative imposed wind speed required')
        if self.factor is None:
            return dict(status='SELF_WEIGHT_BUCKLING', reference_speed_m_s=float(speed_m_s),
                        buckling_utilization=self.buckling_utilization)
        t = self.tree
        u = speed_m_s*self.wind_ratio
        response = np.maximum(t.minimum_drag_fraction,
            np.maximum(1, u/t.reconfiguration_start_m_s)**t.vogel_exponent)
        qc = .5*self.rho*u*u*t.crown_drag_coefficient*self.crown_width*response
        qs = .5*self.rho*u*u*t.stem_drag_coefficient*self.stem_width
        q = qc+qs
        load = distributed_beam_load(self.z, q)
        dofs = np.r_[0., 0., cho_solve(self.factor, load[2:])]
        y, slope = dofs[::2], dofs[1::2]
        ym = (y[:-1]+y[1:])/2+self.dz*(slope[:-1]-slope[1:])/8
        wind_moment = _tail(q*self.dz*self.mid)-self.z*_tail(q*self.dz)
        gravity_moment = _tail(self.weights*ym)-y*self.axial
        moment = wind_moment+gravity_moment
        bending = np.abs(moment)/self.section_modulus
        compression = self.axial/self.area
        stem_util = np.maximum((bending+compression)/t.compression_strength_pa,
                               np.maximum(0, bending-compression)/t.bending_strength_pa)
        crown_force = float(np.sum(qc*self.dz))
        branch_wind_moment = crown_force*self.span/(2*t.branch_count)
        branch_util = math.hypot(branch_wind_moment, self.branch_gravity_moment)/(
            self.branch_section_modulus*t.branch_strength_pa)
        domain = max(float(max(abs(y)))/(.1*t.height_m), float(max(abs(slope)))/.2)
        flags = ['HYPOTHETICAL_TRAITS', 'IMPOSED_STEADY_WIND', 'FIXED_BASE_UNTIL_THRESHOLD',
                 'IDEALIZED_ANCHORAGE', 'INTACT_GEOMETRY_NO_POST_DAMAGE_RESPONSE']
        if domain > 1:
            flags.append('SMALL_DISPLACEMENT_DOMAIN_EXCEEDED')
        if branch_util >= 1:
            flags.append('CROWN_DAMAGE_ONSET_EXCEEDED')
        result = dict(status='CONDITIONAL_STATIC_SCREEN', reference_speed_m_s=float(speed_m_s),
            stem_utilization=float(max(stem_util)), root_utilization=float(moment[0]/self.anchorage['capacity_nm']),
            crown_utilization=branch_util, domain_utilization=domain,
            buckling_utilization=self.buckling_utilization,
            base_wind_moment_nm=float(wind_moment[0]), base_gravity_moment_nm=float(gravity_moment[0]),
            base_total_moment_nm=float(moment[0]), crown_force_n=crown_force,
            total_drag_n=float(np.sum(q*self.dz)), tip_displacement_m=float(y[-1]),
            max_rotation_rad=float(max(abs(slope))), critical_stem_height_m=float(self.z[np.argmax(stem_util)]),
            stem_mass_kg=self.stem_mass, crown_mass_kg=self.crown_mass,
            anchorage_capacity_nm=self.anchorage['capacity_nm'], flags=flags)
        if include_profile:
            result['profile'] = dict(z_m=self.z, displacement_m=y, rotation_rad=slope,
                bending_moment_nm=moment, axial_force_n=self.axial,
                stem_utilization=stem_util, wind_sample_height_m=self.mid, wind_speed_m_s=u)
        return result

    def envelope(self, max_speed_m_s=100.0):
        if not np.isfinite(max_speed_m_s) or max_speed_m_s <= 0:
            raise ValueError('Positive wind search bound required')
        if self.factor is None:
            return dict(status='SELF_WEIGHT_BUCKLING', buckling_utilization=self.buckling_utilization)
        cache = {}
        def at(u):
            if u not in cache:
                cache[u] = self.evaluate(u)
            return cache[u]
        result = dict(status='CONDITIONAL_STATIC_ENVELOPE', maximum_search_speed_m_s=max_speed_m_s,
                      buckling_utilization=self.buckling_utilization)
        for mode in ('stem', 'root', 'crown', 'domain'):
            key = mode+'_utilization'
            if at(0)[key] >= 1:
                value = 0.0
            elif at(max_speed_m_s)[key] < 1:
                value = None
            else:
                value = float(brentq(lambda u: at(u)[key]-1, 0, max_speed_m_s, xtol=1e-7))
            result[mode+'_threshold_m_s'] = value
        domain = result['domain_threshold_m_s']
        crossings = [(result[m+'_threshold_m_s'], m) for m in ('stem','root','crown') if result[m+'_threshold_m_s'] is not None]
        speed, mode = min(crossings) if crossings else (None, 'right_censored')
        result.update(first_threshold_m_s=speed, first_mode=mode,
            first_threshold_within_displacement_domain=bool(speed is not None and (domain is None or speed <= domain)),
            interpretation='thresholds after first damage retain the original intact geometry')
        for mode in ('stem','root','crown'):
            u = result[mode+'_threshold_m_s']
            result[mode+'_threshold_within_domain'] = bool(u is not None and (domain is None or u <= domain))
        return result
