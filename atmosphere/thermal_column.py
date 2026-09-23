"""Steady spherical conduction/advection/Jeans column: an explicitly molecular limit.

The lower atmosphere and homopause temperature are boundary inputs. The upper
profile, exobase and molecular losses are solved, with no imposed exobase floor.
Heating is *net deposited power per lunar surface area*, not incident sunlight.
A valid spectrum-to-heating conversion must be supplied separately. This is NOT
a photochemical or kinetic closure: prescribed well-mixed N2/O2, no ions/atoms,
no explicit IR cooling, tides, eddy transport or momentum acceleration.

Equations and domain flags are documented in research/findings.md.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

from shared.constants import AVOGADRO as NA, BOLTZMANN as KB, MOON_GM as GM, MOON_RADIUS as R

AREA = 4 * math.pi * R**2
MOLAR = np.array([0.0280134, 0.031998])
RS_SPECIES = KB * NA / MOLAR

@dataclass(frozen=True)
class ColumnConfig:
    surface_pressure_pa: float = 121590.0
    surface_temperature_k: float = 288.0
    lower_temperature_k: float = 180.0
    lower_pressure_pa: float = 0.1
    oxygen_mole_fraction: float = 0.175
    lower_cp_j_kg_k: float = 1005.0
    collision_cross_section_m2: float = 5.5e-19  # inherited CO proxy, sensitivity input
    conductivity_coefficient: float = 9.37e-5  # kappa=c*T W/(m K); N2 proxy
    conductivity_multiplier: float = 1.0
    heating_shape: str = 'middle'  # deposited power distribution in log-pressure coordinate
    mesh_points: int = 81
    tolerance: float = 1e-6
    max_nodes: int = 5000

    def validate(self):
        values = [self.surface_pressure_pa, self.surface_temperature_k,
                  self.lower_temperature_k, self.lower_pressure_pa, self.lower_cp_j_kg_k,
                  self.collision_cross_section_m2, self.conductivity_coefficient,
                  self.conductivity_multiplier, self.tolerance]
        if not all(math.isfinite(v) and v > 0 for v in values):
            raise ValueError('Positive finite physical inputs required')
        if not 0 < self.oxygen_mole_fraction < 1:
            raise ValueError('Molecular mixture requires 0 < oxygen fraction < 1')
        if self.lower_pressure_pa >= self.surface_pressure_pa:
            raise ValueError('Column base must be above surface')
        if self.lower_temperature_k > self.surface_temperature_k:
            raise ValueError('This lower-profile prescription requires Tb <= Ts')
        if self.heating_shape not in ('middle', 'low', 'high'):
            raise ValueError('Unknown deposition shape')
        if self.mesh_points < 21 or self.max_nodes < self.mesh_points:
            raise ValueError('Invalid BVP grid')


def heating_cdf(s, shape):
    s = np.asarray(s)
    if shape == 'middle': return 3*s*s - 2*s**3
    if shape == 'low': return 1-(1-s)**3
    if shape == 'high': return s**3
    raise ValueError('Unknown heating shape')


def lower_boundary(cfg: ColumnConfig):
    cfg.validate()
    fractions = np.array([1-cfg.oxygen_mole_fraction, cfg.oxygen_mole_fraction])
    mean_molar = float(fractions @ MOLAR)
    rs = KB*NA/mean_molar
    exponent = rs/cfg.lower_cp_j_kg_k
    x = math.log(cfg.surface_pressure_pa/cfg.lower_pressure_pa)
    xc = math.log(cfg.surface_temperature_k/cfg.lower_temperature_k)/exponent
    integral = cfg.surface_temperature_k/exponent * (-math.expm1(-exponent*min(x,xc)))
    integral += cfg.lower_temperature_k*max(0, x-xc)
    inverse_radius_scaled = 1 - R*rs/GM*integral
    if inverse_radius_scaled <= 0: raise ValueError('Lower profile has no finite radius')
    return inverse_radius_scaled, rs, fractions, mean_molar


def jeans_boundary(u, temperature, log_pressure, fractions):
    """Mass fluxes per *lunar surface area*, and energy at infinity per kg.

    Escape-weighted translational energy at infinity is Rs*T*(lambda+2)/(lambda+1).
    Add Rs*T for the two active rotational degrees of a cold diatomic molecule.
    """
    lam = GM*u/(R*RS_SPECIES*temperature)
    log_j = (log_pressure+np.log(fractions)-.5*np.log(RS_SPECIES*temperature)
             -.5*math.log(2*math.pi)+np.log1p(lam)-lam-2*math.log(u))
    e_inf = RS_SPECIES*temperature*((lam+2)/(lam+1)+1)
    return log_j, e_inf, lam


def solve_column(deposited_heat_w_m2: float, cfg: ColumnConfig = ColumnConfig(),
                 previous=None):
    """Return summary, profile and warm-start solution; raise on failed convergence.

    'previous' is only a numerical warm start. Each returned case is independently
    checked for boundary residuals and flagged for physical approximation limits.
    """
    if not math.isfinite(deposited_heat_w_m2) or deposited_heat_w_m2 < 0:
        raise ValueError('Deposited heat must be nonnegative and finite')
    ub, rs, fractions, mean_molar = lower_boundary(cfg)
    tb, pb, q = cfg.lower_temperature_k, cfg.lower_pressure_pa, deposited_heat_w_m2
    log_kn_const = math.log(GM*(mean_molar/NA)/(cfg.collision_cross_section_m2*pb*R**2))
    c = R*rs*tb/GM
    xmax = ub/c-2
    isothermal_kn = lambda x: log_kn_const+x+2*math.log(ub-c*x)
    if xmax <= 0 or isothermal_kn(xmax) <= 0:
        raise ValueError('No cold molecular exobase for numerical initialization')
    L0 = brentq(isothermal_kn,0,xmax)
    u1 = ub-c*L0
    lj, ei, _ = jeans_boundary(u1,tb,math.log(pb)-L0,fractions)
    scale = max(q,1e-8)
    mesh = np.linspace(0,1,cfg.mesh_points)
    init = np.vstack([ub-c*L0*mesh, np.ones_like(mesh)])
    par = np.r_[lj,math.log(L0),(np.exp(lj)@ei-q)/scale]
    if previous is not None:
        init = previous.sol(mesh)
        par = previous.p.copy()
        par[3] *= previous.heat_scale/scale

    def rhs(s,y,p):
        # Trial-iterate guards prevent NaN; final solution is rejected if nonphysical.
        u = np.maximum(y[0],1e-6)
        t = np.maximum(y[1]*tb,20.)
        length = np.exp(np.clip(p[2],-10,6))
        j = np.exp(np.clip(p[:2],-250,0))
        energy = p[3]*scale+q*heating_cdf(s,cfg.heating_shape)
        advected = np.sum(j[:,None]*(3.5*RS_SPECIES[:,None]*t-GM*u/R),axis=0)
        kappa = cfg.conductivity_coefficient*cfg.conductivity_multiplier*t
        return np.vstack([-R*rs*t*length/GM,
            -R**2*(energy-advected)*rs*t*length/(GM*kappa*tb)])

    def boundaries(a,b,p):
        u=max(b[0],1e-6); t=max(b[1]*tb,20.)
        length=np.exp(np.clip(p[2],-10,6)); j=np.exp(np.clip(p[:2],-250,0))
        lje,ei,_=jeans_boundary(u,t,math.log(pb)-length,fractions)
        return np.r_[a[0]-ub,a[1]-1,(log_kn_const+length+2*math.log(u))/10,
                     (p[:2]-lje)/10,(p[3]*scale+q-j@ei)/scale]

    sol=solve_bvp(rhs,boundaries,mesh,init,p=par,tol=cfg.tolerance,max_nodes=cfg.max_nodes)
    sol.heat_scale=scale
    if not sol.success: raise RuntimeError(f'BVP failed: {sol.message}')
    if np.any(sol.y[0]<=0) or np.any(sol.y[1]<=0): raise RuntimeError('Nonphysical BVP iterate')
    grid=np.linspace(0,1,801)
    u,tn=sol.sol(grid); t=tn*tb; length=math.exp(sol.p[2]); j=np.exp(sol.p[:2])
    p=pb*np.exp(-length*grid); r=R/u
    lje,ei,lam=jeans_boundary(u[-1],t[-1],math.log(p[-1]),fractions)
    energy=sol.p[3]*scale+q*heating_cdf(grid,cfg.heating_shape)
    adv=np.sum(j[:,None]*(3.5*RS_SPECIES[:,None]*t-GM*u/R),axis=0)
    conductive=energy-adv
    density=p/(rs*t)
    speed=j.sum()*u*u/density
    sound=np.sqrt(1.4*rs*t)
    kn=GM*(mean_molar/NA)/(cfg.collision_cross_section_m2*p*r*r)
    first_kn=int(np.argmax(kn>=1-1e-5)) if np.any(kn>=1-1e-5) else len(kn)-1
    flags=['MOLECULAR_LIMIT_ONLY','PRESCRIBED_LOWER_ATMOSPHERE','NO_IR_COOLING',
           'NET_HEATING_INPUT_NOT_SPECTRAL_PREDICTION','FIXED_MIXING_RATIO']
    if min(lam)<15: flags.append('KINETIC_ESCAPE_SENSITIVITY_REQUIRED')
    if min(lam)<3: flags.append('HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED')
    if max(speed/sound)>.1: flags.append('INERTIA_MAY_INVALIDATE_HYDROSTATIC')
    if max(r)>0.1*61_500_000: flags.append('EARTH_TIDES_REQUIRE_REVIEW')
    if first_kn < len(kn)-3: flags.append('EARLIER_EXOBASE_CROSSING')
    summary=dict(deposited_heat_w_m2=q,lower_temperature_k=tb,lower_pressure_pa=pb,
        lower_radius_R=1/ub,oxygen_mole_fraction=cfg.oxygen_mole_fraction,
        conductivity_multiplier=cfg.conductivity_multiplier,heating_shape=cfg.heating_shape,
        exobase_temperature_k=float(t[-1]),max_temperature_k=float(t.max()),
        exobase_radius_R=float(1/u[-1]),exobase_pressure_pa=float(p[-1]),
        jeans_lambda_N2=float(lam[0]),jeans_lambda_O2=float(lam[1]),
        N2_loss_kg_s=float(j[0]*AREA),O2_loss_kg_s=float(j[1]*AREA),
        molecular_loss_kg_s=float(j.sum()*AREA),max_mach=float(max(speed/sound)),
        lower_conductive_flux_W_m2=float(conductive[0]),
        upper_energy_to_infinity_W_m2=float(j@ei),
        net_energy_boundary_residual_W_m2=float(energy[-1]-j@ei),
        max_scaled_BC_residual=float(np.max(abs(boundaries(sol.y[:,0],sol.y[:,-1],sol.p)))),
        max_collocation_residual=float(max(sol.rms_residuals)),mesh_nodes=int(len(sol.x)),
        domain_flags=flags)
    return summary,dict(s=grid,radius_R=r/R,temperature_K=t,pressure_Pa=p,
        conductive_flux_per_surface_W_m2=conductive,knudsen=kn),sol
