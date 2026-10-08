"""Steady spherical conduction/advection/Jeans column: an explicitly molecular limit.

The lower atmosphere and homopause temperature are boundary inputs. The upper
profile, exobase and molecular losses are solved, with no imposed exobase floor.
Heating is *net deposited power per lunar surface area*, not incident sunlight.
A valid spectrum-to-heating conversion must be supplied separately. This is NOT
a photochemical or kinetic closure: prescribed well-mixed N2/O2, no ions/atoms,
tides, eddy transport or momentum acceleration.

Infrared cooling is optional. Given the mole fractions of CO2, atomic oxygen and
NO against log pressure above the base, and CO2's band escape against its column,
the column radiates CO2's 15-um bending band and NO's 5.3-um band (two-level, out
of local thermodynamic equilibrium) and oxygen's 63- and 147-um fine-structure
lines (in equilibrium), each in excess of what air at the base temperature emits:
the base holds the middle atmosphere's radiative balance, so air left at the base
temperature neither cools nor warms by radiation, as without cooling. A photon
leaves with the mean of its escape probabilities up to space and down to the base.
Without the fractions the column has no infrared cooling.

Equations and domain flags are documented in research/findings.md.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
from scipy.special import expn

from shared.constants import (AVOGADRO as NA, BOLTZMANN as KB, MOON_GM as GM, MOON_RADIUS as R, PLANCK,
                              SPEED_OF_LIGHT)

AREA = 4 * math.pi * R**2
MOLAR = np.array([0.0280134, 0.031998])
RS_SPECIES = KB * NA / MOLAR

# Infrared cooling. hc/k in cm K turns wavenumbers into temperatures.
HC_K = PLANCK * SPEED_OF_LIGHT / KB * 100.0
# CO2's bending fundamental at 667.38 cm^-1 as a two-level band (01101 over the ground state, degeneracy 2), with the
# middle atmosphere's Einstein A and deactivation rates (radiative_convective/ck.NonLTE; Lopez-Puertas and Taylor 2001):
# by N2 and O2 7e-17 T^0.5 + 6.7e-10 exp(-83.8 T^-1/3) cm^3/s, by O 3e-12 (T/300)^0.5, uncertain by about two.
CO2_BAND = dict(theta_k=HC_K * 667.38, einstein_a=1.5, degeneracy=2.0)
# NO's fundamental at 1876 cm^-1 (Kockarts 1980, A 13.3 s^-1), excited by atomic oxygen at 4.2e-11 cm^3/s (Hwang et al.
# 2003); excitation by N2 and O2 is three orders of magnitude slower and left out.
NO_BAND = dict(theta_k=HC_K * 1876.0, einstein_a=13.3, o_rate_cm3_s=4.2e-11)
# Atomic oxygen's ground term: levels 3P2, 3P1, 3P0 (energies cm^-1, degeneracies) and its two magnetic-dipole lines,
# 63.2 um (3P1-3P2, A 8.91e-5 s^-1) and 145.5 um (3P0-3P1, A 1.75e-5 s^-1), NIST. Collisions thermalise the levels
# above about 1e6 cm^-3, nearly to the exobase.
O_LEVELS_CM = np.array([0.0, 158.265, 226.977])
O_DEGENERACY = np.array([5.0, 3.0, 1.0])
O_LINES = ((1, 0, 8.91e-5), (2, 1, 1.75e-5))
O_MASS_KG = 0.0159994 / NA
CO2_MOLAR = 0.0440095

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
    # With heating_shape 'traced': the share of the heat deposited below each log pressure above the base,
    # ln(lower_pressure_pa/p), rising from 0 at the base to 1 (a traced deposition profile, renormalised to the heat
    # that lies below the exobase).
    heating_log_pressure: tuple = ()
    heating_fraction: tuple = ()
    # Infrared cooling (see the module's docstring): mole fractions of CO2, O and NO on log pressures above the base,
    # ln(lower_pressure_pa/p), rising from 0; and CO2's band escape, the share of a photon's chances to leave one way
    # through a CO2 column (cm^-2) at the base temperature. Empty: no infrared cooling. The multiplier scales it.
    coolant_log_pressure: tuple = ()
    co2_mole_fraction: tuple = ()
    o_mole_fraction: tuple = ()
    no_mole_fraction: tuple = ()
    co2_escape_column_cm2: tuple = ()
    co2_escape: tuple = ()
    cooling_multiplier: float = 1.0
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
        if self.heating_shape not in ('middle', 'low', 'high', 'traced'):
            raise ValueError('Unknown deposition shape')
        if self.heating_shape == 'traced':
            x = np.asarray(self.heating_log_pressure, float)
            f = np.asarray(self.heating_fraction, float)
            if (x.size < 2 or x.size != f.size or x[0] != 0 or np.any(np.diff(x) <= 0) or f[0] != 0
                    or np.any(np.diff(f) < 0) or abs(f[-1] - 1) > 1e-9):
                raise ValueError('A traced heating shape needs log pressures rising from 0 and shares from 0 to 1')
        if self.coolant_log_pressure:
            x = np.asarray(self.coolant_log_pressure, float)
            fractions = [np.asarray(getattr(self, k), float)
                         for k in ('co2_mole_fraction', 'o_mole_fraction', 'no_mole_fraction')]
            n, e = np.asarray(self.co2_escape_column_cm2, float), np.asarray(self.co2_escape, float)
            if (x.size < 2 or x[0] != 0 or np.any(np.diff(x) <= 0)
                    or any(f.size != x.size or np.any(f < 0) or np.any(f >= 1) for f in fractions)):
                raise ValueError('Infrared cooling needs mole fractions on log pressures rising from 0')
            if (n.size < 2 or n.size != e.size or n[0] <= 0 or np.any(np.diff(n) <= 0) or np.any(e <= 0)
                    or np.any(e > 1) or np.any(np.diff(e) > 0)):
                raise ValueError('CO2 band escape must fall from at most 1 over rising columns')
            if not (math.isfinite(self.cooling_multiplier) and self.cooling_multiplier >= 0):
                raise ValueError('The cooling multiplier must be nonnegative')
        if self.mesh_points < 21 or self.max_nodes < self.mesh_points:
            raise ValueError('Invalid BVP grid')


def heating_cdf(s, shape):
    s = np.asarray(s)
    if shape == 'middle': return 3*s*s - 2*s**3
    if shape == 'low': return 1-(1-s)**3
    if shape == 'high': return s**3
    raise ValueError('Unknown heating shape')


def deposited_share(s, cfg: ColumnConfig, length):
    """Share of the heat deposited below each normalised log pressure s of a column whose exobase lies `length`
    e-folds of pressure above its base."""
    if cfg.heating_shape != 'traced':
        return heating_cdf(s, cfg.heating_shape)
    x, f = np.asarray(cfg.heating_log_pressure), np.asarray(cfg.heating_fraction)
    top = float(np.interp(length, x, f))
    if top <= 0:
        raise ValueError('The traced heat lies above the exobase')
    return np.interp(np.asarray(s) * length, x, f) / top


def doppler_escape(tau):
    """Chance that a photon of a Doppler-broadened line, emitted into the hemisphere toward a boundary at line-centre
    optical depth tau (vertical), leaves through it: the integral over the profile of E2(tau phi(x)/phi(0))."""
    tau = np.asarray(tau, float)
    x = np.linspace(0.0, 6.0, 241)
    weight = np.exp(-x * x) / math.sqrt(math.pi) * 2.0 * np.gradient(x)
    weight[[0, -1]] *= 0.5
    flat = np.maximum(tau.reshape(-1, 1), 0.0) * np.exp(-x * x)
    return (expn(2, flat) @ weight / weight.sum()).reshape(tau.shape)


_ESCAPE_LOG_TAU = np.linspace(-14.0, 14.0, 561)
_ESCAPE_TABLE = []


def doppler_escape_table(tau):
    """doppler_escape interpolated in log tau from a table of it (1 below e^-14, its last value above e^14)."""
    if not _ESCAPE_TABLE:
        _ESCAPE_TABLE.append(doppler_escape(np.exp(_ESCAPE_LOG_TAU)))
    return np.interp(np.log(np.maximum(tau, 1e-300)), _ESCAPE_LOG_TAU, _ESCAPE_TABLE[0], left=1.0)


def _planck_excess(theta, t, tb):
    """1 - B(T_b)/B(T) at a transition of energy theta (K): the share of the emission at T that air at T_b would not
    also give."""
    return 1.0 - np.expm1(theta / t) / np.expm1(theta / tb)


def infrared_cooling(p, t, x, r, cfg: ColumnConfig, mean_molar, parts=False):
    """Net infrared cooling (W/m^3) at pressures p (Pa), temperatures t (K), log pressures x above the base and radii
    r (m); with parts=True, a dict of CO2's, O's and NO's."""
    tb = cfg.lower_temperature_k
    grid = np.asarray(cfg.coolant_log_pressure)
    x_co2, x_o, x_no = (np.interp(x, grid, np.asarray(getattr(cfg, k)))
                        for k in ('co2_mole_fraction', 'o_mole_fraction', 'no_mole_fraction'))
    n = p / (KB * t) * 1e-6                                      # cm^-3
    n_m, n_o = (1.0 - x_co2 - x_o - x_no) * n, x_o * n
    column = p / (mean_molar / NA * GM / r**2) * 1e-4                # molecules cm^-2 above, each way of the base
    below = np.maximum(cfg.lower_pressure_pa / (mean_molar / NA * GM / r**2) * 1e-4 - column, 0.0)
    doppler = np.sqrt(tb / t)
    escape_n, escape_e = np.log(np.asarray(cfg.co2_escape_column_cm2)), np.asarray(cfg.co2_escape)

    def band_escape(c):
        return np.interp(np.log(np.maximum(c * doppler, 1e-30)), escape_n, escape_e, left=1.0)

    # CO2: two-level, n_l h nu A g (e^-theta/T) (1 - B_b/B) * beta eps / (beta + eps), eps the collisional share.
    band = CO2_BAND
    k_m = 7e-17 * np.sqrt(t) + 6.7e-10 * np.exp(-83.8 * t ** (-1.0 / 3.0))
    k_o = 3e-12 * np.sqrt(t / 300.0)
    eps = (k_m * n_m + k_o * n_o) / band['einstein_a']
    beta = 0.5 * (band_escape(x_co2 * column) + band_escape(x_co2 * below))
    co2 = (x_co2 * n * 1e6 * KB * band['theta_k'] * band['einstein_a'] * band['degeneracy']
           * np.exp(-band['theta_k'] / t) * _planck_excess(band['theta_k'], t, tb) * beta * eps / (beta + eps))
    # NO: excited by O, optically thin.
    band = NO_BAND
    eps = band['o_rate_cm3_s'] * n_o / band['einstein_a']
    no = (x_no * n * 1e6 * KB * band['theta_k'] * band['einstein_a'] * np.exp(-band['theta_k'] / t)
          * _planck_excess(band['theta_k'], t, tb) * eps / (1.0 + eps))
    # O: levels in equilibrium at T, each line escaping with its Doppler line's mean chance up and down.
    levels = O_DEGENERACY[:, None] * np.exp(-HC_K * O_LEVELS_CM[:, None] / t)
    share = levels / levels.sum(axis=0)
    o = np.zeros_like(t)
    for upper, lower, a in O_LINES:
        theta = HC_K * (O_LEVELS_CM[upper] - O_LEVELS_CM[lower])
        wavelength_cm = 1.0 / (O_LEVELS_CM[upper] - O_LEVELS_CM[lower])
        width_hz = SPEED_OF_LIGHT * 100.0 / wavelength_cm * np.sqrt(2 * KB * t / O_MASS_KG) / SPEED_OF_LIGHT
        sigma = (share[lower] * O_DEGENERACY[upper] / O_DEGENERACY[lower] * a * wavelength_cm**2 / (8 * math.pi)
                 * -np.expm1(-theta / t) / (math.sqrt(math.pi) * width_hz))
        beta = 0.5 * (doppler_escape_table(sigma * x_o * column) + doppler_escape_table(sigma * x_o * below))
        o += x_o * n * 1e6 * share[upper] * a * KB * theta * _planck_excess(theta, t, tb) * beta
    total = cfg.cooling_multiplier * (co2 + no + o)
    if parts:
        return dict(co2=cfg.cooling_multiplier * co2, o=cfg.cooling_multiplier * o, no=cfg.cooling_multiplier * no)
    return total


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
    cooled = bool(cfg.coolant_log_pressure)
    mesh = np.linspace(0,1,cfg.mesh_points)
    init = np.vstack([ub-c*L0*mesh, np.ones_like(mesh)] + ([np.zeros_like(mesh)] if cooled else []))
    par = np.r_[lj,math.log(L0),(np.exp(lj)@ei-q)/scale]
    if previous is not None and previous.y.shape[0] == init.shape[0]:
        init = previous.sol(mesh)
        par = previous.p.copy()
        par[3] *= previous.heat_scale/scale
        if cooled:
            init[2] *= previous.heat_scale/scale

    def radiated(s,u,t,length):
        """Infrared cooling below each level, per lunar surface area and unit of s (W/m^2)."""
        r = R/u
        return infrared_cooling(pb*np.exp(-length*s),t,length*s,r,cfg,mean_molar)*r**4*rs*t*length/(GM*R**2)

    def rhs(s,y,p):
        # Trial-iterate guards prevent NaN; final solution is rejected if nonphysical.
        u = np.maximum(y[0],1e-6)
        t = np.maximum(y[1]*tb,20.)
        length = np.exp(np.clip(p[2],-10,6))
        j = np.exp(np.clip(p[:2],-250,0))
        energy = p[3]*scale+q*deposited_share(s,cfg,length)-(y[2]*scale if cooled else 0.)
        advected = np.sum(j[:,None]*(3.5*RS_SPECIES[:,None]*t-GM*u/R),axis=0)
        kappa = cfg.conductivity_coefficient*cfg.conductivity_multiplier*t
        out = [-R*rs*t*length/GM, -R**2*(energy-advected)*rs*t*length/(GM*kappa*tb)]
        if cooled:
            out.append(radiated(s,u,t,length)/scale)
        return np.vstack(out)

    def boundaries(a,b,p):
        u=max(b[0],1e-6); t=max(b[1]*tb,20.)
        length=np.exp(np.clip(p[2],-10,6)); j=np.exp(np.clip(p[:2],-250,0))
        lje,ei,_=jeans_boundary(u,t,math.log(pb)-length,fractions)
        lost = b[2]*scale if cooled else 0.
        return np.r_[a[0]-ub,a[1]-1,(log_kn_const+length+2*math.log(u))/10,
                     (p[:2]-lje)/10,(p[3]*scale+q-lost-j@ei)/scale,[a[2]] if cooled else []]

    sol=solve_bvp(rhs,boundaries,mesh,init,p=par,tol=cfg.tolerance,max_nodes=cfg.max_nodes)
    sol.heat_scale=scale
    if not sol.success: raise RuntimeError(f'BVP failed: {sol.message}')
    if np.any(sol.y[0]<=0) or np.any(sol.y[1]<=0): raise RuntimeError('Nonphysical BVP iterate')
    grid=np.linspace(0,1,801)
    state=sol.sol(grid); u,tn=state[:2]; t=tn*tb; length=math.exp(sol.p[2]); j=np.exp(sol.p[:2])
    p=pb*np.exp(-length*grid); r=R/u
    lje,ei,lam=jeans_boundary(u[-1],t[-1],math.log(p[-1]),fractions)
    lost=state[2]*scale if cooled else np.zeros_like(grid)
    energy=sol.p[3]*scale+q*deposited_share(grid,cfg,length)-lost
    adv=np.sum(j[:,None]*(3.5*RS_SPECIES[:,None]*t-GM*u/R),axis=0)
    conductive=energy-adv
    density=p/(rs*t)
    speed=j.sum()*u*u/density
    sound=np.sqrt(1.4*rs*t)
    kn=GM*(mean_molar/NA)/(cfg.collision_cross_section_m2*p*r*r)
    first_kn=int(np.argmax(kn>=1-1e-5)) if np.any(kn>=1-1e-5) else len(kn)-1
    flags=['MOLECULAR_LIMIT_ONLY','PRESCRIBED_LOWER_ATMOSPHERE',
           'IR_COOLING_PRESCRIBED_COMPOSITION' if cooled else 'NO_IR_COOLING',
           'NET_HEATING_INPUT_NOT_SPECTRAL_PREDICTION','FIXED_MIXING_RATIO']
    cooling={}
    if cooled:
        parts=infrared_cooling(p,t,length*grid,r,cfg,mean_molar,parts=True)
        weight=r**4*rs*t*length/(GM*R**2)
        cooling={f'{k}_cooling_W_m2':float(np.trapezoid(v*weight,grid)) for k,v in parts.items()}
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
        **({'infrared_cooling_W_m2': float(lost[-1]), **cooling} if cooled else {}),
        net_energy_boundary_residual_W_m2=float(energy[-1]-j@ei),
        max_scaled_BC_residual=float(np.max(abs(boundaries(sol.y[:,0],sol.y[:,-1],sol.p)))),
        max_collocation_residual=float(max(sol.rms_residuals)),mesh_nodes=int(len(sol.x)),
        domain_flags=flags)
    return summary,dict(s=grid,radius_R=r/R,temperature_K=t,pressure_Pa=p,
        conductive_flux_per_surface_W_m2=conductive,knudsen=kn,
        **({'infrared_cooling_below_W_m2': lost} if cooled else {})),sol
