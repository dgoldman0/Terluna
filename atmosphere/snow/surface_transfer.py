"""Conditional adjoining-facet capture on a prescribed hollow cylinder.

Nelson & Swanson (2019), Appendix B, motivates the mechanism and its length
scale. The finite-volume capture/redistribution below is our own diagnostic,
not their model or a measured kinetic closure. It assumes an absorbing rough
rim, exponential survival along adjoining prism surfaces, and supplied edge
transmission. The bulk field and total mass respond to the added uptake.
"""
from __future__ import annotations
import math
import numpy as np
from atmosphere.snow.flakes import T0
from atmosphere.snow.spatial_crystal_growth import MaskDomain
from atmosphere.snow.facet_growth import ICE_DENSITY
from shared.constants import BOLTZMANN, GAS_CONSTANT


def capture_average(distance_near, distance_far, migration_length, transmission=1.):
    """Cell average of transmission*exp(-distance/xs), with physical lengths.

    This physical collection distance has no grid-pixel multiplier. It is an
    assumed interfacet survival law, distinct from the 2015 vicinal-step rule.
    """
    lo, hi = np.broadcast_arrays(np.asarray(distance_near, float), np.asarray(distance_far, float))
    if np.any(~np.isfinite(lo)) or np.any(~np.isfinite(hi)) or np.any(lo < 0) or np.any(hi <= lo):
        raise ValueError('Finite ordered nonnegative distance intervals required')
    if not math.isfinite(migration_length) or migration_length < 0 or not 0 <= transmission <= 1:
        raise ValueError('Nonnegative migration length and transmission in [0,1] required')
    if migration_length == 0 or transmission == 0:
        return np.zeros_like(lo)
    width = hi-lo
    return transmission*np.exp(-lo/migration_length)*(-np.expm1(-width/migration_length))/(width/migration_length)


class CaptureBoundary:
    """Partition incident uptake between local growth and transfer to the rim."""
    def __init__(self, law, capture_probability, rough, recipient_offset=0.):
        self.law = law
        self.capture = np.asarray(capture_probability, float)
        self.rough = np.asarray(rough, bool)
        self.recipient_offset = recipient_offset
        if not math.isfinite(recipient_offset) or recipient_offset < 0:
            raise ValueError('Nonnegative recipient equilibrium shift required')
        if np.any(~np.isfinite(self.capture)) or np.any(self.capture < 0) or np.any(self.capture > 1):
            raise ValueError('Capture probabilities must lie in [0,1]')

    def flux(self, s):
        s = np.asarray(s, float)
        p, derivative = self.law.flux(s)
        alpha = np.divide(p, s, out=derivative.copy(), where=s > 0)
        dalpha = np.divide(derivative-alpha, s, out=np.zeros_like(s), where=s > 0)
        drive = np.maximum(s-self.recipient_offset, 0.)
        positive = s > self.recipient_offset
        # Capture consumes P*(s-offset)+ of the incident molecular excess;
        # uncaptured excess retains the local attachment probability alpha.
        # At offset=0 this is the independent-probability uptake rule.
        total = p+self.capture*(1-alpha)*drive
        slope = derivative+self.capture*((1-alpha)*positive-dalpha*drive)
        return np.where(self.rough, drive, total), np.where(self.rough, positive.astype(float), slope)


def transfer_rates(domain, air, radius, temperature, supersaturation, laws,
                   migration_length, transmission=1., inner_donor=True,
                   ventilation=1., heating=True, tolerance=1e-9, capillary_length=0.):
    """Frozen-profile rates; supplied rim/cavity and absorbing edge are hypotheses.

    Donor uptake leaves the gas at its actual location and is deposited over the
    top annulus. Isothermal heat uses the full absorbed mass. This conserves
    water but does not solve surface diffusion, edge kinetics or profile motion.
    """
    if not all(math.isfinite(x) and x > 0 for x in (radius, temperature, ventilation)):
        raise ValueError('Positive physical radius, temperature and ventilation required')
    if not math.isfinite(capillary_length) or capillary_length < 0:
        raise ValueError('Nonnegative capillary length required')
    if domain.shape != 'cylinder':
        raise ValueError('Prescribed cylindrical geometry required')
    if not 0 < domain.cavity_radius < 1:
        if capillary_length:
            raise ValueError('Curvature sensitivity is defined only for a hollow cylindrical rim')
    if len(laws) != 2:
        raise ValueError('Prism and basal local laws required')
    d = MaskDomain(domain.rf, domain.zf, domain.solid)
    orientation = d.ice[:, 3].astype(int)
    phi = domain.aspect
    tip = (orientation == 1) & np.isclose(d.ice[:, 5], phi, rtol=0, atol=1e-12)
    outer = (orientation == 0) & np.isclose(d.ice[:, 4], 1., rtol=0, atol=1e-12)
    donor = (orientation == 0) if inner_donor else outer
    probability = np.zeros(len(d.ice))
    j = np.searchsorted(d.z, d.ice[donor, 5])
    probability[donor] = capture_average((phi-d.zf[j+1])*radius,
        (phi-d.zf[j])*radius, migration_length, transmission)
    # Supplied outer-terrace proxy, kappa=2/rim_width. Both direct rim
    # uptake and captured donor excess refer to the same recipient equilibrium.
    # This growth-only sensitivity omits sublimation and a surface-field solve.
    recipient_offset = 2*capillary_length/(radius*(1-domain.cavity_radius))
    combined = [CaptureBoundary(law, probability[orientation == k], tip[orientation == k], recipient_offset)
        for k, law in enumerate(laws)]
    rhoeq = air['ice_saturation_pressure_Pa']/(air['vapor_gas_constant_J_kg_K']*temperature)
    vth = air['vapor_molecular_speed_m_s']/4
    D = air['diffusivity_m2_s']*ventilation
    K = air['conductivity_W_m_K']*ventilation
    L = air['latent_sublimation_J_kg']
    beta = (L/(air['vapor_gas_constant_J_kg_K']*temperature)-1)/temperature
    heat = beta*L*rhoeq*D/K if heating else 0.
    lam = radius*vth/D
    answer = d.coupled(supersaturation, lam, combined, heat_parameter=heat, tolerance=tolerance)
    area = d.ice[:, 1]
    uptake = answer['face_flux']
    captured = lam*probability*np.maximum(answer['surface']-recipient_offset, 0.)
    # The capture term is an actual part of the solved gas-side flux.
    local = uptake-captured
    if np.min(local) < -1e-10*max(np.max(uptake), 1e-30):
        raise ArithmeticError('Capture exceeded the solved vapour supply')
    local = np.maximum(local, 0.)
    captured_total = float(np.dot(area, captured))
    tip_area = float(np.sum(area[tip]))
    deposition = local.copy()
    deposition[tip] += captured_total/tip_area
    deposited_total = float(np.dot(area, deposition))
    received = float(np.dot(area[tip], deposition[tip]-local[tip]))
    Q = answer['half_flux']
    scale = rhoeq*D/(radius*ICE_DENSITY)
    outer_area = float(np.sum(area[outer]))
    tip_rate = scale*float(np.dot(area[tip], deposition[tip]))/tip_area
    prism_rate = scale*float(np.dot(area[outer], deposition[outer]))/outer_area
    return dict(mass_rate_kg_s=2*rhoeq*D*radius*Q,
        tip_mean_normal_rate_m_s=tip_rate, outer_prism_mean_normal_rate_m_s=prism_rate,
        outer_prism_max_normal_rate_m_s=scale*float(np.max(deposition[outer])),
        outer_prism_min_normal_rate_m_s=scale*float(np.min(deposition[outer])),
        tip_min_normal_rate_m_s=scale*float(np.min(deposition[tip])),
        tip_max_normal_rate_m_s=scale*float(np.max(deposition[tip])),
        captured_fraction_of_total=captured_total/Q if Q else 0.,
        tip_fraction_of_deposition=float(np.dot(area[tip], deposition[tip]))/Q if Q else 0.,
        local_outer_prism_fraction=float(np.dot(area[outer], local[outer]))/Q if Q else 0.,
        transfer_mass_residual=abs(deposited_total-Q)/max(Q, 1e-30),
        transfer_recipient_residual=abs(received-captured_total)/max(Q, 1e-30),
        gas_water_residual=answer['relative_mass_residual'],
        nonlinear_residual=answer['relative_nonlinear_residual'], heat_shift=answer['heat_shift'],
        capture_area_m2=2*radius**2*float(np.dot(area, probability)),
        tip_area_m2=2*radius**2*tip_area,
        minimum_surface_supersaturation=float(np.min(answer['surface'])),
        recipient_equilibrium_offset=recipient_offset, requested_nonlinear_tolerance=tolerance,
        cells=d.cells, faces=len(d.ice))


def prism_spiral_rate(temperature, supersaturation, ledge_energy, activation_energy=44e3):
    """Miyamoto et al. (2022) low-supersaturation prism spiral diagnostic.

    Parameters in SI. Reported source values: a=0.4 nm, attempt frequency
    4e13/s, molecular area 8.30e-20 m2, step height 0.39 nm, BCF spacing=19*rc.
    Kappa is supplied within the measured whole-temperature envelope; no
    temperature dependence is inferred by interpolating that envelope.
    """
    if not T0-25 <= temperature <= T0-2.6:
        raise ValueError('Temperature outside the reported -25 to -2.6 C interval')
    if not math.isfinite(supersaturation) or supersaturation <= 0:
        raise ValueError('Positive local surface supersaturation required')
    if not math.isfinite(ledge_energy) or ledge_energy <= 0 or not math.isfinite(activation_energy) or activation_energy <= 0:
        raise ValueError('Positive ledge and activation energies required')
    step_coefficient = .4e-9*4e13*math.exp(-activation_energy/(GAS_CONSTANT*temperature))
    step_velocity = step_coefficient*supersaturation
    spacing = 19*8.30e-20*ledge_energy/(BOLTZMANN*temperature*math.log1p(supersaturation))
    return dict(step_coefficient_m_s=step_coefficient, lateral_step_velocity_m_s=step_velocity,
        spacing_m=spacing, normal_rate_m_s=step_velocity*.39e-9/spacing)
