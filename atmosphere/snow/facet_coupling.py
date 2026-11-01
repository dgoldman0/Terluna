"""Implicit spheroidal surface humidity, attachment and latent-heat balance.

Solves Zhang & Harrington (2014), Eqs. 9, 11, 13--16, with the
Harrington et al. (2021) characteristic supersaturations and facet axis law.
Ventilation multiplies continuum conductance only. Supplying those factors is
a separate physical model; no gravity, drag or habit correlation is implicit.
The earlier semianalytic prototype and external source bytes remain unchanged.
"""
from __future__ import annotations

import math
from scipy.optimize import brentq


def _positive(**values):
    for name, value in values.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be positive and finite')


def efficiency(surface_supersaturation, characteristic, mechanism):
    """Nelson--Baker alpha without a fitted floor; stable over small/large ratios."""
    _positive(characteristic=characteristic, mechanism=mechanism)
    if not 1 <= mechanism <= 15:
        raise ValueError('Tested mechanism interval is 1 <= m <= 15')
    if not math.isfinite(surface_supersaturation) or surface_supersaturation < 0:
        raise ValueError('Growth-only surface supersaturation must be finite')
    if surface_supersaturation == 0:
        return 0.
    log_power = mechanism * (math.log(surface_supersaturation) - math.log(characteristic))
    if log_power < -4:
        # tanh(exp(-log_power)) rounds to unity here.
        return math.exp(log_power)
    inverse = math.exp(-log_power)
    if inverse < 1e-4:
        return 1 - inverse * inverse / 3
    return math.tanh(inverse) / inverse


def surface_state(driving, geometric_resistance, characteristic, mechanism, tolerance=1e-12):
    """Unique root of s_surface[1+L alpha(s_surface)] = driving.

    Bracket in log supersaturation between the perfect-attachment and
    zero-attachment limits. This avoids a source-imposed minimum alpha.
    """
    _positive(geometric_resistance=geometric_resistance, characteristic=characteristic,
              tolerance=tolerance)
    if not math.isfinite(driving) or driving < 0:
        raise ValueError('Invalid growth driving supersaturation')
    if driving == 0:
        return 0., 0.
    log_high = math.log(driving)
    log_low = -math.log1p(geometric_resistance)

    def residual(log_fraction):
        alpha = efficiency(math.exp(log_high + log_fraction), characteristic, mechanism)
        return log_fraction + math.log1p(alpha * geometric_resistance)

    root = brentq(residual, log_low, 0., xtol=tolerance, rtol=tolerance)
    supersaturation = math.exp(log_high + root)
    return supersaturation, efficiency(supersaturation, characteristic, mechanism)


def solve(state, a, c, temperature, ice_supersaturation, mechanism=10., *,
          vapor_ventilation=1., heat_ventilation=1., tolerance=1e-12):
    """Solve both facet humidities and the common linearized thermal balance.

    Source conventions (air properties, alpha_T and critical supersaturations)
    are explicit inputs. Supersaturations use the source's linearized thermal
    convention; a resolved flow field and nonlinear surface temperature are
    outside this calculation.
    """
    _positive(a=a, c=c, temperature=temperature, vapor_ventilation=vapor_ventilation,
              heat_ventilation=heat_ventilation, tolerance=tolerance)
    if not 1 <= mechanism <= 15 or not math.isfinite(mechanism):
        raise ValueError('Tested mechanism interval is 1 <= m <= 15')
    if not math.isfinite(ice_supersaturation) or ice_supersaturation < 0:
        raise ValueError('Growth-only ambient supersaturation must be finite')
    fields = ('capacitance_m', 'boundary_capacitance_m', 'diffusivity_m2_s',
              'conductivity_W_m_K', 'vapor_molecular_speed_m_s', 'thermal_jump_length_m',
              'vapor_gas_constant_J_kg_K', 'ice_saturation_pressure_Pa',
              'latent_sublimation_J_kg', 'characteristic_supersaturation_a',
              'characteristic_supersaturation_c')
    _positive(**{name: state[name] for name in fields})
    cap = state['capacitance_m']
    boundary = state['boundary_capacitance_m']
    if boundary < cap:
        raise ValueError('Enlarging the crystal cannot reduce capacitance')
    dv, ka = state['diffusivity_m2_s'], state['conductivity_W_m_K']
    speed = state['vapor_molecular_speed_m_s']
    rv, latent = state['vapor_gas_constant_J_kg_K'], state['latent_sublimation_J_kg']
    rhoeq = state['ice_saturation_pressure_Pa'] / (rv * temperature)
    thermal = rhoeq * latent / temperature * (latent / (rv * temperature) - 1)
    if thermal < 0:
        raise ValueError('Negative latent-heat resistance outside this ice-growth model')
    areas = (a * c, a * a)
    critical = (state['characteristic_supersaturation_a'], state['characteristic_supersaturation_c'])
    ratio = cap / boundary
    # Multiplying continuum conductance leaves molecular attachment/jump
    # resistances untouched, including at large prescribed ventilation.
    vapor_continuum = ratio / vapor_ventilation
    heat_continuum = ratio / heat_ventilation
    kinetic = [4 * dv * cap / (speed * area) for area in areas]
    lengths = [vapor_continuum / beta for beta in kinetic]
    conductivities = [ka / (heat_continuum + state['thermal_jump_length_m'] * cap / area)
                      for area in areas]
    kstar = (2 * conductivities[0] + conductivities[1]) / 3

    def at_driving(driving):
        facets = [surface_state(driving, length, scrit, mechanism, tolerance)
                  for length, scrit in zip(lengths, critical)]
        alpha = [item[1] for item in facets]
        dx = [aa * dv / (aa * vapor_continuum + beta) for aa, beta in zip(alpha, kinetic)]
        dstar = (2 * dx[0] + dx[1]) / 3
        return facets, dx, dstar

    if ice_supersaturation == 0:
        fraction = 1.
    else:
        # Normalize this root to [0,1]; tiny ambient supersaturation remains
        # resolved without a dimensional absolute root tolerance.
        def heat_residual(fraction):
            _, _, dstar = at_driving(fraction * ice_supersaturation)
            return fraction * (1 + thermal * dstar / kstar) - 1
        fraction = brentq(heat_residual, 0., 1., xtol=tolerance, rtol=tolerance)
    driving = fraction * ice_supersaturation
    facets, dx, dstar = at_driving(driving)
    dm = 4 * math.pi * cap * rhoeq * dstar * driving
    delta_t = latent * dm / (4 * math.pi * cap * kstar)
    beta_t = (latent / (rv * temperature) - 1) / temperature
    residual_scale = max(ice_supersaturation, 1e-300)
    facet_residuals = [(ss * (1 + aa * length) - driving) / residual_scale
                       for (ss, aa), length in zip(facets, lengths)]
    heat_residual = (driving + beta_t * delta_t - ice_supersaturation) / residual_scale
    return dict(alpha_a=facets[0][1], alpha_c=facets[1][1],
        surface_supersaturation_a=facets[0][0], surface_supersaturation_c=facets[1][0],
        thermally_reduced_supersaturation=driving, latent_temperature_increment_K=delta_t,
        effective_diffusivity_m2_s=dstar, effective_conductivity_W_m_K=kstar,
        axis_diffusivities_m2_s=dx, mass_rate_kg_s=dm,
        relative_surface_balance_residuals=facet_residuals,
        relative_thermal_balance_residual=heat_residual,
        vapor_ventilation=vapor_ventilation, heat_ventilation=heat_ventilation,
        root_tolerance=tolerance)


def axis_rates(solution, a, c, deposition_density):
    """2021 facet dc/da=alpha_c/alpha_a with an explicit deposition density."""
    _positive(a=a, c=c, deposition_density=deposition_density)
    dm = solution['mass_rate_kg_s']
    if dm == 0:
        return dict(equatorial_rate_m_s=0., polar_rate_m_s=0., facet_rate_ratio=None)
    aa, ac = solution['alpha_a'], solution['alpha_c']
    _positive(alpha_a=aa, alpha_c=ac, mass_rate=dm)
    # Avoid division by a vanishing alpha while preserving the facet ratio.
    denominator = deposition_density * (4 * math.pi / 3) * (2 * a * c * aa + a * a * ac)
    da, dc = dm * aa / denominator, dm * ac / denominator
    return dict(equatorial_rate_m_s=da, polar_rate_m_s=dc, facet_rate_ratio=ac / aa)
