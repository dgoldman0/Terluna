"""Explicit cylindrical-face transport, with exposed morphology failures.

Own implementation of Nelson & Baker (1996), Eqs. 23--26 and Appendix C,
using the pinned 2026 geometry functions and 2021 air/attachment inputs.
This is a continuum diagnostic, not the unretrieved complete 2026 growth model.
Face-centre checks are necessary conditions for uniform growth, not a resolved
surface field. A thin basal rim is a prescribed asymptotic sensitivity only.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq
from atmosphere.snow.facet_coupling import efficiency, surface_state, _positive
from atmosphere.snow.facet_growth import ICE_DENSITY


def coupled_surfaces(matrix, supersaturation, characteristic, mechanisms, tolerance=1e-13):
    """Solve s_i + sum_j L_ij alpha_j(s_j) s_j = S without alpha floors.

    Eliminate the first surface by a bracketed scalar solve, then bracket the
    second. Nonnegative determinant guarantees a monotone outer residual.
    All supersaturations are relative to ice equilibrium at ambient temperature.
    """
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (2, 2) or not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise ValueError('A finite nonnegative 2 by 2 resistance matrix is required')
    if any(matrix[i, i] <= 0 for i in (0, 1)):
        raise ValueError('Positive self resistances required')
    if matrix[0, 0]*matrix[1, 1] < matrix[0, 1]*matrix[1, 0]*(1-1e-12):
        raise ValueError('Source locations lack the required monotone transport matrix')
    if not math.isfinite(supersaturation) or supersaturation < 0:
        raise ValueError('Finite nonnegative growth supersaturation required')
    _positive(tolerance=tolerance)
    if len(characteristic) != 2 or len(mechanisms) != 2:
        raise ValueError('Two surface characteristics and mechanisms required')
    for critical, mechanism in zip(characteristic, mechanisms):
        efficiency(0., critical, mechanism)  # Validate even at zero forcing.
    if supersaturation == 0:
        return dict(surface_supersaturations=[0., 0.], efficiencies=[0., 0.], residuals=[0., 0.])
    laa, lac = matrix[0]; lca, lcc = matrix[1]
    ca, cc = characteristic; ma, mc = mechanisms

    def first(fraction):
        sc = supersaturation*fraction
        ac = efficiency(sc, cc, mc)
        driving = supersaturation*(1-lac*ac*fraction)
        if driving < -1e-11*supersaturation:
            raise ValueError('Negative remaining vapor supply')
        sa, aa = surface_state(max(0., driving), laa, ca, ma, tolerance)
        return sa, aa, sc, ac

    high = surface_state(supersaturation, lcc, cc, mc, tolerance)[0]/supersaturation
    if lac*efficiency(supersaturation*high, cc, mc)*high > 1:
        high = brentq(lambda x: lac*efficiency(supersaturation*x, cc, mc)*x-1,
                      0., high, xtol=tolerance, rtol=tolerance)

    def residual(fraction):
        sa, aa, sc, ac = first(fraction)
        return (sc*(1+lcc*ac)+lca*aa*sa)/supersaturation-1

    low_value, high_value = residual(0.), residual(high)
    if low_value > 1e-11 or high_value < -1e-11:
        raise ValueError('No nonnegative surface solution for selected step sources')
    # An independently solved diagonal face can put its exact endpoint a few
    # ulps on the wrong side. Accept that endpoint only within the residual
    # bound above; never require a spurious sign change in the decoupled limit.
    if low_value >= 0:
        fraction = 0.
    elif high_value <= 0:
        fraction = high
    else:
        fraction = brentq(residual, 0., high, xtol=tolerance, rtol=tolerance)
    sa, aa, sc, ac = first(fraction)
    ss = np.array([sa, sc]); alpha = np.array([aa, ac])
    errors = (ss + matrix @ (ss*alpha))/supersaturation-1
    if np.max(np.abs(errors)) > 1e-8:
        raise ArithmeticError('Surface balance residual exceeds diagnostic tolerance')
    return dict(surface_supersaturations=ss.tolist(), efficiencies=alpha.tolist(),
                residuals=errors.tolist())


def geometry(shapes, a, c, sources=('edge', 'edge'), thin_basal_rim=None):
    """Cylinder a is radius, c half length; face fluxes are normal mass fluxes."""
    _positive(a=a, c=c)
    if len(sources) != 2 or any(x not in ('edge', 'middle') for x in sources):
        raise ValueError('Specify edge or middle for each of prism and basal faces')
    h = shapes.evaluate(c/a, 'cylinder')
    # Each row is an observation/source position; columns are prism/basal sinks.
    positions = {
        'edge': [a*h['prism_edge'], c*h['basal_edge']],
        'prism_middle': [a*h['prism_middle'], c*h['basal_at_prism_middle']],
        'basal_middle': [a*h['prism_at_basal_middle'], c*h['basal_middle']],
    }
    fraction = 1.
    if thin_basal_rim is not None:
        if tuple(sources) != ('edge', 'edge'):
            raise ValueError('Published thin-rim subtraction is used only at the common edge')
        if not math.isfinite(thin_basal_rim) or not 0 < thin_basal_rim <= .1:
            raise ValueError('Declared thin-rim sensitivity requires 0 < width/a <= 0.1')
        x = 1-thin_basal_rim
        inner = shapes.evaluate(c/(a*x), 'cylinder')
        # NB96 C1 uses c as the basal-flux normalization length on both
        # cylinders. Multiplying the subtracted term by x again is incorrect.
        positions['edge'][1] = c*(h['basal_edge']-inner['basal_edge'])
        fraction = 1-x*x
        if positions['edge'][1] <= 0:
            raise ValueError('Nonpositive thin-ring transport coefficient')
    source_names = [sources[0] if sources[0] == 'edge' else 'prism_middle',
                    sources[1] if sources[1] == 'edge' else 'basal_middle']
    return dict(source_lengths_m=np.array([positions[x] for x in source_names]),
                position_lengths_m=positions, active_face_areas_m2=[4*math.pi*a*c, 2*math.pi*a*a*fraction],
                active_basal_area_fraction=fraction, thin_basal_rim=thin_basal_rim)


def solve(state, shapes, a, c, temperature, ice_supersaturation, *,
          mechanisms=(10., 10.), sources=('edge', 'edge'), thin_basal_rim=None,
          vapor_ventilation=1., heat_ventilation=1., ice_density=ICE_DENSITY, tolerance=1e-13):
    """Local continuum heat/vapor coupling plus independent face kinetics.

    Latent heat follows the same local geometry as vapor, as in NB96; internal
    crystal heat transport and substrate conduction are omitted. Density and
    vapor pressure use the later source's ambient-equilibrium linearization.
    Optional ventilation is an imposed common conductance sensitivity, not a
    calculated flow field. Molecular jumps are omitted and Kn is reported.
    """
    _positive(temperature=temperature, vapor_ventilation=vapor_ventilation,
              heat_ventilation=heat_ventilation, ice_density=ice_density)
    fields = ('diffusivity_m2_s', 'conductivity_W_m_K', 'vapor_molecular_speed_m_s',
              'vapor_gas_constant_J_kg_K', 'ice_saturation_pressure_Pa',
              'latent_sublimation_J_kg', 'vapor_mean_free_path_m',
              'characteristic_supersaturation_a', 'characteristic_supersaturation_c')
    _positive(**{key: state[key] for key in fields})
    g = geometry(shapes, a, c, sources, thin_basal_rim)
    dv = state['diffusivity_m2_s']*vapor_ventilation
    ka = state['conductivity_W_m_K']*heat_ventilation
    kinetic = state['vapor_molecular_speed_m_s']/4
    latent = state['latent_sublimation_J_kg']; rv = state['vapor_gas_constant_J_kg_K']
    rhoeq = state['ice_saturation_pressure_Pa']/(rv*temperature)
    beta = (latent/(rv*temperature)-1)/temperature
    if beta <= 0:
        raise ValueError('Temperature outside the linearized ice-growth regime')
    resistance = kinetic*(1/dv + rhoeq*latent*beta/ka)
    critical = [state['characteristic_supersaturation_a'], state['characteristic_supersaturation_c']]
    solution = coupled_surfaces(g['source_lengths_m']*resistance, ice_supersaturation,
                               critical, mechanisms, tolerance)
    ss = np.array(solution['surface_supersaturations']); alpha = np.array(solution['efficiencies'])
    flux = rhoeq*kinetic*ss*alpha
    areas = np.array(g['active_face_areas_m2'])
    position_states = {}
    for name, lengths in g['position_lengths_m'].items():
        if thin_basal_rim is not None and name != 'edge':
            # The edge asymptote does not provide an off-edge ring field.
            continue
        depletion = float(np.dot(lengths, flux))/(rhoeq*dv)
        dt = latent*float(np.dot(lengths, flux))/ka
        position_states[name] = dict(surface_supersaturation=ice_supersaturation-depletion-beta*dt,
            vapor_depletion=depletion, latent_temperature_increment_K=dt)
    checks = []
    if thin_basal_rim is None:
        for i, name in enumerate(('prism_middle', 'basal_middle')):
            center = position_states[name]['surface_supersaturation']
            required = float(flux[i]/(rhoeq*kinetic*center)) if center > 0 else None
            edge = position_states['edge']['surface_supersaturation']
            required_edge = float(flux[i]/(rhoeq*kinetic*edge)) if edge > 0 else None
            # At zero forcing both the source and centre are saturated.
            possible = ice_supersaturation == 0 or all(
                x is not None and x <= 1+1e-10 for x in (required, required_edge))
            checks.append(dict(face=('prism', 'basal')[i], centre_supersaturation=center,
                required_centre_efficiency=required, required_edge_efficiency=required_edge,
                passes_necessary_uniform_flux_condition=possible,
                centre_exceeds_edge=center > position_states['edge']['surface_supersaturation']+1e-12))
    return dict(surface_supersaturations=solution['surface_supersaturations'],
        efficiencies=solution['efficiencies'], relative_balance_residuals=solution['residuals'],
        normal_mass_fluxes_kg_m2_s=flux.tolist(), equatorial_rate_m_s=float(flux[0]/ice_density),
        polar_rate_m_s=float(flux[1]/ice_density),
        facet_rate_ratio=float(flux[1]/flux[0]) if flux[0] > 0 else None,
        mass_rate_kg_s=float(np.dot(areas, flux)), active_face_areas_m2=areas.tolist(),
        active_basal_area_fraction=g['active_basal_area_fraction'], position_states=position_states,
        uniform_face_checks=checks,
        passes_sampled_uniform_growth_conditions=all(x['passes_necessary_uniform_flux_condition'] for x in checks) if checks else None,
        maximum_mean_free_path_over_semiaxis=state['vapor_mean_free_path_m']/min(a, c),
        mechanisms=list(mechanisms), sources=list(sources), thin_basal_rim=thin_basal_rim,
        vapor_ventilation=vapor_ventilation, heat_ventilation=heat_ventilation,
        geometry='cylinder', root_tolerance=tolerance)
