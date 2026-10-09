"""Conditional Sun-following energy and wind-vector requirements.

No route, control organ, aerodynamic shape or accessible wind layer is established.
Power is per horizontal projected collecting area. The airspeed is relative to
the local air, never inferred from ground speed or a scalar wind percentile.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
INPUT_FILES = ('research/studies/sky_ships/results/sky_ships.json',
               'research/studies/floater_viability/navigation_sources.json')


def relative_velocity(desired_east_m_s, desired_north_m_s,
                      wind_east_m_s, wind_north_m_s):
    """Body ground velocity minus local air velocity; east/north positive."""
    values = (desired_east_m_s, desired_north_m_s, wind_east_m_s, wind_north_m_s)
    if not all(math.isfinite(x) for x in values):
        raise ValueError('Finite velocity components required')
    east = desired_east_m_s - wind_east_m_s
    north = desired_north_m_s - wind_north_m_s
    return dict(air_east_m_s=east, air_north_m_s=north,
                airspeed_m_s=math.hypot(east, north))


def drag_power(airspeed_m_s, air_density_kg_m3, drag_coefficient=.47,
               efficiency=.25, drag_area_per_projected_area=1.):
    """Quasisteady drag work; efficiency bundles actuator and propulsive losses.

    Cd uses frontal drag-reference area. Ratio converts to collecting footprint.
    Cd=.47 and .05 are shape sensitivities, not validated giant-colony values.
    No acceleration, added mass, turns, lift generation or altitude-control cost.
    """
    values = (airspeed_m_s, air_density_kg_m3, drag_coefficient,
              efficiency, drag_area_per_projected_area)
    if (not all(math.isfinite(x) for x in values) or airspeed_m_s < 0
        or min(air_density_kg_m3, drag_coefficient, drag_area_per_projected_area) <= 0
        or not 0 < efficiency <= 1):
        raise ValueError('Physical speed, density, drag area and efficiency required')
    drag = .5 * air_density_kg_m3 * drag_coefficient * drag_area_per_projected_area * airspeed_m_s**2
    mechanical = drag * airspeed_m_s
    return dict(airspeed_m_s=airspeed_m_s, drag_n_m2=drag,
                mechanical_w_m2=mechanical, input_w_m2=mechanical/efficiency,
                efficiency=efficiency, drag_coefficient=drag_coefficient,
                drag_area_per_projected_area=drag_area_per_projected_area)


def wind_match_tolerance(available_input_w_m2, air_density_kg_m3,
                         drag_coefficient=.47, efficiency=.25,
                         drag_area_per_projected_area=1.):
    """Largest continuous vector mismatch spending an already-net margin.

    It is a radius in east/north wind-vector space, not a directional prediction.
    Chemical margin at biological efficiency; electrical margin at engineering
    efficiency. Those energy sources must not be equated without conversion.
    """
    if not math.isfinite(available_input_w_m2) or available_input_w_m2 < 0:
        raise ValueError('Nonnegative available input power required')
    coefficient = drag_power(1., air_density_kg_m3, drag_coefficient, efficiency,
                             drag_area_per_projected_area)['input_w_m2']
    return (available_input_w_m2/coefficient)**(1/3)


def motion_budget(available_input_w_m2, airspeed_m_s, air_density_kg_m3,
                  drag_coefficient=.47, efficiency=.25, duty_fraction=1.,
                  drag_area_per_projected_area=1., carbon_energy_j_kg_c=40e6):
    """Annual mean input expenditure and remaining margin after locomotion.

    The default 40 MJ/kgC is 18 MJ/kg dry biomass divided by .45 carbon fraction,
    an accounting conversion, not a measured muscle/substrate yield. Supply a
    chemically consistent alternative if the same ledger uses another substrate.
    Duty is imposed; it does not establish tracking between powered intervals.
    """
    if (not math.isfinite(available_input_w_m2) or not 0 <= duty_fraction <= 1
        or not math.isfinite(carbon_energy_j_kg_c) or carbon_energy_j_kg_c <= 0):
        raise ValueError('Finite margin, valid duty and positive energy basis required')
    power = drag_power(airspeed_m_s, air_density_kg_m3, drag_coefficient,
                       efficiency, drag_area_per_projected_area)
    mean = power['input_w_m2'] * duty_fraction
    return dict(**power, duty_fraction=duty_fraction,
                mean_input_w_m2=mean, available_input_w_m2=available_input_w_m2,
                remaining_input_w_m2=available_input_w_m2-mean,
                margin_nonnegative=available_input_w_m2 >= mean,
                carbon_energy_j_kg_c=carbon_energy_j_kg_c,
                annual_motion_carbon_equivalent_kg_m2=mean*YEAR_S/carbon_energy_j_kg_c)


def solar_follow_case(latitude_deg, altitude_m, air_density_kg_m3,
                      wind_east_m_s, wind_north_m_s=0., available_input_w_m2=1.,
                      drag_coefficient=.47, efficiency=.25,
                      drag_area_per_projected_area=1.):
    """Follow mean solar longitude at constant latitude in an imposed wind.

    The photoperiod module owns the synodic geometry. This is not an ephemeris,
    path integrator, optical illumination guarantee or meteorological forecast.
    """
    from research.studies.floater_viability.photoperiod import sun_follow_speed
    west = sun_follow_speed(latitude_deg, altitude_m)
    relative = relative_velocity(-west, 0., wind_east_m_s, wind_north_m_s)
    return dict(latitude_deg=latitude_deg, altitude_m=altitude_m,
                desired_ground_east_m_s=-west, desired_ground_north_m_s=0.,
                wind_east_m_s=wind_east_m_s, wind_north_m_s=wind_north_m_s,
                air_east_m_s=relative['air_east_m_s'],
                air_north_m_s=relative['air_north_m_s'],
                **motion_budget(available_input_w_m2, relative['airspeed_m_s'],
                                air_density_kg_m3, drag_coefficient, efficiency,
                                drag_area_per_projected_area=drag_area_per_projected_area))


def two_layer_advection(desired_east_m_s, desired_north_m_s,
                         low_east_m_s, low_north_m_s,
                         high_east_m_s, high_north_m_s):
    """Nearest time-averaged wind in the line segment joining two wind layers.

    Optimistic kinematic bound: instant/free altitude transfer and constant local
    vectors. Residual is the mean velocity shortfall, not a trajectory solution
    or propulsive cost. Zero residual does not mean free vertical steering.
    """
    desired = (desired_east_m_s, desired_north_m_s)
    low = (low_east_m_s, low_north_m_s)
    delta = (high_east_m_s-low_east_m_s, high_north_m_s-low_north_m_s)
    if not all(math.isfinite(x) for x in (*desired, *low, *delta)):
        raise ValueError('Finite wind vectors required')
    norm2 = sum(x*x for x in delta)
    high_fraction = (sum((t-l)*d for t,l,d in zip(desired,low,delta))/norm2
                     if norm2 else 0.)
    high_fraction = min(1., max(0., high_fraction))
    mean = tuple(l+high_fraction*d for l,d in zip(low,delta))
    residual = relative_velocity(*desired, *mean)
    return dict(high_layer_time_fraction=high_fraction,
                low_layer_time_fraction=1-high_fraction,
                mean_wind_east_m_s=mean[0], mean_wind_north_m_s=mean[1],
                minimum_mean_vector_mismatch_m_s=residual['airspeed_m_s'],
                residual_east_m_s=residual['air_east_m_s'],
                residual_north_m_s=residual['air_north_m_s'],
                vertical_control_cost_accounted=False)


def intermittent_correction(mean_velocity_error_m_s, duty_fraction,
                             air_density_kg_m3, drag_coefficient=.47,
                             efficiency=.25):
    """Ideal pulsed correction of constant mean mismatch in constant wind.

    Coasting is comoving; instantaneous acceleration/deceleration ignored. To
    supply mean correction u at duty d needs active airspeed u/d. Mean cubic drag
    work is continuous drag work/d². This idealization does not grant Sun lock
    during coast periods and omits real startup costs.
    """
    if not 0 < duty_fraction <= 1:
        raise ValueError('Correction duty must be in (0,1]')
    power = drag_power(mean_velocity_error_m_s/duty_fraction,
                       air_density_kg_m3, drag_coefficient, efficiency)
    return dict(mean_velocity_error_m_s=mean_velocity_error_m_s,
                active_airspeed_m_s=power['airspeed_m_s'],
                duty_fraction=duty_fraction,
                mean_input_w_m2=power['input_w_m2']*duty_fraction,
                multiplier_over_continuous=1/duty_fraction**2)


def sunlight_advantage(extra_fixation_w_m2, motion_input_w_m2,
                       extra_other_cost_w_m2=0.):
    """Net extra chemical production from a stated, already-achieved light gain.

    Fixation and costs must share chemical energy units and averaging interval.
    No conversion of incident PAR to fixation is supplied by this function.
    Extra maintenance, gas renewal, cooling/water collection can be passed as
    other costs. Existing maintenance is paid in the baseline budget once.
    """
    if (not all(math.isfinite(x) for x in (extra_fixation_w_m2,
                                          motion_input_w_m2, extra_other_cost_w_m2))
        or min(motion_input_w_m2, extra_other_cost_w_m2) < 0):
        raise ValueError('Finite advantage and nonnegative extra costs required')
    net = extra_fixation_w_m2-motion_input_w_m2-extra_other_cost_w_m2
    return dict(extra_fixation_w_m2=extra_fixation_w_m2,
                motion_input_w_m2=motion_input_w_m2,
                extra_other_cost_w_m2=extra_other_cost_w_m2,
                net_advantage_w_m2=net, energetically_advantageous=net>0)


def evaluate():
    """Small algebraic sensitivities on the inherited atmosphere; no simulation."""
    from research.studies.floater_viability.photoperiod import sun_follow_speed
    ships = json.loads((ROOT/INPUT_FILES[0]).read_text())
    if ships.get('schema') != 'terluna.research.sky-ships/1':
        raise ValueError('Unexpected sky-ships schema')
    ten = next(a for a in ships['air'] if a['height_km']==10)
    rho = ten['density_kg_m3']
    speed = sun_follow_speed(0., 10000.)
    shapes = (('sphere_sensitivity',.47,.25), ('streamlined_hypothesis',.05,.25),
              ('sphere_engineering_input',.47,.7), ('streamlined_engineering_input',.05,.7))
    return dict(
        schema='terluna.research.floater-navigation/1',
        evidence='Conditional vector kinematics and quasisteady propulsion requirements. No accessible wind corridor, propulsor, biological navigation system or solar-locked trajectory is demonstrated.',
        reading_rule='W per horizontal projected collecting m². Biological cases use chemical input; engineering comparator uses electrical input. Cd references frontal area. Duty cases alone do not establish Sun tracking.',
        input_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()[:16] for p in INPUT_FILES},
        air_10km=ten,
        source_assumptions=dict(drag_area_per_projected_area=1.,
            biological_overall_efficiency=.25, engineering_overall_efficiency=.7,
            carbon_energy_j_kg_c=40e6,
            reference_gpp_c_kg_m2_year=2.,
            reference_gpp_chemical_equivalent_w_m2=2.*40e6/YEAR_S),
        calm_cases=[dict(shape=name, **solar_follow_case(lat,10000.,rho,0.,
                          drag_coefficient=cd,efficiency=eta))
                    for name,cd,eta in shapes for lat in (0.,30.,60.,70.,80.,85.)],
        imposed_vector_cases=[dict(name=name, **solar_follow_case(0.,10000.,rho,east,north))
             for name,east,north in (('exact_westward_match',-speed,0.),
                 ('westward_with_0.5m_s_shortfall',-speed+.5,0.),
                 ('right_zonal_speed_1m_s_crosswind',-speed,1.),
                 ('equal_eastward_speed',speed,0.))],
        relative_speed_cases=[dict(shape=name, **motion_budget(1.,v,rho,cd,eta))
                              for name,cd,eta in shapes for v in (0.,.25,.5,1.,2.,4.,5.)],
        wind_match_tolerances=[dict(shape=name,available_input_w_m2=q,
                 maximum_continuous_mismatch_m_s=wind_match_tolerance(q,rho,cd,eta))
                 for name,cd,eta in shapes for q in (.1,.5,1.,2.)],
        illustrative_layers=[dict(name='westward_two_layer_match',
                    **two_layer_advection(-speed,0.,-2.,0.,-6.,0.)),
                 dict(name='crosswind_uncorrected',
                    **two_layer_advection(-speed,0.,-2.,1.,-6.,1.))],
        occasional_motion=[motion_budget(1.,v,rho,duty_fraction=d)
                           for v in (1.,2.,5.,10.) for d in (.01,.05,.1)],
        pulse_tracking=[intermittent_correction(.5,d,rho) for d in (1.,.5,.1)],
        sunlight_trade=[sunlight_advantage(gain,drag_power(v,rho)['input_w_m2'])
                        for gain in (.1,.5,1.,2.) for v in (0.,.25,.5,1.,2.)],
        unresolved=['Time-resolved east/north winds and vertical shear along reachable routes; scalar inherited wind percentiles cannot resolve them.',
            'Propulsor/actuator and sensory tissue mass, power, durability, control authority and useful efficiency.',
            'Altitude-transfer energy, mass/trim changes, water supply, freeze/heat stress and storm/collision escape.',
            'Joint time-resolved photon, carbon, hydrogen and water budgets under moving local solar time.',
            'Finite daytime corridor, weather avoidance and biological rhythm need not require exact solar longitude.'])
