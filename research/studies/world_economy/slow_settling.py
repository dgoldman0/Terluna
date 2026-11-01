"""Slow-settling process screens and conditional cost headroom.

Ideal dilute transport, empirical correlation extrapolations and assumed costs
are kept separate. No factory yield, market share or world product is inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from numbers import Real

from shared.constants import BOLTZMANN, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import ROOT, constants_used
from .partial_gravity import stokes_speed
from .mature_economics import component_cost

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'slow_settling_inputs.json'


def number(name, value, positive=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f'{name} must be numeric')
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f'{name} must be finite and nonnegative, or positive when required')
    return value


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH, root=ROOT):
    s = json.loads(Path(path).read_text())
    if s['schema'] != 'terluna.research.slow-settling-inputs/1':
        raise ValueError('Unexpected slow-settling input schema')
    products = {}
    for name, spec in s['products'].items():
        p = Path(root) / name
        if digest(p) != spec['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        d = json.loads(p.read_text())
        if d['schema'] != spec['schema']:
            raise ValueError(f'Unexpected product schema: {name}')
        products[name] = d
    return s, products


def terminal_speed(rho_particle, rho_fluid, diameter, viscosity, gravity):
    """Isolated rigid sphere, continuum fluid, Schiller-Naumann for Re<1000.

    Returns magnitude for either settling or creaming. No wall, concentration,
    cohesion, non-Newtonian, rarefaction or mobile-bubble correction is included.
    """
    rp, rf, d, mu = [number(n, v, True) for n, v in zip(
        ('particle density', 'fluid density', 'diameter', 'viscosity'),
        (rho_particle, rho_fluid, diameter, viscosity))]
    g = number('gravity', gravity)
    stokes = stokes_speed(abs(rp-rf), d/2, mu, g)
    lo, hi = 0., stokes
    for _ in range(90):
        v = (lo+hi)/2
        re = rf*v*d/mu
        if v*(1+.15*re**.687) < stokes:
            lo = v
        else:
            hi = v
    v = (lo+hi)/2
    if stokes == 0:
        v = 0.
    re = rf*v*d/mu
    if re >= 1000:
        raise ValueError('Selected drag branch requires particle Re below 1000')
    return dict(speed_m_s=v, stokes_speed_m_s=stokes, particle_reynolds=re,
                stokes_within_five_percent=(v >= .95*stokes),
                stokes_overprediction_fraction=(stokes/v-1 if v else 0.),
                direction='down' if rp > rf else ('up' if rp < rf else 'none'))


def column_screen(case, gravity, temperature, contrast_limit):
    """Zero-flux exponential concentration profile with constant diffusivity.

    K is a prescribed effective vertical diffusivity, not an inferred mixer
    power. Reflecting boundaries exclude deposition/adhesion and throughflow.
    Brownian motion remains at zero gravity. Drift times describe an isolated
    particle, not the settling front of a concentrated or flocculating slurry.
    """
    if case['rho_particle'] < case['rho_fluid']:
        raise ValueError('This column comparison assumes downward sedimentation')
    values = {k:case[k] for k in ('rho_particle', 'rho_fluid', 'diameter', 'viscosity')}
    r = terminal_speed(**values, gravity=gravity)
    t = number('temperature', temperature, True)
    height = number('height', case['height_m'], True)
    drift = number('drift limit', case['drift_limit_m'], True)
    hold = number('hold time', case['hold_s'])
    contrast = number('concentration contrast', contrast_limit, True)
    if contrast <= 1:
        raise ValueError('Concentration contrast must exceed one')
    brownian = BOLTZMANN*t/(3*math.pi*case['viscosity']*case['diameter'])
    speed = r['speed_m_s']
    required = speed*height/math.log(contrast)
    return dict(**r, brownian_diffusivity_m2_s=brownian,
        gravitational_scale_height_m=brownian/speed if speed else None,
        brownian_only_log_bottom_top_contrast=speed*height/brownian,
        brownian_alone_meets_contrast=(brownian >= required),
        required_total_vertical_diffusivity_m2_s=required,
        required_additional_vertical_diffusivity_m2_s=max(0., required-brownian),
        unbounded_hold_drift_m=speed*hold, reaches_boundary_during_hold=(speed*hold >= height),
        hold_within_drift_limit=(speed*hold <= drift),
        time_to_drift_limit_s=drift/speed if speed else None,
        time_to_cross_height_s=height/speed if speed else None)


def suspension_power(gravity_ratio, other_speed_floor, exponent=.45):
    """Extrapolated Njs~g**exponent and turbulent P~N**3 at fixed geometry.

    Floor is an independently imposed minimum speed relative to Earth's Njs.
    Off-bottom suspension is not homogeneous mixing. Constant power number,
    matched rheology and the same feasible flow regime must be established.
    """
    g = number('gravity ratio', gravity_ratio)
    floor = number('other speed floor', other_speed_floor)
    a = number('gravity exponent', exponent, True)
    earth = max(1., floor)
    selected = max(g**a, floor)
    ratio = selected/earth
    return dict(gravity_ratio=g, other_speed_floor_over_earth_Njs=floor,
        gravity_exponent=a, unconstrained_Njs_ratio=g**a,
        operating_speed_ratio=ratio, turbulent_power_ratio=ratio**3,
        same_viscosity_turbulent_shear_scale_ratio=ratio**1.5,
        interpretation='Correlation extrapolation and a speed-floor sensitivity; no empirical lunar mixer or biological yield prediction.')


def batch_cycle(process_s, other_s, earth_settle_s, gravity_ratio):
    process = number('process time', process_s)
    other = number('other time', other_s)
    settle = number('Earth settling time', earth_settle_s, True)
    g = number('gravity ratio', gravity_ratio)
    earth = process+other+settle
    target = process+other+settle/g if g else None
    return dict(passive_settling_s=settle/g if g else None,
        total_cycle_s=target, throughput_relative_to_earth=earth/target if target else 0.,
        additional_critical_path_s=target-earth if target else None,
        interpretation='Stokes-settling harvest on the critical path. Forced filtration, retention or centrifugation is an alternative requiring its own costs; no passive recovery at zero gravity.')


def fluidized_bed(case, gravity, density_ratio, gas_floor_m_s):
    """Ergun incipient-fluidization force balance plus isolated entrainment screen.

    Spherical noncohesive particles, prescribed onset voidage, continuum gas.
    Pressure work includes a specified distributor/loop loss. The transport
    ceiling is a chosen fraction of isolated terminal speed, not a validated
    elutriation law. Expanded bed height/voidage and gas reaction are unsolved.
    """
    g = number('gravity', gravity)
    rr = number('gas density ratio', density_ratio, True)
    rf = number('gas density', case['rho_gas'], True)*rr
    rp = number('particle density', case['rho_particle'], True)
    if rp <= rf:
        raise ValueError('Bed particles must be denser than gas')
    d = number('diameter', case['diameter_m'], True)
    mu = number('viscosity', case['viscosity_Pa_s'], True)
    h = number('height', case['onset_height_m'], True)
    e = number('voidage', case['onset_voidage'], True)
    if e >= 1:
        raise ValueError('Voidage must be below one')
    factor = number('fluidizing margin', case['operating_multiple_of_Umf'], True)
    if factor < 1:
        raise ValueError('Fluidizing multiple must be at least one')
    margin = number('entrainment ceiling', case['terminal_speed_fraction'], True)
    if margin >= 1:
        raise ValueError('Entrainment fraction must be below one')
    other_dp = number('loop pressure loss', case['other_pressure_loss_Pa'])
    floor = number('process gas flow floor', gas_floor_m_s)
    a = 150*mu*(1-e)/(e**3*d**2)
    b = 1.75*rf/(e**3*d)
    weight = (rp-rf)*g
    umf = 2*weight/(a+math.sqrt(a*a+4*b*weight))
    terminal = terminal_speed(rp, rf, d, mu, g)
    required = max(factor*umf, floor)
    ceiling = margin*terminal['speed_m_s']
    dp = weight*(1-e)*h
    return dict(gravity_m_s2=g, gas_density_kg_m3=rf,
        minimum_fluidization_speed_m_s=umf,
        isolated_terminal_speed_m_s=terminal['speed_m_s'],
        terminal_particle_reynolds=terminal['particle_reynolds'],
        required_superficial_speed_m_s=required,
        assumed_retention_ceiling_m_s=ceiling,
        feasible_flow_interval=(g > 0 and required <= ceiling),
        bed_pressure_drop_Pa=dp,
        flow_pressure_power_W_m2=required*(dp+other_dp),
        interpretation='Pressure-work diagnostic only, before blower efficiency. No stable conventional gravity bed at zero gravity; magnetic/rotating beds are separate competitors. At failed retention, the power is not a viable operating saving.')


def cost_screen(settings, gravity_ratio):
    rows = []
    for share in settings['mixing_electricity_cost_shares']:
        for floor in settings['other_speed_floors']:
            lunar_power = suspension_power(gravity_ratio, floor)['turbulent_power_ratio']
            free_power = suspension_power(0, floor)['turbulent_power_ratio']
            shares = {'mixing_electricity':share, 'other':1-share}
            earth = component_cost(100., shares, {}, settings['earth_accepted_yield'], 1.)
            free = component_cost(100., shares, {'mixing_electricity':free_power}, settings['earth_accepted_yield'], 1.)
            for lunar_yield in settings['moon_accepted_yield_sensitivities']:
                moon = component_cost(100., shares, {'mixing_electricity':lunar_power}, lunar_yield, 1.)
                rotating = component_cost(100., shares, {'mixing_electricity':lunar_power}, lunar_yield, 1.)
                e, m, f, a = [x['factory_cost_per_qualified_unit_U'] for x in (earth, moon, free, rotating)]
                rows.append(dict(mixing_electricity_cost_share=share, other_speed_floor=floor,
                    moon_yield=lunar_yield, rotating_array_yield=lunar_yield,
                    earth_and_freefall_yield=settings['earth_accepted_yield'],
                    earth_cost_U_good=e, moon_cost_U_good=m,
                    array_freefall_before_extra_service_U_good=f,
                    array_rotating_before_extra_service_U_good=a,
                    allowable_added_moon_cost_vs_earth_U_good=e-m,
                    minimum_array_freefall_extra_cost_for_moon_parity_U_good=max(0.,m-f),
                    minimum_array_rotating_extra_cost_for_moon_parity_U_good=max(0.,m-a),
                    lunar_cost_saving_fraction_vs_earth=(e-m)/e,
                    reference_electricity_saving_fraction_of_total=share*(1-lunar_power)))
    return rows



def cycle_cost(mixing_share, time_occupancy_share, power_ratio, throughput_ratio):
    """Equal-output cost with longer vessel occupancy and a fixed recipe.

    Earth costs 100 units per attempt. Time occupancy includes its own capital
    service and fixed overhead. Restored throughput requires proportional vessel
    capacity; no scale economies or extra harvest hardware is assumed. Other
    energy, inputs and yield are held fixed. The ratio is not a growth forecast.
    """
    mix = number('mixing cost share', mixing_share)
    time = number('time occupancy share', time_occupancy_share)
    power = number('power ratio', power_ratio)
    q = number('throughput ratio', throughput_ratio, True)
    if mix+time > 1:
        raise ValueError('Disjoint cost shares must sum to at most one')
    electricity_saving = 100*mix*(1-power)
    occupancy_cost = 100*time*(1/q-1)
    total = 100-electricity_saving+occupancy_cost
    return dict(reference_attempt_cost_U=100., moon_attempt_cost_U=total,
        mixing_electricity_saving_U=electricity_saving,
        added_vessel_occupancy_cost_U=occupancy_cost,
        net_cost_saving_fraction=(100-total)/100,
        incremental_cost_allowance_to_restore_earth_cycle_U=electricity_saving,
        interpretation='Illustrative equal-output factory account. Positive savings are before freight, extra process control and quality changes. Forced recovery can replace the occupancy penalty, but its incremental cost must be counted.')


def results():
    s, products = load_inputs()
    g = MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    gravities = {'earth':STANDARD_GRAVITY, 'moon':MOON_SURFACE_GRAVITY, 'array_freefall':0.}
    columns = [{**c, 'worlds':{w:column_screen(c, a, s['temperature_K'], s['maximum_bottom_top_concentration_ratio']) for w,a in gravities.items()}} for c in s['particles']]
    mixing = [suspension_power(g, f, a) for a in s['mixing']['gravity_exponents'] for f in s['mixing']['other_speed_floors']]
    beds = [{ 'world':w, 'process_flow_floor_m_s':f, **fluidized_bed(s['fluid_bed'], a, 1.2 if w=='moon' else 1., f)} for w,a in gravities.items() for f in s['fluid_bed']['process_flow_floors_m_s']]
    cycles = [{'id':c['id'], 'worlds':{w:batch_cycle(c['process_s'], c['other_s'], c['earth_settle_s'], a/STANDARD_GRAVITY) for w,a in gravities.items()}} for c in s['cycles']]
    coupled = [dict(cycle=c['id'], mixing_cost_share=share, **cycle_cost(share, s['costs']['time_occupancy_cost_share'], suspension_power(g,.25)['turbulent_power_ratio'], c['worlds']['moon']['throughput_relative_to_earth'])) for c in cycles for share in s['costs']['mixing_electricity_cost_shares']]
    files = ('slow_settling.py','slow_settling_inputs.json','slow_settling_sources.json','partial_gravity.py','mature_economics.py')
    return dict(schema='terluna.research.slow-settling/1',
        producer=dict(files={str((HERE/f).relative_to(ROOT)):digest(HERE/f) for f in files}, constants=constants_used([__file__]),input_products=s['products']),
        evidence=s['evidence'], gravity_ratio=g,
        stokes_time_multiplier=1/g, same_stokes_speed_diameter_multiplier=g**-.5,
        same_Brownian_gravity_scale_height_diameter_multiplier=g**(-1/3),
        matched_stokes_viscosity_earth_over_moon=1/g,
        casting_reference=products['research/studies/world_economy/results/partial_gravity.json']['casting'],
        particle_columns=columns,mixing_power=mixing,fluidized_beds=beds,batch_cycles=cycles,
        normalized_costs=cost_screen(s['costs'],g), cycle_adjusted_costs=coupled,
        mixing_cost_share_needed_for_25pct_factory_saving=.25/(1-g**1.35),
        assumptions=s['assumptions'],boundaries=s['boundaries'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'research/runs/world_economy/slow_settling.json')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
