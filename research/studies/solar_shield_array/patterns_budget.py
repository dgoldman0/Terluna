"""Separate service expenditure from unsolved recurring local-pattern cost."""
import json

from shared.provenance import constants_used, constants_changed
from protection.dynamics.coverage_pilot import conditional_cost
from .fleet_run import digest
from .cycling_search import HERE, ROOT, SCENARIO


def main():
    input_path = HERE/'results/patterns.json'
    p = json.loads(input_path.read_text())
    for path, h in p['producer']['source_hashes'].items():
        if digest(ROOT/path) != h:
            raise ValueError('Pattern producer changed: '+path)
    if constants_changed(p['producer']['constants']):
        raise ValueError('Pattern constants changed')
    r = p.get('replay', p['coupled'][0])
    prop = SCENARIO['propulsion']
    import numpy as np
    cant = np.cos(np.deg2rad(prop['cant_deg']))
    mass = r['optical_mass_kg']
    preparation = r['preparation_delta_v_reference_mean_m_s']
    terminal = r['endpoint_reference_velocity_defect_mean_m_s']
    capacity = json.loads((HERE/'results/pilot_validation.json').read_text())
    inventory = capacity['constraint_generation']['refined_solution']['inventory_tiles']
    sources = {'research/studies/solar_shield_array/patterns_budget.py': digest(__file__),
        'protection/dynamics/coverage_pilot.py': digest(ROOT/'protection/dynamics/coverage_pilot.py')}
    out = dict(schema='terluna.research.natural-pattern-budget/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'patterns.json': digest(input_path),
                'pilot_validation.json': digest(HERE/'results/pilot_validation.json')}),
        local_count=r['count'], local_optical_mass_kg=mass,
        service=dict(duration_s=r['duration_s'], electric_delta_v_m_s=0., propulsion_energy_J=0.),
        preparation=dict(reference='Instantaneous common Sun-frame velocity at the exact initial positions; no acquisition propagation.',
            mean_delta_v_m_s=preparation, optical_mass_impulse_N_s=mass*preparation,
            electric_energy_J=mass*preparation*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant),
            propellant_kg=mass*preparation/(prop['exhaust_velocity_m_s']*cant),
            hardware_mass_included=False, maneuver_executed=False),
        endpoint=dict(reference='Velocity mismatch to the final Sun-frame lattice; spatial mismatch remains.',
            mean_velocity_mismatch_m_s=terminal,
            maximum_position_mismatch_m=r['endpoint_reference_position_defect_max_m'],
            closure_maneuver_executed=False),
        recurrence_diagnostic=dict(velocity_only_two_impulse_proxy_m_s=preparation+terminal,
            interpretation='A kinematic comparison, neither a solved return nor a bound on recurring impulse; it leaves kilometre-scale position closure, handovers and recovery unpaid.'),
        inventory_impulse_trade=[dict(inventory_tiles=n, cases=[
            conditional_cost(n, dv, prop, hardware_multiplier=1.2, duty=.1)
            for dv in [1., 3., 5.]]) for n in [inventory, 25e6, 35e6]],
        cost_scope='Conditional recurring impulse inputs, with 20% extra fixed mass, 10% thrust duty, power-hardware/7-day-propellant feedback. These masses require new sail-dynamics propagation.',
        held_screen_benchmark_TW=238.35, measured_recurring_delta_v_m_s_per_day=None,
        acquisition_events_executed=0, service_passages_independently_replayed=int('replay' in p),
        handovers_executed=0, return_legs_executed=0, recovery_events_executed=0,
        recurring_maneuver_frequency_per_day=None, collection_capacity_W=None,
        delivered_power_W=None, accepted_fleet=False)
    (HERE/'results/patterns_budget.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
