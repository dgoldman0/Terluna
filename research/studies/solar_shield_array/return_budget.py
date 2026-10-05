"""Account for the stopped return prefix without inventing recurring service."""
import json

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation, RotationSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.parallel_shadow import minimum_clearance
from protection.dynamics.pattern_return import crossing_witnesses
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO


def checked(name):
    p = json.loads((HERE/'results'/name).read_text())
    for path, h in p['producer']['source_hashes'].items():
        if digest(ROOT/path) != h:
            raise ValueError('Changed source: '+path)
    if constants_changed(p['producer']['constants']):
        raise ValueError('Changed constants')
    return p


def main():
    screen = checked('return_screen.json'); prefix = checked('return_prefix.json')
    fine = prefix['runs'][-1]
    for p in [screen, fine]:
        if digest(ROOT/p['raw_path']) != p['raw_sha256']:
            raise ValueError('Changed raw trace')
    a = np.load(ROOT/screen['raw_path']); b = np.load(ROOT/fine['raw_path'])
    state = b['state'][-1]; t = float(b['t'][-1])
    u, edge_a, edge_b = solar_frame(environment(3.5).at(t)).T
    distance, pair, _ = minimum_clearance(state, np.column_stack([edge_a, -edge_b, -u]))
    i, j = fine['limiting_pair']
    centre_distance = float(length(state[i, :3]-state[j, :3]))
    predicted = CubicHermiteSpline(a['times'], a['state'][:, :, :3], a['state'][:, :, 3:])(b['t'])
    difference = float(length(predicted-b['state'][:, :, :3]).max())
    interpolation = CubicHermiteSpline(b['t'], b['state'][:, :, :3], b['state'][:, :, 3:])
    command = RotationSpline(a['times'], Rotation.from_matrix(a['frames']))
    def state_at(ti, ids):
        return np.c_[interpolation(ti)[ids], interpolation(ti, 1)[ids]]
    crossings = crossing_witnesses(b['t'], b['state'], command, state_at, retain=8)
    mass = len(state)*SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    prop = SCENARIO['propulsion']; cant = np.cos(np.deg2rad(prop['cant_deg']))
    dv = fine['mean_executed_delta_v_m_s']
    sources = {'research/studies/solar_shield_array/return_budget.py': digest(__file__),
        'research/studies/solar_shield_array/fleet_run.py': digest(HERE/'fleet_run.py'),
        'research/studies/solar_shield_array/scenario.json': digest(HERE/'scenario.json'),
        'protection/dynamics/active_formation.py': digest(ROOT/'protection/dynamics/active_formation.py'),
        'protection/dynamics/parallel_shadow.py': digest(ROOT/'protection/dynamics/parallel_shadow.py')}
    sources['protection/dynamics/pattern_return.py'] = digest(ROOT/'protection/dynamics/pattern_return.py')
    out = dict(schema='terluna.research.stopped-return-budget/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={n: digest(HERE/'results'/n) for n in ['return_screen.json', 'return_prefix.json']}),
        frozen_state_diagnostic=dict(at_s=t, actual_command_clearance_m=fine['minimum_clearance_m'],
            sun_facing_plane_minimum_clearance_m=distance,
            sun_facing_limiting_pair=None if pair is None else pair.tolist(),
            failed_pair_centre_distance_m=centre_distance,
            interpretation='Same saved positions with a changed common plane; no alternative force history was propagated. A coplanar pair closer than one side length cannot be separated by changing square roll alone.'),
        linear_prediction_prefix_discrepancy_m=difference,
        prefix_between_sample_crossing_audit=crossings,
        prefix_crossing_audit_scope='Coplanarity roots on a 30 s Hermite reconstruction of the tighter coupled prefix, with finite-square tests; no continuous certificate.',
        linear_comparison_scope='120 s Hermite reconstruction of the saved variational position/velocity proposal; diagnostic only.',
        executed_prefix=dict(members=len(state), optical_mass_kg=mass,
            duration_s=fine['seconds_after_service'], mean_delta_v_m_s=dv,
            maximum_delta_v_m_s=fine['max_executed_delta_v_m_s'],
            optical_mass_impulse_N_s=mass*dv,
            propulsion_energy_J=mass*dv*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant),
            propellant_kg=mass*dv/(prop['exhaust_velocity_m_s']*cant),
            thrust_arcs_started=1, thrust_arcs_completed=0,
            mass_scope='50 g/m2 optical mass; hardware, attitude actuation and propellant feedback unclosed.'),
        endpoint_only_proposal=dict(mean_delta_v_m_s=screen['proposal_mean_delta_v_m_s'],
            calendar_normalized_m_s_per_day=screen['proposal_mean_delta_v_m_s_per_day'],
            interpretation='Geometrically rejected variational proposal, not an operating-cost measurement or universal lower bound.'),
        measured_recurring_delta_v_m_s_per_day=None, fleet_propulsion_W=None,
        collection_capacity_W=None, delivered_power_W=None,
        full_return_executed=False, handovers_executed=0, recovery_executed=False,
        accepted_fleet=False, held_screen_comparison_TW=238.35)
    (HERE/'results/return_budget.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
