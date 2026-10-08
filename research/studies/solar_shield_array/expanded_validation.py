"""Exact trial diagnostics and independently tightened, event-stopped audit."""
import json
import time
import numpy as np
from shared import constants as K
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.pattern_return import beam_and_rate_screen
from protection.dynamics.coupled_control import contact_values, contact_set
from protection.dynamics.optical import length
from .coupled_validate import audit_clearance
from .joint_common import spline, state_at
from .coupled_common import SPECIFIC
from .expanded_model import DIM, PACKAGE


def assess_trial(ex, folder, case):
    from .expanded_run import raw_read, feature, write, energy_parts
    base = raw_read(folder/'base.npz'); model = raw_read(folder/'model.npz')
    raw = raw_read(folder/(case+'_replay.npz'))
    row = json.loads((folder/(case+'_replay.json')).read_text())
    proposal = json.loads((folder/(case+'_proposal.json')).read_text())
    x = np.array(proposal['x']); d = np.array(proposal['step']); cols = model['columns']
    f = feature(ex, raw, x, base, model['cuts'])
    fresh = contact_set(f['state'], f['frames'], reserve=1000.)
    old_keys = set(map(tuple, model['cuts']))
    new_count = sum(tuple(c) not in old_keys for c in fresh)
    predicted = model['terminal']+model['terminal_jac']@d
    defect = (f['terminal']-predicted).reshape(3, 361, 6)
    # All terms use the SAME retained contact set for predicted/actual merit.
    def merit(energy, terminal, gaps, service):
        return energy[-1]/1e14+10*np.mean(terminal**2)+2*np.mean(service**2)+20*max(0., (150.-gaps.min())/1000.)+20*max(abs(terminal))
    before = merit(model['energy'], model['terminal'], model['gaps'], model['service'])
    pred = merit(np.array(proposal['predicted_energy_J']), predicted,
                 model['gaps']+model['gaps_jac']@d, model['service']+model['service_jac']@d)
    actual = merit(f['energy'], f['terminal'], f['gaps'], f['service'])
    ratio = (before-actual)/(before-pred) if before > pred else -1.
    baseline_event = json.loads((folder/'base.json').read_text())['clearance_events_s'][0]
    actual_event = row['clearance_events_s'][0] if row['clearance_events_s'] else row['end_s']
    error = float(length(defect[:, :, :3]*1e5).max())
    physical_caps = bool(f['energy'][0] <= 51.859e12 and f['energy'][1] <= 100e12 and
                         row['overloaded_members'] == 0 and row['maximum_translation_m_s2'] <= .001 and
                         row['maximum_attitude_rate_deg_s'] <= .1)
    accepted_step = bool(ratio >= .1 and actual < before and error <= 150 and
                         actual_event >= baseline_event-.1 and physical_caps)
    coverage = []
    for t in np.linspace(ex.dates(x)[0], ex.dates(x)[1], 13):
        state = state_at(spline(raw), t)
        coverage.append(local_coverage(state, ex.env.at(t), state.mean(axis=0),
                                       limb_count=16, interior=8, rotation=.137))
    out = dict(actual_predicted_merit_ratio=float(ratio), before_merit=float(before),
               predicted_merit=float(pred), actual_merit=float(actual),
               maximum_terminal_prediction_error_m=error,
               maximum_terminal_velocity_prediction_error_m_s=float(length(defect[:, :, 3:]*5.).max()),
               new_contact_rows_exposed=new_count,
               maximum_service_coverage_prediction_error=float(np.max(abs(f['service']-model['service']-model['service_jac']@d))),
               prediction_error_by_stage_m=length(defect[:, :, :3]*1e5).max(axis=1),
               actual_energy_J=f['energy'], predicted_energy_J=proposal['predicted_energy_J'],
               actual_terminal_rms=float(np.sqrt(np.mean(f['terminal']**2))),
               actual_clearance_restoration_slack_m=max(0., 150.-float(f['gaps'].min())),
               terminal_slack=max(abs(f['terminal'])), first_event_hour=actual_event/3600.,
               actual_resource_caps_pass=physical_caps, accepted_restoration_step=accepted_step,
               next_radius=proposal['radius'] if accepted_step else proposal['radius']/2.,
               accepted_cycle=False, next_service_credited=False,
               diagnostic_next_service=dict(dates=13, sources_per_date=25,
                  minimum_source_coverage=min(r['minimum_source_coverage'] for r in coverage),
                  receiver_inside_four_radii=all(r['window_centre_distance_m']+10000 <= 4*K.MOON_RADIUS for r in coverage)),
               cycle=ex.summary(raw, x),
               reason='Post-contact motion is diagnostic. A restoration step is never accepted physical service or recurrence.')
    write(folder/(case+'_assessment.json'), out)
    print('trial', case, 'ratio', ratio, 'prediction error m', error,
          'event h', actual_event/3600, 'step accepted', accepted_step, flush=True)
    return out


def validate(ex, folder, case):
    from .expanded_run import raw_read, write, save_case
    from .expanded_model import check_budget
    proposal = json.loads((folder/(case+'_proposal.json')).read_text()); x = np.array(proposal['x'])
    cmd = ex.command(x); coarse = raw_read(folder/(case+'_replay.npz'))
    coarse_row = json.loads((folder/(case+'_replay.json')).read_text())
    # Independent original-epoch replay: does not inherit the coarse hour-6 state.
    row, fine, sol = ex.evaluate(x, ex.arrays['initial_epoch_state'], 0., fine=True,
                                 stop=True, output_step=60.)
    save_case(case+'_validated', row, fine, folder)
    original_dense = sol.sol
    def checked_dense(t):
        check_budget()
        return original_dense(t)
    sol.sol = checked_dense
    finish = row['end_s']; end = min(finish, coarse['t'][-1])
    times = np.r_[np.arange(21600., end, 30.), end]
    delta = sol.sol(times).T.reshape(-1, 361, 6)-state_at(spline(coarse), times)
    error = float(length(delta[:, :, :3]).max())+row['reconstruction_error_m']+coarse_row['reconstruction_error_m']
    service_guard = audit_clearance(sol, cmd, 0., 21600., error)
    ends = [43200., max(21600., finish-120.)]
    prefixes = [audit_clearance(sol, cmd, 21600., e, error) for e in sorted(set(ends)) if e <= finish]
    stop_guard = audit_clearance(sol, cmd, max(21600., finish-120.), finish, error)
    coverage = []
    for t in np.linspace(0., 21600., 49):
        check_budget()
        state = sol.sol(t).reshape(361, 6)
        coverage.append(local_coverage(state, ex.env.at(t), state.mean(axis=0),
                                       limb_count=16, interior=8, rotation=.137))
    t = fine['t']; p = fine['power_W']; light = fine['first_intercept_W']; att = fine['attitude_acceleration']
    energy = float(np.trapezoid(p.sum(axis=1), t)); mask = t <= 43200.
    endpoint = fine['state'][-1]
    out = dict(schema='terluna.expanded-cycle-validation/1', case=case,
        fine_settings=row, maximum_position_disagreement_with_reconstruction_m=error,
        maximum_velocity_disagreement_m_s=float(length(delta[:, :, 3:]).max()),
        event_time_disagreement_s=abs(finish-coarse_row['clearance_events_s'][0]) if coarse_row['clearance_events_s'] else None,
        service_clearance=service_guard, prefix_clearance=prefixes, stopping_clearance=stop_guard,
        initial_service=dict(dates=49, sources_per_date=25, clear_side_m=9890.,
            minimum_source_coverage=min(r['minimum_source_coverage'] for r in coverage),
            receiver_inside_four_radii=all(r['window_centre_distance_m']+10000 <= 4*K.MOON_RADIUS for r in coverage),
            circle_sagitta_allowance_m=10000*(1-np.cos(np.pi/512)),
            receiver_reference='Actual fine mean; initial controls identical to preserved service. Coarse/fine centre error is included in the positional allowance.'),
        beam_and_rate=beam_and_rate_screen(ex.env, t, fine['state'], cmd),
        account=dict(executed_end_hour=finish/3600., electrical_energy_J=energy,
            twelve_hour_energy_J=float(np.trapezoid(p[mask].sum(axis=1), t[mask])) if finish>=43200 else None,
            attitude_energy_J=float(np.trapezoid((att*ex.mass*SPECIFIC).sum(axis=1), t)),
            translation_energy_J=float(np.trapezoid((length(fine['translation'])*ex.mass*SPECIFIC).sum(axis=1), t)),
            intercepted_light_J=float(np.trapezoid(light.sum(axis=1), t)),
            redirected_light_J=float(np.trapezoid(light.sum(axis=1), t))*ex.env.central_fraction,
            maximum_individual_power_W=float(p.max()), installed_to_peak_minimum=float(np.min(ex.capacity/p.max(axis=0))),
            overloaded_members=int(np.count_nonzero(p.max(axis=0)>ex.capacity)),
            per_member_peak_W=p.max(axis=0), per_member_installed_W=ex.capacity,
            maximum_torque_N_m=float(length(fine['torque_N_m']).max()),
            ideal_exhaust_kg=energy*1.4/30000**2, propulsion_waste_heat_J=.3*energy,
            loaded_mass_kg=float(ex.mass.sum()), hardware_resized=False,
            generation=None, storage=None, unmet_demand=None,
            torque_limit='Shared installed electric-couple capacity, including translation demand; independent axis ratings and placement are unspecified.'),
        terminal=dict(time_s=finish, state_m_m_s=endpoint, attitude_xyzw=cmd(finish).as_quat(),
                      angular_velocity_rad_s=cmd(finish, 1)),
        complete_finite_sequence=False, repeatability_passed=False, subsequent_cycles_executed=0,
        next_service_executed=False, subsequent_departure_executed=False,
        accepted_cycle=False, supply_validated=False, handover_validated=False,
        reading_rule='Measured executed quantities stop at first 100 m event. Post-contact optimization tails earn no service or recurrence credit.')
    write(folder/'validation.json', out)
    print('validated', case, finish/3600, 'error', error, 'TJ', energy/1e12, flush=True)
    return out
