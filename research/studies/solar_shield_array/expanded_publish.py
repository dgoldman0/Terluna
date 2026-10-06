"""Preserve compact inputs, local optimizer state and executed chronologies."""
import importlib.metadata
import hashlib
import io
import json
import shutil
import sys
import time
from pathlib import Path
import numpy as np
from .expanded_model import PACKAGE, HERE, ROOT, DIM, ExpandedExperiment
from .expanded_run import RUN, write, raw_read
from .fleet_run import digest


def main():
    before = time.process_time()
    old_result = HERE/'results/expanded_cycle.json'
    old = json.loads(old_result.read_text()) if old_result.exists() else {}
    exports = old.get('artifact_export_cpu_s', [old['final_export_cpu_s']] if 'final_export_cpu_s' in old else [])
    inputs = json.loads((PACKAGE/'inputs.json').read_text())
    names = ['numpy', 'scipy', 'shapely', 'jplephem', 'pytest']
    versions = {name: importlib.metadata.version(name) for name in names}
    (PACKAGE/'requirements.txt').write_text(''.join(f'{k}=={v}\n' for k, v in versions.items()))
    trial = json.loads((RUN/'actual_proposal.json').read_text())
    shutil.copyfile(RUN/'actual_proposal.json', PACKAGE/'last_trial.json')
    model = raw_read(RUN/'model.npz')
    model['steps'] = np.array([json.loads((RUN/f'direction_{j:02d}.json').read_text())['h'] for j in model['columns']])
    np.savez_compressed(PACKAGE/'optimizer_model.npz', **model)
    # Verification of saved evaluations only; no new trajectory is propagated.
    ex = ExpandedExperiment(); small = raw_read(RUN/'small_step_replay.npz')
    small_proposal = json.loads((RUN/'small_step_proposal.json').read_text())
    scale = np.array([1e5]*3+[5.]*3)
    actual = (ex.terminal(small, small['x'])/scale).ravel()
    predicted = model['terminal']+model['terminal_jac']@np.array(small_proposal['step'])
    error = (actual-predicted).reshape(3, 361, 6)*scale
    position_error = np.linalg.norm(error[:, :, :3], axis=-1).max(axis=1)
    prior_error = json.loads((RUN/'actual_assessment.json').read_text())['maximum_terminal_prediction_error_m']
    write(RUN/'direction_check.json', dict(step_fraction=.1,
        prediction_error_by_stage_m=position_error,
        maximum_terminal_prediction_error_m=float(position_error.max()),
        maximum_terminal_velocity_prediction_error_m_s=float(np.linalg.norm(error[:, :, 3:], axis=-1).max()),
        small_to_full_prediction_error_ratio=float(position_error.max()/prior_error),
        criterion='A tenfold smaller step distinguishes roughly linear stencil bias from quadratic step curvature only locally; two amplitudes cannot isolate shadow-topology and stencil effects.',
        new_trajectory_propagated=False, accepted_cycle=False))
    validation = json.loads((RUN/'validation.json').read_text())
    case = validation['case']; raw = raw_read(RUN/(case+'_validated.npz'))
    # Complete chronology for every member, with no redundant redirected
    # history: the preserved spectral fraction determines it exactly.
    keep = ['t', 'power_W', 'attitude_acceleration', 'translation', 'first_intercept_W', 'torque_N_m']
    history_buffer = io.BytesIO()
    np.savez_compressed(history_buffer, **{k: raw[k] for k in keep})
    history_bytes = history_buffer.getvalue()
    history_parts = []
    for start in range(0, len(history_bytes), 8*1024*1024):
        name = f'executed_history.npz.part{len(history_parts):02d}'
        (PACKAGE/name).write_bytes(history_bytes[start:start+8*1024*1024])
        history_parts.append(name)
    products = {}
    for path in sorted(RUN.glob('*.json')):
        if path.name != 'ledger.json': products[path.stem] = json.loads(path.read_text())
    ledger = json.loads((RUN/'ledger.json').read_text())
    files = set()
    for name in ['expanded_model.py', 'expanded_run.py', 'expanded_validation.py', 'expanded_publish.py',
                 'coupled_common.py', 'coupled_restart.py', 'coupled_cycle.py', 'coupled_stable.py',
                 'joint_common.py', 'fleet_run.py', 'cycling_search.py', 'coupled_validate.py']:
        files.add(HERE/name)
    files.update((ROOT/'protection/dynamics').glob('*.py'))
    files.update((ROOT/'protection/dynamics').glob('*.cpp'))
    files.update([ROOT/'shared/constants.json', HERE/'scenario.json',
                  ROOT/'protection/spectra/stack_short_wave.json',
                  ROOT/'protection/transmission/results/uv_transmission.json'])
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sorted(files)}
    bundle = dict(schema='terluna.expanded-cycle-restart/1',
        physical_selected_x=[0.]*DIM, last_trial_x=trial['x'],
        optimizer_anchor_x=[0.]*DIM,
        next_trust_radius=.1,
        trust_update='Full trial rejected at radius 1; automatic suggestion 0.5; retained diagnostic scaled to one tenth of its step. Rebuild smaller, signed derivative stencils before another substantive step.',
        measured_columns=products['rank']['measured_columns'],
        unmeasured_columns=products['rank']['unmeasured_columns'],
        finite_difference_step=.2, random_seed=20261006, random_draws_used=0,
        source_hashes=source_hashes,
        input_hashes={str(p.relative_to(ROOT)): digest(p) for p in [
            HERE/'results/joint_seed.json', HERE/'results/coupled_restart.json',
            PACKAGE/'inputs.json', PACKAGE/'ephemeris.npz', PACKAGE/'basis.npz']},
        payload_hashes={name:digest(PACKAGE/name) for name in ['last_trial.json', 'optimizer_model.npz']+history_parts},
        executed_history=dict(parts=history_parts, size_bytes=len(history_bytes),
            sha256=hashlib.sha256(history_bytes).hexdigest(),
            encoding='Concatenate parts in listed order to recover the original NumPy NPZ bytes'),
        dependencies=versions, python=sys.version, compiler='g++ -O3 -std=c++17 -shared -fPIC',
        ephemeris=inputs['ephemeris_manifest'], epoch_tdb=inputs['epoch_tdb'],
        accepted_cycle=False, supply_validated=False, physical_hardware_changed=False,
        remaining_constraints=inputs['limits'],
        criterion=inputs['recurrence'],
        roles='Physical selected state remains original slower turn. Last trial and saved Jacobian/contact model include post-contact diagnostic motion. No local campaign was launched.')
    if Path('/tmp/expanded_pre_guard_hashes.json').exists():
        bundle['execution_source_hashes_before_budget_guard_repair'] = json.loads(Path('/tmp/expanded_pre_guard_hashes.json').read_text())
    bundle['budget_guard_repair'] = 'Final source adds explicit process-clock checks after the timer-only validation stage overshot and repairs custom replay output directories. These changes affect stopping/output behavior, not force/control equations or saved numerical values. The executed pre-repair source is reconstructible with budget_guard.patch reversed.'
    write(PACKAGE/'restart.json', bundle)
    out = dict(schema='terluna.research.expanded-cycle/1', completed=True,
        conclusion='Unresolved search with quantified local control, model and budget limitations',
        accepted_cycle=False, accepted_return=False, accepted_fleet=False,
        held_array_comparison_TW=238.35, operating_target_TW=23.835,
        unvalidated_capacity_estimate_tiles=17.29e6,
        primary_magnetic_integration=False, hardware_mass_changed=False,
        source_hashes=source_hashes, products=products, resources=ledger,
        numerical_budget_s=1800., numerical_budget_overrun_s=max(0., ledger['charged_cpu_s']-1800.),
        further_research_stopped=True,
        final_export_cpu_s=time.process_time()-before,
        reading_rule='Exact complete sequence and all physical resource/service gates are required for acceptance; all search tails after first failed clearance remain diagnostic.')
    out['artifact_export_cpu_s'] = exports+[out['final_export_cpu_s']]
    write(HERE/'results/expanded_cycle.json', out)
    print('Published compact continuation; dependencies', versions, flush=True)


if __name__ == '__main__': main()
