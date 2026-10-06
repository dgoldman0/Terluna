"""Checks on control accounting, temporal independence and evidence gates."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.coupled_control import contact_values
from .expanded_model import ExpandedExperiment, DIM, PACKAGE
from .fleet_run import digest


def test_partial_control_energy_and_peak_match_actual_sine_squared_pulse():
    window = np.array([[10., 110.]])
    for f in [.1, .35, .75, 1.]:
        numerical = quad(lambda t: smooth_arc(t, window)[0], 10., 10.+100*f)[0]
        exact = f-np.sin(2*np.pi*f)/(2*np.pi)
        np.testing.assert_allclose(numerical, exact, atol=1e-12)
    assert smooth_arc(60., window)[0] == 2/100.


def test_expanded_controls_include_in_plane_and_distinct_acquisition_times():
    ex = ExpandedExperiment(); x = np.zeros(DIM)
    f = ex.command(x)(13.8*3600).as_matrix()
    assert ex.early.shape == (361, 3, 10)
    assert np.max(abs(ex.early[:, :, 6]@f[:, 2])) < 1e-12
    assert np.max(abs(ex.early[:, :, 7]@f[:, 2])) < 1e-12
    assert np.linalg.norm(ex.early[:, :, 6]) > .1
    x[10] = 1.; x[13] = 1.
    u, windows = ex.controls(x)
    assert windows[1, 1] == windows[2, 0]
    assert np.linalg.norm(u[:, 1]) > 1. and np.linalg.norm(u[:, 2]) > 1.
    assert ex.onset(10) < ex.onset(13) < ex.dates(x)[0]
    np.testing.assert_allclose(np.sum(ex.early*ex.mass[:, None, None], axis=0), 0., atol=1e-7)


def test_compact_ephemeris_and_basis_are_checked_inputs():
    d = json.loads((PACKAGE/'inputs.json').read_text())
    assert digest(PACKAGE/'ephemeris.npz') == d['compact_ephemeris_sha256']
    assert digest(PACKAGE/'basis.npz') == d['basis_sha256']
    ex = ExpandedExperiment()
    assert ex.env.limits[1] > 2*ex.arrival0+46800
    assert ex.initial.shape == (361, 6)
    assert d['supply'] == dict(generation=None, storage=None, unmet_demand=None)
    assert not d['hardware_resized']


def test_new_in_plane_control_can_change_edge_clearance_without_normal_motion():
    states = np.zeros((1, 2, 6)); states[0, 1, 0] = 10050.
    cuts = np.array([[0, 0, 1, 0, -1.]])
    frames = np.eye(3)[None]
    a = contact_values(states, frames, cuts)[0]
    states[0, 1, 0] += 200.
    b = contact_values(states, frames, cuts)[0]
    assert a == 50. and b == 250.


def test_explicit_budget_guard_stops_before_starting_a_coupled_evaluation(monkeypatch):
    import pytest
    from . import expanded_model
    ex = ExpandedExperiment()
    def expired(): raise TimeoutError('test deadline')
    monkeypatch.setattr(expanded_model, 'check_budget', expired)
    with pytest.raises(TimeoutError, match='test deadline'):
        ex.evaluate(np.zeros(DIM), end=21601., account=False)


def test_saved_candidate_replay_starts_without_any_ignored_work_inputs(tmp_path):
    """A 360-second representative replay is a repository verification gate."""
    import os, subprocess, sys
    if not (PACKAGE/'last_trial.json').exists():
        raise AssertionError('Portable candidate has not been published')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    result = subprocess.run([sys.executable, '-m', 'research.studies.solar_shield_array.expanded_run',
        'replay', '--run-dir', str(tmp_path), '--candidate', str(PACKAGE/'last_trial.json'),
        '--case', 'smoke', '--end-hour', '6.1', '--cpu-budget', '60'],
        cwd=PACKAGE.parents[3], env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    row = json.loads((tmp_path/'smoke_replay.json').read_text())
    assert row['completed_horizon'] and row['end_s'] == 6.1*3600.
    assert not row['clearance_events_s'] and not row['accepted_cycle']
    ledger = json.loads((tmp_path/'ledger.json').read_text())
    assert ledger['charged_cpu_s'] < 60.
