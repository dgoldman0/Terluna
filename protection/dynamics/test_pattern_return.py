from types import SimpleNamespace
import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.spatial.transform import Rotation, RotationSpline

from .pattern_return import (smooth_arc, shooting_maps, minimum_component_impulse,
    proposal_state, encounter_scan, crossing_witnesses)


def test_smooth_arc_normalization_and_peak():
    windows = np.array([[100., 1100.], [2100., 4100.]])
    for j, (a, b) in enumerate(windows):
        assert quad(lambda t: smooth_arc(t, windows)[j], a, b)[0] == pytest.approx(1.)
        assert smooth_arc((a+b)/2, windows)[j] == pytest.approx(2/(b-a))
    np.testing.assert_array_equal(smooth_arc(1500., windows), [0., 0.])


def test_variational_shooting_matches_double_integrator_minimum_impulse():
    env = SimpleNamespace(gravity=lambda t, q: np.zeros_like(q))
    reference = SimpleNamespace(sol=lambda t: np.zeros((6,)+np.shape(t)))
    maps, windows = shooting_maps(env, reference, 0., 30000., burn=2000.)
    response = maps.y[36:, -1].reshape(6, 9)
    midpoint = windows.mean(axis=1)
    expected = np.r_[np.kron((30000.-midpoint)[None, :], np.eye(3)), np.tile(np.eye(3), (1, 3))]
    np.testing.assert_allclose(response, expected, atol=2e-5)
    # HiGHS owns a process-global scheduler. Earlier capacity tests initialize
    # it with their default thread count; the return runner starts a fresh
    # process with one solver thread. Exercise that actual launch contract.
    code = '''
import json
import numpy as np
from types import SimpleNamespace
from protection.dynamics.pattern_return import shooting_maps, minimum_component_impulse, proposal_state
env=SimpleNamespace(gravity=lambda t,q: np.zeros_like(q))
ref=SimpleNamespace(sol=lambda t: np.zeros((6,)+np.shape(t)))
maps,windows=shooting_maps(env,ref,0.,30000.,burn=2000.)
initial=np.zeros((2,6)); target=initial.copy(); target[:,0]=[1000.,-1000.]
control,fit=minimum_component_impulse(maps,ref,initial,target,windows)
assert fit['success'],fit
final=proposal_state(ref,maps,initial,control,[30000.])[0]
print(json.dumps(dict(fit=fit,maximum_error=float(np.max(abs(final-target))))))
'''
    result = json.loads(subprocess.check_output([sys.executable, '-c', code], text=True))
    assert result['fit']['mean_delta_v_m_s'] == pytest.approx(2000./(midpoint[-1]-midpoint[0]), rel=1e-7)
    assert result['maximum_error'] < 1e-6


def test_between_sample_coplanarity_detects_a_missed_finite_square_impact():
    times = np.array([0., 10.])
    command = RotationSpline(times, Rotation.from_matrix(np.tile(np.eye(3), (2, 1, 1))))
    def at(t, ids):
        state = np.zeros((2, 6)); state[0, 2] = -1000.+200*t; state[0, 5] = 200.
        return state[ids]
    states = np.array([at(t, np.arange(2)) for t in times])
    assert encounter_scan(times, states, command)['sampled_pair_violations'] == 0
    result = crossing_witnesses(times, states, command, at)
    assert result['physical_intersections'] == 1
    assert result['witnesses'][0]['time_s'] == pytest.approx(5.)
