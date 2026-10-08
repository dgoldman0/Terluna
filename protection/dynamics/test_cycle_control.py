"""Analytic control bound and smooth arrival boundary checks."""
import json
import os
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
from scipy.spatial.transform import Rotation

from .cycle_control import return_command


def test_arrival_matches_feathered_and_service_frames_with_bounded_rate(monkeypatch):
    monkeypatch.setitem(return_command.__globals__,'sun_frame',lambda env,t:np.eye(3))
    c=return_command(SimpleNamespace(),0.,86400.,28800.)
    np.testing.assert_allclose(c(0.).as_matrix(),Rotation.from_rotvec([-np.pi/2,0,0]).as_matrix(),atol=1e-12)
    np.testing.assert_allclose(c(86400.).as_matrix(),np.eye(3),atol=1e-12)
    t=np.linspace(0.,86400.,2001)
    assert np.linalg.norm(c(t,1),axis=1).max()<np.deg2rad(.006)
    assert np.linalg.norm(c(86400.,1))<1e-8


def test_four_arc_endpoint_lp_matches_analytic_free_particle_minimum():
    script='''
import json
import numpy as np
from types import SimpleNamespace
from protection.dynamics.packing_control import response_maps
from protection.dynamics.cycle_control import cycle_windows,endpoint_controls
env=SimpleNamespace(gravity=lambda t,q:np.zeros_like(q))
end=115200.;windows=cycle_windows(0.,end)
maps=response_maps(env,lambda t:np.zeros(6),0.,end,windows,np.eye(3))
base=lambda t:np.zeros((len(np.atleast_1d(t)),1,6))
target=np.array([[(windows[-1].mean()-windows[0].mean())*.01, (windows[-1].mean()-windows[0].mean())*.02,0,0,0,0]])
u,r=endpoint_controls(base,maps,target,windows)
print(json.dumps(r))
'''
    r=subprocess.run([sys.executable,'-c',script],check=True,capture_output=True,text=True,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    fit=json.loads(r.stdout)
    assert fit['success'],fit
    assert abs(fit['component_L1_mean_m_s']-.06)<1e-8
    assert abs(fit['mean_delta_v_m_s']-.02*np.sqrt(5))<1e-8
    assert fit['restricted_linear_euclidean_lower_bound_m_s']<=fit['mean_delta_v_m_s']
    assert fit['endpoint_position_residual_m']<1e-6
    assert fit['endpoint_velocity_residual_m_s']<1e-9
