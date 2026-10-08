from types import SimpleNamespace
import json
import os
import subprocess
import sys

import numpy as np

from .packing_control import packing_offsets, response_maps, matrices, measured_base
from .active_formation import lattice


def test_reordered_packing_preserves_projected_positions_and_depth_population():
    old=lattice(19,8500.,4000.)[0]
    np.testing.assert_array_equal(packing_offsets(4000.),old)
    new=packing_offsets(8000.,True)
    np.testing.assert_array_equal(new[:,1:],old[:,1:])
    np.testing.assert_array_equal(np.sort(new[:,0]),np.sort(2*old[:,0]))
    assert not np.array_equal(new[:,0],2*old[:,0])


def test_two_arc_maps_recover_free_particle_impulse_response_and_defect_transport():
    env=SimpleNamespace(gravity=lambda t,q:np.zeros_like(q))
    windows=np.array([[0.,200.],[400.,600.]])
    maps=response_maps(env,lambda t:np.zeros(6),0.,800.,windows,np.eye(3))
    end=matrices(maps,[800.])[0]
    np.testing.assert_allclose(end[:3,6:],np.c_[np.eye(3)*700,np.eye(3)*300],atol=1e-6)
    np.testing.assert_allclose(end[3:,6:],np.tile(np.eye(3),(1,2)),atol=1e-8)
    base=lambda t:np.zeros((len(np.atleast_1d(t)),1,6))
    t=np.array([0.,200.,400.]);state=base(t)
    state[:,0,0]=100+.2*t;state[:,0,3]=.2
    corrected,_=measured_base(base,maps,t,state)
    np.testing.assert_allclose(corrected([100.,800.])[:,0,0],[120.,260.],atol=1e-6)


def test_independent_members_can_open_a_gap_with_minimum_l1_impulse():
    script='''
import json
import numpy as np
from protection.dynamics.packing_control import independent_impulses
cut=dict(i=0,j=1,gradient=np.array([10000.,0,0,0,0,0]),rhs=4000.)
u,r=independent_impulses(np.zeros((2,2,3)),[cut],np.array([[0.,7200.],[10800.,18000.]]))
print(json.dumps([u.tolist(),r]))
'''
    run=subprocess.run([sys.executable,'-c',script],check=True,capture_output=True,
        text=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    u,result=json.loads(run.stdout);u=np.array(u)
    assert result['success'],result
    assert abs(u[0,0,0]-u[1,0,0]-.4)<1e-8
    assert abs(result['mean_delta_v_m_s']-.2)<1e-8
    assert result['maximum_thrust_acceleration_m_s2']<=.001
