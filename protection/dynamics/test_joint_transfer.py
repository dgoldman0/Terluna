"""Analytic variable-time endpoints and conservative member ratings."""
import os,subprocess,sys,json
import numpy as np


def test_distinct_arrival_times_and_rating_rejection():
    script='''
import json
import numpy as np
from types import SimpleNamespace
from protection.dynamics.packing_control import response_maps,matrices
from protection.dynamics.joint_transfer import fit_controls
windows=np.array([[0.,100.],[300.,400.]])
env=SimpleNamespace(gravity=lambda t,q:np.zeros_like(q))
maps=response_maps(env,lambda t:np.zeros(6),0.,800.,windows,np.eye(3))
base=lambda t:np.zeros((len(np.atleast_1d(t)),2,6))
target=np.zeros((2,6));target[:,0]=[500,300];target[:,3]=[1,.5]
times=np.array([600.,800.])
u,r=fit_controls(base,maps,target,windows,np.ones((2,2))*.03,times)
b=matrices(maps,times,2)[:,:,6:]
actual=np.einsum('nij,nj->ni',b,u.reshape(2,6))
_,failure=fit_controls(base,maps,target,windows,np.ones((2,2))*.001,times)
print(json.dumps([r,actual.tolist(),failure]))
'''
    p=subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    r,actual,failure=json.loads(p.stdout)
    assert r['success'] and not failure['success']
    expected=np.zeros((2,6));expected[:,0]=[500,300];expected[:,3]=[1,.5]
    np.testing.assert_allclose(actual,expected,atol=1e-7)
