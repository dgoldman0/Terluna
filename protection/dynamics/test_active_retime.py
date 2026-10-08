import numpy as np
import pytest

from shared import constants as K
from .active_retime import command,finite_tile_acceleration,rotate_state,windows
from .test_active_formation import sample


class StaticEnvironment:
    sigma=.05
    central_fraction=.138
    def at(self,t):return sample()
    def gravity(self,t,q):return -K.MOON_GM*q/np.linalg.norm(q,axis=-1)[...,None]**3


def test_batch_finite_extent_forces_match_independent_tiles():
    states=np.array([[15e6,0,0,0,570.,0],[14e6,2e6,1e6,-80.,560.,2.]])
    env=StaticEnvironment()
    np.testing.assert_allclose(finite_tile_acceleration(env,0,states),
        np.array([finite_tile_acceleration(env,0,s) for s in states]),atol=1e-16,rtol=0)
    point=env.gravity(0,states[0,:3])
    _,_,finite_error=finite_tile_acceleration(env,0,states[0],diagnostics=True)
    assert 0<finite_error<1e-8


def test_smooth_arc_impulse_and_batched_commands():
    state=np.array([15e6,0,0,0,570.,0]); params=np.arange(1.,13.)/10
    times=np.linspace(0,72000,2401)
    norms=[np.linalg.norm(command(t,state,params,72000,3600)) for t in times]
    assert np.trapezoid(norms,times)==pytest.approx(np.linalg.norm(params.reshape(4,3),axis=1).sum(),rel=1e-12)
    for t in [0,1500,10000,24000,71000,72000]:
        batch=command(t,np.array([state,state]),np.array([params,params*2]),72000,3600)
        np.testing.assert_allclose(batch[0],command(t,state,params,72000,3600))
        np.testing.assert_allclose(batch[1],batch[0]*2)
    with pytest.raises(ValueError):windows(1000,600,4)


def test_phase_target_preserves_radius_speed_energy_and_angular_momentum():
    state=np.array([15e6,0,0,0,570.,0]);target=rotate_state(state,.1)
    assert np.linalg.norm(target[:3])==pytest.approx(np.linalg.norm(state[:3]))
    assert np.linalg.norm(target[3:])==pytest.approx(np.linalg.norm(state[3:]))
    np.testing.assert_allclose(np.cross(target[:3],target[3:]),np.cross(state[:3],state[3:]),rtol=1e-14)
