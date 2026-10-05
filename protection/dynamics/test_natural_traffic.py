import numpy as np
import pytest
from scipy.integrate import solve_ivp

from .active_global import Traffic,Layout,RayIndex,sun_source
from .natural_traffic import NaturalTraffic,NaturalRayIndex
from .optical import length


def test_released_paths_follow_gravity_and_resolve_integer_members():
    from .test_active_global import make_environment
    env=make_environment()
    # Exercise both radius-order and orbital-sense joins without materializing
    # the production-size optical index inside the repository test suite.
    f=NaturalTraffic(env,1.,Layout(planes=64))
    with pytest.raises(ValueError,match='Prescribed-shell bounds'):
        f.separation_certificate()
    f.propagate(days=1/24,plane_stride=8,phase_nodes=128,max_step=300.)
    j=np.array([0,1,13,20,30,31,32,33,61,63])
    k=np.array([0,100,203,17,31,201,399,3,711,1234])
    initial=f.initial_states(j,2*np.pi*k/f.count[j]+f.phase[j],k%2)
    def rhs(t,y):
        y=y.reshape(-1,6)
        return np.column_stack([y[:,3:],env.gravity(t,y[:,:3])]).ravel()
    sol=solve_ivp(rhs,[0,3600],initial.ravel(),method='DOP853',rtol=2e-12,
        atol=np.tile([1e-5]*3+[1e-9]*3,len(j)),t_eval=[1379,3600],max_step=60.)
    assert sol.success
    for t,truth in zip(sol.t,sol.y.T.reshape(-1,len(j),6)):
        state,acc,_=f.path(t,j,k)
        assert length(state[:,:3]-truth[:,:3]).max()<10.
        np.testing.assert_array_equal(acc,env.gravity(t,state[:,:3]))
    # Both spatial indices find the initial finite-member shadows. The natural
    # index follows propagated curve segments, not the circular phase guess.
    source=sun_source(env.at(0),f.frames(0),[.7,-.3])
    xy=np.array([[0.,0.],[5e6,0.],[-2e6,3e6],[1e6,-6e6]])
    actual=NaturalRayIndex(f,0,source).intercept(xy)
    prescribed=RayIndex(Traffic(env,1.,f.layout),0,source).intercept(xy)
    np.testing.assert_array_equal(actual,prescribed)
