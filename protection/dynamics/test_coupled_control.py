"""Regression checks for topology, coupled responses and trust rejection."""
import numpy as np
from scipy.spatial.transform import Rotation
from .coupled_control import contact_values, contact_jacobian, contact_set, secant_update, trust_decision


def test_contact_derivative_includes_both_centres_and_rotating_square_frame():
    states=np.zeros((2,3,6));states[:,1,:3]=[9000.,5000.,450.];states[:,2,:3]=[18000.,200.,1500.]
    frames=Rotation.from_rotvec([[0.,0.,0.],[0.,.03,0.]]).as_matrix()
    cuts=np.array([[0,0,1,2,-1],[1,1,2,0,-1]],float)
    rng=np.random.default_rng(21);response=rng.normal(size=(2,3,6,2))*10
    eps=1e-5;frame_response=[]
    for axis in ['x','z']:
        changed=frames@Rotation.from_euler(axis,eps).as_matrix()
        frame_response.append((changed-frames)/eps)
    fj=np.stack(frame_response,axis=-1)
    actual=contact_jacobian(states,frames,response,fj,cuts)
    base=contact_values(states,frames,cuts)
    for j,axis in enumerate(['x','z']):
        finite=(contact_values(states+eps*response[:,:,:,j],
            frames@Rotation.from_euler(axis,eps).as_matrix(),cuts)-base)/eps
        np.testing.assert_allclose(actual[:,j],finite,atol=3e-4)


def test_normal_side_is_retained_across_contact_not_service_attitude():
    states=np.zeros((3,2,6));states[:,1,0]=5000
    states[:,1,2]=[800.,50.,-200.]
    frames=np.broadcast_to(np.eye(3),(3,3,3))
    cuts=contact_set(states,frames)
    gaps=contact_values(states,frames,cuts)
    np.testing.assert_allclose(gaps,[800,50,-200])


def test_coupled_secant_corrects_cross_member_prediction_and_trust_rejects_false_gain():
    # A control changes a second member through a shared force. The previous
    # diagonal response misses that component; the measured secant must retain it.
    j=np.array([[1.,0.],[0.,1.]])
    step=np.array([.5,0.]);observed=np.array([.6,.2])
    updated=secant_update(j,step,observed)
    np.testing.assert_allclose(updated@step,observed)
    np.testing.assert_allclose(updated[:,1],j[:,1])
    bad=trust_decision(10.,2.,12.,1.,20.)
    assert not bad['accepted'] and bad['next_radius']==.5
    inaccurate=trust_decision(10.,2.,3.,1.,151.)
    assert not inaccurate['accepted']
    good=trust_decision(10.,2.,2.5,1.,5.)
    assert good['accepted'] and good['next_radius']>1
