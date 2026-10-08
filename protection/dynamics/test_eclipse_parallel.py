import numpy as np
from .eclipse_parallel import moments
from .fast_parallel import moments as clear_moments


def test_full_body_and_unocculted_limits():
    q=np.array([[0.,0.,0.],[1000,1000,12000.]])
    s=np.array([0.,0.,1.5e11]);f=np.eye(3)
    b=np.tile([[0.,0.,1.5e7],[1e9,0.,1e8]],(2,1,1));r=[1.7e6,6.4e6]
    np.testing.assert_allclose(moments(q,f,s,b,r),0,atol=1e-8)
    b[:,0,0]=1e8
    np.testing.assert_allclose(moments(q,f,s,b,r),clear_moments(q,f,s),atol=1e-6)


def test_joint_partial_union_and_refinement():
    # A body limb and a foreground tile block the same half: multiplying
    # independent half visibilities would erroneously give one quarter.
    q=np.array([[0.,0.,0.],[4945.,0.,10000.]])
    s=np.array([0.,0.,1.5e11]);f=np.eye(3);r=1.7e6
    b=np.tile([[r,0.,1.5e7],[1e9,0.,1e8]],(2,1,1))
    a=moments(q,f,s,b,[r,6.4e6],grid=128)
    c=moments(q,f,s,b,[r,6.4e6],grid=256)
    assert .49<a[0,0]/9890**2<.51
    assert abs(a[0,0]-c[0,0])/9890**2<.005
