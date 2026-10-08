import numpy as np
from scipy.spatial.transform import Rotation
from .oriented_tiles import moments,square_distance
from .fast_parallel import moments as parallel_moments


def test_parallel_limit_with_both_illuminated_faces():
    q=np.array([[0.,0.,0.],[3500,4200,12000],[8000,2000,-1000.]])
    f=Rotation.from_euler('xyz',[21,35,8],degrees=True).as_matrix()
    for s in [np.array([0.,0.,1.5e11]),np.array([0.,0.,-1.5e11])]:
        np.testing.assert_allclose(moments(q,f,s,every=True),parallel_moments(q,f,s),rtol=2e-7,atol=.02)


def test_contact_parallel_formula_and_transverse_intersection():
    f=np.eye(3);rng=np.random.default_rng(4298)
    for d in rng.uniform(-18000,18000,(25,3)):
        expected=np.linalg.norm(np.maximum(abs(d)-[10000,10000,0],0))
        assert abs(square_distance(np.zeros(3),f,d,f)-expected)<1e-8
    f2=Rotation.from_euler('x',90,degrees=True).as_matrix()
    assert square_distance(np.zeros(3),f,np.zeros(3),f2)<1e-8
    assert abs(square_distance(np.zeros(3),f,np.array([0,0,7000]),f2)-2000)<1e-8


def test_tilted_shadow_against_direct_independent_ray_intersections():
    # A crossing blocker invalidates centre-only depth ordering. Compare the
    # cone union against direct line-plane intersections at deterministic rays.
    q=np.array([[0.,0.,0.],[1200,0,2000.]])
    f=np.array([np.eye(3),Rotation.from_euler('y',60,degrees=True).as_matrix()])
    source=np.array([0.,0.,1.5e11]);m=moments(q,f,source,every=True)
    x=(np.arange(1200)+.5)/1200*9890-9890/2;xx,yy=np.meshgrid(x,x)
    p=np.c_[xx.ravel(),yy.ravel(),np.zeros(xx.size)]
    ray=source-p;n=f[1,:,2];den=ray@n
    tau=((q[1]-p)@n)/den;hit=p+tau[:,None]*ray
    local=(hit-q[1])@f[1]
    blocked=(tau>0)&(tau<1)&(abs(local[:,0])<4945)&(abs(local[:,1])<4945)
    visible=~blocked
    assert abs(m[0,0]/9890**2-visible.mean())<.002
    sampled=np.mean(p[visible,:2],axis=0)*visible.mean()*9890**2
    np.testing.assert_allclose(m[0,1:],sampled,rtol=.006,atol=1e6)
