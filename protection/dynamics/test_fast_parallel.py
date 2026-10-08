import numpy as np
from scipy.spatial.transform import Rotation
from .fast_parallel import moments
from .tile_torque import illuminated_moments
from .parallel_shadow import source_illumination


def test_independent_all_pair_area_and_moments():
    rng=np.random.default_rng(845213)
    q=rng.uniform(-17000,17000,(25,3))
    for frame in [np.eye(3),Rotation.from_euler('xyz',[24,62,-34],degrees=True).as_matrix()]:
        for source in [np.array([0.,0.,1.5e11]),np.array([1.5e11,2e8,-3e8])]:
            old=illuminated_moments(q,frame,source)
            new=moments(q,frame,source)
            np.testing.assert_allclose(new[:,0],old[:,0],rtol=1e-8,atol=1.)
            np.testing.assert_allclose(new[:,1:],old[:,1:],rtol=2e-7,atol=5000.)
            np.testing.assert_allclose(new[:,0]/9890**2,source_illumination(q,frame,source,every=True),atol=1e-8)


def test_stack_no_double_count_and_single_off_centre_shadow():
    q=np.array([[0.,0.,0.],[0,0,12000],[0,0,24000.]])
    m=moments(q,np.eye(3),np.array([0.,0.,1.5e11]))
    np.testing.assert_allclose(m[:,0],[0,0,9890**2],atol=1e-6)
    q=np.array([[0.,0.,0.],[7000.,0.,12000.]])
    m=moments(q,np.eye(3),np.array([0.,0.,1.5e11]))
    assert m[0,1]<0 and abs(m[0,2])<1e-3
