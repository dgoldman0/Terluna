"""Compare the compact union against all-pair geometry and the old force law."""
from types import SimpleNamespace

import numpy as np
from scipy.spatial.transform import Rotation

from shared import constants as K
from .active_formation import lattice
from .compact_shadow import compact_illumination,CompactMassivePattern
from .parallel_shadow import source_illumination
from .mass_feedback import MassiveParallelPattern
from .test_cycling import static_sample
from .optical import length


def test_compaction_matches_unpruned_all_pair_unions_including_full_occlusion():
    source=np.array([0.,0.,K.AU])
    for pitch in [8500.,1000.]:
        offsets=lattice(5,pitch,4000.)[0]
        q=offsets[:,[1,2,0]]+np.array([1e6,2e6,15e6])
        for angle in [0.,.8,1.565,1.575,2.2,np.pi]:
            frame=Rotation.from_rotvec([angle,.12,.4]).as_matrix()
            np.testing.assert_allclose(compact_illumination(q,frame,source),
                source_illumination(q,frame,source,every=True),atol=2e-12,rtol=0.)
    q=np.array([[0.,0.,0.],[0.,0.,10000.],[0.,0.,20000.]])
    np.testing.assert_allclose(compact_illumination(q,np.eye(3),source),[0.,0.,1.],atol=1e-12)


def test_compact_force_preserves_each_member_momentum_and_mass_scaling():
    sample=static_sample();env=SimpleNamespace(at=lambda t:sample,sigma=.05,central_fraction=.138,
        gravity=lambda t,q:-K.MOON_GM*q/length(q)[...,None]**3)
    q=lattice(5,8500.,12000.)[0]+np.array([15e6,15e6,15e6]);state=np.c_[q,np.zeros_like(q)]
    ratio=np.linspace(1.2,1.3,25)
    for angle in [.2,1.565,1.575,2.2]:
        command=lambda t,a=angle:Rotation.from_rotvec([a,.12,.4])
        old=MassiveParallelPattern(env,command,ratio,suns=8).acceleration(0.,state,diagnostics=True)
        new=CompactMassivePattern(env,command,ratio,suns=8).acceleration(0.,state,diagnostics=True)
        np.testing.assert_allclose(new[0],old[0],atol=2e-17,rtol=1e-13)
        np.testing.assert_allclose(new[1],old[1],atol=2e-12,rtol=0.)
        np.testing.assert_allclose(new[2],old[2],atol=1e-17,rtol=1e-12)
