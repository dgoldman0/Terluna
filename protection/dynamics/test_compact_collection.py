from types import SimpleNamespace

import numpy as np
from scipy.spatial.transform import Rotation

from .active_formation import lattice
from .test_cycling import static_sample
from .collection import uneclipsed_collection
from .compact_collection import compact_collection


def test_compact_ledger_matches_original_on_both_faces_and_near_feathering():
    env=SimpleNamespace(at=lambda t:static_sample(),central_fraction=.138)
    q=lattice(5,8500.,12000.)[0]+[15e6,15e6,15e6]
    state=np.c_[q,np.zeros_like(q)]
    for angle in [.2,1.565,1.575,2.2]:
        frame=Rotation.from_rotvec([angle,.12,.4]).as_matrix()
        old=uneclipsed_collection(env,0.,state,frame,suns=8)
        new=compact_collection(env,0.,state,frame,suns=8)
        for key in old:
            np.testing.assert_allclose(old[key],new[key],atol=1e-4,rtol=1e-12)
