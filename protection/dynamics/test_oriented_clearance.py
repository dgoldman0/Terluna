import numpy as np
from scipy.spatial.transform import Rotation
from .oriented_clearance import closest_oriented,projection_lower_bounds
from .oriented_tiles import square_distance


def test_projection_pruning_preserves_global_oriented_minimum():
    rng=np.random.default_rng(701)
    for _ in range(5):
        q=rng.uniform(-15000,15000,(6,3));f=Rotation.random(6,random_state=rng).as_matrix()
        exact=min(square_distance(q[i],f[i],q[j],f[j]) for i in range(6) for j in range(i))
        actual,_=closest_oriented(q,f)
        assert abs(actual-exact)<1e-7
        i,j=np.triu_indices(6,1)
        bounds=projection_lower_bounds(q[i]-q[j],f[i],f[j])
        distances=np.array([square_distance(q[a],f[a],q[b],f[b]) for a,b in zip(i,j)])
        assert np.all(bounds<=distances+1e-7)
