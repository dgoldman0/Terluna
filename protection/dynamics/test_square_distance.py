import numpy as np
from scipy.spatial.transform import Rotation,RotationSpline
from shapely.geometry import box

from .square_distance import closest_squares,distance_guard


def test_candidate_search_matches_all_pairs_independent_polygon_distances():
    rng=np.random.default_rng(721);q=rng.normal(size=(30,3))*12000.
    frame=Rotation.from_rotvec([.2,-.6,.9]).as_matrix();xyz=q@frame
    polygons=[box(x-5000,y-5000,x+5000,y+5000) for x,y,z in xyz]
    values=[(np.hypot(polygons[i].distance(polygons[j]),xyz[i,2]-xyz[j,2]),(i,j))
        for i in range(len(q)) for j in range(i+1,len(q))]
    expected,pair=min(values)
    actual,found,_=closest_squares(q,frame)
    assert abs(actual-expected)<1e-8
    assert tuple(found)==pair


def test_interval_guard_detects_a_plane_crossing_between_clear_endpoints():
    t=np.array([0.,1000.]);states=np.zeros((2,2,6))
    states[:,1,2]=[-3000.,3000.];states[:,1,5]=6.
    command=RotationSpline(t,Rotation.from_matrix(np.tile(np.eye(3),(2,1,1))))
    result=distance_guard(t,states,command,acceleration_bound=0.,rate_bound=0.)
    assert result['minimum_sampled_surface_distance_m']==3000.
    assert result['minimum_conditional_all_pair_clearance_m']==0.
    assert result['intervals_below_100m']==1
