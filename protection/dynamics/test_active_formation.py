import numpy as np
import pytest

from shared import constants as K
from .active_formation import lattice,solar_frame,illumination,swept_square_conflicts,coverage_window
from .fleet import sun_points
from .optical import arriving_ray_vectors


def sample():
    zero=np.zeros(3)
    return dict(positions={'sun':np.array([K.AU,0.,0.]),'earth':np.array([0.,K.EARTH_MOON_DISTANCE,0.])},
        sun_v=zero,earth_v=zero,moon_v=zero,sun_a=zero,earth_a=zero,moon_a=zero)


def test_mutual_force_stencil_matches_all_pair_polygon_unions():
    shapely=pytest.importorskip('shapely')
    q,neighbour,valid,subsets=lattice(5)
    rng=np.random.default_rng(334)
    q+=rng.uniform(-100,100,q.shape); q[:,0]+=15e6
    frame=solar_frame(sample()); state=np.zeros((len(q),6)); state[:,:3]=q@frame.T
    actual=illumination(state,frame,sample(),9890.,neighbour,valid,subsets,suns=8)
    sun,_=arriving_ray_vectors(sample()['positions']['sun'],sample()['positions']['earth'],
        sample()['sun_v'],sample()['earth_v'],sample()['sun_a'],sample()['earth_a'])
    h=9890/2; expected=np.zeros(len(q)); receiver=shapely.box(-h,-h,h,h)
    for source in sun_points(sun,8,.317)@frame:
        for i,p in enumerate(q):
            upstream=q[q[:,0]>p[0]]
            mag=(source[0]-p[0])/(source[0]-upstream[:,0])
            projected=source[1:]+(upstream[:,1:]-source[1:])*mag[:,None]-p[1:]
            rectangles=shapely.box(projected[:,0]-h*mag,projected[:,1]-h*mag,
                                   projected[:,0]+h*mag,projected[:,1]+h*mag)
            blocked=shapely.union_all(rectangles).intersection(receiver).area
            expected[i]+=(1-blocked/receiver.area)/8
    np.testing.assert_allclose(actual,expected,atol=2e-10,rtol=0)


def test_depth_stagger_closes_optical_seams_without_material_overlap():
    frame=solar_frame(sample()); q,_,_,_=lattice(25)
    q[:,0]+=15e6
    state=np.zeros((len(q),6)); state[:,:3]=q@frame.T
    centre=np.r_[frame[:,0]*15e6,[0.,0.,0.]]
    cover=coverage_window(state,-frame[:,0],sample(),10000.,centre,suns=16)
    assert cover['ray_coverage']==pytest.approx(1,abs=1e-12)
    assert cover['all_sampled_sun_coverage']==pytest.approx(1,abs=1e-12)
    result=swept_square_conflicts(np.array([0.,60.]),np.array([state,state]),lambda t:frame)
    assert result['possible_conflict_pairs']==0
    assert result['minimum_proven_axis_clearance_m']>800.
    # Removing a neighbouring tile's depth separation creates a physical
    # overlap even though the centre separation exceeds the tile side/2.
    changed=state.copy(); changed[1,:3]-=1000*frame[:,0]
    result=swept_square_conflicts(np.array([0.,60.]),np.array([changed,changed]),lambda t:frame)
    assert result['possible_conflict_pairs']>0


def test_swept_square_guard_catches_an_intersample_crossing():
    states=np.zeros((2,2,6)); states[:,1,0]=[20000.,-20000.]
    result=swept_square_conflicts(np.array([0.,60.]),states,lambda t:np.eye(3))
    assert result['possible_conflict_pairs']==1
