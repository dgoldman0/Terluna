"""Conservation and geometry checks on the compiled native outcome routine."""
from pathlib import Path
import hashlib
import numpy as np
import pytest
from atmosphere.snow import scale_separation as g


@pytest.fixture
def native_source():
    path=g.CACHE/'lunar/work'/g.SOURCE
    if not path.exists():pytest.skip('Optional pinned SCALE source cache')
    return path.read_text()


def test_patch_rejects_unidentified_source():
    with pytest.raises(ValueError):g.patch('unidentified source',.5)


@pytest.mark.parametrize('value',[-.01,1.01,float('nan')])
def test_invalid_separation(value):
    with pytest.raises(ValueError):g.wrapper(value)


def test_native_zero_limit_mass_conservation_and_planar_span(tmp_path,native_source):
    rng=np.random.default_rng(51022)
    rows=[]
    for _ in range(160):
        a1,a2=np.exp(rng.uniform(np.log(5e-5),np.log(5e-3),2))
        p1,p2=rng.uniform(.02,.3,2);r1,r2=rng.uniform(50,800,2)
        rows.append([1,a1,a2,a1*p1,a2*p2,r1,r2])
    rows=np.array(rows)
    mass0=4*np.pi/3*(rows[:,1]**2*rows[:,3]*rows[:,5]+rows[:,2]**2*rows[:,4]*rows[:,6])
    for s in (0.,.2,.5,.65):
        out=g.pair_outcomes(rows,s,tmp_path/str(s))
        if s==0:np.testing.assert_array_equal(out[:,:3],out[:,3:])
        expected=np.maximum.reduce([rows[:,1],rows[:,2],(1+s)*(rows[:,1]+rows[:,2])/2])
        np.testing.assert_allclose(out[:,3],expected,rtol=2e-15)
        mass1=4*np.pi/3*out[:,3]**2*out[:,4]*out[:,5]
        np.testing.assert_allclose(mass1,mass0,rtol=3e-15)
        assert np.all(out[:,3:] > 0)
        assert np.all(out[:,5] <= 916.8)  # pinned native ice density, not a new lunar constant
        assert np.all(2*out[:,3] <= 2*(rows[:,1]+rows[:,2])*(1+1e-15))


def test_columnar_encounters_preserve_native_outcome(tmp_path,native_source):
    rows=[[n,.0005,.0002,.001,.00002,200,300] for n in (1,2,7)]
    rows += [[n,.0002,.0005,.00002,.001,300,200] for n in (1,2,7)]
    out=g.pair_outcomes(rows,.5,tmp_path)
    np.testing.assert_array_equal(out[:,:3],out[:,3:])


def test_multiple_planar_events_match_sequential_single_events(tmp_path,native_source):
    row=[5,.0005,.0004,.00005,.00004,300,250]
    batch=g.pair_outcomes([row],.5,tmp_path/'batch')[0,3:]
    state=[row[2],row[4],row[6]]
    for event in range(row[0]):
        one=[1,row[1],state[0],row[3],state[1],row[5],state[2]]
        state=g.pair_outcomes([one],.5,tmp_path/str(event))[0,3:]
    np.testing.assert_array_equal(batch,state)
    expected=4*np.pi/3*(row[0]*row[1]**2*row[3]*row[5]+row[2]**2*row[4]*row[6])
    assert 4*np.pi/3*batch[0]**2*batch[1]*batch[2]==pytest.approx(expected,rel=3e-15)


def test_patch_preserves_original_outcome_and_installs_counters(native_source):
    modified=g.patch(native_source,.5)
    original=g.native_routine(native_source)
    assert original.replace(g.ROUTINE,g.ROUTINE+'_native') in modified
    assert modified.count('subroutine '+g.ROUTINE+'(')==1
    assert 'separation_observer.txt' in modified
    assert hashlib.sha256(native_source.encode()).hexdigest()==g.SOURCE_SHA
