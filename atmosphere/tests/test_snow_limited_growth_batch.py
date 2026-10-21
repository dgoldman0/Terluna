"""Bitwise parity of batched reconstruction with frozen per-layer operators."""
import numpy as np
import pytest

from atmosphere.snow.limited_growth_batch import BatchLimitedGrowth
from atmosphere.snow.limited_growth import LimitedGrowth


@pytest.mark.parametrize('bins',[2,3,17,96,144])
def test_batch_matches_scalar_flux_water_and_cfl_bitwise(bins):
    rng=np.random.default_rng(74+bins)
    masses=np.exp(np.cumsum(rng.uniform(.01,.5,bins)))
    growth=rng.uniform(0.,3.,(9,bins))
    growth[0]=0.
    ops=[LimitedGrowth(masses,g,cfl=.1+.05*k) for k,g in enumerate(growth)]
    n=np.exp(rng.uniform(-40.,10.,(9,bins)))
    n[0]=0.;n[1,::2]=0.;n[2]=1e-310;n[3,-1]=np.nextafter(0.,1.)
    production=np.array([0.,0.,1e-312,1e-100,.001,1.,10.,1e8,1e20])
    expected=[o.evaluate(row,p) for o,row,p in zip(ops,n,production)]
    flux,water,limit=BatchLimitedGrowth(ops).evaluate(n,production)
    np.testing.assert_array_equal(flux.view(np.uint64),np.stack([x[0] for x in expected]).view(np.uint64))
    np.testing.assert_array_equal(water.view(np.uint64),np.array([x[1] for x in expected]).view(np.uint64))
    assert limit.hex()==float(min(x[2] for x in expected)).hex()
    assert sum(map(float,water))==sum(x[1] for x in expected)


def test_scalar_production_broadcast_and_empty_population_match():
    masses=np.geomspace(1e-10,1.,96)
    ops=[LimitedGrowth(masses,np.full(96,.1)) for _ in range(4)]
    n=np.zeros((4,96));n[2,20:24]=[1.,2.,3.,4.]
    batch=BatchLimitedGrowth(ops)
    a=batch.evaluate(n,.02);b=batch.evaluate(n,np.full(4,.02))
    for left,right in zip(a[:2],b[:2]):np.testing.assert_array_equal(left,right)
    assert a[2]==b[2]
    zeros,water,limit=batch.evaluate(np.zeros_like(n),.02)
    assert np.count_nonzero(zeros)==0 and np.count_nonzero(water)==0
    assert np.isinf(limit)


def test_batch_rejects_incompatible_grids_and_shapes():
    a=LimitedGrowth([1.,2.],[1.,1.]);b=LimitedGrowth([1.,3.],[1.,1.])
    with pytest.raises(ValueError,match='same mass grid'):BatchLimitedGrowth([a,b])
    with pytest.raises(ValueError,match='nonempty'):BatchLimitedGrowth([])
    batch=BatchLimitedGrowth([a])
    with pytest.raises(ValueError,match='shape'):batch.evaluate([1.,2.],1.)
    with pytest.raises(ValueError,match='production'):batch.evaluate([[1.,2.]],[-1.])
