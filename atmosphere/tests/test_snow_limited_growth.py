"""Convergence, positivity and two conservation laws for limited mass growth."""
import json

import numpy as np
import pytest

from atmosphere.snow import limited_growth as lg
from shared.provenance import constants_changed
from shared.baselines import recorded_file


def _primitive(x):
    # Integral of a unit-normalized compact sin-squared pulse on [2,4].
    x=np.clip(np.asarray(x)-2.,0.,2.)
    return x/2.-np.sin(np.pi*x)/(2*np.pi)


@pytest.mark.parametrize('kind', ['uniform','geometric'])
def test_smooth_translation_approaches_second_order_on_both_grids(kind):
    errors=[]; moment_errors=[]
    for bins in (64,128,256):
        masses=np.linspace(1.,8.,bins) if kind=='uniform' else np.geomspace(1.,8.,bins)
        grid=lg.prepare(masses)
        initial=_primitive(grid['edges'][1:])-_primitive(grid['edges'][:-1])
        velocity=.5
        result=lg.integrate(lg.LimitedGrowth(masses,np.full(bins,velocity)),initial,[0.,1.],100.,
                            max_step_s=.15*grid['widths'].min()/velocity)
        exact=_primitive(grid['edges'][1:]-.5)-_primitive(grid['edges'][:-1]-.5)
        errors.append(np.abs(result['numbers'][-1]-exact).sum())
        moment_errors.append(abs(result['deposited_mass'][-1]-.5))
        assert (result['numbers']>=0).all()
    assert np.log2(errors[-2]/errors[-1])>1.8
    assert errors[-1]<.0015
    assert moment_errors[-1]<moment_errors[0]/20


def test_nonuniform_cells_reconstruct_about_true_centroids_and_bound_faces():
    masses=np.geomspace(1.,100.,31)
    grid=lg.prepare(masses)
    assert not np.allclose(grid['centroids'],masses)
    operator=lg.LimitedGrowth(masses,np.ones_like(masses))
    # A linear density should reconstruct exactly where slopes have neighbors.
    density=1.+.02*grid['centroids']
    faces=operator.face_density(density*grid['widths'])
    np.testing.assert_allclose(faces[1:-1],1.+.02*grid['face_masses'][1:-1],rtol=1e-14)
    pulse=np.zeros_like(masses);pulse[10:13]=[1.,4.,1.]
    faces=operator.face_density(pulse)
    assert (faces>=0).all()
    assert (faces[pulse==0]==0).all()


def test_flux_cap_closes_number_and_discrete_mass_ledger_including_ghost():
    masses=np.geomspace(1.,32.,10)
    operator=lg.LimitedGrowth(masses,np.linspace(.1,2.,10))
    number=np.linspace(.1,1.,10)
    flux,deposited,limit=operator.evaluate(number,.01)
    dn=-flux.copy();dn[1:]+=flux[:-1]
    assert deposited==pytest.approx(.01,rel=1e-14)
    assert dn.sum()+flux[-1]==pytest.approx(0.,abs=1e-16)
    assert np.dot(masses,dn)+operator.upper_ghost_mass*flux[-1]==pytest.approx(deposited,rel=1e-13)
    assert (number+limit*dn>=0).all()
    dry,water,_=operator.evaluate(number,0.)
    np.testing.assert_array_equal(dry,np.zeros(10));assert water==0


def test_sharp_distribution_remains_positive_without_population_clipping():
    masses=np.geomspace(1.,100.,64)
    grid=lg.prepare(masses); initial=np.zeros(64);initial[12]=1.;initial[60]=1e-100
    result=lg.integrate(lg.LimitedGrowth(masses,.3*grid['face_masses']),initial,[0.,2.],10.,max_step_s=1.)
    final=result['numbers'][-1]
    assert (final>=0).all()
    assert initial.sum()-final.sum()-result['upper_exit_number'][-1]==pytest.approx(0.,abs=2e-14)
    source=np.dot(initial,masses)+result['deposited_mass'][-1]
    assert source-np.dot(final,masses)-result['upper_exit_mass'][-1]==pytest.approx(0.,abs=2e-13)
    assert final[-1]>0


def test_rejects_invalid_flux_inputs():
    with pytest.raises(ValueError,match='face growth'):
        lg.LimitedGrowth([1.,2.],[-1.,1.])
    operator=lg.LimitedGrowth([1.,2.],[1.,1.])
    with pytest.raises(ValueError,match='population'):
        operator.evaluate([1.,-1.],1.)
    with pytest.raises(ValueError,match='production'):
        operator.evaluate([1.,1.],-1.)


def test_limited_product_meets_transport_gates_and_binds_inputs():
    if not lg.OUT.exists():
        pytest.skip('limited comparison has not been generated')
    product=json.loads(lg.OUT.read_text())
    assert product['schema']==lg.SCHEMA
    for kind in ('files','inputs'):
        for path,expected in product['producer'][kind].items():
            assert lg.gt.sha(recorded_file(path,expected))==expected,path
    assert not constants_changed(product['producer']['constants'])
    assert all(product['acceptance'].values())
    for run in product['runs']:
        assert (np.asarray(run['final_numbers_m3'])>=0).all()
        for row in run['snapshots']:
            ledger=row['water_ledger']; source=ledger['initial']+ledger['deposited']
            residual=source-ledger['remaining']-ledger['upper_exit']
            assert abs(residual)/source<1e-10
            assert abs(row['number_residual'])/row['number_m3']<1e-10
