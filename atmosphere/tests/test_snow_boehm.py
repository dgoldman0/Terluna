"""Independent McSnow source parity and bounded adapter checks."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import boehm

HERE = Path(boehm.__file__).resolve().parent
REFERENCE = json.loads((HERE/'results/boehm_reference.json').read_text())


def _fixture_arrays():
    rows = REFERENCE['cases']
    keys = ('D1','m1','q1','sigma1','D2','m2','q2','sigma2','g','rho','eta')
    return {key:np.array([row[key] for row in rows]) for key in keys}


def _fixtures():
    a = _fixture_arrays()
    state = dict(rho=a['rho'], eta=a['eta'])
    left = boehm.particle_properties(a['D1'],a['m1'],a['g'],state,a['q1'],sigma=a['sigma1'])
    right = boehm.particle_properties(a['D2'],a['m2'],a['g'],state,a['q2'],sigma=a['sigma2'])
    return left, right, state


def test_reference_identity():
    assert REFERENCE['source_revision'] == boehm.SOURCE_REVISION
    for name, expected in REFERENCE['source_slices'].items():
        assert hashlib.sha256((HERE/'references'/name).read_bytes()).hexdigest() == expected
    for name, expected in REFERENCE['additional_inputs_sha256'].items():
        assert hashlib.sha256((HERE/'references'/name).read_bytes()).hexdigest() == expected
    assert hashlib.sha256((HERE/'references/build_boehm_reference.py').read_bytes()).hexdigest() == REFERENCE['generation_script_sha256']


def test_broadcast_particle_and_pair_results_match_independent_fortran():
    left, right, state = _fixtures()
    for suffix, properties in [('1',left),('2',right)]:
        for key, tolerance in [('v',2e-12),('Re',2e-12),('delta',2e-12),('CDI',5e-10),('alpha_e',2e-14)]:
            expected = np.array([row['expected'][key+suffix] for row in REFERENCE['cases']])
            np.testing.assert_allclose(properties[key],expected,rtol=tolerance,atol=1e-13)
    expected = np.array([row['expected']['Ec'] for row in REFERENCE['cases']])
    actual = boehm.collision_efficiency(left,right,state,diagnostics=True)
    np.testing.assert_allclose(actual['Ec'],expected,rtol=2e-9,atol=5e-13)
    # The independent fixtures exercise all piecewise physical approximations.
    assert (actual['F_over_H'] < 10).any() and (actual['F_over_H'] >= 10).any()
    assert (actual['b_cy'] < .5).any()
    assert ((actual['b_cy'] > .5) & (actual['b_cy'] < 1)).any()
    assert (actual['b_cy'] >= 1).any()


def test_scalar_interface_matches_fortran_and_array_interface():
    left, right, state = _fixtures()
    for index in (0,1,16,47,81,len(REFERENCE['cases'])-1):
        row = REFERENCE['cases'][index]
        air = dict(rho=row['rho'],eta=row['eta'])
        one = boehm.particle_properties(row['D1'],row['m1'],row['g'],air,row['q1'],sigma=row['sigma1'])
        two = boehm.particle_properties(row['D2'],row['m2'],row['g'],air,row['q2'],sigma=row['sigma2'])
        assert all(isinstance(v,np.ndarray) and v.ndim == 0 for v in one.values())
        assert boehm.collision_efficiency(one,two,air) == pytest.approx(row['expected']['Ec'],rel=2e-9,abs=5e-13)
        for key in one:
            assert one[key] == pytest.approx(left[key][index],rel=1e-14)


def test_pair_matrix_is_symmetric_and_preserves_source_anomalies():
    state = dict(rho=.5,eta=1.7e-5)
    D = np.geomspace(2e-4,.6,40)
    properties = boehm.particle_properties(D,.0121*D**1.9,1.62,state,.35)
    left = {key:value[:,None] for key,value in properties.items()}
    right = {key:value[None,:] for key,value in properties.items()}
    diagnostics = boehm.collision_efficiency(left,right,state,diagnostics=True)
    assert diagnostics['Ec'].shape == (40,40)
    np.testing.assert_array_equal(diagnostics['Ec'],diagnostics['Ec'].T)
    assert diagnostics['Ec'].max() > 1
    np.testing.assert_array_equal(diagnostics['outside_unit_interval'],
                                  (diagnostics['Ec']<0)|(diagnostics['Ec']>1))
    kernel = np.pi/4*(D[:,None]+D[None,:])**2 * abs(properties['v'][:,None]-properties['v'][None,:]) * diagnostics['Ec']
    np.testing.assert_array_equal(kernel.diagonal(),np.zeros(40))
    np.testing.assert_array_equal(kernel,kernel.T)


def test_equal_speed_sigma_zero_retains_fortran_sign_convention():
    state = dict(rho=1.,eta=1.7e-5)
    particle = boehm.particle_properties(.001,2e-8,1.62,state,.35,sigma=0)
    result = boehm.collision_efficiency(particle,particle,state,diagnostics=True)
    assert result['vstar'] == pytest.approx(2*particle['v'])
    broad = boehm.particle_properties(.001,2e-8,1.62,state,.35,sigma=.25)
    result_broad = boehm.collision_efficiency(broad,broad,state,diagnostics=True)
    assert result_broad['vstar'] < result['vstar']


def test_archived_earth_graupel_ground_speeds_under_source_default_surface_state():
    # Original run namelist is unavailable. This checks the explicitly stated
    # source-default saturated surface against five rounded archived outputs.
    reference = json.loads((HERE/'references/mcsnow_ground_graupel.json').read_text())
    rows = reference['rows']
    a = {k:np.array([r[k] for r in rows]) for k in rows[0]}
    state = dict(rho=a['rho_kg_m3'],eta=a['eta_Pa_s'])
    actual = boehm.particle_properties(a['D_m'],a['m_kg'],9.81,state,1)
    np.testing.assert_allclose(actual['v'],a['v_author_m_s'],rtol=1e-6,atol=0)
    np.testing.assert_allclose(actual['Re'],a['Re_author'],rtol=1e-6,atol=0)


@pytest.mark.parametrize('changes',[
    {'D':0}, {'m':-1}, {'g':0}, {'q':0}, {'q':1.1}, {'alpha':.8},
    {'sigma':-1}, {'sigma':np.nan}, {'state':{'rho':0,'eta':1e-5}},
])
def test_rejects_unsupported_or_invalid_inputs(changes):
    kwargs = dict(D=.001,m=1e-8,g=1.62,state={'rho':1.,'eta':1.7e-5},q=.35)
    kwargs.update(changes)
    with pytest.raises(ValueError):
        boehm.particle_properties(**kwargs)


def test_collision_rejects_mixed_air_states():
    state = dict(rho=1.,eta=1.7e-5)
    particle = boehm.particle_properties(.001,1e-8,1.62,state,.35)
    with pytest.raises(ValueError,match='density'):
        boehm.collision_efficiency(particle,particle,dict(rho=.5,eta=state['eta']))
