"""Independent author-source fixtures and scope checks for Earth aggregation."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import welss_column as w
from atmosphere.snow import welss_aggregation_physics as a

ROOT=Path(__file__).resolve().parents[2]
FIXTURE=ROOT/'atmosphere/snow/references/welss_aggregation_reference.json'


def air(row):
    t,p,rh=row['temperature_k'],row['pressure_pa'],row['relative_humidity_ice']
    ice,water=w.saturation(t);pv=rh*ice;c=w.SOURCE_CONSTANTS
    return dict(T=t,p=p,pv=pv,psati=ice,psatw=water,
                rho=(p-pv)/(c['Rd']*t)+pv/(c['Rv']*t),eta=w.viscosity(t))


def test_independent_author_components_and_combined_kernel():
    fixture=json.loads(FIXTURE.read_text())
    assert len(fixture['cases'])==182
    for path,expected in fixture['source_files'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    for row in fixture['cases']:
        state=air(row);expected=row['expected'];actual={};velocities=[]
        for side in ('left','right'):
            result=a.microphysics(row[f'mass_{side}_kg'],row[f'monomers_{side}']<1.01,state)
            velocities.append(result['v'])
            for key,unit,out in [('diameter','m','D'),('capacitance','m','C'),
                                 ('deposition','kg_s','growth'),('speed','m_s','v')]:
                actual[f'{key}_{side}_{unit}']=result[out]
            actual[f'Re_{side}']=result['Re']
        result=a.kernel(row['mass_left_kg'],row['mass_right_kg'],state)
        actual.update(interception=result['Ec'],sticking=result['sticking'])
        for key,value in actual.items():
            np.testing.assert_allclose(value,expected[key],rtol=3e-12,atol=0,err_msg=key)
        assert expected['SIG_left']==expected['SIG_right']==.25
        # The outer velocity subtraction is ill-conditioned for nearly equal
        # speeds. Bound its absolute rounding error by the unsubtracted scale.
        scale=result['encounter_area_m2']*result['Ec']*result['sticking']*sum(velocities)
        bound=32*np.finfo(float).eps*scale+3e-12*abs(expected['kernel_m3_s'])
        assert abs(result['K']-expected['kernel_m3_s'])<=bound


def test_monomers_preserve_the_frozen_deposition_reference():
    masses=np.geomspace(1e-14,1e-5,31);state=w.atmosphere(1500.)
    actual=a.microphysics(masses,True,state);baseline=w.microphysics(masses,state)
    for key in baseline:np.testing.assert_array_equal(actual[key],baseline[key])


def test_aggregate_history_changes_growth_without_changing_diagnostic_geometry():
    masses=np.array([1e-12,1e-9,1e-7]);state=w.atmosphere(1800.)
    mono=a.microphysics(masses,True,state);agg=a.microphysics(masses,False,state)
    for key in ('D','A','q','void','v','Re'):np.testing.assert_array_equal(mono[key],agg[key])
    np.testing.assert_allclose(agg['growth']/mono['growth'],agg['C']/mono['C'],rtol=1e-15)
    porous=mono['void']>0
    np.testing.assert_allclose(agg['C'][porous]/mono['C'][porous],w.PI/4,rtol=2e-16)


def test_collision_kernel_is_reciprocal_and_zero_for_equal_terminal_speeds():
    masses=np.geomspace(1e-11,1e-5,9);state=w.atmosphere(1200.)
    result=a.kernel(masses[:,None],masses[None,:],state)
    np.testing.assert_allclose(result['K'],result['K'].T,rtol=1e-14,atol=0)
    np.testing.assert_array_equal(np.diag(result['K']),0.)
    assert np.isfinite(result['K']).all() and (result['K']>=0).all()
    assert (result['Ec']>1).any() # source cross-section ratio, never a clipped probability


def test_scope_rejects_ambiguous_history_and_invalid_inputs():
    with pytest.raises(ValueError,match='boolean'):a.geometry(1e-9,2)
    with pytest.raises(ValueError,match='positive'):a.geometry(-1.,True)
    with pytest.raises(ValueError,match='temperature'):a.sticking(np.nan)
    with pytest.raises(ValueError,match='temperature'):a.sticking(0.)
