"""Independent source parity and analytical controls for the Earth benchmark."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import welss_column as w

HERE = Path(w.__file__).resolve().parent
REFERENCE = json.loads((HERE/'references/welss_microphysics_reference.json').read_text())


def test_independent_fixture_identity():
    assert REFERENCE['source_revision'] == w.REVISION
    for name, digest in REFERENCE['source_files'].items():
        assert hashlib.sha256((HERE/'references'/name).read_bytes()).hexdigest() == digest
    assert hashlib.sha256((HERE/'references/build_welss_reference.py').read_bytes()).hexdigest() == REFERENCE['generation_script_sha256']


@pytest.mark.parametrize('case',REFERENCE['cases'])
def test_independent_fortran_microphysics(case):
    T,p = case['temperature_k'],case['pressure_pa']
    ice,water = w.saturation(T)
    pv = case['relative_humidity_ice']*ice
    state = dict(T=T,p=p,pv=pv,psati=ice,psatw=water,eta=w.viscosity(T),
                 rho=(p-pv)/(w.SOURCE_CONSTANTS['Rd']*T)+pv/(w.SOURCE_CONSTANTS['Rv']*T))
    expected=case['expected']
    for actual,key in [(ice,'ice_saturation_pressure_pa'),(water,'water_saturation_pressure_pa'),
                       (state['rho'],'air_density_kg_m3'),(state['eta'],'dynamic_viscosity_pa_s')]:
        np.testing.assert_allclose(actual,expected[key],rtol=3e-14,atol=0)
    result=w.microphysics(case['mass_kg'],state)
    mapping=dict(D='diameter_m',A='projected_area_m2',C='capacitance_m',v='terminal_speed_m_s',
                 Re='Re',growth='deposition_kg_s',Dv='vapour_diffusivity_m2_s',Ka='thermal_conductivity_w_m_k',
                 L='latent_heat_j_kg',fv='vapour_ventilation_factor',fh='heat_ventilation_factor')
    for key,refkey in mapping.items():
        np.testing.assert_allclose(result[key],expected[refkey],rtol=2e-9,atol=1e-30,err_msg=key)


def test_gamma_and_source_normalization():
    m,mw,z,zw=w.seed_quadrature(12,3)
    assert mw.sum() == pytest.approx(1.)
    assert mw@m == pytest.approx(w.SEED_MEAN,rel=1e-13,abs=0)
    assert mw@(m*m) == pytest.approx(1.5*w.SEED_MEAN**2,rel=1e-13,abs=0)
    assert zw.sum() == pytest.approx(300.)
    assert zw@z/zw.sum() == pytest.approx(2850.)
    assert zw.sum()*w.SOURCE_RATE == pytest.approx(4950.)


def test_current_source_atmosphere_conditions():
    z=np.array([0.,40.,590.,600.,610.,1000.,2700.,3000.,5000.])
    state=w.atmosphere(z)
    np.testing.assert_allclose(state['T'],273.15-.0062*z,rtol=0,atol=6e-14)
    assert np.all(state['qc']>=0)
    assert np.all(state['pv'] <= state['psatw']*(1+2e-15))
    np.testing.assert_allclose(state['pv'][z>=1000]/state['psati'][z>=1000],1.05,rtol=4e-5)
    np.testing.assert_allclose(state['pv'][z<=600]/state['psatw'][z<=600],1.,rtol=2e-15)


def test_arrival_windows_are_distinct():
    ages=np.array([0.,18000.,27000.,34201.,36000.,36001.])
    np.testing.assert_allclose(w.arrival_weight(ages,'last5h'),[1,1,.5,1799/18000,0,0])
    np.testing.assert_allclose(w.arrival_weight(ages,'four_snapshots'),[1,1,1,.75,.25,0])


def test_characteristics_against_constant_growth_and_settling():
    masses=np.array([1e-12,2e-12])
    G,v=3e-13,2.
    sol=w.trajectory(200.,masses,rates=lambda z,m:dict(v=np.full_like(m,v),growth=np.full_like(m,G)),
                     max_step_m=50.,rtol=1e-10)
    z=np.array([200.,150.,0.]);state=sol.sol(z)
    age=(200-z)/v
    np.testing.assert_allclose(state[2:],np.broadcast_to(age,(2,3)),rtol=2e-13,atol=1e-12)
    np.testing.assert_allclose(np.exp(state[:2]),masses[:,None]+G*age,rtol=1e-9)


def test_binned_source_moments_against_exact_steady_solution(monkeypatch):
    v=2.
    def rates(m,state):
        return dict(v=np.full_like(m,v),growth=np.zeros_like(m))
    monkeypatch.setattr(w,'microphysics',rates)
    result=w.integrate(mass_nodes=2,height_nodes_per_cell=1,bin_nodes=2,max_step_m=3000.,wall_limit_s=5.)
    z=result['height']
    flux=w.SOURCE_RATE*np.clip(w.SOURCE_TOP-(z+w.CELL_WIDTH/2),0.,w.SOURCE_TOP-w.SOURCE_BOTTOM)
    expected=np.asarray([flux,flux*w.SEED_MEAN,flux/v,flux*w.SEED_MEAN/v,np.zeros_like(flux)])
    for profile in result['profiles'].values():
        np.testing.assert_allclose(profile,expected,rtol=3e-14,atol=1e-25)


def test_filter_masks_edges_and_preserves_constant_interior():
    x=np.ones(20)
    result=w.lanczos5(x)
    assert np.isnan(result[:2]).all() and np.isnan(result[-2:]).all()
    np.testing.assert_allclose(result[2:-2],1.,rtol=1e-7)
    x[10]=np.nan
    result=w.lanczos5(x)
    assert np.isnan(result[8:13]).all()


def test_reference_hash_mismatch_is_rejected(tmp_path):
    path=tmp_path/'ob_hp0.nc';path.write_bytes(b'not the author reference')
    with pytest.raises(ValueError,match='hash mismatch'):
        w.read_archive(path)


FILTER_REFERENCE = json.loads((HERE/'references/welss_filter_reference.json').read_text())


def test_filter_fixture_identity():
    actual = hashlib.sha256((HERE/'references/build_welss_filter_reference.py').read_bytes()).hexdigest()
    assert actual == FILTER_REFERENCE['generation_script_sha256']


@pytest.mark.parametrize('case',FILTER_REFERENCE['cases'],ids=lambda c:c['name'])
def test_independently_compiled_ncl_filter(case):
    values=np.asarray([np.nan if x is None else x for x in case['input']])
    expected=np.asarray([np.nan if x is None else x for x in case['expected']])
    actual=w.lanczos5(values)
    np.testing.assert_array_equal(np.isnan(actual),np.isnan(expected))
    np.testing.assert_allclose(actual,expected,rtol=5e-14,atol=1e-15,equal_nan=True)


def test_arrival_primitives_across_window_boundaries():
    age=np.array([0.,18000.,27000.,34201.,36000.,36001.])
    np.testing.assert_allclose(w.arrival_primitive(age,'last5h'),
                               [0.,18000.,24750.,27000.-1799.**2/36000.,27000.,27000.])
    np.testing.assert_allclose(w.arrival_primitive(age,'four_snapshots'),
                               [0.,18000.,27000.,34200.75,35100.,35100.])
    np.testing.assert_allclose(w.arrival_primitive(age,'final'),
                               [0.,18000.,27000.,34201.,36000.,36000.])


def test_exact_number_density_for_slow_variable_speed(monkeypatch):
    # An analytically integrable nonuniform speed spans all finite-age windows;
    # one spatial quadrature node must suffice for the endpoint-age Nd formula.
    z0=193.;v0=.001;scale=100.
    masses=np.array([w.SEED_MEAN/2,1.5*w.SEED_MEAN])
    monkeypatch.setattr(w,'seed_quadrature',lambda *args:(masses,np.array([.5,.5]),np.array([z0]),np.array([1.])))
    def rates(m,state):
        z=(w.SOURCE_CONSTANTS['T_3']-state['T'])/w.SOURCE_CONSTANTS['tlr']
        return dict(v=np.broadcast_to(v0*(1+(z/scale)**2),m.shape),growth=np.zeros_like(m))
    monkeypatch.setattr(w,'microphysics',rates)
    result=w.integrate(mass_nodes=2,height_nodes_per_cell=1,bin_nodes=1,max_step_m=5.,rtol=1e-10,wall_limit_s=5.)
    lower=np.arange(0.,z0,w.CELL_WIDTH);upper=np.minimum(lower+w.CELL_WIDTH,z0)
    age=lambda z:scale/v0*(np.arctan(z0/scale)-np.arctan(z/scale))
    for mode,profile in result['profiles'].items():
        expected=w.SOURCE_RATE/w.CELL_WIDTH*(w.arrival_primitive(age(lower),mode)-w.arrival_primitive(age(upper),mode))
        np.testing.assert_allclose(profile[2,:len(lower)],expected,rtol=2e-10,atol=1e-9)
        assert np.all(profile[2,len(lower):]==0.)
