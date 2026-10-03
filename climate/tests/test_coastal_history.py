"""Continuous spectral boundaries, incoming-energy transfer and lunar controls."""
from datetime import timedelta
import json
import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from climate.waves.coastal_checks import iter_nesting_spectra, nesting_spectra
from climate.waves.coastal_history import write_boundary
from climate.waves.model import ORIGIN, sha256
from climate.waves.shore import SITES, stationary_boundary
from climate.waves.shore_checks import group_velocity
from climate.waves import runup
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def test_streamed_boundary_keeps_the_clock_dry_points_and_variance(tmp_path):
    xy=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    frequency=np.array([.05,.1,.2,.4])
    pattern=np.array([0.,2.,1.,0.])[:,None]*np.array([1.,0.,.5,0.])[None,:]
    source=dict(locations=xy,frequency=frequency,direction=np.arange(4)*90.,
                available=np.ones((1,4),bool))
    def records():
        for hour in (0,1,2):
            yield dict(source,dates=[ORIGIN+timedelta(hours=hour)],
                       spectra=np.repeat(pattern[None,None],4,axis=1)*hour)
    path=tmp_path/'continuous.nest'
    audit=write_boundary(path,records(),np.array([[.5,.5],[2.,2.]]),np.array([True,False]),1.)
    stream=list(iter_nesting_spectra(path))
    assert audit['records']==3 and len(stream)==3
    assert stream[-1]['dates']==[ORIGIN+timedelta(hours=2)]
    np.testing.assert_allclose(stream[-1]['spectra'][0,0],2*pattern,atol=3e-8)
    assert not stream[-1]['available'][0,1]
    np.testing.assert_array_equal(nesting_spectra(path)['spectra'],np.concatenate([s['spectra'] for s in stream]))
    with pytest.raises(ValueError,match='absent'):
        list(iter_nesting_spectra(path,{ORIGIN+timedelta(hours=3)}))
    path.write_text(path.read_text().rsplit('NODATA',1)[0])
    with pytest.raises(ValueError,match='Truncated'):
        list(iter_nesting_spectra(path))


def test_frequency_integral_is_exact_for_a_triangle_and_zero_outside():
    actual=runup.spectral_integral([1.,2.,3.],[0.,2.,0.],[-1.,1.,1.5,2.,2.5,3.,4.])
    np.testing.assert_allclose(actual,[0.,0.,.25,1.,1.75,2.,2.],atol=1e-14)


def test_swan_to_swash_transfer_conserves_incoming_power_and_rejects_changed_bytes(tmp_path,monkeypatch):
    monkeypatch.setattr(runup,'ROOT',tmp_path)
    profiles=tmp_path/'profiles.json'
    monkeypatch.setattr(runup,'PROFILES',profiles)
    profiles.write_text(json.dumps(dict(schema='terluna.geography.shore-profiles/1',
        profiles={'west_face':dict(normal_east_north=[1.,0.],station_depth_m=5.)})))
    frequency=np.array([.05,.1,.2,.4])
    pattern=np.array([0.,2.,1.,0.])[:,None]*np.array([1.,0.,3.,0.])[None,:]
    source=dict(locations=np.array(list(SITES.values())),frequency=frequency,direction=np.arange(4)*90.,
                spectra=np.repeat(pattern[None,None],4,axis=1),available=np.ones((1,4),bool))
    path=tmp_path/'sites.spc'
    path.write_text(stationary_boundary(source,0))
    settings=tmp_path/'INPUT'
    settings.write_text(f'SET GRAV {MOON_SURFACE_GRAVITY} RHO {runup.WATER_DENSITY}\n')
    manifest=dict(schema='terluna.climate.shore-case/1',run=dict(
        identity=dict(inputs={'INPUT':sha256(settings)}),output_sha256={'sites.spc':sha256(path)}))
    (tmp_path/'product.json').write_text(json.dumps(manifest))
    text,record=runup.incoming_spectrum(path,'west_face',layers=4)
    rows=np.array([[float(v) for v in line.split()] for line in text.splitlines()[1:]])
    df=record['boundary_frequency_step_hz']
    expected_variance=np.trapezoid(pattern[:,0]*90,frequency)
    # The SWAN interchange rounds its spectral integers to eight digits.
    assert np.sum(rows[:,1])*df==pytest.approx(expected_variance,rel=5e-8)
    cg=group_velocity(rows[:,0],5.,MOON_SURFACE_GRAVITY)
    assert np.sum(rows[:,1]*cg)*df==pytest.approx(record['represented_variance_flux_m3_s'],rel=1e-8)
    # Check the reconstructed spectrum against independent quadrature. The
    # deliberately sparse four-bin source retains a separate quadrature error.
    integrand=lambda f:float(group_velocity(np.array([f]),5.)[0])*np.interp(f,frequency,pattern[:,0]*90)
    exact=sum(quad(integrand,a,b,epsabs=1e-10)[0] for a,b in zip(frequency[:-1],frequency[1:]))
    assert record['represented_variance_flux_m3_s']==pytest.approx(exact,rel=1e-5)
    assert record['transport_relative_difference']>.1
    # A deliberately restricted band must be explicit, and identical at both
    # layer counts for a controlled comparison of their boundary formulations.
    with pytest.raises(ValueError,match='discard'):
        runup.incoming_spectrum(path,'west_face',frequency_ceiling=.3)
    band4,restricted=runup.incoming_spectrum(path,'west_face',layers=4,
        frequency_ceiling=.3,max_removed_variance=.2)
    band6,_=runup.incoming_spectrum(path,'west_face',layers=6,
        frequency_ceiling=.3,max_removed_variance=.2)
    assert band4==band6
    assert .01<restricted['removed_variance_fraction']<.2
    with pytest.raises(ValueError,match='zero resolved'):
        runup.incoming_spectrum(path,'west_face',frequency_ceiling=1e-9)
    settings.write_text(f'SET GRAV {STANDARD_GRAVITY} RHO {runup.WATER_DENSITY}\n')
    manifest['run']['identity']['inputs']['INPUT']=sha256(settings)
    (tmp_path/'product.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='gravity or water density'):
        runup.incoming_spectrum(path,'west_face')
    path.write_text(path.read_text()+'\n')
    with pytest.raises(ValueError,match='differs'):
        runup.incoming_spectrum(path,'west_face')


def test_swash_dispersion_and_gravity_similarity(tmp_path,monkeypatch):
    from climate.waves import runup_controls
    path=os.environ.get('TERLUNA_SWASH_EXECUTABLE')
    if not path:
        pytest.skip('Set TERLUNA_SWASH_EXECUTABLE for independent wave-resolving controls')
    monkeypatch.setattr(runup,'RUNS',tmp_path)
    cases=runup.controls(Path(path))
    assert all(case['passed'] for case in cases.values())
    record=json.loads((tmp_path/'controls.json').read_text())
    assert record['gravity_similarity']['passed']
    monkeypatch.setattr(runup_controls,'RUNS',tmp_path)
    assert runup_controls.reflection_control(Path(path))['passed']


def test_swan_pulse_travel_time_and_file_backed_cache(tmp_path,monkeypatch):
    path=os.environ.get('TERLUNA_SWAN_EXECUTABLE')
    if not path:
        pytest.skip('Set TERLUNA_SWAN_EXECUTABLE for analytic group-delay controls')
    from climate.waves import history_controls
    from climate.waves.pilot import execute
    monkeypatch.setattr(history_controls,'RUNS',tmp_path)
    record=history_controls.run(Path(path))
    assert all(case['passed'] for case in record['cases'].values())
    files,_=history_controls.pulse_files()
    spectrum=tmp_path/'input.spc'
    spectrum.write_text(files['pulse.spc'])
    files['pulse.spc']=spectrum
    cached=execute(Path(path),tmp_path/'dt900',files,['pulse.tbl'],timeout_s=300,threads=1)
    assert cached==record['cases']['900']['run']
    spectrum.write_text(spectrum.read_text()+'\n')
    with pytest.raises(ValueError,match='identity differs'):
        execute(Path(path),tmp_path/'dt900',files,['pulse.tbl'],timeout_s=300,threads=1)
