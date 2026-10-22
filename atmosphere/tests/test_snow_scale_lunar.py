"""Check the matched forcing and mass statistics, independently of cloud skill."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from atmosphere.snow import scale_lunar as l


def test_profile_and_influx_reproduce_the_original_lunar_case():
    w,b,s,states=l.specification()
    assert s['depth_m']==pytest.approx(w.height_of(258.15)-w.z_freeze,abs=1e-9)
    assert s['total_time_s']==pytest.approx(3*s['depth_m']/s['injection_reference_velocity_m_s'])
    assert s['injected_number_flux_m2_s']==pytest.approx(1e5*s['injection_reference_velocity_m_s'])
    assert [v['t'] for v in states]==pytest.approx(273.15-(np.arange(30)+.5)*.5)
    for z,v in zip(np.arange(30)+.5,states):
        assert v==w.state(float(w.z_freeze+z*s['depth_m']/30))
    assert s['snowfall_mm_h']/3600*s['total_time_s']==pytest.approx(397.7678398340825)


def test_lunar_source_patch_rejects_unknown_upstream():
    with pytest.raises(ValueError,match='Unidentified'):l.patch('unknown')


def test_sticking_sensitivity_changes_only_efficiency_after_passive_observer():
    source=l.CACHE/'work'/l.observer.SOURCE_REL
    if not source.exists():pytest.skip('Optional upstream source cache')
    original=source.read_text()
    native=l.collision_patch(original)
    changed=l.collision_patch(original,'inhouse')
    old='             c_rate(tc) = 0.1_RP'
    new='             k=floor(sd_rk(icptc))+1; i=floor(sd_ri(icptc))+1; j=floor(sd_rj(icptc))+1\n             c_rate(tc) = exp(0.09_RP*(t_scale(k,i,j)-t0))'
    assert changed.replace(new,old)==native


def test_arrival_analysis_uses_mass_weights_and_separate_area_and_time(tmp_path):
    # Two time blocks, two different columns, unequal particle weights. Native
    # diameter and the same mass under the old law are deliberately different.
    settings=dict(duration_s=120.,gather_start_s=60.,horizontal_cells=8,dx_m=1.,layers=1,depth_m=1.,mode='prescribed',profile=[dict(e=0.,t=270.)])
    l.boxes.save(tmp_path/'settings.json',settings)
    l.boxes.save(tmp_path/'execution.json',[dict(exit_code=0),dict(exit_code=0)])
    ledger=np.zeros((1,len(l.LEDGER_FIELDS)));ledger[0,0]=120.;ledger[0,1]=1.
    np.savetxt(tmp_path/'lunar_ledger.txt',ledger)
    rows=np.array([[30.,0,10,1e-6,.001,500,1,1,10],
                   [65.,0,2,1e-6,.002,500,1,2,10],
                   [110.,1,1,6e-6,.003,500,1,3,10]],dtype='>f8')
    rows.tofile(tmp_path/'lunar_arrivals.bin')
    r=l.analyze(tmp_path)
    assert r['arrivals']['packets']==2
    assert r['arrivals']['mass_kg_m2']==pytest.approx(1e-6)
    assert r['arrivals']['dimension_mass_quantiles_m'][1]==.003
    assert r['arrivals']['snowfall_mm_h']==pytest.approx(6e-5)
    assert r['horizontal_columns'][0]['snowfall_mm_h']==pytest.approx(12e-5)
    assert sum(x['mass_kg_m2'] for x in r['gathering_time_blocks'])==pytest.approx(1e-6)
    assert sum(x['precipitation_mm_h'] for x in r['gathering_time_blocks'])/4==pytest.approx(6e-5)
    assert r['arrivals']['inherited_geometry_dimension_mass_quantiles_m'][1]>.003
    rows[:,8]=1
    rows.tofile(tmp_path/'lunar_arrivals.bin')
    rain=l.analyze(tmp_path)['arrivals']
    assert rain['snowfall_mm_h']==0
    assert rain['precipitation_mm_h']==pytest.approx(6e-5)
    assert rain['ice_dimension_mass_quantiles_m']==[None]*4


@pytest.mark.parametrize('kw',[{'mode':'unknown'},{'packets':0},{'dt':0.},{'seed':-1},{'dx':0}])
def test_invalid_forcing_does_not_create_directories(tmp_path,kw):
    with pytest.raises(ValueError):l.configure('invalid',cache=tmp_path,**kw)
    assert not list(tmp_path.iterdir())
