"""Controls for author-build integrity, settings and population-profile handling."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from atmosphere.snow import welss_aggregation as w
from atmosphere.snow import welss_aggregation_build as b


def test_seed_patch_changes_only_seed_acquisition():
    original='  INTEGER :: values(1:8), k\n  seed(:) = values(8)\n'
    patched=b.deterministic_main(original)
    assert patched.replace('  CHARACTER(len=32) :: deterministic_seed\n','').replace('  CALL get_command_argument(1, deterministic_seed)\n','').replace('  READ(deterministic_seed, *) values(8)\n','')==original
    with pytest.raises(ValueError):b.deterministic_main('unexpected source')


def test_manifest_blob_integrity_rejects_corruption():
    data=b'example\n';item=dict(path='fixture',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),git_blob_sha1=hashlib.sha1(b'blob 8\0'+data).hexdigest())
    b.verify(data,item)
    with pytest.raises(ValueError):b.verify(data+b'!',item)


def test_published_case_settings_and_source_normalization():
    growth=w.configuration(False);agg=w.configuration(True)
    assert growth['sp_nml']['lwc0']==.0002
    assert agg['sp_nml']['lwc0']==0
    assert growth['sp_nml']['sp_kernel_id']==0
    assert agg['sp_nml']['sp_kernel_id']==2
    assert agg['sp_nml']['sp_stick_effi']==1
    assert agg['general_nml']['rime_type']==0
    assert agg['sp_nml']['sp_kern_sig']==0
    assert agg['sp_nml']['sp_frag']==0
    # This checks only the continuous parameter product; floor admission is below.
    for xi in (825,1650,3300):
        c=w.configuration(True,xi=xi)
        assert c['sp_nml']['nclrate']*xi==pytest.approx(16.5)
    text=w.namelist(agg)
    assert 'interpol_atmos = .TRUE.' in text
    assert "atmofile = 'ecmwf_profile.txt'" in text


@pytest.mark.parametrize('xi,dt,births,physical_rate',[
    (1650,5.,75,16.5),
    (3300,5.,37,16.28),
    (1650,2.5,37,16.28),
    (825,2.5,75,16.5),
])
def test_actual_source_zone_floor_and_refinement_mismatch(xi,dt,births,physical_rate):
    # mo_nucleation:45 floors the entire source zone, then samples uniform z.
    c=w.configuration(True,xi=xi,dtc=dt)
    area=c['grid_nml']['box_area']
    depth=c['atmo_nml']['nh2']-c['atmo_nml']['nh1']
    dz=(c['grid_nml']['dom_top']-c['grid_nml']['dom_bottom'])/c['grid_nml']['nz']
    cell_count=int(depth/dz)
    volume=cell_count*dz*area
    assert cell_count==15
    assert volume==1500.
    actual=int(np.floor(c['sp_nml']['nclrate']*dt*volume))
    assert actual==births
    assert actual*xi/(dt*volume)==pytest.approx(physical_rate)
    if physical_rate!=16.5:
        assert actual*xi/(dt*volume)<c['sp_nml']['nclrate']*xi


def test_ascii_parser_keeps_cell_moments_and_time_units(tmp_path):
    a=np.zeros((251,32));a[:,0]=np.arange(251)*20;a[:,2]=2;a[:,3]=3;a[:,4]=5;a[:,5]=7
    a[:,7]=11;a[:,8]=13;a[:,9]=17;a[:,10]=19
    p=tmp_path/'profiles.dat'
    with p.open('w') as f:
        f.write('#250\n\n#"300.0 min"\n');np.savetxt(f,a)
        f.write('#"600.0 min"\n');np.savetxt(f,2*a*np.r_[.5,np.ones(31)])
    t,frames=w.parse_profiles(p)
    np.testing.assert_array_equal(t,[18000,36000])
    avg=w.means(frames)
    for key,value in [('Nd',3),('Md',4.5),('Fn',7.5),('Fm',10.5),('Nd_mm1',16.5),('Fm_mm1',28.5)]:
        np.testing.assert_array_equal(avg[key],np.full(251,value))


def test_zero_to_missing_is_separate_from_physical_mean():
    a=np.zeros((2,251,32));a[1,:,2]=4
    np.testing.assert_array_equal(w.means(a)['Nd'],np.full(251,2))
    np.testing.assert_array_equal(w.means(a,positive_only=True)['Nd'],np.full(251,4))
    assert np.isnan(w.means(a,positive_only=True)['Md']).all()


def test_metrics_keep_source_region_discrepancy():
    z=np.array([0,20,2700,2980,3000]);reference=np.array([10,10,2,2,np.nan]);model=np.array([10,10,4,4,0])
    r=w.metrics(model,reference,z)
    assert r['below_nucleation']['relative_L1']==0
    assert r['nucleation_belt']['relative_L1']==1
    assert r['all_populated']['relative_L1']==pytest.approx(1/6)


@pytest.mark.parametrize('name',['welss_author_growth.json','welss_author_aggregation.json'])
def test_saved_author_result_identity_if_available(name):
    path=w.HERE/'results'/name
    if not path.exists():pytest.skip('author run not yet recorded')
    result=json.loads(path.read_text())
    assert result['producer']==w.identity()
    assert result['run']['producer']==w.identity()
    assert result['run']['build']['producer']=={k:w.identity()[k] for k in ['welss_aggregation_build.py','welss_aggregation_sources.json']}
    assert result['profiles']['last5h_31']['snapshot_count']==31
    assert result['profiles']['final4']['snapshot_count']==4
    checks=result['run']['pathway_checks']
    for k in ['maximum_saved_riming_rate_kg_s','maximum_saved_graupel_liquid_rimed_category_moment','final_maximum_rime_mass_kg','final_maximum_total_minus_pristine_minus_rime_mass_kg']:
        assert checks[k]==0
    for mode in result['profiles'].values():
        for values in mode['raw_profiles'].values():assert np.isfinite(values).all()
