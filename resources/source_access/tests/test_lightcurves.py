"""Pinned-source identity, qualified joins and shared source protection."""
from copy import deepcopy
import gzip
import hashlib
import io
import json
from pathlib import Path
import zipfile
import pytest
from resources.source_access import lightcurves as lc, neo


def row(number='123',name='Synthetic',designation='2000 AA',period='12.5',flag='',quality='3',notes='',binary=''):
    chars=[' ']*217
    for first,last,value in [(1,7,str(number).rjust(7)),(11,40,name),(42,61,designation),
                             (144,144,flag),(146,158,period),(188,189,quality),
                             (191,195,notes),(197,199,binary)]:
        chars[first-1:last]=list(value.ljust(last-first+1))
    return ''.join(chars)


def summary(*rows):
    return ('ASTEROID LIGHTCURVE DATABASE (LCDB) SUMMARY (Complete)\nGENERATED: synthetic\n\ncolumns\n---\n'+'\n'.join(rows)+'\n').encode()


def index(*rows):
    return lc.parse_summary(summary(*rows),['2-','2','2+','3-','3'])


def test_limits_quality_and_shape_notes_survive_without_resource_inference():
    result=index(row(flag='>',quality='',notes='T?',binary='?'))
    data,method=lc.match({'full_name':'123 Synthetic (2000 AA)'},result)
    assert method=='MPC_number' and data['period_hours']==12.5 and data['period_flag']=='>'
    assert data['period_quality']=='' and not data['period_quality_meets_LCDB_statistical_threshold']
    assert data['notes']=='T?' and data['binary_flag']=='?'
    assert all(k not in data for k in ('albedo','diameter','taxonomy','water_grade','contact_binary'))


def test_number_zero_joins_separate_unnumbered_designations():
    result=index(row('0','2000 AA','2000 AA'),row('0','2000 BB','2000 BB',period='20'))
    a,method=lc.match({'full_name':'(2000 AA)'},result)
    b,_=lc.match({'full_name':'(2000 BB)'},result)
    assert method=='exact_designation' and a['period_hours']==12.5 and b['period_hours']==20
    assert result['by_number']=={}


def test_missing_match_is_unknown_and_approximate_name_is_not_joined():
    result=index(row())
    assert lc.match({'full_name':'999 Other (2000 AB)'},result)==(None,None)


def test_duplicate_number_and_ambiguous_designation_are_rejected():
    with pytest.raises(ValueError,match='Duplicate LCDB'):index(row(),row())
    result=index(row('0','Object A','2000 AA'),row('0','Object B','2000 AA'))
    with pytest.raises(ValueError,match='Ambiguous exact'):lc.match({'full_name':'(2000 AA)'},result)


def test_nonfinite_physical_value_fails_closed():
    with pytest.raises(ValueError,match='Invalid LCDB physical'):index(row(period='nan'))


def test_shared_pools_protect_all_source_sets_and_do_not_authorize():
    config=json.loads((lc.HERE/'access_inputs.json').read_text())
    common={'diameter':'1','albedo':'.05','a':'1','e':'.1','i':'5'}
    records=[dict(common,spkid='20101955',full_name='101955 Bennu (1999 RQ36)'),
             dict(common,spkid='123',full_name='123 Synthetic (2000 AA)')]
    screened=neo.pools(records,config)
    assert all(len(v)==1 for v in screened['source_sets'].values())
    candidates=lc.join(screened['source_sets'],index(row()),config['density_reference_kg_m3'])
    assert len(candidates)==1 and candidates[0]['source_id']=='123'
    assert set(candidates[0]['groups'])==set(screened['source_sets'])
    assert candidates[0]['extraction_authorized'] is False


def test_zip_and_member_identities_are_independently_checked(tmp_path):
    data=b'payload';buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,'w') as z:z.writestr('member.txt',data)
    raw=buffer.getvalue();(tmp_path/'archive.zip').write_bytes(raw)
    config={'archive':{'filename':'archive.zip','bytes':len(raw),'sha256':lc.sha(raw)},
            'members':[{'name':'member.txt','bytes':len(data),'sha256':lc.sha(data)}]}
    assert lc.read_members(config,tmp_path)=={'member.txt':data}
    bad=deepcopy(config);bad['members'][0]['sha256']='0'*64
    with pytest.raises(ValueError,match='identity mismatch'):lc.read_members(bad,tmp_path)
    bad=deepcopy(config);bad['archive']['sha256']='0'*64
    with pytest.raises(ValueError,match='identity mismatch'):lc.read_members(bad,tmp_path)


def test_duplicate_zip_member_is_not_silently_selected(tmp_path):
    buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,'w') as z:
        z.writestr('member.txt',b'a')
        with pytest.warns(UserWarning):z.writestr('member.txt',b'b')
    raw=buffer.getvalue();(tmp_path/'archive.zip').write_bytes(raw)
    config={'archive':{'filename':'archive.zip','bytes':len(raw),'sha256':lc.sha(raw)},
            'members':[{'name':'member.txt','bytes':1,'sha256':lc.sha(b'a')}]}
    with pytest.raises(ValueError,match='Missing or duplicate'):lc.read_members(config,tmp_path)


def test_deterministic_compressed_output(tmp_path):
    a=tmp_path/'a.json.gz';b=tmp_path/'b.json.gz'
    result={'b':2,'a':1}
    digest=lc.write_product(a,result);lc.write_product(b,result)
    assert a.read_bytes()==b.read_bytes()
    assert digest==hashlib.sha256(gzip.decompress(a.read_bytes())).hexdigest()


def test_pinned_observed_product_and_counts():
    config=lc.load_inputs();cache=lc.ROOT/config['archive_directory']/config['archive']['filename']
    if not cache.exists():pytest.skip('Pinned LCDB archive not restored in ignored input cache')
    result=lc.product()
    assert result['groups']['narrow']['candidates']==17
    assert result['groups']['narrow']['lcdb_matches']==2
    assert result['groups']['wider']['candidates']==131
    assert result['groups']['wider']['lcdb_matches']==34
    assert result['groups']['dark_catalogue']['candidates']==496
    assert result['groups']['dark_catalogue']['lcdb_matches']==109
    byid={r['source_id']:r for r in result['candidates']}
    assert byid['20256412']['lcdb']['period_flag']=='>'
    assert byid['20256412']['lcdb']['period_quality']==''
    assert byid['20065679']['lcdb']['period_hours']==7.746
    assert byid['20464798']['lcdb']['notes']=='T'
    assert not {'20101955','20162173','20025143','20000433'} & byid.keys()
    frozen=json.loads(gzip.decompress((lc.HERE/'results/lightcurves.json.gz').read_bytes()))
    assert frozen==result
