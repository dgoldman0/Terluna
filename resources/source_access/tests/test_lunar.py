"""Controls for archive identity, spherical weighting and formal-error units."""
import gzip
import hashlib
import json
import math
import xml.etree.ElementTree as ET

import pytest
from resources.source_access import lunar


def synthetic_table(rows):
    """Independent fixed-width PDS fixture, with the real coordinate order."""
    names=list(lunar.FIELD_NAMES)+[f'E[{i},{j}]' for i in range(9) for j in range(i,9)]
    ns='http://pds.nasa.gov/pds4/pds/v1'
    root=ET.Element(f'{{{ns}}}Product_Observational')
    table=ET.SubElement(root,f'{{{ns}}}Table_Character')
    ET.SubElement(table,f'{{{ns}}}records').text=str(len(rows))
    record=ET.SubElement(table,f'{{{ns}}}Record_Character')
    ET.SubElement(record,f'{{{ns}}}record_length').text=str(16*61+2)
    for i,name in enumerate(names):
        field=ET.SubElement(record,f'{{{ns}}}Field_Character')
        for tag,value in [('name',name),('field_location',16*i+1),('field_length',16)]:
            ET.SubElement(field,f'{{{ns}}}{tag}').text=str(value)
    return b''.join((''.join(f'{v:16.8e}' for v in row)+'\r\n').encode() for row in rows),ET.tostring(root)


def row(index=0,bounds=(-90,90,-180,180)):
    values=[index,*bounds,22,1e-6,.08,.25,.45,.15,.005,.05,500,2,.54]
    C=[[0.0]*9 for _ in range(9)]
    C[1][1]=.09;C[5][5]=.04;C[1][5]=C[5][1]=.01
    C[4][4]=.0001;C[6][6]=100;C[1][6]=C[6][1]=-.4
    return values+[C[i][j] for i in range(9) for j in range(i,9)]


def test_fixed_width_units_full_covariance_and_nonclosure():
    p=lunar.parse_pds_table(*synthetic_table([row()]))[0]
    assert p['abundances']['FeO']==.05
    assert p['abundances']['K']==500 # ppm, not a mass fraction
    assert p['six_oxide_sum']==pytest.approx(.985) # original nonclosure retained
    assert p['covariance'][1][5]==p['covariance'][5][1]==.01
    assert p['covariance'][1][6]==p['covariance'][6][1]==-.4
    assert len(p['covariance'])==9 and all(len(x)==9 for x in p['covariance'])
    assert p['solid_angle_sr']==pytest.approx(4*math.pi)


def test_linear_propagation_includes_cross_terms_and_explicit_units():
    p=lunar.parse_pds_table(*synthetic_table([row()]))[0]
    result=lunar.linear_estimate(p,{'FeO':2,'Al2O3':-3},2)
    assert result['mean']==pytest.approx(-.65)
    assert result['formal_sigma']==pytest.approx(math.sqrt(.85))
    assert result['formal_interval']==pytest.approx([-.65-2*math.sqrt(.85),-.65+2*math.sqrt(.85)])
    # Converting ppm into a mass fraction must scale the variance twice.
    assert lunar.linear_estimate(p,{'K':1e-6})['formal_sigma']==pytest.approx(1e-5)
    p['abundances']['TiO2']=0
    assert lunar.linear_estimate(p,{'TiO2':1})['formal_interval']==pytest.approx([-.01,.01])


@pytest.mark.parametrize('mutation',[lambda r:r.__setitem__(11,-.1),lambda r:r.__setitem__(12,float('nan')),lambda r:r.__setitem__(0,1)])
def test_bad_assays_or_indices_are_rejected(mutation):
    r=row();mutation(r)
    with pytest.raises(ValueError):lunar.parse_pds_table(*synthetic_table([r]))


def test_truncated_or_mislabelled_archive_rejected():
    raw,label=synthetic_table([row()])
    with pytest.raises(ValueError,match='record size'):lunar.parse_pds_table(raw[:-1],label)
    with pytest.raises(ValueError,match='field order'):lunar.parse_pds_table(raw,label.replace(b'W_FeO',b'wrong'))


def test_area_weighting_and_latitude_partition_do_not_create_samples():
    p=lunar.parse_pds_table(*synthetic_table([row()]))
    pieces=[lunar.latitude_band(p,s,s+30) for s in range(-90,90,30)]
    assert sum(lunar.summarize(x)['solid_angle_sr'] for x in pieces)==pytest.approx(4*math.pi)
    south=lunar.summarize(pieces[0])
    assert south['fraction_of_lunar_surface']==pytest.approx((1-math.sqrt(3)/2)/2)
    assert pieces[0][0]['covariance'] is p[0]['covariance']
    assert south['aggregate_assay_uncertainty'] is None
    # Latitude clipping gives new weights, not a fictitious higher-resolution assay.
    assert all(x[0]['abundances']==p[0]['abundances'] for x in pieces)


def test_weighted_quantile_convention_and_zero_weight():
    assert lunar.weighted_quantiles([10,0],[1,1],[0,.5,.75,1])==[0,0,5,10]
    assert lunar.weighted_quantiles([100,10,0],[0,1,1],[.75])==[5]
    for args in [([],[]),([1],[0]),([1],[-1]),([float('nan')],[1])]:
        with pytest.raises(ValueError):lunar.weighted_quantiles(*args)


def test_identity_check_fails_before_parsing_or_substitution(tmp_path):
    raw=b'original'
    config={'files':{'table':{'filename':'table','url':'https://example.invalid/not-used','bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}}}
    (tmp_path/'table').write_bytes(b'changed!')
    with pytest.raises(ValueError,match='hash or length mismatch'):
        lunar.load_inputs(config,archive=tmp_path)


def test_actual_archive_and_saved_product_reproduce():
    config=json.loads((lunar.HERE/'lunar_inputs.json').read_text())
    folder=lunar.ROOT/config['archive_directory']
    if not all((folder/x['filename']).exists() for x in config['files'].values()):
        pytest.skip('Exact external PDS inputs absent; restore with python -m resources.source_access.lunar --download')
    result=lunar.product()
    saved=json.loads(gzip.decompress((lunar.HERE/'results/lunar.json.gz').read_bytes()))
    assert result==saved
    assert result['global_summary']['solid_angle_sr']==pytest.approx(4*math.pi,abs=1e-12)
    assert result['validity']['invalid_or_negative_pixels']==0
    assert result['depth_extrapolation'] is None and result['qualified_reserve'] is None
    p=result['pixels'][0]
    assert p['abundances']['FeO']==.049428
    assert p['abundances']['TiO2']==.0016727
    assert p['abundances']['K']==433.89
    assert p['covariance'][0][0]==6.6536e-6
    assert p['covariance'][0][3]==-1.8852e-6
    assert p['covariance'][0][6]==-.02231
    assert p['covariance'][6][0]==-.02231
    for i in range(9):
        for j in range(9):
            assert p['covariance'][i][j]==p['covariance'][j][i]
    assert result['global_summary']['six_oxide_sum']['minimum']==pytest.approx(.886852)
    assert result['global_summary']['six_oxide_sum']['maximum']==pytest.approx(1.0708517)
    assert math.fsum(b['summary']['solid_angle_sr'] for b in result['latitude_bands'])==pytest.approx(4*math.pi)
    for oxide in lunar.OXIDES:
        weighted=math.fsum(b['summary']['composition'][oxide]['area_weighted_mean']*b['summary']['solid_angle_sr'] for b in result['latitude_bands'])/(4*math.pi)
        assert weighted==pytest.approx(result['global_summary']['composition'][oxide]['area_weighted_mean'])
