import copy,gzip,json,math
from pathlib import Path
import pytest
from resources.source_access import neo,costs,run
HERE=Path(__file__).resolve().parents[1]

def config():return json.loads((HERE/'access_inputs.json').read_text())

def row(d=1,sigma=.1):return {'spkid':'999','full_name':'Test','diameter':str(d),'diameter_sigma':sigma,'spec_B':None,'spec_T':None}

def test_cubic_volume_and_density():
    assert neo.mass(row(2),1000)==pytest.approx(8*neo.mass(row(1),1000))
    assert neo.mass(row(),2000)==pytest.approx(2*neo.mass(row(),1000))

def test_formal_error_excludes_missing_and_does_not_shrink_common_scale():
    c=config();rows=[row(1,.1),{**row(2,None),'spkid':'998','full_name':'Other'}]
    s=neo.summarize(rows,c,{'count':0,'data':[]})
    assert s['reported_diameter_error']['missing_sigma_count']==1
    assert s['reported_diameter_error']['linearized_independent_sigma_kg']==pytest.approx(.3*neo.mass(rows[0],1500))
    assert s['structural_envelopes']['0.2']['mass_kg'][0]==pytest.approx(s['mass_reference_kg']*.8**3*1000/1500)

def test_bad_catalogue_and_ambiguous_join_fail():
    with pytest.raises(ValueError):neo.rows({'fields':['spkid'],'count':2,'data':[['1']]})
    with pytest.raises(ValueError):neo.summarize([row()],config(),{'count':2,'data':[{'fullname':'A'},{'fullname':' A '}]})

def test_source_identity_mismatch(tmp_path,monkeypatch):
    c={'archive_directory':'.','external_files':[{'path':'x','bytes':2,'sha256':'wrong'}]};(tmp_path/'x').write_text('{}');monkeypatch.setattr(neo,'ROOT',tmp_path)
    with pytest.raises(ValueError,match='identity'):neo.load_archive(c)

def test_cost_scaling_and_no_imaginary_yield():
    c=config();coeff=json.loads((HERE/c['cost_coefficients']['path']).read_text());case=c['costs']['cases']['reference']
    a=costs.unit_feed_cost(coeff,case,.03);b=costs.unit_feed_cost(coeff,{**case,'water_grade':case['water_grade']/2},.03)
    assert b['source_water_USD2041_kg']==pytest.approx(2*a['source_water_USD2041_kg'])
    assert b['source_USD2041_kg_feed']==a['source_USD2041_kg_feed']
    assert sum(a['source_feed_components'].values())==pytest.approx(a['source_USD2041_kg_feed'])
    assert sum(a['lunar_bulk_components'].values())==pytest.approx(a['lunar_bulk_USD2041_kg_product'])
    with pytest.raises(ValueError):costs.unit_feed_cost(coeff,{**case,'water_recovery':0},.03)
    assert costs.crf(0,10)==.1

def test_demand_includes_seed_once():
    groups={k:{'mass_reference_kg':100,'structural_envelopes':{'0.2':{'mass_kg':[50,200]}}} for k in ['narrow','wider']}
    case={'physical_scenario':'base','seed':{'source_feed_kg':3,'lunar_raw_kg':5},'annual':[{'year':100,'cumulative_source_feed_kg':7,'cumulative_lunar_raw_kg':11}]}
    result=run.demand_summary({'cases':[case]*81},[100],groups)[0]
    assert result['cumulative_source_feed_kg_range']==[10,10]
    assert result['cumulative_lunar_feed_kg_range']==[16,16]

def test_committed_model_does_not_qualify_unknown_supply():
    x=json.loads(gzip.decompress((HERE/'results/access.json.gz').read_bytes()))
    assert x['qualified_supply_kg'] is None
    assert all(p['release_status']=='unresolved' and not p['extraction_authorized'] for p in x['protection_examples'].values())
    assert x['neo']['groups']['wider']['count']==131
    assert x['neo']['groups']['narrow']['count']==17
    excluded={'20101955','20162173','20025143','20000433'}
    for g in x['neo']['groups'].values():assert not excluded & {r['spkid'] for r in g['objects']}
    assert x['costs']['lunar_oxide_feed_burden']['TiO2']['zero_grade_surface_fraction']>0
    assert x['costs']['lunar_oxide_feed_burden']['TiO2']['formal_1_96sigma_lower_nonpositive_surface_fraction']>.3


def test_integrated_product_reproduction_and_identity():
    archive=neo.ROOT/config()['archive_directory']
    if not (archive/'nea_catalogue.json').exists():pytest.skip('Exact retained source archive is not installed')
    current=run.product()
    saved=json.loads(gzip.decompress((HERE/'results/access.json.gz').read_bytes()))
    assert current==saved
    joined=current['integrated_candidate_evidence']['objects']
    assert len(joined)==len({r['source_id'] for r in joined})==131
    assert all(isinstance(r['source_id'],str) and r['protection_status']=='unresolved' for r in joined)
    assert sum(r['spectral_family']!='no_observation' for r in joined)==17
    assert sum(r['has_period_record'] for r in joined)==34


def test_disjoint_composition_allocation_conserves_count_and_mass():
    x=json.loads(gzip.decompress((HERE/'results/access.json.gz').read_bytes()))
    stratified=x['composition_stratified_estimate']
    for measure,estimate in stratified['conditional_estimates'].items():
        assert sum(estimate['family_measure'].values())+estimate['unallocated_measure']==pytest.approx(estimate['total_known_measure'])
        assert estimate['total_known_measure']==pytest.approx(x['composition_pools']['wider_geometry_all_albedos']['integrated_classifications'][measure]['total_known_measure'])
    assert stratified['strata']['dark']==x['composition_pools']['wider_dark_screen']['integrated_classifications']
