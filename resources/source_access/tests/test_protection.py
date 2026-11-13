"""Protocol status, independent archive gates and consistent source exclusions."""
from copy import deepcopy
import json
import pytest
from resources.source_access import protection as p


@pytest.fixture(scope='module')
def config():
    return p.load_inputs()


def subject(region='near_earth'):
    return dict(source_id='synthetic-ordinary',name='Synthetic ordinary test source',region=region)


def all_attestations(config):
    return {g['id']:{'status':'passed','evidence':['synthetic-test-record:'+g['id']]} for g in config['gates']}


def test_unreviewed_source_is_unknown_not_failed_or_approved(config):
    result=p.evaluate_release(subject(),config=config)
    assert result['release_status']=='unresolved'
    assert result['unknown_conditions'] and not result['failed_conditions']
    assert not result['extraction_authorized']


def test_supplied_complete_records_do_not_authorize_extraction(config):
    result=p.evaluate_release(subject(),all_attestations(config),config)
    assert result['release_status']=='conditions_satisfied'
    assert result['conditions_satisfied'] and not result['extraction_authorized']


@pytest.mark.parametrize('status',['failed','hold'])
def test_explicit_negative_or_hold_is_distinct_from_unknown(config,status):
    records=all_attestations(config);records['pre_disturbance_survey']={'status':status,'evidence':['synthetic-negative-record']}
    result=p.evaluate_release(subject(),records,config)
    assert result['release_status']=='blocked'
    assert result[('failed_conditions' if status=='failed' else 'held_conditions')]==['pre_disturbance_survey']
    assert result['unknown_conditions']==[]


@pytest.mark.parametrize('gate',['native_sample_conditions','unopened_samples','multiple_archive_locations','institutional_continuity','recorded_provenance','open_scientific_review','open_record_access','archive_before_disturbance'])
def test_physical_clearance_does_not_substitute_for_each_data_condition(config,gate):
    records=all_attestations(config);records.pop(gate)
    result=p.evaluate_release(subject(),records,config)
    assert result['release_status']=='unresolved' and result['unknown_conditions']==[gate]


def test_named_hold_survives_all_synthetic_pass_attestations(config):
    s=dict(subject(),source_id='20101955',name='Renamed catalogue display')
    result=p.evaluate_release(s,all_attestations(config),config)
    assert result['release_status']=='blocked' and result['source_holds'][0]['name']=='Bennu'
    assert not result['extraction_authorized']


def test_excluded_body_cannot_return_in_broad_minus_narrow(config):
    bennu=dict(source_id='20101955',name='Bennu');ryugu=dict(spkid='20162173',full_name='162173 Ryugu (1999 JU3)')
    a=dict(source_id='a',name='Ordinary A');b=dict(source_id='b',name='Ordinary B')
    result=p.source_difference([a,b,bennu,ryugu],[a],config)
    assert result['sources']==[b]
    assert len(result['excluded']['broad'])==2
    filtered=p.filter_source_sets({'narrow':[a,bennu],'wide':[a,b,ryugu],'broad':[a,b,bennu,ryugu]},config)
    assert all(not p.exclusion_reasons(r,config) for rows in filtered['source_sets'].values() for r in rows)


def test_protected_categories_and_explicit_holds_are_not_filtered_by_fraction(config):
    s=dict(subject(),protected_category='ocean_world')
    assert p.evaluate_release(s,config=config)['release_status']=='blocked'
    s=dict(subject(),protected=True,protection_record='synthetic-site-designation')
    assert p.filter_source_sets({'one':[s]},config)['source_sets']['one']==[]


def test_unsupported_pass_and_unknown_gate_cannot_satisfy_protocol(config):
    with pytest.raises(ValueError,match='supporting record'):p.evaluate_release(subject(),{'pre_disturbance_survey':{'status':'passed'}},config)
    with pytest.raises(ValueError,match='Unknown protection'):p.evaluate_release(subject(),{'money_paid':{'status':'passed','evidence':['budget']}},config)


def test_lunar_cold_trap_gate_is_scoped_without_importing_open_moon_budgets(config):
    lunar=p.evaluate_release(subject('lunar'),config=config)
    nea=p.evaluate_release(subject(),config=config)
    assert 'polar_record_protection' in lunar['unknown_conditions']
    assert 'polar_record_protection' not in nea['unknown_conditions']
    assert 'ordinary_source_selection' in nea['unknown_conditions']
    assert all(g['id'] not in ('S1','S2','S3','S5','S6','S7') for g in lunar['conditions'])


def test_source_hash_change_fails_closed(config,tmp_path):
    changed=deepcopy(config);changed['sources']['conservation']['sha256']='0'*64
    path=tmp_path/'bad.json';path.write_text(json.dumps(changed))
    with pytest.raises(ValueError,match='Changed protection'):p.load_inputs(path)


def test_protected_id_alias_cannot_be_hidden_by_another_field(config):
    s=dict(subject(),spkid='20101955',name='Renamed catalogue display')
    assert p.exclusion_reasons(s,config)[0]['name']=='Bennu'
    assert p.source_difference([s],[],config)['sources']==[]


def test_conflicting_unprotected_ids_fail_set_subtraction(config):
    s=dict(subject(),spkid='different-catalogue-key')
    with pytest.raises(ValueError,match='Conflicting source IDs'):
        p.source_difference([s],[],config)
