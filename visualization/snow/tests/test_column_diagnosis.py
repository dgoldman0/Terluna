"""The diagnostic figure must preserve the recorded water and size comparisons."""
import copy
import json
import pytest
from visualization.snow import column_diagnosis as view

@pytest.fixture(scope='module')
def product():return json.loads((view.ROOT/view.PRODUCT).read_text())


def test_complete_product_is_admitted(product):
    verified=view.validate(product)
    assert 'atmosphere/snow/column_diagnosis.cpp' in verified


@pytest.mark.parametrize('field',['water','empty','quantile','cdf','decision','producer','incomplete'])
def test_changed_scientific_fields_are_rejected(product,field):
    p=copy.deepcopy(product)
    if field=='water':p['groups']['full_v1']['mean_process_summary']['gather_empty_unused_mm_h']+=.1
    if field=='empty':p['groups']['full_v1']['mean_process_profiles']['empty_fraction'][0]+=.1
    if field=='quantile':p['groups']['full_v1']['mean_summary']['mass_d99_m']*=2
    if field=='cdf':p['groups']['full_v1']['mean_mass_cdf'][50]+=.1
    if field=='decision':p['comparisons'][0]['passed']=not p['comparisons'][0]['passed']
    if field=='producer':p['producer']['files']['atmosphere/snow/column_diagnosis.cpp']='0'*64
    if field=='incomplete':p['completed_cases']-=1
    with pytest.raises(ValueError):view.validate(p)


def test_sticking_only_followup_is_admitted_and_really_removes_both_processes():
    p=json.loads((view.ROOT/view.FOLLOWUP).read_text())
    assert view.validate(p,expected_cases=8)
    for row in p['cases']:
        r=row['result']
        assert r['settings']['growth'] is False and r['settings']['breakup'] is False
        assert r['settings']['collisions'] is True
        assert r['mass_ledger_kg_m2']['deposited']==0
        assert r['number_ledger_m2']['chipping_production']==0
        assert r['number_ledger_m2']['coagulation_loss']>0
