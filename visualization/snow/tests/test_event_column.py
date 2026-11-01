"""A scientific figure must reject changed accounting, distributions and gates."""
from copy import deepcopy
import json
import pytest
from visualization.snow import event_column as view


@pytest.fixture(scope='module')
def product():
    path=view.ROOT/view.PRODUCT
    if not path.exists():pytest.skip('Full-column product has not been produced')
    return json.loads(path.read_text())


def test_complete_product_matches_native_quantities(product):
    verified=view.validate(product)
    assert 'atmosphere/snow/event_column.cpp' in verified


@pytest.mark.parametrize('change',['water','quantile','cdf','missing','acceptance','producer'])
def test_reject_changed_scientific_artifact(product,change):
    value=deepcopy(product)
    if change=='water':value['cases'][0]['result']['mass_ledger_kg_m2']['deposited']*=1.01
    if change=='quantile':value['cases'][0]['result']['summary']['mass_d50_m']*=10
    if change=='cdf':value['groups']['moon_v4']['mean_mass_cdf'][10]+=.01
    if change=='missing':value['cases'].pop()
    if change=='acceptance':value['acceptance']['passed']=not value['acceptance']['passed']
    if change=='producer':value['cases'][0]['producer']['files']['atmosphere/snow/event_column.cpp']='0'*64
    with pytest.raises(ValueError):view.validate(value)
