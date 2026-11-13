import json
import math
import pytest
from research.studies.vacuum_buildout import run as m


@pytest.fixture(scope='module')
def inputs():
    return m.inputs()


def test_commissioning_dates_and_closed_cap(inputs):
    q,_,_=inputs
    assert m.population(8,q)==q['seed_population']
    assert m.population(10,q)==q['seed_population']
    assert m.population(200,q)==pytest.approx(q['cap_population'])
    assert m.population(1000,q)==m.population(200,q)
    starts=sum(m.population(t+2,q)-m.population(t+1,q) for t in range(8,199))
    assert starts+q['seed_population']==pytest.approx(q['cap_population'])
    assert all(m.population(t+1,q)>=m.population(t,q) for t in range(201))
    assert m.population(11,q)-m.population(10,q)<m.population(100,q)-m.population(99,q)
    assert m.population(200,q)-m.population(199,q)<m.population(160,q)-m.population(159,q)


def test_finite_inventory_and_area_conversion(inputs):
    q,e,inv=inputs;rows=m.trajectory(q,e,'services')
    last=rows[-1]
    assert last['cumulative_source_feed_kg']==pytest.approx(sum(r['source_feed_kg_year'] for r in rows))
    assert last['cumulative_lunar_raw_kg']==pytest.approx(sum(r['lunar_raw_kg_year'] for r in rows))
    assert last['excavated_area_km2']*1e6*q['regolith_density_kg_m3']*q['worked_depth_m']==pytest.approx(last['cumulative_lunar_raw_kg'])
    assert last['prospective_area_km2']==last['excavated_area_km2']/q['accessible_area_fraction']


def test_stress_starts_at_reference_and_retains_constraints(inputs):
    q,e,_=inputs
    for stress in ('lunar_stress','source_stress'):
        base=m.trajectory(q,e,'services');poor=m.trajectory(q,e,'services',stress)
        assert poor[0]['source_feed_kg_year']==base[0]['source_feed_kg_year']
        key='source_feed_kg_year' if stress=='source_stress' else 'lunar_raw_kg_year'
        assert poor[-1][key]>base[-1][key]
        assert poor[-1]['poorer_feed_output_share']==1
        assert poor[-1]['feasible']==(poor[-1]['max_capacity_use']<=1+1e-10)


def test_no_extrapolation_beyond_physical_grid(inputs):
    _,e,_=inputs
    with pytest.raises(ValueError):m.response(e,'services',.5)


def test_saved_product_and_all_portfolio_capacity(inputs):
    result=m.product(); saved=json.loads((m.HERE/'results/buildout.json').read_text())
    assert result==saved
    base=[c for c in result['cases'] if c['scenario']=='base']
    assert len(base)==9
    assert all(c['first_capacity_exceedance'] is None for c in base)
    assert all(c['annual'][-2]['population']==pytest.approx(result['schedule']['cap_population']) for c in base)
    for c in base:
        assert c['annual'][-1]['construction_people']==0
        assert c['annual'][-1]['source_feed_kg_year']>0
