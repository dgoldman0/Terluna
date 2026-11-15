"""Stock sufficiency is distinct from access and from product recovery."""
import numpy as np
import pytest
from resources.source_access import potential as m


def demand():
    return {'name':'case','baseline_case':'base','variant':'base','horizon':220,'source_feed_kg':100.,'source_water_equivalent_kg':10.,'source_virgin_metal_kg':2.,'physical_capacity_qualified':True,'mature_annual_mean':{'source_feed_kg':2.,'source_water_equivalent_kg':.5}}

def test_stock_frequency_and_finite_monte_carlo_resolution():
    x=m.stock_test([1.,2.,3.,4.],2.5)
    assert x['probability_sufficient_under_model']==.5
    assert x['monte_carlo_frequency_interval_95percent'][0]<.5<x['monte_carlo_frequency_interval_95percent'][1]
    x=m.stock_test(np.ones(2048)*100,1)
    assert x['probability_sufficient_under_model']==1
    assert x['monte_carlo_frequency_interval_95percent'][0]<1
    for stock,d in [([-1,2],1),([1,2],-1),([1,float('nan')],1)]:
        with pytest.raises(ValueError):m.stock_test(stock,d)


def test_body_feed_and_product_water_use_different_units():
    stock={'body_mass_kg':[1000.,2000.],'water_equivalent_kg':[100.,200.],'native_iron_kg':[20.,40.]}
    x=m.comparison(stock,demand(),.1,.5,.5)
    assert x['feed_mass']['probability_sufficient_under_model']==1
    assert x['water_product']['probability_sufficient_under_model']==.5
    assert x['joint_feed_and_water_probability']==.5
    assert x['native_iron_mass_proxy']['probability_sufficient_under_model']==.5
    # Neither the water recovery nor grade is applied again to the raw feed.
    y=m.comparison(stock,demand(),.1,1.,.5)
    assert y['feed_mass']==x['feed_mass']
    assert y['water_product']['probability_sufficient_under_model']==1
    assert x['remaining_feed_and_water_stock_years_p05_p50_p95']==[0.,0.,0.]
    with pytest.raises(ValueError):m.comparison({**stock,'water_equivalent_kg':[1001,2001]},demand(),.1,.5,.5)


def test_required_share_and_growth_free_sustaining_years():
    stock={'body_mass_kg':[1000.,1000.],'water_equivalent_kg':[100.,100.],'native_iron_kg':[20.,20.]}
    x=m.comparison(stock,demand(),.5,.5,.5)
    # 400 kg feed remains at2/year, but15kg water at0.5/year binds at30years.
    assert x['remaining_feed_and_water_stock_years_p05_p50_p95']==[30.,30.,30.]
    assert m.minimum_share([100,100],10,.5)['share_for_95percent_model_stock_exceedance']==.2
    assert m.minimum_share([0,0],10,.5)['share_at_median_stock'] is None
    with pytest.raises(ValueError):m.comparison(stock,demand(),1.1,.5,.5)


def test_dated_seed_is_counted_once_and_undated_rows_rejected():
    rows=[{'year':y,'source_feed_kg':10.,'lunar_raw_kg':5.,'source_water_equivalent_kg':1.,'source_nitrogen_kg':.01,'source_virgin_metal_kg':.2} for y in range(221)]
    c={'name':'case','baseline_case':'base','variant':'base','annual':rows,'physical_case':{'first_capacity_shortfall':None}}
    x=m.cumulative_audit(c,200)
    assert x['source_feed_kg']==2010 and x['source_water_equivalent_kg']==201
    assert x['mature_annual_mean']['source_feed_kg']==10
    with pytest.raises(ValueError):m.cumulative_audit({**c,'annual':rows[1:]},200)


def test_representative_pool_grade_is_not_selected_ore_grade():
    import json
    coeff=json.loads((m.HERE/'inputs/cost_coefficients.json').read_text())
    case={'equipment_cost_multiplier':1,'availability_multiplier':1,'electricity_USD2041_kWh':.15,'water_grade':.1,'water_recovery':.7,'ore_energy_multiplier':1}
    stock={'body_mass_kg':[100,100],'water_equivalent_kg':[2,2]}
    result=m.representative_grade_cost(stock,coeff,{'reference':case},.03,.7)
    assert result['feed_multiplier_relative_to_10percent_grade_at_same_recovery_p05_p50_p95']==[5,5,5]
    base=m.costs.unit_feed_cost(coeff,case,.03)['source_water_USD2041_kg']
    assert np.allclose(result['partial_source_water_USD2041_kg_by_cost_scenario']['reference'],base*5)


def test_saved_comparison_source_identities_and_exact_demand():
    import gzip,json
    c,products=m.load()
    result=json.loads(gzip.decompress((m.HERE/'results/potential.json.gz').read_bytes()))
    for path,digest in result['dependencies'].items():assert m.sha(m.HERE/path)==digest
    assert result['producer']['sha256']==m.sha(m.__file__)
    assert len(products['reclamation']['cases'])==36
    assert len(result['demand_checkpoints'])==108
    lookup={(r['name'],r['horizon']):r for r in result['demand_checkpoints']}
    for case in products['reclamation']['cases']:
        for h in c['checkpoints']:assert m.cumulative_audit(case,h)==lookup[case['name'],h]
    assert all(row['case_count']==81 for row in result['baseline_demand_ranges'])
    assert result['scientific_protection_release'] is None
    assert result['extraction_authorized'] is False


def test_saved_access_monotonicity_and_conditional_draw_accounting():
    import gzip,json
    c,products=m.load()
    result=json.loads(gzip.decompress((m.HERE/'results/potential.json.gz').read_bytes()))
    groups={}
    for r in result['comparisons']:
        key=tuple(r[k] for k in ('population_scenario','pool','name','horizon'))
        groups.setdefault(key,[]).append(r)
        assert r['joint_feed_and_water_probability']<=min(r['feed_mass']['probability_sufficient_under_model'],r['water_product']['probability_sufficient_under_model'])
    assert len(groups)==7*4*36*3
    for rows in groups.values():
        rows.sort(key=lambda r:r['available_share'])
        assert [r['available_share'] for r in rows]==c['available_resource_shares']
        assert np.all(np.diff([r['joint_feed_and_water_probability'] for r in rows])>=0)
    d=next(r for r in result['demand_checkpoints'] if r['name']=='services:breakthrough:reference:combined' and r['horizon']==200)
    stock=products['population']['scenarios']['central']['pools']['wider']['total']
    r=next(r for r in result['comparisons'] if r['population_scenario']=='central' and r['pool']=='wider' and r['name']==d['name'] and r['horizon']==200 and r['available_share']==.5)
    expected=m.comparison(stock,d,.5,c['water_recovery'],c['native_iron_recovery'])
    assert r=={'population_scenario':'central','pool':'wider',**expected}


def test_joint_required_share_uses_paired_draws():
    import gzip,json
    c,products=m.load()
    result=json.loads(gzip.decompress((m.HERE/'results/potential.json.gz').read_bytes()))
    d=next(r for r in result['demand_checkpoints'] if r['name']=='services:breakthrough:reference:combined' and r['horizon']==200)
    stock=products['population']['scenarios']['central']['pools']['wider']['total']
    required=np.maximum(d['source_feed_kg']/np.asarray(stock['body_mass_kg']),d['source_water_equivalent_kg']/(np.asarray(stock['water_equivalent_kg'])*c['water_recovery']))
    row=next(r for r in result['required_shares'] if r['scenario']=='central' and r['pool']=='wider' and r['name']==d['name'] and r['horizon']==200)
    assert np.allclose(row['joint_required_availability_share_p05_p50_p95'],np.quantile(required,[.05,.5,.95]))
    q=row['joint_required_availability_share_p05_p50_p95'][2]
    assert .949<=np.mean(required<=q)<=.951
