"""Geometric, probabilistic and cost-unit checks for wider lunar potential."""
import gzip,json
from copy import deepcopy
import numpy as np
import pytest
from resources.source_access import lunar_potential as m


def test_column_mass_integrates_density_and_respects_limits():
    for depth in (.01,1.,5.,10.):
        z=np.linspace(0,depth,100001)
        rho=1800-700*np.exp(-z/.06)
        calculated=m.column_mass(depth,1100,1800,.06)
        assert calculated==pytest.approx(np.trapezoid(rho,z),rel=1e-7)
        assert 1100*depth<calculated<1800*depth
    assert m.column_mass(5,1800,1800,.06)==9000
    for args in ((0,1100,1800,.06),(5,1100,1000,.06),(5,1100,1800,0)):
        with pytest.raises(ValueError):m.column_mass(*args)


def test_ranked_prospect_is_fractional_area_not_extra_inventory():
    pixels=[{'abundances':{'FeO':g},'solid_angle_sr':w} for g,w in ((.1,1),(.2,2),(.3,1))]
    before=deepcopy(pixels)
    result=m.ranked_pool(pixels,{'oxide':'FeO','top_area_fraction':.5})
    assert sum(p['solid_angle_sr'] for p in result)==2
    assert [p['abundances']['FeO'] for p in result]==[.3,.2]
    assert pixels==before
    assert m.ranked_pool(pixels,None) is pixels


def test_joint_draws_keep_covariance_and_unnormalized_assays():
    cov=np.eye(9)*1e-6;cov[0,1]=cov[1,0]=.7e-6
    pixel={'abundances':dict(zip(m.lunar.OXIDES,[.2]*6)),'covariance':cov.tolist(),'solid_angle_sr':1}
    draws,diagnostic=m.regional_draws([pixel],20000,np.random.default_rng(42))
    assert np.mean(draws.sum(axis=1))==pytest.approx(1.2,abs=.0001)
    assert np.corrcoef(draws[:,0],draws[:,1])[0,1]==pytest.approx(.7,abs=.025)
    assert diagnostic['rejected_proposals']==0
    pixel['abundances']['TiO2']=0
    draws,diagnostic=m.regional_draws([pixel],2000,np.random.default_rng(42))
    assert np.all(draws>=0) and diagnostic['rejected_proposals']>0


def test_actual_snapshot_and_potential_ordering():
    config=json.loads((m.HERE/'lunar_potential_inputs.json').read_text());config['draws']=4096
    x=m.product(config)
    total=[row['nonpolar_bulk_mass_kg_p05_p50_p95'][1] for row in x['bulk_inventory']]
    assert total[0]<total[1]<total[2]
    assert x['nonpolar_area_km2']==pytest.approx(37_356_050.802,rel=1e-8)
    pools=x['regional_pools'];assert pools['iron_prospect']['area_km2']/x['nonpolar_area_km2']==pytest.approx(.25)
    assert pools['iron_prospect']['composition']['FeO']['grade_p05_p50_p95'][1]>pools['nonpolar']['composition']['FeO']['grade_p05_p50_p95'][1]
    assert pools['titanium_prospect']['composition']['TiO2']['handling_USD2041_per_kg_oxide_equivalent_p05_p50_p95'][1]<pools['nonpolar']['composition']['TiO2']['handling_USD2041_per_kg_oxide_equivalent_p05_p50_p95'][1]
    trials=pools['iron_prospect']['composition']['FeO']['depth_transfer_scenarios']
    assert trials[2]['at_5m_per_1000km2_recovered_oxide_equivalent_kg_p05_p50_p95'][1]==pytest.approx(3*trials[0]['at_5m_per_1000km2_recovered_oxide_equivalent_kg_p05_p50_p95'][1])
    assert x['extraction_authorized'] is False and x['site_specific_exclusion_area_km2'] is None
    for name,expected in x['dependencies'].items():assert m.sha(m.HERE/name)==expected
