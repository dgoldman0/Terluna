"""Conservation, independently integrated purge and product-identity controls."""
import json

import pytest

from shared.constants import GAS_CONSTANT
from research.studies.world_economy import industrial_combustion as m


def test_complete_product_reproduction():
    assert m.results() == json.loads((m.HERE/'results/industrial_combustion.json').read_text())


def test_changed_upstream_identity_is_rejected(tmp_path):
    spec,_=m.load_inputs()
    next(iter(spec['products'].values()))['sha256']='0'*64
    path=tmp_path/'inputs.json';path.write_text(json.dumps(spec))
    with pytest.raises(ValueError,match='Input changed'):m.load_inputs(path)


def test_equal_partial_pressure_preserves_oxygen_inventory():
    r=m.results()['air'];a,b=r['reference'],r['locked']
    assert a['oxygen_mol_m3']==pytest.approx(b['oxygen_mol_m3'])
    for x in (a,b):
        assert x['density_kg_m3']*x['oxygen_mass_fraction']==pytest.approx(x['oxygen_mol_m3']*m.MOLAR_MASS['O2'])
    assert r['density_ratio']*r['oxygen_per_mass_ratio']==pytest.approx(1)
    assert 1.19<r['heat_capacity_per_volume_ratio']<1.21


def test_fixed_fraction_pressure_control_adds_oxygen():
    a=m.air(101325,301,.175);b=m.air(121590,301,.175)
    assert b['oxygen_mol_m3']/a['oxygen_mol_m3']==pytest.approx(1.2)
    assert a['oxygen_mass_fraction']==b['oxygen_mass_fraction']


@pytest.mark.parametrize('oxygen',[.175,.21])
@pytest.mark.parametrize('excess',[1.,1.1,2.])
def test_complete_combustion_conserves_each_element(oxygen,excess):
    r=m.methane(121590,600,oxygen,excess);p=r['products_mol_per_mol_CH4']
    assert p['CO2']==1
    assert 2*p['H2O']==4
    assert 2*p['CO2']+p['H2O']+2*p['O2']==pytest.approx(2*r['oxidizer_mol_per_mol_CH4']*oxygen)
    assert p['N2']==pytest.approx(r['oxidizer_mol_per_mol_CH4']*(1-oxygen))


def test_volume_and_molar_exhaust_changes_remain_distinct():
    r=m.results()['methane'][0]
    assert r['oxidizer_volume_ratio']==pytest.approx(1)
    assert r['nitrogen_per_fuel_ratio']==pytest.approx(.99/.79)
    assert r['wet_product_molar_ratio']>1
    assert r['wet_product_volume_ratio_at_equal_temperature']<1
    assert r['wet_product_volume_ratio_at_equal_temperature']*1.2==pytest.approx(r['wet_product_molar_ratio'])
    a=m.methane(101325,900,.21);b=m.methane(121590,900,.175)
    assert b['wet_product_volume_reference_m3_per_mol_CH4']/a['wet_product_volume_reference_m3_per_mol_CH4']==pytest.approx(r['wet_product_volume_ratio_at_equal_temperature'])


def test_equal_excess_gives_different_dry_exhaust_readings():
    r=m.results()['methane'][1]
    assert r['reference']['excess_oxygen_ratio']==r['locked']['excess_oxygen_ratio']
    assert r['reference']['dry_oxygen_mole_fraction']==pytest.approx(.021105527638190968)
    assert r['locked']['dry_oxygen_mole_fraction']==pytest.approx(.017283950617283963)


def test_mixed_purge_matches_independent_euler_mass_balance():
    p,t,v,y0,yf=121590,301,1.,.175,.01
    r=m.mixed_purge(p,t,v,y0,yf)
    inventory=p*v/(GAS_CONSTANT*t)
    oxygen=y0*inventory
    steps=50000
    feed=r['nitrogen_mol']/steps
    for _ in range(steps):
        oxygen-=feed*oxygen/inventory
    assert oxygen/inventory==pytest.approx(yf,rel=1e-4)


def test_purge_endpoint_changes_sign_of_gas_saving():
    shallow,deep,partial=m.results()['purges']
    assert shallow['nitrogen_moles_ratio']<1
    assert deep['nitrogen_moles_ratio']>1
    for row in (shallow,deep):
        assert row['locked']['final_oxygen_partial_pressure_Pa']/row['reference']['final_oxygen_partial_pressure_Pa']==pytest.approx(1.2)
    assert partial['locked']['final_oxygen_partial_pressure_Pa']==pytest.approx(partial['reference']['final_oxygen_partial_pressure_Pa'])
    assert partial['nitrogen_moles_ratio']==pytest.approx(1.2)


def test_purge_zero_change_and_volume_extensivity():
    assert m.mixed_purge(101325,301,1.,.21,.21)['nitrogen_mol']==0
    one=m.mixed_purge(101325,301,1.,.21,.1)
    two=m.mixed_purge(101325,301,2.,.21,.1)
    assert two['nitrogen_mol']==pytest.approx(2*one['nitrogen_mol'])
    assert one['volume_exchanges']==two['volume_exchanges']


@pytest.mark.parametrize('bad',[0.,-1.,True,float('nan'),float('inf')])
def test_invalid_pressure_rejected(bad):
    with pytest.raises(ValueError):m.air(bad,301,.175)


def test_invalid_chemistry_and_purge_directions_rejected():
    with pytest.raises(ValueError):m.methane(121590,301,.175,.9)
    with pytest.raises(ValueError):m.air(121590,301,1.)
    with pytest.raises(ValueError):m.mixed_purge(121590,301,1.,.175,.21)
