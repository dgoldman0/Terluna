"""Conservation, scaling and accounting controls for the power comparison."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from provisioning.enclosed_power import model as m

CFG = json.loads((Path(m.__file__).parent / 'enclosed_power_inputs.json').read_text())
ENG = m.read_pinned(CFG['inputs']['engineering'])
CONSTANTS = m.read_pinned(CFG['inputs']['constants'])
TECH = CFG['technologies']['reference_assumptions']
PERIOD = CONSTANTS['moon']['synodic_month_days'] * 24


def case(name, **kwargs):
    return m.architecture(name, TECH, CFG, ENG, CONSTANTS, **kwargs)


def priced(item, imported=1, cfg=CFG):
    return m.cost(item, cfg, CONSTANTS, cfg['common']['freight_scenarios']['reference'], imported)


def test_ordinary_dispatch_matches_closed_form():
    sun, dt = m.profile(PERIOD, 4096)
    d = m.dispatch(sun, dt, np.ones(len(sun)), 0, .75, .95, .4675)
    assert d['PV_nameplate_kW'] == pytest.approx((1 + 1/.4675)/(.75*.95))
    assert d['storage_output_kWh'] == pytest.approx(PERIOD/2)
    assert d['solar_to_bus_kWh'] == pytest.approx(d['cycle_load_kWh'] + d['storage_loss_kWh'])
    assert abs(d['cycle_residual_kWh']) < 1e-8


@pytest.mark.parametrize('site', ['ordinary', 'polar_sensitivity'])
def test_load_shift_preserves_energy_and_reduces_storage(site):
    sun, dt = m.profile(PERIOD, 4096, site)
    original = m.dispatch(sun, dt, m.shifted_load(sun, 0), 0, .75, .95, .4675)
    shifted = m.dispatch(sun, dt, m.shifted_load(sun, .5), 0, .75, .95, .4675)
    assert shifted['cycle_load_kWh'] == pytest.approx(PERIOD)
    assert shifted['storage_output_kWh'] == pytest.approx(original['storage_output_kWh']/2)


def test_hybrid_conserves_energy_with_half_storage():
    a, b = case('solar_RFC')['details'], case('hybrid_50')['details']
    assert b['storage_output_kWh'] == pytest.approx(a['storage_output_kWh']/2)
    assert b['solar_to_bus_kWh'] + b['nuclear_to_bus_kWh'] == pytest.approx(b['cycle_load_kWh']+b['storage_loss_kWh'])


def test_rfc_counts_one_conserved_inventory_and_all_tanks():
    a = case('solar_RFC'); d = a['details']
    assert d['oxygen_kg'] == pytest.approx(8*d['hydrogen_kg'])
    assert d['water_equivalent_inventory_kg'] == pytest.approx(9*d['hydrogen_kg'])
    assert d['hydrogen_kg']*39.4*TECH['FC_HHV_efficiency'] == pytest.approx(d['storage_output_kWh'])
    assert d['tank_mass_kg'] == pytest.approx(15*d['hydrogen_kg'] + d['oxygen_kg'] + .1*d['water_equivalent_inventory_kg'])
    inventory = [x for x in a['components_per_mean_kW'] if x['name']=='RFC_inventory']
    assert len(inventory)==1 and inventory[0]['life_years'] is None


def test_thermal_inventory_closes_energy_and_stays_local():
    a = case('solar_thermal'); d = a['details']; c = priced(a)
    recovered = d['thermal_media_kg']*TECH['thermal_store_output_kWhth_kg']*TECH['thermal_store_retention']*TECH['thermal_engine_efficiency']*CFG['common']['grid_efficiency']
    assert recovered == pytest.approx(d['storage_output_kWh'])
    assert d['solar_to_bus_kWh'] == pytest.approx(d['cycle_load_kWh']+d['storage_loss_kWh'])
    assert c['always_local_bulk_kg_per_mean_kW'] == pytest.approx(d['thermal_media_kg'])
    assert c['imported_hardware_and_fluid_kg_per_mean_kW']+c['always_local_bulk_kg_per_mean_kW'] == pytest.approx(a['mass_kg_per_mean_kW'])
    zero = priced(a, 0)
    expected_difference = c['imported_hardware_and_fluid_kg_per_mean_kW']*(1000-20)
    assert c['installed_supply_and_delivery_USD2026_per_mean_kW']-zero['installed_supply_and_delivery_USD2026_per_mean_kW'] == pytest.approx(expected_difference)


def test_beam_utilization_changes_orbital_capacity_only():
    a, b = case('beamed_legacy', beam_useful=.2), case('beamed_legacy', beam_useful=.4)
    for x, y in zip(a['components_per_mean_kW'], b['components_per_mean_kW']):
        assert x['mass_kg'] == pytest.approx(y['mass_kg']*(2 if x['destination']=='orbit' else 1))


def test_small_beam_geometry_does_not_shrink_with_power():
    a, b = case('beamed_legacy', mean_MW=.01), case('beamed_legacy', mean_MW=1)
    for item in [a, b]:
        assert item['details']['aperture_floor_binding']
    masses = [sum(x['mass_kg'] for x in z['components_per_mean_kW'] if x['destination']=='orbit')*z['mean_MW']*1000 for z in [a,b]]
    assert masses[0] == pytest.approx(masses[1])


def test_permanent_inventory_not_rebought_at_plant_renewal():
    cfg = deepcopy(CFG)
    cfg['common']['annual_upkeep_fraction'] = 0
    cfg['common']['annual_thermal_media_replacement_fraction'] = 0
    cfg['common']['annual_RFC_inventory_loss_fraction'] = 0
    a = {'components_per_mean_kW':[{'name':'thermal_media', 'mass_kg':10., 'supply_USD2026':10., 'life_years':None,'destination':'local_bulk'}], 'supply_USD2026_per_mean_kW':10., 'nuclear_mean_share':0}
    result = priced(a, cfg=cfg)
    assert all(v['cost_PV_USD2026_per_mean_kW']==12 for v in result['lifecycle'].values())


def test_replacement_at_horizon_end_is_not_charged():
    cfg = deepcopy(CFG); cfg['common']['annual_upkeep_fraction']=0
    a = {'components_per_mean_kW':[{'name':'hardware', 'mass_kg':0., 'supply_USD2026':100., 'life_years':20,'destination':'surface'}], 'supply_USD2026_per_mean_kW':100., 'nuclear_mean_share':0}
    assert priced(a, cfg=cfg)['lifecycle']['20']['cost_PV_USD2026_per_mean_kW'] == 100


def test_polar_resolution_changes_storage_by_less_than_one_percent():
    values=[]
    for n in [4096,8192]:
        sun, dt = m.profile(PERIOD,n,'polar_sensitivity')
        values.append(m.dispatch(sun,dt,np.ones(n),0,.75,.95,.4675))
    for key in ['storage_output_kWh','PV_nameplate_kW']:
        assert abs(values[0][key]/values[1][key]-1)<.01


def test_product_integrity_and_mapped_mean_load():
    p = json.loads((Path(m.__file__).parent/'results/power_comparison.json').read_text())
    assert p['producer']['sha256']==m.digest(m.__file__)
    for dep in p['dependencies'].values():
        assert m.digest(m.ROOT/dep['path'])==dep['sha256']
    assert len([x for x in p['physical_load_mappings'] if x['year']==200])==81
    for row in p['physical_load_mappings']:
        occupied=ENG['demand']['occupied_mean_W']*row['population']/ENG['normalization']['population']/1e6
        process=row['mean_MW']*row['maximum_shift_if_all_current_lunar_processing_moves']
        assert occupied+process==pytest.approx(row['mean_MW'])
        expected_firm=ENG['demand']['occupied_W']*row['population']/ENG['normalization']['population']/1e6+process*1.25/.8
        assert expected_firm==pytest.approx(row['firm_MW'])


def test_reference_fission_radiator_fits_within_complete_mass():
    d=case('fission')['details']
    assert d['nuclear_radiator_mass_within_system_kg'] < d['nuclear_rated_kW']*TECH['nuclear_kg_kW']


def test_invalid_option_and_negative_imports_fail():
    with pytest.raises(ValueError):case('unpriced_fusion')
    with pytest.raises(ValueError):priced(case('fission'),-1)
