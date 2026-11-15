"""Mass, energy, renewal and provenance controls for the route comparison."""
from copy import deepcopy
import json
import math
import pytest
from . import transport as model

@pytest.fixture(scope='module')
def data():
    q,s,src=model.read_inputs();r=model.run();return q,s,src,r


def test_original_baseline_is_preserved_and_reproduced(data):
    q,s,src,r=data
    old=json.loads((model.HERE/'results/industrial_system.json').read_text())['cases']['reference']
    got=deepcopy(r['cases']['earth_chemical']);got.pop('launch_split')
    assert got==old
    assert r==json.loads((model.HERE/'results/transport.json').read_text())


def test_input_identity_refuses_stale_baseline(data):
    q=deepcopy(data[0]);q['inputs']['industrial_result']['sha256']='0'*64
    with pytest.raises(ValueError,match='hash mismatch'):model.read_inputs(q)


@pytest.mark.parametrize('dry,dv,isp',[(20000,3800,430),(160000,8000,2000),(20000,0,430)])
def test_backhauled_fuel_pays_for_both_burns(dry,dv,isp):
    cargo=200000;q=model.depot_fuelled_burns(cargo,dry,dv,isp)
    ve=isp*model.STANDARD_GRAVITY
    assert ve*math.log((dry+cargo+q['ascent_kg'])/(dry+cargo))==pytest.approx(dv)
    assert ve*math.log((dry+q['backhaul_kg']+q['descent_kg'])/(dry+q['backhaul_kg']))==pytest.approx(dv)
    assert q['supplied_kg']==pytest.approx(q['ascent_kg']+q['descent_kg'])
    assert q['handling_loss_kg']==pytest.approx(0)


def test_loss_and_margin_increase_source_demand():
    a=model.depot_fuelled_burns(200000,80000,3800,430)
    b=model.depot_fuelled_burns(200000,80000,3800,430,1.03,.98)
    assert b['supplied_kg']>a['supplied_kg']
    assert b['supplied_kg']==pytest.approx(b['ascent_kg']+b['descent_kg']+b['handling_loss_kg'])
    with pytest.raises(ValueError):model.depot_fuelled_burns(1,1,1,430,1,0)


@pytest.mark.parametrize('name',['earth_chemical','split_chemical','space_chemical','electric_earth_leg','electric_both_legs'])
def test_launch_account_reconciles_propellant_and_goods(data,name):
    q,s,src,r=data;c=r['cases'][name]
    for phase in ['peak_flows','steady_flows']:
        f=c[phase];a=c['launch_split'][phase]
        hardware=f['earth_deep_kg_year']+f['payload_destination_kg_year']['leo']
        assert a['manufactured_cargo_and_hardware_kg_year']==pytest.approx(hardware)
        assert a['LEO_kg_year']==pytest.approx(hardware+a['transfer_propellant_from_Earth_kg_year'])
        assert a['launches_at_100t']==math.ceil(a['LEO_kg_year']/s['transport']['launch_payload_kg'])
        if name not in ('earth_chemical','split_chemical'):assert a['transfer_propellant_from_Earth_kg_year']==0
    if name=='earth_chemical':
        a=c['launch_split']['peak_flows'];assert a['propellant_fraction']==pytest.approx(.6692989391874657)
        assert a['transfer_propellant_from_Earth_kg_year']==pytest.approx(1710051945.1630747)


@pytest.mark.parametrize('name',['split_chemical','space_chemical','electric_earth_leg','electric_both_legs'])
def test_source_water_and_raw_material_conserve_mass(data,name):
    q,s,src,r=data;c=r['cases'][name];z=s['source']
    for phase in ['peak_flows','steady_flows']:
        f=c[phase]
        assert f['source_water_kg_year']==pytest.approx(f['moon_water_kg_year']+f['depot_water_for_propellant_kg_year']+f['Earth_route_source_water_kg_year']+f['Source_route_water_kg_year'])
        assert f['source_ore_kg_year']*z['water_grade']*z['water_recovery']>=f['source_water_kg_year']*(1-1e-10)
        assert f['source_ore_kg_year']==pytest.approx(f['source_water_kg_year']+f['source_nitrogen_kg_year']+f['source_virgin_metal_kg_year']+f['source_retained_tails_kg_year'])
        assert f['Earth_route_source_water_kg_year']==pytest.approx(f['Earth_route_space_fuel_consumed_kg_year']+f['Earth_route_handling_loss_kg_year']+f['Earth_route_inventory_fill_kg_year'])


@pytest.mark.parametrize('name',['split_chemical','space_chemical','electric_earth_leg','electric_both_legs'])
def test_replacement_imports_and_installed_capacity_close(data,name):
    q,s,src,r=data;c=r['cases'][name];v=model.scenario(q,name)
    new,_=model.design_assets(s,c['demand'],c['peak_flows'],v)
    assert c['convergence']['relative_mass_residual']<1e-10
    for k,a in new.items():assert a['mass_kg']==pytest.approx(c['assets'][k]['mass_kg'],rel=1e-9,abs=.01)
    for dest in ['leo','surface','depot','source']:
        mass=sum(a['mass_kg']*a.get('import_fraction',1)/a['life_years'] for a in c['assets'].values() if a['destination']==dest)
        assert c['steady_flows']['payload_destination_kg_year'][dest]==pytest.approx(mass)
    for prefix,key,flow in [('Earth','Earth_vehicle','earth_deep_kg_year'),('Source','Source_vehicle','nea_payload_kg_year')]:
        if key not in c['peak_flows'] or (name=='split_chemical' and prefix=='Earth'):continue
        ship=c['peak_flows'][key];count=c['engineering']['counts'][prefix+'_route_ships']
        available=count*ship['cargo_kg']*model.JULIAN_YEAR_DAYS/ship['cycle_days']*s['transport']['fleet_availability']
        assert available>=c['peak_flows'][flow]*s['economics']['capacity_margin']
        assert sum(a['mass_kg'] for a in ship['parts'].values())==pytest.approx(ship['dry_kg'],rel=1e-9)


def test_sustaining_source_ship_keeps_installed_mass(data):
    c=data[3]['cases']['electric_both_legs']
    assert c['steady_flows']['Earth_vehicle']['dry_kg']==pytest.approx(c['peak_flows']['Earth_vehicle']['dry_kg'])
    assert c['steady_flows']['Source_vehicle']['dry_kg']==pytest.approx(c['peak_flows']['Source_vehicle']['dry_kg'])


def test_electric_energy_heat_and_renewal_are_charged(data):
    q,s,src,r=data;c=r['cases']['electric_both_legs'];eta=q['electric']['total_efficiency']
    flux=s['power']['emissivity']*model.STEFAN_BOLTZMANN*s['power']['radiator_K']**4
    for key in ['Earth_vehicle','Source_vehicle']:
        v=c['peak_flows'][key];b=v['burns'];prop=b['ascent_kg']+b['descent_kg']
        electric_J=v['electric_bus_W']*v['powered_days']*model.JULIAN_DAY
        assert electric_J*eta==pytest.approx(prop*(v['Isp_s']*model.STANDARD_GRAVITY)**2/2)
        assert v['radiator_emitting_m2']*flux==pytest.approx(v['electric_bus_W']*(1-eta))
        assert v['parts']['drive']['life_years']<v['parts']['hull_tanks']['life_years']
    assert c['assets']['Earth_solar']['life_years']==10
    assert c['assets']['Source_solar']['life_years']==20


@pytest.mark.parametrize('name',['split_chemical','space_chemical','electric_earth_leg','electric_both_legs'])
def test_budget_separates_capital_seed_and_annual_cost(data,name):
    q,s,src,r=data;c=r['cases'][name];b=c['cost'];I=b['inflation_factor'];margin=1+s['economics']['programme_contingency']
    direct=sum(b['groups_USD2041'].values())+b['preservation_USD2041']+b['research_USD2041']+b['seed_USD2041']+b['additional_route_seed_USD2041']
    assert b['capital_direct_USD2041']==pytest.approx(direct)
    assert b['capital_with_contingency_USD2041']==pytest.approx(direct*margin)
    renewal=sum(b['installed_assets_USD2026'][k]/a['life_years'] for k,a in c['assets'].items())*I
    assert b['annual_supplier_renewal_USD2041']==pytest.approx(renewal)
    annual=sum(b[k] for k in ['annual_supplier_operation_USD2041','annual_supplier_renewal_USD2041','annual_settlement_import_manufacture_USD2041','annual_settlement_import_freight_USD2041','annual_local_restoration_USD2041'])
    assert b['annual_supply_direct_USD2041']==pytest.approx(annual)
    assert b['annual_supply_with_contingency_USD2041']==pytest.approx(annual*margin)
    assert b['additional_route_seed_LEO_kg']>=0
    assert 'U_calibration' not in b


def test_low_efficiency_failure_and_sensitivity_direction(data):
    r=data[3]
    assert r['sensitivities']['electric_efficiency_15pct']['status']=='infeasible_in_this_sizing_envelope'
    improved=r['sensitivities']['electric_efficiency_50pct']['case']
    assert improved['cost']['annual_supply_with_contingency_USD2041']<r['cases']['electric_both_legs']['cost']['annual_supply_with_contingency_USD2041']
    # A different flight duration could restore closure; this rejects only the
    # stated mass/power/flight-time combination, not electric propulsion itself.


def test_localization_keeps_refining_and_factory_demands(data):
    r=data[3]
    for name,c in r['source_structure_localization'].items():
        assert c['peak_flows']['source_metal_kg_year']>0
        assert c['assets']['source_metal_trains']['mass_kg']>0
        assert c['assets']['source_fabrication_cells']['mass_kg']>0
        assert c['assets']['asteroid_tugs']['import_fraction']==pytest.approx(.335)
        assert c['demand']==r['cases'][name]['demand']


def test_route_inventory_is_filled_from_space_during_construction(data):
    q,s,src,r=data
    for name,c in r['cases'].items():
        if name=='earth_chemical':continue
        fill=c['peak_flows']['Earth_route_inventory_fill_kg_year']*s['schedule']['delivery_equivalent_years']
        assert fill==pytest.approx(c['engineering']['additional_route_storage_kg'])
        assert c['steady_flows']['Earth_route_inventory_fill_kg_year']==0
        assert c['cost']['additional_route_seed_LEO_kg']<fill


def test_split_refuelling_uses_depot_return_fuel_and_retains_earth_hull(data):
    s=data[1];r=data[3];c=r['cases']['split_chemical']
    assert c['peak_flows']['Earth_vehicle']['dry_kg']==s['transport']['earth_payload_kg']*s['transport']['earth_dry_to_payload']
    assert c['peak_flows']['Earth_route_fuel_backhaul_kg_year']==0
    assert c['peak_flows']['Earth_route_source_water_kg_year']>0
    assert c['launch_split']['peak_flows']['propellant_fraction']>0
    assert c['launch_split']['peak_flows']['LEO_kg_year']<r['cases']['earth_chemical']['launch_split']['peak_flows']['LEO_kg_year']
