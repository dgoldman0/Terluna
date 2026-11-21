"""Numerical and accounting controls; no technology qualification claims."""
from copy import deepcopy
import json
import math
from pathlib import Path
import pytest
from . import run as model

@pytest.fixture(scope='module')
def data():
    s,src=model.read_inputs();r=model.run();return s,src,r


def test_bound_input_identity_is_enforced(data):
    s,_,_=data;q=deepcopy(s);q['inputs']['settlement']['sha256']='0'*64
    with pytest.raises(ValueError,match='hash mismatch'):model.read_inputs(q)


def test_committed_result_reproduces(data):
    assert data[2]==json.loads((model.HERE/'results/industrial_system.json').read_text())


@pytest.mark.parametrize('dv,isp,dry,back',[(0,430,.4,0),(2800,430,.4,0),(2500,430,.1,.06)])
def test_each_rocket_burn_conserves_mass_and_recovers_delta_v(dv,isp,dry,back):
    q=model.round_trip(dv,isp,dry,back);ve=isp*model.STANDARD_GRAVITY
    out_end=1+dry+q['returning'];out_start=out_end+q['outbound']
    assert ve*math.log(out_start/out_end)==pytest.approx(dv)
    assert ve*math.log((dry+back+q['returning'])/(dry+back))==pytest.approx(dv)
    assert q['total']==q['outbound']+q['returning']


def test_water_hydrogen_and_oxygen_have_no_unassigned_coproduct(data):
    s,_,r=data
    for c in r['cases'].values():
        f=c['peak_flows'];H=f['depot_hydrogen_kg_year'];O=f['depot_oxygen_to_fuel_kg_year']
        airO=f['moon_air_kg_year']*(1-s['source']['air_nitrogen_mass_fraction'])
        assert H*9==pytest.approx(f['depot_water_for_propellant_kg_year'])
        assert H*8==pytest.approx(O+airO)
        assert H+O==pytest.approx(f['lunar_propellant_kg_year'])
        assert f['depot_excess_oxygen_kg_year']==0
        assert f['source_water_kg_year']==pytest.approx(f['moon_water_kg_year']+f['depot_water_for_propellant_kg_year']+f['nea_propellant_kg_year'])


def test_ore_supply_meets_both_chemical_requirements(data):
    s,_,r=data
    for name,c in r['cases'].items():
        q=model.with_case(s,name)['source']; f=c['peak_flows'];ore=f['source_ore_kg_year']
        assert ore*q['water_grade']*q['water_recovery']>=f['source_water_kg_year']*(1-1e-12)
        assert ore*q['nitrogen_grade']*q['nitrogen_recovery']>=f['source_nitrogen_kg_year']*(1-1e-12)
        assert ore==pytest.approx(f['source_water_kg_year']+f['source_nitrogen_kg_year']+f['source_virgin_metal_kg_year']+f['source_retained_tails_kg_year'])


def test_supplier_mass_iteration_closes_at_delivered_assets(data):
    s,_,r=data
    for name,c in r['cases'].items():
        q=model.with_case(s,name);new,_=model.design_assets(q,c['demand'],c['peak_flows'])
        assert c['convergence']['relative_mass_residual']<1e-10
        for k,a in new.items():assert a['mass_kg']==pytest.approx(c['assets'][k]['mass_kg'],rel=1e-9,abs=.01)
        assert c['peak_flows']['nea_return_hardware_kg_year']<c['peak_flows']['nea_payload_kg_year']


def test_installed_fleets_have_required_capacity(data):
    s,_,r=data
    for name,c in r['cases'].items():
        q=model.with_case(s,name);t=q['transport'];n=c['engineering']['counts'];f=c['peak_flows'];m=q['economics']['capacity_margin']
        for unit,key,payload,cycles in [('lunar_landers','moon_payload_kg_year',t['lander_payload_kg'],t['lander_cycles_year']),('asteroid_tugs','nea_payload_kg_year',t['nea_payload_kg'],1/t['nea_cycle_years']),('Earth_transfer_tugs','earth_deep_kg_year',t['earth_payload_kg'],t['earth_cycles_year'])]:
            assert n[unit]*payload*cycles*t['fleet_availability']>=f[key]*m


def test_shared_plant_and_power_are_counted_once(data):
    s,src,r=data;c=r['cases']['reference'];d=c['demand'];e=c['engineering']
    assert d['earth_special']+s['lunar']['old_embedded_plant_kg']==src['settlement']['cases']['reference']['resources']['special_imports_kg']
    assert c['assets']['construction_and_service_robots']['mass_kg']==s['lunar']['old_embedded_plant_kg']
    assert e['total_ground_firm_W']==pytest.approx(d['old_total_W']+max(0,e['local_industry_firm_W']-d['old_construction_W']))


def test_capital_and_renewal_budget_sum_without_internal_sales(data):
    s,_,r=data
    for name,c in r['cases'].items():
        q=model.with_case(s,name);b=c['cost'];margin=1+q['economics']['programme_contingency']
        direct=sum(b['groups_USD2041'].values())+b['preservation_USD2041']+b['research_USD2041']+b['seed_USD2041']
        assert direct==pytest.approx(b['capital_direct_USD2041'])
        assert direct*margin==pytest.approx(b['capital_with_contingency_USD2041'])
        campaign=direct+sum(b[k] for k in ('initial_stock_manufacture_USD2041','initial_stock_Earth_freight_USD2041','initial_stock_restoration_USD2041','campaign_supplier_operation_and_renewal_USD2041','construction_period_settlement_makeup_USD2041'))
        assert campaign*margin==pytest.approx(b['industrial_and_initial_material_campaign_USD2041'])
        assert b['U_calibration']['basket_USD2041']/b['U_calibration']['basket_U']==pytest.approx(b['U_calibration']['USD2041_per_U'])


def test_inflation_changes_prices_only(data):
    s,src,r=data;q=deepcopy(s);q['economics']['inflation']=0;c=model.case(q,src);b=r['cases']['reference']
    assert b['peak_flows']==c['peak_flows']
    assert b['cost']['capital_with_contingency_USD2041']==pytest.approx(c['cost']['capital_with_contingency_USD2041']*1.03**15)


def test_resource_and_launch_sensitivities_are_distinct(data):
    _,_,r=data;ref=r['cases']['reference'];ss=r['sensitivities']
    assert ss['launch_50']['source_ore_kg_year']==ref['peak_flows']['source_ore_kg_year']
    assert ss['launch_50']['campaign_USD2041']<ref['cost']['industrial_and_initial_material_campaign_USD2041']
    assert ss['water_makeup_quarter']['source_ore_kg_year']<ref['peak_flows']['source_ore_kg_year']
    assert ss['stock_delivery_10_years']['capital_USD2041']>ss['stock_delivery_40_years']['capital_USD2041']


def test_beam_geometry_and_eclipse_reserve_screen(data):
    s,_,r=data;e=r['cases']['reference']['engineering']
    assert e['sources_after_one_failure_and_one_lunar_eclipse']>=2
    assert e['full_occupied_load_reserve_hours']>=s['power']['common_eclipse_hours']
    # Bound states only the ideal geometry; it does not prove orbital stability.
    assert e['beam_converter_radiator_m2']*s['power']['radiator_kg_m2']<r['cases']['reference']['assets']['lunar_beam_constellation']['mass_kg']


def test_electric_tug_counts_its_power_and_radiators(data):
    s,_,r=data;a,b=r['electric_tug_comparison']
    assert a['dry_kg']>s['electric_tug_comparison']['structural_dry_to_payload']*a['cargo_kg']
    assert a['power_W']>b['power_W'] and a['dry_kg']>b['dry_kg']
    assert a['propellant_water_kg']<r['cases']['reference']['peak_flows']['nea_round_trip']['total']*a['cargo_kg']


def test_source_fabrication_counts_plant_energy_metal_and_import_credit(data):
    s,_,r=data; a=r['cases']['reference'];b=r['cases']['local_fleet_fabrication']
    assert b['peak_flows']['source_metal_kg_year']>0
    assert b['assets']['source_metal_trains']['mass_kg']>0
    assert b['assets']['source_fabrication_cells']['mass_kg']>0
    assert b['peak_flows']['source_virgin_metal_kg_year']<b['peak_flows']['source_metal_kg_year']
    assert b['cost']['capital_with_contingency_USD2041']<a['cost']['capital_with_contingency_USD2041']
    assert b['assets']['asteroid_tugs']['import_fraction']==pytest.approx(.335)
