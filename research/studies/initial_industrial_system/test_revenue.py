"""Independent accounting, capacity, return-flight and source-stock controls."""
from copy import deepcopy
import json
import math
import pytest
from . import revenue as model
from . import resource_inventory as inventory


@pytest.fixture(scope='module')
def data():
    q,src=model.read_inputs()
    settle=src['settlement']['cases']['reference']['resources']
    q['_operating_settlement_Earth_kg']=sum(v for k,v in settle['operating_import_kg_year'].items() if k not in ('water','air'))
    q['_renewal_settlement_Earth_kg']=sum(settle['renewal_import_kg_year'].values())
    q['_renewal_settlement_pressure_kg']=sum(v for k,v in settle['renewal_import_kg_year'].items() if 'pressure' in k)
    return q,src,model.run()


def test_product_reproduces_and_refuses_changed_input(data):
    q,_,r=data
    assert r==json.loads((model.HERE/'results/revenue.json').read_text())
    bad=deepcopy(q);bad['inputs']['transport']['sha256']='0'*64
    with pytest.raises(ValueError,match='hash mismatch'):model.read_inputs(bad)


@pytest.mark.parametrize('mode',['earth_chemical','split_chemical','space_chemical','electric'])
def test_independent_rocket_equations_with_return_goods(mode):
    cargo,goods,dry,dv,isp=200000,20000,80000,3800,430
    a=model.trip_burns(cargo,goods,dry,dv,isp,1,mode)
    ve=isp*model.STANDARD_GRAVITY
    outbound_final=cargo+dry+(a['back'] if mode=='earth_chemical' else 0)
    return_final=goods+dry+(a['out'] if mode not in ('earth_chemical','split_chemical') else 0)
    assert ve*math.log((outbound_final+a['out'])/outbound_final)==pytest.approx(dv)
    assert ve*math.log((return_final+a['back'])/return_final)==pytest.approx(dv)
    empty=model.trip_burns(cargo,0,dry,dv,isp,1,mode)
    assert a['burned']>empty['burned']
    assert a['Earth_fuel']+a['space_supply']==pytest.approx(a['burned'])


def test_backhaul_handling_losses_are_not_free():
    a=model.trip_burns(200000,10000,80000,3800,430,1.03,'space_chemical',.98)
    assert a['space_supply']>a['burned']
    assert a['Earth_fuel']==0


def test_return_payload_reduces_loading_when_tank_is_full():
    kwargs=dict(payload=100000,dry=40000,dv=2800,isp=430,margin=1.03,mode='earth_chemical',retention=1)
    tank=model.trip_burns(100000,0,40000,2800,430,1.03,'earth_chemical')['tank']
    a=model.ferry(1000000,0,tank_limit=tank,**kwargs)
    b=model.ferry(1000000,100000,tank_limit=tank,**kwargs)
    assert b['trips_year']>a['trips_year']
    assert b['outbound_kg_trip']<100000
    assert b['tank_kg_trip']<=tank*(1+1e-10)


@pytest.mark.parametrize('name',['earth_chemical','space_chemical','electric_earth_leg','electric_both_legs','split_chemical','local_structures_split_chemical','local_structures_space_chemical','local_structures_electric_earth_leg','electric_development','favourable_earth_chemical'])
def test_each_supply_case_closes_finance_mass_and_reserved_capacity(data,name):
    q,src,r=data;s,base,tc=model.supply_cases(src)[name];z=s['source']
    for c in r['cases'][name].values():
        assert c['gross_receipts_USD2041']==pytest.approx(sum(c['receipts_by_activity_USD2041'].values()))
        assert c['external_operating_cost_USD2041']==pytest.approx(sum(c['operating_cost_components_USD2041'].values()))
        assert c['external_renewal_USD2041']==pytest.approx(sum(c['renewal_components_USD2041'].values()))
        assert c['operating_cash_balance_USD2041']==pytest.approx(c['gross_receipts_USD2041']-c['external_operating_cost_USD2041'])
        assert c['balance_after_renewal_and_resource_provision_USD2041']==pytest.approx(c['operating_cash_balance_USD2041']-c['external_renewal_USD2041']-c['resource_transition_provision_USD2041'])
        assert c['debt_service_USD2041'] is None
        assert c['annual_procurement_contingency_stress_USD2041']==pytest.approx(.35*(c['external_operating_cost_USD2041']+c['external_renewal_USD2041']))
        for f in (c['physical'],c['commercial_commissioning']['annual_combined_flows']):
            assert max(f['capacity_use_fraction_of_rated'].values())<=1+1e-8
            assert f['source_ore_kg_year']==pytest.approx(f['source_water_kg_year']+f['source_nitrogen_kg_year']+f['source_virgin_metal_kg_year']+f['source_retained_tails_kg_year'])
            assert f['source_ore_kg_year']*z['water_grade']*z['water_recovery']>=f['source_water_kg_year']*(1-1e-10)
            for route in ('Earth','Moon','source'):
                v=f[route]
                assert v['tank_kg_trip']<=v['tank_limit_kg']*(1+1e-10)
                assert v['trips_year']*v['outbound_kg_trip']==pytest.approx(f[{'Earth':'Earth_cargo_kg_year','Moon':'lunar_cargo_kg_year','source':'source_cargo_kg_year'}[route]])
            assert f['Moon']['Earth_fuel_kg_year']==f['source']['Earth_fuel_kg_year']==0
            if tc['earth_route'] not in ('earth_chemical','split_chemical'):assert f['Earth']['Earth_fuel_kg_year']==0
        assert c['resident_export_and_industrial_jobs']<q['boundary']['population']


def test_idle_commerce_does_not_add_phantom_freight(data):
    q,src,r=data
    b={k:0 for k in ['operating_surface_kg','operating_LEO_kg','return_packed_kg','guest_water_kg_year','guest_air_kg_year','propellant_customer_kg','initial_metal_kg','initial_bulk_kg','local_export_metal_kg','added_surface_mean_W']}
    b['renewal_mass_kg']={'surface':0,'leo':0}
    for name,(s,c,tc) in model.supply_cases(src).items():
        f=model.logistics(s,c,tc,b)
        assert f['LEO_kg_year']==pytest.approx(c['steady_flows']['LEO_kg_year'])
        assert f['source_water_kg_year']==pytest.approx(c['steady_flows']['source_water_kg_year'],rel=1e-10)
        assert f['raw_lunar_kg_year']==pytest.approx(c['steady_flows']['raw_lunar_kg_year'])


def test_capsule_fleet_shields_relaunch_and_expected_life_are_distinct(data):
    q,src,r=data;c=r['cases']['local_structures_space_chemical']['reference'];b=c['commercial'];rr=q['return']
    s,_,_=model.supply_cases(src)['local_structures_space_chemical']
    assert b['goods_delivered_saleable_kg']==pytest.approx(15280650)
    assert b['return_flights_year']==pytest.approx(b['return_packed_kg']/rr['capsule_payload_kg'])
    expected_life=sum((1-rr['loss_probability_flight'])**n for n in range(rr['durable_life_flights']))
    assert b['expected_capsule_life_flights']==pytest.approx(expected_life)
    bus=b['assets']['return_capsule_buses'];kg_bus=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*(1-rr['shield_fraction_of_dry'])
    assert bus['mass_kg']/bus['life_years']==pytest.approx(b['return_flights_year']*kg_bus/expected_life)
    kg_shield=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*rr['shield_fraction_of_dry']
    prop=(rr['capsule_payload_kg']+kg_bus+kg_shield)*math.expm1(rr['deorbit_dv_m_s']/(rr['deorbit_Isp_s']*model.STANDARD_GRAVITY))
    assert b['operating_LEO_kg']==pytest.approx(b['return_flights_year']*(kg_bus+kg_shield+prop))
    assert b['operating_purchases_USD2026']['entry_shields']==pytest.approx(b['return_flights_year']*kg_shield*s['prices']['spacecraft_USD2026_kg'])
    # Durable bus manufacture is in renewal; its launch already occurs in sorties.
    assert b['renewal_mass_kg']['leo']==pytest.approx(sum(a['mass_kg']/a['life_years'] for k,a in b['assets'].items() if a['destination']=='leo' and k!='return_capsule_buses'))


def test_factory_material_yield_and_new_guest_food_are_explicit(data):
    q,src,r=data;c=r['cases']['local_structures_space_chemical']['reference'];b=c['commercial']
    assert b['factory_input_kg']==pytest.approx(b['goods_made_kg']+b['factory_rejects_kg'])
    assert b['local_export_metal_kg']==pytest.approx(b['factory_input_kg']*.8+b['goods_made_kg']*.05)
    assert b['guest_years']==pytest.approx(10000*14/model.JULIAN_YEAR_DAYS)
    s,base,tc=model.supply_cases(src)['local_structures_space_chemical']
    assert b['guest_water_kg_year']==pytest.approx(base['demand']['yearly_water']*b['guest_years']/250000)
    assert b['assets']['visitor_life_support_and_food']['mass_kg']==2000*b['visitor_beds']
    assert b['initial_air_kg']>0 and b['initial_water_kg']>0
    ramp=c['commercial_commissioning'];f=ramp['annual_combined_flows']
    assert f['lunar_cargo_kg_year']-c['physical']['lunar_cargo_kg_year']==pytest.approx((b['initial_mass_kg']['surface']+b['initial_air_kg']+b['initial_water_kg'])/ramp['years'])
    assert ramp['additional_source_feed_kg']>ramp['additional_source_water_kg']>0
    assert ramp['initial_resource_provision_USD2041']>0


def test_batch_computing_pays_for_own_power_heat_rejection_and_replacement(data):
    q,src,r=data;c=r['cases']['local_structures_space_chemical']['reference'];b=c['commercial'];a=b['assets'];s,_,_=model.supply_cases(src)['local_structures_space_chemical'];p=s['power']
    assert a['batch_solar']['mass_kg']*p['space_system_W_kg']*p['end_of_life_factor']==pytest.approx(b['compute_peak_W'])
    area=a['batch_radiators']['mass_kg']/p['radiator_kg_m2']
    assert area*p['emissivity']*model.STEFAN_BOLTZMANN*p['radiator_K']**4==pytest.approx(b['compute_peak_W'])
    assert a['batch_IT']['life_years']==5
    assert b['receipts_USD2026']['batch_computing']==pytest.approx(200e6*.6*30)


def test_losses_are_suppressed_and_no_own_fuel_sale_is_revenue(data):
    q,src,r=data
    rejected=r['cases']['earth_chemical']['reference']
    assert rejected['commercial']['goods_made_kg']==0
    assert any(a['activity']=='manufacturing' and a['incremental_balance_after_renewal_and_resource_provision_USD2041']<0 for a in rejected['rejected_activities'])
    focus=r['cases']['local_structures_space_chemical']['reference'];b=focus['commercial']
    assert b['receipts_USD2026']['Earth_mission_propellant']==25e6*300
    assert b['propellant_customer_kg']<focus['physical']['source_water_kg_year']
    assert not q['boundary']['internal_transactions_are_revenue']
    assert b['propellant_sale_USD2026_kg']<=.8*b['Earth_fuel_alternative_USD2026_kg']


def test_normalization_preserves_USD_and_inherits_only_scoped_U_prices(data):
    q,src,r=data;old=src['early_exports'];v=q['visitors'];f=q['manufacturing']
    assert f['nameplate_finished_kg_per_resident_year']==old['factory_nameplate_finished_kg_per_resident_year']
    assert f['retained_U_price_per_kg']==old['factory_finished_price_U_per_kg']
    assert v['retained_U_per_night']==old['stages'][0]['nightly_destination_spend_U']
    assert r['U_calibration']['destination_service_USD2041_per_U_by_market_price_case']['reference']==pytest.approx(20*1.03**15)
    changed=deepcopy(q);changed['normalization']['proposed_USD2041_per_U']=60
    twice=model.run(changed)
    for name,groups in r['cases'].items():
        for market,c in groups.items():
            other=twice['cases'][name][market]
            assert c['balance_after_renewal_and_resource_provision_USD2041']==other['balance_after_renewal_and_resource_provision_USD2041']
            assert c['proposed_normalization_U']['gross_receipts_U']==pytest.approx(other['proposed_normalization_U']['gross_receipts_U']*2)


def test_resource_horizons_are_conditional_inverse_flow_and_mass_not_reserves(data):
    q,src,r=data;c=r['cases']['local_structures_space_chemical']['reference'];h=c['resource_horizons'];flow=c['physical']['source_ore_kg_year']
    for pool,years in h['years_for_specified_usable_feed_pools'].items():assert years*flow==pytest.approx(float(pool))
    for years,req in h['requirements'].items():
        assert req['raw_feed_kg']==pytest.approx(flow*int(years))
        assert math.pi/6*req['equivalent_source_diameter_m']**3*q['resources']['equivalent_body_density_kg_m3']==pytest.approx(req['raw_feed_kg'])
    assert 'not certified' in h['rule']
    assert src['inventory']['inventory']['without_diameter']>src['inventory']['inventory']['with_diameter']


def test_catalogue_filter_keeps_missing_values_and_named_targets_separate():
    fields=['spkid','full_name','diameter','albedo','a','e','i']
    d={'count':5,'fields':fields,'data':[['1','ordinary','1','.05','1','.1','5'],['2','Bennu','2','.04','1','.1','5'],['3','unknown',None,'.05','1','.1','5'],['4','bright','1','.3','1','.1','5'],['5','distant','1','.05','2','.7','30']]}
    a=inventory.summarize(d,1500);b=inventory.summarize(d,3000)
    assert a['with_diameter']==4 and a['without_diameter']==1
    assert a['groups']['dark_geometry_filter']['count']==2
    assert a['groups']['dark_geometry_excluding_named_science_targets']['count']==1
    assert a['groups']['all_with_diameter']['volume_m3']==pytest.approx(math.pi/6*1e9*11)
    assert b['groups']['all_with_diameter']['mass_at_assumed_density_kg']==2*a['groups']['all_with_diameter']['mass_at_assumed_density_kg']
    d['count']=6
    with pytest.raises(ValueError,match='Incomplete'):inventory.summarize(d)


def test_archived_catalogue_reproduces_and_hash_mismatch_stays_an_error(tmp_path):
    if not (inventory.ARCHIVE/'nea_catalogue.json').exists():pytest.skip('Pinned external JPL snapshot is not installed')
    assert inventory.run()==json.loads((model.HERE/'results/resource_inventory.json').read_text())
    bad=tmp_path/'wrong.json';bad.write_text('{}')
    with pytest.raises(ValueError,match='snapshot hash mismatch'):inventory.run(bad)


def test_external_renewal_uses_manufacture_and_freight_once(data):
    q,src,r=data
    for name in ['earth_chemical','local_structures_space_chemical','electric_development']:
        s,base,tc=model.supply_cases(src)[name];c=r['cases'][name]['reference'];I=base['cost']['inflation_factor'];t=s['transport']
        old=base['steady_flows'];prop=old.get('earth_launch_propellant_kg_year',old['earth_transfer_propellant_kg_year'])/old['earth_deep_kg_year']
        freight=(1+prop)*t['LEO_USD2026_kg']+prop*t['earth_propellant_USD2026_kg']
        replacement=0
        for key,a in base['assets'].items():
            manufacture=a['manufacture_USD2026']
            if key=='additional_ground_radiator_controls':manufacture*=q['finance']['ground_radiator_assembly_Earth_fraction']
            # Local tug manufacture was already excluded in its asset price;
            # only the remaining mass is launched. Do not discount twice.
            transport=a['mass_kg']*a.get('import_fraction',1)*(t['LEO_USD2026_kg'] if a['destination']=='leo' else freight)
            replacement+=(manufacture+transport)/a['life_years']
        assert c['renewal_components_USD2041']['baseline_supplier_replacement']==pytest.approx(replacement*I,rel=1e-10)


def test_static_renewal_diagnostic_does_not_claim_a_redesign(data):
    r=data[2];weak=r['cases']['local_structures_space_chemical']['reference'];strong=r['cases']['favourable_earth_chemical']['strong']
    assert not weak['renewal_affordability']['renewal_reduction_alone_can_close_gap']
    assert strong['renewal_affordability']['renewal_reduction_alone_can_close_gap']
    allowed=strong['operating_cash_balance_USD2041']-strong['resource_transition_provision_USD2041']
    assert strong['renewal_affordability']['fraction_of_external_renewal_affordable']==pytest.approx(allowed/strong['external_renewal_USD2041'])
    assert strong['campaign_cash_payback_years_if_constant_positive_net'] is None
