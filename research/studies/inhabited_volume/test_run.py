"""Conservation, physical limits and input provenance of the habitation screen."""
import copy
import json
import math
import pytest
from shared.provenance import constants_changed
from research.studies.inhabited_volume import run as m

@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def test_stored_product_and_inputs(product):
    assert product['schema'] == m.SCHEMA
    assert product == m.results()
    for p,h in {**product['producer']['files'],**product['producer']['inputs']}.items():
        assert m.digest(m.ROOT/p)==h,p
    assert not constants_changed(product['producer']['constants'])


def test_rejects_wrong_schema(tmp_path):
    p=tmp_path/'wrong.json'; p.write_text('{"schema":"wrong"}')
    with pytest.raises(ValueError): m.load_product(p,m.SCHEMA)


def test_population_conserved_and_bad_sum_rejected(product):
    s=m.load_product(m.ROOT/'shared/scenarios/population.json','terluna.scenario.population/1')
    m.validate_cases(s)
    bad=copy.deepcopy(s);bad['cases']['total25']['lunar_aerial']+=1
    with pytest.raises(ValueError): m.validate_cases(bad)
    for c in product['cases'].values():
        p=c['population_billion']
        assert c['lunar_population_billion']==p['lunar_surface']+p['lunar_aerial']
        assert c['wild_aerial_food_credited_people']==0


def test_archimedes_and_mass_allowance_close(product):
    for c in product['lift_cases']:
        air=next(a for a in product['air'] if a['height_km']==c['height_km'])
        displaced=air['density_kg_m3']*c['volume_m3']/1e3
        assert math.isclose(c['gross_supported_t']+c['gas_mass_t'],displaced)
        assert math.isclose(c['fleet_lifting_species_gt']+c['fleet_gas_entrained_air_gt'],c['fleet_gas_gt'])
        assert math.isclose(c['district_payload_t']+c['other_and_reserve_t'],c['gross_supported_t'])
        assert math.isclose(c['residents']*c['district_t_per_person'],c['district_payload_t'])
        # At nominal load,15% air ballonet permits unused lift capacity without overbuoyancy.
        non_gas_loaded=.85*c['gross_supported_t']
        assert math.isclose(non_gas_loaded+c['neutral_trim_gas_t']+c['neutral_trim_ballonet_air_t'],displaced)
        assert math.isclose(c['fleet_equivalent_count']*c['residents'],c['people'])
        assert math.isclose(c['fleet_other_and_reserve_gt'],c['fleet_envelope_control_hardware_gt']+c['fleet_unallocated_reserve_gt'])
        assert c['fleet_count_rounded_up']*c['residents']>=c['people']
        assert (c['fleet_count_rounded_up']-1)*c['residents']<c['people']


def test_share_and_district_budget(product):
    assert math.isclose(sum(m.LIFT_SHARES.values()),1.)
    for case,total in zip(m.DISTRICTS.values(),(20,50,100)):
        assert sum(case.values())==total
    central=product['cases']['total25']['aerial_reference']
    # An independent mass-per-area result: required A = M / (f * rho_lift * 4r/3).
    area=(5e9*50e3)/(.5*1.098*(4*500/3))/1e6
    assert math.isclose(central['sum_projected_area_km2'],area)
    assert math.isclose(central['fleet_district_gt'],250.)
    assert math.isclose(central['fleet_unallocated_reserve_gt'],75.)


def test_radius_and_mass_scaling(product):
    small,medium,large=product['radius_sensitivity']
    assert math.isclose(small['fleet_equivalent_count']/medium['fleet_equivalent_count'],8.)
    assert math.isclose(small['sum_projected_area_km2']/medium['sum_projected_area_km2'],2.)
    assert math.isclose(medium['fleet_gas_gt'],large['fleet_gas_gt'])
    a=product['air'][1]
    one=m.envelope(a,1e9); two=m.envelope(a,2e9)
    assert math.isclose(two['fleet_district_gt']/one['fleet_district_gt'],2.)
    assert m.envelope(a,0)['fleet_count_rounded_up']==0


def test_shear_reference_and_drag_power(product):
    a=product['air'][1]
    zero=m.drag(a,500,0);five=m.drag(a,500,5);ten=m.drag(a,500,10)
    assert zero['input_power_w']==0
    assert math.isclose(ten['drag_n']/five['drag_n'],4.)
    assert math.isclose(ten['input_power_w']/five['input_power_w'],8.)
    assert math.isclose(ten['input_power_w']*.7,ten['drag_n']*10.)
    assert ten['input_power_w']>=ten['thrust_power_w']


def test_gravity_integral_and_water_mass_balance():
    row=m.water(5e9,150,.9,20)
    # Midpoint quadrature of g(h), independent of the closed form.
    n=2000;dh=20000/n
    integrated=sum(m.MOON_GM/(m.MOON_RADIUS+(i+.5)*dh)**2*dh for i in range(n))
    assert math.isclose(row['ideal_lift_j_kg'],integrated,rel_tol=1e-10)
    assert math.isclose(row['makeup_km3_year'],.1*row['gross_km3_year'])
    assert row['makeup_lift_tw']*1e12>=row['makeup_kg_s']*integrated
    assert m.water(5e9,150,1.,20)['makeup_lift_tw']==0
    assert m.water(5e9,150,.9,0)['makeup_lift_tw']==0
    assert m.water(0,150,.9,20)['makeup_kg_s']==0
    assert math.isclose(m.water(5e9,150,.99,20)['makeup_kg_s'],row['makeup_kg_s']/10)


def test_altitude_reduces_mass_lift(product):
    h=[m.envelope(a,1e9) for a in product['air']]
    assert all(a['residents']>b['residents'] for a,b in zip(h,h[1:]))
    for a in product['air']:
        assert m.envelope(a)['residents']>m.envelope(a,gas='helium')['residents']


@pytest.mark.parametrize('args',[(1,150,1.01,10),(1,150,-.1,10),(1,150,.9,-1),(-1,150,.9,10)])
def test_nonphysical_water_rejected(args):
    with pytest.raises(ValueError): m.water(*args)


def test_zero_aerial_share_conserves_people_and_floor(product):
    for c in product['cases'].values():
        base=c['aerial_share_sensitivity'][0]
        assert base['aerial_share']==0 and base['reference']['fleet_district_gt']==0
        for row in c['aerial_share_sensitivity']:
            assert math.isclose(row['reference']['people'],c['lunar_population_billion']*1e9*row['aerial_share'])
            assert math.isclose(row['aerial_floor_km2'],c['lunar_total_floor_km2']*row['aerial_share'])


def test_roaming_policy_and_duty_account(product):
    assert product['assumptions']['habitation_direction']['source']=='Author, 2026-10-09'
    for row in product['roaming_control_sensitivity']:
        continuous=next(r for r in product['held_city_drag'] if r['relative_air_speed_m_s']==row['relative_air_speed_m_s'])
        assert math.isclose(row['fleet_mean_input_tw'],continuous['fleet_input_tw']*row['powered_duty_fraction'])
        assert math.isclose(row['one_event_energy_mwh']['24'],row['one_event_energy_mwh']['1']*24)
        assert row['fleet_mean_input_tw']<continuous['fleet_input_tw']
