import json
import math
import pytest
from research.studies.world_economy import mature_industry as m


def test_product_reproduces():
    assert m.results()==json.loads((m.HERE/'results/mature_industry.json').read_text())


def test_changed_product_rejected(tmp_path):
    s,_=m.load_inputs();s['products'][next(iter(s['products']))]['sha256']='0'*64
    p=tmp_path/'changed.json';p.write_text(json.dumps(s))
    with pytest.raises(ValueError,match='Input changed'):m.load_inputs(p)


def test_pressure_cost_terms_are_distinct_and_fixed_internal_throughput_is_not_scaled():
    r=m.environment_ratios(101325,121590,.1,150,1,1)
    assert r['pressure_load_ratio']==pytest.approx(1.2,rel=1e-6)
    assert r['ideal_pump_work_ratio']==pytest.approx(1.2,rel=1e-5)
    assert 1.01<r['constant_speed_pumpdown_ratio']<1.02
    assert r['ideal_buckling_thickness_ratio']==pytest.approx(math.sqrt(1.2),rel=1e-6)
    assert r['internal_process_load_minimum_inlet_speed_m3_s']==10
    assert r['added_ideal_fixed_throughput_compression_power_W']==pytest.approx(math.log(1.2))


def test_equal_surrounding_pressure_still_requires_pump_work():
    r=m.environment_ratios(101325,101325,.1,150,1,1)
    assert r['ideal_pump_work_ratio']==r['constant_speed_pumpdown_ratio']==1
    assert r['moon_ideal_pump_work_J']>0
    assert r['buoyant_supported_mass_per_volume_ratio']==1


def test_equal_bond_and_rayleigh_size_scalings_preserve_dimensionless_numbers():
    r=m.results()['environment'];g=r['gravity_ratio']
    assert g*r['capillary_length_and_same_Bond_size_ratio']**2==pytest.approx(1)
    assert g*r['same_liquid_Rayleigh_number_length_ratio']**3==pytest.approx(1)
    assert g*r['self_weight_beam_same_absolute_sag_span_ratio']**4==pytest.approx(1)


def test_lift_and_hover_ratios_follow_weight_not_mass_change():
    r=m.results()['environment'];g=r['gravity_ratio'];rho=r['pressure_ratio']
    assert rho*r['same_mass_wing_area_lift_coefficient_speed_ratio']**2==pytest.approx(g)
    assert r['ideal_hover_induced_power_same_mass_disk_area_ratio']**2*rho==pytest.approx(g**3)
    assert r['wind_pressure_same_speed_ratio']==1.2


def test_centrifuge_scale_and_zero_coriolis_for_zero_relative_speed():
    r=m.centrifuge_screen(100,.01,0,.01)
    assert r['radius_for_radial_gradient_m']==10000
    omega=r['rotation_rpm']*2*math.pi/60
    assert omega**2*10000==pytest.approx(m.MOON_SURFACE_GRAVITY)
    assert r['maximum_Coriolis_to_gravity_ratio']==0
    assert r['lunar_central_field_radial_change_over_same_height_fraction']<.0002


def test_reference_costs_depend_on_destination_and_orbital_service():
    d=m.results();rows=[r for r in d['cost_cases'] if r['recipe']=='vacuum_recipe' and r['assumed_gravity_step_saving']==.2 and r['array_site_U_per_good']==5]
    selected={r['destination']:r['comparison'] for r in rows}
    assert selected['moon']['moon_strictly_cheaper_than_all_alternatives']
    assert not selected['array']['moon_strictly_cheaper_than_all_alternatives']
    assert not selected['earth']['moon_strictly_cheaper_than_all_alternatives']


def test_rejects_can_remove_an_otherwise_small_delivered_advantage():
    rows=[r for r in m.results()['yield_sensitivities'] if r['recipe']=='ambient_or_controlled_gas' and r['destination']=='array']
    assert any(r['comparison']['moon_strictly_cheaper_than_all_alternatives'] for r in rows)
    assert any(not r['comparison']['moon_strictly_cheaper_than_all_alternatives'] for r in rows)


def test_mature_renewal_mass_is_distinct_from_awarded_sales():
    x=m.results()['recurring_demand']['aerial']
    assert x['total_annual_renewal_units']==pytest.approx(3.5e12)
    assert x['realized_annual_units']==pytest.approx(35e9)
    assert x['annual_customer_spending_U']==pytest.approx(105e9)


@pytest.mark.parametrize('args',[(1,1,2,1,1,0),(1,1,.1,0,1,0),(1,1,.1,1,0,0)])
def test_bad_pressure_conditions(args):
    with pytest.raises(ValueError):m.environment_ratios(*args)
