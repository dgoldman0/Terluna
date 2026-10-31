import json
import math

import pytest

from shared.constants import GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from research.studies.world_economy import pressure_advantages as m


def test_complete_result_reproduction():
    assert m.results()==json.loads((m.HERE/'results/pressure_advantages.json').read_text())


def test_changed_product_is_rejected(tmp_path):
    s,_=m.load_inputs()
    s['products'][next(iter(s['products']))]['sha256']='0'*64
    p=tmp_path/'inputs.json';p.write_text(json.dumps(s))
    with pytest.raises(ValueError,match='Input changed'):
        m.load_inputs(p)


def test_no_pressure_change_has_no_pressure_only_gain():
    assert all(v==1 for v in m.pressure_scalings(1).values())


@pytest.mark.parametrize('bad',[0,-1,True,float('nan'),float('inf')])
def test_invalid_pressure_ratio_rejected(bad):
    with pytest.raises(ValueError):m.pressure_scalings(bad)


def test_air_density_and_flight_comparisons_hold_required_lift():
    r=m.pressure_scalings(1.2)
    assert 1.2*r['wing_speed_at_same_lift']**2==pytest.approx(1)
    assert 1.2*r['volume_for_same_buoyant_mass']==pytest.approx(1)
    assert r['similar_envelope_area_for_same_buoyant_mass']**1.5==pytest.approx(1/1.2)


def test_composition_preserves_oxygen_without_hiding_density_change():
    r=m.results();p=r['gas_partial_pressures']
    assert p['dry_oxygen_ratio']==pytest.approx(1)
    assert 1<p['saturated_water_interface_oxygen_ratio']<1.01
    assert p['locked_trace_CO2_Pa']/p['reference_trace_CO2_Pa']==pytest.approx(1.2)
    assert 1.19<r['environment']['binary_composition_corrected_density_ratio']<1.2


def test_heat_carriage_and_mass_flux_remain_equal():
    r=m.results();c=r['assumptions']['cooling'];row=r['cooling'][0]
    for name in ('one_atm','locked_same_heat_transport'):
        x=row[name]
        assert x['mass_flow_kg_s']*r['environment']['mass_cp_J_kg_K']*c['air_temperature_rise_K']==pytest.approx(c['heat_transport_W'])
        assert x['volume_flow_m3_s']*x['density_kg_m3']==pytest.approx(x['mass_flow_kg_s'])
    a,b=row['one_atm'],row['locked_same_heat_transport']
    assert a['density_kg_m3']*a['speed_m_s']==pytest.approx(b['density_kg_m3']*b['speed_m_s'])


def test_fixed_heat_and_fixed_volume_have_opposite_fan_cost_changes():
    r=m.results()['cooling'][0]
    assert r['same_heat_transport_power_ratio']==pytest.approx(1/1.2**2)
    assert r['same_volume_flow_power_ratio']==pytest.approx(1.2)
    for k in ('one_atm','locked_same_heat_transport','locked_same_volume_flow'):
        x=r[k]
        assert x['electric_power_W']*.65==pytest.approx(x['volume_flow_m3_s']*x['total_pressure_loss_Pa'])


def test_fixed_pressure_loss_reduces_cooling_saving():
    rows=m.results()['cooling']
    assert rows[0]['same_heat_transport_power_ratio']<rows[1]['same_heat_transport_power_ratio']<1/1.2


def test_fan_cost_saving_applies_only_to_the_affected_cost_share():
    r=m.results()['conditional_whole_cost_savings'][1]
    assert r['fan_electricity_cost_share']==.05
    assert r['whole_cost_saving_fraction']==pytest.approx(.05*(1-1/1.2**2))


def test_clamping_preserves_force_with_less_vacuum_and_less_ideal_work():
    r=m.results()['vacuum_assisted'];a,b=r['fixed_force_cases']
    assert a['clamping_pressure_Pa']==b['clamping_pressure_Pa']==80000
    assert b['final_absolute_Pa']-a['final_absolute_Pa']==pytest.approx(20265)
    assert r['fixed_force_pumpdown_ratio']<1
    assert r['fixed_force_ideal_work_ratio']<1
    assert r['locked_clamping_Pa']-r['one_atm_clamping_Pa']==pytest.approx(20265)


def test_same_absolute_vacuum_has_the_original_extra_pumping_cost():
    a=m.clamp_case(101325,101325-10000,.001,.001)
    b=m.clamp_case(121590,121590-10000,.001,.001)
    assert b['constant_speed_no_load_pumpdown_s']>a['constant_speed_no_load_pumpdown_s']
    assert b['ideal_work_J']>a['ideal_work_J']


def test_impossible_fixed_force_is_rejected():
    with pytest.raises(ValueError):m.clamp_case(101325,101325,1,1)


def test_boiling_law_matches_iapws_reference_and_inverts_pressure():
    assert float(m.water.saturation_pressure(373.1243))==pytest.approx(101325,abs=.03)
    assert float(m.water.saturation_pressure(273.16))==pytest.approx(611.657,abs=.001)
    for p in (101325,121590):
        assert float(m.water.saturation_pressure(m.boiling_temperature(p)))==pytest.approx(p,rel=1e-12)
    assert 378.2<m.boiling_temperature(121590)<378.4


def test_liquid_head_uses_lunar_gravity_without_reducing_pressure_requirement():
    e=m.liquid_margin(101325,301,1000,STANDARD_GRAVITY)
    a=m.liquid_margin(101325,301,1000,MOON_SURFACE_GRAVITY)
    b=m.liquid_margin(121590,301,1000,MOON_SURFACE_GRAVITY)
    assert e['pressure_above_vapour_Pa']==a['pressure_above_vapour_Pa']
    assert a['ideal_suction_column_limit_m']/e['ideal_suction_column_limit_m']==pytest.approx(STANDARD_GRAVITY/MOON_SURFACE_GRAVITY)
    assert (b['ideal_suction_column_limit_m']-a['ideal_suction_column_limit_m'])*1000*MOON_SURFACE_GRAVITY==pytest.approx(20265)


def test_hot_open_water_has_a_real_threshold_not_a_universal_quality_gain():
    a=m.liquid_margin(101325,378.15,1000,MOON_SURFACE_GRAVITY)
    b=m.liquid_margin(121590,378.15,1000,MOON_SURFACE_GRAVITY)
    assert not a['unboiled_surface_possible'] and a['ideal_suction_column_limit_m'] is None
    assert b['unboiled_surface_possible'] and b['pressure_above_vapour_Pa']>0


def test_compression_saving_per_mole_is_independent_of_final_pressure():
    rows=m.results()['compression']
    for r in rows:
        assert r['one_atm']['reversible_isothermal_J_mol']-r['locked']['reversible_isothermal_J_mol']==pytest.approx(GAS_CONSTANT*301*math.log(1.2))
    assert rows[0]['isothermal_work_saving_fraction']>rows[1]['isothermal_work_saving_fraction']>rows[2]['isothermal_work_saving_fraction']
    assert m.compressor(101325,101325,301,1.4)['ideal_isentropic_J_mol']==0


def test_isothermal_height_satisfies_spherical_hydrostatic_pressure():
    r=m.results();z=r['ideal_isothermal_height_of_one_atm_m'];e=r['environment']
    p=e['locked_Pa']*math.exp(-e['molar_mass_kg_mol']*MOON_SURFACE_GRAVITY*MOON_RADIUS*z/((MOON_RADIUS+z)*GAS_CONSTANT*301))
    assert p==pytest.approx(e['reference_Pa'])
    assert 9000<z<11000


@pytest.mark.parametrize('call',[
    lambda:m.boiling_temperature(1),
    lambda:m.compressor(2,1,301,1.4),
    lambda:m.compressor(1,2,301,1),
    lambda:m.air_circuit(1,1,1,1,0,1.1),
    lambda:m.pressure_level_height(1,2,301,.029)])
def test_outside_model_boundaries_rejected(call):
    with pytest.raises(ValueError):call()
