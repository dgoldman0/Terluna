"""Independent mass, pressure and cycle conservation checks."""
import math

import pytest

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from research.studies.floater_viability.trim_cycle import (
    ballonet_air_trim, evaluate, hydrogen_refill_budget,
    net_lift_per_hydrogen_kg, retained_hydrogen_compression,
    trim_ledger, vent_for_neutrality, water_inventory_band,
)


RHO_AIR = 1.204
RHO_MIX = .106
PRESSURE_PA = .988*101325


def test_vented_mixture_and_added_ambient_air_close_actual_mass():
    b = net_lift_per_hydrogen_kg(RHO_AIR,RHO_MIX)
    step = trim_ledger(b,evaporation_kg=1)
    control = vent_for_neutrality(step['final_upward_residual_kg'],b)
    h2_vent = control['additional_hydrogen_vent_kg']
    h2_species_density = RHO_MIX-.02*RHO_AIR
    removed_volume = h2_vent/h2_species_density
    # Fixed outer hull: lost payload + lost mixture + admitted air = zero.
    total_mass_change = -1-RHO_MIX*removed_volume+RHO_AIR*removed_volume
    assert total_mass_change == pytest.approx(0)
    assert h2_vent == pytest.approx(.074608378870674)
    assert control['vent_net_lift_loss_kg'] == pytest.approx(1)


def test_opposing_water_flows_and_leak_are_netted_before_control():
    b = 10.
    step = trim_ledger(b,humidity_water_uptake_kg=.3,rain_water_kg=.4,
                      evaporation_kg=1.,hydrogen_leak_kg=.02)
    assert step['water_mass_change_kg'] == pytest.approx(-.3)
    assert step['net_gas_lift_change_kg'] == pytest.approx(-.2)
    assert step['final_upward_residual_kg'] == pytest.approx(.1)
    control = vent_for_neutrality(step['final_upward_residual_kg'],b)
    assert control['additional_hydrogen_vent_kg'] == pytest.approx(.01)
    assert control['additional_hydrogen_refill_kg'] == 0


def test_closed_payload_cycle_replaces_leak_plus_actual_vent():
    b = 10.
    out = trim_ledger(b,evaporation_kg=5,hydrogen_leak_kg=.1)
    out_control = vent_for_neutrality(out['final_upward_residual_kg'],b)
    back = trim_ledger(b,rain_water_kg=5,hydrogen_leak_kg=.1)
    back_control = vent_for_neutrality(back['final_upward_residual_kg'],b)
    assert out_control['additional_hydrogen_vent_kg'] == pytest.approx(.4)
    assert back_control['additional_hydrogen_refill_kg'] == pytest.approx(.6)
    assert back_control['additional_hydrogen_refill_kg'] == pytest.approx(
        .2+out_control['additional_hydrogen_vent_kg'])
    assert out['ledger_residual_kg'] == pytest.approx(0)
    assert back['ledger_residual_kg'] == pytest.approx(0)


def test_night_water_dump_only_offsets_the_stated_gas_loss():
    b = 13.4
    isolated = trim_ledger(b,hydrogen_leak_kg=.1,water_dump_kg=1.34)
    assert isolated['final_upward_residual_kg'] == pytest.approx(0)
    with_evaporation = trim_ledger(b,hydrogen_leak_kg=.1,water_dump_kg=1.34,
                                  evaporation_kg=.2)
    assert with_evaporation['final_upward_residual_kg'] == pytest.approx(.2)
    # Starch respiration: 162 g dry lost, 90 g product water retained.
    retained_water = trim_ledger(b,dry_mass_change_kg=-.162,metabolic_water_change_kg=.09)
    escaped_water = trim_ledger(b,dry_mass_change_kg=-.162)
    assert retained_water['supported_mass_change_kg'] == pytest.approx(-.072)
    assert escaped_water['supported_mass_change_kg'] == pytest.approx(-.162)


def test_water_capacity_is_not_actual_dumpable_inventory():
    case = water_inventory_band(5,1,20,outflow_kg_day=1,interval_days=7)
    assert case['unfilled_water_capacity_kg'] == 15
    assert case['actual_usable_water_for_dump_kg'] == 4
    assert case['days_to_minimum'] == 4
    assert case['required_extra_water_kg'] == 3
    assert not case['inventory_band_feasible']
    fill = water_inventory_band(5,1,9,inflow_kg_day=3,interval_days=2)
    assert fill['required_extra_drainage_kg'] == 2
    assert fill['unconstrained_final_free_water_kg'] == 11


def test_refill_routes_have_consistent_units_and_molar_feed():
    refill = hydrogen_refill_budget(.002016*2.1,1,projected_area_m2=2)
    assert refill['fermentation_glucose_feed_kg'] == pytest.approx(.180156)
    assert refill['fermentation_feed_carbon_kg'] == pytest.approx(.180156*.4)
    assert refill['mean_hydrogen_chemical_w_m2'] == pytest.approx(
        .002016*2.1*120e6/(2*JULIAN_DAY))
    assert refill['repeated_cycle_fermentation_carbon_kg_m2_year'] == pytest.approx(
        .180156*.4*JULIAN_YEAR_DAYS/2)


def test_retained_gas_compression_conserves_gas_not_constant_density():
    initial_volume = 20/(RHO_AIR-RHO_MIX)
    case = retained_hydrogen_compression(1,initial_volume,RHO_AIR,
                                         PRESSURE_PA,PRESSURE_PA,400)
    vf = case['final_gas_volume_m3']
    pf = case['final_gas_pressure_pa']
    assert PRESSURE_PA*initial_volume == pytest.approx(pf*vf)
    gas_mass = RHO_MIX*initial_volume
    assert RHO_AIR*initial_volume-gas_mass == pytest.approx(20)
    assert RHO_AIR*vf-gas_mass == pytest.approx(19)
    assert case['added_pressure_pa'] > 4700
    assert not case['extra_pressure_allowance_sufficient']
    assert case['minimum_work_above_ambient_j'] > 0
    assert case['minimum_work_above_ambient_j'] < case['compression_work_on_gas_j']
    assert case['required_design_pressure_pa'] == pytest.approx(400+case['added_pressure_pa'])


def test_equal_pressure_ballonet_adds_real_air_mass_and_raises_pressure():
    outer,gas = 30.,20.
    case = ballonet_air_trim(1,outer,gas,RHO_AIR,PRESSURE_PA,PRESSURE_PA,400)
    pf = case['final_internal_pressure_pa']
    vf = case['final_gas_volume_m3']
    rho_air_final = RHO_AIR*pf/PRESSURE_PA
    old_air_mass = RHO_AIR*(outer-gas)
    new_air_mass = rho_air_final*(outer-vf)
    assert new_air_mass-old_air_mass == pytest.approx(1)
    assert pf*vf == pytest.approx(PRESSURE_PA*gas)
    assert case['minimum_isothermal_pumping_work_j'] > 0
    assert not case['extra_pressure_allowance_sufficient']
    # Independent numerical quadrature of R_air*T*ln(p/p_ambient) per kg admitted.
    n = 1000
    numerical = sum((PRESSURE_PA/RHO_AIR)*math.log(
        1+(i+.5)/n/(RHO_AIR*outer)) for i in range(n))/n
    assert case['minimum_isothermal_pumping_work_j'] == pytest.approx(numerical,rel=1e-6)


def test_zero_mass_change_requires_neither_vent_nor_compression():
    control = vent_for_neutrality(0,13.4,available_hydrogen_kg=0)
    assert control['additional_hydrogen_vent_kg'] == 0
    assert control['additional_hydrogen_refill_kg'] == 0
    case = retained_hydrogen_compression(0,20,RHO_AIR,PRESSURE_PA,PRESSURE_PA,400)
    assert case['added_pressure_pa'] == 0
    assert case['minimum_work_above_ambient_j'] == 0
    assert case['extra_pressure_allowance_sufficient']


def test_impossible_controls_are_not_certified_and_small_evaluation_runs():
    assert not vent_for_neutrality(20,10,available_hydrogen_kg=1)['vent_inventory_feasible']
    impossible = retained_hydrogen_compression(30,20,RHO_AIR,PRESSURE_PA,PRESSURE_PA,400)
    assert not impossible['geometry_feasible']
    assert impossible['final_gas_pressure_pa'] is None
    with pytest.raises(ValueError):
        trim_ledger(10,evaporation_kg=-1)
    with pytest.raises(ValueError):
        water_inventory_band(5,6,10)
    with pytest.raises(ValueError):
        ballonet_air_trim(1,10,20,RHO_AIR,PRESSURE_PA,PRESSURE_PA,400)
    assert len(evaluate(RHO_AIR,RHO_MIX,PRESSURE_PA)['daily_evaporation_then_refill']) == 4
