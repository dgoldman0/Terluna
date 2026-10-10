"""Dimensional, conservation and independent limiting checks for water and heat."""
import json
import math

import pytest

from research.studies.aerophytes import environment as e


def test_saturation_values_and_physical_range():
    assert e.saturation_pa(0)==pytest.approx(611.21)
    assert e.saturation_pa(20)==pytest.approx(2338.3,rel=2e-4)
    assert e.saturation_pa(30)>e.saturation_pa(20)>e.saturation_pa(10)
    with pytest.raises(ValueError): e.saturation_pa(100)


def test_diffusion_includes_boundary_and_closure_is_not_no_loss():
    assert e.series_conductance(.1,.5)==pytest.approx(1/12)
    case=e.transpiration(20,20,.5,100000,.1)
    # 0.5*2338.3Pa deficit /100kPa, multiplied by molar conductance and water mass.
    expected=(1/12)*(.5*2338.3/100000)*.01801528*86400
    assert case['water_kg_m2_day']==pytest.approx(expected,rel=2e-4)
    assert case['latent_w_m2']==pytest.approx(case['water_kg_m2_day']/86400*e.LATENT_J_KG)
    assert e.transpiration(20,20,.5,100000,.0016)['water_kg_m2_day']>0
    assert e.transpiration(20,20,1,100000,.1)['water_kg_m2_day']==0
    assert e.transpiration(20,5,.9,100000,.1)['condensation_possible']
    assert e.transpiration(20,5,.9,100000,.1)['water_kg_m2_day']==0


def test_stomatal_carbon_water_relation_independent_molar_budget():
    out=e.water_for_assimilation(120.11,20,20,.5,100000,co2_ppm=400)
    # 120.11g C =10mol; DeltaCO2=12Pa. Need1.6*DeltaH2O/12 mol water/mol C.
    expected=10*(1.6*.5*e.saturation_pa(20)/12)*e.WATER_MOLAR_KG
    assert out['water_kg_m2_year']==pytest.approx(expected)
    enriched=e.water_for_assimilation(120.11,20,20,.5,100000,co2_ppm=2000)
    assert enriched['water_kg_m2_year']==pytest.approx(out['water_kg_m2_year']/5)
    twice=e.water_for_assimilation(240.22,20,20,.5,100000)
    assert twice['water_kg_m2_year']==pytest.approx(2*out['water_kg_m2_year'])


def test_storage_spends_only_mobile_fraction_and_counts_leaf_area():
    row=e.water_cycle(2,0,3,rainless_days=10,water_store_kg_m2=100,usable_share=.2)
    assert row['demand_kg_m2']==30
    assert row['usable_store_kg_m2']==20
    assert row['autonomy_days']==pytest.approx(20/3)
    assert row['required_total_store_kg_m2']==150
    assert e.water_cycle(0,0)['autonomy_days'] is None


def test_fog_requires_relative_flow_and_declared_droplets():
    case=e.fog_collection(.1,3,.3,collector_area_ratio=2,fog_duty=.1)
    assert case['in_fog_kg_m2_hour']==pytest.approx(.648)
    assert case['mean_kg_m2_day']==pytest.approx(.648*24*.1)
    assert e.fog_collection(.1,0,.3)['mean_kg_m2_day']==0
    assert e.fog_collection(0,3,.3)['mean_kg_m2_day']==0


def test_dew_requires_both_cooling_and_subdewpoint_temperature():
    no_dew=e.dew_collection(50,20,20,.6)
    assert not no_dew['below_dewpoint'] and no_dew['energy_ceiling_kg_m2_day']==0
    row=e.dew_collection(50,5,20,.6)
    assert row['energy_ceiling_kg_m2_day']*e.LATENT_J_KG==pytest.approx(50*86400)


def test_rain_and_drainage_account_conserves_arrival():
    assert e.rain_load(60,.25)['accumulated_kg_m2']==15
    row=e.rain_load(60,.25,.5,10)
    assert row['accumulated_kg_m2']==5
    assert e.rain_load(60,1,.5,30)['accumulated_kg_m2']==0


def test_leaf_energy_closes_and_transpiration_cools():
    zero=e.leaf_equilibrium(20,20,0,1,100000,0)
    assert zero['leaf_c']==pytest.approx(20)
    closed=e.leaf_equilibrium(20,20,300,.6,100000,.0016)
    open_=e.leaf_equilibrium(20,20,300,.6,100000,.1)
    assert open_['leaf_c']<closed['leaf_c']
    assert open_['water_kg_m2_day']>closed['water_kg_m2_day']
    for row in (closed,open_):
        assert row['radiative_loss_w_m2']+row['sensible_loss_w_m2']+row['latent_w_m2']==pytest.approx(300)
        assert abs(row['energy_residual_w_m2'])<1e-8


def test_heat_storage_does_not_imply_night_survival():
    row=e.thermal_store(9,10,1)
    assert row['autonomy_days']==pytest.approx(9*4180*10/86400)
    assert not row['covers_interval']
    assert e.thermal_store(9,10,0)['covers_interval']


def test_high_air_heating_inverts_leaf_energy_balance():
    required=e.leaf_heat_requirement(-13.9,-13.9,5,.6,55000,.0016,10)
    leaf=e.leaf_equilibrium(-13.9,-13.9,required['net_heat_required_w_m2'],.6,55000,.0016,10)
    assert leaf['leaf_c']==pytest.approx(5)
    assert required['sensible_loss_w_m2']==pytest.approx(189)


def test_evaluation_reads_actual_air_and_never_assigns_cloud_residence():
    result=e.evaluate()
    assert json.loads(json.dumps(result,allow_nan=False))==result
    air={v['height_km']:v for v in result['air']}
    assert air[10]['temperature_c']==14.2
    assert air[40]['temperature_c']<0
    assert result['mean_freezing_height_km']==25.7
    assert result['co2']['inherited_canopy_partial_pa']['moon']==48.636
    assert all('relative_humidity' not in row for row in result['air'])
    assert {row['co2_ppm'] for row in result['gross_to_net_leaf_water']}=={400,1000,2000}
