"""Independent unit, capacity, conservation and lifecycle checks of the study."""
from copy import deepcopy
import json
import math

import pytest

from shared.constants import STANDARD_ATMOSPHERE
from shared.provenance import constants_changed
from . import buildout


@pytest.fixture(scope='module')
def inputs():
    return buildout.read_inputs()


@pytest.fixture(scope='module')
def reference(inputs):
    return buildout.scenario(*inputs)


def test_population_food_and_shared_production_scope(reference):
    c=reference;f=c['food']
    assert c['city_population']==200000
    assert c['town_population']==50000
    assert c['population']==f['population']==250000
    assert c['geometry']['floor_demand_m2']['ordinary']==14e6
    assert c['resources']['everyday_goods_kg_year']==125e6
    assert f['plant_intake_kcal_person_day']+f['fish_intake_kcal_person_day']+f['microbial_counted_kcal_person_day']==pytest.approx(2500)
    assert f['diet_protein_g_person_day']>94.2


def test_fish_feed_and_real_spare_capacity(reference):
    f=reference['food']
    assert f['live_fish_kg_year']==pytest.approx(6_391_875)
    assert f['feed_kg_year']==pytest.approx(9_587_812.5)
    assert f['feed_crop_kcal_year']>0.1*f['crop_kcal_year']
    assert f['installed_crop_bed_m2']*0.75==pytest.approx(f['active_crop_bed_m2'])
    assert f['one_city_out_crop_capacity_fraction']==pytest.approx(1)
    assert f['stored_food_kg']*3500/250000/2500==pytest.approx(198)


def test_photon_energy_and_crop_area_units(reference):
    f=reference['food']
    # 40 mol/day at 3.5 micromol/J requires 3.174603 kWh/m2/day.
    energy_kWh=f['lighting_W']/1000*24/f['active_crop_bed_m2']
    assert energy_kWh==pytest.approx(3.1746031746031744)
    assert f['active_crop_bed_m2']/250000==pytest.approx(80.2479288)
    assert f['farm_floor_m2']==pytest.approx(8_693_525.62)


def test_enclosure_capacity_and_independent_pressure_loads(inputs,reference):
    g=reference['geometry']
    assert g['available_floor_m2']>=g['required_floor_m2']
    assert g['modules_per_city']==80
    assert g['equivalent_cave_length_m_per_city']==11200
    assert g['cross_section_containment_bound']<1
    assert g['module_end_force_N']==pytest.approx(STANDARD_ATMOSPHERE*math.pi*30**2)
    assert g['overburden_ratio_to_atmosphere']>4
    s,src=inputs;weaker=deepcopy(s);weaker['geometry']['effective_allowable_Pa']=400e6
    other=buildout.scenario(weaker,src)['geometry']
    assert other['pressure_reinforcement_kg']==pytest.approx(1.5*g['pressure_reinforcement_kg'])
    assert other['pressure_liner_kg']==g['pressure_liner_kg']
    bad=deepcopy(s);bad['geometry']['cave_width_m']=180
    with pytest.raises(ValueError,match='do not fit'):
        buildout.geometry(bad,reference['food'])


def test_one_utility_budget_and_explicit_town_credits(inputs,reference):
    s,src=inputs;c=reference['construction'];m=reference['maintenance']
    old=src['towns']['cases']['500']['costs']['reference']
    assert sum(old['direct_U'].values())-sum(c['town_direct_credits_U'].values())==pytest.approx(c['retained_town_direct_U'])
    assert c['retained_town_direct_U']+sum(c['additions_direct_U'].values())==pytest.approx(c['direct_U'])
    assert c['initial_U']==pytest.approx(c['direct_U']+c['contingency_U']+c['transition_reserve_U'])
    assert not {'town_firm_power','town_emergency_storage','town_radiators'} & m['assets_U'].keys()
    assert 'town_utility_operation' not in m['variable_operations_U_year']
    assert c['town_direct_credits_U']['embedded_production_capital']==6.25e9


def test_food_water_and_nutrient_imports_are_not_free(reference):
    f=reference['food'];r=reference['resources']
    assert f['crop_water_makeup_kg_year']==pytest.approx(146_552_779.971)
    assert f['ras_water_makeup_kg_year']==pytest.approx(24_289_125)
    assert f['nutrient_feed_delivered_kg_year']>sum(f['nutrient_element_makeup_kg_year'].values())
    assert f['initial_working_feed_kg']>90e6
    assert r['operating_import_kg_year']['water']>f['crop_water_makeup_kg_year']+f['ras_water_makeup_kg_year']
    assert r['total_annual_import_kg']==pytest.approx(sum(r['operating_import_kg_year'].values())+sum(r['renewal_import_kg_year'].values()))
    assert r['local_goods_makeup_kg_year']+r['operating_import_kg_year']['everyday_carbon_feed']==pytest.approx(r['virgin_goods_kg_year'])


def test_feed_mass_and_capacity_account_includes_old_foundations(inputs,reference):
    s,src=inputs;q=src['towns']['cases']['500']['quantities'];r=reference['resources'];d=reference['schedule']['production_demand']
    assert d['local_bulk']+d['local_metal']+q['shielding_kg']==pytest.approx(r['local_finished_kg'])
    assert d['handling']>r['raw_feed_kg']+r['additional_grading_kg']
    assert d['pressure_components']+d['other_imports']==pytest.approx(r['special_imports_kg'])
    assert reference['utilities']['construction_firm_W']>=s['production']['local_bulk_kg_year']*s['resources']['bulk_MJ_kg']/3.6/8766*1000+s['production']['local_metal_kg_year']*s['resources']['metal_kWh_kg']/8766*1000+s['production']['handling_kg_year']*s['resources']['handling_MJ_kg']/3.6/8766*1000


def test_expected_renewal_is_separate_from_horizon_cash(reference):
    m=reference['maintenance'];initial=reference['construction']['initial_U']
    assert m['total_U_year']==pytest.approx(m['recurring_U_year']+sum(m['renewal_U_year'].values()))
    assert buildout.replacement_count(50,25)==1
    assert buildout.replacement_count(51,25)==2
    assert buildout.replacement_count(25,25)==0
    with pytest.raises(ValueError):buildout.replacement_count(0,25)
    for year,h in m['horizons'].items():
        assert h['initial_plus_vacuum_cash_U']==pytest.approx(initial+float(year)*m['recurring_U_year']+h['actual_replacements_U'])


def test_lower_crop_yield_and_weaker_recovery_have_material_effects(inputs,reference):
    s,src=inputs;changed=deepcopy(s);changed['food']['crop_productivity_kcal_m2_day']/=2
    q=buildout.scenario(changed,src)
    assert q['food']['active_crop_bed_m2']==2*reference['food']['active_crop_bed_m2']
    assert q['geometry']['equivalent_cave_length_m_per_city']>reference['geometry']['equivalent_cave_length_m_per_city']
    assert q['construction']['initial_U']>reference['construction']['initial_U']
    changed=deepcopy(s);changed['food']['crop_water_recovery']=.98;q=buildout.scenario(changed,src)
    assert q['food']['crop_water_makeup_kg_year']==pytest.approx(4*reference['food']['crop_water_makeup_kg_year'])
    assert q['resources']['total_annual_import_kg']>reference['resources']['total_annual_import_kg']


def test_changed_input_is_rejected(inputs):
    s,_=inputs;wrong=deepcopy(s);wrong['inputs']['towns']['sha256']='0'*64
    with pytest.raises(ValueError,match='Input hash mismatch'):
        buildout.read_inputs(wrong)
    wrong=deepcopy(s);wrong['inputs']['food']['schema']='wrong'
    with pytest.raises(ValueError,match='Input schema mismatch'):
        buildout.read_inputs(wrong)


def test_saved_product_and_constants_reproduce():
    saved=json.loads((buildout.HERE/'results/buildout.json').read_text())
    assert not constants_changed(saved['producer']['constants'])
    assert saved==buildout.run()
