"""Conservation, temporal allocation and reproducible integrated products."""
from dataclasses import replace
import hashlib
import json
import math

import pytest

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from research.studies.floater_viability import run as model


def air():
    data=json.loads((model.ROOT/'research/studies/sky_ships/results/sky_ships.json').read_text())
    return next(x for x in data['air'] if x['height_km']==10)


def case(**kw):
    return model.coupled_case(20.,air(),replace(model.Traits(),**kw))


def test_supported_mass_includes_all_stocks_and_exact_neutral_trim():
    c=case();m=c['mass'];a=math.pi*20**2
    explicit=m['structure_kg']+a*(m['living_and_community_wet_kg_m2']+m['free_water_kg_m2']+m['stored_reserve_wet_kg_m2'])
    assert explicit==pytest.approx(m['actual_supported_mass_kg'])
    assert m['total_mass_closure_kg']==pytest.approx(0,abs=1e-8)
    assert c['hydrogen']['inventory_kg']==m['actual_hydrogen_mass_kg']


def test_annual_carbon_is_allocated_once_and_budget_stage_is_explicit():
    c=case()['carbon']
    sinks=sum(c[k] for k in ('host_maintenance_c_kg_m2_year','night_h2_gross_substrate_c_kg_m2_year',
        'consumer_food_c_kg_m2_year','turnover_construction_c_kg_m2_year','remaining_assimilate_c_kg_m2_year'))
    assert c['gpp_after_h2_c_kg_m2_year']==pytest.approx(sinks)
    assert 'Post-maintenance' in c['stage']


def test_continuous_productive_light_removes_dark_feed_and_reserve():
    c=case(annual_dark_fraction=0,longest_dark_days=0)
    assert c['carbon']['night_h2_gross_substrate_c_kg_m2_year']==0
    assert c['storage']['minimum_reserve_dry_kg_m2']==0
    assert c['solar_only_night_alternative']['required_water_ballast_release_kg_m2']==0
    assert c['carbon']['host_maintenance_c_kg_m2_year']==pytest.approx(.003*.4*JULIAN_YEAR_DAYS)


def test_longest_dark_interval_changes_storage_without_recharging_annual_food():
    a=case(longest_dark_days=2);b=case(longest_dark_days=10)
    assert b['storage']['minimum_reserve_dry_kg_m2']==pytest.approx(5*a['storage']['minimum_reserve_dry_kg_m2'])
    assert a['carbon']==b['carbon']


def test_photons_for_dark_leakage_are_paid_next_day_with_water_ballast():
    f=case();b=case(dark_hydrogen_strategy='water_ballast')
    assert b['carbon']['night_h2_gross_substrate_c_kg_m2_year']==0
    assert b['hydrogen']['required_hydrogen_area_time_fraction']==pytest.approx(2*f['hydrogen']['required_hydrogen_area_time_fraction'])
    assert b['storage']['minimum_reserve_dry_kg_m2']<f['storage']['minimum_reserve_dry_kg_m2']
    assert b['water']['ballast_replenishment_kg_m2_year']>0


def test_budget_rejects_high_permeability_without_hiding_shortfall():
    c=case(permeability_barrer=100)
    assert c['carbon']['remaining_assimilate_c_kg_m2_year']<0
    assert not c['gates']['solar_h2_allocation']
    assert c['hydrogen']['unmet_hydrogen_w_m2']>0


def test_parent_funded_reproduction_closes_carbon_and_exclusive_area():
    c=case();r=c['reproduction']
    assert r['allocations_close']
    assert 0 < r['total_gas_area_fraction'] < 1
    assert r['carbon_residual_kg_m2']==pytest.approx(0,abs=1e-8)
    assert r['parent_funded_years']>r['full_surface_gas_fill_years']


def test_nonpositive_surplus_never_creates_a_reproduction_time():
    r=model.parent_funded_reproduction(-1,2,1,1e8,2,.1)
    assert r['parent_funded_years'] is None
    assert not r['allocations_close']


def test_simultaneous_water_and_gas_flows_close_without_vent_tax():
    trim=case()['trim_comparisons']
    assert trim['balanced_daily_water_and_gas_flow']['final_upward_residual_kg']==pytest.approx(0)
    assert trim['one_kg_water_loss_and_unreplaced_daily_gas_loss']['final_upward_residual_kg']>0
    assert not trim['retained_gas_one_kg_water_loss']['extra_pressure_allowance_sufficient']
    assert not trim['fixed_hull_air_ballonet_one_kg_water_loss']['extra_pressure_allowance_sufficient']


def test_contracting_gas_solution_limits_and_fixed_area_bound():
    a=model.contracting_pure_hydrogen_bladder(10,0,1e5,280,1e6)
    assert a['inventory_fraction_remaining']==1
    a=model.contracting_pure_hydrogen_bladder(10,1,1e5,280,1e9)
    assert a['inventory_fraction_remaining']==0
    flux=1e-6;t=10*JULIAN_DAY;r=10.;p=1e5;temp=280
    a=model.contracting_pure_hydrogen_bladder(r,flux,p,temp,t)
    linear=3*flux*model.GAS_CONSTANT*temp*t/(p*r)
    assert 0 < 1-a['inventory_fraction_remaining'] <= linear


def test_invalid_zero_dark_regime_is_rejected():
    with pytest.raises(ValueError):case(annual_dark_fraction=0,longest_dark_days=1)


def test_saved_product_binds_all_producers_inputs_and_named_constants():
    data=json.loads(model.OUT.read_text())
    assert data['schema']==model.SCHEMA
    for category in ('files','inputs'):
        for path,digest in data['producer'][category].items():
            assert hashlib.sha256((model.ROOT/path).read_bytes()).hexdigest()==digest,path
    files=[model.ROOT/p for p in data['producer']['files']]
    assert data['producer']['constants']==constants_used(files)


def test_saved_product_regenerates_exactly(tmp_path,monkeypatch):
    before=model.OUT.read_bytes()
    target=tmp_path/'floater_viability.json'
    monkeypatch.setattr(model,'OUT',target)
    model.run()
    assert target.read_bytes()==before
