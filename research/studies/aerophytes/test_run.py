"""Conservation, temporal allocation and reproducible integrated products."""
from dataclasses import replace
import hashlib
import json
import math

import pytest

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from research.studies.aerophytes import run as model


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
    target=tmp_path/'aerophytes.json'
    monkeypatch.setattr(model,'OUT',target)
    model.run()
    assert target.read_bytes()==before


# Schema-1 reference bodies (product of 2026-10-09 before the giants round):
# diameter -> supported mass per m², lifting-mix share, remaining assimilate.
FIRST_ROUND = {40.: (18.70278700256628, 0.6387563866996678, 1.0292073468154646),
               100.: (24.589244728636576, 0.33591864383383296, 0.7346490022029066)}


def product():
    return json.loads(model.OUT.read_text())


# Surviving ranges at 10 km: family -> (smallest, largest, gate that fails just above the range).
EXPECTED_RANGES_10_KM={'round_collagen':(100.,200.,'mass'),'round':(50.,700.,'carbon'),
                       'round_strong':(50.,1000.,'carbon'),'round_strong_gpp5':(50.,3000.,'carbon'),
                       'colony_module':(60.,150.,'carbon')}


def test_schema_binds_every_component_and_its_inputs():
    data=product()
    assert data['schema']=='terluna.research.aerophytes/1'
    files={p.rsplit('/',1)[-1] for p in data['producer']['files']}
    for name in ('run.py','mechanics.py','envelope.py','biology.py','environment.py','navigation.py',
                 'photoperiod.py','trim_cycle.py','structure.py','growth.py','water.py','storms.py',
                 'sailing.py','gas_biology.py','resources.py'):
        assert name in files
    inputs=set(data['producer']['inputs'])
    study='research/studies/aerophytes/'
    # This product owns its producer components' sources. The later sail-biology
    # requirement study has an independent product and regeneration/hash test.
    for component in files:
        source=component.removesuffix('.py')+'_sources.json'
        if (model.ROOT/study/source).exists():
            assert study+source in inputs
    for path in ('climate/results/gcm/global_winds_A28_dim5_moon.json','climate/gcm/products/climatology_A28_dim5_moon.npz',
                 'climate/results/crm/ring_ring_equator.json','climate/results/crm/ring_ring_80s_lsw.json',
                 'research/studies/megaforest_wind/results/reference_cases.csv','biosphere/canopy/results/plant.json'):
        assert path in inputs


def test_first_round_bodies_keep_their_numbers():
    rows={c['diameter_m']:c for c in product()['cases'][:4]}
    for d,(supported,share,surplus) in FIRST_ROUND.items():
        m=rows[d]['mass']
        assert m['supported_kg_m2']==pytest.approx(supported,rel=1e-12)
        assert m['required_lifting_mix_volume_fraction']==pytest.approx(share,rel=1e-12)
        assert rows[d]['carbon']['remaining_assimilate_c_kg_m2_year']==pytest.approx(surplus,rel=1e-12)
        assert m['tendon_dry_kg_m2']==m['fittings_kg']==m['seams_kg']==m['network_kg']==m['skin_water_kg']==0


def test_fittings_seams_and_ties_enter_the_mass_once():
    lobed=dict(architecture='lobed',tendon_allowable_mpa=60.)
    a=case(**lobed);b=case(**lobed,fitting_share=.1,seam_mass_share=.05,network_kg_m2=.3)
    tendon=a['skin']['tendon_mass_kg']
    assert tendon>0
    assert b['mass']['fittings_kg']==pytest.approx(.1*tendon)
    assert b['mass']['seams_kg']==pytest.approx(.05*(a['skin']['total_structure_mass_kg']-tendon))
    assert b['mass']['network_kg']==pytest.approx(.3*math.pi*20**2)
    added=b['mass']['fittings_kg']+b['mass']['seams_kg']+b['mass']['network_kg']
    assert b['mass']['structure_kg']-a['mass']['structure_kg']==pytest.approx(added)
    assert b['mass']['total_mass_closure_kg']==pytest.approx(0,abs=1e-8)


def test_tendons_renew_at_their_own_rate():
    lobed=dict(architecture='lobed',inert_turnover_per_year=.1)
    unset=case(**lobed);slow=case(**lobed,tendon_turnover_per_year=.01)
    assert unset['carbon']==case(**lobed,tendon_turnover_per_year=.1)['carbon']
    tendon=slow['mass']['tendon_dry_kg_m2']
    saved=unset['carbon']['turnover_dry_kg_m2_year']-slow['carbon']['turnover_dry_kg_m2_year']
    assert tendon>0 and saved==pytest.approx(tendon*(.1-.01))


def test_size_limits_end_on_the_lift_reserve_below_and_on_carbon_above():
    ten=product()['size_limits']['10_km']
    # family: (smallest, largest, gate that fails just above the range)
    expected=EXPECTED_RANGES_10_KM
    for name,(small,large,gate) in expected.items():
        r=ten['ranges'][name]
        assert (r['smallest_m'],r['largest_m'])==(small,large),name
        rows={x['diameter_m']:x for x in ten['rows'][name]}
        below=[d for d in rows if d<small];above=[d for d in rows if d>large]
        if below:
            assert not rows[max(below)]['gates']['mass'],name
        if above:
            assert not rows[min(above)]['gates'][gate],name
        for d in r['passing_m']:
            assert rows[d]['surplus_c_kg_m2_year']>0 and rows[d]['gas_fraction']<=.7


def test_carbon_stop_lies_between_the_rows_it_interpolates():
    for height in ('10_km','20_km'):
        ten=product()['size_limits'][height]
        for name,r in ten['ranges'].items():
            stop=r['carbon_stop_m']
            if stop is None:
                assert all(x['surplus_c_kg_m2_year']>0 for x in ten['rows'][name]) or ten['rows'][name][0]['surplus_c_kg_m2_year']<=0
                continue
            rows=sorted(ten['rows'][name],key=lambda x:x['diameter_m'])
            lo=max(x['diameter_m'] for x in rows if x['diameter_m']<=stop)
            hi=min(x['diameter_m'] for x in rows if x['diameter_m']>=stop)
            by={x['diameter_m']:x['surplus_c_kg_m2_year'] for x in rows}
            assert by[lo]>0>=by[hi] or lo==hi


def test_round_cells_and_skin_water_answer_damage_and_fire():
    for name,rows in product()['size_limits']['10_km']['rows'].items():
        for r in rows:
            if r['hydrogen_kg_m2'] is None:
                continue
            assert r['skin_water_carried_kg_m2']>=model.SKIN_WATER_KG_M2
            assert r['skin_water_carried_kg_m2']>=r['fire_skin_water_kg_m2']*(1-1e-6)
            if not name.startswith('colony'):
                assert r['compartments']>=math.ceil(r['supported_kg_m2']/model.DUMPABLE_WATER_KG_M2)
                assert r['compartments']>=model.COMPARTMENTS_ROUND and r['compartments']%2==0


def test_settled_cell_count_is_the_fewest_that_answers_a_loss():
    structural=product()['structure']
    fam=model.design_families(10,structural)
    for d in (200.,1000.,3000.):
        c,t,n=model.settle_case(d,air(),fam['round_strong'])
        assert n>=math.ceil(c['mass']['supported_kg_m2']/model.DUMPABLE_WATER_KG_M2)
        if n>model.COMPARTMENTS_ROUND:
            fewer,_,_=model.settle_case(d,air(),replace(fam['round_strong'],partition_count=1+(n-2)/2),cells=n-2)
            assert math.ceil(fewer['mass']['supported_kg_m2']/model.DUMPABLE_WATER_KG_M2)>n-2


def test_founding_module_is_the_smallest_that_passes():
    rows={r['diameter_m']:r for r in product()['size_limits']['10_km']['rows']['colony_module']}
    d=2*model.MODULE_RADIUS_M
    assert rows[d]['all_gates']
    assert not any(rows[x]['all_gates'] for x in rows if x<d)


def test_carbon_stops_are_where_the_surplus_reaches_zero():
    ten=product()['size_limits']['10_km']
    structural=product()['structure']
    fam=model.design_families(10,structural)
    stop=ten['ranges']['round']['carbon_stop_m']
    above,_,_=model.settle_case(stop*1.002,air(),fam['round'])
    below,_,_=model.settle_case(stop/1.002,air(),fam['round'])
    assert below['carbon']['remaining_assimilate_c_kg_m2_year']>0>=above['carbon']['remaining_assimilate_c_kg_m2_year']


def test_hulls_carry_their_spare_superpressure():
    structural=json.loads(model.OUT.read_text())['structure']
    fam=model.design_families(10,structural)
    assert fam['round'].bottom_overpressure_pa==pytest.approx(model.Traits().bottom_overpressure_pa+model.SPARE_TRIM_PA)
    module=model.colony_module(25.,fam)
    assert module.bottom_overpressure_pa==pytest.approx(model.Traits().bottom_overpressure_pa+model.SPARE_TRIM_PA_COLONY)
    assert module.tendon_turnover_per_year==pytest.approx(.01)
    ten=product()['size_limits']['10_km']['rows']
    for r in ten['round']:
        assert r['spare_superpressure_pa']==model.SPARE_TRIM_PA and r['design_pressure_pa']>model.SPARE_TRIM_PA


def test_walls_between_gas_cells_lose_no_hydrogen():
    one=case();cells=case(partition_count=5.)
    assert cells['hydrogen']['hydrogen_loss_kg_day']==pytest.approx(one['hydrogen']['hydrogen_loss_kg_day'])
    assert cells['mass']['partitions_kg']==pytest.approx(5*one['mass']['partitions_kg'])
    air_side=case(partition_count=5.,gas_air_partitions=2.)
    assert air_side['hydrogen']['hydrogen_loss_kg_day']>one['hydrogen']['hydrogen_loss_kg_day']
    with pytest.raises(ValueError):case(partition_count=1.,gas_air_partitions=2.)


def test_skin_water_is_carried_once():
    dry=case();wet=case(skin_water_kg_m2=.5)
    exterior=dry['skin']['area_m2']
    assert wet['mass']['skin_water_kg']==pytest.approx(.5*exterior)
    assert wet['mass']['actual_supported_mass_kg']-dry['mass']['actual_supported_mass_kg']==pytest.approx(.5*exterior,rel=1e-6)
    assert wet['mass']['structure_kg']==pytest.approx(dry['mass']['structure_kg'])


def test_ages_rise_with_size_and_stop_where_growth_stops():
    ages=product()['ages']
    for row in ages['round']:
        years=list(row['years_to_radius'].values())
        reached=[y for y in years if y is not None]
        assert reached==sorted(reached)
        if None in years:
            assert all(y is None for y in years[years.index(None):])
        stop=row['growth_stops_at_radius_m']
        for radius,y in row['years_to_radius'].items():
            if stop is not None and float(radius)>stop:
                assert y is None
    pairs={}
    for row in ages['round']:
        pairs.setdefault((row['fibre'],row['gpp_c_kg_m2_year'],row['hydrogen_route']),{})[row['growth_share']]=row
    for both in pairs.values():
        for radius,y in both[1.]['years_to_radius'].items():
            if y is not None:
                assert both[.5]['years_to_radius'][radius]==pytest.approx(2*y)


def test_colonies_double_at_construction_over_surplus():
    for row in product()['ages']['colony']:
        c=row['construction']['total_c_kg_m2'];s=row['surplus_c_kg_m2_year']*row['growth_share']
        assert row['doubling_years']==pytest.approx(c*math.log(2)/s)
        for area,y in row['years_to_area'].items():
            assert y==pytest.approx(row['doubling_years']*math.log2(float(area)/row['start_area_m2']))


def test_nursery_hydrogen_is_made_once_from_host_sugar_and_own_light():
    from research.studies.aerophytes import growth
    nursery=product()['canopy_nursery']
    for r in nursery['juveniles']:
        made=r['years']*(r['diverted_c_kg_year']/growth.fermentation_carbon_per_hydrogen(r['mol_h2_per_mol_glucose'])
                         +r['own_hydrogen_kg_year']-r['gas_loss_kg_year'])
        assert made==pytest.approx(r['hydrogen_kg'])
        assert r['gas_loss_kg_year']>0 and r['years']<r['years_on_own_light_only']
    for r in nursery['two_hundred_metre_gas']:
        assert r['sugar_c_kg']==pytest.approx(r['hydrogen_kg']*growth.fermentation_carbon_per_hydrogen(r['mol_h2_per_mol_glucose']))
        assert r['hectare_years']==pytest.approx(r['sugar_c_kg']/(r['npp_c_kg_m2_year']*1e4))


def test_rain_closure_matches_the_zonal_rain_it_reads():
    from research.studies.aerophytes import water
    out=product()['water']
    zonal=out['zonal_rain']
    for r in out['rain_closure']:
        edge=water.rain_closing_latitude(zonal,r['need_mm_day'],r['catch_share'],r['height_factor'])
        assert (r['north_deg'],r['south_deg'])==(edge['north_deg'],edge['south_deg'])
        rows=[(l,p*r['catch_share']*r['height_factor']) for l,p in zip(zonal['latitude_deg'],zonal['rain_mm_day'])]
        met=[l for l,p in rows if p>=r['need_mm_day'] and l>0]
        if met and r['north_deg'] is not None:
            assert max(met)<=r['north_deg']
        assert 0<=r['moon_share']<=1


def test_dry_gaps_use_the_sailing_periods_and_daylight_only():
    data=product()
    days={(x['latitude_deg'],x['height_km']):x['relative_period_days'] for x in data['sailing']['day_length']}
    rows={x['path']:x for x in data['water']['dry_intervals']}
    assert rows['equator_riding_10_km_mean_wind']['relative_period_days']==pytest.approx(days[(0.,10.)])
    assert rows['equator_riding_20_km_mean_wind']['relative_period_days']==pytest.approx(days[(0.,20.)])
    need={n['name']:n['water_kg_m2_day'] for n in data['water']['needs']}
    for r in rows.values():
        assert r['daylight_dry_days']<r['dry_days']
        assert r['storage_kg_m2']['reference']==pytest.approx(2*need['reference']*r['daylight_dry_days'])


def test_lifespans_follow_their_strike_rates():
    for r in product()['lifespans']:
        assert r['survive_300_years_if_fire_spreads']==pytest.approx(math.exp(-300*r['strikes_per_year']))
        assert r['survive_1000_years_if_fire_spreads']==pytest.approx(math.exp(-1000*r['strikes_per_year']))
        assert r['strikes_in_300_years']==pytest.approx(300*r['strikes_per_year'])
        assert r['strikes_per_century']==pytest.approx(100*r['strikes_per_year'])
