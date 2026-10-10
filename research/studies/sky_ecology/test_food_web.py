"""Carbon and energy closures, allometries and limits of the food-web component."""
import itertools
import json
import math

import pytest

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from research.studies.sky_ecology import food_web as fw
from research.studies.sky_ecology import run as m
from research.studies.sky_ecology import stoichiometry as st


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())['food_web']


@pytest.fixture(scope='module')
def ships():
    return json.loads((m.ROOT/'research/studies/sky_ships/results/sky_ships.json').read_text())


@pytest.fixture(scope='module')
def air10(ships):
    return next(a for a in ships['air'] if a['height_km'] == 10)


def test_production_partition_closes_the_carbon_ledger(product):
    for r in product['production']:
        assert abs(r['ledger_residual_c']) < 1e-12
        assert math.isclose(r['renewal_biomass_c'], r['green_tissue_c']+r['structure_c'])
        assert math.isclose(r['biomass_c'], r['renewal_biomass_c']+r['growth_biomass_c'])
        assert math.isclose(r['surplus_c'], r['growth_biomass_c']+r['growth_respiration_c']+r['growth_gas_and_reserve_c'])
        for h, mature in zip(fw.HAZARD_PER_YEAR, r['mature_biomass_c']):
            assert math.isclose(mature, r['renewal_biomass_c']+r['standing_host_c_kg_m2']*h)
        assert 0 < r['biomass_share_of_growth'] < 1
        assert math.isclose(r['growth_rate_per_year'], r['growth_biomass_c']/r['standing_host_c_kg_m2'])


def test_biomass_leaves_the_tenants_sugar_out(product):
    sky = product['sky']
    for state in ('growing', 'mature'):
        lo, hi = sky['kg_c_m2_year'][state]
        assert math.isclose(sky['mt_c_year_at_1pct'][state][0], fw.sky_totals(lo, .01))
        assert lo < hi
    rows = [r for r in product['production'] if r['all_gates'] and r['gpp_c_kg_m2_year'] == 2.]
    assert min(r['biomass_c'] for r in rows) == sky['kg_c_m2_year']['growing'][0]
    assert all(r['food_to_tenants_c'] == sky['kg_c_m2_year']['food_to_tenants'][0] for r in rows)


def test_consumer_level_closes():
    lv = fw.consumer_level(1., .5, .4, .9, .3)
    assert math.isclose(lv['consumption_c'], lv['egestion_c']+lv['assimilation_c'])
    assert math.isclose(lv['assimilation_c'], lv['production_c']+lv['respiration_c'])
    assert math.isclose(lv['assimilation_c'], 1.*.4+.5*.9)
    with pytest.raises(ValueError):
        fw.consumer_level(1., 0., 1.2, .9, .3)


def test_trophic_chain_litter_accounts_for_every_flow(product):
    for key, ch in product['chains'].items():
        for edge in ('low', 'high'):
            t = ch[edge]
            offer = fw.tree_offer(ch['production'], edge) if key != 'eaten_base_10m' else fw.eaten_offer(ch['production'], edge)
            assert math.isclose(t['litter_c'], offer['grazable_c']-t['grazed_c']+offer['shed_c']+t['level2']['egestion_c'])
            assert t['level3_production_c'] < t['level2']['production_c'] < t['level2']['consumption_c']
            assert t['level2_production_c'] <= t['level2']['production_c']+1e-15
            assert t['great_flyer_food_c'] <= t['level3_production_c']
        assert ch['low']['level2_production_c'] < ch['high']['level2_production_c']


def test_edges_bound_every_corner_of_the_phosphorus_limit(product):
    for key in ('sky_reef_60m', 'eaten_base_10m'):
        ch = product['chains'][key]
        offer = fw.tree_offer(ch['production']) if key != 'eaten_base_10m' else fw.eaten_offer(ch['production'])
        produced = []
        for gs, a, pe, fp, pa, cp in itertools.product(offer['grazed_share_range'], fw.RANGES['herbivore_assimilation'],
                                                       fw.RANGES['poikilotherm_production_efficiency'], offer['food_p_g_kg'],
                                                       st.CONSUMER_P_ASSIMILATION, st.TENANT_P_G_KG):
            lv = fw.consumer_level(gs*offer['grazable_c'], offer['food_to_tenants_c'], a, fw.SUGAR_ASSIMILATION, pe)
            produced.append(st.consumer_p_balance(gs*offer['grazable_c'], lv['production_c'], fp, pa, cp)['p_limited_production_c'])
        assert ch['low']['level2_production_c'] <= min(produced)+1e-15
        assert ch['high']['level2_production_c'] >= max(produced)-1e-15
        lo, hi = ch['p_corners']['share_range']
        assert 0 < lo < hi <= 1 and ch['p_corners']['corners'] == 64


def test_eaten_base_feeds_on_living_tissue_only(product):
    per_living = {r['edible_c_kg_m2']/r['living_dry_kg_m2'] for r in product['eaten_base']}
    assert max(per_living)-min(per_living) < 1e-12 and min(per_living) >= fw.C_PER_DRY   # living tissue with its resident community
    for r in product['eaten_base']:
        shares = r['construction_shares']
        assert math.isclose(sum(shares.values()), 1.)
        assert math.isclose(r['edible_share_of_construction'], shares['edible_tissue'])
        assert shares['edible_tissue'] < shares['envelope'] < shares['gas']
        assert math.isclose(r['turnover_per_year'], r['surplus_c_kg_m2_year']/r['construction_c_kg_m2'])
        assert math.isclose(r['edible_production_c_kg_m2_year'], r['turnover_per_year']*r['edible_c_kg_m2'])
        assert math.isclose(r['envelope_production_c_kg_m2_year'], r['turnover_per_year']*r['envelope_c_kg_m2'])
    small = [r for r in product['eaten_base'] if r['gas_route'] == 'photolysis' and r['gpp_c_kg_m2_year'] == 2.]
    assert all(a['turnover_per_year'] > b['turnover_per_year'] for a, b in zip(small, small[1:]))
    offer = fw.eaten_offer(product['chains']['eaten_base_10m']['production'])
    assert offer['shed_c'] == product['chains']['eaten_base_10m']['production']['envelope_production_c_kg_m2_year']


def test_invertebrate_pb_follows_banse_and_mosher():
    one = fw.invertebrate_pb(1.)
    ten = fw.invertebrate_pb(10.)
    assert 4.5 < one < 5.3 and 1.9 < ten < 2.3
    assert math.isclose(one/ten, 10**.37, rel_tol=1e-9)


def test_allometry_lies_between_the_similarity_terms(product):
    """The allometric ratio is bound to lie between the flight and upkeep factors; the agreement is consistency."""
    people = json.loads((m.ROOT/'biosphere/ecology/results/people.json').read_text())['lift']['flyers']
    assert math.isclose(fw.bird_fmr_j_day(1.), 10.5*1000**.681*1e3)
    for f in product['flyers']:
        src = people[f['flyer']]
        i = 0 if f['edge'] == 'low' else 1
        k = src['moon_span_m'][i]/src['earth_span_m'][i]
        assert math.isclose(f['length_factor'], k)
        lo, hi = f['similarity_bounds']
        assert lo < f['allometric_ratio'] < hi
        assert all(lo <= s <= hi for s in f['similarity_ratio'])
        assert math.isclose(f['allometric_ratio'], (f['mass_kg']/f['earth_mass_kg'])**fw.BIRD_FMR[1])


def test_flight_per_km_is_a_sixth_of_earths():
    moon = fw.flight_j_per_km(1000., 20., .25)
    earth = fw.flight_j_per_km(1000., 20., .25, STANDARD_GRAVITY)
    assert math.isclose(moon/earth, MOON_SURFACE_GRAVITY/STANDARD_GRAVITY)
    assert math.isclose(moon, 1000.*MOON_SURFACE_GRAVITY*1e3/20./.25)


def test_glide_ratios_lie_within_the_sky_fleets_gliders(product):
    best = list(product['sky_fleet_best_glide'].values())
    lo, hi = fw.RANGES['glide_ratio']
    assert min(best) < lo < hi < max(best)


def test_flyers_fed_inverts_fmr_and_carries_the_residual(product):
    fmr = 2e8
    food = fmr*365.25/(fw.FOOD_J_PER_KG_C*.5)
    assert math.isclose(fw.flyers_fed(food, .5, fmr), 1.)
    for row in product['flyers_fed']:
        assert math.isclose(row['as_grazers_over_fmr_residual'][0], row['as_grazers'][0]/fw.BIRD_FMR_RESIDUAL[1])
        assert math.isclose(row['as_grazers_over_fmr_residual'][1], row['as_grazers'][1]/fw.BIRD_FMR_RESIDUAL[0])


def test_small_body_sits_on_the_gas_gate(air10):
    from research.studies.aerophytes import run as aero, structure
    s = structure.evaluate(round_spare_pa=aero.SPARE_TRIM_PA, module_spare_pa=aero.SPARE_TRIM_PA_COLONY)
    c, t = fw.small_body(10., air10, s)
    assert abs(c['mass']['required_lifting_mix_volume_fraction']-.7) < 1e-9
    assert c['traits']['skin_water_kg_m2'] == 0.
    bigger, tb = fw.small_body(20., air10, s)
    assert tb.living_dry_kg_m2 > t.living_dry_kg_m2


def test_reference_ten_metre_check(product, air10):
    chk = product['reference_small_check']['with_ballonet_partition']
    assert math.isclose(chk['gross_lift_kg_m2'], 2/3*10*air10['lift_hydrogen_kg_m3'], rel_tol=1e-9)
    assert abs(chk['gas_fraction']-.7) < 1e-9
    assert .33 < chk['living_share_of_reference'] < .4


def test_filter_return_and_open_sky(product):
    assert fw.filter_return(4.46, 4.46) == 1.
    assert math.isclose(fw.filter_return(2., 4.), .5)
    lo, hi = product['filters']['open_sky_aeroplankton_mg_m3']
    assert math.isclose(lo, 1e-5, rel_tol=1e-9) and math.isclose(hi, 1e-4, rel_tol=1e-9)
    assert max(max(v) for v in product['filters']['energy_return'].values()) < 1e-3


def test_fall_speed_follows_the_air_and_drift_integrates_the_wind(ships):
    v = fw.fall_speed(1.2, areal_density_kg_m2=.1, gravity=1.6)
    assert math.isclose(v, math.sqrt(2*.1*1.6/(1.2*1.2)))
    assert fw.fall_speed(1.0, areal_density_kg_m2=.1) > fw.fall_speed(1.4, areal_density_kg_m2=.1)
    flat = [dict(height_km=0., density_kg_m3=1.2), dict(height_km=40., density_kg_m3=1.2)]
    out = fw.fall(10., flat, [0., 20.], [2., 2.], areal_density_kg_m2=.1)
    speed = fw.fall_speed(1.2, areal_density_kg_m2=.1)
    assert math.isclose(out['hours'], 1e4/speed/3600, rel_tol=1e-9)
    assert math.isclose(out['drift_km'], 2.*1e4/speed/1e3, rel_tol=1e-9)
    for a in ships['air']:
        assert math.isclose(fw.density_at(a['height_km'], ships['air']), a['density_kg_m3'], rel_tol=1e-12)
    with pytest.raises(ValueError):
        fw.fall_speed(0., areal_density_kg_m2=.1)


def test_sky_snow_reads_the_wind_product(product):
    snow = product['sky_snow']
    winds, _ = m.cross_branch_json('climate/results/gcm/zonal_winds_A28_dim5_moon.json')
    if winds is None:
        pytest.skip('zonal-wind product unavailable')
    assert snow['model_cell_km'] == winds['trajectories']['cell_km']
    for start in ('from_10_km', 'from_20_km'):
        lat = snow[start]['slowest_flakes_by_latitude']
        assert lat['largest_northward_km'] < snow['model_cell_km']


def test_sky_totals_match_the_aerial_study():
    aerial = json.loads((m.ROOT/'research/studies/aerial_ecology/results/aerial_ecology.json').read_text())
    assert math.isclose(fw.sky_totals(1., .01)*1e9, aerial['central']['projected_area_m2'], rel_tol=1e-9)


def test_grazing_tolerance_scales_with_share(product):
    tol = product['chains']['sky_reef_60m']['tolerance']
    assert tol[0]['regrowth_c'] < tol[1]['regrowth_c']
    assert all(0 < t['share_of_surplus'] < .2 for t in tol)
    assert [t['grazed_share'] for t in tol] == list(fw.RANGES['tree_grazed_share'])


@pytest.mark.parametrize('food_p', [(0., 0.), (1., 3.), (100., 100.)])
def test_phosphorus_limited_ledger_keeps_unbuilt_carbon_explicit(food_p):
    offer = dict(grazable_c=.225, grazed_share=.12, food_to_tenants_c=.1,
                 food_p_g_kg=food_p, shed_c=.4)
    for edge in ('low', 'high'):
        row = fw.trophic_chain(offer, edge)
        level = row['level2']
        assert level['production_c'] == row['level2_production_c']
        assert level['unallocated_assimilate_c'] >= 0
        assert math.isclose(level['consumption_c'], level['egestion_c']+
                            level['respiration_c']+level['production_c']+
                            level['unallocated_assimilate_c'])
        assert math.isclose(row['carbon_limited']['level2_production_c'],
                            level['production_c']+level['unallocated_assimilate_c'])
        if food_p == (0., 0.):
            assert row['level2_production_c'] == row['level3_production_c'] == 0
            assert level['unallocated_assimilate_c'] > 0
        if food_p == (100., 100.):
            assert level['unallocated_assimilate_c'] == 0
