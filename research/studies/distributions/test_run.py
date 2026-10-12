"""Bindings, conservation and rules of the distributions product."""
import json
import math
import re

import pytest

from shared.provenance import constants_changed
from research.studies.distributions import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def missing_inputs():
    missing = [p for p in m.INPUTS if not (m.ROOT / p).exists()]
    for key in m.BOUND:
        try:
            m.bound_bytes(key)
        except m.MissingInput as error:
            missing.append(str(error))
    return missing


@pytest.fixture(scope='module')
def fresh():
    missing = missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    return json.loads(m.dumps(m.results()))


@pytest.fixture(scope='module')
def basis():
    missing = missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    d = {k: m.load_json(m.KEYS[k], m.INPUTS[m.KEYS[k]]) for k in m.KEYS}
    b = {k: m.load_bound(k) for k in m.BOUND if m.BOUND[k]['schema']}
    return m.build_basis(d, b), d


def close(a, b, rel=1e-6, abs_=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


# ---- bindings ---------------------------------------------------------------------------------------------------

def test_product_regenerates(product, fresh):
    assert product['schema'] == m.SCHEMA
    assert product == fresh


def test_producer_and_input_hashes(product):
    producer = product['producer']
    for path, h in {**producer['files'], **producer['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    assert set(producer['inputs']) == set(m.INPUTS)
    assert producer['sources'] == m.digest(m.HERE / 'sources.json')
    assert not constants_changed(producer['constants'])


def test_bound_inputs_pinned(product):
    bound = product['producer']['bound_inputs']
    assert set(bound) == set(m.BOUND)
    for key, b in bound.items():
        pinned = m.BOUND[key]
        assert (b['path_at_integration'], b['branch'], b['commit'], b['sha256'], b['blob']) == \
            (pinned['path'], pinned['branch'], pinned['commit'], pinned['sha256'], pinned['blob'])
        try:
            raw = m.bound_bytes(key)
        except m.MissingInput:
            pytest.skip(f'{key} is not readable here')
        if pinned['schema']:
            assert json.loads(raw)['schema'] == pinned['schema']


def test_rejects_wrong_schema_and_missing(tmp_path, monkeypatch):
    (tmp_path / 'wrong.json').write_text('{"schema": "wrong"}')
    monkeypatch.setattr(m, 'ROOT', tmp_path)
    with pytest.raises(ValueError):
        m.load_json('wrong.json', 'terluna.research.food-systems/1')
    with pytest.raises(m.MissingInput):
        m.load_json('absent.json', 'terluna.research.food-systems/1')


def test_transcribed_phrases_occur_in_the_bound_note(product):
    try:
        text = m.bound_bytes('industry_md').decode('utf-8')
    except m.MissingInput:
        pytest.skip('the industry note is not readable here')
    for key, entry in m.TRANSCRIBED.items():
        assert entry['phrase'] in text, key
    assert all(product['checks']['transcribed_phrases'].values())
    for name in m.PLACES_AU:                       # each place's row carries its distance and its sunlight
        row = next(line for line in text.splitlines() if line.startswith(m.TRANSCRIBED[name]['phrase']))
        assert m.TRANSCRIBED[f'{name}_sunlight']['phrase'] in row


def test_every_published_figure_has_a_source():
    sources = {s['key'] for s in json.loads((m.HERE / 'sources.json').read_text())['sources']}
    used = {src for entry in m.LIT.values() for src in entry['source']}
    assert used == sources


def test_design_guesses_are_ranges(product):
    for key, g in m.GUESSES.items():
        lo, hi = g['value']
        assert lo < hi, key
        assert g['unit'] and g['note'], key
    assert product['assumptions']['design_guesses'].keys() == m.GUESSES.keys()
    assert len(m.COMBOS) == 8 and len(set(m.COMBOS)) == 8


def test_every_input_feeds_the_basis(basis):
    B, _ = basis
    assert B['sky']['powered_m_s'] > 0 and B['sky']['duty'][0] < B['sky']['duty'][1]      # inhabited volume
    assert B['sky']['breathes_like_earth_m'] > 0 and B['earth']['transfer_days']          # daily life
    assert B['earth']['shared_cases_billion'][0] <= B['earth']['shared_cases_billion'][1]  # shared cases


# ---- the Moon's placement and its accounts ----------------------------------------------------------------------

def test_people_are_conserved_at_every_combination(basis):
    B, _ = basis
    for moon in m.MOON_BILLION:
        for aloft in m.ALOFT_BILLION:
            ground = moon - aloft
            if ground < 1:
                continue
            for combo in m.COMBOS:
                for rule in ('spread', 'twilit'):
                    for dry_rule in ('land', 'poleward'):
                        place = m.placement(B, ground, combo, rule, dry_rule, 0.5)
                        assert close(sum(place['people'].values()), ground, rel=1e-12)
                        assert close(sum(place['scattered_by_home'].values()), place['people']['scattered'],
                                     rel=1e-12)
                        assert close(sum(place['dry_people_by_band'].values()), place['dry_highland'], rel=1e-12)
                        assert close(place['dry_highland'],
                                     place['highland_share'] * ground - B['metropolis']['billion'])


def test_settled_areas_follow_densities(basis):
    B, _ = basis
    place = m.placement(B, 3.0, (1, 1, 0))
    a, p = place['areas_km2'], place['people']
    assert close(a['tower_districts'], p['tower_districts'] * 1e9 / B['towers']['people_per_km2'])
    assert close(a['highland_towns'], p['highland_towns'] * 1e9 / B['towns']['people_per_km2'])
    assert close(a['scattered'], p['scattered'] * 1e9 * B['scattered']['m2_per_person'][0] / 1e6)
    assert close(a['metropolis'], B['metropolis']['dry_land_km2'])
    assert close(sum(place['dry_area_by_band'].values()), a['tower_districts'] + a['highland_towns'])
    least, most = m.placement(B, 3.0, (1, 1, 0)), m.placement(B, 3.0, (0, 0, 1))
    assert least['settled_km2'] < most['settled_km2']


def test_the_guesses_are_untied(product):
    """The dry high ground's settled share spans every combination of the two guesses (round 1, R4)."""
    land = product['moon']['cases']['6_3']['kept_30']['full']['land']
    lo, hi = land['dry_high_ground_settled_share']
    assert hi / lo > 1.5


def test_food_land_replaces_the_food_studys_towns(basis):
    B, d = basis
    food = d['food']
    for moon in m.MOON_BILLION:
        for kept in m.KEPT:
            fb = m.food_block(food, m.CHOSEN_MIX, moon, kept)
            blk = fb['block']
            for e in (0, 1):
                towns = m.food_town_km2(B, blk['people_billion']['surface'][e], e, fb['reach'])
                f = blk['fields']
                parts = (sum(f['field_km2_by_zone'][z][e] for z in f['field_km2_by_zone'])
                         + f['fog_desert_watered_km2'][e] + blk['land_km2']['glass'][e])
                assert close(parts + towns, blk['land_km2']['total'][e], rel=5e-3)


def test_chosen_mix_fits_the_food_studys_zones(basis):
    B, d = basis
    for moon in m.MOON_BILLION:
        for kept in m.KEPT:
            fb = m.food_block(d['food'], m.CHOSEN_MIX, moon, kept)
            assert fb is not None and fb['both_ends']


def test_land_bands_partition_the_dry_land_and_flag_overfill(product, basis):
    B, d = basis
    lb = B['land_bands']
    total = sum(lb['dry_land_km2'][k] for k in m.FOOD_BANDS)
    assert close(total, sum(lb['zones_km2'].values()), rel=5e-3)
    for kept in m.KEPT:
        c = m.moon_corner(B, d['food'], 8, 4, kept, 'fusion', (1, 0, 0), 1)
        assert close(sum(c['bands']['capacity_km2'].values()), total * (1 - kept), rel=5e-3)
        assert c['buildable'] == all(v <= 1.0 for v in c['bands']['overfill_km2'].values())
        assert c['buildable'] == all(c['bands']['use_km2'][b] <= c['bands']['capacity_km2'][b] + 1.0
                                     for b in m.LAND_BANDS)
    for case in product['moon']['cases'].values():
        for entry in case.values():
            assert entry['buildable_corners'] + len(entry['overfilling_corners']) == entry['corners']
            if entry['feasible'] is None:
                assert entry['buildable_corners'] == entry['corners']
            else:
                f, p = entry['full']['land'], entry['feasible']['land']
                assert f['wild_share_outside_kept'][0] <= p['wild_share_outside_kept'][0]
                assert p['wild_share_outside_kept'][1] <= f['wild_share_outside_kept'][1]
                assert entry['remedy'] is not None


def test_the_fields_reach_the_dry_high_ground_by_the_food_rule(basis):
    """At 8 billion with a kept half and the high-cost end the fill reaches the 30-45 degree highlands (round 1, R1)."""
    B, d = basis
    c = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 1)
    hb = B['land_bands']['highland_km2']
    fb = m.food_block(d['food'], m.CHOSEN_MIX, 8, 0.5)
    towns = m.food_town_km2(B, fb['block']['people_billion']['surface'][1], 1, fb['reach'])
    expect = hb['30-45'] * 0.5 - towns * hb['30-45'] / sum(hb.values())
    assert close(c['bands']['fields_on_dry_high_ground_km2'], expect, rel=2e-2)
    assert c['land']['dry_high_ground_used_share'] > c['land']['dry_high_ground_settled_share']
    low = m.moon_corner(B, d['food'], 8, 4, 0.3, 'fusion', (0, 1, 0), 1)
    assert low['bands']['fields_on_dry_high_ground_km2'] == 0.0


def test_the_middle_wet_zone_is_its_own_band(basis):
    """Rain-fed fields fill the middle-wet zone at 8 billion with a kept half and the high-cost end; the scattered
    communities' middle-band land covers the zone by area unless they keep to the fog-desert side (round 2)."""
    B, d = basis
    lb = B['land_bands']
    assert close(lb['middle_wet_share_of_middle_homes'], lb['zones_km2']['middle_wet'] / lb['middle_band_homes_km2'])
    assert 0.0 < lb['middle_wet_share_of_middle_homes'] < 1.0
    c = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 1)
    bands = c['bands']
    assert close(bands['fields_on_middle_wet_km2'], lb['zones_km2']['middle_wet'] * 0.5, rel=1e-3)
    assert close(bands['scattered_on_middle_wet_km2'],
                 bands['scattered_middle_band_km2'] * lb['middle_wet_share_of_middle_homes'])
    assert not c['buildable'] and bands['overfill_km2']['middle_wet_zone'] > 1.0
    dry = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 1, wet_rule='dry_side')
    assert dry['buildable'] and dry['bands']['scattered_on_middle_wet_km2'] <= 1.0
    assert close(dry['land']['used_km2'], c['land']['used_km2'])
    assert dry['bands']['middle_homes_dry_side_used_share'] < 1.0
    low = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 0, wet_rule='dry_side')
    ref = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 0)
    assert close(low['bands']['scattered_on_middle_wet_km2'], ref['bands']['scattered_on_middle_wet_km2'])


def test_the_remedy_moves_the_dry_highland_towns_poleward(product, basis):
    B, d = basis
    base = m.moon_corner(B, d['food'], 8, 2, 0.5, 'fusion', (1, 0, 1), 1)
    two = m.moon_corner(B, d['food'], 8, 2, 0.5, 'fusion', (1, 0, 1), 1, rule='twilit', dry_rule='poleward')
    rem = m.moon_corner(B, d['food'], 8, 2, 0.5, 'fusion', (1, 0, 1), 1, rule='twilit', dry_rule='poleward',
                        wet_rule='dry_side')
    assert not base['buildable'] and not two['buildable'] and rem['buildable']
    assert rem['dry_people_by_band']['30-45'] < base['dry_people_by_band']['30-45']
    assert rem['water']['lift_gw'] > base['water']['lift_gw']
    assert rem['scattered_share_0_30'] < base['scattered_share_0_30']
    for case in product['moon']['cases'].values():
        r = case['kept_50']['remedy']
        if r:
            assert r['buildable_corners'] == r['buildable_corners_by_parts']['all_three'] == 16
            assert r['buildable_corners_by_parts']['poleward'] < 16


def test_energy_bases_and_heat_components(basis):
    B, d = basis
    c = m.moon_corner(B, d['food'], 6, 3, 0.3, 'fusion', (0, 0, 1), 1)
    land, en = c['land'], c['energy']
    assert close(land['wild_share_outside_kept'], 1 - land['used_km2'] / (B['moon']['dry_land_km2'] * 0.7))
    assert close(land['used_km2'], land['food_km2'] + land['grazing_km2'] + land['settled_km2'])
    assert close(en['electricity_tw'], en['food_tw'] + en['sky_tw'] + en['water_lift_tw'])
    assert close(en['primary_tw'], en['resident_primary_tw']
                 + en['electricity_tw'] * B['energy']['heat_per_w']['fusion'][1])
    assert close(en['heat_tw'], sum(en['heat_by_component_tw'].values()))
    assert close(sum(en['heat_share'].values()), 1.0)
    assert close(en['heat_w_m2'], en['heat_tw'] * B['energy']['w_m2_per_tw'])
    assert close(en['added_dimming_percent'], en['heat_tw'] * B['energy']['dimming_percent_per_tw'][1])
    assert close(en['sky_tw'], sum(en['sky_parts_tw'].values()))
    # at the high-cost end the hydrogen remade gives the largest share of the heat (round 1, R8)
    assert max(en['heat_share'], key=en['heat_share'].get) == 'sky_hydrogen_remade'
    a = m.moon_corner(B, d['food'], 6, 3, 0.3, 'array', (0, 0, 1), 1)
    assert a['energy']['heat_tw'] < en['heat_tw']


def test_phosphorus_follows_the_need(basis):
    B, d = basis
    c = m.moon_corner(B, d['food'], 8, 4, 0.3, 'fusion', (0, 0, 0), 1)
    p, ph = c['phosphorus'], B['phosphorus']
    fb = m.food_block(d['food'], m.CHOSEN_MIX, 8, 0.3)
    assert close(p['need_kt'], fb['block']['phosphorus_mt_year']['lost'][1] * 1e3)
    assert close(p['need_kt'], 8 * ph['design_kg_per_person'][1] * 1e9 / 1e6, rel=1e-2)
    assert close(p['dredged_sediment_gt'], p['need_kt'] * 1e6 * ph['sediment_kg_per_kg'][1] / 1e12)
    assert close(p['share_of_delta_catch'], p['need_kt'] / (ph['delta_particles_kg_yr'][0] / 1e6))


def test_lit_sky_by_design(product, basis):
    B, _ = basis
    rule = product['checks']['lit_sky_rule']
    for k in rule['bound']:
        for a, b in zip(rule['bound'][k], rule['recomputed'][k]):
            assert close(a, b, rel=5e-3)
    place = m.placement(B, 4.0, (0, 0, 1))
    for e in (0, 1):
        for des in m.DESIGNS:
            lit = m.lit_sky(B, place, 4.0, e, des)
            assert lit['scattered'] <= sum(B['scattered']['homes_km2'].values())
            # each design takes its own metropolis figure (round 1, R2)
            assert close(lit['metropolis_lit_as_a_town'],
                         B['metropolis']['billion'] * 1e9 * B['light']['metropolis_town_lit'][des][e] / 1e6)
    twilit = m.placement(B, 4.0, (0, 0, 1), 'twilit')
    assert m.lit_sky(B, twilit, 4.0, 1, 'shielded')['scattered'] < m.lit_sky(B, place, 4.0, 1, 'shielded')['scattered']


# ---- the sky -------------------------------------------------------------------------------------------------------

def test_spacing_and_encounter_rules_reproduce_their_products(product):
    s = product['checks']['spacing_rule']
    assert close(s['spacing_km'], s['ways_of_living_km'], rel=5e-3)
    e = product['checks']['encounter_rule']
    for a, b in zip(e['recomputed'], e['sky_ecology']):
        assert close(a, b, rel=5e-3)


def test_traffic_scales_and_fleet_totals_use_the_met_density(basis):
    B, _ = basis
    two, four = m.traffic(B, 2.0, 0.3), m.traffic(B, 4.0, 0.3)
    assert close(four['towns'], 2 * two['towns'])
    for a, b in zip(two['encounters_per_town_day_busiest_band']['drifting'],
                    four['encounters_per_town_day_busiest_band']['drifting']):
        assert close(b, 2 * a)
    for a, b in zip(two['encounters_fleet_day']['powered'], four['encounters_fleet_day']['powered']):
        assert close(b, 4 * a)
    # the fleet total takes the occupancy-weighted density, below the busiest band's (round 1, R6)
    ratio = two['met_density_over_even'] / two['busiest_band_gathering']
    for busy, mean in zip(two['encounters_per_town_day_busiest_band']['drifting'],
                          two['encounters_per_town_day_mean']['drifting']):
        assert close(mean, busy * ratio)
    for mean, fleet in zip(two['encounters_per_town_day_mean']['drifting'], two['encounters_fleet_day']['drifting']):
        assert close(fleet, two['towns'] * mean / 2)
    half = m.traffic(B, 2.0, 0.5)
    assert half['spacing_even_km'] < two['spacing_even_km']
    assert close(two['spacing_even_km'] / two['spacing_busiest_band_km'], math.sqrt(two['busiest_band_gathering']))


def test_wind_spread_is_positive(basis):
    B, _ = basis
    w = B['sky']['wind_10_km']
    assert w['sigma_u_m_s'] > 0 and w['sigma_v_m_s'] > 0
    assert close(w['uncorrelated_relative_m_s'], math.sqrt(2 * (w['sigma_u_m_s'] ** 2 + w['sigma_v_m_s'] ** 2)))


# ---- tubes, synergies, Earth, the array and the colonies ----------------------------------------------------------

def test_great_circle():
    quarter = m.gc_km(0.0, 0.0, 0.0, 90.0)
    assert close(quarter, math.pi / 2 * m.MOON_RADIUS / 1e3)
    assert close(m.gc_km(10.0, 20.0, 10.0, 20.0), 0.0, abs_=1e-6)


def test_tube_pairs_and_ties(product):
    t = product['tubes_and_undersea']
    assert t['undersea_domes'] == len(t['pairs']) == 13
    joined = t['joined']
    assert set(joined['sites']) == {p['site'] for p in t['pairs'] if p['km'] <= t['join_km']}
    assert all(p['km'] > 0 for p in t['pairs'])
    assert close(joined['domes_million'], len(joined['sites']) * t['dome_residents'] / 1e6, rel=5e-3)
    # tied candidates all join (round 1, R5): Lunokhod 1 lies equally far from rilles D and E1
    tie = next(p for p in t['pairs'] if p['site'].startswith('Lunokhod 1'))
    assert tie['tied'] and len(tie['tubes']) == 2
    assert set(tie['tube_keys']) <= set(joined['tubes'])


def test_earth_block_and_sensitivity(product, basis):
    B, _ = basis
    low, high = m.earth_block(B, 9.0), m.earth_block(B, 12.0)
    assert low['within_boundaries'] and low['off_land_food_tw'] == [0.0, 0.0]
    assert not high['within_boundaries']
    beyond = 12.0 - B['earth']['within_boundaries_billion']
    assert close(high['off_land_food_tw'][1], beyond * B['food']['power_made_kw'][1])
    assert close(high['farmland_km2'][0], B['earth']['within_boundaries_billion']
                 * B['earth']['farmland_km2_per_billion'][0])
    assert close(high['energy_tw'][1], 12.0 * B['energy']['band_kw'][1]
                 + high['off_land_food_tw'][1] * B['energy']['heat_per_w']['fusion'][1])
    sens = product['earth_sensitivity']
    assert set(sens) == {'8', '9', '11'}
    assert sens['8']['off_land_food_tw'] == [0.0, 0.0] and sens['11']['beyond_boundaries_billion'] > 0
    for k, v in sens['11']['totals_billion'].items():
        assert close(v, product['distributions'][k]['people_billion']['total'] + 2.0)


def test_colonies_credit_no_capacity_and_count_the_mirrors(product):
    c = product['colonies']
    assert c['capacity_credited'] is False
    k = m.LIT['kalpana_one']['value']
    assert close(c['habitat_gt']['kalpana_one'], k['total_mass_t'] / k['residents'], rel=5e-3)
    farm = c['food']['sunlit_farm_m2_per_person_at_1_au']
    for place in c['food']['sunlit_farm_mirrors'].values():
        assert close(place['mirrors_beyond_the_farm_m2'], place['collecting_area_m2'] - farm, rel=5e-3)


def test_distributions_respect_the_working_ranges(product):
    dists = product['distributions']
    assert 3 <= len(dists) <= 5
    for key, dist in dists.items():
        p = dist['people_billion']
        assert 20.0 <= p['total'] <= 30.0, key
        assert close(p['total'], p['earth'] + p['moon'] + p['array'] + p['colonies'])
        assert 4.0 <= p['moon'] <= 8.0 and 2.0 <= p['moon_aloft'] <= 4.0, key
        assert p['earth'] < 10.0, key
        assert dist['assumptions'] and dist['uncertainties']
        e = dist['energy_tw']
        for i in (0, 1):
            assert close(e['total'][i], e['earth'][i] + e['moon'][i] + e['array'][i] + e['colonies'][i], rel=2e-2)
    four = [d for d in dists.values() if d['people_billion']['moon'] == 4.0]
    assert four
    space = max(d['people_billion']['array'] + d['people_billion']['colonies'] for d in dists.values()
                if d['people_billion']['total'] == four[0]['people_billion']['total'])
    assert four[0]['people_billion']['array'] + four[0]['people_billion']['colonies'] == space


def test_the_authors_distribution_keeps_cs_totals(product):
    c, cp = product['distributions']['C'], product['distributions']["C'"]
    assert cp['people_billion'] == c['people_billion']
    assert cp['moon']['case'] == 'cprime'
    for kept in ('kept_30', 'kept_50', 'wild_30_hybrid_20', 'wild_50_hybrid_20'):
        assert product['moon']['cprime'][kept]['corners'] == 32


def test_cprime_moves_the_middle_latitudes_and_keeps_homesteads(basis):
    """Homesteads at the guessed density among the middle-wet fields; the rest of the 30-48.6 degree lowland people
    placed in the ranked order, each destination within its spare land; everyone conserved (C')."""
    B, d = basis
    order = ['high_latitude_lowland', 'dry_60_90', 'dry_45_60', 'dry_30_45', 'lowland_0_30']
    for kept, e, h in ((0.3, 1, 1), (0.5, 0, 0), (0.5, 1, 0), (0.5, 1, 1)):
        for combo in ((0, 0, 0), (1, 1, 1), (0, 1, 1)):
            base = m.moon_corner(B, d['food'], 8, 4, kept, 'fusion', combo, e)
            c = m.moon_corner(B, d['food'], 8, 4, kept, 'fusion', combo, e, middle=dict(order=order, h_end=h))
            mv = c['moved']
            density = m.GUESSES['homestead_people_per_km2']['value'][h]
            place = m.placement(B, 4.0, combo, kept=kept)
            middle = place['people']['scattered'] * place['scattered_share_middle']
            assert close(mv['middle_billion'], middle)
            assert close(mv['homesteads_billion'], min(density * mv['middle_wet_fields_km2'] / 1e9, middle))
            assert close(mv['moved_billion'], mv['middle_billion'] - mv['homesteads_billion'])
            assert close(sum(mv['placed_billion'].values()) + mv['unplaced_billion'], mv['moved_billion'])
            assert mv['unplaced_billion'] == 0.0
            assert abs(c['people_unplaced_billion']) < 1e-12
            assert c['bands']['scattered_middle_band_km2'] == 0.0
            assert close(c['land']['settled_km2'], base['land']['settled_km2'])
            if mv['middle_wet_fields_km2'] == 0.0:
                assert mv['homesteads_billion'] == 0.0


def test_the_destination_ranking_follows_its_measures(product, basis):
    B, _ = basis
    spare = {k: 1.0 for k in m.DESTINATIONS}
    r = m.rank_destinations(B, spare)
    assert all(v['spare_land'] == 3.0 for v in r['ranks'].values())          # ties share the mean rank
    for kept, rk in product['moon']['destination_ranking'].items():
        order = rk['order']
        means = [rk['ranks'][k]['mean'] for k in order]
        assert means == sorted(means)
        assert set(order) == set(m.DESTINATIONS)
        assert product['moon']['cprime'][kept]['destination_order'] == order
    D = product['basis']['destinations']
    assert D['dry_30_45']['water_gw_per_billion_nothing_recovered'][0] < D['dry_60_90'][
        'water_gw_per_billion_nothing_recovered'][0]
    assert D['high_latitude_lowland']['lit_m2_per_person'] == [0.0, 0.0]


def test_cprime_band_check_at_the_costly_corners(product, basis):
    """Outside a kept half the homesteads' homes overfill the full middle-wet zone at the high-cost end when they sit
    beside the fields; with the remedy's poleward and twilit placements and their homes counted within the fields
    they work, every corner fits (C')."""
    B, d = basis
    order = product['moon']['cprime']['kept_50']['destination_order']
    c = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 1, middle=dict(order=order, h_end=1))
    bands = c['bands']
    assert close(bands['use_km2']['middle_wet_zone'],
                 bands['fields_on_middle_wet_km2'] + bands['homestead_land_on_middle_wet_km2'])
    assert close(bands['overfill_km2']['middle_wet_zone'], c['moved']['homestead_land_km2'], rel=1e-6)
    k30 = product['moon']['cprime']['kept_30']
    assert k30['buildable_corners'] == k30['corners']
    v = product['moon']['cprime']['kept_50']['variants']
    assert v['homes_within_the_fields']['buildable_corners'] == 32
    assert v['poleward_and_twilit']['buildable_corners'] < 32
    within = m.moon_corner(B, d['food'], 8, 4, 0.5, 'fusion', (0, 1, 0), 1, rule='twilit', dry_rule='poleward',
                           middle=dict(order=order, h_end=1, homes_on='within_fields'))
    share = within['bands']['homestead_share_of_middle_wet_fields']
    assert close(share, within['moved']['homestead_land_km2'] / within['bands']['fields_on_middle_wet_km2'])
    assert within['buildable'] and 0.0 < share < 0.01
    crews = B['crews']
    assert crews['km_to_45_degrees'] > crews['sky_boat_range_with_climb_km']
    assert crews['km_to_45_degrees'] < crews['regional_ship_route_km']


def test_land_rules_match_the_food_study(product):
    rules = product['assumptions']['land_rules']
    assert rules == product['checks']['land_rules']['food_study'] and product['checks']['land_rules']['equal']
    half, comp = (next(k for k in m.KEPT if m.kept_key(k) == key) for key in ('kept_50', 'wild_30_hybrid_20'))
    assert half != comp and half == 0.5 and comp != 0.5 and len({half, comp}) == 2
    with pytest.raises(ValueError):
        m.LandRule(0.5, 0.2)
    assert set(rules) == {m.kept_key(k) for k in m.KEPT}
    for k in m.KEPT:
        rule = rules[m.kept_key(k)]
        assert close(float(k), rule['closed']) and close(m.hybrid_share(k), rule['hybrid'])


def test_cprime_under_the_authors_compromise(product, basis):
    """30% wild with 20% hybrid closes the same half to farms and towns as a kept half, and its hybrid fifth is the
    author's parks and gardens, so the results for C' match the kept half's and the parks add 20% of the dry land; every
    corner fits with the remedy's poleward and twilit placements once the homesteads' homes count within their
    fields (C')."""
    B, d = basis
    cp = product['moon']['cprime']
    fits = {k: (cp[k]['buildable_corners'], *(cp[k]['variants'][v]['buildable_corners']
                                              for v in ('poleward_and_twilit', 'homes_within_the_fields',
                                                        'to_the_dry_highlands')))
            for k in ('kept_30', 'kept_50', 'wild_30_hybrid_20', 'wild_50_hybrid_20')}
    assert fits == {'kept_30': (32, 32, 32, 32), 'kept_50': (16, 16, 32, 16),
                    'wild_30_hybrid_20': (16, 16, 32, 16), 'wild_50_hybrid_20': (16, 16, 32, 16)}
    half, comp = cp['kept_50'], cp['wild_30_hybrid_20']
    for key in ('full', 'feasible'):
        a, b = half[key], comp[key]
        assert {k: v for k, v in a.items() if k != 'land'} == {k: v for k, v in b.items() if k != 'land'}
        assert ({k: v for k, v in a['land'].items() if not k.startswith('parks')}
                == {k: v for k, v in b['land'].items() if not k.startswith('parks')})
    dry = B['moon']['dry_land_km2']
    assert close(comp['full']['land']['parks_km2'][0], 0.2 * dry, rel=5e-3)
    assert half['full']['land']['parks_km2'] == [0.0, 0.0]


def test_the_dark_sky_holds_with_every_park_lit(product):
    """Under the compromise the towns, cities and sky towns light at most 17.2% of the dry land, so even with every
    park lit more than half keeps a dark sky, above the preferred half and the compromise's 30% (C')."""
    comp = product['moon']['cprime']['wild_30_hybrid_20']
    lit = max(comp['full']['light']['shielded']['share_of_dry_land_metropolis_at_earth_area'])
    parks = comp['full']['land']['parks_share_of_dry_land'][1]
    assert lit < 0.175 and 1.0 - lit - parks > 0.5


def test_homesteads_homes_stay_on_the_farm_half_under_a_hybrid_rule(product, basis):
    """Beside their fields the homesteads' homes count against the middle-wet zone's farm half, which the fields fill at
    the high-cost end, so they overfill it by their own land until they count within the fields (C')."""
    B, d = basis
    rule = next(k for k in m.KEPT if m.kept_key(k) == 'wild_30_hybrid_20')
    order = product['moon']['cprime']['wild_30_hybrid_20']['destination_order']
    c = m.moon_corner(B, d['food'], 8, 4, rule, 'fusion', (0, 1, 0), 1, rule='twilit', dry_rule='poleward',
                      middle=dict(order=order, h_end=1))
    homes = c['moved']['homestead_land_km2']
    bands = c['bands']
    assert homes > 0.0 and close(bands['homestead_land_on_middle_wet_km2'], homes)
    assert close(bands['overfill_km2']['middle_wet_zone'], homes, rel=1e-6)
    within = m.moon_corner(B, d['food'], 8, 4, rule, 'fusion', (0, 1, 0), 1, rule='twilit', dry_rule='poleward',
                           middle=dict(order=order, h_end=1, homes_on='within_fields'))
    assert within['buildable']


def test_middle_wet_fields_are_held_to_the_fills_room_alone(basis):
    """The cap on the middle-wet fields removes only the food product's rounding: at most 0.5% of its field area."""
    B, d = basis
    for kept in m.KEPT:
        for moon in m.MOON_BILLION:
            blk = m.food_block(d['food'], m.CHOSEN_MIX, moon, kept)['block']
            lb = B['land_bands']
            hb, zones = lb['highland_km2'], lb['zones_km2']
            room = (min(lb['dry_land_km2']['30-45'], zones['middle_wet'] + hb['30-45']) - hb['30-45']) * (1.0 - kept)
            for e in (0, 1):
                listed = max(blk['fields']['field_km2_by_zone']['middle_wet'][e], 0.0)
                held = m.middle_wet_fields(B, blk, e, kept)
                assert held <= room * (1 + 1e-12) and held <= listed, (m.kept_key(kept), moon, e)
                assert listed - held <= 5e-3 * listed + 1e-9, (m.kept_key(kept), moon, e)


def test_the_hybrid_part_is_parkland(basis):
    """The hybrid part holds the author's parks and gardens: grazing and the homesteads' homes keep off it, the fog
    band's grazing counts against the land outside the closed share, and the parks are 20% of the dry land."""
    B, d = basis
    zones = B['land_bands']['zones_km2']
    area = zones['fog_desert'] + zones['polar']
    dry = B['moon']['dry_land_km2']
    for kept in m.KEPT:
        c = m.moon_corner(B, d['food'], 8, 4, kept, 'fusion', (0, 0, 0), 1, mix='rain_fed_first')
        bands = c['bands']
        assert c['land']['grazing_km2'] > 0.0
        assert close(bands['use_share']['fog_desert_and_polar'],
                     bands['use_km2']['fog_desert_and_polar'] / (area * (1.0 - kept)))
        assert close(c['land']['parks_km2'], dry * m.hybrid_share(kept))
        assert 'homestead_land_on_hybrid_km2' not in bands


def test_wild_shares_count_against_the_land_outside_the_closed_share(basis):
    B, d = basis
    dry = B['moon']['dry_land_km2']
    for kept in m.KEPT:
        c = m.moon_corner(B, d['food'], 8, 4, kept, 'fusion', (1, 0, 1), 1)
        assert close(c['land']['wild_share_outside_kept'], 1.0 - c['land']['used_km2'] / (dry * (1.0 - kept)))
        s = c['sky']
        assert close(s['wild_sky_outside_kept_by_grid'], 1.0 - s['grid_share_of_moon'] / (1.0 - kept))
        assert close(s['towns_grid_share_of_their_sky'], s['grid_share_of_moon'] / (1.0 - kept))


def test_the_four_distributions_under_the_authors_compromise_and_the_comparison(product):
    """Under the compromise A and B fit as placed and C and D with the remedy's three placements; with half wild and
    20% hybrid B needs the remedy and D fits 10 of its 16 corners with it."""
    cases = product['moon']['cases']
    comp = {k: cases[k]['wild_30_hybrid_20'] for k in ('4_2', '6_3', '8_4', '8_2')}
    assert comp['4_2']['buildable_corners'] == comp['6_3']['buildable_corners'] == 16
    assert comp['8_4']['remedy']['buildable_corners'] == comp['8_2']['remedy']['buildable_corners'] == 16
    for k in comp:
        assert comp[k]['buildable_corners'] == cases[k]['kept_50']['buildable_corners']
    wide = {k: cases[k]['wild_50_hybrid_20'] for k in ('6_3', '8_2')}
    assert wide['6_3']['buildable_corners'] < 16 and wide['6_3']['remedy']['buildable_corners'] == 16
    assert wide['8_2']['remedy']['buildable_corners'] == 10


def test_dumps_keeps_value_lists_on_one_line():
    text = m.dumps({'a': [1.0, 2.0], 'b': {'c': [3.0]}})
    assert '[1.0, 2.0]' in text and '[3.0]' in text


# ---- writing -----------------------------------------------------------------------------------------------------

BANNED = (r'\bthin\b', r'\bregolith\b', r'\bfloaters?\b', r'\bleaks?\b', r'\bleakage\b')


def test_the_writing_keeps_the_projects_words():
    for name in ('README.md',):
        text = (m.HERE / name).read_text().lower()
        for pattern in BANNED:
            assert not re.search(pattern, text), (name, pattern)
        fenced = False
        for i, line in enumerate(text.splitlines(), 1):
            if line.startswith('```'):
                fenced = not fenced
                continue
            if fenced or line.startswith('|'):
                continue                       # prose lines only: tables and code keep their width
            assert len(line) <= 120, (name, i)
