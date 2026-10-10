"""Phosphorus budgets close, the resources check reproduces the published case, and the routes conserve flux."""
import io
import json
import math

import numpy as np
import pytest

from research.studies.sky_ecology import phosphorus as ph
from research.studies.sky_ecology import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())['phosphorus']


@pytest.fixture(scope='module')
def food_web():
    return json.loads(m.OUT.read_text())['food_web']


@pytest.fixture(scope='module')
def winds():
    data, source = m.cross_branch_json('climate/results/gcm/zonal_winds_A28_dim5_moon.json')
    if data is None:
        pytest.skip(f'zonal-wind product unavailable ({source})')
    return data


def test_stocks_sum_by_tissue(product):
    for s in product['stocks'].values():
        for i in (0, 1):
            parts = sum(s[k][i] for k in ('green_g_m2', 'envelope_g_m2', 'tendon_g_m2', 'reserve_g_m2', 'tenants_g_m2',
                                          'pools_g_m2'))
            assert math.isclose(parts, s['total_g_m2'][i])
        assert s['total_g_m2'][1] < s['redfield_g_m2']       # one Redfield ratio overstates every tissue mix


def test_throughput_in_both_states(product):
    for key, f in product['throughput'].items():
        for i in (0, 1):
            assert math.isclose(f['growing']['total_g_m2_year'][i], f['renewal_g_m2_year'][i]+f['growth_g_m2_year'][i])
            assert math.isclose(f['mature']['total_g_m2_year'][i], f['renewal_g_m2_year'][i]+f['replacement_g_m2_year'][i])
            assert math.isclose(f['replacement_g_m2_year'][i],
                                product['stocks'][key]['total_g_m2'][i]*ph.HAZARD_PER_YEAR[i])
        for state in ph.STATES:
            assert f[state]['share_of_redfield'][1] < .2


def test_loops_close_in_both_states(product):
    for l in product['loops']:
        assert math.isclose(l['losses_g_m2_year'], l['litter_g_m2_year']+l['falls_g_m2_year']
                            + l['harvest_not_returned_g_m2_year']+l['tenants_leaving_g_m2_year'])
        rest = l['renewal_g_m2_year']-l['harvest_g_m2_year']-l['tenants_leaving_g_m2_year']
        assert math.isclose(rest, l['recycled_g_m2_year']+l['litter_g_m2_year'])
        assert math.isclose(l['need_g_m2_year'], l['losses_g_m2_year']+l['stock_change_g_m2_year'])
        assert math.isclose(l['frass_g_m2_year']+l['grazers_built_g_m2_year'], l['grazed_g_m2_year'])
        if l['state'] == 'mature':
            assert l['stock_change_g_m2_year'] == 0. and l['new_bodies_g_m2_year'] == l['falls_g_m2_year']
        else:
            assert math.isclose(l['stock_change_g_m2_year'], l['new_bodies_g_m2_year']-l['falls_g_m2_year'])
        assert math.isclose(l['recycled_share'], 1-(l['need_g_m2_year']-l['harvest_not_returned_g_m2_year'])
                            / l['throughput_g_m2_year'])
        if l['retention'] is not None:
            assert math.isclose(l['retained_share_of_rest'], l['retention'])
        assert l['net_need_g_m2_year'] <= l['need_g_m2_year']


def test_organism_alone_keeps_least_and_growth_needs_more(product):
    for body in ('sky_reef_60m', 'round_giant_1km', 'lesser_100m'):
        for edge in ('low', 'high'):
            by = {(l['state'], l['retention_case']): l for l in product['loops'] if l['body'] == body and l['edge'] == edge}
            for state in ph.STATES:
                r = {c: by[(state, c)]['recycled_share'] for c in ph.RETENTION_CASES}
                assert r['organism_only'] < r['earth_forest'] < r['designed_pools'] < .999
            for c in ph.RETENTION_CASES:
                g, mt = by[('growing', c)], by[('mature', c)]
                assert math.isclose(g['need_g_m2_year']-mt['need_g_m2_year'],
                                    g['new_bodies_g_m2_year']-mt['falls_g_m2_year'])


def test_sky_need_weights_the_classes(product):
    for row in product['sky_need']:
        need = sum(ph.CLASS_SHARES[cls]*next(l['net_need_g_m2_year'] for l in product['loops'] if l['body'] == b
                                             and l['state'] == row['state'] and l['retention_case'] == row['retention_case']
                                             and l['edge'] == row['edge']) for cls, b in ph.CLASS_BODIES.items())
        assert math.isclose(row['need_t_year_at_1pct'], need*.01*ph.MOON_AREA_M2/1e6)
    assert math.isclose(sum(ph.CLASS_SHARES.values()), 1.)


def test_grazers_cannot_build_more_p_than_they_eat(product):
    reef = next(l for l in product['loops'] if l['body'] == 'sky_reef_60m')
    prod = dict(green_tissue_dry_kg=1., structure_dry_kg=1., renewal_biomass_c=1., standing_host_c_kg_m2=1.,
                living_dry_kg_m2=1., structural_dry_kg_m2=1.)
    stock = dict(total_g_m2=[1., 1.])
    with pytest.raises(ValueError):
        ph.sky_loop(prod, stock, .1, grazer_p=1.)
    assert reef['grazers_built_g_m2_year'] <= reef['grazed_g_m2_year']


def test_consumer_p_balance(product, food_web):
    for key in ('sky_reef_60m', 'eaten_base_10m'):
        for edge, c in product['consumers'][key].items():
            chain = food_web['chains'][key][edge]
            assert c['p_limited_production_c'] <= chain['level2']['production_c']+1e-15
            assert c['buildable_p_g_m2_year'] <= c['food_p_g_m2_year']
            assert math.isclose(c['production_p_g_m2_year'], min(c['buildable_p_g_m2_year'],
                                                                 c['carbon_limited_need_p_g_m2_year']))
            if c['subsidy_for_carbon_limit_g_m2_year'] == 0:
                assert c['p_limited_share_of_carbon_limit'] == 1.


def test_resources_budget_formula():
    b = ph.resources_budget(1e6, .05, .999, .99, 100.)
    assert math.isclose(b['export_kg'], (.95*.001+.05*.01)*1e6)
    assert math.isclose(b['additional_return_needed_kg'], b['export_kg']-100.)
    assert ph.resources_budget(1e6, .05, .999, .99, 1e9)['additional_return_needed_kg'] == 0.
    for s in (.3, .9, .99):
        r = ph.equivalent_recycling(s)
        assert math.isclose((1-ph.RESOURCES_HARVEST)*(1-r), 1-s)


def test_resources_check_reproduces_the_published_case(product):
    c = product['resources_check']
    assert c['reproduces']
    assert math.isclose(c['reproduced']['additional_return_needed_kg'], 10_725_773.750154935, rel_tol=1e-9)
    assert math.isclose(c['redfield_throughput_g_m2_year'], 1000/(106*12.011/30.974), rel_tol=1e-9)
    loops = {(l['state'], l['retention_case'], l['edge']): l for l in product['loops'] if l['body'] == 'sky_reef_60m'}
    for case in c['cases']:
        loop = loops[(case['state'], case['retention_case'], case['edge'])]
        assert math.isclose(case['p_per_kg_c_g'], loop['throughput_g_m2_year']/loop['biomass_c_kg_m2_year'])
        natural = (1-ph.RESOURCES_HARVEST)*(1-case['natural_recycling'])*case['throughput_kg']
        assert math.isclose(natural, (1-case['recycled_share'])*case['throughput_kg'])
        assert case['throughput_share_of_redfield'] < .2


def test_resources_check_against_the_bound_product():
    data, source = m.cross_branch_json('resources/results/phosphorus.json')
    if data is None:
        pytest.skip(f'resources product unavailable ({source})')
    row = next(r for r in data['budget_grid'] if r['food_return'] == .99 and r['natural_recycling'] == .999
               and r['upward_from_surface_p_kg_yr'] < 3e6)
    b = ph.resources_budget(data['central']['throughput_p_kg_yr'], .05, .999, .99, row['upward_from_surface_p_kg_yr'])
    assert math.isclose(b['additional_return_needed_kg'], row['additional_return_needed_p_kg_yr'], rel_tol=1e-9)


def test_occupancy_uses_exact_band_areas_and_deposition_conserves_the_flux(product):
    for h, occ in product['occupancy'].items():
        area = ph.band_areas(occ['lat_edges_deg'])
        assert np.allclose(occ['area_share'], area, rtol=1e-12)
        assert math.isclose(area.sum(), 1., rel_tol=1e-12)
        assert np.allclose(np.array(occ['concentration'])*area, occ['time_share'], rtol=1e-12)
        assert abs(sum(occ['time_share'])-1) < 2e-3
        rows = product['deposition'][h]['per_kt']
        total = sum(r['g_m2_year']*ph.MOON_AREA_M2*a for r, a in zip(rows, area))
        assert math.isclose(total, 1e9*sum(occ['time_share']), rel_tol=1e-9)


def test_band_areas_are_exact():
    assert np.allclose(ph.band_areas([-90., 0., 90.]), [.5, .5])
    assert math.isclose(ph.band_areas([80., 90.])[0], (1-math.sin(math.radians(80.)))/2)


def _atlas(water):
    lat = np.linspace(89.875, -89.875, 720)
    lon = np.linspace(.125, 359.875, 1440)
    return dict(water_label=water.astype(int), lat_deg=lat, lon_deg=lon)


def test_land_share_on_synthetic_atlases():
    lat = np.linspace(89.875, -89.875, 720)
    lon = np.linspace(.125, 359.875, 1440)
    lons = np.array([[10., 200., 359.99], [90., 270., 0.01]])
    lats = np.array([[5., -40., 60.], [-5., 40., -60.]])
    w = np.ones_like(lons)
    assert ph.land_share_along(lons, lats, w, _atlas(np.ones((720, 1440))))['all'] == 0.
    assert ph.land_share_along(lons, lats, w, _atlas(np.zeros((720, 1440))))['all'] == 1.
    north_land = _atlas(np.repeat((lat < 0)[:, None], 1440, axis=1))          # land north of the equator
    out = ph.land_share_along(lons, lats, w, north_land)
    assert math.isclose(out['all'], .5) and out['by_record'] == [2/3, 1/3]
    east_land = _atlas(np.repeat((lon >= 180.)[None, :], 720, axis=0))        # land on longitudes 0-180
    out = ph.land_share_along(lons, lats, w, east_land)
    assert out['by_record'] == [1/3, 2/3]                         # 359.99 and 0.01 take the last and first columns
    sea = ph.band_sea_share(north_land, [-90., 0., 90.])
    assert np.allclose(sea, [1., 0.])
    occ = dict(time_share=[.25, .75])
    assert math.isclose(ph.band_land_share(occ, sea), .75)


def test_subsolar_longitude_follows_the_wind_codes_frame(winds):
    sun = winds['sun']
    # The wind code centres output n at (n + 1/2) output intervals; its track steps west by one interval's turn.
    assert math.isclose(sun['track_step_deg'], 360.*sun['output_interval_s']/(sun['synodic_month_days']*86400.),
                        rel_tol=1e-5)
    assert math.isclose(ph.subsolar_longitude(0., sun), (sun['track_start_lon_deg']+sun['track_step_deg']/2) % 360.)
    half = sun['output_interval_s']/2/86400.
    assert math.isclose(ph.subsolar_longitude(half, sun), sun['track_start_lon_deg'] % 360., abs_tol=1e-4)
    assert math.isclose((ph.subsolar_longitude(1., sun)-ph.subsolar_longitude(0., sun)) % 360.,
                        360.-360./sun['synodic_month_days'])


def test_route_positions_add_the_hour_to_the_sun(winds):
    raw, _ = m.cross_branch_bytes('climate/results/gcm/zonal_winds_A28_dim5_moon.npz')
    if raw is None:
        pytest.skip('zonal-wind positions unavailable')
    npz = np.load(io.BytesIO(raw))
    t = winds['trajectories']
    lon, lat, w = ph.route_positions(npz, winds, 10., first_record=20)
    per = len(t['release_lat_deg'])*len(t['release_hour_deg'])
    k = t['heights_km'].index(10.)
    raw_pos = np.asarray(npz['long_run_positions'], float)[20, k*per, :]/100.
    expect = (ph.subsolar_longitude(20*t['position_every_days'], winds['sun'])+raw_pos[0]) % 360.
    assert math.isclose(lon[0, 0], expect) and math.isclose(lat[0, 0], raw_pos[1])
    assert lon.min() >= 0 and lon.max() < 360 and abs(lat).max() <= 90
    assert lon.shape == lat.shape == w.shape and lon.shape[1] == per


def test_land_share_used_is_the_band_weighted_one(product):
    used = product['land_share_used']
    if product['band_sea_share'] is None:
        assert 'absent' in used['basis']
        return
    assert math.isclose(used['value'], ph.band_land_share(product['occupancy']['10.0'], product['band_sea_share']))
    snap = product['land_share_of_routes']['10.0']
    assert snap['snapshot_range'][0] < snap['snapshots'] < snap['snapshot_range'][1]


def test_falls(product):
    stock = dict(sky_reef_2km='sky_reef_60m', round_giant_1km='round_giant_1km', lesser_150m='lesser_150m')
    for f in product['falls']:
        s = product['stocks'][stock[f['body']]]['total_g_m2']
        assert math.isclose(f['p_t'][0], s[0]*f['area_m2']/1e6) and math.isclose(f['p_t'][1], s[1]*f['area_m2']/1e6)
        assert math.isclose(f['g_m2_on_ground'][0], s[0]/ph.DRAPE_FACTOR[1])
        assert f['g_m2_on_ground'][0] < f['g_m2_on_ground'][1]
        assert f['years_of_land_loss'][0] < f['years_of_land_loss'][1]
        assert all(math.isclose(n, f['bodies_at_cover']*h) for n, h in zip(f['falls_per_year'], ph.HAZARD_PER_YEAR))
        assert f['falls_per_year_at_growth_rate'][0] > f['falls_per_year'][0]
    assert all(x > 0 for x in product['land_loss_g_m2_year'])
    assert .05 < product['land_loss_g_m2_year'][0] < .07 and .17 < product['land_loss_g_m2_year'][1] < .19


def test_sea_return_arithmetic(product):
    s = product['sea_to_land']
    assert math.isclose(s['earth_seabirds_share_of_loss'][0], 99e6/4.16e9, rel_tol=1e-9)
    assert math.isclose(s['earth_seabirds_share_of_dissolved'][1], 99e6/s['dissolved_loss_kg_year'][0], rel_tol=1e-9)
    n = ph.flyers_for_return(1000., 50., .5, .8)
    assert math.isclose(n*50*.5*.8, 1000.)
    with pytest.raises(ValueError):
        ph.flyers_for_return(1000., 0., .5)
