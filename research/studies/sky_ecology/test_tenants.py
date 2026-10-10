"""Light geometry, real estate, capacity, domatia and contact of the tenants component."""
import json
import math

import numpy as np
import pytest

from research.studies.sky_ecology import run as m
from research.studies.sky_ecology import tenants as tn


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())['tenants']


@pytest.fixture(scope='module')
def surface():
    return json.loads((m.ROOT/'illumination/surface_light/results/surface_light.json').read_text())


def test_layer_reflectance_limits(surface):
    r = np.array(surface['columns']['moon_1.2atm']['reflectance_from_below'])
    assert np.allclose(tn.layer_reflectance(r, 1.), r)
    assert np.allclose(tn.layer_reflectance(r, 0.), 0.)
    assert np.all(tn.layer_reflectance(r, .2) <= tn.layer_reflectance(r, .4))


def test_underside_factor_limits():
    assert math.isclose(tn.underside_factor(.1, 0.), .1)
    assert math.isclose(tn.underside_factor(0., .3), .3)
    assert tn.underside_factor(.1, .3) < .1+.3


def test_two_stream_matches_the_product_and_bounds_its_absorption(surface, product):
    lam = np.array(surface['centres_nm'])
    col = surface['columns']['moon_1.2atm']
    diff = np.array(col['reflectance_from_below'])-tn.two_stream_column(col['rayleigh_optical_depth_550nm'], lam)
    par = (lam >= 400) & (lam <= 700)
    blue_green = (lam >= 400) & (lam <= 560)
    assert diff[par].max() < .005                         # absorption only lowers the reflectance
    assert np.median(np.abs(diff[blue_green])) < .005
    chk = product['light']['check']
    assert math.isclose(chk['max_abs_difference_400_560_nm'], np.abs(diff[blue_green]).max())
    assert math.isclose(chk['share_of_400_560_nm_bins_within_0_01'], (np.abs(diff[blue_green]) <= .01).mean())


def test_height_factor_is_two_stream_adding():
    for r in (0., .1, .3):
        for a in (0., .1, .5):
            f = tn.height_factor(r, a)
            ground = 1./f                                       # what reaches the ground of a unit arriving at height
            assert math.isclose(ground, (1-r)/(1-a*r))
            assert f >= 1.
    assert math.isclose(tn.height_factor(0., .3), 1.)


def test_light_at_height_exceeds_the_ground(product):
    for z in product['light']['cases']:
        assert 1. < z['height_gain'] < 1.15
        assert math.isclose(z['level_top_par_umol_m2_s'], z['ground_par_umol_m2_s']*z['height_gain'])


def test_zone_areas():
    s = tn.sphere_zone_areas()
    assert all(math.isclose(v, 1.) for v in s.values())
    r = tn.raft_zone_areas()
    assert math.isclose(sum(r.values()), 4*math.pi/(2*math.sqrt(3)))


def test_lone_sphere_view_factors_are_analytic():
    bands = [dict(sun_deg=a, direct=0., diffuse=100., up=0.) for a in (0., 90.)]
    z = tn.zone_light(bands, raft=False, n_dir=4000)
    for name, (lo, hi) in tn.ZONES.items():
        # mean of (1 + cos b)/2 over the zone's area, weighted by sin b
        b = np.radians(np.linspace(lo, hi, 2001))
        w = np.sin(b)
        expect = 100.*((1+np.cos(b))/2*w).sum()/w.sum()
        assert abs(z['zones'][name]['sky']-expect) < 1., name


def test_neighbours_only_take_light_away(product):
    for albedo in tn.ALBEDOS:
        lone = next(z for z in product['light']['cases'] if not z['raft'] and z['albedo'] == albedo)['zones']
        raft = next(z for z in product['light']['cases'] if z['raft'] and z['albedo'] == albedo)['zones']
        for k in lone:
            assert raft[k]['par_umol_m2_s'] <= lone[k]['par_umol_m2_s']+1e-9, (albedo, k)
        assert lone['under_cap']['share_of_level_top'] < lone['top_cap']['share_of_level_top']


def test_pocket_volume_by_integration():
    r = 1.
    rng = np.random.default_rng(1)
    pts = rng.uniform([-1., -1., 0.], [1., 1., 1.], size=(400_000, 3))*[math.sqrt(3), 1.5, 1.]
    # One cell of the hexagonal lattice: two triangles of side 2r; test the triangle with corners
    # (0,0), (2,0), (1, sqrt3) and its upper half-space prism outside the three spheres.
    tri = pts.copy()
    tri[:, 0] = rng.uniform(0, 2, len(pts))
    tri[:, 1] = rng.uniform(0, math.sqrt(3), len(pts))
    inside = (tri[:, 1] <= math.sqrt(3)*tri[:, 0]) & (tri[:, 1] <= math.sqrt(3)*(2-tri[:, 0]))
    centres = np.array([[0, 0, 0], [2, 0, 0], [1, math.sqrt(3), 0]])
    out = np.ones(len(tri), bool)
    for c in centres:
        out &= ((tri-c)**2).sum(axis=1) > r*r
    box = 2*math.sqrt(3)*1.
    vol = box*(inside & out).mean()
    assert abs(vol-(math.sqrt(3)-math.pi/3)) < .01
    assert math.isclose(tn.raft_pocket_volume_m3_per_m2(30.), (math.sqrt(3)-math.pi/3)*30/math.sqrt(3))


def test_capacity_rows_follow_the_product(product):
    aero = json.loads((m.ROOT/'research/studies/aerophytes/results/aerophytes.json').read_text())
    rows = {r['diameter_m']: r for r in aero['size_limits']['10_km']['rows']['round']}
    for c in product['capacity']['round']:
        r = rows[c['diameter_m']]
        gross = r['supported_kg_m2']/r['gas_fraction']
        assert math.isclose(c['steady_tenants_kg_m2'], (.7-r['gas_fraction'])*gross)
        assert math.isclose(c['sudden_trim_kg_m2'], r['water_swing_kg_m2_at_spare'])
        assert c['century_budget_tenants_kg_m2'][0] <= c['century_budget_tenants_kg_m2'][1] <= c['steady_tenants_kg_m2']+1e-9
    big = {c['diameter_m']: c for c in product['capacity']['round_strong']}
    assert big[1000.]['sudden_trim_t'] > 1000 > big[100.]['sudden_trim_t']


def test_domatia_carry_more_than_green_modules(product):
    green = {c['diameter_m']: c for c in product['capacity']['colony_module']}
    for d in product['domatia']:
        assert d['mass_gate']
        assert d['steady_tenants_kg_m2'] > green[d['diameter_m']]['steady_tenants_kg_m2']
        assert d['sudden_trim_kg_m2'] > green[d['diameter_m']]['sudden_trim_kg_m2']
        assert d['carbon_balance_c_kg_m2_year'] < 0


def test_spores_from_the_bodies_height(product):
    first = json.loads((m.ROOT/'biosphere/ecology/results/first_screen.json').read_text())['settling']
    out = product['contact']['fall_days_from_body_height']
    for name, p in first['particles'].items():
        assert math.isclose(out['days'][name], out['start_km']*1e3/p['fall_m_s']['moon']/86400.)


def test_encounter_rate():
    rate = tn.encounter_rate_per_day(100., .01, .5, 50.)
    density = .01/(math.pi*50**2)
    assert math.isclose(rate, density*.5*2*150*86400)
    assert math.isclose(tn.encounter_rate_per_day(100., .02, .5, 50.), 2*rate)


def test_landing_counts(product):
    for kind in ('round', 'round_strong', 'colony_module'):
        for row in product['landing'][kind]:
            assert 0 < row['at_once_by_trim'] < row['at_once_dropping_water']
            assert row['roosting_steadily'] > 0
