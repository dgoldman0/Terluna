"""Checks of the twilight ledger's pure pieces (its channels, its sphere integral, its share among the GCM's layers)
and of its product, when present."""
import json

import numpy as np
import pytest

from climate.gcm import twilight as tw
from climate.gcm import plasim_shortwave as ps
from climate.tests.test_plasim_shortwave import SIGMA


def test_plasims_default_levels_are_the_runs():
    assert np.allclose(tw.plasim_sigma(), SIGMA, atol=1e-4)


def lines(n=4000, seed=1):
    """A band of strong and weak lines over a smooth continuum, on three layers."""
    rng = np.random.default_rng(seed)
    nu = np.linspace(0.0, 1.0, n)
    strength = np.zeros(n)
    for centre, height in zip(rng.uniform(0, 1, 40), 10.0 ** rng.uniform(-2, 3, 40)):
        strength += height / (1.0 + ((nu - centre) / 0.002) ** 2)
    absorption = np.outer([0.2, 0.3, 0.5], strength + 0.01)
    return absorption, np.full((3, n), 0.05), np.full(n, 1.0 / n), 500.0 + 10.0 * nu


def test_faithful_channels_keep_the_energy_and_the_absorbed_beam():
    absorption, scattering, energy, wavelength = lines()
    edges = np.array([500.0, 510.0 + 1e-9])
    ch = tw.faithful_channels(absorption, scattering, energy, wavelength, edges, 0.01)
    coarse = tw.faithful_channels(absorption, scattering, energy, wavelength, edges, 1e9)
    assert ch['energy'].sum() == pytest.approx(1.0) and coarse['energy'].size == 1
    for m in (0.1, 1.0, 10.0):
        exact = energy @ (1.0 - np.exp(-m * absorption.sum(0)))
        kept = ch['energy'] @ (1.0 - np.exp(-m * ch['absorption'].sum(0)))
        one = coarse['energy'] @ (1.0 - np.exp(-m * coarse['absorption'].sum(0)))
        assert kept == pytest.approx(exact, abs=0.01)
        assert one >= exact                            # an averaged optical depth absorbs more than its members
        if m <= 1.0:
            assert one > exact + 0.05


def test_faithful_channels_with_heights_keep_where_the_beam_is_lost():
    # Two kinds of member with the same whole column, one taking its light in the top layer and one in the bottom.
    n = 400
    absorption = np.zeros((3, n))
    absorption[2, :n // 2] = 4.0                                        # layers bottom-up
    absorption[0, n // 2:] = 4.0
    scattering, energy, wavelength = np.zeros((3, n)), np.full(n, 1.0 / n), np.linspace(1300.0, 1310.0, n)
    edges = np.array([1300.0, 1310.0 + 1e-9])
    flat = tw.faithful_channels(absorption, scattering, energy, wavelength, edges, 0.01)
    tall = tw.faithful_channels(absorption, scattering, energy, wavelength, edges, 0.01, heights=True)
    assert flat['energy'].size == 1 and tall['energy'].size >= 2      # the whole columns agree; the heights do not
    above = lambda a, e: e @ (1.0 - np.exp(-a[1:].sum(0)))              # the beam's loss above the bottom layer
    exact = above(absorption, energy)
    assert abs(above(tall['absorption'], tall['energy']) - exact) <= 0.01
    assert above(flat['absorption'], flat['energy']) - exact > 0.3       # one channel takes the light too high


def test_the_sphere_takes_the_disks_power_from_a_black_ground():
    a = np.radians(np.linspace(-90.0, 90.0, 721))
    assert tw.sphere_share(np.maximum(np.sin(a), 0.0), a) == pytest.approx(1.0, rel=1e-5)
    parts = [tw.class_share(np.maximum(np.sin(a), 0.0), a, lo, hi) for _, lo, hi in tw.SUN_CLASSES]
    assert sum(parts) == pytest.approx(1.0, rel=1e-5) and parts[0] == 0.0


def test_the_shells_share_out_among_the_gcm_layers_by_pressure():
    p = 1.0e5 * np.exp(-np.linspace(0.0, 12.0, 101))
    air = np.ones((100, 3))
    layers = tw.onto_gcm_layers(air, p)
    assert layers.sum(0) == pytest.approx(air.sum(0))
    half = np.r_[0.0, ps.half_levels(tw.plasim_sigma())]
    s = p / p[0]
    inside = (s[:-1] <= half[-1]) & (s[1:] >= half[-2])               # the shells within the lowest layer
    assert layers[-1, 0] == pytest.approx(inside.sum(), abs=1.0)


def test_the_slant_path_follows_the_ray_through_the_shells():
    R, H = 1737.4e3, 50e3
    edges = R + H * np.linspace(0.0, 10.0, 101)
    density = np.exp(-0.5 * (edges[1:] + edges[:-1] - 2 * R) / H)
    elevations = np.radians([90.0, 30.0, 0.0])
    m = tw.air_mass(edges, density, elevations)
    s = np.linspace(0.0, 6e6, 600_001)
    for a, got in zip(elevations, m):
        shell = np.searchsorted(edges, np.sqrt(R * R + s * s + 2 * R * s * np.sin(a))) - 1
        along = np.where(shell < 100, density[np.clip(shell, 0, 99)], 0.0)
        assert got == pytest.approx(np.trapezoid(along, s) / (density @ np.diff(edges)), rel=1e-3)
    assert m[0] == pytest.approx(1.0) and m[1] < 1.9                  # the curve shortens even a 30 degree path
    assert m[2] == pytest.approx(np.sqrt(np.pi * R / (2 * H)), rel=0.01)


def test_plasims_layers_add_up_to_its_air():
    from climate.tests.test_plasim_shortwave import moist_column
    net = ps.net_flux(moist_column(), 0.4, 'scale', 1.356)
    s = ps.clear_sky(moist_column(), 0.4, 'scale', 1.356)
    assert np.sum(net[:-1] - net[1:]) == pytest.approx(s['air']) and net[-1] == pytest.approx(s['ground'])


def test_the_301_k_column_is_built_between_the_stored_cases_that_bracket_it():
    from illumination.surface_light import model as surface
    assert tw.design_column() is surface.COLUMNS[tw.WORLD]
    assert tw.design_column(surface.DESIGN_TS) is surface.COLUMNS[tw.WORLD]
    column = tw.design_column(301.0)
    assert column.name == 'moon_1.2atm_ts301' and column.surface_temperature_k == 301.0
    (lo, wlo), (hi, whi) = column.stored
    assert (lo, hi) == ('moon_1.2atm_titania_stack_ts298', 'moon_1.2atm_titania_stack_ts308')
    assert (wlo, whi) == pytest.approx((0.7, 0.3))
    # The design column's own weights come out of the same rule.
    again = tw.design_column(294.9000001)
    assert [w for _, w in again.stored] == pytest.approx([w for _, w in surface.COLUMNS[tw.WORLD].stored], abs=1e-6)
    with pytest.raises(ValueError):
        tw.design_column(312.0)
    assert tw.product_path(column).name == 'twilight_301k.json' and tw.product_path(tw.design_column()) == tw.OUT
    assert tw.registered(column) == column.name and surface.COLUMNS[column.name] == column


LEDGERS = [tw.OUT, tw.RESULTS / 'twilight_301k.json']


@pytest.fixture(scope='module', params=LEDGERS, ids=[p.name for p in LEDGERS])
def product(request):
    if not request.param.exists():
        pytest.skip(f'{request.param.name} not built')
    return json.loads(request.param.read_text())


def test_the_301_k_ledger_records_its_column_and_its_solved_orders():
    path = tw.RESULTS / 'twilight_301k.json'
    if not path.exists():
        pytest.skip('twilight_301k.json not built')
    product = json.loads(path.read_text())
    assert product['schema'] == tw.SCHEMA and product['column_name'] == 'moon_1.2atm_ts301'
    source = product['column']['source']
    assert source['surface_temperature_k'] == 301.0
    assert source['stored_cases'] == {'moon_1.2atm_titania_stack_ts298': 0.7, 'moon_1.2atm_titania_stack_ts308': 0.3}
    channels = product['channels']
    assert channels['orders_converged'] is True                      # the sky model's own rule (solved_sky.py)
    for row in channels['continued_orders']:
        assert row['final_order'] > channels['order_cap']
        assert max(row['ground_increment_fraction'], row['maximum_increment_fraction']) < channels['increment_tolerance']


def test_the_product_integrates_back_and_its_channels_keep_the_energy(product):
    assert product['schema'] in (tw.SCHEMA, 'terluna.climate.gcm-twilight/1')
    by = product['by_elevation']
    a = np.radians(by['elevation_deg'])
    flux = product['sunlight']['covered_w_m2']
    for name in ('sphere', 'sphere_shallow', 'pseudo_spherical', 'flat', 'plasim'):
        whole = product['global_mean_w_m2']['whole_spectrum'][name]
        ground = flux / 4.0 * tw.sphere_share(np.array(by[name]['ground']), a)
        air = flux / 4.0 * tw.sphere_share(np.sum(by[name]['layers'], 0), a)
        assert ground == pytest.approx(whole['ground'], rel=1e-4) and air == pytest.approx(whole['air'], rel=1e-4)
        assert sum(x['air'] for x in whole['by_sun'].values()) == pytest.approx(whole['air'], rel=1e-4)
    for row in product['reduction_check']:
        for key in ('air', 'ground'):
            assert row['channels'][key] == pytest.approx(row['fine'][key], rel=0.01)
        # and layer by layer, within 1% of the band's absorbed sunlight
        assert np.allclose(row['channels']['air_by_layer'], row['fine']['air_by_layer'], rtol=0.0,
                           atol=0.01 * row['fine']['air'])


def test_the_pseudo_spherical_column_casts_no_light_below_the_horizon(product):
    for name in ('pseudo_spherical', 'flat', 'plasim'):
        below = product['global_mean_w_m2']['whole_spectrum'][name]['by_sun']['below_horizon']
        assert below['ground'] == 0.0 and below['air'] == 0.0
    sphere = product['global_mean_w_m2']['whole_spectrum']['sphere']['by_sun']['below_horizon']
    assert sphere['ground'] > 0.0 and sphere['air'] > 0.0


def test_the_cone_and_the_curved_light_make_up_the_geometry(product):
    total, cone, light = (product[k] for k in ('sphere_minus_pseudo_spherical_w_m2', 'cone_w_m2', 'curved_light_w_m2'))
    for key in ('ground', 'air'):
        assert cone[key] + light[key] == pytest.approx(total[key], rel=1e-5, abs=1e-6)   # six figures kept
    assert cone['ground'] == 0.0 and cone['air'] > 0.0                  # the widening adds air above, not ground
    term = product['gcm_term']
    assert term['global_mean_w_m2']['air'] == pytest.approx(light['air'], rel=0.02)
    assert term['global_mean_w_m2']['ground'] == pytest.approx(light['ground'], rel=0.02, abs=0.01)
