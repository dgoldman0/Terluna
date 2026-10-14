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


@pytest.fixture(scope='module')
def product():
    if not tw.OUT.exists():
        pytest.skip('twilight.json not built')
    return json.loads(tw.OUT.read_text())


def test_the_product_integrates_back_and_its_channels_keep_the_energy(product):
    assert product['schema'] == tw.SCHEMA
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


def test_the_pseudo_spherical_column_casts_no_light_below_the_horizon(product):
    for name in ('pseudo_spherical', 'flat', 'plasim'):
        below = product['global_mean_w_m2']['whole_spectrum'][name]['by_sun']['below_horizon']
        assert below['ground'] == 0.0 and below['air'] == 0.0
    sphere = product['global_mean_w_m2']['whole_spectrum']['sphere']['by_sun']['below_horizon']
    assert sphere['ground'] > 0.0 and sphere['air'] > 0.0


def test_the_cone_and_the_curved_light_make_up_the_geometry(product):
    total, cone, light = (product[k] for k in ('sphere_minus_pseudo_spherical_w_m2', 'cone_w_m2', 'curved_light_w_m2'))
    for key in ('ground', 'air'):
        assert cone[key] + light[key] == pytest.approx(total[key], abs=1e-6)
    assert cone['ground'] == 0.0 and cone['air'] > 0.0                  # the widening adds air above, not ground
    term = product['gcm_term']
    assert term['global_mean_w_m2']['air'] == pytest.approx(light['air'], rel=0.02)
    assert term['global_mean_w_m2']['ground'] == pytest.approx(light['ground'], rel=0.02, abs=0.01)
