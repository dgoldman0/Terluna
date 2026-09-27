"""Checks on the canopy model: leaf optics, light in the canopy, the leaf model, cycle integrals and stored results."""
import json
import math

import numpy as np
import pytest
from scipy.special import expn

from biosphere.canopy import fetch_inputs, leaf, leaf_optics as lo, model, radiation as rad

HAVE_PROSPECT = fetch_inputs.restore()['status'] == 'PASS'
STORED = model.RESULTS / 'canopy.json'


def test_plate_transmissivity_is_twice_the_third_exponential_integral():
    k = np.array([1e-4, 0.01, 0.3, 1.0, 4.0])
    assert lo.plate_transmissivity(k) == pytest.approx(2 * expn(3, k), rel=1e-10)
    assert lo.plate_transmissivity(np.array([0.0]))[0] == 1.0


def test_interface_transmissivity_limits():
    # isotropic light over the whole hemisphere and light near the normal (Fresnel) for n = 1.5
    assert lo.interface_transmissivity(90.0, np.array([1.5]))[0] == pytest.approx(0.908, abs=2e-3)
    assert lo.interface_transmissivity(1.0, np.array([1.5]))[0] == pytest.approx(1 - (0.5 / 2.5) ** 2, abs=1e-3)


@pytest.mark.skipif(not HAVE_PROSPECT, reason='PROSPECT-D coefficients not restored (fetch_inputs --download)')
def test_leaf_optics_match_the_independent_reference():
    # The default leaf, as the prosail package's PROSPECT-D gives it (agreement to 1e-15 on 2026-09-26).
    r, t = lo.optics()
    at = lambda a, nm: a[lo.WAVELENGTH_NM == nm][0]
    for nm, rr, tt in ((450, 0.041, 0.001), (550, 0.151, 0.150), (670, 0.036, 0.006), (700, 0.127, 0.135),
                       (730, 0.368, 0.395), (800, 0.443, 0.475)):
        assert at(r, nm) == pytest.approx(rr, abs=6e-4) and at(t, nm) == pytest.approx(tt, abs=6e-4)
    bare = lo.optics(lo.Leaf(chlorophyll=0, carotenoids=0, water=0, dry_matter=0))
    assert bare[0] + bare[1] == pytest.approx(np.ones(lo.WAVELENGTH_NM.size), abs=1e-12)
    pale = lo.optics(lo.Leaf(chlorophyll=20, carotenoids=4))
    red = lo.WAVELENGTH_NM == 670
    assert (pale[0] + pale[1])[red] > (r + t)[red]


def test_canopy_of_black_leaves_transmits_exactly():
    zero = np.zeros(3)
    canopy = rad.Canopy(lai=3.0)
    direct = rad.solve(canopy, zero, zero, zero, mu0=0.5)
    assert direct['to_soil'] == pytest.approx(np.full(3, math.exp(-0.5 * 3.0 / 0.5)), rel=1e-12)
    diffuse = rad.solve(canopy, zero, zero, zero)
    assert diffuse['to_soil'] == pytest.approx(np.full(3, 2 * expn(3, 1.5)), rel=1e-5)


@pytest.mark.parametrize('mu0', [None, 1.0, 0.3, 0.05])
def test_canopy_conserves_energy_and_counts_sunlit_leaves(mu0):
    rng = np.random.default_rng(1)
    rho, tau, soil = rng.uniform(0.02, 0.5, 6), rng.uniform(0.0, 0.45, 6), rng.uniform(0.0, 0.4, 6)
    canopy = rad.Canopy(lai=4.0, clumping=0.8)
    out = rad.solve(canopy, rho, tau, soil, mu0=mu0)
    assert out['reflected'] + out['leaf_absorbed'] + out['soil_absorbed'] == pytest.approx(np.ones(6), abs=1e-7)
    n, dl = canopy.layers()
    sunlit_lai = out['sunlit'].sum() * dl
    expected = 0.0 if mu0 is None else mu0 / rad.G * (1 - math.exp(-rad.G * 0.8 * 4.0 / mu0))
    assert sunlit_lai == pytest.approx(expected, rel=1e-12, abs=1e-15)


def test_thinner_layers_change_little():
    rng = np.random.default_rng(2)
    rho, tau, soil = rng.uniform(0.02, 0.5, 4), rng.uniform(0.0, 0.45, 4), rng.uniform(0.0, 0.4, 4)
    coarse = rad.solve(rad.Canopy(lai=5.0), rho, tau, soil, mu0=0.5)
    fine = rad.solve(rad.Canopy(lai=5.0, layer_lai=0.00625), rho, tau, soil, mu0=0.5)
    assert np.max(np.abs(coarse['leaf_absorbed'] - fine['leaf_absorbed'])) < 3e-3
    assert np.max(np.abs(coarse['reflected'] - fine['reflected'])) < 3e-3


def test_leaf_model_kinetics_and_limits():
    gamma, kc, ko = leaf.kinetics(25.0)
    # Bernacchi et al. (2001) at 25 C: 42.75 umol/mol, 404.9 umol/mol and 278.4 mmol/mol, near 100 kPa
    assert (gamma, kc, ko) == pytest.approx((4.275, 40.49, 27840.0), rel=0.01)
    t = leaf.Traits()
    ci = t.ci_ca * 41.0
    saturated = leaf.gross(1e9, 60.0, 25.0, 41.0, 21227.0)
    rubisco = 60.0 * (ci - gamma) / (ci + kc * (1 + 21227.0 / ko))
    electron = 100.0 * (ci - gamma) / (4 * ci + 8 * gamma)
    assert saturated == pytest.approx(min(rubisco, electron), rel=1e-6)
    light = np.array([0.0, 0.1, 10.0, 100.0, 1000.0])
    a = leaf.gross(light, 60.0, 25.0, 41.0, 21227.0)
    assert a[0] == 0.0 and np.all(np.diff(a) > 0)
    # initial slope: I2 = I (1 - f) / 2 electrons, four or more per CO2
    assert a[1] == pytest.approx(0.1 * 0.425 * (ci - gamma) / (4 * ci + 8 * gamma), rel=1e-3)


def test_cycle_mean_is_exact_for_a_rate_linear_in_the_sine_of_the_sun():
    heights = (90.0, 45.0, 10.0, 1.0)
    for lat in (0.0, 30.0, 60.0):
        values = [10.0 * math.sin(math.radians(h)) for h in heights]
        assert model.cycle_mean(heights, values, lat) == pytest.approx(10.0 * math.cos(math.radians(lat)) / math.pi,
                                                                       rel=1e-6)


@pytest.fixture(scope='module')
def product():
    if not STORED.is_file():
        pytest.skip('results/canopy.json has not been generated')
    return json.loads(STORED.read_text())


def test_stored_product(product):
    assert product['schema'] == model.SCHEMA
    for world in model.WORLDS:
        rows = product['by_sun_height'][world]
        gpp = [r['gpp_umol_m2_s'] for r in rows]
        assert np.all(np.diff(gpp) < 0)                                  # heights run from 90 down to 1
        for r in rows:
            assert 0 < r['absorbed_umol_m2_s'] < r['incident_par_umol_m2_s']
            assert 0 < r['canopy_albedo_par'] < 0.1
        lai = product['lai'][world]
        by_lai = [lai[k]['gpp_mean_umol_m2_s'] for k in sorted(lai, key=float)]
        assert np.all(np.diff(by_lai) > 0)
        cycle = product['summaries'][world]['0.0']['cycle']
        assert cycle['balance_closes']
        assert cycle['store_for_the_night_g_c_m2'] >= cycle['night_respiration_g_c_m2'] * 0.99
        assert product['variants']['extended_par'][world]['gpp_mean_umol_m2_s'] > \
            product['summaries'][world]['0.0']['gpp_mean_umol_m2_s']
    steps = product['decomposition']
    assert steps['moon']['0.0'] == pytest.approx(product['summaries']['moon']['0.0']['gpp_mean_umol_m2_s'], rel=1e-6)
    assert steps['earth']['0.0'] == pytest.approx(product['summaries']['earth']['0.0']['gpp_mean_umol_m2_s'], rel=1e-6)
