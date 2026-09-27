"""Checks on the surface-light product: the ground relation, bins and bands, day integrals and stored results."""
import json
import math

import numpy as np
import pytest

from shared.constants import AVOGADRO, PLANCK, SPEED_OF_LIGHT
from atmosphere.radiative_convective import shortwave as sw
from illumination.surface_light import model as m

STORED = m.RESULTS / 'surface_light.json'


def synthetic(layers=20, points=60, seed=3):
    """Layers top-down with Rayleigh-like scattering and patchy absorption, and their level radii."""
    rng = np.random.default_rng(seed)
    tau_sca = np.outer(np.linspace(0.01, 0.2, layers), np.linspace(3.0, 0.2, points))
    tau_abs = rng.uniform(0.0, 0.3, (layers, points)) * (rng.uniform(size=(layers, points)) < 0.4)
    radius = 1737.4e3 + np.linspace(400e3, 0.0, layers + 1)
    return tau_abs, tau_sca, radius


@pytest.mark.parametrize('mu0', [1.0, 0.4, math.sin(math.radians(2))])
def test_ground_relation_matches_direct_solutions(mu0):
    tau_abs, tau_sca, radius = synthetic()
    g = np.zeros_like(tau_sca)
    black = sw.column_fluxes(tau_abs, tau_sca, g, 0.0, mu0, radius)
    refl, trans = m.from_below(tau_abs, tau_sca)
    for albedo in (0.2, 0.8):
        direct, diffuse, up = sw.column_fluxes(tau_abs, tau_sca, g, albedo, mu0, radius)
        got = m.with_ground(black[0][-1], black[1][-1], black[2][0], refl, trans, albedo)
        for a, b in zip(got, (direct[-1], diffuse[-1], up[0])):
            assert a == pytest.approx(b, rel=1e-10, abs=1e-15)


def test_atmosphere_seen_from_below_conserves_energy_without_absorption():
    _, tau_sca, _ = synthetic()
    refl, trans = m.from_below(np.zeros_like(tau_sca), tau_sca)
    assert np.all((refl > 0) & (refl < 1))
    assert refl + trans == pytest.approx(np.ones_like(refl), abs=1e-6)


def test_bins_and_bands():
    assert m.CENTRES_NM[0] == 202.5 and m.CENTRES_NM[-1] == 999.5 and m.CENTRES_NM.size == 798
    assert m.N_UV + m.N_LBL == m.CENTRES_NM.size
    band = lambda name: m._band(*m.BANDS_NM[name])
    assert np.array_equal(band('par'), band('blue') | band('green') | band('red'))
    assert not np.any(band('blue') & band('green'))
    assert band('red_660').sum() == 10 and band('far_red_730').sum() == 10


def test_day_integral_is_exact_for_a_flux_linear_in_the_sine_of_the_sun():
    heights = (90.0, 45.0, 10.0, 1.0)
    flux = [100.0 * math.sin(math.radians(h)) for h in heights]
    for lat in (0.0, 30.0, 60.0):
        expected = 100.0 * math.cos(math.radians(lat)) * 86400.0 / math.pi
        assert m.day_integral(heights, flux, lat, 86400.0) == pytest.approx(expected, rel=1e-6)


@pytest.fixture(scope='module')
def product():
    if not STORED.is_file():
        pytest.skip('results/surface_light.json has not been generated')
    return json.loads(STORED.read_text())


def test_stored_product_reproduces_the_middle_atmosphere_uv_index(product):
    assert product['schema'] == m.SCHEMA
    check = product['checks']['earth_subsolar_uv_index']
    assert check['this_model'] == pytest.approx(check['middle_atmosphere'], rel=0.03)


def test_stored_light_behaves(product):
    summ = lambda case, h, a: product['cases'][case]['summaries'][h][a]
    for case in product['cases']:
        for h in product['cases'][case]['summaries']:
            par = [summ(case, h, a)['photons_umol_m2_s']['par'] for a in map(str, m.ALBEDOS)]
            share = [summ(case, h, a)['diffuse_share']['par'] for a in map(str, m.ALBEDOS)]
            assert np.all(np.diff(par) > 0) and np.all(np.diff(share) > 0)
        noon = [summ(case, str(h), '0.1')['photons_umol_m2_s']['par'] for h in m.SUN_DEG]
        assert np.all(np.diff(noon) < 0)
    for h in map(str, m.SUN_DEG):
        design, bare = summ('moon_1.2atm/design', h, '0.1'), summ('moon_1.2atm/unfiltered', h, '0.1')
        assert design['photons_umol_m2_s']['uv_b'] < 1e-9 < bare['photons_umol_m2_s']['uv_b']
        assert design['photons_umol_m2_s']['par'] < bare['photons_umol_m2_s']['par']


def test_reading_rule_reproduces_the_exact_summaries(product):
    """The stored black-ground spectra and binned R give the band photons of a bright ground within 1%."""
    lam = np.array(product['centres_nm'])
    to_umol = lam * 1e-9 / (PLANCK * SPEED_OF_LIGHT) * 1e6 / AVOGADRO
    for case, data in product['cases'].items():
        refl = np.array(product['columns'][data['column']]['reflectance_from_below'])
        for h in ('90.0', '30.0', '5.0'):
            spec = data['spectra'][h]
            direct, black = np.array(spec['direct_w_m2_nm']), np.array(spec['diffuse_black_w_m2_nm'])
            total = (direct + black) / (1 - 0.8 * refl) * to_umol
            exact = data['summaries'][h]['0.8']['photons_umol_m2_s']
            for name in ('blue', 'par', 'far_red', 'near_infrared'):
                assert total[m._band(*m.BANDS_NM[name])].sum() == pytest.approx(exact[name], rel=0.01)
