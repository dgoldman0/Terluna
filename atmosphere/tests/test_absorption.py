"""Checks on the shield's protected radius and the exosphere's ion production."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from shared.constants import BOLTZMANN, MOON_GM, MOON_RADIUS
from atmosphere.loss_response import absorption as ab


def direct_partition(lam, lam_c, n=100001):
    """The share of the exobase Maxwellian found at radius r, by brute force over velocity space."""
    u = np.linspace(1e-6, 8.0, n)
    s = (lam / lam_c) ** 2 * (1 + (lam_c - lam) / u ** 2)     # sin^2 of the widest angle that still meets the exobase
    cone = 1 - np.sqrt(np.clip(1 - s, 0, None))
    w = 4 / math.sqrt(math.pi) * u ** 2 * np.exp(-u ** 2)
    return (np.trapezoid(np.where(u ** 2 < lam, w * cone, 0), u),
            np.trapezoid(np.where(u ** 2 >= lam, w * cone / 2, 0), u))


def isothermal_profile(t=250.0, r_base=2.0, r_top=3.0, n_base=1e16, points=4000):
    k = MOON_GM * ab.AIR_KG / (BOLTZMANN * t * MOON_RADIUS)
    r = np.linspace(r_base, r_top, points)
    n = n_base * np.exp(k / r - k / r_base)
    density = lambda x: n_base * math.exp(k / x - k / r_base)
    return dict(radius_R=r, temperature_K=np.full_like(r, t), pressure_Pa=n * BOLTZMANN * t), density, k


@pytest.mark.parametrize('lam_c', [6.0, 11.5])
@pytest.mark.parametrize('ratio', [1.0, 2.0, 10.0])
def test_partition_matches_direct_integral(lam_c, ratio):
    ballistic, escaping = ab.partition(np.array([lam_c / ratio]), lam_c)
    b, e = direct_partition(lam_c / ratio, lam_c)
    assert ballistic[0] == pytest.approx(b, rel=1e-3)
    assert escaping[0] == pytest.approx(e, rel=1e-3)


def test_grazing_column_matches_quadrature_and_chapman():
    profile, density, k = isothermal_profile()
    for b in (2.2, 2.6, 2.9):
        span = math.sqrt(3.0 ** 2 - b ** 2)
        exact = 2 * quad(lambda s: density(math.hypot(b, s)), 0, span, limit=200)[0] * MOON_RADIUS
        column = ab.thermosphere_columns(profile, [b])[0]
        assert column == pytest.approx(exact, rel=1e-3)
    # far below the top, the grazing column is n(b) sqrt(2 pi b H), up to a correction of order 1/lambda
    b = 2.2
    scale_height = b * MOON_RADIUS / (k / b)
    assert ab.thermosphere_columns(profile, [b])[0] == pytest.approx(
        density(b) * math.sqrt(2 * math.pi * b * MOON_RADIUS * scale_height), rel=0.1)


def test_thin_limb_deposit_counts_the_molecules_outside_the_shadow():
    profile, density, _ = isothermal_profile(n_base=1e10)
    weights = dict(xuv=1.0, lyman_alpha=0.0, schumann_runge=0.0)
    b, heat, _ = ab.limb_heating_curve(profile, weights)
    radius = 2.3
    outside = quad(lambda r: density(r) * 4 * math.pi * r * r * math.sqrt(1 - (radius / r) ** 2), radius, 3.0, limit=200)[0]
    expected = ab.CROSS_SECTIONS_M2['xuv'] * outside * MOON_RADIUS / math.pi
    assert float(np.interp(radius, b, heat)) == pytest.approx(expected, rel=0.01)
    assert heat[-1] == 0.0 and np.all(np.diff(heat) <= 0)


def test_ion_production_falls_with_the_shadow_and_ends_at_the_outer_radius():
    profile, _, _ = isothermal_profile()
    radii = np.array([0.0, 1.0, 3.0, 4.0, 6.0, 10.0, 12.0])
    production = ab.ion_production(profile, 'quiet', radii, ab.OUTER['third_hill'])
    assert np.all(np.diff(production) < 0)
    assert production[-1] == 0.0
    # with no shadow, the production is every ballistic molecule above the exobase times m J
    n_c, t_c, r_c = ab.exobase_state(profile)
    count = quad(lambda r: ab.exosphere_density(np.array([r]), n_c, t_c, r_c)[0][0] * 4 * math.pi * r * r,
                 r_c, ab.OUTER['third_hill'], limit=400)[0] * MOON_RADIUS ** 3
    assert production[0] == pytest.approx(ab.AIR_KG * ab.ionization_rate('quiet') * count, rel=1e-3)


def test_ionization_rate_is_of_the_tabulated_order():
    assert 3e-7 < ab.ionization_rate('quiet') < 1e-6
    assert ab.ionization_rate('solar_maximum') > ab.ionization_rate('quiet')


def test_september_aperture_is_the_reference():
    assert ab.aperture(3.0)['area_over_september'] == pytest.approx(1.0)
    assert ab.aperture(3.0)['aperture_radius_km'] == pytest.approx(3 * MOON_RADIUS / 1000 + 362.7, abs=0.5)


def test_design_case_needs_a_radius_just_under_the_exobase():
    row = ab.case('titania_stack', 'lte', 'quiet', 0.003, 10, ab.band_weights())
    radius = row['heating']['0.1']['radius_R']
    assert row['exobase_radius_R'] - 0.2 < radius < row['exobase_radius_R']
    assert row['tau_one_radius_R']['xuv'] < radius
    production = row['exosphere_ion_production_kg_s']['third_hill']
    assert production['heating_radius'] > production['6'] > production['10']
