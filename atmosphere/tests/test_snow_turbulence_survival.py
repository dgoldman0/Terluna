"""Survival by load equivalence: the settling scale factor and the turbulent increment."""
import math

import pytest

from atmosphere.snow import flakes
from atmosphere.snow import turbulence_survival as survival


def test_settling_only_equivalence_is_the_gravity_ratio_to_the_inverse_mass_exponent():
    moon, _ = flakes.moon_column()
    earth = flakes.earth_column()
    for law, (_, b) in survival.MASS_LAWS.items():
        span = survival.equivalent_span(moon, earth, .04, 1e-4, 1e-4, law, 'weight_n')
        assert span / .04 == pytest.approx((earth.gravity / moon.gravity) ** (1 / b), rel=1e-9)


def test_velocity_increment_joins_the_viscous_and_inertial_ranges():
    nu, eps = 2e-5, 1e-3
    eta = (nu ** 3 / eps) ** .25
    below = survival.velocity_increment(eta * .999999, eps, nu)
    above = survival.velocity_increment(eta * 1.000001, eps, nu)
    assert below == pytest.approx(above, rel=1e-5)
    assert survival.velocity_increment(.1, eps, nu) == pytest.approx((eps * .1) ** (1 / 3))


def test_turbulence_only_shortens_the_equivalent_span():
    moon, _ = flakes.moon_column()
    earth = flakes.earth_column()
    law = 'mitchell1996_mixed_aggregates'
    calm = survival.equivalent_span(moon, earth, .06, 1e-4, 1e-4, law, 'weight_n')
    rough = survival.equivalent_span(moon, earth, .06, 1e-4, 1e-3, law)
    assert rough < calm
    assert math.isfinite(rough) and rough > .06
