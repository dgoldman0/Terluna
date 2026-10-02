"""Checks of the small-ion balance, mobility scaling and charge relaxation."""
import numpy as np
import pytest

from atmosphere.electricity import conductivity as cd


def test_the_ion_balance_has_its_two_limits_and_satisfies_the_equation():
    q, alpha = 2.0e6, 1.6e-12                                            # 2 ion pairs per cm3 per s; 1.6e-6 cm3/s
    assert cd.ion_density(q, alpha) == pytest.approx(np.sqrt(q / alpha))   # recombination only
    assert cd.ion_density(q, 0.0, 0.01) == pytest.approx(q / 0.01)         # attachment only
    n = cd.ion_density(q, alpha, 3.0e-3)
    assert alpha * n ** 2 + 3.0e-3 * n == pytest.approx(q)


def test_mobility_scales_inversely_with_density_and_conductivity_and_relaxation_follow():
    mu = cd.mobility(1.4e-4, np.array([2.5e25, 1.25e25]), 2.5e25)
    assert mu == pytest.approx([1.4e-4, 2.8e-4])
    sigma = cd.conductivity(5.0e8, 1.4e-4, 1.9e-4)
    assert sigma == pytest.approx(1.602176634e-19 * 5.0e8 * 3.3e-4)
    assert cd.relaxation_time(1.0e-14) == pytest.approx(885.4, rel=1e-3)   # about 15 minutes at 1e-14 S/m
