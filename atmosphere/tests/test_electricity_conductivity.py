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


def test_the_ion_chain_reproduces_earth_s_measured_fair_weather_conductivity():
    from atmosphere.electricity import fetch_inputs
    try:
        fetch_inputs.path('CRII_tables.zip')
    except FileNotFoundError:
        pytest.skip('CRAC:CRII tables not restored (python -m atmosphere.electricity.fetch_inputs --download)')
    z = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
    q, sigma, positive = cd.earth_column(z, pc_gv=3.0, phi_mv=400.0)          # Germany's cutoff, solar minimum
    assert positive == pytest.approx(cd.gringel_positive(z), rel=0.10)
    assert 10e6 < q[1] < 50e6                                                  # tens of ion pairs per cm3 per s at 10 km


def test_the_crac_crii_reader_returns_the_published_values():
    from atmosphere.electricity import fetch_inputs
    try:
        depths, y = cd.crac_crii(fetch_inputs.path('CRII_tables.zip'), 0.1, 400.0)
    except FileNotFoundError:
        pytest.skip('CRAC:CRII tables not restored')
    assert y[np.argmin(abs(depths - 1025.0))] == pytest.approx(1500.0)
    assert y[np.argmin(abs(depths - 705.0))] == pytest.approx(5182.0)


def test_droplets_and_aerosol_take_up_ions_as_their_coefficients_say():
    assert cd.aerosol_attachment(0.1, 1.0e9) == pytest.approx((4.36e-6 - 9.2e-8) * 1e-6 * 1e9)
    mu = 2.0e-4
    d = mu * cd.BOLTZMANN * 270.0 / cd.ELEMENTARY_CHARGE
    assert cd.droplet_attachment(270.0, mu, 1e8, 7e-6) == pytest.approx(4 * np.pi * d * 1e8 * 7e-6)
    assert cd.recombination(300.0, 2.5e25) == pytest.approx((6e-8 + 6e-26 * 2.5e19) * 1e-6)
