"""The hydrogen escape's diffusion limit: Hunten's flux, its independence of the homopause's height, the carriers and
methane's share broken by the glow."""
import math

import pytest

from shared.constants import MOON_RADIUS
from atmosphere.loss_response import hydrogen as hy


def test_the_limiting_flux_gives_earths_textbook_value():
    # On Earth b f (m_a - m) g / (k T) for H2 near 200 K is about 2.5e13 f cm^-2 s^-1 (Hunten and Strobel 1974).
    from shared.constants import BOLTZMANN
    b = hy.diffusion_parameter('H2', 200.0) * 1e-2                      # cm^-1 s^-1
    per_cm = (28.96 - 2.016) * hy.AMU_KG * 1e3 * 980.0 / (BOLTZMANN * 1e7 * 200.0)
    assert 2.0e13 < b * per_cm < 2.6e13


def test_the_moons_loss_does_not_depend_on_where_the_homopause_lies():
    for r in (1.0, 1.3, 2.0):
        radius = r * MOON_RADIUS
        whole = hy.limiting_flux_per_m2('H2', 1e-6, 170.0, radius) * 4 * math.pi * radius ** 2
        assert whole == pytest.approx(hy.limiting_rate('H2', 1e-6, 170.0), rel=1e-12)


def test_methane_unbroken_diffuses_slower_than_its_hydrogen_as_h2():
    m = hy.air_mass_amu()
    unbroken = hy.hydrogen_kg_s({'CH4': hy.limiting_rate('CH4', 1e-6, 170.0, m)})
    broken = hy.hydrogen_kg_s({'H2': hy.limiting_rate('H2', 2e-6, 170.0, m)})
    assert 0.1 < unbroken / broken < 0.2
    kg = hy.escape_rates(170.0, 3e-8, 0.25)
    assert kg['total_unbroken'] < kg['total'] < kg['total_broken']
    assert kg['total'] == pytest.approx(kg['without_methane'] + 0.75 * kg['methane_unbroken'] + 0.25 * kg['methane_broken'])


def test_the_glow_breaks_methane_on_a_scale_of_years_above_the_base():
    share = hy.broken_share(170.0, 1.2 * MOON_RADIUS)
    assert 300.0 < share['glow_lifetime_days'] < 3000.0
    assert 0.0 < share['share'] < 1.0 and share['homopause_pa'] < hy.BASE_PA


def test_the_product_checks_the_oxygen_steps_hydrogen_and_counts_methane():
    import json
    if not hy.OUT.exists():
        pytest.skip('the hydrogen escape has not been run')
    out = json.loads(hy.OUT.read_text())
    assert all(0.8 < r < 1.2 for r in out['check']['ratio'])            # the full transport sits at the limit
    for tr in out['treatments'].values():
        for row in tr['grounds'].values():
            assert row['methane'] > row['without_methane']               # methane carries most of the hydrogen
            assert row['cycles']['waters_years'] > 1e9
