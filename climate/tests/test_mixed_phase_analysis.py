"""Checks of the mixed-phase analysis: Morrison's size distributions and fall speeds at the host gravity, the
charging zone, patches on a periodic domain and the levels where the air crosses a temperature."""
from math import gamma

import numpy as np
import pytest

from climate.crm import mixed_phase_analysis as mpa

GRAUPEL, ICE = mpa.SPECIES['graupel'], mpa.SPECIES['ice']


def test_fall_speeds_scale_with_gravity_as_the_patch_does():
    g = 1.624
    assert mpa.gravity_factor(GRAUPEL, g) == pytest.approx((g / 9.81) ** (1.37 / 3.0))
    assert mpa.gravity_factor(GRAUPEL, g) == pytest.approx(0.44, abs=0.005)
    assert mpa.gravity_factor(ICE, g) == pytest.approx(0.30, abs=0.005)
    assert mpa.gravity_factor(GRAUPEL, 9.81) == 1.0


def test_the_slope_follows_mass_and_number_and_stays_within_the_scheme_s_bounds():
    q, n = 1.0e-3, 1.0e3                                                  # 1 g/kg in 1000 particles per kg
    lam, n_adj = mpa.slope(q, n, GRAUPEL)
    assert lam == pytest.approx((np.pi * 400.0 * n / q) ** (1.0 / 3.0))
    assert n_adj == pytest.approx(n)
    lam, n_adj = mpa.slope(1.0e-3, 1.0e-3, GRAUPEL)                       # too few particles: held at 2 mm
    assert lam == pytest.approx(1.0 / 2000.0e-6) and n_adj == pytest.approx(lam ** 3 * 1.0e-3 / (np.pi * 400.0))
    lam, n_adj = mpa.slope(0.0, 10.0, GRAUPEL)
    assert np.isnan(lam) and n_adj == 0.0


def test_the_mass_weighted_speed_is_the_distribution_s_and_the_cap_scales_with_gravity():
    lam, rho = 1.0 / 1.0e-3, mpa.RHOSU
    # integrate V(D) D^3 exp(-lam D) over D against the closed form
    d = np.linspace(1.0e-7, 0.05, 400001)
    weight = d ** 3 * np.exp(-lam * d)
    numeric = np.trapezoid(mpa.speed_at(d, rho, GRAUPEL, 1.624) * weight, d) / np.trapezoid(weight, d)
    assert mpa.fall_speed(lam, rho, GRAUPEL, 1.624) == pytest.approx(numeric, rel=1e-4)
    assert mpa.fall_speed(lam, rho, GRAUPEL, 1.624, 'number') == pytest.approx(
        mpa.gravity_factor(GRAUPEL, 1.624) * 19.3 * gamma(1.37) / lam ** 0.37)
    assert mpa.fall_speed(1.0, rho, GRAUPEL, 1.624) == pytest.approx(20.0 * mpa.gravity_factor(GRAUPEL, 1.624))


def test_the_charging_zone_needs_all_three_between_0_and_minus_40_c():
    t = np.array([-5.0, -5.0, -5.0, 2.0, -45.0, -15.0])
    qc = np.array([2e-5, 0.0, 2e-5, 2e-5, 2e-5, 2e-5])
    qg = np.array([2e-5, 2e-5, 2e-5, 2e-5, 2e-5, 2e-5])
    qi = np.array([1e-5, 1e-5, 0.0, 1e-5, 1e-5, 0.4e-5])
    qs = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.7e-5])                       # ice and snow count together
    assert list(mpa.zone_mask(t, qc, qg, qi, qs, 1e-5)) == [True, False, False, False, False, True]


def test_patches_join_across_the_periodic_edges():
    mask = np.zeros((5, 6), bool)
    mask[2, 0] = mask[2, 5] = True                                        # one patch split by the x edge
    mask[0, 3] = mask[4, 3] = True                                        # one split by the y edge
    labels = mpa.periodic_labels(mask)
    assert labels[2, 0] == labels[2, 5] != 0 and labels[0, 3] == labels[4, 3] != 0
    assert len(np.unique(labels[labels > 0])) == 2
    ring = mpa.periodic_labels(np.array([True, False, True, True, False, True]))
    assert len(np.unique(ring[ring > 0])) == 2                            # the ring's ends are one patch


def test_the_crossing_height_interpolates_and_ignores_warming_aloft():
    zh = np.array([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
    t = np.array([10.0, 4.0, -2.0, -8.0, 5.0])                            # warming again above, like a stratosphere
    assert mpa.crossing_height(zh, t, 0.0) == pytest.approx(1000.0 + 1000.0 * 4.0 / 6.0)
    assert np.isnan(mpa.crossing_height(zh, t, -40.0))


def test_the_charging_pass_matches_a_direct_pair_integral_and_earth_speeds_charge_faster():
    from atmosphere.electricity import charging as ch
    law = lambda t_c, lwc, vg, other: ch.power_law(1.0e-15, 100.0e-6, 1.0, 2.0, 2.5)
    rho = np.array([0.8])
    cells = dict(rho=rho, t_c=np.array([-15.0]), lwc_g_m3=np.array([0.2]),
                 n={'graupel': np.array([800.0]), 'ice': np.array([3.0e3]), 'snow': np.array([2.0e3])},
                 lam={'graupel': np.array([1.0 / 1.5e-3]), 'ice': np.array([1.0 / 0.1e-3]), 'snow': np.array([1.0 / 1.0e-3])})
    out = mpa.cell_charging(cells, 1.624, {'test': law})
    speed = lambda key, g: (lambda d: mpa.speed_at(d, rho[:, None, None], mpa.SPECIES[key], g))
    direct = sum(ch.pair_rate(800.0, 1.0 / 1.5e-3, cells['n'][o], cells['lam'][o], speed('graupel', 1.624),
                              speed(o, 1.624), law(None, None, None, o), mpa.separation_efficiency(o, -15.0))
                 for o in ('ice', 'snow'))
    assert out[('test', 'run')] == pytest.approx(direct)
    assert out[('test', 'earth_speeds')][0] > 10.0 * out[('test', 'run')][0]


def test_separation_follows_the_nssl_settings():
    assert mpa.separation_efficiency('ice', -10.0) == pytest.approx(0.5 * (1.0 - 0.1 * np.exp(-1.0)))
    snow = mpa.separation_efficiency('snow', -10.0)
    stick = 0.1 * np.exp(-1.0) * (400.0 - 300.0) / 300.0
    assert snow(np.array([30e-6, 95e-6, 400e-6])) == pytest.approx([0.0, 0.125 * (1 - stick), 0.25 * (1 - stick)])
