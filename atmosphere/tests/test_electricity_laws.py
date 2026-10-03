"""Checks of the published charging laws as WRF-ELEC implements them."""
import numpy as np
import pytest

from atmosphere.electricity import charging as ch


def test_the_critical_rime_rate_follows_brooks_above_minus_15_and_saunders_peck_below():
    t = np.array([-5.0, -10.0, -15.0, -20.0, -23.7, -30.0, -34.0])
    rarc = ch.rar_critical(t)
    assert rarc[0] == pytest.approx(0.1)                                # Brooks gives none above -7.35 C; floored
    assert rarc[1] == pytest.approx(-1.47 + 2.0)
    assert rarc[2:5] == pytest.approx([1.65, 2.76, 3.40], abs=0.01)    # s(T) of Saunders and Peck (1998)
    assert rarc[5] == pytest.approx(1.80, abs=0.01) and rarc[6] == pytest.approx(0.1)


def test_the_charge_function_is_positive_above_the_critical_rate_and_a_parabola_below():
    t = -20.0
    rarc = float(ch.rar_critical(t))
    assert ch.saunders_peck_charge_fc(rarc + 1.0, t) == pytest.approx(6.74)
    mid = (rarc + 0.1) / 2.0
    assert ch.saunders_peck_charge_fc(mid, t) == pytest.approx(-6.5 * min(1.0, 0.6 * (rarc - 0.1)))
    assert ch.saunders_peck_charge_fc(np.array([0.05, 0.1]), t) == pytest.approx([0.0, 0.0])
    assert ch.saunders_peck_charge_fc(5.0, -33.0) == 0.0                # no charging below -32.47 C


def test_the_size_factor_joins_across_its_size_breaks():
    for d, positive in ((155e-6, True), (452e-6, True), (253e-6, False)):
        below, above = ch.size_speed_factor(d * (1 - 1e-6), 1.0, positive), ch.size_speed_factor(d * (1 + 1e-6), 1.0, positive)
        assert below == pytest.approx(above, rel=0.02)
    assert ch.size_speed_factor(100e-6, 2.0, True) == pytest.approx(ch.size_speed_factor(100e-6, 1.0, True) * 2 ** 2.5)


def test_the_saunders_peck_law_stops_charging_at_low_rime_rates():
    law = ch.saunders_peck(np.array([-10.0]), np.array([0.1]), lambda d: 1.0 + 0.0 * d)
    dq = law(np.array([2e-3]), np.array([300e-6]), np.array([1.0]))   # RAR = 0.7 x 0.1 x 1 = 0.07 < 0.1
    assert dq == pytest.approx([0.0])
    law = ch.saunders_peck(np.array([-10.0]), np.array([1.0]), lambda d: 3.0 + 0.0 * d)
    dq = law(np.array([2e-3]), np.array([300e-6]), np.array([3.0]))   # RAR 2.1 > 0.53: positive
    assert dq[0] == pytest.approx(4.0e6 * (300e-6) ** 1.9 * 3.0 ** 2.5 * 6.74 * (2.1 - 0.53) * 1e-15)


def test_takahashi_s_table_reads_and_interpolates_as_the_taka_routine():
    from atmosphere.electricity import fetch_inputs
    try:
        table = ch.takahashi_table(fetch_inputs.path('takahashi.txt'))
    except FileNotFoundError:
        pytest.skip('Takahashi table not restored (python -m atmosphere.electricity.fetch_inputs --download)')
    cwc, charge = table
    assert len(cwc) == 30 and charge.shape == (30, 31) and cwc[0] == pytest.approx(0.01) and cwc[-1] == 30.0
    # entries quoted in the literature note: 0.9 g/m3 at -5 C, 2.0 g/m3 at -14 C, 0.1 g/m3 at -10 C
    assert ch.takahashi_charge_fc(np.array([0.9, 2.0, 0.1]), np.array([-5.0, -14.0, -10.0]), table) == \
        pytest.approx([66.0, -33.0, 19.14])
    mid = ch.takahashi_charge_fc(0.95, -5.5, table)
    corners = [ch.takahashi_charge_fc(c, t, table) for c in (0.9, 1.0) for t in (-5.0, -6.0)]
    assert mid == pytest.approx(np.mean(corners))
    assert ch.takahashi_charge_fc(0.005, -10.0, table) == 0.0 and ch.takahashi_charge_fc(0.5, -40.0, table) == 0.0
    dq = ch.takahashi(table)(-5.0, 0.9)(None, 100e-6, 8.0)
    assert dq == pytest.approx(5.0 * 66.0e-15)                          # alpha = 5 at 100 um and 8 m/s


def test_the_low_speed_brackets_hold_their_measured_ranges():
    q = ch.low_speed_kumar(np.array([-10.0, -5.0, -20.0]), np.array([0.5, 0.5, 0.5]))(None, None, np.ones(3))
    assert q == pytest.approx([(1.392 * 0.25 - 3.149 * 0.5 + 8.108) * 1e-15, 0.0, 0.0])
    q = ch.low_speed_avila(np.array([-10.0, -15.0]), None)(None, None, np.ones(2))
    assert q == pytest.approx([0.2e-15, 0.0])
