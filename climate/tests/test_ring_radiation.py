import numpy as np

from climate.crm import ring_radiation as rr


def test_vapour_mole_fraction_from_the_mixing_ratio():
    assert rr.mole_fraction(0.0, 0.02897) == 0.0
    x = rr.mole_fraction(0.01, 0.02897)
    assert abs(x - 0.01 * 0.02897 / 0.01801528 / (1 + 0.01 * 0.02897 / 0.01801528)) < 1e-12


def test_levels_carry_the_layers_to_the_interfaces_and_continue_above():
    zw = np.array([0.0, 1000.0, 2000.0, 3000.0])
    zh = 0.5 * (zw[1:] + zw[:-1])
    prs = 1.0e5 * np.exp(-zh / 8000.0)
    temp = 290.0 - 0.0065 * zh
    p, t, q = rr.column_levels(zw, zh, prs, temp, np.full(3, 0.005), psfc=1.0e5, top_pa=1.0, extra=4)
    assert p.size == zw.size + 4 and np.isclose(p[0], 1.0e5) and np.isclose(p[-1], 1.0)
    assert np.allclose(p[1:3], 1.0e5 * np.exp(-zw[1:3] / 8000.0))
    assert np.allclose(t[:4], 290.0 - 0.0065 * zw) and np.all(t[4:] == t[3])
    assert np.all(np.diff(p) < 0) and np.allclose(q, 0.005)
