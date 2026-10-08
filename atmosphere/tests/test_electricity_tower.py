"""The tower screen's field solver solves the electrified CM1's difference equations."""
import numpy as np

from atmosphere.electricity import tower
from climate.tests.test_cm1_elec_field import direct_potential, stretched


def test_potential_matches_the_direct_solve_of_the_same_equations():
    rng = np.random.default_rng(3)
    zh, zf = stretched(7, 100.0, 1.3, 20000.0)
    q = rng.normal(size=(5, 4, 7)) * 1e-9                       # ni, nj, nk as the direct solve takes them
    ref = direct_potential(q, zh, zf, 3000.0, 2500.0)
    phi = tower.potential(np.transpose(q, (2, 1, 0)), zh, zf, 3000.0, 2500.0)   # nk, nj, ni
    assert np.allclose(np.transpose(phi, (2, 1, 0)), ref, rtol=1e-9, atol=1e-6 * np.abs(ref).max())


def test_a_charged_layer_between_grounded_plates():
    zh, zf = stretched(60, 50.0, 1.05, 30000.0)
    q = np.zeros((zh.size, 3, 3)); k = np.argmin(np.abs(zh - 10000.0)); q[k] = 1e-9
    phi = tower.potential(q, zh, zf, 1000.0, 1000.0)[:, 1, 1]
    sigma, z0, top = 1e-9 * (zf[k + 1] - zf[k]), zh[k], zf[-1]
    exact = sigma / tower.EPS0 * z0 * (top - z0) / top          # sheet charge between two grounded plates
    assert abs(phi[k] - exact) / exact < 0.02


def test_rizk_potential_falls_with_density():
    assert abs(tower.rizk_mv(1e6, 1.0) - 1.556) < 1e-3
    assert tower.rizk_mv(24000.0, 0.5) < tower.rizk_mv(24000.0, 1.0)
