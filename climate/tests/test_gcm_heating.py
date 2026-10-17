import numpy as np

from climate.crm import gcm_heating as gh


def test_gcm_levels_run_from_the_ground_through_the_half_levels_and_on_to_the_top():
    sigma = np.array([0.0383, 0.1191, 0.2108, 0.3168, 0.4368, 0.5668, 0.6994, 0.8233, 0.9241, 0.9833])
    ta = np.linspace(220.0, 293.0, 10)
    p, t, q = gh.gcm_column_levels(sigma, ta, np.linspace(1e-6, 0.01, 10), 121590.0, top_pa=1.0, extra=8)
    assert p.size == 10 + 9 and np.isclose(p[0], 121590.0) and np.isclose(p[-1], 1.0)
    assert np.all(np.diff(p) < 0)
    assert np.isclose(p[9], 0.0766 * 121590.0, rtol=5e-3) and np.isclose(p[10], 0.0383 * 121590.0)
    assert t[0] > ta[-1] and np.all(t[10:] == ta[0])                    # warmer at the ground, isothermal above
