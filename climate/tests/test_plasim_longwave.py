import numpy as np

from climate.gcm import plasim_longwave as pl

SIGMA = np.array([0.0383, 0.1191, 0.2108, 0.3168, 0.4368, 0.5668, 0.6994, 0.8233, 0.9241, 0.9833])
T = np.append(np.linspace(215.0, 292.0, 10), 296.0)
Q = np.geomspace(2e-6, 0.012, 10)


def column(**kw):
    return pl.lwr(SIGMA, 121590.0, T, Q, 1.0, 1.6242, 1008.3, **kw)


def test_half_levels_come_back_from_midpoints():
    half = np.array([0.0, 0.1, 0.3, 0.6, 1.0])
    assert np.allclose(pl.half_sigma((half[:-1] + half[1:]) / 2.0), half)


def test_a_layers_correction_takes_its_escape_from_that_layer_and_gives_it_to_space():
    tref = np.append(T[:-1], T[-1])
    zero = dict(land=np.zeros(11), sea=np.zeros(11), tref=tref)
    one = dict(land=np.where(np.arange(11) == 8, 5.0, 0.0), sea=np.zeros(11), tref=tref)
    net0, h0, f0 = column()
    netz, hz, fz = column(correction=zero)
    net1, h1, f1 = column(correction=one)
    assert np.allclose(hz, h0) and np.isclose(fz['olr'], f0['olr'])
    cap = 1008.3 * 121590.0 * np.diff(pl.half_sigma(SIGMA)) / 1.6242 * (1.0 + pl.ADV * Q)   # J/m2/K, moist
    lost = (h0 - h1) * cap                                                      # W/m2 each layer gave up
    assert np.isclose(lost[8], 5.0) and np.allclose(np.delete(lost, 8), 0.0, atol=1e-9)
    assert np.isclose(f1['olr'] - f0['olr'], 5.0) and np.isclose(f1['ground_down'], f0['ground_down'])


def test_the_sea_takes_the_sea_correction():
    tref = T.copy()
    corr = dict(land=np.zeros(11), sea=np.full(11, 1.0), tref=tref)
    _, _, land = pl.lwr(SIGMA, 121590.0, T, Q, 1.0, 1.6242, 1008.3, correction=corr)
    _, _, sea = pl.lwr(SIGMA, 121590.0, T, Q, 0.0, 1.6242, 1008.3, correction=corr)
    _, _, bare = pl.lwr(SIGMA, 121590.0, T, Q, 0.0, 1.6242, 1008.3)
    _, _, bare_land = pl.lwr(SIGMA, 121590.0, T, Q, 1.0, 1.6242, 1008.3)
    assert np.isclose(land['olr'], bare_land['olr']) and np.isclose(sea['olr'] - bare['olr'], 11.0)
