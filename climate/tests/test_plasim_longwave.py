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


def test_the_error_summary_reads_the_layers_from_the_ground():
    errs = np.array([[-2.0, 1.0, 0.5, 3.0], [-4.0, 1.0, -0.5, 1.0]])      # three layers, then the outgoing longwave
    out = pl.error_summary(errs, np.array([True, True]))
    assert out['columns'] == 2 and list(out['mean']) == ['layer_1', 'layer_2', 'layer_3', 'olr']
    assert out['mean']['layer_1'] == -3.0 and out['rms']['layer_3'] == 0.5
    assert out['largest_layer_mean'] == 3.0                            # the layers only, not the outgoing longwave


def test_the_check_near_301_k_binds_the_corrections_it_checks():
    import hashlib
    import json
    import pytest
    if not pl.CHECK.exists():
        pytest.skip('the longwave check has not been run')
    out = json.loads(pl.CHECK.read_text())
    assert out['schema'] == pl.CHECK_SCHEMA and out['holds_w_m2'] == pl.HOLDS_W_M2
    assert set(out['corrections']) == set(pl.CORRECTIONS)
    for name, path in pl.CORRECTIONS.items():
        assert out['corrections'][name]['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    near = [c for c in out['cases'] if abs(c['global_ground_k'] - 301.0) < 0.75]
    assert [f"{c['run']} year {c['year']}" for c in near] == out['near_301k']['cases'] and near
    for name in pl.CORRECTIONS:
        worst = max(s[name][g]['largest_layer_mean'] for c in near for s in c['samples'].values()
                    if name in s for g in ('land', 'sea') if g in s[name])
        assert worst == out['near_301k'][name]['largest_layer_mean_w_m2']
        assert out['near_301k'][name]['holds'] == (worst <= pl.HOLDS_W_M2)
    if 'all_cases' in out:                                              # products since 2026-10-21
        assert out['all_cases']['cases'] == [f"{c['run']} year {c['year']}" for c in out['cases']]
    # Each correction is marked on the sample it was fitted on, and only there.
    marked = [(c['run'], c['year'], k) for c in out['cases'] for k, s in c['samples'].items()
              for name in pl.CORRECTIONS if name in s and s[name]['fitted_on_these_columns']]
    assert ('A28_corrected_970_lw_ice', 2, 'cloud_free') in marked and len(marked) <= len(pl.CORRECTIONS)
