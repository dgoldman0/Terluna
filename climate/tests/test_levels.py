"""The 14- and 15-level sets and the conservative remap that carries a restart onto them (climate/gcm/levels.py)."""
import struct

import numpy as np
import pytest

from climate.gcm import levels as lv


def test_plasims_own_grid_is_its_formula_with_full_levels_midway():
    half = lv.half_levels(lv.DEFAULT)
    full = lv.full_levels(half)
    assert np.allclose(full, [0.0383, 0.1191, 0.21085, 0.31685, 0.4368, 0.5668, 0.69935, 0.82335, 0.9241, 0.9833])
    assert half[-1] == 1.0 and lv.namelist_entries(lv.DEFAULT) == []


def test_the_14_level_set_keeps_the_upper_six_levels_and_the_lowest_and_lands_the_new_ones():
    ten, l14 = lv.half_levels('10'), lv.half_levels('l14')
    f10, f14 = lv.full_levels(ten), lv.full_levels(l14)
    assert l14.size == 14 and l14[-1] == 1.0 and np.all(np.diff(l14) > 0)
    assert np.array_equal(l14[:6], ten[:6]) and np.array_equal(f14[:6], f10[:6])    # the same numbers
    assert l14[12] == ten[8] and f14[13] == f10[9]                                   # the lowest layer
    z = lv.heights_km(f14)
    for target in lv.L14_TARGETS_KM:
        assert np.min(np.abs(z - target)) < 0.05, target
    assert abs(z[6] - 19.24) < 0.01 and abs(z[8] - 10.24) < 0.01 and abs(z[11] - 4.04) < 0.01
    dsig = np.diff(np.r_[0.0, l14])
    assert np.all(dsig[1:] / dsig[:-1] > 0.5) and np.all(dsig[1:] / dsig[:-1] < 1.35)


def test_the_15_level_set_is_the_14_with_its_lowest_layer_split_near_0_6_km():
    l14, l15 = lv.half_levels('l14'), lv.half_levels('l15')
    f14, f15 = lv.full_levels(l14), lv.full_levels(l15)
    assert l15.size == 15 and l15[-1] == 1.0 and np.all(np.diff(l15) > 0)
    assert np.array_equal(l15[:13], l14[:13]) and np.array_equal(f15[:13], f14[:13])  # every layer above 1.73 km
    assert l15[13] == 0.9883 and f15[14] == 0.5 * (0.9883 + 1.0)                     # midway in sigma
    z = lv.heights_km(f15)
    assert abs(z[14] - lv.L15_LOWEST_KM) < 0.005 and abs(lv.heights_km(l15[13]) - 0.60) < 0.005
    assert abs(z[13] - 1.16) < 0.01
    dsig = np.diff(np.r_[0.0, l15])
    assert np.all(dsig[1:] / dsig[:-1] > 0.5) and np.all(dsig[1:] / dsig[:-1] < 1.35)
    assert sorted(lv.shared_layers(l14, l15)) == list(range(13))
    assert sorted(lv.shared_layers(lv.half_levels('10'), l15)) == [0, 1, 2, 3, 4, 5]


@pytest.mark.parametrize('name', ['l14', 'l15'])
def test_a_run_on_the_set_writes_its_half_levels_and_plasims_top_friction(name):
    entries = dict(lv.namelist_entries(name))
    assert entries['NEQSIG'] == '-1'
    half = [float(x) for x in entries['SIGH'].split(',')]
    assert half == list(lv.half_levels(name))                        # the values round-trip exactly
    friction = [float(x) for x in entries['TFRC'].split(',')]
    assert friction[:2] == [20.0 * 86400.0, 100.0 * 86400.0] and friction[2:] == [0.0] * (len(half) - 2)
    with pytest.raises(ValueError):
        lv.half_levels('l99')


def test_the_remap_keeps_column_integrals_reproduces_linear_profiles_and_copies_shared_layers():
    a, b = lv.half_levels('10'), lv.half_levels('l14')
    da, db = np.diff(np.r_[0.0, a]), np.diff(np.r_[0.0, b])
    rng = np.random.default_rng(3)
    x = rng.normal(size=(10, 7)) * 50.0 + 250.0
    y = lv.remap(x, a, b)
    assert y.shape == (14, 7)
    assert np.allclose(db @ y, da @ x, rtol=0, atol=1e-11)
    for k, j in lv.shared_layers(a, b).items():
        assert np.array_equal(y[k], x[j])
    assert sorted(lv.shared_layers(a, b)) == [0, 1, 2, 3, 4, 5, 13]
    fa, fb = lv.full_levels(a), lv.full_levels(b)
    assert np.allclose(lv.remap(2.0 - 3.0 * fa, a, b), 2.0 - 3.0 * fb, atol=1e-13)
    # back onto the same grid is the identity
    assert np.array_equal(lv.remap(x, a, a), x)


@pytest.mark.parametrize('source', ['10', 'l14'])
def test_the_remap_onto_15_levels_keeps_column_integrals_and_splits_only_what_it_must(source):
    a, b = lv.half_levels(source), lv.half_levels('l15')
    da, db = np.diff(np.r_[0.0, a]), np.diff(np.r_[0.0, b])
    rng = np.random.default_rng(11)
    x = rng.normal(size=(a.size, 9)) * 40.0 + 260.0
    y = lv.remap(x, a, b)
    assert y.shape == (15, 9)
    assert np.allclose(db @ y, da @ x, rtol=0, atol=1e-11)
    shared = lv.shared_layers(a, b)
    for k, j in shared.items():
        assert np.array_equal(y[k], x[j])
    # both grids' lowest layer starts at sigma 0.9666, and its two parts on 15 levels hold its mass-weighted mean
    assert a[-2] == b[-3]
    assert np.allclose((db[13] * y[13] + db[14] * y[14]) / (db[13] + db[14]), x[-1], atol=1e-11)
    fa, fb = lv.full_levels(a), lv.full_levels(b)
    assert np.allclose(lv.remap(2.0 - 3.0 * fa, a, b), 2.0 - 3.0 * fb, atol=1e-13)


def _restart(path, nlev, records):
    rows = [('nstep', struct.pack('<i', 12345)), ('nlev', struct.pack('<i', nlev))]
    rows += [(name, np.ascontiguousarray(v, '<f8').tobytes()) for name, v in records]
    lv.write_records(path, rows)


def test_a_restart_is_carried_onto_the_set_with_its_other_records_untouched(tmp_path):
    a, b = lv.half_levels('10'), lv.half_levels('l14')
    rng = np.random.default_rng(5)
    st = rng.normal(size=(10, 506))
    aadt = rng.normal(size=(11, 2048))
    dcc = rng.uniform(0.0, 1.0, size=(10, 2048))
    soil = rng.normal(size=(5, 2048))
    src, dst = tmp_path / 'rest10', tmp_path / 'rest14'
    _restart(src, 10, [('st', st), ('aadt', aadt), ('dcc', dcc), ('dsoilt', soil)])
    report = lv.carry_restart(src, dst, a, b)
    got = {name.strip(): data for name, data in lv.read_records(dst)}
    assert struct.unpack('<i', got['nlev'])[0] == 14 and struct.unpack('<i', got['nstep'])[0] == 12345
    st14 = np.frombuffer(got['st'], '<f8').reshape(14, 506)
    assert np.allclose(np.diff(np.r_[0.0, b]) @ st14, np.diff(np.r_[0.0, a]) @ st, atol=1e-12)
    assert np.array_equal(st14[:6], st[:6]) and np.array_equal(st14[13], st[9])
    aadt14 = np.frombuffer(got['aadt'], '<f8').reshape(15, 2048)
    assert np.array_equal(aadt14[14], aadt[10])                        # the ground's value kept
    dcc14 = np.frombuffer(got['dcc'], '<f8').reshape(14, 2048)
    assert dcc14.min() >= 0.0 and dcc14.max() <= 1.0
    assert got['dsoilt'] == np.ascontiguousarray(soil, '<f8').tobytes()
    assert set(report['fields']) == {'st', 'aadt', 'dcc'}
    # an unknown field held by level is refused, and so is a restart on another grid
    _restart(tmp_path / 'odd', 10, [('zzz', rng.normal(size=(10, 2048)))])
    with pytest.raises(ValueError):
        lv.carry_restart(tmp_path / 'odd', tmp_path / 'odd14', a, b)
    with pytest.raises(ValueError):
        lv.carry_restart(dst, tmp_path / 'again', a, b)


def test_the_15_level_inputs_split_the_lowest_layers_terms_and_keep_their_sums():
    import json
    from climate.gcm import exoplasim_run as er
    # the twilight terms: the thirteen upper layers as on 14 levels, the two lowest summing to its lowest layer's
    t14, t15 = (json.loads(er.levels_product(er.TWILIGHT_LEDGERS['301k'], name).read_text()) for name in ('l14', 'l15'))
    assert np.allclose(t15['levels']['half_levels'], lv.half_levels('l15'), rtol=1e-6, atol=0)   # six figures kept
    for key in ('gcm_term', 'gcm_term_cone'):
        a, b = np.array(t15[key]['layers']), np.array(t14[key]['layers'])
        assert a.shape == (15, b.shape[1])
        np.testing.assert_allclose(a[:13], b[:13], rtol=0, atol=1e-12)
        np.testing.assert_allclose(a[13] + a[14], b[13], rtol=1e-4, atol=1e-8)
    for key in ('curved', 'cone'):
        a, b = t15['longwave']['gcm_layers'][key], t14['longwave']['gcm_layers'][key]
        assert len(a) == 15 and abs(sum(a) - sum(b)) < 1e-4 and abs(a[13] + a[14] - b[13]) < 1e-4
    # the 301 K refit carried from ten levels: a value for each layer and the ground, the split layers' sum kept
    carried = json.loads(er.levels_product(er.LONGWAVE_FITS['301k'], 'l15').read_text())
    for key in ('land', 'sea'):
        assert len(carried['correction'][key]) == 16
        assert abs(carried['check'][key]['sums'][0] - carried['check'][key]['sums'][1]) < 1e-3


def test_the_runner_reads_the_provisional_refit_as_it_is_on_15_levels_and_refuses_it_elsewhere():
    import json
    from climate.gcm import exoplasim_run as er
    path = er.longwave_product('l15_provisional', 'l15', 15)
    assert path == er.LONGWAVE_FITS['l15_provisional']
    refit = json.loads(path.read_text())
    assert refit['levels']['name'] == 'l15' and refit['levels']['half_levels'] == list(lv.half_levels('l15'))
    assert all(len(refit['correction'][k]) == 16 for k in ('land', 'sea', 'tref'))
    assert refit['correction']['land'][-1] == 0.0 and refit['correction']['sea'][-1] == 0.0   # the ground's
    assert er.longwave_product('301k', 'l15', 15).name == 'longwave_A28_corrected_970_lw_ice_l15.json'
    assert er.longwave_product('301k', lv.DEFAULT, 10) == er.LONGWAVE_FITS['301k']
    for name, layers in (('l14', 14), (lv.DEFAULT, 10)):
        with pytest.raises(ValueError):
            er.longwave_product('l15_provisional', name, layers)
    with pytest.raises(ValueError):
        er.longwave_product('301k', 'l15', 14)                     # a layer count other than the correction's
