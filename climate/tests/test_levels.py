"""The 14-level set and the conservative remap that carries a restart onto it (climate/gcm/levels.py)."""
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


def test_a_run_on_the_set_writes_its_half_levels_and_plasims_top_friction():
    entries = dict(lv.namelist_entries('l14'))
    assert entries['NEQSIG'] == '-1'
    half = [float(x) for x in entries['SIGH'].split(',')]
    assert half == list(lv.half_levels('l14'))                       # the values round-trip exactly
    friction = [float(x) for x in entries['TFRC'].split(',')]
    assert friction[:2] == [20.0 * 86400.0, 100.0 * 86400.0] and friction[2:] == [0.0] * 12
    with pytest.raises(ValueError):
        lv.half_levels('l15')


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
