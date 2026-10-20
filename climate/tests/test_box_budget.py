"""The box budget's reading of CM1's domain diagnostics and of the prescribed day-night heating (no CM1 run)."""
import numpy as np
import pytest

from climate.crm import box_budget as bb
from climate.crm import cm1_output as co


def test_the_day_night_heating_follows_the_patch_round_the_day():
    z = np.array([0.0, 1000.0, 2000.0])
    th = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0], [4.0, 4.0, 4.0]])   # four bins of 90 degrees
    zh = np.array([500.0, 1500.0, 2500.0])
    # var18 = -180 puts the hour angle at -180 degrees at time 0: halfway between the last bin's centre and the first's
    assert np.allclose(bb.day_night_heating(z, th, zh, 0.0, -180.0, 360.0), [2.5, 2.5, 0.0])
    # 45 degrees later the hour angle stands at the first bin's centre
    assert np.allclose(bb.day_night_heating(z, th, zh, 45.0, -180.0, 360.0), [1.0, 1.0, 0.0])
    # a quarter of the way on from there, between the first and second bins
    assert np.allclose(bb.day_night_heating(z, th, zh, 67.5, -180.0, 360.0), [1.25, 1.25, 0.0])


def test_layer_sums_weigh_each_level_by_its_density_and_depth():
    zw = np.array([0.0, 1000.0, 3000.0, 14000.0, 30000.0])
    rho = np.array([1.4, 1.3, 1.0, 0.5])
    rate = np.array([1.0, 2.0, -1.0, 0.5])                               # K/day
    out = bb.layer_sums(rate, rho, zw)
    assert out['0-1.77 km'] == pytest.approx(1.0 * 1.4 * co.CP * 1000.0 / bb.DAY)
    assert out['1.77-6.53 km'] == pytest.approx(2.0 * 1.3 * co.CP * 2000.0 / bb.DAY)
    assert out['column'] == pytest.approx(sum(r * p * co.CP * d for r, p, d in zip(rate, rho, np.diff(zw))) / bb.DAY)


def test_a_diagnostic_output_reads_by_its_descriptor(tmp_path):
    (tmp_path / 'cm1out_diag_s.ctl').write_text('dset ^cm1out_diag_00%y4_s.dat\nxdef 1 linear 0 1\nydef 1 linear 0 1\n'
                                                'zdef 3 levels\n0.05\n0.15\n0.25\nvars 2\nmtime 0 99 time\n'
                                                'th 3 99 theta\nendvars\n')
    path = tmp_path / 'cm1out_diag_000002_s.dat'
    np.array([10800.0, 300.0, 301.0, 302.0], dtype='<f4').tofile(path)
    d = bb.read_diag(tmp_path, path)
    assert d['mtime'] == 10800.0 and np.allclose(d['th'], [300.0, 301.0, 302.0])
    np.array([10800.0, 300.0, 301.0], dtype='<f4').tofile(path)
    with pytest.raises(ValueError):
        bb.read_diag(tmp_path, path)
