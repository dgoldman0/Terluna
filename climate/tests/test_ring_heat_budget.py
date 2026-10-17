import numpy as np
import pytest

from climate.crm import ring_heat_budget as hb


def test_the_nudging_weight_rises_from_its_floor_to_full_strength():
    weight = hb.nudge_weight([0.0, 8000.0, 12000.0, 16000.0, 20000.0], 8000.0, 16000.0)
    assert weight.tolist() == [0.0, 0.0, 0.5, 1.0, 1.0]


def test_layer_sums_and_means_weigh_each_level_by_its_depth_and_density():
    zw = np.array([0.0, 100.0, 300.0, 600.0, 1000.0])
    rho = np.array([1.2, 1.1, 1.0, 0.9])
    heating = np.full(4, 2.0)
    assert hb.layer_integral(heating, rho, zw, 0.0, 1000.0) == pytest.approx(2.0 * hb.CP * (120 + 220 + 300 + 360))
    mean = hb.layer_mean(np.array([10.0, 20.0, 30.0, 40.0]), rho, zw, 0.0, 300.0)
    assert mean == pytest.approx((10 * 120 + 20 * 220) / 340)


def test_the_forcing_files_read_as_cm1_writes_them(tmp_path):
    (tmp_path / 'wls.txt').write_text('3\n0.000 0.0\n1000.000 -1.0e-03\n150000.000 0.0\n')
    z, w = hb.read_wls(tmp_path / 'wls.txt')
    assert z.tolist() == [0.0, 1000.0, 150000.0] and w[1] == -1.0e-3
    (tmp_path / 'ref.dat').write_text(' * h\n 0.0 1.0e30\n * p\n * z theta\n'
                                      ' 0.0 278.5 9.2 0.0 0.0\n 50.0 278.6 9.2 0.0 0.0\n')
    z, th = hb.read_reference(tmp_path / 'ref.dat')
    assert z.tolist() == [0.0, 50.0] and th.tolist() == [278.5, 278.6]
    assert hb.namelist_number(' lsnudge_tau        =    259200.0,\n var14     =   8000.0,\n', 'var14') == 8000.0


def test_block_means_average_consecutive_outputs():
    a = np.arange(12.0).reshape(12, 1)
    assert np.allclose(hb.block_means(a, 3)[:, 0], [1.0, 4.0, 7.0, 10.0])
    assert np.allclose(hb.block_means(a, 5)[:, 0], [2.0, 7.0])           # the last two outputs make no block
    assert hb.block_means(a, 1) is a


def test_coarser_means_miss_the_eddies_faster_than_them():
    # A wave that turns over every output: each output holds its covariance, and pairs of outputs average it away.
    times, lon = 22, 16
    phase = np.arange(times)[:, None] * np.pi + np.linspace(0.0, 2.0 * np.pi, lon, endpoint=False)[None, :]
    v = np.cos(phase)[:, None, None, :] * np.ones((1, 2, 3, 1))
    t = 2.0 * v
    fine, n = hb.zonal_eddy_sums(v, t, np.zeros_like(v))
    coarse, m = hb.zonal_eddy_sums(v, t, np.zeros_like(v), block=2)
    assert (n, m) == (22, 11)
    assert np.allclose(fine['eddy'] / n, 1.0)
    assert np.allclose(coarse['eddy'] / m, 0.0)
    assert np.allclose(fine['t'] / n, coarse['t'] / m)                  # the zonal means are the same either way
