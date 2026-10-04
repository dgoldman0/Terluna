"""Wave-slope identities, dispersion limits and file-format rejection."""
import numpy as np
import pytest

from .model import fraction_above, iter_swan, slope_moments, wavenumber


def test_dispersion_deep_and_shallow_limits():
    f = np.array([.02, .1, .4])
    g = 2.0
    assert wavenumber(f, 1e7, g) == pytest.approx((2 * np.pi * f) ** 2 / g, rel=1e-12)
    assert wavenumber(f, 1e-8, g) == pytest.approx(2 * np.pi * f / np.sqrt(g * 1e-8), rel=1e-7)


def test_isotropic_and_single_direction_slope_covariance():
    f, angle = np.linspace(.01, .3, 101), np.arange(0, 360, 10.)
    e = np.ones((len(f), len(angle)))
    r = slope_moments(f, angle, e, 100, 2)
    assert r['variance_m2'] == pytest.approx(.29 * 360)
    assert r['covariance'][0, 0] == pytest.approx(r['covariance'][1, 1])
    assert r['covariance'][0, 1] == pytest.approx(0, abs=1e-12)
    assert np.trace(r['covariance']) == pytest.approx(np.trapezoid(r['slope_density'], f))
    assert r['dispersion_relative_residual'] < 1e-14
    e[:, 1:] = 0
    r = slope_moments(f, angle, e, 100, 2)
    assert r['covariance'][1, 1] == 0
    assert r['covariance'][0, 0] == pytest.approx(np.trapezoid(r['slope_density'], f))


def test_sinusoid_slope_variance_and_rotation():
    f, angle = np.array([.05, .1, .15]), np.arange(0, 360, 45.)
    # One spectral line with variance a^2/2; central trapezoidal weight = .05 Hz.
    amplitude = .7
    e = np.zeros((3, 8)); e[1, 1] = amplitude ** 2 / 2 / (.05 * 45)
    r = slope_moments(f, angle, e, 25, 1.6)
    k = wavenumber(np.array([.1]), 25, 1.6)[0]
    assert r['variance_m2'] == pytest.approx(amplitude ** 2 / 2)
    assert r['covariance'] == pytest.approx(np.full((2, 2), (amplitude * k) ** 2 / 4))


def test_exact_cutoff_fraction():
    assert fraction_above([1, 2, 3], [1, 2, 3], 1.5) == pytest.approx((9 - 2.25) / 8)
    assert fraction_above([1, 2], [0, 0], 1.5) is None


def swan_text(body):
    return ('SWAN 1\n$ comment\nTIME\n1\nLONLAT\n2\n10 20\n30 40\nAFREQ\n2\n.1\n.2\n'
            'CDIR\n4\n0\n90\n180\n270\nQUANT\n1\nVaDens\nm2/Hz/degr\n-99\n' + body)


def test_swan_units_factor_and_missing_record(tmp_path):
    p = tmp_path / 'waves.spc'
    p.write_text(swan_text('20000101.000000\nFACTOR\n.01\n1 2 3 4\n5 6 7 8\nNODATA\n'
                           '20000101.010000\nZERO\nZERO\n'))
    rows = list(iter_swan(p))
    assert rows[0]['variance'][0] == pytest.approx(np.arange(1, 9).reshape(2, 4) / 100)
    assert rows[0]['available'].tolist() == [True, False]
    assert rows[1]['available'].all()
    assert rows[1]['variance'].sum() == 0
    p.write_text(swan_text('20000101.000000\nFACTOR\n.01\n1 2 3 4\n'))
    with pytest.raises(ValueError, match='Truncated'):
        list(iter_swan(p))


def test_invalid_spectrum_rejected():
    with pytest.raises(ValueError):
        slope_moments([.1, .2], [0, 30, 60], np.ones((2, 3)), 10, 2)
    with pytest.raises(ValueError):
        wavenumber([0], 10, 2)
