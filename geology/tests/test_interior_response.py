import math

import numpy as np
import pytest
from scipy.special import lpmv

from geology import interior_response as ir


def test_harmonics_recover_a_known_degree():
    """The sag's harmonic analysis on a regular grid returns a field's degree-3 part and drops its degree-6 part."""
    step = 1.0
    lat = np.arange(-90 + step / 2, 90, step)
    lon = np.arange(-180 + step / 2, 180, step)
    x = np.cos(np.radians(90.0 - lat))[:, None]
    phi = np.radians(lon)[None, :]
    three = lpmv(2, 3, x) * np.cos(2 * phi)
    six = lpmv(4, 6, x) * np.sin(4 * phi)
    _, synth = ir.real_harmonics(three + 1e-3 * six, lat, lon, 8)
    part = synth({3: 1.0})
    assert np.max(np.abs(part - three)) / np.max(np.abs(three)) < 2e-3
    assert np.max(np.abs(synth({6: 1.0}) - 1e-3 * six)) / np.max(np.abs(1e-3 * six)) < 2e-2


def test_piecewise_load_history():
    a, b = np.eye(3), 3 * np.eye(3)
    load = ir.piecewise({10.0: a, 30.0: b})
    assert np.allclose(load(0.0), 0) and np.allclose(load(5.0), 0.5 * a)
    assert np.allclose(load(20.0), 2 * np.eye(3)) and np.allclose(load(1e9), b)


def test_equivalent_lid_interpolates_the_lid_curve():
    curve = [(50, -0.63), (100, -0.47), (300, -0.21)]
    te = ir.equivalent_te(np.array([-0.47, -0.34, -0.9]), curve)
    assert te[0] == pytest.approx(100.0) and 100 < te[1] < 300 and math.isnan(te[2])
