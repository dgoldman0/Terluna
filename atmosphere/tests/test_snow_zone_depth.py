"""Snow-zone depth above the freezing level: a synthetic column with known answers."""
import numpy as np
import pytest

from atmosphere.snow import snow_zone_depth as zone
from atmosphere.snow import turbulence as t
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT


def synthetic(nz=60, dz=1000.):
    interfaces = np.arange(nz + 1) * dz
    z = (interfaces[:-1] + interfaces[1:]) / 2
    T = 290. - 1e-3 * z                      # crosses 273.15 K at 16.85 km
    p = 1e5 * np.exp(-z / 30000.)
    th = T * (t.CM1_P00 / p) ** (CM1_DRY_AIR_GAS_CONSTANT / t.CM1_CP)
    qs = np.where((z > 18000.) & (z < 30000.), 2e-5, 0.)
    qi = np.where((z > 25000.) & (z < 40000.), 2e-5, 0.)
    frame = dict(th=th[:, None], prs=p[:, None], qs=qs[:, None], qi=qi[:, None])
    return dict(z=z, dz=np.diff(interfaces)), frame


def test_freezing_level_snow_top_and_ice_top():
    item, frame = synthetic()
    rows = zone.column_measures(item, frame, 1e-5)
    assert rows.shape == (1, 7)
    z0, top, top_c, half, ninety, mass, ice_top = rows[0]
    assert z0 == pytest.approx(16850., abs=1.)
    assert top == pytest.approx(29500. - z0)
    assert top_c == pytest.approx(290. - 29.5 - 273.15)
    assert 18500. - z0 < half < 29500. - z0 and half < ninety
    assert ice_top == pytest.approx(39500. - z0)


def test_columns_without_snow_above_freezing_are_skipped():
    item, frame = synthetic()
    frame['qs'][:] = 0.
    assert zone.column_measures(item, frame, 1e-5).shape == (0, 7)


def test_weighted_quantiles_follow_the_weights():
    # Equal weights give cumulative fractions 1/3, 2/3 and 1: the median lies halfway from 1 to 2.
    assert zone.quantiles([1., 2., 3.], [1., 1., 1.])['p50'] == pytest.approx(1.5)
    assert zone.quantiles([3., 1., 2.])['p50'] == pytest.approx(2.)
    assert zone.quantiles([np.nan]) is None
