"""RRTMG's table edge: the held-edge edit, the twins inside the tables and the offline driver (no CM1 run)."""
import shutil

import numpy as np
import pytest

from climate.crm import rrtmg_edge as re


def column(ground_pa=121_400.0, layers=60):
    """A lunar column on fine layers near the ground, as the ring's are."""
    from atmosphere.radiative_convective import thermodynamics as td
    air = td.earthlike_air(121590.0, 400.0)
    p = ground_pa * np.exp(-np.linspace(0.0, 9.0, layers + 1) ** 1.5 / 3.0)
    t = np.maximum(300.0 - 6.0e-3 * 51_000.0 * np.log(p[0] / p), 190.0)
    x = np.maximum(0.02 * (p / p[0]) ** 4, 1e-6)
    return td.column_from_levels(td.MOON, air, p, t, x, constant_gravity=True), (p, t, x), air


def test_the_held_edge_changes_only_the_interpolation_fraction():
    old, new = re.EDGE_EDIT
    assert '5._rb *(preflog(jp(lay)) - plog)' in old and 'max(0._rb, 5._rb *(preflog(jp(lay)) - plog))' in new
    source = re.BUILD / re.MODULE
    if not source.exists():
        pytest.skip('the lunar CM1 build is not at hand')
    assert source.read_text().count(old) == 1


def test_the_twins_keep_every_layer_inside_rrtmgs_tables():
    col, levels, air = column()
    assert col['layer_p_pa'][0] / 100.0 > re.TABLE_EDGE_HPA          # the column's lowest layers lie beyond the edge
    twin = re.twins(dict(levels=[levels], air=air))[0]
    assert twin['p_pa'][0] / 100.0 == pytest.approx(re.TWIN_GROUND_HPA)
    assert (twin['layer_p_pa'] / 100.0 < re.TABLE_EDGE_HPA).all()
    assert np.allclose(twin['layer_t_k'], col['layer_t_k']) and np.allclose(twin['layer_x_h2o'], col['layer_x_h2o'])


def test_the_driver_reads_cm1s_rrtmg_tables_and_cools_the_air(tmp_path):
    if shutil.which('gfortran') is None or not (re.BUILD / 'module_ra_rrtmg_lw.o').exists():
        pytest.skip('no gfortran or no lunar CM1 build')
    from atmosphere.radiative_convective import thermodynamics as td
    col, _, _ = column()
    rate = re.rrtmg(re.build('as_built', tmp_path), [col], [300.0], 1004.64, td.MOON.surface_gravity, tmp_path)[0]
    assert np.isfinite(rate).all() and rate.size == col['layer_p_pa'].size
    assert rate[:20].mean() < 0.0                                     # the lower air loses heat to space
