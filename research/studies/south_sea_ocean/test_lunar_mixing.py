import numpy as np
import pytest

from research.studies.south_sea_ocean import circulation as C
from research.studies.south_sea_ocean import lunar_mixing as M


def test_the_inertial_band_keeps_the_motion_turning_the_inertial_way():
    n, dt = M.INERTIAL_SAMPLES, C.YEAR_S / M.INERTIAL_SAMPLES
    f = np.array([[-2.0 * np.pi * 20 / (n * dt)]])          # a southern sea's f, on a frequency the year resolves
    t = np.arange(n)[:, None, None] * dt
    inertial = 0.05 * np.exp(-1j * f * t)                   # free inertial motion turns as exp(-i f t)
    contrary = 0.03 * np.exp(1j * f * t)                    # the same period turning the other way
    slow = 0.02 * np.exp(-1j * 0.3 * f * t)
    band = M.inertial_band(inertial + contrary + slow + 0.01, dt, f, M.INERTIAL_BAND)
    assert np.allclose(band, inertial, atol=1e-12)
    # the wind's work on that band: a stress turning with it does work, one turning against it does none
    stress = 0.01 * np.exp(-1j * f * t)
    assert (stress.real * band.real + stress.imag * band.imag).mean() == pytest.approx(0.01 * 0.05)


def test_each_tidal_line_is_fitted_apart():
    t = np.arange(6.5, 35.5, 1.0) * 86400.0
    month, fortnight = (2 * np.pi / (p * 86400.0) for p in M.TIDE_LINES_DAYS)
    u = 0.002 * np.cos(month * t) + 0.001 * np.sin(fortnight * t) + 0.0005
    v = 0.002 * np.sin(month * t)
    lines = M.line_mean_squares(t, u[:, None], v[:, None])
    assert lines[0][0] == pytest.approx(0.002 ** 2, rel=1e-9)         # half of each component's amplitude squared
    assert lines[1][0] == pytest.approx(0.5 * 0.001 ** 2, rel=1e-9)


def test_a_tidal_line_below_the_inertial_frequency_makes_no_internal_waves():
    month = 2 * np.pi / (27.32 * 86400.0)
    f = np.array([0.0, 0.5 * month, 1.5 * month])
    e = M.line_power(1025.0, 1.0e4, 1.0e-3, 1.0e-6, f, 27.32)
    assert e[0] == pytest.approx(0.5 * 1025.0 * M.KAPPA_PER_M * 1.0e4 * 1.0e-3 * 1.0e-6)
    assert e[1] == pytest.approx(e[0] * np.sqrt(0.75))
    assert e[2] == 0.0


def test_floor_roughness_is_the_spread_about_each_cells_plane():
    rng = np.random.default_rng(1)
    cell = np.repeat([0, 1, 2], 400)
    x, y = rng.uniform(-15e3, 15e3, cell.size), rng.uniform(-15e3, 15e3, cell.size)
    noise = rng.normal(0.0, 80.0, cell.size)
    z = 3000.0 + 0.02 * x - 0.05 * y + noise                 # a sloping floor with 80 m of roughness
    z[cell == 1] = 2000.0 + 0.1 * x[cell == 1]               # a smooth slope
    var = M.plane_residual_variance(cell, x, y, z, 4)
    assert var[0] == pytest.approx(np.var(noise[cell == 0]), rel=0.02)
    assert var[1] == pytest.approx(0.0, abs=1e-6)
    assert var[3] == 0.0                                      # an empty cell
    few = M.plane_residual_variance(np.zeros(5, int), x[:5], y[:5], z[:5], 1)
    assert few[0] == 0.0


def test_oxygen_scales_with_the_export_and_is_bounded_by_the_dye():
    used, new = np.array([1.0, 3.0, 9.0]), np.array([0.0, 0.5, 1.0])
    out = M.oxygen_at_exports(6.0, used, new, M.EXPORTS_G_C[0])
    reference = out[M.EXPORTS_G_C[0]]
    assert reference['high'].tolist() == [5.0, 3.0, 0.0]                    # none where more is spent than held
    assert reference['low'].tolist() == [5.0, 0.0, 0.0]                     # water not renewed since the start
    half = M.oxygen_at_exports(6.0, used, new, 2.0 * M.EXPORTS_G_C[0])[M.EXPORTS_G_C[0]]
    assert half['high'].tolist() == [5.5, 4.5, 1.5]


def test_the_paused_runs_reading_shows_a_renewed_top_and_a_still_layer():
    import json
    out = json.loads((M.HERE / 'results' / 'lunar_mixing_1deg.json').read_text())
    assert out['schema'] == 'terluna.research.south-sea-lunar-mixing/1' and out['model_years'] == len(out['years']) == 3
    assert out['setting']['tracers'] == list(M.RUN_TRACERS) and [y['year'] for y in out['years']] == [51.0, 52.0, 53.0]
    mid = np.array(out['level_mid_m'])
    at = lambda depth: int(np.abs(mid - depth).argmin())
    dye = np.array([[np.nan if v is None else v for v in y['new_water_by_depth']] for y in out['years']])
    used = np.array([[np.nan if v is None else v for v in y['used_by_depth_mg_l']] for y in out['years']])
    k = np.array([np.nan if v is None else v for v in out['years'][-1]['diffusivity_by_depth_m2_s']])
    for depth in (117.0, 179.0, 269.0):
        assert (np.diff(dye[:, at(depth)]) < 0).all()                       # the upper zone's water reaches the surface
    assert dye[-1, at(117.0)] < 0.4 < 0.8 < dye[-1, at(269.0)]
    assert (dye[-1, mid > 580.0][:6] > 0.9998).all()                        # under 0.02% from 585 m down in three years
    assert k[at(328.0)] > 10.0 * k[at(399.0)] and (k[(mid > 390.0) & (mid < 900.0)] < 1e-6).all()
    assert k[at(1493.0)] > 1e-3                                             # the deep pool, stirred within itself
    assert (np.diff(used[:, at(483.0)]) > 0).all() and used[-1, at(117.0)] < 0.5
    bounds = out['oxygen']['34 g C']
    assert bounds['high']['eels_layer_mg_l'] > 6.0 and bounds['low']['eels_layer_mg_l'] < 0.5
