"""Growth and aggregation together: declared patches, column and parcel configurations, forcing tables."""
import json
import math

import numpy as np
import pytest

from atmosphere.snow import mcsnow_aggregation as agg
from atmosphere.snow import mcsnow_growth_aggregation as grow


def test_seed_is_a_four_micrometre_ice_sphere_of_mcsnow_density():
    assert grow.SEED_MASS == pytest.approx(919. * math.pi / 6 * (4e-6) ** 3)


def test_water_limited_column_spreads_the_supply_and_saturates_over_water():
    cfg = grow.column_config(3., .1)
    sp, grid, atmo = cfg['sp_nml'], cfg['grid_nml'], cfg['atmo_nml']
    assert sp['water_supply'] == pytest.approx(3. / 3600.)
    assert atmo['h1'] == 0. and atmo['h2'] == grid['dom_top'] and atmo['IGF_id'] == 2
    assert sp['sp_habit'] == 1 and sp['plate_base'] and not sp['box_array'] and sp['sp_kernel_id'] == 2
    per_step = math.floor(sp['nclrate'] * sp['dtc'] * (atmo['nh2'] - atmo['nh1']) * grid['box_area'])
    flux = per_step * sp['xi0'] / (sp['dtc'] * grid['box_area'])
    assert flux == pytest.approx(.1 * grow.REFERENCE_FLUX, rel=1e-3)
    assert cfg['general_nml']['dt_1dprof'] == cfg['general_nml']['time_end']
    saturated = grow.column_config(0., 1., saturated=True)
    assert saturated['sp_nml']['water_supply'] < 0


def test_reference_column_keeps_the_aggregation_runs_depth_and_the_deep_column_reaches_minus_30_c():
    reference = grow.column_config(3., 1.)['grid_nml']
    assert reference['dom_top'] == agg.configuration('moon', monomer_mass=1e-8, water_mm_h=1.)['grid_nml']['dom_top']
    deep = grow.column_config(.3, 1., top_C=-30.)
    assert deep['grid_nml']['dom_top'] == pytest.approx(30. / agg.lunar_column()['tlr'], abs=20.)
    with pytest.raises(ValueError):
        grow.column_config(.3, 1., top_C=-35.)


def test_source_patches_apply_and_leave_every_new_option_off_by_default():
    if not (grow.SOURCE / 'src/mo_sp_nml.f90').exists():
        pytest.skip('pinned McSnow source copy is not installed')
    src = grow.patched_sources('earth')
    assert {'mo_depo.f90', 'mo_atmo.f90', 'mo_loop_array.f90', 'mo_nucleation.f90'} <= set(src)
    nml = src['mo_sp_nml.f90']
    for default in ('water_supply   = -1.0_wp', 'depo_max_fraction = 0.0_wp', 'band_mass      = 1.0_wp',
                    'plate_base     = .FALSE.', 'box_array      = .FALSE.', 'dt_snap        = 0.0_wp'):
        assert default in nml
    depo = src['mo_depo.f90']
    assert 'IF ( .NOT. plate_base ) fv      = fv + c2 / phi * X_equiv**(3./2.)' in depo
    assert 'C = capacitance( sp, Re ) * band_mass' in depo
    assert 'IF ( depo_max_fraction <= 0._wp ) THEN\n    CALL sp_depo(' in depo
    loop = src['mo_loop_array.f90']
    assert 'IF ( .NOT. box_array .OR. box_settle ) CALL sp_sedi( sp, atmo, dtc )' in loop
    assert 'IF ( water_supply >= 0.0_wp .AND. diffusion_type == 1 ) THEN' in loop


def test_band_factors_are_the_error_band_edges():
    if not grow.EARTH_BAND.exists():
        pytest.skip('error band product absent')
    band = json.loads(grow.EARTH_BAND.read_text())['column']['allowance_F']
    assert grow.band_factors(None) == (1., 1., 1.)
    assert grow.band_factors('combined_small') == tuple(band['combined_small']['factors'])


def test_path_columns_seed_each_release_cell_once_and_stop_at_the_window_end():
    cfg = grow.path_column_config(3, [3600., 7200., 1800.], 1e5)
    sp, grid, atmo = cfg['sp_nml'], cfg['grid_nml'], cfg['atmo_nml']
    assert grid['nz'] == 3 * grow.N_SUB and grid['dom_top'] == 3 * grow.N_SUB * grow.COLUMN_DZ
    per_step = math.floor(sp['nclrate'] * sp['dtc'] * (atmo['nh2'] - atmo['nh1']) * grid['box_area'])
    assert per_step == 3 * 64 and sp['ncl_until'] == sp['dtc']
    assert sp['box_array'] and sp['box_settle'] and sp['box_sub'] == grow.N_SUB and sp['skip_bulk']
    assert sp['depo_residual'] == pytest.approx(grow.RESIDUAL_MASS) and sp['depo_censor'] == grow.CENSOR_M
    assert sp['xi0'] * 64 / (grow.COLUMN_DZ * grow.PATH_AREA) == pytest.approx(1e5, rel=1e-6)
    assert cfg['general_nml']['time_end'] >= 7200. and cfg['general_nml']['time_end'] % 600 == 0


def test_residual_and_censor_are_off_unless_set():
    if not (grow.SOURCE / 'src/mo_sp_nml.f90').exists():
        pytest.skip('pinned McSnow source copy is not installed')
    src = grow.patched_sources('moon')
    assert 'depo_residual  = 0.0_wp' in src['mo_sp_nml.f90'] and 'depo_censor    = 0.0_wp' in src['mo_sp_nml.f90']
    assert 'IF ( m1 < depo_residual ) THEN' in src['mo_depo.f90']
    assert "IF ( depo_censor > 0._wp .AND. sp%mm < 1.01_wp ) THEN" in src['mo_depo.f90']
    assert 'IF ( .NOT. skip_bulk ) THEN' in src['mo_loop_array.f90']
    column = grow.column_config(3., 1.)['sp_nml']
    assert 'depo_residual' not in column and 'depo_censor' not in column and 'box_sub' not in column
