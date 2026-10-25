"""McSnow aggregation runs: observed relations, run configuration and declared source patches."""
import math

import numpy as np
import pytest

from atmosphere.snow import mcsnow_aggregation as agg


def test_observed_median_volume_diameters_follow_the_exponential_fits():
    d0 = agg.observed_d0_mm(1.)
    assert d0['gunn_marshall_1958'] == pytest.approx(3.67 / 2.55)
    assert d0['sekhon_srivastava_1970'] == pytest.approx(3.67 / 2.29)
    assert agg.observed_d0_mm(3.)['gunn_marshall_1958'] == pytest.approx(3.67 / (2.55 * 3 ** -.48))


def test_weighted_quantile_matches_a_hand_case():
    # Flux weights are always positive (multiplicity x mass x speed), so the cumulative weight rises strictly.
    values = np.array([1., 2., 3., 4.])
    assert agg.weighted_quantile(values, np.array([1., 1., 1., 1.]), .5) == pytest.approx(2.)
    # Cumulative weights 1, 2, 3, 8: the median weight 4 lies a fifth of the way from 3 to 8.
    assert agg.weighted_quantile(values, np.array([1., 1., 1., 5.]), .5) == pytest.approx(3.2)


def test_configuration_delivers_the_requested_water_through_integer_superparticles():
    cfg = agg.configuration('earth', monomer_mass=agg.EARTH_MONOMER_MASS, water_mm_h=1., stick=2)
    sp, grid, atmo = cfg['sp_nml'], cfg['grid_nml'], cfg['atmo_nml']
    per_step = math.floor(sp['nclrate'] * sp['dtc'] * (atmo['nh2'] - atmo['nh1']) * grid['box_area'])
    flux = per_step * sp['xi0'] / (sp['dtc'] * grid['box_area'])
    assert flux * sp['nclmass'] * 3600 == pytest.approx(1., rel=.03)
    assert atmo['ssat'] == 0. and atmo['h2'] > atmo['h1'] and atmo['nh2'] < grid['dom_top']
    assert sp['sp_kernel_id'] == 2 and sp['sp_habit'] == 0 and cfg['general_nml']['agggeo_type'] == 2


def test_source_patches_change_only_what_they_declare():
    if not (agg.SOURCE / 'src/mo_sp_nml.f90').exists():
        pytest.skip('pinned McSnow source copy is not installed')
    earth = agg.patched_sources('earth')
    assert set(earth) == {'mo_sp_nml.f90', 'mo_velocity.f90', 'mo_mass2diam.f90'}
    assert '  vt = vt * band_speed\nEND FUNCTION vterm_snow' in earth['mo_velocity.f90']
    assert 'mono_alf       = 0.7453_wp' in earth['mo_sp_nml.f90']
    moon = agg.patched_sources('moon')
    assert 'grav  = 1.62421887656' in moon['mo_atmo_types.f90']
