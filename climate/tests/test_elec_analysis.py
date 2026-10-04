"""Checks of the electrified storms' analysis: logs taken up again after a restart, flashes summed by kind and
sub-step, WRF-ELEC's charging totals read from CM1's log by step, and the charge regions found in height and
temperature."""
import numpy as np
import pytest

from climate.crm import elec_analysis as ea


def test_a_restart_replaces_the_rows_written_after_the_time_it_takes_up_from(tmp_path):
    log = tmp_path / 'terluna_field.txt'
    rows = lambda times: ''.join(f'{t:.1f} {int(t / 6)} 1.0e3 0 0 0 1.0 -1.0 5.0 1e-10 1e-11\n' for t in times)
    log.write_text('# time_s step e_max_V_m x_m y_m z_m positive_C negative_C energy_J noninductive_max_C_m3_s '
                   'inductive_max_C_m3_s\n'
                   '# run from 0.0\n' + rows([60, 120, 180, 240]) + '# run from 120.0\n' + rows([180, 240, 300]))
    out = ea.read_log(log, ea.FIELD_COLUMNS)
    assert list(out['time_s']) == [60, 120, 180, 240, 300] and list(out['step']) == [10, 20, 30, 40, 50]
    assert ea.read_log(tmp_path / 'missing.txt', ea.FIELD_COLUMNS)['time_s'].size == 0


def test_flashes_are_summed_by_kind_and_the_field_is_read_before_each_sub_step_s_first(tmp_path):
    log = tmp_path / 'terluna_flashes.txt'
    row = ('{t:.1f} {step} {sub} {kind} 1000.0 2000.0 8000.0 1.2e5 1.0e5 {over} 3.0e7 1 {qp} {qn} 2.0e9 1.0e9 '
           '5000.0 11000.0 40 {nox}\n')
    log.write_text(
        '# time_s step substep kind x_m y_m z_m e_V_m e_break_V_m e_over_break_max extent_m2 regions positive_C '
        'negative_C energy_before_J energy_after_J z_low_m z_high_m points nox_mol\n# run from 0.0\n' +
        row.format(t=60.0, step=10, sub=1, kind=1, over=2.0, qp=10.0, qn=10.0, nox=5.0) +
        row.format(t=60.0, step=10, sub=1, kind=2, over=1.5, qp=0.0, qn=20.0, nox=8.0) +
        row.format(t=60.5, step=10, sub=2, kind=1, over=1.2, qp=12.0, qn=12.0, nox=6.0) +
        row.format(t=66.0, step=11, sub=1, kind=3, over=1.1, qp=30.0, qn=0.0, nox=0.0))
    f = ea.read_log(log, ea.FLASH_COLUMNS)
    assert list(f['kind']) == [1, 2, 1, 3] and list(f['nox_mol']) == [5, 8, 6, 0]
    edges = np.arange(0.0, 400.0, 300.0)
    ic = ea.flash_summary({k: v[f['kind'] == 1] for k, v in f.items()}, edges)
    assert ic['count'] == 2 and ic['per_5_min'] == [2.0] and ic['nox_total_mol'] == 11.0
    assert ic['positive_c']['mean'] == 11.0 and ic['energy_dissipated_j']['max'] == 1.0e9 and 'regions' not in ic
    first = np.r_[True, (np.diff(f['step']) != 0) | (np.diff(f['substep']) != 0)]
    assert list(first) == [True, False, True, True]


def test_the_ice_above_the_radiation_s_cap_is_summed_by_mass(tmp_path):
    log = tmp_path / 'terluna_ice_optics.txt'
    log.write_text('# header\n# run from 0.0\n'
                   '800.0 10 0.0 0.0 0.0 0.0 0\n'
                   '1600.0 20 1.0e6 0.0 0.0 0.0 0\n'
                   '2400.0 30 2.0e6 5.0e5 175.0 190.0 12\n'
                   '3200.0 40 1.0e6 5.0e5 154.0 160.0 4\n')
    s = ea.ice_cap_summary(ea.read_log(log, ea.ICE_COLUMNS))
    assert s['rows'] == 3 and s['share_above_max'] == pytest.approx(0.5) and s['share_above_overall'] == pytest.approx(0.25)
    assert s['mean_radius_above_um'] == pytest.approx(164.5) and s['max_radius_um'] == 190.0
    assert s['optical_depth_overstated'] == pytest.approx(164.5 / 140.0)
    assert ea.ice_cap_summary(ea.read_log(tmp_path / 'missing.txt', ea.ICE_COLUMNS)) is None


def test_charging_totals_go_to_the_step_whose_time_line_follows_them(tmp_path):
    block = ('Integrated pos/neg charging rates:\n'
             'ctswin,ctswip = -1.00000E-02,  2.00000E-02\n'
             'ctghsn,ctghsp = -3.00000E+00,  4.00000E+00\n'
             'ctghin,ctghip = -5.00000E+00,  6.00000E+00\n'
             'ctghwn,ctghwp = {w:.5E},  0.00000E+00\n')
    (tmp_path / 'cm1.log').write_text(block.format(w=-1.0) + '   turbfac =    3.3E-03\n'
                                      '             1              0.100000 min \n' + block.format(w=-2.0) +
                                      '             2              0.200000 min \n')
    out = ea.charging(tmp_path)
    assert list(out['step']) == [1, 2] and list(out['time_s']) == pytest.approx([6.0, 12.0])
    assert out['ctghi'][0].tolist() == [-5.0, 6.0] and out['ctghw'][:, 0].tolist() == [-1.0, -2.0]
    assert np.isnan(out['cthi']).all()


def test_the_main_negative_region_and_the_positive_ones_around_it():
    zh = np.arange(0.5, 15.0, 1.0)                                       # km
    t = 300.0 - 6.5 * zh[:, None] * np.ones((1, 4))                      # K
    net = np.zeros((zh.size, 4))
    net[7, 1] = -2.0e-9                                                  # 7.5 km, -18.75 C
    net[10, 1] = 1.0e-9                                                  # 10.5 km
    net[3, 2] = 0.5e-9                                                   # 3.5 km
    graupel = np.where(net < 0, net, 0.0)
    s = dict(net=net, t=t, volume=np.full(net.shape, 1.0e9), charges=dict(graupel=graupel, cloud_ice=net - graupel))
    out = ea.structure(s, zh)
    assert out['main_negative']['z_km'] == 7.5 and out['main_negative']['t_c'] == pytest.approx(300.0 - 6.5 * 7.5 - 273.15)
    assert out['main_negative']['charge_c'] == pytest.approx(-2.0)
    assert out['main_negative']['carried_c']['graupel'] == pytest.approx(-2.0)
    assert out['upper_positive']['z_km'] == 10.5 and out['lower_positive']['z_km'] == 3.5
    assert out['positive_c'] == pytest.approx(1.5) and out['min_nc_m3'] == pytest.approx(-2.0)
    assert ea.structure(dict(s, net=np.abs(net)), zh).get('main_negative') is None
