"""Checks of the electrified storms' analysis: logs taken up again after a restart, WRF-ELEC's charging totals read
from CM1's log by step, and the charge regions found in height and temperature."""
import numpy as np
import pytest

from climate.crm import elec_analysis as ea


def test_a_restart_replaces_the_rows_written_after_the_time_it_takes_up_from(tmp_path):
    log = tmp_path / 'terluna_field.txt'
    rows = lambda times: ''.join(f'{t:.1f} {int(t / 6)} 1.0e3 0 0 0 1.0 -1.0 5.0 0 1e-10 1e-11\n' for t in times)
    log.write_text('# time_s step e_max_V_m x_m y_m z_m positive_C negative_C energy_J discharges\n'
                   '# run from 0.0\n' + rows([60, 120, 180, 240]) + '# run from 120.0\n' + rows([180, 240, 300]))
    out = ea.read_log(log, ea.FIELD_COLUMNS)
    assert list(out['time_s']) == [60, 120, 180, 240, 300] and list(out['step']) == [10, 20, 30, 40, 50]
    assert ea.read_log(tmp_path / 'missing.txt', ea.FIELD_COLUMNS)['time_s'].size == 0


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
