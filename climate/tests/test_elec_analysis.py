"""Checks of the electrified storms' analysis: logs taken up again after a restart, flashes summed by kind and
sub-step, WRF-ELEC's charging totals read from CM1's log by step, the charge regions found in height and temperature,
and a run summarized over a window of its days."""
import json

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


def test_the_field_at_the_ground_and_its_point_discharge_are_read_by_interval(tmp_path):
    log = tmp_path / 'terluna_ground.txt'
    log.write_text('# header\n# run from 0.0\n'
                   '      400.0       10  2.0000E+00 -1.0000E+00    12000.0     6000.0      30       4\n'
                   '      800.0       20  5.0000E+00  0.0000E+00     4000.0     3000.0       2       0\n')
    s = ea.ground_summary(ea.read_log(log, ea.GROUND_COLUMNS))
    assert s['positive_c'] == pytest.approx(7.0) and s['negative_c'] == pytest.approx(-1.0) and s['columns_max'] == 30
    assert s['e_ground_max_kv_m'] == pytest.approx(12.0) and s['e_ground_median_kv_m'] == pytest.approx(8.0)
    assert s['e_ground_dry_max_kv_m'] == pytest.approx(6.0) and s['columns_wet_max'] == 4
    assert ea.ground_summary(ea.read_log(tmp_path / 'missing.txt', ea.GROUND_COLUMNS)) is None


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


def test_a_window_takes_only_its_own_flashes_charging_and_segments(tmp_path, monkeypatch):
    """A run summarized over a window of its days: the segments overlapping it, its flashes by hour and kind, and the
    charge each collision type separated over it, the non-inductive pairs counted once."""
    case = tmp_path / 'run'
    case.mkdir()
    monkeypatch.setattr(ea.ra, 'RUNS', tmp_path)
    (case / 'case.json').write_text(json.dumps(dict(
        configuration=dict(output_s=10800.0), grid=dict(dx_m=1000.0), build=dict(executable_sha256='ab' * 32),
        electricity=dict(isaund=1, lightning=4, ground_m=-1.0, corona_v_m=0.0, leakage=0, substep_s=6.8))))
    (case / 'cm1out_s.ctl').write_text('zdef 2 linear 0.5 1.0\n')
    (case / 'input_grid_z').write_text('0.0\n1000.0\n2000.0\n')
    (case / 'progress.json').write_text(json.dumps(dict(start_s=3600.0, segments=[
        dict(model_s=7200.0, log='cm1_segment_001.log', executable='e1'),
        dict(model_s=10800.0, log='cm1_segment_002.log', executable='e2'),
        dict(model_s=14400.0, log='cm1_segment_003.log', executable='e2')])))
    block = ('ctghsn,ctghsp = -1.00000E+00,  2.00000E+00\n'
             'ctgsn,ctgsp = -1.00000E+00,  2.00000E+00\n'
             'ctghwn,ctghwp = -5.00000E-01,  0.00000E+00\n')
    for k, steps in ((1, (1, 2)), (2, (3, 4)), (3, (5, 6))):
        (case / f'cm1_segment_00{k}.log').write_text(''.join(
            block + f'             {s}              {3600.0 + 1800.0 * s:.6f} sec \n' for s in steps))
    row = '{t:.1f} {s} 1 {kind} 0 0 8000.0 1e5 1e5 1.1 1e7 1 {qp} {qn} 2e9 1e9 5000.0 9000.0 4 1.0\n'
    (case / 'terluna_flashes.txt').write_text('# header\n# run from 3600.0\n' + ''.join(
        row.format(t=t, s=int(t // 60), kind=kind, qp=qp, qn=qn) for t, kind, qp, qn in (
            (6000.0, 1, 10.0, 10.0), (8000.0, 1, 20.0, 20.0), (8100.0, 2, 0.0, 50.0), (12000.0, 1, 5.0, 5.0))))
    (case / 'terluna_field.txt').write_text('# header\n# run from 3600.0\n' + ''.join(
        f'{t:.1f} {int(t // 60)} {e:.1f} 0 0 0 10.0 -8.0 1e9 1e-12 1e-13\n' for t, e in ((7500.0, 1.2e5), (9000.0, 0.9e5))))
    (case / 'terluna_ground.txt').write_text('# header\n# run from 3600.0\n' + ''.join(
        f'{t:.1f} {int(t // 60)} 0.0 0.0 {e:.1f} {e:.1f} 0 0\n' for t, e in ((7000.0, 9.0e3), (8000.0, 4.0e3), (9000.0, 2.0e3))))
    assert [s['executable'] for s in ea.window_segments(case, 7200.0, 10800.0)] == ['e2']
    w = ea.window('run', 7200.0, 10800.0)
    assert w['law'] == 'takahashi' and w['executables'] == ['e2'] and w['snapshots'] == []
    assert w['lightning']['count'] == 2 and w['lightning']['hours_with_flashes'] == 1
    assert w['lightning']['first_h'] == pytest.approx(800.0 / 3600.0)
    assert w['lightning']['in_cloud']['count'] == 1 and w['lightning']['negative_to_ground']['negative_c']['max'] == 50.0
    # steps 3 and 4 (9000 and 10800 s) are inside: the first counts no time, the second 1800 s
    assert w['charge']['noninductive_c'] == pytest.approx(3.0 * 1800.0)                   # ctghs only, not ctgs again
    assert w['charge']['inductive_c'] == pytest.approx(0.5 * 1800.0)
    assert w['field']['e_max_kv_m'] == pytest.approx(120.0) and w['field']['negative_c_max'] == pytest.approx(8.0)
    assert w['ground']['rows'] == 2 and w['ground']['e_ground_dry_max_kv_m'] == pytest.approx(4.0)
