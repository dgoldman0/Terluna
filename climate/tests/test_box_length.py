"""The boundary-layer length test's reading (climate/crm/box_length.py) on made-up CM1 files: the verdict set before the
run, the outputs of CM1's second output configuration (suffix .02, numbered from 1) found by their metadata beside a
source's copied metadata, the statistics joined over configurations, the rain by quarter from the running total, and
the cloud tops of a snapshot."""
import json

import numpy as np
import pytest

from climate.crm import box_length as bln
from climate.crm import box_lid as bl


def test_the_verdict_follows_the_criteria_set_before_the_run():
    assert bln.verdict(6.2, 0.5) == 'robust'                       # within 1 km and 1 g/kg
    assert bln.verdict(5.48, 1.35) == 'robust'                     # on both edges
    assert bln.verdict(4.9, 0.4) == 'dependent on the length'      # the depth falls by more than 1.5 km
    assert bln.verdict(6.5, 1.5) == 'dependent on the length'      # the drop grows by more than 1 g/kg
    assert bln.verdict(5.2, 0.4) == 'ambiguous'                    # a fall of 1.28 km
    assert bln.verdict(7.6, 0.3) == 'ambiguous'                    # deeper by more than 1 km


def write_ctl(path, names, nx=1):
    path.write_text(f'xdef {nx} linear 0 1\nydef 1 linear 0 1\nzdef 1 levels\n1\nvars {len(names)}\n'
                    + ''.join(f'{n} 0 99 {n}\n' for n in names) + 'endvars\n')


def test_a_second_output_configuration_is_found_by_its_metadata_beside_the_sources_copy(tmp_path):
    names = ['mtime', 'nstep', 'dt', 'adt', 'nwrite']
    for suffix, rows in (('', [(0.0, 1.0), (10800.0, 2.0)]),                 # the source's, its files not held
                         ('.02', [(21600.0, 1.0), (43200.0, 2.0)])):
        write_ctl(tmp_path / f'cm1out_metadata{suffix}.ctl', names)
        np.array([[t, 0.0, 1.0, 1.0, n] for t, n in rows], '<f4').tofile(tmp_path / f'cm1out_metadata{suffix}.dat')
    write_ctl(tmp_path / 'cm1out_s.02.ctl', ['qc', 'qi'], nx=2)
    for n in (1, 2):
        np.arange(4, dtype='<f4').tofile(tmp_path / f'cm1out_t{n:06d}_s.02.dat')
    found = bln.output_times(tmp_path)
    assert sorted(found) == [21600, 43200] and found[43200][0].name == 'cm1out_t000002_s.02.dat'
    d = bln.read_snapshot(*found[21600])
    assert d['qc'].tolist() == [0.0, 1.0] and d['qi'].tolist() == [2.0, 3.0]
    np.arange(5, dtype='<f4').tofile(tmp_path / 'cm1out_t000002_s.02.dat')
    with pytest.raises(ValueError):
        bln.read_snapshot(*found[43200])


def test_the_statistics_join_over_configurations_and_give_the_rain_by_quarter(tmp_path, monkeypatch):
    monkeypatch.setattr(bln.ra, 'RUNS', tmp_path)
    case = tmp_path / 'box'
    case.mkdir()
    (case / 'case.json').write_text(json.dumps(dict(configuration=dict(start_hour_angle_deg=0.0, site=dict(lon_deg=0.0)),
                                                    grid=dict(nx=1, ny=1, dx_m=1.0))))
    names = ['mtime', 'train']
    day = 4.0 * 86400.0                                                # a made-up solar day of four Earth days
    monkeypatch.setattr('climate.crm.cm1_run.solar_day_s', lambda: day)
    write_ctl(case / 'cm1out_stats.ctl', names)                        # the source's layout, its file not held
    write_ctl(case / 'cm1out_stats.02.ctl', names)
    t = np.arange(0.0, 4.01, 0.25) * 86400.0                            # every 6 hours over one solar day
    rate = np.where((t > 0.0) & (t <= 86400.0), 2.0, 0.0)              # 2 mm/day over the first quarter alone
    train = np.cumsum(rate * 0.25)
    rows = np.vstack([np.r_[0.0, t], np.r_[0.0, train]]).T             # a leading record the restart left unwritten
    rows.astype('<f4').tofile(case / 'cm1out_stats.02.dat')
    s = bln.stats_rows(case)
    assert s['mtime'].size == t.size - 1 and s['mtime'][0] == t[1]      # unwritten records (time nil) are gone
    assert (np.diff(s['mtime']) > 0).all()
    rain = bln.rain_by_quarter('box', (0.0, 4.0))
    assert rain['noon-sunset'] == pytest.approx(2.0) and rain['sunset-midnight'] == pytest.approx(0.0)
    assert rain['midnight-sunrise'] == pytest.approx(0.0) and rain['sunrise-noon'] == pytest.approx(0.0)


def test_the_cloud_tops_count_each_column_s_highest_cloud_below_40_km():
    zh = np.array([1000.0, 5000.0, 9000.0, 20000.0, 45000.0])
    qc, qi = np.zeros((5, 4)), np.zeros((5, 4))
    qc[0:2, 0] = 2e-5                                                   # cloud to 5 km
    qc[1:3, 1] = 2e-5                                                   # to 9 km
    qi[3, 2] = 2e-5                                                     # at 20 km alone
    qi[4, 3] = 2e-5                                                     # above 40 km only: no cloud counted below it
    out = bln.cloud_tops(qc, qi, zh)
    assert out['cloudy_share'] == pytest.approx(0.75) and out['top_median_m'] == 9000.0
    assert out['above_16km_share'] == pytest.approx(0.25)


def test_the_diagnostics_of_both_configurations_are_read_each_with_its_own_layout(tmp_path):
    for suffix in ('', '.02'):
        for k in 'sw':
            write_ctl(tmp_path / f'cm1out_diag_{k}{suffix}.ctl', ['mtime', 'x'])
    for name in ('cm1out_diag_000474_s.dat', 'cm1out_diag_000474_w.dat', 'cm1out_diag_000001_s.02.dat',
                 'cm1out_diag_000001_w.02.dat'):
        np.zeros(2, '<f4').tofile(tmp_path / name)
    files = bl.diag_files(tmp_path)
    assert sorted((s.name, w.name) for s, w, *_ in files) == [
        ('cm1out_diag_000001_s.02.dat', 'cm1out_diag_000001_w.02.dat'),
        ('cm1out_diag_000474_s.dat', 'cm1out_diag_000474_w.dat')]
    assert all(lay_s == [('mtime', 1), ('x', 1)] for _, _, lay_s, _ in files)
