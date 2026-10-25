"""The resolution test's reading of the held-night boxes (box_resolution.py) on made-up outputs: the highest cloud from
CM1's fractions of columns, the rain from its running total, a day's means, and the control's check and the reading set
before the runs."""
import pytest

from climate.crm import box_resolution as br


def test_the_highest_cloud_is_the_highest_level_any_column_holds_water_or_ice():
    zh = [500.0, 5000.0, 20000.0, 45000.0]
    assert br.cloud_top_m([0.1, 0.02, 0.0, 0.0], [0.0, 0.0, 0.001, 0.0], zh) == 20000.0
    assert br.cloud_top_m([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], zh) is None


def test_cloud_aloft_that_stands_apart_is_told_from_deep_cloud():
    import numpy as np
    zh = np.array([2000.0, 8000.0, 20000.0, 30000.0, 50000.0])
    qc = np.zeros((5, 3))
    qi = np.zeros((5, 3))
    qc[:2, 0] = 1e-4                                   # a shallow cloud to 8 km
    qi[:, 1] = 1e-4                                    # a deep column, cloud from 2 to 50 km
    qc[0, 2], qi[4, 2] = 1e-4, 1e-4                    # low cloud, and ice at 50 km with nothing between
    out = br.cloud_layers(qc, qi, zh)
    assert out['top_below_40km_highest_m'] == 30000.0 and out['top_below_40km_mean_m'] == pytest.approx(40000.0 / 3)
    assert out['cloud_16_40km_share'] == pytest.approx(1 / 3) and out['apart_above_40km_share'] == pytest.approx(1 / 3)


def test_the_rain_of_a_span_comes_from_the_running_total():
    rows = [dict(hour=h, rain_total_mm=total) for h, total in ((0.0, 0.0), (12.0, 0.1), (24.0, 0.2), (48.0, 1.2))]
    assert br.rain_over(rows, 0.0, 24.0) == pytest.approx(0.2)                # mm in a day
    assert br.rain_over(rows, 24.0, 48.0) == pytest.approx(1.0)
    assert br.rain_over(rows, 6.0, 18.0) == pytest.approx(0.2)                 # 0.1 mm in half a day, linear between
    assert br.rain_over(rows, 24.0, 60.0) is None


def made_up_day(rain, w, p999, deep, hours=(24.0, 48.0)):
    """Half-hourly diagnostics and 6-hourly snapshots over two days, the second day's as given."""
    rows, snaps, total = [], [], 0.0
    for k in range(97):
        h = 0.5 * k
        late = h > 24.0
        total += (rain if late else 0.0) / 48.0 if k else 0.0
        rows.append(dict(hour=h, rain_total_mm=total, rain_mm_day=rain if late else 0.0, cm1_cape_j_kg=1300.0,
                         cm1_cin_j_kg=0.5, hpbl_m=6000.0, w_max_below_16km_m_s=w, w_max_below_60km_m_s=w + 0.5 * late,
                         cloud_top_m=8000.0 if k % 2 else None, parcel_cape_j_kg=1400.0, parcel_cin_j_kg=0.4,
                         parcel_lfc_m=6200.0, parcel_cape_eps8_j_kg=10.0))
    for k in range(9):
        h = 6.0 * k
        late = h > 24.0
        snaps.append(dict(hour=h, w_max_m_s=w, w_max_height_m=5000.0, w_p999_below_16km_m_s=p999 if late else 0.5,
                          w_p999_16_60km_m_s=0.3, cloud_top_max_m=30000.0 if late else None,
                          cloud_tops_share=dict(below_16km=0.2, from_16_to_40km=0.01, above_40km=deep if late else 0.0),
                          cloud_below_16km_share=0.2, columns=dict(cape_mean_j_kg=1300.0, hpbl_mean_m=6000.0),
                          clouds=dict(top_below_40km_highest_m=9000.0 if late else None, top_below_40km_mean_m=7000.0,
                                      cloud_16_40km_share=0.0, apart_above_40km_share=deep if late else 0.0)))
    return br.day_summary(rows, snaps, *hours)


def test_a_day_takes_its_rain_from_the_total_and_its_means_after_its_start():
    day = made_up_day(rain=0.24, w=1.5, p999=1.1, deep=0.01)
    assert day['outputs'] == 48 and day['snapshots'] == 4 and day['rain_mm_day'] == pytest.approx(0.24)
    assert day['w_max_below_60km_m_s'] == pytest.approx(2.0) and day['w_max_below_16km_m_s'] == pytest.approx(1.5)
    assert day['cloud_top_m'] == 8000.0 and day['hpbl_m'] == 6000.0
    assert day['snapshot']['w_p999_below_16km_m_s'] == pytest.approx(1.1)
    assert day['snapshot']['cloud_tops_share']['above_40km'] == pytest.approx(0.01)
    assert day['snapshot']['cloud_top_below_40km_highest_m'] == 9000.0
    first = made_up_day(rain=0.24, w=1.5, p999=1.1, deep=0.01, hours=(0.0, 24.0))
    assert first['rain_mm_day'] == pytest.approx(0.0) and first['snapshot']['w_p999_below_16km_m_s'] == 0.5


def test_the_control_check_and_the_reading_follow_the_rules_set_before_the_runs():
    night = made_up_day(rain=0.1, w=1.5, p999=1.1, deep=0.01)
    check = br.quiet(night)
    assert check['quiet'] and not br.convects(night)
    assert check['rain_mm_day'] == dict(value=pytest.approx(0.1), below=0.3, holds=True)
    storms = made_up_day(rain=1.2, w=9.0, p999=3.0, deep=0.3)
    assert not br.quiet(storms)['quiet'] and br.convects(storms)
    drizzle = made_up_day(rain=0.6, w=2.5, p999=1.5, deep=0.0)                  # wet but shallow: neither
    assert not br.quiet(drizzle)['quiet'] and not br.convects(drizzle)
    assert br.reading(True, {'3 km': (False, True), '2 km': (False, True)}).startswith('storms start')
    assert br.reading(True, {'3 km': (True, False), '2 km': (True, False)}).startswith('none of the three')
    assert br.reading(True, {'3 km': (True, False), '2 km': (False, True)}) == \
        'in between: 3 km stays quiet; 2 km convects'
    assert br.reading(False, {'3 km': (False, True), '2 km': (False, True)}).startswith('the control does not')


def test_the_refill_counts_the_microphysics_vapour_where_the_air_holds_no_water(tmp_path, monkeypatch):
    """Two levels a day apart: the lower holds vapour, the upper none, and there the microphysics puts back 0.2 mm/day
    in its 500 kg/m2; the lower level's microphysics term is not counted."""
    import json
    import numpy as np
    monkeypatch.setattr(br.ra, 'RUNS', tmp_path)
    case = tmp_path / 'box'
    case.mkdir()
    (case / 'case.json').write_text(json.dumps(dict(grid=dict(nx=1, dx_m=1000.0, length_m=1000.0),
                                                    configuration=dict(start_hour_angle_deg=0.0))))
    (case / 'input_grid_z').write_text('0\n1000\n2000\n')
    (case / 'terluna_surface.txt').write_text('1\n0 1000 2 16 300 300\n')
    names = ('qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'rho', 'qvb_mp')
    (case / 'cm1out_diag_s.ctl').write_text('xdef 1 linear 0 1\nydef 1 linear 0 1\nzdef 2 levels\n500\n1500\n'
                                            f'vars {1 + len(names)}\nmtime 0 99 t\n'
                                            + ''.join(f'{n} 2 99 x\n' for n in names) + 'endvars\n')
    per = 1.0 / (br.DAY * 500.0)                                          # 1 mm/day in the upper level's 500 kg/m2
    for n, day in enumerate((1.0, 2.0, 3.0), 1):
        fields = dict(qv=(0.01, 0.0), qc=(1e-4, 0.0), qr=(0.0, 0.0), qi=(0.0, 0.0), qs=(0.0, 0.0), qg=(0.0, 0.0),
                      rho=(1.0, 0.5), qvb_mp=(-3.0 * per, 0.2 * per))
        np.array([day * br.DAY] + [v for name in names for v in fields[name]], dtype='<f4').tofile(
            case / f'cm1out_diag_{n:06d}_s.dat')
    assert br.refilled('box', 1.0, 2.0) == pytest.approx(0.2, rel=1e-5)
