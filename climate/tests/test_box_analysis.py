"""Checks of the box comparison's pure pieces: the air as the GCM measures it and the variables it adds."""
import numpy as np
import pytest

from climate.crm import box_analysis as ba


def column(theta_of_p, pressures, surface=1.0e5):
    p = np.asarray(pressures, dtype=float)[:, None]
    return dict(prs=p, th=theta_of_p(p), psfc=np.array([surface]))


def test_the_gcm_layer_of_an_isothermal_column_is_carried_down_a_dry_adiabat():
    kappa = 287.04 / 1004.64
    d = column(lambda p: 280.0 * (1.0e5 / p) ** kappa, [99000.0, 97000.0, 95000.0, 90000.0, 80000.0])
    got = ba.gcm_layer_air_c(d, 0.9833)
    assert got[0] == pytest.approx(280.0 * 0.9833 ** -kappa - 273.15, abs=1e-6)


def test_the_gcm_layer_of_a_well_mixed_column_gives_the_air_at_the_ground():
    levels = 1.0e5 - np.arange(0.5, 40.0) * 250.0                         # 250-Pa layers through the bottom 10 kPa
    d = column(lambda p: np.full_like(p, 300.0), levels)                   # constant potential temperature
    assert ba.gcm_layer_air_c(d, 0.9833)[0] == pytest.approx(300.0 - 273.15, abs=0.01)


def test_the_box_adds_the_gcm_measure_and_rain_cover():
    assert ba.KEYS[-2:] == ('air_gcm_layer_c', 'raining') and 'air_c' in ba.KEYS


def test_the_microphysics_budget_takes_the_running_totals_over_the_span(tmp_path):
    import json
    names = ('mtime', 'dt', 'tcond', 'tevac', 'tevar', 'train')
    (tmp_path / 'cm1out_stats.ctl').write_text('dset ^cm1out_stats.dat\nvars 6\n' + ''.join(f'{n} 1 99 {n}\n' for n in names)
                                               + 'endvars\n')
    days = np.array([0.0, 1.0, 2.0])
    area = 2 * 2 * 1000.0 ** 2                                         # 2 x 2 columns 1 km wide: 1 mm is 4e6 kg
    rows = [[d * 86400.0, 60.0, 4e6 * 30 * d, 4e6 * 5 * d, 4e6 * 20 * d, 4e6 * 5 * d] for d in days]
    np.array(rows, dtype='<f4').tofile(tmp_path / 'cm1out_stats.dat')
    (tmp_path / 'case.json').write_text(json.dumps(dict(grid=dict(nx=2, ny=2, dx_m=1000.0))))
    out = ba.microphysics_budget(tmp_path, 1.0)
    assert out['span_days'] == [1.0, 2.0]
    assert (out['condensed_mm'], out['evaporated_from_cloud_mm']) == pytest.approx((30.0, 5.0))
    assert (out['evaporated_on_the_way_down_mm'], out['reached_the_ground_mm']) == pytest.approx((20.0, 5.0))
    assert out['share_of_falling_water_evaporated'] == pytest.approx(0.8)
