"""Checks of the flat rings' comfort analysis: the terrain twin's height classes, the correction drawn from them, the
nearest land cell and the counts of comfortable hours."""
import numpy as np
import pytest

from climate.crm import ring_comfort as rcf


def test_the_class_table_takes_each_class_mean_by_day_and_by_night_over_land_only():
    heights = np.array([0.0, 300.0, 700.0, 900.0, 2500.0, 2500.0])
    land = np.array([True, True, True, True, True, False])
    day = np.array([[True] * 6, [False] * 6])
    d_air = np.array([-heights / 1000.0, -2.0 * heights / 1000.0])    # 1 C cooler per km by day, 2 by night
    d_dew = np.array([np.full(6, -1.0), np.full(6, -3.0)])
    table = rcf.class_table(heights, land, d_air, d_dew, day, classes=((0, 500), (500, 1000), (2000, 3000)))
    assert [r['columns'] for r in table] == [2, 2, 1]                  # the lake at 2.5 km is left out
    assert [r['mean_height_m'] for r in table] == pytest.approx([150.0, 800.0, 2500.0])
    assert [r['air_day_c'] for r in table] == pytest.approx([-0.15, -0.8, -2.5])
    assert [r['air_night_c'] for r in table] == pytest.approx([-0.3, -1.6, -5.0])
    assert {r['dewpoint_day_c'] for r in table} == {-1.0} and {r['dewpoint_night_c'] for r in table} == {-3.0}
    assert rcf.per_km(table, 'air_night_c') == pytest.approx(-2.0)


def test_the_correction_is_nothing_at_sea_level_interpolates_between_classes_and_holds_above_the_highest():
    table = [dict(mean_height_m=1000.0, air_day_c=-1.0, dewpoint_day_c=-2.0, air_night_c=-0.5, dewpoint_night_c=-1.5),
             dict(mean_height_m=3000.0, air_day_c=-3.0, dewpoint_day_c=-4.0, air_night_c=-2.5, dewpoint_night_c=-3.5)]
    air, dew = rcf.correction(table, np.array([0.0, 500.0, 2000.0, 5000.0]), np.array([True, True, True, False]))
    assert air == pytest.approx([0.0, -0.5, -2.0, -2.5])
    assert dew == pytest.approx([0.0, -1.0, -3.0, -3.5])


def test_the_nearest_land_cell_passes_over_nearer_water():
    land = np.array([[True, False, False, True]])
    rows, cols = rcf.nearest_land_cells(np.array([0.0, 0.0]), np.array([12.0, 22.0]), np.array([0.0]),
                                        np.array([0.0, 10.0, 20.0, 30.0]), land)
    assert list(rows) == [0, 0] and list(cols) == [0, 3]


def test_hours_are_shares_of_the_lunar_day():
    from climate.crm.cm1_run import solar_day_s
    assert rcf.hours(0.5) == pytest.approx(solar_day_s() / 7200.0)
    assert rcf.hours(np.array([0.0, 1.0])) == pytest.approx([0.0, solar_day_s() / 3600.0])


def ring(lat, height, strict, gcm_strict, land=None):
    n = len(height)
    return dict(land=np.ones(n, bool) if land is None else np.asarray(land), lat=np.full(n, lat), height=np.asarray(height, float),
                strict=np.asarray(strict, float), loose=np.asarray(strict, float), strict_flat=np.zeros(n),
                air=np.full(n, 20.0), dew=np.full(n, 12.0), gcm_strict=np.asarray(gcm_strict, float),
                gcm_loose=np.asarray(gcm_strict, float), gcm_height=np.full(n, 1500.0), gcm_air=np.full(n, 21.0),
                gcm_dew=np.full(n, 11.0))


def test_the_rings_land_is_pooled_by_latitude_and_height_and_counted_in_shares():
    a = ring(10.0, [100, 100, 3000, 3000], [0.0, 0.0, 0.2, 0.4], [0.5] * 4)
    b = ring(50.0, [3000, 3000], [0.3, 0.3], [0.1, 0.1], land=[True, False])
    rows = rcf.by_latitude_and_height([a, b])
    assert [(r['from_deg'], r['from_m'], r['columns']) for r in rows] == [(0, 0, 2), (0, 2000, 2), (30, 2000, 1)]
    assert rows[1]['comfortable_hours']['strict'] == pytest.approx(rcf.hours(0.3))
    shares = rcf.land_shares([a, b])                                    # hours 0, 0, 142, 283, 213 against 354 (x4) and 71
    assert shares['columns'] == 5
    assert shares['at_least']['100'] == pytest.approx(dict(corrected=0.6, gcm=0.8))
    assert shares['at_least']['350'] == pytest.approx(dict(corrected=0.0, gcm=0.8))


def test_column_wetness_reads_the_class_each_column_was_given(tmp_path):
    (tmp_path / 'terluna_surface.txt').write_text('3\n-1e9 2.0 2.0 16 300.0 300.0\n2.0 4.0 1.0 20 295.0 295.0\n'
                                                  '4.0 1e9 1.0 30 295.0 295.0\n')
    (tmp_path / 'LANDUSE.TBL').write_text("USGS\n33,2, 'ALBD   SLMO'\nSUMMER\n16,  8., 1.0, .98,'Water Bodies'\n"
                                          "20,  20., .05, .95,'Terluna placeholder land, wetness 0.05'\n"
                                          "30,  20., .90, .95,'Terluna placeholder land, wetness 0.90'\nWINTER\n"
                                          "20,  20., .07, .95,'winter row'\n")
    wet = rcf.column_wetness(tmp_path, np.array([1.0, 3.0, 5.0]))
    assert np.isnan(wet[0]) and list(wet[1:]) == [0.05, 0.90]           # water, then the summer rows
