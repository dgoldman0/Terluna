"""The climate carried onto a changed land mask: the land-water contrast fitted on switched cells, the land's
evaporation along its curve, and the estimate's closure of the surface's water."""
import math

import numpy as np
import pytest

from climate.gcm import surface_change as sc

LAT = np.array([60.0, 20.0, -20.0, -60.0])
LON = np.arange(6) * 60.0
LAND = np.array([[1, 1, 0, 0, 0, 0],
                 [1, 1, 1, 0, 0, 0],
                 [0, 1, 1, 1, 0, 0],
                 [0, 0, 0, 1, 0, 0]], float)


def window(mask=LAND, rain_land=4.0, rain_water=1.5, evaporation_land=2.0, evaporation_water=3.0):
    land = mask > 0.5
    rows = np.arange(4)[:, None] * np.ones(6)
    return dict(lat=LAT, lon=LON, lsm=mask.copy(), ts_k=np.full((4, 6), 300.0),
                pr_mm_day=np.where(land, rain_land + rows, rain_water),            # the land's rain grows southward
                evap_mm_day=np.where(land, evaporation_land, evaporation_water + 0.1 * rows))


def test_the_line_is_fitted_over_cells_switched_either_way():
    a = window()
    other = LAND.copy()
    other[0, 0] = other[1, 2] = other[2, 3] = 0.0                    # three cells land in a, water in b
    other[3, 4] = other[0, 3] = 1.0                                  # two cells water in a, land in b
    b = window(other)
    land_a, land_b = a['lsm'] > 0.5, b['lsm'] > 0.5                  # as water each cell rains 1.2 + 0.1 x its rain as land
    b['pr_mm_day'] = np.where(land_a & ~land_b, 1.2 + 0.1 * a['pr_mm_day'], b['pr_mm_day'])
    a['pr_mm_day'] = np.where(land_b & ~land_a, 1.2 + 0.1 * b['pr_mm_day'], a['pr_mm_day'])
    cells = sc.switched(a, b)
    assert cells['rain_land'].size == 5
    fit = sc.water_rain_fit(cells)
    assert fit['a_mm_day'] == pytest.approx(1.2) and fit['b'] == pytest.approx(0.1) and fit['spread_mm_day'] < 1e-9


def test_the_land_evaporates_along_its_own_curve():
    n = 3.0
    rain, evaporation = np.array([2.0, 1.0, 6.0]), np.array([1.5, 1.0, 3.0])
    em = sc.ceilings(rain, evaporation, n)
    assert sc.curve(rain[0], em[0], n) == pytest.approx(1.5) and math.isinf(em[1])     # the cell lies on its curve
    more = sc.land_evaporation(rain * 1.2, rain, evaporation, n)
    assert 1.5 < more[0] < 1.5 + 0.4                                 # part of the added rain evaporates
    assert more[1] == pytest.approx(1.2)                             # land that evaporated all its rain still does
    assert more[2] - 3.0 < 0.05 * 1.2                                # wet land sheds nearly all of it
    assert sc.land_evaporation(rain * 100.0, rain, evaporation, n)[0] == pytest.approx(em[0], rel=1e-3)   # the ceiling


def test_the_exponent_is_found_on_land_that_lies_on_one_curve():
    f = window()
    land = f['lsm'] > 0.5
    f['pr_mm_day'] = np.where(land, np.linspace(0.5, 9.0, 24).reshape(4, 6), 1.5)
    f['evap_mm_day'] = np.where(land, sc.curve(f['pr_mm_day'], 3.5, 2.5), 3.0)
    fit = sc.fit_exponent(f, within=90.0)
    assert fit['n'] == pytest.approx(2.5, abs=0.051) and fit['ceiling_mm_day'] == pytest.approx(3.5, abs=0.051)


WET = np.where(LAND > 0.5, 0.1, 0.9)                                # the surface the run stood on: a tenth of each land
#                                                                     cell water, a tenth of each water cell land


def test_more_water_in_a_cell_moves_its_rain_by_the_contrast_and_keeps_the_surface_s_water():
    base = window()
    wet = WET.copy()
    wet[1, 2] = 0.6                                                  # a land cell raining 5 mm/day: half of it turns to water
    wet[2, 4] = 0.4                                                  # a water cell beside land: half of it turns to land
    fit = dict(a_mm_day=1.2, b=0.1)
    out, s = sc.carry(base, WET, wet, fit, 3.0)
    w = sc.weights(LAT, LON)
    factor = s['land_rain_factor']
    # the land cell: its rain moves half way from 5 to the line's 1.7, then its land share takes the factor
    assert out['pr_mm_day'][1, 2] == pytest.approx(5.0 + 0.5 * (1.7 - 5.0) + 0.4 * 5.0 * (factor - 1.0))
    # the water cell: half way from its 1.5 to the 6.5 of the land beside it, and its land share takes the factor
    assert out['pr_mm_day'][2, 4] == pytest.approx(1.5 + 0.5 * (6.5 - 1.5) + 0.6 * 6.5 * (factor - 1.0))
    assert out['evap_mm_day'][2, 4] == base['evap_mm_day'][2, 4]     # a cell the run held as water keeps its evaporation
    held_water = LAND < 0.5
    assert np.array_equal(out['evap_mm_day'][held_water], base['evap_mm_day'][held_water])
    assert np.array_equal(out['lsm'], base['lsm'])                   # the mask is the run's
    held = ((base['pr_mm_day'] - base['evap_mm_day']) * w).sum()
    assert ((out['pr_mm_day'] - out['evap_mm_day']) * w).sum() == pytest.approx(held, abs=1e-9)
    assert abs(s['rain_less_evaporation_left_mm_day']) < 1e-9
    assert s['cells'] == dict(more_water_by_a_tenth=1, less_water_by_a_tenth=1)


def test_the_run_s_own_surface_leaves_the_window_as_it_was():
    base = window()
    out, s = sc.carry(base, WET, WET, dict(a_mm_day=1.2, b=0.1), 3.0)
    assert s['land_rain_factor'] == pytest.approx(1.0, abs=1e-9)
    assert out['pr_mm_day'] == pytest.approx(base['pr_mm_day']) and out['evap_mm_day'] == pytest.approx(base['evap_mm_day'])


def test_wet_land_turned_to_water_raises_the_land_s_rain_and_small_changes_move_it_in_proportion():
    base = window(rain_land=6.0)
    fit = dict(a_mm_day=1.2, b=0.1)
    factors = []
    for more in (0.2, 0.4):
        wet = WET.copy()
        wet[2, 2] += more
        _, s = sc.carry(base, WET, wet, fit, 3.0)
        factors.append(s['land_rain_factor'])
        assert s['more_water']['rain_mm_day'][1] < s['more_water']['rain_mm_day'][0]
    assert 1.0 < factors[0] < factors[1]
    assert (factors[1] - 1.0) / (factors[0] - 1.0) == pytest.approx(2.0, rel=0.15)   # twice the change, about twice the rain moved


def test_windows_are_named_with_or_without_a_label():
    assert sc.parse_window('RUN_a:45-59') == ('RUN_a', 45, 59, None)
    assert sc.parse_window('RUN_a:45-59:scaled300.31K') == ('RUN_a', 45, 59, 'scaled300.31K')
    assert sc.product_path('climatology', 'R', 1, 2, 'x').name == 'climatology_R_y1-2_x.npz'
