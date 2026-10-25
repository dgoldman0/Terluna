"""A box's lunar days one after another (box_days.py), its pure pieces on made-up outputs: a lunar day takes the outputs
between its restart times, its rain falls in the quarters of the site's hour angle, and its temperature's change runs
from its first output to its last."""
import numpy as np
import pytest

from climate.crm import box_days as bd


def made_up_rows(days=(0.0, 59.0), step=0.125, lunar_day=29.5):
    """Outputs every 3 hours: the hour angle from sunrise round the lunar day, rain only in the afternoon, the air at
    40 km cooling 0.1 K a day."""
    rows = []
    for day in np.arange(days[0], days[1] + 1e-9, step):
        hour = (-90.0 + 360.0 * day / lunar_day + 180.0) % 360.0 - 180.0
        rows.append(dict(day=float(day), hour=float(hour), rain=2.0 if 0.0 <= hour < 90.0 else 0.0, cape=100.0,
                         cin=3.0, case='box', t_k=[260.0 - 0.1 * day, 230.0]))
    return rows


def test_a_lunar_day_takes_the_outputs_between_its_restart_times():
    rows = made_up_rows()
    second = bd.span_summary(rows, (29.5, 59.0))
    assert second['span_days'] == [29.5, 59.0] and second['outputs'] == 237
    assert second['rain_by_quarter_mm_day'][1] == pytest.approx(2.0) and second['rain_by_quarter_mm_day'][3] == 0.0
    assert second['rain_mm_day'] == pytest.approx(0.5, abs=0.02)
    assert second['temperature_k']['40 km']['change_k_day'] == pytest.approx(-0.1)
    assert second['temperature_k']['60 km']['mean'] == pytest.approx(230.0)
    assert bd.span_summary(rows, (59.0, 88.5)) is None                  # one output is not a day


def test_a_quarter_without_outputs_holds_none():
    rows = [r for r in made_up_rows() if r['hour'] < 0.0]
    assert bd.quarters(rows, 'rain')[1:3] == [None, None]
