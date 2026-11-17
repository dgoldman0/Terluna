"""Checks of the neck's hydraulics and the near-side record: the level search agrees with the network study's tables,
critical flow over a rectangular throat is the textbook B g^0.5 (2E/3)^1.5, friction alone sets a small head's flow,
the two-part balance conserves water, and the record holds together once computed."""
import json
import math

import numpy as np
import pytest

from research.studies.atlas_round import near_side_neck as N
from research.studies.drainage_network import network as NW


def test_the_level_search_agrees_with_the_network_tables():
    rng = np.random.default_rng(1)
    n = 4000
    joined = np.sort(rng.uniform(0.0, 100.0, n))
    table = NW.Table(joined, np.ones(n), rng.normal(1.0, 3.0, n), joined.copy(), 1000.0)
    for supply in (-5.0, 0.0, 10.0, 300.0, 2000.0, 1e6):
        for below in (20.0, 60.0, 101.0):
            assert N.level_for(table, supply, below) == table.level_for(supply, below)


def rectangle(width=1000.0, crest=0.0, sill=1.0):
    """A throat with vertical sides: its section is width times the water over the crest."""
    levels = np.linspace(crest, crest + 50.0, 501)
    return N.Neck(crest, levels, width * (levels - crest), np.full(levels.size, width), sill)


def test_critical_flow_over_a_rectangular_throat():
    neck = rectangle()
    for e in (0.5, 2.0, 8.0):
        textbook = 1000.0 * math.sqrt(N.G) * (2.0 * e / 3.0) ** 1.5
        assert neck.critical(e) == pytest.approx(textbook, rel=0.01)
        assert neck.flow(e, -100.0) == pytest.approx(textbook, rel=0.01)       # the far side fallen away


def test_friction_sets_a_small_heads_flow():
    neck = rectangle(sill=10000.0)
    up, down = 10.0, 9.99
    a, r = 1000.0 * 9.995, 9.995
    expected = a * math.sqrt((up - down) / (N.DRAG * 10000.0 / (N.G * r) + N.EXIT_HEADS / (2 * N.G)))
    assert neck.flow(up, down) == pytest.approx(expected, rel=1e-3)
    assert neck.flow(down, up) == pytest.approx(-expected, rel=1e-3)
    assert neck.flow(-1.0, -2.0) == 0.0                                           # both under the crest


class Flat:
    """A part whose water covers `area` (m2) at every level and loses `rate` (m/s) over it, with the network.Table
    interface the balance reads."""
    def __init__(self, area, rate, bottom=-500.0, top=500.0, n=20001):
        self.c = np.linspace(bottom, top, n)
        step_area = area / n
        self.area = np.cumsum(np.full(n, step_area))
        self.loss = np.cumsum(np.full(n, step_area * rate))

    def at(self, level):
        i = int(np.searchsorted(self.c, level))
        return (0.0, 0.0, 0.0) if i == 0 else (float(self.area[i - 1]), float(self.loss[i - 1]), 0.0)


def test_the_two_part_balance_conserves_water():
    tab = dict(main=Flat(4e12, 2e-8, n=200001), south=Flat(1e12, 2e-8, n=200001))
    neck = rectangle(width=200.0, crest=200.0, sill=6000.0)       # about 30 m of water over it, a head of metres
    supply = dict(main=70000.0, south=3000.0)
    s = N.balance(tab, supply, neck)
    loss_m = tab['main'].at(s['main_level'])[1]
    loss_s = tab['south'].at(s['south_level'])[1]
    assert loss_m == pytest.approx(supply['main'] - s['flow_m3s'], rel=2e-3)
    assert loss_s == pytest.approx(supply['south'] + s['flow_m3s'], rel=2e-3)
    assert s['main_level'] > s['south_level'] > neck.crest
    assert neck.flow(s['main_level'], s['south_level']) == pytest.approx(s['flow_m3s'], rel=0.02)
    assert not s['choked']


@pytest.fixture(scope='module')
def record():
    if not N.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.near_side_neck first')
    return json.loads(N.OUT.read_text())


def test_the_record_holds_together(record):
    crest = record['crest']['four_connected_m']
    levels = [t['level_m'] for t in record['throat']['by_level']]
    sections = [t['section_m2'] for t in record['throat']['by_level']]
    assert levels == sorted(levels) and sections == sorted(sections) and levels[0] > crest
    for label, b in record['budgets'].items():
        if b['near_side_would_settle']['reaches_its_rim']:
            continue
        steady = b['natural_neck']['steady']
        assert steady['south_level_m'] <= steady['main_level_m'] + 1e-6
        assert steady['main_level_m'] >= crest - 1e-6
        rows = b['natural_neck']['with_extra_loss']
        south = [r['south_level_m'] for r in rows if 'south_level_m' in r]
        assert all(a >= b_ - 0.01 for a, b_ in zip(south, south[1:]))           # more loss, a lower south part
        for s in b['straits']:
            assert s['bed_width_m'] > 0 and s['volume_km3'] > 0 and s['bed_m'] < s['shared_level_m']
        for t in b['tide'].values():
            assert t['south_monthly_range_m'] <= t['south_open_neck_monthly_range_m'] + 0.05
