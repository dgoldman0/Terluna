"""Checks of the drainage on the climate run on the final surface: the climatology carried onto today's ground is
turned_climatology's carrying with the final surface's rotation, the loader reads it as the network study's loader
reads a climatology, and the record holds together once computed."""
import json

import numpy as np
import pytest

from climate.gcm import turned_climatology as TC
from research.studies.atlas_round import final_drainage as FD


@pytest.fixture(scope='module')
def carried():
    path = FD.ROOT / FD.CLIMATOLOGY
    if not path.exists() or not (FD.ROOT / FD.TURNED_CLIMATOLOGY).exists():
        pytest.skip('the run\'s climatology is not on this machine')
    return path


def test_the_climatology_is_carried_with_the_final_surfaces_rotation(carried):
    z = np.load(carried)
    meta = json.loads(str(z['metadata']))
    src = np.load(FD.ROOT / FD.TURNED_CLIMATOLOGY)
    clim = {k: src[k] for k in src.files if k != 'metadata'}
    _, frame = FD.frame_record()
    expected = TC.to_today(clim, np.asarray(frame['frame']['rotation_today_to_turned'], float))
    for k in ('pr_mm_day', 'evap_mm_day', 'lsm'):
        np.testing.assert_allclose(z[k], expected[k], rtol=0, atol=1e-6)
    assert meta['frame'] == 'today' and meta['run'] == FD.RUN and list(meta['years']) == FD.YEARS
    assert meta['carried']['climatology']['sha256'] == FD.sha256(FD.ROOT / FD.TURNED_CLIMATOLOGY)


def test_the_loader_keeps_the_network_studys_balance(carried):
    runoff, lake_net, precip, evap = FD.load(1.0)
    _, lake_net_closed, _, _ = FD.load(1.2)
    assert runoff.shape == (2880, 5760)
    assert float(runoff.min()) >= 0.0
    np.testing.assert_allclose(runoff, np.maximum(precip - evap, 0.0))
    assert float((lake_net_closed - lake_net).max()) <= 1e-6                    # more evaporation, less balance


@pytest.fixture(scope='module')
def record():
    if not FD.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.final_drainage first')
    return json.loads(FD.OUT.read_text())


def test_the_record_holds_together(record):
    b = record['budget']
    assert b['surplus_m3s'] == b['delivered_to_the_seas_m3s'] - b['seas_loss_m3s']
    assert b['open_water_evaporation_factor_closing_it'] == pytest.approx(
        1.0 + b['surplus_m3s'] / b['open_water_evaporation_m3s'], abs=6e-4)
    for label in FD.BUDGETS:
        case = record['cases'][label]
        seas = case['seas']
        for held in ('south_pole_aitken', 'humboldtianum'):
            assert abs(seas[held]['balance_m3s']) <= 2                             # the weirs meet each held sea's need
        assert sum(case['released_m3s'].values()) == pytest.approx(case['moscoviense_surplus_m3s'], abs=3)
        settle = seas['near_side']['would_settle']
        assert settle['reaches_its_rim'] == (seas['near_side']['balance_m3s'] > 0)
    closed = record['cases']['closed_budget']['evaporation']
    assert closed == record['budget']['open_water_evaporation_factor_closing_it']
