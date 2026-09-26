"""Checks on the protection design point."""
import json

import pytest

from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import model as lr
from research.studies.protection_architecture import run


def test_cycle_time_of_one_kilogram_a_second_is_about_a_hundred_billion_years():
    assert run.cycle_time_years(1.0, run.atmosphere_mass_kg()) == pytest.approx(9.95e10, rel=0.01)


def test_design_point_holds_together():
    wind = lr.solar_wind_losses()
    p = run.point('titania_stack', 'lte', 'quiet', 2e-4, ab.band_weights(), wind, run.atmosphere_mass_kg())
    assert p['status'] == 'thermal_column'
    assert 0.0 < p['glow_share_of_heating'] < 1.0
    assert p['heating_as_equivalent_transmission'] > p['transmission']
    assert p['total_loss_kg_s']['low'] == pytest.approx(p['uv_driven_loss_kg_s'] + wind['low']['total_kg_s'])
    assert p['cycle_time_years']['high'] < p['cycle_time_years']['low']
    assert p['exobase_radius_R'] - 0.2 < p['protected_radius_R'] <= p['exobase_radius_R']


def test_results_cover_every_level_and_shield():
    out = json.loads((run.HERE / 'results' / 'design_point.json').read_text())
    assert out['schema'] == run.SCHEMA
    assert set(out['summary']) == {'tight', 'standard', 'relaxed'}
    for level in out['summary'].values():
        assert set(level) == set(lr.SHIELDS)
