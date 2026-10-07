"""Checks on the protection design point."""
import json

import pytest

from atmosphere.loss_response import model as lr
from research.studies.protection_architecture import run


def test_cycle_time_of_one_kilogram_a_second_is_about_a_hundred_billion_years():
    assert run.cycle_time_years(1.0, run.atmosphere_mass_kg()) == pytest.approx(9.95e10, rel=0.01)


def test_design_point_holds_together():
    wind = lr.solar_wind_losses()
    p = run.point('titania_stack', 'lte', 'quiet', 2e-4, wind, run.atmosphere_mass_kg())
    assert p['status'] == 'thermal_column' and p['covers_heating']
    assert 0.0 < p['glow_share_of_heating'] < 1.0
    # traced, the whole heat counts as more transmission than the swarm passes
    assert p['heating_as_disk_count_transmission'] > p['transmission']
    none, wake, september = (p['total_loss_kg_s'][s] for s in run.SCENARIOS)
    assert none['low'] == pytest.approx(p['uv_driven_loss_kg_s'] + wind['low']['proton_sputtering_kg_s']
                                        + p['exosphere']['loss_kg_s']['no_magnetosphere']['low'])
    assert september['central'] == pytest.approx(p['uv_driven_loss_kg_s']
                                                 + p['exosphere']['loss_kg_s']['september_magnetosphere']['central'])
    assert september['central'] < none['central']
    assert wake['central'] <= none['central']
    for scenario in run.SCENARIOS:
        assert p['cycle_time_years'][scenario]['high'] < p['cycle_time_years'][scenario]['low']
    assert p['heating_radius_R'] < p['protected_radius_R']
    # the total at the ring fleet's radius is the headline's central total
    at = p['total_at_protected_radius_kg_s']
    assert at['4']['september_magnetosphere'] == pytest.approx(september['central'])


def test_results_cover_every_level_and_shield():
    out = json.loads((run.HERE / 'results' / 'design_point.json').read_text())
    assert out['schema'] == run.SCHEMA
    assert set(out['summary']) == {'tight', 'standard', 'relaxed'}
    for level in out['summary'].values():
        assert set(level) == set(lr.SHIELDS)
