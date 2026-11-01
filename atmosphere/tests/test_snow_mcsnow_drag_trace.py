"""Fixed-state interface and archive-interpretation checks; no empirical validation."""
from copy import deepcopy
import json

import pytest

from atmosphere.snow import mcsnow_drag_trace as trace


@pytest.fixture
def product():
    return json.loads(trace.PRODUCT.read_text())


def test_primary_trace_preserves_archived_state_and_prior_native_result(product):
    state = next(s for s in product['states'] if s['variant'] == 'IGF2plus' and
                 s['elapsed_s'] == 600 and s['temperature_K'] == 258.65)
    row = next(r for r in product['traces'] if r['case_id'] == state['case_id'] and r['control'] == 'archived')
    assert row['input_mass_kg'] == row['drag_mass_kg'] == 4.206150e-9
    assert row['ice_volume_m3'] == 1.3383342e-11
    assert row['aspect_ratio'] == .01470726
    # Independent expectation: the preceding archived-state execution, before this trace caller existed.
    assert row['terminal_speed_m_s'] == pytest.approx(.1461897908353587, rel=2e-14)
    assert row['archive_speed_m_s'] == .2597
    assert row['author_helper_speed_m_s'] == pytest.approx(row['terminal_speed_m_s'], rel=2e-14)


def test_native_geometry_and_three_call_paths_are_consistent(product):
    columns = ('case_id', 'control_id', 'temperature_K', 'pressure_Pa', 'input_mass_kg', 'ice_volume_m3',
               'aspect_ratio', 'speed_from_reynolds_m_s', 'speed_boehm_m_s', 'terminal_speed_m_s',
               'equatorial_dimension_m', 'polar_dimension_m', 'diameter_m', 'ice_density_kg_m3',
               'drag_mass_kg', 'projected_area_m2', 'native_projected_area_m2')
    raw = [{k: r[k] for k in columns} for r in product['traces']]
    trace.validate_trace(raw, product['states'])
    assert product['maximum_helper_native_relative_difference'] < 1.1e-7
    assert not product['build']['physics_patches']


@pytest.mark.parametrize('damage', ['missing', 'duplicate', 'mass', 'axis'])
def test_incomplete_or_stale_native_trace_is_rejected(product, damage):
    rows = [{k: v for k, v in r.items() if isinstance(v, (int, float))} for r in product['traces']]
    if damage == 'missing':
        rows.pop()
    elif damage == 'duplicate':
        rows.append(deepcopy(rows[0]))
    elif damage == 'mass':
        rows[0]['input_mass_kg'] *= 2
    else:
        rows[0]['equatorial_dimension_m'] /= 2
    with pytest.raises(ValueError):
        trace.validate_trace(rows, product['states'])


def test_unknown_author_code_is_rejected_before_execution(tmp_path, monkeypatch):
    (tmp_path/'McSnow_paper_plotscripts.zip').write_bytes(b'not the identified author archive')
    monkeypatch.setattr(trace.model, 'CACHE', tmp_path)
    with pytest.raises(ValueError, match='checksum'):
        trace.author_helper()


def test_speed_disagreement_does_not_reject_matching_plotted_fields():
    units = {k: .5 for k in trace.model.QUANTITIES}
    units['terminal_speed_m_s'] = 1000
    row = dict(group='grid', variant='IGF2plus', mode='native', elapsed_s=600,
               temperature_K=258.65, comparison=dict(rounding_units=units))
    result = trace.figure_fields_summary(dict(comparisons=[row]))[0]
    assert result['matches_plotted_mass_axes_and_aspect_ratio'] == 1
    assert result['matches_all_fields_except_speed'] == 1
    units['ice_volume_m3'] = 3
    result = trace.figure_fields_summary(dict(comparisons=[row]))[0]
    assert result['matches_plotted_mass_axes_and_aspect_ratio'] == 1
    assert result['matches_all_fields_except_speed'] == 0
    units['equatorial_dimension_m'] = 3
    result = trace.figure_fields_summary(dict(comparisons=[row]))[0]
    assert result['matches_plotted_mass_axes_and_aspect_ratio'] == 0


def test_solid_envelope_is_not_an_archive_match_or_universal_repair(product):
    selected = [s['case_id'] for s in product['states'] if s['elapsed_s'] == 600 and s['temperature_K'] == 258.65]
    dense = [r for r in product['traces'] if r['case_id'] in selected and r['control'] == 'solid_envelope']
    assert len(dense) == 3
    assert all(r['drag_mass_kg'] > 2*r['input_mass_kg'] for r in dense)
    assert all(not r['within_archived_speed_rounding'] for r in dense)
    assert any(r['native_archive_speed_ratio'] > 1.8 for r in dense)
