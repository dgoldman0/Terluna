"""Check source-trace timing and conditional inference, without claiming validation."""
from copy import deepcopy
import json

import pytest

from atmosphere.snow import mcsnow_growth_trace as trace


@pytest.fixture
def product():
    return json.loads(trace.PRODUCT.read_text())


def test_matching_mass_interpolates_time_and_preserves_bounds():
    rows = [dict(elapsed_s=0.,mass_kg=2.,equatorial_dimension_m=3.),
            dict(elapsed_s=4.,mass_kg=6.,equatorial_dimension_m=9.)]
    result = trace.matched_mass(rows,3.)
    assert result['state'] == dict(elapsed_s=1.,mass_kg=3.,equatorial_dimension_m=4.5)
    assert result['bracket_elapsed_s'] == [0.,4.]
    assert trace.matched_mass(rows,1.) is None
    assert trace.matched_mass(rows,7.) is None
    with pytest.raises(ValueError,match='increase strictly'):
        trace.matched_mass(rows[::-1],3.)


def test_branch_mass_recovers_known_two_density_history():
    # Grow 2 kg of compact ice, then deposit 3 kg at one-quarter ice density.
    volume = 2/1000+3/250
    assert trace.branch_mass(5.,volume,1000.,.25) == pytest.approx(2.)
    with pytest.raises(ValueError,match='branching plate'):
        trace.branch_mass(5.,volume,1000.,1.)


def test_unknown_source_rejected_before_instrumentation():
    with pytest.raises(ValueError,match='Unrecognized'):
        trace.source_patch('unidentified Fortran','native')
    with pytest.raises(ValueError,match='Unknown'):
        trace.source_patch('unidentified Fortran','fitted')


def actual_single_step(product):
    selected = lambda r:r['variant']=='IGF2plus' and r['mode']=='native' and r['dt_s']==.015625 and r['age_s']==600
    before = deepcopy(next(r['trace'] for r in product['source_traces'] if selected(r)))
    after = deepcopy(next(r['state'] for r in product['comparisons'] if selected(r)))
    # Shift one actually evaluated step to elapsed zero; retain all physical fields.
    before.update(evaluation_elapsed_s=0.,state_elapsed_s=.015625)
    after['elapsed_s'] = .015625
    return [before],[dict(elapsed_s=0.),after],dict(dt_s=.015625,duration_s=.015625)


def test_actual_source_trace_predicts_saved_mass_volume_and_habit(product):
    trace.validate_trace(*actual_single_step(product))


@pytest.mark.parametrize('damage',['time','mass','volume','habit','ventilation','missing'])
def test_stale_or_incomplete_source_trace_rejected(product,damage):
    rows,states,config = actual_single_step(product)
    if damage=='time':
        rows[0]['evaluation_elapsed_s'] = .015625
    elif damage=='mass':
        rows[0]['mass_kg'] *= 1.1
    elif damage=='volume':
        rows[0]['equatorial_dimension_m'] *= 1.1
    elif damage=='habit':
        rows[0]['effective_growth_factor'] *= 2
    elif damage=='ventilation':
        rows[0]['ventilation_vapor'] *= 1.1
    else:
        rows.clear()
    with pytest.raises(ValueError):
        trace.validate_trace(rows,states,config)


def test_trace_preserves_preceding_trajectories_and_output_cadence(product):
    assert len(product['parity']) == 16
    assert all(r['exact_common_snapshot_and_air_equality'] for r in product['parity'])
    assert all(r['common_snapshots_exactly_equal'] for r in product['interpolation'])
    for row in product['branch_mass']:
        assert row['native_final_inferred_branch_mass_kg'] == pytest.approx(
            row['first_branch_evaluation']['mass_kg'],rel=1e-10)


def test_rate_intervention_is_distinct_from_directional_intervention(product):
    rows = {r['mode']:r['trace'] for r in product['fixed_state']
            if r['variant']=='IGF2plus' and r['author_age_s']==600}
    native,base,radius = (rows[m] for m in ('native','plate_base','sphere_radius'))
    for key in ('mass_kg','equatorial_dimension_m','polar_dimension_m','terminal_speed_m_s',
                'capacitance_m','deposition_density_kg_m3'):
        assert native[key] == base[key] == radius[key]
    assert native['effective_growth_factor'] == base['effective_growth_factor']
    assert base['deposition_rate_kg_s'] < native['deposition_rate_kg_s']
    assert radius['deposition_rate_kg_s'] == native['deposition_rate_kg_s']
    assert radius['effective_growth_factor'] < native['effective_growth_factor']


def test_better_endpoint_does_not_get_promoted_to_reproduction(product):
    primary = [r for r in product['comparisons'] if r['variant']=='IGF2plus'
               and r['age_s']==600 and r['dt_s']==.015625]
    assert all(not r['comparison']['all_within_printed_rounding'] for r in primary)
    assert all(r['matched'] is None for r in product['matched_mass']
               if r['mode']=='historical' and r['author_age_s']==600)
