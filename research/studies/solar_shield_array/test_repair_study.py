"""Retain failed maneuvers, measured energy and hypothetical fleet assumptions."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/name).read_text())


def test_repair_products_pin_sources_parents_and_unclosed_cycle():
    paths=list((HERE/'results').glob('repair_*.json'))
    assert len(paths)>=7
    for path in paths:
        p=json.loads(path.read_text())
        assert p['completed'] and not p['accepted_cycle'] and not p['accepted_fleet']
        for source,h in p['producer']['source_hashes'].items():assert digest(ROOT/source)==h,source
        for source,h in p['producer']['inputs'].items():assert digest(HERE/'results'/source)==h,source
        assert not constants_changed(p['producer']['constants'])


def test_proposal_failures_and_solver_limit_remain_distinct():
    p=read('repair_screen.json')
    assert len(p['cases'])==11
    assert all(not r['early_proposal_pass'] for r in p['cases'])
    assert p['actual_encounter'][-1]['displacement_body_m'][2]==pytest.approx(100.,abs=1e-5)
    local=read('repair_local_screen.json')
    assert local['iterations'][-1]['fit']['status']==1  # time limit, not proven infeasibility
    refined=read('repair_local_refined.json')
    assert all(r['fit']['success'] for r in refined['iterations'])
    assert not refined['early_proposal_pass']


def test_partial_energy_sums_members_and_keeps_full_cycle_unknown():
    p=read('repair_budget.json')
    assert p['full_cycle_energy_J'] is None and p['full_cycle_average_W'] is None
    for trial in p['trials']:
        assert trial['start_s']==43200.
        for e in trial['evaluations']:
            cost=e['cost'];assert len(cost['per_member_energy_J'])==361
            assert cost['propulsion_energy_J']==pytest.approx(sum(cost['per_member_energy_J']))
            opt=e['optical'];a=opt['first_intercept_bolometric_equivalent_W']
            assert a['energy_J']==pytest.approx(sum(a['per_member_energy_J']))
            assert a['energy_J']==pytest.approx(opt['front_first_intercept_W']['energy_J']+opt['back_first_intercept_W']['energy_J'])
        full=trial['whole_executed_prefix']
        assert full['duration_s']==trial['end_s']
        assert full['first_intercept_mean_W']==pytest.approx(full['first_intercept_energy_J']/full['duration_s'])


def test_static_and_moving_optical_comparison_uses_explicit_duty_and_band():
    p=read('repair_trade.json');r=p['one_AU_face_on_source_comparison']
    assert r[0]['extra_five_percent_dimming_W']==pytest.approx(585.5707346646982e12,rel=.002)
    assert r[2]['combined_filter_and_dimming_band_W']>30e15
    assert p['broad_static_allocation']['optically_available_W']>200e15
    moving=p['capacity_relaxation_sensitivity']
    local=p['measured_twelve_hour_passage']['first_interception_fraction_of_same_clear_film']
    assert .58<local<.60
    assert moving['passage_only_fraction_if_repeated_every_two_days']==pytest.approx(local/4)
    assert p['actual_fleet_collection_W'] is None and not p['new_target_adopted']
    assert p['original_target_W']==23.835e12
