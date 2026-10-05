"""Preserve the energy search's rejected trajectories and explicit cost scope."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/name).read_text())


def test_energy_products_bind_executed_sources_and_parent_products():
    for name in ['energy_screen.json','energy.json','energy_correction.json','energy_replay.json','energy_budget.json','energy_collection.json']:
        p=read(name)
        for path,h in p['producer']['source_hashes'].items():assert digest(ROOT/path)==h,path
        for path,h in p['producer']['inputs'].items():assert digest(HERE/'results'/path)==h,path
        assert not constants_changed(p['producer']['constants'])
        assert not p['accepted_fleet']


def test_lower_nominal_costs_never_override_coupled_clearance_failures():
    screen=read('energy_screen.json');failed=read('energy.json')
    correction=read('energy_correction.json');last=read('energy_replay.json')
    assert len(screen['cases'])==45
    assert screen['cases'][screen['ranked_proposals'][0]]['nominal_sum_m_s']<1.
    assert not failed['accepted_departure'] and len(failed['runs'])==2
    assert not correction['proposal_pass']
    assert not correction['iterations'][-1]['optimization']['signed_half_space_LP_feasible']
    assert last['candidate_origin']=='next_screened_family_after_infeasible_measured_defect_correction'
    assert not last['accepted_departure'] and len(last['runs'])==1
    for r in failed['runs']+last['runs']:
        assert r['clearance_event'] and not r['completed']
        assert r['geometry']['minimum_sampled_clearance_m']==pytest.approx(100.,abs=1e-6)
        assert not r['local_gates_pass']


def test_torque_account_preserves_rejected_prefix_and_does_not_claim_savings():
    p=read('energy_budget.json')
    assert p['completed'] and not p['accepted_new_departure']
    assert p['same_model_fractional_energy_reduction'] is None
    assert p['conditional_two_day_comparisons']==[]
    assert not p['datasets']['new_departure']['completed_intended_interval']
    assert 'new' not in p['same_model_totals']
    for data in p['datasets'].values():
        assert data['agreements']['maximum_collection_relative_error']<1e-10
        assert data['agreements']['maximum_bound_excess_N_m_kg']<=1e-10
        for r in data['evaluations']:
            c=r['cost']
            assert c['translation_m_s']+c['attitude_equivalent_m_s']==pytest.approx(c['total_equivalent_m_s'])
            assert r['optical']['redirected_band_optical_W']['energy_J']>0
    for field in ['measured_recurring_impulse_m_s_day','fleet_propulsion_W','collection_capacity_W','delivered_power_W']:
        assert p[field] is None


def test_each_rejected_trial_retains_its_own_optical_inventory():
    p=read('energy_collection.json')
    runs=read('energy.json')['runs']+read('energy_replay.json')['runs']
    assert p['completed'] and not p['accepted_fleet']
    assert {r['case'] for r in p['trials']}=={24,20,41}
    assert p['collection_capacity_W'] is None and p['delivered_power_W'] is None
    for trial,run in zip(p['trials'],runs):
        assert not trial['accepted_departure']
        assert trial['duration_s']==pytest.approx(run['end_s']-trial['start_s'])
        light=trial['optical']
        first=light['first_intercept_bolometric_equivalent_W']['energy_J']
        assert 0<light['redirected_band_optical_W']['energy_J']<first
        assert first==pytest.approx(light['front_first_intercept_W']['energy_J']+light['back_first_intercept_W']['energy_J'])
        assert light['unshadowed_aperture_W']['energy_J']==pytest.approx(first+light['mutual_shadow_equivalent_loss_W']['energy_J'])
