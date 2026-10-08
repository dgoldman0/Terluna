"""Keep completed passages, a failed return and conditional costs distinct."""
import json
from pathlib import Path

import numpy as np
import pytest

from shared.provenance import constants_changed
from .closure_budget import cost_account
from .cycling_search import SCENARIO
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/name).read_text())


def test_mass_weighted_energy_counts_each_member_and_partial_burn():
    times=np.arange(0.,3601.,60.);zero=np.zeros((len(times),2,3))
    command=lambda t,order:np.zeros((len(np.atleast_1d(t)),3))
    controls=np.array([[[1.,0,0]],[[2.,0,0]]]);ratio=np.array([1.2,1.5])
    result=cost_account(times,zero,zero,command,controls,np.array([[0.,7200.]]),ratio)
    assert result['mean_applied_translation_m_s']==pytest.approx(.75)
    assert result['optical_mass_equivalent_m_s']==pytest.approx(1.05)
    assert result['mass_weighted_equivalent_m_s']==pytest.approx(2.1/2.7)
    p=SCENARIO['propulsion'];mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    expected=2.1*mass*p['exhaust_velocity_m_s']/(2*p['efficiency']*np.cos(np.deg2rad(p['cant_deg'])))
    assert result['propulsion_energy_J']==pytest.approx(expected)
    assert sum(result['per_member_energy_J'])==pytest.approx(expected)
    assert result['ideal_exhaust_propellant_kg']==pytest.approx(2*p['efficiency']*expected/p['exhaust_velocity_m_s']**2)


def test_closure_products_bind_sources_constants_and_parent_results():
    for name in ['closure_screen','closure_hardware','closure_sizing','closure_return','closure_validation','closure_budget','closure_audit']:
        p=read(name+'.json')
        assert p['completed'] and not p['accepted_fleet']
        for path,h in p['producer']['source_hashes'].items():assert digest(ROOT/path)==h,path
        for path,h in p['producer']['inputs'].items():assert digest(HERE/'results'/path)==h,path
        assert not constants_changed(p['producer']['constants'])


def test_failed_return_branches_remain_local_evidence():
    p=read('closure_screen.json')
    assert len(p['cases'])==12 and len(p['refinements'])==3
    assert not p['accepted_return']
    for c in p['cases']:
        fit=c['shooting']
        if fit['success']:
            assert fit['restricted_linear_euclidean_lower_bound_m_s']<=fit['mean_delta_v_m_s']
            assert fit['endpoint_position_residual_m']<1e-5
            assert c['encounters']['sampled_pair_violations']>0
    for r in p['refinements']:
        assert not r['proposal_geometry_pass']
        assert r['steps'][0]['crossings']['physical_intersections']>0


def test_tighter_replay_starts_at_epoch_and_preserves_hardware_mass():
    p=read('closure_validation.json');sizing=read('closure_sizing.json')
    assert p['accepted_hardware_passage'] and not p['accepted_return']
    assert len(p['runs'])==3
    for r in p['runs']:
        assert r['comparison']['position_difference_plus_reconstruction_m']<50.
        assert r['beam_margin_deg']>0 and r['rate_deg_s']<.1
    assert p['runs'][0]['start_s']==0 and p['runs'][0]['end_s']==21600.
    assert p['runs'][1]['end_s']==p['runs'][2]['start_s']==43200.
    assert p['runs'][0]['coverage']['passed']
    assert p['runs'][2]['clearance_event'] and not p['runs'][2]['completed']
    assert p['runs'][2]['minimum_sampled_surface_distance_m']==pytest.approx(100.,abs=1e-5)
    assert sizing['total_allocated_mass_kg']>sizing['optical_mass_kg']*1.2
    audit=read('closure_audit.json')
    assert audit['passed'] and audit['total_probes']>audit['grid_probes']
    assert not audit['continuous_certificate']


def test_executed_collection_and_energy_do_not_create_a_recurring_cycle():
    p=read('closure_budget.json')
    assert not p['accepted_return'] and not p['recurring_cost_measured']
    for field in ['fleet_propulsion_W','collection_capacity_W','delivered_power_W']:assert p[field] is None
    for label,entry in p['datasets'].items():
        for e in entry['evaluations']:
            assert len(e['cost']['per_member_energy_J'])==361
            assert e['cost']['propulsion_energy_J']==pytest.approx(sum(e['cost']['per_member_energy_J']))
            optical=e['optical'];first=optical['first_intercept_bolometric_equivalent_W']['energy_J']
            assert first==pytest.approx(optical['front_first_intercept_W']['energy_J']+optical['back_first_intercept_W']['energy_J'])
            assert 0.<optical['redirected_band_optical_W']['energy_J']<first
    sizing=read('closure_sizing.json');prop=SCENARIO['propulsion']
    ratio=np.array(sizing['per_member_total_to_optical_mass_ratio'])
    capacity=(ratio-1.2)*SCENARIO['baseline_areal_mass_kg_m2']*10000.**2*prop['specific_power_W_kg']
    peaks=np.maximum.reduce([e['cost']['per_member_peak_W']
        for label in ['hardware_service','hardware_departure','rejected_return']
        for e in p['datasets'][label]['evaluations']])
    assert p['hardware_account']['all_members_fit_installed_capacity']==bool(np.all(peaks<=capacity))
    assert p['hardware_account']['full_margin_target_met']==bool(np.all(prop['peak_margin_factor']*peaks<=capacity))
