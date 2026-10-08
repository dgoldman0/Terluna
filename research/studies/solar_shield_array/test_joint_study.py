"""Preserve the distinction between executed geometry and hardware feasibility."""
import json
import hashlib
from pathlib import Path
import numpy as np
from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/(name+'.json')).read_text())


def test_joint_sources_parents_and_constants_remain_bound():
    for path in (HERE/'results').glob('joint_*.json'):
        data=json.loads(path.read_text())
        assert data.get('completed'),path
        producer=data.get('producer')
        if producer is None:continue
        for name,h in producer['source_hashes'].items():assert digest(ROOT/name)==h,(path,name)
        for name,h in producer['inputs'].items():assert digest(HERE/'results'/name)==h,(path,name)
        assert not constants_changed(producer['constants'])


def test_geometric_return_is_separate_from_power_feasible_cycle():
    fine=read('joint_corridor_108_fine');account=read('joint_account');audit=read('joint_audit')
    assert len(fine['runs'])==3 and all(r['completed'] for r in fine['runs'])
    assert fine['runs'][-1]['coverage']['passed']
    assert audit['geometric_return_plus_service_replayed']
    assert not fine['accepted_cycle'] and not fine['actuator_and_supply_feasible']
    assert not account['accepted_cycle'] and not account['accepted_fleet']
    assert account['actual_generation_W'] is None and account['full_operating_cycle_cost_J'] is None
    for row in fine['comparison']:
        assert row['maximum_position_with_reconstruction_m']<50
        assert row['conditional_clearance_after_replay_allowance_m']>100
    measured=account['measured_sequence']
    assert measured['overloaded_members']>0
    assert measured['minimum_installed_to_peak']<1
    assert not measured['resized_actuator_mass_propagated']
    supply=account['local_supply_upper_bound']
    assert supply['buffer_60s_unserved_J']>0 and supply['first_unserved_s'] is not None
    assert not supply['assumptions']['collector_mass_and_photon_forces_replayed']
    assert audit['handovers_executed']==0


def test_whole_executed_cost_and_light_are_counted_once():
    account=read('joint_account');stages=account['stages'];whole=account['measured_sequence']
    assert len(stages)==3
    assert [s['start_s'] for s in stages][1:]==[s['end_s'] for s in stages][:-1]
    energies=[s['evaluations'][-1]['propulsion_energy_J'] for s in stages]
    np.testing.assert_allclose(sum(energies),whole['propulsion_energy_J'],rtol=1e-14)
    intercepted=sum(s['evaluations'][-1]['first_intercept_J'] for s in stages)
    np.testing.assert_allclose(intercepted,whole['first_intercept_J'],rtol=1e-14)
    np.testing.assert_allclose(whole['ideal_exhaust_kg'],2*.7*whole['propulsion_energy_J']/30000**2)
    assert 0<whole['redirected_band_J']<whole['first_intercept_J']
    # Extra collector light is an uninstalled sensitivity, never added here.
    assert account['local_supply_upper_bound']['available_bus_upper_energy_J']>0


def test_original_regressions_and_failed_proposals_survive():
    regression=read('joint_regression')
    assert regression['runs'][0]['coverage']['passed']
    assert regression['runs'][1]['completed']
    assert regression['runs'][2]['clearance_event']
    assert regression['runs'][2]['witness']['pair']==[325,347]
    assert abs(regression['runs'][2]['end_s']/3600-13.893535)<.001
    assert len(read('joint_screen')['cases'])==6
    assert len(read('joint_slots')['cases'])==3
    assert len(read('joint_bridge')['cases'])==2
    assert not read('joint_grouped')['accepted_cycle']


def test_exported_seed_round_trips_exact_loaded_states():
    seed=read('joint_seed')
    for name,values in seed['arrays'].items():
        array=np.asarray(values,dtype='<f8')
        assert array.shape==((361,) if name=='mass_ratio' else (361,6))
        assert np.isfinite(array).all()
        assert hashlib.sha256(array.tobytes()).hexdigest()==seed['float64_little_endian_hashes'][name]
    assert np.asarray(seed['arrays']['mass_ratio']).min()>1.2


def test_maintenance_preserves_reconstructible_executed_source():
    checks=json.loads((HERE/'joint_checks.json').read_text())
    maintenance=checks['post_run_provenance_maintenance']
    source=(HERE/'joint_screen.py').read_text()
    assert hashlib.sha256(source.encode()).hexdigest()==maintenance['current_source_sha256']
    original=source.replace('from shared import constants as K\n','').replace('4*K.MOON_RADIUS','4*1737400')
    assert hashlib.sha256(original.encode()).hexdigest()==maintenance['executed_source_sha256']
    for row in checks['products']:
        assert digest(ROOT/row['path'])==row['sha256']
