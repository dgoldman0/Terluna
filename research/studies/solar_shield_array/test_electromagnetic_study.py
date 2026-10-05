"""Scientific boundaries and conserved accounts for the bounded EM study."""
import json
from pathlib import Path
import numpy as np
import pytest
from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/name).read_text())


def test_em_products_bind_inputs_and_do_not_accept_a_trajectory():
    for name in ['electromagnetic.json','electromagnetic_audit.json']:
        p=read(name)
        assert p['completed'] and not p['accepted_trajectory'] and not p['accepted_cycle'] and not p['accepted_fleet']
        for f,h in p['producer']['source_hashes'].items():assert digest(ROOT/f)==h,f
        for f,h in p['producer']['inputs'].items():assert digest(HERE/'results'/f)==h,f
        assert not constants_changed(p['producer']['constants'])
    optical=read('electromagnetic_audit.json')['producer']['optical_product']
    assert digest(ROOT/optical['path'])==optical['sha256']


def test_inherited_failed_return_and_energy_are_preserved():
    p=read('electromagnetic.json');r=p['stages'][-1]
    assert p['original_target_W']==23.835e12 and not p['new_target_adopted']
    assert p['reference']['overloaded_members']==160
    assert r['net_energy_deficit_J']==pytest.approx(76.0063e12,rel=1e-5)
    assert r['peak_pooled_deficit_W']>24e9
    assert all(s['propulsion_energy_J']==pytest.approx(s['translation_energy_J']+s['attitude_energy_J']) for s in p['stages'])
    assert r['sum_local_positive_deficit_J']>=r['pooled_positive_deficit_J']>=r['net_energy_deficit_J']


def test_actual_square_interactions_converge_and_conserve_momentum():
    p=read('electromagnetic.json')
    assert max(r['relative_refinement'] for r in p['pair_samples'])<1e-4
    assert max(r['action_reaction_relative'] for r in p['pair_samples'])<1e-10
    assert max(r['angular_momentum_relative'] for r in p['pair_samples'])<1e-9
    assert max(r['maximum_pair_force_relative_refinement'] for r in p['group_fields'])<1e-4
    for r in read('electromagnetic_audit.json')['current_allocations']:
        assert np.linalg.norm(r['net_internal_force_N'])<1e-6
        assert np.linalg.norm(r['internal_total_torque_about_origin_N_m'])<1e-3
        assert r['inductance_min_eigenvalue_H']>0
        assert len(r['all_starts'])==9


def test_storage_chronology_and_beam_crossing_assumptions_are_explicit():
    p=read('electromagnetic.json');a=read('electromagnetic_audit.json')
    conservative=p['storage_cases'][1]['battery_mass_kg']
    assert 140e6<a['chronological_storage']['battery_mass_at_200Wh_kg']<=conservative
    assert any(r['clear_links']==0 for r in a['beam_geometry'])
    for r in a['beam_geometry']:
        for link in r['links']:
            assert link['transmission_if_each_crossing_99percent']==pytest.approx(.99**len(link['blockers']))
    assert a['film_crossings']['laser_normal_incidence_T1064']>.95
    assert a['film_crossings']['microwave_lossless_slab_T']>.999


def test_pair_scaling_and_power_source_mass_are_not_free():
    p=read('electromagnetic.json')
    assert all(x['source_bus_W']>x['delivered_W'] for x in p['link_cases'])
    assert all(x['generation']['installed_mass_kg']>0 for x in p['hub_power'])
    rows=[r for r in p['coil_scaling'] if r['required_pair_force_N']==100]
    assert all(a['ampere_turns']<b['ampere_turns'] for a,b in zip(rows,rows[1:]))
    assert rows[-1]['installed_mass_kg']>5e6
