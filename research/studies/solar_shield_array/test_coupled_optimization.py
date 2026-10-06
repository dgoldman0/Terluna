"""Scientific evidence gates and the optimizer's near-duplicate-knot regression."""
import json
from pathlib import Path
import numpy as np
from shared.provenance import constants_changed
from .fleet_run import digest
from . import coupled_stable

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):
    return json.loads((HERE/'results'/('coupled_'+name+'.json')).read_text())


def test_nearly_grid_aligned_timing_has_bounded_continuous_attitude_load(monkeypatch):
    monkeypatch.setattr(coupled_stable,'sun_frame',lambda env,t:np.eye(3))
    ex=coupled_stable.StableExperiment.__new__(coupled_stable.StableExperiment)
    ex.env=None;ex.arrival=166099.3073618387;ex.end=ex.arrival+46800
    exact=np.zeros(8);exact[6]=-1
    close=exact.copy();close[6]+=5e-14
    a=ex.command(exact);b=ex.command(close)
    assert np.diff(b.times).min()>.001
    dates=np.r_[np.linspace(0,ex.end,1000),39600+np.linspace(-120,120,101)]
    np.testing.assert_allclose(a(dates).as_matrix(),b(dates).as_matrix(),atol=1e-12)
    np.testing.assert_allclose(a(dates,2),b(dates,2),atol=1e-12)
    assert np.linalg.norm(b(dates,2),axis=1).max()<1e-6


def test_result_export_normalizes_numpy_scalars_without_changing_values():
    from .coupled_finalize import native
    data=native({'passed':np.bool_(False),'distance':np.float64(100.),'count':np.int64(361)})
    assert json.loads(json.dumps(data))=={'passed':False,'distance':100.,'count':361}


def test_coupled_derivative_is_repeatable_and_improves_the_same_baseline_map():
    b=read('benchmark');s=read('optimization')
    assert b['repeatability']['maximum_position_m']==0
    assert b['directional_convergence']['half_step_prediction_error_m']<.1
    r=s['response_comparison']
    assert r['coupled_maximum_position_error_m']<r['gravity_only_maximum_position_error_m']/50
    assert s['derivative_reuse_check']['maximum_frame_element_difference']<1e-12
    assert len(s['basis']['names'])==8


def test_failed_branches_and_unpaid_cycle_obligations_are_retained():
    old=read('optimization_attempt');s=read('optimization');v=read('validation')
    assert not old['completed']
    assert old['iterations'][0]['replay']['overloaded_members']>0
    assert not old['iterations'][0]['decision']['accepted']
    assert s['interrupted_attempt_charge']['additional_partial_allowance_s']==90
    assert s['final_cycle_debt']['dates_s'][-1]>s['final_cycle_debt']['dates_s'][1]
    for data in [b for b in [old,s,v]]:
        assert not data['accepted_return'] and not data['accepted_cycle']
    assert v['initial_service_coverage']['passed']
    assert not v['next_service_executed'] and not v['subsequent_departure_executed']
    assert v['replacement_service_required_from_s']==21600
    assert v['held_array_comparison_TW']==238.35 and v['target_TW']==23.835


def test_executed_account_does_not_credit_restoration_tail_or_unspecified_supply():
    v=read('validation');a=v['account']
    np.testing.assert_allclose(a['electrical_energy_J'],a['attitude_energy_J']+a['translation_energy_J'],rtol=1e-12)
    np.testing.assert_allclose(a['electrical_energy_J'],sum(r['electrical_energy_J'] for r in v['runs']),rtol=1e-12)
    assert a['duration_s']==v['terminal']['time_s']
    assert np.asarray(v['terminal']['state_m_m_s']).shape==(361,6)
    assert a['overloaded_members']==0 and a['minimum_installed_to_peak']>1
    for k in ['actual_generation_W','actual_bus_delivery_W','actual_storage_state_J',
              'actual_delivery_losses_J','actual_unmet_demand_J','full_cycle_energy_J']:
        assert a[k] is None
    assert not v['stopping_clearance']['passed']


def test_saved_products_bind_sources_and_restart_survives_without_raw_runs():
    manifest=json.loads((HERE/'coupled_checks.json').read_text())
    for entry in manifest['products']:
        path=ROOT/entry['path'];assert digest(path)==entry['sha256']
        data=json.loads(path.read_text())
        for name,expected in data['producer']['source_hashes'].items():assert digest(ROOT/name)==expected
        for name,expected in data['producer']['inputs'].items():assert digest(HERE/'results'/name)==expected
        assert not constants_changed(data['producer']['constants'])
    from .coupled_restart import restore
    state,arrays=restore()
    assert arrays['initial_optimization_state'].shape==(361,6)
    assert arrays['endpoint_state_jac'].shape==(361,6,8)
    np.testing.assert_allclose(arrays['x'],read('optimization')['selected_x'],atol=0,rtol=0)
    budget=manifest['budget'];assert budget['charged_cpu_s']<=1800
    assert budget['new_raw_MiB']<512 and budget['peak_rss_MiB']<2048
