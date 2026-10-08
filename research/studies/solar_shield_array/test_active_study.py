import json

import numpy as np
import pytest

from shared.provenance import constants_changed,constants_used
from .fleet_run import HERE,ROOT,digest


def test_active_products_bind_sources_and_keep_global_coverage_open():
    for path in (HERE/'results').glob('active_*.json'):
        p=json.loads(path.read_text())
        if 'producer' not in p:continue  # Historical rejected search trial.
        sources=p['producer']['source_hashes']
        for name,expected in sources.items():assert digest(ROOT/name)==expected,(path.name,name)
        assert not constants_changed(p['producer']['constants'])
        assert p['producer']['constants']==constants_used(sources)
        for name,expected in p['producer'].get('input_products',{}).items():
            assert digest(HERE/'results'/name)==expected
        accepted=p.get('acceptance',p)
        assert accepted.get('global_coverage_demonstrated',accepted.get('global_four_lunar_radius_coverage_demonstrated')) is False
        assert accepted.get('delivered_power_W',accepted.get('delivered_electrical_power_demonstrated_W')) is None


def test_control_accounting_uses_integrated_thrust_and_reports_failed_trials():
    path=HERE/'results/active_retime.json'
    if not path.exists():pytest.skip('Retiming calculation has not run')
    p=json.loads(path.read_text());prop=p['propulsion_assumptions']
    cant=np.cos(np.deg2rad(prop['cant_deg']))
    for c in p['cases'].values():
        dv=np.linalg.norm(np.array(c['parameters_m_s']).reshape(-1,3),axis=1).sum()
        assert c['delta_v_m_s']==pytest.approx(dv)
        assert c['propellant_kg']==pytest.approx(c['optical_mass_kg']*dv/(prop['exhaust_velocity_m_s']*cant))
        assert c['ideal_electric_energy_J']==pytest.approx(c['propellant_kg']*prop['exhaust_velocity_m_s']**2/(2*prop['efficiency']))
        assert c['replay_agreement_within_50m']==(c['independent_replay_position_error_m']<50)
        assert c['minimum_earth_beam_margin_deg']>0
    failed=json.loads((HERE/'results/active_search.json').read_text())['rejected_trials']
    assert any(c['independent_replay_position_error_m']>50 for c in failed)


def test_local_acquisition_preserves_its_initial_gaps_and_collision_check():
    path=HERE/'results/active_acquisition.json'
    if not path.exists():pytest.skip('Acquisition calculation has not run')
    p=json.loads(path.read_text())
    assert p['coverage'][0]['ray_coverage']<.99
    assert p['coverage'][-1]['ray_coverage']==pytest.approx(1,abs=1e-10)
    assert p['first_sample_with_full_sun_coverage_minutes']>0
    assert p['collisions']['possible_conflict_pairs']==0
    assert p['collisions']['minimum_proven_axis_clearance_m']>p['collisions']['required_clearance_m']
