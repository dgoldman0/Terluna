"""Published evidence boundaries, source identity and capacity optimization."""
import json
from pathlib import Path
import numpy as np
import pytest
from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def test_pilot_products_bind_current_sources_and_retain_rejected_gates():
    for name in ['pilot_24.json','pilot_48.json','pilot_validation.json','pilot_budget.json']:
        p=json.loads((HERE/'results'/name).read_text())
        for path,expected in p['producer']['source_hashes'].items():
            assert digest(ROOT/path)==expected,path
        assert not constants_changed(p['producer']['constants'])
    p=json.loads((HERE/'results/pilot_validation.json').read_text())
    assert p['accepted_fleet'] is False
    assert p['recurring_correction_measured'] is False
    for name,expected in p['producer']['inputs'].items():
        assert digest(HERE/'results'/name)==expected
    assert p['independent_replay']['passes_50m']
    assert p['constraint_generation']['minimum_old_solution_capacity'] < 1
    assert not p['attitude_witness']['slew_path_safety_certified']
    budget=json.loads((HERE/'results/pilot_budget.json').read_text())
    assert budget['producer']['input_sha256']==digest(HERE/'results/pilot_validation.json')
    assert budget['collection_capacity_W'] is None and budget['delivered_power_W'] is None


def test_phase_population_relaxation_is_nested_and_lp_duals_close():
    p=json.loads((HERE/'results/pilot_validation.json').read_text())
    models=p['phase_population_optimization']
    uniform=models['uniform_solution'];free=models['free_phase_solution'];gated=models['sampled_attitude_gated_solution']
    assert free['inventory_tiles']<=uniform['inventory_tiles']+1
    if gated['success']:
        assert gated['inventory_tiles']>=free['inventory_tiles']-1
    for result in [uniform,free,gated,p['constraint_generation']['refined_solution']]:
        if result['success']:
            assert result['minimum_capacity']>=1-1e-6
            assert result['dual_max_column_price']<=1+1e-6
            assert result['inventory_tiles']/1e6==pytest.approx(result['dual_lower_bound_million_tiles'],abs=1e-6)
