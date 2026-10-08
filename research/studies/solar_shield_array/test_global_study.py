import json
import numpy as np
import pytest

from shared.provenance import constants_used,constants_changed
from .fleet_run import HERE,ROOT,digest


def test_rejected_schedule_does_not_claim_its_interrupted_checks():
    p=json.loads((HERE/'results/active_global_rejected.json').read_text())
    src=p['producer']['source_hashes']
    assert all(digest(ROOT/name)==h for name,h in src.items())
    assert p['producer']['constants']==constants_used(src)
    assert p['producer']['input_products']['holding.json']==digest(HERE/'results/holding.json')
    assert p['status']=='REJECTED_HIGH_RECURRING_COST'
    assert not p['interrupted_global_coverage']['completed']
    assert p['interrupted_global_coverage']['last_day']<30
    assert p['compact_rejected']['return_pole']['actual_square_intersections']>0
    assert p['separation']['passes']
    assert not p['acceptance']['low_energy_active_gap_repair_demonstrated']
    low=p['preliminary_budget']['mean_power_bounds_W'][0]
    assert low>p['held_comparison']['mean_power_TW']*1e12


def test_natural_release_keeps_missing_repairs_out_of_its_power_account():
    path=HERE/'results/natural_traffic.json'
    if not path.exists():pytest.skip('Natural release has not completed')
    p=json.loads(path.read_text())
    src=p['producer']['source_hashes']
    assert all(digest(ROOT/name)==h for name,h in src.items())
    assert p['producer']['constants']==constants_used(src)
    assert not constants_changed(p['producer']['constants'])
    assert p['independent_replay']['interpolation_within_50m']
    assert p['independent_replay']['geometry_prefix_maximum_position_difference_m']<50
    assert p['coverage']['horizon_days']==p['independent_replay']['geometry_horizon_days']
    assert max(r['day'] for r in p['coverage']['series'])<=p['independent_replay']['geometry_horizon_days']
    assert not p['acceptance']['low_energy_gap_repair_demonstrated']
    assert p['acceptance']['delivered_electrical_power_demonstrated_W'] is None
    assert 'not a complete fleet operating budget' in p['limitations'][0]


def test_natural_collision_replay_confirms_identified_members_without_grid_interpolation():
    path=HERE/'results/natural_collision_replay.json'
    p=json.loads(path.read_text())
    src=p['producer']['source_hashes']
    assert all(digest(ROOT/name)==h for name,h in src.items())
    assert p['producer']['constants']==constants_used(src)
    assert p['producer']['input_product_sha256']==digest(HERE/'results/natural_traffic.json')
    assert sum(x['confirmed_square_intersections'] for x in p['series'])>0
    assert all(min(x['minimum_beam_margins_deg'])>0 for x in p['series'])
