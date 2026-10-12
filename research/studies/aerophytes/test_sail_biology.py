"""Invariant and reading-boundary tests for sail requirements, not viability."""
import math
import json

import pytest

from research.studies.aerophytes import sail_biology as model


def test_local_signal_loop_can_beat_a_fast_remote_loop():
    local = model.propagation(1., .3)
    remote = model.propagation(5000., 100.)
    assert local['round_trip_s'] < remote['one_way_s']
    assert remote['round_trip_s'] == 2*remote['one_way_s']


def test_zero_distance_does_not_create_a_detection_time():
    assert model.propagation(0., 1.)['round_trip_s'] == 0.


def test_hydraulic_stages_keep_the_total_head():
    g = dict(gm_m3_s2=1., radius_m=100.)
    whole = model.liquid_head(10., 20., 1000., **g)['static_head_pa']
    upper = model.liquid_head(4., 20., 1000., **g)['static_head_pa']
    lower = model.liquid_head(6., 16., 1000., **g)['static_head_pa']
    assert whole == pytest.approx(upper+lower)


def test_retrieval_speed_changes_time_not_lifting_energy():
    args = dict(depth_m=5., paid_line_m=6., payload_kg=20., top_height_m=10.,
                gm_m3_s2=1., radius_m=100.)
    slow, fast = (model.retrieval(reel_speed_m_s=v, **args) for v in (1., 2.))
    assert slow['payload_gravitational_energy_j'] == fast['payload_gravitational_energy_j']
    assert slow['travel_s'] == 2*fast['travel_s']
    assert fast['payload_mean_lifting_power_w'] == 2*slow['payload_mean_lifting_power_w']


def test_a_better_hinge_arm_reduces_tendon_force_not_mechanical_work():
    args = dict(load_n_m2=2., area_m2=10., force_arm_m=.5, angle_deg=30., seconds=10.)
    short = model.panel_actuation(tendon_arm_m=.1, **args)
    long = model.panel_actuation(tendon_arm_m=.2, **args)
    assert short['tendon_force_n'] == 2*long['tendon_force_n']
    assert short['mechanical_work_j'] == long['mechanical_work_j']


def test_split_tethers_preserve_total_section_but_increase_width():
    result = model.equal_tether_split(16)
    assert 16*result['diameter_ratio_each']**2 == 1.
    assert result['total_projected_width_ratio'] == 4.


def test_recycling_changes_export_not_new_tissue_production():
    args = dict(area_m2=10., dry_kg_m2=.1, replacements_per_year=2.,
                carbon_fraction=.5, phosphorus_fraction=.001)
    no_recovery, full_recovery = (model.renewal(phosphorus_recovery=r, **args) for r in (0., 1.))
    assert no_recovery['incorporated_carbon_kg_year'] == full_recovery['incorporated_carbon_kg_year']
    assert no_recovery['lost_phosphorus_kg_year'] > 0.
    assert full_recovery['lost_phosphorus_kg_year'] == 0.


@pytest.mark.parametrize('speed', [0., -1., math.nan, math.inf])
def test_invalid_propagation_rate_is_rejected(speed):
    with pytest.raises(ValueError):
        model.propagation(1., speed)


@pytest.mark.parametrize('count', [0, -1, 1.5, True])
def test_invalid_tether_count_is_rejected(count):
    with pytest.raises(ValueError):
        model.equal_tether_split(count)


def test_material_fractions_cannot_exceed_mass():
    with pytest.raises(ValueError):
        model.renewal(10., .1, 1., .9, .2, .8)


def test_no_below_sea_level_hydraulic_endpoint():
    with pytest.raises(ValueError):
        model.potential_difference(100., 10., 1., 100.)


def test_line_cannot_be_shorter_than_drop():
    with pytest.raises(ValueError):
        model.retrieval(10., 5., 1., 1., 10., 1., 100.)


def test_parent_product_change_requires_explicit_review(tmp_path):
    parent = tmp_path/model.PARENT
    parent.parent.mkdir(parents=True)
    parent.write_text('{}\n')
    with pytest.raises(ValueError, match='parent product changed'):
        model.evaluate(tmp_path)


def test_stored_product_regenerates_exactly():
    saved = json.loads((model.ROOT/model.STUDY/'results/sail_biology.json').read_text())
    assert model.evaluate() == saved


@pytest.mark.parametrize('water', [-0.1, math.nan, math.inf, -math.inf])
def test_retained_water_cannot_remove_mass_or_be_nonfinite(water):
    with pytest.raises(ValueError):
        model.wet_load(314159., water)
