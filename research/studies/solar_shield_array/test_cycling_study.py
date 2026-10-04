import json
from pathlib import Path

import pytest

from shared import constants as K
from shared.provenance import constants_changed, constants_used
from .cycling_analysis import digest
from .cycling_search import ROOT, HERE, SCENARIO


def product():
    return json.loads((HERE/"results/cycling.json").read_text())


def test_cycling_product_binds_current_sources_and_comparison_input():
    p = product()
    for name, expected in p["producer"]["source_hashes"].items():
        assert digest(ROOT/name) == expected, name
    assert not constants_changed(p["producer"]["constants"])
    assert p["producer"]["constants"] == constants_used(p["producer"]["source_hashes"])
    assert p["producer"]["holding_product_sha256"] == digest(HERE/"results/holding.json")
    config = json.loads((HERE/"cycling_scenario.json").read_text())
    assert p["areal_mass_kg_m2"] == config["areal_mass_kg_m2"] == SCENARIO["baseline_areal_mass_kg_m2"]
    assert config["protected_lunar_radii"] == SCENARIO["protected_radii"]


def test_bounded_motion_and_independent_return_are_separate_from_formation_control():
    p = product()
    for name in ("inner_filter", "middle_filter", "middle_full", "inner_filter_three_year"):
        case = p["cases"][name]
        assert case["completed"]
        assert case["moon_range_km"][0]*1000 > 4*K.MOON_RADIUS
        assert case["minimum_reflected_beam_clearance_deg"] > 0
        assert case["maximum_sunward_photon_acceleration_m_s2"] < 1e-15
    assert p["cases"]["inner_filter_three_year"]["duration_days"] == 3*K.JULIAN_YEAR_DAYS
    arc = p["independently_replayed_return_arc"]
    assert arc["closure_position_m"] < 50
    assert arc["closure_velocity_m_s"] < .001
    assert arc["gravity_only_counterfactual_closure_position_m"] > 1e6
    assert p["convergence"]["useful_fraction_relative_difference"] < .005
    assert p["convergence"]["maximum_perturbation_separation_m"] > 1000
    assert p["fleet_coverage_demonstrated"] is False
    assert not p["cases"]["outer_full"]["completed"]


def test_inventory_is_an_area_comparison_at_fixed_mass_not_a_solved_fleet():
    p = product()
    mass = p["reference_holding"]["base_optical_mass_kg"]
    for name, inventory in p["reference_area_inventory"].items():
        duty = p["cases"][name]["useful_projected_area_fraction"]
        assert inventory["reference_aperture_area_multiplier"] == pytest.approx(1/duty)
        assert inventory["base_optical_inventory_kg"] == pytest.approx(mass/duty)
    assert .15 < p["cases"]["inner_filter"]["useful_projected_area_fraction"] < .17
