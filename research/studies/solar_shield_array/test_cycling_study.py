import json
from pathlib import Path

import pytest
import numpy as np

from shared import constants as K
from shared.provenance import constants_changed, constants_used
from .cycling_analysis import digest, time_average
from .cycling_search import ROOT, HERE, SCENARIO


def product():
    return json.loads((HERE/"results/cycling.json").read_text())


@pytest.mark.parametrize("values, expected", [
    ([True, True, True, True], 1.),
    ([False, False, False, False], 0.),
    ([False, True, True, False], 5/8),
    ([True, False, False, True], 3/8),
    ([0., .2, .6, 1.], 4.9/8),
])
def test_time_average_counts_boolean_plateaus_on_nonuniform_grid(values, expected):
    # Three intervals of 1, 2 and 5 seconds. Boolean addition must count
    # adjacent True samples twice, just as a constant indicator integrates to 1.
    assert time_average(np.asarray(values), np.array([0., 1., 3., 8.])) == pytest.approx(expected)


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


def test_raw_service_time_bounds_the_weighted_area_statistic():
    for case in product()["cases"].values():
        if "raw_service_time_fraction" in case:
            # The threshold can discard at most 0.001 of weighted area.
            assert 0 <= case["raw_service_time_fraction"] <= 1
            assert case["useful_projected_area_fraction"] <= case["raw_service_time_fraction"] + .001
    assert product()["cases"]["inner_filter"]["raw_service_time_fraction"] == pytest.approx(.15823446649935358)
