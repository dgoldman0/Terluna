"""Numerical, limiting-case and provenance controls; no engineering validation."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

from research.studies.world_economy import industrial_physics as model


def loading(start=0, thresholds=(.1, 100, 1000)):
    return model.loading_screen(model.MOON_RADIUS, model.MOON_SURFACE_GRAVITY,
        3.1413494805269903e18, 121590, 500, start, thresholds,
        model.JULIAN_DAY * model.JULIAN_YEAR_DAYS)


def test_loading_hydrostatic_bound_and_reference_identity():
    r = loading()
    assert r["thin_equivalent_final_pressure_Pa"] == pytest.approx(134508.98955579198)
    assert r["spherical_mass_factor"] == pytest.approx(1.1062504281256023)
    for item in r["thresholds"]:
        thin = item["thin_column"]
        calibrated = item["reference_column_scaled"]
        assert thin["retained_mass_kg"] < calibrated["retained_mass_kg"]
        assert calibrated["retained_mass_kg"] / r["final_mass_kg"] == pytest.approx(item["threshold_Pa"] / 121590)


def test_vacuum_threshold_is_tiny_fraction_not_centuries():
    row = loading()["thresholds"][0]["reference_column_scaled"]
    assert row["fraction_of_final_mass"] == pytest.approx(8.224360555966774e-7)
    assert row["elapsed_loading_seconds"] / 3600 == pytest.approx(3.604737231680237)


def test_delay_preserves_end_and_accelerates_subsequent_loading():
    immediate, delayed = loading(), loading(100)
    assert delayed["end_year"] == immediate["end_year"]
    assert delayed["retained_loading_kg_s"] == pytest.approx(immediate["retained_loading_kg_s"] * 1.25)
    for first, second in zip(immediate["thresholds"], delayed["thresholds"]):
        a, b = first["reference_column_scaled"], second["reference_column_scaled"]
        assert b["elapsed_loading_years"] == pytest.approx(a["elapsed_loading_years"] * .8)
        assert b["crossing_year"] == pytest.approx(100 + b["elapsed_loading_years"])


def test_threshold_above_final_not_reached():
    row = loading(thresholds=(200000,))["thresholds"][0]
    for method in ("thin_column", "reference_column_scaled"):
        assert row[method]["crossing_year"] is None
        assert row[method]["elapsed_loading_seconds"] is None
        assert row[method]["reached_by_end"] is False


@pytest.mark.parametrize("start", [500, 501, -1, math.inf, True])
def test_invalid_loading_start_rejected(start):
    with pytest.raises(ValueError):
        loading(start)


@pytest.mark.parametrize("thresholds", [(), (.1, .1), (-1,), (0,), (math.nan,)])
def test_invalid_thresholds_rejected(thresholds):
    with pytest.raises(ValueError):
        loading(thresholds=thresholds)


def test_impossible_calibrated_pressure_rejected():
    with pytest.raises(ValueError):
        model.loading_screen(1, 1, 1, 1, 500, 0, [.1], 1)


def test_sphere_membrane_force_balance_and_mass_scaling():
    first = model.tensile_shell_screen(25, 121590, 0, 2700, 100e6)
    second = model.tensile_shell_screen(100, 121590, 0, 2700, 100e6)
    t = first["sphere"]["thickness_m"]
    assert 2 * math.pi * 25 * t * 100e6 == pytest.approx(math.pi * 25**2 * 121590)
    assert first["sphere"]["shell_mass_per_volume_kg_m3"] == pytest.approx(4.924395)
    assert second["sphere"]["shell_mass_kg"] / first["sphere"]["shell_mass_kg"] == pytest.approx(64)
    assert second["sphere"]["shell_mass_per_volume_kg_m3"] == first["sphere"]["shell_mass_per_volume_kg_m3"]
    assert first["long_cylinder_barrel"]["thickness_m"] == 2*t


def test_pressure_difference_controls_tension_and_equal_pressure_limit():
    a = model.tensile_shell_screen(25, 121590, 20265, 2700, 100e6)
    b = model.tensile_shell_screen(25, 101325, 0, 2700, 100e6)
    c = model.tensile_shell_screen(25, 121590, 121590, 2700, 100e6)
    assert a["sphere"] == b["sphere"]
    assert c["sphere"]["shell_mass_kg"] == 0
    with pytest.raises(ValueError):
        model.tensile_shell_screen(25, 0, 121590, 2700, 100e6)


def test_buckling_sensitivity_is_square_root_not_linear():
    r = model.external_sphere_screen(25, 121590, 2700, 70e9, .33, .2)
    perfect = r["perfect_shell_thickness_m"]
    reduced = r["assumed_reduced_capacity_thickness_m"]
    assert reduced / perfect == pytest.approx(math.sqrt(5))
    recovered = .2 * 2 * 70e9 / math.sqrt(3*(1-.33**2)) * (reduced/25)**2
    assert recovered == pytest.approx(121590)
    zero = model.external_sphere_screen(25, 0, 2700, 70e9, .33, .2)
    assert zero["perfect_shell_thickness_m"] == 0


@pytest.mark.parametrize("nu,k", [(.5,.2), (-.1,.2), (.33,0), (.33,1.1)])
def test_invalid_buckling_parameters(nu,k):
    with pytest.raises(ValueError):
        model.external_sphere_screen(25, 121590, 2700, 70e9, nu, k)


def test_handling_keeps_inertia_while_reducing_weight():
    r = model.handling_screen(1000, model.MOON_SURFACE_GRAVITY,
        model.STANDARD_GRAVITY, 10, 1, 1, .5)
    assert r["lunar_weight_N"] == pytest.approx(1624.2188765608128)
    assert r["reference_weight_N"] == pytest.approx(9806.65)
    assert r["inertial_force_both_N"] == 1000
    assert r["kinetic_energy_both_J"] == 500
    assert r["lunar_lift_energy_J"] / r["reference_lift_energy_J"] == pytest.approx(r["weight_ratio"])


def test_orbital_energy_includes_final_kinetic_energy():
    r = model.ring_energy_screen(model.MOON_GM, model.MOON_RADIUS, 20e6)
    assert r["mechanical_energy_increase_J_kg"] == pytest.approx(
        r["potential_energy_increase_J_kg"] + r["final_circular_kinetic_energy_J_kg"])
    assert r["mechanical_energy_increase_J_kg"] < r["surface_escape_energy_J_kg"]
    assert r["final_circular_kinetic_energy_J_kg"] == pytest.approx(122570.00295)
    with pytest.raises(ValueError):
        model.ring_energy_screen(model.MOON_GM, model.MOON_RADIUS, model.MOON_RADIUS)


def test_beam_spans_preserve_the_named_constraint_only():
    r = model.beam_span_screen(model.MOON_SURFACE_GRAVITY, model.STANDARD_GRAVITY)
    gravity_ratio = model.MOON_SURFACE_GRAVITY/model.STANDARD_GRAVITY
    assert gravity_ratio*r["equal_absolute_sag_span_ratio"]**4 == pytest.approx(1)
    assert gravity_ratio*r["equal_bending_stress_span_ratio"]**2 == pytest.approx(1)
    assert r["equal_absolute_sag_span_ratio"] == pytest.approx(1.5675, abs=5e-5)
    assert r["equal_bending_stress_span_ratio"] == pytest.approx(2.4572, abs=5e-5)


def test_capillary_length_and_bond_number_reproduce_selected_example():
    moon = model.capillary_screen(2400, .9, model.MOON_SURFACE_GRAVITY, [.003,.020], .2)
    earth = model.capillary_screen(2400, .9, model.STANDARD_GRAVITY, [.003,.020], .2)
    assert moon["capillary_length_m"] == pytest.approx(.015195, abs=1e-6)
    assert earth["capillary_length_m"] == pytest.approx(.006184, abs=1e-6)
    assert [r["bond_number"] for r in moon["bond_numbers"]] == pytest.approx([.038981253,1.732500135], rel=1e-7)
    assert [r["bond_number"] for r in earth["bond_numbers"]] == pytest.approx([.2353596,10.460426667], rel=1e-7)
    for r in moon["bond_numbers"]:
        assert r["bond_number"] == pytest.approx((r["length_m"]/moon["capillary_length_m"])**2)


def test_reversible_pumping_work_limits_and_near_equal_pressures():
    r = model.ideal_pump_work(150,101325,.1)
    assert r["ideal_reversible_work_J"] == pytest.approx(15.1985e6, abs=30)
    assert model.ideal_pump_work(150,101325,0)["ideal_reversible_work_J"] == 15198750
    assert model.ideal_pump_work(150,101325,101325)["ideal_reversible_work_J"] == 0
    assert model.ideal_pump_work(150,101325,101325-1e-8)["ideal_reversible_work_J"] > 0
    with pytest.raises(ValueError):
        model.ideal_pump_work(150,101325,101326)


def test_batch_pumpdown_capacity_ratio_and_zero_limit():
    values = [model.batch_pumpdown_screen(t,p)["ideal_capacity_ratio_if_all_pumpdown_removed"]
              for t in [120,480] for p in [70,16]]
    assert values == pytest.approx([190/120,136/120,550/480,496/480])
    assert model.batch_pumpdown_screen(120,0)["ideal_capacity_ratio_if_all_pumpdown_removed"] == 1


def test_historical_head_rates_are_volume_conversions_and_integer_capacity():
    seconds = model.JULIAN_DAY*model.JULIAN_YEAR_DAYS
    low = model.deposition_head_screen(.00033,2700,.7,1e9,seconds)
    high = model.deposition_head_screen(.0025,2700,.7,1e9,seconds)
    assert low["mass_rate_kg_h"] == pytest.approx(.891)
    assert high["mass_rate_kg_h"] == pytest.approx(6.75)
    assert high["deposited_mass_per_head_kg_year"] == pytest.approx(41419.35)
    assert high["minimum_integer_heads"] == 24144
    for r in [low,high]:
        annual = r["deposited_mass_per_head_kg_year"]
        count = r["minimum_integer_heads"]
        assert annual*count >= 1e9
        assert annual*(count-1) < 1e9
    with pytest.raises(ValueError):
        model.deposition_head_screen(.0025,2700,1.1,1e9,seconds)


def test_accelerator_length_does_not_shrink_with_lunar_weight():
    cases = [model.accelerator_screen(2400,a) for a in [10,100,1000]]
    assert [r["constant_acceleration_length_m"] for r in cases] == [288000,28800,2880]
    assert [r["acceleration_duration_s"] for r in cases] == [240,24,2.4]


def test_product_identity_and_rejection_of_changed_bytes(tmp_path):
    data, products = model.load_inputs()
    assert products["historical_air"]["values"]["inventory.mass_kg"] == products["air_inventory"]["values"]["inventories.air_kg"]
    data["products"]["air_inventory"]["sha256"] = "0"*64
    path = tmp_path / "changed.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Input changed"):
        model.load_inputs(path)


def test_results_are_finite_and_keep_all_four_release_scenarios():
    r = model.results()
    assert r["schema"] == model.SCHEMA
    assert [x["start_year"] for x in r["loading"]] == [0,10,50,100]
    assert len(r["enclosure_cases"]) == 4
    assert r["producer"]["constants"]["MOON_SURFACE_GRAVITY"] == model.MOON_SURFACE_GRAVITY
    json.dumps(r, allow_nan=False)
