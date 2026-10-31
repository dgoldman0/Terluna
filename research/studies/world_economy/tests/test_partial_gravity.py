"""Dimensional, control and provenance checks for the conditional process screen."""
import copy
import hashlib
import json
import math
from pathlib import Path

import pytest

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_changed

from research.studies.world_economy import partial_gravity as model


@pytest.fixture
def inputs():
    return model.load_inputs()


def test_radius_conversion_and_fixed_interval_window(inputs):
    material = inputs["paper_tabulated_material_inputs"]["values"]
    process = inputs["illustrative_casting_process"]["values"]
    window = model.gravity_window(material, process)
    assert window["minimum_for_bubble_arrival_m_s2"] == pytest.approx(1.1153846153846154)
    assert window["maximum_for_particle_uniformity_m_s2"] == pytest.approx(2.047058823529412)
    assert window["interval_nonempty"]
    # Doubling particle diameter quadruples drift and divides its gravity bound by four.
    larger = dict(material, particle_diameter_m=2 * material["particle_diameter_m"])
    assert model.gravity_window(larger, process)["maximum_for_particle_uniformity_m_s2"] == pytest.approx(
        window["maximum_for_particle_uniformity_m_s2"] / 4)


def test_shared_gravity_controls_classification(inputs):
    result = model.calculate(inputs)
    moon, earth, zero = (result["casting"][k] for k in ("moon", "earth", "ideal_microgravity"))
    assert moon["gravity_m_s2"] == MOON_SURFACE_GRAVITY
    assert earth["gravity_m_s2"] == STANDARD_GRAVITY
    assert moon["both_constraints_pass"]
    assert not earth["particle_drift_within_limit"]
    assert earth["bubble_reaches_surface_before_immobilization"]
    assert zero["particle_drift_within_limit"]
    assert not zero["both_constraints_pass"]
    assert zero["bubble_surface_arrival_time_s"] is None
    assert zero["particle_drift_limit_time_s"] is None
    assert zero["bubble_displacement_m"] == zero["particle_displacement_m"] == 0


def test_dimensionless_numbers_use_diameter_and_dynamic_viscosity(inputs):
    result = model.calculate(inputs)["casting"]["moon"]
    assert result["bubble_reynolds"] == pytest.approx(
        1170 * result["bubble_speed_m_s"] * 0.0002 / 0.029)
    assert result["bubble_bond_diameter"] == pytest.approx(
        1170 * MOON_SURFACE_GRAVITY * 0.0002**2 / 0.03)
    assert result["particle_reynolds"] < result["bubble_reynolds"] < 0.01


def test_earth_time_and_viscosity_controls_match_both_displacements(inputs):
    result = model.calculate(inputs)
    moon = result["casting"]["moon"]
    controls = result["earth_matching_controls"]
    for name in ("time_control", "viscosity_control"):
        for quantity in ("bubble_displacement_m", "particle_displacement_m"):
            assert controls[name][quantity] == pytest.approx(moon[quantity])
        assert controls[name]["both_constraints_pass"]
    assert controls["same_viscosity_immobilize_after_s"] < 100
    assert controls["same_duration_dynamic_viscosity_Pa_s"] > 0.029


def test_uniform_viscosity_change_scales_both_gravity_bounds(inputs):
    material = inputs["paper_tabulated_material_inputs"]["values"]
    process = inputs["illustrative_casting_process"]["values"]
    base = model.gravity_window(material, process)
    changed = model.gravity_window(dict(material, dynamic_viscosity_Pa_s=0.058), process)
    for key in ("minimum_for_bubble_arrival_m_s2", "maximum_for_particle_uniformity_m_s2"):
        assert changed[key] == pytest.approx(2 * base[key])
    assert changed["bubble_to_particle_speed_ratio"] == pytest.approx(base["bubble_to_particle_speed_ratio"])


def test_long_constant_viscosity_interval_fails_three_tested_environments(inputs):
    result = model.calculate(inputs)
    long = result["long_interval_sensitivity"]
    assert long["duration_s"] == 2700
    assert all(not row["both_constraints_pass"] for row in long["cases"].values())
    assert long["cases"]["moon"]["particle_displacement_m"] == pytest.approx(
        27 * result["casting"]["moon"]["particle_displacement_m"])
    # An interval remains at other small, nonzero accelerations; do not infer a universal failure.
    assert long["gravity_window"]["interval_nonempty"]


def test_bubble_birth_deadline_exposes_late_gas_limit(inputs):
    result = model.calculate(inputs)["casting"]["moon"]
    assert result["latest_bubble_birth_time_s"] + result["bubble_surface_arrival_time_s"] == pytest.approx(100)
    assert 0 < result["latest_bubble_birth_time_s"] < 50


def test_spin_shape_and_inertial_load_scaling(inputs):
    result = model.calculate(inputs)["spin_paraboloid"]
    moon, earth = (result["cases"][k] for k in ("moon", "earth"))
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    assert moon["angular_speed_rad_s"] / earth["angular_speed_rad_s"] == pytest.approx(math.sqrt(ratio))
    assert moon["rotation_load_proxy_m2_s2"] / earth["rotation_load_proxy_m2_s2"] == pytest.approx(ratio)
    assert moon["edge_sag_m"] == earth["edge_sag_m"] == pytest.approx(5**2 / 80)
    assert not result["cases"]["ideal_microgravity"]["hydrostatic_paraboloid_defined"]
    assert result["cases"]["ideal_microgravity"]["angular_speed_rad_s"] is None


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True])
def test_nonphysical_gravity_rejected(value):
    with pytest.raises(ValueError):
        model.stokes_speed(1000, 0.0001, 1, value)


def test_material_density_and_process_duration_rejected(inputs):
    material = inputs["paper_tabulated_material_inputs"]["values"]
    process = inputs["illustrative_casting_process"]["values"]
    with pytest.raises(ValueError):
        model.casting_screen(dict(material, particle_density_kg_m3=1000), process, 1)
    with pytest.raises(ValueError):
        model.casting_screen(material, dict(process, liquid_interval_s=0), 1)


def test_product_is_reproducible_and_binds_transcription_model_constants():
    first = model.results()
    second = json.loads((model.HERE / "results/partial_gravity.json").read_text())
    assert first == second
    encoded = json.dumps(first, allow_nan=False)
    assert json.loads(encoded) == first
    producer = first["producer"]
    assert not constants_changed(producer["constants"])
    assert producer["constants"] == {
        "MOON_SURFACE_GRAVITY": MOON_SURFACE_GRAVITY, "STANDARD_GRAVITY": STANDARD_GRAVITY}
    for label, digest in producer["files"].items():
        path = Path(model.__file__).with_name(Path(label).name)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_wrong_input_schema_rejected(tmp_path, inputs):
    changed = copy.deepcopy(inputs)
    changed["schema"] = "wrong"
    path = tmp_path / "inputs.json"
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError):
        model.load_inputs(path)
