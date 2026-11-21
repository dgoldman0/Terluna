"""Independent limiting cases and integrity checks for the engineering screen."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("enclosed_feasibility_screen", HERE / "feasibility.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_encircled_power_origin_and_airy_first_zero():
    assert m.encircled_power(0) == 0
    # Standard first J1 zero: 83.78% of uniformly illuminated aperture power.
    assert m.encircled_power(3.8317059702075125) == pytest.approx(0.8377848692, abs=1e-10)


def test_capture_inversion_and_geometric_scaling():
    radius = m.capture_radius(0.8, 0.12, 1e7, 1200)
    x = math.pi * 1200 * radius / (0.12 * 1e7)
    assert m.encircled_power(x) == pytest.approx(0.8, abs=1e-12)
    assert m.capture_radius(0.8, 0.12, 2e7, 1200) == pytest.approx(2 * radius)
    assert m.capture_radius(0.8, 0.12, 1e7, 2400) == pytest.approx(radius / 2)


def test_capture_rejects_nonphysical_inputs():
    with pytest.raises(ValueError):
        m.capture_radius(1, 0.12, 1e7, 1200)


def test_radiator_heat_balance_and_temperature_scaling():
    sigma = 5.670374419e-8  # Independent SI reference used only for this control.
    area = m.radiator_area(1e9, 350, 0.9, sigma, 1.25)
    assert area * 0.9 * sigma * 350 ** 4 == pytest.approx(1.25e9)
    assert m.radiator_area(1e9, 700, 0.9, sigma, 1.25) == pytest.approx(area / 16)
    assert m.radiator_area(1e9, 350, 0.9, sigma, 1.25, 0.5) == pytest.approx(2 * area)


def test_traffic_probability_complements_and_small_limit():
    assert m.loss_probability(0, 0.1) == 0
    assert m.loss_probability(2, 0.1) == pytest.approx(0.19)
    assert m.loss_probability(10000, 1e-12) == pytest.approx(1e-8, rel=1e-8)


def test_refurbishment_no_repair_and_mass_break_even():
    assert m.refurbishment_purchase_ratio(0, 0.2, 0.5) == 1
    assert m.refurbishment_purchase_ratio(0.48, 0.2, 0.2) == pytest.approx(1)
    assert m.refurbishment_purchase_ratio(0.48, 0.2, 1) == pytest.approx(0.616)
    assert m.refurbishment_purchase_ratio(0.48, 0.2, 0) > 1


def test_input_hash_mismatch_is_an_error(tmp_path):
    path = tmp_path / "fixture.json"
    path.write_text('{"fixture": true}')
    entry = {"path": path.name, "sha256": m.sha256(path)}
    assert m.read_pinned(tmp_path, entry)["fixture"]
    path.write_text('{"fixture": false}')
    with pytest.raises(ValueError, match="hash mismatch"):
        m.read_pinned(tmp_path, entry)


def test_product_provenance_and_unchanged_case_admission():
    product = json.loads((HERE / "results/feasibility.json").read_text())
    config = json.loads((HERE / "feasibility_inputs.json").read_text())
    assert product["producer"]["sha256"] == m.sha256(HERE / "feasibility.py")
    for entry in product["dependencies"].values():
        assert m.sha256(m.ROOT / entry["path"]) == entry["sha256"]
    assert product["admission"] == {"combined": {"total": 81, "passed": 81},
                                     "combined_poor_source": {"total": 6, "passed": 6},
                                     "combined_rollout": {"total": 6, "passed": 4}}
    assert len(product["portfolio_ranges_year200"]) == 9
    assert len(product["requirements"]) == 93
    for case in product["requirements"]:
        assert [r["year"] for r in case["checkpoints"]] == config["checkpoints"]


def test_mission_counts_energy_storage_and_pressure_conservation():
    product = json.loads((HERE / "results/feasibility.json").read_text())
    cfg = json.loads((HERE / "feasibility_inputs.json").read_text())
    data = m.read_pinned(m.ROOT, cfg["inputs"]["engineering"])
    p = data["industrial"]["power"]
    ref = next(c for c in product["requirements"] if c["name"] == product["reference_case"])
    for row, screen in zip(ref["checkpoints"], product["checkpoint_screens"]):
        assert screen["LEO_100t_launches_year"] * data["industrial"]["transport"]["launch_payload_kg"] == pytest.approx(row["LEO_kg_year"])
        delivered = screen["beam_installed_DC_W"] * 2 / p["stations"] * p["end_of_life_factor"] * p["RF_efficiency"] * p["beam_capture"] * p["receiver_efficiency"]
        assert delivered == pytest.approx(row["ground_firm_W"])
        assert screen["beam_station_RF_W_at_end_of_life"] + screen["beam_station_RF_conversion_heat_W"] == pytest.approx(screen["beam_installed_DC_W"] / p["stations"] * p["end_of_life_factor"])
        assert screen["usable_firm_load_hours"] * row["ground_firm_W"] == pytest.approx(screen["reserve_Wh"] * product["storage_basis"]["dispatch_fraction_after_margin"])
        assert screen["inherited_critical_load_hours"] == pytest.approx(product["storage_basis"]["target_hours"])
    for beam in product["beam"]["ideal_focus_cases"]:
        for site in beam["site_comparisons"]:
            assert site["flux_sized_area_m2"] * p["receiver_mean_flux_W_m2"] == pytest.approx(site["delivered_MW"] * 1e6)
    for dome in product["crater_domes"]:
        assert dome["shield_3m_weight_N"] / dome["uplift_N"] == pytest.approx(dome["shield_3m_uplift_fraction"])
    assert product["pressure_hall"]["end_force_N"] == pytest.approx(2.864896880624872e8)


def test_reference_opening_does_not_hide_seed_or_operating_peaks():
    product = json.loads((HERE / "results/feasibility.json").read_text())
    ref = next(c for c in product["requirements"] if c["name"] == product["reference_case"])
    assert product["prestart"]["seed_feed_to_opening_annual_feed_ratio"] > 10
    terminal = next(r for r in ref["checkpoints"] if r["year"] == 200)
    assert ref["operating_peaks"]["lunar_raw_kg_year"]["value"] > 30 * terminal["lunar_raw_kg_year"]
