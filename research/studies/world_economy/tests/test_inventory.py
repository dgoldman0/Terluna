import json
import pytest
from research.studies.world_economy import inventory


def test_committed_inventory_reproduces_and_population_closes():
    actual = inventory.results()
    assert actual == json.loads((inventory.HERE / "results/inventory.json").read_text())
    p = actual["population_billion"]
    assert p["moon_ground"] + p["moon_aloft"] == p["moon"]
    assert sum(p[k] for k in ("earth", "moon", "array", "colonies")) == p["total"]
    assert actual["equal_per_resident_benchmark"]["moon_per_resident_premium_for_equal_total"] == .125


def test_lift_capacity_is_not_counted_as_manufactured_stock():
    d = inventory.results()["derived"]
    assert d["aerial_district_gt"] == 200
    assert d["aerial_envelope_flight_allowance_gt"] == 100
    assert d["aerial_control_allowance_gt"] == 40
    assert d["aerial_unmanufactured_reserve_gt"] == 60
    assert d["array_habitat_gt"] == 5025
    assert d["aerial_hydrogen_gt"] == pytest.approx(29.84)


def test_changed_input_fails_before_silent_substitution(tmp_path):
    manifest = json.loads((inventory.HERE / "inputs.json").read_text())
    name = next(iter(manifest["products"]))
    manifest["products"] = {name:manifest["products"][name]}
    source = tmp_path / name
    source.parent.mkdir(parents=True)
    source.write_text("{}")
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Input changed"):
        inventory.load_inputs(tmp_path, path)
