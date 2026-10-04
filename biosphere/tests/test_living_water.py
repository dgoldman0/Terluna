"""The living-water design guesses are complete and each best guess lies inside its range."""
import json
from pathlib import Path

WATERS = Path(__file__).resolve().parents[1] / "living_water/waters.json"


def test_every_guess_lies_in_its_range_and_names_a_known_soil():
    record = json.loads(WATERS.read_text())
    assert record["schema"] == "terluna.biosphere.living-water/1"
    for water in record["waters"].values():
        for key in ("chlorophyll_mg_m3", "dissolved_organic_440_per_m", "fines_g_m3"):
            low, high = water[key]["range"]
            assert low <= water[key]["best"] <= high and low < high
        assert water["soil"] in record["soils"] and water["basis"]
    cycle = record["light_cycle"]["chlorophyll_dusk_over_dawn"]
    assert cycle["range"][0] <= cycle["best"] <= cycle["range"][1]
