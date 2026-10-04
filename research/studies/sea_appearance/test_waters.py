"""Checks of the water-colour product: provenance, the design guesses it used and what it reports."""
import hashlib
import json

from research.studies.sea_appearance import waters as W
from shared.provenance import constants_changed


def test_the_committed_colours_are_current():
    product = json.loads(W.RESULT.read_text())
    assert product["schema"] == "terluna.research.sea-water-colours/1"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((W.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    assert product["inputs"]["biosphere/living_water/waters.json"] == hashlib.sha256(W.WATERS.read_bytes()).hexdigest()


def test_clearer_water_lets_light_deeper_and_dusk_is_greener():
    product = json.loads(W.RESULT.read_text())
    mean = {(c["water"], c["soil"]): c for c in product["cases"] if c["phase"] == "mean"}
    depth = [mean[(w, "12001")]["moon"]["depth_one_percent_photons_m"] for w in ("open_sea", "productive_coast", "river_mouth")]
    assert depth[0] > depth[1] > depth[2]
    open_sea = {c["phase"]: c for c in product["cases"] if c["water"] == "open_sea"}
    assert open_sea["dusk"]["chlorophyll_mg_m3"] > open_sea["mean"]["chlorophyll_mg_m3"] > open_sea["dawn"]["chlorophyll_mg_m3"]
    assert open_sea["dusk"]["moon"]["chromaticity_xy"][1] > open_sea["dawn"]["moon"]["chromaticity_xy"][1]
    coast = {s: mean[("productive_coast", s)]["moon"]["luminance_cd_m2"] for s in ("10084", "12001", "62231")}
    assert coast["10084"] < coast["12001"] < coast["62231"]
