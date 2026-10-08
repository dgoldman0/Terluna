"""Checks of the regime product: provenance, the moments it chose and what its scenes report."""
import hashlib
import json

import numpy as np
import pytest

from research.studies.sea_appearance import regimes as G
from shared.provenance import constants_changed


def test_the_committed_regimes_are_current():
    product = json.loads(G.RESULT.read_text())
    assert product["schema"] == "terluna.research.sea-regimes/1"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((G.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    for path, expected in product["inputs"].items():
        if (G.ROOT / path).exists():
            assert hashlib.sha256((G.ROOT / path).read_bytes()).hexdigest() == expected, path
    panoramas = G.ROOT / product["panoramas"]["path"]
    if panoramas.exists():
        assert hashlib.sha256(panoramas.read_bytes()).hexdigest() == product["panoramas"]["sha256"]


def test_each_coast_runs_from_noon_to_its_darkest_hour():
    product = json.loads(G.RESULT.read_text())
    assert len(product["scenes"]) == len(G.COASTS) * len(G.MOMENTS)
    for coast in G.COASTS:
        scenes = {s["moment"]: s for s in product["scenes"] if s["coast"] == coast}
        light = [scenes[m]["light_on_ground_lux"] for m, _ in G.MOMENTS]
        assert light[0] == max(light) and light[-1] == min(light)
        for name, target in G.MOMENTS[1:-1]:
            assert scenes[name]["sun"]["elevation_deg"] == pytest.approx(target, abs=0.6)


def test_the_moment_finder_follows_the_sinking_sun():
    sun = np.r_[np.linspace(-80, 80, 50), np.linspace(80, -80, 50)]
    series = dict(sun_elevation_deg=sun)
    i = G.moment(series, "low Sun", 4.0, np.zeros(100))
    assert sun[i - 1] >= 4.0 > sun[i] and i > 50
