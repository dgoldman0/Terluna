"""Checks of the scene descriptions' product: provenance, the frames' geometry and what each frame reports."""
import hashlib
import json
import math

import numpy as np
import pytest

from research.studies.sea_appearance import scenes as S
from shared.provenance import constants_changed


def test_the_committed_scenes_are_current():
    product = json.loads(S.RESULT.read_text())
    assert product["schema"] == "terluna.research.sea-scenes/1"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    for path, expected in product["inputs"].items():
        if (S.ROOT / path).exists():
            assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == expected, path


def test_every_regime_moment_has_a_frame_and_its_sources_where_the_regimes_put_them():
    product = json.loads(S.RESULT.read_text())
    regimes = json.loads(S.REGIMES.read_text())["scenes"]
    assert len(product["scenes"]) == len(regimes) == len(S.VIEWS) == 24
    for scene, regime in zip(product["scenes"], regimes):
        assert (scene["coast"], scene["moment"]) == (regime["coast"], regime["moment"])
        shares = scene["light"]["share_of_frame"]
        assert sum(shares.values()) == pytest.approx(1.0, abs=0.002)
        # Cameras moved at most 5.5 km from the regime points: the Sun and the Earth shift by under 0.2 degrees.
        assert scene["sun"]["elevation_deg"] == pytest.approx(regime["sun"]["elevation_deg"], abs=0.2)
        assert scene["earth"]["elevation_deg"] == pytest.approx(regime["earth"]["elevation_deg"], abs=0.2)
        if scene["sun"]["in_frame"] and scene["sun"]["elevation_deg"] > 0:
            assert scene["sun"]["y"] < scene["camera"]["horizon_row_at_centre"]


def test_the_frame_projects_its_own_rays_back_to_their_pixels():
    frame = S.Frame(262.0, 1.0)
    x = np.array([0.0, 400.0, 960.0, 1919.0])
    y = np.array([0.0, 700.0, 540.0, 1079.0])
    px, py = frame.project(frame.direction(x, y))
    np.testing.assert_allclose(px, x, atol=1e-6)
    np.testing.assert_allclose(py, y, atol=1e-6)
    centre = frame.direction(np.array([960.0]), np.array([540.0]))[0]
    assert math.degrees(math.asin(centre[2])) == pytest.approx(1.0, abs=0.05)


def test_colour_words_follow_the_display_colour():
    assert S.colour_name("#ffffff") == "white"
    assert S.colour_name("#ff9b00").endswith("orange")
    assert "aqua" in S.colour_name("#cdece7") or "cyan" in S.colour_name("#cdece7")
    assert S.to_hex(np.array([0.9505, 1.0, 1.089]) * 0.5, 1.0) == S.to_hex(np.array([0.9505, 1.0, 1.089]) * 0.5, 1.0)
