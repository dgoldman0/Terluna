"""Checks of the scene descriptions' product: provenance, the frames' geometry, that each frame's colours come
from its own sea state and sources, and that scenes.md is the template filled from the product."""
import hashlib
import json
import math

import numpy as np
import pytest

from research.studies.sea_appearance import scenes as S
from shared.provenance import constants_changed


@pytest.fixture(scope="module")
def product():
    return json.loads(S.RESULT.read_text())


def test_the_committed_scenes_are_current(product):
    assert product["schema"] == "terluna.research.sea-scenes/2"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    for path, expected in product["inputs"].items():
        if (S.ROOT / path).exists():
            assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == expected, path


def test_every_regime_moment_has_a_frame_and_its_sources_where_the_regimes_put_them(product):
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


def test_each_frames_colours_come_from_its_own_waves_and_sources(product):
    for scene in product["scenes"]:
        basis = scene["colour_basis"]
        assert basis["slope_covariance_east_north"] == scene["waves"]["slope_covariance_east_north"]
        assert np.trace(basis["slope_covariance_east_north"]) == pytest.approx(
            scene["waves"]["mean_square_slope"]["total"], abs=1e-5)
        assert basis["sun"][0] == pytest.approx(scene["sun"]["elevation_deg"], abs=0.01)
        assert basis["earth"][0] == pytest.approx(scene["earth"]["elevation_deg"], abs=0.01)
        assert basis["camera_lon_lat_deg"] == scene["camera"]["lon_lat_deg"]


def test_disk_sizes_follow_the_rectilinear_lens(product):
    focal = S.WIDTH / 2 / math.tan(math.radians(S.FOV_DEG) / 2)
    assert math.degrees(1 / focal) == pytest.approx(0.0380, abs=1e-4)
    for scene in product["scenes"]:
        for body in ("sun", "earth"):
            d = scene[body]
            centre = 2 * math.tan(math.radians(d["diameter_deg"] / 2)) * focal
            assert d["diameter_px_at_centre"] == pytest.approx(centre, abs=0.06)
            if d["in_frame"]:
                # Away from the axis a rectilinear lens stretches, by at most 1/cos^2 of the angle off the axis.
                off = math.atan(math.hypot(d["x"] - S.WIDTH / 2, d["y"] - S.HEIGHT / 2) / focal)
                for size in (d["width_px"], d["height_px"]):
                    assert centre * 0.99 <= size <= centre / math.cos(off) ** 2 * 1.01
    sun = next(s["sun"] for s in product["scenes"] if s["coast"] == "Smythii headland" and s["moment"] == "low Sun")
    assert sun["width_px"] == pytest.approx(14.2, abs=0.2)


def test_recorded_stars_are_listed_in_full(product):
    for scene in product["scenes"]:
        stars = scene["stars"]
        assert len(stars["recorded"]) == stars["long_exposure_shows"]
        assert all(s["long_exposure_shows_it"] for s in stars["recorded"])


def test_the_description_is_its_template_filled_from_the_product(product):
    assert S.TEXT.read_text() == S.render(product, S.TEMPLATE.read_text())


def test_the_frame_projects_its_own_rays_back_to_their_pixels():
    frame = S.Frame(262.0, 1.0)
    x = np.array([0.0, 400.0, 960.0, 1919.0])
    y = np.array([0.0, 700.0, 540.0, 1079.0])
    px, py = frame.project(frame.direction(x, y))
    np.testing.assert_allclose(px, x, atol=1e-6)
    np.testing.assert_allclose(py, y, atol=1e-6)
    centre = frame.direction(np.array([960.0]), np.array([540.0]))[0]
    assert math.degrees(math.asin(centre[2])) == pytest.approx(1.0, abs=0.05)
    corner = frame.direction(np.array([0.0]), np.array([0.0]))
    assert float(frame.pixel_solid_angle(corner)[0]) < float(frame.pixel_solid_angle(centre)[0])


def test_colour_words_follow_the_display_colour():
    assert S.colour_name("#ffffff") == "white"
    assert S.colour_name("#ff9b00").endswith("orange")
    assert "aqua" in S.colour_name("#cdece7") or "cyan" in S.colour_name("#cdece7")
