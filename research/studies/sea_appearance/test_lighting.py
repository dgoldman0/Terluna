"""Checks of the lighting calendar: the local frame, the phase law, the modes of vision, the sky's readout and the
committed product."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from research.studies.sea_appearance import lighting as L
from shared.provenance import constants_changed


def test_directions_overhead_and_on_the_eastern_horizon():
    frame = L.local_frame(30.0, 20.0)
    el, _ = L.horizontal(L.unit(30.0, 20.0), frame)
    assert el == pytest.approx(90.0)
    el, az = L.horizontal(frame[0], frame)
    assert el == pytest.approx(0.0, abs=1e-12) and az == pytest.approx(90.0)
    el, az = L.horizontal(frame[1], frame)
    assert az == pytest.approx(0.0) or az == pytest.approx(360.0)


def test_the_earth_stands_high_over_the_sub_earth_point():
    jd = np.array([L.JD_START + 10.0])
    g = L.geometry(jd, 0.0, 0.0)
    assert g["earth_elevation"][0] > 80
    assert 0 <= g["earth_phase_angle"][0] <= 180 and g["earth_sunlight_factor"][0] == pytest.approx(1.0, abs=0.04)


def test_full_earth_gives_the_observed_visual_earthlight_above_the_air():
    """Robinson et al. (2025): I/F = 0.23/pi at full phase in the visual band; photopically a little less, as
    Glenar et al.'s Earth is bluer than the Sun."""
    if not (L.SPHERICAL / "moon_1.2atm_standard.npz").exists():
        pytest.skip("The cached scattering solution lives on the research drive")
    sky = L.Sky()
    light = L.Earthlight(sky)
    weights = light.weights(np.array([0.0, 60.0]), np.full(2, L.EARTH_MOON_DISTANCE), np.ones(2))
    lux = weights @ sky.channel_y
    visual = 0.23 * (6_371_000.0 / L.EARTH_MOON_DISTANCE) ** 2 * sky.unfiltered_above_air
    assert lux[0] == pytest.approx(visual, rel=0.1) and lux[0] < visual
    assert lux[1] / lux[0] == pytest.approx(0.44, abs=0.02)


def test_modes_of_vision_for_an_eighteen_percent_grey_surface():
    day, night = 5.0 * math.pi / L.GREY, 0.005 * math.pi / L.GREY
    assert list(L.regime([day * 1.01, day * 0.99, night * 1.01, night * 0.99])) == \
        ["photopic", "mesopic", "mesopic", "scotopic"]


def test_the_nearest_coast_is_water_of_the_body_beside_land():
    labels = np.zeros((8, 16), int)
    labels[2:6, 4:12] = 7
    lon, lat = L.nearest_coast(labels, 7, 0.0, 88.0, ppd=0.08)
    row, col = int(round((90 - lat) * 0.08 - 0.5)), int(round(((lon % 360) * 0.08) - 0.5))
    assert labels[row, col] == 7


def test_the_sky_reproduces_its_stored_samples():
    if not (L.SPHERICAL / "moon_1.2atm_standard.npz").exists():
        pytest.skip("The cached scattering solution lives on the research drive")
    sky = L.Sky()
    assert sky.readout_check < 1e-9
    assert sky.last_order_fraction < 1e-3
    e = np.arange(-90, 90.5, 0.5)
    assert np.all(np.diff(sky.ground(e)) > 0)


def test_the_committed_calendar_is_current():
    product = json.loads(L.RESULT.read_text())
    assert product["schema"] == "terluna.research.sea-lighting-calendar/1"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((L.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    assert product["earth_over_seas"]["sha256"] == hashlib.sha256(L.MAP.read_bytes()).hexdigest()
    assert len(product["sites"]) == 6
    for site in product["sites"]:
        assert sum(site["hours_by_mode"].values()) == len(product["time_hours"])
