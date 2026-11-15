"""Arithmetic, units and binding checks for the sea-salt supply screen."""
import json

import pytest

from resources import sea_salt as s
from shared.provenance import constants_changed

PRODUCT = json.loads(s.OUT.read_text())
VALUES = json.loads(s.SOURCES.read_text())["values"]


def test_product_binds_inputs_and_constants():
    assert PRODUCT["schema"] == s.SCHEMA
    for path, sha in {**PRODUCT["producer"]["files"], **PRODUCT["producer"]["inputs"]}.items():
        assert s.digest(s.ROOT / path) == sha, path
    assert not constants_changed(PRODUCT["producer"]["constants"])
    assert PRODUCT == json.loads(json.dumps(s.run()))


def test_seawaters_ions_stay_within_its_salinity():
    sea = VALUES["seawater_g_per_kg"]
    major = sum(sea[k] for k in ("cl", "na", "so4", "mg", "ca", "k"))
    assert 0.0 < sea["total"] - major < 0.3            # bicarbonate, bromide, strontium, boron and fluoride
    assert sea["cl"] / sea["total"] == pytest.approx(0.55, abs=0.005)


def test_the_need_matches_the_sea_chemistrys_salt():
    """The sea chemistry study gives 174,000 Gt of salt for the South Pole-Aitken Sea's 4.84 million km3."""
    own = s.need(4.84e6, VALUES)
    assert own["salt_kg"] == pytest.approx(1.74e17, rel=5e-3)
    assert own["cl_kg"] / own["salt_kg"] == pytest.approx(VALUES["seawater_g_per_kg"]["cl"] / 35.165)
    half = s.need(4.84e6, VALUES, 0.5)
    assert half["cl_kg"] == pytest.approx(0.5 * own["cl_kg"])


def test_riding_chlorine_follows_grade_times_mass():
    screen = json.loads(s.SCREEN.read_text())
    riding = s.riding_chlorine(VALUES, screen)
    water = screen["inventories"]["water"]
    standing = water["seas_kg"] + water["lakes_kg"]
    ice = riding["water_from_comet_like_ice"]
    assert ice["cl_kg"]["central"] == pytest.approx(2.8e-4 * standing)
    rock = riding["water_from_ci_like_rock"]
    assert rock["water_share_of_rock"] == pytest.approx(0.166, abs=0.001)       # 18,600 ppm of hydrogen as water
    assert rock["cl_kg"]["central"] == pytest.approx(rock["rock_kg"] * 717e-6)
    n2 = screen["inventories"]["n2_kg"]
    air = riding["nitrogen_from_ci_like_rock"]
    assert air["cl_kg"]["central"] == pytest.approx(717.0 / 2500.0 * n2)
    assert air["cl_kg"]["low"] < air["cl_kg"]["central"] < air["cl_kg"]["high"]
    assert air["cl_kg_at_the_disputed_low_grade"] == pytest.approx(115.0 / 2500.0 * n2)
    comet = riding["nitrogen_from_comet_like_ice"]
    assert comet["cl_kg_floor"][0] < comet["cl_kg_floor"][1] < comet["cl_kg_gas_over_gas"]


def test_rock_for_the_chloride_brings_more_sodium_than_the_sea_needs():
    own = PRODUCT["needs"][s.SEA]
    rock = s.rock_for_the_chloride(own["cl_kg"], own["na_kg"], VALUES)
    assert rock["ci_like_rock_kg"]["at_717_ppm"] == pytest.approx(own["cl_kg"] / 717e-6)
    assert rock["sodium_over_the_seas_need"] == pytest.approx(
        rock["na_over_cl_in_ci_rock"] / rock["na_over_cl_in_seawater"], rel=1e-3)
    assert rock["sodium_over_the_seas_need"] > 10.0
    low, high = rock["lunar_ground_for_the_chloride_share_of_the_moon"]
    assert low == pytest.approx(own["cl_kg"] / 70e-6 / s.MOON_MASS) and high == pytest.approx(10.0 * low)
    assert s.MOON_MASS == pytest.approx(7.346e22, rel=1e-3)


def test_near_earth_asteroids_scale_with_their_mass_and_stay_far_below_the_need():
    asteroids = PRODUCT["near_earth_asteroids"]
    for pool in s.ASTEROID_POOLS:
        row = asteroids["pools"][pool]
        median = row["body_mass_kg_p05_p50_p95"][1]
        assert row["at_the_median"]["cl_kg"]["ci"] == pytest.approx(median * 717e-6)
        assert row["at_the_median"]["share_of_the_seas_need"]["ci"] < 1e-3
    feed = asteroids["feed_through_year_200"]
    assert feed["low"]["cl_kg"]["ci"] < feed["high"]["cl_kg"]["ci"]
    assert feed["high"]["seawater_km3"]["ci"] == pytest.approx(
        feed["high"]["cl_kg"]["ci"] / asteroids["cl_kg_per_km3_of_seawater"])


def test_the_salt_is_a_small_share_of_the_delivery():
    d = PRODUCT["delivery"]
    assert 0.0 < d["salt_share_of_the_delivery"][0] < d["salt_share_of_the_delivery"][1] < 0.02
    salt = PRODUCT["needs"][s.SEA]["salt_kg"]
    assert d["years_of_the_stream"][1] == pytest.approx(salt / d["delivery_flow_kg_s"][0] / s.YEAR_S)


def test_every_value_names_its_evidence():
    def walk(node, path):
        if isinstance(node, dict):
            if any(isinstance(v, (int, float)) for v in node.values()):
                assert "evidence" in node, path
            else:
                for key, value in node.items():
                    walk(value, f"{path}/{key}")
    for name, entry in VALUES.items():
        walk(entry, name)
    keys = {src["key"] for src in json.loads(s.SOURCES.read_text())["sources"]}
    assert {"LoddersFegley2023", "Rubin2019", "Gargano2020", "LunarSourcebook"} <= keys
