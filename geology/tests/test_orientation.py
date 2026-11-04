import json
import math
from pathlib import Path

import numpy as np
import pytest

from geology import orientation as ori
from shared.constants import MOON_EQUATOR_TO_ECLIPTIC_DEG, MOON_RADIUS

PRODUCT = Path(ori.__file__).parent / "results" / "moon_orientation.json"
GRIDS = Path("/media/projectspace/terluna-research/crustal-loading/crustal_loading_nmax359.npz")


def coarse_grid(step=1.0):
    lat = np.arange(-90 + step / 2, 90, step)
    lon = np.arange(-180 + step / 2, 180, step)
    return lat, lon


def test_moments_match_the_report():
    assert ori.C_MOMENT - ori.A_MOMENT == pytest.approx(5.50e31, rel=2e-3)
    assert ori.B_MOMENT - ori.A_MOMENT == pytest.approx(1.99e31, rel=3e-3)
    J2 = (ori.C_MOMENT - (ori.A_MOMENT + ori.B_MOMENT) / 2) / ori.MR2
    C22 = (ori.B_MOMENT - ori.A_MOMENT) / (4 * ori.MR2)
    assert J2 == pytest.approx(2.0321e-4, rel=1e-3)          # GRAIL GL0420A, as the report checks
    assert C22 == pytest.approx(2.2380e-5, rel=1e-3)


def test_cassini_relation_reproduces_todays_tilt():
    assert ori.cassini_tilt(ori.C_MOMENT, ori.A_MOMENT) == pytest.approx(MOON_EQUATOR_TO_ECLIPTIC_DEG, abs=0.005)


def test_darwin_radau():
    assert ori.fluid_love(0.4) == pytest.approx(1.5, abs=1e-12)
    assert ori.fluid_love(ori.C_MR2) == pytest.approx(1.436, abs=1e-3)


def test_cap_inertia_matches_point_mass():
    lat, lon = coarse_grid(0.25)
    weight, units = ori.grid_geometry(lat, lon)
    assert weight.sum() / (4 * np.pi) == pytest.approx(1.0, abs=1e-9)
    n0 = ori.unit(40.0, -30.0)
    cap = np.where(np.tensordot(n0, units, axes=1) > math.cos(math.radians(1.0)), 1.0e6, 0.0)
    dI, m = ori.inertia_of_load(cap, lat, lon)
    exact = m * MOON_RADIUS**2 * (np.eye(3) - np.outer(n0, n0))
    exact -= np.eye(3) * np.trace(exact) / 3
    assert np.max(np.abs(dI - exact)) / np.max(np.abs(exact)) < 2e-3
    # a load's own geoid (degree 2) carries the same inertia
    _, Q = ori.degree2(cap, lat, lon)
    rho = 3 * ori.MASS / (4 * np.pi * MOON_RADIUS**3)
    N2 = 3 / (5 * rho) * np.einsum("ij,i...,j...->...", Q, units, units)
    assert np.allclose(ori.inertia_of_geoid(N2, lat, lon), dI, rtol=1e-6, atol=1e-6 * np.abs(dI).max())


def test_frames_round_trip():
    frame = ori.principal(np.array([[0, 0, -7e30], [0, 0, 0], [-7e30, 0, 0]], float))[1]
    lat, lon = np.array([45.0, -53.0, 0.0, 89.0]), np.array([-36.1, -169.0, 0.0, 10.0])
    tl, tn = ori.to_turned(lat, lon, frame)
    bl, bn = ori.to_today(tl, tn, frame)
    assert np.allclose(bl, lat, atol=1e-9) and np.allclose((bn - lon + 180) % 360 - 180, 0, atol=1e-9)
    assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-12) and np.linalg.det(frame) == pytest.approx(1.0)


def test_evolve_limits():
    """A constant load relaxing through one mode moves from the elastic state to the relaxed one."""
    load = np.array([[0, -1e30, -7e30], [-1e30, 0, 0], [-7e30, 0, 0]], float)
    rates = np.array([1e-3])
    times = np.r_[0.0, np.logspace(0, 6, 25)]
    rows = ori.evolve(times, [(lambda t: load, (-0.02, np.array([-0.5])))], (0.0, np.array([0.0])), rates, bulge=False)
    assert rows[0]["pole_deg"] == pytest.approx(ori.describe(0.98 * load)["pole_deg"], rel=1e-9)
    assert rows[-1]["pole_deg"] == pytest.approx(ori.describe(0.48 * load)["pole_deg"], rel=1e-6)
    # the bulge re-forming about the moved axes enlarges the long-term turn
    full = ori.evolve(times, [(lambda t: load, (-0.02, np.array([-0.5])))], (0.0, np.array([ori.fluid_love(ori.C_MR2)])), rates)
    assert full[-1]["pole_deg"] > rows[-1]["pole_deg"]


@pytest.mark.skipif(not PRODUCT.exists(), reason="results/moon_orientation.json not written")
def test_product_reproduces_the_reports_elastic_state():
    product = json.loads(PRODUCT.read_text())
    assert product["schema"] == "terluna.geology.moon-orientation/1"
    state = ori.describe(np.array(product["load_tensors_kg_m2"]["elastic_share"]))
    assert state["pole_deg"] == pytest.approx(8.05, abs=0.01)
    assert state["pole_toward_lon_deg"] == pytest.approx(-167.1, abs=0.1)
    assert state["sub_earth_lat_deg"] == pytest.approx(8.00, abs=0.01)
    assert state["sub_earth_lon_deg"] == pytest.approx(6.09, abs=0.01)
    frame = np.array(product["elastic_frame"]["rotation_today_to_turned"])
    assert np.allclose(frame, state["frame"], atol=1e-9)
    assert product["elastic_frame"]["north_pole"]["lat_deg"] == pytest.approx(81.95, abs=0.01)
    assert product["elastic_frame"]["sub_earth_point"]["lat_deg"] == pytest.approx(8.0, abs=0.05)
    lat, lon = ori.to_turned(45.0, -36.1, frame)
    assert (round(float(lat), 1), round(float(lon), 1)) == (39.4, -37.0)
    back = ori.to_today(lat, lon, frame)
    assert np.allclose([float(back[0]), float(back[1])], [45.0, -36.1], atol=1e-9)


@pytest.mark.skipif(not (PRODUCT.exists() and GRIDS.exists()), reason="crustal-loading grids or product absent")
def test_product_tensors_come_from_the_grids():
    product = json.loads(PRODUCT.read_text())
    grids = np.load(GRIDS)
    dI = ori.inertia_of_geoid(grids["elastic__share__N"].astype(float), grids["lat_deg"], grids["lon_deg"])
    assert np.allclose(dI, np.array(product["load_tensors_kg_m2"]["elastic_share"]), rtol=1e-9, atol=1e20)
