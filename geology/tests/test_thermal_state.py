import json
import math
from pathlib import Path

import numpy as np
import pytest

from geology import thermal_state as ts
from shared.constants import MOON_RADIUS
from shared.provenance import constants_changed

RESULT = Path(ts.__file__).parent / "results" / "interior_response.json"


def test_profiles_pass_through_their_nodes_and_warm_with_depth():
    for name, spec in ts.SPEC["temperature_profiles"].items():
        if name == "about":
            continue
        z = np.array(spec["depth_km"]) * 1e3
        assert np.allclose(ts.temperature_k(z, name), np.array(spec["temperature_c"]) + 273.15)
        deep = np.linspace(50e3, 1237e3, 200)
        assert np.all(np.diff(ts.temperature_k(deep, name)) > 0)
    mid = np.array([100e3, 200e3, 400e3])                   # the central profile joins the warm one below 700 km
    assert np.all(ts.temperature_k(mid, "cool") < ts.temperature_k(mid, "central"))
    assert np.all(ts.temperature_k(mid, "central") < ts.temperature_k(mid, "warm"))


def test_flow_law_units():
    """Dry diffusion creep at 1,600 K, no pressure, 1 mm grains: sigma / (3 strain rate) by hand."""
    eta = ts.creep_viscosity(np.array([1600.0]), np.array([0.0]), 0.001, 1.0, False)["diffusion"][0]
    rate_per_mpa = 1.5e9 * (1000.0) ** -3 * math.exp(-375e3 / (ts.GAS_CONSTANT * 1600.0))
    assert eta == pytest.approx(1e6 / (3 * rate_per_mpa), rel=1e-9)
    wet = ts.creep_viscosity(np.array([1600.0]), np.array([0.0]), 0.001, 1.0, True)["diffusion"][0]
    assert wet < eta / 10                                   # 100 ppm water weakens diffusion creep tens of times


def test_pressure_and_viscosity_profile():
    assert ts.pressure_pa(1237.4e3) == pytest.approx(4.4e9, rel=0.05)
    z = np.linspace(60e3, 1230e3, 100)
    eta = ts.viscosity_pa_s(z)
    assert np.all(np.diff(np.log10(eta)) < 0)
    assert math.isinf(ts.viscosity_pa_s(np.array([20e3]))[0])
    assert ts.viscosity_pa_s(np.array([1300e3]))[0] == ts.LVZ_VISCOSITY_PA_S[0]
    assert ts.lid_depth_m(1e7, ts.STIFF) > ts.lid_depth_m(1e7) > ts.lid_depth_m(1e7, ts.SOFT)


def test_shells_cover_the_moon():
    stack = ts.shells()
    tops = [s[0] for s in stack]
    assert tops[-1] == MOON_RADIUS and np.all(np.diff(tops) > 0)
    assert stack[0][2] == 0.0 and stack[1][3] == ts.LVZ_VISCOSITY_PA_S[0] and math.isinf(stack[-1][3])


@pytest.mark.skipif(not RESULT.exists(), reason="results/interior_response.json not written")
def test_product_header_and_constants():
    product = json.loads(RESULT.read_text())
    assert product["schema"] == "terluna.geology.viscoelastic-response/1"
    assert not constants_changed(product["producer"]["constants"])
    assert product["response"]["central"]["spectrum_misfit"] < 1e-2
    k2 = product["response"]["central"]["by_degree"]["2"]["k_load"]["at_years"]
    values = [k2[y] for y in ("1", "1000", "1e+06", "1e+09")]
    assert all(b < a for a, b in zip(values, values[1:]))       # the load sinks in as time passes
