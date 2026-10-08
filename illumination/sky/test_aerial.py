"""The air between an observer and a surface, against the solved sky it comes from."""
import json
from pathlib import Path

import numpy as np
import pytest

from illumination.sky import aerial
from illumination.sky.solved_transport import evaluate

SOLUTION = Path(__file__).resolve().parents[2] / "research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz"
pytestmark = pytest.mark.skipif(not SOLUTION.exists(), reason="the cached Moon sky solution is not restored")


@pytest.fixture(scope="module")
def solved():
    with np.load(SOLUTION, allow_pickle=False) as z:
        solution = {k: z[k] for k in z.files if k != "meta"}
        meta = json.loads(str(z["meta"]))
    return solution, meta


def test_whole_rays_are_the_solved_sky(solved):
    solution, meta = solved
    elevations, azimuths = np.array([0.5, 4.0, 20.0, 60.0]), np.array([0.0, 50.0, 180.0])
    for sun in (30.0, -15.0):
        inscatter, transmittance, length = aerial.view(solution, sun, elevations, azimuths, [1e9])
        reference = evaluate(solution, meta, np.array([sun]), elevations, azimuths)[0] @ solution["xyz"]
        np.testing.assert_allclose(inscatter[:, :, 0], reference, rtol=1e-9)


def test_partial_rays_grow_toward_the_sky_and_transmit_less(solved):
    solution, _ = solved
    distances = np.array([0.0, 100.0, 1e3, 1e4, 1e5, 1e9])
    inscatter, transmittance, length = aerial.view(solution, 20.0, np.array([0.2, -0.05]), np.array([90.0]),
                                                   distances, height_m=2.0)
    y = inscatter[:, 0, :, 1]
    assert np.all(y[:, 0] == 0) and np.all(np.diff(y, axis=1) >= 0)
    assert np.all(transmittance[:, 0] == 1) and np.all(np.diff(transmittance, axis=1) <= 1e-15)
    # Short paths near the ground: in-scatter grows in proportion to distance, transmittance falls as Beer-Lambert.
    chi = solution["scattering"][0] + solution["absorption"][0]
    np.testing.assert_allclose(transmittance[0, 1], np.exp(-chi * 100.0), rtol=1e-6)
    assert y[0, 2] / y[0, 1] == pytest.approx(10.0, rel=0.02)


def test_rays_from_a_height_dip_below_the_level_to_the_horizon(solved):
    solution, _ = solved
    dip = np.degrees(np.sqrt(2 * 2.0 / solution["edges"][0]))
    _, _, length = aerial.view(solution, 20.0, np.array([-0.5 * dip, -1.5 * dip]), np.array([0.0]), [1e9],
                               height_m=2.0)
    assert length[0] > 1e5                                  # skims past the horizon into the sky
    assert length[1] == pytest.approx(2.0 / np.sin(np.radians(1.5 * dip)), rel=0.2)   # meets the sea
