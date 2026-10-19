"""The propagation's controls: geometry, delays, spreading and area by independent calculation, the thunder model's
absorption in the target's air, and a field covariance that is a valid one and converged."""
import numpy as np
import pytest

import propagation
import scenes
from atmosphere.electricity import thunder

PRODUCT = scenes.load_arrivals()
RATE = propagation.SAMPLE_RATE
STILL_AIR = dict(sound_speed_m_s=343.0)


def render(position, area=1.0, spacing=0.2):
    impulse = np.zeros(RATE // 2)
    impulse[100] = 1.0
    return propagation.propagate_mono(impulse, np.asarray(position), STILL_AIR, spacing_m=spacing,
                                      listener_height_m=1.6, source_area_m2=area, absorption=False)


def test_absorption_is_the_thunder_models_in_the_targets_air():
    f = np.array([125.0, 1000.0, 4000.0, 10000.0])
    for target in scenes.TARGETS:
        for world in scenes.WORLDS:
            air = scenes.air(PRODUCT, target, world)
            expected = thunder.absorption_db_per_m(f, air['temperature_k'], air['pressure_pa'], air['vapour_percent'],
                                                   air['x_o2'], air['x_n2'])
            assert np.array_equal(propagation.air_absorption_db_per_m(f, air), expected)
        # The Moon's air holds less oxygen and more nitrogen than Earth's: its nitrogen relaxation absorbs more at
        # low frequencies, and its oxygen relaxation less at high ones.
        moon = scenes.air(PRODUCT, target, 'moon')
        ratio = (propagation.air_absorption_db_per_m([125.0, 4000.0, 8000.0], moon)
                 / propagation.air_absorption_db_per_m([125.0, 4000.0, 8000.0],
                                                       dict(moon, x_o2=thunder.EARTH_O2, x_n2=thunder.EARTH_N2)))
        assert 1.0 < ratio[0] < 1.05 and ratio[2] < ratio[1] < 0.97


def test_mirrored_and_central_sources():
    x, _ = render([3.0, 2.0, 0.0])
    mirrored, _ = render([-3.0, 2.0, 0.0])
    assert np.max(np.abs(x - mirrored[:, ::-1])) < 1e-10
    centred, _ = render([0.0, 3.0, 0.0])
    assert np.max(np.abs(centred[:, 0] - centred[:, 1])) < 1e-10


def test_paths_delays_and_gains():
    x, meta = render([3.0, 2.0, 0.0])
    receivers = np.array([[-0.1, 0.0, 1.6], [0.1, 0.0, 1.6]])
    distances = np.linalg.norm(np.array([3.0, 2.0, 0.0])[None, :] - receivers, axis=1)
    assert np.max(np.abs(np.array(meta['distances_m']) - distances)) < 1e-12
    lag = (np.argmax(np.abs(x[:, 0])) - np.argmax(np.abs(x[:, 1]))) / RATE
    assert abs(lag - (distances[0] - distances[1]) / 343.0) <= 1.1 / RATE       # the right microphone hears it first
    assert np.max(np.abs(x.sum(axis=0) - 1 / distances)) < 1e-8
    near, _ = render([0.0, 0.0, 3.6], spacing=0.0)
    far, _ = render([0.0, 0.0, 5.6], spacing=0.0)
    assert np.sum(far[:, 0]) / np.sum(near[:, 0]) == pytest.approx(0.5, abs=1e-8)
    quadruple, _ = render([3.0, 2.0, 0.0], area=4.0)
    assert np.max(np.abs(quadruple - 2 * x)) < 1e-10


def test_ring_cells_keep_the_annulus_power():
    edges, height = np.array([0.4, 2.0, 6.0, 16.0, 40.0]), 1.6
    patches = propagation.ring_patches(edges, sectors=8, listener_height_m=height)
    for ring, (inner, outer) in enumerate(zip(edges[:-1], edges[1:])):
        modeled = sum(p['source_area_m2'] / np.sum((np.array(p['source_xyz_m']) - [0, 0, height])**2)
                      for p in patches if p['ring'] == ring)
        assert modeled == pytest.approx(np.pi * np.log((outer**2 + height**2) / (inner**2 + height**2)), rel=1e-12)


def test_the_field_covariance_is_valid_and_converged():
    air = scenes.air(PRODUCT, 'resort', 'moon')
    f = np.r_[0.0, np.geomspace(20, 20000, 160)]
    field = propagation.field_covariance(f, air)
    c = propagation.covariance_matrix(field)
    assert np.max(np.abs(c - c.conj().transpose(0, 2, 1))) < 1e-12 * np.max(np.abs(c))
    assert np.min(np.linalg.eigvalsh(c)) >= -1e-12 * np.max(np.abs(c))
    h = propagation.covariance_factor(field)
    assert np.max(np.abs(h @ h.conj().transpose(0, 2, 1) - c)) < 1e-10 * np.max(np.abs(c))
    assert np.max(np.abs(field['power_left'] / field['power_right'] - 1)) < 1e-12        # a symmetric sheet
    coincident = propagation.field_covariance(f, air, spacing_m=0.0)
    assert np.max(np.abs(coincident['cross_rl'] / coincident['power_left'] - 1)) < 1e-12
    expected = np.pi * np.log((40**2 + 1.6**2) / (0.2**2 + 1.6**2))
    assert coincident['power_left'][0] == pytest.approx(expected, rel=1e-12)
    dense = propagation.field_covariance(f, air, azimuth_count=256)
    assert np.max(np.abs(field['cross_rl'] / field['power_left'] - dense['cross_rl'] / dense['power_left'])) < 1e-6


def test_fractional_delays():
    x = np.zeros(1000)
    x[100] = 1.0
    assert np.array_equal(propagation.fractional_delay(x, 25.0)[125:126], [1.0])
    half = propagation.fractional_delay(x, 25.5)
    assert np.argmax(half) in (125, 126) and half[125] == pytest.approx(half[126], rel=2e-3)
