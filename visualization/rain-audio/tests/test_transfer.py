"""The impact-force transfer scales as its pulses say it must, keeps Earth unchanged, and filters causally."""
import math

import numpy as np
import pytest

import scenes
import transfer

PRODUCT = scenes.load_arrivals()


def drops(diameter_mm, flux, speed, label='test'):
    return scenes.Population(label, np.array([diameter_mm * 1e-3]), np.array([float(flux)]), np.array([float(speed)]))


def low_frequency_limit(p):
    mass = p.liquid_density_kg_m3 * math.pi / 6 * p.diameter_m**3
    return 2 * float(np.sum(p.number_flux_m2_s * (mass * p.normal_speed_m_s)**2))


def test_the_spectrum_starts_at_twice_the_impulse_power():
    p = scenes.direct_rain(PRODUCT, 'resort', 'moon', 5.0)
    assert transfer.force_psd(p, [0.0])[0] == pytest.approx(low_frequency_limit(p), rel=1e-12)


def test_speed_scalings_at_low_and_high_frequency():
    slow, fast = drops(2.0, 100, 2.0), drops(2.0, 100, 4.0)
    f = np.array([0.0, 1e8])
    ratio = transfer.force_psd(fast, f) / transfer.force_psd(slow, f)
    assert ratio[0] == pytest.approx(4.0, rel=1e-12)          # J^2 grows as U^2
    assert ratio[1] == pytest.approx(32.0, rel=1e-6)          # and the pulse's tail as (U / D)^3


def test_fewer_larger_drops_at_fixed_rain_force_more_at_low_frequency():
    rain = 5.0 / scenes.MM_H_PER_M_S
    small = drops(1.0, rain / (math.pi / 6 * 1e-9), 3.0)
    large = drops(2.0, rain / (math.pi / 6 * 8e-9), 3.0)
    assert small.rain_mm_h == pytest.approx(large.rain_mm_h)
    assert transfer.force_psd(large, [0.0])[0] / transfer.force_psd(small, [0.0])[0] == pytest.approx(8.0)


def test_flux_and_pulse_width_scalings():
    p = drops(2.0, 100, 3.0)
    f = np.array([0.0, 2e3, 1e8])
    assert np.allclose(transfer.force_psd(drops(2.0, 400, 3.0), f), 4 * transfer.force_psd(p, f))
    wide = transfer.force_psd(p, f, transfer.PUDDLE_WIDTH_SCALE) / transfer.force_psd(p, f)
    assert wide[0] == pytest.approx(1.0) and wide[-1] == pytest.approx(transfer.PUDDLE_WIDTH_SCALE**-3, rel=1e-6)


def test_identical_scenes_leave_the_texture_untouched():
    earth = scenes.scene(PRODUCT, 'lock', 'earth', 'concrete')
    magnitude, meta = transfer.make_transfer(earth, earth)
    assert meta['exact_identity'] and np.all(magnitude == 1.0)
    audio = np.random.default_rng(0).standard_normal(4800)
    out, detail = transfer.apply_transfer(audio, transfer.FREQUENCIES, magnitude)
    assert detail['identity'] and np.array_equal(out, audio)


@pytest.mark.parametrize('surface', ('concrete', 'puddles', 'grass-leaves'))
def test_the_moons_rain_is_quieter_and_duller_in_every_target(surface):
    for target in scenes.TARGETS:
        earth = scenes.scene(PRODUCT, target, 'earth', surface)
        moon = scenes.scene(PRODUCT, target, 'moon', surface)
        magnitude, meta = transfer.make_transfer(earth, moon)
        assert np.all(magnitude < 1.0)
        assert meta['gain_db'][1] - meta['gain_db'][5] > 8          # 100 Hz against 10 kHz
        assert np.all(np.diff(magnitude[transfer.FREQUENCIES > 50]) <= 1e-12)


def test_the_filter_is_causal_and_keeps_the_magnitude():
    earth = scenes.scene(PRODUCT, 'resort', 'earth', 'grass-leaves')
    moon = scenes.scene(PRODUCT, 'resort', 'moon', 'grass-leaves')
    magnitude, _ = transfer.make_transfer(earth, moon)
    kernel, detail = transfer._minimum_phase_kernel(transfer.FREQUENCIES, magnitude, transfer.SAMPLE_RATE)
    assert detail['magnitude_error_db_rms'] < 0.01 and detail['discarded_tail_energy_fraction'] <= 1e-8
    impulse = np.zeros(20000)
    impulse[5000] = 1.0
    out, _ = transfer.apply_transfer(impulse, transfer.FREQUENCIES, magnitude)
    assert np.max(np.abs(out[:5000])) < 1e-12 * np.max(np.abs(out)) and len(out) == len(impulse)


def test_large_drops_arrive_only_within_their_events():
    earth = scenes.scene(PRODUCT, 'resort', 'earth', 'concrete')
    heavy = scenes.scene(PRODUCT, 'resort', 'moon', 'concrete', tail_fraction=0.05)
    bed = np.random.default_rng(1).standard_normal(48000 * 20) * 1e-3
    stem, meta = transfer.large_drop_events(bed, earth, heavy, 20.0, 99, area_m2=20.0)
    assert meta['count'] > 0 and meta['expected_count'] == pytest.approx(
        heavy.tail.arrivals_m2_s * 20.0 * 20.0)
    support = np.zeros(len(stem), dtype=bool)
    for e in meta['events']:
        support[e['sample_start']:e['sample_end']] = True
    assert np.all(stem[~support] == 0) and np.any(stem[support] != 0)
