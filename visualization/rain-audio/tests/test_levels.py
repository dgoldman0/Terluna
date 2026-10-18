"""Loudness follows ITU-R BS.1770-4's reference readings, and WAV files round-trip."""
import numpy as np
import pytest

import levels

RATE = levels.SAMPLE_RATE


def sine(f, seconds, amplitude):
    t = np.arange(round(seconds * RATE)) / RATE
    return amplitude * np.sin(2 * np.pi * f * t)


def test_reference_readings():
    # A full-scale sine near 1 kHz in one channel reads -3.01 LKFS (BS.1770-4, section 2.1).
    assert levels.loudness(sine(997.0, 10, 1.0)) == pytest.approx(-3.01, abs=0.02)
    # EBU Tech 3341 case 1: a stereo 1 kHz sine at -23 dBFS reads -23.0 LUFS.
    x = sine(1000.0, 20, 10**(-23 / 20))
    assert levels.loudness(np.column_stack((x, x))) == pytest.approx(-23.0, abs=0.1)
    # The gates leave out silence and near-silence.
    quiet = np.r_[x, sine(1000.0, 20, 1e-5)]
    assert levels.loudness(np.column_stack((quiet, quiet))) == pytest.approx(-23.0, abs=0.1)


def test_gain_to_a_target():
    x = sine(440.0, 5, 0.01) + 1e-3 * np.random.default_rng(0).standard_normal(5 * RATE)
    gain = levels.gain_to(x, -35.4)
    assert levels.loudness(x * 10**(gain / 20)) == pytest.approx(-35.4, abs=1e-6)


def test_wav_round_trip_and_no_clipping(tmp_path):
    x = np.column_stack((sine(300.0, 1, 0.5), sine(500.0, 1, -0.25)))
    for bits, tolerance in ((24, 1 / 2**23), (16, 1 / 2**15), (32, 1e-7)):
        path = levels.write_wav(tmp_path / f'x{bits}.wav', x, bits=bits)
        rate, y = levels.read_wav(path)
        assert rate == RATE and y.shape == x.shape and np.max(np.abs(y - x)) <= tolerance
    with pytest.raises(ValueError):
        levels.write_wav(tmp_path / 'loud.wav', 3.0 * x)
