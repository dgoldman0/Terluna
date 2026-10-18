"""Digital levels and audio files for the auditions: loudness by ITU-R BS.1770-4, peaks, fades and WAV writing.

Loudness is the gated integrated loudness of BS.1770-4 (LUFS): each channel K-weighted by the standard's two biquads
at 48 kHz, mean squares over 400 ms blocks overlapping by 75%, an absolute gate at -70 LUFS and a relative gate 10 LU
below, channels summed with weight one. The rain work matched its comparison passages to -35.4 LUFS with ffmpeg's
meter of the same standard; this module measures them without ffmpeg. All levels are digital: they say nothing about
sound pressure in the world.
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
from scipy import signal

SAMPLE_RATE = 48000
TARGET_LUFS = -35.4          # the rain work's listening level for comparison passages
GAP_SECONDS = 0.8
FADE_SECONDS = 0.04
# K-weighting at 48 kHz (BS.1770-4, tables 1 and 2): a high-frequency shelf, then the RLB high-pass.
SHELF = ([1.53512485958697, -2.69169618940638, 1.19839281085285], [1.0, -1.69065929318241, 0.73248077421585])
HIGHPASS = ([1.0, -2.0, 1.0], [1.0, -1.99004745483398, 0.99007225036621])


def _channels(x):
    x = np.asarray(x, dtype=np.float64)
    return x[:, None] if x.ndim == 1 else x


def loudness(x, fs=SAMPLE_RATE):
    """Integrated loudness (LUFS) of mono or channels-last audio; None if every block is below the absolute gate."""
    if fs != SAMPLE_RATE:
        raise ValueError('The K-weighting coefficients are those for 48 kHz')
    x = _channels(x)
    weighted = signal.lfilter(*HIGHPASS, signal.lfilter(*SHELF, x, axis=0), axis=0)
    block, step = round(0.4 * fs), round(0.1 * fs)
    if len(weighted) < block:
        raise ValueError('Loudness needs at least 400 ms of audio')
    starts = np.arange(0, len(weighted) - block + 1, step)
    power = np.cumsum(np.r_[np.zeros((1, x.shape[1])), weighted**2], axis=0)
    z = (power[starts + block] - power[starts]) / block           # (blocks, channels) mean squares
    blocks = -0.691 + 10 * np.log10(np.maximum(z.sum(axis=1), 1e-300))
    gated = blocks > -70.0
    if not np.any(gated):
        return None
    relative = -0.691 + 10 * np.log10(z[gated].sum(axis=1).mean()) - 10.0
    keep = gated & (blocks > relative)
    return float(-0.691 + 10 * np.log10(z[keep].sum(axis=1).mean()))


def rms_dbfs(x):
    value = float(np.sqrt(np.mean(np.asarray(x, dtype=np.float64)**2)))
    return float(20 * np.log10(value)) if value > 0 else None


def true_peak(x, chunk=SAMPLE_RATE * 5, overlap=256):
    """The peak after fourfold oversampling, as a fraction of full scale (BS.1770-4, annex 2, in outline), taken in
    overlapping five-second pieces to bound memory."""
    x = _channels(x)
    peak = 0.0
    for start in range(0, len(x), chunk):
        piece = x[max(0, start - overlap):start + chunk + overlap]
        peak = max(peak, float(np.max(np.abs(signal.resample_poly(piece, 4, 1, axis=0)))))
    return peak


def fade(x, seconds=FADE_SECONDS, fs=SAMPLE_RATE):
    """Sine-squared fades at both ends."""
    x = np.array(x, dtype=np.float64, copy=True)
    n = min(round(seconds * fs), len(x) // 2)
    if n:
        w = np.sin(np.linspace(0., np.pi / 2., n))**2
        shape = (n,) + (1,) * (x.ndim - 1)
        x[:n] *= w.reshape(shape)
        x[-n:] *= w[::-1].reshape(shape)
    return x


def gain_to(x, target_lufs=TARGET_LUFS):
    """The scalar gain (dB) that brings x to the target loudness. A quiet clip is measured on a prescaled copy, so
    the absolute gate never decides it."""
    prescale_db = -30.0 - rms_dbfs(x)
    measured = loudness(np.asarray(x) * 10**(prescale_db / 20))
    return prescale_db + target_lufs - measured


def sequence(passages, gap_seconds=GAP_SECONDS, fs=SAMPLE_RATE):
    """Passages (each already faded and scaled) joined with silent gaps, and each passage's start and end in s."""
    channels = max(_channels(p).shape[1] for p in passages)
    gap = np.zeros((round(gap_seconds * fs), channels))
    parts, spans, t = [], [], 0.0
    for i, p in enumerate(passages):
        p = _channels(p)
        if p.shape[1] != channels:
            p = np.repeat(p, channels, axis=1)
        if i:
            parts.append(gap)
            t += len(gap) / fs
        parts.append(p)
        spans.append([t, t + len(p) / fs])
        t += len(p) / fs
    out = np.concatenate(parts)
    return (out[:, 0] if channels == 1 else out), spans


def write_wav(path, x, fs=SAMPLE_RATE, bits=24):
    """Write PCM WAV of 16 or 24 bits, or 32-bit float with bits=32. Samples beyond full scale stop the write rather
    than clip."""
    x = _channels(x)
    if not np.all(np.isfinite(x)):
        raise ValueError('Non-finite audio')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    channels = x.shape[1]
    if bits == 32:
        data, fmt = x.astype('<f4').tobytes(), 3
    else:
        if np.max(np.abs(x)) > 1.0:
            raise ValueError(f'{path.name}: samples exceed full scale; lower the gain')
        full = 2**(bits - 1) - 1
        q = np.round(x * full).astype('<i4')
        data = q.astype('<i2').tobytes() if bits == 16 else q.view(np.uint8).reshape(-1, 4)[:, :3].tobytes()
        fmt = 1
    width = bits // 8
    header = b'RIFF' + struct.pack('<I', 36 + len(data)) + b'WAVE' + b'fmt ' + struct.pack(
        '<IHHIIHH', 16, fmt, channels, fs, fs * channels * width, channels * width, bits) + b'data' + struct.pack(
        '<I', len(data))
    path.write_bytes(header + data)
    return path


def read_wav(path):
    """A WAV written here, or any PCM or float WAV scipy reads, as float64 channels-last and its rate."""
    from scipy.io import wavfile
    rate, raw = wavfile.read(path)
    if np.issubdtype(raw.dtype, np.integer):
        x = raw.astype(np.float64) / float(-np.iinfo(raw.dtype).min)
    else:
        x = raw.astype(np.float64)
    return rate, x


def describe(x, fs=SAMPLE_RATE) -> dict:
    x = _channels(x)
    return dict(seconds=len(x) / fs, channels=int(x.shape[1]), loudness_lufs=loudness(x, fs), rms_dbfs=rms_dbfs(x),
                sample_peak=float(np.max(np.abs(x))), true_peak=true_peak(x))
