"""Thunder heard at the ground, rendered from the thunder study's product.

The atmospheric-electricity study traces each flash's sound to the ground (research/studies/atmospheric_electricity,
results/thunder.json, from atmosphere/electricity/thunder.py) and gives, for a listener below the flash, the delay of
the first sound, the span holding 90% of the audible energy, and each octave band's loudest one-second level in dB re
20 uPa. This layer renders one such event: in each octave band whose level reaches the threshold of hearing (ISO 226,
as the thunder model carries it), and in the 16 Hz band of infrasound, band-limited noise follows a common envelope and
is scaled so that its loudest second has the product's level. The envelope's fine shape, a rise, a decay, four
arrivals of constructive wavefronts and a slow irregular swell, stretched to the product's span, and the noise itself
are synthesized: the product gives no channel waveform. The result is pressure in pascals at the listener, ground
reflection included as the product includes it.

The layer is the storm listening study's of 11 October 2026 (render_storm.py, make_thunder), brought into the
repository on 17 October. The product it read then is the one the repository holds now (the same sha256). It renders
the study's flash in the box's own air, which is the old design climate's: no target climate's thunder is computed.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import fft, signal

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from atmosphere.electricity import thunder  # noqa: E402

PRODUCT = ROOT / 'research' / 'studies' / 'atmospheric_electricity' / 'results' / 'thunder.json'
SCHEMA = 'terluna.research.atmospheric-electricity-thunder/1'
SAMPLE_RATE = 48000
P_REF = thunder.P_REF                 # 20 uPa
# The events this lab renders: (label, path into the product). Each is a listener below the flash.
EVENTS = {
    'moon-cloud-flash': ('lunar', 'typical_eff0.0018', '45'),
    'moon-ground-strike': ('lunar', 'ground_typical_eff0.0018', '45'),
    'earth-cloud-flash': ('earth_benchmark_flash', '0.0018'),
    'earth-ground-strike': ('earth_benchmark_ground_strike', '0.0018'),
}


def load(path=PRODUCT) -> dict:
    product = json.loads(Path(path).read_text())
    if product.get('schema') != SCHEMA:
        raise ValueError(f'Expected {SCHEMA}, found {product.get("schema")}')
    return product


def sha256(path=PRODUCT) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def event(product, name) -> dict:
    """The product's entry below the flash for one of EVENTS: first_s, span_s and peak_db_by_band."""
    path = EVENTS[name]
    if path[0] == 'lunar':
        return product['lunar'][path[1]][path[2]]['below']
    return product[path[0]]['results'][path[1]]['below']


def rendered_bands(levels: dict) -> list:
    """The bands that are rendered: the 16 Hz band, and each audible band at or above the threshold of hearing."""
    out = []
    for band, threshold in zip(thunder.BANDS, thunder.THRESHOLD_DB):
        level = levels.get(f'{band:.1f}')
        if level is None:
            continue
        if math.isnan(threshold) or level >= threshold:
            out.append((float(band), float(level)))
    return out


def band_noise(rng, n, lo, hi):
    """Unit-variance noise confined to lo-hi Hz by a tapered spectral window."""
    count = fft.next_fast_len(n)
    spectrum = fft.rfft(rng.standard_normal(count))
    freq = fft.rfftfreq(count, 1 / SAMPLE_RATE)
    weight = np.zeros_like(freq)
    inside = (freq >= lo) & (freq <= hi)
    weight[inside] = signal.windows.tukey(int(inside.sum()), .24)
    x = fft.irfft(spectrum * weight, count)[:n]
    return x / max(float(np.std(x)), 1e-15)


def smooth_random(rng, duration, dt, sigma):
    """A smooth random curve: Gaussian knots every dt seconds, lightly smoothed and joined linearly."""
    knots = np.arange(-dt, duration + 2 * dt, dt)
    values = np.convolve(rng.normal(0, sigma, len(knots)), [.25, .5, .25], mode='same')
    return np.interp(np.arange(int(round(duration * SAMPLE_RATE))) / SAMPLE_RATE, knots, values).astype(np.float32)


def loudest_second(x):
    """The root-mean-square pressure of the loudest one-second window, in quarter-second steps."""
    ss = np.r_[0, np.cumsum(np.asarray(x, dtype=np.float64)**2)]
    starts = np.arange(0, max(1, len(x) - SAMPLE_RATE + 1), SAMPLE_RATE // 4)
    return math.sqrt(float(np.max((ss[starts + SAMPLE_RATE] - ss[starts]) / SAMPLE_RATE)))


def envelope(rng, span_s):
    """The event's envelope, stretched so that 90% of its energy falls within span_s."""
    seconds = 44.0
    tn = np.arange(round(seconds * SAMPLE_RATE)) / SAMPLE_RATE
    env = ((1 - np.exp(-tn / .70)) * np.exp(-np.maximum(tn - 4, 0) / 13.2)
           * np.exp(smooth_random(rng, seconds, 1.25, .34)))
    env *= 1 / (1 + np.exp((tn - 39) * 1.3))
    ramp = int(.03 * SAMPLE_RATE)
    env[:ramp] *= np.linspace(0, 1, ramp)
    for mu, amp, width in ((3.8, .6, 1.0), (8.7, .75, 1.8), (16.2, .9, 2.4), (25.5, .8, 2.8)):
        env += amp * np.exp(-.5 * ((tn - mu) / width)**2) * (1 - np.exp(-tn / .6))
    cumulative = np.cumsum(env**2)
    q5, q95 = np.interp(np.array([.05, .95]) * cumulative[-1], cumulative, tn)
    stretch = span_s / (q95 - q5)
    t2 = np.arange(round(seconds * stretch * SAMPLE_RATE)) / SAMPLE_RATE
    return np.interp(t2 / stretch, tn, env, right=0)


def render(product, name, seed=20261011):
    """One event's pressure (Pa) at the listener, from its first sound, and a record of its bands' levels."""
    levels = event(product, name)
    rng = np.random.default_rng(seed)
    env = envelope(rng, levels['span_s'])
    pressure = np.zeros(len(env))
    bands = {}
    for f, level in rendered_bands(levels['peak_db_by_band']):
        x = band_noise(rng, len(env), max(9.0, f / math.sqrt(2)), f * math.sqrt(2))
        x *= env * np.exp(smooth_random(rng, len(env) / SAMPLE_RATE, .55, .12))
        x *= P_REF * 10**(level / 20) / max(loudest_second(x), 1e-12)
        pressure += x
        bands[f'{f:g}'] = dict(product_db=level, rendered_db=20 * math.log10(loudest_second(x) / P_REF))
    return pressure, dict(event=name, path=list(EVENTS[name]), first_sound_s=levels['first_s'],
                          span_s=levels['span_s'], peak_dba=levels.get('peak_dba'), bands=bands, seed=seed,
                          synthesized=['envelope shape within the span', 'band noise', 'slow swells'])


def stereo(pressure, azimuth_deg, sound_speed_m_s, spacing_m=0.08):
    """A distant event at an azimuth (clockwise from straight ahead) on the lab's spaced pair: a plane wave whose
    only difference between channels is its delay, spacing sin(azimuth) / c."""
    from propagation import fractional_delay
    lag = spacing_m * math.sin(math.radians(azimuth_deg)) / sound_speed_m_s * SAMPLE_RATE
    pad = int(math.ceil(abs(lag))) + 1
    x = np.r_[pressure, np.zeros(pad)]
    left = fractional_delay(x, max(lag, 0.0))
    right = fractional_delay(x, max(-lag, 0.0))
    return np.column_stack((left, right))
