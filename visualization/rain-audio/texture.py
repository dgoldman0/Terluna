"""Earth rain textures: fit compact statistics from a recording, and render new rain from those statistics alone.

A texture model describes how a recording's sound fluctuates in 64 overlapping frequency bands spaced evenly on the
ERB-number scale (Glasberg and Moore 1990), as amplitude envelopes measured every 0.25 ms: each band's 257 marginal
quantiles, and the smoothed complex cross-spectra of the bands' modulations, which keep how the bands rise and fall
together and drop the order in which they did. Rendering draws fresh Gaussian innovations with that cross-spectral
structure, maps each band through its quantiles, corrects the overall intensity distribution with a monotone quantile
map calibrated on a 120-second pilot render (version 3 of the renderer), and imposes the envelopes on fresh noise
carriers. Rendering reads no recording, envelope trace or recorded event.

The method is the rain work's of 11 October 2026, brought into the repository on 17 October. With the same model, seed
and duration the port reproduces the archive's raw masters bit for bit (tests/test_texture.py). The three models are
statistics of CC0 recordings (models/models.json): textures, with no mapping to drop size, impact speed or rain rate.
The audio sample rate is 48 kHz throughout.

    import texture
    model = texture.load('concrete')
    audio = texture.synthesize(model, seconds=10, seed=1001)
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import warnings
from pathlib import Path

import numpy as np
from scipy import fft, linalg, special
from scipy.io import wavfile

HERE = Path(__file__).resolve().parent
MODELS = HERE / 'models'
SURFACES = ('concrete', 'puddles', 'grass-leaves')
SAMPLE_RATE = 48000
GUARD_SECONDS = 2.0          # audio rendered beyond each end of a clip and cropped, so filters start settled
LATENT_KNOTS = 257           # quantile knots per band, evenly spaced in the Gaussian latent from -4.25 to 4.25
PILOT = dict(seconds=120.0, seed=17011)   # the pilot render behind the intensity correction
MAX_SECONDS = 300.0


def load(surface_or_path) -> dict:
    """A texture model as a dict of arrays, from a surface's name or an NPZ path."""
    path = Path(surface_or_path)
    if not path.suffix:
        if str(surface_or_path) not in SURFACES:
            raise ValueError(f'surface must be one of {SURFACES}')
        path = MODELS / f'{surface_or_path}.npz'
    with np.load(path, allow_pickle=False) as data:
        return dict(data)


def model_path(surface: str) -> Path:
    return MODELS / f'{surface}.npz'


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _scalar(model, name):
    return np.asarray(model[name]).item()


# --- bands, envelopes and carriers ------------------------------------------------------------------------------

def band_weights(n, bands):
    """The bands' amplitude weights on an n-sample real FFT grid: half-cosines one band-spacing wide on each side of
    centres spaced evenly in ERB number from 0 Hz to the Nyquist frequency."""
    if not isinstance(bands, (int, np.integer)) or bands < 2:
        raise ValueError('Use at least two audio bands.')
    f = fft.rfftfreq(n, 1 / SAMPLE_RATE)
    u = 21.4 * np.log10(1 + .00437 * f)
    grid = np.linspace(0, 21.4 * np.log10(1 + .00437 * SAMPLE_RATE / 2), bands)
    for center in grid:
        d = (u - center) / (grid[1] - grid[0])
        yield np.where(np.abs(d) < 1, np.cos(d * np.pi / 2), 0).astype(np.float32)


def analytic(spectrum, n):
    """The analytic signal of a one-sided spectrum on an n-sample grid."""
    z = np.zeros(n, dtype=np.complex64)
    z[:len(spectrum)] = spectrum
    z[1:len(spectrum) - (1 if n % 2 == 0 else 0)] *= 2
    return fft.ifft(z, workers=2)


def extract_envelopes(wave, bands=64, envelope_rate=4000):
    """Each band's root-mean-square amplitude over frames of 1 / envelope_rate s, shape (bands, frames)."""
    if not isinstance(envelope_rate, (int, np.integer)) or envelope_rate <= 0 or SAMPLE_RATE % envelope_rate:
        raise ValueError('Envelope rate must be a positive integer divisor of 48000.')
    if not isinstance(bands, (int, np.integer)) or bands < 2:
        raise ValueError('Use at least two audio bands.')
    hop = SAMPLE_RATE // envelope_rate
    frames = len(wave) // hop
    n = fft.next_fast_len(len(wave))
    spectrum = fft.rfft(wave, n=n, workers=2)
    result = np.empty((bands, frames), dtype=np.float32)
    for b, h in enumerate(band_weights(n, bands)):
        a = analytic(spectrum * h, n)[:frames * hop]
        result[b] = np.sqrt(np.mean((np.abs(a)**2).reshape(frames, hop), axis=1))
    return result


def render_envelopes(env, envelope_rate, seed):
    """Audio whose bands follow the given envelopes on fresh noise carriers, one independent carrier per band."""
    bands, frames = env.shape
    duration = frames / float(envelope_rate)
    samples = round(duration * SAMPLE_RATE)
    n = fft.next_fast_len(samples)
    out = np.zeros(n // 2 + 1, dtype=np.complex64)
    partition = np.zeros(n // 2 + 1, dtype=np.float32)
    rng = np.random.default_rng(seed)
    # A measured envelope belongs to the centre of its integration frame.
    frame_time = (np.arange(frames) + .5) / envelope_rate
    time_axis = np.arange(n) / SAMPLE_RATE
    for b, h in enumerate(band_weights(n, bands)):
        noise = rng.standard_normal(n, dtype=np.float32)
        carrier = analytic(fft.rfft(noise, workers=2) * h, n)
        a = np.abs(carrier)
        floor = max(float(np.sqrt(np.mean(a * a))) * .05, 1e-12)
        unit = carrier / np.maximum(a, floor)
        modulation = np.interp(time_axis, frame_time, env[b]).astype(np.float32)
        out += fft.rfft(np.real(unit) * modulation, workers=2) * h
        partition += h**4
    out /= np.sqrt(np.maximum(partition, 1e-8))
    return fft.irfft(out, n, workers=2)[:samples]


# --- the statistical model ----------------------------------------------------------------------------------------

def _spectral_factors(z, rate, frequency_bins=48, shrinkage=.04):
    """Smoothed one-sided cross-spectral density of the bands' latent modulations, as Cholesky factors.

    Welch windows average over the recording; logarithmic groups average neighbouring modulation frequencies;
    shrinkage towards the diagonal keeps every estimate positive definite. The factors are normalized so that each
    band has unit variance after synthesis, before its marginal map.
    """
    bands, samples = z.shape
    nperseg = min(16384, 2**int(np.floor(np.log2(samples))))
    if nperseg < 512:
        raise ValueError('At least 512 envelope samples are required.')
    step = nperseg // 2
    starts = np.arange(0, samples - nperseg + 1, step, dtype=int)
    if starts[-1] != samples - nperseg:
        starts = np.append(starts, samples - nperseg)
    window = np.hanning(nperseg).astype(np.float32)
    norm = np.sqrt(2.0 / (rate * float(np.sum(window**2))))
    bins = nperseg // 2 + 1
    edges = np.unique(np.rint(np.geomspace(1, bins, frequency_bins + 1)).astype(int))
    groups = list(zip(edges[:-1], edges[1:]))
    covariances = np.zeros((len(groups), bands, bands), dtype=np.complex128)
    for start in starts:
        spectrum = fft.rfft(z[:, start:start + nperseg] * window, axis=1, workers=2) * norm
        # DC is left out; the Nyquist bin is real and has no doubled one-sided term.
        spectrum[:, -1] /= np.sqrt(2.0)
        for i, (lo, hi) in enumerate(groups):
            part = spectrum[:, lo:hi]
            covariances[i] += (part @ part.conj().T) / (hi - lo)
    covariances /= len(starts)
    frequencies, factors = [], []
    for (lo, hi), cov in zip(groups, covariances):
        cov = (cov + cov.conj().T) * .5
        diagonal = np.maximum(np.real(np.diag(cov)), 1e-20)
        cov *= 1.0 - shrinkage
        cov[np.diag_indices(bands)] += shrinkage * diagonal + 1e-10 * diagonal.mean()
        frequencies.append((lo + hi - 1) * .5 * rate / nperseg)
        factors.append(linalg.cholesky(cov, lower=True, check_finite=False))
    frequencies = np.asarray([0.0, *frequencies, rate / 2], dtype=np.float64)
    factors = np.asarray([factors[0], *factors, factors[-1]], dtype=np.complex128)
    # Linear interpolation between factors changes each band's variance; account for it analytically.
    energy = np.zeros(bands)
    for i, width in enumerate(np.diff(frequencies)):
        left, right = factors[i], factors[i + 1]
        energy += width / 3.0 * (np.sum(np.abs(left)**2, axis=1) + np.sum(np.abs(right)**2, axis=1)
                                 + np.real(np.sum(left * right.conj(), axis=1)))
    factors /= np.sqrt(np.maximum(energy, 1e-20))[None, :, None]
    return (frequencies.astype(np.float32), factors.astype(np.complex64),
            {'welch_window_samples': int(nperseg), 'welch_windows': int(len(starts)),
             'welch_resolution_hz': float(rate / nperseg), 'frequency_groups': len(groups),
             'diagonal_shrinkage': float(shrinkage)})


def fit(envelopes, envelope_rate=2000.0) -> dict:
    """The Gaussian-copula statistics of non-negative band envelopes, shape (bands, frames): each band's marginal
    quantiles and the smoothed complex cross-spectra of the bands' latent modulations. No frame sequence is kept."""
    a = np.asarray(envelopes, dtype=np.float32)
    rate = float(envelope_rate)
    if a.ndim != 2 or a.shape[0] < 1:
        raise ValueError('envelopes must have shape (bands, time).')
    if not np.isfinite(rate) or rate <= 0 or a.shape[1] < max(512, rate):
        raise ValueError('Use a positive envelope rate and at least one second of data.')
    if not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError('Envelopes must contain finite nonnegative amplitudes.')
    latent_knots = np.linspace(-4.25, 4.25, LATENT_KNOTS, dtype=np.float64)
    probabilities = special.ndtr(latent_knots)
    quantiles = np.quantile(a, probabilities, axis=1).T
    z = np.empty_like(a)
    constant = np.ptp(a, axis=1) <= np.maximum(np.max(a, axis=1) * 1e-7, 1e-20)
    for b in range(a.shape[0]):
        # A repeated quantile, as in a silent band, sits at the middle of its probability mass in the latent.
        values, inverse = np.unique(quantiles[b], return_inverse=True)
        centers = np.bincount(inverse, weights=latent_knots) / np.bincount(inverse)
        if constant[b]:
            z[b] = 0.0
        else:
            z[b] = np.interp(a[b], values, centers)
            z[b] -= z[b].mean(dtype=np.float64)
            z[b] /= max(float(z[b].std(dtype=np.float64)), 1e-8)
    frequencies, factors, details = _spectral_factors(z, rate)
    details.update({
        'method': 'Gaussian envelope copula with smoothed complex modulation cross-spectra',
        'source_traces_stored': False,
        'source_rows_sampled_at_render': False,
        'physical_drop_density_calibrated': False,
        'limitations': [
            'Univariate rise/decay asymmetry and higher-order event sequences are not fitted.',
            'Rare extremes beyond fitted quantile knots are bounded.',
            'Stationary fit does not represent evolving storm regimes.',
            'Low-frequency estimates depend on training duration and Welch resolution.',
            'No perceptual or physical validation is implied by fitted statistics.'
        ],
    })
    return {
        'model_version': np.asarray(1, dtype=np.int32),
        'envelope_rate': np.asarray(rate, dtype=np.float64),
        'bands': np.asarray(a.shape[0], dtype=np.int32),
        'training_seconds': np.asarray(a.shape[1] / rate, dtype=np.float64),
        'latent_knots': latent_knots.astype(np.float32),
        'amplitude_quantiles': quantiles.astype(np.float32),
        'modulation_frequency_hz': frequencies,
        'modulation_cholesky': factors,
        'constant_bands': constant,
        'mean_amplitude': a.mean(axis=1, dtype=np.float64).astype(np.float32),
        'rms_amplitude': np.sqrt(np.mean(a.astype(np.float64)**2, axis=1)).astype(np.float32),
        'metadata_json': np.asarray(json.dumps(details, sort_keys=True)),
    }


def generate(model, seconds, seed):
    """Fresh positive band envelopes, shape (bands, round(seconds * envelope_rate)), drawn from a model's statistics."""
    rate = float(_scalar(model, 'envelope_rate'))
    duration = float(seconds)
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError('seconds must be positive and finite.')
    requested = int(round(duration * rate))
    if requested < 1:
        raise ValueError('Duration must produce at least one sample.')
    bands = int(_scalar(model, 'bands'))
    nodes = np.asarray(model['modulation_frequency_hz'], dtype=np.float64)
    factors = np.asarray(model['modulation_cholesky'], dtype=np.complex64)
    # Eight seconds of padding keep the periodic FFT boundaries out of the clip and leave room for sub-hertz variation.
    padding = int(round(8.0 * rate))
    n = fft.next_fast_len(requested + 2 * padding)
    frequency = fft.rfftfreq(n, 1 / rate)
    result_spectrum = np.zeros((bands, len(frequency)), dtype=np.complex64)
    rng = np.random.default_rng(seed)
    scale = np.sqrt(rate * n / 2.0)
    for i in range(len(nodes) - 1):
        lo = max(1, int(np.searchsorted(frequency, nodes[i], side='left')))
        hi = int(np.searchsorted(frequency, nodes[i + 1], side='left'))
        hi = min(hi, len(frequency) - 1)
        if hi <= lo:
            continue
        weight = ((frequency[lo:hi] - nodes[i]) / max(nodes[i + 1] - nodes[i], 1e-20)).astype(np.float32)
        innovation = (rng.standard_normal((bands, hi - lo), dtype=np.float32)
                      + 1j * rng.standard_normal((bands, hi - lo), dtype=np.float32)) / np.sqrt(2.0)
        result_spectrum[:, lo:hi] = scale * ((factors[i] @ innovation) * (1.0 - weight)
                                             + (factors[i + 1] @ innovation) * weight)
    # A zero DC bin holds the latent mean at zero; the single Nyquist bin is left at zero as well.
    latent = fft.irfft(result_spectrum, n=n, axis=1, workers=2)
    latent = latent[:, padding:padding + requested]
    knots = np.asarray(model['latent_knots'], dtype=np.float64)
    quantiles = np.asarray(model['amplitude_quantiles'], dtype=np.float64)
    out = np.empty((bands, requested), dtype=np.float32)
    tiny = np.finfo(np.float32).tiny
    for b in range(bands):
        out[b] = np.maximum(np.interp(latent[b], knots, quantiles[b]), tiny)
    if not np.all(np.isfinite(out)):
        raise RuntimeError('Nonfinite generated envelopes.')
    return out


def correct_intensity(env, model):
    """Version 3's monotone map of the summed band intensity onto the recording's, keeping spectral proportions."""
    radius = np.sqrt(np.sum(env * env, axis=0))
    desired = np.interp(radius, model['radial_pilot_quantiles'], model['radial_reference_quantiles'])
    env *= np.power(desired / np.maximum(radius, 1e-20), float(model['radial_correction_strength']))[None, :]
    return env


def seeds(seed):
    """The envelope and carrier seeds a render's seed spawns."""
    envelope_seed, carrier_seed = np.random.SeedSequence(seed).generate_state(2)
    return int(envelope_seed), int(carrier_seed)


def synthesize_guarded(model, seconds, seed):
    """A new render of a model with GUARD_SECONDS of audio before and after the clip, and the seeds it used."""
    if not np.isfinite(seconds) or seconds < 1 or seconds > MAX_SECONDS:
        raise ValueError('Use 1-300 seconds per offline render.')
    if str(np.asarray(model.get('envelope_coordinates', ''))) != 'band_amplitudes_with_radial_quantile_correction':
        raise ValueError('Use a version 3 model, with its intensity correction.')
    envelope_seed, carrier_seed = seeds(seed)
    env = correct_intensity(generate(model, seconds + 2 * GUARD_SECONDS, envelope_seed), model)
    audio = render_envelopes(env, int(np.asarray(model['envelope_rate'])), carrier_seed)
    return audio, dict(seed=int(seed), envelope_seed=envelope_seed, carrier_seed=carrier_seed,
                       source_audio_read=False, recorded_timing_replayed=False)


def crop(audio, seconds):
    """The clip inside a guarded render."""
    return audio[round(GUARD_SECONDS * SAMPLE_RATE):round((GUARD_SECONDS + seconds) * SAMPLE_RATE)].copy()


def synthesize(model, seconds, seed):
    """A new render of a model, seconds long (float32 at 48 kHz, digital level as the model gives it)."""
    audio, _ = synthesize_guarded(model, seconds, seed)
    return crop(audio, seconds)


def master(model, seconds, seed, gain_db=0.0):
    """The renderer's master as the rain work wrote it: a fixed gain, a whole-file attenuation to keep the peak at
    0.85 of full scale if it would exceed it, and 50 ms raised-cosine fades. No compression or limiting."""
    if not np.isfinite(gain_db) or gain_db < -120 or gain_db > 60:
        raise ValueError('Playback gain must be finite and between -120 and +60 dB.')
    x = synthesize(model, seconds, seed)
    x *= 10**(gain_db / 20)
    peak = float(np.max(np.abs(x)))
    x *= min(1., .85 / max(peak, 1e-20))
    fade_len = min(round(.05 * SAMPLE_RATE), len(x) // 2)
    fade = .5 - .5 * np.cos(np.linspace(0, np.pi, fade_len))
    x[:fade_len] *= fade
    x[-fade_len:] *= fade[::-1]
    return x.astype(np.float32)


def carrier_pair(model, seconds, seed, envelope_correlation=.85):
    """Two guarded renders for a stereo field: the first exactly synthesize_guarded's; the second from envelopes whose
    latent innovations correlate with the first's by envelope_correlation, on independent carriers. Each keeps the
    model's marginal statistics."""
    if not 0 <= envelope_correlation <= 1:
        raise ValueError('Invalid envelope correlation.')
    envelope_seed, carrier_seed = seeds(seed)
    other_envelope_seed, other_carrier_seed = seeds(seed + 190003)
    a = generate(model, seconds + 2 * GUARD_SECONDS, envelope_seed)
    fresh = generate(model, seconds + 2 * GUARD_SECONDS, other_envelope_seed)
    b = np.empty_like(a)
    knots = model['latent_knots']
    for i, q in enumerate(model['amplitude_quantiles']):
        values, indices = np.unique(q, return_index=True)
        if len(values) == 1:
            b[i] = values[0]
            continue
        z_a = np.interp(a[i], values, knots[indices])
        z_b = np.interp(fresh[i], values, knots[indices])
        b[i] = np.interp(envelope_correlation * z_a + np.sqrt(1 - envelope_correlation**2) * z_b, knots, q)
    del fresh
    a, b = correct_intensity(a, model), correct_intensity(b, model)
    rate = int(model['envelope_rate'])
    pair = np.column_stack((render_envelopes(a, rate, carrier_seed), render_envelopes(b, rate, other_carrier_seed)))
    return pair, dict(seed=int(seed), envelope_seeds=[envelope_seed, other_envelope_seed],
                      carrier_seeds=[carrier_seed, other_carrier_seed],
                      latent_envelope_correlation=envelope_correlation, source_audio_read=False,
                      recorded_timing_replayed=False)


# --- fitting from a recording -------------------------------------------------------------------------------------

def decode(path) -> np.ndarray:
    """A recording's first channel at 48 kHz as float32. PCM WAV is read directly, scaled as ffmpeg scales it to
    32-bit float; other formats (the grass recording is an MP3) need ffmpeg."""
    path = Path(path)
    if path.suffix.lower() == '.wav':
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', wavfile.WavFileWarning)     # broadcast-wave metadata chunks
            rate, raw = wavfile.read(path)
        if rate != SAMPLE_RATE or not np.issubdtype(raw.dtype, np.signedinteger):
            raise ValueError('Expected signed integer PCM at 48 kHz')
        first = raw if raw.ndim == 1 else raw[:, 0]
        return (first.astype(np.float64) / float(-np.iinfo(raw.dtype).min)).astype(np.float32)
    if shutil.which('ffmpeg') is None:
        raise RuntimeError(f'Decoding {path.name} needs ffmpeg')
    result = subprocess.run(['ffmpeg', '-v', 'error', '-i', str(path), '-af', 'pan=mono|c0=c0', '-ar',
                             str(SAMPLE_RATE), '-f', 'f32le', '-'], capture_output=True, check=True)
    return np.frombuffer(result.stdout, dtype='<f4').copy()


def train(source, start=2., seconds=36., bands=64, envelope_rate=4000) -> dict:
    """Fit a version 3 model to seconds of a recording from start, with two seconds of guard on each side."""
    if not np.isfinite(start) or start < 2 or not np.isfinite(seconds) or seconds < 1:
        raise ValueError('Use a finite training start >=2 s and duration >=1 s.')
    x = decode(source)
    pad = 2.
    lo, hi = round((start - pad) * SAMPLE_RATE), round((start + seconds + pad) * SAMPLE_RATE)
    if lo < 0 or hi > len(x):
        raise ValueError('Training interval requires two seconds of source guard at each end.')
    env = extract_envelopes(x[lo:hi], bands, envelope_rate)
    env = env[:, round(pad * envelope_rate):round((pad + seconds) * envelope_rate)]
    total = np.sqrt(np.sum(env * env, axis=0))
    model = fit(env, envelope_rate)
    pilot = generate(model, PILOT['seconds'], PILOT['seed'])
    pilot_total = np.sqrt(np.sum(pilot * pilot, axis=0))
    probabilities = special.ndtr(model['latent_knots'])
    model['radial_reference_quantiles'] = np.quantile(total, probabilities).astype(np.float32)
    model['radial_pilot_quantiles'] = np.quantile(pilot_total, probabilities).astype(np.float32)
    model['radial_correction_strength'] = np.asarray(1., dtype=np.float32)
    model['radial_calibration_seed'] = np.asarray(PILOT['seed'])
    model['radial_calibration_seconds'] = np.asarray(PILOT['seconds'])
    wave = x[round(start * SAMPLE_RATE):round((start + seconds) * SAMPLE_RATE)]
    model.update({
        'audio_sample_rate': np.asarray(SAMPLE_RATE),
        'reference_rms': np.asarray(float(np.sqrt(np.mean(wave.astype(np.float64)**2)))),
        'source_sha256': np.asarray(sha256(source)),
        'source_file': np.asarray(Path(source).name),
        'training_start_seconds': np.asarray(start),
        'training_duration_seconds': np.asarray(seconds),
        'renderer_version': np.asarray(3),
        'audio_bands': np.asarray(bands),
        'envelope_coordinates': np.asarray('band_amplitudes_with_radial_quantile_correction'),
    })
    return model
