"""The lunar rain's excitation relative to Earth's: an impact-force transfer applied to the Earth texture.

Each impact delivers the drop's incident normal momentum J = m U as a force pulse of duration tau = a D / U, shaped as
a gamma pulse of order k = 1.5: a square-root rise and an exponential tail (after Soto et al. 2014, who measure a peak
force near pi rho (D/2)^2 U^2 lasting about D / U, and Mitchell et al. 2019, who measure the early square-root rise and
the late exponential decay). Its spectrum is |F(f)|^2 = J^2 (1 + (2 pi f tau)^2)^-k. Three equally likely widths,
a = 0.24, 0.32 and 0.44, stand for the spread of real pulses, and on puddles the pulses are 1.6 times wider. For
independent arrivals at number flux j per class the force on a unit of surface has the one-sided spectrum

    Q(f) = 2 sum_classes j J^2 mean_widths (1 + (2 pi f tau)^2)^-k.

The same surface, the same path to the microphone and a linear response would leave the ratio of two such spectra as
the ratio of the sounds' powers, so the Moon's rain is the Earth texture filtered by sqrt(Q_moon / Q_earth): a causal
minimum-phase filter that keeps the texture's timing and its material's resonances, with no global slowing, pitch
shift or invented bubble. The ratio is relative: incident momentum fixes no absolute sound pressure.

Drops of 8-12 mm, when a scene carries them, are rare enough to be heard one by one, so each is scheduled as an event
of its own (large_drop_events) over the filtered texture.

The model is the rain work's of 13 October 2026 (acoustic_transfer.py and the large-drop scheduler of renderer.py),
brought into the repository on 17 October and driven by the rain domain's arriving drops.
"""
from __future__ import annotations

import numpy as np
from scipy import fft, signal
from scipy.ndimage import gaussian_filter1d

PULSE_SHAPE_K = 1.5
WIDTH_COEFFICIENTS = np.asarray([.24, .32, .44], dtype=np.float64)
PUDDLE_WIDTH_SCALE = 1.6
SAMPLE_RATE = 48000
# Frequencies on which transfers are evaluated: zero and 640 points from 5 Hz to the Nyquist frequency, logarithmic.
FREQUENCIES = np.r_[0., np.geomspace(5., SAMPLE_RATE / 2, 640)]


def tau_scale(surface):
    return PUDDLE_WIDTH_SCALE if surface == 'puddles' else 1.0


def _frequencies(frequency_hz):
    f = np.asarray(frequency_hz, dtype=np.float64)
    if f.ndim != 1 or len(f) == 0 or not np.isfinite(f).all() or np.any(f < 0):
        raise ValueError('frequency_hz must be a nonempty finite nonnegative vector')
    return f


def force_psd(population, frequency_hz, width_scale=1.0):
    """The one-sided continuous spectrum of the impacts' force on one square metre (N^2/Hz per m^2), the mean force at
    zero frequency excluded. A population with no arrivals gives zero."""
    f = _frequencies(frequency_hz)
    scale = float(width_scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('width_scale must be positive and finite')
    d, j, u = population.diameter_m, population.number_flux_m2_s, population.normal_speed_m_s
    active = j > 0
    out = np.zeros_like(f)
    if not np.any(active):
        return out
    d, j, u = d[active], j[active], u[active]
    mass = population.liquid_density_kg_m3 * (np.pi / 6.) * d**3
    impulse_power = j * (mass * u)**2
    for lo in range(0, len(f), 16384):
        ff = f[lo:lo + 16384]
        shape = np.zeros((len(d), len(ff)), dtype=np.float64)
        for coefficient in WIDTH_COEFFICIENTS:
            tau = coefficient * scale * d / u
            shape += (1. + (2. * np.pi * tau[:, None] * ff[None, :])**2)**(-PULSE_SHAPE_K)
        out[lo:lo + len(ff)] = 2. * (impulse_power @ shape) / len(WIDTH_COEFFICIENTS)
    return out


def source_spectrum(scene, frequency_hz):
    """A scene's stationary force spectrum: its direct rain and drips, each counted once; the tail is scheduled."""
    width = tau_scale(scene.surface)
    total = np.zeros(len(_frequencies(frequency_hz)))
    for part in (scene.direct, scene.drips):
        if part is not None:
            total = total + force_psd(part, frequency_hz, width)
    return total


def make_transfer(earth, moon, frequency_hz=FREQUENCIES):
    """The amplitude ratio sqrt(Q_moon / Q_earth) between two scenes on the same surface, and what it assumed."""
    f = _frequencies(frequency_hz)
    if np.any(np.diff(f) <= 0):
        raise ValueError('Transfer frequency grid must be strictly increasing')
    if earth.surface != moon.surface:
        raise ValueError('A transfer compares two scenes on the same surface')
    pe, pm = source_spectrum(earth, f), source_spectrum(moon, f)
    if not np.all(pe > 0):
        raise ValueError('A relative transfer needs Earth forcing at every frequency')
    magnitude = np.sqrt(pm / pe)
    if not np.isfinite(magnitude).all():
        raise RuntimeError('Nonfinite relative source magnitude')
    identity = bool(np.array_equal(pe, pm))
    if identity:
        magnitude = np.ones_like(magnitude)
    probes = [x for x in (0., 100., 500., 1000., 4000., 10000., 16000.) if f[0] <= x <= f[-1]]
    gains = np.interp(probes, f, magnitude)
    return magnitude, dict(
        method='sqrt of the ratio of the scenes\' impact-force spectra',
        surface=moon.surface, pulse_shape_k=PULSE_SHAPE_K, width_coefficients=WIDTH_COEFFICIENTS.tolist(),
        width_scale=tau_scale(moon.surface), exact_identity=identity, frequency_hz=probes,
        amplitude_ratio=gains.tolist(), gain_db=[float(20 * np.log10(g)) if g > 0 else None for g in gains],
        earth=earth.describe(), moon=moon.describe())


def _minimum_phase_kernel(frequency_hz, magnitude, fs):
    """A finite causal FIR with a smooth positive magnitude, by the real cepstrum, with its tail audited."""
    largest = float(np.max(magnitude))
    floor = max(largest * 1e-12, np.finfo(float).tiny)
    tail_target = 1e-10
    accepted = None
    # Design FFT sizes, independent of the clip's length; the FIR is applied by linear convolution.
    for n in (32768, 65536, 131072, 262144, 524288):
        frequencies = fft.rfftfreq(n, 1. / fs)
        target = np.interp(frequencies, frequency_hz, magnitude)
        cepstrum = fft.irfft(np.log(np.maximum(target, floor)), n=n)
        causal = np.zeros(n, dtype=np.float64)
        causal[0] = cepstrum[0]
        causal[1:n // 2] = 2. * cepstrum[1:n // 2]
        causal[n // 2] = cepstrum[n // 2]
        full = fft.irfft(np.exp(fft.rfft(causal)), n=n)
        energy = full * full
        total = float(np.sum(energy))
        if not np.isfinite(total) or total <= 0:
            raise RuntimeError('Invalid minimum-phase impulse energy')
        cumulative = np.cumsum(energy) / total
        essential = int(np.searchsorted(cumulative, 1. - tail_target)) + 1
        if essential <= n // 2 - 256:
            length = max(64, essential + 128)
            kernel = full[:length].copy()
            fade = min(128, length // 4)
            kernel[-fade:] *= np.cos(np.linspace(0., np.pi / 2., fade))**2
            discarded = float(np.sum(energy[length:]) / total)
            late_half = float(np.sum(energy[n // 2:]) / total)
            if discarded <= 1e-8 and late_half <= 1e-8:
                accepted = (kernel, n, target, discarded, late_half)
                break
    if accepted is None:
        raise RuntimeError('Relative filter has unresolved impulse tail; use a smoother frequency response')
    kernel, n, target, discarded, late_half = accepted
    actual = np.abs(fft.rfft(kernel, n=n))
    compare = target > max(float(target.max()) * 1e-9, floor)
    error_db = 20. * np.log10(np.maximum(actual[compare], floor) / target[compare])
    return kernel, dict(minimum_phase_method='real cepstrum of the relative magnitude; finite causal FIR',
                        design_fft_samples=n, impulse_samples=len(kernel), impulse_seconds=len(kernel) / float(fs),
                        discarded_tail_energy_fraction=discarded, design_late_half_energy_fraction=late_half,
                        magnitude_error_db_max=float(np.max(np.abs(error_db))),
                        magnitude_error_db_rms=float(np.sqrt(np.mean(error_db**2))))


def apply_transfer(audio, frequency_hz, magnitude, fs=SAMPLE_RATE):
    """Audio (mono, or channels last) filtered causally by the magnitude, keeping its length. The grid must run from
    0 Hz to the Nyquist frequency. A magnitude of exactly one returns a copy; render with guard audio and crop after."""
    x = np.asarray(audio)
    f = _frequencies(frequency_hz)
    mag = np.asarray(magnitude, dtype=np.float64)
    fs = float(fs)
    if x.ndim not in (1, 2) or len(x) == 0 or not np.isfinite(x).all():
        raise ValueError('Audio must be a nonempty finite mono or channel-last array')
    if mag.shape != f.shape or not np.isfinite(mag).all() or np.any(mag < 0):
        raise ValueError('Magnitude must be a finite nonnegative array matching the grid')
    if len(f) < 2 or np.any(np.diff(f) <= 0) or f[0] != 0 or f[-1] < fs / 2:
        raise ValueError('Frequency grid must rise strictly from 0 Hz to the Nyquist frequency')
    if np.all(mag == 1.):
        return x.copy(), dict(identity=True, impulse_samples=1)
    if np.all(mag == mag[0]):
        return x * mag[0], dict(identity=False, constant_gain=float(mag[0]), impulse_samples=1)
    kernel, details = _minimum_phase_kernel(f, mag, fs)
    if x.ndim == 1:
        result = signal.fftconvolve(x.astype(np.float64), kernel, mode='full')[:len(x)]
    else:
        result = signal.fftconvolve(x.astype(np.float64), kernel[:, None], mode='full', axes=0)[:len(x)]
    if not np.isfinite(result).all():
        raise RuntimeError('Nonfinite filtered output')
    return result, dict(identity=False, **details)


# Each surface's assumed decay of a large drop's response (s), as the rain work set it.
LARGE_DROP_DECAY_S = {'concrete': .012, 'grass-leaves': .024, 'puddles': .040}


def large_drop_events(earth_audio, earth, moon, seconds, seed, area_m2=4.0):
    """The scene's 8-12 mm drops as separately scheduled events on an effective receiving area (a stem the length of
    the clip). Arrivals are Poisson at the tail's flux; each event is fresh broadband sound whose spectrum is the
    Earth texture's response to force (its power spectrum over the Earth scene's force spectrum) times the event's own
    force spectrum, with a random colouring of contact, random phases, an attack of tau and the surface's decay
    lengthened by the drop's capillary time. Puddle events carry a weaker secondary packet after a capillary and
    gravity delay. Area, decay, colouring and the secondary packet are declared estimates; the event energy is relative
    to the texture, like the transfer."""
    surface = moon.surface
    stem = np.zeros(round(seconds * SAMPLE_RATE), dtype=np.float64)
    if moon.tail is None:
        return stem, dict(events=[], count=0, expected_count=0.0, area_m2=area_m2)
    pop = moon.tail
    diam, flux, vel = pop.diameter_m, pop.number_flux_m2_s, pop.normal_speed_m_s
    total_rate = float(flux.sum() * area_m2)
    rng = np.random.default_rng(seed)
    count = int(rng.poisson(total_rate * seconds))
    times = np.sort(rng.uniform(0, seconds, count))
    f, p = signal.welch(earth_audio.astype(np.float64), fs=SAMPLE_RATE, nperseg=8192, noverlap=4096)
    reference_psd = np.exp(gaussian_filter1d(np.interp(FREQUENCIES, f, np.log(np.maximum(p, 1e-30))), 2.))
    reference_force = area_m2 * source_spectrum(earth, FREQUENCIES)
    response = np.sqrt(reference_psd / np.maximum(reference_force, 1e-40))
    # The event response keeps 45 Hz to 19 kHz: below, a rain texture's energy is ambience; above, the band's edge.
    response *= np.where(FREQUENCIES > 0, 1 / np.sqrt(1 + (45 / np.maximum(FREQUENCIES, 1e-20))**8), 0.)
    response *= 1 / np.sqrt(1 + (FREQUENCIES / 19000.)**12)
    rho = pop.liquid_density_kg_m3
    sigma = moon.air['surface_tension_n_m']
    gravity = moon.air['gravity_m_s2']
    base_decay = LARGE_DROP_DECAY_S[surface]
    events, volume = [], 0.
    for t0 in times:
        i = int(rng.choice(len(diam), p=flux / flux.sum()))
        d, u = float(diam[i]), float(vel[i])
        impulse = rho * np.pi * d**3 / 6 * u
        tau = float(rng.choice(WIDTH_COEFFICIENTS)) * d / u * tau_scale(surface)
        capillary = np.sqrt(rho * d**3 / sigma)
        decay = base_decay * float(rng.uniform(.75, 1.35)) + .15 * capillary
        support = min(.9, 7 * decay)
        length = round(support * SAMPLE_RATE)
        n = fft.next_fast_len(max(length, 1024))
        ff = fft.rfftfreq(n, 1 / SAMPLE_RATE)
        target = np.interp(ff, FREQUENCIES, response) * impulse / (1 + (2 * np.pi * ff * tau)**2)**.75
        knots = np.linspace(0, np.log1p(SAMPLE_RATE / 2), 14)
        target *= np.exp(np.interp(np.log1p(ff), knots, rng.normal(0, .18, len(knots))))
        target[0] = target[-1] = 0.
        event_energy = float(2 * np.trapezoid(target * target, ff))
        carrier = fft.irfft(target * np.exp(1j * rng.uniform(-np.pi, np.pi, len(ff))), n)[:length]
        tt = np.arange(length) / SAMPLE_RATE
        envelope = (1 - np.exp(-tt / max(tau, .00025))) * np.exp(-tt / decay)
        if surface == 'puddles':
            delay = float(.35 * capillary + .25 * np.sqrt(d / gravity))
            late = np.maximum(tt - delay, 0)
            envelope += .35 * (1 - np.exp(-late / .003)) * np.exp(-late / (.55 * decay)) * (tt >= delay)
        edge = min(round(.01 * SAMPLE_RATE), length // 3)
        envelope[-edge:] *= np.linspace(1, 0, edge)**2
        wave = carrier * envelope
        wave -= envelope * wave.sum() / max(envelope.sum(), 1e-20)
        wave *= np.sqrt(event_energy * SAMPLE_RATE / max(float(np.dot(wave, wave)), 1e-40))
        start = round(t0 * SAMPLE_RATE)
        stop = min(len(stem), start + len(wave))
        stem[start:stop] += wave[:stop - start]
        volume += np.pi * d**3 / 6
        events.append(dict(time_s=float(t0), diameter_mm=d * 1e3, impact_speed_m_s=u, impulse_kg_m_s=float(impulse),
                           force_time_s=tau, decay_s=decay, sample_start=start, sample_end=stop))
    return stem, dict(events=events, count=count, expected_count=total_rate * seconds, area_m2=area_m2,
                      realized_tail_rain_mm_h=volume / area_m2 / seconds * 3.6e6,
                      expected_tail_rain_mm_h=pop.rain_mm_h)
