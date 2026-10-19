"""Stereo and distance: the rain's sound carried from the ground to a pair of microphones through the target's air.

The receivers are two omnidirectional microphones, 8 cm apart and 1.6 m above flat ground, by default: x runs to the
right, y forward and z up, in metres. From a source the pressure falls as 1/r, arrives after r/c, and loses
frequency-dependent energy to the air's molecular relaxation, which the thunder model computes by ISO 9613-1 with its
oxygen and nitrogen terms scaled to the air's own fractions (atmosphere/electricity/thunder.py, absorption_db_per_m);
c is the sound speed in the same air. An area of independent sources radiates power in proportion to its size.

The rain falling around the listener is a sheet of independent sources from 0.2 to 40 m. Its sound at the two
microphones is a stationary 2 x 2 cross-spectral covariance, integrated over 4,096 cells of the sheet
(field_covariance), and the stereo rain is two fresh texture renders mixed by that covariance's square root. One rain
patch heard from a place (propagate_mono) carries its own absolute delay. Ground reflection, refraction by wind and
temperature, scattering, the radiation pattern of an impact, and a listener's head and ears are left out: the
receivers are a spaced pair of microphones.

The model is the rain work's of 13 October 2026 (propagation.py, field_covariance.py and render_spatial.py), brought
into the repository on 17 October with the air's absorption read from the thunder model. Its 8 cm spacing is an
effective value chosen there by matching the grass recording's coherence between channels in three bands from 250 Hz
to 12 kHz (models/models.json); the recorder's real geometry is unknown.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import fft
from scipy.signal import oaconvolve

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from atmosphere.electricity import thunder  # noqa: E402

SAMPLE_RATE = 48000
SPACING_M = 0.08
LISTENER_HEIGHT_M = 1.6
RADIUS_EDGES_M = (0.2, 2.0, 6.0, 16.0, 40.0)
AZIMUTHS = 128
RADIAL_SUBDIVISIONS = 8
FIELD_FREQUENCIES = np.linspace(0, SAMPLE_RATE / 2, 2401)


def air_absorption_db_per_m(frequency_hz, air):
    """The air's absorption (dB/m) at each frequency, from the thunder model."""
    f = np.asarray(frequency_hz, dtype=float)
    if np.any(~np.isfinite(f)) or np.any(f < 0):
        raise ValueError('Frequencies must be finite and nonnegative')
    return thunder.absorption_db_per_m(f, air['temperature_k'], air['pressure_pa'], air['vapour_percent'],
                                       air['x_o2'], air['x_n2'])


def stereo_geometry(source_xyz, sound_speed_m_s, spacing_m=SPACING_M, listener_height_m=LISTENER_HEIGHT_M,
                    reference_distance_m=1.0):
    """The direct paths from a point to the two receivers, left then right: distances, delays and 1/r gains."""
    source = np.asarray(source_xyz, dtype=float)
    if source.shape != (3,) or not np.all(np.isfinite(source)):
        raise ValueError('source_xyz must contain three finite metre coordinates')
    if sound_speed_m_s <= 0 or spacing_m < 0 or reference_distance_m <= 0:
        raise ValueError('Sound speed and reference distance must be positive, the spacing non-negative')
    receivers = np.array([[-spacing_m / 2, 0, listener_height_m], [spacing_m / 2, 0, listener_height_m]], dtype=float)
    distances = np.linalg.norm(source[None, :] - receivers, axis=1)
    if np.min(distances) < 0.05:
        raise ValueError('Source must be at least 5 cm from both receivers')
    return dict(source_xyz_m=source.tolist(), receiver_xyz_m=receivers.tolist(), channel_order=['left', 'right'],
                distances_m=distances.tolist(), delays_s=(distances / sound_speed_m_s).tolist(),
                pressure_gains=(reference_distance_m / distances).tolist(),
                interchannel_delay_right_minus_left_s=float((distances[1] - distances[0]) / sound_speed_m_s),
                sound_speed_m_s=float(sound_speed_m_s), spacing_m=float(spacing_m))


def fractional_delay(audio, delay_samples, half_length=64):
    """A positive delay of the same length with zeros shifted in: exact for whole samples, and a 129-tap
    Kaiser-windowed sinc for the fraction, its own latency removed."""
    x = np.asarray(audio, dtype=float)
    if x.ndim != 1 or not np.all(np.isfinite(x)) or not np.isfinite(delay_samples) or delay_samples < 0:
        raise ValueError('Expected finite mono audio and a non-negative delay')
    n = len(x)
    integer = int(math.floor(delay_samples))
    fraction = float(delay_samples - integer)
    if abs(fraction) < 1e-12:
        out = np.zeros(n, dtype=float)
        if integer < n:
            out[integer:] = x[:n - integer]
        return out
    k = np.arange(-half_length, half_length + 1, dtype=float)
    kernel = np.sinc(k - fraction) * np.kaiser(len(k), 8.6)
    kernel /= kernel.sum()
    filtered = oaconvolve(x, kernel, mode='full')
    out = np.zeros(n, dtype=float)
    source_start = max(0, half_length - integer)
    output_start = max(0, integer - half_length)
    count = max(0, min(n - output_start, len(filtered) - source_start))
    out[output_start:output_start + count] = filtered[source_start:source_start + count]
    return out


def absorption_kernel(distance_m, air, fs=SAMPLE_RATE, n_fft=8192, taps=1024):
    """A causal minimum-phase FIR with the air's attenuation over a distance; the time of flight is applied apart."""
    if distance_m < 0 or taps >= n_fft // 2:
        raise ValueError('Invalid absorption distance or FIR length')
    frequencies = np.fft.rfftfreq(n_fft, 1 / fs)
    logmag = -np.log(10) / 20 * air_absorption_db_per_m(frequencies, air) * distance_m
    cepstrum = np.fft.irfft(logmag, n_fft)
    causal = np.zeros(n_fft)
    causal[0] = cepstrum[0]
    causal[1:n_fft // 2] = 2 * cepstrum[1:n_fft // 2]
    causal[n_fft // 2] = cepstrum[n_fft // 2]
    full = np.fft.irfft(np.exp(np.fft.rfft(causal)), n_fft)
    kernel = full[:taps].copy()
    kernel[-128:] *= np.cos(np.linspace(0, np.pi / 2, 128))**2
    kernel /= kernel.sum()
    obtained = abs(np.fft.rfft(kernel, n_fft))
    band = (frequencies >= 20) & (frequencies <= min(20000, fs * 0.45))
    error = 20 * np.log10(np.maximum(obtained[band], 1e-20) / np.exp(logmag)[band])
    return kernel, dict(fir_taps=int(taps), maximum_magnitude_error_db_20hz_20khz=float(np.max(np.abs(error))),
                        discarded_impulse_energy_fraction=float(np.sum(full[taps:]**2) / np.sum(full**2)),
                        absorption_db_at_1_4_10_khz=(air_absorption_db_per_m(np.array([1000., 4000., 10000.]), air)
                                                     * distance_m).tolist())


def propagate_mono(audio, source_xyz, air, fs=SAMPLE_RATE, spacing_m=SPACING_M, listener_height_m=LISTENER_HEIGHT_M,
                   reference_distance_m=1.0, source_area_m2=1.0, absorption=True):
    """One independent patch's sound at the two receivers, (samples, 2), with its geometry. The patch's audio is its
    pressure at the reference distance per square root of area; its area weighs the amplitude by sqrt(area)."""
    x = np.asarray(audio, dtype=float)
    if x.ndim != 1 or not np.all(np.isfinite(x)):
        raise ValueError('Expected finite mono audio')
    if not np.isfinite(source_area_m2) or source_area_m2 < 0:
        raise ValueError('Source area must be finite and nonnegative')
    geometry = stereo_geometry(source_xyz, air['sound_speed_m_s'], spacing_m, listener_height_m, reference_distance_m)
    stereo = np.zeros((len(x), 2), dtype=float)
    filters = []
    for channel, (distance, delay, gain) in enumerate(zip(geometry['distances_m'], geometry['delays_s'],
                                                          geometry['pressure_gains'])):
        if absorption:
            kernel, detail = absorption_kernel(distance, air, fs)
            filtered = oaconvolve(x, kernel, mode='full')[:len(x)]
        else:
            filtered, detail = x, dict(disabled=True)
        stereo[:, channel] = fractional_delay(filtered, delay * fs) * gain * math.sqrt(source_area_m2)
        filters.append(detail)
    return stereo, dict(geometry, source_area_m2=float(source_area_m2), absorption_filters=filters)


def ring_patches(edges_m=(0.4, 2, 6, 16, 40), sectors=8, listener_height_m=LISTENER_HEIGHT_M, source_height_m=0.0,
                 rotation_rad=0.0):
    """Annulus-sector cells of a uniform sheet, each placed at the radius that keeps the annulus's unabsorbed power
    at the central receiver, the integral of dA / (r^2 + h^2)."""
    edges = np.asarray(edges_m, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or edges[0] < 0 or np.any(np.diff(edges) <= 0) or sectors < 2:
        raise ValueError('Expected increasing nonnegative ring edges and at least two sectors')
    h2 = (listener_height_m - source_height_m)**2
    if h2 == 0 and edges[0] == 0:
        raise ValueError('Coplanar receiver and source sheet require an excluded central radius')
    result = []
    for ring, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        radius = np.sqrt((b * b - a * a) / np.log((b * b + h2) / (a * a + h2)) - h2)
        area = np.pi * (b * b - a * a) / sectors
        for sector in range(sectors):
            angle = rotation_rad + 2 * np.pi * (sector + .5) / sectors
            result.append(dict(ring=int(ring), sector=int(sector), inner_radius_m=float(a), outer_radius_m=float(b),
                               source_area_m2=float(area), azimuth_rad=float(angle),
                               source_xyz_m=[float(radius * np.sin(angle)), float(radius * np.cos(angle)),
                                             float(source_height_m)]))
    return result


def field_covariance(frequency_hz, air, spacing_m=SPACING_M, listener_height_m=LISTENER_HEIGHT_M,
                     radius_edges_m=RADIUS_EDGES_M, azimuth_count=AZIMUTHS, radial_subdivisions=RADIAL_SUBDIVISIONS):
    """The two receivers' powers and cross-spectrum E[R conj(L)] from a uniform sheet of independent sources of unit
    power per square metre at 1 m, area weights included."""
    frequency = np.asarray(frequency_hz, dtype=float)
    if frequency.ndim != 1 or frequency.size == 0 or np.any(~np.isfinite(frequency)) or np.any(frequency < 0):
        raise ValueError('frequency_hz must be a nonempty vector of finite nonnegative frequencies')
    if not isinstance(azimuth_count, (int, np.integer)) or azimuth_count < 4:
        raise ValueError('azimuth_count must be an integer >= 4')
    if not isinstance(radial_subdivisions, (int, np.integer)) or radial_subdivisions < 1:
        raise ValueError('radial_subdivisions must be a positive integer')
    speed = float(air['sound_speed_m_s'])
    edges = np.asarray(radius_edges_m, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(~np.isfinite(edges)) or edges[0] < 0 or np.any(np.diff(edges) <= 0):
        raise ValueError('radius_edges_m must be finite increasing nonnegative boundaries')
    # Eight radial cells per ring keep a broad annulus from acting as one coherent radius at high frequencies.
    fine_edges = []
    for a, b in zip(edges[:-1], edges[1:]):
        cells = np.geomspace(a, b, radial_subdivisions + 1) if a > 0 else np.linspace(a, b, radial_subdivisions + 1)
        fine_edges.extend(cells[:-1])
    fine_edges.append(edges[-1])
    patches = ring_patches(fine_edges, sectors=int(azimuth_count), listener_height_m=listener_height_m)
    alpha = air_absorption_db_per_m(frequency, air)
    power_left = np.zeros_like(frequency)
    power_right = np.zeros_like(frequency)
    cross_rl = np.zeros(len(frequency), dtype=complex)
    for patch in patches:
        # A different azimuthal offset for each radial cell keeps the spokes of the quadrature apart.
        theta = patch['azimuth_rad'] + patch['ring'] * np.pi / azimuth_count
        radius = float(np.hypot(*patch['source_xyz_m'][:2]))
        geometry = stereo_geometry([radius * np.sin(theta), radius * np.cos(theta), 0.0], speed, spacing_m,
                                   listener_height_m)
        left, right = geometry['distances_m']
        area = patch['source_area_m2']
        power_left += area / left**2 * np.power(10.0, -alpha * left / 10.0)
        power_right += area / right**2 * np.power(10.0, -alpha * right / 10.0)
        cross_rl += (area / (left * right) * np.power(10.0, -alpha * (left + right) / 20.0)
                     * np.exp(-2j * np.pi * frequency * (right - left) / speed))
    return dict(power_left=power_left, power_right=power_right, cross_rl=cross_rl,
                metadata=dict(spacing_m=float(spacing_m), listener_height_m=float(listener_height_m),
                              radius_edges_m=edges.tolist(), azimuth_count=int(azimuth_count),
                              radial_subdivisions=int(radial_subdivisions), cells=len(patches),
                              field_area_m2=float(np.pi * (edges[-1]**2 - edges[0]**2)), sound_speed_m_s=speed))


def covariance_matrix(cov):
    c = np.empty((len(cov['power_left']), 2, 2), complex)
    c[:, 0, 0] = cov['power_left']
    c[:, 1, 1] = cov['power_right']
    c[:, 1, 0] = cov['cross_rl']
    c[:, 0, 1] = cov['cross_rl'].conj()
    return c


def covariance_factor(cov, calibration=None):
    """The positive square root H of each frequency's 2 x 2 covariance, C = H H*, times an optional source
    calibration."""
    c = covariance_matrix(cov)
    s = np.sqrt(np.maximum(np.real(c[:, 0, 0] * c[:, 1, 1]) - abs(c[:, 0, 1])**2, 0))
    den = np.sqrt(np.maximum(np.real(c[:, 0, 0] + c[:, 1, 1]) + 2 * s, 1e-30))
    h = c.copy()
    h[:, 0, 0] += s
    h[:, 1, 1] += s
    h /= den[:, None, None]
    if calibration is not None:
        h *= np.asarray(calibration)[:, None, None]
    return h


def apply_matrix(pair, freq, h):
    """Two input channels mixed by a frequency-dependent 2 x 2 matrix, by one FFT over the guarded clip."""
    n = fft.next_fast_len(len(pair))
    ff = fft.rfftfreq(n, 1 / SAMPLE_RATE)
    spec = fft.rfft(pair, n, axis=0, workers=2)
    out = np.zeros_like(spec, dtype=complex)
    for i in range(2):
        for j in range(2):
            hij = np.interp(ff, freq, h[:, i, j].real) + 1j * np.interp(ff, freq, h[:, i, j].imag)
            hij[0] = hij[0].real
            hij[-1] = hij[-1].real
            out[:, i] += spec[:, j] * hij
    return fft.irfft(out, n, axis=0, workers=2)[:len(pair)]
