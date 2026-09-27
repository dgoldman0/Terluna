"""Wind power at a tower: the resource in a speed distribution, what devices of several kinds take from it, the
thrust they put on the structure, and the chance that a flier crossing a rotor is struck.

A device is described by a power coefficient on its frontal area (the rotor's swept area, or a bladeless
mast's diameter times its height), a cut-in, rated and cut-out speed, a thrust coefficient while it works and a
drag coefficient when it is stopped in a storm. Its mean output is the integral of its power curve over a Weibull
distribution of wind speed. The collision probability is the single-transit probability of Band's collision risk
model for a flier crossing a rotor at its own airspeed, averaged over the rotor disc.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import gamma, log
import numpy as np


def weibull_from_percentiles(p50: float, p90: float) -> tuple:
    """Weibull shape k and scale c (m/s) with the given median and 90th percentile."""
    k = log(log(10.0) / log(2.0)) / log(p90 / p50)
    return k, p50 / log(2.0) ** (1.0 / k)


def weibull_quantile(k: float, c: float, q: float) -> float:
    return c * (-log(1.0 - q)) ** (1.0 / k)


def weibull_moment(k: float, c: float, n: float) -> float:
    """Mean of the speed to the power n."""
    return c ** n * gamma(1.0 + n / k)


def resource_w_m2(k: float, c: float, density_kg_m3: float) -> float:
    """Mean kinetic power flux of the wind, 0.5 rho <V^3>, in W per m2 facing it."""
    return 0.5 * density_kg_m3 * weibull_moment(k, c, 3.0)


@dataclass(frozen=True)
class Device:
    name: str
    power_coefficient: float        # share of the wind's power taken on the frontal area, up to rated speed
    cut_in_m_s: float
    rated_m_s: float
    cut_out_m_s: float
    thrust_coefficient: float       # on the frontal area, while working
    parked_drag_coefficient: float  # on the frontal area, stopped for a storm
    tip_speed_ratio: float = 0.0    # blade tip speed / wind speed; 0 for a device without a rotor
    blades: int = 0
    chord_share: float = 0.0        # blade chord / rotor radius
    pitch_deg: float = 0.0


def mean_output_w_m2(device: Device, k: float, c: float, density_kg_m3: float, screen_factor: float = 1.0) -> float:
    """Mean electrical-shaft output per m2 of the device's frontal area over a Weibull wind distribution."""
    v = np.linspace(0.0, max(3.0 * c, device.cut_out_m_s) + 5.0, 4001)
    pdf = (k / c) * (v / c) ** (k - 1.0) * np.exp(-(v / c) ** k)
    on = (v >= device.cut_in_m_s) & (v <= device.cut_out_m_s)
    power = np.where(on, 0.5 * density_kg_m3 * device.power_coefficient * np.minimum(v, device.rated_m_s) ** 3, 0.0)
    return float(np.trapezoid(power * pdf, v) * screen_factor)


def working_share(device: Device, k: float, c: float) -> float:
    """Share of the time the device turns (wind between cut-in and cut-out)."""
    cdf = lambda x: 1.0 - np.exp(-(x / c) ** k)
    return float(cdf(device.cut_out_m_s) - cdf(device.cut_in_m_s))


def band_collision_probability(device: Device, rotor_radius_m: float, wind_m_s: float, flier_length_m: float,
                               flier_span_m: float, flier_speed_m_s: float, flapping: bool = True,
                               points: int = 400) -> float:
    """Chance that a flier crossing the rotor disc once, at its airspeed plus the wind, is struck by a blade.

    Band (2012): at radius r the probability is (B omega / (2 pi v)) (|c sin(gamma) + alpha c cos(gamma)| + max(L,
    W alpha F)), with alpha = v / (r omega), B blades of chord c and pitch gamma, flier length L, span W, F = 1 for
    flapping flight, and v the flier's speed through the disc; the disc average weights radius by area. The flier
    here crosses downwind, so v is its airspeed plus the wind.
    """
    if device.blades == 0:
        return 0.0
    omega = device.tip_speed_ratio * wind_m_s / rotor_radius_m
    v = flier_speed_m_s + wind_m_s
    r = (np.arange(points) + 0.5) / points * rotor_radius_m
    chord = device.chord_share * rotor_radius_m
    alpha = v / (r * omega)
    gam = np.radians(device.pitch_deg)
    body = np.maximum(flier_length_m, flier_span_m * alpha * (1.0 if flapping else 2.0 / np.pi))
    p = device.blades * omega / (2.0 * np.pi * v) * (np.abs(chord * np.sin(gam) + alpha * chord * np.cos(gam)) + body)
    p = np.minimum(p, 1.0)
    return float((p * 2.0 * r).sum() / points / rotor_radius_m)
