"""Rain as numbers of drops: size distributions at the cloud base, what reaches the ground, and arrival fluxes.

A cloud-base spectrum is a boundary condition chosen for a scenario: a gamma distribution of drop concentration,
normalized to a water flux entering the column. drops.fall follows each size class to the ground; a surviving class
keeps its number flux J (drops per m2 per second) while it shrinks, so its concentration below is J / u at its speed u.
The ledger reports the surface rain, the counts and the sizes and speeds weighted by number and by water separately.
Drops here fall independently; the encounter diagnostic measures, along their paths, how often they would meet. The
arrival flux carries a fixed rain rate and a fixed share of impacts by size, which fix the arrival rate whatever the
fall speed; speed sets the concentration needed to supply it.

    from atmosphere.rain import distributions
    d, width = distributions.diameter_grid(0.1, 6.0, 64)
    spectrum = distributions.gamma_flux(d, width, speeds, dm_mm=1.5, shape_mu=3, cloud_base_rain_mm_h=5)

Written for the Open Moon's rain on 11 October 2026 (the cohort ledger and encounter diagnostic) and 10-11 October
(the arrival flux), and brought into the repository on 17 October. Sizes are diameters in metres; rain is in mm/h.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MM_H_PER_M_S = 3.6e6


def _positive(value, name, zero=False):
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'non-negative' if zero else 'positive'}")
    return value


def _vector(value, name, zero=False):
    value = np.asarray(value, dtype=float)
    if value.ndim != 1 or not value.size or not np.all(np.isfinite(value)):
        raise ValueError(f'{name} must be a non-empty finite vector')
    if np.any(value < 0 if zero else value <= 0):
        raise ValueError(f'{name} has an entry of the wrong sign')
    return value


def drop_volume(diameter_m):
    """Volume (m3) of a drop of each volume-equivalent diameter."""
    return np.pi / 6.0 * _vector(diameter_m, 'diameter_m')**3


def diameter_grid(min_mm=0.1, max_mm=8.0, bins=64):
    """Logarithmic size classes: each class's geometric-mean diameter and its width (m)."""
    lo, hi = _positive(min_mm, 'min_mm'), _positive(max_mm, 'max_mm')
    if hi <= lo or isinstance(bins, bool) or int(bins) != bins or bins < 2:
        raise ValueError('Need max_mm > min_mm and a whole number of bins, at least 2')
    edges = np.geomspace(lo * 1e-3, hi * 1e-3, int(bins) + 1)
    return np.sqrt(edges[:-1] * edges[1:]), np.diff(edges)


def gamma_flux(diameter_m, width_m, downward_speed_m_s, *, dm_mm=1.5, shape_mu=0.0, cloud_base_rain_mm_h=10.0):
    """A gamma concentration spectrum n(D) ~ D^mu exp(-(mu + 4) D / Dm), normalized to the cloud-base water flux.

    dm_mm is the mass-weighted mean diameter M4/M3 of the untruncated distribution; the truncated grid's own value is
    returned beside it. flux_m2_s is each class's number flux, already integrated over its width. Every class must
    enter the column downward (speed over the ground positive); classes held up by an updraft need another boundary.
    """
    d, w, v = (_vector(x, n) for x, n in ((diameter_m, 'diameter_m'), (width_m, 'width_m'),
                                           (downward_speed_m_s, 'downward_speed_m_s')))
    if not d.shape == w.shape == v.shape:
        raise ValueError('Diameter, width and speed must have one entry per class')
    dm = _positive(dm_mm, 'dm_mm') * 1e-3
    mu = float(shape_mu)
    if not np.isfinite(mu) or mu <= -1:
        raise ValueError('shape_mu must be finite and above -1')
    q = _positive(cloud_base_rain_mm_h, 'cloud_base_rain_mm_h', zero=True) / MM_H_PER_M_S
    log_shape = mu * np.log(d / dm) - (mu + 4) * d / dm
    shape = np.exp(log_shape - np.max(log_shape))
    unscaled = shape * w * v
    scale = q / np.sum(np.pi * d**3 / 6 * unscaled)
    concentration = shape * w * scale
    realized = np.sum(concentration * d**4) / np.sum(concentration * d**3) if q else None
    return dict(diameter_m=d, width_m=w, concentration_m3=concentration, flux_m2_s=unscaled * scale,
                cloud_base_rain_mm_h=q * MM_H_PER_M_S, shape_mu=mu, untruncated_dm_mm=dm * 1e3,
                realized_dm_mm=None if realized is None else realized * 1e3, lambda_per_mm=(mu + 4) / (dm * 1e3))


def _weighted_quantiles(value, weight):
    """p10, p50, p90 and p99, interpolated through each class's weighted mid-rank: an approximation to the
    continuous distribution, which refining the classes sharpens."""
    order = np.argsort(value)
    x, p = value[order], weight[order]
    keep = p > 0
    x, p = x[keep], p[keep]
    if not p.size:
        return None
    cdf = (np.cumsum(p) - 0.5 * p) / p.sum()
    return dict(zip(('p10', 'p50', 'p90', 'p99'), map(float, np.interp([0.1, 0.5, 0.9, 0.99], cdf, x))))


def surface_summary(initial_diameter_m, cloud_base_flux_m2_s, surface_diameter_m, surface_speed_m_s, survived, *,
                    water_density_kg_m3=998.2):
    """The ledger at the ground: water, counts, concentration, and sizes and speeds by number and by water.

    Classes that did not survive may carry None or NaN; survivors need a positive size and downward speed. The
    momentum and kinetic-energy fluxes are what the rain brings to a horizontal surface, before any impact physics.
    """
    d0 = _vector(initial_diameter_m, 'initial_diameter_m')
    j0 = _vector(cloud_base_flux_m2_s, 'cloud_base_flux_m2_s', zero=True)
    d = np.asarray(surface_diameter_m, dtype=float)
    v = np.asarray(surface_speed_m_s, dtype=float)
    alive = np.asarray(survived, dtype=bool)
    if not d0.shape == j0.shape == d.shape == v.shape == alive.shape:
        raise ValueError('Every class array must have the same shape')
    if np.any(~np.isfinite(d[alive]) | (d[alive] <= 0) | ~np.isfinite(v[alive]) | (v[alive] <= 0)):
        raise ValueError('Surviving classes need a positive finite size and downward speed')
    rho = _positive(water_density_kg_m3, 'water_density_kg_m3')
    q0 = float(np.sum(j0 * np.pi * d0**3 / 6))
    d, v, j = d[alive], v[alive], j0[alive]
    water = np.pi * d**3 / 6 * j
    n = j / v
    held = n * np.pi * d**3 / 6
    q = float(water.sum())
    mm = d * 1e3
    return dict(
        cloud_base_rain_mm_h=q0 * MM_H_PER_M_S, surface_rain_mm_h=q * MM_H_PER_M_S,
        water_survival_fraction=q / q0 if q0 else None,
        number_survival_fraction=float(j.sum() / j0.sum()) if j0.sum() else None,
        surface_arrivals_m2_s=float(j.sum()), surface_concentration_m3=float(n.sum()),
        surface_liquid_water_g_m3=float(held.sum() * rho * 1e3),
        arrival_diameter_mm_by_number=_weighted_quantiles(mm, j),
        arrival_diameter_mm_by_water=_weighted_quantiles(mm, water),
        concentration_diameter_mm_by_number=_weighted_quantiles(mm, n),
        concentration_dm_mm=float(np.sum(mm * held) / held.sum()) if held.sum() else None,
        arrival_speed_m_s_by_number=_weighted_quantiles(v, j),
        arrival_speed_m_s_by_water=_weighted_quantiles(v, water),
        momentum_flux_pa=float(np.sum(rho * water * v)),
        kinetic_energy_flux_w_m2=float(np.sum(0.5 * rho * water * v**2)))


def encounters(altitude_m, diameter_m, downward_speed_m_s, terminal_speed_m_s, cloud_base_flux_m2_s, *,
               collision_efficiency=1.0):
    """How often each class would meet others while falling, along the unchanged paths of independent drops.

    altitude_m has shape (Z,), the other arrays (Z, N), the flux (N,); NaN or None marks a class no longer present.
    The kernel is E pi ((Di + Dj) / 2)^2 |ui - uj| on the actual speeds (a common air speed cancels), each class's
    concentration J / u, and the exposure tau_i = integral of sum_j K_ij n_j dz / u_i. E = 1 is the straight-path
    geometric reference; turbulence, wakes and charge change it, and it is no upper bound. The volume exposure weighs
    each encounter by the other drop's volume over the collector's. Both are warnings: a large value means
    collisions would reshape the spectrum that produced it, so the independent-drop spectrum stops being a prediction.
    """
    z = _vector(altitude_m, 'altitude_m', zero=True)
    if z.size < 2 or not (np.all(np.diff(z) > 0) or np.all(np.diff(z) < 0)):
        raise ValueError('Need at least two strictly monotone heights')
    d, u, vt = (np.asarray(x, dtype=float) for x in (diameter_m, downward_speed_m_s, terminal_speed_m_s))
    flux = _vector(cloud_base_flux_m2_s, 'cloud_base_flux_m2_s', zero=True)
    if not d.shape == u.shape == vt.shape == (z.size, flux.size):
        raise ValueError('Path arrays must have shape (heights, classes)')
    efficiency = _positive(collision_efficiency, 'collision_efficiency', zero=True)
    if efficiency > 1:
        raise ValueError('The collision efficiency lies between 0 and 1')
    valid = np.isfinite(d) & (d > 0) & np.isfinite(u) & (u > 0) & np.isfinite(vt) & (vt >= 0)
    rate = np.zeros_like(d)
    volume_rate = np.zeros_like(d)
    deviation = []
    for k in range(z.size):
        ids = np.flatnonzero(valid[k])
        size, motion, settling = d[k, ids], u[k, ids], vt[k, ids]
        relative = np.abs(motion[:, None] - motion[None, :])
        if ids.size:
            deviation.append(float(np.max(np.abs(relative - np.abs(settling[:, None] - settling[None, :])))))
        kernel = np.pi * ((size[:, None] + size[None, :]) / 2)**2 * relative * efficiency
        n = flux[ids] / motion
        rate[k, ids] = kernel @ n / motion
        volume = np.pi * size**3 / 6
        volume_rate[k, ids] = kernel @ (n * volume) / volume / motion
    dz = np.abs(np.diff(z))[:, None]
    tau = np.sum((rate[:-1] + rate[1:]) * 0.5 * dz, axis=0)
    volume_tau = np.sum((volume_rate[:-1] + volume_rate[1:]) * 0.5 * dz, axis=0)
    # The last partial step of a class that evaporates is a trapezoid; refine the heights near evaporation fronts.
    ground = int(np.argmin(z))
    survives = valid[ground]
    water = np.where(survives, np.pi * np.nan_to_num(d[ground])**3 / 6 * flux, 0.0)
    arriving = np.where(survives, flux, 0.0)

    def mean(weights):
        return float(np.sum(tau * weights) / np.sum(weights)) if np.sum(weights) else None

    share = (lambda mask: float(water[mask].sum() / water.sum()) if water.sum() else None)
    return dict(tau=tau.tolist(), probability_of_an_encounter=(-np.expm1(-tau)).tolist(),
                volume_exposure=volume_tau.tolist(),
                volume_exposure_by_surface_water=float(np.sum(volume_tau * water) / water.sum()) if water.sum()
                else None,
                surface_water_share_volume_exposure_above_0_1=share(volume_tau > 0.1),
                surface_water_share_volume_exposure_above_1=share(volume_tau > 1),
                largest_gap_from_terminal_relative_speed_m_s=max(deviation, default=0.0),
                tau_by_cloud_base_number=mean(flux), tau_by_surface_number=mean(arriving),
                tau_by_surface_water=mean(water), surface_water_share_tau_above_0_1=share(tau > 0.1),
                surface_water_share_tau_above_1=share(tau > 1), collision_efficiency=efficiency)


@dataclass(frozen=True)
class ArrivalFlux:
    """Number flux (drops per m2 per second) by size class across a receiving plane of area_m2.

    For vertical rain on sloping ground, take the horizontal catch area with the vertical speed, or the sloping area
    with the speed normal to it, but never project twice. The stored arrays are read-only copies.
    """
    diameter_m: np.ndarray
    number_flux_m2_s: np.ndarray
    area_m2: float

    def __post_init__(self):
        diameter = _vector(self.diameter_m, 'diameter_m').copy()
        flux = _vector(self.number_flux_m2_s, 'number_flux_m2_s', zero=True).copy()
        if diameter.shape != flux.shape:
            raise ValueError('Diameter and number flux must have one entry per class')
        diameter.flags.writeable = False
        flux.flags.writeable = False
        object.__setattr__(self, 'diameter_m', diameter)
        object.__setattr__(self, 'number_flux_m2_s', flux)
        object.__setattr__(self, 'area_m2', _positive(self.area_m2, 'area_m2', zero=True))

    @property
    def rain_mm_h(self):
        return float(np.dot(drop_volume(self.diameter_m), self.number_flux_m2_s)) * MM_H_PER_M_S

    @property
    def arrivals_per_second(self):
        return float(self.area_m2 * self.number_flux_m2_s.sum())

    @property
    def water_m3_per_second(self):
        return self.area_m2 * self.rain_mm_h / MM_H_PER_M_S

    def concentration_m3(self, normal_speed_m_s):
        """Drops per m3 in each class, J / v, at the speeds normal to the plane; a class with no flux and no speed
        holds none, and a class with flux needs a speed."""
        speed = _vector(normal_speed_m_s, 'normal_speed_m_s', zero=True)
        if speed.shape != self.number_flux_m2_s.shape:
            raise ValueError('Need one normal speed per class')
        if np.any((speed == 0) & (self.number_flux_m2_s > 0)):
            raise ValueError('A class arriving at all needs a positive normal speed')
        return np.divide(self.number_flux_m2_s, speed, out=np.zeros_like(speed), where=speed > 0)


def flux_from_impact_shares(*, rain_mm_h, diameter_m, impact_number_share, area_m2):
    """An arrival flux from a rain rate and the share of impacts in each class, which must already sum to one."""
    rain = _positive(rain_mm_h, 'rain_mm_h', zero=True)
    volume = drop_volume(diameter_m)
    share = _vector(impact_number_share, 'impact_number_share', zero=True)
    if volume.shape != share.shape:
        raise ValueError('Need one impact share per class')
    if not np.isclose(share.sum(), 1.0, rtol=0, atol=1e-12):
        raise ValueError('The impact shares must sum to one')
    total = rain / MM_H_PER_M_S / np.dot(share, volume)
    return ArrivalFlux(np.asarray(diameter_m, dtype=float), total * share, area_m2)


def flux_from_concentration(*, diameter_m, concentration_m3, normal_speed_m_s, area_m2):
    """An arrival flux J = n v from the drops in the air and their speeds normal to the plane."""
    diameter = _vector(diameter_m, 'diameter_m')
    concentration = _vector(concentration_m3, 'concentration_m3', zero=True)
    speed = _vector(normal_speed_m_s, 'normal_speed_m_s', zero=True)
    if not diameter.shape == concentration.shape == speed.shape:
        raise ValueError('Every class array must have the same shape')
    return ArrivalFlux(diameter, concentration * speed, area_m2)
