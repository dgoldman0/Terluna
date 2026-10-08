"""Thunder at the ground: geometric acoustics through a horizontally uniform atmosphere with wind.

Rays leave points along a lightning channel in every direction of one vertical plane, an azimuth. In a layered
atmosphere the ray parameter a = cos(theta) / c_eff is conserved along a ray, theta its elevation and c_eff the sound
speed plus the wind along the azimuth (the effective-sound-speed approximation, good to the Mach number squared), so
in a layer where c_eff changes linearly with height a ray follows a circular arc of radius 1 / (a |dc_eff/dz|). A ray
turns where a c_eff reaches 1, reflects at a hard, flat ground and leaves through the top. Each carries the energy of
its share of the sphere around its source point, and the energy reaching the ground is gathered by range and arrival
time, which keeps the caustics of geometric acoustics finite.

Absorption follows ISO 9613-1 by octave band along each ray, its oxygen and nitrogen relaxation scaled to the air's
mole fractions of each. A flash's acoustic energy is a share of the electrostatic energy it releases (Holmes et al.
1971: 0.18 %), spread evenly along its channel, in a spectrum x^2 / (1 + x^4) of x = f / f_m, peaked at Few's (1969)
f_m = 0.63 c0 (p0 / E_l)^(1/2) for E_l the energy per metre of channel. At the ground the sound adds to its own
reflection: pressure doubles at 125 Hz and below, where the reflection arrives in phase at a listener's ear, and
energy doubles above.
"""
from __future__ import annotations

import numpy as np

R_GAS = 8.314462618
GAMMA = 1.4
P_REF = 2.0e-5                                  # Pa
ISO_P = 101325.0
ISO_T0 = 293.15
ISO_T01 = 273.16
EARTH_O2, EARTH_N2 = 0.209, 0.781               # the mole fractions ISO 9613-1's coefficients carry
BANDS = np.array([16.0, 31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0])          # octave centres, Hz
AUDIBLE = BANDS >= 31.5                         # the 16-Hz band is infrasound
THRESHOLD_DB = np.array([np.nan, 59.5, 37.5, 22.1, 11.4, 4.4, 2.4, -1.3, -5.4])            # ISO 226:2003, hearing
A_WEIGHT_DB = np.array([-56.7, -39.4, -26.2, -16.1, -8.6, -3.2, 0.0, 1.2, 1.0])
COHERENT_HZ = 125.0                             # pressure doubling at the ground up to this band
HOLMES_EFFICIENCY = 0.0018                      # acoustic share of a flash's energy (Holmes et al. 1971)


def sound_speed(t_k, qv=0.0, molar_mass=0.028964):
    """The speed of sound (m/s) in moist air of temperature t_k (K) and vapour mixing ratio qv (kg/kg)."""
    eps = 0.018015 / molar_mass
    tv = np.asarray(t_k) * (1.0 + np.asarray(qv) / eps) / (1.0 + np.asarray(qv))
    return np.sqrt(GAMMA * R_GAS * tv / molar_mass)


def vapour_percent(qv, molar_mass=0.028964):
    """Molar concentration of water vapour (%) for a mixing ratio qv (kg/kg)."""
    eps = 0.018015 / molar_mass
    return 100.0 * np.asarray(qv) / (np.asarray(qv) + eps)


def absorption_db_per_m(f, t_k, p_pa, h_percent, x_o2=EARTH_O2, x_n2=EARTH_N2):
    """Atmospheric absorption (dB/m) at frequency f (Hz) by ISO 9613-1, its relaxation terms scaled to the air's
    oxygen and nitrogen."""
    f = np.asarray(f, dtype=float)
    pr = np.asarray(p_pa) / ISO_P
    tr = np.asarray(t_k) / ISO_T0
    h = np.asarray(h_percent)
    fr_o = pr * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    fr_n = pr * tr ** -0.5 * (9.0 + 280.0 * h * np.exp(-4.170 * (tr ** (-1.0 / 3.0) - 1.0)))
    return 8.686 * f ** 2 * (1.84e-11 / pr * tr ** 0.5 + tr ** -2.5 * (
        0.01275 * (x_o2 / EARTH_O2) * np.exp(-2239.1 / np.asarray(t_k)) / (fr_o + f ** 2 / fr_o)
        + 0.1068 * (x_n2 / EARTH_N2) * np.exp(-3352.0 / np.asarray(t_k)) / (fr_n + f ** 2 / fr_n)))


def few_peak_hz(energy_per_m, p_pa, c_m_s):
    """Few's (1969) peak frequency of thunder from a channel of energy_per_m (J/m) in air at p_pa and c_m_s."""
    return 0.63 * c_m_s * np.sqrt(p_pa / energy_per_m)


def band_shares(f_peak):
    """The share of the acoustic energy in each octave band for a spectrum x^2 / (1 + x^4), x = f / f_peak, with the
    16-Hz band taking everything below 22 Hz and the top band everything above it."""
    x = np.geomspace(1.0e-4, 1.0e4, 20001)
    s = x ** 2 / (1.0 + x ** 4)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(x))])
    total = np.pi / (2.0 * np.sqrt(2.0))
    edges = np.concatenate([[0.0], BANDS[:-1] * np.sqrt(2.0), [np.inf]]) / f_peak
    at = np.interp(np.clip(edges, x[0], x[-1]), x, cum)
    at[0], at[-1] = 0.0, total
    return np.diff(at) / total


class Atmosphere:
    """A layered atmosphere on edges z (m, the ground at 0): temperature (K), pressure (Pa), vapour mixing ratio and
    wind (u eastward, v northward, m/s) at the edges, and the air's molar mass and oxygen and nitrogen fractions."""

    def __init__(self, z, t_k, p_pa, qv, u, v, molar_mass=0.028964, x_o2=EARTH_O2, x_n2=EARTH_N2):
        self.z = np.asarray(z, dtype=float)
        self.t, self.p, self.qv = (np.asarray(a, dtype=float) for a in (t_k, p_pa, qv))
        self.u, self.v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
        self.molar_mass, self.x_o2, self.x_n2 = molar_mass, x_o2, x_n2
        self.c = sound_speed(self.t, self.qv, molar_mass)
        mid = lambda a: 0.5 * (a[1:] + a[:-1])
        h = vapour_percent(mid(self.qv), molar_mass)
        self.alpha = np.array([absorption_db_per_m(f, mid(self.t), mid(self.p), h, x_o2, x_n2) for f in BANDS]).T
        self.rho_c_ground = self.p[0] * molar_mass / (R_GAS * self.t[0]) * self.c[0]

    def c_eff(self, azimuth_deg):
        a = np.radians(azimuth_deg)
        return self.c + self.u * np.sin(a) + self.v * np.cos(a)          # azimuth clockwise from north


def trace(atm: Atmosphere, z_source, azimuth_deg=90.0, n_rays=3601, x_max=500.0e3, max_bounces=8):
    """Rays from a point at height z_source (on one of the atmosphere's edges) launched every 180/(n_rays-1) degrees
    of elevation in one azimuth. Returns the ground hits: range (m), arrival time (s), absorption by band (dB), the sine
    of the grazing angle, the launch weight cos(theta0) dtheta0 and the bounce number."""
    z, c = atm.z, atm.c_eff(azimuth_deg)
    k0 = int(np.argmin(abs(z - z_source)))
    if abs(z[k0] - z_source) > 1.0:
        raise ValueError('the source must lie on an edge of the atmosphere')
    top = len(z) - 1
    cmax_above = np.maximum.accumulate(c[::-1])[::-1]
    grad = np.diff(c) / np.diff(z)
    theta = np.linspace(-90.0, 90.0, n_rays)[1:-1]
    theta = theta[theta != 0.0]
    dtheta = np.radians(180.0 / (n_rays - 1))
    th = np.radians(theta)
    a = np.cos(th) / c[k0]
    k = np.full(th.size, k0)
    d = np.where(th > 0.0, 1, -1)
    x = np.zeros(th.size)
    t = np.zeros(th.size)
    absorb = np.zeros((th.size, BANDS.size))
    bounces = np.zeros(th.size, dtype=int)
    alive = np.ones(th.size, dtype=bool)
    weight = np.cos(th) * dtheta
    hits = {key: [] for key in ('x', 't', 'absorb', 'sin', 'weight', 'bounce')}
    while alive.any():
        ground = alive & (d < 0) & (k == 0)
        if ground.any():
            g = np.nonzero(ground)[0]
            hits['x'].append(x[g]); hits['t'].append(t[g]); hits['absorb'].append(absorb[g].copy())
            hits['sin'].append(np.sqrt(np.clip(1.0 - (a[g] * c[0]) ** 2, 0.0, 1.0)))
            hits['weight'].append(weight[g]); hits['bounce'].append(bounces[g].copy())
            d[g] = 1
            bounces[g] += 1
            alive[g] &= bounces[g] < max_bounces
        escape = alive & (d > 0) & ((k >= top) | (a * cmax_above[np.minimum(k, top)] < 1.0))
        alive &= ~escape
        move = np.nonzero(alive & ~ground)[0]
        if move.size == 0:
            continue
        kk, dd, aa = k[move], d[move], a[move]
        kn = kk + dd
        layer = np.minimum(kk, kn)
        gl = grad[layer] * dd                                  # the change of c_eff along the ray's way, per metre
        u1 = np.clip(aa * c[kk], 0.0, 1.0)
        u2 = aa * c[kn]
        s1 = np.sqrt(1.0 - u1 ** 2)
        dz = abs(z[kn] - z[kk])
        dx = np.empty(move.size); dt = np.empty(move.size); ds = np.empty(move.size)
        cross = u2 < 1.0
        flat = np.abs(gl) < 1.0e-9
        # crossing a layer with a straight ray
        m = cross & flat
        dx[m] = dz[m] * u1[m] / s1[m]; ds[m] = dz[m] / s1[m]; dt[m] = ds[m] * aa[m] / u1[m]
        # crossing along an arc
        m = cross & ~flat
        s2 = np.sqrt(np.clip(1.0 - u2[m] ** 2, 0.0, 1.0))
        dx[m] = np.abs(s1[m] - s2) / (aa[m] * np.abs(gl[m]))
        dt[m] = np.abs(np.log(u2[m] * (1.0 + s1[m]) / (u1[m] * (1.0 + s2))) / gl[m])
        ds[m] = np.abs(np.arccos(u1[m]) - np.arccos(np.clip(u2[m], 0.0, 1.0))) / (aa[m] * np.abs(gl[m]))
        # turning inside the layer: up to the turning height and back
        m = ~cross
        dx[m] = 2.0 * s1[m] / (aa[m] * np.abs(gl[m]))
        dt[m] = 2.0 * np.abs(np.log((1.0 + s1[m]) / u1[m]) / gl[m])
        ds[m] = 2.0 * np.arccos(u1[m]) / (aa[m] * np.abs(gl[m]))
        x[move] += dx
        t[move] += dt
        absorb[move] += atm.alpha[layer] * ds[:, None]
        k[move] = np.where(cross, kn, kk)
        d[move] = np.where(cross, dd, -dd)
        alive[move] &= x[move] < x_max
    return {key: (np.concatenate(v) if v else np.array([])) for key, v in hits.items()}


def ground_levels(hits_by_source, energy_by_source, f_peak, atm: Atmosphere, x_edges, time_bin=1.0, min_sin=0.035):
    """Sound at the ground from point sources of acoustic energy (J) whose hits trace() gave, at the range bins
    x_edges: the sound exposure (Pa^2 s) and the loudest one-second level (dB) by band, and the arrival of the first
    sound and the span holding 90 % of the energy. Grazing angles below asin(min_sin), about 2 degrees, count as 2
    degrees, where geometric acoustics fails."""
    shares = band_shares(f_peak)
    nx = len(x_edges) - 1
    centres = 0.5 * (x_edges[1:] + x_edges[:-1])
    area = 0.5 * (x_edges[1:] ** 2 - x_edges[:-1] ** 2)       # per radian of azimuth
    doubling = np.where(BANDS <= COHERENT_HZ, 4.0, 2.0)
    t_all = np.concatenate([h['t'] for h in hits_by_source if h['t'].size]) if hits_by_source else np.array([0.0])
    t_edges = np.arange(0.0, t_all.max() + 2.0 * time_bin, time_bin)
    exposure = np.zeros((nx, t_edges.size - 1, BANDS.size))
    for hits, energy in zip(hits_by_source, energy_by_source):
        if not hits['x'].size:
            continue
        ix = np.searchsorted(x_edges, hits['x'], side='right') - 1
        ok = (ix >= 0) & (ix < nx)
        it = np.clip(np.searchsorted(t_edges, hits['t'], side='right') - 1, 0, t_edges.size - 2)
        flux = (energy * shares[None, :] * hits['weight'][:, None] / (4.0 * np.pi) * 10.0 ** (-hits['absorb'] / 10.0)
                / np.maximum(hits['sin'], min_sin)[:, None])
        np.add.at(exposure, (ix[ok], it[ok]), flux[ok] / area[ix[ok], None])
    exposure *= atm.rho_c_ground * doubling[None, None, :]
    total = exposure.sum(axis=1)
    with np.errstate(divide='ignore'):
        peak = 10.0 * np.log10(exposure.max(axis=1) / (time_bin * P_REF ** 2))
        a_series = (exposure[:, :, AUDIBLE] * 10.0 ** (A_WEIGHT_DB[AUDIBLE] / 10.0)).sum(axis=2)
        peak_a = 10.0 * np.log10(a_series.max(axis=1) / (time_bin * P_REF ** 2))
    energy_t = exposure[:, :, AUDIBLE].sum(axis=2)
    cum = np.cumsum(energy_t, axis=1)
    first = np.full(nx, np.nan); span = np.full(nx, np.nan)
    for i in range(nx):
        if cum[i, -1] > 0.0:
            lo = np.searchsorted(cum[i], 0.05 * cum[i, -1]); hi = np.searchsorted(cum[i], 0.95 * cum[i, -1])
            first[i] = t_edges[np.argmax(energy_t[i] > 0.0)]
            span[i] = (hi - lo + 1) * time_bin
    return dict(x_m=centres, exposure_pa2s=total, peak_db=peak, peak_dba=peak_a, first_s=first, span_s=span)


def audible_range(levels, quiet_dba=None):
    """The farthest range bin where some audible band's loudest second reaches the threshold of hearing, or, given
    quiet_dba, where the loudest second's A-weighted level reaches that background."""
    if quiet_dba is None:
        ok = np.any(levels['peak_db'][:, AUDIBLE] >= THRESHOLD_DB[AUDIBLE], axis=1)
    else:
        ok = levels['peak_dba'] >= quiet_dba
    return float(levels['x_m'][ok].max()) if ok.any() else 0.0
