"""Solar fluxes of a column: delta-scaled two-stream layers combined by adding.

Each layer is a homogeneous two-stream slab with the practical improved flux
method (PIFM) coefficients of Zdunkowski et al. (1980), delta-scaled for forward
scattering (Joseph et al. 1976), which keeps reflectance non-negative for
strongly absorbing layers. Layers are combined by the adding method over a
Lambertian surface. The direct beam is attenuated along straight rays through
spherical shells (pseudo-spherical), which matters at low Sun for the Moon's
deep atmosphere; each layer's effective cosine reproduces that attenuation, so
the layer solutions stay energy-conserving. Planetary means integrate over the
sunlit disk with Gauss nodes in mu0 weighted by 2 mu0.
"""
from __future__ import annotations
import numpy as np

from shared.constants import AU, SOLAR_CONSTANT, SUN_RADIUS
from atmosphere.radiative_convective import fetch_inputs
from atmosphere.radiative_convective.longwave import planck


def solar_spectrum(nu, solar_constant=SOLAR_CONSTANT, tail_temperature=5772.0, shield=None):
    """Top-of-atmosphere spectral irradiance (W m^-2 per cm^-1) and the flux below nu[0].

    TSIS-1 HSRS for 202-2729.9 nm; beyond 2729.9 nm a blackbody tail joined to it
    continuously; the whole spectrum is scaled so that its total, including the
    part longward of nu[0], equals the solar constant. Irradiance shortward of
    202 nm (about 0.1 W/m^2) is omitted. `shield`, if given, is a transmission
    function of wavelength in nm applied to the sunlight (its value at 5000 nm
    is used longward of the grid).
    """
    data = np.genfromtxt(fetch_inputs.path('TSIS1_HSRS_stride100.csv'), delimiter=',', skip_header=1)
    lam_nm, f_lam = data[:, 0], data[:, 1]
    wn = 1e7 / lam_nm[::-1]
    f_wn = (f_lam * lam_nm**2 / 1e7)[::-1]                # per cm^-1
    edge = wn[0]
    dilution = np.pi * (SUN_RADIUS / AU) ** 2

    def tail(x):
        return dilution * planck(x, tail_temperature)
    tail_scale = f_wn[0] / tail(edge)
    fine = np.linspace(1e-3, edge, 200001)
    below_edge = np.trapezoid(tail_scale * tail(fine), fine)
    hsrs_total = np.trapezoid(f_wn, wn)
    scale = solar_constant / (below_edge + hsrs_total)
    nu = np.asarray(nu, dtype=float)
    spectrum = np.where(nu < edge, tail_scale * tail(nu), np.interp(nu, wn, f_wn, right=0.0)) * scale
    lo = np.linspace(1e-3, nu[0], 20001)
    longward = scale * np.trapezoid(tail_scale * tail(lo), lo) if nu[0] < edge else None
    if shield is not None:
        spectrum = spectrum * shield(1e7 / nu)
        if longward is not None:
            longward = longward * float(shield(np.array([5000.0]))[0])
    return spectrum, longward


def layer_properties(tau, omega, g, mu):
    """Delta-scaled PIFM two-stream reflectance and transmittance of a homogeneous layer.

    Returns diffuse reflectance and transmittance, the direct-beam reflectance and
    diffuse transmittance (per unit direct flux on a horizontal surface) and the
    direct transmittance exp(-tau/mu).
    """
    f = g * g
    tau_s = (1 - omega * f) * tau
    om = np.minimum((1 - f) * omega / (1 - omega * f), 1 - 1e-9)
    gs = (g - f) / (1 - f)
    g1 = (8 - om * (5 + 3 * gs)) / 4
    g2 = 3 * om * (1 - gs) / 4
    g3 = (2 - 3 * gs * mu) / 4
    g4 = 1 - g3
    lam = np.sqrt(np.maximum(g1 * g1 - g2 * g2, 1e-30))
    # Avoid the resonance lam*mu = 1 of the particular solution.
    mu = np.where(np.abs(1 - (lam * mu) ** 2) < 1e-6, mu * (1 + 1e-4), mu)
    gam = g2 / (g1 + lam)
    e = np.exp(-lam * tau_s)
    eb = np.exp(-tau_s / mu)
    denom = 1 - gam * gam * e * e
    r_dif = gam * (1 - e * e) / denom
    t_dif = e * (1 - gam * gam) / denom
    det = 1 / mu**2 - lam * lam
    a, b = g1 + 1 / mu, g1 - 1 / mu
    up_p = om * (-b * g3 - g2 * g4) / det
    dn_p = -om * (a * g4 + g2 * g3) / det
    a1 = (gam * e * up_p * eb - dn_p) / denom
    a2 = -gam * e * a1 - up_p * eb
    r_dir = (gam * a1 + a2 * e + up_p) / mu
    t_dir = (a1 * e + gam * a2 + dn_p * eb) / mu
    return r_dif, t_dif, np.maximum(r_dir, 0.0), np.maximum(t_dir, 0.0), eb


def slant_factors(radius_levels, mu0):
    """Path-length factors P[k, j]: slant optical depth to level k = sum_j P[k, j] tau_j.

    Levels are ordered top-down; layer j lies between levels j and j+1. Straight
    rays reach each level along the same local vertical with zenith cosine mu0.
    """
    r = np.asarray(radius_levels, dtype=float)
    n = len(r) - 1
    thick = r[:-1] - r[1:]
    sin2 = 1 - mu0 * mu0
    p = np.zeros((n + 1, n))
    for k in range(1, n + 1):
        b2 = r[k] ** 2 * sin2
        top = np.sqrt(np.maximum(r[:k] ** 2 - b2, 0.0))
        bottom = np.sqrt(np.maximum(r[1:k + 1] ** 2 - b2, 0.0))
        # (top - bottom) written without cancellation.
        p[k, :k] = (r[:k] + r[1:k + 1]) / np.maximum(top + bottom, 1e-300)
    return p


def column_fluxes(tau_abs, tau_sca, g_sca, surface_albedo, mu0, radius_levels=None):
    """Spectral fluxes for unit incident flux normal to the beam.

    tau_abs, tau_sca, g_sca: (layers, nu), layers top-down. Returns downward
    direct, downward diffuse and upward diffuse fluxes at levels (levels, nu).
    """
    tau = tau_abs + tau_sca
    omega = np.where(tau > 0, tau_sca / np.maximum(tau, 1e-300), 0.0)
    nlay, nnu = tau.shape
    if radius_levels is None:
        mu_eff = np.full_like(tau, mu0)
    else:
        slant = slant_factors(radius_levels, mu0) @ tau
        slant = np.maximum.accumulate(slant, axis=0)
        dslant = np.diff(slant, axis=0)
        mu_eff = np.where(dslant > 0, tau / np.maximum(dslant, 1e-300), 1.0)
        mu_eff = np.clip(mu_eff, mu0, 1.0)
    r_dif, t_dif, r_dir, t_dir, tdd = layer_properties(tau, omega, g_sca, mu_eff)
    rb_dif = np.empty((nlay + 1, nnu)); rb_dir = np.empty((nlay + 1, nnu))
    rb_dif[nlay] = surface_albedo; rb_dir[nlay] = surface_albedo
    for k in range(nlay - 1, -1, -1):
        d = 1 - r_dif[k] * rb_dif[k + 1]
        rb_dif[k] = r_dif[k] + t_dif[k] ** 2 * rb_dif[k + 1] / d
        rb_dir[k] = r_dir[k] + t_dif[k] * (tdd[k] * rb_dir[k + 1] + t_dir[k] * rb_dif[k + 1]) / d
    direct = np.empty((nlay + 1, nnu)); diffuse = np.empty((nlay + 1, nnu)); up = np.empty((nlay + 1, nnu))
    direct[0] = mu0; diffuse[0] = 0.0
    up[0] = rb_dir[0] * mu0
    for k in range(nlay):
        direct[k + 1] = direct[k] * tdd[k]
        d = 1 - r_dif[k] * rb_dif[k + 1]
        diffuse[k + 1] = (t_dir[k] * direct[k] + t_dif[k] * diffuse[k] + r_dif[k] * rb_dir[k + 1] * direct[k + 1]) / d
        up[k + 1] = rb_dir[k + 1] * direct[k + 1] + rb_dif[k + 1] * diffuse[k + 1]
    return direct, diffuse, up


def disk_nodes(n=6):
    """Gauss nodes mu0 on (0, 1] and weights w (sum 1) for disk averages 2 * integral(f mu0 dmu0)."""
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


def planetary_mean(nu, spectrum, tau_abs, tau_sca, g_sca, surface_albedo, radius_levels=None,
                   n_mu0=6, chunk=40000):
    """Global-mean solar budget of the column (W/m^2), all levels top-down.

    Returns incident (S/4 over the spectral grid), reflected at the top, net
    downward flux at each level (so layer absorption is its decrease), and the
    spectral reflectance for diagnostics.
    """
    mu, w = disk_nodes(n_mu0)
    step = nu[1] - nu[0]
    weights = np.full(nu.size, step); weights[0] *= 0.5; weights[-1] *= 0.5
    nlev = tau_abs.shape[0] + 1
    net = np.zeros(nlev)
    reflected = 0.0
    spectral_up = np.zeros(nu.size)
    spectral_surface = np.zeros(nu.size)
    for a in range(0, nu.size, chunk):
        sl = slice(a, a + chunk)
        f = spectrum[sl] * weights[sl]
        for m, wt in zip(mu, w):
            direct, diffuse, up = column_fluxes(tau_abs[:, sl], tau_sca[:, sl], g_sca[:, sl],
                                                surface_albedo, m, radius_levels)
            # Planetary mean of a per-area flux F(mu0) = 0.5 * integral F dmu0 (unit F0 -> S).
            scale = 0.5 * wt
            net += scale * ((direct + diffuse - up) @ f)
            reflected += scale * float(up[0] @ f)
            spectral_up[sl] += scale * up[0] * spectrum[sl]
            spectral_surface[sl] += scale * (direct[-1] + diffuse[-1] - up[-1]) * spectrum[sl]
    incident = 0.25 * float(spectrum @ weights)
    return dict(incident=incident, reflected=reflected, net=net, spectral_reflected=spectral_up,
                spectral_surface=spectral_surface, weights=weights)
