"""Leaf reflectance and transmittance from PROSPECT-D (Feret et al. 2017).

A leaf is a pile of N absorbing plates (Jacquemoud and Baret 1990). The first
surface receives light within alpha = 40 degrees of its normal and the inner
plates receive isotropic light. A plate's transmissivity follows from the
specific absorption of its chlorophyll, carotenoids, anthocyanins, brown pigments,
water and dry matter (Allen et al. 1969); the interface transmissivities are the
averages over incidence angles of Stern (1964) and Allen (1973); the pile is
combined with Stokes' (1862) equations. The coefficients are the published
PROSPECT-D table, restored by fetch_inputs.py.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import exp1

from biosphere.canopy import fetch_inputs

WAVELENGTH_NM = np.arange(400, 2501)
CONSTITUENTS = ('chlorophyll', 'carotenoids', 'anthocyanins', 'brown', 'water', 'dry_matter')
FIRST_SURFACE_DEG = 40.0


@dataclass(frozen=True)
class Leaf:
    """PROSPECT-D traits. The defaults are a typical green broadleaf leaf."""
    structure: float = 1.5        # N, plates
    chlorophyll: float = 40.0     # Cab, ug/cm2
    carotenoids: float = 8.0      # Car, ug/cm2
    anthocyanins: float = 0.0     # Ant, ug/cm2
    brown: float = 0.0            # Cbrown, arbitrary units
    water: float = 0.01           # Cw, cm
    dry_matter: float = 0.009     # Cm, g/cm2


def coefficients():
    """Refractive index ('nr') and specific absorption coefficients at WAVELENGTH_NM."""
    data = np.loadtxt(fetch_inputs.path('prospect_d_spectra.txt'), comments='#')
    if data.shape != (WAVELENGTH_NM.size, 8) or not np.array_equal(data[:, 0], WAVELENGTH_NM):
        raise ValueError('Unexpected PROSPECT-D table')
    return dict(nr=data[:, 1], **{name: data[:, i + 2] for i, name in enumerate(CONSTITUENTS)})


def interface_transmissivity(alpha_deg, n):
    """Mean transmissivity of a plane dielectric interface of index n for isotropic light within
    alpha_deg of the normal (Stern 1964; Allen 1973)."""
    n2 = n * n
    n_plus, n_minus = n2 + 1, n2 - 1
    a = (n + 1) ** 2 / 2
    k = -n_minus ** 2 / 4
    s2 = np.sin(np.radians(alpha_deg)) ** 2
    root = np.sqrt((s2 - n_plus / 2) ** 2 + k) if alpha_deg != 90 else 0.0
    b = root - (s2 - n_plus / 2)
    ts = (k ** 2 / (6 * b ** 3) + k / b - b / 2) - (k ** 2 / (6 * a ** 3) + k / a - a / 2)
    tp = (-2 * n2 * (b - a) / n_plus ** 2
          - 2 * n2 * n_plus * np.log(b / a) / n_minus ** 2
          + n2 * (1 / b - 1 / a) / 2
          + 16 * n2 ** 2 * (n2 ** 2 + 1) * np.log((2 * n_plus * b - n_minus ** 2) / (2 * n_plus * a - n_minus ** 2))
          / (n_plus ** 3 * n_minus ** 2)
          + 16 * n2 ** 3 * (1 / (2 * n_plus * b - n_minus ** 2) - 1 / (2 * n_plus * a - n_minus ** 2)) / n_plus ** 3)
    return (ts + tp) / (2 * s2)


def plate_transmissivity(k):
    """Transmissivity of an absorbing plate to isotropic light: (1 - k) exp(-k) + k^2 E1(k)."""
    k = np.asarray(k, dtype=float)
    out = np.ones_like(k)
    pos = k > 0
    out[pos] = (1 - k[pos]) * np.exp(-k[pos]) + k[pos] ** 2 * exp1(k[pos])
    return out


def optics(leaf: Leaf = Leaf(), coeffs=None):
    """Reflectance and transmittance of the leaf at WAVELENGTH_NM."""
    c = coefficients() if coeffs is None else coeffs
    k = sum(getattr(leaf, name) * c[name] for name in CONSTITUENTS) / leaf.structure
    tau = plate_transmissivity(k)
    n = c['nr']
    t_first = interface_transmissivity(FIRST_SURFACE_DEG, n)
    t_in = interface_transmissivity(90.0, n)
    t_out = t_in / n ** 2
    r_out = 1 - t_out
    d = 1 - (r_out * tau) ** 2
    # The first plate, lit within 40 degrees of its normal, and an inner plate, lit isotropically.
    t_a = t_first * tau * t_out / d
    r_a = (1 - t_first) + r_out * tau * t_a
    t = t_in * tau * t_out / d
    r = (1 - t_in) + r_out * tau * t
    # The remaining N - 1 plates (Stokes), with the lossless limit taken separately.
    m = leaf.structure - 1
    lossless = r + t >= 1.0
    rr, tt = np.where(lossless, 0.5, r), np.where(lossless, 0.25, t)
    disc = np.sqrt(np.maximum((1 + rr + tt) * (1 + rr - tt) * (1 - rr + tt) * (1 - rr - tt), 0.0))
    a = (1 + rr ** 2 - tt ** 2 + disc) / (2 * rr)
    b = (1 - rr ** 2 + tt ** 2 + disc) / (2 * tt)
    bm = b ** m
    denom = a ** 2 * bm ** 2 - 1
    r_sub = a * (bm ** 2 - 1) / denom
    t_sub = bm * (a ** 2 - 1) / denom
    t_sub = np.where(lossless, t / (t + (1 - t) * m), t_sub)
    r_sub = np.where(lossless, 1 - t_sub, r_sub)
    join = 1 - r_sub * r
    return r_a + t_a * r_sub * t / join, t_a * t_sub / join
