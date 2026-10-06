"""Full-range equilibrium water spectra; separate from pinned short-wave studies."""
from __future__ import annotations
import numpy as np
from illumination.water_surface.short_waves import curvature, spreading, _peak

def full_curvature(k, u10, u_star, omega, gravity):
    """Full-range spectrum with the documented correction to ECKV equation 41.

    Mobley's Light & Water supplement 6 gives F_m = L_PM J_p exp(...).
    The short-wave asymptote in curvature() omits that low-k suppression and
    must not be integrated down toward k=0 as an elevation spectrum.
    https://misclab.umeoce.maine.edu/education/Light&Water/D:/PAPER/supp6.pdf
    Existing high-wavenumber-only callers retain their published baseline.
    """
    k = np.asarray(k, float)
    if np.any(k <= 0) or not np.isfinite(k).all():
        raise ValueError('Positive finite wavenumbers required')
    long, short = curvature(k, u10, u_star, omega, gravity)
    kp, _ = _peak(u10, omega, gravity)
    gamma = 1.7 if omega <= 1 else 1.7 + 6 * np.log10(omega)
    sigma = .08 * (1 + 4 * omega ** -3)
    peak = np.exp(-(np.sqrt(k/kp)-1)**2 / (2*sigma*sigma))
    return long, short * np.exp(-1.25*(kp/k)**2) * gamma**peak


def unified_density(u10, u_star, omega, toward_deg, gravity):
    """Conditional equilibrium ECKV density across gravity and capillary waves.

    Uses the documented full-range short-wave correction. Wave age and wind
    are inputs; this is not a fetch-, duration- or terrain-resolving wave solver.
    """
    wind = np.radians(toward_deg)
    def density(kx, ky):
        k = np.hypot(kx, ky)
        out = np.zeros_like(k)
        valid = k > 0
        kk = k[valid]
        long, short = full_curvature(kk, u10, u_star, omega, gravity)
        delta = spreading(kk, u10, u_star, omega, gravity)
        phi = np.arctan2(ky[valid], kx[valid])
        out[valid] = (long+short)/kk**4 * (1+delta*np.cos(2*(phi-wind)))/(2*np.pi)
        return out
    return density
