"""Light inside a plant canopy, by wavelength: sunlit and shaded leaves, layer by layer.

A horizontally uniform canopy of leaf area index L is cut into thin layers. Leaves
are spread with clumping index Omega and a spherical angle distribution (mean
projection G = 1/2), and scatter bi-Lambertianly with reflectance rho and
transmittance tau (Ross 1981). The direct beam falls off as exp(-G Omega L / mu0).
Diffuse light travels in streams at Gauss-Legendre cosines of the zenith angle,
downward and upward; a layer absorbs 1 - rho - tau of what it intercepts and
scatters the rest, a share beta backward and the rest forward, isotropically within
each hemisphere, with omega beta = [rho + tau + (rho - tau) cos^2 theta] / 2 and
cos^2 theta = 1/4 for spherical leaves (Sellers 1985, as in CLM5). The soil is
Lambertian. Scattering is iterated until no layer's scattered light changes by more
than 1e-11 of the incident light.

The fields are linear in the incident light, so solve() returns them for unit direct
light at one Sun height or for unit isotropic diffuse light; any incident spectrum
combines the two.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.signal import lfilter

G = 0.5
COS2_THETA = 0.25


@dataclass(frozen=True)
class Canopy:
    lai: float = 5.0
    clumping: float = 1.0
    layer_lai: float = 0.025
    streams: int = 8

    def layers(self):
        n = max(1, int(math.ceil(self.lai / self.layer_lai - 1e-9)))
        return n, self.lai / n

    def depth(self):
        """Cumulative leaf area above the middle of each layer."""
        n, dl = self.layers()
        return (np.arange(n) + 0.5) * dl


def gauss_streams(m):
    """Cosines and weights (summing to 1) of m streams on (0, 1)."""
    x, w = np.polynomial.legendre.leggauss(m)
    return 0.5 * (x + 1), 0.5 * w


def backscatter(rho, tau):
    """Backward share of the light a spherical-leaf layer scatters (Sellers 1985)."""
    omega = rho + tau
    return np.where(omega > 0, 0.5 * (omega + (rho - tau) * COS2_THETA) / np.maximum(omega, 1e-300), 0.5)


def solve(canopy: Canopy, rho, tau, soil, mu0=None, tol=1e-11, max_iter=20000):
    """Response of the canopy to unit incident light.

    mu0 given: unit horizontal flux of direct sunlight at that zenith cosine; otherwise unit
    isotropic diffuse light. rho and tau are (wavelengths,) or (layers, wavelengths); soil is
    (wavelengths,). Returns per layer the diffuse light (sky and scattered) absorbed per unit
    leaf area, the sunlit share of the layer's leaves and the direct light absorbed per unit
    area of a leaf facing the Sun; and per wavelength the light reflected upward at the top,
    reaching the soil, and absorbed by leaves and by soil.
    """
    n, dl = canopy.layers()
    rho = np.broadcast_to(np.asarray(rho, dtype=float), (n, np.shape(rho)[-1]))
    tau = np.broadcast_to(np.asarray(tau, dtype=float), rho.shape)
    soil = np.asarray(soil, dtype=float)
    nw = rho.shape[1]
    omega = rho + tau
    beta = backscatter(rho, tau)
    mu, w = gauss_streams(canopy.streams)
    share = 2 * mu * w                                        # flux share of each stream in isotropic light
    extinction = G * canopy.clumping * dl
    t = np.exp(-extinction / mu)
    if mu0 is not None:
        beam = np.exp(-extinction / mu0 * np.arange(n + 1))    # horizontal direct flux at the levels
        sky = 0.0
    else:
        beam = np.zeros(n + 1)
        sky = 1.0
    beam_hit = (beam[:-1] - beam[1:])[:, None]                # direct light intercepted by each layer
    e_up = np.zeros((n, nw))
    e_down = np.zeros((n, nw))

    def sweep(e_up, e_down):
        down = np.empty((n + 1, canopy.streams, nw))
        up = np.empty((n + 1, canopy.streams, nw))
        for j in range(canopy.streams):
            x = np.vstack([np.full((1, nw), share[j] * sky), share[j] * e_down])
            down[:, j] = lfilter([1.0], [1.0, -t[j]], x, axis=0)
        to_soil = down[n].sum(axis=0) + beam[n]
        for j in range(canopy.streams):
            x = np.vstack([share[j] * soil * to_soil, share[j] * e_up[::-1]])[None]
            up[::-1, j] = lfilter([1.0], [1.0, -t[j]], x[0], axis=0)
        hit_down = np.einsum('kjw,j->kw', down[:-1], 1 - t)
        hit_up = np.einsum('kjw,j->kw', up[1:], 1 - t)
        return down, up, to_soil, hit_down, hit_up

    for _ in range(max_iter):
        down, up, to_soil, hit_down, hit_up = sweep(e_up, e_down)
        from_above = hit_down + beam_hit
        new_up = omega * (beta * from_above + (1 - beta) * hit_up)
        new_down = omega * ((1 - beta) * from_above + beta * hit_up)
        change = max(np.max(np.abs(new_up - e_up)), np.max(np.abs(new_down - e_down)))
        e_up, e_down = new_up, new_down
        if change < tol:
            break
    else:
        raise RuntimeError('Canopy scattering did not converge')
    down, up, to_soil, hit_down, hit_up = sweep(e_up, e_down)
    absorbed_diffuse = (1 - omega) * (hit_down + hit_up)
    absorbed_direct = (1 - omega) * beam_hit
    sunlit = (mu0 * (beam[:-1] - beam[1:]) / (G * dl)) if mu0 is not None else np.zeros(n)
    return dict(absorbed_diffuse=absorbed_diffuse / dl,
                sunlit=sunlit,
                direct_on_facing_leaf=(1 - omega) / mu0 if mu0 is not None else np.zeros_like(omega),
                reflected=up[0].sum(axis=0),
                to_soil=to_soil,
                leaf_absorbed=(absorbed_diffuse + absorbed_direct).sum(axis=0),
                soil_absorbed=(1 - soil) * to_soil)
