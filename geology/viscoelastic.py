"""Load and tidal Love numbers of a layered, self-gravitating Maxwell Moon through time.

The static spheroidal equations of a self-gravitating sphere (Takeuchi & Saito 1972) hold in the Laplace domain
for a viscoelastic body once each shell's rigidity is replaced by its transform: for a Maxwell solid
mu(s) = mu s / (s + mu / eta), the correspondence principle. At each real s the equations are integrated from the
fluid core to the surface exactly as for an elastic Moon, and the boundary conditions give the load Love numbers
h'_n(s), k'_n(s) and the tidal h_n(s), k_n(s). Shells whose transformed rigidity has fallen to a negligible
fraction of rho g a at that s join the fluid region above the core, where every density step lies on an
equipotential.

Each shell is crossed in steps of the fourth-order Magnus integrator (the matrix exponential of the two-point
Gauss average and its commutator), batched over s, with QR re-orthonormalisation between steps; steps shrink with
degree and with the stiffness that a weak shell's small rigidity brings. The time response to a load applied at
t = 0 is the inverse transform of k(s) / s. It is found as a relaxation spectrum on fixed rates spaced evenly in
log s, k(s) = k_e + sum_j a_j r_j / (s + r_j), fitted by non-negative least squares to the computed k(s) between the
elastic and the fluid limits; in time this is k(t) = k_e + sum_j a_j (1 - exp(-r_j t)), a discretised sum of
normal modes. Gaver-Stehfest inversion of the same k(s) checks it independently.

The time response uses incompressible shells, as the interior's source does (Briaud et al. 2023, whose ALMA code
is incompressible): a compressible Maxwell shell of uniform density is unstably stratified in its fluid limit.
Shells may carry a bulk modulus for the elastic check against the crustal-loading study's compressible numbers.

Conventions follow the crustal-loading study (geography/love_numbers.py on domain/hydrosphere): the potential is
positive over a mass excess, y6 = dy5/dr - 4 pi G rho y1, and the equations are written in units of the surface
radius, the model's mean density and its surface gravity, where 4 pi G = 3. A surface load sigma_n has potential
psi_n = 4 pi G a sigma_n / (2n + 1); it moves the ground by h'_n psi_n / g and the equipotential by
(1 + k'_n) psi_n / g. Times are in years and s in 1/yr.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import expm
from scipy.optimize import nnls

from shared.constants import GRAVITATIONAL_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS

FPG = 3.0                       # 4 pi G in the equations' units
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
FLUID_FRACTION = 1e-9           # a shell is fluid where mu(s) (2n + 1) / (rho g a) falls below this
START_DECADES = 25.0            # high degrees start where (r/a)^n = 1e-25, below which nothing reaches the surface
MAX_EXPONENT = 20.0             # largest growth allowed within one step (e^20)
RATES_PER_YEAR = np.logspace(-16.0, 7.0, 185)         # relaxation spectrum, 8 rates a decade
FIT_S_PER_YEAR = np.logspace(-13.0, 8.0, 127)         # where k(s) is fitted, 6 a decade
GAUSS = math.sqrt(3.0) / 6.0


@dataclass(frozen=True)
class Shell:
    """A homogeneous shell: outer radius (m), density (kg/m^3), rigidity (Pa; 0 for an inviscid fluid), viscosity
    (Pa s; infinite for an elastic solid) and bulk modulus (Pa; None for an incompressible solid)."""
    top_m: float
    density: float
    rigidity: float
    viscosity: float = math.inf
    bulk: float | None = None


class Moon:
    """A stack of shells, centre outward, in the units of the equations."""

    def __init__(self, shells):
        self.shells = list(shells)
        self.a = self.shells[-1].top_m
        top = np.array([s.top_m for s in self.shells])
        bottom = np.r_[0.0, top[:-1]]
        rho = np.array([s.density for s in self.shells])
        self.mass = float(4 * math.pi / 3 * np.sum(rho * (top**3 - bottom**3)))
        self.mean_density = self.mass / (4 * math.pi / 3 * self.a**3)
        self.g = GRAVITATIONAL_CONSTANT * self.mass / self.a**2
        self.stress = self.mean_density * self.g * self.a
        self.top, self.bottom = top / self.a, bottom / self.a
        self.rho = rho / self.mean_density
        self.mu = np.array([s.rigidity for s in self.shells]) / self.stress
        self.relax = np.array([0.0 if (s.rigidity == 0 or math.isinf(s.viscosity)) else
                               s.rigidity / s.viscosity * YEAR_S for s in self.shells])     # mu / eta, 1/yr
        self.kappa = np.array([np.inf if s.bulk is None else s.bulk / self.stress for s in self.shells])
        self.below = np.r_[0.0, np.cumsum(self.rho * (self.top**3 - self.bottom**3))[:-1]]
        if self.mu[0] > 0:
            self.core_top = 0
        else:
            solid = np.flatnonzero(self.mu > 0)
            if np.any(self.mu[solid[0]:] == 0):
                raise ValueError("Inviscid fluid shells must lie below every solid shell")
            self.core_top = int(solid[0])           # first shell above the inviscid core

    def gravity(self, k, r):
        return (self.below[k] + self.rho[k] * (r**3 - self.bottom[k]**3)) / r**2

    def rigidity(self, k, s):
        """mu(s) of shell k at Laplace variable s (1/yr; inf for the elastic limit)."""
        s = np.asarray(s)
        if self.relax[k] == 0.0:
            return np.full(s.shape, self.mu[k], dtype=np.result_type(s, float))
        with np.errstate(invalid="ignore"):
            out = self.mu[k] * s / (s + self.relax[k])
        return np.where(np.isinf(s), self.mu[k], out)

    def matrices(self, k, r, n, mu, b):
        """The 6 x 6 systems dy/dr = A y of shell k at radius r, one per s (mu, b arrays)."""
        L = n * (n + 1.0)
        rho, g = self.rho[k], self.gravity(k, r)
        lb = 1.0 - 2.0 * mu * b                 # lambda / (lambda + 2 mu)
        c = 3.0 * mu - 4.0 * mu * mu * b        # mu (3 lambda + 2 mu) / (lambda + 2 mu)
        A = np.zeros(mu.shape + (6, 6), dtype=mu.dtype)
        A[:, 0, 0], A[:, 0, 1], A[:, 0, 2] = -2.0 * lb / r, b, L * lb / r
        A[:, 1, 0] = -4.0 * rho * g / r + 4.0 * c / r**2
        A[:, 1, 1] = -4.0 * mu * b / r
        A[:, 1, 2] = rho * g * L / r - 2.0 * c * L / r**2
        A[:, 1, 3], A[:, 1, 5] = L / r, -rho
        A[:, 2, 0], A[:, 2, 2], A[:, 2, 3] = -1.0 / r, 1.0 / r, 1.0 / mu
        A[:, 3, 0] = rho * g / r - 2.0 * c / r**2
        A[:, 3, 1] = -lb / r
        A[:, 3, 2] = 2.0 * mu * ((2.0 * L - 1.0) - 2.0 * mu * b * L) / r**2
        A[:, 3, 3], A[:, 3, 4] = -3.0 / r, -rho / r
        A[:, 4, 0], A[:, 4, 5] = FPG * rho, 1.0
        A[:, 5, 2], A[:, 5, 4], A[:, 5, 5] = -FPG * rho * L / r, L / r**2, -2.0 / r
        return A

    def compliance(self, k, mu):
        """1 / (lambda + 2 mu) of shell k for transformed rigidity mu: 0 if incompressible."""
        if np.isinf(self.kappa[k]):
            return np.zeros_like(mu)
        return 1.0 / (self.kappa[k] + 4.0 * mu / 3.0)

    def fluid_start(self, n, j):
        """Potential and its gradient at the top of the fluid region of shells 0..j, interfaces on equipotentials,
        and the three solutions there: the potential, a free radial displacement and free slip."""
        psi, dpsi = 1.0, n / self.top[0]
        for k in range(1, j + 1):
            r0 = self.bottom[k]
            u = psi / self.gravity(k - 1, r0)
            dpsi -= FPG * (self.rho[k - 1] - self.rho[k]) * u
            grow = ((n + 1) * psi + dpsi * r0) / (2 * n + 1)      # coefficient of (r/r0)^n
            x = self.top[k] / r0
            psi, dpsi = grow * x**n + (psi - grow) * x**(-n - 1), (n * grow * x**n - (n + 1) * (psi - grow) * x**(-n - 1)) / self.top[k]
            scale = abs(psi) or 1.0
            psi, dpsi = psi / scale, dpsi / scale
        r = self.top[j]
        rho, g = self.rho[j], self.gravity(j, r)
        start = np.zeros((6, 3))
        start[:, 0] = [0.0, -rho * psi, 0.0, 0.0, psi, dpsi]
        start[:, 1] = [1.0, rho * g, 0.0, 0.0, 0.0, -FPG * rho]
        start[:, 2] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        return start

    def surface(self, n, s, step=0.01):
        """Three independent solutions at the surface for each s, as an array (len(s), 6, 3)."""
        s = np.atleast_1d(np.asarray(s))
        dtype = np.result_type(s, float)
        ns, nk = len(s), len(self.shells)
        scale = (2 * n + 1)
        weak = np.array([np.abs(self.rigidity(k, s)) * scale < FLUID_FRACTION * self.rho[k] * self.gravity(k, self.top[k]) * self.top[k]
                         if k >= self.core_top else np.ones(ns, bool) for k in range(nk)])
        # The fluid region for each s: the inviscid core and the weak shells contiguous above it.
        fluid_top = np.full(ns, self.core_top - 1)
        for k in range(self.core_top, nk):
            extend = weak[k] & (fluid_top == k - 1)
            fluid_top[extend] = k
        deep = 10 ** (-START_DECADES / n)
        Y = np.zeros((ns, 6, 3), dtype=dtype)
        started = np.zeros(ns, bool)
        for k in range(self.core_top, nk):
            if deep >= self.top[k]:
                continue                        # below where the regular solutions reach the surface
            # s whose solid part begins in this shell: at `deep` inside it (any three independent vectors grow into
            # the regular solutions), at the centre of a solid sphere, or on the fluid region below
            begin = (~started) & (fluid_top < k)
            r_lo = self.bottom[k]
            if np.any(begin):
                if deep > r_lo or k == 0:
                    r_lo = max(deep, r_lo, 1e-4)
                    Y[begin] = np.eye(6)[:, [0, 2, 4]]
                else:
                    Y[begin] = self.fluid_start(n, k - 1)
                started |= begin
            if not np.any(started):
                continue
            idx = np.flatnonzero(started)
            mu = self.rigidity(k, s[idx])
            b = self.compliance(k, mu)
            # s that need about the same number of steps cross the shell together
            thickness = self.top[k] - r_lo
            stiff = (self.rho[k] * max(self.gravity(k, r_lo), 1e-12) * n * (n + 1) / (np.abs(mu) * r_lo**2)) ** (1 / 3)
            count = np.maximum(thickness * (stiff + (n + 1) / r_lo) / MAX_EXPONENT, thickness / min(step, 0.25 * r_lo / (n + 1)))
            group = np.ceil(np.log2(np.maximum(count, 1.0))).astype(int)
            for key in np.unique(group):
                sel = group == key
                Y[idx[sel]] = self.cross(k, n, r_lo, Y[idx[sel]], mu[sel], b[sel], step)
        # A body fluid to its surface: the free-slip column carries no stress, so stand a unit tangential stress
        # in its place; the boundary condition y4 = 0 then gives it no weight.
        whole = fluid_top == nk - 1
        if np.any(whole):
            start = self.fluid_start(n, nk - 1)
            start[:, 2] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            Y[whole] = start
        return Y

    def cross(self, k, n, r_lo, y, mu, b, step):
        """Carry the solutions y (len(mu), 6, 3) across shell k from r_lo to its top."""
        weakest = float(np.abs(mu).min())
        r = r_lo
        while r < self.top[k] * (1 - 1e-12):
            # steps grow with radius, shrink with degree and with a weak shell's stiffness
            stiff = (self.rho[k] * max(self.gravity(k, r), 1e-12) * n * (n + 1) / (weakest * r * r)) ** (1 / 3) + (n + 1) / r
            h = min(step, 0.25 * r / (n + 1), MAX_EXPONENT / stiff, self.top[k] - r)
            if self.top[k] - (r + h) < 1e-3 * h:
                h = self.top[k] - r
            A1 = self.matrices(k, r + (0.5 - GAUSS) * h, n, mu, b)
            A2 = self.matrices(k, r + (0.5 + GAUSS) * h, n, mu, b)
            omega = 0.5 * h * (A1 + A2) + (math.sqrt(3.0) / 12.0) * h * h * (A2 @ A1 - A1 @ A2)
            y = expm(omega) @ y
            y, _ = np.linalg.qr(y)
            r += h
        return y

    def love(self, n, s, step=0.01):
        """Load (h', k') and tidal (h, k) Love numbers of degree n >= 2 at each s, as arrays."""
        Y = self.surface(n, s, step)
        B = np.stack([Y[:, 1], Y[:, 3], Y[:, 5] + (n + 1) * Y[:, 4]], axis=1)
        out = {}
        for name, rhs in (("load", [-(2 * n + 1) / 3.0, 0.0, 2 * n + 1.0]), ("tide", [0.0, 0.0, 2 * n + 1.0])):
            c = np.linalg.solve(B, np.broadcast_to(np.array(rhs, dtype=B.dtype), B.shape[:2])[..., None])
            out[name] = (Y @ c)[..., 0]
        return dict(h_load=out["load"][:, 0], k_load=out["load"][:, 4] - 1.0, h=out["tide"][:, 0], k=out["tide"][:, 4] - 1.0)


# ------------------------------------------------------------- the time domain

def fit_spectrum(s, values, elastic, rates=RATES_PER_YEAR):
    """Non-negative amplitudes a_j (signed as the relaxation) with values - elastic = sum_j a_j r_j / (s + r_j).
    Returns the signed amplitudes and the largest misfit relative to the relaxation's size."""
    target = np.asarray(values, float) - elastic
    sign = 1.0 if target[np.argmin(s)] >= 0 else -1.0
    G = rates[None, :] / (s[:, None] + rates[None, :])
    amp, _ = nnls(G, sign * target, maxiter=50 * len(rates))
    misfit = float(np.max(np.abs(G @ amp - sign * target)) / max(np.max(np.abs(target)), 1e-30))
    return sign * amp, misfit


def step_response(t_years, elastic, amplitudes, rates=RATES_PER_YEAR):
    """k(t) after a load applied at t = 0: elastic + sum_j a_j (1 - exp(-r_j t))."""
    t = np.atleast_1d(np.asarray(t_years, float))
    return elastic + (1.0 - np.exp(-np.outer(t, rates))) @ amplitudes


def stehfest_weights(N=14):
    V = np.zeros(N)
    half = N // 2
    for k in range(1, N + 1):
        total = 0.0
        for j in range((k + 1) // 2, min(k, half) + 1):
            total += (j**half * math.factorial(2 * j)) / (math.factorial(half - j) * math.factorial(j) * math.factorial(j - 1)
                                                         * math.factorial(k - j) * math.factorial(2 * j - k))
        V[k - 1] = (-1) ** (k + half) * total
    return V


def stehfest_points(t_years, N=14):
    """The s at which Gaver-Stehfest needs k(s) for a step response at time t."""
    return np.arange(1, N + 1) * math.log(2.0) / t_years


def stehfest_step(values_at_points, N=14):
    """Step response from k(s) at stehfest_points: f(t) = sum_k V_k / k * k(s_k), since F(s) = k(s) / s."""
    return float(np.sum(stehfest_weights(N) / np.arange(1, N + 1) * np.asarray(values_at_points)))


def response(moon, n, times_years, kinds=("h_load", "k_load", "h", "k")):
    """Elastic and fluid-side limits, the fitted spectrum and the step responses of degree n."""
    s = FIT_S_PER_YEAR
    values = moon.love(n, np.r_[np.inf, s])
    out = {}
    for kind in kinds:
        elastic = float(values[kind][0].real)
        amp, misfit = fit_spectrum(s, values[kind][1:].real, elastic)
        out[kind] = dict(elastic=elastic, slowest=float(values[kind][1 + int(np.argmin(s))].real),
                         amplitudes=amp, misfit=misfit, series=step_response(times_years, elastic, amp))
    return out


def homogeneous_maxwell(n, rigidity_ratio, maxwell_years):
    """Closed forms of a homogeneous incompressible Maxwell sphere, rigidity_ratio = mu / (rho g a): load and tidal
    elastic and fluid Love numbers and the single relaxation time (years)."""
    m = (2 * n**2 + 4 * n + 3) / n * rigidity_ratio
    return dict(k_load=(-1 / (1 + m), -1.0), h_load=(-(2 * n + 1) / 3 / (1 + m), -(2 * n + 1) / 3),
                k=(1.5 / (n - 1) / (1 + m), 1.5 / (n - 1)), h=((2 * n + 1) / (2 * (n - 1)) / (1 + m), (2 * n + 1) / (2 * (n - 1))),
                tau_years=maxwell_years * (1 + m))
