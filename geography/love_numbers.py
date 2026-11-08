"""Load and tidal Love numbers of a layered, self-gravitating Moon.

Two states of the same interior answer how the crust moves under a surface load:

- Elastic: the static spheroidal equations of a layered, self-gravitating,
  elastic sphere (Takeuchi & Saito 1972) integrated from the fluid core to the
  surface, compressible or incompressible. Gives the tidal Love numbers h_n and
  k_n, checked against the Moon's measured k2 and h2, and the load Love numbers
  h'_n and k'_n of the immediate response.
- Relaxed: an elastic lithosphere, a thin spherical shell of effective thickness
  Te (Turcotte et al. 1981), over an interior that has flowed to hydrostatic
  equilibrium. The crust's column deflects with the shell and displaces mantle;
  the core and the interior's other density steps settle on their equipotentials.
  Self-gravitation enters through the mass sheets of the load, the deflected
  surface and Moho and the settled interfaces.

Conventions: the potential perturbation is positive over a mass excess,
y6 = dy5/dr - 4 pi G rho y1, and the equations are written in units of the
surface radius, the model's mean density and its surface gravity, where
4 pi G = 3. A surface load of density sigma_n (kg/m^2) has potential
psi_n = 4 pi G a sigma_n / (2n + 1); it moves the ground by h'_n psi_n / g and
the equipotential by (1 + k'_n) psi_n / g. Degree 1 is given in the frame of the
solid Moon's centre of mass (k'_1 = 0), where a rigid Moon has h'_1 = 0.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from shared.constants import GRAVITATIONAL_CONSTANT, MOON_RADIUS

HERE = Path(__file__).resolve().parent
SOURCES = HERE / "crustal_loading_sources.json"
FPG = 3.0                    # 4 pi G in units of the radius, mean density and surface gravity
START_DECADES = 25.0         # high degrees start where (r/a)^n = 10^-25


@dataclass(frozen=True)
class Layer:
    """A homogeneous shell: outer radius (m), density (kg/m^3), rigidity (Pa; 0 for a fluid),
    first Lame parameter (Pa; None for an incompressible solid)."""
    top_m: float
    density: float
    rigidity: float
    lame: float | None = None


def interior_model(path=SOURCES, compressible=True, poisson=0.25):
    """The published interior from the sources file, the crust carried to the Moon's radius.

    The source gives density and rigidity only (its code treats the Moon as incompressible);
    a compressible solid takes Poisson's ratio `poisson` (lambda = 2 mu nu / (1 - 2 nu))."""
    spec = json.loads(Path(path).read_text())["interior_model"]
    layers = []
    for row in spec["layers"]:
        mu = float(row["rigidity_pa"])
        lam = 2 * mu * poisson / (1 - 2 * poisson) if (compressible and mu > 0) else None
        layers.append(Layer(float(row["top_km"]) * 1e3, float(row["density_kg_m3"]), mu, lam))
    layers[-1] = Layer(MOON_RADIUS, layers[-1].density, layers[-1].rigidity, layers[-1].lame)
    return layers


def mass_and_inertia(layers):
    """Total mass (kg) and C / (M a^2) of a stack of homogeneous shells."""
    mass = inertia = 0.0
    bottom = 0.0
    for layer in layers:
        mass += 4 * math.pi / 3 * layer.density * (layer.top_m**3 - bottom**3)
        inertia += 8 * math.pi / 15 * layer.density * (layer.top_m**5 - bottom**5)
        bottom = layer.top_m
    return mass, inertia / (mass * layers[-1].top_m**2)


class Interior:
    """A layer stack in the units of the equations: radius, density, gravity and moduli."""

    def __init__(self, layers):
        self.a = layers[-1].top_m
        mass, _ = mass_and_inertia(layers)
        self.mean_density = mass / (4 * math.pi / 3 * self.a**3)
        self.g = GRAVITATIONAL_CONSTANT * mass / self.a**2
        stress = self.mean_density * self.g * self.a
        self.top = np.array([layer.top_m / self.a for layer in layers])
        self.bottom = np.r_[0.0, self.top[:-1]]
        self.rho = np.array([layer.density / self.mean_density for layer in layers])
        self.mu = np.array([layer.rigidity / stress for layer in layers])
        self.b = np.array([0.0 if layer.lame is None else stress / (layer.lame + 2 * layer.rigidity)
                           for layer in layers])
        self.below = np.r_[0.0, np.cumsum(self.rho * (self.top**3 - self.bottom**3))[:-1]]

    def gravity(self, k, r):
        return (self.below[k] + self.rho[k] * (r**3 - self.bottom[k]**3)) / r**2

    def matrix(self, k, r, n):
        """The 6 x 6 system dy/dr = A y in solid layer k (Takeuchi & Saito 1972, static)."""
        L = n * (n + 1.0)
        rho, mu, b, g = self.rho[k], self.mu[k], self.b[k], self.gravity(k, r)
        lb = 1.0 - 2.0 * mu * b                 # lambda / (lambda + 2 mu)
        c = 3.0 * mu - 4.0 * mu * mu * b        # mu (3 lambda + 2 mu) / (lambda + 2 mu)
        A = np.zeros((6, 6))
        A[0, 0], A[0, 1], A[0, 2] = -2.0 * lb / r, b, L * lb / r
        A[1, 0] = -4.0 * rho * g / r + 4.0 * c / r**2
        A[1, 1] = -4.0 * mu * b / r
        A[1, 2] = rho * g * L / r - 2.0 * c * L / r**2
        A[1, 3], A[1, 5] = L / r, -rho
        A[2, 0], A[2, 2], A[2, 3] = -1.0 / r, 1.0 / r, 1.0 / mu
        A[3, 0] = rho * g / r - 2.0 * c / r**2
        A[3, 1] = -lb / r
        A[3, 2] = 2.0 * mu * ((2.0 * L - 1.0) - 2.0 * mu * b * L) / r**2
        A[3, 3], A[3, 4] = -3.0 / r, -rho / r
        A[4, 0], A[4, 5] = FPG * rho, 1.0
        A[5, 2], A[5, 4], A[5, 5] = -FPG * rho * L / r, L / r**2, -2.0 / r
        return A


def _fluid_start(interior, n, top_index):
    """Three solutions at the top of the fluid region (layers 0..top_index) for degree n.

    The fluid's interfaces lie on equipotentials; its top meets the solid with any radial
    displacement, its stress set by hydrostatic equilibrium, and free tangential slip."""
    psi, dpsi = 1.0, n / interior.top[0]
    for k in range(1, top_index + 1):
        r0, r1 = interior.bottom[k], interior.top[k]
        u = psi / interior.gravity(k - 1, r0)
        dpsi -= FPG * (interior.rho[k - 1] - interior.rho[k]) * u
        a = ((n + 1) + dpsi * r0 / psi) / (2 * n + 1)
        x = r1 / r0
        grow, decay = x**n, x**(-n - 1)
        psi, dpsi = psi * (a * grow + (1 - a) * decay), psi * (n * a * grow - (n + 1) * (1 - a) * decay) / r1
        scale = abs(psi)
        psi, dpsi = psi / scale, dpsi / scale
    r = interior.top[top_index]
    rho, g = interior.rho[top_index], interior.gravity(top_index, r)
    start = np.zeros((6, 3))
    start[:, 0] = [0.0, -rho * psi, 0.0, 0.0, psi, dpsi]
    start[:, 1] = [1.0, rho * g, 0.0, 0.0, 0.0, -FPG * rho]
    start[:, 2] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    return start, r


def surface_solutions(interior, n, rtol=1e-10):
    """Three independent solutions at the surface, regular at depth, as a 6 x 3 matrix."""
    solid = np.flatnonzero(interior.mu > 0)
    first = int(solid[0])
    if np.any(interior.mu[first:] == 0):
        raise ValueError("Fluid layers must lie below every solid layer")
    fluid_top = interior.top[first - 1] if first > 0 else 0.0
    # With a solid centre, starting at r = 1e-4 leaves irregular parts of order 1e-4^(2n - 1).
    deep = max(10 ** (-START_DECADES / max(n, 1)), 1e-4 if first == 0 else 0.0)
    if deep > fluid_top:
        # Below `deep` the regular solutions are negligible at the surface; any three independent
        # vectors grow into them, and the irregular ones die out on the way up.
        Y, r = np.eye(6)[:, [0, 2, 4]].copy(), deep
    else:
        Y, r = _fluid_start(interior, n, first - 1)
    for k in range(first, len(interior.top)):
        if interior.top[k] <= r:
            continue
        r = max(r, interior.bottom[k])
        while r < interior.top[k] - 1e-15:
            r_next = min(interior.top[k], r * (1.0 + 4.0 / (n + 1)))
            sol = solve_ivp(lambda x, y: (interior.matrix(k, x, n) @ y.reshape(6, 3)).ravel(),
                            (r, r_next), Y.ravel(), method="DOP853", rtol=rtol, atol=1e-14 * np.abs(Y).max())
            Y, _ = np.linalg.qr(sol.y[:, -1].reshape(6, 3))
            r = r_next
    return Y


def elastic_love_numbers(interior, n):
    """Tidal (h, k) and load (h', k') Love numbers of degree n; degree 1 load only, in the
    frame of the solid Moon's centre of mass."""
    Y = surface_solutions(interior, n)
    if n == 1:
        # Centre of mass of Moon and load fixed: no external dipole, then shift to the solid Moon's frame.
        x = np.linalg.solve(np.vstack([Y[1], Y[3], Y[4]]), [-1.0, 0.0, 0.0])
        y = Y @ x
        return dict(n=1, h=None, k=None, h_load=1.0 + y[0], k_load=0.0, gauss_residual=float(y[5] + 2 * y[4] - 3.0))
    B = np.vstack([Y[1], Y[3], Y[5] + (n + 1) * Y[4]])
    tide = Y @ np.linalg.solve(B, [0.0, 0.0, 2 * n + 1.0])
    load = Y @ np.linalg.solve(B, [-(2 * n + 1) / 3.0, 0.0, 2 * n + 1.0])
    return dict(n=n, h=float(tide[0]), k=float(tide[4] - 1.0), h_load=float(load[0]), k_load=float(load[4] - 1.0))


def elastic_spectrum(interior, nmax, exact_up_to=40, samples=24):
    """h'_n and k'_n for n = 0..nmax: computed for n <= exact_up_to and at `samples` degrees
    spaced in log n above it, interpolated in log n between them (both vary smoothly there)."""
    exact = list(range(1, min(exact_up_to, nmax) + 1))
    sparse = sorted(set(int(n) for n in np.geomspace(exact_up_to + 1, max(nmax, exact_up_to + 2), samples).round()) - set(exact))
    rows = {n: elastic_love_numbers(interior, n) for n in exact + [n for n in sparse if n > exact_up_to]}
    known = np.array(sorted(rows))
    h = np.zeros(nmax + 1)
    k = np.zeros(nmax + 1)
    degrees = np.arange(1, nmax + 1)
    h[1:] = np.interp(np.log(degrees), np.log(known), [rows[n]["h_load"] for n in known])
    k[1:] = np.interp(np.log(degrees), np.log(known), [rows[n]["k_load"] * n for n in known]) / degrees
    return h, k, rows


# ------------------------------------------------------------- relaxed interior

def _sheet_potential(r, r_sheet, n):
    """Potential at radius r of a unit-density degree-n sheet at r_sheet, over 4 pi G r_sheet / (2n + 1)."""
    return (r_sheet / r) ** (n + 1) if r >= r_sheet else (r / r_sheet) ** n


def relaxed_love_numbers(layers, nmax, te_m, young_pa, poisson, moho_index=None):
    """Load Love numbers h'_n, k'_n (n = 0..nmax) of an elastic shell over a relaxed interior.

    The shell's flexural rigidity is E Te^3 / 12 (1 - nu^2) and its membrane stiffness E Te
    (Turcotte et al. 1981; Te the published effective elastic thickness, with that literature's
    E and nu). The crust's column deflects with the shell and displaces mantle at the base of
    the shell or at the Moho, whichever is deeper; the interfaces below settle on equipotentials.
    Degree 1 holds the centre of mass of Moon and load, then shifts to the solid Moon's frame."""
    a = layers[-1].top_m
    mass, _ = mass_and_inertia(layers)
    gravity = lambda r: GRAVITATIONAL_CONSTANT * _mass_within(layers, r) / r**2
    g = GRAVITATIONAL_CONSTANT * mass / a**2
    moho_index = len(layers) - 2 if moho_index is None else moho_index
    r_moho = layers[moho_index].top_m
    rho_crust, rho_mantle = layers[-1].density, layers[moho_index].density
    r_base = min(a - te_m, r_moho)
    g_base = gravity(r_base)
    # Fluid interfaces below the deflecting column: radius and density step (below minus above).
    fluid = [(layers[i].top_m, layers[i].density - layers[i + 1].density)
             for i in range(moho_index) if layers[i].top_m < r_base]
    D = young_pa * te_m**3 / (12 * (1 - poisson**2))
    h = np.zeros(nmax + 1)
    k = np.zeros(nmax + 1)
    for n in range(1, nmax + 1):
        L = n * (n + 1.0)
        shell = 0.0 if n == 1 else (D * (L**3 - 4 * L**2) + young_pa * te_m * a**2 * (L - 2)) / (a**4 * (L - 1 + poisson))
        # Sheets: the load at a (unit potential), the column's surface and Moho, then the fluid interfaces.
        radii = [a, a, r_moho] + [r for r, _ in fluid]
        unit = lambda r, rs: 4 * math.pi * GRAVITATIONAL_CONSTANT * rs * _sheet_potential(r, rs, n) / (2 * n + 1)
        sigma = (2 * n + 1) / (4 * math.pi * GRAVITATIONAL_CONSTANT * a)      # load with psi_L = 1
        m = 1 + len(fluid)
        A = np.zeros((m, m))
        rhs = np.zeros(m)
        # Unknowns: w (column deflection), u_j (fluid interfaces). Equipotential height N(r) = V(r) / g(r).
        def potential_row(r):
            row = np.zeros(m)
            row[0] = rho_crust * unit(r, a) + (rho_mantle - rho_crust) * unit(r, r_moho)
            for j, (rj, dr) in enumerate(fluid):
                row[1 + j] = dr * unit(r, rj)
            return row, sigma * unit(r, a)
        row, load = potential_row(r_base)
        A[0] = -rho_mantle * row
        A[0, 0] += shell + rho_mantle * g_base
        rhs[0] = -g * sigma + rho_mantle * load
        for j, (rj, _) in enumerate(fluid):
            row, load = potential_row(rj)
            A[1 + j] = -row / gravity(rj)
            A[1 + j, 1 + j] += 1.0
            rhs[1 + j] = load / gravity(rj)
        if n == 1:
            moment = np.r_[rho_crust * a**3 + (rho_mantle - rho_crust) * r_moho**3, [dr * rj**3 for rj, dr in fluid]]
            x, *_ = np.linalg.lstsq(np.vstack([A, moment]), np.r_[rhs, -sigma * a**3], rcond=None)
        else:
            x = np.linalg.solve(A, rhs)
        row, load = potential_row(a)
        h[n] = g * x[0]
        k[n] = row @ x
        if n == 1:
            # Into the solid Moon's frame: no external dipole from the Moon itself.
            h[1], k[1] = h[1] + 1.0, 0.0
    return h, k


def _mass_within(layers, r):
    mass, bottom = 0.0, 0.0
    for layer in layers:
        top = min(layer.top_m, r)
        if top > bottom:
            mass += 4 * math.pi / 3 * layer.density * (top**3 - bottom**3)
        bottom = layer.top_m
        if bottom >= r:
            break
    return mass


def relaxed_shell_ode(layers, n, te_m, young_pa, poisson):
    """Cross-check of the thin shell: the same shell as a thick elastic layer (uniform E, nu) over
    the relaxed interior, integrated as in the elastic case. Load (h', k') of degree n >= 2."""
    a = layers[-1].top_m
    lam = young_pa * poisson / ((1 + poisson) * (1 - 2 * poisson))
    mu = young_pa / (2 * (1 + poisson))
    stack = []
    for layer in layers:
        bottom = stack[-1].top_m if stack else 0.0
        if layer.top_m <= a - te_m:
            stack.append(Layer(layer.top_m, layer.density, 0.0, None))
        elif bottom < a - te_m:
            stack.append(Layer(a - te_m, layer.density, 0.0, None))
            stack.append(Layer(layer.top_m, layer.density, mu, lam))
        else:
            stack.append(Layer(layer.top_m, layer.density, mu, lam))
    row = elastic_love_numbers(Interior(stack), n)
    return row["h_load"], row["k_load"]


def relaxation_time_years(n, viscosity_pa_s, te_m, young_pa, poisson, mantle_density, gravity, radius=MOON_RADIUS):
    """e-folding time of degree n's approach to the relaxed state: a Newtonian half-space under a
    thin elastic plate, tau = 2 eta k / (rho g + D k^4), k = (n + 1/2) / a."""
    k = (n + 0.5) / radius
    D = young_pa * te_m**3 / (12 * (1 - poisson**2))
    return 2 * viscosity_pa_s * k / (mantle_density * gravity + D * k**4) / (365.25 * 86400.0)


def homogeneous_reference(n, rigidity_ratio):
    """Closed forms for a homogeneous incompressible elastic sphere, rigidity_ratio = mu / (rho g a):
    tidal h_n, k_n and load h'_n, k'_n."""
    m = (2 * n**2 + 4 * n + 3) / n * rigidity_ratio
    return dict(h=(2 * n + 1) / (2 * (n - 1)) / (1 + m), k=1.5 / (n - 1) / (1 + m),
                h_load=-(2 * n + 1) / 3 / (1 + m), k_load=-1 / (1 + m))
