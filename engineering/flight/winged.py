"""Winged and rotor flight on any world: the speed a wing needs to carry a mass, the gusts it meets, what its wing
weighs as it grows, and the power to fly or hover.

**Lift and power.** A wing of area S carrying mass m flies where m g = 0.5 rho V^2 S C_L, so at a given lift
coefficient its speed goes as sqrt(g (m/S) / rho). The power to fly at a lift-to-drag ratio, m g V / (L/D), goes as
g at a given speed and as g^1.5 at a given lift coefficient; a hovering rotor's ideal induced power,
(m g)^1.5 / sqrt(2 rho A), as g^1.5 for a given disc; a kilometre of flight costs m g / (L/D), in proportion to g.

**Gusts.** The limit load factor is the larger of the manoeuvre factor and 1 + Δn from a gust U at speed V, with
Pratt's alleviation factor as 14 CFR 23.341 gives it: Δn = rho U V a K_g / (2 g m/S) with K_g = 0.88 mu / (5.3 + mu)
and mu = 2 (m/S) / (rho c a). Gravity cancels in mu. The gust's extra lift per square metre,
Δn g m/S, is the same on any world at the same speed and density, so a flyer carrying little weight per square metre
takes gusts as larger load factors.

**Wing weight.** The spar caps carry the bending of the ultimate lift, n_ult m g, spread elliptically across the
span, on a box whose depth is a share of the local thickness. For a straight-tapered wing their ideal mass, as a
share of the flyer's mass, is n_ult g b AR (rho_m / sigma) times a planform factor over (t/c): it grows with span and
with gravity, so under a sixth of the gravity a wing of the same proportions and load factor carries its weight with
the same share at six times the span. A non-optimum factor (webs, ribs, skins, joints, control surfaces, fuel
sealing) scales the ideal caps, and an areal term per square metre of wing covers what no load sizes. The sky-ship
study fits both to the wing groups of ten transports. This is first-order sizing, for comparing sizes and worlds.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from engineering.flight.buoyant import Material

# The spar caps carry the ultimate load at this stress; the fitted non-optimum factor absorbs the difference from the
# alloys' real allowables.
ALUMINIUM = Material('aluminium alloy wing box (7000-series)', 2800.0, 300e6)


def speed_for_lift(mass_per_area_kg_m2, gravity_m_s2, density_kg_m3, lift_coefficient):
    """Airspeed (m/s) at which a wing carries its mass per square metre at a lift coefficient."""
    return np.sqrt(2.0 * gravity_m_s2 * np.asarray(mass_per_area_kg_m2) / (density_kg_m3 * lift_coefficient))


def gust_load_factor(mass_per_area_kg_m2, chord_m, speed_m_s, gust_m_s, density_kg_m3, gravity_m_s2,
                     lift_slope_per_rad=5.0):
    """Pratt's increment in load factor from a sharp-edged gust, with the gust alleviation factor."""
    ms = np.asarray(mass_per_area_kg_m2, dtype=float)
    mu = 2.0 * ms / (density_kg_m3 * chord_m * lift_slope_per_rad)
    kg = 0.88 * mu / (5.3 + mu)
    return density_kg_m3 * gust_m_s * speed_m_s * lift_slope_per_rad * kg / (2.0 * gravity_m_s2 * ms)


def planform_factor(taper: float, n: int = 2000) -> float:
    """∫ mu(eta) / (1 - (1 - taper) eta) d eta over the semi-span, with mu the elliptic loading's bending moment per
    (total lift x span) at eta = 2y/b; with a root chord 2S / (b (1 + taper)), cap mass = n_ult m g b AR (rho/sigma)
    (1 + taper) x this / (depth share x t/c)."""
    eta = (np.arange(n) + 0.5) / n
    # Elliptic lift per unit span, per (total lift / semi-span): (4/pi) sqrt(1 - eta^2); moment at eta, per
    # (lift x span): integral of the load outboard of eta times its arm, over the semi-span.
    s = np.sqrt(np.clip(1.0 - eta ** 2, 0.0, None))
    load = 4.0 / np.pi * s / n                    # share of the half-wing's lift in each strip (sums to 1)
    mu = np.array([0.5 * (0.5 * (load[i:] * (eta[i:] - eta[i])).sum()) for i in range(n)])  # (L/2)(b/2) arm
    return float((mu / (1.0 - (1.0 - taper) * eta)).sum() / n)


@dataclass(frozen=True)
class Wing:
    span_m: float
    aspect_ratio: float = 9.0
    taper: float = 0.3
    thickness: float = 0.12        # thickness / chord at the root
    depth_share: float = 0.8       # effective depth of the spar caps / thickness

    @property
    def area_m2(self):
        return self.span_m ** 2 / self.aspect_ratio

    @property
    def mean_chord_m(self):
        return self.span_m / self.aspect_ratio


def wing_mass_kg(wing: Wing, mass_kg, gravity_m_s2, ultimate_load_factor, material: Material = ALUMINIUM,
                 non_optimum: float = 2.0, areal_kg_m2: float = 10.0):
    """Wing structural mass: the ideal spar caps times the non-optimum factor, plus an areal term."""
    caps = (ultimate_load_factor * np.asarray(mass_kg, dtype=float) * gravity_m_s2 * wing.span_m * wing.aspect_ratio
            * material.density_kg_m3 / material.allowable_pa * (1.0 + wing.taper) * planform_factor(wing.taper)
            / (wing.depth_share * wing.thickness))
    return non_optimum * caps + areal_kg_m2 * wing.area_m2


def power_w(mass_kg, gravity_m_s2, speed_m_s, lift_to_drag, efficiency=0.8):
    """Power to fly level at a lift-to-drag ratio."""
    return np.asarray(mass_kg) * gravity_m_s2 * speed_m_s / (lift_to_drag * efficiency)


def hover_power_w(mass_kg, gravity_m_s2, density_kg_m3, disc_area_m2, figure_of_merit=0.7):
    """Power to hover: momentum theory's induced power over the figure of merit."""
    thrust = np.asarray(mass_kg) * gravity_m_s2
    return thrust ** 1.5 / np.sqrt(2.0 * density_kg_m3 * disc_area_m2) / figure_of_merit
