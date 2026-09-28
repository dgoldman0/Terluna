"""Winged and rotor flight on any world: the speed a wing needs to carry a mass, the gusts it meets, what its wing
weighs as it grows, the power to fly or hover, a flyer's drag polar (power, sink and glide against speed), the
descent under a canopy, and the energy and mass of an electric flyer that takes off vertically.

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

**The drag polar.** A flyer of span b and parasite drag area f (the wing's profile drag and every other part's) has
D = q f + (m g)^2 / (q pi e b^2). Its best glide ratio, 0.5 sqrt(pi e b^2 / f), depends on neither gravity, air nor
mass; the speeds of least drag and least power go as sqrt(g), the least power as g^1.5 and the least sink as
sqrt(g). So under a sixth of the gravity the same glider glides as far, 2.46 times more slowly, and sinks 2.46 times
more slowly; a flyer of the same shape with a sixth of the wing area and drag area flies at Earth's speed and sink.

**Canopies.** Under drag alone a mass descends at sqrt(2 m g / (rho C_D A)).

**Vertical take-off.** A lift-and-cruise flyer hovers at the momentum-theory power over the figure of merit,
climbs against gravity and cruises on its wing; an electric one carries its trip's energy in batteries and drive
for its hover power, and its airframe as a share of its mass.
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


# ---------------------------------------------------------------------------------------------------------------
# A flyer's drag polar: the power to fly, the sink in a glide, turning

@dataclass(frozen=True)
class Polar:
    """A parabolic drag polar, D = q f + (m g)^2 / (q pi e b^2): q the dynamic pressure, f the parasite drag area
    (the wing's profile drag and every other part's), e the span efficiency. The lift coefficient stays under its
    maximum."""
    span_m: float
    area_m2: float
    drag_area_m2: float
    span_efficiency: float = 0.9
    max_lift_coefficient: float = 1.4

    @property
    def aspect_ratio(self):
        return self.span_m ** 2 / self.area_m2

    @property
    def best_glide_ratio(self):
        """(L/D)max = 0.5 sqrt(pi e b^2 / f), whatever the gravity, the air and the mass."""
        return 0.5 * np.sqrt(np.pi * self.span_efficiency * self.span_m ** 2 / self.drag_area_m2)

    def drag_n(self, weight_n, density_kg_m3, speed_m_s):
        q = 0.5 * density_kg_m3 * np.asarray(speed_m_s, dtype=float) ** 2
        induced = np.asarray(weight_n, dtype=float) ** 2 / (q * np.pi * self.span_efficiency * self.span_m ** 2)
        return q * self.drag_area_m2 + induced

    def stall_speed_m_s(self, weight_n, density_kg_m3):
        return np.sqrt(2.0 * np.asarray(weight_n, dtype=float) / (density_kg_m3 * self.area_m2 *
                                                                  self.max_lift_coefficient))


def flight_point(polar: Polar, mass_kg, gravity_m_s2, density_kg_m3, speed_m_s, efficiency=1.0) -> dict:
    """Level flight (or a steady glide) at a speed: drag, lift-to-drag ratio, shaft power over the efficiency, and
    the sink rate the same flyer has gliding at that speed."""
    w = np.asarray(mass_kg, dtype=float) * gravity_m_s2
    v = np.asarray(speed_m_s, dtype=float)
    d = polar.drag_n(w, density_kg_m3, v)
    return dict(speed_m_s=v, drag_n=d, lift_to_drag=w / d, power_w=d * v / efficiency, sink_m_s=d * v / w)


def _polar_terms(polar: Polar, weight_n, density_kg_m3):
    """D = A V^2 + B / V^2."""
    a = 0.5 * density_kg_m3 * polar.drag_area_m2
    b = np.asarray(weight_n, dtype=float) ** 2 / (0.5 * density_kg_m3 * np.pi * polar.span_efficiency *
                                                  polar.span_m ** 2)
    return a, b


def min_power(polar: Polar, mass_kg, gravity_m_s2, density_kg_m3, efficiency=1.0, stall_margin=1.2) -> dict:
    """Flight at the speed of least power, V^4 = B / 3A, or at `stall_margin` times the stall speed where that is
    faster: the least power to stay up, and the least sink in a glide."""
    w = np.asarray(mass_kg, dtype=float) * gravity_m_s2
    a, b = _polar_terms(polar, w, density_kg_m3)
    v = np.maximum((b / (3.0 * a)) ** 0.25, stall_margin * polar.stall_speed_m_s(w, density_kg_m3))
    return flight_point(polar, mass_kg, gravity_m_s2, density_kg_m3, v, efficiency)


def best_glide(polar: Polar, mass_kg, gravity_m_s2, density_kg_m3, stall_margin=1.2) -> dict:
    """Flight at the speed of least drag, V^4 = B / A, or at `stall_margin` times the stall speed where faster:
    the flattest glide, and the least energy per kilometre in powered flight."""
    w = np.asarray(mass_kg, dtype=float) * gravity_m_s2
    a, b = _polar_terms(polar, w, density_kg_m3)
    v = np.maximum((b / a) ** 0.25, stall_margin * polar.stall_speed_m_s(w, density_kg_m3))
    return flight_point(polar, mass_kg, gravity_m_s2, density_kg_m3, v)


def turn_radius_m(speed_m_s, gravity_m_s2, bank_deg):
    """The radius of a level, coordinated turn: V^2 / (g tan(bank))."""
    return np.asarray(speed_m_s, dtype=float) ** 2 / (gravity_m_s2 * np.tan(np.radians(bank_deg)))


# ---------------------------------------------------------------------------------------------------------------
# Canopies

def descent_speed_m_s(mass_kg, gravity_m_s2, density_kg_m3, drag_area_m2):
    """Steady descent under drag alone, from m g = 0.5 rho v^2 (C_D A)."""
    return np.sqrt(2.0 * np.asarray(mass_kg, dtype=float) * gravity_m_s2 / (density_kg_m3 * drag_area_m2))


def canopy_diameter_m(mass_kg, gravity_m_s2, density_kg_m3, speed_m_s, drag_coefficient):
    """The round canopy, with its drag coefficient on its circular area, that lands a mass at a speed."""
    area = 2.0 * np.asarray(mass_kg, dtype=float) * gravity_m_s2 / (density_kg_m3 * drag_coefficient * speed_m_s ** 2)
    return np.sqrt(4.0 * area / np.pi)


# ---------------------------------------------------------------------------------------------------------------
# Vertical take-off: a trip's energy and an electric flyer's mass

def vtol_trip(mass_kg, gravity_m_s2, density_kg_m3, disc_area_m2, lift_to_drag, cruise_m_s, distance_m,
              hover_s=90.0, climb_m=0.0, figure_of_merit=0.7, efficiency=0.8) -> dict:
    """A trip of a lift-and-cruise flyer: hovering (take-off, landing, waiting) at the momentum-theory power over
    the figure of merit, climbing against gravity, and cruising on its wing at a lift-to-drag ratio. Energy in J."""
    m = np.asarray(mass_kg, dtype=float)
    hover = hover_power_w(m, gravity_m_s2, density_kg_m3, disc_area_m2, figure_of_merit)
    cruise = power_w(m, gravity_m_s2, cruise_m_s, lift_to_drag, efficiency)
    parts = dict(hover_j=hover * hover_s, climb_j=m * gravity_m_s2 * climb_m / efficiency,
                 cruise_j=cruise * distance_m / cruise_m_s)
    return dict(hover_w=hover, cruise_w=cruise, energy_j=sum(parts.values()), **parts)


@dataclass(frozen=True)
class Electric:
    """An electric lift-and-cruise flyer's technology: the airframe (structure, cabin, systems) as a share of its
    mass, the rotors' disc loading in kg of take-off mass per m2, and the drive and batteries."""
    airframe_share: float = 0.40
    disc_loading_kg_m2: float = 50.0
    lift_to_drag: float = 12.0
    figure_of_merit: float = 0.7
    efficiency: float = 0.8            # propeller, motor and controller in cruise
    drive_w_kg: float = 2000.0         # motors, controllers and rotors per watt of hover power
    battery_wh_kg: float = 200.0       # pack
    usable: float = 0.8                # share of the pack's energy a trip may use
    reserve: float = 1.3               # energy carried over the trip's


def electric_vtol(payload_kg, gravity_m_s2, density_kg_m3, cruise_m_s, range_m, hover_s=90.0, climb_m=0.0,
                  tech: Electric = Electric(), iterations: int = 100) -> dict:
    """Size an electric lift-and-cruise flyer for a payload and a trip: take-off mass = payload + airframe (a share
    of the take-off mass) + drive (for the hover power) + batteries (the trip's energy with the reserve, over the
    usable share of the pack). Returns None where no mass closes."""
    m = payload_kg / (1.0 - tech.airframe_share)
    for _ in range(iterations):
        trip = vtol_trip(m, gravity_m_s2, density_kg_m3, m / tech.disc_loading_kg_m2, tech.lift_to_drag, cruise_m_s,
                         range_m, hover_s, climb_m, tech.figure_of_merit, tech.efficiency)
        battery = trip['energy_j'] * tech.reserve / (tech.usable * tech.battery_wh_kg * 3600.0)
        drive = trip['hover_w'] / tech.drive_w_kg
        new = (payload_kg + battery + drive) / (1.0 - tech.airframe_share)
        if not np.isfinite(new) or new > 1e3 * payload_kg:
            return None
        if abs(new - m) < 1e-9 * m:
            m = new
            break
        m = new
    trip = vtol_trip(m, gravity_m_s2, density_kg_m3, m / tech.disc_loading_kg_m2, tech.lift_to_drag, cruise_m_s,
                     range_m, hover_s, climb_m, tech.figure_of_merit, tech.efficiency)
    return dict(mass_kg=m, payload_kg=payload_kg, airframe_kg=tech.airframe_share * m,
                battery_kg=trip['energy_j'] * tech.reserve / (tech.usable * tech.battery_wh_kg * 3600.0),
                drive_kg=trip['hover_w'] / tech.drive_w_kg, disc_area_m2=m / tech.disc_loading_kg_m2, **trip)
