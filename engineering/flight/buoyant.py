"""Buoyant ships on any world: gross lift, the hull, and the loads that decide how large a ship can be.

**Lift.** A buoyant ship's gross lift is the mass of air its gas displaces, less the gas's own mass. Per cubic metre
of gas cell that is the air's density less the density of the gas mixture in the cell (the lifting gas at its
purity, the rest air, at the gas's temperature). Gravity cancels: a cubic metre of hydrogen lifts the same mass on
the Moon as on Earth in air of the same density.

**Loads.** Gravity returns in the loads, which fall in three families that scale differently with size:

- *Weight-driven.* Bending from an uneven spread of lift and weight (the classic design case is a deflated gas cell,
  which leaves a gap in the lift), the gas's hydrostatic pressure on the hull's frames, and concentrated weights. In
  hulls of one shape the material these need, as a share of the lift, grows as g L rho_m / sigma.
- *Aerodynamic.* A gust U across a hull flying at V bends it with at most 0.10 (U/V) q Vol^(2/3) L between 0.5 and
  0.65 of its length: Woodward's envelope (1976), built from the Guggenheim Airship Institute's tank tests on
  free-flying models of the Akron with their fins and from Calligeros and McDavitt's gust theory, which covers
  the other aerodynamic cases measured or calculated on the rigid airships (circling, rudder reversal, flight
  through rough air). The material this needs, as a share of the lift, depends on rho V U and the lift per
  volume; size and gravity leave it unchanged.
- *Areal.* The outer cover and the gas cells weigh what their fabric weighs per square metre, set by handling and
  gas-tightness more than by load, so as a share of the lift they fall as 1 / D.

Hence the rule of similarity: in air of the same density, at the same speed and through the same gusts, a hull
scaled up by g_earth / g in every length needs the same share of its lift for its weight-driven and aerodynamic
loads, and a smaller share for its cover and cells.

**Sizing.** `rigid` sizes a Zeppelin-type hull by the groups of the airship weight statements:
- *longitudinals* for the ultimate bending moment at midship, as a thin-walled tube of the hull's radius;
- *frames* (main and intermediate rings, wiring, gangways, reinforcement), in two parts: one in proportion to the
  longitudinals, for the hull's bending and shear, and one in proportion to the ring force the netting brings them
  from the gas pressure;
- *empennage* for the fins' own gust load;
- outer cover and gas cells by area, engines for a top speed, and fixed equipment as a share of the gross lift.
Each structural group's ideal amount carries its own factor; the sky-ship study fits them to the weight statements
of four rigid airships (LZ 129 Hindenburg, LZ 127 Graf Zeppelin, ZR-3 Los Angeles and ZRS-5 Macon).
`pressure_hull` sizes a non-rigid envelope: its pressure must exceed the dynamic pressure by a margin and keep the
hull from wrinkling under the bending moment, and its hoop tension, the pressure at its crown times its radius, must
stay within the fabric's allowable. This is first-order sizing for comparing concepts: joints, fatigue, ground
handling, dynamic response and the details of the fins enter only through the fitted factors.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

GAS_CONSTANT = 8.314462618


@dataclass(frozen=True)
class LiftingGas:
    name: str
    molar_mass_kg_mol: float


HYDROGEN = LiftingGas('hydrogen', 2.01588e-3)
HELIUM = LiftingGas('helium', 4.002602e-3)


@dataclass(frozen=True)
class Air:
    """The air a ship flies in: gravity, density, temperature and the air's molar mass."""
    gravity_m_s2: float
    density_kg_m3: float
    temperature_k: float = 288.15
    molar_mass_kg_mol: float = 0.0289644

    @property
    def pressure_pa(self) -> float:
        return self.density_kg_m3 * GAS_CONSTANT * self.temperature_k / self.molar_mass_kg_mol


def gas_density(air: Air, gas: LiftingGas = HYDROGEN, purity: float = 0.98, superheat_k: float = 0.0):
    """Density of the cell's gas at the air's pressure: the lifting gas at `purity` (by volume), the rest air, at the
    air's temperature plus the superheat."""
    mix = purity * gas.molar_mass_kg_mol / air.molar_mass_kg_mol + (1.0 - purity)
    return air.density_kg_m3 * mix * air.temperature_k / (air.temperature_k + superheat_k)


def net_lift_kg_m3(air: Air, gas: LiftingGas = HYDROGEN, purity: float = 0.98, superheat_k: float = 0.0):
    """Gross lift per cubic metre of cell volume, in kg: displaced air less the cell's gas."""
    return air.density_kg_m3 - gas_density(air, gas, purity, superheat_k)


def lamb_coefficients(fineness):
    """Lamb's apparent-mass coefficients (k1 axial, k2 transverse) of a prolate spheroid of this fineness ratio."""
    f = np.asarray(fineness, dtype=float)
    e = np.sqrt(1.0 - 1.0 / f ** 2)
    ln = np.log((1.0 + e) / (1.0 - e))
    a0 = 2.0 * (1.0 - e ** 2) / e ** 3 * (0.5 * ln - e)
    b0 = 1.0 / e ** 2 - (1.0 - e ** 2) / (2.0 * e ** 3) * ln
    return a0 / (2.0 - a0), b0 / (2.0 - b0)


@dataclass(frozen=True)
class Hull:
    """A streamlined body of revolution. The coefficients fit the great rigid airships: V = 0.475 L D^2 gives
    LZ 129's 200,000 m3 from 245 m and 41.2 m, and the area follows a spheroid of the same fineness."""
    length_m: float
    fineness: float = 6.0                 # length / greatest diameter
    volume_coefficient: float = 0.475     # V = c L D^2 (a spheroid has pi/6 = 0.524)
    area_coefficient: float = 2.45        # A = c L D (a spheroid of fineness 4-6 has 2.50-2.53)

    @property
    def diameter_m(self):
        return self.length_m / self.fineness

    @property
    def radius_m(self):
        return 0.5 * self.diameter_m

    @property
    def volume_m3(self):
        return self.volume_coefficient * self.length_m * self.diameter_m ** 2

    @property
    def area_m2(self):
        return self.area_coefficient * self.length_m * self.diameter_m

    @property
    def reference_area_m2(self):
        """Volume to the two-thirds, the reference area of airship drag coefficients."""
        return self.volume_m3 ** (2.0 / 3.0)

    @classmethod
    def of_volume(cls, volume_m3, fineness=6.0, **kw):
        c = kw.get('volume_coefficient', cls.volume_coefficient)
        return cls((np.asarray(volume_m3, dtype=float) * fineness ** 2 / c) ** (1.0 / 3.0), fineness, **kw)


@dataclass(frozen=True)
class Material:
    name: str
    density_kg_m3: float
    allowable_pa: float      # design stress at limit load in girders, with their buckling in it


# The 1930s rigid airships' duralumin (17S / 2017): tensile strength about 400 MPa, yield about 250 MPa; their
# open lattice girders buckled locally well below yield, so the design stress is taken at 140 MPa.
DURALUMIN = Material('duralumin girders (1930s)', 2790.0, 140e6)
# A carbon-fibre composite tube truss: a laminate of about 75 GPa held to the 0.004 design strain common in aerospace
# structure (Tenney and Dexter 1985), at the 1.57 g/cm3 of IM7/8552.
CARBON = Material('carbon-fibre composite tube truss', 1570.0, 300e6)


@dataclass(frozen=True)
class Fabric:
    name: str
    areal_density_kg_m2: float
    strength_n_m: float       # ultimate tensile strength (warp), N per metre of width


@dataclass(frozen=True)
class Design:
    """Design conditions and the non-structural parts."""
    speed_m_s: float = 36.0            # top speed, for strength and power
    gust_m_s: float = 10.7             # the largest transverse gust (35 ft/s, the Guggenheim Institute's maximum)
    gust_moment: float = 0.10          # Woodward's peak of bending moment / ((U/V) q Vol^(2/3) L)
    safety_factor: float = 1.5         # ultimate / limit
    static_bending: float = 0.016      # static bending moment / (gross lift weight x length), one cell deflated
    ring_pressure_share: float = 0.5   # mean gas-pressure head on the rings, as a share of Δρ g D
    girder_factor: float = 1.0         # actual / ideal longitudinals
    frame_factor: float = 0.5          # frames, wiring, gangways and reinforcement: per kg of ideal longitudinals
    frame_gas_factor: float = 7.5      # and per kg of ideal rings for the gas pressure
    fin_factor: float = 20.0           # actual / ideal fin structure
    frame: Material = CARBON
    cover_kg_m2: float = 0.25          # outer cover, on the hull and both faces of the fins
    cell_kg_m2: float = 0.20           # gas cells with their nettings and valves, on the cells' whole area
    cells: int = 16
    fin_area_share: float = 0.25       # area of the four fins / Vol^(2/3)
    fin_aspect_ratio: float = 0.5      # span^2 / area of one fin (long, low fins)
    fin_lift_slope: float = 2.5        # normal force per radian on one fin, with the hull's interference
    fin_thickness: float = 0.10        # depth of the fin's spar box / chord
    drag_coefficient: float = 0.025    # on Vol^(2/3), hull, fins, cars and engines
    propulsive_efficiency: float = 0.7
    power_plant_kg_w: float = 2.0e-3   # engines, propellers and power supply per watt at top speed
    fixed_share: float = 0.05          # controls, electrical, ballast, water recovery, crew spaces: share of gross lift


def midship_moments(hull: Hull, air: Air, lift_kg: float, design: Design) -> dict:
    """Limit bending moments at midship (N m): the gust's, by Woodward's envelope, and the static one with a cell
    deflated; with the dynamic pressure and the gust's angle for the fins."""
    q = 0.5 * air.density_kg_m3 * design.speed_m_s ** 2
    alpha = np.arctan(design.gust_m_s / design.speed_m_s)
    gust = design.gust_moment * design.gust_m_s / design.speed_m_s * q * hull.reference_area_m2 * hull.length_m
    static = design.static_bending * lift_kg * air.gravity_m_s2 * hull.length_m
    return dict(dynamic_pressure_pa=q, gust_angle_rad=alpha, gust=gust, static=static)


def munk_moment(hull: Hull, air: Air, speed_m_s: float, angle_rad: float):
    """Munk's free moment on a hull at an angle of attack in potential flow, q Vol (k2 - k1) sin 2 alpha."""
    k1, k2 = lamb_coefficients(hull.fineness)
    return 0.5 * air.density_kg_m3 * speed_m_s ** 2 * hull.volume_m3 * (k2 - k1) * np.sin(2.0 * angle_rad)


def ideal_structure(hull: Hull, air: Air, lift_per_m3: float, design: Design) -> dict:
    """The ideal amounts (kg) behind each structural group: longitudinal material for the gust's and the static
    moment, the rings' share of the gas pressure, and the fins' spar caps."""
    V, R, L, D = hull.volume_m3, hull.radius_m, hull.length_m, hull.diameter_m
    m = midship_moments(hull, air, lift_per_m3 * V, design)
    sf, rho_m, sigma = design.safety_factor, design.frame.density_kg_m3, design.frame.allowable_pa
    # Longitudinals: a thin tube of radius R with section modulus (area x R / 2); the moment tapers to the ends, so
    # they average 0.6 of the midship section over the length.
    girders = lambda moment: 2.0 * sf * moment / (sigma * R) * L * 0.6 * rho_m
    # Rings: the netting brings the cells' pressure head (a share of Δρ g D) over each bay to the rings, carried as
    # ring force p s R over the ring's circumference; summed over the rings that is 2 pi R^2 L p / sigma.
    head = design.ring_pressure_share * lift_per_m3 * air.gravity_m_s2 * D
    fins, fin_area = fin_structure(hull, m, design)
    return dict(gust=girders(m['gust']), static=girders(m['static']),
                rings=sf * 2.0 * np.pi * R ** 2 * L * head / sigma * rho_m, fins=fins, fin_area_m2=fin_area,
                moments=m)


def rigid(hull: Hull, air: Air, lift_per_m3: float, design: Design = Design()) -> dict:
    """Size a rigid hull and return its masses (kg) by group, its loads and its useful lift. Arrays broadcast
    through the hull.

    The shares of the gross lift split the structure by what sizes it: weight-driven (the static moment's part of
    the longitudinals and frames, and the rings' gas pressure), aerodynamic (the gust's part, and the empennage) and
    areal (cover and cells)."""
    V, A, R, L, D = hull.volume_m3, hull.area_m2, hull.radius_m, hull.length_m, hull.diameter_m
    lift = lift_per_m3 * V
    i = ideal_structure(hull, air, lift_per_m3, design)
    kg, kf, kp = design.girder_factor, design.frame_factor, design.frame_gas_factor
    parts = dict(longitudinals=kg * (i['gust'] + i['static']), frames=kf * (i['gust'] + i['static']) + kp * i['rings'],
                 empennage=design.fin_factor * i['fins'],
                 cover=design.cover_kg_m2 * (A + 2.0 * i['fin_area_m2']),
                 cells=design.cell_kg_m2 * (A + 2.0 * design.cells * np.pi * R ** 2 * 0.75))
    power = power_w(hull, air, design.speed_m_s, design)
    parts.update(power_plant=design.power_plant_kg_w * power, fixed=design.fixed_share * lift)
    empty = sum(parts.values())
    weight_driven = (kg + kf) * i['static'] + kp * i['rings']
    aerodynamic = (kg + kf) * i['gust'] + parts['empennage']
    areal = parts['cover'] + parts['cells']
    return dict(volume_m3=V, area_m2=A, length_m=L, diameter_m=D, gross_lift_kg=lift, parts_kg=parts,
                empty_kg=empty, useful_kg=lift - empty, useful_share=(lift - empty) / lift,
                structure_kg=parts['longitudinals'] + parts['frames'] + parts['empennage'],
                weight_driven_share=weight_driven / lift, aerodynamic_share=aerodynamic / lift, areal_share=areal / lift,
                moments_n_m=i['moments'], power_w=power, gas_head_pa=lift_per_m3 * air.gravity_m_s2 * D)


def fin_structure(hull: Hull, moments: dict, design: Design):
    """Ideal spar caps of the four fins (kg) and their total area: each fin takes the gust's normal force,
    q S a alpha, whose root moment acts at 0.45 of its span on a box `fin_thickness` of its chord deep; the caps taper
    to the tip, averaging half the root section."""
    area = design.fin_area_share * hull.reference_area_m2
    one = area / 4.0
    span = np.sqrt(design.fin_aspect_ratio * one)
    chord = one / span
    force = moments['dynamic_pressure_pa'] * one * design.fin_lift_slope * moments['gust_angle_rad']
    root = design.safety_factor * force * 0.45 * span
    caps = 2.0 * root / (design.frame.allowable_pa * design.fin_thickness * chord) * span * 0.5 * design.frame.density_kg_m3
    return 4.0 * caps, area


def power_w(hull: Hull, air: Air, speed_m_s, design: Design = Design()):
    """Shaft power to fly at a speed: drag on Vol^(2/3) over the propulsive efficiency."""
    q = 0.5 * air.density_kg_m3 * np.asarray(speed_m_s, dtype=float) ** 2
    return q * design.drag_coefficient * hull.reference_area_m2 * speed_m_s / design.propulsive_efficiency


def pressure_hull(hull: Hull, air: Air, lift_per_m3: float, fabric: Fabric, design: Design = Design(),
                  pressure_margin: float = 1.25, fabric_safety: float = 4.0) -> dict:
    """A non-rigid envelope: the pressure it needs, its hoop tension and whether the fabric carries it.

    The pressure is the larger of `pressure_margin` times the dynamic pressure at top speed and what keeps the
    hull's longitudinal tension, Δp R / 2, above the bending stress M / (pi R^2) under the limit moment. The crown,
    a diameter above the bottom, holds the gas head Δρ g D more. The hoop tension there, times the fabric's safety
    factor, must stay within the fabric's strength; the envelope weighs the fabric's areal density, or more where the
    tension needs a stronger fabric of the same kind."""
    R, D = hull.radius_m, hull.diameter_m
    lift = lift_per_m3 * hull.volume_m3
    m = midship_moments(hull, air, lift, design)
    moment = m['gust'] + m['static']
    dp = np.maximum(pressure_margin * m['dynamic_pressure_pa'], 2.0 * moment / (np.pi * R ** 3))
    crown = dp + lift_per_m3 * air.gravity_m_s2 * D
    hoop = crown * R
    needed = hoop * fabric_safety
    areal = fabric.areal_density_kg_m2 * np.maximum(1.0, needed / fabric.strength_n_m)
    return dict(pressure_pa=dp, crown_pressure_pa=crown, hoop_n_m=hoop, fabric_use=needed / fabric.strength_n_m,
                envelope_kg=areal * hull.area_m2, gross_lift_kg=lift)

