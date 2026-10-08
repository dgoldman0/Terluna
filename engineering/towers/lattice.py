"""A tapered lattice tower sized for gravity and wind on any world: how its mass, base and stiffness grow with height.

The tower is a square space frame. Its four corner legs carry the weight and the overturning moment; bracing in
its four faces carries the shear. Floors are enclosed bands of storeys spaced up the frame, each with the same
floor area, and the open frame between them lets the wind through. The frame's width follows

    b(z) = b_top + (b_base - b_top) (1 - z/H)^p,

so p = 1 is a straight taper and larger p flares the base. Marching down from the top, each leg takes a quarter
of the weight above and the moment of the wind above (wind along a diagonal, the worst case), and its area is set
so that its stress stays at the allowable stress; the bracing carries the wind's shear, or a set fraction of the
legs' mass if that is more. From the sized frame the model finds

- the sway under a wind, which is proportional to the wind speed squared for a given frame;
- the first natural frequency, by Rayleigh's method with the static wind deflection as the mode shape;
- the factor by which the frame's gravity load could grow before the frame buckles as a whole, by the
  Rayleigh-Ritz energy method with the same shape;
- the load and bearing area of each footing.

This is first-order sizing of the kind used to compare concepts. Member buckling enters as one reduction factor,
the dynamic response to gusts as one amplification factor on the peak-gust pressure, and joints, fatigue,
construction, ice and the ground's own mechanics are left out.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

STRUT_SHEAR_FACTOR = 5.2   # bracing mass per unit height = factor x shear x density / allowable stress (X-bracing, 30% for struts)


@dataclass(frozen=True)
class Material:
    name: str
    density_kg_m3: float
    strength_pa: float      # yield stress, or design compressive strength for a composite
    modulus_pa: float


# EN 10025-6 S690: yield 690 MPa for thin plate; E and density are the usual structural-steel values.
STEEL = Material('high-strength structural steel (S690)', 7850.0, 690e6, 210e9)
# A unidirectional carbon-fibre composite in compression, roughly: longitudinal compressive strength about 1 GPa.
CARBON = Material('carbon-fibre composite (unidirectional, in compression)', 1600.0, 1.0e9, 130e9)


@dataclass(frozen=True)
class Frame:
    material: Material = STEEL
    buckling_factor: float = 0.75   # members buckle between bracing points at this share of the strength
    safety_factor: float = 1.6      # on stress, for gravity and wind together
    bracing_fraction: float = 0.35  # bracing mass as a share of the legs', unless the shear needs more
    solidity: float = 0.12          # projected member area / gross area of one face
    round_members: bool = True
    top_width_m: float = 60.0
    top_mass_kg: float = 0.0        # a crown at the top: spire, beacon or docking head
    min_leg_area_m2: float = 0.01   # the smallest leg section, where the loads alone would need less

    @property
    def allowable_pa(self) -> float:
        return self.material.strength_pa * self.buckling_factor / self.safety_factor

    def force_coefficient(self) -> float:
        """Force coefficient of a square lattice tower on the projected member area of one face, for wind along a
        diagonal (TIA-222-G, section 2.6.9.1.1 and Table 2-6): Cf = 4.0 e^2 - 5.9 e + 4.0 from the solidity e, times
        0.57 - 0.14 e + 0.86 e^2 - 0.24 e^3 (at most 1) for round members in subcritical flow, the larger of the
        standard's two round-member factors, times 1 + 0.75 e (at most 1.2) for the diagonal."""
        e = self.solidity
        cf = 4.0 * e * e - 5.9 * e + 4.0
        if self.round_members:
            cf *= min(1.0, 0.57 - 0.14 * e + 0.86 * e * e - 0.24 * e ** 3)
        return cf * min(1.2, 1.0 + 0.75 * e)


@dataclass(frozen=True)
class Floors:
    """Enclosed layers, a band of storeys every `spacing_m`, set within the frame: each storey has the programme's
    floor area, or `plan_share` of the frame's plan where the frame is too narrow for it. The band's facade spans
    the frame's width and takes the wind."""
    spacing_m: float = 200.0
    storeys: int = 4
    storey_height_m: float = 4.0
    area_per_storey_m2: float = 50_000.0
    plan_share: float = 0.5
    mass_per_area_kg_m2: float = 500.0   # floor structure, fit-out, facade share and occupants, per m2 of floor
    drag_coefficient: float = 1.4        # of the band's facade, on its width times its height

    @property
    def band_height_m(self) -> float:
        return self.storeys * self.storey_height_m

    def area_per_m(self, width):
        """Floor area per metre of height in a frame of the given width."""
        return self.storeys * np.minimum(self.area_per_storey_m2, self.plan_share * np.asarray(width) ** 2) / self.spacing_m

    @property
    def full_width_m(self) -> float:
        """The narrowest frame that holds the full programme."""
        return float(np.sqrt(self.area_per_storey_m2 / self.plan_share))


@dataclass(frozen=True)
class World:
    name: str
    gravity_m_s2: float
    air_density: object          # callable: kg/m3 at a height (m) above the tower's base


def size(height_m, base_width_m, exponent, design_wind_m_s, world: World, frame: Frame = Frame(),
         floors: Floors | None = Floors(), dynamic_factor: float = 1.2, steps: int = 600,
         added_drag: float = 0.0) -> dict:
    """Size towers (arrays broadcast together) and return their masses, stiffness and checks.

    design_wind_m_s is the peak gust used for strength, taken as uniform up the tower; the pressure is
    0.5 rho V^2 times dynamic_factor. Sway at any other wind V is sway_per_wind2_m * V^2. added_drag is a
    drag coefficient times the share of the frame's face it covers (wind devices in its openings, or
    cladding), added to the open frame and the floor bands.
    """
    H, b0, p, V = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in (height_m, base_width_m, exponent,
                                                                             design_wind_m_s)))
    H, b0, p, V = (a.ravel() for a in (H, b0, p, V))
    n = H.size
    bt = np.minimum(frame.top_width_m, b0)
    dz = H / steps
    sig, rho_s, E, g = frame.allowable_pa, frame.material.density_kg_m3, frame.material.modulus_pa, world.gravity_m_s2
    cf_phi = frame.force_coefficient() * frame.solidity
    band = (0.0 if floors is None else floors.drag_coefficient * floors.band_height_m / floors.spacing_m) + added_drag

    keep = {k: np.empty((steps, n)) for k in ('z', 'b', 'area', 'm', 'w1', 'm1', 'N')}
    mass_above = np.full(n, frame.top_mass_kg)
    shear1 = np.zeros(n)     # wind shear and moment per (m/s)^2 of design wind
    moment1 = np.zeros(n)
    frame_mass = np.zeros(n)
    floor_mass = np.zeros(n)
    floor_area = np.zeros(n)
    for k in range(steps):
        s = (k + 0.5) / steps
        z = H * (1.0 - s)
        b = bt + (b0 - bt) * s ** p
        db_dz = (b0 - bt) * p * s ** (p - 1.0) / H
        area_floor = np.zeros(n) if floors is None else floors.area_per_m(b)
        m_floor = area_floor * (0.0 if floors is None else floors.mass_per_area_kg_m2)
        w1 = 0.5 * world.air_density(z) * dynamic_factor * (cf_phi + band) * b
        shear_new = shear1 + w1 * dz
        moment_new = moment1 + shear1 * dz + 0.5 * w1 * dz * dz
        lean = np.sqrt(1.0 + 0.5 * db_dz ** 2)   # a corner leg moves (db/2) sqrt(2) sideways per unit height
        # Leg force at the bottom of this step, counting the step's own frame weight implicitly.
        known = g * (mass_above + m_floor * dz) / 4.0 + moment_new * V * V / (np.sqrt(2.0) * b)
        own = g * rho_s * lean * (1.0 + frame.bracing_fraction) * dz / sig
        area = np.maximum(known / sig / (1.0 - own), frame.min_leg_area_m2)
        leg_mass = 4.0 * area * rho_s * lean
        brace_mass = np.maximum(frame.bracing_fraction * leg_mass, STRUT_SHEAR_FACTOR * shear_new * V * V * rho_s / sig)
        m_frame = leg_mass + brace_mass
        mass_above = mass_above + (m_frame + m_floor) * dz
        frame_mass += m_frame * dz
        floor_mass += m_floor * dz
        floor_area += area_floor * dz
        for key, val in (('z', z), ('b', b), ('area', area), ('m', m_frame + m_floor), ('w1', w1), ('m1', moment_new),
                         ('N', g * mass_above)):
            keep[key][k] = val
        shear1, moment1 = shear_new, moment_new

    up = {k: v[::-1] for k, v in keep.items()}
    sway, omega2, buckling = response(dz, E * up['area'] * up['b'] ** 2, up['m'], up['w1'], up['m1'], up['N'])
    base_leg = keep['N'][-1] / 4.0 + keep['m1'][-1] * V * V / (np.sqrt(2.0) * b0)
    width_2_3 = bt + (b0 - bt) * (1.0 / 3.0) ** p
    return dict(height_m=H, base_width_m=b0, exponent=p, design_wind_m_s=V, frame_mass_kg=frame_mass,
                floor_mass_kg=floor_mass, floor_area_m2=floor_area, total_mass_kg=frame_mass + floor_mass + frame.top_mass_kg,
                sway_per_wind2_m=sway, frequency_hz=np.sqrt(omega2) / (2.0 * np.pi), buckling_factor=buckling,
                base_leg_force_n=base_leg, base_shear_n=shear1 * V * V, base_moment_nm=moment1 * V * V,
                width_two_thirds_m=width_2_3, frontal_area_m2=(up['b'] * dz).sum(axis=0))


def response(dz, EI, mass_per_m, load_per_m, moment, axial):
    """Sway, Rayleigh frequency squared and whole-frame buckling factor of a cantilever, from arrays running from
    the base up (steps x designs): bending stiffness, mass and lateral load per metre, the bending moment of that
    load, and the axial compression from gravity. Sway is the top deflection under the load; the frequency uses
    the static deflection as the mode shape, omega^2 = sum(w u) / sum(m u^2); the buckling factor is the
    Rayleigh-Ritz ratio sum(EI u''^2) / sum(N u'^2) for the same shape. Both are NaN where there is no lateral
    load to give the shape (still air)."""
    kappa = moment / EI
    theta = np.cumsum(kappa * dz, axis=0) - 0.5 * kappa * dz
    sway = np.cumsum(theta * dz, axis=0) - 0.5 * theta * dz
    with np.errstate(invalid='ignore', divide='ignore'):
        omega2 = (load_per_m * sway).sum(axis=0) / (mass_per_m * sway ** 2).sum(axis=0)
        buckling = (EI * kappa ** 2).sum(axis=0) / (axial * theta ** 2).sum(axis=0)
    return sway[-1] + 0.5 * theta[-1] * dz, omega2, buckling


def lightest(result: dict, groups, service_wind_m_s, sway_limit=1.0 / 500.0, min_buckling=3.0,
             max_base_share=0.5) -> list:
    """For each group of designs (an array of labels), the lightest frame meeting the checks: sway under the
    service wind within sway_limit x height, whole-frame buckling factor at least min_buckling, and a base no
    wider than max_base_share x height. Returns (label, index or None) pairs."""
    groups = np.asarray(groups)
    ok = ((result['sway_per_wind2_m'] * service_wind_m_s ** 2 <= sway_limit * result['height_m'])
          & (result['buckling_factor'] >= min_buckling)
          & (result['base_width_m'] <= max_base_share * result['height_m']))
    out = []
    for label in dict.fromkeys(groups.tolist()):
        idx = np.flatnonzero((groups == label) & ok)
        out.append((label, None if idx.size == 0 else int(idx[np.argmin(result['frame_mass_kg'][idx])])))
    return out
