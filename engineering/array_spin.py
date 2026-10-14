"""Necessary-condition screen for a ballasted circular spinning membrane.

This is an ideal, axisymmetric continuum with distributed independent rim
ballast. Radial traction has a finite inertial reaction, rather than a massless
fixed edge. The static sag solution assumes its boundary remains planar; edge
warp, spin-flexible dynamics, deployment and the controller are not solved.
It is not an accepted sail or a complete spacecraft mass model.

The supplied command frame follows the membrane normal but excludes its added
spin. Uniform mass and isotropic in-plane elasticity are assumed. No empirical
material properties or astronomical constants are introduced here.
"""
from __future__ import annotations

from math import hypot, isfinite, pi

import numpy as np

from engineering.array_support import normal_acceleration_coefficients


def _positive(**values):
    if any(not isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f"Expected finite positive values: {values}")


def disk_grid_coverage(radius_m: float, pitch_x_m: float, pitch_y_m: float) -> dict:
    """Coplanar identical disks on a rectangular lattice, normal projection.

    Every point is covered iff the cell corner is within a disk: 2R >=
    hypot(px, py). Uniformly shrinking both pitches to meet that condition
    gives the quoted count factor at fixed protected area. This does not
    supply three-dimensional clearance, rolled projections or support shadows.
    """
    _positive(radius_m=radius_m, pitch_x_m=pitch_x_m, pitch_y_m=pitch_y_m)
    required = hypot(pitch_x_m, pitch_y_m)/2
    scale = min(1., radius_m/required)
    cell = pitch_x_m*pitch_y_m
    return dict(radius_m=radius_m, pitch_x_m=pitch_x_m, pitch_y_m=pitch_y_m,
                minimum_covering_radius_m=required,
                original_grid_fully_covered=radius_m >= required,
                pitch_scale_to_cover=scale,
                count_multiplier_at_fixed_protected_area=1/scale**2,
                gross_film_area_per_original_cell=pi*radius_m**2/cell,
                gross_film_area_per_covered_cell=pi*radius_m**2/(cell*scale**2),
                evidence="exact planar rectangular-lattice covering geometry; no formation or clearance model")


def ballasted_disk_screen(*, radius_m: float, areal_mass_kg_m2: float,
                         spin_rad_s: float, solar_pressure_Pa,
                         allowable_sag_m: float, load_bearing_thickness_m: float,
                         admissible_stress_Pa: float, times_s,
                         omega_rad_s, alpha_rad_s2, radial_body,
                         mean_motion_squared_s2, poisson_ratio: float = .3,
                         dynamic_factor: float = 2.) -> dict:
    """Size edge prestress and ideal ballast for a prescribed rolling history.

    `omega` and `alpha` are expressed in the command frame. They exclude
    added constant positive spin s about its z axis. With c the nonspinning
    normal-acceleration coefficients (including subtraction of lunar tidal
    acceleration), the spinning disk coefficients are c+2*s*omega_xy.
    Their maximum over the disk is R*norm(c). Adding full unrelieved solar
    pressure is a deliberately conservative residual-load assumption.

    T_required(t)=dynamic_factor*q(t)*R^2/(4*allowable_sag) uses the uniform-
    tension disk bound. Actual axial spin is v(t)=s+omega_z(t), so independent
    rim ballast obeys Tedge(t)=M_rim*v(t)^2/(2*pi). Size one finite ballast
    mass to satisfy T_required at every sampled time. It has no separately
    counted hoop-force support. For an isotropic membrane,
    Tr(r)=Tedge+A*(R^2-r^2), A=(3+nu)*sigma*v^2/8. The exact axisymmetric
    uniform-load static sag is q/(4*A)*log(1+A*R^2/Tedge). This also bounds
    signed static loads no greater in magnitude under the same planar edge.

    All moving film/substrate/trim mass belongs in `areal_mass_kg_m2`; only
    the effective thickness that actually carries tension belongs in
    `load_bearing_thickness_m`. Their compatibility must be established
    separately. Static stress excludes roll-induced strain and local anchors.
    No equipment, motors, energy storage or assembly margin is included.
    """
    _positive(radius_m=radius_m, areal_mass_kg_m2=areal_mass_kg_m2,
              spin_rad_s=spin_rad_s, allowable_sag_m=allowable_sag_m,
              load_bearing_thickness_m=load_bearing_thickness_m,
              admissible_stress_Pa=admissible_stress_Pa)
    if not isfinite(poisson_ratio) or not 0 <= poisson_ratio < .5:
        raise ValueError("Expected ordinary isotropic plane-stress Poisson ratio in [0, 0.5)")
    if not isfinite(dynamic_factor) or dynamic_factor < 1:
        raise ValueError("dynamic_factor must be finite and at least one")
    t = np.asarray(times_s, dtype=float)
    if t.ndim != 1 or len(t) < 2 or not np.isfinite(t).all() or not (np.diff(t) > 0).all():
        raise ValueError("times_s must be a finite strictly increasing vector with at least two samples")
    w = np.asarray(omega_rad_s, dtype=float)
    a = np.asarray(alpha_rad_s2, dtype=float)
    u = np.asarray(radial_body, dtype=float)
    if any(v.shape != (len(t), 3) or not np.isfinite(v).all() for v in (w, a, u)):
        raise ValueError("omega, alpha and radial_body must be finite arrays of shape (len(times_s), 3)")
    n2 = np.asarray(mean_motion_squared_s2, dtype=float)
    if n2.ndim != 0 or not np.isfinite(n2) or n2 < 0:
        raise ValueError("mean_motion_squared_s2 must be a finite nonnegative scalar")
    n2 = float(n2)
    pressure = np.asarray(solar_pressure_Pa, dtype=float)
    if pressure.ndim == 0:
        pressure = np.full(len(t), float(pressure))
    if pressure.shape != t.shape or not np.isfinite(pressure).all() or (pressure < 0).any():
        raise ValueError("solar_pressure_Pa must be a finite nonnegative scalar or time vector")
    base = normal_acceleration_coefficients(w, a, u, n2)
    coefficients = base+2*spin_rad_s*w[:, :2]
    normal_acceleration = radius_m*np.linalg.norm(coefficients, axis=1)
    load = pressure+areal_mass_kg_m2*normal_acceleration
    q = float(load.max())
    _positive(maximum_effective_pressure_Pa=q)
    actual_spin = spin_rad_s+w[:, 2]
    if (actual_spin <= 0).any():
        raise ValueError("Command plus added spin must remain positive; no zero-spin ballast support is supplied")
    required_tension = dynamic_factor*load*radius_m**2/(4*allowable_sag_m)
    area = pi*radius_m**2
    film_mass = areal_mass_kg_m2*area
    rim_mass = float(np.max(2*pi*required_tension/actual_spin**2))
    tension = rim_mass*actual_spin**2/(2*pi)
    spin_tension_coefficient = (3+poisson_ratio)*areal_mass_kg_m2*actual_spin**2/8
    spin_centre_tension = spin_tension_coefficient*radius_m**2
    peak_tension = float(np.max(tension+spin_centre_tension))
    sag = float(np.max(load/(4*spin_tension_coefficient)*np.log1p(spin_centre_tension/tension)))
    constant_tension_sag = float(np.max(load*radius_m**2/(4*tension)))
    peak_stress = peak_tension/load_bearing_thickness_m
    axial_inertia = (.5*film_mass+rim_mass)*radius_m**2
    transverse_inertia = axial_inertia/2
    diagonal = np.array([transverse_inertia, transverse_inertia, axial_inertia])
    relative_spin_momentum = axial_inertia*spin_rad_s
    actual_spin_momentum = axial_inertia*actual_spin
    momentum = w*diagonal
    momentum[:, 2] += relative_spin_momentum
    external_torque = a*diagonal+np.cross(w, momentum)
    gravity_torque = 3*n2*np.cross(u, u*diagonal)
    nongravitational_torque = external_torque-gravity_torque
    precession_rate = np.linalg.norm(w[:, :2], axis=1)
    _positive(film_mass_kg=film_mass, rim_mass_kg=rim_mass,
              axial_inertia_kg_m2=axial_inertia, peak_stress_Pa=peak_stress,
              exact_axisymmetric_static_sag_m=sag)
    if not np.isfinite(external_torque).all() or not np.isfinite(nongravitational_torque).all():
        raise ValueError("Parameters overflowed the rigid torque diagnostic")
    return dict(
        geometry="filled circular membrane with ideal distributed independent rim ballast",
        radius_m=radius_m, film_area_m2=area, areal_mass_kg_m2=areal_mass_kg_m2,
        spin_rad_s=spin_rad_s, relative_spin_period_s=2*pi/spin_rad_s,
        minimum_actual_axial_spin_rad_s=float(actual_spin.min()),
        maximum_actual_axial_spin_rad_s=float(actual_spin.max()),
        maximum_rim_speed_m_s=float(radius_m*actual_spin.max()), poisson_ratio=poisson_ratio,
        load_bearing_thickness_m=load_bearing_thickness_m,
        admissible_stress_Pa=admissible_stress_Pa,
        allowable_sag_m=allowable_sag_m, selected_dynamic_factor=dynamic_factor,
        peak_unrelieved_solar_pressure_Pa=float(pressure.max()),
        peak_normal_acceleration_m_s2=float(normal_acceleration.max()),
        peak_effective_pressure_Pa=q,
        pressure_peak_time_s=float(t[np.argmax(load)]),
        peak_normal_precession_rate_rad_s=float(precession_rate.max()),
        edge_prestress_N_m=float(tension.max()),
        minimum_edge_prestress_N_m=float(tension.min()),
        edge_prestress_rule="one fixed ballast mass; edge traction follows actual axial spin squared",
        spin_centre_tension_increment_N_m=float(spin_centre_tension.max()),
        peak_tension_resultant_N_m=peak_tension,
        peak_membrane_stress_Pa=peak_stress,
        stress_utilization=peak_stress/admissible_stress_Pa,
        static_stress_screen_passes=peak_stress <= admissible_stress_Pa,
        exact_axisymmetric_static_sag_m=sag,
        constant_tension_static_sag_bound_m=constant_tension_sag,
        selected_multiplier_times_static_sag_m=dynamic_factor*sag,
        axisymmetric_edge_slope_rad=float(np.max(load*radius_m/(2*tension))),
        small_slope_screen_passes=bool(np.max(load*radius_m/(2*tension)) <= .1),
        film_mass_kg=film_mass, ideal_rim_ballast_mass_kg=rim_mass,
        ideal_film_plus_ballast_mass_kg=film_mass+rim_mass,
        ideal_film_plus_ballast_areal_mass_kg_m2=(film_mass+rim_mass)/area,
        rim_line_density_kg_m=rim_mass/(2*pi*radius_m),
        locked_axisymmetric_inertia_kg_m2=np.diag(diagonal).tolist(),
        spin_angular_momentum_N_m_s=float(actual_spin_momentum.max()),
        minimum_spin_angular_momentum_N_m_s=float(actual_spin_momentum.min()),
        spin_kinetic_energy_J=float(.5*axial_inertia*np.max(actual_spin**2)),
        peak_spin_precession_torque_N_m=float(np.max(actual_spin_momentum*precession_rate)),
        peak_full_rigid_external_torque_N_m=float(np.linalg.norm(external_torque, axis=1).max()),
        peak_full_rigid_nongravitational_torque_N_m=float(np.linalg.norm(nongravitational_torque, axis=1).max()),
        peak_axis_nongravitational_torque_N_m=np.max(np.abs(nongravitational_torque), axis=0).tolist(),
        samples=len(t), duration_s=float(t[-1]-t[0]),
        accepted_complete_design=False,
        evidence="necessary-condition mass, prestress, static shape and rigid torque diagnostic; no validated operating sail",
        reading_rule="Compare ideal film-plus-ballast masses only to equally incomplete ledgers; a passing static stress screen is not a feasible spacecraft or continuously covered formation.",
        unresolved=[
            "planar rim and membrane shape under nonaxisymmetric rolling loads",
            "spin-flexible modes, Coriolis effects and resonant forcing; selected multiplier is not a dynamic bound",
            "rim attachment mass, local stress, nonuniform ballast and film/substrate strain compatibility",
            "external torque allocation, periodic control and loss-of-spin behavior",
            "motors, power, energy storage, hub, deployment hardware and manufacturing margin",
            "optical spectrum, projection, overlap and three-dimensional neighbour clearance"],
    )
