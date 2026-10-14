"""Conditional finite-load-path supports for scalloped, prestressed membranes.

Unlike a free cable net, the edge cables here react against four compression
booms and a finite square hub.  Cable equilibrium and retained planar area are
analytic.  Truss mass is a strength/stiffness/buckling sizing screen with finite
chords, webs, posts and joint allowance, not a qualified assembled structure.
All dimensions and material properties are explicit SI inputs.  No new
empirical material properties are supplied by this module.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import lru_cache
from math import atan, ceil, cos, isfinite, pi, sin, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss

from engineering.array_tile import Material


class SupportInfeasible(ValueError):
    """The assumed geometry and self-load do not yield a converged mass screen."""


def _positive(**values):
    if any(not isfinite(value) or value <= 0 for value in values.values()):
        raise ValueError(f"Expected finite positive values: {values}")


def _nonnegative(**values):
    if any(not isfinite(value) or value < 0 for value in values.values()):
        raise ValueError(f"Expected finite nonnegative values: {values}")


def normal_acceleration_coefficients(omega, alpha, radial_body,
                                     lunar_mean_motion_squared: float):
    """Normal rigid-body acceleration minus linear lunar differential gravity.

    Inputs are body-frame angular velocity, angular acceleration and a unit
    outward radial vector, each with final dimension three. Leading dimensions
    must broadcast. The result has final dimension two, giving cx,cy such that
    a_normal(x,y)=cx*x+cy*y for a flat membrane point (x,y,0).

    The expression is the z component of alpha×r + omega×(omega×r) minus
    n²*(3*u*(u·r)-r), with n²=GM/R³. It is a local spherical-gravity linear
    expansion. It excludes translation, radiation-pressure imbalance, flexible
    accelerations and higher lunar gravity terms; it is not an all-load bound.
    """
    if np.ndim(lunar_mean_motion_squared) != 0:
        raise ValueError("lunar_mean_motion_squared must be a finite nonnegative scalar")
    try:
        mean_motion_squared = float(lunar_mean_motion_squared)
    except (TypeError, ValueError) as exc:
        raise ValueError("lunar_mean_motion_squared must be a finite nonnegative scalar") from exc
    _nonnegative(lunar_mean_motion_squared=mean_motion_squared)
    vectors = [np.asarray(value, dtype=float) for value in (omega, alpha, radial_body)]
    if any(value.ndim == 0 or value.shape[-1] != 3 or value.size == 0 or
           not np.isfinite(value).all() for value in vectors):
        raise ValueError("omega, alpha and radial_body must be finite nonempty arrays with final dimension three")
    try:
        angular_velocity, angular_acceleration, radial = np.broadcast_arrays(*vectors)
    except ValueError as exc:
        raise ValueError("omega, alpha and radial_body leading dimensions must broadcast") from exc
    if not np.allclose(np.linalg.norm(radial, axis=-1), 1., rtol=1e-9, atol=1e-10):
        raise ValueError("radial_body must contain unit vectors")
    cx = (-angular_acceleration[..., 1]+angular_velocity[..., 2]*angular_velocity[..., 0]
          -3*mean_motion_squared*radial[..., 2]*radial[..., 0])
    cy = (angular_acceleration[..., 0]+angular_velocity[..., 2]*angular_velocity[..., 1]
          -3*mean_motion_squared*radial[..., 2]*radial[..., 1])
    return np.stack((cx, cy), axis=-1)


def square_normal_acceleration_envelope(side_m: float, coefficients) -> float:
    """Exact spatial maximum for a square; maximum over supplied time samples.

    This returns L/2*max(|cx|+|cy|), since extrema of the linear normal
    acceleration occur at corners. A time history needs separate sampling
    convergence; the function does not bound unobserved times.
    """
    _positive(side_m=side_m)
    values = np.asarray(coefficients, dtype=float)
    if (values.ndim == 0 or values.shape[-1] != 2 or values.size == 0 or
            not np.isfinite(values).all()):
        raise ValueError("coefficients must be finite nonempty arrays with final dimension two")
    return float(side_m/2*np.max(np.sum(np.abs(values), axis=-1)))


def box_normal_acceleration_envelope(half_spans_m, omega, alpha, radial_body,
                                    lunar_mean_motion_squared: float) -> float:
    """Maximum normal acceleration over a centred box and supplied samples.

    Uses the same rigid acceleration and spherical linear gravity gradient as
    normal_acceleration_coefficients. For finite z, the additional coefficient
    is cz=-(omega_x²+omega_y²)-n²*(3*u_z²-1). The spatial maximum of the
    resulting linear field over a box is |cx|*hx+|cy|*hy+|cz|*hz. This is
    an exact maximum of that *linear field*, not a nonlinear gravity or
    flexible-structure bound; unobserved times remain a sampling question.
    """
    half_spans = np.asarray(half_spans_m, dtype=float)
    if (half_spans.shape != (3,) or not np.isfinite(half_spans).all() or
            (half_spans < 0).any()):
        raise ValueError("half_spans_m must be three finite nonnegative values")
    planar = normal_acceleration_coefficients(omega, alpha, radial_body,
                                              lunar_mean_motion_squared)
    velocity, _, radial = np.broadcast_arrays(np.asarray(omega, dtype=float),
                                              np.asarray(alpha, dtype=float),
                                              np.asarray(radial_body, dtype=float))
    cz = (-(velocity[..., 0]**2+velocity[..., 1]**2)-
          float(lunar_mean_motion_squared)*(3*radial[..., 2]**2-1))
    maximum = (np.sum(np.abs(planar)*half_spans[:2], axis=-1)+
               np.abs(cz)*half_spans[2])
    return float(np.max(maximum))


@lru_cache(maxsize=8)
def _quadrature(order):
    return leggauss(order)


def scalloped_geometry(side_m: float, edge_retreat_m: float, *,
                       quadrature_order: int = 64) -> dict:
    """Four inward circular arcs connecting a square's corners.

    For retreat f < (sqrt(2)-1)*L/2, adjacent circular segments do not overlap. An
    isotropic membrane stress resultant T loads each cable normally, giving
    constant cable tension T*R.  This is equilibrium of a *prescribed* uniform
    stress field; compatible cut pattern, attachment strains and deployment
    are not established.  The centred square of side L-2*f lies inside the
    membrane domain, before any shadows from finite supporting members.
    """
    _positive(side_m=side_m, edge_retreat_m=edge_retreat_m)
    if isinstance(quadrature_order, bool) or not isinstance(quadrature_order, int) or quadrature_order < 8:
        raise ValueError("quadrature_order must be an integer of at least eight")
    if edge_retreat_m >= side_m*(sqrt(2)-1)/2:
        raise ValueError("edge retreat makes neighbouring circular segments overlap")
    radius = (side_m**2 + 4*edge_retreat_m**2)/(8*edge_retreat_m)
    theta = 2*atan(2*edge_retreat_m/side_m)
    segment_factor = (theta-sin(theta)*cos(theta) if theta >= 1e-3 else
                      2*theta**3/3-2*theta**5/15+4*theta**7/315)
    segment = radius**2*segment_factor
    retained = side_m**2-4*segment
    nodes, weights = _quadrature(quadrature_order)
    a = side_m/2
    x = a*nodes
    # h=f-R+sqrt(R²-x²), rationalized to avoid loss of significance.
    removed_height = edge_retreat_m-x*x/(radius+np.sqrt(radius*radius-x*x))
    segment_quadrature = float(a*np.dot(weights, removed_height))
    top_x_second = float(a*np.dot(weights, x*x*removed_height))
    top_y_second = float(a*np.dot(weights, a*a*removed_height-
                                  a*removed_height**2+removed_height**3/3))
    area_second = side_m**4/12-2*(top_x_second+top_y_second)
    per_mass = area_second/retained
    membrane_inertia = [[per_mass, 0., 0.], [0., per_mass, 0.], [0., 0., 2*per_mass]]
    return dict(side_m=side_m, edge_retreat_m=edge_retreat_m,
                cable_radius_m=radius, cable_half_angle_rad=theta,
                cable_length_per_edge_m=2*radius*theta,
                lost_area_m2=4*segment, retained_area_m2=retained,
                retained_area_fraction=retained/side_m**2,
                membrane_area_inertia_per_mass_m2=membrane_inertia,
                retained_area_second_moment_x_m4=area_second,
                retained_area_second_moment_y_m4=area_second,
                quadrature_order=quadrature_order,
                retained_area_quadrature_m2=side_m**2-4*segment_quadrature,
                quadrature_area_relative_error=abs(4*(segment-segment_quadrature))/retained,
                guaranteed_centred_square_side_m=side_m-2*edge_retreat_m,
                guaranteed_centred_circle_radius_m=side_m/2-edge_retreat_m,
                area_rule="planar scalloped membrane; finite-member shadows excluded")


def _line_inertia(a, b, mass):
    """Exact inertia of uniform straight line about the tile origin."""
    second = [[mass*((a[i]*a[j]+b[i]*b[j])/3+
                     (a[i]*b[j]+b[i]*a[j])/6) for j in range(3)] for i in range(3)]
    trace = sum(second[i][i] for i in range(3))
    return [[(trace if i == j else 0)-second[i][j] for j in range(3)] for i in range(3)]


def _add_inertia(target, value):
    for i in range(3):
        for j in range(3):
            target[i][j] += value[i][j]


def _tube_area(force, length, radius, material, buckling_factor):
    limit = min(material.admissible_stress_Pa,
                pi**2*material.youngs_modulus_Pa*radius**2/
                (2*length**2*buckling_factor))
    return force/limit


def _truss_screen(length_m, axial_force_N, depth_m, bay_m, chord_radius_m,
                  material, brace_material, allowable_shortening_m,
                  allowable_lateral_deflection_m, tip_force_N, tip_moment_N_m,
                  line_load_N_m, structural_acceleration_m_s2,
                  joint_mass_fraction, buckling_factor, imperfection_fraction,
                  minimum_wall_thickness_m):
    """Square-section truss, conservatively sized with constant member sections.

    Four chords give I=A_chord*d^2. The same root-fixed, tip-free lateral
    boundary is used for bending and Euler buckling (effective length 2L).
    Euler and web shear receive equal
    compliance allocations so Pcrit=1/(1/PE+1/GA) >= buckling_factor*P.
    Web GA=2*E*A_web*b*d^2/l_diagonal^3 has one active tension diagonal on
    each of two faces; all four faces have tension-only X braces. An explicit
    notional lateral shear sizes braces/posts even without external bending.

    Lateral response is a root-fixed cantilever envelope, with beam-column
    amplification bounded by factor/(factor-1). Half the displacement allowance
    goes to bending and half to shear. This does not demonstrate actual hub
    rotational restraint or validate the Engesser-type truss approximation.
    """
    bays = ceil(length_m/bay_m)
    bay = length_m/bays
    diagonal = sqrt(bay**2+depth_m**2)
    rho, modulus = material.density_kg_m3, material.youngs_modulus_Pa
    er, rhor = brace_material.youngs_modulus_Pa, brace_material.density_kg_m3
    amplification = buckling_factor/(buckling_factor-1)
    minimum_area = 2*pi*chord_radius_m*minimum_wall_thickness_m
    mass = 0.
    converged = False
    for iteration in range(300):
        distributed = line_load_N_m+structural_acceleration_m_s2*mass/length_m
        root_force = tip_force_N+distributed*length_m
        root_moment = tip_moment_N_m+tip_force_N*length_m+distributed*length_m**2/2
        max_chord_force = axial_force_N/4+amplification*root_moment/(2*depth_m)
        chord_area = max(
            minimum_area,
            _tube_area(max_chord_force, bay, chord_radius_m, material, buckling_factor),
            8*buckling_factor*axial_force_N*length_m**2/(pi**2*modulus*depth_m**2),
            axial_force_N*length_m/(4*modulus*allowable_shortening_m),
            amplification*(tip_force_N*length_m**3/3+tip_moment_N_m*length_m**2/2+
                           distributed*length_m**4/8)/
            (modulus*depth_m**2*(allowable_lateral_deflection_m/2)))
        shear = amplification*root_force+imperfection_fraction*axial_force_N
        brace_force = shear*diagonal/(2*depth_m)
        required_ga = max(2*buckling_factor*axial_force_N,
                          amplification*(tip_force_N*length_m+distributed*length_m**2/2)/
                          (allowable_lateral_deflection_m/2))
        brace_area = max(brace_force/brace_material.admissible_stress_Pa,
                         required_ga*diagonal**3/(2*er*bay*depth_m**2))
        post_force = shear
        post_area = max(minimum_area,
                        _tube_area(post_force, depth_m, chord_radius_m,
                                   material, buckling_factor))
        components = dict(chords=4*length_m*chord_area*rho,
                          x_braces=8*bays*diagonal*brace_area*rhor,
                          posts=4*(bays+1)*depth_m*post_area*rho)
        components["joints"] = joint_mass_fraction*sum(components.values())
        updated_mass = sum(components.values())
        if not isfinite(updated_mass) or updated_mass > 1e25:
            raise SupportInfeasible("Self-load sizing diverged for this truss geometry")
        if abs(updated_mass-mass) <= 1e-10*max(1., updated_mass):
            mass = updated_mass
            converged = True
            break
        mass = updated_mass
    if not converged:
        raise SupportInfeasible("Self-load sizing did not converge in 300 iterations")
    inertia = chord_area*depth_m**2
    ga = 2*er*brace_area*bay*depth_m**2/diagonal**3
    pe = pi**2*modulus*inertia/(4*length_m**2)
    pcrit = 1/(1/pe+1/ga)
    bending = amplification*(tip_force_N*length_m**3/3+tip_moment_N_m*length_m**2/2+
                             distributed*length_m**4/8)/(modulus*inertia)
    shear_displacement = amplification*(tip_force_N*length_m+distributed*length_m**2/2)/ga
    return dict(mass_kg=mass, components_kg=components, length_m=length_m,
                depth_m=depth_m, bays=bays, bay_m=bay, diagonal_length_m=diagonal,
                chord_area_m2=chord_area, brace_area_m2=brace_area,
                post_area_m2=post_area, chord_radius_m=chord_radius_m,
                chord_thin_wall_ratio=chord_area/(2*pi*chord_radius_m**2),
                post_thin_wall_ratio=post_area/(2*pi*chord_radius_m**2),
                chord_wall_thickness_m=chord_area/(2*pi*chord_radius_m),
                axial_force_N=axial_force_N, maximum_chord_compression_N=max_chord_force,
                axial_shortening_m=axial_force_N*length_m/(4*modulus*chord_area),
                ideal_global_buckling_force_N=pcrit, global_buckling_utilization=axial_force_N/pcrit,
                euler_effective_length_m=2*length_m,
                root_transverse_force_N=root_force, root_bending_moment_N_m=root_moment,
                amplified_root_bending_moment_N_m=amplification*root_moment,
                lateral_bending_displacement_m=bending,
                lateral_shear_displacement_m=shear_displacement,
                lateral_displacement_m=bending+shear_displacement,
                beam_column_amplification_bound=amplification,
                distributed_load_including_self_N_m=distributed,
                self_load_iterations=iteration+1)


def _truss_inertia(start, end, screen):
    length = screen["length_m"]
    axis = [(end[i]-start[i])/length for i in range(3)]
    transverse = [-axis[1], axis[0], 0.]
    depth = screen["depth_m"]
    n = screen["bays"]
    signs = ((-1, -1), (1, -1), (1, 1), (-1, 1))
    def node(k, corner):
        su, sz = signs[corner]
        return [start[i]+axis[i]*length*k/n+su*depth/2*transverse[i]+
                (sz*depth/2 if i == 2 else 0.) for i in range(3)]
    result = [[0.]*3 for _ in range(3)]
    masses = screen["components_kg"]
    for corner in range(4):
        _add_inertia(result, _line_inertia(node(0, corner), node(n, corner), masses["chords"]/4))
    for k in range(n):
        for corner in range(4):
            nxt = (corner+1)%4
            _add_inertia(result, _line_inertia(node(k, corner), node(k+1, nxt), masses["x_braces"]/(8*n)))
            _add_inertia(result, _line_inertia(node(k, nxt), node(k+1, corner), masses["x_braces"]/(8*n)))
    for k in range(n+1):
        for corner in range(4):
            _add_inertia(result, _line_inertia(node(k, corner), node(k, (corner+1)%4),
                                             masses["posts"]/(4*(n+1))))
            point = node(k, corner)
            _add_inertia(result, _line_inertia(point, point, masses["joints"]/(4*(n+1))))
    return result


def scalloped_support(side_m: float, tension_N_m: float, edge_retreat_m: float,
                      compression_material: Material, cable_material: Material, *,
                      boom_depth_m: float = 100., bay_m: float = 50.,
                      chord_radius_m: float = 2., hub_side_m: float = 1000.,
                      allowable_anchor_retreat_m: float = 5.,
                      allowable_lateral_deflection_m: float = 5.,
                      tip_transverse_force_N: float = 0., tip_moment_N_m: float = 0.,
                      membrane_normal_pressure_Pa: float = 0.,
                      distributed_transverse_load_N_m: float = 0.,
                      structural_acceleration_m_s2: float = 0.,
                      joint_mass_fraction: float = .2, buckling_factor: float = 2.,
                      imperfection_fraction: float = .01,
                      minimum_wall_thickness_m: float = 1e-4) -> dict:
    """Circular edge cables + four radial trusses + four finite hub trusses.

    At each outer corner the edge cables require radial compression
    P=sqrt(2)*T*R*(cos(theta)+sin(theta)). Four hub sides each take P/sqrt(2).
    No support is massless. The finite hub is also sized with each radial
    boom's root transverse-force/moment envelope, using a cantilever-member
    screen; this is deliberately conservative member loading, not a derived
    stiffness solution for the assembled hub. Coupled joint moments, overall
    elastic stability, equipment reactions and deployment remain unvalidated.

    Tip load and moment inputs apply to EACH radial boom. Distributed load
    applies per metre of each radial boom; structural acceleration applies
    to its own iterated mass and separately to hub-member mass. Members that
    leave the thin-wall regime are flagged, never silently called qualified.
    Optional membrane normal pressure contributes q*retained_area/4 to each
    tip load, in addition to equipment tip forces. This transfers the full
    membrane resultant through the supports as a conservative cantilever
    member-loading case; it does not solve cable normal displacement or the
    actual distributed edge reactions. Do not also include that same membrane
    load in tip_transverse_force_N.
    """
    _positive(tension_N_m=tension_N_m, boom_depth_m=boom_depth_m, bay_m=bay_m,
              chord_radius_m=chord_radius_m, hub_side_m=hub_side_m,
              allowable_anchor_retreat_m=allowable_anchor_retreat_m,
              allowable_lateral_deflection_m=allowable_lateral_deflection_m)
    _nonnegative(tip_transverse_force_N=tip_transverse_force_N, tip_moment_N_m=tip_moment_N_m,
                 distributed_transverse_load_N_m=distributed_transverse_load_N_m,
                 structural_acceleration_m_s2=structural_acceleration_m_s2,
                 joint_mass_fraction=joint_mass_fraction, imperfection_fraction=imperfection_fraction,
                 minimum_wall_thickness_m=minimum_wall_thickness_m)
    _nonnegative(membrane_normal_pressure_Pa=membrane_normal_pressure_Pa)
    if not isfinite(buckling_factor) or buckling_factor <= 1:
        raise ValueError("buckling_factor must exceed one")
    geometry = scalloped_geometry(side_m, edge_retreat_m)
    length = (side_m-hub_side_m)/sqrt(2)
    if hub_side_m >= side_m or boom_depth_m >= min(length, hub_side_m):
        raise ValueError("Hub or truss depth is inconsistent with the tile")
    if 4*chord_radius_m >= min(bay_m, boom_depth_m):
        raise ValueError("Chord radius is inconsistent with the truss cross section")
    theta, radius = geometry["cable_half_angle_rad"], geometry["cable_radius_m"]
    cable_force = tension_N_m*radius
    cable_area = cable_force/cable_material.admissible_stress_Pa
    cable_length = geometry["cable_length_per_edge_m"]
    cable_mass = 4*cable_length*cable_area*cable_material.density_kg_m3
    boom_force = sqrt(2)*cable_force*(cos(theta)+sin(theta))
    normal_membrane_force = membrane_normal_pressure_Pa*geometry["retained_area_m2"]
    combined_tip_force = tip_transverse_force_N+normal_membrane_force/4
    common = dict(depth_m=boom_depth_m, bay_m=bay_m, chord_radius_m=chord_radius_m,
                  material=compression_material, brace_material=cable_material,
                  allowable_lateral_deflection_m=allowable_lateral_deflection_m,
                  structural_acceleration_m_s2=structural_acceleration_m_s2,
                  joint_mass_fraction=joint_mass_fraction, buckling_factor=buckling_factor,
                  imperfection_fraction=imperfection_fraction,
                  minimum_wall_thickness_m=minimum_wall_thickness_m)
    boom = _truss_screen(length, boom_force,
                         allowable_shortening_m=allowable_anchor_retreat_m/2,
                         tip_force_N=combined_tip_force, tip_moment_N_m=tip_moment_N_m,
                         line_load_N_m=distributed_transverse_load_N_m, **common)
    hub = _truss_screen(hub_side_m, boom_force/sqrt(2),
                       allowable_shortening_m=allowable_anchor_retreat_m/sqrt(2),
                       tip_force_N=boom["root_transverse_force_N"],
                       tip_moment_N_m=boom["root_bending_moment_N_m"],
                       line_load_N_m=0., **common)
    inertia = [[0.]*3 for _ in range(3)]
    corners = ((-1, -1), (1, -1), (1, 1), (-1, 1))
    for index, (sx, sy) in enumerate(corners):
        start = [sx*hub_side_m/2, sy*hub_side_m/2, 0.]
        end = [sx*side_m/2, sy*side_m/2, 0.]
        _add_inertia(inertia, _truss_inertia(start, end, boom))
        nx, ny = corners[(index+1)%4]
        _add_inertia(inertia, _truss_inertia(start, [nx*hub_side_m/2, ny*hub_side_m/2, 0.], hub))
    # Use local sag coordinates and series to avoid subtracting enormous
    # nearly equal circle-centre coordinates for shallow edge retreat.
    if theta < .01:
        x_factor = theta**2/3-theta**4/15+2*theta**6/315-theta**8/2835
        one_minus_sinc = theta**2/6-theta**4/120+theta**6/5040-theta**8/362880
        cos_deficit_squared = theta**4/20-theta**6/168+theta**8/2880
    else:
        x_factor = .5-sin(2*theta)/(4*theta)
        one_minus_sinc = 1-sin(theta)/theta
        cos_deficit_squared = 1.5-2*sin(theta)/theta+sin(2*theta)/(4*theta)
    y0 = side_m/2-edge_retreat_m
    x2 = radius**2*x_factor
    y2 = y0**2+2*y0*radius*one_minus_sinc+radius**2*cos_deficit_squared
    qxy = cable_mass/2*(x2+y2)
    _add_inertia(inertia, [[qxy, 0., 0.], [0., qxy, 0.], [0., 0., 2*qxy]])
    components = dict(edge_cables=cable_mass, radial_booms=4*boom["mass_kg"],
                      finite_hub=4*hub["mass_kg"])
    cable_diameter = sqrt(4*cable_area/pi)
    # depth is the separation of chord centre lines, not the outer material
    # thickness. Tube radii are the mean radii used by the thin-wall screens.
    # Include half a wall thickness and the solid-equivalent brace radius.
    member_radius = max(
        chord_radius_m+member[area_key]/(4*pi*chord_radius_m)
        for member in (boom, hub) for area_key in ("chord_area_m2", "post_area_m2"))
    member_radius = max(member_radius,
                        *(sqrt(member["brace_area_m2"]/pi) for member in (boom, hub)))
    enclosing_width = boom_depth_m+2*member_radius
    # Union area is bounded above by the sum of member-envelope strips.
    structure_shadow = min(geometry["retained_area_m2"],
                           (4*length+4*hub_side_m)*enclosing_width+8*enclosing_width**2)
    cable_shadow = min(geometry["retained_area_m2"], 4*cable_length*cable_diameter)
    shadow_upper = min(geometry["retained_area_m2"], structure_shadow+cable_shadow)
    thin_wall = max(boom["chord_thin_wall_ratio"], hub["chord_thin_wall_ratio"],
                    boom["post_thin_wall_ratio"], hub["post_thin_wall_ratio"])
    return dict(geometry=geometry, mass_kg=sum(components.values()), components_kg=components,
                support_areal_mass_per_gross_kg_m2=sum(components.values())/side_m**2,
                compression_material=asdict(compression_material), cable_material=asdict(cable_material),
                cable_tension_N=cable_force, cable_area_m2=cable_area,
                cable_diameter_m=cable_diameter,
                cable_axial_strain=cable_force/(cable_material.youngs_modulus_Pa*cable_area),
                cable_equilibrium_rule="isotropic membrane normal traction; cable tension T*R",
                membrane_normal_pressure_Pa=membrane_normal_pressure_Pa,
                membrane_normal_resultant_N=normal_membrane_force,
                membrane_normal_force_per_boom_N=normal_membrane_force/4,
                membrane_normal_force_transfer_residual_N=normal_membrane_force-4*(normal_membrane_force/4),
                normal_boundary_rule="full membrane normal resultant charged at four boom tips; actual cable normal displacement and distributed attachment reactions unresolved",
                radial_boom=boom, hub_member=hub, inertia_kg_m2=inertia,
                diagonal_anchor_retreat_m=boom["axial_shortening_m"]+hub["axial_shortening_m"]/sqrt(2),
                projected_support_obscuration_upper_bound_m2=shadow_upper,
                unobscured_membrane_area_lower_bound_m2=geometry["retained_area_m2"]-shadow_upper,
                obscuration_rule="sum of enclosing truss strips plus cable-width strips; overlaps overcounted; opaque support shadows remove transmitted light, not UV protection",
                maximum_member_outer_radius_m=member_radius,
                enclosing_truss_width_m=enclosing_width,
                outside_span_m=max(side_m+boom_depth_m/sqrt(2)+2*member_radius,
                                   hub_side_m+boom_depth_m+2*member_radius,
                                   side_m+cable_diameter),
                half_depth_m=max(boom_depth_m/2+member_radius, cable_diameter/2),
                material_envelope_rule="encloses straight tube/brace material and planar edge cables; chord depth is centre-line separation; joint protrusions and deformations remain unresolved",
                tip_station_positions_m=[[sx*side_m/2, sy*side_m/2, 0.] for sx, sy in corners],
                thin_wall_screen_passes=thin_wall <= .1,
                evidence="analytic planar cable equilibrium; conditional truss member sizing; unvalidated assembled structure",
                unresolved=["compatible uniformly prestressed cut pattern and edge attachments",
                            "tube shell/local buckling and composite orthotropy",
                            "assembled hub restraint, joint moment transfer and global eigenmodes",
                            "out-of-plane cable response, flexible modes and dynamic amplification",
                            "whole-body load equilibrium and equipment attachment reactions",
                            "deployment, fatigue, thermal distortion and damage tolerance"])


def cellular_support(side_m: float, cells_per_side: int, tension_N_m: float,
                     retreat_fraction: float, compression_material: Material,
                     cable_material: Material, **cell_options) -> dict:
    """Area-equivalent inventory of n² separately supported smaller membranes.

    This adds no massless structural connections between cells. Each is an
    independent unit with the supplied hub/boom dimensions. Their individual
    equipment, control and spacing costs must be added by the caller. The
    array inertia below is an instantaneous geometric inventory, not inertia
    of a rigidly coupled 10 km spacecraft. This wrapper provides no rigidity.
    """
    if isinstance(cells_per_side, bool) or not isinstance(cells_per_side, int) or cells_per_side < 1:
        raise ValueError("cells_per_side must be a positive integer")
    _positive(side_m=side_m, retreat_fraction=retreat_fraction)
    count = cells_per_side**2
    cell_side = side_m/cells_per_side
    cell = scalloped_support(cell_side, tension_N_m, cell_side*retreat_fraction,
                             compression_material, cable_material, **cell_options)
    mass = cell["mass_kg"]
    # Sum of squared offsets of an n×n uniform grid, exact in closed form.
    offset_second = side_m**2*(cells_per_side**2-1)/(12*cells_per_side**2)
    tensor = [[count*cell["inertia_kg_m2"][i][j] for j in range(3)] for i in range(3)]
    tensor[0][0] += count*mass*offset_second
    tensor[1][1] += count*mass*offset_second
    tensor[2][2] += 2*count*mass*offset_second
    return dict(cells=count, cells_per_side=cells_per_side, cell_side_m=cell_side,
                mass_kg=count*mass, support_areal_mass_per_gross_kg_m2=count*mass/side_m**2,
                retained_area_m2=count*cell["geometry"]["retained_area_m2"],
                instantaneous_inventory_inertia_kg_m2=tensor, cell=cell,
                evidence="independent-cell structural inventory; no shared control or inter-cell structure supplied")
