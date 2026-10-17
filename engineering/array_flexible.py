"""Conditional energy-consistent normal dynamics for a supported membrane.

This new reduced model leaves the earlier support inventory unchanged. It
discretizes prescribed, uniform prestress with linear triangles, tensioned
edge cables and finite normal support compliances. It is a small-slope energy
model, not a nonlinear membrane/boom finite-element qualification. In-plane
strain, wrinkling, changing prestress, thermal distortion and failure are absent.

Three affine normal motions (translation and two infinitesimal tilts) have
exactly zero elastic energy. All elastic modes are mass orthogonal to them.
The returned normal modes may be used in a separate floating-frame dynamics
model; this module does not supply or approximate its inertial coupling terms.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt

import numpy as np
from scipy import linalg
from scipy.sparse import coo_matrix, csr_matrix, eye, hstack
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class FlexibleModel:
    nodes_m: np.ndarray
    triangles: np.ndarray
    nodal_masses_kg: np.ndarray
    modes: np.ndarray
    omega_squared_s2: np.ndarray
    stiffness_N_m: csr_matrix
    modal_mass: np.ndarray
    modal_stiffness: np.ndarray
    membrane_node_count: int
    corner_indices: np.ndarray
    hub_indices: np.ndarray
    rigid_normal_modes: np.ndarray
    metadata: dict

    @property
    def scalar_modes(self):
        return self.modes[:, 2, :]

    def displaced_nodes(self, modal_coordinates):
        values = np.asarray(modal_coordinates, dtype=float)
        if values.shape != self.omega_squared_s2.shape or not np.isfinite(values).all():
            raise ValueError("modal_coordinates must match the finite retained modal vector")
        return self.nodes_m+np.einsum('nik,k->ni', self.modes, values)


def _positive(**values):
    if any(not isfinite(value) or value <= 0 for value in values.values()):
        raise ValueError(f"Expected finite positive values: {values}")


def _divisions(value):
    if isinstance(value, bool) or not isinstance(value, int) or not 2 <= value <= 40:
        raise ValueError("mesh_divisions must be an integer in [2, 40] for the bounded dense modal solve")


def scalloped_mesh(side_m: float, edge_retreat_m: float, mesh_divisions: int):
    """Coons map from a structured square to four circular boundary arcs.

    Nodes run first along x, then y. Triangles are counterclockwise from +z.
    Returned corners and boundary loops start lower-left and proceed CCW.
    The zero-retreat square is supplied for independent finite-element checks.
    Positive triangle areas are verified rather than assumed from the mapping.
    """
    _positive(side_m=side_m)
    _divisions(mesh_divisions)
    if not isfinite(edge_retreat_m) or not 0 <= edge_retreat_m < side_m*(sqrt(2)-1)/2:
        raise ValueError("edge_retreat_m must preserve nonoverlapping circular scallops")
    n = mesh_divisions
    s, t = np.meshgrid(np.linspace(0., 1., n+1), np.linspace(0., 1., n+1), indexing='ij')
    if edge_retreat_m:
        radius = (side_m**2+4*edge_retreat_m**2)/(8*edge_retreat_m)
        def height(u):
            x = side_m*(u-.5)
            return edge_retreat_m-x*x/(radius+np.sqrt(radius*radius-x*x))
        x = side_m*(s-.5)+(1-2*s)*height(t)
        y = side_m*(t-.5)+(1-2*t)*height(s)
    else:
        x, y = side_m*(s-.5), side_m*(t-.5)
    nodes = np.column_stack((x.ravel(), y.ravel(), np.zeros((n+1)**2)))
    triangles = []
    index = lambda i, j: i*(n+1)+j
    for i in range(n):
        for j in range(n):
            a, b, c, d = index(i, j), index(i+1, j), index(i+1, j+1), index(i, j+1)
            triangles.extend(((a, b, c), (a, c, d)) if (i+j)%2 == 0 else ((a, b, d), (b, c, d)))
    triangles = np.asarray(triangles, dtype=int)
    corners = np.asarray([index(0, 0), index(n, 0), index(n, n), index(0, n)], dtype=int)
    boundaries = [np.asarray([index(i, 0) for i in range(n+1)]),
                  np.asarray([index(n, j) for j in range(n+1)]),
                  np.asarray([index(i, n) for i in range(n, -1, -1)]),
                  np.asarray([index(0, j) for j in range(n, -1, -1)])]
    triangle_areas(nodes, triangles)
    return nodes, triangles, corners, boundaries


def triangle_areas(nodes, triangles):
    xy = np.asarray(nodes, dtype=float)[np.asarray(triangles), :2]
    first, second = xy[:, 1]-xy[:, 0], xy[:, 2]-xy[:, 0]
    area = .5*(first[:, 0]*second[:, 1]-first[:, 1]*second[:, 0])
    if not np.isfinite(area).all() or (area <= 0).any():
        raise ValueError("Mesh contains degenerate or inverted triangles")
    return area


def _membrane_matrices(nodes, triangles, tension_N_m, areal_mass_kg_m2):
    area = triangle_areas(nodes, triangles)
    xy = nodes[triangles, :2]
    gx = np.column_stack((xy[:, 1, 1]-xy[:, 2, 1], xy[:, 2, 1]-xy[:, 0, 1], xy[:, 0, 1]-xy[:, 1, 1]))/(2*area[:, None])
    gy = np.column_stack((xy[:, 2, 0]-xy[:, 1, 0], xy[:, 0, 0]-xy[:, 2, 0], xy[:, 1, 0]-xy[:, 0, 0]))/(2*area[:, None])
    element = tension_N_m*area[:, None, None]*(gx[:, :, None]*gx[:, None, :]+gy[:, :, None]*gy[:, None, :])
    rows = np.repeat(triangles, 3, axis=1).ravel()
    cols = np.tile(triangles, (1, 3)).ravel()
    stiffness = coo_matrix((element.ravel(), (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()
    mass = np.bincount(triangles.ravel(), weights=np.repeat(area*areal_mass_kg_m2/3, 3), minlength=len(nodes))
    return stiffness, mass, area


def _tip_compliance(member, compression_material, cable_material):
    """Finite inherited cantilever bending + web shear compliance.

    Includes the prior model's declared beam-column amplification. This is
    an imposed effective stiffness, not a newly derived assembled-structure
    tangent including geometric stiffness and hub joint constraints.
    """
    length, depth = member['length_m'], member['depth_m']
    modulus = compression_material['youngs_modulus_Pa']
    brace_modulus = cable_material['youngs_modulus_Pa']
    second_moment = member['chord_area_m2']*depth**2
    ga = (2*brace_modulus*member['brace_area_m2']*member['bay_m']*depth**2/
          member['diagonal_length_m']**3)
    amplification = member['beam_column_amplification_bound']
    compliance = amplification*(length**3/(3*modulus*second_moment)+length/ga)
    _positive(tip_compliance_m_N=compliance)
    return compliance


def _elastic_modes(stiffness, masses, nodes, retained_modes):
    n = len(masses)
    if isinstance(retained_modes, bool) or not isinstance(retained_modes, int) or not 1 <= retained_modes <= n-3:
        raise ValueError("retained_modes must lie between one and node count minus three")
    root_mass = np.sqrt(masses)
    inverse_root = 1/root_mass
    weighted = stiffness.multiply(inverse_root[:, None]).multiply(inverse_root[None, :]).toarray()
    weighted = (weighted+weighted.T)/2
    raw_values, raw_vectors = linalg.eigh(weighted, subset_by_index=(0, retained_modes+2), check_finite=False)
    scale = max(float(np.linalg.norm(weighted, ord=np.inf)), 1e-300)
    null_tolerance = scale*1e-10
    if np.max(np.abs(raw_values[:3])) > null_tolerance or raw_values[3] <= null_tolerance:
        raise ValueError("Model does not have exactly three separated affine normal rigid modes")
    rigid = np.column_stack((np.ones(n), nodes[:, 0], nodes[:, 1]))
    orthogonal_rigid, _ = np.linalg.qr(root_mass[:, None]*rigid, mode='reduced')
    elastic = raw_vectors[:, 3:]
    elastic -= orthogonal_rigid@(orthogonal_rigid.T@elastic)
    elastic, _ = np.linalg.qr(elastic, mode='reduced')
    # Small Rayleigh-Ritz re-diagonalization after exact rigid-mode removal.
    reduced = elastic.T@weighted@elastic
    values, vectors = linalg.eigh((reduced+reduced.T)/2, check_finite=False)
    shapes = inverse_root[:, None]*(elastic@vectors)
    mass_matrix = shapes.T@(masses[:, None]*shapes)
    stiffness_matrix = shapes.T@(stiffness@shapes)
    residual = stiffness@shapes-masses[:, None]*shapes*values
    relative_residual = np.linalg.norm(residual)/max(np.linalg.norm(stiffness@shapes), 1e-300)
    return values, shapes, rigid, mass_matrix, stiffness_matrix, dict(
        removed_rigid_normal_modes=3, raw_null_eigenvalues_s2=raw_values[:3].tolist(),
        generalized_eigen_residual_relative=float(relative_residual),
        mass_orthonormal_residual=float(np.linalg.norm(mass_matrix-np.eye(retained_modes))),
        rigid_mass_orthogonality_residual=float(np.linalg.norm(orthogonal_rigid.T@(root_mass[:, None]*shapes))))


def model_from_support(support: dict, sheet_areal_mass_kg_m2: float, *,
                       mesh_divisions: int = 12, retained_modes: int = 12,
                       hub_equipment_mass_kg: float = 0.) -> FlexibleModel:
    """Construct a free, finite-compliance normal model from a support product.

    Let h be the four hub displacements and A h the best-fit hub plane on
    the membrane. Membrane/cable energy uses r=w-Ah. A radial boom's residual
    is its outer displacement minus the extrapolated hub plane and its own
    hub corner's twist displacement. Hub twist energy is a finite quadratic
    penalty on (I-P)h, where P projects hub values onto affine planes.

    The hub twist coefficient is four times one inherited member's reciprocal
    tip compliance, an explicit effective ring restraint assumption. It has
    not been calibrated to a beam/frame model. Boom and hub support masses
    are endpoint lumped; mass and centre of mass are preserved, while detailed
    distributed-frame inertia is not. Hub equipment is split equally among
    the four finite hub points, not treated as a massless attitude anchor.

    Modes are mass normalized: modal coordinates have units m*sqrt(kg), and
    physical displacement is modes @ coordinates. The normal-only kinematics
    cannot resolve in-plane cable/film motions or changing prestress.
    """
    _positive(sheet_areal_mass_kg_m2=sheet_areal_mass_kg_m2)
    if not isfinite(hub_equipment_mass_kg) or hub_equipment_mass_kg < 0:
        raise ValueError("hub_equipment_mass_kg must be finite and nonnegative")
    geometry = support['geometry']
    side, retreat = geometry['side_m'], geometry['edge_retreat_m']
    membrane_nodes, triangles, corners, boundaries = scalloped_mesh(side, retreat, mesh_divisions)
    nm = len(membrane_nodes)
    hub_side = support['hub_member']['length_m']
    hub_xy = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])*hub_side/2
    hub_nodes = np.column_stack((hub_xy, np.zeros(4)))
    nodes = np.vstack((membrane_nodes, hub_nodes))
    hub_indices = np.arange(nm, nm+4)
    tension = support['cable_tension_N']/geometry['cable_radius_m']
    _positive(tension_N_m=tension)
    membrane_k, film_mass, areas = _membrane_matrices(membrane_nodes, triangles, tension, sheet_areal_mass_kg_m2)
    # The polygon tends to the circular boundary. Preserve the producer's
    # exact physical sheet mass while recording the geometry discretization.
    mesh_area = float(areas.sum())
    exact_area = geometry['retained_area_m2']
    film_mass *= exact_area/mesh_area
    membrane_affine = np.column_stack((np.ones(nm), membrane_nodes[:, :2]/side))
    hub_affine = np.column_stack((np.ones(4), hub_xy/side))
    hub_fit = np.linalg.pinv(hub_affine)
    plane_at_membrane = membrane_affine@hub_fit
    plane_at_hub = hub_affine@hub_fit
    twist_projector = np.eye(4)-plane_at_hub
    residual_map = hstack((eye(nm, format='csr'), csr_matrix(-plane_at_membrane)), format='csr')
    # Cable geometric stiffness and lumped mass on polygon edge segments.
    edge_pairs = np.vstack([np.column_stack((edge[:-1], edge[1:])) for edge in boundaries])
    edge_lengths = np.linalg.norm(membrane_nodes[edge_pairs[:, 1]]-membrane_nodes[edge_pairs[:, 0]], axis=1)
    cable_stiffness = support['cable_tension_N']/edge_lengths
    erows = np.column_stack((edge_pairs[:, 0], edge_pairs[:, 0], edge_pairs[:, 1], edge_pairs[:, 1])).ravel()
    ecols = np.column_stack((edge_pairs[:, 0], edge_pairs[:, 1], edge_pairs[:, 0], edge_pairs[:, 1])).ravel()
    evals = (cable_stiffness[:, None]*np.array([1., -1., -1., 1.])).ravel()
    cable_k = coo_matrix((evals, (erows, ecols)), shape=(nm, nm)).tocsr()
    stiffness = (residual_map.T@(membrane_k+cable_k)@residual_map).tocsr()
    masses = np.concatenate((film_mass, np.zeros(4)))
    cable_mass = support['components_kg']['edge_cables']
    cable_segment_mass = cable_mass*edge_lengths/edge_lengths.sum()
    np.add.at(masses, edge_pairs.ravel(), np.repeat(cable_segment_mass/2, 2))
    boom_mass = support['components_kg']['radial_booms']/4
    masses[corners] += boom_mass/2
    masses[hub_indices] += boom_mass/2+support['components_kg']['finite_hub']/4+hub_equipment_mass_kg/4
    boom_compliance = _tip_compliance(support['radial_boom'], support['compression_material'], support['cable_material'])
    hub_compliance = _tip_compliance(support['hub_member'], support['compression_material'], support['cable_material'])
    boom_map = np.zeros((4, nm+4))
    boom_map[np.arange(4), corners] = 1.
    boom_map[:, hub_indices] = -plane_at_membrane[corners]-twist_projector
    hub_map = np.zeros((4, nm+4))
    hub_map[:, hub_indices] = twist_projector
    stiffness += csr_matrix(boom_map.T@boom_map/boom_compliance+4*hub_map.T@hub_map/hub_compliance)
    stiffness = ((stiffness+stiffness.T)/2).tocsr()
    if not np.isfinite(masses).all() or (masses <= 0).any():
        raise ValueError("Every generalized normal node must carry finite positive mass")
    values, shapes, rigid, modal_mass, modal_stiffness, checks = _elastic_modes(stiffness, masses, nodes, retained_modes)
    modes = np.zeros((len(nodes), 3, retained_modes))
    modes[:, 2, :] = shapes
    second_moment = nodes.T@(masses[:, None]*nodes)
    nodal_inertia = np.trace(second_moment)*np.eye(3)-second_moment
    reference_inertia = (np.asarray(support['inertia_kg_m2'])+
                         exact_area*sheet_areal_mass_kg_m2*
                         np.asarray(geometry['membrane_area_inertia_per_mass_m2']))
    reference_inertia += hub_equipment_mass_kg*np.diag([hub_side**2/4, hub_side**2/4, hub_side**2/2])
    metadata = dict(
        schema='terluna.engineering.flexible-array-model/1', mesh_divisions=mesh_divisions,
        retained_modes=retained_modes, membrane_tension_N_m=tension,
        cable_tension_N=support['cable_tension_N'],
        exact_membrane_area_m2=exact_area, polygon_membrane_area_m2=mesh_area,
        area_relative_error=(mesh_area-exact_area)/exact_area,
        membrane_mass_kg=exact_area*sheet_areal_mass_kg_m2,
        edge_cable_mass_kg=cable_mass,
        radial_boom_mass_kg=support['components_kg']['radial_booms'],
        hub_structure_mass_kg=support['components_kg']['finite_hub'],
        hub_equipment_mass_kg=hub_equipment_mass_kg,
        nodal_reference_inertia_kg_m2=nodal_inertia.tolist(),
        distributed_support_and_sheet_reference_inertia_kg_m2=reference_inertia.tolist(),
        nodal_inertia_relative_frobenius_error=float(np.linalg.norm(nodal_inertia-reference_inertia)/np.linalg.norm(reference_inertia)),
        boom_tip_compliance_m_N=boom_compliance, hub_member_tip_compliance_m_N=hub_compliance,
        hub_twist_stiffness_N_m=4/hub_compliance,
        modal_coordinate_units='m sqrt(kg); shapes have units kg^(-1/2)',
        mass_rule='triangle-lumped film rescaled to exact curved area; cable mass distributed by polygon edge length; boom and hub endpoint lumping',
        support_rule='finite effective cantilever compliance and hub twist penalty; no calibrated assembled-frame stiffness',
        geometry_rule='planar circular-arc Coons mesh with piecewise-linear triangle boundary; area refinement required',
        checks=checks,
        evidence='implemented linear energetic reduced model; finite compliant boundary; no empirical structural validation',
        unresolved=['in-plane strain and prestress compatibility', 'wrinkling and nonlinear stress stiffening',
                    'assembled boom/hub modes and equipment attachment flexibility',
                    'endpoint-lumped support inertia versus distributed hardware',
                    'thermal deformation, damping, damage and deployment',
                    'nonlinear finite-slope geometry and follower-load stability'])
    return FlexibleModel(nodes, triangles, masses, modes, values, stiffness,
                         modal_mass, modal_stiffness, nm, corners, hub_indices, rigid, metadata)


def fixed_square_membrane(*, side_m: float, tension_N_m: float,
                          areal_mass_kg_m2: float, mesh_divisions: int,
                          pressure_Pa: float = 0., mode_count: int = 3) -> dict:
    """Independent fixed-edge square benchmark using the same triangle kernel.

    Analytic eigenvalues are (T/sigma)*pi²*(m²+n²)/L². Uniform pressure
    also permits comparison with the independent sine-series sag solution.
    This benchmark does not validate the model's prescribed support stiffness.
    """
    _positive(side_m=side_m, tension_N_m=tension_N_m, areal_mass_kg_m2=areal_mass_kg_m2)
    if not isfinite(pressure_Pa):
        raise ValueError("pressure_Pa must be finite")
    nodes, triangles, _, edges = scalloped_mesh(side_m, 0., mesh_divisions)
    stiffness, masses, areas = _membrane_matrices(nodes, triangles, tension_N_m, areal_mass_kg_m2)
    boundary = np.unique(np.concatenate(edges))
    interior = np.setdiff1d(np.arange(len(nodes)), boundary)
    if isinstance(mode_count, bool) or not isinstance(mode_count, int) or not 1 <= mode_count <= len(interior):
        raise ValueError("mode_count must fit the nonempty square interior")
    reduced_k = stiffness[interior][:, interior]
    inverse_root = 1/np.sqrt(masses[interior])
    weighted = reduced_k.multiply(inverse_root[:, None]).multiply(inverse_root[None, :]).toarray()
    values = linalg.eigh(weighted, subset_by_index=(0, mode_count-1), eigvals_only=True)
    force = pressure_Pa*masses/areal_mass_kg_m2
    displacement = np.zeros(len(nodes))
    displacement[interior] = spsolve(reduced_k, force[interior])
    centre = np.argmin(np.linalg.norm(nodes[:, :2], axis=1))
    return dict(nodes_m=nodes, triangles=triangles, nodal_masses_kg=masses,
                area_m2=float(areas.sum()), omega_squared_s2=values,
                centre_displacement_m=float(displacement[centre]),
                displacement_m=displacement, potential_energy_J=float(.5*displacement@(stiffness@displacement)))
