"""The summit port's chosen form, designed and sized by hand calculation (research/decisions.md, the summit port and
metropolis, 2026-09-27 and 2026-09-28):

    OPENBLAS_NUM_THREADS=1 python -m research.studies.summit_tower.form    # results/summit_tower_form.json

The form: a round diagrid of 24 members each way, gathered into six splayed legs through a deep transfer ring on its
fourth node level, with an arch between each pair of legs; six wide rings in the lower tower whose decks cantilever
from the frame; seven disks from about 12 km up, each a spoked wheel round a central core; the core; and the
long-haul terminal with its docking arms in the flight band. Each system is designed below (its members, spans,
depths and spacings, chosen with their reasons) and then sized for lunar gravity and the study's design gust, with
the study's steel (S690 at 323 MPa allowable) and cables of 1,770 MPa wire at 800 MPa.

Method and limits: first-order design, as at concept stage. Forces come from statics, member areas from the allowable
stress, and member buckling from a slenderness limit that keeps the frame model's buckling factor valid. The ring
cantilevers are pre-cambered for their dead load. The disks are analysed as radial cable trusses without the deck's
hoop stiffness, which overstates the forces near the core. The frame's sway, period and whole-frame buckling use the
lattice model's Rayleigh routine. Joints, fatigue, construction stages, gust dynamics beyond one factor, ice and the
ground's own mechanics are left out. Evidence state: a proposed design sized by hand calculation, not checked with a
structural analysis model.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from atmosphere.radiative_convective.thermodynamics import MOON
from engineering.towers import lattice as lt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STUDY = HERE / 'results' / 'summit_tower.json'
GCM_WINDS = ROOT / 'climate' / 'results' / 'gcm' / 'site_winds_A28_dim5_summit.json'
SCHEMA = 'terluna.research.summit-tower-form/1'

G = MOON.surface_gravity
SIGMA = lt.Frame().allowable_pa          # S690: 690 MPa x 0.75 (member buckling) / 1.6 (safety), as the study
SIGMA_CABLE = 800e6                      # locked-coil cable of 1,770 MPa wire at a factor of 2.2
TAU = SIGMA / math.sqrt(3.0)
RHO = lt.STEEL.density_kg_m3
E = lt.STEEL.modulus_pa
DYNAMIC = 1.2                            # the study's gust amplification on the peak-gust pressure
BEARING_PA = 1.0e6                       # the study's footing pressure

# ---- the frame: the study's outline and the diagrid chosen on 2026-09-27
H, B0, BT, P = 24000.0, 8193.0, 320.0, 1.5
N = 24                                   # members each way
C = math.tan(math.radians(35.0)) / (B0 / 2)
NODE_DZ = math.pi / (N * C)              # node levels (the members cross) every 765.8 m
Z_T = 4 * NODE_DZ                        # the transfer ring
SOLIDITY = 0.12                          # the study's frame solidity: it sets the drag
MEMBER_SHARE = 0.055                     # each member's breadth as a share of the spacing between members
MIN_MEMBER = 6.0                         # m: near the top the solidity rule alone would give 3 m members
SLENDER = 20.0                           # unbraced length at most 20 breadths (slenderness about 40 for four
#                                          chords at the corners), so the 0.75 buckling factor holds; rings between
#                                          the node levels brace the members where the nodes are too far apart
BRACE_SHARE = 0.015                      # a ring's radial bracing force, as a share of each member's force
LEG_FLARE, FOOT_Z = 350.0, -350.0        # the legs flare outward 350 m to feet about 350 m below the centre

# ---- floors: masses per m2
FLOOR_KG_M2 = 500.0                      # a storey: structure, fit-out, facade share and occupants (the study's)
ROOF_SLAB_KG_M2 = 300.0                  # the roof slab under a park, with its waterproofing
PARK_KG_M2 = 1000.0                      # 0.6 m of lightweight soil at 1.2 t/m3 wet (720), drainage and paving
#                                          (150), trees, water, snow and people (130)
DISK_COVER, DISK_STOREYS = 0.30, 3       # buildings cover 30% of a disk's deck, three storeys high
HALLS_KG_M2 = 150.0                      # sealed disks: pressurised glass halls over a third of the deck, spread
SEALED_FROM = 15000.0

# Rings: node level, width (m), storey heights (m). Widest where the open air is mildest; fewer storeys under the
# wide parks. Ring 0 carries the base programme (6 km2) so the ground between the legs stays open.
RINGS = [(4, 150.0, (9.0, 6.0)), (6, 120.0, (6.0, 6.0)), (8, 100.0, (6.0, 6.0)),
         (10, 80.0, (5.0, 5.0, 5.0)), (12, 60.0, (5.0, 5.0, 5.0)), (14, 50.0, (5.0, 5.0, 5.0))]
TRUSS_SPACING = 15.0                     # m between radial cantilever trusses at the frame (two facade bays)
ROOT_DEPTH_RATIO = 1.0 / 6.0             # the truss below the floors: depth at the frame / cantilever length
TIP_DEPTH = 2.0
SECONDARY = 0.30                         # posts, bracing and deck beams, as a share of the chords and diagonals

# Disks: node levels every second node (1,532 m apart), from about 12 km to 21.4 km.
DISKS = [16, 18, 20, 22, 24, 26, 28]
SPOKES = 144                             # radial cable trusses, every 2.5 degrees
LENS = 0.20                              # depth of the cable at the core / the disk's radius
POSTS = 8                                # struts from each cable up to the deck
FILL = 0.90                              # a locked-coil cable's steel share of its circle

# The core: from the lowest disk's lower hub to the terminal; lift shafts, stairs and services in a lattice tube.
CORE_DIAMETER = 120.0
CORE_KG_M = 40e3                         # lattice 12 t/m, 24 shafts with rails and cabs 20 t/m, services 8 t/m

# The terminal and its docking arms in the flight band (the sealed disks hold the rest of the long-haul programme).
TERMINAL_FLOOR_M2 = 1.0e6
ARMS = [(900.0 + 180.0 * (j % 3)) for j in range(8)]
GANGWAY_KG_M = 5e3


def width(z):
    return BT + (B0 - BT) * max(0.0, 1.0 - z / H) ** P


def radius(z):
    return width(z) / 2


def breadth(z):
    return max(MIN_MEMBER, MEMBER_SHARE * 2 * math.pi * radius(z) / N)


def ring_outer(z):
    """A floor's outer edge: just inside the members' node hubs."""
    return radius(z) - 0.65 * breadth(z) - 2.0


def lean(z):
    """A member's angle from the vertical: its helix turns C radians per metre of height at radius r."""
    return math.atan(C * radius(z))


def air_density(gcm: dict, ground_above_sea_m: float):
    z_sea = np.asarray(gcm['height_m']) + gcm['ground_height_m']
    log_rho = np.log(np.asarray(gcm['density_kg_m3']))
    return lambda z: float(np.exp(np.interp(ground_above_sea_m + z, z_sea, log_rho)))


# ------------------------------------------------------------------------------------------------ the wide rings
def cantilever(L, q_pa, spacing=TRUSS_SPACING):
    """A radial truss cantilevered L metres from the frame, carrying q_pa over its strip, its depth tapering from
    L/6 at the frame to 2 m at the tip. Chords and diagonals at the allowable stress; pre-cambered for dead load, so
    the stress-sized chords are enough. Returns the truss's steel per m2 of deck, its root moment and depth."""
    d_root = max(ROOT_DEPTH_RATIO * L, 8.0)
    x = np.linspace(0.0, L, 401)
    depth = d_root + (TIP_DEPTH - d_root) * x / L
    m = q_pa * spacing * (L - x) ** 2 / 2
    v = q_pa * spacing * (L - x)
    chords = 2 * m / (depth * SIGMA)                  # top (tension) and bottom (compression)
    diagonals = 2 * v / SIGMA                         # 45-degree Warren diagonals: area v sqrt2, length sqrt2 per m
    kg = RHO * (1 + SECONDARY) * np.trapezoid(chords + diagonals, x)
    chord_root = float(m[0] / (d_root * SIGMA))
    return dict(steel_kg_m2=kg / (spacing * L), root_moment_nm=float(m[0]), root_depth_m=d_root,
                chord_area_root_m2=chord_root, diagonal_area_root_m2=float(v[0] * math.sqrt(2) / SIGMA))


def root_girder(r_o, q_line, m_line, depth):
    """The box girder along the frame at a ring's root: it carries the deck's reaction q_line (N/m) and its moment
    m_line (N m/m, as torsion) to the frame. It rests on the nodes, on a post from the node below at each mid-span
    (the diamond's lower corner lies under it) and on posts from the diagonal members half a node level below at the
    quarter points, so it spans a quarter of the node spacing. Returns kg per metre of ring and the span."""
    span = 2 * math.pi * r_o / N / 4
    bend = q_line * span ** 2 / 10
    chords = 2 * bend / (depth * SIGMA)
    width_ = depth / 2
    torque = m_line * span / 2
    walls = 2 * (depth + width_) * torque / (2 * depth * width_ * TAU)
    return RHO * (1 + SECONDARY) * (chords + walls), span


def frame_band(r_o, m_line, z):
    """The frame takes each ring's root moment as hoop forces a node level above and below, and its members bend
    to carry it there: the steel added to the tie rings and to the four member segments at each of the 24 nodes."""
    hoop = m_line * r_o / NODE_DZ
    rings = 2 * 2 * math.pi * r_o * hoop / SIGMA * RHO
    per_segment = m_line * 2 * math.pi * r_o / N / 4
    chord = per_segment / breadth(z) / SIGMA
    members = N * 4 * 2 * chord * (NODE_DZ / math.cos(lean(z))) * 0.5 * RHO
    return rings + members


def rings():
    out = []
    for level, W, storeys in RINGS:
        z = level * NODE_DZ
        r_o = ring_outer(z)
        r_c = r_o - W / 2
        area = 2 * math.pi * r_c * W
        dead = len(storeys) * FLOOR_KG_M2 + ROOF_SLAB_KG_M2 + PARK_KG_M2
        truss = cantilever(W, G * (dead + 150.0))
        for _ in range(3):                            # the truss's own weight, to convergence
            truss = cantilever(W, G * (dead + truss['steel_kg_m2']))
        q = G * (dead + truss['steel_kg_m2'])
        girder_kg_m, span = root_girder(r_o, q * W, truss['root_moment_nm'] / TRUSS_SPACING, truss['root_depth_m'])
        band = frame_band(r_o, truss['root_moment_nm'] / TRUSS_SPACING, z)
        structure = truss['steel_kg_m2'] * area + girder_kg_m * 2 * math.pi * r_o
        out.append(dict(
            index=len(out), node_level=level, z_m=z, outer_radius_m=r_o, width_m=W, storey_heights_m=list(storeys),
            floor_m2=area * len(storeys), park_m2=area, dead_kg_m2=dead,
            truss=dict(spacing_m=TRUSS_SPACING, root_depth_m=truss['root_depth_m'], tip_depth_m=TIP_DEPTH,
                       chord_area_root_m2=truss['chord_area_root_m2'],
                       diagonal_area_root_m2=truss['diagonal_area_root_m2'], steel_kg_m2=truss['steel_kg_m2']),
            root_girder=dict(depth_m=truss['root_depth_m'], width_m=truss['root_depth_m'] / 2, span_m=span,
                             steel_kg_m=girder_kg_m),
            mass_kg=(dead * area + structure), structure_kg=structure, frame_band_kg=band,
            sealed=z >= SEALED_FROM))
    return out


# ----------------------------------------------------------------------------------------------------- the disks
def disk(level):
    """A disk as a spoked wheel: the deck is the compression chord, 144 radial cables the tension chord, hanging
    from the rim at the frame down to a tension hub on the core LENS x R below the deck, with posts from each cable
    up to the deck. Without hoop stiffness, each radian carries M(r) = q (R^3 - r^3) / 6, so the cable hangs
    h (1 - (r/R)^3) below the deck and its horizontal force per radian is F = q R^3 / (6 h). The chords' horizontal
    forces balance at the rim and at the hubs, so the frame takes only the deck's weight."""
    z = level * NODE_DZ
    R = ring_outer(z)
    sealed = z >= SEALED_FROM
    dead = PARK_KG_M2 + ROOF_SLAB_KG_M2 + DISK_COVER * DISK_STOREYS * FLOOR_KG_M2 + (HALLS_KG_M2 if sealed else 0.0)
    h = LENS * R
    r_hub = CORE_DIAMETER / 2 + 10.0
    structure_kg_m2 = 250.0
    for _ in range(5):
        q = G * (dead + structure_kg_m2)
        F = q * R ** 3 / (6 * h)
        cable_area = 2 * math.pi * F / SIGMA_CABLE
        rr = np.linspace(r_hub, R, 400)
        sag = h * (1 - (rr / R) ** 3)
        slope = 3 * h * rr ** 2 / R ** 3
        cable_len = float(np.trapezoid(np.sqrt(1 + slope ** 2), rr))
        cables = cable_area * cable_len * RHO
        deck = 2 * math.pi * F * (R - r_hub) / SIGMA * RHO          # the deck's radial compression
        hubs = 2 * (2 * math.pi * r_hub) * F / SIGMA * RHO          # compression ring on the deck, tension ring below
        # posts: each carries its share of the deck down to the cable, as a lattice strut 1/30 of its length wide
        post_r = r_hub + (R - r_hub) * (np.arange(POSTS) + 0.5) / POSTS
        post_len = h * (1 - (post_r / R) ** 3)
        post_load = q * math.pi * (R ** 2 - r_hub ** 2) / (SPOKES * POSTS)
        posts = SPOKES * float(np.sum(post_load / SIGMA * post_len)) * RHO * 2.0   # x2 for buckling and bracing
        rim_kg_m, span = root_girder(R, q * R / 2, 0.0, max(20.0, R / 40))
        rim = rim_kg_m * 2 * math.pi * R
        total = (cables + deck + hubs + posts) * (1 + SECONDARY) + rim
        structure_kg_m2 = total / (math.pi * R ** 2)
    cable_d = math.sqrt(4 * cable_area / SPOKES / (math.pi * FILL))
    return dict(node_level=level, z_m=z, radius_m=R, deck_m2=math.pi * R ** 2, sealed=sealed,
                open_winter=(not sealed) and z > 13000.0, dead_kg_m2=dead, lens_depth_m=h, hub_radius_m=r_hub,
                horizontal_force_per_radian_n=F, spokes=SPOKES, cable_diameter_m=cable_d, cable_length_m=cable_len,
                hub_ring_area_m2=F / SIGMA, posts_per_spoke=POSTS,
                post_radii_m=[float(v) for v in post_r], post_lengths_m=[float(v) for v in post_len],
                rim_girder=dict(depth_m=max(20.0, R / 40), span_m=span, steel_kg_m=rim_kg_m),
                floor_m2=math.pi * R ** 2 * DISK_COVER * DISK_STOREYS,
                structure_kg=total, structure_kg_m2=structure_kg_m2, mass_kg=(dead + structure_kg_m2) * math.pi * R ** 2)


# ------------------------------------------------------------------------------------------------- the frame
def frame(floors: list, core: dict, crown: dict, rho_air, v_design: float) -> dict:
    """March down the diagrid from the top: each member carries a 48th of the weight above and its share of the
    wind's overturning moment (2M / (48 r) at the extreme member), along its helix. Node rings brace the members at
    every node level, and intermediate rings wherever the node levels are more than 20 breadths apart along a
    member. Floors, the core and the terminal hang on the frame at their levels."""
    dz = 10.0
    zs = np.arange(H - dz / 2, Z_T, -dz)
    drag_frame = lt.Frame().force_coefficient() * SOLIDITY
    point = {}                                        # level -> (mass, frontal area x drag coefficient)
    for f in floors:
        point.setdefault(round(f['z_m'] / dz), []).append((f['mass_kg'] + f.get('frame_band_kg', 0.0), f['cd_area']))
    point.setdefault(round(crown['z_m'] / dz), []).append((crown['mass_kg'], crown['cd_area']))
    mass_above, shear, moment = 0.0, 0.0, 0.0
    rows, node_rings, mid_rings = [], [], []
    next_node = math.floor((H - Z_T) / NODE_DZ) * NODE_DZ + Z_T
    member_kg = ring_kg = 0.0
    for z in zs:
        r = radius(z)
        a = lean(z)
        b = width(z)
        for m, cda in point.get(round(z / dz), []):
            mass_above += m
            shear += 0.5 * rho_air(z) * DYNAMIC * v_design ** 2 * cda
        if core['z0_m'] <= z <= core['z1_m']:
            mass_above += CORE_KG_M * dz
        w = 0.5 * rho_air(z) * DYNAMIC * v_design ** 2 * drag_frame * b
        moment += shear * dz + 0.5 * w * dz * dz
        shear += w * dz
        force = (G * mass_above / (2 * N) + 2 * moment / (2 * N * r)) / math.cos(a)
        own = G * RHO * dz / SIGMA                     # the step's own member weight, implicitly
        area = max(force / SIGMA / (1 - own), 0.02)
        seg_kg = 2 * N * area * RHO * dz / math.cos(a)
        member_kg += seg_kg
        mass_above += seg_kg
        if z <= next_node + 1e-6:                      # a node level: its ring, and the rings between it and the next
            ring_area = max(BRACE_SHARE * force * 2 * N / (2 * math.pi) / SIGMA, 0.05)
            kg = 2 * math.pi * r * ring_area * RHO
            node_rings.append(dict(z_m=next_node, radius_m=r, ring_area_m2=ring_area))
            span = NODE_DZ / math.cos(a)
            n_mid = max(0, math.ceil(span / (SLENDER * breadth(z))) - 1)
            for k in range(1, n_mid + 1):
                mid_rings.append(dict(z_m=next_node + k * NODE_DZ / (n_mid + 1), ring_area_m2=ring_area * 0.6))
            kg *= 1 + 0.6 * n_mid
            ring_kg += kg
            mass_above += kg
            next_node -= NODE_DZ
        rows.append(dict(z_m=float(z), radius_m=r, lean_deg=math.degrees(a), breadth_m=breadth(z),
                         member_force_n=force, member_area_m2=area, mass_above_kg=mass_above, moment_nm=moment,
                         shear_n=shear))
    return dict(rows=rows[::-1], member_kg=member_kg, ring_kg=ring_kg, node_rings=node_rings[::-1],
                intermediate_rings=sorted(mid_rings, key=lambda x: x['z_m']), weight_at_transfer_n=G * mass_above,
                moment_at_transfer_nm=moment, shear_at_transfer_n=shear)


def base(fr: dict, rho_air, v_design: float) -> dict:
    """Six legs, six arches, the transfer ring and a buried tie between the feet. The diagrid lands on 24 nodes at
    the transfer ring: six on the legs and three between each pair of legs, which the arch below picks up at its
    crown and its two spandrel posts. So each arch carries three node loads to the feet it springs from, each leg
    one node load and its share of the wind's overturning moment, and the ring ties the legs' tops, which lean
    inward under load. The feet push outward (the legs' lean and the arches' thrust); a tie between the feet along
    the hexagon's edges holds them."""
    W = fr['weight_at_transfer_n']
    node = W / 24
    r_top = radius(Z_T)
    r_foot = radius(0.0) + LEG_FLARE
    leg_len = math.hypot(r_foot - r_top, Z_T - FOOT_Z)
    beta = math.atan2(r_foot - r_top, Z_T - FOOT_Z)
    # the wind below the transfer ring, on the six legs as lattice columns about 500 m wide
    below = 0.5 * rho_air(Z_T / 2) * DYNAMIC * v_design ** 2 * 2.0 * 0.25 * 500 * 6 * (Z_T - FOOT_Z)
    m_base = fr['moment_at_transfer_nm'] + fr['shear_at_transfer_n'] * (Z_T - FOOT_Z) + below * (Z_T - FOOT_Z) / 2
    leg_force = (node + m_base / (3 * (r_top + r_foot) / 2)) / math.cos(beta)
    leg_area = leg_force / SIGMA
    leg_kg = leg_area * leg_len * RHO * 1.35
    leg_panel = 0.45 * 515.0                              # panels about 0.45 of the leg's mean width (650 to 380 m)
    leg_chord_breadth = leg_panel / SLENDER               # each chord a lattice box, braced at every panel point
    # arches: parabolic, from foot to foot (the hexagon's edge) rising to the ring's underside
    span = 2 * r_foot * math.sin(math.pi / 6)
    rise = Z_T - 280.0 - FOOT_Z
    thrust = 3 * node * span / (8 * rise)
    arch_force = math.hypot(thrust, 1.5 * node)          # at the springing; at the crown it is the thrust alone
    arch_area = arch_force / SIGMA
    s_ = 4 * rise / span
    arch_len = span / 2 * (math.sqrt(1 + s_ * s_) + math.asinh(s_) / s_)
    # axial force H / cos(slope) along the parabola: steel = H / sigma x integral of (1 + y'^2) dx
    arch_kg = thrust / SIGMA * span * (1 + 16 * rise ** 2 / (3 * span ** 2)) * RHO * 1.35
    arch_rib = arch_len / 4 / SLENDER                      # braced at the springings, the posts and the crown
    # transfer ring: hoop compression from the legs' inward push at their tops. The node loads land on the legs and
    # on the arch's crown and posts, so the ring's own bending is only ring 0's, which ring 0's root girder carries;
    # a lattice 280 m deep and 120 m wide, with 30% for its bracing and the wind's shear
    inward = leg_force * math.sin(beta)
    ring_force = inward / (2 * math.sin(math.pi / 6))
    ring_area = ring_force / SIGMA
    ring_kg = ring_area * 2 * math.pi * r_top * RHO * 1.30
    # feet: vertical load, footing area, and the outward push on the buried tie
    leg_weight = G * (leg_kg + arch_kg)
    foot_v = W / 6 + leg_weight / 6 + m_base / (3 * r_foot)
    outward = leg_force * math.sin(beta) + thrust
    tie_force = outward / (2 * math.sin(math.pi / 6))
    tie_area = tie_force / SIGMA
    tie_kg = tie_area * span * 6 * RHO
    return dict(weight_at_transfer_n=W, node_load_n=node, base_moment_nm=m_base,
                leg=dict(count=6, length_m=leg_len, lean_deg=math.degrees(beta), force_n=leg_force, area_m2=leg_area,
                         chord_area_m2=leg_area / 4, chord_breadth_m=leg_chord_breadth, kg=leg_kg * 6),
                arch=dict(count=6, span_m=span, rise_m=rise, length_m=arch_len, thrust_n=thrust, force_n=arch_force,
                          area_m2=arch_area,
                          chord_area_m2=arch_area / 4, rib_breadth_m=arch_rib, kg=arch_kg * 6),
                transfer_ring=dict(z_m=Z_T, radius_m=r_top, depth_m=280.0, force_n=ring_force, area_m2=ring_area,
                                   kg=ring_kg),
                feet=dict(vertical_load_n=foot_v, footing_area_m2=foot_v / BEARING_PA, outward_push_n=outward),
                tie=dict(force_n=tie_force, area_m2=tie_area, kg=tie_kg))


def crown_terminal(ground_above_sea_m: float) -> dict:
    """The long-haul terminal: eight storeys of 5 m round the frame's top, a docking arm per berth group; the arms
    as first drawn, braced from the frame by two struts each."""
    floors_kg = TERMINAL_FLOOR_M2 * (FLOOR_KG_M2 + 100.0)
    arms_kg = sum(L * (GANGWAY_KG_M + 12e3) for L in ARMS)       # a lattice arm of about 12 t/m and its gangway
    return dict(z_m=35000.0 - ground_above_sea_m + 30.0, floor_m2=TERMINAL_FLOOR_M2, storeys=8, storey_m=5.0,
                mass_kg=floors_kg + arms_kg, arms=len(ARMS), arm_lengths_m=ARMS,
                cd_area=1.4 * 40.0 * 420.0 + sum(1.2 * 0.35 * 60.0 * L for L in ARMS))


def dynamics(fr: dict, v_service: float, rho_air) -> dict:
    """Sway under the service gust, first period and whole-frame buckling, from the members' bending stiffness
    (48 members on a circle: I = 48 A r^2 / 2, their axes leaning by the helix angle)."""
    rows = fr['rows']
    z = np.array([r['z_m'] for r in rows])
    dz = float(np.mean(np.diff(z)))
    area = np.array([r['member_area_m2'] for r in rows])
    rad = np.array([r['radius_m'] for r in rows])
    cos3 = np.cos(np.radians([r['lean_deg'] for r in rows])) ** 3
    EI = E * 2 * N * area * rad ** 2 / 2 * cos3
    mass = np.gradient(-np.array([r['mass_above_kg'] for r in rows]), z)
    mass = np.maximum(mass, 1.0)
    ratio = (v_service / DESIGN_GUST) ** 2
    load = np.gradient(np.array([r['shear_n'] for r in rows]), -z) * ratio
    moment = np.array([r['moment_nm'] for r in rows]) * ratio
    axial = G * np.array([r['mass_above_kg'] for r in rows])
    sway, omega2, buckling = lt.response(dz, EI[:, None], mass[:, None], load[:, None], moment[:, None], axial[:, None])
    return dict(sway_at_service_m=float(sway[0]), period_s=float(2 * math.pi / math.sqrt(omega2[0])),
                buckling_factor=float(buckling[0]))


DESIGN_GUST = 27.3


def results() -> dict:
    global DESIGN_GUST
    study = json.loads(STUDY.read_text())
    gcm = json.loads(GCM_WINDS.read_text())
    ground = study['site']['ground_above_sea_m']
    rho_air = air_density(gcm, ground)
    DESIGN_GUST = study['design_winds']['model_based_design_gust_m_s']
    v_service = study['design_winds']['service_gust_m_s']
    ring_list = rings()
    disk_list = [disk(level) for level in DISKS]
    for r in ring_list:
        depth = sum(r['storey_heights_m']) + 2.0 + r['truss']['root_depth_m'] / 2
        r['cd_area'] = 1.4 * depth * 2 * (r['outer_radius_m'] + 10.0)
    for d in disk_list:
        d['cd_area'] = 1.4 * (12.0 + DISK_COVER * 15.0) * 2 * d['radius_m'] + 1.2 * 0.1 * d['lens_depth_m'] * d['radius_m']
    crown = crown_terminal(ground)
    core = dict(z0_m=disk_list[0]['z_m'] - disk_list[0]['lens_depth_m'], z1_m=crown['z_m'] - 30.0,
                diameter_m=CORE_DIAMETER, kg_m=CORE_KG_M)
    core['kg'] = CORE_KG_M * (core['z1_m'] - core['z0_m'])
    fr = frame(ring_list + disk_list, core, crown, rho_air, DESIGN_GUST)
    bs = base(fr, rho_air, DESIGN_GUST)
    dyn = dynamics(fr, v_service, rho_air)
    steel = dict(members=fr['member_kg'], node_and_intermediate_rings=fr['ring_kg'],
                 ring_frame_bands=sum(r['frame_band_kg'] for r in ring_list),
                 legs=bs['leg']['kg'], arches=bs['arch']['kg'], transfer_ring=bs['transfer_ring']['kg'],
                 buried_tie=bs['tie']['kg'],
                 ring_structures=sum(r['structure_kg'] for r in ring_list),
                 disk_structures=sum(d['structure_kg'] for d in disk_list), core=core['kg'])
    floors_mass = (sum(r['mass_kg'] - r['structure_kg'] for r in ring_list)
                   + sum(d['mass_kg'] - d['structure_kg'] for d in disk_list) + crown['mass_kg'])
    floor_area = (sum(r['floor_m2'] for r in ring_list) + sum(d['floor_m2'] for d in disk_list)
                  + crown['floor_m2'])
    rows = fr['rows']
    sample = [min(rows, key=lambda r: abs(r['z_m'] - z)) for z in (Z_T + 50, 6000, 9000, 12000, 15000, 18000,
                                                                       21000, 23000)]
    return dict(
        schema=SCHEMA,
        evidence=('A proposed design sized by first-order hand calculation (statics, allowable stress, slenderness '
                  'limits, Rayleigh dynamics); not checked with a structural analysis model.'),
        reading_rule=('Masses in kg, forces in N, lengths in m; z is height above the summit ground at the tower; '
                      'areas are steel cross-sections; chord areas are one of four chords.'),
        basis=dict(gravity_m_s2=G, steel_allowable_pa=SIGMA, cable_allowable_pa=SIGMA_CABLE,
                   design_gust_m_s=DESIGN_GUST, service_gust_m_s=v_service, dynamic_factor=DYNAMIC,
                   ground_above_sea_m=ground, node_spacing_m=NODE_DZ, transfer_ring_m=Z_T,
                   masses_kg_m2=dict(floor=FLOOR_KG_M2, roof_slab=ROOF_SLAB_KG_M2, park=PARK_KG_M2,
                                     sealed_halls=HALLS_KG_M2)),
        frame=dict(members_each_way=N, member_breadth_rule=f'max({MIN_MEMBER} m, {MEMBER_SHARE} x spacing)',
                   slenderness_limit_breadths=SLENDER,
                   samples=[dict(z_m=r['z_m'], radius_m=r['radius_m'], breadth_m=r['breadth_m'],
                                 lean_deg=r['lean_deg'], member_force_n=r['member_force_n'],
                                 member_area_m2=r['member_area_m2'], chord_area_m2=r['member_area_m2'] / 4)
                            for r in sample],
                   profile=[dict(z_m=r['z_m'], area_m2=r['member_area_m2'], breadth_m=r['breadth_m'],
                                 member_force_n=r['member_force_n'])
                            for r in rows[::50]],
                   node_rings=fr['node_rings'], intermediate_rings=fr['intermediate_rings']),
        base=bs, rings=ring_list, disks=disk_list, core=core, terminal=crown, dynamics=dyn,
        steel_kg=steel, steel_total_kg=sum(steel.values()), floors_and_parks_kg=floors_mass,
        floor_area_m2=floor_area,
        study_square_lattice=dict(frame_mass_kg=study['port']['design']['frame_mass_mt'] * 1e9,
                                  floor_mass_kg=study['port']['design']['floor_mass_mt'] * 1e9,
                                  period_s=study['port']['design']['period_s'],
                                  sway_at_service_m=study['port']['design']['sway_at_service_m']))


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main() -> int:
    product = results()
    product['producer'] = dict(study='summit_tower', files={'research/studies/summit_tower/form.py': digest(__file__)},
                               inputs={'research/studies/summit_tower/results/summit_tower.json': digest(STUDY),
                                       'climate/results/gcm/site_winds_A28_dim5_summit.json': digest(GCM_WINDS)})
    out = HERE / 'results' / 'summit_tower_form.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    s = product['steel_kg']
    print(f"steel {product['steel_total_kg'] / 1e9:.1f} Mt: " + ', '.join(f'{k} {v / 1e9:.2f}' for k, v in s.items()))
    print(f"floors and parks {product['floors_and_parks_kg'] / 1e9:.1f} Mt; floor {product['floor_area_m2'] / 1e6:.1f} km2")
    print(f"dynamics {product['dynamics']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
