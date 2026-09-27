"""An architectural model of the summit port (hexagonal frame) in Blender, from the summit tower study's results.

    blender -b --factory-startup --python visualization/summit-port/model.py -- visualization/summit-port/build/renders [VIEW ...]

Frame height, width profile, band spacing, storey area and zones come from the tower model. The architecture on it is
a design proposal: six leg towers (lattice legs with lift cores), a four-storey wing every 200 m wrapped round two
opposite legs and stepping one leg round per band, berth fingers every 60 m on each wing's outer rim, roof and inner
gardens, gate interchanges at the leg feet, and a sealed terminal drum with docking arms at the crown. Everything is
built through the data API (no operators, so no undo copies), with ships, people, trees and berths instanced.
"""
import json
import math
import sys
from pathlib import Path

import bmesh
import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Matrix, Vector

ARGS = sys.argv[sys.argv.index('--') + 1:]
OUT = Path(ARGS[0])
VIEWS = set(ARGS[1:])
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(__file__).resolve().parents[2]
TOWER = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'
R_ = json.loads(TOWER.read_text())
D = R_['port']['design']
H = D['height_km'] * 1e3
B0, BT, P = D['base_width_m'], 320.0, D['flare_exponent']
FLIGHT_LO = 35000.0 - R_['site']['ground_above_sea_m']
BAND, STOREY_H, STOREYS, DEPTH, CAP, SHARE = 200.0, 4.0, 4, 40.0, 100_000.0, 0.02
WING_H = STOREY_H * STOREYS
BERTH = 60.0

# Where the model is drawn in full (slabs, columns, berths, ships, people); elsewhere wings are solid blocks.
DETAIL = [(0, 0.0, 1400.0), (0, 5500.0, 8000.0), (3, 5500.0, 8000.0), (1, 5500.0, 8000.0), (5, 5500.0, 8000.0)]
PEOPLE_ZONES = [(0, 6500.0, 6700.0), (0, 0.0, 700.0)]


def width(z):
    return BT + (B0 - BT) * max(0.0, 1.0 - z / H) ** P


def radius(z):
    return width(z) / 2


def storey(z):
    return min(CAP, SHARE * width(z) ** 2)


def leg_width(z):
    return 0.018 * width(z) + 24.0


def corner(j, z):
    a = math.radians(60 * j)
    return Vector((radius(z) * math.cos(a), radius(z) * math.sin(a), z))


def detailed(j, z):
    return any(j == jj and lo <= z <= hi for jj, lo, hi in DETAIL) or z > 22500


def srgb(h, a=1.0):
    h = h.lstrip('#')
    c = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    return (*c, a)


C = dict(chord='#8b9299', brace='#a4aaaf', core='#ece8df', slab='#f4f2ec', column='#d9d5cc', glass='#b7d3e2',
         roof='#86b36b', tree='#5b8f45', dock='#d7d2c8', hull='#fbfbf9', taxi='#e8b53b', regional='#d9774b',
         longhaul='#4f7ea6', person='#2f3337', ground='#cfc8b8', plaza='#e8e2d6', rail='#b7b0a4', lobby='#f0d9a8',
         water='#7fb8d6', sealed='#a9cfe0', open_slab='#efe9dc')

# ------------------------------------------------------------------------------------------------------ scene
bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.context.preferences.edit.undo_steps = 0
scene = bpy.context.scene
scene.render.engine = 'BLENDER_WORKBENCH'
sh = scene.display.shading
sh.light = 'STUDIO'
sh.studio_light = 'outdoor.sl'
sh.use_world_space_lighting = True
sh.studiolight_rotate_z = math.radians(35)
sh.color_type = 'OBJECT'
sh.show_cavity = True
sh.cavity_type = 'BOTH'
sh.cavity_ridge_factor = 1.0
sh.cavity_valley_factor = 1.0
sh.show_object_outline = True
sh.object_outline_color = (0.18, 0.2, 0.22)
sh.show_shadows = False
sh.shadow_intensity = 0.45
scene.display.light_direction = (0.45, -0.35, 0.82)
scene.display.shadow_shift = 0.05
scene.display.shadow_focus = 0.0
scene.display.render_aa = '16'
scene.view_settings.view_transform = 'Standard'
world = bpy.data.worlds.new('world')
world.color = (0.86, 0.89, 0.92)
scene.world = world


def obj(name, me, col):
    o = bpy.data.objects.new(name, me)
    o.color = srgb(C[col])
    scene.collection.objects.link(o)
    return o


def mesh_from(bm, name):
    me = bpy.data.meshes.new(name)
    bm.normal_update()
    bm.to_mesh(me)
    bm.free()
    return me


def box(bm, centre, size, rot=0.0, tilt=None):
    """A box of size (length along heading, width, height) centred at `centre`, turned `rot` about z."""
    lx, ly, lz = (s / 2 for s in size)
    m = Matrix.Rotation(rot, 3, 'Z')
    vs = []
    for sx, sy, sz in ((-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)):
        vs.append(bm.verts.new(Vector(centre) + m @ Vector((sx * lx, sy * ly, sz * lz))))
    for f in ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
        bm.faces.new([vs[i] for i in f])


def quad_prism(bm, quad, z0, z1):
    """A vertical prism over a plan quadrilateral (four xy points)."""
    lo = [bm.verts.new((x, y, z0)) for x, y in quad]
    hi = [bm.verts.new((x, y, z1)) for x, y in quad]
    bm.faces.new(lo[::-1])
    bm.faces.new(hi)
    for i in range(4):
        bm.faces.new((lo[i], lo[(i + 1) % 4], hi[(i + 1) % 4], hi[i]))


def tubes(name, members, col, res=1):
    cu = bpy.data.curves.new(name, 'CURVE')
    cu.dimensions = '3D'
    cu.bevel_depth = 1.0
    cu.bevel_resolution = res
    cu.use_fill_caps = True
    for pts, radii in members:
        sp = cu.splines.new('POLY')
        sp.points.add(len(pts) - 1)
        for i, (p, r) in enumerate(zip(pts, radii)):
            sp.points[i].co = (p[0], p[1], p[2], 1.0)
            sp.points[i].radius = r
    return obj(name, cu, col)


def instance(name, me, col, placements):
    """Linked copies of one mesh: placements are (location, rotation about z, scale)."""
    for i, (loc, rot, scale) in enumerate(placements):
        o = bpy.data.objects.new(f'{name}{i}', me)
        o.location = loc
        o.rotation_euler = (0, 0, rot)
        o.scale = (scale, scale, scale)
        o.color = srgb(C[col])
        scene.collection.objects.link(o)


# ------------------------------------------------------------------------------------------------- the frame
ZZ = [H * (i / 240) ** 1.0 for i in range(241)]


def build_frame():
    chords, braces = [], []
    for j in range(6):
        a = math.radians(60 * j)
        radial, tang = Vector((math.cos(a), math.sin(a), 0)), Vector((-math.sin(a), math.cos(a), 0))
        # four chords of the leg column
        for sr, st in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            pts = [corner(j, z) + (radial * sr + tang * st) * leg_width(z) / 2 for z in ZZ]
            chords.append((pts, [0.05 * leg_width(z) for z in ZZ]))
        # the leg's own bracing: an X in each face, one panel per leg width
        z = 0.0
        while z < H - 1:
            w = leg_width(z)
            z2 = min(H, z + w)
            quad0 = [corner(j, z) + (radial * sr + tang * st) * w / 2 for sr, st in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
            quad1 = [corner(j, z2) + (radial * sr + tang * st) * leg_width(z2) / 2 for sr, st in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
            rb = 0.022 * w
            for i in range(4):
                braces.append(([quad0[i], quad1[(i + 1) % 4]], [rb, rb]))
                braces.append(([quad0[(i + 1) % 4], quad1[i]], [rb, rb]))
                braces.append(([quad1[i], quad1[(i + 1) % 4]], [rb * 0.8, rb * 0.8]))
            z = z2
    # the frame's faces: an X between neighbouring legs, one panel per 0.9 of the face's length, and ring girders
    faces = []
    z = 0.0
    while z < H - 1:
        z2 = min(H, z + radius(z) * 0.9)
        for j in range(6):
            p0, p1 = corner(j, z), corner(j + 1, z)
            q0, q1 = corner(j, z2), corner(j + 1, z2)
            r0, r1 = 0.28 * 0.05 * leg_width(z) * 4, 0.28 * 0.05 * leg_width(z2) * 4
            faces.append(([p0, q1], [r0, r1]))
            faces.append(([p1, q0], [r0, r1]))
            faces.append(([q0, q1], [r1 * 1.2, r1 * 1.2]))
        z = z2
    tubes('leg_chords', chords, 'chord', res=1)
    tubes('leg_braces', braces, 'brace', res=0)
    tubes('face_braces', faces, 'chord', res=1)
    # inclined lift-and-stair cores inside each leg: a 36 m square tube swept up the leg's centre line
    bm = bmesh.new()
    for j in range(6):
        a = math.radians(60 * j)
        radial, tang = Vector((math.cos(a), math.sin(a), 0)), Vector((-math.sin(a), math.cos(a), 0))
        rings = []
        for z in ZZ[:-2]:
            c = corner(j, z)
            rings.append([bm.verts.new(c + (radial * sr + tang * st) * 18.0)
                          for sr, st in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
        for r0, r1 in zip(rings[:-1], rings[1:]):
            for i in range(4):
                bm.faces.new((r0[i], r0[(i + 1) % 4], r1[(i + 1) % 4], r1[i]))
        bm.faces.new(rings[-1])
    obj('cores', mesh_from(bm, 'cores'), 'core')


# -------------------------------------------------------------------------------------------------- the wings
def wing_arms(j, z):
    """The two arms of a wing on leg j at height z: each runs from the leg's corner along one face, straddling it."""
    area = storey(z) / 2                         # two wings share a band's storey
    arm = area / DEPTH / 2
    c = corner(j, z).to_2d()
    arms = []
    for nb in (j - 1, j + 1):
        d = (corner(nb, z).to_2d() - c).normalized()
        n = Vector((d.y, -d.x))
        if n.dot(c) < 0:                          # outward normal of this face
            n = -n
        arms.append((c, d, n, arm))
    return arms


def arm_quad(c, d, n, length, offset=0.0):
    p0, p1 = c + d * offset, c + d * (offset + length)
    return [p0 - n * DEPTH / 2, p1 - n * DEPTH / 2, p1 + n * DEPTH / 2, p0 + n * DEPTH / 2]


def zone(z):
    return 'open' if z < 3000 else ('sealed' if z >= 15000 else 'enclosed')


anchors_world = {}


def build_wings():
    solid, slabs, columns, glass, roof, dock_bm, lobby = (bmesh.new() for _ in range(7))
    berths, ships_taxi, ships_regional, trees, people = [], [], [], [], []
    for k in range(1, int(H // BAND)):
        z = k * BAND
        for j in (k % 3, k % 3 + 3):
            arms = wing_arms(j, z)
            zn = zone(z)
            detail = detailed(j, z)
            # transfer lobby round the core where the wing meets the leg
            cc = corner(j, z)
            box(lobby, (cc.x, cc.y, z + WING_H / 2 + 2), (leg_width(z) * 0.8, leg_width(z) * 0.8, WING_H + 8),
                rot=math.radians(60 * j))
            for c, d, n, length in arms:
                q = arm_quad(c, d, n, length, offset=leg_width(z) * 0.3)
                if not detail:
                    quad_prism(solid if zn != 'sealed' else glass, [tuple(v) for v in q], z, z + WING_H)
                    quad_prism(roof, [tuple(v) for v in arm_quad(c, d, n * 0.8, length, leg_width(z) * 0.3)],
                               z + WING_H, z + WING_H + 1.2)
                    continue
                # slabs for each storey and the roof, perimeter columns every 10 m, glazing for sealed wings
                for s in range(STOREYS + 1):
                    quad_prism(slabs, [tuple(v) for v in q], z + s * STOREY_H, z + s * STOREY_H + 0.6)
                heading = math.atan2(d.y, d.x)
                steps = int(length // 10)
                for i in range(steps + 1):
                    p = c + d * (leg_width(z) * 0.3 + i * length / steps)
                    for side in (-1, 1):
                        cp = p + n * side * (DEPTH / 2 - 1.0)
                        box(columns, (cp.x, cp.y, z + WING_H / 2), (0.8, 0.8, WING_H), rot=heading)
                if zn != 'open':
                    for s_ in range(STOREYS):
                        quad_prism(glass, [tuple(v) for v in arm_quad(c, d, n * 0.95, length, leg_width(z) * 0.3)],
                                   z + s_ * STOREY_H + 0.6, z + (s_ + 1) * STOREY_H)
                # roof garden in two beds with a walk between (open and enclosed zones)
                for side in ((-1, 1) if zn != 'sealed' else ()):
                    bed = [c + d * (leg_width(z) * 0.3) + n * side * 4, c + d * (leg_width(z) * 0.3 + length) + n * side * 4,
                           c + d * (leg_width(z) * 0.3 + length) + n * side * 17, c + d * (leg_width(z) * 0.3) + n * side * 17]
                    quad_prism(roof, [tuple(v) for v in bed], z + WING_H + 0.6, z + WING_H + 1.4)
                # berth fingers every 60 m on the outer rim, at the top storey, with ships at some
                nb = int(length // BERTH)
                for b in range(nb):
                    p = c + d * (leg_width(z) * 0.3 + (b + 0.5) * BERTH) + n * (DEPTH / 2 + 22)
                    box(dock_bm, (p.x, p.y, z + 3 * STOREY_H - 0.6), (8.0, 44.0, 1.2), rot=heading)
                    heading_ship = heading
                    if b % 2 == 0:
                        big = zn == 'enclosed' and z >= 10000
                        sp = p + n * (30 if not big else 55) + d * (0 if b % 4 == 0 else 6)
                        (ships_regional if (big or b % 6 == 4) else ships_taxi).append(
                            ((sp.x, sp.y, z + 3 * STOREY_H + 6), heading_ship, 1.0))
                # trees along the beds and people on the roof walk, where the camera comes close
                for i in range(0, int(length) if zn != 'sealed' else 0, 9):
                    for side in (-1, 1):
                        tp = c + d * (leg_width(z) * 0.3 + i + 4) + n * side * 10.5
                        trees.append(((tp.x, tp.y, z + WING_H + 1.4), 0.0, 0.8 + 0.4 * ((i * 7 + side) % 5) / 5))
                if any(j == jj and lo <= z <= hi for jj, lo, hi in PEOPLE_ZONES):
                    for i in range(0, int(length), 5):
                        for lane in (-1.5, 1.0):
                            if (i * 13 + int(lane * 10)) % 3:
                                continue
                            pp = c + d * (leg_width(z) * 0.3 + i + 2) + n * lane
                            people.append(((pp.x, pp.y, z + WING_H + 0.6), 0.0, 1.0))
                if j == 0 and k == 33:
                    anchors_world['v3_wing'] = tuple(c + d * (leg_width(z) * 0.3 + 150))
    obj('wings_solid', mesh_from(solid, 'wings_solid'), 'slab')
    obj('wing_slabs', mesh_from(slabs, 'wing_slabs'), 'slab')
    obj('wing_columns', mesh_from(columns, 'wing_columns'), 'column')
    obj('wing_glass', mesh_from(glass, 'wing_glass'), 'glass')
    obj('roof_gardens', mesh_from(roof, 'roof_gardens'), 'roof')
    obj('berths', mesh_from(dock_bm, 'berths'), 'dock')
    obj('lobbies', mesh_from(lobby, 'lobbies'), 'lobby')
    return ships_taxi, ships_regional, trees, people


# ----------------------------------------------------------------------------------------- ships, trees, people
def ship_mesh(length, name):
    """A sky ship: a lifting hull with a gondola and tail fins, `length` metres long, pointing along +x."""
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=24, v_segments=12, radius=1.0)
    bmesh.ops.scale(bm, vec=(length / 2, length * 0.12, length * 0.1), verts=bm.verts)
    box(bm, (length * 0.05, 0, -length * 0.12), (length * 0.35, length * 0.07, length * 0.05))
    box(bm, (-length * 0.42, 0, length * 0.08), (length * 0.12, length * 0.01, length * 0.12))
    box(bm, (-length * 0.42, 0, 0), (length * 0.12, length * 0.24, length * 0.01))
    return mesh_from(bm, name)


def tree_mesh():
    bm = bmesh.new()
    bmesh.ops.create_cone(bm, cap_ends=True, segments=8, radius1=0.4, radius2=0.3, depth=3.0,
                          matrix=Matrix.Translation((0, 0, 1.5)))
    bmesh.ops.create_icosphere(bm, subdivisions=1, radius=3.2, matrix=Matrix.Translation((0, 0, 5.5)))
    return mesh_from(bm, 'tree')


def person_mesh():
    bm = bmesh.new()
    bmesh.ops.create_cone(bm, cap_ends=True, segments=8, radius1=0.24, radius2=0.2, depth=1.45,
                          matrix=Matrix.Translation((0, 0, 0.72)))
    bmesh.ops.create_icosphere(bm, subdivisions=1, radius=0.13, matrix=Matrix.Translation((0, 0, 1.6)))
    return mesh_from(bm, 'person')


# --------------------------------------------------------------------------------------- the base and the crown
def build_base():
    ground = bmesh.new()
    box(ground, (0, 0, -2), (160000, 160000, 4))
    obj('ground', mesh_from(ground, 'ground'), 'ground')
    plaza, gates, rail = bmesh.new(), bmesh.new(), bmesh.new()
    for j in range(6):
        a = math.radians(60 * j)
        c = corner(j, 0)
        # the footing, then a gate interchange wrapped round the leg foot: a ring of halls with a plaza inside
        box(plaza, (c.x, c.y, 0.5), (620, 620, 1), rot=a)
        for i in range(4):
            ang = a + math.pi / 2 * i
            off = Vector((math.cos(ang), math.sin(ang), 0)) * 210
            box(gates, (c.x + off.x, c.y + off.y, 14), (140, 330, 28), rot=ang)
        # the rail viaduct round the base, joining the six gates
        c2 = corner(j + 1, 0)
        mid, d = (c + c2) / 2, (c2 - c)
        box(rail, (mid.x, mid.y, 14), (d.length - 700, 16, 4), rot=math.atan2(d.y, d.x))
        for t in range(1, 40):
            p = c + d * (t / 40)
            box(rail, (p.x, p.y, 6), (4, 4, 12))
    obj('plazas', mesh_from(plaza, 'plazas'), 'plaza')
    obj('gates', mesh_from(gates, 'gates'), 'core')
    obj('rail', mesh_from(rail, 'rail'), 'rail')


def build_crown():
    top_z = FLIGHT_LO - 120
    drum, arms, glass = bmesh.new(), bmesh.new(), bmesh.new()
    # the sealed long-haul terminal: a hexagonal drum round the frame's top, six storeys
    r = radius(top_z) + 40
    hexq = [(r * math.cos(math.radians(60 * i)), r * math.sin(math.radians(60 * i))) for i in range(6)]
    for s in range(7):
        zz = top_z + s * 5
        lo = [drum.verts.new((x, y, zz)) for x, y in hexq]
        hi = [drum.verts.new((x, y, zz + 0.8)) for x, y in hexq]
        drum.faces.new(lo[::-1])
        drum.faces.new(hi)
        for i in range(6):
            drum.faces.new((lo[i], lo[(i + 1) % 6], hi[(i + 1) % 6], hi[i]))
    gl = [glass.verts.new((x * 0.97, y * 0.97, top_z)) for x, y in hexq]
    gh = [glass.verts.new((x * 0.97, y * 0.97, top_z + 30)) for x, y in hexq]
    for i in range(6):
        glass.faces.new((gl[i], gl[(i + 1) % 6], gh[(i + 1) % 6], gh[i]))
    longhaul = []
    for jn in range(8):
        z = FLIGHT_LO + 60 + jn * 65
        a = math.radians(jn * 137.50776 + 20)
        length = 900 + 180 * (jn % 3)
        d = Vector((math.cos(a), math.sin(a), 0))
        nrm = Vector((-d.y, d.x, 0))
        r0 = radius(z)
        mid = d * (r0 + length / 2)
        box(arms, (mid.x, mid.y, z), (length, 36, 14), rot=a)
        box(glass, (mid.x, mid.y, z + 9), (length, 30, 4), rot=a)
        side = nrm if nrm.x >= 0 else -nrm
        for b in range(3):
            t = 0.3 + 0.3 * b
            p = d * (r0 + t * length) + side * 90
            longhaul.append(((p.x, p.y, z + 4), a, 1.0))
            g = d * (r0 + t * length) + side * 40
            box(arms, (g.x, g.y, z + 2), (8, 44, 3), rot=a + math.pi / 2)
        anchors_world[f'arm{jn}'] = tuple(d * (r0 + length))[:2] + (z,)
    obj('terminal_drum', mesh_from(drum, 'drum'), 'slab')
    obj('crown_arms', mesh_from(arms, 'arms'), 'dock')
    obj('crown_glass', mesh_from(glass, 'crown_glass'), 'glass')
    return longhaul


# ------------------------------------------------------------------------------------------------------ build
build_frame()
taxi, regional, trees, people = build_wings()
build_base()
longhaul = build_crown()
import random
random.seed(7)
for i in range(70):
    zf = random.uniform(5900.0, 7600.0)
    j = random.choice((0, 0, 0, 1, 5))
    c = corner(j, zf)
    ang = math.radians(60 * j) + random.uniform(-0.9, 0.9)
    dist = random.uniform(80.0, 900.0) * random.choice((1, -1))
    p = c + Vector((math.cos(ang), math.sin(ang), 0)) * dist
    (taxi if i % 5 else regional).append(((p.x, p.y, zf), random.uniform(0, 2 * math.pi), 1.0))
instance('taxi', ship_mesh(34.0, 'taxi'), 'taxi', taxi)
instance('regional', ship_mesh(95.0, 'regional'), 'regional', regional)
instance('longhaul', ship_mesh(230.0, 'longhaul'), 'longhaul', longhaul)
instance('tree', tree_mesh(), 'tree', trees)
instance('person', person_mesh(), 'person', people)
# the Burj Khalifa beside the base for scale
burj = bmesh.new()
bmesh.ops.create_cone(burj, cap_ends=True, segments=12, radius1=45, radius2=3, depth=828,
                      matrix=Matrix.Translation((radius(0) + 60, -520, 414)))
obj('burj', mesh_from(burj, 'burj'), 'regional')
print('BUILT', len(bpy.data.objects), 'objects')

cam_data = bpy.data.cameras.new('cam')
cam_data.clip_start = 1.0
cam_data.clip_end = 300000.0
cam = bpy.data.objects.new('cam', cam_data)
scene.collection.objects.link(cam)
scene.camera = cam
labels = json.loads((OUT / 'labels.json').read_text()) if (OUT / 'labels.json').is_file() else {}


def aim(loc, target, lens, clip_start=1.0, clip_end=60000.0):
    cam.location = Vector(loc)
    cam.rotation_euler = (Vector(target) - Vector(loc)).to_track_quat('-Z', 'Y').to_euler()
    cam_data.lens = lens
    cam_data.clip_start = clip_start
    cam_data.clip_end = clip_end


def render(name, rx, ry, anchors):
    if VIEWS and name not in VIEWS:
        return
    scene.render.resolution_x, scene.render.resolution_y = rx, ry
    scene.render.filepath = str(OUT / f'{name}.png')
    bpy.ops.render.render(write_still=True)
    out = {}
    for k, p in anchors.items():
        v = world_to_camera_view(scene, cam, Vector(p))
        out[k] = [round(v.x * rx, 1), round((1 - v.y) * ry, 1), round(v.z, 1)]
    labels[name] = dict(size=[rx, ry], anchors=out)
    print('RENDERED', name)


c0 = corner(0, 0)
# 1: the whole tower
aim((-30000, -52000, 6000), (0, 0, 11000), 70, clip_start=50, clip_end=300000)
render('tower', 1500, 1700, {
    'z0': tuple(corner(3, 0)), 'z3': tuple(corner(3, 3000)), 'z10': tuple(corner(3, 10000)),
    'z15': tuple(corner(3, 15000)), 'z24': tuple(corner(3, H)),
    'segment': tuple(corner(0, 6800)), 'base': tuple(corner(0, 300)), 'crown': (0, 0, FLIGHT_LO + 200),
    'burj': (radius(0) + 60, -520, 828)})

# 2: a segment of leg 0 at 6-7.5 km: three of its wings, its neighbours' wings across the faces
L = corner(0, 6600)
aim((L.x + 1150, L.y - 1250, 6420), (L.x - 150, L.y + 260, 6660), 30, clip_start=2)
arm1 = (corner(1, 6600).to_2d() - corner(0, 6600).to_2d()).normalized()
arm5 = (corner(5, 6600).to_2d() - corner(0, 6600).to_2d()).normalized()
render('segment', 1600, 1100, {
    'wing_6600': tuple(corner(0, 6600).to_2d() + arm5 * 300) + (6616,),
    'wing_7200': tuple(corner(0, 7200).to_2d() + arm1 * 250) + (7216,),
    'lobby': tuple(corner(0, 6630)), 'leg': tuple(corner(0, 6950)), 'core': tuple(corner(0, 6300)),
    'berths': tuple(corner(0, 6600).to_2d() + arm1 * 420 + Vector((arm1.y, -arm1.x)) * 45) + (6611,)})

# 3: human scale: along a wing's roof and berths at 6.6 km
W = Vector(anchors_world['v3_wing'])
arm_dir = (corner(1, 6600).to_2d() - corner(0, 6600).to_2d()).normalized()
out_n = Vector((arm_dir.y, -arm_dir.x))
if out_n.dot(corner(0, 6600).to_2d()) < 0:
    out_n = -out_n
eye = W + out_n * 70 + arm_dir * 95
aim((eye.x, eye.y, 6600 + WING_H + 24), (W.x - arm_dir.x * 40, W.y - arm_dir.y * 40, 6600 + WING_H - 2), 24,
    clip_start=0.5)
render('wing', 1600, 1000, {'roof': (W.x, W.y, 6600 + WING_H + 2), 'berth': tuple(W + out_n * 42) + (6600 + 11,),
                            'core': tuple(corner(0, 6616))})

# 4: the base at leg 0: gate interchange, footing, the first wings above
aim((c0.x + 1500, c0.y - 1450, 330), (c0.x - 150, c0.y + 120, 620), 28, clip_start=2)
render('base', 1600, 1000, {'gate': (c0.x, c0.y, 30), 'wing_600': tuple(corner(0, 616)),
                            'wing_1200': tuple(corner(0, 1216)), 'rail': tuple((corner(0, 0) + corner(1, 0)) / 2)})

# 5: the crown in the flight band
aim((1900, -2600, FLIGHT_LO + 260), (0, 0, FLIGHT_LO + 120), 30, clip_start=5)
render('crown', 1600, 1000, {'drum': (radius(FLIGHT_LO) + 40, 0, FLIGHT_LO - 100),
                             **{k: v for k, v in anchors_world.items() if k.startswith('arm')}})

# 6: looking up inside the frame from the plaza of the base
aim((0, 0, 60), (900, 300, 9000), 16, clip_start=2)
render('inside', 1200, 1500, {'z3': (radius(3000) * 0.8, 0, 3000)})

(OUT / 'labels.json').write_text(json.dumps(labels, indent=1))
print('DONE')
