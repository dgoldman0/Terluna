"""The long-haul ships at the summit port's crown, to scale: a massing model in Blender (Eevee), one scene per ship
class, from the sky-fleet study's berth layouts and the summit tower study's frame.

    blender -b --factory-startup --python visualization/sky-fleet/crown.py -- visualization/sky-fleet/build/renders [NAME ...]

Read from the products: the frame's height and width profile, the flight band's lower edge, the disks' heights, each
class's hull length and diameter, fin span, berth count, and the arms, levels and pitch of its berths. Drawn as a
design proposal at plausible sizes: the round diagrid's members (two opposite spirals at 62 degrees), the sealed
disks and their core, the glazed terminal drum in the flight band, the arms as deep trusses with an enclosed
gangway, and each ship's 36-sided hull (as the great rigid airships were built), fins, passenger decks and engine
pods. Nothing here is structural design. The views are the same for every class, so the sizes compare directly: a
perspective from the south-south-east, a plan, and a side elevation of all classes beside the Hindenburg.
"""
import json
import math
import sys
from pathlib import Path

import bmesh
import bpy
from mathutils import Matrix, Vector

ARGS = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
OUT = Path(ARGS[0]) if ARGS else Path(__file__).resolve().parent / 'build' / 'renders'
ONLY = {a for a in ARGS[1:] if not a.startswith('--')}
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(__file__).resolve().parents[2]
FLEET = json.loads((ROOT / 'research' / 'studies' / 'sky_fleet' / 'results' / 'sky_fleet.json').read_text())
FORM = json.loads((ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower_form.json').read_text())
SHAPE = FORM['width_profile']        # the chosen form's frame
CROWN = FLEET['long_haul']['crown']
H = SHAPE['height_m']
B0, BT, P = SHAPE['base_width_m'], CROWN['top_width_m'], SHAPE['flare_exponent']
BAND = CROWN['band_bottom_above_summit_m']
Z0 = H - 5300.0                      # the lowest height drawn
LIFT = H - 400.0                     # the scene's origin, 400 m below the top, keeps coordinates small
HAZE_M = 30000.0                     # aerial perspective: the share of the haze colour is 1 - exp(-distance / HAZE_M)
HAZE_COLOUR = '#dde4e9'
SKY = ('#edf0f1', '#d3dde4', '#93a9bb')
SMALL = '--small' in sys.argv
HINDENBURG = dict(length_m=245.0, diameter_m=41.2, fin_span_m=15.0)


def width(z):
    return BT + (B0 - BT) * max(0.0, 1.0 - z / H) ** P


def radius(z):
    return width(z) / 2.0


def lin(h):
    h = h.lstrip('#')
    return tuple((int(h[i:i + 2], 16) / 255.0) ** 2.2 for i in (0, 2, 4)) + (1.0,)


# ------------------------------------------------------------------------------------------------- materials

def principled(name, colour, rough=0.6, metal=0.0, spec=0.5):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']
    b.inputs['Base Color'].default_value = lin(colour)
    b.inputs['Roughness'].default_value = rough
    b.inputs['Metallic'].default_value = metal
    if 'Specular IOR Level' in b.inputs:
        b.inputs['Specular IOR Level'].default_value = spec
    m.diffuse_color = lin(colour)
    return m, b


def banded(name, a, b, axis, period, share, rough=(0.25, 0.6)):
    """A facade: colour `a` (glass) broken by bands of colour `b` (slab edges or mullions) every `period` metres
    along object axis `axis`, the bands taking `share` of it."""
    m, bsdf = principled(name, a, rough[0])
    nt = m.node_tree
    tc = nt.nodes.new('ShaderNodeTexCoord')
    sep = nt.nodes.new('ShaderNodeSeparateXYZ')
    nt.links.new(tc.outputs['Object'], sep.inputs['Vector'])
    div = nt.nodes.new('ShaderNodeMath')
    div.operation = 'DIVIDE'
    div.inputs[1].default_value = period
    nt.links.new(sep.outputs[axis], div.inputs[0])
    frac = nt.nodes.new('ShaderNodeMath')
    frac.operation = 'FRACT'
    nt.links.new(div.outputs[0], frac.inputs[0])
    step = nt.nodes.new('ShaderNodeMath')
    step.operation = 'LESS_THAN'
    step.inputs[1].default_value = share
    nt.links.new(frac.outputs[0], step.inputs[0])
    mix = nt.nodes.new('ShaderNodeMix')
    mix.data_type = 'RGBA'
    mix.inputs[6].default_value = lin(a)
    mix.inputs[7].default_value = lin(b)
    nt.links.new(step.outputs[0], mix.inputs[0])
    nt.links.new(mix.outputs[2], bsdf.inputs['Base Color'])
    rmix = nt.nodes.new('ShaderNodeMath')
    rmix.operation = 'MULTIPLY_ADD'
    rmix.inputs[1].default_value = rough[1] - rough[0]
    rmix.inputs[2].default_value = rough[0]
    nt.links.new(step.outputs[0], rmix.inputs[0])
    nt.links.new(rmix.outputs[0], bsdf.inputs['Roughness'])
    return m


def hazed(m):
    """Fade a material toward the haze colour with its distance from the camera."""
    nt = m.node_tree
    out = nt.nodes['Material Output']
    surface = out.inputs['Surface'].links[0].from_socket
    cam = nt.nodes.new('ShaderNodeCameraData')
    k = nt.nodes.new('ShaderNodeMath')
    k.operation = 'MULTIPLY'
    k.inputs[1].default_value = -1.0 / HAZE_M
    nt.links.new(cam.outputs['View Distance'], k.inputs[0])
    e = nt.nodes.new('ShaderNodeMath')
    e.operation = 'EXPONENT'
    nt.links.new(k.outputs[0], e.inputs[0])
    f = nt.nodes.new('ShaderNodeMath')
    f.operation = 'SUBTRACT'
    f.inputs[0].default_value = 1.0
    nt.links.new(e.outputs[0], f.inputs[1])
    glow = nt.nodes.new('ShaderNodeEmission')
    glow.inputs['Color'].default_value = lin(HAZE_COLOUR)
    glow.inputs['Strength'].default_value = 1.0
    mix = nt.nodes.new('ShaderNodeMixShader')
    nt.links.new(f.outputs[0], mix.inputs['Fac'])
    nt.links.new(surface, mix.inputs[1])
    nt.links.new(glow.outputs['Emission'], mix.inputs[2])
    nt.links.new(mix.outputs['Shader'], out.inputs['Surface'])
    return m


def materials():
    mats = {}
    mats['frame'] = principled('frame', '#c9c6be', 0.55)[0]
    mats['core'] = principled('core', '#d9d4ca', 0.7)[0]
    mats['slab'] = principled('slab', '#ebe7df', 0.7)[0]
    mats['disk'] = banded('disk', '#40566a', '#e4e0d7', 'Z', 4.5, 0.22)
    mats['drum'] = banded('drum', '#3f5669', '#e6e2d9', 'Z', 5.0, 0.2)
    mats['arm'] = principled('arm', '#d3cfc6', 0.55)[0]
    mats['gangway'] = banded('gangway', '#4a6173', '#eeeae2', 'Y', 6.0, 0.25)
    mats['hull'] = principled('hull', '#dcdfe2', 0.36, 0.35, 0.6)[0]
    mats['fin'] = principled('fin', '#bcc3c9', 0.45, 0.2)[0]
    mats['deck'] = banded('deck', '#27333d', '#cfd4d8', 'X', 7.0, 0.18, rough=(0.15, 0.4))
    mats['pod'] = principled('pod', '#7d868d', 0.45, 0.4)[0]
    mats['dock'] = principled('dock', '#6f777d', 0.5, 0.3)[0]
    return {k: hazed(m) for k, m in mats.items()}


def link(name, data, mat, collection):
    o = bpy.data.objects.new(name, data)
    if mat is not None:
        o.data.materials.append(mat)
    collection.objects.link(o)
    return o


def mesh_from(bm, name, smooth=False, crease_deg=None):
    me = bpy.data.meshes.new(name)
    bm.normal_update()
    bm.to_mesh(me)
    bm.free()
    if smooth:
        for p in me.polygons:
            p.use_smooth = True
        if crease_deg is not None and hasattr(me, 'use_auto_smooth'):
            me.use_auto_smooth = True
            me.auto_smooth_angle = math.radians(crease_deg)
    return me


def tubes(name, members, mat, collection, res=2):
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
    return link(name, cu, mat, collection)


def cylinder(bm, r, z0, z1, n=128, centre=(0.0, 0.0)):
    lo = [bm.verts.new((centre[0] + r * math.cos(2 * math.pi * i / n), centre[1] + r * math.sin(2 * math.pi * i / n), z0))
          for i in range(n)]
    hi = [bm.verts.new((v.co.x, v.co.y, z1)) for v in lo]
    bm.faces.new(lo[::-1])
    bm.faces.new(hi)
    for i in range(n):
        bm.faces.new((lo[i], lo[(i + 1) % n], hi[(i + 1) % n], hi[i]))


def box(bm, lo, hi):
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    v = [bm.verts.new(p) for p in ((x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                                   (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1))]
    for f in ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
        bm.faces.new([v[i] for i in f])


# ------------------------------------------------------------------------------------------------------ scene

def setup_scene():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    bpy.context.preferences.edit.undo_steps = 0
    scene = bpy.context.scene
    # Eevee on the graphics card (this Blender build has no GPU path tracing and no denoiser).
    scene.render.engine = 'BLENDER_EEVEE'
    ee = scene.eevee
    ee.taa_render_samples = 96
    ee.use_gtao = True
    ee.gtao_distance = 80.0
    ee.gtao_factor = 1.0
    ee.use_soft_shadows = True
    ee.shadow_cascade_size = '4096'
    ee.use_ssr = True
    scene.render.filter_size = 1.1
    scene.view_settings.view_transform = 'AgX'
    scene.view_settings.look = 'AgX - Medium High Contrast'
    scene.view_settings.exposure = 0.0

    world = bpy.data.worlds.new('world')
    world.use_nodes = True
    nt = world.node_tree
    for n in list(nt.nodes):
        if n.type != 'OUTPUT_WORLD':
            nt.nodes.remove(n)
    out = nt.nodes['World Output']
    # A drawing's sky behind the model: pale at the foot of the frame, a soft cool grey-blue at its top. The model is
    # lit by an even, dimmer fill of the same tone and by the sun.
    tc = nt.nodes.new('ShaderNodeTexCoord')
    sep = nt.nodes.new('ShaderNodeSeparateXYZ')
    nt.links.new(tc.outputs['Window'], sep.inputs['Vector'])
    ramp = nt.nodes.new('ShaderNodeValToRGB')
    ramp.color_ramp.interpolation = 'EASE'
    els = ramp.color_ramp.elements
    els[0].position, els[0].color = 0.0, lin(SKY[0])
    els[1].position, els[1].color = 1.0, lin(SKY[2])
    mid = els.new(0.45)
    mid.color = lin(SKY[1])
    nt.links.new(sep.outputs['Y'], ramp.inputs['Fac'])
    shown = nt.nodes.new('ShaderNodeBackground')
    nt.links.new(ramp.outputs['Color'], shown.inputs['Color'])
    fill = nt.nodes.new('ShaderNodeBackground')
    fill.inputs['Color'].default_value = lin('#d4dde4')
    fill.inputs['Strength'].default_value = 0.36
    path = nt.nodes.new('ShaderNodeLightPath')
    mix = nt.nodes.new('ShaderNodeMixShader')
    nt.links.new(path.outputs['Is Camera Ray'], mix.inputs['Fac'])
    nt.links.new(fill.outputs['Background'], mix.inputs[1])
    nt.links.new(shown.outputs['Background'], mix.inputs[2])
    nt.links.new(mix.outputs['Shader'], out.inputs['Surface'])
    scene.world = world

    sun = bpy.data.lights.new('sun', 'SUN')
    sun.energy = 4.6
    sun.color = (1.0, 0.965, 0.92)
    sun.angle = math.radians(1.2)
    sun.shadow_cascade_count = 4
    sun.shadow_cascade_max_distance = 9000.0
    sun.shadow_cascade_exponent = 0.7
    sun.shadow_cascade_fade = 0.1
    sun_o = bpy.data.objects.new('sun', sun)
    # From the south-west, 36 degrees up.
    sun_o.rotation_euler = (math.radians(54), 0.0, math.radians(-40))
    scene.collection.objects.link(sun_o)
    return scene


# ----------------------------------------------------------------------------------------------- the crown

def interp(z, xs, ys):
    """Linear interpolation, held constant beyond the ends."""
    if z <= xs[0]:
        return ys[0]
    for (x0, y0), (x1, y1) in zip(zip(xs[:-1], ys[:-1]), zip(xs[1:], ys[1:])):
        if z <= x1:
            return y0 + (y1 - y0) * (z - x0) / (x1 - x0)
    return ys[-1]


def build_tower(col, M):
    """The tower's top from the summit tower study's sized form: 24 members each way on helices whose lean from
    vertical and breadth follow the form, its node rings and the rings between them, the 120 m core, the sealed
    disks as spoked wheels hanging their cables to hubs on the core, and the eight-storey terminal at the band."""
    frame = FORM['frame']
    lean_z = [x['z_m'] for x in frame['samples']]
    lean = [x['lean_deg'] for x in frame['samples']]
    breadth_z = [x['z_m'] for x in frame['profile']]
    breadth = [x['breadth_m'] for x in frame['profile']]
    n, dz = frame['members_each_way'], 20.0
    zs = [Z0 + i * dz for i in range(int((H - Z0) / dz) + 1)]
    turn = [0.0]
    for a, b in zip(zs[:-1], zs[1:]):
        zm = 0.5 * (a + b)
        turn.append(turn[-1] + dz * math.tan(math.radians(interp(zm, lean_z, lean))) / radius(zm))
    members = []
    for sense in (1.0, -1.0):
        for i in range(n):
            t0 = 2 * math.pi * (i + (0.5 if sense < 0 else 0.0)) / n
            pts = [(radius(z) * math.cos(t0 + sense * t), radius(z) * math.sin(t0 + sense * t), z - LIFT)
                   for z, t in zip(zs, turn)]
            members.append((pts, [0.5 * interp(z, breadth_z, breadth) for z in zs]))
    tubes('diagrid', members, M['frame'], col)
    rings = []
    heights = [(r['z_m'], 3.2) for r in frame['node_rings']] + [(r['z_m'], 2.2) for r in frame['intermediate_rings']]
    for z, rr in heights + [(H - 1.0, 4.0)]:
        if Z0 <= z <= H:
            m = 160
            rad = radius(z)
            rings.append(([(rad * math.cos(2 * math.pi * k / m), rad * math.sin(2 * math.pi * k / m), z - LIFT)
                           for k in range(m + 1)], [rr] * (m + 1)))
    tubes('rings', rings, M['frame'], col, res=1)
    core = FORM['core']
    bm = bmesh.new()
    cylinder(bm, core['diameter_m'] / 2, max(Z0, core['z0_m']) - LIFT, core['z1_m'] - LIFT, 72)
    link('core', mesh_from(bm, 'core', True), M['core'], col)
    for d in FORM['disks']:
        z, R = d['z_m'], d['radius_m']
        if not (Z0 < z < H):
            continue
        bm = bmesh.new()
        cylinder(bm, R, z - LIFT - 3.0, z - LIFT, 160)
        link(f"deck{z:.0f}", mesh_from(bm, 'deck', True), M['slab'], col)
        bm = bmesh.new()
        inner = R * math.sqrt(0.7)            # the sealed block of three storeys over 30% of the deck, at its rim
        m = 160
        lo_o = [bm.verts.new((R * math.cos(2 * math.pi * k / m), R * math.sin(2 * math.pi * k / m), z - LIFT)) for k in range(m)]
        hi_o = [bm.verts.new((v.co.x, v.co.y, z - LIFT + 15.0)) for v in lo_o]
        lo_i = [bm.verts.new((inner * math.cos(2 * math.pi * k / m), inner * math.sin(2 * math.pi * k / m), z - LIFT))
                for k in range(m)]
        hi_i = [bm.verts.new((v.co.x, v.co.y, z - LIFT + 15.0)) for v in lo_i]
        for k in range(m):
            j = (k + 1) % m
            bm.faces.new((lo_o[k], lo_o[j], hi_o[j], hi_o[k]))
            bm.faces.new((hi_i[k], hi_i[j], lo_i[j], lo_i[k]))
            bm.faces.new((hi_o[k], hi_o[j], hi_i[j], hi_i[k]))
        link(f"block{z:.0f}", mesh_from(bm, 'block', True), M['disk'], col)
        cables = []
        hub = d['hub_radius_m']
        for k in range(d['spokes']):
            a = 2 * math.pi * (k + 0.5) / d['spokes']
            cables.append(([(R * math.cos(a), R * math.sin(a), z - LIFT - 3.0),
                            (hub * math.cos(a), hub * math.sin(a), z - LIFT - d['lens_depth_m'])], [0.6, 0.6]))
        tubes(f"cables{z:.0f}", cables, M['arm'], col, res=0)
        bm = bmesh.new()
        cylinder(bm, hub + 4.0, z - LIFT - d['lens_depth_m'] - 6.0, z - LIFT - d['lens_depth_m'] + 6.0, 48)
        link(f"hub{z:.0f}", mesh_from(bm, 'hub', True), M['dock'], col)
    t = FORM['terminal']
    r_t = math.sqrt(t['floor_m2'] / t['storeys'] / math.pi)
    z0, z1 = t['z_m'], t['z_m'] + t['storeys'] * t['storey_m']
    bm = bmesh.new()
    cylinder(bm, r_t, z0 - LIFT, z1 - LIFT)
    link('terminal', mesh_from(bm, 'terminal', True), M['drum'], col)
    bm = bmesh.new()
    cylinder(bm, r_t + 3.0, z1 - LIFT, z1 - LIFT + 2.5)
    cylinder(bm, r_t + 3.0, z0 - LIFT - 2.5, z0 - LIFT)
    link('terminal_slabs', mesh_from(bm, 'terminal_slabs', True), M['slab'], col)
    return r_t


def truss(col, M, name, y0, y1, zc, depth, width, x0=0.0):
    """A deep truss arm along y from y0 to y1 at height zc: four chords, Warren diagonals on its faces, cross
    members, and an enclosed, glazed gangway through its middle."""
    sgn = 1.0 if y1 > y0 else -1.0
    length = abs(y1 - y0)
    panels = max(2, int(round(length / depth)))
    chords, webs = [], []
    rc, rw = max(2.2, depth / 22.0), max(1.3, depth / 45.0)
    for x in (x0 - width / 2, x0 + width / 2):
        for z in (zc - depth / 2, zc + depth / 2):
            chords.append(([(x, y0, z), (x, y1, z)], [rc, rc]))
        for k in range(panels):
            ya, yb = y0 + sgn * length * k / panels, y0 + sgn * length * (k + 1) / panels
            lo, hi = zc - depth / 2, zc + depth / 2
            webs.append(([(x, ya, lo if k % 2 == 0 else hi), (x, yb, hi if k % 2 == 0 else lo)], [rw, rw]))
    for k in range(panels + 1):
        yk = y0 + sgn * length * k / panels
        for z in (zc - depth / 2, zc + depth / 2):
            webs.append(([(x0 - width / 2, yk, z), (x0 + width / 2, yk, z)], [rw * 0.8, rw * 0.8]))
        webs.append(([(x0 - width / 2, yk, zc - depth / 2), (x0 - width / 2, yk, zc + depth / 2)], [rw * 0.7] * 2))
        webs.append(([(x0 + width / 2, yk, zc - depth / 2), (x0 + width / 2, yk, zc + depth / 2)], [rw * 0.7] * 2))
    tubes(name + '_chords', chords, M['arm'], col)
    tubes(name + '_webs', webs, M['arm'], col, res=1)
    bm = bmesh.new()
    g = 0.5 * width
    box(bm, (x0 - g / 2, min(y0, y1), zc - 6.0), (x0 + g / 2, max(y0, y1), zc + 6.0))
    link(name + '_gangway', mesh_from(bm, name + '_gangway'), M['gangway'], col)


# ------------------------------------------------------------------------------------------------- the ships

PEAK = 0.45 / 1.25
NORM = PEAK ** 0.45 * (1 - PEAK) ** 0.8


def profile(s):
    """The hull's radius over its greatest radius at s along its length: s^0.45 (1 - s)^0.8, normalised; it holds
    the flight models' volume, 0.475 L D^2."""
    return (max(s, 0.0) ** 0.45) * (max(1.0 - s, 0.0) ** 0.8) / NORM


SIDES = 36


def hull_mesh(length, diameter):
    """A hull of 36 flat panels round (the great rigid airships' polygon), smooth along its length."""
    R = diameter / 2
    bm = bmesh.new()
    ss = [0.5 * (1 - math.cos(math.pi * i / 110)) for i in range(111)]
    rings = []
    for s in ss[1:-1]:
        r = R * profile(s) / math.cos(math.pi / SIDES)
        rings.append([bm.verts.new((s * length, r * math.cos(2 * math.pi * (j + 0.5) / SIDES),
                                    r * math.sin(2 * math.pi * (j + 0.5) / SIDES))) for j in range(SIDES)])
    nose, tail = bm.verts.new((0.0, 0.0, 0.0)), bm.verts.new((length, 0.0, 0.0))
    for j in range(SIDES):
        bm.faces.new((nose, rings[0][(j + 1) % SIDES], rings[0][j]))
        bm.faces.new((tail, rings[-1][j], rings[-1][(j + 1) % SIDES]))
    for a, b in zip(rings[:-1], rings[1:]):
        for j in range(SIDES):
            bm.faces.new((a[j], a[(j + 1) % SIDES], b[(j + 1) % SIDES], b[j]))
    return mesh_from(bm, 'hull', True, crease_deg=6.0)


def fin_mesh(length, diameter, fin_span, fin_area):
    """Four fins in a cross, each of the model's area: root chord 4/3 and tip chord 2/3 of the mean, the root on the
    hull, the outer edge parallel to the axis a span out from the root's mean radius, the trailing edge upright."""
    R = diameter / 2
    mean = fin_area / fin_span
    root, tip = 4.0 * mean / 3.0, 2.0 * mean / 3.0
    x_te = 0.965 * length
    x_le = x_te - root
    r_mid = R * profile((x_le + x_te) / 2 / length)
    r_tip = r_mid + fin_span
    t = max(0.8, 0.03 * root)
    bm = bmesh.new()
    n = 12
    for k in range(4):
        rot = Matrix.Rotation(math.pi / 2 * k, 3, 'X')
        root_pts = [(x_le + root * i / n, R * profile((x_le + root * i / n) / length) * 0.97) for i in range(n + 1)]
        outline = root_pts + [(x_te, r_tip), (x_te - tip, r_tip)]
        lo = [bm.verts.new(rot @ Vector((x, -t / 2, r))) for x, r in outline]
        hi = [bm.verts.new(rot @ Vector((x, t / 2, r))) for x, r in outline]
        bm.faces.new(lo[::-1])
        bm.faces.new(hi)
        for i in range(len(outline)):
            j = (i + 1) % len(outline)
            bm.faces.new((lo[i], lo[j], hi[j], hi[i]))
    return mesh_from(bm, 'fins')


def deck_mesh(length, diameter):
    """The passenger decks: a band of windows on each lower flank, from 0.2 to 0.62 of the length."""
    R = diameter / 2
    bm = bmesh.new()
    ss = [0.20 + 0.42 * i / 80 for i in range(81)]
    angles = [math.radians(-9.0 - 22.0 * i / 6) for i in range(7)]
    for side in (1.0, -1.0):
        rows = []
        for a in angles:
            rows.append([bm.verts.new((s * length, side * (R * profile(s) + 0.4) * math.cos(a),
                                       (R * profile(s) + 0.4) * math.sin(a))) for s in ss])
        for ra, rb in zip(rows[:-1], rows[1:]):
            for i in range(len(ss) - 1):
                f = (ra[i], ra[i + 1], rb[i + 1], rb[i]) if side > 0 else (ra[i], rb[i], rb[i + 1], ra[i + 1])
                bm.faces.new(f)
    return mesh_from(bm, 'decks', True)


def pod_mesh(length, diameter):
    """Eight engine pods on short pylons along the lower flanks, each with a spinner and a propeller disc."""
    R = diameter / 2
    bm = bmesh.new()
    pl, pr = 0.032 * length, 0.0070 * length
    for s in (0.30, 0.44, 0.58, 0.72):
        for side in (1.0, -1.0):
            a = math.radians(-50.0)
            r = R * profile(s)
            off = r + 2.2 * pr
            c = Vector((s * length, side * off * math.cos(a), off * math.sin(a)))
            rings = []
            for u, scale in ((-0.5, 0.55), (-0.35, 1.0), (0.3, 1.0), (0.5, 0.35)):
                rings.append([bm.verts.new(c + Vector((u * pl, scale * pr * math.cos(t), scale * pr * math.sin(t))))
                              for t in [2 * math.pi * i / 20 for i in range(20)]])
            bm.faces.new(rings[0][::-1])
            bm.faces.new(rings[-1])
            for ra, rb in zip(rings[:-1], rings[1:]):
                for i in range(20):
                    bm.faces.new((ra[i], ra[(i + 1) % 20], rb[(i + 1) % 20], rb[i]))
            # a four-bladed propeller behind the pod
            hub = c + Vector((0.56 * pl, 0.0, 0.0))
            blade, chord, thick = 2.4 * pr, 0.45 * pr, 0.12 * pr
            for q in range(4):
                t = math.pi / 4 + math.pi / 2 * q
                d = Vector((0.0, math.cos(t), math.sin(t)))
                n = Vector((0.0, -math.sin(t), math.cos(t)))
                quad = [hub + d * 0.3 * pr - n * chord / 2, hub + d * blade - n * chord / 3,
                        hub + d * blade + n * chord / 3, hub + d * 0.3 * pr + n * chord / 2]
                lo = [bm.verts.new(v - Vector((thick / 2, 0.0, 0.0))) for v in quad]
                hi = [bm.verts.new(v + Vector((thick / 2, 0.0, 0.0))) for v in quad]
                bm.faces.new(lo[::-1])
                bm.faces.new(hi)
                for i in range(4):
                    bm.faces.new((lo[i], lo[(i + 1) % 4], hi[(i + 1) % 4], hi[i]))
            ro = Vector((s * length, side * r * math.cos(a), r * math.sin(a)))
            w = 0.35 * pl
            box(bm, (ro.x - w / 2, min(ro.y, c.y) - 0.5, min(ro.z, c.z) - 0.5), (ro.x + w / 2, max(ro.y, c.y) + 0.5,
                                                                              max(ro.z, c.z) + 0.5))
    return mesh_from(bm, 'pods', True)


def ship_parts(row):
    L, Dm = row['length_m'], row['diameter_m']
    vol = 0.475 * L * Dm ** 2
    fin_area = 0.25 * vol ** (2.0 / 3.0) / 4.0
    fin_span = math.sqrt(0.5 * fin_area)
    return dict(hull=hull_mesh(L, Dm), fins=fin_mesh(L, Dm, fin_span, fin_area), decks=deck_mesh(L, Dm),
                pods=pod_mesh(L, Dm))


def place_ship(col, M, parts, name, nose, heading=0.0):
    for key, mat in (('hull', 'hull'), ('fins', 'fin'), ('decks', 'deck'), ('pods', 'pod')):
        o = bpy.data.objects.new(f'{name}_{key}', parts[key])
        if not o.data.materials:
            o.data.materials.append(M[mat])
        col.objects.link(o)
        o.location = nose
        o.rotation_euler = (0.0, 0.0, heading)


# ------------------------------------------------------------------------------------------------- the classes

def build_berths(col, M, option):
    """The option's arms and ships, from its berth layout in the study's product."""
    lay = option['crown']
    across = option['diameter_m'] + 2 * option['fin_span_m']
    r0 = radius(BAND) - 14.0
    parts = ship_parts(option)
    placed = 0
    for level, zc in enumerate(lay['levels_above_summit_m']):
        depth = max(30.0, min(0.15 * lay['arm_length_m'], 0.5 * across))
        for side in (-1.0, 1.0):
            if placed >= lay['berths']:
                break
            ships_here = min(lay['ships_per_arm'], lay['berths'] - placed)
            arm = lay['arm_length_m'] if ships_here == lay['ships_per_arm'] else \
                lay['pitch_m'] * (ships_here - 0.5) + FLEET['long_haul']['crown']['arm_beyond'] * across
            truss(col, M, f'arm{level}{int(side)}', side * (radius(zc) - 12.0), side * (r0 + arm), zc - LIFT, depth,
                  26.0)
            for k in range(ships_here):
                y = side * (r0 + lay['pitch_m'] * (k + 0.5))
                bm = bmesh.new()
                box(bm, (13.0, y - 5.0, zc - LIFT - 5.0), (36.0, y + 5.0, zc - LIFT + 5.0))
                link(f'boom{level}{k}{int(side)}', mesh_from(bm, 'boom'), M['dock'], col)
                place_ship(col, M, parts, f'ship{level}_{int(side)}_{k}', (36.0, y, zc - LIFT))
                placed += 1


def camera(scene, name, location, target, lens=50.0, ortho=None):
    cam = bpy.data.cameras.new(name)
    cam.clip_start, cam.clip_end = 20.0, 80000.0
    if ortho:
        cam.type = 'ORTHO'
        cam.ortho_scale = ortho
    else:
        cam.lens = lens
    o = bpy.data.objects.new(name, cam)
    o.location = location
    o.rotation_euler = (Vector(target) - Vector(location)).to_track_quat('-Z', 'Y').to_euler()
    scene.collection.objects.link(o)
    return o


def extent(options):
    """Points on every class's ships and arms and on the crown, in scene coordinates, for framing."""
    r0 = radius(BAND) - 14.0
    pts = [Vector((sx * r0, sy * r0, z)) for sx in (-1, 1) for sy in (-1, 1) for z in (BAND - LIFT - 40.0, H - LIFT + 10.0)]
    for o in options:
        lay = o['crown']
        across = o['diameter_m'] + 2 * o['fin_span_m']
        placed = 0
        for level, zl in enumerate(lay['levels_above_summit_m']):
            zc = zl - LIFT
            for side in (-1.0, 1.0):
                n = min(lay['ships_per_arm'], lay['berths'] - placed)
                for k in range(max(n, 0)):
                    y = side * (r0 + lay['pitch_m'] * (k + 0.5))
                    pts += [Vector((36.0, y, zc)), Vector((36.0 + o['length_m'], y, zc)),
                            Vector((36.0 + 0.36 * o['length_m'], y, zc + 0.5 * across)),
                            Vector((36.0 + 0.36 * o['length_m'], y, zc - 0.5 * across))]
                placed += max(n, 0)
    return pts


def framed_camera(scene, options, zc, az_deg=-40.0, el_deg=25.0, lens=70.0, fill=0.92):
    """One camera from the south-east for every class, drawn back until all their berths fit the frame."""
    from bpy_extras.object_utils import world_to_camera_view
    pts = extent(options)
    target = sum(pts, Vector()) / len(pts)
    target.z -= 60.0
    az, el = math.radians(az_deg), math.radians(el_deg)
    look = Vector((math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)))
    cam = camera(scene, 'view', target + 4000.0 * look, target, lens)
    scene.render.resolution_x, scene.render.resolution_y = 2400, 1500
    for _ in range(12):
        bpy.context.view_layer.update()
        uv = [world_to_camera_view(scene, cam, p) for p in pts]
        span = max(max(abs(u.x - 0.5) for u in uv) / 0.5, max(abs(u.y - 0.5) for u in uv) / 0.5)
        d = (cam.location - target).length * span / fill
        cam.location = target + d * look
    return cam


def render(scene, name, cam, size=(2400, 1500)):
    if ONLY and name not in ONLY and name.split('_')[0] not in ONLY:
        return
    scene.camera = cam
    scene.render.resolution_x, scene.render.resolution_y = size
    scene.render.resolution_percentage = 25 if SMALL else 100
    scene.render.filepath = str(OUT / f'{name}.png')
    bpy.ops.render.render(write_still=True)


def main():
    scene = setup_scene()
    M = materials()
    options = FLEET['long_haul']['options']
    tower = bpy.data.collections.new('tower')
    scene.collection.children.link(tower)
    build_tower(tower, M)
    zc = BAND + 250.0 - LIFT
    views = dict(view=framed_camera(scene, options, zc))
    layouts = {}
    for option in options:
        col = bpy.data.collections.new(f"berths_{option['length_m']:.0f}")
        scene.collection.children.link(col)
        build_berths(col, M, option)
        layouts[option['length_m']] = col
    for L, col in layouts.items():
        for other in layouts.values():
            other.hide_render = other is not col
        for view, cam in views.items():
            render(scene, f'{view}_{L:.0f}', cam)
    for other in layouts.values():
        other.hide_render = True
    tower.hide_render = True
    # The side elevation: every class beside the Hindenburg, noses aligned, each on its own row.
    line = bpy.data.collections.new('lineup')
    scene.collection.children.link(line)
    rows = [HINDENBURG] + options
    tallest = max(r['diameter_m'] for r in rows)
    z, centres = 0.0, []
    for i, row in enumerate(rows):
        parts = ship_parts(row)
        z -= row['diameter_m'] / 2
        place_ship(line, M, parts, f'line{i}', (0.0, 0.0, z))
        centres.append(z)
        z -= row['diameter_m'] / 2 + 0.3 * tallest
    span = max(r['length_m'] for r in rows)
    height = -z
    scale = span * 1.42
    size = (2400, int(round(2400 * height * 1.08 / scale / 10.0)) * 10)
    cx = span / 2 - 0.11 * scale
    cam = camera(scene, 'lineup', (cx, -2500.0, z / 2), (cx, 0.0, z / 2), ortho=scale)
    render(scene, 'lineup', cam, size)
    from bpy_extras.object_utils import world_to_camera_view
    scene.render.resolution_x, scene.render.resolution_y = size
    bpy.context.view_layer.update()
    labels = []
    for row, zc in zip(rows, centres):
        uv = world_to_camera_view(scene, cam, Vector((0.0, 0.0, zc)))
        labels.append(dict(length_m=row['length_m'], x=uv.x * size[0], y=(1.0 - uv.y) * size[1]))
    (OUT / 'lineup_labels.json').write_text(json.dumps(dict(size=size, rows=labels)))


main()
