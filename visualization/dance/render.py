"""Blender side of the dance renderer: a skinned dancer on a lit stage, driven by a motion file.

    blender -b --threads 3 --python render.py -- build  --blend DANCER.blend --report-dir DIR
    blender -b --threads 3 --python render.py -- render --blend DANCER.blend --motion M.npz --out OUT.mp4 \
            [--garment G.npz ...] [--preset test|final] [--stills 2.0,4.4] [--save-scene S.blend] \
            [--bodice R,G,B [--bodice-look matte|satin]] [--yoke R,G,B,A]

Use run_render.py, which supplies these arguments. Blender 4.0.2, EEVEE (legacy).

build: imports Blender Studio's CC0 realistic female base mesh (Human Base Meshes v1.4.1), turns it into the dance
frame (x forward, y left, z up, soles on z = 0), measures its joint landmarks and compares them with body.json
(landmark_check.json, and proposal_landmarks.json when a difference exceeds 2 cm), weights it to a fit armature
placed at its own joints (Blender's heat weights, then rule-based clean-up), applies two multires levels, reposes it
into body.json's rest pose with the limbs carried onto body.json's joint centres, and binds it to an armature built
from body.json alone. Costume zones (leotard, slippers and the sheer yoke's zone), hair, eyes and materials are added
and the rigged dancer is saved.

render: opens the rigged dancer, builds the stage (black Marley floor, cyclorama with a graded wash, ground row,
booms, soft key, rim, top light, haze), keys every bone from the motion's segment poses sampled at the render rate,
plays optional cloth sequences through Mesh Cache (PC2), sets a smoothly following camera and renders stills and/or an
H.264 MP4. For design stills it can dress the leotard zone as a bodice (set_bodice) and lay a sheer chiffon yoke over
the skin above it (set_yoke).
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from pathlib import Path

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BODY_JSON = ROOT / 'humanities' / 'dance' / 'body' / 'body.json'
ASSET = Path('/media/projectspace/terluna-research/humanities/assets/human-base-meshes-bundle-v1.4.1/'
             'human_base_meshes_bundle.blend')
BODY_NAME = 'GEO-body_female_realistic'
MULTIRES_LEVEL = 2          # 10.6 k base vertices -> about 170 k; level 3 (677 k) adds skin detail a full-body shot cannot show
RENDER_FPS = 30
# draft is the working size for whole videos until a complete video is near; stills check detail at final
RESOLUTION = dict(draft=(854, 480), test=(1280, 720), detail=(1280, 720), video=(1920, 1080), final=(1920, 1080))

# Face sets of the base mesh (left side; the right side is the mirrored set). The artist's sets end at the joints.
FS = dict(torso=1, head=17, abdomen=19, pelvis=18,
          upper_arm_L=21, upper_arm_R=20, forearm_L=12, forearm_R=11, hand_L=9, hand_R=10,
          thigh_L=24, thigh_R=23, shank_L=15, shank_R=16, foot_L=14, foot_R=13, ear_L=4, ear_R=5)
# The artist's head set (17) also covers the neck, and the jaw and mouth are a set of their own (22), so the head and
# neck are weighted by height up the neck instead: planes z = a + b x on the base mesh (x forward, metres) through
# the neck's base (the notch in front at 1.335 m, the seventh vertebra behind at 1.395 m) and through the jaw line
# (under the chin at 1.40 m, the skull's base behind at 1.46 m).
NECK_BASE = (1.366, -0.484)
JAW_LINE = (1.443, -0.43)


def neck_height(V):
    """Height up the neck: 0 at its base, 1 at the jaw line and the skull's base; above 1 is the head."""
    zn = NECK_BASE[0] + NECK_BASE[1] * V[:, 0]
    zb = JAW_LINE[0] + JAW_LINE[1] * V[:, 0]
    return (V[:, 2] - zn) / np.maximum(zb - zn, 1e-3)


# ----------------------------------------------------------------------------------------------- small helpers
def log(*a):
    print('[dance-render]', *a, flush=True)


def smoothstep(e0, e1, x):
    t = np.clip((np.asarray(x, float) - e0) / (e1 - e0), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def binned(values, coord, centres, half, func):
    """func(values in each bin of coord), NaN for empty bins."""
    out = np.full(len(centres), np.nan)
    for i, c in enumerate(centres):
        m = np.abs(coord - c) < half
        if m.any():
            out[i] = func(values[m])
    return out


def body_spec():
    return json.loads(BODY_JSON.read_text())


def seg_targets(spec):
    """Head and tail of every bone from body.json: tails run to the child joint, or along the segment's tip."""
    seg = {s['name']: s for s in spec['segments']}
    org = {n: np.array(s['origin'], float) for n, s in seg.items()}
    tail_child = dict(pelvis='abdomen', abdomen='thorax', thorax='neck', neck='head', L_clavicle='L_upper_arm',
                      R_clavicle='R_upper_arm', L_upper_arm='L_forearm', R_upper_arm='R_forearm',
                      L_forearm='L_hand', R_forearm='R_hand', L_thigh='L_shank', R_thigh='R_shank',
                      L_shank='L_foot', R_shank='R_foot', L_foot='L_toes', R_foot='R_toes')
    heads, tails = {}, {}
    for n, s in seg.items():
        heads[n] = org[n]
        tails[n] = org[tail_child[n]] if n in tail_child else org[n] + np.array(s['tip'], float)
    return heads, tails


def quat_to_mat(q):
    """(..., 4) w,x,y,z -> (..., 3, 3)."""
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    w, x, y, z = np.moveaxis(q, -1, 0)
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], -1),
                     np.stack([2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], -1),
                     np.stack([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1)], -2)


def quat_mul(p, q):
    w1, x1, y1, z1 = np.moveaxis(p, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(q, -1, 0)
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], -1)


def rot_between(a, b):
    """Smallest rotation (3x3) taking direction a onto direction b."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(a @ b)
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * (1 / (1 + c))


def set_active(obj, select=()):
    bpy.ops.object.select_all(action='DESELECT')
    for o in select:
        o.select_set(True)
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj


def mesh_coords(obj, evaluated=False):
    if evaluated:
        dg = bpy.context.evaluated_depsgraph_get()
        ev = obj.evaluated_get(dg)
        me = ev.to_mesh()
        co = np.empty(len(me.vertices) * 3)
        me.vertices.foreach_get('co', co)
        ev.to_mesh_clear()
        return co.reshape(-1, 3)
    co = np.empty(len(obj.data.vertices) * 3)
    obj.data.vertices.foreach_get('co', co)
    return co.reshape(-1, 3)


def weights_matrix(obj, names):
    """(V, B) dense weights for the named vertex groups."""
    idx = {g.index: i for i, g in enumerate(obj.vertex_groups) if g.name in names}
    order = {obj.vertex_groups[n].index: names.index(n) for n in names if n in obj.vertex_groups}
    W = np.zeros((len(obj.data.vertices), len(names)))
    for v in obj.data.vertices:
        for g in v.groups:
            j = order.get(g.group)
            if j is not None:
                W[v.index, j] = g.weight
    return W


def write_weights(obj, names, W, eps=1e-3):
    for n in names:
        if n in obj.vertex_groups:
            obj.vertex_groups.remove(obj.vertex_groups[n])
    for j, n in enumerate(names):
        g = obj.vertex_groups.new(name=n)
        nz = np.nonzero(W[:, j] > eps)[0]
        vals = W[nz, j]
        # group by rounded weight to keep the number of add() calls small
        r = np.round(vals, 3)
        for val in np.unique(r):
            g.add(nz[r == val].tolist(), float(val), 'REPLACE')


# ----------------------------------------------------------------------------------------------- build: import
def import_body():
    with bpy.data.libraries.load(str(ASSET), link=False) as (src, dst):
        dst.collections = ['Body Female - Realistic']
    coll = dst.collections[0]
    bpy.context.scene.collection.children.link(coll)
    bpy.context.view_layer.update()
    body = bpy.data.objects[BODY_NAME]
    eyes = [bpy.data.objects[BODY_NAME + '.eye.L'], bpy.data.objects[BODY_NAME + '.eye.R']]
    for e in eyes:
        mw = e.matrix_world.copy()
        e.parent = None
        e.matrix_world = mw
    bpy.context.view_layer.update()
    # Into the dance frame: the file's body faces -y with its left on +x, centred at x = -1.3401.
    zmin = float(mesh_coords(body)[:, 2].min())
    M = Matrix.Rotation(math.pi / 2, 4, 'Z') @ Matrix.Translation((1.3401, 0.0, -zmin))
    for o in [body] + eyes:
        o.matrix_world = M @ o.matrix_world
        bpy.context.view_layer.update()
        o.data.transform(o.matrix_world)
        o.matrix_world = Matrix.Identity(4)
    return body, eyes


def face_sets(body):
    a = body.data.attributes.get('.sculpt_face_set')
    fs = np.zeros(len(body.data.polygons), int)
    if a is not None:
        a.data.foreach_get('value', fs)
    return fs


def polys(body):
    return [list(p.vertices) for p in body.data.polygons]


# ----------------------------------------------------------------------------------------------- build: landmarks
def measure_landmarks(body):
    """Joint landmarks of the A-posed mesh in the dance frame (left side measured, right side mirrored).

    Rules (each an estimate; uncertainty about 1-2 cm):
      elbow, wrist, knee, ankle: centroid of the artist's face-set boundary loop at that joint (the closest
        approach of the upper-arm and forearm axes and the patella apex are recorded as cross-checks);
      shoulder (glenohumeral centre): 4 cm below the acromion (body.json's own rule) and 2.5 cm medial to it, the
        acromion taken as the lateral end of the shoulder's top contour where it turns steeper than 45 deg;
      hip joint centre: at the height of the centroid of the artist's leg-hole loop (the greatest hip breadth lies at
        the upper thighs on this body, so the trochanteric rule fails), a quarter of the greatest hip breadth from the
        midline, and at the front-back centre of the thigh section just below the crotch;
      MTP (ball of the foot): where the foot is broadest in its front half, at 45 % of the foot's height there;
      foot length: heel to longest toe along the foot's axis; sole: the lowest vertex (z = 0 after import).
    """
    V0 = mesh_coords(body)
    P = polys(body)
    fs = face_sets(body)
    vsets = {}
    for f, s in zip(P, fs):
        vsets.setdefault(int(s), set()).update(f)
    mod = body.modifiers['Multires']
    mod.levels = MULTIRES_LEVEL
    V2 = mesh_coords(body, evaluated=True)

    def loop(a, b):
        s = sorted(vsets[a] & vsets[b])
        return V0[s].mean(0)

    def axis(Q, a, b, lo, hi, nb=10):
        d = (b - a) / np.linalg.norm(b - a)
        L = np.linalg.norm(b - a)
        t = (Q - a) @ d
        cs = []
        for s in np.linspace(lo, hi, nb):
            m = np.abs(t - s * L) < 0.012 * L + 0.003
            if m.sum() > 4:
                cs.append(Q[m].mean(0))
        cs = np.array(cs)
        c = cs.mean(0)
        dd = np.linalg.svd(cs - c)[2][0]
        return c, (dd if dd @ d > 0 else -dd)

    lm, notes = {}, {}
    # --- arm
    sh_loop = loop(FS['upper_arm_L'], FS['torso'])
    el_loop = loop(FS['upper_arm_L'], FS['forearm_L'])
    wr_loop = loop(FS['forearm_L'], FS['hand_L'])
    ua_pts = V0[sorted(vsets[FS['upper_arm_L']])]
    fa_pts = V0[sorted(vsets[FS['forearm_L']])]
    cu, du = axis(ua_pts, sh_loop + np.array([0, 0.04, 0]), el_loop, 0.35, 0.9)
    cf, df = axis(fa_pts, el_loop, wr_loop, 0.15, 0.85)
    # acromion: top contour of the left shoulder on the level-2 surface
    S = V2[(V2[:, 1] > 0.08) & (V2[:, 1] < 0.27) & (V2[:, 2] > 1.15) & (np.abs(V2[:, 0]) < 0.09)]
    log('level-%d surface: %d vertices' % (MULTIRES_LEVEL, len(V2)))
    ys = np.arange(0.10, 0.26, 0.005)
    ztop = np.array([S[np.abs(S[:, 1] - y) < 0.004][:, 2].max() if (np.abs(S[:, 1] - y) < 0.004).any() else np.nan
                     for y in ys])
    ok = ~np.isnan(ztop)
    ys, ztop = ys[ok], ztop[ok]
    slope = np.diff(ztop) / np.diff(ys)
    i = int(np.argmax((slope < -1.0) & (ys[:-1] > 0.15)))
    acro = np.array([0.0, ys[i], ztop[i]])
    on_axis = cu + du * ((acro[2] - 0.04 - cu[2]) / du[2])
    gh = np.array([on_axis[0], acro[1] - 0.025, acro[2] - 0.04])
    # closest approach of the two arm axes
    w0 = cu - cf
    a, b, c_, d_, e_ = du @ du, du @ df, df @ df, du @ w0, df @ w0
    den = a * c_ - b * b
    bend = 0.5 * ((cu + du * (b * e_ - c_ * d_) / den) + (cf + df * (a * e_ - b * d_) / den))
    elbow = el_loop
    hand_pts = V0[sorted(vsets[FS['hand_L']] | set().union(*[vsets[k] for k in range(64, 84) if k in vsets]))]
    hand_tip = wr_loop + df * float(((hand_pts - wr_loop) @ df).max())
    lm.update(L_shoulder=gh, L_elbow=elbow, L_wrist=wr_loop, L_hand_tip=hand_tip, L_acromion=acro)
    notes['L_elbow'] = dict(loop=el_loop.round(4).tolist(), axes_bend_crosscheck=bend.round(4).tolist(),
                            spread_m=round(float(np.linalg.norm(el_loop - bend)), 4))
    notes['L_acromion'] = acro.round(4).tolist()
    notes['L_shoulder_on_upper_arm_axis_crosscheck'] = on_axis.round(4).tolist()
    notes['arm_abduction_deg'] = round(float(np.degrees(np.arccos(-du[2]))), 1)
    # --- hip
    pel = V0[sorted(vsets[FS['pelvis']])]
    crotch = float(pel[np.abs(pel[:, 1]) < 0.012][:, 2].min())
    H = V2[(V2[:, 2] > crotch) & (V2[:, 2] < crotch + 0.16) & (np.abs(V2[:, 1]) < 0.20)]
    zs = np.arange(crotch, crotch + 0.16, 0.004)
    br = binned(np.abs(H[:, 1]), H[:, 2], zs, 0.003, np.max)
    j = int(np.nanargmax(br))
    th = V2[(V2[:, 1] > 0.02) & (V2[:, 1] < 0.26) & (np.abs(V2[:, 2] - (crotch - 0.03)) < 0.003)]
    leg_hole = loop(FS['thigh_L'], FS['pelvis'])
    hip = np.array([th[:, 0].mean(), br[j] / 2, leg_hole[2]])
    lm.update(L_hip=hip)
    notes['L_hip_leg_hole_loop'] = leg_hole.round(4).tolist()
    notes['greatest_hip_breadth_z'] = round(float(zs[j]), 4)
    notes['crotch_z'] = round(crotch, 4)
    notes['hip_breadth_m'] = round(float(2 * br[j]), 4)
    # --- knee, ankle
    lm['L_knee'] = loop(FS['thigh_L'], FS['shank_L'])
    lm['L_ankle'] = loop(FS['shank_L'], FS['foot_L'])
    # patella apex as a cross-check of the knee loop
    K = V2[(np.abs(V2[:, 1] - lm['L_knee'][1]) < 0.025) & (V2[:, 2] > 0.38) & (V2[:, 2] < 0.58)]
    log('knee band points', len(K), 'knee loop', lm['L_knee'].round(4))
    zk = np.arange(0.385, 0.575, 0.005)
    xf = binned(K[:, 0], K[:, 2], zk, 0.004, np.max)
    ok = ~np.isnan(xf)
    if ok.sum() > 5:
        trend = np.polyval(np.polyfit(zk[ok], xf[ok], 1), zk)
        notes['patella_apex_z'] = round(float(zk[int(np.nanargmax(xf - trend))]), 4)
    # --- foot
    F = V2[(V2[:, 2] < 0.05) & (V2[:, 1] > 0.03) & (V2[:, 1] < 0.26)]
    c2 = F[:, :2].mean(0)
    fd = np.linalg.svd(F[:, :2] - c2)[2][0]
    fd = fd if fd[0] > 0 else -fd
    fl = np.array([-fd[1], fd[0]])
    t = (F[:, :2] - c2) @ fd
    heel, tip = t.min(), t.max()
    tb = np.arange(heel + 0.4 * (tip - heel), tip - 0.03, 0.004)
    wid = binned((F[:, :2] - c2) @ fl, t, tb, 0.003, np.ptp)
    k = int(np.nanargmax(wid))
    sec = F[np.abs(t - tb[k]) < 0.002]
    top = V2[(np.abs((V2[:, :2] - c2) @ fd - tb[k]) < 0.003) & (np.abs((V2[:, :2] - c2) @ fl) < 0.05)
             & (V2[:, 2] < 0.12)][:, 2].max()
    mtp_xy = c2 + fd * tb[k] + fl * float(((sec[:, :2] - c2) @ fl).mean())
    lm['L_mtp'] = np.array([mtp_xy[0], mtp_xy[1], 0.45 * top])
    toe_xy = c2 + fd * tip + fl * float(((F[t > tip - 0.01][:, :2] - c2) @ fl).mean())
    lm['L_toe_tip'] = np.array([toe_xy[0], toe_xy[1], 0.012])
    lm['L_heel_back'] = np.r_[c2 + fd * heel, 0.03]
    notes['foot_length_m'] = round(float(tip - heel), 4)
    notes['foot_axis_turnout_deg'] = round(float(np.degrees(np.arctan2(fd[1], fd[0]))), 1)
    notes['sole_z'] = round(float(V2[:, 2].min()), 4)
    # --- head and trunk
    ear = V0[sorted(vsets[FS['ear_L']])].mean(0)
    lm['ear'] = np.array([ear[0], 0.0, ear[2]])
    lm['head_top'] = V2[np.argmax(V2[:, 2])]
    mid = V2[np.abs(V2[:, 1]) < 0.01]

    def depth(z):
        m = mid[np.abs(mid[:, 2] - z) < 0.004]
        return float(m[:, 0].min()), float(m[:, 0].max())
    lm['_depth'] = depth
    lm['stature'] = float(V2[:, 2].max() - V2[:, 2].min())
    mod.levels = 0
    return lm, notes


def mirror(p):
    return np.array([p[0], -p[1], p[2]])


def mesh_skeleton(lm, spec):
    """Bone heads and tails at the mesh's own joints (A-pose), for weighting and reposing."""
    heads_j, tails_j = seg_targets(spec)
    depth = lm['_depth']
    h, t = {}, {}
    hipL = lm['L_hip']
    h['pelvis'] = np.array([hipL[0], 0.0, hipL[2]])
    for n in ('abdomen', 'thorax', 'neck'):
        z = heads_j[n][2]
        xb, xf = depth(z)
        h[n] = np.array([xb + 0.40 * (xf - xb), 0.0, z])     # vertebral line at 40 % of the trunk's depth
    oc = lm['ear'] + np.array([0.0, 0.0, -0.025])
    h['head'] = oc
    t['head'] = np.array([oc[0], 0.0, lm['head_top'][2] - 0.02])
    xb, xf = depth(heads_j['L_clavicle'][2])
    for side, f in (('L', lambda p: p), ('R', mirror)):
        h[f'{side}_clavicle'] = f(np.array([xf - 0.02, 0.02, heads_j['L_clavicle'][2]]))
        h[f'{side}_upper_arm'] = f(lm['L_shoulder'])
        h[f'{side}_forearm'] = f(lm['L_elbow'])
        h[f'{side}_hand'] = f(lm['L_wrist'])
        t[f'{side}_hand'] = f(lm['L_hand_tip'])
        h[f'{side}_thigh'] = f(lm['L_hip'])
        h[f'{side}_shank'] = f(lm['L_knee'])
        h[f'{side}_foot'] = f(lm['L_ankle'])
        h[f'{side}_toes'] = f(lm['L_mtp'])
        t[f'{side}_toes'] = f(lm['L_toe_tip'])
    chain = dict(pelvis='abdomen', abdomen='thorax', thorax='neck', neck='head')
    for side in 'LR':
        chain.update({f'{side}_clavicle': f'{side}_upper_arm', f'{side}_upper_arm': f'{side}_forearm',
                      f'{side}_forearm': f'{side}_hand', f'{side}_thigh': f'{side}_shank',
                      f'{side}_shank': f'{side}_foot', f'{side}_foot': f'{side}_toes'})
    for a, b in chain.items():
        t[a] = h[b]
    return h, t


def rest_equivalent(lm, spec):
    """Where the mesh's joints would sit in body.json's rest pose: the mesh's root joints and segment lengths carried
    along body.json's rest directions."""
    H, T = seg_targets(spec)
    out = {}
    for side, f in (('L', lambda p: p), ('R', mirror)):
        sh = f(lm['L_shoulder'])
        lu = np.linalg.norm(lm['L_elbow'] - lm['L_shoulder'])
        lf = np.linalg.norm(lm['L_wrist'] - lm['L_elbow'])
        du = (H[f'{side}_forearm'] - H[f'{side}_upper_arm'])
        dfr = (H[f'{side}_hand'] - H[f'{side}_forearm'])
        el = sh + lu * du / np.linalg.norm(du)
        wr = el + lf * dfr / np.linalg.norm(dfr)
        hp = f(lm['L_hip'])
        lt = np.linalg.norm(lm['L_knee'] - lm['L_hip'])
        ls = np.linalg.norm(lm['L_ankle'] - lm['L_knee'])
        dt = H[f'{side}_shank'] - H[f'{side}_thigh']
        ds = H[f'{side}_foot'] - H[f'{side}_shank']
        kn = hp + lt * dt / np.linalg.norm(dt)
        an = kn + ls * ds / np.linalg.norm(ds)
        out.update({f'{side}_shoulder': (sh, H[f'{side}_upper_arm']), f'{side}_elbow': (el, H[f'{side}_forearm']),
                    f'{side}_wrist': (wr, H[f'{side}_hand']), f'{side}_hip': (hp, H[f'{side}_thigh']),
                    f'{side}_knee': (kn, H[f'{side}_shank']), f'{side}_ankle': (an, H[f'{side}_foot'])})
    return out


def landmark_report(lm, notes, spec, report_dir: Path):
    H, T = seg_targets(spec)
    rows = {}
    for name, (m, j) in rest_equivalent(lm, spec).items():
        d = m - j
        rows[name] = dict(mesh=m.round(4).tolist(), body_json=j.round(4).tolist(), diff=d.round(4).tolist(),
                          dist_m=round(float(np.linalg.norm(d)), 4))
    seglen = {}
    for nm, a, b, ja, jb in (('upper_arm', 'L_shoulder', 'L_elbow', 'L_upper_arm', 'L_forearm'),
                             ('forearm', 'L_elbow', 'L_wrist', 'L_forearm', 'L_hand'),
                             ('thigh', 'L_hip', 'L_knee', 'L_thigh', 'L_shank'),
                             ('shank', 'L_knee', 'L_ankle', 'L_shank', 'L_foot'),
                             ('foot_ankle_to_mtp', 'L_ankle', 'L_mtp', 'L_foot', 'L_toes')):
        m = float(np.linalg.norm(lm[b] - lm[a]))
        j = float(np.linalg.norm(H[jb] - H[ja]))
        seglen[nm] = dict(mesh=round(m, 4), body_json=round(j, 4), diff=round(m - j, 4))
    m = float(np.linalg.norm(lm['L_hand_tip'] - lm['L_wrist']))
    seglen['hand_wrist_to_fingertip'] = dict(mesh=round(m, 4), body_json=round(float(np.linalg.norm(T['L_hand'] - H['L_hand'])), 4))
    seglen['hand_wrist_to_fingertip']['diff'] = round(m - seglen['hand_wrist_to_fingertip']['body_json'], 4)
    j_toe = T['L_toes']
    heel_j = H['L_foot'] + np.array([-0.04, 0, -0.0696])     # body.json heel contact point (skin behind it)
    foot = dict(mesh_length_heel_to_toe=notes['foot_length_m'],
                body_json_heel_contact_to_toe_tip=round(float(j_toe[0] - heel_j[0]), 4),
                note='body.json gives contact points, not the heel skin; its toe tip is the toes segment tip')
    sole = dict(mesh_ankle_height_above_sole=round(float(lm['L_ankle'][2] - notes['sole_z']), 4),
                body_json_ankle_height_above_sole=round(float(H['L_foot'][2] - (H['L_foot'][2] - 0.069555)), 4),
                mesh_mtp_height=round(float(lm['L_mtp'][2]), 4), body_json_mtp_height=round(float(H['L_toes'][2]), 4))
    depth = lm['_depth']
    spine = {}
    for n in ('abdomen', 'thorax', 'neck', 'head', 'L_clavicle'):
        p = H[n]
        xb, xf = depth(p[2])
        spine[n] = dict(body_json=p.round(4).tolist(), mesh_back_skin_x=round(xb, 4), mesh_front_skin_x=round(xf, 4),
                        depth_fraction_from_back=round((p[0] - xb) / (xf - xb), 3))
    head_oc = lm['ear'] + np.array([0, 0, -0.025])
    spine['head']['mesh_occipital_condyles_estimate'] = head_oc.round(4).tolist()
    spine['head']['note'] = 'estimate: 2.5 cm below the centre of the ear'
    out = dict(
        schema='terluna.humanities.dance-landmark-check/1',
        mesh='Blender Studio Human Base Meshes v1.4.1, Body Female - Realistic (CC0), multires level %d' % MULTIRES_LEVEL,
        frame='dance frame (x forward, y left, z up), soles on z = 0, mesh x/y as imported (bbox-centred)',
        method=measure_landmarks.__doc__.strip(),
        comparison='mesh joints carried into body.json rest pose: the mesh root joints (shoulder, hip) with the mesh '
                   'segment lengths laid along body.json rest directions',
        stature=dict(mesh=round(lm['stature'], 4), body_json=spec['stature_m']),
        joints=rows, segment_lengths=seglen, foot=foot, sole=sole, trunk=spine, notes=notes,
        threshold_m=0.02,
        exceeding=[k for k, r in rows.items() if r['dist_m'] > 0.02]
        + [f'segment:{k}' for k, r in seglen.items() if abs(r['diff']) > 0.02]
        + [f'trunk:{k}' for k, r in spine.items() if r['depth_fraction_from_back'] < 0.15])
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / 'landmark_check.json').write_text(json.dumps(out, indent=1))
    log('landmark check exceeding 2 cm:', out['exceeding'])
    return out


# ----------------------------------------------------------------------------------------------- build: armatures
def make_armature(name, heads, tails, parents=None):
    arm = bpy.data.armatures.new(name)
    obj = bpy.data.objects.new(name, arm)
    bpy.context.scene.collection.objects.link(obj)
    set_active(obj)
    bpy.ops.object.mode_set(mode='EDIT')
    for n in heads:
        b = arm.edit_bones.new(n)
        b.head = Vector(heads[n])
        b.tail = Vector(tails[n])
        b.roll = 0.0
        b.use_deform = True
    if parents:
        for n, p in parents.items():
            if p:
                arm.edit_bones[n].parent = arm.edit_bones[p]
                arm.edit_bones[n].use_connect = False
    bpy.ops.object.mode_set(mode='OBJECT')
    arm.display_type = 'STICK'
    return obj


def auto_weights(body, arm):
    set_active(arm, select=[body])
    bpy.ops.object.parent_set(type='ARMATURE_AUTO')
    body.parent = None
    body.matrix_world = Matrix.Identity(4)
    for m in list(body.modifiers):
        if m.type == 'ARMATURE':
            body.modifiers.remove(m)


def fix_weights(body, names, mesh_heads):
    """Rule-based clean-up of the heat weights on the A-posed base mesh."""
    V = mesh_coords(body)
    fs_face = face_sets(body)
    P = polys(body)
    vfs = np.zeros(len(V), int)
    for f, s in zip(P, fs_face):
        vfs[f] = s
    W = weights_matrix(body, names)
    zero = W.sum(1) < 1e-6
    log('vertices without heat weight:', int(zero.sum()))
    if zero.any():   # nearest bone segment
        for i in np.nonzero(zero)[0]:
            best, bj = 1e9, 0
            for j, n in enumerate(names):
                a = mesh_heads[0][n]
                b = mesh_heads[1][n]
                d = b - a
                t = np.clip((V[i] - a) @ d / (d @ d), 0, 1)
                dist = np.linalg.norm(V[i] - (a + t * d))
                if dist < best:
                    best, bj = dist, j
            W[i, bj] = 1.0
    y, z = V[:, 1], V[:, 2]
    col = {n: j for j, n in enumerate(names)}
    for n, j in col.items():
        if n.startswith('L_'):
            W[:, j] *= smoothstep(-0.012, 0.004, y)
        if n.startswith('R_'):
            W[:, j] *= smoothstep(-0.012, 0.004, -y)
    arm_sets = {FS[k] for k in ('upper_arm_L', 'upper_arm_R', 'forearm_L', 'forearm_R', 'hand_L', 'hand_R')} \
        | set(range(64, 104))
    leg_sets = {FS[k] for k in ('thigh_L', 'thigh_R', 'shank_L', 'shank_R', 'foot_L', 'foot_R')} | set(range(25, 64))
    for s in 'LR':
        # arms own nothing on the legs and lower trunk; the upper arm reaches the torso only at the shoulder
        for part in ('upper_arm', 'forearm', 'hand'):
            j = col[f'{s}_{part}']
            W[np.isin(vfs, list(leg_sets | {FS['pelvis'], FS['abdomen'], FS['head']})), j] = 0
        j = col[f'{s}_upper_arm']
        tors = vfs == FS['torso']
        W[tors, j] *= smoothstep(0.135, 0.175, np.abs(y[tors])) * smoothstep(1.17, 1.26, z[tors])
        for part in ('forearm', 'hand'):
            j = col[f'{s}_{part}']
            W[~np.isin(vfs, list(arm_sets)), j] = 0
        # legs own nothing above the waist or on the arms
        for part in ('thigh', 'shank', 'foot', 'toes'):
            j = col[f'{s}_{part}']
            W[np.isin(vfs, list(arm_sets)), j] = 0
        W[:, col[f'{s}_thigh']] *= 1 - smoothstep(0.90, 0.97, z)
        W[:, col[f'{s}_shank']] *= 1 - smoothstep(0.56, 0.62, z)
        W[:, col[f'{s}_foot']] *= 1 - smoothstep(0.15, 0.20, z)
        W[:, col[f'{s}_toes']] *= 1 - smoothstep(0.05, 0.07, z)
        W[:, col[f'{s}_clavicle']] *= 1 - smoothstep(0.0, 0.05, 1.20 - z)
    W[np.isin(vfs, list(leg_sets | arm_sets)), col['head']] = 0
    W[:, col['head']] *= smoothstep(1.36, 1.42, z)
    W[:, col['pelvis']] *= 1 - smoothstep(1.02, 1.08, z)
    s = W.sum(1, keepdims=True)
    W = np.where(s > 1e-9, W / np.maximum(s, 1e-9), W)
    # head and neck: the whole face, jaw and mouth go with the skull; the neck's bend spreads along its length, from
    # the thorax at its base through the neck bone to the head at the jaw line, over a column 6-8.5 cm either side
    u = neck_height(V)
    limbs = np.isin(vfs, list(leg_sets | arm_sets))
    head_zone = (u >= 1.0) & ~limbs
    th, hd = 1 - smoothstep(0.0, 0.45, u), smoothstep(0.55, 1.0, u)
    N = np.zeros_like(W)
    N[:, col['thorax']], N[:, col['head']] = th, hd
    N[:, col['neck']] = 1 - th - hd
    f = (1 - smoothstep(0.055, 0.085, np.abs(y))) * smoothstep(-0.25, 0.0, u)
    f = np.where(head_zone, 1.0, np.where(limbs, 0.0, f))[:, None]
    W = (1 - f) * W + f * N
    write_weights(body, names, W)
    # Blender's own smoothing across the joints, then limit and normalise
    set_active(body)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.object.vertex_group_smooth(group_select_mode='ALL', factor=0.5, repeat=4, expand=0.0)
    bpy.ops.object.vertex_group_clean(group_select_mode='ALL', limit=0.01)
    bpy.ops.object.vertex_group_limit_total(group_select_mode='ALL', limit=4)
    bpy.ops.object.vertex_group_normalize_all(group_select_mode='ALL', lock_active=False)
    bpy.ops.object.mode_set(mode='OBJECT')
    # the head stays rigid above the jaw line
    W = weights_matrix(body, names)
    W[head_zone] = 0
    W[head_zone, col['head']] = 1
    write_weights(body, names, W)


def repose(body, eyes, fit, mesh_ht, target_ht, names):
    """Carry the A-posed mesh into body.json's rest pose: limbs rotated onto body.json's directions, stretched to its
    lengths along their axes, and placed on its joint centres; the trunk, neck and head keep the mesh's shape."""
    mh, mt = mesh_ht
    jh, jt = target_ht
    trunk = ('pelvis', 'abdomen', 'thorax', 'neck', 'head')
    set_active(fit)
    bpy.ops.object.mode_set(mode='POSE')
    for n in names:
        pb = fit.pose.bones[n]
        pb.rotation_mode = 'QUATERNION'
        if n in trunk:
            pb.matrix_basis = Matrix.Identity(4)
            continue
        a0, b0 = mh[n], mt[n]
        if n.endswith('clavicle'):
            a1, b1 = a0, jh[n.replace('clavicle', 'upper_arm')]
        else:
            a1, b1 = jh[n], jt[n]
        R = rot_between(b0 - a0, b1 - a1)
        s = np.linalg.norm(b1 - a1) / np.linalg.norm(b0 - a0)
        rest = np.array(fit.data.bones[n].matrix_local)
        Rr = rest[:3, :3]
        M = np.eye(4)
        M[:3, :3] = R @ Rr @ np.diag([1.0, s, 1.0])
        M[:3, 3] = a1
        pb.matrix_basis = Matrix((np.linalg.inv(rest) @ M).tolist())
    bpy.ops.object.mode_set(mode='OBJECT')
    for o in [body] + eyes:
        if o is not body:
            g = o.vertex_groups.new(name='head')
            g.add(list(range(len(o.data.vertices))), 1.0, 'REPLACE')
        m = o.modifiers.new('repose', 'ARMATURE')
        m.object = fit
        set_active(o)
        bpy.ops.object.modifier_apply(modifier=m.name)


def leotard_and_slipper_masks(body, names):
    """Costume zones as smooth vertex fields on the rest-posed mesh (rest coordinates, metres)."""
    V = mesh_coords(body)
    W = weights_matrix(body, names)
    col = {n: j for j, n in enumerate(names)}
    x, y, z = V.T
    ay = np.abs(y)
    arm = sum(W[:, col[f'{s}_{p}']] for s in 'LR' for p in ('upper_arm', 'forearm', 'hand'))
    headneck = W[:, col['head']] + W[:, col['neck']]
    crotch = 0.752
    front = smoothstep(-0.03, 0.03, x)
    k = 0.30 + (0.95 - 0.30) * front          # high cut in front, covering the seat behind
    zcut = crotch + k * np.maximum(ay - 0.022, 0)
    neck_front = 1.255 + 0.85 * ay
    neck_back = 1.215 + 1.1 * ay
    zneck = neck_back + (neck_front - neck_back) * front
    d = np.minimum.reduce([z - zcut, zneck - z, (0.32 - arm) * 0.08, (0.30 - headneck) * 0.08,
                           np.maximum(0.128 - ay, 1.16 - z)])
    E = mesh_edges(body)
    leo = smoothstep(0.35, 0.65, smooth_field(smoothstep(-0.0015, 0.0015, d), E, 6))
    foot = sum(W[:, col[f'{s}_{p}']] for s in 'LR' for p in ('foot', 'toes'))
    slip = smoothstep(0.35, 0.65, smooth_field(smoothstep(-0.0015, 0.0015, np.minimum(0.055 - z, (foot - 0.4) * 0.05)),
                                             E, 4))
    for nm, val in (('leotard', leo), ('slipper', slip)):
        a = body.data.attributes.get(nm) or body.data.attributes.new(nm, 'FLOAT', 'POINT')
        a.data.foreach_set('value', val.astype(np.float32))
    return leo, slip


# The sheer yoke's neckline follows the base of the neck. The neck leans forward, so it is taken as an elliptic
# cylinder about the axis x = 0.012 + 0.4 (z - 1.36) m (through the neck's base centre, tilted with the neck), with
# semi-axes 6.0 cm front to back and 7.0 cm side to side: on the rest mesh it meets the skin at the jugular notch
# (z 1.337 m), at the sides of the neck's base (1.39 m) and just below the seventh vertebra behind (1.39 m). The plane
# of neck_height alone would cut the shoulder tops 11 cm from the midline and leave the trapezius bare beside the neck.
NECK_AXIS = (0.012, 1.36, 0.4)        # x at z0, z0, dx/dz
NECK_SEMI = (0.060, 0.070)            # front-back, side-to-side (m)


def yoke_mask(body, names):
    """The sheer yoke as smooth vertex fields on the rest-posed mesh (rest coordinates, metres).

    'yoke' covers the skin from under the leotard up to a jewel neckline at the base of the neck, front and back,
    over the shoulder tops out to the shoulder point, and down to the leotard at the sides, its armholes on the
    leotard's own armhole line (summed arm weight 0.32, which on the shoulder top lies at |y| = 0.195 m, where the
    contour turns down the arm). Below its visible part it runs on under the leotard (from z = 1.12 m, the leotard's
    sides reaching 1.16 m), so the leotard's edge is the yoke's edge there. 'yoke_hem' is a 4 to 6 mm band inside the
    yoke's neckline and armholes, where the fabric is turned back on itself. The render's chiffon look reads both
    (set_yoke); without --yoke the skin shows as before."""
    V = mesh_coords(body)
    W = weights_matrix(body, names)
    col = {n: j for j, n in enumerate(names)}
    x, y, z = V.T
    arm = sum(W[:, col[f'{s}_{p}']] for s in 'LR' for p in ('upper_arm', 'forearm', 'hand'))
    x0, z0, tilt = NECK_AXIS
    a, b = NECK_SEMI
    r = np.sqrt(((x - x0 - tilt * (z - z0)) / a) ** 2 + (y / b) ** 2)
    nh = neck_height(V)
    d_neck = np.maximum((r - 1.0) * min(a, b), -(nh + 0.6) * 0.08)    # outside the neck, or well below its base
    d_head = (0.6 - nh) * 0.08                                         # never up the neck past its lower part
    d_arm = (0.32 - arm) * 0.08                                        # the leotard's armhole line
    d = np.minimum.reduce([d_neck, d_head, d_arm, z - 1.12])
    E = mesh_edges(body)
    yoke = smoothstep(0.35, 0.65, smooth_field(smoothstep(-0.0015, 0.0015, d), E, 6))
    hem = yoke * (1 - smoothstep(0.004, 0.006, np.minimum.reduce([d_neck, d_head, d_arm])))
    for nm, val in (('yoke', yoke), ('yoke_hem', hem)):
        at = body.data.attributes.get(nm) or body.data.attributes.new(nm, 'FLOAT', 'POINT')
        at.data.foreach_set('value', np.clip(val, 0, 1).astype(np.float32))
    log('yoke: %d vertices (%d of them clear of the leotard), hem %d' % (
        (yoke > 0.5).sum(), ((yoke > 0.5) & (body_leotard(body) < 0.5)).sum(), (hem > 0.5).sum()))
    return yoke


def body_leotard(body):
    a = body.data.attributes['leotard']
    v = np.empty(len(a.data))
    a.data.foreach_get('value', v)
    return v


def mesh_edges(obj):
    e = np.empty(len(obj.data.edges) * 2, int)
    obj.data.edges.foreach_get('vertices', e)
    return e.reshape(-1, 2)


def smooth_field(F, E, iters, mask=None, lam=0.5):
    """Laplacian smoothing of a per-vertex field over the mesh graph (only `mask` vertices move)."""
    n = len(F)
    deg = np.maximum(np.bincount(E.ravel(), minlength=n), 1).astype(float)
    S = np.array(F, float)
    one = S.ndim == 1
    if one:
        S = S[:, None]
    for _ in range(iters):
        avg = np.stack([np.bincount(E[:, 0], S[E[:, 1], c], n) + np.bincount(E[:, 1], S[E[:, 0], c], n)
                        for c in range(S.shape[1])], 1) / deg[:, None]
        if mask is None:
            S += lam * (avg - S)
        else:
            S[mask] += lam * (avg[mask] - S[mask])
    return S[:, 0] if one else S


def fabric_fill(body, leo, iters=60):
    """Stretch fabric bridges hollows: under the leotard, move the surface outward (never inward) toward a smoothed
    copy of itself, so the cleavage, navel and spine groove fill as knit fabric would."""
    V = mesh_coords(body)
    E = mesh_edges(body)
    N = np.empty(len(V) * 3)
    body.data.vertices.foreach_get('normal', N)
    N = N.reshape(-1, 3)
    region = leo > 0.02
    S = smooth_field(V, E, iters, mask=region)
    D = (S - V) * (leo * smoothstep(-0.002, 0.004, ((S - V) * N).sum(1)))[:, None]
    mag = np.linalg.norm(D, axis=1)
    D *= (0.006 * np.tanh(mag / 0.006) / np.maximum(mag, 1e-9))[:, None]     # saturate at about 6 mm
    D = smooth_field(D, E, 12)
    body.data.vertices.foreach_set('co', (V + D).ravel())
    body.data.update()
    d = np.linalg.norm(D, axis=1)
    log('fabric fill: max %.1f mm, mean over leotard %.2f mm' % (1e3 * d.max(), 1e3 * d[leo > 0.5].mean()))


def face_paint(body, eyes):
    """Lips (the artist's lip face sets) and brows (an arch above each eye) as soft vertex fields."""
    V = mesh_coords(body)
    E = mesh_edges(body)
    fs = face_sets(body)
    lips = np.zeros(len(V))
    for poly, s_ in zip(body.data.polygons, fs):
        if s_ in (7, 8):
            lips[list(poly.vertices)] = 1.0
    lips = smooth_field(lips, E, 3)
    brow = np.zeros(len(V))
    for e in eyes:
        c = mesh_coords(e).mean(0)
        dy = (V[:, 1] - c[1]) * np.sign(c[1])          # outward from the nose
        zc = c[2] + 0.021 + 0.005 * (1 - np.clip(((dy - 0.004) / 0.022) ** 2, 0, 1))
        band = smoothstep(0.0045, 0.002, np.abs(V[:, 2] - zc)) * smoothstep(-0.018, -0.012, dy) \
            * smoothstep(0.03, 0.024, dy) * smoothstep(c[0] - 0.02, c[0] - 0.008, V[:, 0])
        brow = np.maximum(brow, band)
    brow = smooth_field(brow, E, 2)
    for nm, val in (('lips', lips), ('brow', brow)):
        a = body.data.attributes.get(nm) or body.data.attributes.new(nm, 'FLOAT', 'POINT')
        a.data.foreach_set('value', np.clip(val, 0, 1).astype(np.float32))


def eye_iris_attribute(eye):
    V = mesh_coords(eye)
    c = V.mean(0)
    r = np.linalg.norm(V - c, axis=1).max()
    d = (V - c) / np.maximum(np.linalg.norm(V - c, axis=1, keepdims=True), 1e-9)
    ang = np.degrees(np.arccos(np.clip(d[:, 0], -1, 1)))      # angle from the forward (+x) axis
    a = eye.data.attributes.new('iris', 'FLOAT', 'POINT')
    a.data.foreach_set('value', ang.astype(np.float32))
    return c, r


def build_hair(body, names):
    """A hair cap from the scalp of the rest mesh (offset 3.5 mm, thickened) and a classical bun at the crown."""
    V = mesh_coords(body)
    W = weights_matrix(body, ['head'])
    x, y, z = V.T
    ay = np.abs(y)
    plane = 1.468 + 0.82 * (x + 0.045)                  # nape hairline at the back, forehead hairline in front
    temple = np.where(x > 0.035, 1.588 - 0.9 * np.maximum(ay - 0.03, 0), 0)
    ears = np.where((ay > 0.052) & (x > -0.03), 1.538, 0)
    hair_v = (W[:, 0] > 0.98) & (z > np.maximum.reduce([plane, temple, ears]))
    me = body.data
    me.calc_normals_split() if hasattr(me, 'calc_normals_split') else None
    N = np.empty(len(V) * 3)
    me.vertices.foreach_get('normal', N)
    N = N.reshape(-1, 3)
    keep = [p for p in me.polygons if all(hair_v[i] for i in p.vertices)]
    used = sorted({i for p in keep for i in p.vertices})
    remap = {o: n for n, o in enumerate(used)}
    faces = [[remap[i] for i in p.vertices] for p in keep]
    # the cap thins to 1.5 mm at the hairline over about 2 cm (hops from the boundary)
    from collections import Counter, deque
    ec = Counter(tuple(sorted(e)) for f in faces for e in zip(f, f[1:] + f[:1]))
    adj = [[] for _ in used]
    for (a_, b_) in ec:
        adj[a_].append(b_)
        adj[b_].append(a_)
    hops = np.full(len(used), 99)
    q = deque()
    for (a_, b_), c_ in ec.items():
        if c_ == 1:
            for v_ in (a_, b_):
                if hops[v_]:
                    hops[v_] = 0
                    q.append(v_)
    while q:
        v_ = q.popleft()
        for w_ in adj[v_]:
            if hops[w_] > hops[v_] + 1:
                hops[w_] = hops[v_] + 1
                q.append(w_)
    off = 0.0015 + 0.0045 * smoothstep(0, 7, hops)
    verts = (V[used] + off[:, None] * N[used]).tolist()
    edge = smoothstep(0, 6, hops)
    hm = bpy.data.meshes.new('hair_cap')
    hm.from_pydata(verts, [], faces)
    hm.update()
    cap = bpy.data.objects.new('hair_cap', hm)
    bpy.context.scene.collection.objects.link(cap)
    for p in hm.polygons:
        p.use_smooth = True
    sm = cap.modifiers.new('smooth', 'SMOOTH')
    sm.factor, sm.iterations = 0.6, 6
    so = cap.modifiers.new('thick', 'SOLIDIFY')
    so.thickness, so.offset = 0.004, 1.0
    sub = cap.modifiers.new('sub', 'SUBSURF')
    sub.levels = sub.render_levels = 1
    # bun: behind the crown
    back = V[(np.abs(y) < 0.01) & (np.abs(z - 1.575) < 0.01)][:, 0].min()
    bpy.ops.mesh.primitive_uv_sphere_add(segments=48, ring_count=24, radius=1.0,
                                         location=(back - 0.016, 0.0, 1.592))
    bun = bpy.context.active_object
    bun.name = 'hair_bun'
    bun.scale = (0.034, 0.048, 0.046)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    bpy.ops.object.shade_smooth()
    # comb lines: azimuth about the axis from the cranium's centre to the bun (meridians converge on the bun)
    C = np.array(verts)
    hc = 0.5 * (C.min(0) + C.max(0)) - np.array([0.0, 0.0, 0.02])
    bc = np.array(bun.location)
    ax_ = (bc - hc) / np.linalg.norm(bc - hc)
    e1 = np.array([0.0, 0.0, -1.0]) - ax_ * (-ax_[2])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(ax_, e1)
    d_ = C - hc
    phi = np.arctan2(d_ @ e2, d_ @ e1)
    at = hm.attributes.new('hair_phi', 'FLOAT', 'POINT')
    at.data.foreach_set('value', phi.astype(np.float32))
    ae = hm.attributes.new('hair_edge', 'FLOAT', 'POINT')
    ae.data.foreach_set('value', edge.astype(np.float32))
    for o in (cap, bun):
        g = o.vertex_groups.new(name='head')
        g.add(list(range(len(o.data.vertices))), 1.0, 'REPLACE')
    strands = build_strands(C, faces, edge, bc, np.array(bun.scale) if tuple(bun.scale) != (1, 1, 1)
                            else np.array((0.034, 0.048, 0.046)), ax_)
    return cap, bun, strands


def build_strands(Vc, Fc, edge, bc, brad, bun_axis, n_strands=16000, n_rings=2600, npts=10, seed=7):
    """Hair as curves: strands rooted over the cap, combed along the scalp to the bun (each step moves toward the
    bun in the scalp's tangent plane and is snapped back onto the cap at its own height), plus strands wound round
    the bun. Rooting thins toward the hairline. Returns a hair-curves object, or None if this Blender lacks the API."""
    from mathutils.bvhtree import BVHTree
    if not hasattr(bpy.data, 'hair_curves'):
        log('no hair curves API; keeping the cap only')
        return None
    rng = np.random.default_rng(seed)
    tree = BVHTree.FromPolygons([tuple(v) for v in Vc], [tuple(f) for f in Fc])
    tris = np.array([(f[0], f[k], f[k + 1]) for f in Fc for k in range(1, len(f) - 1)])
    A, B, Cc = Vc[tris[:, 0]], Vc[tris[:, 1]], Vc[tris[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(B - A, Cc - A), axis=1)
    w = area * (0.12 + 0.88 * edge[tris].mean(1))
    idx = rng.choice(len(tris), n_strands, p=w / w.sum())
    u = rng.random((n_strands, 2))
    flip = u.sum(1) > 1
    u[flip] = 1 - u[flip]
    roots = A[idx] + u[:, :1] * (B[idx] - A[idx]) + u[:, 1:] * (Cc[idx] - A[idx])
    lift = 0.0006 + 0.0032 * rng.random(n_strands) ** 1.5
    bias = rng.normal(0.0, 0.10, n_strands)
    rb = float(np.mean(brad))
    bcv = Vector(bc)
    paths = []
    for i in range(n_strands):
        loc, nor, _, _ = tree.find_nearest(Vector(roots[i]))
        pt = loc + nor * lift[i]
        pts = [pt.copy()]
        for _ in range(40):
            loc, nor, _, _ = tree.find_nearest(pt)
            d = bcv - pt
            if d.length < rb + 0.004:
                break
            t = d - nor * d.dot(nor)
            if t.length < 1e-6:
                break
            t.normalize()
            side = nor.cross(t)
            step = t + side * (bias[i] + rng.normal(0.0, 0.05))
            step.normalize()
            pt = pt + step * 0.006
            loc, nor, _, _ = tree.find_nearest(pt)
            pt = loc + nor * lift[i]
            pts.append(pt.copy())
        if len(pts) < 3:
            continue
        P = np.array([tuple(q) for q in pts])
        L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        tt = np.linspace(0, L[-1], npts)
        paths.append(np.stack([np.interp(tt, L, P[:, c]) for c in range(3)], 1))
    # strands wound round the bun: spirals on its ellipsoid surface
    a = bun_axis / np.linalg.norm(bun_axis)
    e1 = np.cross(a, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    ring_pts = 16
    for i in range(n_rings):
        th0 = np.arccos(1 - 2 * rng.random() * 0.95)
        ph = rng.random() * 2 * np.pi + np.linspace(0, 1.6 * np.pi, ring_pts)
        th = th0 + np.linspace(0, 0.35, ring_pts) * rng.choice([-1, 1])
        dirs = (np.cos(th)[:, None] * a + np.sin(th)[:, None] * (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2))
        r = 1.0 / np.sqrt(((dirs / brad) ** 2).sum(1))
        P = bc + dirs * (r * (1.02 + 0.04 * rng.random()))[:, None]
        L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        tt = np.linspace(0, L[-1], npts)
        paths.append(np.stack([np.interp(tt, L, P[:, c]) for c in range(3)], 1))
    paths = np.array(paths)
    cv = bpy.data.hair_curves.new('hair_strands')
    cv.add_curves([npts] * len(paths))
    cv.attributes['position'].data.foreach_set('vector', paths.reshape(-1).astype(np.float32))
    rad = cv.attributes.get('radius') or cv.attributes.new('radius', 'FLOAT', 'POINT')
    rr = np.tile(np.linspace(0.00011, 0.00006, npts), len(paths))
    rad.data.foreach_set('value', rr.astype(np.float32))
    ob = bpy.data.objects.new('hair_strands', cv)
    bpy.context.scene.collection.objects.link(ob)
    log('hair strands: %d combed, %d round the bun' % (len(paths) - n_rings, n_rings))
    return ob


# ----------------------------------------------------------------------------------------------- materials
def node_mat(name):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    for n in list(nt.nodes):
        nt.nodes.remove(n)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    return m, nt, out


def principled(nt, **kw):
    p = nt.nodes.new('ShaderNodeBsdfPrincipled')
    for k, v in kw.items():
        p.inputs[k].default_value = v
    return p


def srgb(*c):
    return tuple(((v + 0.055) / 1.055) ** 2.4 if v > 0.04045 else v / 12.92 for v in c) + (1.0,)


SKIN_RAMP = ((0.43, 0.26, 0.175), (0.49, 0.27, 0.19))     # the skin's sRGB colour range (blotch noise between them)


def body_material():
    """Skin (subsurface), the leotard and the slippers as zones of one material, mixed by vertex fields."""
    m, nt, out = node_mat('dancer_body')
    L = nt.links
    tex = nt.nodes.new('ShaderNodeTexCoord')
    # skin: warm medium-brown, subsurface with a red-weighted radius, faint blotch and pore detail
    noise = nt.nodes.new('ShaderNodeTexNoise')
    noise.inputs['Scale'].default_value = 9.0
    noise.inputs['Detail'].default_value = 3.0
    L.new(tex.outputs['Object'], noise.inputs['Vector'])
    ramp = nt.nodes.new('ShaderNodeValToRGB')
    ramp.color_ramp.elements[0].color = srgb(*SKIN_RAMP[0])
    ramp.color_ramp.elements[1].color = srgb(*SKIN_RAMP[1])
    L.new(noise.outputs['Fac'], ramp.inputs['Fac'])
    skin = principled(nt, **{'Roughness': 0.55, 'Subsurface Weight': 1.0, 'Subsurface Scale': 0.0035,
                             'Subsurface Radius': (1.0, 0.30, 0.15), 'Specular IOR Level': 0.38})
    a_lip = nt.nodes.new('ShaderNodeAttribute')
    a_lip.attribute_name = 'lips'
    a_brow = nt.nodes.new('ShaderNodeAttribute')
    a_brow.attribute_name = 'brow'
    mlip = nt.nodes.new('ShaderNodeMixRGB')
    mlip.inputs['Color2'].default_value = srgb(0.40, 0.16, 0.15)
    flip = nt.nodes.new('ShaderNodeMath')
    flip.operation = 'MULTIPLY'
    flip.inputs[1].default_value = 0.8
    L.new(a_lip.outputs['Fac'], flip.inputs[0])
    L.new(flip.outputs['Value'], mlip.inputs['Fac'])
    L.new(ramp.outputs['Color'], mlip.inputs['Color1'])
    mbrow = nt.nodes.new('ShaderNodeMixRGB')
    mbrow.inputs['Color2'].default_value = srgb(0.07, 0.045, 0.035)
    fbrow = nt.nodes.new('ShaderNodeMath')
    fbrow.operation = 'MULTIPLY'
    fbrow.inputs[1].default_value = 0.85
    L.new(a_brow.outputs['Fac'], fbrow.inputs[0])
    L.new(fbrow.outputs['Value'], mbrow.inputs['Fac'])
    L.new(mlip.outputs['Color'], mbrow.inputs['Color1'])
    L.new(mbrow.outputs['Color'], skin.inputs['Base Color'])
    pores = nt.nodes.new('ShaderNodeTexNoise')
    pores.inputs['Scale'].default_value = 900.0
    L.new(tex.outputs['Object'], pores.inputs['Vector'])
    bump = nt.nodes.new('ShaderNodeBump')
    bump.inputs['Strength'].default_value = 0.08
    bump.inputs['Distance'].default_value = 0.0004
    L.new(pores.outputs['Fac'], bump.inputs['Height'])
    L.new(bump.outputs['Normal'], skin.inputs['Normal'])
    # leotard: deep garnet knit with sheen
    leo = principled(nt, **{'Base Color': srgb(0.17, 0.016, 0.042), 'Roughness': 0.7,
                            'Sheen Weight': 0.3, 'Sheen Roughness': 0.45, 'Sheen Tint': srgb(0.6, 0.1, 0.18),
                            'Specular IOR Level': 0.15})
    knit = nt.nodes.new('ShaderNodeTexWave')
    knit.inputs['Scale'].default_value = 900.0
    knit.bands_direction = 'Z'
    L.new(tex.outputs['Object'], knit.inputs['Vector'])
    kb = nt.nodes.new('ShaderNodeBump')
    kb.inputs['Strength'].default_value = 0.05
    L.new(knit.outputs['Fac'], kb.inputs['Height'])
    L.new(kb.outputs['Normal'], leo.inputs['Normal'])
    # slippers: pale pink satin
    sat = principled(nt, **{'Base Color': srgb(0.95, 0.70, 0.70), 'Roughness': 0.28, 'Anisotropic': 0.5,
                            'Sheen Weight': 0.8, 'Sheen Tint': srgb(1.0, 0.85, 0.85), 'Specular IOR Level': 0.7})
    a_leo = nt.nodes.new('ShaderNodeAttribute')
    a_leo.attribute_name = 'leotard'
    a_sl = nt.nodes.new('ShaderNodeAttribute')
    a_sl.attribute_name = 'slipper'
    mix1 = nt.nodes.new('ShaderNodeMixShader')
    L.new(a_leo.outputs['Fac'], mix1.inputs['Fac'])
    L.new(skin.outputs['BSDF'], mix1.inputs[1])
    L.new(leo.outputs['BSDF'], mix1.inputs[2])
    mix2 = nt.nodes.new('ShaderNodeMixShader')
    L.new(a_sl.outputs['Fac'], mix2.inputs['Fac'])
    L.new(mix1.outputs['Shader'], mix2.inputs[1])
    L.new(sat.outputs['BSDF'], mix2.inputs[2])
    L.new(mix2.outputs['Shader'], out.inputs['Surface'])
    if hasattr(m, 'use_sss_translucency'):
        m.use_sss_translucency = True
    return m


def hair_material(comb=False):
    m, nt, out = node_mat('hair_cap' if comb else 'hair_bun')
    L = nt.links
    tex = nt.nodes.new('ShaderNodeTexCoord')
    if comb:   # sin(140 phi + noise): continuous across the azimuth seam because 140 is whole
        at = nt.nodes.new('ShaderNodeAttribute')
        at.attribute_name = 'hair_phi'
        nz = nt.nodes.new('ShaderNodeTexNoise')
        nz.inputs['Scale'].default_value = 40.0
        L.new(tex.outputs['Object'], nz.inputs['Vector'])
        m1 = nt.nodes.new('ShaderNodeMath')
        m1.operation = 'MULTIPLY_ADD'
        m1.inputs[1].default_value = 260.0
        L.new(at.outputs['Fac'], m1.inputs[0])
        L.new(nz.outputs['Fac'], m1.inputs[2])
        wave = nt.nodes.new('ShaderNodeMath')
        wave.operation = 'SINE'
        L.new(m1.outputs['Value'], wave.inputs[0])
        wave_out = wave.outputs['Value']
    else:
        wave = nt.nodes.new('ShaderNodeTexWave')
        wave.wave_type = 'RINGS'
        wave.inputs['Scale'].default_value = 60.0
        wave.inputs['Distortion'].default_value = 6.0
        wave.inputs['Detail'].default_value = 4.0
        L.new(tex.outputs['Object'], wave.inputs['Vector'])
        wave_out = wave.outputs['Fac']
    p = principled(nt, **{'Base Color': srgb(0.05, 0.032, 0.026), 'Roughness': 0.42, 'Specular IOR Level': 0.4,
                          'Sheen Weight': 0.5, 'Sheen Tint': srgb(0.45, 0.32, 0.24)})
    b = nt.nodes.new('ShaderNodeBump')
    b.inputs['Strength'].default_value = 0.18 if comb else 0.35
    b.inputs['Distance'].default_value = 0.0004
    L.new(wave_out, b.inputs['Height'])
    L.new(b.outputs['Normal'], p.inputs['Normal'])
    L.new(p.outputs['BSDF'], out.inputs['Surface'])
    if comb:   # the hairline thins out: alpha from the distance to the cap's edge, broken up by fine noise
        ae = nt.nodes.new('ShaderNodeAttribute')
        ae.attribute_name = 'hair_edge'
        fn = nt.nodes.new('ShaderNodeTexNoise')
        fn.inputs['Scale'].default_value = 900.0
        L.new(tex.outputs['Object'], fn.inputs['Vector'])
        mx = nt.nodes.new('ShaderNodeMath')
        mx.operation = 'MULTIPLY_ADD'
        mx.inputs[1].default_value = 1.6
        mx.inputs[2].default_value = -0.55
        L.new(ae.outputs['Fac'], mx.inputs[0])
        ad = nt.nodes.new('ShaderNodeMath')
        ad.operation = 'ADD'
        ad.use_clamp = True
        L.new(mx.outputs['Value'], ad.inputs[0])
        L.new(fn.outputs['Fac'], ad.inputs[1])
        L.new(ad.outputs['Value'], p.inputs['Alpha'])
        m.blend_method = 'HASHED'
        m.shadow_method = 'HASHED'
    return m


def eye_material():
    m, nt, out = node_mat('eye')
    L = nt.links
    a = nt.nodes.new('ShaderNodeAttribute')
    a.attribute_name = 'iris'
    mp = nt.nodes.new('ShaderNodeMapRange')
    mp.inputs['From Min'].default_value = 0.0
    mp.inputs['From Max'].default_value = 60.0
    L.new(a.outputs['Fac'], mp.inputs['Value'])
    r = nt.nodes.new('ShaderNodeValToRGB')
    el = r.color_ramp.elements
    el[0].position, el[0].color = 0.0, (0.004, 0.003, 0.003, 1)
    el[1].position, el[1].color = 0.17, srgb(0.20, 0.11, 0.06)
    e = el.new(0.43)
    e.color = srgb(0.30, 0.17, 0.09)
    e = el.new(0.47)
    e.color = srgb(0.86, 0.83, 0.80)
    e = el.new(1.0)
    e.color = srgb(0.80, 0.74, 0.72)
    L.new(mp.outputs['Result'], r.inputs['Fac'])
    p = principled(nt, **{'Roughness': 0.05, 'Coat Weight': 1.0, 'Coat Roughness': 0.02, 'Specular IOR Level': 0.6})
    L.new(r.outputs['Color'], p.inputs['Base Color'])
    L.new(p.outputs['BSDF'], out.inputs['Surface'])
    return m


# ----------------------------------------------------------------------------------------------- build entry
def build(args):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    spec = body_spec()
    names = [s['name'] for s in spec['segments']]
    report_dir = Path(args.report_dir)
    body, eyes = import_body()
    lm, notes = measure_landmarks(body)
    check = landmark_report(lm, notes, spec, report_dir)
    mh, mt = mesh_skeleton(lm, spec)
    jh, jt = seg_targets(spec)
    if check['exceeding']:
        prop = dict(schema='terluna.humanities.dance-landmark-proposal/1',
                    status='proposal for the body model; body.json is unchanged and the renderer fits the mesh to it',
                    items={})
        for k in check['exceeding']:
            if k in check['joints']:
                r = check['joints'][k]
                prop['items'][k] = dict(body_json=r['body_json'], mesh_estimate=r['mesh'], diff=r['diff'])
            elif k.startswith('segment:'):
                prop['items'][k] = check['segment_lengths'][k.split(':', 1)[1]]
            elif k.startswith('trunk:'):
                n = k.split(':', 1)[1]
                r = check['trunk'][n]
                xb, xf = r['mesh_back_skin_x'], r['mesh_front_skin_x']
                prop['items'][k] = dict(body_json=r['body_json'], depth_fraction_from_back=r['depth_fraction_from_back'],
                                        proposed_x_at_40pct_depth=round(xb + 0.40 * (xf - xb), 4))
        (report_dir / 'proposal_landmarks.json').write_text(json.dumps(prop, indent=1))
    # weights on the base mesh, at the mesh's own joints
    fit = make_armature('fit_rig', mh, mt)
    auto_weights(body, fit)
    fix_weights(body, names, (mh, mt))
    # detail: apply two multires levels (vertex groups are interpolated by the subdivision)
    body.modifiers['Multires'].levels = MULTIRES_LEVEL
    set_active(body)
    bpy.ops.object.modifier_apply(modifier='Multires')
    log('body vertices after multires', len(body.data.vertices))
    repose(body, eyes, fit, (mh, mt), (jh, jt), names)
    bpy.data.objects.remove(fit)
    # final armature from body.json alone; bones are parentless and each follows its segment's world pose
    rig = make_armature('dancer_rig', jh, jt)
    rig.data.display_type = 'OCTAHEDRAL'
    # costume zones, fabric, face, hair, eyes
    leo, _slip = leotard_and_slipper_masks(body, names)
    yoke_mask(body, names)
    fabric_fill(body, leo)
    face_paint(body, eyes)
    cap, bun, strands = build_hair(body, names)
    for e in eyes:
        eye_iris_attribute(e)
    # skinning: dual-quaternion armature (volume at shoulders and hips) and corrective smooth
    Vr = mesh_coords(body)
    dq = body.vertex_groups.new(name='dq_mix')
    wdq = np.round(smoothstep(0.50, 0.66, Vr[:, 2]) * (1 - smoothstep(1.02, 1.14, Vr[:, 2])), 2)
    for val in np.unique(wdq):
        if val > 0:
            dq.add(np.nonzero(wdq == val)[0].tolist(), float(val), 'REPLACE')
    leo_vg = body.vertex_groups.new(name='leotard_vg')
    lv = np.round(leo, 2)
    for val in np.unique(lv):
        if val > 0:
            leo_vg.add(np.nonzero(lv == val)[0].tolist(), float(val), 'REPLACE')
    disp = body.modifiers.new('leotard_thickness', 'DISPLACE')
    disp.direction = 'NORMAL'
    disp.strength = 0.0018
    disp.mid_level = 0.0
    disp.vertex_group = 'leotard_vg'
    for o in [body] + eyes + [cap, bun]:
        m = o.modifiers.new('skin', 'ARMATURE')
        m.object = rig
        m.use_deform_preserve_volume = False
        o.parent = rig
        if o is cap:   # keep the hair's own modifiers after the armature
            while o.modifiers[0].name != 'skin':
                bpy.ops.object.select_all(action='DESELECT')
                set_active(o)
                bpy.ops.object.modifier_move_to_index(modifier='skin', index=0)
    m2 = body.modifiers.new('skin_dq', 'ARMATURE')
    m2.object = rig
    m2.use_deform_preserve_volume = True
    m2.use_multi_modifier = True
    m2.vertex_group = 'dq_mix'
    cs = body.modifiers.new('corrective', 'CORRECTIVE_SMOOTH')
    cs.smooth_type = 'LENGTH_WEIGHTED'
    cs.iterations = 8
    cs.factor = 0.5
    cs.rest_source = 'ORCO'
    cs.use_pin_boundary = False
    # joints-only mask for the corrective smooth: everything but the head, hands and feet
    W = weights_matrix(body, names)
    col = {n: j for j, n in enumerate(names)}
    keep = W[:, col['head']] + sum(W[:, col[f'{s}_{p}']] for s in 'LR' for p in ('hand', 'toes'))
    g = body.vertex_groups.new(name='corrective_mask')
    vals = np.round(np.clip(1.0 - keep, 0, 1), 2)
    for val in np.unique(vals):
        if val > 0:
            g.add(np.nonzero(vals == val)[0].tolist(), float(val), 'REPLACE')
    cs.vertex_group = 'corrective_mask'
    for p in body.data.polygons:
        p.use_smooth = True
    body.data.materials.clear()
    body.data.materials.append(body_material())
    cap.data.materials.append(hair_material(comb=True))
    bun.data.materials.append(hair_material())
    if strands is not None:   # rigid with the head: bone-parented, at rest where it was built
        sm_, snt, sout = node_mat('hair_strands')
        sp = principled(snt, **{'Base Color': srgb(0.032, 0.021, 0.016), 'Roughness': 0.55,
                                'Specular IOR Level': 0.12})
        snt.links.new(sp.outputs['BSDF'], sout.inputs['Surface'])
        strands.data.materials.append(sm_)
        hb = rig.data.bones['head']
        strands.parent = rig
        strands.parent_type = 'BONE'
        strands.parent_bone = 'head'
        strands.matrix_parent_inverse = (rig.matrix_world @ hb.matrix_local
                                         @ Matrix.Translation((0.0, hb.length, 0.0))).inverted()
    em = eye_material()
    for e in eyes:
        e.data.materials.clear()
        e.data.materials.append(em)
        for p in e.data.polygons:
            p.use_smooth = True
    for o in [body] + eyes + [cap, bun]:
        o.name = o.name if o.name.startswith(('hair', 'GEO')) else o.name
    body.name = 'dancer_body'
    coll = bpy.data.collections.new('dancer')
    bpy.context.scene.collection.children.link(coll)
    for o in [rig, body, cap, bun] + eyes + ([strands] if strands is not None else []):
        for c in list(o.users_collection):
            c.objects.unlink(o)
        coll.objects.link(o)
    for c in list(bpy.data.collections):
        if c.name == 'Body Female - Realistic':
            bpy.data.collections.remove(c)
    meta = dict(source=str(ASSET), collection='Body Female - Realistic', licence='CC0 (Blender Studio)',
                multires_level=MULTIRES_LEVEL, body_json=str(BODY_JSON), producer='visualization/dance/render.py build')
    bpy.context.scene['dance_dancer_meta'] = json.dumps(meta)
    Path(args.blend).parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.blend), compress=True)
    log('saved', args.blend)


# ----------------------------------------------------------------------------------------------- render: motion
def load_motion(path):
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    return d


def sample_motion(mot, fps_out):
    """Segment poses at the render rate: positions interpolated linearly, orientations by normalised lerp."""
    fps = float(mot['fps'])
    P, Q = mot['seg_pos'], mot['seg_quat'].copy()
    for i in range(1, len(Q)):   # sign continuity
        flip = (Q[i] * Q[i - 1]).sum(-1) < 0
        Q[i][flip] *= -1
    dur = (len(P) - 1) / fps
    t = np.arange(0, dur + 1e-9, 1 / fps_out)
    f = t * fps
    i0 = np.clip(np.floor(f).astype(int), 0, len(P) - 1)
    i1 = np.clip(i0 + 1, 0, len(P) - 1)
    a = (f - i0)[:, None, None]
    pos = P[i0] * (1 - a) + P[i1] * a
    q = Q[i0] * (1 - a) + Q[i1] * a
    q /= np.linalg.norm(q, axis=-1, keepdims=True)
    com = mot['com'][i0] * (1 - a[:, :, 0]) + mot['com'][i1] * a[:, :, 0]
    return t, pos, q, com


def key_armature(rig, names, pos, quat, frame0=1):
    """Key each bone so that it follows its segment's world pose: basis = rest^-1 * [R(q) R_rest | p]."""
    act = bpy.data.actions.new('dance_motion')
    rig.animation_data_create()
    rig.animation_data.action = act
    n = len(pos)
    frames = np.arange(frame0, frame0 + n, dtype=float)
    for j, name in enumerate(names):
        pb = rig.pose.bones[name]
        pb.rotation_mode = 'QUATERNION'
        rest = np.array(rig.data.bones[name].matrix_local)
        Rr = rest[:3, :3]
        qr = np.array(Matrix(Rr.tolist()).to_quaternion())
        qri = qr * np.array([1, -1, -1, -1])
        qb = quat_mul(quat_mul(np.broadcast_to(qri, (n, 4)), quat[:, j]), np.broadcast_to(qr, (n, 4)))
        for i in range(1, n):
            if qb[i] @ qb[i - 1] < 0:
                qb[i] = -qb[i]
        loc = (pos[:, j] - rest[:3, 3]) @ Rr          # R_rest^T (p - head)
        for path, arr in ((f'pose.bones["{name}"].rotation_quaternion', qb),
                          (f'pose.bones["{name}"].location', loc)):
            for c in range(arr.shape[1]):
                fc = act.fcurves.new(path, index=c, action_group=name)
                fc.keyframe_points.add(n)
                co = np.empty(2 * n)
                co[0::2] = frames
                co[1::2] = arr[:, c]
                fc.keyframe_points.foreach_set('co', co)
                fc.keyframe_points.foreach_set('interpolation', [1] * n)    # LINEAR
                fc.update()
    return frames


# ----------------------------------------------------------------------------------------------- render: garments
def write_pc2(path, verts, fps):
    verts = np.asarray(verts, np.float32)
    with open(path, 'wb') as f:
        f.write(struct.pack('<12siiffi', b'POINTCACHE2\0', 1, verts.shape[1], 0.0, 1.0, verts.shape[0]))
        f.write(verts.tobytes())


# How each fabric of the garment library looks: sheer silks translucent, satins opaque with a silk satin's soft
# directional sheen, wool, velvet and jersey matte and opaque. Sleeves keep the pale rose of the first renders; skirts take a pale gold so that layers read
# apart. (fabric -> kind)
SHEER_METHOD = {'hashed': 'HASHED', 'blend': 'BLEND'}
SAMPLES = 0                         # render samples; 0 keeps the preset's
FABRIC_LOOK = {'silk_chiffon': 'sheer', 'silk_georgette': 'sheer', 'silk_organza': 'sheer', 'nylon_tulle': 'sheer',
               'silk_habotai': 'silk', 'silk_crepe_de_chine': 'silk', 'silk_taffeta': 'silk',
               'silk_satin_soft': 'satin', 'silk_satin_heavy': 'satin',
               'wool_plain': 'matte', 'cotton_plain': 'matte', 'silk_velvet': 'velvet', 'viscose_jersey': 'matte'}


def garment_fabric(z):
    try:
        meta = json.loads(str(z['meta']))
        fab = meta.get('fabric')
        return fab.get('name') if isinstance(fab, dict) else str(fab)
    except Exception:
        return None


def garment_render(z):
    """A garment file's own render settings, if it carries them (meta 'render': look, sRGB color, alpha)."""
    try:
        return json.loads(str(z['meta'])).get('render') or {}
    except Exception:
        return {}


def add_garment(npz_path, motion_fps, t_render, pc2_dir: Path, index, sheer='hashed'):
    with np.load(npz_path, allow_pickle=False) as z:
        faces = z['faces'].astype(int)
        verts = z['verts'].astype(float)
        gfps = float(z['fps']) if 'fps' in z.files else motion_fps
        fabric = garment_fabric(z)
        own = garment_render(z)
    f = np.clip(t_render * gfps, 0, len(verts) - 1)
    i0 = np.floor(f).astype(int)
    i1 = np.minimum(i0 + 1, len(verts) - 1)
    a = (f - i0)[:, None, None]
    vs = verts[i0] * (1 - a) + verts[i1] * a
    pc2_dir.mkdir(parents=True, exist_ok=True)
    pc2 = pc2_dir / f'{Path(npz_path).stem}_{index}.pc2'
    write_pc2(pc2, vs, RENDER_FPS)
    me = bpy.data.meshes.new(f'garment_{index}')
    me.from_pydata(vs[0].tolist(), [], faces.tolist())
    me.update()
    for p in me.polygons:
        p.use_smooth = True
    ob = bpy.data.objects.new(f'garment_{index}', me)
    bpy.context.scene.collection.objects.link(ob)
    mc = ob.modifiers.new('cloth_cache', 'MESH_CACHE')
    mc.cache_format = 'PC2'
    mc.filepath = str(pc2)
    mc.time_mode = 'FRAME'
    mc.play_mode = 'SCENE'
    mc.frame_start = 1.0
    mc.frame_scale = 1.0
    sm = ob.modifiers.new('smooth', 'SUBSURF')    # one level of smoothing: folds read as cloth, not facets
    sm.levels = sm.render_levels = 1
    look = own.get('look') or FABRIC_LOOK.get(fabric, 'sheer')
    if look != 'sheer':
        # opaque cloth is given its 1.5 mm thickness; a sheer layer stays one surface seen from both sides (backface
        # culling off), since a solidified one is two veils 1.5 mm apart and stacked layers then read as opaque white
        so = ob.modifiers.new('thickness', 'SOLIDIFY')
        so.thickness, so.offset = 0.0015, 0.0
    m, nt, out = node_mat(f'garment_{index}')
    sleeves = 'sleeve' in Path(npz_path).stem
    if look != 'sheer':
        # opaque cloth: a pale gold silk, satin or wool; satin with a silk satin's sheen, silk soft, wool and jersey
        # matte, velvet deep. Satin was glossier until 2 November (roughness 0.22, specular 0.7, anisotropy 0.6), which
        # read as rubber or plastic in the design stills.
        base = {'satin': srgb(0.88, 0.78, 0.58), 'silk': srgb(0.90, 0.80, 0.62), 'matte': srgb(0.80, 0.66, 0.46),
                'velvet': srgb(0.55, 0.16, 0.24)}[look]
        if own.get('color'):
            base = srgb(*own['color'])
        rough = {'satin': 0.38, 'silk': 0.38, 'matte': 0.8, 'velvet': 0.7}[look]
        p = principled(nt, **{'Base Color': base, 'Roughness': rough,
                              'Sheen Weight': {'satin': 0.2, 'silk': 0.5, 'matte': 0.6, 'velvet': 1.0}[look],
                              'Sheen Roughness': 0.4, 'Specular IOR Level': {'satin': 0.45, 'silk': 0.45}.get(look, 0.25),
                              'Anisotropic': 0.3 if look == 'satin' else 0.0})
        nt.links.new(p.outputs['BSDF'], out.inputs['Surface'])
        m.use_backface_culling = False
        ob.data.materials.append(m)
        return ob
    # sheer silk (chiffon, georgette): pale rose, translucent, partly see-through (alpha hashed)
    tint = own.get('color') or (0.93, 0.74, 0.74)
    p = principled(nt, **{'Base Color': srgb(*tint), 'Roughness': 0.5, 'Sheen Weight': 0.9,
                          'Sheen Roughness': 0.35, 'Sheen Tint': srgb(1.0, 0.9, 0.9), 'Specular IOR Level': 0.3,
                          'Alpha': own.get('alpha', 0.72 if sleeves else 0.6)})
    tr = nt.nodes.new('ShaderNodeBsdfTranslucent')
    tr.inputs['Color'].default_value = srgb(*[min(1.0, 1.02 * c) for c in tint]) if own.get('color') \
        else srgb(0.95, 0.70, 0.66)
    mx = nt.nodes.new('ShaderNodeMixShader')
    mx.inputs['Fac'].default_value = 0.4
    tp = nt.nodes.new('ShaderNodeBsdfTransparent')
    mx2 = nt.nodes.new('ShaderNodeMixShader')
    mx2.inputs['Fac'].default_value = 0.28 if sleeves else 0.4
    nt.links.new(p.outputs['BSDF'], mx.inputs[1])
    nt.links.new(tr.outputs['BSDF'], mx.inputs[2])
    nt.links.new(mx.outputs['Shader'], mx2.inputs[1])
    nt.links.new(tp.outputs['BSDF'], mx2.inputs[2])
    nt.links.new(mx2.outputs['Shader'], out.inputs['Surface'])
    # hashed transparency is dithered and smooths out only with many samples; blended transparency is smooth at any
    # sample count but sorts by object, so overlapping layers of one garment can draw out of order
    m.blend_method = SHEER_METHOD[sheer]
    m.show_transparent_back = True
    m.shadow_method = 'HASHED'
    m.use_backface_culling = False
    ob.data.materials.append(m)
    return ob


# ----------------------------------------------------------------------------------------------- render: stage
def add_light(name, kind, loc, target, energy, color, size=None, spot_deg=None, blend=0.3, soft=0.1):
    ld = bpy.data.lights.new(name, kind)
    ld.energy = energy
    ld.color = color
    if kind == 'SPOT':
        ld.spot_size = math.radians(spot_deg)
        ld.spot_blend = blend
        ld.shadow_soft_size = soft
    elif kind == 'AREA':
        ld.shape = 'RECTANGLE'
        ld.size, ld.size_y = size
    else:
        ld.shadow_soft_size = soft
    ld.use_shadow = True
    if hasattr(ld, 'use_contact_shadow'):
        ld.use_contact_shadow = True
        ld.contact_shadow_distance = 0.15
    ob = bpy.data.objects.new(name, ld)
    bpy.context.scene.collection.objects.link(ob)
    ob.location = loc
    d = Vector(target) - Vector(loc)
    ob.rotation_euler = d.to_track_quat('-Z', 'Y').to_euler()
    return ob


# The stage's lighting (--light). stage: side light from booms, amber from stage left and blue from stage right
# (2.2-2.7 kW each), a 4.2 kW rim from behind and a soft key in front several times dimmer than a boom, graded with
# AgX's medium-high contrast; the dancer's skin turns orange or blue-grey with the side it faces and her face dark
# when it turns from the booms to the camera. soft: the booms paler and at 0.8 of their power, a key five times as
# strong, a soft fill from the camera's side, half the rim and AgX's base contrast (a proposal of 3 November).
LIGHT_RIGS = ('stage', 'soft')


def build_stage(cam_az_deg=30.0, centre=(0.0, 0.0), R=7.5, light='stage'):
    sc = bpy.context.scene
    # floor: black Marley, matte with a soft sheen and faint seams between the strips
    bpy.ops.mesh.primitive_plane_add(size=40.0, location=(0, 0, 0))
    floor = bpy.context.active_object
    floor.name = 'marley'
    m, nt, out = node_mat('marley')
    L = nt.links
    tc = nt.nodes.new('ShaderNodeTexCoord')
    n1 = nt.nodes.new('ShaderNodeTexNoise')
    n1.inputs['Scale'].default_value = 0.6
    n1.inputs['Detail'].default_value = 6.0
    L.new(tc.outputs['Object'], n1.inputs['Vector'])
    rr = nt.nodes.new('ShaderNodeMapRange')
    rr.inputs['To Min'].default_value = 0.52
    rr.inputs['To Max'].default_value = 0.68
    L.new(n1.outputs['Fac'], rr.inputs['Value'])
    p = principled(nt, **{'Base Color': srgb(0.028, 0.028, 0.030), 'Specular IOR Level': 0.5})
    L.new(rr.outputs['Result'], p.inputs['Roughness'])
    L.new(p.outputs['BSDF'], out.inputs['Surface'])
    floor.data.materials.append(m)
    # cyclorama: a curved cloth upstage of the camera's view, graded from a warm horizon glow to deep indigo
    az = math.radians(cam_az_deg)
    H = 9.0
    cx, cy = centre
    verts, faces = [], []
    na = 64
    for i in range(na + 1):
        a = math.pi + az + math.radians(-80 + 160 * i / na)
        for z in (0.0, H):
            verts.append((cx + R * math.cos(a), cy + R * math.sin(a), z))
    for i in range(na):
        faces.append([2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1])
    me = bpy.data.meshes.new('cyc')
    me.from_pydata(verts, [], faces)
    cyc = bpy.data.objects.new('cyclorama', me)
    sc.collection.objects.link(cyc)
    m, nt, out = node_mat('cyc')
    L = nt.links
    tc = nt.nodes.new('ShaderNodeTexCoord')
    sep = nt.nodes.new('ShaderNodeSeparateXYZ')
    L.new(tc.outputs['Object'], sep.inputs['Vector'])
    mr = nt.nodes.new('ShaderNodeMapRange')
    mr.inputs['From Min'].default_value = 0.4
    mr.inputs['From Max'].default_value = 6.5
    L.new(sep.outputs['Z'], mr.inputs['Value'])
    ramp = nt.nodes.new('ShaderNodeValToRGB')
    el = ramp.color_ramp.elements
    el[0].position, el[0].color = 0.0, srgb(1.0, 0.60, 0.32)
    el[1].position, el[1].color = 1.0, srgb(0.006, 0.008, 0.03)
    e = el.new(0.07)
    e.color = srgb(0.90, 0.34, 0.36)
    e = el.new(0.20)
    e.color = srgb(0.34, 0.16, 0.50)
    e = el.new(0.42)
    e.color = srgb(0.06, 0.07, 0.24)
    ramp.color_ramp.interpolation = 'B_SPLINE'
    L.new(mr.outputs['Result'], ramp.inputs['Fac'])
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Strength'].default_value = 0.75
    L.new(ramp.outputs['Color'], em.inputs['Color'])
    L.new(em.outputs['Emission'], out.inputs['Surface'])
    cyc.data.materials.append(m)
    # ground row: low black masking in front of the cyc's foot
    gv, gf = [], []
    for i in range(na + 1):
        a = math.pi + az + math.radians(-80 + 160 * i / na)
        for z in (0.0, 0.42):
            gv.append((cx + (R - 0.6) * math.cos(a), cy + (R - 0.6) * math.sin(a), z))
    for i in range(na):
        gf.append([2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1])
    gm = bpy.data.meshes.new('groundrow')
    gm.from_pydata(gv, [], gf)
    gr = bpy.data.objects.new('groundrow', gm)
    sc.collection.objects.link(gr)
    mb, ntb, outb = node_mat('masking_black')
    pb = principled(ntb, **{'Base Color': (0.004, 0.004, 0.004, 1), 'Roughness': 0.9, 'Specular IOR Level': 0.1})
    ntb.links.new(pb.outputs['BSDF'], outb.inputs['Surface'])
    gr.data.materials.append(mb)
    # world: near black, with a light haze for the beams
    w = bpy.data.worlds.new('stage')
    sc.world = w
    w.use_nodes = True
    w.node_tree.nodes['Background'].inputs['Color'].default_value = (0.0006, 0.0007, 0.001, 1)
    vol = w.node_tree.nodes.new('ShaderNodeVolumePrincipled')
    vol.inputs['Density'].default_value = 0.008
    vol.inputs['Anisotropy'].default_value = 0.45
    w.node_tree.links.new(vol.outputs['Volume'], w.node_tree.nodes['World Output'].inputs['Volume'])
    # lights: side light from booms (warm stage left, cool stage right), soft key, rim, top pool. The 'soft' rig
    # (LIGHT_RIGS) evens the skin out: the booms' colours paler and their power 0.8, a key five times as strong, a
    # soft fill from the camera's side and half the rim.
    soft = light == 'soft'
    warm, cool = ((1.0, 0.86, 0.72), (0.78, 0.85, 1.0)) if soft else ((1.0, 0.70, 0.42), (0.50, 0.64, 1.0))
    lights = []
    for side, colr, y in (('SL', warm, 4.2), ('SR', cool, -4.2)):
        for lvl, z, e in (('shin', 0.35, 1400), ('mid', 1.45, 1700), ('head', 2.45, 1500)):
            lights.append(add_light(f'boom_{side}_{lvl}', 'SPOT', (0.0, y, z), (0.0, 0.0, max(z * 0.6, 0.5)),
                                    e * 1.6 * (0.8 if soft else 1.0), colr, spot_deg=36, blend=0.45, soft=0.16))
    lights.append(add_light('key_soft', 'AREA', (5.2 * math.cos(math.radians(30) + 0.25), 5.2 * math.sin(math.radians(30) + 0.25), 4.2),
                            (0, 0, 1.1), 950 if soft else 190, (1.0, 0.93, 0.86), size=(2.4, 2.4)))
    if soft:
        lights.append(add_light('fill', 'AREA', (4.0 * math.cos(math.radians(30)), 4.0 * math.sin(math.radians(30)), 1.5),
                                (0, 0, 1.1), 300, (1.0, 0.96, 0.92), size=(3.0, 2.0)))
    lights.append(add_light('rim', 'SPOT', (-3.6, -1.0, 5.4), (0, 0, 1.3), 2100 if soft else 4200, (0.85, 0.9, 1.0),
                            spot_deg=24, blend=0.5, soft=0.15))
    lights.append(add_light('top_pool', 'SPOT', (0.0, 0.0, 7.5), (0, 0, 0), 210, (1.0, 0.92, 0.82),
                            spot_deg=21, blend=0.8, soft=0.25))
    # the rig was designed for a camera bearing of 30 deg; turn it with the camera and let it travel with the dancer
    rig_ = bpy.data.objects.new('light_rig', None)
    sc.collection.objects.link(rig_)
    rig_.rotation_euler = (0.0, 0.0, math.radians(cam_az_deg - 30.0))
    for ob in lights:
        ob.parent = rig_
    return floor, cyc, lights, rig_


def setup_render(preset, out_path, n_frames, light='stage'):
    sc = bpy.context.scene
    sc.render.hair_type = 'STRAND'
    sc.render.hair_subdiv = 1
    sc.render.engine = 'BLENDER_EEVEE'
    ee = sc.eevee
    final = preset in ('final', 'detail')        # detail: the final quality at 720p
    video = preset == 'video'
    sc.render.resolution_x, sc.render.resolution_y = RESOLUTION[preset]
    sc.render.resolution_percentage = 100
    sc.render.fps = RENDER_FPS
    ee.taa_render_samples = 80 if final else (32 if video else 24)
    if SAMPLES:
        ee.taa_render_samples = SAMPLES
    ee.use_gtao = True
    ee.gtao_distance = 0.35
    ee.gtao_factor = 1.0
    ee.use_bloom = True
    ee.bloom_threshold = 1.0
    ee.bloom_intensity = 0.04
    ee.bloom_radius = 5.5
    ee.use_ssr = True
    ee.ssr_quality = 0.5
    ee.ssr_thickness = 0.05
    ee.use_soft_shadows = not video
    ee.shadow_cube_size = '1024' if final else '512'
    ee.shadow_cascade_size = '1024'
    ee.use_shadow_high_bitdepth = True
    ee.sss_samples = 9
    ee.light_threshold = 0.005
    ee.use_volumetric_lights = True
    ee.use_volumetric_shadows = True
    ee.volumetric_tile_size = '16'      # an 8-pixel tile with 64 samples returned black 1080p frames (4.0.2, integrated GPU)
    ee.volumetric_samples = 32 if final else 16
    ee.volumetric_start = 0.5
    ee.volumetric_end = 25.0
    ee.use_motion_blur = False          # EEVEE's motion blur returned black frames in background renders (4.0.2, integrated GPU)
    ee.motion_blur_shutter = 0.4
    vs = sc.view_settings
    try:
        vs.view_transform = 'AgX'
        vs.look = 'AgX - Base Contrast' if light == 'soft' else 'AgX - Medium High Contrast'
    except TypeError:
        vs.view_transform = 'Filmic'
        vs.look = 'Medium High Contrast'
    vs.exposure = 0.0
    sc.render.film_transparent = False
    sc.frame_start, sc.frame_end = 1, n_frames
    if out_path.suffix.lower() == '.mp4':
        video_output(sc, out_path)


def video_output(sc, out_path):
    """H.264 in an MP4 container through Blender's own FFmpeg, written to exactly out_path."""
    sc.render.image_settings.file_format = 'FFMPEG'
    ff = sc.render.ffmpeg
    ff.format = 'MPEG4'
    ff.codec = 'H264'
    ff.constant_rate_factor = 'HIGH'
    ff.ffmpeg_preset = 'GOOD'
    ff.gopsize = 30
    ff.audio_codec = 'NONE'
    sc.render.use_file_extension = False
    sc.render.filepath = str(out_path)


def gauss_smooth(x, dt, sigma_s):
    k = np.arange(-int(3 * sigma_s / dt), int(3 * sigma_s / dt) + 1) * dt
    g = np.exp(-0.5 * (k / sigma_s) ** 2)
    g /= g.sum()
    pad = len(k) // 2
    xp = np.pad(x, ((pad, pad),) + ((0, 0),) * (x.ndim - 1), mode='edge')
    if x.ndim == 1:
        return np.convolve(xp, g, mode='valid')
    return np.stack([np.convolve(xp[:, c], g, mode='valid') for c in range(x.shape[1])], 1)


def camera_plan(com, t, pelvis_quat, az_arg, dist_arg):
    """Bearing and distance: a fixed 30 deg bearing for phrases in place; for travelling phrases the camera stands
    square to the travel (the dancer's line in profile), on the side the dancer's front turns toward, and further back."""
    dt = t[1] - t[0] if len(t) > 1 else 1 / RENDER_FPS
    sm = gauss_smooth(com[:, :2], dt, 0.8)
    travel = sm[-1] - sm[0]
    if az_arg is not None:
        az = az_arg
    elif np.linalg.norm(travel) < 1.0:
        az = 30.0
    else:
        d = travel / np.linalg.norm(travel)
        perp = np.array([-d[1], d[0]])
        fwd = quat_to_mat(pelvis_quat[len(pelvis_quat) // 2])[:2, 0]
        if perp @ fwd < 0:
            perp = -perp
        az = math.degrees(math.atan2(perp[1], perp[0]))
    dist = dist_arg if dist_arg else (5.5 if np.linalg.norm(travel) < 1.0 else 6.8)
    return az, dist, sm


def follow_camera(com, t, az_deg=30.0, dist=5.5, height=1.1, aim_dz=0.06, sigma_s=0.8, rig=None, fstop=2.8):
    """A camera at a fixed bearing that follows the smoothed centre of mass: across the floor (sigma 0.8 s) and, more
    closely, in height (sigma 0.35 s), so a leap stays framed. It aims `aim_dz` above the centre of mass and rises
    with its aim above the default 6 cm, so a close still of the bodice is taken level."""
    sc = bpy.context.scene
    cd = bpy.data.cameras.new('camera')
    cd.lens = 50.0
    cd.sensor_width = 36.0
    cd.dof.use_dof = True
    cd.dof.aperture_fstop = fstop
    height += aim_dz - 0.06
    cam = bpy.data.objects.new('camera', cd)
    sc.collection.objects.link(cam)
    sc.camera = cam
    tgt = bpy.data.objects.new('camera_target', None)
    sc.collection.objects.link(tgt)
    dt = t[1] - t[0] if len(t) > 1 else 1 / RENDER_FPS
    sm = gauss_smooth(com[:, :2], dt, sigma_s)
    zs = gauss_smooth(com[:, 2], dt, 0.35)
    zref = float(np.median(com[:, 2]))
    az = math.radians(az_deg)
    off = np.array([math.cos(az), math.sin(az)]) * dist
    ct = cam.constraints.new('TRACK_TO')
    ct.target = tgt
    ct.track_axis = 'TRACK_NEGATIVE_Z'
    ct.up_axis = 'UP_Y'
    cd.dof.focus_object = tgt
    for i, p in enumerate(sm):
        f = i + 1
        cam.location = (p[0] + off[0], p[1] + off[1], height + 0.8 * (zs[i] - zref))
        cam.keyframe_insert('location', frame=f)
        tgt.location = (p[0], p[1], zs[i] + aim_dz)
        tgt.keyframe_insert('location', frame=f)
        if rig is not None:
            rig.location = (p[0], p[1], 0.0)
            rig.keyframe_insert('location', frame=f)
    return cam, sm


# The bodice's looks (--bodice-look). matte: a matte silk to sit with sheer chiffon skirts, its light spread wide
# (roughness 0.5) with a soft sheen tinted by its own colour, low specular and no anisotropy, the leotard's knit ribbing
# all but gone. satin: the glossy satin of the first design stills (until 2 November the only look), which read as a
# rubber catsuit beside chiffon; its sheen keeps the leotard's garnet tint.
BODICE_LOOK = {'matte': dict(roughness=0.5, sheen=0.4, sheen_roughness=0.5, specular=0.3, anisotropic=0.0, bump=0.004),
               'satin': dict(roughness=0.3, sheen=0.25, sheen_roughness=None, specular=0.6, anisotropic=0.4, bump=0.05)}


def set_bodice(body, rgb, look='matte'):
    """Dress the leotard zone as a bodice of the given sRGB colour in one of BODICE_LOOK's looks. The leotard's node is
    the principled BSDF built with the deep garnet base colour."""
    nt = body.data.materials[0].node_tree
    garnet = srgb(0.17, 0.016, 0.042)
    s = BODICE_LOOK[look]
    for n in nt.nodes:
        if n.type == 'BSDF_PRINCIPLED' and all(abs(a - b) < 1e-3 for a, b in zip(n.inputs['Base Color'].default_value,
                                                                                  garnet)):
            n.inputs['Base Color'].default_value = srgb(*rgb)
            n.inputs['Roughness'].default_value = s['roughness']
            n.inputs['Sheen Weight'].default_value = s['sheen']
            n.inputs['Specular IOR Level'].default_value = s['specular']
            n.inputs['Anisotropic'].default_value = s['anisotropic']
            if s['sheen_roughness'] is not None:
                n.inputs['Sheen Roughness'].default_value = s['sheen_roughness']
                n.inputs['Sheen Tint'].default_value = srgb(*[0.5 + 0.5 * c for c in rgb])
            for link in n.inputs['Normal'].links:
                if link.from_node.type == 'BUMP':
                    link.from_node.inputs['Strength'].default_value = s['bump']
            log('bodice: %s, sRGB %s' % (look, ','.join('%.2f' % c for c in rgb)))
            return True
    log('bodice: leotard node not found')
    return False


def set_yoke(body, rgba):
    """Lay a sheer chiffon yoke over the skin in the build's 'yoke' zone (yoke_mask), at render time.

    rgba: the chiffon's sRGB colour and its strength a (0 to 1), the share of the skin it veils seen face on. Toward
    the silhouette the light crosses the fabric obliquely, so the veil thickens as a sheet's would: transmission
    (1 - a) ** (1 / cos theta), with cos theta from the Layer Weight node's facing output, floored at 0.3; the hem band
    ('yoke_hem') is two layers, transmission squared. The veil mixes a chiffon BSDF over the skin BSDF, so the skin
    keeps its subsurface colour beneath, lightened toward the fabric. The leotard is mixed in above the yoke, so its
    edge is the yoke's edge at the bodice. The material built by `build` is unchanged; the dancer must carry the
    'yoke' attributes (dancers built before 2 November do not: rebuild them).

    Where threads cover the skin they reflect part of the light in their own colour and pass part on to the skin and
    back, so the chiffon BSDF's colour is 0.35 x the fabric's plus 0.25 x the skin's mean albedo (linear; a single pass
    through translucent silk), with roughness 0.5, a soft sheen tinted toward white and low specular. Ivory over this
    skin reads as a pale grey veil at a = 0.5, and a = 0.3 shows the skin clearly.

    Tried and set aside: with the fabric's full colour and a cos theta floor of 0.2, an ivory yoke at a = 0.5 read as
    an opaque white top. Linking the skin's colour nodes into the chiffon BSDF made Blender abort with bad_alloc under
    the 14 GB memory cap, so the skin's albedo enters as a constant. A weave bump (threads 0.5 mm apart, below a pixel
    at every framing used) made 80-sample stills come back black with the GPU free while 24-sample stills rendered:
    summed from three band textures it blackened full-figure stills, and as one band texture, built as the leotard's
    knit is, it blackened close views at 2.4 m. Without it the full-figure 80-sample still rendered. The cause is not
    known, and the yoke has no weave bump."""
    me = body.data
    if 'yoke' not in me.attributes or 'yoke_hem' not in me.attributes:
        log('yoke: this dancer has no yoke zone; rebuild it (--rebuild)')
        return False
    nt = me.materials[0].node_tree
    L = nt.links
    mix_leo = next((n for n in nt.nodes if n.type == 'MIX_SHADER' and n.inputs['Fac'].links
                    and getattr(n.inputs['Fac'].links[0].from_node, 'attribute_name', '') == 'leotard'), None)
    if mix_leo is None or not mix_leo.inputs[1].links:
        log('yoke: the leotard mix was not found')
        return False
    skin_out = mix_leo.inputs[1].links[0].from_socket
    r, g, b, a = rgba
    a = min(max(a, 0.0), 0.999)

    def math(op, x0, x1=None, clamp=False):
        m = nt.nodes.new('ShaderNodeMath')
        m.operation = op
        m.use_clamp = clamp
        for k, v in enumerate((x0, x1)):
            if v is None:
                continue
            if isinstance(v, (int, float)):
                m.inputs[k].default_value = v
            else:
                L.new(v, m.inputs[k])
        return m.outputs['Value']

    thread = srgb(r, g, b)
    chiffon = principled(nt, **{'Roughness': 0.5, 'Sheen Weight': 0.5,
                                'Sheen Roughness': 0.35, 'Sheen Tint': srgb(*[0.5 + 0.5 * c for c in (r, g, b)]),
                                'Specular IOR Level': 0.3})
    skin = [0.5 * (u + v) for u, v in zip(srgb(*SKIN_RAMP[0])[:3], srgb(*SKIN_RAMP[1])[:3])]
    chiffon.inputs['Base Color'].default_value = tuple(0.35 * f + 0.25 * k for f, k in zip(thread[:3], skin)) + (1.0,)
    # veil = yoke * (1 - (1 - a) ** ((1 + hem) / max(cos theta, 0.3)))
    lw = nt.nodes.new('ShaderNodeLayerWeight')
    lw.inputs['Blend'].default_value = 0.5               # facing = 1 - |cos theta|
    cos_t = math('MAXIMUM', math('SUBTRACT', 1.0, lw.outputs['Facing']), 0.3)
    a_y = nt.nodes.new('ShaderNodeAttribute')
    a_y.attribute_name = 'yoke'
    a_h = nt.nodes.new('ShaderNodeAttribute')
    a_h.attribute_name = 'yoke_hem'
    path = math('DIVIDE', math('ADD', a_h.outputs['Fac'], 1.0), cos_t)
    veil = math('MULTIPLY', math('SUBTRACT', 1.0, math('POWER', 1.0 - a, path)), a_y.outputs['Fac'], clamp=True)
    mix_y = nt.nodes.new('ShaderNodeMixShader')
    mix_y.name = 'yoke_mix'
    L.new(veil, mix_y.inputs['Fac'])
    L.new(skin_out, mix_y.inputs[1])
    L.new(chiffon.outputs['BSDF'], mix_y.inputs[2])
    L.new(mix_y.outputs['Shader'], mix_leo.inputs[1])
    log('yoke: chiffon sRGB %.2f,%.2f,%.2f, strength %.2f' % (r, g, b, a))
    return True


# Dressed as a colour on the leotard zone, a bodice follows every contour of the skin and reads as moulded on the body.
# bodice_shell makes it a garment: a shell over the torso that spans the body's hollows, lifted 3 mm off it, given
# cloth's thickness and skinned as the body is, with a neckline cut by the shader, a surface style and a piped edge in
# the trim colour.
BODICE_STYLES = ('structured', 'ruched', 'pleated')
BODICE_NECKLINES = ('sweetheart', 'vee', 'straight')
BODICE_BOTTOM = 1.0                    # m, rest height of the lower edge: under the sash, 2 cm below the natural waist
BUST_APEX = (0.158, 0.076, 1.183)      # forward, to the side and up (m): the rest mesh's bust apex (measured 3 November)


def neckline_z(x, ay, kind):
    """Rest height (m) of the bodice's top edge at a point forward x and |y| to the side: the front's shape by kind,
    straight across the back below the shoulder blades (1.205 m), the two blended across the sides as the leotard's
    neckline is. The sweetheart dips between the breasts (1.212 m) and arches over each (1.252 m at the apex's |y|);
    the vee runs from 1.165 m at the centre front; the straight line is level at 1.235 m, 5 cm above the apex."""
    if kind == 'sweetheart':
        zf = 1.212 + 0.040 * np.sin(np.pi * np.clip(ay / (2 * BUST_APEX[1]), 0.0, 1.0))
    elif kind == 'vee':
        zf = np.minimum(1.165 + 0.95 * ay, 1.27)
    else:
        zf = np.full(np.shape(ay), 1.235)
    front = smoothstep(-0.03, 0.03, x)
    return 1.205 + (zf - 1.205) * front


def bodice_shell(body, names, style='structured', neckline='sweetheart', rgb=(0.93, 0.90, 0.84),
                 trim=(0.70, 0.52, 0.25)):
    """Build the bodice shell from the body (rest mesh, skinned by the body's own armature modifiers) and clear the
    body's leotard zone above its neckline, so that the skin or the sheer yoke shows there."""
    import bmesh
    if style not in BODICE_STYLES or neckline not in BODICE_NECKLINES:
        raise SystemExit(f'bodice shell: style in {BODICE_STYLES}, neckline in {BODICE_NECKLINES}')
    V = mesh_coords(body)
    W = weights_matrix(body, names)
    col = {n: j for j, n in enumerate(names)}
    arm = sum(W[:, col[f'{s}_{q}']] for s in 'LR' for q in ('upper_arm', 'forearm', 'hand'))
    headneck = W[:, col['head']] + W[:, col['neck']]
    x, y, z = V.T
    nz = neckline_z(x, np.abs(y), neckline)
    keep = (arm < 0.3) & (headneck < 0.3) & (z > BODICE_BOTTOM) & (z < nz + 0.025)
    # clear the body's own bodice above the new neckline (2 mm below it, under the shell's edge)
    leo = np.empty(len(V))
    body.data.attributes['leotard'].data.foreach_get('value', leo)
    leo *= 1.0 - smoothstep(-0.006, -0.002, z - nz)
    body.data.attributes['leotard'].data.foreach_set('value', leo.astype(np.float32))
    shell = body.copy()
    shell.data = body.data.copy()
    shell.name = shell.data.name = 'bodice_shell'
    bpy.context.scene.collection.objects.link(shell)
    bm = bmesh.new()
    bm.from_mesh(shell.data)
    bm.verts.ensure_lookup_table()
    bmesh.ops.delete(bm, geom=[bm.verts[i] for i in np.where(~keep)[0]], context='VERTS')
    bm.verts.ensure_lookup_table()
    boundary = np.array([v.is_boundary for v in bm.verts])
    bm.to_mesh(shell.data)
    bm.free()
    shell.data.update()
    V = mesh_coords(shell)
    E = mesh_edges(shell)
    N = np.empty(len(V) * 3)
    shell.data.vertices.foreach_get('normal', N)
    N = N.reshape(-1, 3)
    rest = V.copy()
    # as boning and a lined shell span the abdomen's relief, the cleavage and the spine's groove: the torso with its
    # relief a few centimetres across taken out by Taubin's smoothing (steps of 0.5 and -0.53, which keep the torso's
    # and the bust's size), the shell at most 1.5 mm inside it and up to 15 mm out, then filled outward only until no
    # hollow is left, and lifted 3 mm
    S = V.copy()
    for _ in range(60):
        S = smooth_field(smooth_field(S, E, 1, mask=~boundary, lam=0.5), E, 1, mask=~boundary, lam=-0.53)
    dn = np.clip(((S - V) * N).sum(1), -0.0015, 0.015)
    for _ in range(300):
        X = V + N * dn[:, None]
        rise = ((smooth_field(X, E, 1, mask=~boundary) - X) * N).sum(1)
        dn = np.minimum(dn + np.maximum(rise, 0.0), 0.015)
    D = smooth_field(N * (dn + 0.003)[:, None], E, 8)
    shell.data.vertices.foreach_set('co', (V + D).ravel())
    shell.data.update()
    # the copy brings the body's custom split normals, which shade the shell as the body itself, its abdominal muscles
    # showing through whatever the shell's shape (the first shells, 3 November); the shell shades by its own surface
    with bpy.context.temp_override(object=shell, active_object=shell, selected_objects=[shell],
                                   selected_editable_objects=[shell]):
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
    for nm, typ, val in (('rest_co', 'FLOAT_VECTOR', rest),
                         ('below_neck', 'FLOAT', neckline_z(rest[:, 0], np.abs(rest[:, 1]), neckline) - rest[:, 2])):
        a = shell.data.attributes.get(nm) or shell.data.attributes.new(nm, typ, 'POINT')
        a.data.foreach_set('vector' if typ == 'FLOAT_VECTOR' else 'value', val.astype(np.float32).ravel())
    for m in list(shell.modifiers):
        if m.type == 'DISPLACE':
            shell.modifiers.remove(m)
    so = shell.modifiers.new('thickness', 'SOLIDIFY')
    so.thickness, so.offset = 0.0015, 1.0
    sm = shell.modifiers.new('smooth', 'SUBSURF')
    sm.levels = sm.render_levels = 1
    shell.data.materials.clear()
    shell.data.materials.append(bodice_material(style, rgb, trim))
    log('bodice shell: %s, %s neckline, %d vertices; lift %.1f-%.1f mm' % (
        style, neckline, len(V), 1e3 * np.linalg.norm(D, axis=1).min(), 1e3 * np.linalg.norm(D, axis=1).max()))
    return shell


def bodice_material(style, rgb, trim):
    """Matte silk in rgb with an 8 mm piped edge in satin of the trim colour along the neckline (from the shell's
    below_neck field), cut at the neckline by alpha clipping, and the style's surface as a bump from the rest
    position: princess seams (structured), gathers running down to the centre front (ruched), or pleats radiating from
    below the waist (pleated)."""
    m, nt, out = node_mat(f'bodice_{style}')
    L = nt.links.new
    rest = nt.nodes.new('ShaderNodeAttribute')
    rest.attribute_name = 'rest_co'
    sep = nt.nodes.new('ShaderNodeSeparateXYZ')
    L(rest.outputs['Vector'], sep.inputs[0])
    X, Y, Z = sep.outputs[0], sep.outputs[1], sep.outputs[2]

    def op(kind, a, b=None, c=None):
        n = nt.nodes.new('ShaderNodeMath')
        n.operation = kind
        for k, v in enumerate((a, b, c)):
            if v is None:
                continue
            if isinstance(v, (int, float)):
                n.inputs[k].default_value = float(v)
            else:
                L(v, n.inputs[k])
        return n.outputs[0]

    def ramp(v, a, b, lo, hi):                     # smoothstep from (a -> lo) to (b -> hi)
        n = nt.nodes.new('ShaderNodeMapRange')
        n.interpolation_type = 'SMOOTHSTEP'
        L(v, n.inputs['Value'])
        n.inputs['From Min'].default_value, n.inputs['From Max'].default_value = a, b
        n.inputs['To Min'].default_value, n.inputs['To Max'].default_value = lo, hi
        return n.outputs['Result']

    AY = op('ABSOLUTE', Y)
    if style == 'structured':                      # princess seams over the bust apex, and down the back
        ys_front = op('MULTIPLY_ADD', op('SUBTRACT', Z, BUST_APEX[2]), 0.12, BUST_APEX[1])
        front = op('GREATER_THAN', X, 0.0)
        ys = op('ADD', op('MULTIPLY', front, op('SUBTRACT', ys_front, 0.065)), 0.065)
        d = op('ABSOLUTE', op('SUBTRACT', AY, ys))
        height, strength, dist = ramp(d, 0.0015, 0.005, 0.0, 1.0), 0.9, 0.004
    elif style == 'ruched':                        # gathers rising from the centre front's waist to the sides
        # irregular as gathered cloth lies: two noise fields shift the folds (about 17 mm apart) by up to a third of a
        # fold and vary their depth, and the crests are broad; a regular 11 mm sine read as stripes (3 November)
        def noise(scale):
            n = nt.nodes.new('ShaderNodeTexNoise')
            n.inputs['Scale'].default_value = scale
            n.inputs['Detail'].default_value = 2.0
            L(rest.outputs['Vector'], n.inputs['Vector'])
            return n.outputs['Fac']
        line = op('MULTIPLY_ADD', AY, 0.55, op('SUBTRACT', Z, BODICE_BOTTOM))
        phase = op('MULTIPLY_ADD', op('SUBTRACT', noise(14.0), 0.5), 4.0, op('MULTIPLY', line, 2 * math.pi / 0.017))
        crest = op('POWER', op('MULTIPLY_ADD', op('SINE', phase), 0.5, 0.5), 0.6)
        height, strength, dist = op('MULTIPLY', crest, op('MULTIPLY_ADD', noise(9.0), 1.2, 0.4)), 0.4, 0.003
    else:                                          # pleats fanning from 7 cm below the waist
        theta = op('ARCTAN2', Y, op('SUBTRACT', Z, 0.93))
        height, strength, dist = op('FRACT', op('MULTIPLY', theta, 1.0 / 0.055)), 0.5, 0.0025
    bump = nt.nodes.new('ShaderNodeBump')
    bump.inputs['Strength'].default_value = strength
    bump.inputs['Distance'].default_value = dist
    L(height, bump.inputs['Height'])
    neck = nt.nodes.new('ShaderNodeAttribute')
    neck.attribute_name = 'below_neck'
    below = neck.outputs['Fac']
    alpha = op('GREATER_THAN', below, 0.0)
    piping = ramp(below, 0.007, 0.009, 1.0, 0.0)
    silk = principled(nt, **{'Base Color': srgb(*rgb), 'Roughness': 0.5, 'Sheen Weight': 0.4, 'Sheen Roughness': 0.5,
                             'Specular IOR Level': 0.3})
    satin = principled(nt, **{'Base Color': srgb(*trim), 'Roughness': 0.35, 'Specular IOR Level': 0.5,
                              'Anisotropic': 0.3})
    for b in (silk, satin):
        L(bump.outputs['Normal'], b.inputs['Normal'])
        L(alpha, b.inputs['Alpha'])
    mix = nt.nodes.new('ShaderNodeMixShader')
    L(piping, mix.inputs['Fac'])
    L(silk.outputs['BSDF'], mix.inputs[1])
    L(satin.outputs['BSDF'], mix.inputs[2])
    L(mix.outputs['Shader'], out.inputs['Surface'])
    m.blend_method = 'CLIP'
    m.alpha_threshold = 0.5
    m.shadow_method = 'CLIP'
    m.use_backface_culling = False
    return m


def render(args):
    bpy.ops.wm.open_mainfile(filepath=str(args.blend))
    spec = body_spec()
    names = [s['name'] for s in spec['segments']]
    mot = load_motion(args.motion)
    assert [str(s) for s in mot['segments']] == names, 'motion segments differ from body.json'
    t, pos, quat, com = sample_motion(mot, RENDER_FPS)
    rig = bpy.data.objects['dancer_rig']
    key_armature(rig, names, pos, quat)
    if args.bodice:
        set_bodice(bpy.data.objects['dancer_body'], [float(v) for v in args.bodice.split(',')], args.bodice_look)
    if args.yoke:
        rgba = [float(v) for v in args.yoke.split(',')]
        if len(rgba) != 4:
            raise SystemExit('--yoke takes r,g,b,a')
        set_yoke(bpy.data.objects['dancer_body'], rgba)
    if args.bodice_shell:
        style, neckline = (args.bodice_shell.split(',') + ['sweetheart'])[:2]
        rgb = [float(v) for v in args.bodice.split(',')] if args.bodice else [0.93, 0.90, 0.84]
        bodice_shell(bpy.data.objects['dancer_body'], names, style, neckline, rgb,
                     [float(v) for v in args.trim.split(',')])
    az, dist, path = camera_plan(com, t, quat[:, names.index('pelvis')], args.cam_az, args.cam_dist)
    centre = 0.5 * (path.min(0) + path.max(0))
    radius = 7.5 + 0.5 * float(np.linalg.norm(path.max(0) - path.min(0)))
    _f, _c, _l, light_rig = build_stage(az, tuple(centre), radius, args.light)
    follow_camera(com, t, az_deg=az, dist=dist, rig=light_rig, aim_dz=args.cam_aim, fstop=args.cam_fstop)
    log('camera bearing %.1f deg, distance %.1f m; cyclorama centre (%.2f, %.2f), radius %.1f m'
        % (az, dist, centre[0], centre[1], radius))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for i, g in enumerate(args.garment or []):
        add_garment(g, float(mot['fps']), t, out.parent / 'pc2', i, args.sheer)
    setup_render(args.stills_preset if args.stills else args.preset, out, len(t), args.light)
    if args.save_scene:
        bpy.ops.wm.save_as_mainfile(filepath=str(args.save_scene), compress=True)
    sc = bpy.context.scene
    if args.stills:
        stem = out.with_suffix('')
        img = sc.render.image_settings
        fmt = img.file_format
        for ts in [float(s) for s in args.stills.split(',') if s]:
            f = int(round(ts * RENDER_FPS)) + 1
            sc.frame_set(f)
            img.file_format = 'PNG'
            sc.render.filepath = f'{stem}_t{ts:05.2f}.png'
            bpy.ops.render.render(write_still=True)
            log('still', sc.render.filepath)
        img.file_format = fmt
    if out.suffix.lower() == '.mp4' and not args.stills_only:
        setup_render(args.preset, out, len(t))      # the video's own size and output
        if args.frames:
            a, b = (int(x) for x in args.frames.split(':'))
            sc.frame_start, sc.frame_end = a, b
        bpy.ops.render.render(animation=True)
        log('video', out)


def main():
    argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['build', 'render'])
    ap.add_argument('--blend', required=True)
    ap.add_argument('--report-dir', default='.')
    ap.add_argument('--motion')
    ap.add_argument('--out')
    ap.add_argument('--garment', action='append')
    ap.add_argument('--preset', default='draft', choices=list(RESOLUTION))
    ap.add_argument('--stills-preset', default='final', choices=list(RESOLUTION))
    ap.add_argument('--bodice', default='', help='sRGB r,g,b: dress the leotard as a bodice of this colour')
    ap.add_argument('--bodice-look', default='matte', choices=list(BODICE_LOOK),
                    help="the bodice's look: matte silk (default) or glossy satin")
    ap.add_argument('--yoke', default='',
                    help='sRGB r,g,b and strength a (0-1): a sheer chiffon yoke over the skin above the bodice')
    ap.add_argument('--sheer', default='hashed', choices=list(SHEER_METHOD))
    ap.add_argument('--bodice-shell', default='')
    ap.add_argument('--trim', default='0.70,0.52,0.25')
    ap.add_argument('--samples', type=int, default=0)
    ap.add_argument('--stills', default='')
    ap.add_argument('--stills-only', action='store_true')
    ap.add_argument('--frames', default='')
    ap.add_argument('--cam-az', type=float, default=None)
    ap.add_argument('--cam-dist', type=float, default=None)
    ap.add_argument('--cam-aim', type=float, default=0.06)
    ap.add_argument('--light', default='stage', choices=LIGHT_RIGS)
    ap.add_argument('--cam-fstop', type=float, default=2.8)
    ap.add_argument('--save-scene', default='')
    args = ap.parse_args(argv)
    global SAMPLES
    SAMPLES = args.samples
    if args.cmd == 'build':
        build(args)
    else:
        render(args)


if __name__ == '__main__':
    main()
