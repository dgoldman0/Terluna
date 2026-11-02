"""Measure the rigged dancer's skin against body.json's capsules and write the dress form that the cloth of dressed
renders collides with (Blender side, run on the built dancer):

    blender -b DANCER.blend --python visualization/dance/dress_form.py -- --out visualization/dance/dress_form.json

body.json's capsules are the dynamics' body. The rendered skin lies outside them by up to 5 cm at the thighs and
4 cm at the waist, so cloth simulated on them passes through the skin the renderer shows. The form splits each
capsule into three tapered pieces. A limb's pieces reach the 95th percentile of the skin's distance from the axis
near each node, never less than the capsule's radius; the hands keep theirs, since a water sleeve's cuff closes
round the wrist and its tube holds the hand. The trunk's upright capsules (pelvis, abdomen, thorax) become
left and right pairs whose union spans the skin's cross-section at each node (its depth as the radius, the rest of
its width as the pair's spread), with the leotard's thickness added; the trunk's crosswise capsules stay as they are.
The waist's outline at the natural waist (2 cm above the omphalion, as the garments place the waistband) is
written as a polar radius about the abdomen's axis for the skirts' waistband.
"""
import hashlib
import json
import sys
from pathlib import Path

import bpy
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BODY_JSON = ROOT / 'humanities' / 'dance' / 'body' / 'body.json'
SCHEMA = 'terluna.visualization.dance-dress-form/1'
LEOTARD_M = 0.0018                 # the renderer's leotard thickness (render.py build)
WAIST_ABOVE_OMPHALION = 0.02       # as humanities/dance/garment/garments.py
TRUNK = ('pelvis', 'abdomen', 'thorax')
NODES = 4                          # three pieces per capsule
LIMB_PERCENTILE = 95


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seg_dist(P, a, b):
    ab = b - a
    t = np.clip((P - a) @ ab / max(float(ab @ ab), 1e-12), 0.0, 1.0)
    return np.linalg.norm(P - (a + t[:, None] * ab), axis=1), t


def main():
    argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
    out = Path(argv[argv.index('--out') + 1])
    spec = json.loads(BODY_JSON.read_text())
    names = [s['name'] for s in spec['segments']]
    body = bpy.data.objects['dancer_body']
    me = body.data
    V = np.empty(len(me.vertices) * 3)
    me.vertices.foreach_get('co', V)
    V = V.reshape(-1, 3)                       # rest pose, body.json's world frame
    gi = {g.index: names.index(g.name) for g in body.vertex_groups if g.name in names}
    W = np.zeros((len(V), len(names)))
    for v in me.vertices:
        for g in v.groups:
            j = gi.get(g.group)
            if j is not None:
                W[v.index, j] = g.weight
    dom = W.argmax(1)
    trunk_v = np.isin(dom, [names.index(n) for n in TRUNK])
    pieces, report = [], {}
    for i, s in enumerate(spec['segments']):
        o = np.array(s['origin'], float)
        caps = [(np.array(c['a'], float), np.array(c['b'], float), float(c['radius'])) for c in s['capsules']]
        own = dom == i
        if caps:
            G = np.stack([seg_dist(V[own], o + a, o + b)[0] - r for a, b, r in caps], 1)
            nearest = np.zeros(len(V), int) - 1
            nearest[own] = G.argmin(1)
        for k, (a, b, r) in enumerate(caps):
            A, B = o + a, o + b
            upright = s['name'] in TRUNK and abs((B - A)[2]) > 0.8 * np.linalg.norm(B - A)
            ts = np.linspace(0.0, 1.0, NODES)
            if upright:
                # left and right pair spanning the trunk's cross-section at each node's height
                nodes = []
                for t in ts:
                    z = A[2] + t * (B[2] - A[2])
                    m = trunk_v & (np.abs(V[:, 2] - z) < 0.01)
                    if m.sum() < 30:
                        nodes.append((A[0], 0.0, r))
                        continue
                    x, y = V[m, 0], V[m, 1]
                    front, back = np.percentile(x, 98), np.percentile(x, 2)
                    half_w = np.percentile(np.abs(y), 98)
                    rr = max(0.5 * (front - back) + LEOTARD_M, 0.02)
                    nodes.append((0.5 * (front + back), max(half_w + LEOTARD_M - rr, 0.0), rr))
                for j in range(NODES - 1):
                    for sgn in (1.0, -1.0):
                        pa = np.array([nodes[j][0], sgn * nodes[j][1], A[2] + ts[j] * (B[2] - A[2])])
                        pb = np.array([nodes[j + 1][0], sgn * nodes[j + 1][1], A[2] + ts[j + 1] * (B[2] - A[2])])
                        pieces.append(dict(segment=s['name'], a=(pa - o).round(4).tolist(), b=(pb - o).round(4).tolist(),
                                           ra=round(nodes[j][2], 4), rb=round(nodes[j + 1][2], 4)))
                report[f"{s['name']}:{k}"] = dict(kind='pair', capsule_radius=r,
                                                   nodes=[dict(x=round(n[0], 4), spread=round(n[1], 4),
                                                               radius=round(n[2], 4)) for n in nodes])
                continue
            m = nearest == k
            d, t = seg_dist(V, A, B)
            radii = []
            for tj in ts:
                sel = m & (np.abs(t - tj) <= 0.5 / (NODES - 1))
                gap = np.percentile(d[sel], LIMB_PERCENTILE) if sel.sum() >= 20 else r
                radii.append(max(r, float(gap)))
            if s['name'] in TRUNK or s['name'].endswith('hand'):
                radii = [r] * NODES        # crosswise trunk capsules, and the hands inside the water sleeves' cuffs
            for j in range(NODES - 1):
                pa, pb = a + ts[j] * (b - a), a + ts[j + 1] * (b - a)
                pieces.append(dict(segment=s['name'], a=pa.round(4).tolist(), b=pb.round(4).tolist(),
                                   ra=round(radii[j], 4), rb=round(radii[j + 1], 4)))
            report[f"{s['name']}:{k}"] = dict(kind='tapered', capsule_radius=r, node_radii=[round(x, 4) for x in radii])
    # the waist's outline for the waistband
    sp = {s['name']: s for s in spec['segments']}
    ab = sp['abdomen']
    axis_x = ab['origin'][0] + ab['capsules'][0]['a'][0]
    z_w = ab['origin'][2] + WAIST_ABOVE_OMPHALION
    m = trunk_v & (np.abs(V[:, 2] - z_w) < 0.006)
    th = np.mod(np.arctan2(V[m, 1], V[m, 0] - axis_x), 2 * np.pi)
    rad = np.hypot(V[m, 1], V[m, 0] - axis_x)
    nb = 72
    b = np.floor(th / (2 * np.pi) * nb).astype(int) % nb
    R = np.full(nb, np.nan)
    for j in range(nb):
        if (b == j).any():
            R[j] = rad[b == j].max()
    ok = ~np.isnan(R)
    ang = (np.arange(nb) + 0.5) * 2 * np.pi / nb
    R = np.interp(ang, ang[ok], R[ok], period=2 * np.pi)
    R = (np.roll(R, 1) + 2 * R + np.roll(R, -1)) / 4 + LEOTARD_M
    perim = float(np.sum(np.hypot(np.diff(np.r_[R, R[0]] * np.cos(np.r_[ang, ang[0]])),
                                  np.diff(np.r_[R, R[0]] * np.sin(np.r_[ang, ang[0]])))))
    form = dict(
        schema=SCHEMA,
        producer=dict(lane='visualization', file='visualization/dance/dress_form.py',
                      sha256=sha(__file__), dancer_blend=bpy.data.filepath, dancer_sha256=sha(bpy.data.filepath),
                      body_json_sha256=sha(BODY_JSON)),
        evidence='Measured from the renderer\'s fitted CC0 base mesh at body.json\'s rest pose: a collision stand-in '
                 'that keeps dressed renders\' cloth outside the skin they show, until the body model and the '
                 'rendered mesh are fitted as one body (question 1). The dynamics and the garment study keep '
                 'body.json\'s capsules.',
        units='metres, radians',
        reading=dict(pieces='tapered capsules in their segment\'s frame, a and b relative to its rest origin; the '
                            'radius runs linearly from ra at a to rb at b along the axis',
                     waist='polar radius r(theta) of the skin (leotard included) about (axis_x, 0) at height z in the '
                           'rest pose, theta from +x toward +y; linear in theta, periodic'),
        settings=dict(nodes_per_capsule=NODES, limb_percentile=LIMB_PERCENTILE, leotard_m=LEOTARD_M),
        pieces=pieces,
        waist=dict(axis_x=round(float(axis_x), 4), z=round(float(z_w), 4), theta=ang.round(5).tolist(),
                   r=R.round(4).tolist(), perimeter_m=round(perim, 4)),
        capsules=report)
    out.write_text(json.dumps(form, indent=1))
    print('[dress-form] wrote', out, len(pieces), 'pieces; waist perimeter %.3f m, r %.3f-%.3f m'
          % (perim, R.min(), R.max()))
    for k, v in report.items():
        print('[dress-form]', k, v)


main()
