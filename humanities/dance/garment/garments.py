"""Ordinary dance garments cut from flat patterns: circle skirts (a full circle, or cut from part of a circle or from
more than one, with or without a weighted hem), a pair of Chinese water sleeves, and free panels for trailing tests.

Full circle skirt. The pattern is a flat ring whose inner circle is the waist (circumference = the waistband's) and
whose width is the skirt length; centre front lies on the straight grain, so the warp hangs vertical at centre front
and back, the weft at the sides and the bias between them. The waistband is sewn at the natural waist, 2 cm above the
omphalion, just outside the abdomen capsule, and follows the pelvis and abdomen frames equally.

Water sleeves (shuixiu). Each is a tube of silk sewn to the costume's cuff at the wrist and open at its far end, warp
along its length. The cuff ring follows the forearm; the hand sits inside the tube.

Both start flat or straight in a dressing pose (arms raised to the side) and settle before the motion begins.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .body import Attachment, Body, qrot_np
from .patterns import Piece, annulus, develop, gore, grid, rectangle, sector, slit

WAIST_ABOVE_OMPHALION = 0.02       # m: the natural waist, just above the navel (placement)
BAND_CLEARANCE = 0.006             # m: waistband radius beyond the abdomen capsule (or the dress form's waist)


def waist_band(waist, axis_x, z_w, clearance=BAND_CLEARANCE):
    """The waistband on a dress form's waist outline (body.waist), offset by BAND_CLEARANCE: returns its perimeter
    and band(f) -> (body angle, radius) at fraction f of the perimeter from centre back (theta = pi), turning as
    theta increases. Points are spaced by arc length, as the pattern's inner edge is."""
    assert abs(waist['axis_x'] - axis_x) < 1e-3 and abs(waist['z'] - z_w) < 1e-3, 'dress form from another body'
    th0, r0 = np.asarray(waist['theta']), np.asarray(waist['r']) + clearance
    th = np.pi + np.linspace(0.0, 2 * np.pi, 721)
    r = np.interp(np.mod(th, 2 * np.pi), th0, r0, period=2 * np.pi)
    P = np.c_[r * np.cos(th), r * np.sin(th)]
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    perim = float(s[-1])

    def band(f):
        t = np.interp(np.asarray(f) * perim, s, th)
        return t, np.interp(np.mod(t, 2 * np.pi), th0, r0, period=2 * np.pi)
    return perim, band


@dataclass
class Garment:
    name: str
    piece: Piece
    X0: np.ndarray                  # (V, 3) initial positions in the pose at the start of dressing
    pinned: np.ndarray              # (P,) vertex indices
    attachment: Attachment
    params: dict
    extra_mass: np.ndarray | None = None   # (V,) mass sewn on beyond the fabric's own (a weighted hem), kg


def merge(pieces):
    F, UV, rings, vuv, off = [], [], {}, [], 0
    for tag, p in pieces:
        F.append(p.F + off)
        UV.append(p.UV)
        vuv.append(p.vertex_uv)
        for k, v in p.rings.items():
            rings[f'{tag}_{k}'] = [i + off for i in v]
        off += p.n_vertices
    out = Piece(np.concatenate(F), np.concatenate(UV), off, warp=pieces[0][1].warp, rings=rings,
                vertex_uv=np.concatenate(vuv))
    return out


def circle_skirt(body: Body, pose0, length=0.85, h=0.02, seed=0) -> Garment:
    """pose0 = (seg_pos, seg_quat) at the start of dressing. Its band is round: with a dress form, use skirt(cut=2 pi)."""
    assert body.waist is None, 'circle_skirt has a round waistband; skirt() puts the band on a dress form'
    sp = {s['name']: s for s in body.spec['segments']}
    ab = sp['abdomen']
    cap = ab['capsules'][0]
    omph = ab['origin'][2]
    axis_x = ab['origin'][0] + cap['a'][0]
    r_band = cap['radius'] + BAND_CLEARANCE
    z_w = omph + WAIST_ABOVE_OMPHALION
    piece = annulus(r_band, r_band + length, h, seed=seed, warp=(1.0, 0.0))
    inner = np.array(piece.rings['inner'])
    uv = piece.vertex_uv
    centre_rest = np.array([axis_x, 0.0, z_w])
    p_rest = centre_rest + np.c_[uv[inner], np.zeros(len(inner))]
    att = Attachment(body, p_rest, ['pelvis', 'abdomen'], [0.5, 0.5])
    pos, quat = pose0
    tgt = att.targets(pos, quat)
    c0 = tgt.mean(axis=0)
    q = quat[body.index['abdomen']]
    X0 = c0 + qrot_np(np.broadcast_to(q, (piece.n_vertices, 4)), np.c_[uv, np.zeros(piece.n_vertices)])
    X0[inner] = tgt
    params = dict(kind='full circle skirt', length_m=length, waist_radius_m=r_band,
                  waist_circumference_m=2 * np.pi * r_band, hem_circumference_m=2 * np.pi * (r_band + length),
                  natural_waist_height_rest_m=z_w, pattern='flat ring, centre front on the straight grain',
                  attachment='waistband ring pinned to pelvis 0.5 + abdomen 0.5',
                  area_m2=piece.area(), mean_edge_m=piece.edge_length(), vertices=piece.n_vertices,
                  triangles=len(piece.F))
    return Garment('circle_skirt', piece, X0, inner, att, params)


def _cut_slits(piece: Piece, r_in, cut, slits, slit_top, band):
    """Slits up a sector skirt's pattern at body angles `slits` (degrees from centre front toward her left), each along
    the pattern's radial line that hangs from that angle's band point, from the hem up to slit_top down the pattern
    from the waist; they must lie away from the seam at centre back. Returns the cut piece, its 'outer' ring with each
    slit's new hem vertex after the old, and each slit's lips as two vertex lists from its closed end to the hem (the
    end shared), the lip toward centre back first."""
    wrap = lambda a: np.angle(np.exp(1j * np.asarray(a)))
    f = np.linspace(0.0, 1.0, 7201)
    th_f = np.asarray(band(f)[0])
    lips = []
    for deg in slits:
        uv = piece.vertex_uv
        rho = np.hypot(uv[:, 0], uv[:, 1])
        psi = np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi)
        psi_s = f[np.argmin(np.abs(wrap(th_f - math.radians(deg))))] * cut
        path = []
        for r in np.unique(np.round(rho[rho >= r_in + slit_top - 1e-6], 7)):    # one vertex per row, nearest the line
            on = np.where(np.abs(rho - r) < 1e-6)[0]
            path.append(int(on[np.argmin(np.abs(wrap(psi[on] - psi_s)))]))
        on_path = np.zeros(piece.n_vertices, bool)
        on_path[path] = True
        off = ~on_path[piece.F]
        ang = wrap(np.arctan2(piece.UV[..., 1], piece.UV[..., 0]) - psi_s)        # each corner's side of the line
        side = on_path[piece.F].any(1) & ((np.where(off, ang, 0.0).sum(1) / np.maximum(off.sum(1), 1)) > 0)
        piece, copies = slit(piece, path, side)
        outer = list(piece.rings['outer'])
        k = outer.index(path[-1])
        piece.rings['outer'] = outer[:k + 1] + [int(copies[-1])] + outer[k + 1:]
        lips.append((path, [path[0]] + [int(c) for c in copies]))
    return piece, lips


def skirt(body: Body, pose0, cut=2 * np.pi, length=0.85, h=0.02, hem_weight_kg_m=0.0,
          band_clearance=BAND_CLEARANCE, start_clear=None, binding=False, binding_kg_m=0.0023,
          slits=(), slit_top=None) -> Garment:
    """A circle skirt cut from `cut` radians of a ring (pi: a half circle; 2 pi: a full circle; 4 pi: a double
    circle, two full circles sewn at the sides, whose starting placement is not built yet), with the seam at centre
    back and the straight grain at centre front, on the same waistband as circle_skirt. Its pattern's inner radius
    makes the waist: cut x r_in = the waistband's circumference. The skirt starts as the cone it makes when laid
    smooth (flare asin(cut / 2 pi): 30 degrees for a half circle, flat for a full circle), which is isometric to its
    pattern. hem_weight_kg_m adds mass along the hem (a weighted hem tape or chain) that the air does not see.

    On a dress form's waist outline no smooth cone fits the pattern exactly, since the band is not a circle. Each line
    of the pattern starts out from its band point along the outline's outward normal and, over the first 6 cm, turns
    into the pattern's own cone about the outline's centre, so only that blend carries strain (under 16% at the 99th
    percentile). Run out from the abdomen capsule's axis as on a round band, the lines started crowded at the front
    and spread at the back of an outline centred 5 cm forward of that axis, a median strain of 17-20% and up to 180%,
    and the solver's recoil from that start kept a full circle moving through the dressing (2 November). start_clear
    (default: on a dress form) then moves any starting point inside the body out of it, 4 mm clear: a half circle's
    cone, laid from the small of the back, starts partly inside the hips, and on the form that start held the hem up
    at the sides. With binding, the hem's last strip of triangles is a bound edge, listed in params['binding_tris']
    for the binding's own fabric, and the hem carries binding_kg_m more for the fold (as sleeve_tier's hems).

    slits: body angles (degrees from centre front toward her left; 90 and -90 for the sides) of slits cut up the
    skirt from its hem to slit_top down the pattern from the waist, each along the line of the pattern that hangs
    from its band point (_cut_slits). Their lips start closed, as cut, and with binding are bound like the hem: the
    triangles along each lip join params['binding_tris'] and each lip carries binding_kg_m along its length.
    params['slit_lips'] lists each slit's two lips (vertex indices from its closed end to the hem)."""
    sp = {s['name']: s for s in body.spec['segments']}
    ab = sp['abdomen']
    cap = ab['capsules'][0]
    omph = ab['origin'][2]
    axis_x = ab['origin'][0] + cap['a'][0]
    r_band = cap['radius'] + band_clearance
    z_w = omph + WAIST_ABOVE_OMPHALION
    if body.waist is not None:          # the band on the dress form's waist outline; r_band its mean radius
        perim, band = waist_band(body.waist, axis_x, z_w, band_clearance)
        r_band = perim / (2 * np.pi)
    else:
        band = lambda f: (np.pi + 2 * np.pi * np.asarray(f), np.full(np.shape(f), r_band))
    r_in = 2 * np.pi * r_band / cut
    piece = sector(r_in, r_in + length, cut, h)
    lips = []
    if len(slits):
        assert slit_top is not None and 0.0 < slit_top < length, 'a slit runs from slit_top down to the hem'
        piece, lips = _cut_slits(piece, r_in, cut, slits, slit_top, band)
    inner = np.array(piece.rings['inner'])
    hem = np.array(piece.rings['outer'])
    uv = piece.vertex_uv
    rho = np.hypot(uv[:, 0], uv[:, 1])
    psi = np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi)
    sa = min(1.0, cut / (2 * np.pi))
    ca = np.sqrt(1.0 - sa * sa)
    if cut > 2 * np.pi + 1e-9:
        raise NotImplementedError('a skirt cut from more than one circle needs a fluted starting placement')
    th, rb = band(psi / cut)                                # body angle, centre back (the seam) at -x
    d = rho - r_in
    if body.waist is None:
        R = rb + d * sa                                     # rho sa on a round band
        local = np.c_[R * np.cos(th), R * np.sin(th), -d * ca]
    else:                       # out along the outline's normal, turning into the pattern's cone about its centre
        f = psi / cut
        df = 0.02 / perim                                   # the normal smoothed over 2 cm of the outline
        t1, r1 = band(np.mod(f + df, 1.0))
        t0, r0 = band(np.mod(f - df, 1.0))
        T = np.c_[r1 * np.cos(t1) - r0 * np.cos(t0), r1 * np.sin(t1) - r0 * np.sin(t0)]
        N = np.c_[T[:, 1], -T[:, 0]]                        # outward: the outline turns counter-clockwise
        N /= np.linalg.norm(N, axis=1, keepdims=True)
        near = np.c_[rb * np.cos(th), rb * np.sin(th)] + (d * sa)[:, None] * N
        tc, rc = band(np.linspace(0.0, 1.0, 720, endpoint=False))
        centre = np.c_[rc * np.cos(tc), rc * np.sin(tc)].mean(0)
        phi = np.pi + 2 * np.pi * f                         # centre back at f = 0, as on the band
        far = centre + ((r_in + d) * sa)[:, None] * np.c_[np.cos(phi), np.sin(phi)]
        w = np.clip(d / 0.06, 0.0, 1.0)
        w = (w * w * (3 - 2 * w))[:, None]
        local = np.c_[(1 - w) * near + w * far, -d * ca]
    centre_rest = np.array([axis_x, 0.0, z_w])
    th_in, rb_in = th[inner], rb[inner]
    p_rest = centre_rest + np.c_[rb_in * np.cos(th_in), rb_in * np.sin(th_in), np.zeros(len(inner))]
    att = Attachment(body, p_rest, ['pelvis', 'abdomen'], [0.5, 0.5])
    pos, quat = pose0
    tgt = att.targets(pos, quat)
    # the cone is laid about the band's axis, which a band on a waist outline does not have at its mean
    c0 = Attachment(body, centre_rest[None], ['pelvis', 'abdomen'], [0.5, 0.5]).targets(pos, quat)[0]
    q = quat[body.index['abdomen']]
    X0 = c0 + qrot_np(np.broadcast_to(q, (piece.n_vertices, 4)), local)
    X0[inner] = tgt
    if start_clear if start_clear is not None else body.waist is not None:
        X0 = clear_of_body(body, pose0, X0, np.setdiff1d(np.arange(piece.n_vertices), inner))
    extra = None
    if hem_weight_kg_m > 0 or binding:
        extra = np.zeros(piece.n_vertices)
        extra[hem] = (hem_weight_kg_m + (binding_kg_m if binding else 0.0)) * cut * (r_in + length) / len(hem)
    binding_tris = []
    if binding:                                             # the last strip of triangles, between the last two rows
        nr = max(1, int(round(length / (h * np.sqrt(3) / 2))))
        last = rho >= r_in + length * (nr - 1) / nr - 1e-9
        binding_tris = [int(t) for t in np.where(last[piece.F].all(1))[0]]
        if lips:                                            # each lip's triangles, and its fold along its length
            cutv = np.zeros(piece.n_vertices, bool)
            for lip in (l for pair in lips for l in pair):
                cutv[lip[1:]] = True
                extra[lip[1:]] += binding_kg_m * (length - slit_top) / len(lip[1:])
            binding_tris = sorted(set(binding_tris) | {int(t) for t in np.where(cutv[piece.F].any(1))[0]})
    params = dict(kind=f'circle skirt cut from {cut / (2 * np.pi):g} of a circle', cut_rad=cut, length_m=length,
                  waist_radius_m=r_band, waist_circumference_m=2 * np.pi * r_band, pattern_inner_radius_m=r_in,
                  hem_circumference_m=cut * (r_in + length),
                  smooth_cone_flare_deg=float(np.degrees(np.arcsin(sa))),
                  hem_weight_kg_m=hem_weight_kg_m, hem_weight_total_kg=hem_weight_kg_m * cut * (r_in + length),
                  natural_waist_height_rest_m=z_w, bound_hem=binding, binding_tris=binding_tris,
                  slits_deg=[float(s) for s in slits], slit_top_m=slit_top if len(slits) else None, slit_lips=lips,
                  waistband=('on the dress form\'s waist outline + %g m' % band_clearance if body.waist is not None
                             else 'round, the abdomen capsule + %g m' % band_clearance),
                  pattern='ring sector sewn at centre back, centre front on the straight grain',
                  attachment='waistband ring pinned to pelvis 0.5 + abdomen 0.5',
                  area_m2=piece.area(), mean_edge_m=piece.edge_length(), vertices=piece.n_vertices,
                  triangles=len(piece.F))
    return Garment('skirt', piece, X0, inner, att, params, extra)


def clear_of_body(body: Body, pose0, X, free, margin=0.004, rounds=20):
    """Starting positions moved out of the body: each free point inside a body piece (or within `margin` of it) is
    moved along that piece's outward direction to `margin` beyond its surface, repeated while any remain (a point
    pushed out of one piece can land in another)."""
    X = np.array(X, float)
    free = np.asarray(free)
    A, B, r = body.capsules(*pose0)
    ab = B - A
    L2 = np.maximum((ab * ab).sum(1), 1e-12)
    dr = body.cap_rb - body.cap_r
    for _ in range(rounds):
        x = X[free]
        t = np.clip(((x[:, None] - A[None]) * ab[None]).sum(2) / L2[None], 0.0, 1.0)
        c = A[None] + t[..., None] * ab[None]
        dist = np.linalg.norm(x[:, None] - c, axis=2)
        pen = r[None] + t * dr[None] - dist
        k = pen.argmax(1)
        p = pen[np.arange(len(free)), k]
        bad = p > -margin
        if not bad.any():
            break
        n = x[bad] - c[np.where(bad)[0], k[bad]]
        n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-9)
        X[free[bad]] = x[bad] + (p[bad] + margin + 1e-4)[:, None] * n
    return X


def waist_tethers(G: Garment, slack=1.05):
    """Long-range tethers for a skirt: each vertex to the nearest waistband vertex in the flat pattern (across the
    seam too), at `slack` times their distance in the pattern, a limit on how far the cloth below can be pulled away
    from the waist (5% leaves room for a bias that drops)."""
    uv = G.piece.vertex_uv
    inner = np.asarray(G.pinned)
    cut = G.params.get('cut_rad', 2 * np.pi)
    # pattern angle between vertex and band point, the short way round the sewn piece (across the seam if shorter);
    # images of the band turned by the whole cut would overlap the band itself for cuts between a half and a full
    # circle and tether cloth to the wrong part of the waist
    rho = np.hypot(uv[:, 0], uv[:, 1])
    psi = np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi)
    dpsi = np.mod(psi[:, None] - psi[inner][None, :] + cut / 2, cut) - cut / 2
    d = np.sqrt(np.maximum(rho[:, None] ** 2 + rho[inner][None, :] ** 2
                           - 2 * rho[:, None] * rho[inner][None, :] * np.cos(dpsi), 0.0))
    k = d.argmin(1)
    return inner[k], slack * d[np.arange(len(uv)), k]


def cuff_tethers(G: Garment, slack=1.05):
    """Long-range tethers for water sleeves: each vertex to the nearest vertex of its own sleeve's cuff in the flat
    pattern, at `slack` times their distance there, so that a sleeve yanked from the floor by a cast or flung by a
    landing cannot stretch away from the wrist further than silk allows."""
    uv = G.piece.vertex_uv
    half = G.piece.n_vertices // 2
    anchor = np.full(len(uv), -1)
    dist = np.zeros(len(uv))
    for side, lo in (('L', 0), ('R', half)):
        ring = np.array(G.piece.rings[f'{side}_start'])
        idx = np.arange(lo, lo + half)
        d = np.linalg.norm(uv[idx, None, :] - uv[None, ring, :], axis=2)
        k = d.argmin(1)
        anchor[idx] = ring[k]
        dist[idx] = d[np.arange(len(idx)), k]
    return anchor, slack * dist


def water_sleeves(body: Body, pose0, length=1.5, circumference=0.50, h=0.02) -> Garment:
    pieces, p_rest, segs = [], [], []
    for side in ('L', 'R'):
        pieces.append((side, rectangle(circumference, length, h, wrap=True, warp=(0.0, 1.0))))
    piece = merge(pieces)
    uv = piece.vertex_uv
    rc = circumference / (2 * np.pi)
    pos, quat = pose0
    X0 = np.zeros((piece.n_vertices, 3))
    pinned = []
    tgt_all = []
    for side in ('L', 'R'):
        ring = np.array(piece.rings[f'{side}_start'])
        hand = body.spec['segments'][body.index[f'{side}_hand']]
        wrist = np.array(hand['origin'], float)
        phi = 2 * np.pi * uv[ring, 0] / circumference
        pr = wrist + np.c_[rc * np.cos(phi), rc * np.sin(phi), np.zeros(len(ring))]
        pinned.append(ring)
        p_rest.append(pr)
    pinned = np.concatenate(pinned)
    p_rest = np.concatenate(p_rest)
    # one attachment object per side would be tidier; a single one with per-vertex segments is built here
    att = _SideAttachment(body, p_rest, pinned, piece)
    tgt = att.targets(pos, quat)
    half = piece.n_vertices // 2
    for side, lo in (('L', 0), ('R', half)):
        # each vertex starts on the line from its cuff point along the forearm's distal axis
        down = qrot_np(quat[body.index[f'{side}_forearm']], np.array([0.0, 0.0, -1.0]))
        idx = np.arange(lo, lo + half)
        ring = np.array(piece.rings[f'{side}_start'])
        cuff_pts = tgt[np.searchsorted(pinned, ring)]
        X0[idx] = cuff_pts[(idx - lo) % len(ring)] + uv[idx, 1:2] * down
    params = dict(kind='water sleeves (shuixiu), a pair', length_m=length, circumference_m=circumference,
                  cuff_radius_m=rc, attachment='cuff ring pinned to each forearm at the wrist joint centre',
                  area_m2=piece.area(), mean_edge_m=piece.edge_length(), vertices=piece.n_vertices,
                  triangles=len(piece.F))
    return Garment('water_sleeves', piece, X0, pinned, att, params)


class _SideAttachment:
    """Cuff rings: each pinned vertex follows its own side's forearm."""

    def __init__(self, body, p_rest, pinned, piece):
        half = piece.n_vertices // 2
        self.parts = []
        for side, sel in (('L', pinned < half), ('R', pinned >= half)):
            self.parts.append((sel, Attachment(body, p_rest[sel], [f'{side}_forearm'], [1.0])))
        self.n = len(pinned)

    def targets(self, pos, quat):
        out = np.zeros((self.n, 3))
        for sel, a in self.parts:
            out[sel] = a.targets(pos, quat)
        return out


# ----------------------------------------------------------------------------------------------- garments in parts
class CombinedAttachment:
    """The pinned targets of several garments' attachments, in order, for a garment combined from them."""

    def __init__(self, parts):
        self.parts = list(parts)

    def targets(self, pos, quat):
        return np.concatenate([a.targets(pos, quat) for a in self.parts])


def _merge_aligned(pieces, warp=(0.0, 1.0)):
    """Pieces merged into one, each pattern first turned so its warp lies along `warp` (a turn leaves a pattern's
    geometry unchanged), so every piece keeps its own grain in the merged cloth, which has one warp direction; seam
    turns (corner_rot) are carried over, zero for pieces without seams."""
    F, UV, vuv, rot, rings, off = [], [], [], [], {}, 0
    target = math.atan2(warp[1], warp[0])
    for k, p in enumerate(pieces):
        a = target - math.atan2(p.warp[1], p.warp[0])
        Rm = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
        F.append(p.F + off)
        UV.append(p.UV @ Rm.T)
        vuv.append(p.vertex_uv @ Rm.T)
        rot.append(p.corner_rot if p.corner_rot is not None else np.zeros(p.F.shape))
        for key, v in p.rings.items():
            rings[f'p{k}_{key}'] = [i + off for i in v]
        off += p.n_vertices
    any_rot = any(p.corner_rot is not None for p in pieces)
    return Piece(np.concatenate(F), np.concatenate(UV), off, warp=tuple(warp), rings=rings,
                 vertex_uv=np.concatenate(vuv), corner_rot=np.concatenate(rot) if any_rot else None)


def combine(garments, name, params=None):
    """One Garment from several built on the same body and pose (a gown's layers, its petals, its sleeves), simulated
    together so that they collide: pieces merged in order, pins, starting positions and sewn-on mass carried over.
    Returns it with the vertex offset of each part, for per-part fabrics (cloth.merged_tri_material) and tethers."""
    piece = _merge_aligned([G.piece for G in garments])
    offs = np.cumsum([0] + [G.piece.n_vertices for G in garments])
    pinned = np.concatenate([np.asarray(G.pinned) + o for G, o in zip(garments, offs[:-1])])
    extra = None
    if any(G.extra_mass is not None for G in garments):
        extra = np.concatenate([G.extra_mass if G.extra_mass is not None else np.zeros(G.piece.n_vertices)
                                for G in garments])
    params = dict(params or {}, parts=[dict(name=G.name, **G.params) for G in garments])
    G = Garment(name, piece, np.concatenate([G.X0 for G in garments]), pinned,
                CombinedAttachment([G.attachment for G in garments]), params, extra)
    return G, offs


def edge_tethers(G: Garment, ring, slack=1.05):
    """Tethers from every vertex of a piece to the nearest vertex of `ring` in its flat pattern, at `slack` times
    their distance there (a petal to its top edge)."""
    uv = G.piece.vertex_uv
    ring = np.asarray(ring)
    d = np.linalg.norm(uv[:, None, :] - uv[None, ring, :], axis=2)
    k = d.argmin(1)
    return ring[k], slack * d[np.arange(len(uv)), k]


def petals(body: Body, pose0, n=8, top_w=0.165, hem_w=0.40, length=0.80, h=0.015, bias=True, hip_drop=0.10,
           hip_axes=(0.125, 0.17), centre_dx=0.0, hem_weight_kg_m=0.0, start_deg=0.0) -> Garment:
    """A ring of n petals hanging free below the high hip: each an open gore top_w wide at its top edge and hem_w at
    its hem, length long, cut on the bias (warp at 45 degrees to its centre line) or on the straight grain, its top
    edge pinned along an elliptical hip line hip_drop below the natural waist, about the pelvis capsule's axis
    (moved forward by centre_dx), with semi-axes hip_axes (front-back, side to side) chosen to clear the body. Petals
    overlap where n top_w exceeds the hip line; alternate petals are set 4 mm further out so the overlaps start
    clear. Each starts hanging straight down from its pin line. The pin line follows the pelvis frame."""
    sp = {s['name']: s for s in body.spec['segments']}
    pel = sp['pelvis']
    capv = [c for c in pel['capsules'] if abs(c['a'][2] - c['b'][2]) > 0.05][0]
    ab = sp['abdomen']
    z_hip = ab['origin'][2] + WAIST_ABOVE_OMPHALION - hip_drop
    cx = pel['origin'][0] + capv['a'][0] + centre_dx
    angle = (hem_w - top_w) / length
    r_in = top_w / angle
    a_fb, a_lr = hip_axes
    t = np.linspace(0.0, 2 * np.pi, 2001)
    E = np.c_[a_fb * np.cos(t), a_lr * np.sin(t)]
    s_arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(E, axis=0), axis=1))]
    perim = s_arc[-1]
    gs, pinned_rest = [], []
    for k in range(n):
        warp_angle = angle / 2 + (np.pi / 4 if bias else 0.0)
        pc = gore(r_in, r_in + length, angle, h, warp=(np.cos(warp_angle), np.sin(warp_angle)))
        uv = pc.vertex_uv
        rho = np.hypot(uv[:, 0], uv[:, 1])
        psi = np.arctan2(uv[:, 1], uv[:, 0])
        lift = 0.004 * (k % 2)
        # the top edge spans top_w of the hip line, centred at this petal's place round it
        centre_s = (start_deg / 360.0 + k / n) * perim
        s_v = centre_s + (psi - angle / 2) * r_in
        tt = np.interp(np.mod(s_v, perim), s_arc, t)
        base = np.c_[(a_fb + lift) * np.cos(tt), (a_lr + lift) * np.sin(tt)]
        out = base / np.maximum(np.linalg.norm(base, axis=1, keepdims=True), 1e-9)
        drop = rho - r_in
        local = np.c_[base + 0.02 * np.tanh(drop / 0.1)[:, None] * out, -drop]
        p_rest = np.array([cx, 0.0, z_hip]) + local
        inner = np.array(pc.rings['inner'])
        att = Attachment(body, p_rest[inner], ['pelvis'], [1.0])
        pos, quat = pose0
        q = quat[body.index['pelvis']]
        c_now = Attachment(body, np.array([[cx, 0.0, z_hip]]), ['pelvis'], [1.0]).targets(pos, quat)[0]
        X0 = c_now + qrot_np(np.broadcast_to(q, (pc.n_vertices, 4)), local)
        extra = None
        if hem_weight_kg_m > 0:
            hem = np.array(pc.rings['outer'])
            extra = np.zeros(pc.n_vertices)
            extra[hem] = hem_weight_kg_m * hem_w / len(hem)
        gs.append(Garment(f'petal_{k}', pc, X0, inner, att,
                          dict(kind='petal', top_width_m=top_w, hem_width_m=hem_w, length_m=length,
                               grain='bias' if bias else 'straight', vertices=pc.n_vertices), extra))
    G, offs = combine(gs, 'petals', dict(kind=f'{n} petals from the high hip', petals=n, top_width_m=top_w,
                                         hem_width_m=hem_w, length_m=length, grain='bias' if bias else 'straight',
                                         hip_line_m=float(perim), hip_drop_m=hip_drop, hem_weight_kg_m=hem_weight_kg_m))
    anchor = np.full(G.piece.n_vertices, -1)
    length_t = np.zeros(G.piece.n_vertices)
    for g, o in zip(gs, offs[:-1]):
        a, L = edge_tethers(g, g.piece.rings['inner'])
        anchor[o:o + g.piece.n_vertices] = a + o
        length_t[o:o + g.piece.n_vertices] = L
    G.params['tethers'] = (anchor, length_t)
    return G


# ----------------------------------------------------------------------------------------------- a gown's sleeves etc.
def _ring_frame(axis):
    """Two unit vectors at right angles to `axis` and to each other (a ring's plane)."""
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    ref = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(a, ref)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(a, e1)


def arm_sleeves(body: Body, pose0, width=0.50, length=0.58, armhole_r=0.070, armhole_drop=0.04, cuff_r=0.048,
                h=0.015) -> Garment:
    """A pair of full sleeves from the armhole to the wrist: each a tube `width` round and `length` long, pinned round
    the upper arm armhole_drop below the shoulder joint (radius armhole_r) and round the wrist (radius cuff_r), so its
    fullness gathers into both rings and blouses between them (a bishop sleeve); the length leaves room for the
    elbow to bend. The armhole follows the upper arm, the cuff the forearm. Each starts as the straight line of
    cloth from ring to ring."""
    sp = {s['name']: s for s in body.spec['segments']}
    pos, quat = pose0
    gs = []
    for side in ('L', 'R'):
        ua, fa, hand = sp[f'{side}_upper_arm'], sp[f'{side}_forearm'], sp[f'{side}_hand']
        cap = ua['capsules'][0]
        axis = np.array(cap['b'], float) - np.array(cap['a'], float)
        axis /= np.linalg.norm(axis)
        c_top = np.array(ua['origin'], float) + armhole_drop * axis
        c_cuff = np.array(hand['origin'], float)
        fcap = fa['capsules'][0]
        faxis = np.array(fcap['b'], float) - np.array(fcap['a'], float)
        pc = rectangle(width, length, h, wrap=True, warp=(0.0, 1.0))
        uv = pc.vertex_uv
        top, end = np.array(pc.rings['start']), np.array(pc.rings['end'])
        phi = 2 * np.pi * uv[:, 0] / width
        e1, e2 = _ring_frame(axis)
        f1, f2 = _ring_frame(faxis)
        ring_top = c_top + armhole_r * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
        ring_cuff = c_cuff + cuff_r * (np.cos(phi)[:, None] * f1 + np.sin(phi)[:, None] * f2)
        a_top = Attachment(body, ring_top[top], [f'{side}_upper_arm'], [1.0])
        a_cuff = Attachment(body, ring_cuff[end], [f'{side}_forearm'], [1.0])
        A = Attachment(body, ring_top, [f'{side}_upper_arm'], [1.0]).targets(pos, quat)
        C = Attachment(body, ring_cuff, [f'{side}_forearm'], [1.0]).targets(pos, quat)
        s = (uv[:, 1] / length)[:, None]
        X0 = (1 - s) * A + s * C
        pinned = np.concatenate([top, end])
        order = np.argsort(pinned)
        att = CombinedAttachment([a_top, a_cuff])
        if not (order == np.arange(len(order))).all():
            raise ValueError('sleeve rings out of index order')
        gs.append(Garment(f'{side}_arm_sleeve', pc, X0, pinned, att,
                          dict(kind='full sleeve, armhole to wrist', width_m=width, length_m=length,
                               armhole_radius_m=armhole_r, cuff_radius_m=cuff_r, vertices=pc.n_vertices)))
    G, offs = combine(gs, 'arm_sleeves', dict(kind='a pair of full sleeves gathered at armhole and wrist', width_m=width,
                                              length_m=length, armhole_radius_m=armhole_r, cuff_radius_m=cuff_r))
    return G


def cuff_flounces(body: Body, pose0, ring_r=0.078, cut=1.5 * np.pi, length=0.22, h=0.015) -> Garment:
    """A pair of flounces at the wrists: each a circle cut from `cut` radians, sewn at its ends into the cone it makes,
    its inner edge pinned on a ring of radius ring_r round the wrist joint and opening along the forearm toward and
    over the hand, like a small skirt; several, graded in length, make tiers. Each follows its forearm."""
    sp = {s['name']: s for s in body.spec['segments']}
    pos, quat = pose0
    gs = []
    sa = min(1.0, cut / (2 * np.pi))
    ca = np.sqrt(1.0 - sa * sa)
    for side in ('L', 'R'):
        fa, hand = sp[f'{side}_forearm'], sp[f'{side}_hand']
        fcap = fa['capsules'][0]
        faxis = np.array(fcap['b'], float) - np.array(fcap['a'], float)
        faxis /= np.linalg.norm(faxis)
        c = np.array(hand['origin'], float)
        r_in = 2 * np.pi * ring_r / cut
        pc = sector(r_in, r_in + length, cut, h)
        uv = pc.vertex_uv
        rho = np.hypot(uv[:, 0], uv[:, 1])
        th = 2 * np.pi * np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi) / cut
        e1, e2 = _ring_frame(faxis)
        radial = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
        d = (rho - r_in)[:, None]
        rest = c + ring_r * radial + d * (sa * radial + ca * faxis)
        inner = np.array(pc.rings['inner'])
        att = Attachment(body, rest[inner], [f'{side}_forearm'], [1.0])
        X0 = Attachment(body, rest, [f'{side}_forearm'], [1.0]).targets(pos, quat)
        gs.append(Garment(f'{side}_flounce', pc, X0, inner, att,
                          dict(kind='wrist flounce', cut_rad=cut, length_m=length, ring_radius_m=ring_r,
                               vertices=pc.n_vertices)))
    G, offs = combine(gs, 'cuff_flounces', dict(kind='a pair of wrist flounces', cut_rad=cut, length_m=length,
                                                ring_radius_m=ring_r))
    G.params['offsets'] = [int(o) for o in offs]
    return G


def sleeve_tier(body: Body, pose0, length=0.75, top=0.44, ring_r=0.055, flare=0.2, slant=0.0, binding=True,
                h=0.015, binding_kg_m=0.0023) -> Garment:
    """One tier of a tiered sleeve end on each wrist: a sheer tube `length` long, gathered at its top into a ring of
    radius ring_r round the wrist joint (at right angles to the forearm, which it follows), and widening by `flare` of
    its top circumference `top` toward its hem, so that it is a cone cut from a narrow ring sector. Several tiers,
    graded in length and each a little fuller than the one inside it, echo a tiered skirt. slant > 0 cuts the hem on
    a slant, longer by slant x length at the back of the arm and shorter by as much at the front, so that the tier
    falls to a soft point. With binding, the hem's last strip of triangles is a bound edge, listed in
    params['binding_tris'] for the binding's own fabric, and the hem carries binding_kg_m more for the fold. Each tier
    starts as its cone hanging along the forearm, its top gathered into the ring and widening to its own width over at
    least 3 cm, and never faster than half its run down the cloth (a full tier on a small ring cannot widen within 3 cm
    of cloth at all)."""
    sp = {s_['name']: s_ for s_ in body.spec['segments']}
    pos, quat = pose0
    angle = flare * top / length                       # the circumference grows by `angle` per metre of length
    r_in = top / angle
    sa = angle / (2 * np.pi)                           # sine of the cone's half-angle
    ca = np.sqrt(1.0 - sa * sa)
    gs = []
    for side in ('L', 'R'):
        fa, hand = sp[f'{side}_forearm'], sp[f'{side}_hand']
        fcap = fa['capsules'][0]
        faxis = np.array(fcap['b'], float) - np.array(fcap['a'], float)
        faxis /= np.linalg.norm(faxis)
        c = np.array(hand['origin'], float)
        e1, e2 = _ring_frame(faxis)
        back = np.array([-1.0, 0.0, 0.0])
        th_back = np.arctan2(back @ e2, back @ e1)
        hem = None
        if slant > 0:
            hem = lambda u: length * (1.0 + slant * np.cos(2 * np.pi * np.asarray(u) - th_back))
        pc = sector(r_in, r_in + length * (1.0 + slant), angle, h, hem=hem)
        uv = pc.vertex_uv
        rho = np.hypot(uv[:, 0], uv[:, 1])
        u = np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi) / angle
        th = 2 * np.pi * u
        d = rho - r_in                                 # distance down the cloth from the top
        radial = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
        r_free = top / (2 * np.pi)
        d_g = max(0.03, 2.0 * (r_free - ring_r))         # the gather widens at most half as fast as it runs down
        R = ring_r + (r_free - ring_r) * np.clip(d / d_g, 0.0, 1.0) + d * sa
        slope = (r_free - ring_r) / d_g + sa             # the cloth's outward run per metre down it, in the gather
        along = np.where(d < d_g, d * np.sqrt(1.0 - slope ** 2), d_g * np.sqrt(1.0 - slope ** 2) + (d - d_g) * ca)
        rest = c + R[:, None] * radial + along[:, None] * faxis
        inner = np.array(pc.rings['inner'])
        hem_v = np.array(pc.rings['outer'])
        extra = None
        bt = []
        if binding:
            nr = max(1, int(round(length * (1.0 + slant) / (h * np.sqrt(3) / 2))))
            full = length * (1.0 + slant * np.cos(th - th_back)) if slant > 0 else np.full(len(d), length)
            last = d / full >= (nr - 1) / nr - 1e-6
            bt = np.where(last[pc.F].all(1))[0]
            extra = np.zeros(pc.n_vertices)
            hem_len = top * (1.0 + flare)
            extra[hem_v] = binding_kg_m * hem_len / len(hem_v)
        att = Attachment(body, rest[inner], [f'{side}_forearm'], [1.0])
        X0 = Attachment(body, rest, [f'{side}_forearm'], [1.0]).targets(pos, quat)
        gs.append(Garment(f'{side}_sleeve_tier', pc, X0, inner, att,
                          dict(kind='sleeve tier', length_m=length, top_circumference_m=top, flare=flare,
                               slant=slant, ring_radius_m=ring_r, vertices=pc.n_vertices,
                               binding_tris=[int(t) for t in bt]), extra))
    G, offs = combine(gs, 'sleeve_tier', dict(kind='a pair of sleeve tiers at the wrists', length_m=length,
                                              top_circumference_m=top, hem_circumference_m=top * (1.0 + flare),
                                              flare=flare, slant=slant, ring_radius_m=ring_r, bound_hem=binding))
    n_tri = np.cumsum([0] + [len(g.piece.F) for g in gs])
    G.params['binding_tris'] = [int(t) + int(o) for g, o in zip(gs, n_tri[:-1]) for t in g.params['binding_tris']]
    for part in G.params['parts']:
        part.pop('binding_tris', None)
    return G


def butterfly_shoulders(body: Body, pose0, width=0.30, height=0.12, h=0.008, tilt=0.35, pleat=0.03, cap_r=0.075,
                        arc_deg=(30.0, 150.0)) -> Garment:
    """A pair of stiff gathered sleeve caps after the Philippine terno's butterfly sleeves: each a rectangle `width`
    long and `height` high whose lower edge is gathered into an arch over the shoulder joint, from the front of the
    armhole over the top of the deltoid to the back (cap_r from the joint, arc_deg in the front-back plane, about half
    the panel's length), so that the crisp fabric stands up and out from the shoulder like a wing and leans outward by
    `tilt`. The arch is part of the sleeve and follows the upper arm. Each starts folded as its gathers hold it: an
    accordion of folds `pleat` of fabric each, laid on whole mesh columns, isometric to the flat panel (laid flat and
    squeezed, the solver fans a gathered panel out in its own plane and overstretches it)."""
    sp = {s_['name']: s_ for s_ in body.spec['segments']}
    pos, quat = pose0
    gs = []
    th = np.radians(np.asarray(arc_deg, float))
    for side, sgn in (('L', 1.0), ('R', -1.0)):
        c = np.array(sp[f'{side}_upper_arm']['origin'], float) + np.array([0.0, sgn * 0.015, 0.0])
        pc = rectangle(width, height, h, wrap=False, warp=(0.0, 1.0))
        uv = pc.vertex_uv
        nu = max(2, int(round(width / h)))
        pl = max(1, int(round(pleat / (width / nu)))) * width / nu
        arc_len = cap_r * (th[1] - th[0])
        ratio = width / arc_len
        fold_no = np.floor(uv[:, 0] / pl + 1e-9)
        t = uv[:, 0] / pl - fold_no
        ang = th[0] + (fold_no + t) * pl / width * (th[1] - th[0])
        radial = np.c_[np.cos(ang), np.zeros(len(ang)), np.sin(ang)]
        # the panel rises at right angles to the arch, out along its radius and sideways, so that each height's
        # arch has radius cap_r + v lean; it must not outgrow the panel's width (lean <= (ratio - 1) cap_r / height),
        # and within that limit the panel leans out by `tilt`. The folds open as the arch widens upward, their depth
        # set at each height so that the folded cloth keeps its cut length.
        lean = min(1.0 / np.sqrt(1.0 + tilt * tilt), 0.95 * (ratio - 1.0) * cap_r / height)
        rise = lean * radial + np.sqrt(1.0 - lean * lean) * np.array([0.0, sgn, 0.0])
        spread = (cap_r + uv[:, 1] * lean) / cap_r
        amp = 0.5 * pl * np.sqrt(np.clip(1.0 - (spread / ratio) ** 2, 0.0, None))
        w = amp * np.where(fold_no % 2 == 0, 2 * t - 1, 1 - 2 * t)
        tangent = np.c_[-np.sin(ang), np.zeros(len(ang)), np.cos(ang)]
        fold = np.cross(rise, tangent)
        fold /= np.linalg.norm(fold, axis=1, keepdims=True)
        rest = c + cap_r * radial + uv[:, 1:2] * rise + w[:, None] * fold
        start = np.array(pc.rings['start'])
        att = Attachment(body, rest[start], [f'{side}_upper_arm'], [1.0])
        X0 = Attachment(body, rest, [f'{side}_upper_arm'], [1.0]).targets(pos, quat)
        gs.append(Garment(f'{side}_butterfly', pc, X0, start, att,
                          dict(kind='gathered standing sleeve cap', width_m=width, height_m=height,
                               seam_m=float(arc_len), vertices=pc.n_vertices)))
    G, offs = combine(gs, 'butterfly_shoulders', dict(kind='a pair of terno-like sleeve caps', width_m=width,
                                                      height_m=height))
    return G


def sash(body: Body, pose0, clearance=0.048, band_w=0.05, tail_w=0.07, tail_len=0.85, h=0.015) -> Garment:
    """A sash at the natural waist: a band band_w wide held round the waist `clearance` beyond it (outside the skirts'
    bands), every point of it following pelvis and abdomen equally, and two tails tail_w wide and tail_len long hanging
    free from it at the centre back. On a dress form with a waist outline the band follows the outline, so the sash
    fits the rendered waist; otherwise it is round, `clearance` beyond the abdomen capsule. A round band about the
    capsule's axis stands off the rendered waist's sides and back (by 5 and 13 cm at 4.8 cm clearance, 2 November)."""
    sp = {s['name']: s for s in body.spec['segments']}
    ab = sp['abdomen']
    cap = ab['capsules'][0]
    axis_x = ab['origin'][0] + cap['a'][0]
    z_w = ab['origin'][2] + WAIST_ABOVE_OMPHALION
    if body.waist is not None:
        perim, band_at = waist_band(body.waist, axis_x, z_w, clearance)
    else:
        perim = 2 * np.pi * (cap['radius'] + clearance)
        band_at = lambda f: (np.pi + 2 * np.pi * np.asarray(f), np.full(np.shape(f), cap['radius'] + clearance))
    r = perim / (2 * np.pi)                             # the band's mean radius
    pos, quat = pose0
    segs, wts = ['pelvis', 'abdomen'], [0.5, 0.5]
    band = rectangle(perim, band_w, h, wrap=True, warp=(1.0, 0.0))
    uv = band.vertex_uv
    phi, rb = band_at(uv[:, 0] / perim)                 # centre back at u = 0
    rest = np.c_[axis_x + rb * np.cos(phi), rb * np.sin(phi), z_w + band_w / 2 - uv[:, 1]]
    all_v = np.arange(band.n_vertices)
    gs = [Garment('sash_band', band, Attachment(body, rest, segs, wts).targets(pos, quat), all_v,
                  Attachment(body, rest, segs, wts), dict(kind='sash band, held', radius_m=r, width_m=band_w,
                                                         vertices=band.n_vertices))]
    for k, dphi in enumerate((-0.12, 0.12)):
        tail = rectangle(tail_w, tail_len, h, wrap=False, warp=(0.0, 1.0))
        tv = tail.vertex_uv
        phi0, r0 = band_at(np.mod(dphi * r / perim, 1.0))
        out = np.array([np.cos(phi0), np.sin(phi0), 0.0])
        tang = np.array([-np.sin(phi0), np.cos(phi0), 0.0])
        base = np.array([axis_x, 0.0, z_w]) + (r0 + 0.012) * out
        trest = base + (tv[:, 0:1] - tail_w / 2) * tang - tv[:, 1:2] * np.array([0.0, 0.0, 1.0]) + 0.004 * k * out
        start = np.array(tail.rings['start'])
        gs.append(Garment(f'sash_tail_{k}', tail, Attachment(body, trest, segs, wts).targets(pos, quat), start,
                          Attachment(body, trest[start], segs, wts),
                          dict(kind='sash tail', width_m=tail_w, length_m=tail_len, vertices=tail.n_vertices)))
    G, offs = combine(gs, 'sash', dict(kind='sash with two tails at the back', radius_m=r, band_width_m=band_w,
                                       tail_width_m=tail_w, tail_length_m=tail_len))
    G.params['offsets'] = [int(o) for o in offs]
    return G


class _BlendAttachment:
    """Pinned vertices that follow two segments, each vertex with its own share of the first."""

    def __init__(self, body: Body, p_rest, segments, first_share):
        self.parts = [Attachment(body, p_rest, [s], [1.0]) for s in segments]
        self.w = np.asarray(first_share, float)[:, None]

    def targets(self, pos, quat):
        return self.w * self.parts[0].targets(pos, quat) + (1.0 - self.w) * self.parts[1].targets(pos, quat)


def radial_envelope(X, axis_x, n_theta=72, dz=0.01):
    """How far cloth already dressed reaches round the body: a function (theta, z) -> radius about the vertical line
    x = axis_x, y = 0, the largest radius of the points X (rest frame) in that angle's and height's bin and the bins
    next to it, zero where there is no cloth. A garment put on over others is laid outside it."""
    X = np.asarray(X, float)
    th = np.mod(np.arctan2(X[:, 1], X[:, 0] - axis_x), 2 * np.pi)
    r = np.hypot(X[:, 0] - axis_x, X[:, 1])
    z0 = X[:, 2].min() - 2 * dz
    nz = int(np.ceil((X[:, 2].max() - z0) / dz)) + 3
    E = np.zeros((n_theta, nz))
    np.maximum.at(E, ((th / (2 * np.pi) * n_theta).astype(int) % n_theta, ((X[:, 2] - z0) / dz).astype(int)), r)
    D = E.copy()
    for di in (-1, 0, 1):                     # the empty bins at both ends in height make the roll there harmless
        for dj in (-1, 0, 1):
            D = np.maximum(D, np.roll(np.roll(E, di, 0), dj, 1))

    def env(theta, z):
        theta, z = np.broadcast_arrays(np.asarray(theta, float), np.asarray(z, float))
        i = (np.mod(theta, 2 * np.pi) / (2 * np.pi) * n_theta).astype(int) % n_theta
        j = np.floor((z - z0) / dz).astype(int)
        return np.where((j >= 0) & (j < nz), D[i, np.clip(j, 0, nz - 1)], 0.0)
    return env


def _fold_profile(gather, folds, weight, n=1601):
    """Across a band gathered into `folds` folds of depths set by weight(s): the chord fraction s, the outward height d
    (in chords) and the fraction of the cloth's width t at n points, for the profile d = a weight(s) sin^2(pi folds s)
    whose arc is `gather` chords."""
    s = np.linspace(0.0, 1.0, n)
    shape = weight(s) * np.sin(np.pi * folds * s) ** 2

    def arc(a):
        return np.r_[0.0, np.cumsum(np.hypot(np.diff(s), a * np.diff(shape)))]
    lo, hi = 0.0, 1.0
    while arc(hi)[-1] < gather:
        hi *= 2.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if arc(mid)[-1] < gather else (lo, mid)
    L = arc(hi)
    return s, hi * shape, L / L[-1]


def _clear_radially(body: Body, pose, P, axis_x, margin, rounds=40, skip=('upper_arm', 'forearm', 'hand')):
    """Points moved horizontally away from the vertical line x = axis_x, y = 0 until each is `margin` clear of the
    body's pieces other than those of segments named with `skip` (in the rest pose the arms hang beside the waist),
    keeping their heights (clear_of_body moves points along a piece's normal, which near the ends of the trunk's
    pieces turns up or down)."""
    P = np.array(P, float)
    keep = np.array([not any(k in n for k in skip) for n in body.cap_name])
    A, B, r = (x[keep] for x in body.capsules(*pose))
    ab = B - A
    L2 = np.maximum((ab * ab).sum(1), 1e-12)
    dr = (body.cap_rb - body.cap_r)[keep]
    for _ in range(rounds):
        t = np.clip(((P[:, None] - A[None]) * ab[None]).sum(2) / L2[None], 0.0, 1.0)
        dist = np.linalg.norm(P[:, None] - (A[None] + t[..., None] * ab[None]), axis=2)
        pen = (r[None] + t * dr[None] - dist).max(1)
        bad = pen > -margin
        if not bad.any():
            break
        out = np.c_[P[bad, 0] - axis_x, P[bad, 1], np.zeros(bad.sum())]
        out /= np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-9)
        P[bad] += (pen[bad] + margin + 1e-4)[:, None] * out
    return P


def sash_band(body: Body, knot_deg=55.0, top=(0.025, 0.02), bottom=(0.035, 0.06), clearance=0.018):
    """Where a draped sash's band lies (rest frame): its held edges' heights z_top(theta) and z_bot(theta), top[0]
    above and bottom[0] below the natural waist, falling by top[1] and bottom[1] more toward the knot at knot_deg from
    centre front toward the left, and base(theta, z), the points `clearance` off the body at those angles and heights,
    moved straight out from the waist outline. Returns (z_top, z_bot, base, axis_x)."""
    sp = {s['name']: s for s in body.spec['segments']}
    ab = sp['abdomen']
    cap = ab['capsules'][0]
    axis_x = ab['origin'][0] + cap['a'][0]
    z_w = ab['origin'][2] + WAIST_ABOVE_OMPHALION
    rest_pose = (body.origin, np.tile([1.0, 0.0, 0.0, 0.0], (len(body.names), 1)))
    if body.waist is not None:
        th0, r0 = np.asarray(body.waist['theta']), np.asarray(body.waist['r'])
        r_waist = lambda th: np.interp(np.mod(th, 2 * np.pi), th0, r0, period=2 * np.pi)
    else:
        r_waist = lambda th: np.full(np.shape(th), float(cap['radius']))
    tk = math.radians(knot_deg)
    bump = lambda th: ((1.0 + np.cos(np.asarray(th, float) - tk)) / 2.0) ** 2
    z_top = lambda th: z_w + top[0] - top[1] * bump(th)
    z_bot = lambda th: z_w - bottom[0] - bottom[1] * bump(th)

    def base(th, z):
        th, z = np.broadcast_arrays(np.asarray(th, float), np.asarray(z, float))
        rw = r_waist(th)
        P = np.stack([axis_x + rw * np.cos(th), rw * np.sin(th), z], -1)
        return _clear_radially(body, rest_pose, P.reshape(-1, 3), axis_x, clearance).reshape(P.shape)
    return z_top, z_bot, base, axis_x


def draped_sash(body: Body, pose0, knot_deg=55.0, top=(0.025, 0.02), bottom=(0.035, 0.06), gather=1.9, folds=4,
                clearance=0.018, knot_r=0.028, knot_tilt_deg=25.0,
                tails=((0.55, 0.06, 0.13, 0.03), (0.75, 0.06, 0.13, -0.03)), tail_slant_deg=35.0, h=0.006,
                h_tail=0.012, envelope=None, margin=0.018, fold_drift=1.0, seed=0, hold='edges', tack_every=0.08,
                share=None) -> Garment:
    """A draped sash knotted at the hip. The band (sash_band) is a strip of soft satin gathered across its width into
    `folds` folds of unequal depth, which drift across it fold_drift times round the waist (its cut width `gather` times the 6-10 cm it covers), held along both edges
    `clearance` off the body on a slant that falls toward the knot at the left hip, so it covers the skirts' waistbands;
    it is cut to its longest rows (patterns.develop), so it starts compressed where its folds' hollows run shorter and
    nowhere stretched. Both its ends go into the knot, a held bundle of wrapped satin knot_r round with spiral folds,
    its long axis tilted knot_tilt_deg below the waist's line toward the front. From the knot hang tails (length, width
    at the knot, width at the end, offset along the waist), free, cut on a slant of tail_slant_deg at the end. Held
    parts follow the pelvis and abdomen, the abdomen's share falling from half at the top edge to a fifth at the knot.

    The sash goes on over the skirts, as one tied on after the skirt hangs: `envelope(theta, z) -> r` (radial_envelope
    of the dressed skirts, rest frame) tilts each tail out from the hip until it starts `margin` outside their cloth.
    The skirts' cloth under the band is pressed under it by the run that dresses them (it lies within the band's
    reach of the body, sash_band).

    hold='top' holds the band's top edge all round and tacks its bottom edge every tack_every round the waist, so the
    gathered cloth between the tacks can settle into its own folds; 'edges' (the default) holds both edges all round.
    share sets one pelvis share for every held part (band, knot and tails), so that a twist of the torso turns the
    sash as one piece; None keeps the shares above."""
    z_top, z_bot, base, axis_x = sash_band(body, knot_deg, top, bottom, clearance)
    env = envelope if envelope is not None else (lambda th, z: np.zeros(np.broadcast(th, z).shape))
    tk = math.radians(knot_deg)
    segs = ['pelvis', 'abdomen']
    outward = lambda th: np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], -1)

    # the band: columns round the waist from the knot back to it, rows across its gathered width
    r_k = float(np.hypot(*(base(tk, z_top(tk))[:2] - [axis_x, 0.0])))
    gap = 0.8 * knot_r / r_k
    th_f = np.linspace(tk + gap, tk + 2 * np.pi - gap, 721)
    Tf = base(th_f, z_top(th_f))
    n_i = max(8, int(round(np.linalg.norm(np.diff(Tf, axis=0), axis=1).sum() / h)) + 1)
    th_c = np.linspace(tk + gap, tk + 2 * np.pi - gap, n_i)
    T, B = base(th_c, z_top(th_c)), base(th_c, z_bot(th_c))
    chord = np.linalg.norm(B - T, axis=1)
    n_j = max(5, int(round(gather * chord.max() / h)) + 1)
    t_rows = np.linspace(0.0, 1.0, n_j)
    phase = np.random.default_rng(seed).uniform(0.0, 2 * np.pi)
    X = np.zeros((n_i, n_j, 3))
    for i in range(n_i):
        c = B[i] - T[i]
        cdir = c / chord[i]
        out = outward(th_c[i]) - (outward(th_c[i]) @ cdir) * cdir
        out /= np.linalg.norm(out)
        weight = lambda s, i=i: 1.0 + 0.3 * np.sin(2 * np.pi * s + phase + fold_drift * th_c[i])   # unequal folds
        sf, d, t = _fold_profile(gather, folds, weight)
        X[i] = T[i] + np.interp(t_rows, t, sf)[:, None] * c + (chord[i] * np.interp(t_rows, t, d))[:, None] * out
    band = grid(develop(X, common_rows=True), warp=(1.0, 0.0))
    Xb = X.transpose(1, 0, 2).reshape(-1, 3)              # vertex k is grid point (k % n_i, k // n_i)
    assert hold in ('edges', 'top'), hold
    one = share
    share = 0.5 + 0.3 * np.repeat(t_rows, n_i)            # the pelvis's share: half at the top edge, 0.8 at the bottom
    ends = np.r_[band.rings['left'], band.rings['right']]
    share[ends] = 0.8
    if one is not None:
        share[:] = one
    k_tail = 0.8 if one is None else one
    held_bottom = band.rings['end'] if hold == 'edges' else band.rings['end'][::max(1, int(round(tack_every / h)))]
    pins = np.unique(np.r_[band.rings['start'], held_bottom, ends])
    gs = [Garment('sash_band', band, _BlendAttachment(body, Xb, segs, share).targets(*pose0), pins,
                  _BlendAttachment(body, Xb[pins], segs, share[pins]),
                  dict(kind='draped sash band, held along both edges' if hold == 'edges' else
                       f'draped sash band, held along its top edge and tacked every {tack_every:g} m along its bottom',
                       gather=gather, folds=folds,
                       width_on_body_m=[float(chord.min()), float(chord.max())], vertices=band.n_vertices))]

    # the knot: a held bundle over the band's ends, its long axis along the waist tilted down toward the front
    zk = 0.5 * (float(z_top(tk)) + float(z_bot(tk)))
    bk = base(tk, zk)
    out_k, tan_k, down = outward(np.float64(tk)), np.array([-math.sin(tk), math.cos(tk), 0.0]), np.array([0.0, 0.0, -1.0])
    tilt = math.radians(knot_tilt_deg)
    e1 = -math.cos(tilt) * tan_k + math.sin(tilt) * down       # toward the front and down
    e3 = out_k - (out_k @ e1) * e1
    e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3, e1)
    centre = bk + 0.8 * knot_r * e3
    u = np.linspace(0.1, 0.9, 9)
    v = np.linspace(0.0, 2 * np.pi, 17)
    R = knot_r * np.sqrt(np.sin(np.pi * u))[None, :] * (1.0 + 0.15 * np.sin(3.0 * v[:, None] + 3 * np.pi * u[None, :]))
    flat_back = np.where(np.sin(v) < 0, 0.35, 1.0)[:, None, None]        # the side against the band lies flatter
    K = (centre + (2 * u - 1)[None, :, None] * 1.1 * knot_r * e1
         + R[..., None] * (np.cos(v)[:, None, None] * e2 + flat_back * np.sin(v)[:, None, None] * e3))
    knot = grid(develop(K), warp=(1.0, 0.0), wrap=True)
    Xk = K[:-1].transpose(1, 0, 2).reshape(-1, 3)
    held_k = _BlendAttachment(body, Xk, segs, np.full(len(Xk), k_tail))
    gs.append(Garment('sash_knot', knot, held_k.targets(*pose0), np.arange(knot.n_vertices), held_k,
                      dict(kind='sash knot, held', radius_m=knot_r, vertices=knot.n_vertices)))

    # the tails: planar trapezoids hanging from under the knot, tilted out from the hip until they clear the skirts
    slant = math.tan(math.radians(tail_slant_deg))
    for k, (L, w0, w1, side) in enumerate(tails):
        n_u = max(3, int(round(w1 / h_tail)) + 1)
        n_v = max(3, int(round(L / h_tail)) + 1)
        uu, vv = np.linspace(-0.5, 0.5, n_u), np.linspace(0.0, 1.0, n_v)
        flat = np.stack([uu[:, None] * (w0 + (w1 - w0) * vv[None, :]),
                         vv[None, :] * (L + uu[:, None] * w1 * slant)], -1)          # (n_u, n_v, 2)
        top_c = centre - 0.7 * knot_r * e3 + 0.6 * knot_r * down + side * tan_k + (0.4 * knot_r + 0.012 * k) * out_k
        for alpha in np.arange(0.1, 1.2, 0.02):
            dirn = math.cos(alpha) * down + math.sin(alpha) * out_k
            Xt = top_c + flat[..., :1] * tan_k + flat[..., 1:] * dirn
            th_t = np.arctan2(Xt[..., 1], Xt[..., 0] - axis_x)
            r_t = np.hypot(Xt[..., 0] - axis_x, Xt[..., 1])
            low = flat[..., 1] > 2.5 * knot_r                  # below the knot's own reach
            if (r_t[low] >= env(th_t[low], Xt[..., 2][low]) + margin).all():
                break
        tail = grid(flat, warp=(0.0, 1.0))
        Xtv = Xt.transpose(1, 0, 2).reshape(-1, 3)
        top_ids = np.asarray(tail.rings['start'])
        gs.append(Garment(f'sash_tail_{k}', tail,
                          _BlendAttachment(body, Xtv, segs, np.full(len(Xtv), k_tail)).targets(*pose0), top_ids,
                          _BlendAttachment(body, Xtv[top_ids], segs, np.full(len(top_ids), k_tail)),
                          dict(kind='sash tail', length_m=L, width_m=[w0, w1], tilt_deg=float(np.degrees(alpha)),
                               vertices=tail.n_vertices)))
    G, offs = combine(gs, 'draped_sash', dict(kind='draped sash knotted at the hip', knot_deg=knot_deg, gather=gather,
                                              folds=folds, knot_radius_m=knot_r,
                                              over_dressed_cloth=envelope is not None))
    G.params['offsets'] = [int(o) for o in offs]
    return G
