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

from dataclasses import dataclass

import numpy as np

from .body import Attachment, Body, qrot_np
from .patterns import Piece, annulus, rectangle, sector

WAIST_ABOVE_OMPHALION = 0.02       # m: the natural waist, just above the navel (placement)
BAND_CLEARANCE = 0.006             # m: waistband radius beyond the abdomen capsule


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
    """pose0 = (seg_pos, seg_quat) at the start of dressing."""
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


def skirt(body: Body, pose0, cut=2 * np.pi, length=0.85, h=0.02, hem_weight_kg_m=0.0) -> Garment:
    """A circle skirt cut from `cut` radians of a ring (pi: a half circle; 2 pi: a full circle; 4 pi: a double
    circle, two full circles sewn at the sides, whose starting placement is not built yet), with the seam at centre
    back and the straight grain at centre front, on the same waistband as circle_skirt. Its pattern's inner radius
    makes the waist: cut x r_in = the waistband's circumference. The skirt starts as the cone it makes when laid
    smooth (flare asin(cut / 2 pi): 30 degrees for a half circle, flat for a full circle), which is isometric to its
    pattern. hem_weight_kg_m adds mass along the hem (a weighted hem tape or chain) that the air does not see."""
    sp = {s['name']: s for s in body.spec['segments']}
    ab = sp['abdomen']
    cap = ab['capsules'][0]
    omph = ab['origin'][2]
    axis_x = ab['origin'][0] + cap['a'][0]
    r_band = cap['radius'] + BAND_CLEARANCE
    z_w = omph + WAIST_ABOVE_OMPHALION
    r_in = 2 * np.pi * r_band / cut
    piece = sector(r_in, r_in + length, cut, h)
    inner = np.array(piece.rings['inner'])
    hem = np.array(piece.rings['outer'])
    uv = piece.vertex_uv
    rho = np.hypot(uv[:, 0], uv[:, 1])
    psi = np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi)
    sa = min(1.0, cut / (2 * np.pi))
    ca = np.sqrt(1.0 - sa * sa)
    if cut > 2 * np.pi + 1e-9:
        raise NotImplementedError('a skirt cut from more than one circle needs a fluted starting placement')
    th = np.pi + psi * 2 * np.pi / cut                      # body angle, centre back (the seam) at -x
    local = np.c_[rho * sa * np.cos(th), rho * sa * np.sin(th), -(rho - r_in) * ca]
    centre_rest = np.array([axis_x, 0.0, z_w])
    th_in = np.pi + psi[inner] * 2 * np.pi / cut
    p_rest = centre_rest + np.c_[r_band * np.cos(th_in), r_band * np.sin(th_in), np.zeros(len(inner))]
    att = Attachment(body, p_rest, ['pelvis', 'abdomen'], [0.5, 0.5])
    pos, quat = pose0
    tgt = att.targets(pos, quat)
    c0 = tgt.mean(axis=0)
    q = quat[body.index['abdomen']]
    X0 = c0 + qrot_np(np.broadcast_to(q, (piece.n_vertices, 4)), local)
    X0[inner] = tgt
    extra = None
    if hem_weight_kg_m > 0:
        extra = np.zeros(piece.n_vertices)
        extra[hem] = hem_weight_kg_m * cut * (r_in + length) / len(hem)
    params = dict(kind=f'circle skirt cut from {cut / (2 * np.pi):g} of a circle', cut_rad=cut, length_m=length,
                  waist_radius_m=r_band, waist_circumference_m=2 * np.pi * r_band, pattern_inner_radius_m=r_in,
                  hem_circumference_m=cut * (r_in + length),
                  smooth_cone_flare_deg=float(np.degrees(np.arcsin(sa))),
                  hem_weight_kg_m=hem_weight_kg_m, hem_weight_total_kg=hem_weight_kg_m * cut * (r_in + length),
                  natural_waist_height_rest_m=z_w,
                  pattern='ring sector sewn at centre back, centre front on the straight grain',
                  attachment='waistband ring pinned to pelvis 0.5 + abdomen 0.5',
                  area_m2=piece.area(), mean_edge_m=piece.edge_length(), vertices=piece.n_vertices,
                  triangles=len(piece.F))
    return Garment('skirt', piece, X0, inner, att, params, extra)


def waist_tethers(G: Garment, slack=1.05):
    """Long-range tethers for a skirt: each vertex to the nearest waistband vertex in the flat pattern (across the
    seam too), at `slack` times their distance in the pattern, a limit on how far the cloth below can be pulled away
    from the waist (5% leaves room for a bias that drops)."""
    uv = G.piece.vertex_uv
    inner = np.asarray(G.pinned)
    cut = G.params.get('cut_rad', 2 * np.pi)
    c, s = np.cos(cut), np.sin(cut)
    images = [uv[inner], uv[inner] @ np.array([[c, s], [-s, c]]), uv[inner] @ np.array([[c, -s], [s, c]])]
    best = np.full(len(uv), np.inf)
    anchor = np.full(len(uv), -1)
    for im in images:
        d = np.linalg.norm(uv[:, None, :] - im[None, :, :], axis=2)
        k = d.argmin(1)
        dm = d[np.arange(len(uv)), k]
        better = dm < best
        best[better] = dm[better]
        anchor[better] = inner[k[better]]
    return anchor, slack * best


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
