"""Ordinary dance garments cut from flat patterns: a full circle skirt and a pair of Chinese water sleeves.

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
from .patterns import Piece, annulus, rectangle

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
