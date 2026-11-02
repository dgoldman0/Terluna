"""The dancer in a named shape: segment poses by forward kinematics over the body model, and the whole body's inertia
about its centre of mass.

A shape is the vocabulary's combination of regional poses (poses.json: legs, arms, torso, head, feet) for a gesture
side G, laid over the carriage the motion engine uses between shapes. This module places the joints the way the
engine does (a ball joint swings its limb's long axis toward an elevation and azimuth, then twists about it; a hinge
turns about its axis, held within its range), in plain numpy, so that the physics screen (envelope.py) can ask how a
shape turns in the air without the engine's libraries. The engine's tests check the two agree.

Frames: x forward, y to the dancer's left, z up, with the pelvis upright at the origin. A somersault turns about y, a
cartwheel about x and a twist about z.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from humanities.dance.body.model import spec as body_spec

POSES = Path(__file__).resolve().parents[1] / 'vocabulary' / 'poses.json'
CARRIAGE = {'legs': 'turned_out', 'arms': 'bras_bas', 'torso': 'upright', 'head': 'forward'}
# Spine and neck bends and rotations are written toward G; the model's positive sense is toward the left.
SIDED_SIGN = {'lumbar_bend', 'lumbar_rot', 'thoracic_bend', 'thoracic_rot', 'neck_bend', 'neck_rot'}


def regions():
    return json.loads(POSES.read_text())['regions']


def axis_angle(axis, angle):
    """The rotation matrix about a unit axis (Rodrigues)."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * K @ K


def ball_matrix(spec, side_sign):
    """A ball joint's rotation from elev/az/twist: swing the long axis (down at rest) toward elevation elev and azimuth
    az (outward is the limb's own side), then twist about the long axis."""
    elev = math.radians(spec.get('elev', 0.0))
    az = math.radians(spec.get('az', 0.0))
    twist = math.radians(spec.get('twist', 0.0))
    d = np.array([math.cos(az), side_sign * math.sin(az), 0.0])
    axis = np.cross([0.0, 0.0, -1.0], d)
    swing = axis_angle(axis, elev) if np.linalg.norm(axis) > 1e-9 else np.eye(3)
    return swing @ axis_angle([0.0, 0.0, side_sign], twist)


def joint_values(shape, side='R', carriage=True):
    """{joint name: ('hinge', radians) or ('ball', 3x3)} for a shape ({region: pose name}) with gesture side `side`;
    mirrored shapes swap G and S. The carriage fills the regions the shape leaves out, as in the motion engine."""
    G, S = side, {'L': 'R', 'R': 'L'}[side]
    if shape.get('mirror'):
        G, S = S, G
    full = {**CARRIAGE, **shape} if carriage else dict(shape)
    lib = regions()
    out = {}
    for region in ('legs', 'arms', 'torso', 'head', 'feet'):
        name = full.get(region)
        if not name:
            continue
        for jname, val in lib[region][name]['joints'].items():
            absolute = jname.replace('G_', f'{G}_').replace('S_', f'{S}_')
            limb = absolute[0] if absolute[:2] in ('L_', 'R_') else None
            if isinstance(val, dict):
                out[absolute] = ('ball', ball_matrix(val, 1 if limb == 'L' else -1))
            else:
                ang = math.radians(val)
                if jname in SIDED_SIGN and G == 'R':
                    ang = -ang
                out[absolute] = ('hinge', ang)
    return out


def frames(values, sp=None):
    """Each segment's world rotation and origin with the pelvis upright at the origin, the joints set to `values`."""
    sp = sp or body_spec()
    children = {}
    for s in sp['segments']:
        children.setdefault(s['parent'], []).append(s)
    out = {}

    def visit(seg, R_parent, p_parent):
        if seg['parent'] is None:
            R, p = np.eye(3), np.zeros(3)
        else:
            p = p_parent + R_parent @ np.array(seg['offset'], float)
            R = R_parent.copy()
            for j in seg['joints']:
                kind, val = values.get(j['name'], (None, None))
                if j['kind'] == 'hinge':
                    lo, hi = (math.radians(x) for x in j['range_deg'])
                    ang = min(max(val, lo), hi) if kind == 'hinge' else 0.0
                    R = R @ axis_angle(j['axis'], ang)
                elif kind == 'ball':
                    R = R @ val
        out[seg['name']] = (R, p)
        for ch in children.get(seg['name'], []):
            visit(ch, R, p)

    for root in children[None]:
        visit(root, None, None)
    return out


def inertia(shape, side='R', carriage=True):
    """The whole body's mass, centre of mass and inertia tensor about it (kg m^2, the pelvis's frame) in a shape."""
    sp = body_spec()
    fr = frames(joint_values(shape, side, carriage), sp)
    items = []
    for s in sp['segments']:
        R, p = fr[s['name']]
        c = p + R @ np.array(s['com'], float)
        I = R @ np.diag(s['inertia_kg_m2']) @ R.T
        items.append((s['mass_kg'], c, I))
    M = sum(m for m, _, _ in items)
    com = sum(m * c for m, c, _ in items) / M
    total = np.zeros((3, 3))
    for m, c, I in items:
        d = c - com
        total += I + m * (np.dot(d, d) * np.eye(3) - np.outer(d, d))
    return M, com, total
