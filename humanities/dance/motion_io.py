"""The motion file every dance consumer reads: segment poses through time, with the physics that made them.

A motion is an .npz archive (schema terluna.humanities.dance-motion/1):

    time          (N,)        seconds from the phrase's start, uniform at `fps`
    fps           ()          frames per second
    segments      (S,)        segment names, in body.json's order
    seg_pos       (N, S, 3)   world position of each segment frame (its proximal joint centre), metres
    seg_quat      (N, S, 4)   world orientation of each segment frame, (w, x, y, z); identity in the rest pose
    com           (N, 3)      whole-body centre of mass
    contact_names (C,)        the sole and palm points of body.json
    contact_on    (N, C)      the point bears load
    contact_force (N, C, 3)   the floor's force on each point, newtons, world frame
    q             (N, nq)     generalised coordinates in MuJoCo's order (free root, then joints), for replay
    meta          ()          JSON text: gravity, air density, score, producer and evidence

Segment frames are world-aligned at rest (body.json), so a renderer drives each bone by its segment's world rotation
and position alone, and a garment solver reads its attachment and collision frames the same way. Joint conventions
stay inside the body model.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

SCHEMA = 'terluna.humanities.dance-motion/1'
BODY = Path(__file__).resolve().parent / 'body' / 'body.json'


def body_spec() -> dict:
    return json.loads(BODY.read_text())


def save_motion(path, time, seg_pos, seg_quat, com, q, meta: dict, contact_on=None, contact_force=None):
    spec = body_spec()
    names = [s['name'] for s in spec['segments']]
    cnames = list(spec['contacts'])
    n = len(time)
    if contact_on is None:
        contact_on = np.zeros((n, len(cnames)), bool)
    if contact_force is None:
        contact_force = np.zeros((n, len(cnames), 3))
    fps = 1.0 / float(time[1] - time[0]) if n > 1 else 0.0
    np.savez_compressed(path, schema=SCHEMA, time=np.asarray(time, float), fps=fps, segments=np.array(names),
                        seg_pos=np.asarray(seg_pos, float), seg_quat=np.asarray(seg_quat, float),
                        com=np.asarray(com, float), contact_names=np.array(cnames),
                        contact_on=np.asarray(contact_on, bool), contact_force=np.asarray(contact_force, float),
                        q=np.asarray(q, float), meta=json.dumps(meta))


def load_motion(path) -> dict:
    with np.load(path, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files}
    assert str(out['schema']) == SCHEMA, out['schema']
    out['meta'] = json.loads(str(out['meta']))
    return out


def mujoco_frames(model, data, segment_names):
    """Segment positions and orientations from a MuJoCo state whose bodies carry the segment names."""
    import mujoco
    ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, n) for n in segment_names]
    return np.array([data.xpos[i] for i in ids]), np.array([data.xquat[i] for i in ids])
