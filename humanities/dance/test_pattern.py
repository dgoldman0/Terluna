"""A kinematic test pattern for the renderer and the garment solver: poses through the body's range, no physics.

    python -m humanities.dance.test_pattern OUT.npz      # needs MuJoCo (the dance environment)

Twelve seconds at 60 frames a second: the arms through bras bas, first, second and fifth; a plié and a rise to
demi-pointe; a développé à la seconde and an arabesque line; a slow whole-turn about the vertical axis that a skirt
can flare in; and a back bend. The root is lowered each frame so the lowest sole point rests on the floor. It exercises
skinning at extreme poses and garment attachment; it is a test pattern, with no claim to dance or to dynamics.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from humanities.dance.motion_io import body_spec, save_motion, mujoco_frames

HERE = Path(__file__).resolve().parent
FPS = 60
DURATION = 12.0


def smooth(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * x * (10 - 15 * x + 6 * x * x)


def bump(t, t0, t1, t2, t3):
    """0 before t0, rises to 1 by t1, holds, falls back to 0 between t2 and t3."""
    return smooth((t - t0) / (t1 - t0)) * (1 - smooth((t - t2) / (t3 - t2)))


def quat_axis_angle(axis, angle):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    return np.r_[math.cos(angle / 2), math.sin(angle / 2) * a]


def quat_mul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def pose(t):
    """Joint targets at time t: hinge angles (radians) and ball-joint quaternions."""
    d = math.radians
    h, b = {}, {}
    arms_first = bump(t, 0.5, 1.5, 2.0, 2.6)        # rounded forward
    arms_second = bump(t, 2.0, 2.6, 3.4, 4.0)       # open to the side
    arms_fifth = bump(t, 3.4, 4.2, 5.4, 6.2)        # overhead
    plie = bump(t, 1.0, 1.8, 2.2, 3.0)
    releve = bump(t, 3.2, 3.8, 4.6, 5.2)
    devel = bump(t, 5.2, 6.2, 6.8, 7.6)             # right leg à la seconde
    arab = bump(t, 7.6, 8.4, 8.8, 9.6)              # left leg behind, torso forward
    arch = bump(t, 10.0, 10.8, 11.2, 11.9)
    for s, L in ((1, 'L'), (-1, 'R')):
        # Shoulders: forward flexion about -y, abduction about s*x; elbows soften.
        q = quat_axis_angle((0, -1, 0), d(55) * arms_first + d(165) * arms_fifth * 0.0)
        q = quat_mul(quat_axis_angle((s, 0, 0), d(80) * arms_second + d(170) * arms_fifth), q)
        q = quat_mul(q, quat_axis_angle((0, 0, s), d(25) * arms_first))
        b[f'{L}_shoulder'] = q
        h[f'{L}_elbow'] = d(15) + d(40) * arms_first + d(20) * arms_fifth
        h[f'{L}_wrist_flex'] = d(10) * arms_second
        h[f'{L}_clav_elev'] = d(12) * arms_fifth
        # Legs: plié turned out, relevé onto demi-pointe.
        turnout = d(45)
        knee = d(70) * plie
        hq = quat_mul(quat_axis_angle((0, 0, s), turnout), quat_axis_angle((0, -1, 0), d(35) * plie))
        b[f'{L}_hip'] = hq
        h[f'{L}_knee'] = knee
        h[f'{L}_ankle'] = -d(38) * plie + d(55) * releve
        h[f'{L}_mtp'] = d(60) * releve
    # Développé à la seconde, right leg; arabesque, left leg behind.
    b['R_hip'] = quat_mul(quat_axis_angle((-1, 0, 0), d(95) * devel), b['R_hip'])
    h['R_ankle'] = h['R_ankle'] + d(70) * devel
    b['L_hip'] = quat_mul(quat_axis_angle((0, 1, 0), d(80) * arab), b['L_hip'])
    h['L_ankle'] = h['L_ankle'] + d(70) * arab
    h['lumbar_flex'] = d(20) * arab - d(25) * arch
    h['thoracic_flex'] = d(10) * arab - d(25) * arch
    h['neck_flex'] = -d(30) * arch
    h['thoracic_rot'] = d(20) * arms_second
    return h, b


def main(out) -> int:
    import mujoco
    spec = body_spec()
    model = mujoco.MjModel.from_xml_path(str(HERE / 'body' / 'dancer.xml'))
    data = mujoco.MjData(model)
    names = [s['name'] for s in spec['segments']]
    sole = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, n) for n in spec['contacts'] if 'palm' not in n]
    times = np.arange(0, DURATION + 1e-9, 1 / FPS)
    P, Q, C, QP = [], [], [], []
    for t in times:
        data.qpos[:] = model.qpos0
        h, b = pose(t)
        for n, v in h.items():
            j = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, n)
            data.qpos[model.jnt_qposadr[j]] = v
        for n, v in b.items():
            j = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, n)
            a = model.jnt_qposadr[j]
            data.qpos[a:a + 4] = v / np.linalg.norm(v)
        spin = 2 * math.pi * smooth((t - 6.0) / 4.0)                 # one slow turn between 6 and 10 s
        data.qpos[3:7] = quat_axis_angle((0, 0, 1), spin)
        mujoco.mj_kinematics(model, data)
        low = min(data.site_xpos[i][2] for i in sole)
        data.qpos[2] -= low
        mujoco.mj_kinematics(model, data)
        mujoco.mj_comPos(model, data)
        p, q = mujoco_frames(model, data, names)
        P.append(p)
        Q.append(q)
        C.append(data.subtree_com[1].copy())
        QP.append(data.qpos.copy())
    save_motion(out, times, P, Q, C, QP, meta=dict(
        kind='test pattern', gravity=None, air_density=None,
        evidence='Kinematic poses for renderer and garment development; no dynamics, not a dance.',
        producer='humanities/dance/test_pattern.py'))
    print(f'wrote {out}: {len(times)} frames')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else 'test_pattern.npz'))
