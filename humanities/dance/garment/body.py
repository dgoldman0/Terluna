"""The body as the garment sees it: capsules and attachment frames from body.json, posed from a motion each substep.

Segment frames are world-aligned at rest (body.json), so a point fixed to a segment moves as
p(t) = seg_pos(t) + R(seg_quat(t)) (p_rest - origin_rest). Poses between frames are interpolated linearly in
position and by normalised linear interpolation in orientation. Before a motion starts a garment is dressed: the
dancer holds the first frame's pose with the arms raised to the side, the garment settles, and the arms come down to
the first frame, so that the cloth ends outside the arms as it is worn.
"""
from __future__ import annotations

import math

import numpy as np
import torch

from humanities.dance.motion_io import body_spec, load_motion
from .cloth import quat_rotate


def qmul(p, q):
    w1, x1, y1, z1 = np.moveaxis(p, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(q, -1, 0)
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], -1)


def qaxis(axis, ang):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    return np.r_[math.cos(ang / 2), math.sin(ang / 2) * a]


def qrot_np(q, v):
    w, u = q[..., :1], q[..., 1:]
    t = 2 * np.cross(u, v)
    return v + w * t + np.cross(u, t)


class Body:
    def __init__(self, spec=None, form=None):
        """form: a dress form (schema terluna.visualization.dance-dress-form/1) whose tapered pieces replace
        body.json's capsules for the cloth, so that dressed renders' cloth stays outside the skin they show; its
        waist outline places the skirts' waistband. Without it the cloth collides with body.json's capsules."""
        sp = spec or body_spec()
        self.spec = sp
        self.names = [s['name'] for s in sp['segments']]
        self.index = {n: i for i, n in enumerate(self.names)}
        self.origin = np.array([s['origin'] for s in sp['segments']], float)       # world, rest
        caps = []
        if form is not None:
            assert form['schema'] == 'terluna.visualization.dance-dress-form/1', form['schema']
            for p in form['pieces']:
                caps.append((self.index[p['segment']], p['a'], p['b'], p['ra'], p['segment'], p['rb']))
        else:
            for i, s in enumerate(sp['segments']):
                for c in s['capsules']:
                    caps.append((i, c['a'], c['b'], c['radius'], s['name'], c['radius']))
        self.form = form
        self.waist = form['waist'] if form is not None else None
        self.cap_seg = np.array([c[0] for c in caps])
        self.cap_a = np.array([c[1] for c in caps], float)                           # relative to the segment origin
        self.cap_b = np.array([c[2] for c in caps], float)
        self.cap_r = np.array([c[3] for c in caps], float)                           # radius at a
        self.cap_rb = np.array([c[5] for c in caps], float)                          # radius at b (tapered pieces)
        self.cap_name = [c[4] for c in caps]

    def capsules(self, pos, quat):
        """World capsule end points for one pose (numpy)."""
        s = self.cap_seg
        return pos[s] + qrot_np(quat[s], self.cap_a), pos[s] + qrot_np(quat[s], self.cap_b), self.cap_r

    def arms_raised(self, pos, quat, angle):
        """The pose with each arm (upper arm, forearm, hand) turned about its shoulder by `angle` of abduction."""
        pos, quat = pos.copy(), quat.copy()
        forward = qrot_np(quat[self.index['thorax']], np.array([1.0, 0.0, 0.0]))    # abduct about the body's own axis
        for side, sgn in (('L', 1.0), ('R', -1.0)):
            q = qaxis(sgn * forward, angle)
            sh = pos[self.index[f'{side}_upper_arm']].copy()
            for seg in ('upper_arm', 'forearm', 'hand'):
                k = self.index[f'{side}_{seg}']
                pos[k] = sh + qrot_np(q, pos[k] - sh)
                quat[k] = qmul(q, quat[k])
        return pos, quat


class Attachment:
    """Pinned vertices fixed to one or two segments: target = sum_s w_s (pos_s + R_s (p_rest - origin_s))."""

    def __init__(self, body: Body, p_rest: np.ndarray, segments, weights):
        self.body = body
        self.p_rest = np.asarray(p_rest, float)
        self.seg = [body.index[s] for s in segments]
        self.wts = list(weights)
        self.off = [self.p_rest - body.origin[k] for k in self.seg]

    def targets(self, pos, quat):
        out = 0.0
        for k, w, off in zip(self.seg, self.wts, self.off):
            out = out + w * (pos[k] + qrot_np(quat[k], off))
        return out


class MotionPose:
    """Segment poses at any time from a motion file; times before 0 run the dressing sequence on the first frame."""

    def __init__(self, path_or_motion, body: Body, raise_deg=80.0, t_raised=3.0, t_lower=1.5, t_hold=1.5):
        m = load_motion(path_or_motion) if isinstance(path_or_motion, (str, bytes)) or hasattr(path_or_motion, 'parent') else path_or_motion
        self.m = m
        names = list(m['segments'])
        order = [names.index(n) for n in body.names]
        self.pos = np.asarray(m['seg_pos'], float)[:, order]
        self.quat = np.asarray(m['seg_quat'], float)[:, order]
        for k in range(1, len(self.quat)):              # keep successive quaternions in one hemisphere
            flip = (self.quat[k] * self.quat[k - 1]).sum(-1) < 0
            self.quat[k][flip] *= -1
        self.fps = float(m['fps'])
        self.n = len(self.pos)
        self.duration = (self.n - 1) / self.fps
        self.body = body
        self.raise_rad = math.radians(raise_deg)
        self.t_raised, self.t_lower, self.t_hold = t_raised, t_lower, t_hold
        self.pre = t_raised + t_lower + t_hold

    def __call__(self, t):
        if t < 0:
            tau = t + self.pre                       # 0 at the start of dressing
            if tau < self.t_raised:
                f = 1.0
            elif tau < self.t_raised + self.t_lower:
                x = (tau - self.t_raised) / self.t_lower
                f = 1.0 - x * x * x * (10 - 15 * x + 6 * x * x)
            else:
                f = 0.0
            return self.body.arms_raised(self.pos[0], self.quat[0], f * self.raise_rad)
        u = min(t * self.fps, self.n - 1 - 1e-9)
        i = int(math.floor(u))
        a = u - i
        p = (1 - a) * self.pos[i] + a * self.pos[i + 1]
        q = (1 - a) * self.quat[i] + a * self.quat[i + 1]
        q = q / np.linalg.norm(q, axis=-1, keepdims=True)
        return p, q


class SpinPose:
    """A fixed pose turned about the vertical through `axis_xy` at the rate profile omega(t) (rad/s), with a smooth
    spin-up over t_up seconds; used for the steady-turn validation."""

    def __init__(self, pos, quat, omega, t_up=2.0, axis_xy=(0.0, 0.0)):
        self.pos0, self.quat0 = pos, quat
        self.omega, self.t_up = omega, t_up
        self.c = np.array([axis_xy[0], axis_xy[1], 0.0])

    def angle(self, t):
        if t <= 0:
            return 0.0
        if t < self.t_up:                                   # rate rises as smoothstep; angle is its integral
            x = t / self.t_up
            return self.omega * self.t_up * (2.5 * x ** 4 - 3 * x ** 5 + x ** 6)
        return self.omega * self.t_up * 0.5 + self.omega * (t - self.t_up)

    def __call__(self, t):
        q = qaxis((0, 0, 1), self.angle(t))
        pos = self.c + qrot_np(np.broadcast_to(q, (len(self.pos0), 4)), self.pos0 - self.c)
        quat = qmul(np.broadcast_to(q, self.quat0.shape), self.quat0)
        return pos, quat


def run(cloth, pose, attachment, body: Body, t0: float, t1: float, n_sub: int, fps: float, record=None,
        damping=None, callback=None):
    """Advance `cloth` from t0 to t1 with n_sub substeps per 1/fps. The body is posed at frame times on the host;
    capsule end points and pinned targets are interpolated linearly on the device within each frame.
    record(t) -> bool picks the frame times whose vertices are kept; damping(t) -> 1/s sets settling damping."""
    dt = 1.0 / (fps * n_sub)
    dev = cloth.dev
    dtype = cloth.x.dtype
    tt = lambda a: torch.as_tensor(np.asarray(a), dtype=dtype, device=dev)

    def frame(t):
        p, q = pose(t)
        a, b, r = body.capsules(p, q)
        tg = tt(attachment.targets(p, q)) if attachment is not None else None
        return tt(a), tt(b), tg

    if not cloth.has_body:
        cloth.set_body(body.cap_r, body.cap_rb)
    n_frames = int(round((t1 - t0) * fps))
    A0, B0, T0 = frame(t0)
    frames, times = [], []
    for f in range(n_frames):
        t = t0 + (f + 1) / fps
        A1, B1, T1 = frame(t)
        cloth.set_frame(A0, A1, B0, B1, T0, T1, damping=damping(t) if damping is not None else None)
        for k in range(n_sub):
            cloth.substep(dt, k / n_sub, (k + 1) / n_sub)
        A0, B0, T0 = A1, B1, T1
        if record is not None and record(t):
            frames.append(cloth.x.detach().to('cpu', torch.float32).numpy().copy())
            times.append(t)
        if callback is not None:
            callback(t, cloth)
    return (np.array(frames, np.float32) if frames else None), np.array(times)
