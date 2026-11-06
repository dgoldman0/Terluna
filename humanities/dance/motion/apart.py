"""Keeping the body out of itself: the distances between the surfaces of the body's capsules that can trap cloth
between them, as a Crocoddyl residual for the optimal control and as a check on a solved motion.

The solve keeps the soles, knees and hands above the floor, but until 3 November nothing kept one part of the body
out of another: the motions of 1-2 November passed the legs through each other in strides and turns (the shins 82 mm
deep in the leap's run-up, the thighs up to 140 mm in the whirl) and, in the hanging salto's tuck, the arms through
each other and an upper arm into the head (README.md, the body through itself). The pairs kept apart are the two legs
(the thigh, shank and foot of one against those of the other), the two arms (upper arm, forearm and hand) and each
arm against the head; parts that meet at a joint, or rest against the trunk as the arms do, are left out. A barrier
on each pair's distance keeps it at or above zero: the parts may touch, as legs pressed together in fifth do, but not
pass into each other.
"""
from __future__ import annotations

import threading

import numpy as np
import pinocchio as pin
import crocoddyl

from humanities.dance.motion_io import body_spec

LEGS = ('thigh', 'shank', 'foot')
ARMS = ('upper_arm', 'forearm', 'hand')


def pairs():
    """The pairs kept apart: (group, segment, segment)."""
    out = [('legs', f'L_{a}', f'R_{b}') for a in LEGS for b in LEGS]
    out += [('arms', f'L_{a}', f'R_{b}') for a in ARMS for b in ARMS]
    out += [('arm-head', f'{s}_{a}', 'head') for s in 'LR' for a in ARMS]
    return out


def capsules(spec=None):
    """Each segment's first capsule: its end points in the segment's frame (relative to its origin, at rest) and its
    radius. The limbs and the head have one capsule each."""
    spec = spec or body_spec()
    return {s['name']: (np.array(s['capsules'][0]['a'], float), np.array(s['capsules'][0]['b'], float),
                        float(s['capsules'][0]['radius'])) for s in spec['segments'] if s.get('capsules')}


def closest(P1, Q1, P2, Q2):
    """Parameters s, t in [0, 1] of the closest points of segments P1-Q1 and P2-Q2, row by row (arrays (n, 3)), for
    segments of nonzero length (Ericson, Real-Time Collision Detection, 2005, section 5.1.9)."""
    d1, d2, r = Q1 - P1, Q2 - P2, P1 - P2
    a, e = (d1 * d1).sum(1), (d2 * d2).sum(1)
    b, c, f = (d1 * d2).sum(1), (d1 * r).sum(1), (d2 * r).sum(1)
    den = a * e - b * b
    s = np.where(den > 1e-14, np.clip((b * f - c * e) / np.where(den > 1e-14, den, 1.0), 0.0, 1.0), 0.0)
    t = (b * s + f) / e
    s = np.where(t < 0.0, np.clip(-c / a, 0.0, 1.0), np.where(t > 1.0, np.clip((b - c) / a, 0.0, 1.0), s))
    return s, np.clip(t, 0.0, 1.0)


class SelfDistance:
    """The pairs' surface distances on a Pinocchio model built from body.json (pinmodel.build), where each segment's
    frame sits on its last joint with its rest orientation the world's."""

    def __init__(self, model, spec=None):
        self.model, self.data = model, model.createData()
        self.prs = pairs()
        caps = capsules(spec)
        segs = sorted({p[1] for p in self.prs} | {p[2] for p in self.prs})
        self.seg_joint = np.array([model.frames[model.getFrameId(s)].parentJoint for s in segs])
        k = {s: i for i, s in enumerate(segs)}
        self.i1 = np.array([k[p[1]] for p in self.prs])
        self.i2 = np.array([k[p[2]] for p in self.prs])
        self.ab = np.array([[caps[s][0], caps[s][1]] for s in segs])            # (segments, 2, 3)
        self.rad = np.array([caps[s][2] for s in segs])
        self.groups = [p[0] for p in self.prs]

    def _ends(self, q):
        pin.forwardKinematics(self.model, self.data, q)
        R = np.array([self.data.oMi[int(j)].rotation for j in self.seg_joint])
        o = np.array([self.data.oMi[int(j)].translation for j in self.seg_joint])
        return np.einsum('kij,kej->kei', R, self.ab) + o[:, None, :], o            # ends (segments, 2, 3), origins

    def distances(self, q):
        """Each pair's surface distance (m; negative where the capsules overlap), the closest points and the unit
        vector from the second part's point to the first's."""
        E, o = self._ends(q)
        P1, Q1, P2, Q2 = E[self.i1, 0], E[self.i1, 1], E[self.i2, 0], E[self.i2, 1]
        s, t = closest(P1, Q1, P2, Q2)
        c1 = P1 + s[:, None] * (Q1 - P1)
        c2 = P2 + t[:, None] * (Q2 - P2)
        dv = c1 - c2
        dn = np.linalg.norm(dv, axis=1)
        mid = 0.5 * (P1 + Q1) - 0.5 * (P2 + Q2)                # where the axes cross, the midpoints' direction
        n = np.where(dn[:, None] > 1e-9, dv / np.maximum(dn, 1e-9)[:, None],
                     mid / np.maximum(np.linalg.norm(mid, axis=1), 1e-12)[:, None])
        return dn - self.rad[self.i1] - self.rad[self.i2], c1, c2, n, o

    def jacobian(self, q):
        """d(distance)/dq in the tangent space (pairs, nv): n . (J(c1) - J(c2)), each closest point's velocity Jacobian
        from its segment's joint Jacobian. The closest points' own motion along the axes adds nothing at first order."""
        d, c1, c2, n, o = self.distances(q)
        pin.computeJointJacobians(self.model, self.data, q)
        J = np.array([pin.getJointJacobian(self.model, self.data, int(j), pin.LOCAL_WORLD_ALIGNED)
                      for j in self.seg_joint])                               # (segments, 6, nv)

        def point_jac(i, c):
            r = c - o[i]
            Jl, Ja = J[i, :3], J[i, 3:]
            # v = v_o + w x r = v_o - r x w
            return Jl - np.cross(r[:, :, None], Ja, axisa=1, axisb=1, axisc=1)

        return d, np.einsum('pi,piv->pv', n, point_jac(self.i1, c1) - point_jac(self.i2, c2))


_SHARED = {}


def shared_distance(model, spec=None):
    """One SelfDistance, with its Pinocchio data, per model and thread, shared by every node's residual. With one per
    node the whole phrase's windowed solve grew by about 2 MB a node and ran out of its 6 GB at 24 of its 37 seconds
    (6 November); a residual computes everything it reads within each call, so the residuals of one thread can share
    it, and each thread has its own in case the solver evaluates nodes in parallel."""
    key = (model.name, model.nq, model.nv, model.nframes, id(spec), threading.get_ident())
    if key not in _SHARED:
        _SHARED[key] = SelfDistance(model, spec)
    return _SHARED[key]


class ResidualApart(crocoddyl.ResidualModelAbstract):
    """The pairs' surface distances as a residual, for a barrier that keeps them at or above zero
    (ActivationModelQuadraticBarrier): at the node's configuration and, with an interval T > 0 between nodes, at the
    configurations the body passed through on its way there from the node before. The solver sees the body only at its
    nodes, 25 ms apart, and a swinging shin moving at 4-6 m/s slides past a standing one within that time: in the
    leap's run-up the shins stood clear or touching at every node and passed 66 mm through each other between two (3
    November). The integrator steps the configuration by the next node's velocity (symplectic Euler, q' = q + v' T),
    and the motion file interpolates between nodes along that step, so the path from the node before is exactly q
    integrated back along v for g T, at each fraction g of `between`; an impulse at a landing changes v at the node, and
    there the path is approximate. The derivatives are exact: each distance's Jacobian through the integration's
    (pin.dIntegrate)."""

    def __init__(self, state, nu, spec=None, interval=0.0, between=(0.25, 0.5, 0.75)):
        self.m, self.spec = state.pinocchio, spec
        self.n, self.T = len(shared_distance(self.m, spec).prs), float(interval)
        self.back = tuple(between) if self.T > 0 else ()
        crocoddyl.ResidualModelAbstract.__init__(self, state, self.n * (1 + len(self.back)), nu, True,
                                                 bool(self.back), False)
        self.nq, self.nv = self.m.nq, self.m.nv

    @property
    def sd(self):
        return shared_distance(self.m, self.spec)

    def calc(self, data, x, u=None):
        q, v = x[:self.nq], x[self.nq:]
        sd = self.sd
        data.r[:self.n] = sd.distances(q)[0]
        for k, g in enumerate(self.back, 1):
            data.r[k * self.n:(k + 1) * self.n] = sd.distances(pin.integrate(self.m, q, -g * self.T * v))[0]

    def calcDiff(self, data, x, u=None):
        q, v = x[:self.nq], x[self.nq:]
        sd = self.sd
        data.Rx[:, :] = 0.0
        data.Rx[:self.n, :self.nv] = sd.jacobian(q)[1]
        for k, g in enumerate(self.back, 1):
            w = -g * self.T * v
            J = sd.jacobian(pin.integrate(self.m, q, w))[1]
            rows = slice(k * self.n, (k + 1) * self.n)
            data.Rx[rows, :self.nv] = J @ pin.dIntegrate(self.m, q, w, pin.ArgumentPosition.ARG0)
            data.Rx[rows, self.nv:] = (-g * self.T) * (J @ pin.dIntegrate(self.m, q, w, pin.ArgumentPosition.ARG1))

    def createData(self, collector):
        return crocoddyl.ResidualDataAbstract(self, collector)


def apart_cost(state, nu, weight, spec=None, interval=0.0, between=(0.25, 0.5, 0.75), margin=0.0):
    """The barrier cost: zero while every pair is at or above `margin` (and, with an interval between nodes, on the way
    from the node before, at each fraction of `between`), quadratic in the shortfall below it."""
    res = ResidualApart(state, nu, spec, interval, between)
    n = res.nr
    act = crocoddyl.ActivationModelQuadraticBarrier(crocoddyl.ActivationBounds(np.full(n, float(margin)),
                                                                                np.full(n, 1e3)))
    return crocoddyl.CostModelResidual(state, act, res), weight


def scan(model, qs, spec=None):
    """Each group's least surface distance over a motion's configurations qs (n, nq): {group: (distance m, index)}."""
    sd = SelfDistance(model, spec)
    out = {}
    for k, q in enumerate(qs):
        d = sd.distances(q)[0]
        for g in set(sd.groups):
            dm = float(d[[i for i, gg in enumerate(sd.groups) if gg == g]].min())
            if g not in out or dm < out[g][0]:
                out[g] = (dm, k)
    return out
