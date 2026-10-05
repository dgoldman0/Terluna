"""Depth arrangements and sparse individual two-burn encounter corrections.

Gravity response is a proposal model. Measured sail/shadow defects are added
between coupled replays; an LP solution alone never accepts a trajectory.
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, eye, hstack, vstack
from scipy.spatial import cKDTree

from .active_formation import lattice
from .pattern_return import smooth_arc
from .optical import length


def packing_offsets(depth, swapped=False, rows=19, pitch=8500.):
    offsets = lattice(rows, pitch, depth)[0]
    if swapped:
        iy, ix = np.indices((rows, rows))
        offsets[:, 0] = ((2*(ix % 2)+(iy % 2)).ravel()-1.5)*depth
    return offsets


def response_maps(env, reference, start, end, windows, basis):
    """Common gravity STM and responses to independently weighted smooth arcs."""
    arcs = len(windows)
    def rhs(t, flat):
        matrix = flat.reshape(6, 6+3*arcs)
        q = reference(t)[:3]; eps = 10.
        g = env.gravity(t, q+np.r_[np.eye(3), -np.eye(3)]*eps)
        gradient = ((g[:3]-g[3:])/(2*eps)).T
        out = np.empty_like(matrix)
        out[:3] = matrix[3:]; out[3:] = gradient@matrix[:3]
        for k, scale in enumerate(smooth_arc(t, windows)):
            out[3:, 6+3*k:9+3*k] += scale*basis
        return out.ravel()
    initial = np.c_[np.eye(6), np.zeros((6, 3*arcs))]
    sol = solve_ivp(rhs, [start, end], initial.ravel(), method='DOP853',
        rtol=2e-11, atol=1e-9, max_step=120., dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    return sol


def matrices(maps, times, arcs=2):
    return maps.sol(np.atleast_1d(times)).T.reshape(-1, 6, 6+3*arcs)


def controls_state(base, maps, times, controls, original):
    times = np.atleast_1d(times)
    b = matrices(maps, times, controls.shape[1])[:, :, 6:]
    return base(times)+np.einsum('tij,nj->tni', b, (controls-original).reshape(len(controls), -1))


def measured_base(previous, maps, times, states):
    """Add the measured state defect, transporting its last value over the tail."""
    residual = states-previous(times)
    defect = CubicHermiteSpline(times, residual[:, :, :3], residual[:, :, 3:])
    last = float(times[-1]); inverse = np.linalg.inv(matrices(maps, [last])[0, :, :6])
    def evaluate(t):
        t = np.atleast_1d(t); out = previous(t).copy()
        before = t <= last
        out[before] += np.concatenate([defect(t[before]), defect(t[before], 1)], axis=-1)
        if np.any(~before):
            phi = matrices(maps, t[~before])[:, :, :6]
            out[~before] += np.einsum('tij,nj->tni', phi@inverse, residual[-1])
        return out
    return evaluate, float(length(residual[:, :, :3]).max())


def separation_cuts(times, states, command, maps, controls, margin=400., per_date=48):
    """Choose local finite-square branches; constrain closest pairs first."""
    cuts = []
    for t, state, frame, matrix in zip(times, states, command(times).as_matrix(), matrices(maps, times)):
        pairs = cKDTree(state[:, :3]).query_pairs(15000., output_type='ndarray')
        if not len(pairs): continue
        projected = (state[pairs[:, 0], :3]-state[pairs[:, 1], :3])@frame
        clear = abs(projected)-[10000., 10000., 0.]
        scores = clear.max(axis=1)
        for k in np.argsort(scores)[:per_date]:
            i, j = pairs[k]; axis = int(np.argmax(clear[k]))
            sign = 1. if projected[k, axis] >= 0 else -1.
            gradient = sign*frame[:, axis]@matrix[:3, 6:]
            rhs = margin-clear[k, axis]+gradient@(controls[i]-controls[j]).ravel()
            if np.linalg.norm(gradient) < 1e-8 and rhs <= 0: continue
            cuts.append(dict(t=float(t), i=int(i), j=int(j), axis=axis, sign=sign,
                gradient=gradient, rhs=float(rhs)))
    return cuts


def crossing_cuts(witnesses, command, maps, controls, margin=400.):
    cuts = []
    for w in witnesses:
        t, i, j = w['time_s'], w['i'], w['j']
        d = np.array(w['in_plane_projections_m']); axis = int(np.argmax(abs(d)))
        sign = 1. if d[axis] >= 0 else -1.
        gradient = sign*command(t).as_matrix()[:, axis]@matrices(maps, [t])[0, :3, 6:]
        rhs = margin+10000.-abs(d[axis])+gradient@(controls[i]-controls[j]).ravel()
        cuts.append(dict(t=float(t), i=int(i), j=int(j), axis=axis, sign=sign,
            gradient=gradient, rhs=float(rhs)))
    return cuts


def independent_impulses(controls, cuts, windows, trust=.75, cap=.001, wall_s=25.):
    """Sparse L1 fuel surrogate, with inner component bounds on vector thrust."""
    count, arcs, _ = controls.shape; width = arcs*3; variables = count*width
    row = []; col = []; data = []; rhs = []
    for n, cut in enumerate(cuts):
        for member, sign in [(cut['i'], -1.), (cut['j'], 1.)]:
            row.extend([n]*width); col.extend(member*width+np.arange(width))
            data.extend(sign*cut['gradient']/10000.)
        rhs.append(-cut['rhs']/10000.)
    geometry = coo_matrix((data, (row, col)), shape=(len(cuts), variables)).tocsr()
    identity = eye(variables, format='csr')
    a = vstack([hstack([geometry, coo_matrix(geometry.shape)]),
        hstack([identity, -identity]), hstack([-identity, -identity])], format='csr')
    b = np.r_[rhs, np.zeros(2*variables)]
    limits = np.repeat(cap*np.diff(windows, axis=1).ravel()/(2*np.sqrt(3)), 3)
    limits = np.tile(limits, count)
    flat = controls.ravel()
    bounds = list(zip(np.maximum(-limits, flat-trust), np.minimum(limits, flat+trust)))
    bounds += [(0., None)]*variables
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected')
        fit = linprog(np.r_[np.zeros(variables), np.ones(variables)/count],
            A_ub=a, b_ub=b, bounds=bounds, method='highs',
            options={'threads':1, 'parallel':False, 'time_limit':wall_s})
    result = dict(success=bool(fit.success), status=int(fit.status), message=fit.message,
        cuts=len(cuts), independent_vector_coefficients=variables,
        component_trust_m_s=trust, maximum_solver_wall_s=wall_s)
    if not fit.success: return controls.copy(), result
    u = fit.x[:variables].reshape(count, arcs, 3)
    peak = float(np.max(length(u)*2/np.diff(windows, axis=1).ravel()))
    result.update(component_objective_m_s=float(fit.fun),
        mean_delta_v_m_s=float(length(u).sum(axis=1).mean()),
        maximum_tile_delta_v_m_s=float(length(u).sum(axis=1).max()),
        maximum_thrust_acceleration_m_s2=peak,
        minimum_cut_residual_m=float(np.min(np.array(rhs)-geometry@u.ravel())*10000.) if cuts else None,
        burn_count=int(np.count_nonzero(length(u)>1e-8)),
        maximum_control_change_m_s=float(length(u-controls).max()))
    result['success'] = bool(peak <= cap+1e-10)
    return u, result
