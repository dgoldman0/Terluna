"""Convex impulse refinement within a three-mode departure control family.

Each separating half-space is local to a selected side of a finite-square
encounter. Successful optimization is not a general optimum or certificate.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from .active_formation import solar_frame
from .pattern_departure import layer_labels, expansion_state
from .optical import length


def control_modes(env, start, rows=19):
    iy, ix = np.indices((rows, rows)); mid = (rows-1)/2
    u, b, c = solar_frame(env.at(start)).T
    return np.stack([layer_labels(rows)[:, None]*u,
        ((ix.ravel()-mid)/mid)[:, None]*b,
        ((iy.ravel()-mid)/mid)[:, None]*c], axis=-1)


def mode_impulses(modes, parameters):
    return np.einsum('nij,j->ni', modes, parameters)


def local_cuts(maps, initial, modes, command, times, parameters, margin=400., per_date=12):
    impulses = mode_impulses(modes, parameters)
    states = expansion_state(maps, initial, impulses, times)
    cuts = []
    for t, state, frame in zip(times, states, command(times).as_matrix()):
        pairs = cKDTree(state[:, :3]).query_pairs(14500., output_type='ndarray')
        if not len(pairs): continue
        projected = (state[pairs[:, 0], :3]-state[pairs[:, 1], :3])@frame
        separation = abs(projected)-[10000., 10000., 0.]
        scores = separation.max(axis=1)
        response = maps.sol(t)[42:].reshape(6, 3)[:3]
        for k in np.argsort(scores)[:per_date]:
            i, j = pairs[k]
            gradients = frame.T@response@(modes[i]-modes[j])
            cuts.append(choose_cut(t, i, j, projected[k], gradients, parameters, margin))
    return [c for c in cuts if c is not None]


def choose_cut(t, i, j, projected, gradients, parameters, margin, in_plane=False):
    choices = []
    for axis in range(2 if in_plane else 3):
        extent = 10000. if axis < 2 else 0.
        for sign in [-1., 1.]:
            gradient = sign*gradients[axis]
            current = sign*projected[axis]-extent
            norm = np.linalg.norm(gradient)
            if norm < 1e-6:
                continue
            # Prefer an already safe half-space. Otherwise choose the smallest
            # local parameter displacement; the subsequent cap may reject it.
            distance = (margin-current)/norm
            choices.append((distance, axis, sign, gradient, margin-current+gradient@parameters))
    if not choices:
        return None
    _, axis, sign, gradient, rhs = min(choices, key=lambda c: c[0])
    return dict(time_s=float(t), i=int(i), j=int(j), axis=axis, sign=sign,
        gradient=gradient.tolist(), rhs=float(rhs))


def crossing_cuts(maps, modes, command, parameters, witnesses, margin=400.):
    cuts = []
    for w in witnesses:
        t, i, j = w['time_s'], w['i'], w['j']
        frame = command(t).as_matrix()
        response = maps.sol(t)[42:].reshape(6, 3)[:3]
        gradient = frame.T@response@(modes[i]-modes[j])
        projected = np.r_[w['in_plane_projections_m'], 0.]
        cuts.append(choose_cut(t, i, j, projected, gradient, parameters, margin, in_plane=True))
    return [c for c in cuts if c is not None]


def minimize_impulse(modes, start, cuts, burn, cap=.001):
    a = np.array([c['gradient'] for c in cuts])/10000.
    b = np.array([c['rhs'] for c in cuts])/10000.
    def objective(x):
        impulses = mode_impulses(modes, x)
        norm = np.maximum(length(impulses), 1e-12)
        return float(norm.mean()), np.einsum('ni,nij->j', impulses/norm[:, None], modes)/len(modes)
    # One smooth common arc: the largest individual integrated vector sets
    # the exact peak-acceleration cap, not a component-box approximation.
    squared_modes = np.sum(modes*modes, axis=1)
    constraints = [dict(type='ineq', fun=lambda x: a@x-b, jac=lambda x: a),
        dict(type='ineq', fun=lambda x: (cap*burn/2)**2-squared_modes@(x*x),
             jac=lambda x: -2*squared_modes*x)]
    fit = minimize(lambda x: objective(x), start, jac=True, method='SLSQP',
        bounds=[(0., 3.), (0., 4.), (0., 4.)], constraints=constraints,
        options={'ftol': 1e-10, 'maxiter': 150})
    actual_margin = float((a@fit.x-b).min()) if len(cuts) else np.inf
    peak = float(length(mode_impulses(modes, fit.x)).max()*2/burn)
    success = bool(fit.success and actual_margin >= -1e-7 and peak <= cap+1e-10)
    return fit.x, dict(success=success, solver_success=bool(fit.success), status=int(fit.status),
        message=str(fit.message), iterations=int(fit.nit), cuts=len(cuts),
        parameters_m_s=fit.x.tolist(), mean_delta_v_m_s=objective(fit.x)[0],
        maximum_acceleration_m_s2=peak, minimum_cut_residual_m=actual_margin*10000.)
