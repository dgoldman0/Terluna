"""Small coupled-response trust-region SQP steps with explicit restoration.

The caller supplies derivatives measured with the full coupled propagator.
Finite-square half-spaces are frozen only within a step, then rebuilt. The
SQP is a local subproblem; slacks, rejected steps and cycle debt remain data.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree


def contact_set(states, frames, reserve=1500.):
    """All potentially close pairs at sampled dates, preserving approach side.

    Pick a separating face from the last pre-encounter state. Keeping that
    branch prevents a post-crossing positive unsigned distance hiding contact.
    This is a local topological choice, not a global nonintersection model.
    """
    cuts = []; approach_sign = {}
    for k, (y, frame) in enumerate(zip(states, frames)):
        pairs = cKDTree(y[:, :3]).query_pairs(np.sqrt(2)*10000+reserve, output_type='ndarray')
        if not len(pairs):
            continue
        delta = (y[pairs[:, 0], :3]-y[pairs[:, 1], :3])@frame
        gap = abs(delta)-[10000., 10000., 0.]
        for p in np.flatnonzero(np.linalg.norm(np.maximum(gap, 0), axis=1) < reserve):
            i, j = map(int, pairs[p]); axis = int(np.argmax(gap[p])); sign = float(np.sign(delta[p, axis]) or 1.)
            if axis == 2:
                # Retain the side on first entry to the near-contact band.
                # Using the service normal would confuse the 90-degree turn
                # with a forbidden crossing of the later feathered plane.
                sign = approach_sign.setdefault((i,j), sign)
            cuts.append((k, i, j, axis, sign))
    return np.asarray(cuts, dtype=float).reshape(-1, 5)


def contact_values(states, frames, cuts):
    if not len(cuts):
        return np.empty(0)
    k, i, j, a = cuts[:, :4].astype(int).T
    direction = frames[k, :, a]*cuts[:, 4, None]
    return np.sum((states[k, i, :3]-states[k, j, :3])*direction, axis=1)-np.where(a < 2, 10000., 0.)


def contact_jacobian(states, frames, state_jac, frame_jac, cuts):
    k, i, j, a = cuts[:, :4].astype(int).T
    direction = frames[k, :, a]*cuts[:, 4, None]
    first = np.einsum('ri,ric->rc', direction, state_jac[k, i, :3]-state_jac[k, j, :3])
    second = np.einsum('ri,ric->rc', states[k, i, :3]-states[k, j, :3],
                       frame_jac[k, :, a, :]*cuts[:, 4, None, None])
    return first+second


def secant_update(jacobian, step, actual_change):
    """Broyden update from an actually propagated direction; no force fallback."""
    norm = float(step@step)
    if norm <= 1e-20:
        return jacobian.copy()
    defect = actual_change-np.einsum('...c,c->...', jacobian, step)
    return jacobian+defect[..., None]*step/norm


def trust_decision(before, predicted, actual, radius, maximum_state_error_m):
    reduction = before-predicted
    ratio = (before-actual)/reduction if reduction > 1e-12 else -1.
    accept = bool(ratio >= .1 and actual < before and maximum_state_error_m <= 150.)
    if not accept or ratio < .25:
        next_radius = radius/2
    elif ratio > .75 and maximum_state_error_m < 30.:
        next_radius = min(1.5, radius*1.5)
    else:
        next_radius = radius
    return dict(accepted=accept, actual_predicted_reduction_ratio=float(ratio),
                next_radius=float(next_radius), maximum_state_prediction_error_m=float(maximum_state_error_m))


def solve_step(x, radius, gaps, gap_jac, energy, power_margin, terminal,
               terminal_limit, margin=400.):
    """Constrained SQP subproblem: 8 controls and one explicit distance slack.

    Energy/power use supplied coupled attitude responses plus exact thrust
    norms. Terminal is a vector of normalized, cycle-linked residuals. Its
    RMS may not exceed the stated bound. The scalar clearance slack is in km;
    it never relaxes power, energy, or the terminal-debt bound.
    """
    dimension = len(x)
    def objective(z):
        residual = terminal(z[:dimension])
        return energy(z[:dimension])[0]/1e14+.2*np.mean(residual**2)+100*z[-1]
    def constraints(z):
        d=z[:dimension]; total, prefix=energy(d); residual=terminal(d)
        return np.r_[(gaps+gap_jac@d-margin)/1000+z[-1],
            power_margin(d), (1e14-total)/1e14, (51.859e12-prefix)/1e14,
            terminal_limit**2-np.mean(residual**2)]
    lower=np.maximum(-radius, -3.-x); upper=np.minimum(radius, 3.-x)
    # Timing controls are restricted to ±0.5 h; no service-time changes.
    lower[6:]=np.maximum(lower[6:], -1.-x[6:]);upper[6:]=np.minimum(upper[6:], 1.-x[6:])
    start=np.r_[np.zeros(dimension),max(0., (margin-float(gaps.min()))/1000)]
    fit=minimize(objective,start,method='SLSQP',bounds=list(zip(lower,upper))+[(0,None)],
        constraints=[{'type':'ineq','fun':constraints}],options={'maxiter':100,'ftol':1e-9})
    residual=constraints(fit.x)
    return fit.x[:dimension], dict(success=bool(fit.success and residual.min()>=-1e-6),
        status=int(fit.status), message=fit.message, iterations=int(fit.nit),
        clearance_restoration_slack_m=float(fit.x[-1]*1000),
        minimum_scaled_constraint=float(residual.min()), objective=float(fit.fun))
