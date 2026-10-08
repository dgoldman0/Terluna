"""Signed preparation impulses within the same three reusable pattern modes."""
import numpy as np
from scipy.optimize import minimize, linprog
import warnings

from .departure_control import mode_impulses
from .optical import length


def signed_impulse(modes, start, cuts, burn, cap=.001):
    a=np.array([c['gradient'] for c in cuts])/10000.
    b=np.array([c['rhs'] for c in cuts])/10000.
    bounds=[(-3.,3.),(-4.,4.),(-4.,4.)]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        feasible=linprog(np.zeros(3),A_ub=-a,b_ub=-b,bounds=bounds,method='highs',
            options={'threads':1,'parallel':False,'time_limit':20.})
    if not feasible.success:
        return start,dict(success=False,message=feasible.message,status=int(feasible.status),
            signed_half_space_LP_feasible=False,cuts=len(cuts))
    def objective(x):
        impulse=mode_impulses(modes,x);norm=np.maximum(length(impulse),1e-12)
        return float(norm.mean()),np.einsum('ni,nij->j',impulse/norm[:,None],modes)/len(modes)
    squared=np.sum(modes*modes,axis=1)
    fit=minimize(lambda x:objective(x),feasible.x,jac=True,method='SLSQP',bounds=bounds,
        constraints=[dict(type='ineq',fun=lambda x:a@x-b,jac=lambda x:a),
            dict(type='ineq',fun=lambda x:(cap*burn/2)**2-squared@(x*x),jac=lambda x:-2*squared*x)],
        options={'ftol':1e-10,'maxiter':150})
    peak=float(length(mode_impulses(modes,fit.x)).max()*2/burn)
    residual=float((a@fit.x-b).min()*10000.)
    return fit.x,dict(success=bool(fit.success and residual>=-1e-3 and peak<=cap+1e-10),
        solver_success=bool(fit.success),status=int(fit.status),message=fit.message,
        signed_half_space_LP_feasible=True,cuts=len(cuts),parameters_m_s=fit.x.tolist(),
        mean_delta_v_m_s=objective(fit.x)[0],maximum_acceleration_m_s2=peak,minimum_cut_residual_m=residual)
