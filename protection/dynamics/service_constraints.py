"""Sparse service-ray and encounter constraints, with free terminal states.

Only a proposal model: ray assignments, linearized gravity responses and
selected separating branches do not establish nonlinear service or clearance.
"""
import warnings
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,eye,hstack,vstack
from .packing_control import matrices
from .natural_pattern import retarded_sun
from .cycling import sail_basis
from .optical import unit,length
from shared import constants as K


def projection(q,source,u,b,c):
    z=q@u;depth=source@u;source_xy=np.array([source@b,source@c])
    factor=depth/(depth-z)
    xy=source_xy+factor[:,None]*(np.c_[q@b,q@c]-source_xy)
    jac=factor[:,None,None]*np.stack([b,c])[None,:,:]+(xy-source_xy)[:,:,None]*u[None,None,:]/(depth-z)[:,None,None]
    return xy,jac,factor


def service_rows(env,base,maps,arrival,arcs,margin=250.,dates=7,assignment=None):
    """Assign sampled receiver/Sun rays to nearest natural tile, not a lattice."""
    rows=[];theta=np.arange(12)*2*np.pi/12
    receiver=np.r_[np.zeros((1,2)),10000*np.c_[np.cos(theta),np.sin(theta)]]
    for t in np.linspace(arrival,arrival+21600,dates):
        y=base([t])[0];sample=env.at(t);sun=retarded_sun(sample);u=unit(sun);_,b,c=sail_basis(sun)
        centre=y.mean(axis=0)[:3];window=projection(centre[None,:],sun,u,b,c)[0][0]+receiver
        angles=np.arange(8)*2*np.pi/8+.137
        sources=np.r_[sun[None,:],sun+K.SUN_RADIUS*(np.cos(angles)[:,None]*b+np.sin(angles)[:,None]*c)]
        B=matrices(maps,[t],arcs)[0,:3,6:]
        for source in sources:
            xy,jac,mag=projection(y[:,:3],source,u,b,c)
            half=9890/2*mag-margin
            assigned_xy=xy if assignment is None else projection(assignment([t])[0,:,:3],source,u,b,c)[0]
            assignments=np.argmin(np.max(abs(window[:,None,:]-assigned_xy[None,:,:])/half[None,:,None],axis=2),axis=1)
            for target,i in zip(window,assignments):
                for axis in range(2):
                    gradient=jac[i,axis]@B
                    for sign in [-1,1]:
                        rows.append((int(i),sign*gradient,float(half[i]+sign*(target[axis]-xy[i,axis]))))
    return rows


def fit_free_controls(base,maps,windows,mass,caps,service,cuts=(),total_budget_J=1e14,prefix_budget_J=5.1859e13,specific=30000/(1.4*np.cos(np.pi/4)),wall_s=8.):
    n=len(mass);arcs=len(windows);width=arcs*3;size=n*width
    rr=[];cc=[];vv=[];bound=[]
    for i,gradient,rhs in service:
        k=len(bound);rr.extend([k]*width);cc.extend(i*width+np.arange(width));vv.extend(gradient/1e4);bound.append(rhs/1e4)
    if cuts:
        dates=np.unique([c[0] for c in cuts]);B=matrices(maps,dates,arcs)[:,:3,6:];q=base(dates)[:,:,:3];index={float(t):k for k,t in enumerate(dates)}
        for t,i,j,axis,required in cuts:
            ix=index[float(t)];g=axis@B[ix];k=len(bound)
            rr.extend([k]*(2*width));cc.extend(np.r_[i*width+np.arange(width),j*width+np.arange(width)]);vv.extend(np.r_[-g,g]/1e4)
            bound.append((axis@(q[ix,i]-q[ix,j])-required)/1e4)
    A=coo_matrix((vv,(rr,cc)),shape=(len(bound),size)).tocsr();I=eye(size,format='csr')
    rows=[hstack([A,coo_matrix(A.shape)]),hstack([I,-I]),hstack([-I,-I])]
    rhs=[np.array(bound),np.zeros(2*size)]
    weights=np.repeat(mass*specific,width)
    # Integral of the sin² arc by the fixed twelve-hour comparison boundary.
    from .pattern_return import smooth_arc
    dates=np.unique(np.r_[21600.,np.arange(21600.,43200.,30.),43200.])
    fractions=np.trapezoid(np.array([smooth_arc(t,windows) for t in dates]),dates,axis=0)
    prefix_weights=weights*np.tile(np.repeat(fractions,3),n)
    rows.extend([coo_matrix((np.r_[np.zeros(size),weights]/1e14)[None,:]),coo_matrix((np.r_[np.zeros(size),prefix_weights]/1e14)[None,:])])
    rhs.append(np.array([total_budget_J,prefix_budget_J])/1e14)
    component=np.repeat(np.asarray(caps)*np.diff(windows,axis=1).ravel()[None,:]/(2*np.sqrt(3)),3,axis=1).ravel()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        fit=linprog(np.r_[np.zeros(size),weights/1e14],A_ub=vstack(rows,format='csr'),b_ub=np.r_[*rhs],
            bounds=list(zip(-component,component))+[(0,None)]*size,method='highs',options={'threads':1,'parallel':False,'time_limit':wall_s})
    result=dict(success=bool(fit.success),status=int(fit.status),message=fit.message,service_ray_inequalities=len(service),encounter_cuts=len(cuts),
        total_translation_budget_J=total_budget_J,prefix_translation_budget_J=prefix_budget_J)
    if not fit.success:return None,result
    control=fit.x[:size].reshape(n,arcs,3)
    result.update(mean_translation_m_s=float(length(control).sum(axis=1).mean()),weighted_component_energy_J=float(fit.fun*1e14))
    return control,result
