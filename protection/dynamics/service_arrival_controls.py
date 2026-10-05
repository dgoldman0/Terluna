"""Free-state service LP with early AND arrival arcs and exact prefix charging."""
import warnings
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,eye,hstack,vstack
from .packing_control import matrices
from .optical import length


def arc_fraction(time,windows):
    x=np.clip((time-windows[:,0])/np.diff(windows,axis=1).ravel(),0,1)
    return x-np.sin(2*np.pi*x)/(2*np.pi)


def fit_service_arcs(base,maps,windows,mass,caps,service,cuts,total_budget_J,prefix_budget_J,wall_s=8.):
    n=len(mass);arcs=len(windows);width=arcs*3;size=n*width
    rr=[];cc=[];vv=[];bounds=[]
    for i,g,rhs in service:
        k=len(bounds);rr.extend([k]*width);cc.extend(i*width+np.arange(width));vv.extend(g/1e4);bounds.append(rhs/1e4)
    if cuts:
        dates=np.unique([c[0] for c in cuts]);B=matrices(maps,dates,arcs)[:,:3,6:];q=base(dates)[:,:,:3];index={float(t):k for k,t in enumerate(dates)}
        for t,i,j,axis,required in cuts:
            ix=index[float(t)];g=axis@B[ix];k=len(bounds)
            rr.extend([k]*(2*width));cc.extend(np.r_[i*width+np.arange(width),j*width+np.arange(width)]);vv.extend(np.r_[-g,g]/1e4)
            bounds.append((axis@(q[ix,i]-q[ix,j])-required)/1e4)
    A=coo_matrix((vv,(rr,cc)),shape=(len(bounds),size)).tocsr();I=eye(size,format='csr')
    rows=[hstack([A,coo_matrix(A.shape)]),hstack([I,-I]),hstack([-I,-I])];rhs=[np.array(bounds),np.zeros(2*size)]
    weights=np.repeat(mass*30000/(1.4*np.cos(np.pi/4)),width)
    prefix=weights*np.tile(np.repeat(arc_fraction(43200.,windows),3),n)
    rows.extend([coo_matrix((np.r_[np.zeros(size),weights]/1e14)[None,:]),coo_matrix((np.r_[np.zeros(size),prefix]/1e14)[None,:])])
    rhs.append(np.array([total_budget_J,prefix_budget_J])/1e14)
    component=np.repeat(np.asarray(caps)*np.diff(windows,axis=1).ravel()[None,:]/(2*np.sqrt(3)),3,axis=1).ravel()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        fit=linprog(np.r_[np.zeros(size),weights/1e14],A_ub=vstack(rows,format='csr'),b_ub=np.r_[*rhs],
            bounds=list(zip(-component,component))+[(0,None)]*size,method='highs',options={'threads':1,'parallel':False,'time_limit':wall_s})
    result=dict(success=bool(fit.success),status=int(fit.status),message=fit.message,service_ray_inequalities=len(service),encounter_cuts=len(cuts),total_translation_budget_J=total_budget_J,prefix_translation_budget_J=prefix_budget_J)
    if not fit.success:return None,result
    controls=fit.x[:size].reshape(n,arcs,3);cost=length(controls)*mass[:,None]*30000/(1.4*np.cos(np.pi/4))
    result.update(weighted_component_energy_J=float(fit.fun*1e14),translation_euclidean_J=float(cost.sum()),prefix_translation_euclidean_J=float((cost*arc_fraction(43200.,windows)).sum()),mean_translation_m_s=float(length(controls).sum(axis=1).mean()))
    return controls,result
