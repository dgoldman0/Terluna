"""Endpoint-flexible, individually power-capped linear trajectory proposals.

Shooting matrices are gravity linearizations. Every control is charged from
its actual start. A successful program is a proposal, requiring independent
nonlinear shadow-force replay and actual next-service coverage.
"""
import warnings
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix,coo_matrix,eye,hstack,vstack,block_diag
from scipy.spatial import cKDTree
from .packing_control import matrices
from .optical import length


def fit_controls(base,maps,target,windows,caps,endpoint_times,cuts=(),slack=None,wall_s=15.):
    n=len(target);arcs=len(windows);width=3*arcs;size=n*width
    scale=np.array([1e4]*3+[1.]*3)
    B=matrices(maps,endpoint_times,arcs)[:,:,6:]
    bg=base(endpoint_times)[np.arange(n),np.arange(n)]
    rhs=(target-bg)/scale
    A=block_diag([csr_matrix(b/scale[:,None]) for b in B],format='csr')
    A=hstack([A,csr_matrix((6*n,size))],format='csr')
    ident=eye(size,format='csr')
    upper=[hstack([ident,-ident]),hstack([-ident,-ident])];limits=[np.zeros(size*2)]
    equality=A;eq_rhs=rhs.ravel()
    if slack is not None:
        upper.extend([A,-A]);s=np.broadcast_to(slack,(n,6))/scale
        limits.extend([(rhs+s).ravel(),(-rhs+s).ravel()]);equality=None;eq_rhs=None
    if cuts:
        ts=np.unique([c[0] for c in cuts]);bm=matrices(maps,ts,arcs)[:,:3,6:];ys=base(ts)[:,:,:3]
        lookup={float(t):i for i,t in enumerate(ts)};row=[];col=[];value=[];bound=[]
        for k,(t,i,j,axis,required) in enumerate(cuts):
            ix=lookup[float(t)];b=axis@bm[ix]
            row.extend([k]*(2*width));col.extend(np.r_[i*width+np.arange(width),j*width+np.arange(width)])
            value.extend(np.r_[-b,b]/1e4);bound.append((axis@(ys[ix,i]-ys[ix,j])-required)/1e4)
        upper.append(coo_matrix((value,(row,col)),shape=(len(cuts),2*size)).tocsr());limits.append(np.array(bound))
    # Bounds conservatively inscribe the per-member vector acceleration ball.
    component=np.repeat(np.asarray(caps)*np.diff(windows,axis=1).ravel()[None,:]/(2*np.sqrt(3)),3,axis=1).ravel()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        fit=linprog(np.r_[np.zeros(size),np.ones(size)/n],A_ub=vstack(upper,format='csr'),
            b_ub=np.r_[*limits],A_eq=equality,b_eq=eq_rhs,
            bounds=list(zip(-component,component))+[(0.,None)]*size,method='highs',
            options={'threads':1,'parallel':False,'time_limit':wall_s})
    row=dict(success=bool(fit.success),status=int(fit.status),message=fit.message,cuts=len(cuts))
    if not fit.success:return None,row
    controls=fit.x[:size].reshape(n,arcs,3);dv=length(controls).sum(axis=1)
    row.update(mean_translation_m_s=float(dv.mean()),maximum_translation_m_s=float(dv.max()),
        maximum_acceleration_m_s2=float((length(controls)*2/np.diff(windows,axis=1).ravel()).max()))
    return controls,row


def propose(base,maps,controls,times):
    b=matrices(maps,times,controls.shape[1])[:,:,6:]
    return base(times)+np.einsum('tij,nj->tni',b,controls.reshape(len(controls),-1))


def violations(times,states,cmd,margin=500.,retain=2500,history=None):
    records=[];minimum=np.inf;count=0
    for t,y,f in zip(times,states,cmd(times).as_matrix()):
        pairs=cKDTree(y[:,:3]).query_pairs(10000*np.sqrt(2)+margin,output_type='ndarray')
        if not len(pairs):continue
        delta=(y[pairs[:,0],:3]-y[pairs[:,1],:3])@f
        sep=abs(delta)-[10000,10000,0];distance=length(np.maximum(sep,0))
        minimum=min(minimum,float(distance.min()))
        for k in np.flatnonzero(distance<margin):
            count+=1;i,j=map(int,pairs[k]);axis=int(np.argmax(sep[k]));sign=1 if delta[k,axis]>=0 else -1
            # Retain the sign at the last safe state for normal-axis avoidance;
            # this gives earlier controls a consistent escape branch.
            if axis==2 and history is not None:
                sign=1 if (history[i,:3]-history[j,:3])@f[:,2]>=0 else -1
            direction=f[:,axis]*sign;required=(10000 if axis<2 else 0)+margin
            records.append((float(distance[k]),(float(t),i,j,direction,float(required))))
    records.sort(key=lambda x:x[0])
    return [r[1] for r in records[:retain]],dict(violating_samples=count,minimum_sampled_distance_m=None if not np.isfinite(minimum) else minimum)
