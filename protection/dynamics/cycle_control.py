"""Smooth return attitudes and independent endpoint/encounter proposals.

The independent-member force model omits mutual shadows. Gravity response
maps and L1 control costs supply proposals and restricted-model bounds only.
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix, eye, hstack, kron, vstack
from scipy.spatial.transform import Rotation, RotationSpline

from shared import constants as K
from .pattern_departure import sun_frame
from .packing_control import matrices
from .pattern_return import smooth_arc
from .cycling import ray_geometry, visible_sun
from .fleet import sun_points
from .natural_pattern import finite_gravity, retarded_sun
from .optical import unit, length


def return_command(env, start, end, turn_s):
    """Continue the verified feathered frame; make a quintic service turn."""
    if turn_s <= 0 or end-turn_s < start:
        raise ValueError('Arrival turn does not fit')
    times=np.unique(np.r_[np.arange(start,end,120.),end-turn_s,end])
    f=np.clip((times-end+turn_s)/turn_s,0.,1.)
    angle=-np.pi/2*(1-10*f**3+15*f**4-6*f**5)
    frames=np.array([sun_frame(env,t) for t in times])
    frames=frames@Rotation.from_rotvec(np.c_[angle,np.zeros((len(times),2))]).as_matrix()
    return RotationSpline(times,Rotation.from_matrix(frames))


def independent_force(env,t,state,command,mass_ratio=1.):
    """Unshadowed finite-source force for proposal generation, with body fraction."""
    q=state[:,:3];frame=command(t).as_matrix();sample=env.at(t)
    rays=ray_geometry(q,sample)
    visible=visible_sun(rays['sun'],rays['earth'],rays['moon'])
    n=frame[:,2];coefficient=np.zeros(len(q))
    for source in sun_points(retarded_sun(sample),8,.317):
        ray=source-q;cosine=-unit(ray)@n
        coefficient+=cosine*abs(cosine)*(K.AU/length(ray))**2/8
    scale=2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)
    sail=(scale*(9890/10000)**2*coefficient*visible/np.asarray(mass_ratio))[:,None]*n
    return finite_gravity(env,t,q,frame[:,0],frame[:,1])+sail


def independent_path(env,initial,start,end,command,controls=None,windows=None,mass_ratio=1.):
    def rhs(t,flat):
        state=flat.reshape(-1,6);force=independent_force(env,t,state,command,mass_ratio)
        if controls is not None:
            force+=np.einsum('k,nki->ni',smooth_arc(t,windows),controls)
        return np.c_[state[:,3:],force].ravel()
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=300.,
        rtol=2e-11,atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def cycle_windows(start,end):
    middle=(start+end)/2
    windows=np.array([[start,start+7200.],[middle-3600.,middle+3600.],
                      [end-21600.,end-14400.],[end-7200.,end]])
    if np.any(windows[1:,0]<windows[:-1,1]):raise ValueError('Burn windows overlap')
    return windows


def corrected_state(base,maps,controls,times):
    t=np.atleast_1d(times);response=matrices(maps,t,controls.shape[1])[:,:,6:]
    return base(t)+np.einsum('tij,nj->tni',response,controls.reshape(len(controls),-1))


def endpoint_controls(base,maps,target,windows,cuts=(),cap=.001,wall_s=25.):
    """Minimum component L1, explicit Euclidean cost and a relaxation bound."""
    count=len(target);width=3*len(windows);size=count*width
    response=matrices(maps,[maps.t[-1]],len(windows))[0,:,6:]
    scale=np.array([1e4]*3+[1.]*3)
    rhs=target-base([maps.t[-1]])[0]
    equality=hstack([kron(eye(count),csr_matrix(response/scale[:,None])),csr_matrix((count*6,size))],format='csr')
    identity=eye(size,format='csr')
    upper=[hstack([identity,-identity]),hstack([-identity,-identity])]
    limits=[np.zeros(2*size)]
    row=[];col=[];data=[];bounds=[]
    for k,cut in enumerate(cuts):
        t,i,j,axis,required=cut
        b=matrices(maps,[t],len(windows))[0,:3,6:]
        q=base([t])[0,:,:3];a=axis@b
        row.extend([k]*(2*width));col.extend(np.r_[i*width+np.arange(width),j*width+np.arange(width)])
        data.extend(np.r_[-a,a]/1e4);bounds.append((axis@(q[i]-q[j])-required)/1e4)
    if cuts:
        geometry=coo_matrix((data,(row,col)),shape=(len(cuts),2*size)).tocsr()
        upper.append(geometry);limits.append(np.array(bounds))
    component=np.tile(np.repeat(cap*np.diff(windows,axis=1).ravel()/(2*np.sqrt(3)),3),count)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected')
        fit=linprog(np.r_[np.zeros(size),np.ones(size)/count],A_ub=vstack(upper,format='csr'),
            b_ub=np.r_[*limits],A_eq=equality,b_eq=(rhs/scale).ravel(),
            bounds=list(zip(-component,component))+[(0.,None)]*size,method='highs',
            options={'threads':1,'parallel':False,'time_limit':wall_s})
    result=dict(success=bool(fit.success),status=int(fit.status),message=fit.message,
        cuts=len(cuts),control_coefficients=size,acceleration_cap_m_s2=cap)
    if not fit.success:return None,result
    u=fit.x[:size].reshape(count,len(windows),3)
    defect=u.reshape(count,width)@response.T-rhs
    norms=length(u)
    result.update(component_L1_mean_m_s=float(fit.fun),
        restricted_linear_euclidean_lower_bound_m_s=float(fit.fun/np.sqrt(3)),
        mean_delta_v_m_s=float(norms.sum(axis=1).mean()),maximum_delta_v_m_s=float(norms.sum(axis=1).max()),
        maximum_acceleration_m_s2=float((norms*2/np.diff(windows,axis=1).ravel()).max()),
        burn_count=int(np.count_nonzero(norms>1e-8)),
        endpoint_position_residual_m=float(length(defect[:,:3]).max()),
        endpoint_velocity_residual_m_s=float(length(defect[:,3:]).max()))
    return u,result
