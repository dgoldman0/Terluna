"""Finite 10 km tiles, smooth thrust arcs and orbital retiming."""
import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from .fleet import fleet_acceleration,square_vertices
from .optical import unit

GX,GW=np.polynomial.legendre.leggauss(3)
UV=np.array([(x,y) for x in GX for y in GX])
WEIGHTS=np.array([x*y/4 for x in GW for y in GW])


def finite_tile_acceleration(env,t,state,side=10000.,diagnostics=False):
    single=np.asarray(state).ndim==1
    states=np.atleast_2d(state)
    d=fleet_acceleration(env,t,states,4*K.MOON_RADIUS,side,diagnostics=True)
    corners=square_vertices(states[:,:3],d['normal'],d['sun'],side)
    a=(corners[:,1]-corners[:,0])/side; b=(corners[:,3]-corners[:,0])/side
    points=states[:,None,:3]+side/2*(UV[None,:,:1]*a[:,None,:]+UV[None,:,1:]*b[:,None,:])
    gravity=np.einsum('j,ijk->ik',WEIGHTS,env.gravity(t,points))
    sail=d['sail']*((side-110.)/side)**2
    acceleration=gravity+sail
    if diagnostics:
        error=np.linalg.norm(gravity-env.gravity(t,states[:,:3]),axis=-1)
        return (acceleration[0],d,float(error[0])) if single else (acceleration,d,error)
    return acceleration[0] if single else acceleration


def windows(duration,burn,count):
    starts=np.linspace(0,duration-burn,count)
    if count<2 or np.any(np.diff(starts)<burn):raise ValueError('Thrust arcs overlap')
    return [(float(s),float(s+burn)) for s in starts]


def command(t,state,parameters,duration,burn):
    """Sin² profiles, each parameter vector has units of integrated delta-v.

    Components are radial / tangential / orbit-normal in the current state.
    The norm integral per arc is exactly the norm of its parameter vector.
    """
    shape=np.asarray(state).shape[:-1]
    coefficients=np.asarray(parameters).reshape(shape+(-1,3))
    for k,(start,end) in enumerate(windows(duration,burn,coefficients.shape[-2])):
        if start<=t<=end:break
    else:return np.zeros_like(state[...,:3])
    e=unit(state[...,:3]); h=unit(np.cross(state[...,:3],state[...,3:])); tangent=np.cross(h,e)
    v=coefficients[...,k,:]
    return (e*v[...,0,None]+tangent*v[...,1,None]+h*v[...,2,None])*(2/burn)*np.sin(np.pi*(t-start)/burn)**2


def propagate(env,initial,duration,parameters,burn=21600.,rtol=2e-10,max_step=1200.,output_step=None):
    shape=np.asarray(initial).shape
    count=int(np.prod(shape[:-1])) if len(shape)>1 else 1
    def f(t,y):
        state=y.reshape(shape)
        return np.concatenate([state[...,3:],finite_tile_acceleration(env,t,state)+
            command(t,state,parameters,duration,burn)],axis=-1).ravel()
    # Explicit arc boundaries prevent the integrator stepping across an
    # unobserved short thrust interval.
    state=np.array(initial).ravel(); pieces=[]
    boundaries=sorted(set([0.,duration]+[x for arc in windows(duration,burn,np.asarray(parameters).shape[-1]//3) for x in arc]))
    for begin,end in zip(boundaries[:-1],boundaries[1:]):
        if end<=begin:continue
        times=None if output_step is None else np.linspace(begin,end,int(np.ceil((end-begin)/output_step))+1)
        sol=solve_ivp(f,[begin,end],state,method='DOP853',rtol=rtol,
            atol=np.tile([1e-4]*3+[1e-8]*3,count),max_step=max_step,t_eval=times)
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y[:,-1]; pieces.append(sol)
    t=np.concatenate([s.t if i==0 else s.t[1:] for i,s in enumerate(pieces)])
    y=np.concatenate([s.y.T if i==0 else s.y.T[1:] for i,s in enumerate(pieces)]).reshape((-1,)+shape)
    return t,y,sum(s.nfev for s in pieces)


def rotate_state(state,angle):
    axis=unit(np.cross(state[:3],state[3:]))
    def turn(v):return v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
    return np.r_[turn(state[:3]),turn(state[3:])]
