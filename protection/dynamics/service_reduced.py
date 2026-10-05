"""Inexpensive centre-gravity/source dynamics for screening only."""
import numpy as np
from scipy.integrate import solve_ivp
from .cycling import ray_geometry,visible_sun
from .optical import unit,length
from .pattern_return import smooth_arc
from shared import constants as K


def reduced_force(env,t,state,frames,mass_ratio):
    q=state[:,:3];f=np.broadcast_to(frames,(len(q),3,3));n=f[:,:,2]
    rays=ray_geometry(q,env.at(t));direction=unit(rays['sun'])
    cosine=-np.sum(direction*n,axis=1)
    light=visible_sun(rays['sun'],rays['earth'],rays['moon'])
    scale=2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)*(9890/10000)**2
    sail=(scale*cosine*abs(cosine)*light*(K.AU/length(rays['sun']))**2/np.asarray(mass_ratio))[:,None]*n
    return env.gravity(t,q)+sail


def reduced_path(env,initial,start,end,command,controls=None,windows=None,mass_ratio=1.):
    def rhs(t,y):
        state=y.reshape(-1,6);a=reduced_force(env,t,state,command(t).as_matrix(),mass_ratio)
        if controls is not None:a+=np.einsum('k,nki->ni',smooth_arc(t,windows),controls)
        return np.c_[state[:,3:],a].ravel()
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=600.,rtol=2e-10,
        atol=np.tile([1e-4]*3+[1e-8]*3,len(initial)),dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol
