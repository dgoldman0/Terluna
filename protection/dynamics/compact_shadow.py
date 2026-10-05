"""Exact parallel shadows after removing empty receiving-plane rectangles."""
import numpy as np

from shared import constants as K
from .mass_feedback import MassiveParallelPattern
from .parallel_shadow import source_neighbours
from .natural_pattern import rectangle_union,retarded_sun,finite_gravity
from .fleet import sun_points
from .cycling import ray_geometry
from .collection import require_uneclipsed_tiles
from .optical import length,unit


def compact_illumination(q,frame,source,clear_side=9890.):
    neighbour,valid,_=source_neighbours(q,source)
    centre=q.mean(axis=0);basis=frame[:,[2,0,1]]
    xyz=(q-centre)@basis;s=(source-centre)@basis
    qi,qj=xyz[:,None,:],xyz[neighbour];to_source=s[0]-qi[...,0]
    with np.errstate(divide='ignore',invalid='ignore'):
        fraction=(qj[...,0]-qi[...,0])/to_source
        mag=to_source/(s[0]-qj[...,0])
    active=valid&(fraction>0)&(fraction<1);mag=np.where(active,mag,1.)
    shift=s[1:]+(qj[...,1:]-s[1:])*mag[...,None]-qi[...,1:]
    h=clear_side/2
    lo=np.maximum(-h,shift-h*mag[...,None]);hi=np.minimum(h,shift+h*mag[...,None])
    active &= np.all(hi>lo,axis=-1)
    full=np.any(active&np.all(lo==-h,axis=-1)&np.all(hi==h,axis=-1),axis=1)
    active[full]=False
    width=int(active.sum(axis=1).max())
    visible=np.ones(len(q));visible[full]=0.
    if not width:return visible
    index=np.argsort(~active,axis=1,kind='stable')[:,:width]
    low=np.take_along_axis(lo,index[...,None],axis=1)
    high=np.take_along_axis(hi,index[...,None],axis=1)
    selected=np.take_along_axis(active,index,axis=1)
    area=rectangle_union(low,high,selected)
    visible[~full]=np.clip(1-area[~full]/clear_side**2,0.,1.)
    return visible


class CompactMassivePattern(MassiveParallelPattern):
    """Same finite-source force and mass law; algebraically compact shadow unions."""
    def acceleration(self,t,state,diagnostics=False):
        sample=self.env.at(t);frame=self.command(t).as_matrix();q=state[:,:3]
        require_uneclipsed_tiles(ray_geometry(q,sample),self.side)
        sources=sun_points(retarded_sun(sample),self.suns,.317)
        coefficient=np.zeros(len(q));light=[];n=frame[:,2]
        for source in sources:
            visible=compact_illumination(q,frame,source,self.side-110.)
            incident=source-q;cosine=-unit(incident)@n
            coefficient+=visible*cosine*abs(cosine)*(K.AU/length(incident))**2/self.suns
            light.append(visible)
        scale=2*self.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*self.env.sigma)
        sail=(scale*((self.side-110.)/self.side)**2*coefficient/np.broadcast_to(self.mass_ratio,(len(q),)))[:,None]*n
        gravity=finite_gravity(self.env,t,q,frame[:,0],frame[:,1],self.side)
        if diagnostics:return gravity+sail,np.array(light),sail
        return gravity+sail
