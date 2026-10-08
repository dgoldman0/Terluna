"""Compiled exact rectangle unions, independently checked against NumPy unions.

The historical producers remain unchanged. This implementation uses stable
relative coordinates and an all-pair blocker loop, without a neighbour stencil.
"""
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from shared import constants as K
from .natural_pattern import finite_gravity, retarded_sun
from .fleet import sun_points
from .cycling import ray_geometry
from .collection import require_uneclipsed_tiles
from .tile_torque import gravity_torque
from .optical import unit, length

_LIB = None

def library():
    global _LIB
    if _LIB is None:
        source=Path(__file__).with_suffix('.cpp')
        key=hashlib.sha256(source.read_bytes()).hexdigest()[:20]
        target=Path(tempfile.gettempdir())/('terluna_parallel_'+key+'.so')
        if not target.exists():
            subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(source),'-o',str(target)],check=True)
        _LIB=ctypes.CDLL(str(target))
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        _LIB.parallel_moments.argtypes=[ctypes.c_int,array,array,ctypes.c_double,array]
        _LIB.parallel_moments.restype=None
    return _LIB


def moments(q,frame,source,clear_side=9890.):
    centre=np.mean(q,axis=0)
    xyz=np.ascontiguousarray((q-centre)@frame,dtype=np.float64)
    sun=np.ascontiguousarray((source-centre)@frame,dtype=np.float64)
    out=np.zeros((len(q),3))
    library().parallel_moments(len(q),xyz,sun,clear_side,out)
    return out


def load(env,t,state,frame,suns=8,rotation=.317,torque=False):
    q=state[:,:3];sample=env.at(t)
    require_uneclipsed_tiles(ray_geometry(q,sample),10000.)
    first=np.zeros(len(q));unshadowed=first.copy();front=first.copy();back=first.copy()
    force=np.zeros_like(q);moment=np.zeros_like(q)
    for source in sun_points(retarded_sun(sample),suns,rotation):
        area,mx,my=moments(q,frame,source).T
        ray=source-q;cosine=unit(ray)@frame[:,2]
        irradiance=K.SOLAR_CONSTANT*(K.AU/length(ray))**2/suns
        density=irradiance*abs(cosine);unique=density*area
        first+=unique;unshadowed+=density*9890.**2
        front+=np.where(cosine>=0,unique,0);back+=np.where(cosine<0,unique,0)
        pressure=-2*env.central_fraction/K.SPEED_OF_LIGHT*density*cosine
        force+=(pressure*area)[:,None]*frame[:,2]
        moment+=pressure[:,None]*np.c_[my,-mx,np.zeros(len(q))]
    result=dict(reflected_force_N=force,radiation_torque_N_m=moment,
        optical=dict(first_intercept_bolometric_equivalent_W=first,
            unshadowed_aperture_W=unshadowed,front_first_intercept_W=front,
            back_first_intercept_W=back,mutual_shadow_equivalent_loss_W=unshadowed-first,
            redirected_band_optical_W=first*env.central_fraction))
    if torque:result['gravity_torque_per_mass']=gravity_torque(env,t,q,frame)
    return result


class FastParallelPattern:
    def __init__(self,env,command,mass_ratio,suns=8):
        self.env,self.command,self.mass_ratio,self.suns=env,command,np.asarray(mass_ratio),suns
    def acceleration(self,t,state):
        frame=self.command(t).as_matrix()
        light=load(self.env,t,state,frame,self.suns)
        return finite_gravity(self.env,t,state[:,:3],frame[:,0],frame[:,1])+light['reflected_force_N']/(self.env.sigma*10000.**2*self.mass_ratio[:,None])
