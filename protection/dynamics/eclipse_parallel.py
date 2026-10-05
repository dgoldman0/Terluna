"""Joint finite-Sun body/parallel-array visibility, with area refinement.

No product of separately integrated body and tile visibility fractions is
used. At partial body limbs the same area rays test the union of both bodies
and every foreground clear aperture. Force and energy use this same union.
"""
import ctypes,hashlib,subprocess,tempfile
from pathlib import Path
import numpy as np
from shared import constants as K
from .fast_parallel import load as uneclipsed_load
from .natural_pattern import finite_gravity,retarded_sun
from .fleet import sun_points
from .cycling import ray_geometry
from .collection import require_uneclipsed_tiles
from .tile_torque import gravity_torque
from .optical import unit,length
_LIB=None

def library():
    global _LIB
    if _LIB is None:
        src=Path(__file__).with_suffix('.cpp');other=src.with_name('fast_parallel.cpp')
        key=hashlib.sha256(src.read_bytes()+other.read_bytes()).hexdigest()[:20]
        target=Path(tempfile.gettempdir())/('terluna_eclipse_'+key+'.so')
        if not target.exists():subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(src),'-o',str(target)],check=True)
        _LIB=ctypes.CDLL(str(target));a=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        _LIB.masked_moments.argtypes=[ctypes.c_int,a,a,a,a,ctypes.c_double,ctypes.c_int,a]
        _LIB.masked_moments.restype=None
    return _LIB


def moments(q,frame,source,bodies,radii,grid=16,side=9890.):
    centre=q.mean(axis=0);xyz=np.ascontiguousarray((q-centre)@frame)
    sun=np.ascontiguousarray((source-centre)@frame)
    body=np.ascontiguousarray(np.asarray(bodies)@frame);radii=np.ascontiguousarray(radii,dtype=float)
    out=np.zeros((len(q),3));library().masked_moments(len(q),xyz,sun,body,radii,side,grid,out)
    return out


def load(env,t,state,frame,suns=8,rotation=.317,torque=False,grid=16):
    q=state[:,:3];sample=env.at(t);rays=ray_geometry(q,sample)
    try:
        require_uneclipsed_tiles(rays)
        return uneclipsed_load(env,t,state,frame,suns,rotation,torque)
    except ValueError as exc:
        if 'Body eclipse' not in str(exc):raise
    bodies=np.stack([rays['moon'],rays['earth']],axis=1)
    first=np.zeros(len(q));bare=first.copy();front=first.copy();back=first.copy()
    force=np.zeros_like(q);moment=np.zeros_like(q)
    for source in sun_points(retarded_sun(sample),suns,rotation):
        area,mx,my=moments(q,frame,source,bodies,[K.MOON_RADIUS,K.EARTH_RADIUS],grid).T
        ray=source-q;cosine=unit(ray)@frame[:,2]
        density=K.SOLAR_CONSTANT*(K.AU/length(ray))**2/suns*abs(cosine)
        unique=density*area;first+=unique;bare+=density*9890**2
        front+=np.where(cosine>=0,unique,0);back+=np.where(cosine<0,unique,0)
        pressure=-2*env.central_fraction/K.SPEED_OF_LIGHT*density*cosine
        force+=(pressure*area)[:,None]*frame[:,2]
        moment+=pressure[:,None]*np.c_[my,-mx,np.zeros(len(q))]
    out=dict(reflected_force_N=force,radiation_torque_N_m=moment,
        optical=dict(first_intercept_bolometric_equivalent_W=first,unshadowed_aperture_W=bare,
            body_and_mutual_shadow_equivalent_loss_W=bare-first,front_first_intercept_W=front,
            back_first_intercept_W=back,redirected_band_optical_W=first*env.central_fraction))
    if torque:out['gravity_torque_per_mass']=gravity_torque(env,t,q,frame)
    return out


class EclipseParallelPattern:
    def __init__(self,env,command,mass_ratio,suns=8,grid=16):
        self.env,self.command,self.ratio,self.suns,self.grid=env,command,np.asarray(mass_ratio),suns,grid
    def acceleration(self,t,state):
        f=self.command(t).as_matrix();a=load(self.env,t,state,f,self.suns,grid=self.grid)
        return finite_gravity(self.env,t,state[:,:3],f[:,0],f[:,1])+a['reflected_force_N']/(self.env.sigma*1e8*self.ratio[:,None])
