"""Service-led return primitives, without a prescribed terminal lattice."""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline
from scipy.optimize import minimize_scalar
from .pattern_departure import sun_frame
from .natural_pattern import retarded_sun
from .optical import unit, length
from shared import constants as K


def smooth_step(x):
    x=np.clip(x,0,1)
    return x**3*(10+x*(-15+6*x))


def free_command(env,end,tilt=90.,axis=0,arrival=None,turn_start=7.,turn_hours=4.,return_hours=8.):
    """Same first service, charged smooth departure and optional service turn."""
    knots=[0.,end,21600.,turn_start*3600,(turn_start+turn_hours)*3600]
    if arrival is not None:knots.extend([arrival-return_hours*3600,arrival,arrival+21600])
    times=np.unique(np.r_[np.arange(0,end,120.),knots]);times=times[(times>=0)&(times<=end)]
    angle=-np.deg2rad(tilt)*smooth_step((times-turn_start*3600)/(turn_hours*3600))
    if arrival is not None:angle*=1-smooth_step((times-arrival+return_hours*3600)/(return_hours*3600))
    vector=np.zeros((len(times),3));vector[:,axis]=angle
    frames=np.array([sun_frame(env,t) for t in times])@Rotation.from_rotvec(vector).as_matrix()
    return RotationSpline(times,Rotation.from_matrix(frames))


def natural_opportunities(env,path,start,end,step=1800.):
    """Find Sunward transverse-distance minima; these only propose service dates."""
    def metric(t):
        state=path(t);q=state[:3];u=unit(retarded_sun(env.at(t)))
        return float(length(q-(q@u)*u)) if q@u>0 else 1e12
    ts=np.arange(start,end+step/2,step);r=np.array([metric(t) for t in ts]);out=[]
    for k in range(1,len(ts)-1):
        if r[k]<r[k-1] and r[k]<r[k+1] and r[k]<4*K.MOON_RADIUS:
            fit=minimize_scalar(metric,bounds=(ts[k-1],ts[k+1]),method='bounded',options={'xatol':1.})
            arrival=float(fit.x-10800.)
            if arrival>=start and arrival+21600<=end:
                out.append(dict(arrival_s=arrival,midpass_s=float(fit.x),minimum_transverse_m=float(fit.fun)))
    return out


def cadence_account(arrival_s,service_s=21600.,tiles=361,mass_kg=2263298487.05685):
    """Necessary time-capacity inventory only; spatial coverage/handover unproved."""
    multiple=arrival_s/service_s
    return dict(start_to_start_cadence_s=float(arrival_s),credited_service_s=service_s,
        time_capacity_group_lower_bound=float(multiple),integer_groups_if_undivided=int(np.ceil(multiple)),
        time_capacity_tile_lower_bound=float(tiles*multiple),time_capacity_mass_lower_bound_kg=float(mass_kg*multiple),
        spatial_coverage_and_handover_demonstrated=False)
