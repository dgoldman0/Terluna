"""Repeated attitude schedule with millisecond knot coalescing.

Optimization coefficients approaching a grid-aligned turn boundary can create
two timestamps separated by 1e-10 s. Coalesce such knots before constructing
the spline; keep the exact smooth angle law at the surviving timestamps.
"""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline
from protection.dynamics.service_search import smooth_step
from protection.dynamics.pattern_departure import sun_frame
from .coupled_cycle import CycleExperiment


def coalesce(times,tolerance=.001):
    ordered=np.sort(np.asarray(times,dtype=float));keep=[ordered[0]]
    for t in ordered[1:]:
        if t-keep[-1]>tolerance:keep.append(t)
    return np.array(keep)


class StableExperiment(CycleExperiment):
    def command(self,x):
        duration=(6.+.5*x[6])*3600.
        knots=[0.,21600.,21600.+duration,self.arrival-28800.,self.arrival,
               self.arrival+21600.,self.arrival+21600.+duration,self.end]
        times=coalesce(np.r_[np.arange(0.,self.end,120.),knots])
        times=times[(times>=0)&(times<=self.end)]
        first=smooth_step((times-21600.)/duration)*(1-smooth_step((times-self.arrival+28800.)/28800.))
        second=smooth_step((times-self.arrival-21600.)/duration)
        angles=-np.pi/2*(first+second)
        frames=np.array([sun_frame(self.env,t) for t in times])@Rotation.from_rotvec(np.c_[angles,np.zeros((len(times),2))]).as_matrix()
        return RotationSpline(times,Rotation.from_matrix(frames))
