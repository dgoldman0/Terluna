"""Exact global surface distance for parallel squares and conditional time guards."""
import numpy as np
from scipy.spatial import cKDTree

from .optical import length


def closest_squares(q,frame,side=10000.):
    """Global minimum, with a sphere-bound candidate search of proven radius."""
    tree=cKDTree(q);_,neighbour=tree.query(q,k=2)
    delta=abs((q-q[neighbour[:,1]])@frame)
    upper=float(length(np.maximum(delta-[side,side,0.],0.)).min())
    # The difference of two square points has length <= sqrt(2)*side.
    # A farther centre pair cannot improve the existing surface-distance upper bound.
    pairs=tree.query_pairs(np.sqrt(2)*side+upper+1e-6,output_type='ndarray')
    projected=(q[pairs[:,0]]-q[pairs[:,1]])@frame
    values=length(np.maximum(abs(projected)-[side,side,0.],0.))
    k=int(np.argmin(values))
    return float(values[k]),pairs[k],projected[k]


def distance_guard(times,states,command,acceleration_bound=.15,rate_bound=np.deg2rad(.1),point_error=0.):
    """All-pair distance guard, conditional on actual centre/rate bounds.

    Subtract common translational motion, then bound distance to the nearer
    endpoint using relative speed, centre curvature and square-corner motion.
    point_error bounds each reconstructed centre at these evaluation dates.
    """
    distances=[];witnesses=[]
    for t,state,frame in zip(times,states,command(times).as_matrix()):
        value,pair,projection=closest_squares(state[:,:3],frame)
        distances.append(value);witnesses.append((float(t),pair.tolist(),projection.tolist()))
    distances=np.array(distances);dt=np.diff(times)
    relative=states[:,:,3:]-states[:,:,3:].mean(axis=1)[:,None,:]
    speed=length(relative).max(axis=1)
    loss=(np.maximum(speed[:-1],speed[1:])+10000./np.sqrt(2)*rate_bound)*dt
    loss+=acceleration_bound*dt**2/2+2*point_error
    lower=np.minimum(distances[:-1],distances[1:])-loss
    k=int(np.argmin(distances));interval=int(np.argmin(lower))
    return dict(minimum_sampled_surface_distance_m=float(distances[k]),
        sampled_witness=dict(time_s=witnesses[k][0],pair=witnesses[k][1],projection_m=witnesses[k][2]),
        minimum_conditional_all_pair_clearance_m=float(lower[interval]),
        limiting_interval_s=times[interval:interval+2].tolist(),
        intervals_below_100m=int(np.count_nonzero(lower<100.)),
        acceleration_bound_m_s2=acceleration_bound,rate_bound_rad_s=rate_bound,
        reconstruction_error_per_centre_m=point_error,maximum_step_s=float(dt.max()))
