"""Controlled finite 10 km membranes in a moving local orbital formation.

The formation is an explicit set of separate translational states. All panels
share a Sun-facing normal, with four depth levels so adjacent optical footprints
can overlap without coplanar material overlap. Finite-Sun mutual interception
changes each tile's photon force. Electric acceleration is bounded feedback;
attitude torques, flexible dynamics and exhaust transport remain outside this
model. A local covered window is not a complete lunar shield.
"""
from itertools import combinations

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

from shared import constants as K
from .cycling import sail_basis, ray_geometry, visible_sun
from .fleet import square_vertices, sun_points, projected_polygons, beam_margins
from .optical import unit, length, arriving_ray_vectors


def solar_frame(sample):
    sun, _ = arriving_ray_vectors(sample['positions']['sun'], sample['positions']['earth'],
        sample['sun_v'], sample['earth_v'], sample['sun_a'], sample['earth_a'])
    u = unit(sun); _, b, _ = sail_basis(sun)
    return np.column_stack([u, b, np.cross(u, b)])


def lattice(rows, pitch=9750., depth=1000.):
    """Finite checkerboard depth assignment; no panel is area-weighted."""
    iy, ix = np.indices((rows, rows))
    z = ((ix % 2)+2*(iy % 2))*depth
    q = np.column_stack([z.ravel()-1.5*depth,
        (ix.ravel()-(rows-1)/2)*pitch, (iy.ravel()-(rows-1)/2)*pitch])
    offsets = np.array([(x,y) for y in (-1,0,1) for x in (-1,0,1) if x or y])
    neighbour = np.empty((rows*rows,8), dtype=int)
    valid = np.empty_like(neighbour, dtype=bool)
    for k,(dx,dy) in enumerate(offsets):
        nx, ny = ix.ravel()+dx, iy.ravel()+dy
        valid[:,k] = (nx>=0)&(ny>=0)&(nx<rows)&(ny<rows)
        neighbour[:,k] = np.clip(ny,0,rows-1)*rows+np.clip(nx,0,rows-1)
    # Only neighbouring footprints whose lattice offsets fit in one cell can
    # have common interior. The runner verifies the displacement bound needed
    # for this sparse stencil, including uncontrolled failure cases.
    subsets=[]
    for count in (1,2,3):
        subsets += [(ids,(-1)**(count+1)) for ids in combinations(range(8),count)
                    if np.all(np.ptp(offsets[list(ids)],axis=0)<=1)]
    return q, neighbour, valid, subsets


def illumination(state, frame, sample, clear_side, neighbour, valid, subsets, suns=8):
    """Exact axis-aligned square union in each sampled source perspective.

    Shared panel normals make projected rectangles axis aligned. Inclusion /
    exclusion over the sparse near-neighbour stencil resolves their overlap.
    Multiple reflections and realistic spectral/angular transmission are omitted.
    """
    xyz=state[:,:3]@frame
    source=sample['positions']['sun']
    sun,_=arriving_ray_vectors(source,sample['positions']['earth'],sample['sun_v'],
        sample['earth_v'],sample['sun_a'],sample['earth_a'])
    sources=sun_points(sun,suns,.317)@frame
    qi=xyz[:,None,:]; qj=xyz[neighbour]
    upstream=valid & (qj[...,0]>qi[...,0])
    visible=np.zeros(len(state))
    h=clear_side/2
    for s in sources:
        mag=(s[0]-qi[...,0])/(s[0]-qj[...,0])
        centre=s[1:]+(qj[...,1:]-s[1:])*mag[...,None]-qi[...,1:]
        lo=np.maximum(-h,centre-h*mag[...,None])
        hi=np.minimum(h,centre+h*mag[...,None])
        total=np.zeros(len(state))
        for ids,sign in subsets:
            exists=np.all(upstream[:,ids],axis=1)
            widths=np.maximum(0,np.min(hi[:,ids,:],axis=1)-np.max(lo[:,ids,:],axis=1))
            total+=sign*exists*np.prod(widths,axis=1)
        visible+=1-total/clear_side**2
    return np.clip(visible/suns,0,1)


class Formation:
    def __init__(self, env, centre, duration, rows=25, side=10000., pitch=9750., depth=1000.,
                 cap=1e-3, response_s=900., suns=8):
        self.env,self.centre,self.side,self.cap=env,centre,side,cap
        self.offsets,self.neighbour,self.valid,self.subsets=lattice(rows,pitch,depth)
        self.clear_side=side-110.
        self.kp,self.kd=1/response_s**2,2/response_s
        self.suns=suns
        ts=np.linspace(-600,duration+600,int(duration/600)+3)
        self.frames=CubicSpline(ts,np.array([solar_frame(env.at(t)) for t in ts]))

    def target(self,t):
        c=self.centre(t)
        frame=self.frames(t)
        q=c[:3]+self.offsets@frame.T
        v=c[3:]+self.offsets@self.frames(t,1).T
        ca=self.env.gravity(t,c[:3])+self.isolated_sail(t,c[None,:])[0]
        acc=ca+self.offsets@self.frames(t,2).T
        return np.column_stack([q,v]),acc

    def isolated_sail(self,t,state):
        sample=self.env.at(t); rays=ray_geometry(state[:,:3],sample)
        n=-unit(solar_frame(sample)[:,0])
        cosine=np.maximum(np.sum(-unit(rays['sun'])*n,axis=-1),0)
        visible=visible_sun(rays['sun'],rays['earth'],rays['moon'])
        scale=2*self.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*self.env.sigma)
        return (scale*(K.AU/length(rays['sun']))**2*visible*cosine**2*
                (self.clear_side/self.side)**2)[:,None]*n

    def forces(self,t,state):
        target,target_acc=self.target(t)
        frame=self.frames(t); sample=self.env.at(t)
        visible=illumination(state,frame,sample,self.clear_side,self.neighbour,self.valid,
                             self.subsets,self.suns)
        sail=self.isolated_sail(t,state)*visible[:,None]
        natural=self.env.gravity(t,state[:,:3])+sail
        desired=target_acc-natural+self.kp*(target[:,:3]-state[:,:3])+self.kd*(target[:,3:]-state[:,3:])
        control=desired*np.minimum(1,self.cap/np.maximum(length(desired),1e-100))[:,None]
        return natural,control,visible,target

    def derivative(self,t,flat):
        state=flat.reshape(-1,6)
        natural,control,_,_=self.forces(t,state)
        return np.column_stack([state[:,3:],natural+control]).ravel()


def coverage_window(state,normal,sample,side,centre,radius=50000.,suns=32):
    """Finite-Sun shadow union on a moving receiver window in the lunar plane."""
    import shapely
    from shapely.geometry import Point
    sun,_=arriving_ray_vectors(sample['positions']['sun'],sample['positions']['earth'],
        sample['sun_v'],sample['earth_v'],sample['sun_a'],sample['earth_a'])
    u=unit(sun); _,b,c=sail_basis(sun)
    q=state[:,:3]; tau=np.maximum(q@u,0)/K.SPEED_OF_LIGHT
    q=q-(state[:,3:]+sample['moon_v'])*tau[:,None]
    tau_c=max(centre[:3]@u,0)/K.SPEED_OF_LIGHT
    qc=centre[:3]-(centre[3:]+sample['moon_v'])*tau_c
    # The centre's central-ray intercept defines the moving window.
    ratio=(sun@u)/(sun@u-qc@u)
    intercept=sun+(qc-sun)*ratio
    disk=Point(intercept@b,intercept@c).buffer(radius,quad_segs=128)
    vertices=square_vertices(q,np.broadcast_to(normal,q.shape),np.broadcast_to(sun,q.shape),side-110.)
    full=disk; covered=0.
    for source in sun_points(sun,suns,.127):
        polygons=projected_polygons(vertices,source,u,b,c)
        union=shapely.union_all(polygons).intersection(disk)
        covered+=union.area; full=full.intersection(union)
    return dict(ray_coverage=covered/(suns*disk.area),all_sampled_sun_coverage=full.area/disk.area,
                window_centre_distance_m=float(np.hypot(intercept@b,intercept@c)))


def swept_square_conflicts(times,states,frames,side=10000.,clearance=100.,
                           acceleration_bound=.15,normal_rate_bound=1e-6):
    """Swept separating-axis exclusion for parallel finite squares.

    Each interval uses a common midpoint frame. Bounds include curvature of
    both centre histories and the corners' attitude motion. A flagged pair is
    a possible conflict, not necessarily an actual impact.
    """
    pairs_seen=set(); smallest_clearance=np.inf
    for k,dt in enumerate(np.diff(times)):
        q0,q1=states[k,:,:3],states[k+1,:,:3]
        mid=(q0+q1)/2; travel=float(length(q1-q0).max())
        sag=acceleration_bound*dt**2/4
        rotation=side*np.sqrt(2)*normal_rate_bound*dt/2
        pairs=cKDTree(mid).query_pairs(np.sqrt(2)*side+clearance+travel+sag,output_type='ndarray')
        if not len(pairs):continue
        frame=frames((times[k]+times[k+1])/2)
        a=(q0[pairs[:,0]]-q0[pairs[:,1]])@frame
        z=(q1[pairs[:,0]]-q1[pairs[:,1]])@frame
        lower=np.where(a*z<=0,0,np.minimum(np.abs(a),np.abs(z)))
        extent=np.array([0.,side,side])+sag+rotation
        separation=np.max(lower-extent,axis=1)
        smallest_clearance=min(smallest_clearance,float(separation.min()))
        pairs_seen.update(map(tuple,pairs[separation<clearance]))
    return dict(possible_conflict_pairs=len(pairs_seen),
        minimum_proven_axis_clearance_m=float(smallest_clearance),
        acceleration_bound_m_s2=acceleration_bound,normal_rate_bound_rad_s=normal_rate_bound,
        required_clearance_m=clearance)
